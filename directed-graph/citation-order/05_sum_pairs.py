"""Step 5: sum the directed pairs of step 4 by pair and direction.

One row per direction (A -> B and B -> A are separate rows): in how many
articles A is cited before B, and the sum and median of both distances of
step 4: 1/delta and 1/tool_gap (closer citations weigh more). Each row also
carries the numbers of the opposite direction and
`share_forward`, the fraction of the pair's articles that go this way
(0.5 = no direction). Writes everything to a CSV, no threshold on the
number of articles applied.

MAX_DELTA keeps only the article pairs with delta <= MAX_DELTA before summing,
in both directions (None = all). With a limit the output gets a suffix, e.g.
directed_edges_delta1.csv, so the unfiltered file is not overwritten.

Needs the output of 04_directed_pairs.py. Run from BH2026-InSoLiTo/:
    .venv/bin/python directed-graph/citation-order/05_sum_pairs.py
    .venv/bin/python directed-graph/citation-order/05_sum_pairs.py --max-delta 1
    .venv/bin/python directed-graph/citation-order/05_sum_pairs.py --max-delta 1 --dir directed-graph/edges-processed \
        --raw directed-graph/edges-processed/directed_edges_raw-all_citations.parquet
"""

import argparse

import duckdb

OUT_DIR = 'directed-graph/citation-order/data'
MAX_DELTA = None  # e.g. 1 = only tools cited right after each other

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--max-delta', type=int, default=MAX_DELTA, help=f'keep pairs with delta <= this (default: {MAX_DELTA})')
parser.add_argument('--dir', default=OUT_DIR, help=f'folder with directed_edges_raw.parquet, also the output folder (default: {OUT_DIR})')
parser.add_argument('--raw', help='raw pairs to read instead of directed_edges_raw.parquet in --dir')
args = parser.parse_args()

RAW = args.raw or f'{args.dir}/directed_edges_raw.parquet'
OUT_CSV = f'{args.dir}/directed_edges{"" if args.max_delta is None else f"_delta{args.max_delta}"}.csv'
DELTA_FILTER = 'true' if args.max_delta is None else f'delta <= {args.max_delta}'
print(f'Reading {RAW}, {"all pairs" if args.max_delta is None else f"only pairs with delta <= {args.max_delta}"}')

con = duckdb.connect()

con.sql(f"""
    CREATE TEMP TABLE edges AS
    WITH summed AS (
        SELECT source, target,
               count(*) AS n_articles,
               sum(1.0 / delta) AS sum_inv_delta,
               median(delta) AS median_delta,
               sum(1.0 / tool_gap) AS sum_inv_tool_gap,
               median(tool_gap) AS median_tool_gap
        FROM '{RAW}'
        WHERE {DELTA_FILTER}
        GROUP BY source, target
    )
    SELECT f.source, f.target, f.n_articles,
           f.sum_inv_delta, f.median_delta,
           f.sum_inv_tool_gap, f.median_tool_gap,
           coalesce(r.n_articles, 0) AS n_articles_rev,
           coalesce(r.sum_inv_delta, 0) AS sum_inv_delta_rev,
           coalesce(r.sum_inv_tool_gap, 0) AS sum_inv_tool_gap_rev,
           f.n_articles / (f.n_articles + coalesce(r.n_articles, 0)) AS share_forward
    FROM summed f
    LEFT JOIN summed r ON r.source = f.target AND r.target = f.source
""")

print('Directed edges by number of articles:')
con.sql("""
    SELECT count(*) AS edges,
           count(*) FILTER (WHERE n_articles >= 2) AS in_2plus,
           count(*) FILTER (WHERE n_articles >= 5) AS in_5plus,
           count(*) FILTER (WHERE n_articles >= 10) AS in_10plus,
           count(*) FILTER (WHERE n_articles_rev = 0) AS one_way_only
    FROM edges
""").show()

# One row per pair, in its most frequent direction (ties: alphabetical).
print('Most frequent pairs, in their main direction:')
con.sql("""
    SELECT source, target, n_articles, n_articles_rev,
           round(share_forward, 2) AS share_forward,
           round(sum_inv_delta, 1) AS sum_inv_delta,
           round(sum_inv_delta_rev, 1) AS sum_inv_delta_rev,
           median_delta,
           round(sum_inv_tool_gap, 1) AS sum_inv_tool_gap,
           round(sum_inv_tool_gap_rev, 1) AS sum_inv_tool_gap_rev,
           median_tool_gap
    FROM edges
    WHERE n_articles > n_articles_rev
       OR (n_articles = n_articles_rev AND source < target)
    ORDER BY n_articles + n_articles_rev DESC
    LIMIT 20
""").show(max_width=250)

print('Example: edges from maxquant, most frequent first:')
con.sql("""
    SELECT source, target, n_articles, n_articles_rev,
           round(share_forward, 2) AS share_forward, median_delta, median_tool_gap
    FROM edges
    WHERE source = 'maxquant'
    ORDER BY n_articles DESC
    LIMIT 10
""").show()

con.sql(f"""
    COPY (SELECT * FROM edges ORDER BY n_articles DESC, source, target)
    TO '{OUT_CSV}' (HEADER)
""")
print(f'Written {OUT_CSV}')
