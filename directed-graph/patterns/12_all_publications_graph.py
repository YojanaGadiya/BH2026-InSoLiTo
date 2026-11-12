"""Step 12: clean the all-publications directed graph into one node per tool.

The all-publications edge list (citation-order run on every Europe PMC OA article, not
only proteomics and metabolomics) has one node per bio.tools entry. Entries that list the
same paper all get its citations, so Galaxy and EBI wrappers and sub-entries show up as
copies of the tool they wrap (samtools_sort, samtools_flagstat, sam_to_bam ~ SAMtools; a
dozen EBI and Ensembl BLAST entries ~ BLAST).

Each paper's owner comes from step 1 (paper_owner.parquet: title names, only entry,
name in title, main publication; shared platform papers have no owner). Then:
  - an entry that lists papers but owns none is dropped (its edges are other tools');
  - entries owning the same papers are merged, and an entry whose papers are a subset
    of a bigger entry's is merged into it (bigger = more article links in this graph);
  - entries not in paper_owner (papers never cited in the omics corpus) stay as they are.
Edges between merged entries count the same articles, so a tool pair takes the largest
count of its entry pairs, not the sum.

Entries that own their own paper but also list another tool's get that tool's citations too
(bio.tools' quasr and splicemap entries list the HISAT and Bowtie papers). Where that makes
two nodes near-identical copies (same neighbours, same weight) they are merged and named
"a/b"; `citations_mostly_of` gives the owner of their most-cited shared paper. Smaller
leaks of this kind can't be fixed without the article-level data.

Input: directed-graph/citation-order/data/all_publications/directed_edges.parquet
(or .csv). Writes patterns/data/allpub_tools.parquet (tool, entries, strength, citations_mostly_of) and
allpub_edges.parquet (one row per tool pair, in its main direction: source, target,
n_fwd, n_rev, n_total, share_forward, median_delta). Needs step 1.
Run from BH2026-InSoLiTo/:
    python directed-graph/patterns/12_all_publications_graph.py
"""

import re
from collections import defaultdict
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import scipy.sparse as sp

SRC = 'directed-graph/citation-order/data/all_publications/directed_edges'
DATA = 'directed-graph/patterns/data'
TOOL_WORKS = 'directed-graph/citation-order/data/tool_works.parquet'
AFFIXES = {'ebi', 'rest', 'ws', 'api', 'toolshed', 'wrapper', 'galaxy', 'hgv', 'bioconductor',
           'r', 'py', 'web', 'server', 'online', 'biolib', 'app', 'cli', 'plus', 'tool', 'tools'}

con = duckdb.connect()
con.sql("SET preserve_insertion_order = false")
src = f"'{SRC}.parquet'" if Path(f'{SRC}.parquet').exists() else f"read_csv('{SRC}.csv')"
con.sql(f"CREATE TABLE raw AS SELECT source, target, n_articles, n_articles_rev, median_delta FROM {src}")
strength = con.sql("""SELECT t, sum(n) AS n FROM (SELECT source t, n_articles n FROM raw UNION ALL
                      SELECT target, n_articles FROM raw) GROUP BY 1""").df().set_index('t').n.to_dict()

po = con.sql(f"SELECT work_id, owners, listed_by FROM '{DATA}/paper_owner.parquet'").df()
owned, listing = defaultdict(set), set()
for r in po.itertuples(index=False):
    listing.update(r.listed_by)
    for o in r.owners:
        owned[o].add(r.work_id)

# Merge identical paper sets, then subsets into the strongest superset.
by_set = defaultdict(list)
for e, ws in owned.items():
    by_set[frozenset(ws)].append(e)
groups = {}
for ws, members in by_set.items():
    rep = max(members, key=lambda e: (strength.get(e, 0), -len(e)))
    groups[rep] = (set(ws), set(members))
work_groups = defaultdict(set)
for rep, (ws, _) in groups.items():
    for w in ws:
        work_groups[w].add(rep)
for rep in sorted(groups, key=lambda r: strength.get(r, 0)):
    ws, members = groups[rep]
    supers = set.intersection(*(work_groups[w] for w in ws)) - {rep}
    supers = {s for s in supers if s in groups and groups[s][0] > ws}
    if supers:
        target = max(supers, key=lambda s: (strength.get(s, 0), -len(s)))
        groups[target][1].update(members)
        for w in ws:
            work_groups[w].discard(rep)
        del groups[rep]


def display_name(members):
    def affix_count(e):
        return sum(t in AFFIXES for t in re.split(r'[^a-z0-9]+', e.lower()))
    return min(members, key=lambda e: (affix_count(e), -strength.get(e, 0), len(e), e))


entry_tool = {}
for rep, (_, members) in groups.items():
    name = display_name(members)
    for m in members:
        entry_tool[m] = name
nodes = set(strength)
dropped = {e for e in nodes if e in listing and e not in entry_tool}
for e in nodes - set(entry_tool) - dropped:
    entry_tool[e] = e  # never cited in the omics corpus: keep as is



def build_pairs():
    m = pd.DataFrame([{'entry': e, 'tool': t} for e, t in entry_tool.items() if e in nodes])
    con.register('m_df', m)
    con.sql("CREATE OR REPLACE TABLE m AS SELECT * FROM m_df")  # a real table: joins on a registered DataFrame are slow
    con.sql("""
        CREATE OR REPLACE TABLE pairs AS
        SELECT a.tool AS s, b.tool AS t, max(r.n_articles) AS n, max(r.n_articles_rev) AS n_rev,
               arg_max(r.median_delta, r.n_articles) AS median_delta
        FROM raw r JOIN m a ON a.entry = r.source JOIN m b ON b.entry = r.target
        GROUP BY 1, 2
    """)
    # Self-pairs (two entries of one tool) dropped afterwards: as a WHERE or HAVING filter, DuckDB
    # turns a.tool <> b.tool into a join condition and runs out of memory.
    con.sql("DELETE FROM pairs WHERE s = t")
    return m


m = build_pairs()

# Copies the paper rules can't see: entries that list another tool's paper as well as their
# own get that tool's citations (bio.tools' quasr and splicemap entries both list the HISAT
# and Bowtie papers, boltz and chai-1 list ColabFold). They show up as nodes with the same
# neighbours and the same weight. Merge nodes with near-identical link profiles
# (cosine >= COPY_COS, weights within COPY_RATIO of each other).
COPY_COS, COPY_RATIO, COPY_MIN = 0.995, 0.9, 500
und = con.sql("SELECT s, t, n + n_rev AS w FROM pairs").df()
wsum = und.groupby('s').w.sum()
big = wsum[wsum >= COPY_MIN].index
ix = {t: i for i, t in enumerate(big)}
u2 = und[und.s.isin(ix) & und.t.isin(ix)]
M = sp.csr_matrix((u2.w.astype(float), (u2.s.map(ix), u2.t.map(ix))), shape=(len(ix), len(ix)))
M = M + M.T
norm = np.sqrt(M.multiply(M).sum(1)).A1 + 1e-9
Mn = sp.diags(1 / norm) @ M
S = (Mn @ Mn.T).tocoo()
parent = {}
def find(x):
    while parent.get(x, x) != x:
        x = parent[x]
    return x
copies = []
for a, b, c in zip(S.row, S.col, S.data):
    if a < b and c >= COPY_COS:
        ta, tb = big[a], big[b]
        if min(wsum[ta], wsum[tb]) / max(wsum[ta], wsum[tb]) >= COPY_RATIO:
            copies.append((ta, tb, round(float(c), 4)))
            ra, rb = find(ta), find(tb)
            if ra != rb:
                parent[rb] = ra
copy_groups = defaultdict(set)
for x in list(parent) + [p for p in parent.values()]:
    copy_groups[find(x)].add(x)
listed_works = con.sql(f"SELECT biotoolsID, list(work_id) AS w FROM '{TOOL_WORKS}' GROUP BY 1").df().set_index('biotoolsID').w.to_dict()
cites = con.sql(f"SELECT work_id, max(n_citations) AS n FROM '{TOOL_WORKS}' GROUP BY 1").df().set_index('work_id').n.to_dict()
owners_of = dict(zip(po.work_id, po.owners))
copy_rows = []
for members in copy_groups.values():
    entries_ = [e for e, t in entry_tool.items() if t in members]
    shared = set.intersection(*(set(listed_works.get(e, [])) for e in entries_)) if entries_ else set()
    best = max(shared, key=lambda w: cites.get(w, 0), default=None)
    owner = ', '.join(owners_of.get(best, [])) if best is not None else ''
    label = '/'.join(sorted(members))
    for e in entries_:
        entry_tool[e] = label
    copy_rows.append({'tool': label, 'members': sorted(members), 'citations_mostly_of': owner})
if copy_groups:
    m = build_pairs()

# One row per unordered pair in its main direction. n(s->t) is taken from the (s, t) row,
# falling back to the reverse count of the (t, s) row when only that one survived.
con.sql("""
    CREATE TABLE edges AS
    WITH f AS (
        SELECT p.s, p.t, greatest(p.n, coalesce(q.n_rev, 0)) AS n_fwd,
               greatest(p.n_rev, coalesce(q.n, 0)) AS n_rev, p.median_delta
        FROM pairs p LEFT JOIN pairs q ON q.s = p.t AND q.t = p.s
    )
    SELECT s AS source, t AS target, n_fwd, n_rev, n_fwd + n_rev AS n_total,
           n_fwd / (n_fwd + n_rev) AS share_forward, median_delta
    FROM f WHERE n_fwd > n_rev OR (n_fwd = n_rev AND s < t)
""")
tools = m.groupby('tool').entry.apply(lambda x: sorted(x)).reset_index(name='entries')
st = con.sql("""SELECT tool, sum(n) AS strength FROM (SELECT source AS tool, n_total AS n FROM edges UNION ALL
                SELECT target, n_total FROM edges) GROUP BY 1""").df()
tools = tools.merge(st, on='tool', how='inner').sort_values('strength', ascending=False)
cp = pd.DataFrame(copy_rows, columns=['tool', 'members', 'citations_mostly_of'])
tools = tools.merge(cp[['tool', 'citations_mostly_of']], on='tool', how='left')
con.register('tools_df', tools)
con.sql(f"COPY tools_df TO '{DATA}/allpub_tools.parquet' (FORMAT parquet)")
con.sql(f"COPY (SELECT * FROM edges ORDER BY n_total DESC) TO '{DATA}/allpub_edges.parquet' (FORMAT parquet)")

raw_n = con.sql("SELECT count(*) FROM raw").fetchone()[0]
e = con.sql("SELECT count(*), sum(n_total) FROM edges").fetchone()
print(f'Entries: {len(nodes)}; dropped (list papers, own none): {len(dropped)}; '
      f'kept as is (not in the omics corpus): {sum(1 for x in nodes if entry_tool.get(x) == x and x not in owned)}')
print(f'Tools after merging: {len(tools)}; merged groups with 2+ entries: {(tools.entries.str.len() > 1).sum()}')
print(f'Directed rows in: {raw_n:,}; tool pairs out: {e[0]:,} ({e[1]:,} article links)')
print('\nBiggest merges:')
print(tools.assign(k=tools.entries.str.len()).nlargest(8, 'k')[['tool', 'k', 'strength']].to_string(index=False))
print(f'\nNear-identical copies merged ({len(copy_groups)} groups):')
print(cp.merge(tools[['tool', 'strength']], on='tool').sort_values('strength', ascending=False).head(15).to_string(index=False))
print('\nStrongest dropped entries (wrappers / platform listings):')
print(sorted(((strength[d], d) for d in dropped), reverse=True)[:12])
print('\nTop pairs:')
print(con.sql("SELECT source, target, n_fwd, n_rev, round(share_forward, 2) AS fwd_share FROM edges ORDER BY n_total DESC LIMIT 12").df().to_string(index=False))
