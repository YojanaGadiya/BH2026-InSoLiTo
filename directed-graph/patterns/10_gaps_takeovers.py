"""Step 10: gaps in the tool graph, and which tools took over from others.

Uses the EDAM annotations of step 8 (bio.tools operations, input/output data and formats).

1. Compatible but rarely combined: tool A outputs a data type or format that tool B
   takes as input (generic ones such as TSV, HTML or "Data" left out), both are used in
   the same kinds of papers (field mix), yet A then B appears in far fewer articles than
   expected if their use were independent. These are links bio.tools says should work
   that the literature hardly uses.
2. Adapter steps: tools that mostly sit between two others in a methods section (A, then
   the tool, then B), counted over ordered triples in each article. A tool that is mostly
   in the middle, with an EDAM operation such as format conversion or data handling,
   marks a step that exists to move data between tools.
3. Takeovers: as in step 7 (an older tool falling, a newer one rising, the same specific
   neighbours, rarely used together), but now the two must share a specific EDAM
   operation, so only tools doing the same job count.

Writes patterns/data/gaps_compatible.parquet, adapters.parquet, takeovers_edam.parquet.
Needs steps 1, 2, 7 and 8. Run from BH2026-InSoLiTo/:
    python directed-graph/patterns/10_gaps_takeovers.py
"""

from collections import Counter, defaultdict

import duckdb
import numpy as np
import pandas as pd

DATA = 'directed-graph/patterns/data'
EXPORT = 'directed-graph/europepmc/data/export'
KEPT = 'directed-graph/citation-order/data/kept_articles.parquet'

GENERIC_FORMATS = {'html', 'textual format', 'tsv', 'csv', 'xml', 'json', 'image format', 'png', 'pdf', 'txt',
                   'plain text format (unformatted)', 'xlsx', 'jpeg', 'jpg', 'svg', 'binary format', 'zip', 'gzip',
                   'format (by type of data)', 'tabular format', 'yaml', 'r data format', 'rdata', 'xls', 'markdown'}
GENERIC_DATA = {'data', 'report', 'plot', 'image', 'text', 'textual format', 'file', 'identifier', 'data reference',
                'url', 'database entry', 'sequence record', 'annotation track', 'free text'}
GENERIC_OPS = {'visualisation', 'data retrieval', 'query and retrieval', 'database search', 'statistical calculation',
               'data handling', 'deposition', 'annotation', 'analysis', 'plotting', 'formatting', 'editing',
               'standardisation and normalisation', 'data loading', 'filtering', 'sorting',
               # catch-all terms that put unrelated tools together (XDS ~ PISA, ESPript ~ anything)
               'structure visualisation', 'sequence analysis', 'expression analysis', 'phylogenetic tree analysis',
               'protein secondary structure comparison', 'protein secondary structure prediction',
               'sequence comparison', 'protein structure analysis'}
ADAPTER_OPS = {'format conversion', 'formatting', 'data handling', 'file reformatting', 'conversion', 'data conversion',
               'sequence file editing', 'read pre-processing', 'sequence trimming', 'parsing', 'indexing'}

con = duckdb.connect()
edam = con.sql(f"SELECT * FROM '{DATA}/edam_functions.parquet'").df()
edam['t'] = edam.term.str.lower()
tools = con.sql(f"SELECT tool, n_articles FROM '{DATA}/tools.parquet'").df().set_index('tool').n_articles
trends = con.sql(f"SELECT * FROM '{DATA}/tool_trends.parquet'").df().set_index('tool')
pairs = con.sql(f"SELECT source, target, sum(n) AS n FROM '{DATA}/pairs_by_year.parquet' GROUP BY 1, 2").df()
n_articles = con.sql(f"SELECT count(DISTINCT pmcid) FROM '{DATA}/tool_use.parquet'").fetchone()[0]
co = defaultdict(int)
for s, t, n in pairs.itertuples(index=False):
    co[(s, t)] += n
pair_any = defaultdict(int)
for (s, t), n in co.items():
    pair_any[tuple(sorted((s, t)))] += n

# ---------------------------------------------------------------- 1. compatible, rarely combined
outs = edam[edam.kind.isin(['output_format', 'output_data'])]
ins = edam[edam.kind.isin(['input_format', 'input_data'])]
outs = outs[~outs.t.isin(GENERIC_FORMATS | GENERIC_DATA)]
ins = ins[~ins.t.isin(GENERIC_FORMATS | GENERIC_DATA)]
out_by = outs.groupby('tool').t.apply(set)
in_by = ins.groupby('tool').t.apply(set)
lean = trends.field_lean
rows = []
for a, oa in out_by.items():
    for b, ib in in_by.items():
        if a == b:
            continue
        shared = oa & ib
        if not shared or abs(lean.get(a, 0) - lean.get(b, 0)) > 3:  # same kinds of papers
            continue
        expected = tools.get(a, 0) * tools.get(b, 0) / n_articles
        observed = co.get((a, b), 0)
        together = pair_any.get(tuple(sorted((a, b))), 0)
        rows.append({'from_tool': a, 'to_tool': b, 'via': ', '.join(sorted(shared)), 'a_to_b': observed,
                     'together_any_order': together, 'expected_if_independent': round(expected, 1),
                     'ratio': round((together + 0.5) / (expected + 0.5), 2),
                     'from_uses': int(tools.get(a, 0)), 'to_uses': int(tools.get(b, 0))})
compat = pd.DataFrame(rows)
gaps = compat[(compat.expected_if_independent >= 5) & (compat.ratio < 0.5)].sort_values('expected_if_independent', ascending=False)
working = compat[compat.ratio >= 2].sort_values('together_any_order', ascending=False)

# ---------------------------------------------------------------- 2. adapter steps
ranked = con.sql(f"""
    WITH fm AS (
        SELECT m.pmcid, m.ref_id, min(m.mention_idx) AS i FROM '{EXPORT}/mentions.parquet' m JOIN '{KEPT}' USING (pmcid)
        WHERE m.sec_class = 'methods' AND m.in_table IS NOT TRUE AND m.in_caption IS NOT TRUE GROUP BY 1, 2
    )
    SELECT fm.pmcid, p.tool, min(fm.i) AS i
    FROM fm JOIN '{EXPORT}/citations.parquet' c USING (pmcid, ref_id) JOIN '{DATA}/tool_papers.parquet' p ON p.work_id = c.work_id
    GROUP BY 1, 2 ORDER BY 1, 3
""").df()
middle, ends = Counter(), Counter()
flank = defaultdict(Counter)
for _, g in ranked.groupby('pmcid'):
    seq = list(g.tool)
    if len(seq) < 3:
        for t in seq:
            ends[t] += 1
        continue
    for i, t in enumerate(seq):
        if 0 < i < len(seq) - 1:
            middle[t] += 1
            flank[t][(seq[i - 1], seq[i + 1])] += 1
        else:
            ends[t] += 1
ops = edam[edam.kind == 'operation'].groupby('tool').t.apply(set)
ad = []
for t in set(middle) | set(ends):
    m, e = middle[t], ends[t]
    if m + e < 100:
        continue
    adapter_ops = sorted(ops.get(t, set()) & ADAPTER_OPS)
    ad.append({'tool': t, 'in_middle': m, 'at_ends': e, 'middle_share': round(m / (m + e), 3),
               'adapter_operations': ', '.join(adapter_ops),
               'top_flanks': '; '.join(f'{a} > {t} > {b} ({n})' for (a, b), n in flank[t].most_common(3))})
adapters = pd.DataFrame(ad).sort_values('middle_share', ascending=False)

# ---------------------------------------------------------------- 3. takeovers with a shared operation
spec_ops = edam[(edam.kind == 'operation') & ~edam.t.isin(GENERIC_OPS)].groupby('tool').t.apply(set)
by_year = con.sql(f"SELECT * FROM '{DATA}/pairs_by_year.parquet'").df()


def neighbours(lo, hi):
    w = by_year[by_year.year.between(lo, hi)].groupby(['source', 'target']).n.sum()
    vec = defaultdict(dict)
    for (s, t), n in w.items():
        vec[s][t] = vec[s].get(t, 0) + n
        vec[t][s] = vec[t].get(s, 0) + n
    return vec


early, late, allv = neighbours(2012, 2017), neighbours(2021, 2026), neighbours(2012, 2026)
idf = {k: np.log(len(allv) / (1 + sum(1 for v in allv.values() if k in v))) for k in allv}


def cosine(a, b):
    a = {k: np.log1p(v) * idf.get(k, 0) for k, v in a.items()}
    b = {k: np.log1p(v) * idf.get(k, 0) for k, v in b.items()}
    num = sum(a[k] * b[k] for k in set(a) & set(b))
    den = np.sqrt(sum(v * v for v in a.values())) * np.sqrt(sum(v * v for v in b.values()))
    return num / den if den else 0.0


g = trends.growth_2020_2025
use_tot = tools
falling = [t for t in g.index if g[t] < -0.03 and sum(early.get(t, {}).values()) >= 30]
rising = [t for t in g.index if g[t] > 0.03 and sum(late.get(t, {}).values()) >= 30]
tk = []
for a in falling:
    for b in rising:
        shared = spec_ops.get(a, set()) & spec_ops.get(b, set())
        if a == b or not shared:
            continue
        together = pair_any.get(tuple(sorted((a, b))), 0) / max(1, min(use_tot.get(a, 1), use_tot.get(b, 1)))
        if together >= 0.2:
            continue
        sim = cosine({k: v for k, v in early.get(a, {}).items() if k != b}, {k: v for k, v in late.get(b, {}).items() if k != a})
        tk.append({'older': a, 'newer': b, 'shared_operation': ', '.join(sorted(shared)[:3]),
                   'neighbour_similarity': round(sim, 3), 'used_together': round(together, 3),
                   'older_growth': round(float(g[a]), 3), 'newer_growth': round(float(g[b]), 3),
                   'older_peak': int(trends.peak_year.get(a, 0)), 'newer_first_year': trends.first_year.get(b)})
take = pd.DataFrame(tk).sort_values('neighbour_similarity', ascending=False)
take = take[take.neighbour_similarity >= 0.3]
best = take.groupby('older').head(2)

for name, df in (('gaps_compatible', gaps), ('compatible_working', working.head(200)), ('adapters', adapters), ('takeovers_edam', best)):
    con.register('df', df)
    con.sql(f"COPY df TO '{DATA}/{name}.parquet' (FORMAT parquet)")

pd.set_option('display.width', 220)
pd.set_option('display.max_colwidth', 70)
fmt_tools = set(out_by.index) | set(in_by.index)
print(f'Tools with specific input/output types: {len(fmt_tools)} of {len(tools)}')
print(f'Compatible pairs (same kinds of papers): {len(compat)}; used together 2x+ more than chance: {(compat.ratio >= 2).sum()}; '
      f'rarely (< half of chance, >= 5 expected): {len(gaps)}')
print('\nCompatible but rarely combined:')
print(gaps.head(15)[['from_tool', 'to_tool', 'via', 'together_any_order', 'expected_if_independent', 'ratio']].to_string(index=False))
print('\nCompatible and used together most:')
print(working.head(8)[['from_tool', 'to_tool', 'via', 'together_any_order', 'expected_if_independent', 'ratio']].to_string(index=False))
print('\nTools most often in the middle of a sequence (100+ appearances):')
print(adapters.head(15)[['tool', 'in_middle', 'at_ends', 'middle_share', 'adapter_operations', 'top_flanks']].to_string(index=False))
print('\nTakeovers sharing a specific EDAM operation:')
print(best.head(25)[['older', 'newer', 'shared_operation', 'neighbour_similarity', 'used_together', 'older_growth', 'newer_growth']].to_string(index=False))
