"""Step 15: pipeline steps, pipeline position and breadth of use in the all-publications graph.

A. Real steps vs loose mentions. The edge list counts every pair of tools in an article,
   not only neighbours, and keeps sum(1 / tool_gap): tool_gap = 1 when no other tool is
   cited between the two. adjacency = that sum / articles (1 = always next to each other,
   0.2 = typically 5 tools apart). Pairs are compared on counts and on step weight
   (the 1 / tool_gap sum), and the pipeline chains of step 13 are rebuilt on step weight.
B. Pipeline position. A SpringRank score per tool from the step-weighted arrows (each arrow
   pulls its source above its target), turned into a 0-1 position (0 = cited first, 1 =
   last, by rank among tools). Compared with the tools' EDAM operations from step 8: the
   mean position per operation, and tools far from where their operations put them.
C. Breadth. How evenly a tool's links spread over the 14 workflow areas of step 13
   (normalised entropy: 0 = one area, 1 = all areas equally) and in how many areas it has
   at least 5% of its links. For the omics-specific tools: where their links outside omics
   papers go.

Writes patterns/data/allpub_steps.parquet (pairs), allpub_position.parquet (tools),
allpub_position_edam.parquet (operations), allpub_breadth.parquet (tools),
allpub_omics_crossing.parquet, allpub_chains_steps.parquet. Needs steps 8, 12 and 13.
Run from BH2026-InSoLiTo/:
    python directed-graph/patterns/15_allpub_steps.py
"""

from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.sparse.linalg as spla

DATA = 'directed-graph/patterns/data'
RAW = 'directed-graph/citation-order/data/all_publications/directed_edges.parquet'
MIN_ARTICLES = 10
CHAIN_SHARE, CHAIN_MIN_STEP = 0.8, 60   # step weight ~ articles where the two are neighbours
ALPHA = 0.1                             # SpringRank regularisation
AREAS = 14

con = duckdb.connect()
con.sql("SET preserve_insertion_order = false")
tools = con.sql(f"SELECT tool, entries FROM '{DATA}/allpub_tools.parquet'").df()
emap = tools.explode('entries').rename(columns={'entries': 'entry'})
con.register('emap_df', emap)
con.sql("CREATE TABLE emap AS SELECT * FROM emap_df")
# Merged entries count the same articles: per tool pair, keep the entry pair with most articles.
con.sql(f"""
    CREATE TABLE p AS
    SELECT a.tool AS s, b.tool AS t, max(r.n_articles) AS n, arg_max(r.sum_inv_tool_gap, r.n_articles) AS g,
           arg_max(r.sum_inv_delta, r.n_articles) AS d, arg_max(r.n_articles_rev, r.n_articles) AS n_rev,
           arg_max(r.sum_inv_tool_gap_rev, r.n_articles) AS g_rev, arg_max(r.sum_inv_delta_rev, r.n_articles) AS d_rev
    FROM '{RAW}' r JOIN emap a ON a.entry = r.source JOIN emap b ON b.entry = r.target
    GROUP BY 1, 2
""")
con.sql("DELETE FROM p WHERE s = t")  # not as WHERE: see step 12
edges = con.sql(f"SELECT source, target, n_fwd, n_rev, n_total, share_forward FROM '{DATA}/allpub_edges.parquet' WHERE n_total >= {MIN_ARTICLES}").df()
pp = con.sql("SELECT * FROM p").df().set_index(['s', 't'])
fw = pp.reindex(pd.MultiIndex.from_frame(edges[['source', 'target']]))
bw = pp.reindex(pd.MultiIndex.from_arrays([edges.target, edges.source]))
# forward-direction sums from the (s, t) row, or the reverse columns of the (t, s) row
edges['step_fwd'] = np.fmax(fw.g.values, bw.g_rev.values)
edges['step_rev'] = np.fmax(fw.g_rev.values, bw.g.values)
edges['close_sum'] = np.fmax(fw.d.values, bw.d_rev.values) + np.fmax(fw.d_rev.values, bw.d.values)
edges[['step_fwd', 'step_rev', 'close_sum']] = edges[['step_fwd', 'step_rev', 'close_sum']].fillna(0)
edges['step_total'] = edges.step_fwd + edges.step_rev
edges['adjacency'] = (edges.step_total / edges.n_total).clip(upper=1)
edges['closeness'] = (edges.close_sum / edges.n_total).clip(upper=1)
edges['step_share_forward'] = np.where(edges.step_total > 0, edges.step_fwd / edges.step_total.replace(0, np.nan), np.nan)

# ---------------------------------------------------------------- A. steps vs mentions
hubs = con.sql(f"SELECT tool, community, strength, omics_share FROM '{DATA}/allpub_hubs.parquet'").df().set_index('tool')
comm = hubs.community
big = edges[edges.n_total >= 200]
w = big.n_total
sharp_n = float((np.abs(2 * big.share_forward - 1) * w).sum() / w.sum())
sharp_s = float((np.abs(2 * big.step_share_forward - 1) * w).sum() / w.sum())
loose = big[big.adjacency < 0.25].sort_values('n_total', ascending=False)
tight = big[big.adjacency >= 0.6].sort_values('n_total', ascending=False)
strong_n = edges[(edges.share_forward >= 0.8) & (edges.n_total >= 100)]
# chains on step weight (same rule as step 13)
st = edges[(edges.step_share_forward >= CHAIN_SHARE) & (edges.step_fwd >= CHAIN_MIN_STEP)]
succ = st.sort_values('step_fwd', ascending=False).groupby('source').head(1).set_index('source').target.to_dict()
prev2 = st.sort_values('step_fwd', ascending=False).groupby('target').head(2).groupby('target').source.apply(set).to_dict()
stepmap = {a: b for a, b in succ.items() if a in prev2.get(b, set())}
w_of = {(r.source, r.target): r.step_fwd for r in st.itertuples(index=False)}
has_prev = set(stepmap.values())
chains = []
for a in stepmap:
    if a in has_prev:
        continue
    path, seen = [a], {a}
    while path[-1] in stepmap and stepmap[path[-1]] not in seen:
        path.append(stepmap[path[-1]])
        seen.add(path[-1])
    if len(path) >= 3:
        ws = [w_of[(x, y)] for x, y in zip(path, path[1:])]
        chains.append({'chain': ' > '.join(path), 'length': len(path), 'weakest_step': round(min(ws)),
                       'step_weights': ', '.join(str(round(x)) for x in ws)})
chains = pd.DataFrame(chains).sort_values(['length', 'weakest_step'], ascending=False)

# ---------------------------------------------------------------- B. pipeline position
nodes = sorted(set(edges.source) | set(edges.target))
ix = {n: i for i, n in enumerate(nodes)}
N = len(nodes)
A = sp.csr_matrix((np.r_[edges.step_fwd, edges.step_rev],
                   (np.r_[edges.source.map(ix), edges.target.map(ix)], np.r_[edges.target.map(ix), edges.source.map(ix)])),
                  shape=(N, N))
kout, kin = np.asarray(A.sum(1)).ravel(), np.asarray(A.sum(0)).ravel()
L = sp.diags(kout + kin + ALPHA) - (A + A.T)
score = spla.spsolve(L.tocsc(), kout - kin)  # SpringRank: higher = earlier in the pipeline
pos = pd.DataFrame({'tool': nodes, 'springrank': score, 'out_step': kout, 'in_step': kin})
pos['strength'] = pos.tool.map(hubs.strength)
pos = pos[(pos.out_step + pos.in_step) >= 30]
pos['position'] = 1 - pos.springrank.rank(pct=True)  # 0 = first, 1 = last
pos['community'] = pos.tool.map(comm)

edam_file = Path(DATA) / 'edam_functions_allpub.parquet'
if not edam_file.exists():
    edam_file = Path(DATA) / 'edam_functions.parquet'
edam = con.sql(f"SELECT DISTINCT tool, term FROM '{edam_file}' WHERE kind = 'operation'").df()
# step-8 tool names are the omics ones; map their entries to the all-publications tools
if edam_file.name == 'edam_functions.parquet':
    ent = con.sql(f"SELECT DISTINCT tool, entry FROM '{edam_file}'").df()
    e2t = emap.set_index('entry').tool.to_dict()
    ent['allpub_tool'] = ent.entry.map(e2t)
    edam = edam.merge(ent[['tool', 'allpub_tool']].dropna().drop_duplicates(), on='tool')[['allpub_tool', 'term']] \
               .rename(columns={'allpub_tool': 'tool'}).drop_duplicates()
ep = edam.merge(pos[['tool', 'position']], on='tool')
ops = ep.groupby('term').position.agg(['mean', 'count', 'std']).reset_index()
ops = ops[ops['count'] >= 5].sort_values('mean')
# How far each tool sits from the mean position of its own operations
tool_exp = ep.merge(ops[['term', 'mean']], on='term').groupby('tool')['mean'].mean().rename('expected')
pos = pos.join(tool_exp, on='tool')
pos['off_by'] = pos.position - pos.expected
r_ops = float(np.corrcoef(pos.dropna(subset=['expected']).position, pos.dropna(subset=['expected']).expected)[0, 1])

# ---------------------------------------------------------------- C. breadth
om = con.sql(f"SELECT source, target, n_total, omics_total FROM '{DATA}/allpub_omics_pairs.parquet'").df()
both = pd.concat([om.rename(columns={'source': 'tool', 'target': 'nb'}), om.rename(columns={'target': 'tool', 'source': 'nb'})])
both['area'] = both.nb.map(comm).where(lambda c: c < AREAS, AREAS)  # small communities pooled as "other"
area_w = both.groupby(['tool', 'area']).n_total.sum().unstack(fill_value=0)
P = area_w.div(area_w.sum(1), axis=0)
ent_ = -(P * np.log(P.where(P > 0, 1))).sum(1) / np.log(P.shape[1])
breadth = pd.DataFrame({'entropy': ent_, 'areas_5pct': (P >= 0.05).sum(1), 'strength': area_w.sum(1)})
breadth['home'] = breadth.index.map(comm)
breadth['home_share'] = [P.loc[t, comm[t]] if comm.get(t, AREAS) in P.columns and comm[t] < AREAS else np.nan for t in breadth.index]
breadth = breadth.join(hubs.omics_share)
breadth = breadth[breadth.strength >= 500].reset_index(names='tool')
both['non_omics'] = (both.n_total - both.omics_total).clip(lower=0)
spec = breadth[(breadth.omics_share >= 0.5) & (breadth.strength >= 1000)].tool
cross = both[both.tool.isin(spec)].groupby(['tool', 'area'])[['non_omics', 'omics_total']].sum().reset_index()
cross_top = both[both.tool.isin(spec)].sort_values('non_omics', ascending=False).groupby('tool').head(4) \
    .groupby('tool').apply(lambda g: ', '.join(f'{r.nb} ({int(r.non_omics)})' for r in g.itertuples()), include_groups=False)

for name, df in (('allpub_steps', edges), ('allpub_position', pos), ('allpub_position_edam', ops),
                 ('allpub_breadth', breadth), ('allpub_omics_crossing', cross), ('allpub_chains_steps', chains)):
    con.register('df', df)
    con.sql(f"COPY df TO '{DATA}/{name}.parquet' (FORMAT parquet)")

pd.set_option('display.width', 230)
pd.set_option('display.max_colwidth', 120)
print(f'A. {len(edges):,} pairs (>= {MIN_ARTICLES} articles). Median adjacency {edges.adjacency.median():.2f}; '
      f'pairs with 200+ articles: {len(big):,}, of which {(big.adjacency >= 0.6).mean():.0%} tight (>= 0.6) and {(big.adjacency < 0.25).mean():.0%} loose (< 0.25)')
print(f'   Correlation of adjacency with articles (log): {np.corrcoef(np.log(edges.n_total), edges.adjacency)[0, 1]:.2f}')
print(f'   Direction sharpness |2*share-1| (weighted, 200+ articles): counts {sharp_n:.2f}, step weight {sharp_s:.2f}')
print(f'   Strong count arrows (>= 80% one way, 100+ articles): {len(strong_n):,}; of these loose (< 0.25): {(strong_n.adjacency < 0.25).mean():.0%}')
print('\nLoose: used together in many articles, rarely next to each other:')
print(loose.head(12)[['source', 'target', 'n_total', 'adjacency', 'share_forward']].round(2).to_string(index=False))
print('\nTight: usually cited one right after the other:')
print(tight.head(12)[['source', 'target', 'n_total', 'adjacency', 'share_forward', 'step_share_forward']].round(2).to_string(index=False))
print(f'\nChains on step weight ({len(chains)}):')
print(chains.head(25).to_string(index=False))
print(f'\nB. Positions for {len(pos):,} tools. EDAM from {edam_file.name}: {ep.tool.nunique()} tools with operations; '
      f'correlation of position with operation mean: {r_ops:.2f}')
print('Earliest operations:'); print(ops.head(12).round(2).to_string(index=False))
print('Latest operations:'); print(ops.tail(12).round(2).to_string(index=False))
print('First tools:'); print(pos[pos.strength >= 5000].nsmallest(12, 'position')[['tool', 'position', 'out_step', 'in_step']].round(2).to_string(index=False))
print('Last tools:'); print(pos[pos.strength >= 5000].nlargest(12, 'position')[['tool', 'position', 'out_step', 'in_step']].round(2).to_string(index=False))
print('Furthest from their EDAM operations:')
print(pos.dropna(subset=['expected']).assign(a=lambda d: d.off_by.abs()).nlargest(15, 'a')[['tool', 'position', 'expected', 'off_by']].round(2).to_string(index=False))
print(f'\nC. Breadth for {len(breadth):,} tools (500+ links)')
print('Broadest:'); print(breadth[breadth.strength >= 5000].nlargest(12, 'entropy')[['tool', 'entropy', 'areas_5pct', 'strength', 'omics_share']].round(2).to_string(index=False))
print('Narrowest (5000+ links):'); print(breadth[breadth.strength >= 5000].nsmallest(10, 'entropy')[['tool', 'entropy', 'areas_5pct', 'home_share', 'strength']].round(2).to_string(index=False))
print('Omics tools: links outside omics papers');
sp_ = breadth[breadth.tool.isin(spec)].set_index('tool')
print(pd.DataFrame({'omics_share': sp_.omics_share.round(2), 'entropy': sp_.entropy.round(2), 'top_non_omics': cross_top}).sort_values('omics_share').to_string())
