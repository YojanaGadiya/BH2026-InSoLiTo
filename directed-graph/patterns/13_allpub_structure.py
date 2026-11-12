"""Step 13: structure of the all-publications tool graph, and how the omics subset differs.

Graph: the cleaned all-publications tool pairs of step 12 seen in at least MIN_ARTICLES
articles, as arrows in their main direction when at least MIN_SHARE of articles agree.

1. Communities (Louvain on the undirected weighted graph), hubs (PageRank, approximate
   betweenness), bridges (participation across communities x log strength) and order
   (source_share: ~1 used first, ~0 used last), as in step 9.
2. The omics subset: the proteomics/metabolomics edge list (citation-order step 05, same
   method, per bio.tools entry) mapped to the same tools. omics_share = share of a tool's
   or a pair's article links that come from omics papers. Omics-specific tools sit where
   in the full graph; which pairs are much more (or less) common in omics papers than in
   biology overall (lift = pair's omics share / overall omics share, with a floor so tiny
   pairs don't dominate); and whether the arrow points the same way in both.
3. Pipeline backbone: strongly one-way arrows (share >= CHAIN_SHARE, >= CHAIN_MIN
   articles). A step n -> m is kept when m is n's most common next tool and n is one of
   m's two most common previous tools; chains of these steps are the typical pipelines. Pairs used together a lot with no consistent order are listed too.

Writes patterns/data/allpub_hubs.parquet, allpub_community_links.parquet,
allpub_omics_pairs.parquet, allpub_chains.parquet, allpub_unordered.parquet.
Needs step 12. Run from BH2026-InSoLiTo/:
    python directed-graph/patterns/13_allpub_structure.py
"""

import duckdb
import networkx as nx
import numpy as np
import pandas as pd

DATA = 'directed-graph/patterns/data'
OMICS = 'directed-graph/citation-order/data/directed_edges.csv'
MIN_ARTICLES, MIN_SHARE = 10, 0.55
RESOLUTION, SEED, BT_SAMPLES = 2.0, 11, 400  # 1.0 puts proteomics, RNA-seq and enrichment in one community of 1,250 tools
CHAIN_SHARE, CHAIN_MIN = 0.8, 100

con = duckdb.connect()
con.sql("SET preserve_insertion_order = false")
edges = con.sql(f"SELECT * FROM '{DATA}/allpub_edges.parquet' WHERE n_total >= {MIN_ARTICLES}").df()
tools = con.sql(f"SELECT tool, entries, strength, citations_mostly_of FROM '{DATA}/allpub_tools.parquet'").df()

# ------------------------------------------------------------------ 1. communities and hubs
D = nx.DiGraph()
U = nx.Graph()
for r in edges.itertuples(index=False):
    U.add_edge(r.source, r.target, weight=float(r.n_total))
    if r.share_forward >= MIN_SHARE:
        D.add_edge(r.source, r.target, weight=float(r.n_fwd), distance=1.0 / r.n_total)
comms = nx.community.louvain_communities(U, weight='weight', seed=SEED, resolution=RESOLUTION)
comms = sorted(comms, key=lambda c: -sum(d for _, d in U.degree(c, weight='weight')))
comm = {n: i for i, c in enumerate(comms) for n in c}
pr = nx.pagerank(D, weight='weight')
bt = nx.betweenness_centrality(D, k=min(BT_SAMPLES, len(D)), weight='distance', seed=SEED)
rows = []
for n in U.nodes:
    k = {}
    for nb, d in U[n].items():
        k[comm[nb]] = k.get(comm[nb], 0) + d['weight']
    tot = sum(k.values())
    out_w = sum(d['weight'] for _, _, d in D.out_edges(n, data=True)) if n in D else 0
    in_w = sum(d['weight'] for _, _, d in D.in_edges(n, data=True)) if n in D else 0
    rows.append({'tool': n, 'community': comm[n], 'strength': tot, 'degree': U.degree(n),
                 'pagerank': pr.get(n, 0.0), 'betweenness': bt.get(n, 0.0),
                 'participation': 1 - sum((v / tot) ** 2 for v in k.values()),
                 'source_share': out_w / (out_w + in_w) if out_w + in_w else np.nan})
hubs = pd.DataFrame(rows)
hubs['bridge_score'] = hubs.participation * np.log1p(hubs.strength)
names = {c: ', '.join(hubs[hubs.community == c].nlargest(4, 'strength').tool) for c in hubs.community.unique()}
hubs['community_name'] = hubs.community.map(names)
cl = {}
for s, t, d in U.edges(data=True):
    a, b = sorted((comm[s], comm[t]))
    if a != b:
        x = cl.setdefault((a, b), [0.0, []])
        x[0] += d['weight']
        x[1].append((d['weight'], s, t))
links = pd.DataFrame([{'a': names[a], 'b': names[b], 'weight': w,
                       'top_links': '; '.join(f'{s}|{t} ({int(x)})' for x, s, t in sorted(l, reverse=True)[:4])}
                      for (a, b), (w, l) in cl.items()]).sort_values('weight', ascending=False)

# ------------------------------------------------------------------ 2. omics subset
emap = tools[['tool', 'entries']].explode('entries').rename(columns={'entries': 'entry'})
con.register('emap_df', emap)
con.sql("CREATE TABLE emap AS SELECT * FROM emap_df")
con.sql(f"""
    CREATE TABLE om_pairs AS
    SELECT a.tool AS s, b.tool AS t, max(r.n_articles) AS n, max(r.n_articles_rev) AS n_rev
    FROM read_csv('{OMICS}') r JOIN emap a ON a.entry = r.source JOIN emap b ON b.entry = r.target
    GROUP BY 1, 2
""")
con.sql("DELETE FROM om_pairs WHERE s = t")  # see step 12: not as a WHERE filter
om = con.sql("""
    SELECT p.s, p.t, greatest(p.n, coalesce(q.n_rev, 0)) AS n_fwd, greatest(p.n_rev, coalesce(q.n, 0)) AS n_rev
    FROM om_pairs p LEFT JOIN om_pairs q ON q.s = p.t AND q.t = p.s
""").df().set_index(['s', 't'])
e2 = edges.set_index(['source', 'target'])
om_fwd = om.n_fwd.reindex(e2.index).fillna(0)  # same orientation as the all-publications arrow
om_rev = om.n_rev.reindex(e2.index).fillna(0)
e2['omics_fwd'], e2['omics_rev'] = om_fwd.values, om_rev.values
# The omics run is not a strict subset (20 of 6,579 pairs with 20+ omics articles have more
# omics than all-publication articles, e.g. decontam, whose preprint the omics run maps),
# so omics counts are capped at the all-publications count.
scale = np.minimum(1.0, e2.n_total / (e2.omics_fwd + e2.omics_rev).replace(0, np.nan)).fillna(1.0)
e2['omics_fwd'], e2['omics_rev'] = e2.omics_fwd * scale, e2.omics_rev * scale
e2['omics_total'] = e2.omics_fwd + e2.omics_rev
e2 = e2.reset_index()
base = e2.omics_total.sum() / e2.n_total.sum()
e2['omics_share'] = e2.omics_total / e2.n_total
e2['lift'] = (e2.omics_total + 1) / (e2.n_total * base + 1)
e2['omics_share_forward'] = np.where(e2.omics_total > 0, e2.omics_fwd / e2.omics_total.replace(0, np.nan), np.nan)
tool_om = pd.concat([e2[['source', 'omics_total', 'n_total']].rename(columns={'source': 'tool'}),
                     e2[['target', 'omics_total', 'n_total']].rename(columns={'target': 'tool'})]).groupby('tool').sum()
hubs = hubs.join((tool_om.omics_total / tool_om.n_total).rename('omics_share'), on='tool')
hubs['omics_lift'] = hubs.omics_share / base

# ------------------------------------------------------------------ 3. pipeline backbone
strong = edges[(edges.share_forward >= CHAIN_SHARE) & (edges.n_total >= CHAIN_MIN)]
B = nx.DiGraph()
for r in strong.itertuples(index=False):
    B.add_edge(r.source, r.target, w=float(np.log(r.n_fwd)), n=int(r.n_fwd))
# A pipeline step: n -> m where m is n's most common next tool and n is one of m's two most
# common previous tools. Each tool then has at most one next step; chains are the paths.
nxt = {n: max(B.successors(n), key=lambda m: B[n][m]['n']) for n in B if B.out_degree(n)}
top_prev = {m: set(sorted(B.predecessors(m), key=lambda n: -B[n][m]['n'])[:2]) for m in B if B.in_degree(m)}
step = {n: m for n, m in nxt.items() if n in top_prev.get(m, set())}
has_prev = set(step.values())
chains = []
for start in step:
    if start in has_prev and any(step.get(p) == start and p != start for p in step):
        continue  # not a chain start
    path, seen = [start], {start}
    while path[-1] in step and step[path[-1]] not in seen:
        path.append(step[path[-1]])
        seen.add(path[-1])
    if len(path) >= 3:
        links_ = [B[a][b]['n'] for a, b in zip(path, path[1:])]
        chains.append({'chain': ' > '.join(path), 'length': len(path), 'weakest_link': min(links_),
                       'articles_per_step': ', '.join(map(str, links_)),
                       'community': names.get(comm.get(path[0]), '')})
chains = pd.DataFrame(chains).sort_values(['weakest_link'], ascending=False)
unordered = edges[(edges.share_forward < 0.55) & (edges.n_total >= 500)].sort_values('n_total', ascending=False)

for name, df in (('allpub_hubs', hubs), ('allpub_community_links', links), ('allpub_omics_pairs', e2),
                 ('allpub_chains', chains), ('allpub_unordered', unordered)):
    con.register('df', df)
    con.sql(f"COPY df TO '{DATA}/{name}.parquet' (FORMAT parquet)")

pd.set_option('display.width', 230)
pd.set_option('display.max_colwidth', 90)
print(f'Graph: {U.number_of_nodes()} tools, {U.number_of_edges()} pairs (>= {MIN_ARTICLES} articles), '
      f'{D.number_of_edges()} arrows (>= {MIN_SHARE:.0%} one way)')
print(f'{len(comms)} communities; the largest:')
for c in range(min(14, len(comms))):
    sub = hubs[hubs.community == c]
    print(f'  {c:2d}: {len(sub):4d} tools, {sub.strength.sum() / hubs.strength.sum():5.1%} of links, '
          f'omics share {sub.omics_share.mean():.2f}  {", ".join(sub.nlargest(8, "strength").tool)}')
cols = ['tool', 'community', 'pagerank', 'betweenness', 'participation', 'source_share', 'degree', 'omics_share']
print('\nPageRank:'); print(hubs.nlargest(12, 'pagerank')[cols].round(4).to_string(index=False))
print('\nBetweenness:'); print(hubs.nlargest(12, 'betweenness')[cols].round(4).to_string(index=False))
print('\nBridges (20+ links):'); print(hubs[hubs.degree >= 20].nlargest(15, 'bridge_score')[cols].round(3).to_string(index=False))
print(f'\nOmics papers: {base:.1%} of all article links')
print('Most omics-specific tools (strength >= 2000):')
print(hubs[hubs.strength >= 2000].nlargest(25, 'omics_share')[cols + ['strength']].round(3).to_string(index=False))
print('\nPairs most over-represented in omics papers (>= 50 omics articles):')
print(e2[e2.omics_total >= 50].nlargest(15, 'lift')[['source', 'target', 'n_total', 'omics_total', 'omics_share', 'lift']].round(2).to_string(index=False))
print('\nCommon pairs rare in omics (>= 2000 articles):')
print(e2[e2.n_total >= 2000].nsmallest(10, 'lift')[['source', 'target', 'n_total', 'omics_total', 'lift']].round(3).to_string(index=False))
agree = e2[(e2.omics_total >= 30) & (e2.share_forward >= 0.6)]
flip = agree[agree.omics_share_forward < 0.4]
print(f'\nDirection: of {len(agree)} clear arrows with 30+ omics articles, {len(flip)} point the other way in omics papers')
print(flip.nlargest(10, 'omics_total')[['source', 'target', 'n_total', 'share_forward', 'omics_total', 'omics_share_forward']].round(2).to_string(index=False))
print(f'\nBackbone: {B.number_of_nodes()} tools, {B.number_of_edges()} strong arrows; {len(chains)} chains of 3+ steps. Strongest:')
print(chains.head(30)[['chain', 'length', 'articles_per_step']].to_string(index=False))
print('\nUsed together a lot, no fixed order:')
print(unordered.head(12)[['source', 'target', 'n_total', 'share_forward']].round(2).to_string(index=False))
