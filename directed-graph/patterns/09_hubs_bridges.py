"""Step 9: hubs and bridges of the tool graph.

Graph: the top tools (step 1) and every arrow seen in at least MIN_ARTICLES articles,
in its main direction, weighted by articles (all years, step 7). On it:

  pagerank        where workflows flow to (weighted, directed)
  betweenness     how often a tool lies on the shortest paths between others
                  (distance = 1 / articles, so strong links are short)
  source_share    share of the tool's arrows that leave it: ~1 = used first (inputs,
                  aligners, converters), ~0 = used last (statistics, enrichment, plots)
  community       workflow community (Louvain on the undirected weighted graph)
  participation   how evenly a tool's links spread over communities (0 = all in its own;
                  near 1 = spread over many). Bridges = high participation and many links.

Also the strongest links between communities, and how bridging changed between the
early (2012-2017) and late (2021-2026) years.

Writes patterns/data/hubs.parquet, community_links.parquet. Needs steps 1, 2 and 7.
Run from BH2026-InSoLiTo/:
    python directed-graph/patterns/09_hubs_bridges.py
"""

import duckdb
import networkx as nx
import numpy as np
import pandas as pd

DATA = 'directed-graph/patterns/data'
MIN_ARTICLES = 10
MIN_SHARE = 0.55
SEED = 11
RESOLUTION = 1.5  # 1.0 merges RNA-seq, proteomics and enrichment into one community; 1.5 gives 8 that read as workflows

con = duckdb.connect()
pairs = con.sql(f"SELECT * FROM '{DATA}/pairs_by_year.parquet' WHERE year BETWEEN 2012 AND 2026").df()
tools = con.sql(f"SELECT tool, rank, n_articles FROM '{DATA}/tools.parquet'").df().set_index('tool')
trends = con.sql(f"SELECT tool, shape, field_lean FROM '{DATA}/tool_trends.parquet'").df().set_index('tool')


def build(lo, hi):
    w = pairs[pairs.year.between(lo, hi)].groupby(['source', 'target']).n.sum()
    rev = w.copy()
    rev.index = rev.index.swaplevel()
    df = pd.DataFrame({'fwd': w, 'rev': rev.reindex(w.index).fillna(0)})
    df['share'] = df.fwd / (df.fwd + df.rev)
    df = df[(df.fwd + df.rev >= MIN_ARTICLES) & (df.share >= MIN_SHARE)]
    D = nx.DiGraph()
    for (s, t), r in df.iterrows():
        D.add_edge(s, t, weight=float(r.fwd), distance=1.0 / float(r.fwd + r.rev))
    U = nx.Graph()
    for s, t, d in D.edges(data=True):
        prev = U.get_edge_data(s, t, {'weight': 0})['weight']
        U.add_edge(s, t, weight=prev + d['weight'])
    return D, U


def measures(D, U, communities=None):
    if communities is None:
        communities = nx.community.louvain_communities(U, weight='weight', seed=SEED, resolution=RESOLUTION)
    comm = {n: i for i, c in enumerate(sorted(communities, key=len, reverse=True)) for n in c}
    pr = nx.pagerank(D, weight='weight')
    bt = nx.betweenness_centrality(D, weight='distance', normalized=True, seed=SEED)
    rows = []
    for n in D.nodes:
        out_w = sum(d['weight'] for _, _, d in D.out_edges(n, data=True))
        in_w = sum(d['weight'] for _, _, d in D.in_edges(n, data=True))
        k = {}
        for nb, d in U[n].items():
            k[comm.get(nb, -1)] = k.get(comm.get(nb, -1), 0) + d['weight']
        tot = sum(k.values())
        part = 1 - sum((v / tot) ** 2 for v in k.values()) if tot else 0
        rows.append({'tool': n, 'pagerank': pr[n], 'betweenness': bt[n], 'degree': U.degree(n),
                     'strength': tot, 'source_share': out_w / (out_w + in_w) if out_w + in_w else np.nan,
                     'community': comm.get(n, -1), 'participation': part,
                     'n_communities_linked': sum(1 for c, v in k.items() if v / tot >= 0.05) if tot else 0})
    return pd.DataFrame(rows), comm, communities


D, U = build(2012, 2026)
hubs, comm, communities = measures(D, U)
hubs = hubs.join(tools, on='tool').join(trends, on='tool')
hubs['bridge_score'] = hubs.participation * np.log1p(hubs.strength)

# Name each community after its three most-used tools.
names = {c: ', '.join(hubs[hubs.community == c].nsmallest(3, 'rank').tool) for c in hubs.community.unique()}
hubs['community_name'] = hubs.community.map(names)

# Links between communities.
cl = {}
for s, t, d in U.edges(data=True):
    a, b = sorted((comm[s], comm[t]))
    if a != b:
        cl.setdefault((a, b), [0, []])
        cl[(a, b)][0] += d['weight']
        cl[(a, b)][1].append((d['weight'], s, t))
links = pd.DataFrame([{'a': names[a], 'b': names[b], 'weight': w,
                       'top_links': '; '.join(f'{s}|{t} ({int(x)})' for x, s, t in sorted(l, reverse=True)[:4])}
                      for (a, b), (w, l) in cl.items()]).sort_values('weight', ascending=False)

# Bridging early vs late, with the full-period communities kept fixed.
change = []
for lo, hi in ((2012, 2017), (2021, 2026)):
    De, Ue = build(lo, hi)
    m, _, _ = measures(De, Ue, communities=[c & set(De.nodes) for c in communities if c & set(De.nodes)])
    change.append(m.set_index('tool')[['participation', 'pagerank', 'degree']].add_suffix(f'_{lo}'))
change = pd.concat(change, axis=1)
hubs = hubs.join(change, on='tool')

for name, df in (('hubs', hubs), ('community_links', links)):
    con.register('df', df)
    con.sql(f"COPY df TO '{DATA}/{name}.parquet' (FORMAT parquet)")

pd.set_option('display.width', 200)
print(f'Graph: {D.number_of_nodes()} tools, {D.number_of_edges()} arrows (>= {MIN_ARTICLES} articles, >= {MIN_SHARE:.0%} one way)')
print(f'{len(communities)} communities:')
for c in sorted(hubs.community.unique()):
    sub = hubs[hubs.community == c]
    print(f'  {c}: {len(sub):3d} tools  {", ".join(sub.nsmallest(8, "rank").tool)}')
cols = ['tool', 'community_name', 'pagerank', 'betweenness', 'participation', 'source_share', 'degree']
print('\nTop hubs by PageRank:')
print(hubs.nlargest(12, 'pagerank')[cols].round(3).to_string(index=False))
print('\nTop by betweenness:')
print(hubs.nlargest(12, 'betweenness')[cols].round(3).to_string(index=False))
print('\nTop bridges (participation x log strength), at least 10 links:')
print(hubs[hubs.degree >= 10].nlargest(15, 'bridge_score')[cols].round(3).to_string(index=False))
print('\nMost upstream (source_share), at least 10 links:')
print(hubs[hubs.degree >= 10].nlargest(8, 'source_share')[['tool', 'source_share', 'community_name']].round(2).to_string(index=False))
print('Most downstream:')
print(hubs[hubs.degree >= 10].nsmallest(8, 'source_share')[['tool', 'source_share', 'community_name']].round(2).to_string(index=False))
print('\nStrongest links between communities:')
print(links.head(10).to_string(index=False))
d = hubs[(hubs.degree_2012 >= 8) & (hubs.degree_2021 >= 8)]  # present in both periods
d = d.assign(change=d.participation_2021 - d.participation_2012)
print('\nBiggest gains in bridging, 2012-17 -> 2021-26 (tools with 8+ links in both):')
print(d.nlargest(8, 'change')[['tool', 'participation_2012', 'participation_2021']].round(2).to_string(index=False))
print('Biggest losses:')
print(d.nsmallest(6, 'change')[['tool', 'participation_2012', 'participation_2021']].round(2).to_string(index=False))
