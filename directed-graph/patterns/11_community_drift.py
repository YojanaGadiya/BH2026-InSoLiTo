"""Step 11: do tools move between application areas over time?

The communities of step 9 (Louvain on the whole 2012-2026 graph) are kept fixed as
"application areas". For each 3-year window, a tool's community mix is the share of its
workflow links (articles where it is cited right before or after another top tool in the
methods) that go to tools of each community. A tool whose mix shifts from one area to
another is being used in different kinds of workflows, even if its own community label
never changes.

The areas grow at different rates (sequence analysis was 37% of all link ends in
2012-16 and 17% in 2022-26; RNA-seq, statistics and metagenomics grew), so part of every
tool's drift is the field changing around it. `moved` is the raw shift (how the tool's
use actually changed); `moved_adjusted` divides each area's share by that area's share
of all links in the window first, so it counts only shifts beyond the field as a whole.
The adjusted one blows up small areas in the early years (microbiome was 2% of links),
so it is a check, not the main measure.

A small prior (PRIOR articles' worth of links to the tool's own community) keeps tools
with few links in a window from jumping around.

Writes patterns/data/community_mix.parquet (tool, window, n_links, mix_0..mix_7) and
community_movers.parquet (one row per tool: early and late mix, distance moved, main
area early and late). Needs steps 7 and 9. Run from BH2026-InSoLiTo/:
    python directed-graph/patterns/11_community_drift.py
"""

import duckdb
import numpy as np
import pandas as pd

DATA = 'directed-graph/patterns/data'
WINDOWS = [(y, y + 2) for y in range(2012, 2025)]  # same 3-year windows as step 7
PRIOR = 5          # articles' worth of links added to the tool's own community
MIN_LINKS = 30     # links needed in both the early and the late period to rank movers
EARLY, LATE = (2012, 2016), (2022, 2026)
AREA = ['Proteomics, annotation & networks', 'Sequence analysis & phylogenetics', 'RNA-seq', 'Structural biology',
        'Statistics, enrichment & single-cell', 'Metagenome assembly & annotation', 'Alignment, variants & ChIP-seq',
        'Microbiome & metabolomics']

con = duckdb.connect()
hubs = con.sql(f"SELECT tool, community, n_articles FROM '{DATA}/hubs.parquet'").df().set_index('tool')
comm = hubs.community.to_dict()
K = len(AREA)
pairs = con.sql(f"SELECT source, target, year, n FROM '{DATA}/pairs_by_year.parquet' WHERE year BETWEEN 2012 AND 2026").df()
pairs = pairs[pairs.source.isin(comm) & pairs.target.isin(comm) & (pairs.source != pairs.target)]
# Undirected: each article link counts for both tools, toward the other tool's community.
both = pd.concat([pairs.rename(columns={'source': 'tool', 'target': 'other'}),
                  pairs.rename(columns={'target': 'tool', 'source': 'other'})])
both['c'] = both.other.map(comm)


def mix(lo, hi, adjust=False):
    w = both[both.year.between(lo, hi)].groupby(['tool', 'c']).n.sum().unstack(fill_value=0)
    w = w.reindex(index=list(comm), columns=range(K), fill_value=0).astype(float)
    n = w.sum(axis=1)
    g = w.sum(axis=0) / w.values.sum()  # each area's share of all link ends in the window
    for t in w.index:
        w.loc[t, comm[t]] += PRIOR
    w = w.div(w.sum(axis=1), axis=0)
    if adjust:
        w = w / g
        w = w.div(w.sum(axis=1), axis=0)
    return w, n


rows = []
for lo, hi in WINDOWS:
    m, n = mix(lo, hi)
    for t in m.index:
        rows.append({'tool': t, 'window': f'{lo}-{hi}', 'n_links': int(n[t]), **{f'mix_{c}': round(m.loc[t, c], 4) for c in range(K)}})
mix_df = pd.DataFrame(rows)

(me, ne), (ml, nl) = mix(*EARLY), mix(*LATE)
(ae, _), (al, _) = mix(*EARLY, adjust=True), mix(*LATE, adjust=True)
mv = []
for t in comm:
    if ne[t] < MIN_LINKS or nl[t] < MIN_LINKS:
        continue
    a, b = me.loc[t].values, ml.loc[t].values
    gain = b - a
    mv.append({'tool': t, 'home': AREA[comm[t]], 'links_early': int(ne[t]), 'links_late': int(nl[t]),
               'moved': round(0.5 * np.abs(gain).sum(), 3),  # total variation: share of links that changed area
               'moved_adjusted': round(0.5 * np.abs(al.loc[t].values - ae.loc[t].values).sum(), 3),
               'main_early': AREA[int(a.argmax())], 'main_late': AREA[int(b.argmax())],
               'gained_most': AREA[int(gain.argmax())], 'gain': round(float(gain.max()), 3),
               'lost_most': AREA[int(gain.argmin())], 'loss': round(float(gain.min()), 3),
               'home_share_early': round(float(a[comm[t]]), 3), 'home_share_late': round(float(b[comm[t]]), 3)})
movers = pd.DataFrame(mv).sort_values('moved', ascending=False)

for name, df in (('community_mix', mix_df), ('community_movers', movers)):
    con.register('df', df)
    con.sql(f"COPY df TO '{DATA}/{name}.parquet' (FORMAT parquet)")

pd.set_option('display.width', 220)
print(f'{len(movers)} tools with {MIN_LINKS}+ links in both {EARLY[0]}-{EARLY[1]} and {LATE[0]}-{LATE[1]}')
print(f'Median share of links that changed area: {movers.moved.median():.2f} ({movers.moved_adjusted.median():.2f} adjusted for area growth)')
cols = ['tool', 'home', 'moved', 'moved_adjusted', 'main_early', 'main_late', 'gained_most', 'gain', 'lost_most', 'loss']
print('\nMoved most:')
print(movers.head(20)[cols].to_string(index=False))
print('\nMost fixed (among the 40 most used):')
top = movers[movers.tool.isin(hubs.nlargest(40, 'n_articles').index)]
print(top.nsmallest(10, 'moved')[cols].to_string(index=False))
for t in ('cytoscape', 'maxquant'):
    if t in movers.tool.values:
        print('\n', movers[movers.tool == t][cols].to_string(index=False))
