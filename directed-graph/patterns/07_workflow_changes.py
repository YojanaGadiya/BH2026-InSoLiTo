"""Step 7: how workflows changed: the directed tool graph year by year.

Uses the cleaned tools of step 1 (one owner per paper). For each research article, its
top tools are put in the order of their first mention in the methods section (tables
and captions left out), as in citation-order/03, and every pair gives an arrow from the
tool cited first to the one cited later.

Writes to patterns/data/:
  pairs_by_year.parquet   source, target, year, n (articles with that order that year)
  workflow_frames.json    the animation: 3-year windows, a fixed layout, and per window
                          the share of articles using each tool and each arrow
  takeovers.parquet       tool pairs where a newer tool took an older one's place: the
                          same neighbours, rarely used together, one falling, one rising
  edge_turnover.parquet   per window: of the TOP_K arrows (by articles per 1,000), how many
                          were also in the window before and in the first window. A fixed
                          article count would mostly measure the corpus growing tenfold.

Needs steps 1-2. Run from BH2026-InSoLiTo/:
    python directed-graph/patterns/07_workflow_changes.py
"""

import json

import duckdb
import networkx as nx
import numpy as np
import pandas as pd

DATA = 'directed-graph/patterns/data'
EXPORT = 'directed-graph/europepmc/data/export'
KEPT = 'directed-graph/citation-order/data/kept_articles.parquet'
MANIFEST = 'directed-graph/europepmc/data/manifest.tsv'
WINDOWS = [(y - 1, y + 1) for y in range(2013, 2026)]   # 2012-14 ... 2024-26 (2026 partial)
MIN_EDGE = 5           # articles an arrow needs in a window to count there
TOP_K = 200            # arrows compared between windows for turnover
N_NODES = 80           # tools in the animation
N_EDGES = 170          # arrows in the animation
# Mass-spectrometry tools are used less than RNA-seq tools, so they would lose out to them;
# this corpus is about proteomics and metabolomics, so each keeps its 2 strongest arrows.
ALWAYS = ['maxquant', 'perseus', 'skyline', 'msfragger', 'DIA-NN', 'percolator', 'MASCOT', 'proteowizard',
          'comet', 'gnps', 'metaboanalyst', 'ms_dial', 'hmdb', 'canopus']
EARLY, LATE = (2012, 2017), (2021, 2026)

con = duckdb.connect()
con.sql(f"""
    CREATE TEMP TABLE ranked AS
    WITH first_mention AS (
        SELECT m.pmcid, m.ref_id, min(m.mention_idx) AS first_idx
        FROM '{EXPORT}/mentions.parquet' m JOIN '{KEPT}' k USING (pmcid)
        WHERE m.sec_class = 'methods' AND m.in_table IS NOT TRUE AND m.in_caption IS NOT TRUE
        GROUP BY 1, 2
    ), ref_rank AS (
        SELECT pmcid, ref_id, dense_rank() OVER (PARTITION BY pmcid ORDER BY first_idx) AS r FROM first_mention
    )
    SELECT c.pmcid, p.tool, min(rr.r) AS r
    FROM ref_rank rr JOIN '{EXPORT}/citations.parquet' c USING (pmcid, ref_id)
    JOIN '{DATA}/tool_papers.parquet' p ON p.work_id = c.work_id
    GROUP BY 1, 2
""")
con.sql(f"""
    CREATE TEMP TABLE pairs AS
    SELECT a.tool AS source, b.tool AS target, TRY_CAST(m.pub_year AS INTEGER) AS year, count(*) AS n
    FROM ranked a JOIN ranked b ON a.pmcid = b.pmcid AND a.r < b.r
    JOIN read_csv('{MANIFEST}', delim='\t', header=true, all_varchar=true) m ON m.pmcid = a.pmcid
    GROUP BY 1, 2, 3
""")
con.sql(f"COPY pairs TO '{DATA}/pairs_by_year.parquet' (FORMAT parquet)")
pairs = con.sql("SELECT * FROM pairs WHERE year BETWEEN 2012 AND 2026").df()
articles_per_year = con.sql(f"""
    SELECT TRY_CAST(m.pub_year AS INTEGER) AS year, count(*) AS n FROM '{KEPT}' k
    JOIN read_csv('{MANIFEST}', delim='\t', header=true, all_varchar=true) m USING (pmcid) GROUP BY 1
""").df().dropna().astype({'year': int}).set_index('year').n.sort_index()
use = con.sql(f"SELECT tool, year FROM '{DATA}/tool_use.parquet' WHERE year BETWEEN 2012 AND 2026").df()
tools = con.sql(f"SELECT * FROM '{DATA}/tools.parquet'").df()
trends = con.sql(f"SELECT tool, shape, growth_2020_2025, field_lean FROM '{DATA}/tool_trends.parquet'").df().set_index('tool')
families = con.sql(f"SELECT * FROM '{DATA}/tool_families.parquet'").df().set_index('tool').family


def window_sum(df, lo, hi, keys):
    return df[df.year.between(lo, hi)].groupby(keys).n.sum()


# ---------------------------------------------------------------- overall direction
tot = pairs.groupby(['source', 'target']).n.sum()
rev = tot.rename_axis(['target', 'source']).reorder_levels(['source', 'target'])
und = pd.DataFrame({'fwd': tot, 'rev': rev.reindex(tot.index).fillna(0)})
und['share'] = und.fwd / (und.fwd + und.rev)
main_dir = und[(und.share >= 0.5) & ((und.share > 0.5) | (und.index.get_level_values(0) < und.index.get_level_values(1)))]

# ---------------------------------------------------------------- turnover per window
edge_sets, turnover = [], []
for lo, hi in WINDOWS:
    w = window_sum(pairs, lo, hi, ['source', 'target'])
    w = w[w.index.isin(main_dir.index) & (w >= MIN_EDGE)]
    edge_sets.append(set(w.sort_values(ascending=False).head(TOP_K).index))
for i, (lo, hi) in enumerate(WINDOWS):
    prev = edge_sets[i - 1] if i else edge_sets[0]
    cur = edge_sets[i]
    turnover.append({'window': f'{lo}-{hi}', 'top_arrows': len(cur), 'also_in_window_before': len(cur & prev),
                     'also_in_first_window': len(cur & edge_sets[0])})
turnover = pd.DataFrame(turnover)

# ---------------------------------------------------------------- animation selection
use_w = {(lo, hi): use[use.year.between(lo, hi)].groupby('tool').size() / articles_per_year.loc[lo:hi].sum() * 1000
         for lo, hi in WINDOWS}
peak_use = pd.concat(use_w, axis=1).max(axis=1)
# Articles per window for each arrow in its main direction.
wins = []
for lo, hi in WINDOWS:
    wins.append(window_sum(pairs, lo, hi, ['source', 'target']).rename(f'{lo}'))
per_win = pd.concat(wins, axis=1).fillna(0)
per_win = per_win.loc[per_win.index.intersection(main_dir.index)]
cand = per_win[(per_win.max(axis=1) >= 20)]
cand = cand.assign(peak=cand.max(axis=1)).join(main_dir.share)
cand = cand[cand.share >= 0.6]
# Keep the strongest arrows per window, so tools that only mattered early or late still show.
chosen = set()
for col in per_win.columns:
    chosen |= set(cand.sort_values(col, ascending=False).head(N_EDGES // 6).index)
chosen |= set(cand.sort_values('peak', ascending=False).head(N_EDGES).index)
small = per_win[per_win.max(axis=1) >= 8]  # MS tools' arrows are smaller: a lower floor for them
small = small.assign(peak=small.max(axis=1)).join(main_dir.share)
small = small[small.share >= 0.6]
cand = pd.concat([cand, small[~small.index.isin(cand.index)]])
forced = set()
for t in ALWAYS:
    mine = small[[t in e for e in small.index]]
    forced |= set(mine.sort_values('peak', ascending=False).head(2).index)
nodes = pd.Series([n for e in chosen for n in e]).value_counts()
node_set = set(peak_use.reindex(nodes.index).sort_values(ascending=False).head(N_NODES).index)
node_set |= {n for e in forced for n in e}
chosen = [e for e in chosen | forced if e[0] in node_set and e[1] in node_set]
chosen = sorted(chosen, key=lambda e: (e not in forced, -cand.loc[e, 'peak']))[:N_EDGES + len(forced)]
node_set = sorted({n for e in chosen for n in e})

# Fixed layout from the whole period, so changes show as fading in place. It uses every
# link between the shown tools (not only the shown arrows, which fall apart into pieces),
# on a log scale so a few huge links don't dominate.
G = nx.Graph()
G.add_nodes_from(node_set)
shown = set(node_set)
for (s, t), n in pairs.groupby(['source', 'target']).n.sum().items():
    if s in shown and t in shown and n >= 5:
        w = G.get_edge_data(s, t, {'weight': 0})['weight'] + float(n)
        G.add_edge(s, t, weight=w)
for _, _, d in G.edges(data=True):
    d['weight'] = float(np.log1p(d['weight']))
pos = nx.spring_layout(G, weight='weight', seed=7, k=2.2 / np.sqrt(len(G)), iterations=600)
xy = np.array([pos[n] for n in node_set])
lo_, hi_ = np.percentile(xy, 2, axis=0), np.percentile(xy, 98, axis=0)   # robust to stragglers
xy = np.clip((xy - lo_) / (hi_ - lo_), -0.04, 1.04)
xy = (xy + 0.04) / 1.08
norm = {n: [round(float(x), 4), round(float(y), 4)] for n, (x, y) in zip(node_set, xy)}

frames = {
    'windows': [f'{lo}-{hi}' for lo, hi in WINDOWS],
    'partial_last': True,
    'nodes': [{'id': n, 'xy': norm[n], 'family': int(families.get(n, 0)), 'shape': trends['shape'].get(n),
               'use': [round(float(use_w[w].get(n, 0)), 2) for w in WINDOWS]} for n in node_set],
    'edges': [{'s': s, 't': t, 'share': round(float(cand.loc[(s, t), 'share']), 3),
               'n': [int(per_win.loc[(s, t), f'{lo}']) for lo, hi in WINDOWS],
               'per1000': [round(float(per_win.loc[(s, t), f'{lo}'] / articles_per_year.loc[lo:hi].sum() * 1000), 3)
                           for lo, hi in WINDOWS]} for s, t in chosen],
}
with open(f'{DATA}/workflow_frames.json', 'w') as fh:
    json.dump(frames, fh)

# ---------------------------------------------------------------- takeovers
def neighbour_vec(lo, hi):
    w = window_sum(pairs, lo, hi, ['source', 'target'])
    vec = {}
    for (s, t), n in w.items():
        vec.setdefault(s, {})[t] = vec.get(s, {}).get(t, 0) + n
        vec.setdefault(t, {})[s] = vec.get(t, {}).get(s, 0) + n
    return vec


early, late = neighbour_vec(*EARLY), neighbour_vec(*LATE)
co = {}
both = pairs.groupby(['source', 'target']).n.sum()
for (s, t), n in both.items():
    key = tuple(sorted((s, t)))
    co[key] = co.get(key, 0) + n
use_tot = use.groupby('tool').size()


# Hub neighbours (Cytoscape, KEGG...) are shared by almost everyone: weight each neighbour
# by how specific it is (inverse document frequency), so shared specific neighbours count.
all_vec = neighbour_vec(2012, 2026)
n_tools = len(all_vec)
idf = {k: np.log(n_tools / (1 + sum(1 for v in all_vec.values() if k in v))) for k in all_vec}


def cosine(a, b):
    a = {k: np.log1p(v) * idf.get(k, 0) for k, v in a.items()}
    b = {k: np.log1p(v) * idf.get(k, 0) for k, v in b.items()}
    keys = set(a) & set(b)
    num = sum(a[k] * b[k] for k in keys)
    den = np.sqrt(sum(v * v for v in a.values())) * np.sqrt(sum(v * v for v in b.values()))
    return num / den if den else 0.0


g = trends.growth_2020_2025
falling = [t for t in g.index if g[t] < -0.05 and t in early and sum(early[t].values()) >= 50]
rising = [t for t in g.index if g[t] > 0.05 and t in late and sum(late[t].values()) >= 50]
rows = []
for a in falling:
    for b in rising:
        if a == b:
            continue
        ea = {k: v for k, v in early[a].items() if k != b}
        lb = {k: v for k, v in late[b].items() if k != a}
        sim = cosine(ea, lb)
        together = co.get(tuple(sorted((a, b))), 0) / min(use_tot.get(a, 1), use_tot.get(b, 1))
        rows.append({'older': a, 'newer': b, 'neighbour_similarity': round(sim, 3), 'used_together': round(together, 3),
                     'older_growth': round(float(g[a]), 3), 'newer_growth': round(float(g[b]), 3)})
take = pd.DataFrame(rows)
take = take[take.used_together < 0.15].sort_values('neighbour_similarity', ascending=False)
best = take.groupby('older').head(1).head(25)

for name, df in (('takeovers', best), ('edge_turnover', turnover)):
    con.register('df', df)
    con.sql(f"COPY df TO '{DATA}/{name}.parquet' (FORMAT parquet)")

print(f'Top {TOP_K} arrows per window: how many were also in the window before / in the first window:')
print(turnover.to_string(index=False))
e0, e1 = edge_sets[0], edge_sets[-2]
print(f'\nOf the top {len(e1)} arrows in {WINDOWS[-2][0]}-{WINDOWS[-2][1]}, {len(e1 - e0)} ({len(e1 - e0) / len(e1):.0%}) '
      f'were not in the top {TOP_K} of {WINDOWS[0][0]}-{WINDOWS[0][1]}')
print(f'\nAnimation: {len(node_set)} tools, {len(chosen)} arrows, {len(WINDOWS)} windows')
print('\nLikely takeovers (newer tool in the older one\'s place, rarely used together):')
print(best.to_string(index=False))
