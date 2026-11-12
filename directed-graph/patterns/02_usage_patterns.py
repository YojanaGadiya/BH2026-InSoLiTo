"""Step 2: usage patterns of the top tools over time and together.

1. Adoption: for every tool and year, uses per 1,000 research articles of the tool's own
   field mix. The corpus grows over time and its mix shifts (proteomics-only articles
   fell from 83% in 2012 to 47% in 2026, metabolomics rose from 10% to 40%), so a raw
   share would make every proteomics tool look like it is declining. Instead each
   year's denominator weights the proteomics, metabolomics and both-topic articles of
   that year by the tool's own split of uses across them. Each tool's curve is scaled to its own peak and the tools are clustered
   by curve shape (Ward linkage), then each cluster is named from its average curve.
2. Co-usage: for the most-used tools, how much more often two tools are used in the
   same article than chance (normalised pointwise mutual information, -1 to 1),
   clustered into families of tools that travel together.
3. Field lean: how much more a tool is used in the proteomics or the metabolomics part
   of the corpus than the corpus as a whole.

Writes to directed-graph/patterns/data/ (tables) and patterns/figures/ (clustermaps).
Needs step 1. Run from BH2026-InSoLiTo/:
    python directed-graph/patterns/02_usage_patterns.py
"""

from pathlib import Path

import duckdb
import matplotlib
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage

matplotlib.use('Agg')
import seaborn as sns  # noqa: E402

YEARS = list(range(2012, 2027))  # 2026 runs to the search date (28 Sep 2026)
N_SHAPES = 6                     # adoption clusters
N_COUSAGE = 80                   # tools in the co-usage map
N_FAMILIES = 10                  # co-usage clusters
MIN_PAIR = 10                    # articles a pair needs before its NPMI counts

DATA = 'directed-graph/patterns/data'
FIG = 'directed-graph/patterns/figures'
MANIFEST = 'directed-graph/europepmc/data/manifest.tsv'
KEPT = 'directed-graph/citation-order/data/kept_articles.parquet'

Path(FIG).mkdir(parents=True, exist_ok=True)
con = duckdb.connect()
tools = con.sql(f"SELECT * FROM '{DATA}/tools.parquet' ORDER BY rank").df()
use = con.sql(f"SELECT * FROM '{DATA}/tool_use.parquet'").df()
corpus = con.sql(f"""
    SELECT TRY_CAST(a.pub_year AS INTEGER) AS year, a.topics
    FROM '{KEPT}' k JOIN read_csv('{MANIFEST}', delim='\t', header=true, all_varchar=true) a USING (pmcid)
""").df()

# ---------------------------------------------------------------- 1. adoption
GROUPS = ['proteomics', 'metabolomics', 'metabolomics,proteomics']
corpus['group'] = corpus.topics.where(corpus.topics.isin(GROUPS), 'proteomics')
use['group'] = use.topics.where(use.topics.isin(GROUPS), 'proteomics')
per_group_year = corpus.groupby(['group', 'year']).size().unstack(fill_value=0).reindex(index=GROUPS, columns=YEARS, fill_value=0)
mix = use.groupby(['tool', 'group']).size().unstack(fill_value=0).reindex(index=tools.tool, columns=GROUPS, fill_value=0)
mix = mix.div(mix.sum(axis=1), axis=0)                # tool's split of uses across fields
denominator = mix.values @ per_group_year.values       # tool x year: articles of its field mix
counts = (use[use.year.isin(YEARS)].groupby(['tool', 'year']).size()
          .unstack(fill_value=0).reindex(index=tools.tool, columns=YEARS, fill_value=0))
share = counts / pd.DataFrame(denominator, index=counts.index, columns=YEARS) * 1000
smooth = share.T.rolling(3, center=True, min_periods=1).mean().T
shape = smooth.div(smooth.max(axis=1).replace(0, np.nan), axis=0).fillna(0)

Z = linkage(shape.values, method='ward')
cluster = fcluster(Z, N_SHAPES, criterion='maxclust')


# First year with at least 5 uses (all years, not only YEARS): when a tool really arrived.
all_counts = use.dropna(subset=['year']).groupby(['tool', 'year']).size()
first_used = all_counts[all_counts >= 5].reset_index().groupby('tool').year.min().reindex(tools.tool)


def describe(curve: pd.Series, first: int) -> str:
    """Name a cluster from its average curve and its members' median first-use year."""
    peak = int(curve.idxmax())
    last = curve.loc[2024:2026].mean()
    if first >= 2019:
        return f'Steep rise: most took off around {first}'
    if peak >= 2023 and curve.loc[2012:2016].mean() < 0.35:
        return 'Rising: still growing'
    if peak >= 2022:
        return 'Established: high and steady'
    if last < 0.55:
        return f'Declining: peaked around {peak}'
    return f'Plateaued: flat since about {peak}'


centroids = shape.groupby(cluster).mean()
cluster_first = pd.Series(first_used.values, index=shape.index).groupby(cluster).median()
names = {c: describe(centroids.loc[c], int(cluster_first.loc[c])) for c in centroids.index}
log_share = np.log(share.loc[:, 2020:2025] + 0.05)
slope = pd.Series(np.polyfit(range(2020, 2026), log_share.T.values, 1)[0], index=share.index)

# ---------------------------------------------------------------- 3. field lean
def lean(topics: pd.Series) -> pd.Series:
    t = topics.fillna('')
    return pd.Series({'proteomics_only': (t == 'proteomics').mean(),
                      'metabolomics_only': (t == 'metabolomics').mean(),
                      'both': (t == 'metabolomics,proteomics').mean()})


base = lean(corpus.topics)
by_tool = use.groupby('tool').topics.apply(lean).unstack()
field_lean = np.log2((by_tool.proteomics_only + 0.01) / (base.proteomics_only + 0.01)) - \
    np.log2((by_tool.metabolomics_only + 0.01) / (base.metabolomics_only + 0.01))

trends = tools[['tool', 'rank', 'n_articles']].copy()
trends['shape'] = [names[c] for c in cluster]
trends['shape_id'] = cluster
trends['peak_year'] = smooth.idxmax(axis=1).values
trends['first_year'] = first_used.values
trends['growth_2020_2025'] = (np.exp(slope.values) - 1).round(3)  # yearly change in share
trends = trends.join(by_tool, on='tool')
trends['field_lean'] = trends.tool.map(field_lean).round(2)       # >0 proteomics, <0 metabolomics

con.register('trends', trends)
con.sql(f"COPY trends TO '{DATA}/tool_trends.parquet' (FORMAT parquet)")
long = share.stack().rename('per_1000').reset_index()
long['n'] = counts.stack().values
con.register('tool_year', long)
con.sql(f"COPY tool_year TO '{DATA}/tool_year.parquet' (FORMAT parquet)")

order = trends.sort_values(['shape_id', 'peak_year', 'rank']).tool
show = trends.nsmallest(120, 'rank').tool
rows = [t for t in order if t in set(show)]
palette = dict(zip(sorted(set(cluster)), sns.color_palette('Set2', N_SHAPES)))
g = sns.clustermap(shape.loc[rows], row_linkage=None, row_cluster=False, col_cluster=False,
                   cmap='mako_r', figsize=(10, 26), yticklabels=True, xticklabels=True,
                   row_colors=[palette[c] for c in trends.set_index('tool').loc[rows, 'shape_id']],
                   cbar_kws={'label': 'share of peak year'})
g.ax_heatmap.set_ylabel('')
g.figure.suptitle('Adoption of the 120 most-used tools (uses per 1,000 articles of their field, scaled to each tool\'s peak)', y=1.0)
g.savefig(f'{FIG}/adoption_clustermap.png', dpi=130, bbox_inches='tight')

# ---------------------------------------------------------------- 2. co-usage
top = tools.nsmallest(N_COUSAGE, 'rank').tool.tolist()
arts = use[use.tool.isin(top)].groupby('pmcid').tool.apply(set)
n_all = corpus.shape[0]
idx = {t: i for i, t in enumerate(top)}
M = np.zeros((len(top), len(arts)), dtype=np.uint8)
for j, ts in enumerate(arts):
    for t in ts:
        M[idx[t], j] = 1
co = M.astype(np.int32) @ M.T.astype(np.int32)
p = np.diag(co) / n_all
pij = co / n_all
with np.errstate(divide='ignore', invalid='ignore'):
    npmi = np.log(pij / np.outer(p, p)) / -np.log(pij)
npmi[co < MIN_PAIR] = 0
np.fill_diagonal(npmi, 1)
npmi = pd.DataFrame(np.nan_to_num(npmi), index=top, columns=top)

Zc = linkage(npmi.values, method='average', metric='correlation')
fam = fcluster(Zc, N_FAMILIES, criterion='maxclust')
families = pd.DataFrame({'tool': top, 'family': fam})
con.register('families', families)
con.sql(f"COPY families TO '{DATA}/tool_families.parquet' (FORMAT parquet)")
pairs = npmi.where(np.triu(np.ones(npmi.shape, dtype=bool), 1)).stack().rename('npmi').reset_index()
pairs.columns = ['a', 'b', 'npmi']
pairs['n_both'] = [co[idx[a], idx[b]] for a, b in zip(pairs.a, pairs.b)]
con.register('pairs', pairs)
con.sql(f"COPY pairs TO '{DATA}/tool_cousage.parquet' (FORMAT parquet)")

g = sns.clustermap(npmi, row_linkage=Zc, col_linkage=Zc, cmap='vlag', center=0, vmin=-0.3, vmax=0.6,
                   figsize=(18, 18), xticklabels=True, yticklabels=True, cbar_kws={'label': 'NPMI'})
g.figure.suptitle(f'Which of the {N_COUSAGE} most-used tools are used together (NPMI)', y=1.0)
g.savefig(f'{FIG}/cousage_clustermap.png', dpi=120, bbox_inches='tight')

# ---------------------------------------------------------------- report
print('Adoption shapes:')
print(trends.groupby('shape').agg(tools=('tool', 'size'), example=('tool', lambda s: ', '.join(s[:6]))).to_string())
print('\nFastest growing 2020-2025 (yearly change in share), among tools with 300+ articles:')
print(trends[trends.n_articles >= 300].nlargest(12, 'growth_2020_2025')[['tool', 'n_articles', 'growth_2020_2025', 'shape']].to_string(index=False))
print('\nFastest declining:')
print(trends[trends.n_articles >= 300].nsmallest(12, 'growth_2020_2025')[['tool', 'n_articles', 'growth_2020_2025', 'shape']].to_string(index=False))
print('\nCo-usage families:')
for f, grp in families.groupby('family'):
    print(f, ', '.join(grp.tool))
print(f'\nWritten {DATA}/tool_trends.parquet, tool_year.parquet, tool_families.parquet, tool_cousage.parquet')
print(f'Figures in {FIG}/')
