"""Step 6: are tools used more by their own group than by other labs?

For every research article that uses a top tool in its methods (step 1), with authors
from OpenAlex (step 3):
  self       the article shares an author (OpenAlex author id) with any of the tool's papers
  same_org   no shared author, but a shared institution (OpenAlex institution id)
  cross_lab  neither
Articles OpenAlex doesn't know, or that have no authors, are left out.

Big institutions publish a lot, so a shared institution can happen by chance. Each rate
is compared with the chance level: the share of all articles in the same pool (articles
using any top tool) that share an author or an institution with that tool's team. The
ratio is the "lift" (1 = no more than chance).

Also: whether in-group use is higher in a tool's first years than later, and whether
more prominent tools are less in-group.

Writes patterns/data/ingroup_by_tool.parquet and ingroup_by_age.parquet.
Needs steps 1-3 (the articles step fully fetched). Run from BH2026-InSoLiTo/:
    python directed-graph/patterns/06_ingroup.py
"""

from collections import defaultdict

import duckdb
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

DATA = 'directed-graph/patterns/data'
MIN_KNOWN = 50   # uses with known authors a tool needs for its own rates

con = duckdb.connect()
tools = con.sql(f"SELECT tool, rank, n_articles FROM '{DATA}/tools.parquet'").df()
use = con.sql(f"SELECT pmcid, tool, year FROM '{DATA}/tool_use.parquet'").df()
# Read the fetch cache directly, so a partly finished fetch can already be looked at.
ARTICLES = f"read_json('{DATA}/openalex/articles.jsonl', format='newline_delimited')"
arts = con.sql(f"SELECT pmcid, author_ids, institution_ids FROM {ARTICLES} WHERE found AND len(author_ids) > 0").df()
papers = con.sql(f"""
    SELECT p.tool, CAST(o.year AS INTEGER) AS year, o.cited_by_count, o.authors
    FROM '{DATA}/tool_papers.parquet' p JOIN '{DATA}/oa_tool_papers.parquet' o USING (work_id) WHERE o.found
""").df()

fetched = con.sql(f"SELECT count(*) FROM {ARTICLES}").fetchone()[0]
wanted = use.pmcid.nunique()
print(f'Articles fetched from OpenAlex: {fetched:,} of {wanted:,} ({fetched / wanted:.0%}); '
      f'with authors: {len(arts):,}')

# The tool's team: every author and institution on any of its papers.
team_authors, team_insts = defaultdict(set), defaultdict(set)
for r in papers.itertuples():
    for a in r.authors:
        if a['author_id']:
            team_authors[r.tool].add(a['author_id'])
        team_insts[r.tool].update(a['institution_ids'])
release = papers.groupby('tool').year.min()
oa_cites = papers.groupby('tool').cited_by_count.sum()

art_authors = dict(zip(arts.pmcid, map(set, arts.author_ids)))
art_insts = dict(zip(arts.pmcid, map(set, arts.institution_ids)))
pool = list(art_authors)

# Chance level: share of the pool sharing an author / an institution with the team.
by_author, by_inst = defaultdict(set), defaultdict(set)
for p in pool:
    for a in art_authors[p]:
        by_author[a].add(p)
    for i in art_insts[p]:
        by_inst[i].add(p)


def chance(tool):
    sa = set().union(*(by_author.get(a, ()) for a in team_authors[tool])) if team_authors[tool] else set()
    si = set().union(*(by_inst.get(i, ()) for i in team_insts[tool])) if team_insts[tool] else set()
    return len(sa) / len(pool), len(si - sa) / len(pool)


use = use[use.pmcid.isin(art_authors)].copy()
use['kind'] = [('self' if art_authors[p] & team_authors[t]
                else 'same_org' if art_insts[p] & team_insts[t] else 'cross_lab')
               for p, t in zip(use.pmcid, use.tool)]
use['age'] = use.year - use.tool.map(release)

per_tool = use.groupby('tool').kind.value_counts().unstack(fill_value=0).reindex(columns=['self', 'same_org', 'cross_lab'], fill_value=0)
per_tool['known'] = per_tool.sum(axis=1)
for k in ('self', 'same_org', 'cross_lab'):
    per_tool[f'{k}_rate'] = per_tool[k] / per_tool.known
ch = pd.DataFrame([chance(t) for t in per_tool.index], index=per_tool.index, columns=['self_chance', 'same_org_chance'])
per_tool = per_tool.join(ch)
per_tool['self_lift'] = per_tool.self_rate / per_tool.self_chance.replace(0, np.nan)
per_tool['same_org_lift'] = per_tool.same_org_rate / per_tool.same_org_chance.replace(0, np.nan)
per_tool = per_tool.join(tools.set_index('tool')).join(oa_cites.rename('oa_citations')).join(release.rename('release_year'))
per_tool = per_tool.reset_index()

age_bins = pd.cut(use.age, [-100, 2, 5, 10, 100], labels=['0-2 years', '3-5 years', '6-10 years', '11+ years'])
by_age = use.groupby(age_bins, observed=True).kind.value_counts(normalize=True).unstack().reset_index().rename(columns={'age': 'tool_age'})

for name, df in (('ingroup_by_tool', per_tool), ('ingroup_by_age', by_age.astype({'tool_age': str}))):
    con.register('df', df)
    con.sql(f"COPY df TO '{DATA}/{name}.parquet' (FORMAT parquet)")

total = per_tool[['self', 'same_org', 'cross_lab']].sum()
print(f'\nAll {int(total.sum()):,} uses with known authors: self {total.self / total.sum():.1%}, '
      f'same organization {total.same_org / total.sum():.1%}, cross-lab {total.cross_lab / total.sum():.1%}')
ok = per_tool[per_tool.known >= MIN_KNOWN]
print(f'Tools with at least {MIN_KNOWN} known uses: {len(ok)}')
print(f'  median self rate {ok.self_rate.median():.1%} (chance {ok.self_chance.median():.2%}, median lift {ok.self_lift.median():.1f}x)')
print(f'  median same-org rate {ok.same_org_rate.median():.1%} (chance {ok.same_org_chance.median():.1%}, median lift {ok.same_org_lift.median():.1f}x)')
for col in ('self_rate', 'same_org_rate'):
    r = spearmanr(np.log(ok.n_articles), ok[col])
    print(f'  {col} vs prominence (log uses): rho={r.statistic:+.2f}, p={r.pvalue:.3f}')
print('\nBy tool age at the citing article:')
print(by_age.round(3).to_string(index=False))
print('\nMost in-group tools (self + same org), at least 50 known uses:')
ok = ok.assign(ingroup=ok.self_rate + ok.same_org_rate)
print(ok.nlargest(12, 'ingroup')[['tool', 'known', 'self_rate', 'same_org_rate', 'self_lift', 'same_org_lift']].round(3).to_string(index=False))
print('\nMost-used tools:')
print(ok.nsmallest(12, 'rank')[['tool', 'known', 'self_rate', 'same_org_rate', 'self_lift', 'same_org_lift']].round(3).to_string(index=False))
