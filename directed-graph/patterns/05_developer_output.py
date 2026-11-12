"""Step 5: the people behind the tools.

1. Developer output: papers per year of each tool's first author (usually the developer)
   and last author (usually the lab head), from 5 years before to 5 years after the
   tool's first paper, scaled to the same person's average in the 5 years before. Careers
   grow anyway, so the curves are compared between tools that became heavily used (top
   third by citations) and the rest.
2. Who builds well-used tools: tool use against the team behind the first paper:
   team size, the lab head's papers before the tool (seniority), whether the lab head
   has other top tools, and the country and institution. Only the top tools are covered,
   so this compares well-used tools with very well-used ones.

Use = citations of all of the tool's papers in OpenAlex, and research articles using the
tool in methods in our corpus.

Writes patterns/data/developer_trajectories.parquet, tool_teams.parquet and
patterns/figures/developer_output.png. Needs steps 1-3. Run from BH2026-InSoLiTo/:
    python directed-graph/patterns/05_developer_output.py
"""

import duckdb
import matplotlib
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

DATA = 'directed-graph/patterns/data'
FIG = 'directed-graph/patterns/figures'
WINDOW = range(-5, 6)
rng = np.random.default_rng(23)

con = duckdb.connect()
tools = con.sql(f"SELECT tool, rank, n_articles FROM '{DATA}/tools.parquet'").df()
papers = con.sql(f"""
    SELECT p.tool, o.openalex_id, CAST(o.year AS INTEGER) AS year, o.cited_by_count, o.authors
    FROM '{DATA}/tool_papers.parquet' p JOIN '{DATA}/oa_tool_papers.parquet' o USING (work_id)
    WHERE o.found
""").df().drop_duplicates(['tool', 'openalex_id'])
years = con.sql(f"SELECT * FROM '{DATA}/oa_author_years.parquet'").df()
authors = con.sql(f"SELECT * FROM '{DATA}/oa_authors.parquet'").df().set_index('author_id')

cites = papers.groupby('tool').cited_by_count.sum().rename('oa_citations')
first = papers.sort_values(['tool', 'year', 'cited_by_count'], ascending=[True, True, False]).groupby('tool').head(1)
first = first.set_index('tool')
tools = tools.join(cites, on='tool').join(first[['year']].rename(columns={'year': 'release_year'}), on='tool')
tools['success'] = pd.qcut(tools.oa_citations.rank(method='first'), 3, labels=['lower third', 'middle third', 'top third'])

per_year = years.pivot_table(index='author_id', columns='year', values='n_works', fill_value=0)


def role_authors(auths, role):
    if role == 'first':
        return [a['author_id'] for a in auths if a['position'] == 'first' and a['author_id']]
    return [a['author_id'] for a in auths if a['position'] == 'last' and a['author_id']]


# ---------------------------------------------------------------- 1. trajectories
traj = []
for t, row in first.iterrows():
    R = int(row.year)
    for role in ('first', 'last'):
        for aid in role_authors(row.authors, role):
            if aid not in per_year.index:
                continue
            series = per_year.loc[aid]
            pre = np.mean([series.get(R + h, 0) for h in range(-5, 0)])
            if pre < 1:   # fewer than 1 paper a year before: too sparse to scale
                continue
            for h in WINDOW:
                traj.append({'tool': t, 'role': role, 'author_id': aid, 'release_year': R, 'h': h,
                             'relative_output': series.get(R + h, 0) / pre})
traj = pd.DataFrame(traj)
traj = traj[traj.release_year + 5 <= 2025]  # full 5 years after observed
traj = traj.merge(tools[['tool', 'success']], on='tool')


def ci(x):
    x = np.asarray(x)
    b = rng.choice(x, (2000, len(x))).mean(axis=1)
    return x.mean(), np.percentile(b, 2.5), np.percentile(b, 97.5), len(x)


summary = (traj.groupby(['role', 'success', 'h'], observed=True).relative_output
           .apply(lambda s: pd.Series(ci(s.values), index=['mean', 'lo', 'hi', 'n'])).unstack().reset_index())
after = traj[traj.h.between(1, 5)].groupby(['role', 'success', 'author_id', 'tool'], observed=True).relative_output.mean().reset_index()
tests = {}
for role in ('first', 'last'):
    a = after[(after.role == role) & (after.success == 'top third')].relative_output
    b = after[(after.role == role) & (after.success != 'top third')].relative_output
    tests[role] = (a.median(), b.median(), mannwhitneyu(a, b).pvalue, len(a), len(b))

# ---------------------------------------------------------------- 2. teams
lab_heads = {t: (role_authors(r.authors, 'last') or [None])[0] for t, r in first.iterrows()}
head_tools = pd.Series(lab_heads).dropna().value_counts()
team = []
for t, r in first.iterrows():
    R = int(r.year)
    head = lab_heads[t]
    prior = per_year.loc[head, [y for y in per_year.columns if y < R]].sum() if head in per_year.index else np.nan
    insts = [i for a in r.authors for i in a['institution_ids']]
    countries = sorted({c for a in r.authors for c in a['countries']})
    head_auth = next((a for a in r.authors if a['author_id'] == head), None)
    team.append({'tool': t, 'release_year': R, 'team_size': len(r.authors), 'lab_head': head,
                 'lab_head_name': authors.name.get(head) if head in authors.index else None,
                 'lab_head_prior_papers': prior,
                 'lab_head_top_tools': int(head_tools.get(head, 0)) if head else 0,
                 'lab_head_institutions': head_auth['institution_ids'] if head_auth else [],
                 'institutions': sorted(set(insts)), 'countries': countries,
                 'n_countries': len(countries)})
team = pd.DataFrame(team).merge(tools.drop(columns='release_year'), on='tool')
team['multi_tool_lab'] = team.lab_head_top_tools >= 2
# Older tools have had more years to collect citations: compare citations per year since release.
team['citations_per_year'] = team.oa_citations / (2026 - team.release_year).clip(lower=1)

corr = {}
for col in ('team_size', 'lab_head_prior_papers', 'n_countries', 'release_year'):
    d = team[[col, 'citations_per_year', 'n_articles']].dropna()
    r1 = spearmanr(d[col], d.citations_per_year)
    r2 = spearmanr(d[col], d.n_articles)
    corr[col] = (r1.statistic, r1.pvalue, r2.statistic, r2.pvalue, len(d))
mt = team.groupby('multi_tool_lab').agg(tools=('tool', 'size'), median_citations_per_year=('citations_per_year', 'median'),
                                        median_articles=('n_articles', 'median'))
mt_p = mannwhitneyu(team[team.multi_tool_lab].citations_per_year, team[~team.multi_tool_lab].citations_per_year).pvalue
mt_p2 = mannwhitneyu(team[team.multi_tool_lab].n_articles, team[~team.multi_tool_lab].n_articles).pvalue
heads = (team[team.multi_tool_lab].groupby(['lab_head', 'lab_head_name'])
         .agg(tools=('tool', lambda s: ', '.join(sorted(s))), n=('tool', 'size'), citations=('oa_citations', 'sum'))
         .sort_values('n', ascending=False).reset_index())
countries = team.explode('countries').groupby('countries').agg(tools=('tool', 'nunique'), citations=('oa_citations', 'sum')).sort_values('tools', ascending=False)

for name, df in (('developer_trajectories', summary), ('tool_teams', team)):
    con.register('df', df)
    con.sql(f"COPY df TO '{DATA}/{name}.parquet' (FORMAT parquet)")

fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
colors = {'lower third': '#b9c2c8', 'middle third': '#7d8c96', 'top third': '#2f6db5'}
for ax, role in zip(axes, ('first', 'last')):
    for succ, grp in summary[summary.role == role].groupby('success', observed=True):
        ax.plot(grp.h, grp['mean'], 'o-', color=colors[succ], label=f'{succ} of tools by citations')
        ax.fill_between(grp.h, grp.lo, grp.hi, color=colors[succ], alpha=0.15)
    ax.axvline(0, color='0.6', ls=':')
    ax.axhline(1, color='0.4', lw=0.8)
    ax.set_title(f'{role.capitalize()} author of the tool paper')
    ax.set_xlabel('years from the tool\'s first paper')
axes[0].set_ylabel('papers per year, vs own average before')
axes[0].legend(frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(f'{FIG}/developer_output.png', dpi=150)

print(f'Trajectories: {traj.author_id.nunique()} authors of {traj.tool.nunique()} tools (released by 2020)')
print(summary[summary.h.isin([-1, 0, 2, 5])].round(2).to_string(index=False))
for role, (a, b, p, na, nb) in tests.items():
    print(f'{role} author, years +1..+5: median relative output {a:.2f} (top-third tools, n={na}) vs {b:.2f} (others, n={nb}), p={p:.3f}')
print('\nSpearman correlation with OpenAlex citations per year since release, and with corpus use:')
for k, (r1, p1, r2, p2, n) in corr.items():
    print(f'  {k:24s} citations/yr rho={r1:+.2f} p={p1:.3f}   corpus use rho={r2:+.2f} p={p2:.3f}  n={n}')
print('\nLab heads with 2+ top tools vs others:')
print(mt.to_string(), f'\n  Mann-Whitney p: citations/yr {mt_p:.3f}, corpus use {mt_p2:.3f}')
print(heads.head(12).to_string(index=False))
print('\nCountries on the first paper (tools, summed citations):')
print(countries.head(10).to_string())
