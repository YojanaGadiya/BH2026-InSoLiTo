"""Step 4: do update papers boost a tool's use, and does growing use lead to updates?

Use = citations per year to all of a tool's papers in OpenAlex (the whole literature,
not only our corpus), from step 3. Citations grow every year for everyone, so each
tool's yearly growth is compared with the median growth of the top tools that year
("excess growth"; 0 = growing like a typical top tool).

An update event = a tool paper published after the tool's first paper (a new version or
an update article), in 2015-2023 so that three years before and two after are observed
(2026 is left out as a partial year). Several update papers in one year count once.

  1. Event study: average excess growth from 3 years before to 2 years after an update,
     with 95% bootstrap intervals.
  2. Expansion or takeover: after an update, how much of the tool's citations go to the
     new paper, and whether its older papers keep growing.
  3. Does growth lead to updates: the chance of an update paper in a year, by how fast
     the tool grew in the two years before.

Writes patterns/data/update_events.parquet, update_event_study.parquet, update_hazard.parquet
and patterns/figures/update_event_study.png. Run from BH2026-InSoLiTo/:
    python directed-graph/patterns/04_updates_vs_usage.py
"""

import json

import duckdb
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

DATA = 'directed-graph/patterns/data'
FIG = 'directed-graph/patterns/figures'
YEARS = range(2012, 2026)
EVENT_YEARS = range(2015, 2024)
WINDOW = range(-3, 3)
MIN_CITES = 20       # a tool needs this many citations the year before to count
N_BOOT = 2000
rng = np.random.default_rng(17)

con = duckdb.connect()
papers = con.sql(f"""
    SELECT p.tool, o.openalex_id, CAST(o.year AS INTEGER) AS pub_year, o.counts_by_year
    FROM '{DATA}/tool_papers.parquet' p JOIN '{DATA}/oa_tool_papers.parquet' o USING (work_id)
    WHERE o.found
""").df().drop_duplicates(['tool', 'openalex_id'])

rows = [(r.tool, r.openalex_id, r.pub_year, c['year'], c['cited_by_count'])
        for r in papers.itertuples() for c in json.loads(r.counts_by_year or '[]') if c['year'] in YEARS]
cites = pd.DataFrame(rows, columns=['tool', 'paper', 'pub_year', 'year', 'cites'])
T = cites.pivot_table(index='tool', columns='year', values='cites', aggfunc='sum', fill_value=0).reindex(columns=list(YEARS), fill_value=0)

g = np.log(T + 1).diff(axis=1)                       # yearly log growth
eligible = T.shift(1, axis=1) >= MIN_CITES
market = g.where(eligible).median()                   # typical top-tool growth each year
excess = g.sub(market, axis=1).where(eligible)

first = papers.groupby('tool').pub_year.min()
ev = papers.merge(first.rename('first'), on='tool')
ev = ev[(ev.pub_year > ev['first']) & ev.pub_year.isin(EVENT_YEARS)][['tool', 'pub_year']].drop_duplicates()
ev = ev[[bool(eligible.loc[t, y]) if t in eligible.index else False for t, y in zip(ev.tool, ev.pub_year)]]

# ---------------------------------------------------------------- 1. event study
grid = pd.DataFrame([{'tool': t, 'event_year': y, 'h': h, 'excess': excess.loc[t, y + h] if (y + h) in excess.columns else np.nan}
                     for t, y in zip(ev.tool, ev.pub_year) for h in WINDOW])


def boot_mean(x):
    x = np.asarray(x[~np.isnan(x)])
    if len(x) < 5:
        return np.nan, np.nan, np.nan, len(x)
    b = rng.choice(x, (N_BOOT, len(x))).mean(axis=1)
    return x.mean(), np.percentile(b, 2.5), np.percentile(b, 97.5), len(x)


study = pd.DataFrame([(h, *boot_mean(grid[grid.h == h].excess.values)) for h in WINDOW],
                     columns=['h', 'mean_excess', 'lo', 'hi', 'n'])

# Within-tool jump: each event's growth in the update year and the next, minus the same
# tool's growth in the 3 years before (removes "tools that update were growing anyway").
wide = grid.pivot_table(index=['tool', 'event_year'], columns='h', values='excess')
jump = (wide[[0, 1]].mean(axis=1) - wide[[-3, -2, -1]].mean(axis=1)).dropna()
jump_mean, jump_lo, jump_hi, jump_n = boot_mean(jump.values)
jump_share_up = (jump > 0).mean()

# Placebo: tool-years at least 3 years away from any update.
near = {(t, y + d) for t, y in zip(ev.tool, ev.pub_year) for d in range(-3, 4)}
placebo = np.array([excess.loc[t, y] for t in excess.index for y in excess.columns
                    if (t, y) not in near and not np.isnan(excess.loc[t, y])])
pl_mean, pl_lo, pl_hi, pl_n = boot_mean(placebo)

# ---------------------------------------------------------------- 2. expansion or takeover
new_share, old_growth = [], []
for t, y in zip(ev.tool, ev.pub_year):
    c = cites[cites.tool == t]
    tot = c.groupby('year').cites.sum()
    new = c[c.pub_year == y].groupby('year').cites.sum()
    old = c[c.pub_year < y].groupby('year').cites.sum()
    for h in (0, 1, 2):
        if y + h in tot.index and tot[y + h] > 0:
            new_share.append((h, new.get(y + h, 0) / tot[y + h]))
    if all(k in old.index for k in (y - 1, y + 2)) and old[y - 1] > 0:
        # old papers' growth over the 3 years after, minus the typical top tool's
        mk = market.loc[y:y + 2].sum()
        old_growth.append(np.log(old[y + 2] + 1) - np.log(old[y - 1] + 1) - mk)
share_by_h = pd.DataFrame(new_share, columns=['h', 'share']).groupby('h').share.median()

# ---------------------------------------------------------------- 3. does growth lead to updates?
haz = []
upd = set(zip(ev.tool, ev.pub_year))
for t in excess.index:
    for y in EVENT_YEARS:
        prior = excess.loc[t, [y - 2, y - 1]]
        if prior.isna().any() or first.get(t, 9999) >= y:
            continue
        haz.append({'tool': t, 'year': y, 'prior_excess': prior.mean(), 'update': (t, y) in upd})
haz = pd.DataFrame(haz)
haz['prior_growth'] = pd.qcut(haz.prior_excess, 4, labels=['slowest 25%', '25-50%', '50-75%', 'fastest 25%'])
from scipy.stats import chi2_contingency  # noqa: E402
hazard = haz.groupby('prior_growth', observed=True).agg(tool_years=('update', 'size'), updates=('update', 'sum'),
                                                        rate=('update', 'mean')).reset_index()
chi2_p = chi2_contingency(np.column_stack([hazard.updates, hazard.tool_years - hazard.updates]))[1]

# ---------------------------------------------------------------- save + report
con.register('ev', ev.rename(columns={'pub_year': 'event_year'}))
con.sql(f"COPY ev TO '{DATA}/update_events.parquet' (FORMAT parquet)")
study['placebo_mean'], study['placebo_lo'], study['placebo_hi'] = pl_mean, pl_lo, pl_hi
con.register('study', study)
con.sql(f"COPY study TO '{DATA}/update_event_study.parquet' (FORMAT parquet)")
hz = hazard.astype({'prior_growth': str})
con.register('hz', hz)
con.sql(f"COPY hz TO '{DATA}/update_hazard.parquet' (FORMAT parquet)")

fig, ax = plt.subplots(figsize=(7, 4))
ax.axhspan(pl_lo, pl_hi, color='0.85', label='years away from any update (95% CI)')
ax.axhline(0, color='0.4', lw=0.8)
ax.axvline(-0.5, color='0.6', ls=':', lw=0.8)
ax.errorbar(study.h, study.mean_excess, yerr=[study.mean_excess - study.lo, study.hi - study.mean_excess],
            fmt='o-', color='#2f6db5', capsize=3, label=f'around {len(ev)} update papers')
ax.set_xlabel('years from the update paper (0 = year it is published)')
ax.set_ylabel('excess citation growth (log)')
ax.set_title('Tool citation growth around an update paper, vs a typical top tool')
ax.legend(frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(f'{FIG}/update_event_study.png', dpi=150)

print(f'{len(ev)} update events in {ev.tool.nunique()} tools ({EVENT_YEARS.start}-{EVENT_YEARS.stop - 1})')
print('\nExcess growth by year from the update (log points; 0.10 ~ 10% faster than a typical top tool):')
print(study.round(3).to_string(index=False))
print(f'\nYears away from any update: mean {pl_mean:.3f} [{pl_lo:.3f}, {pl_hi:.3f}], n={pl_n}')
print(f'Within-tool jump (years 0-1 minus years -3..-1): {jump_mean:.3f} [{jump_lo:.3f}, {jump_hi:.3f}], '
      f'n={jump_n}; {jump_share_up:.0%} of updates are followed by faster growth than the same tool had before')
print('\nShare of the tool\'s citations going to the new paper (median):', share_by_h.round(2).to_dict())
og = np.array(old_growth)
print(f'Older papers after an update: 3-year excess growth median {np.median(og):.2f}, '
      f'{(og < 0).mean():.0%} grow slower than a typical top tool (n={len(og)})')
print('\nChance of an update paper in a year, by growth over the two years before:')
print(hazard.to_string(index=False))
print(f'Chi-squared test, update rate vs prior growth quartile: p = {chi2_p:.3f}')
