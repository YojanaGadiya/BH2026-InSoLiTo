"""Step 16: registered workflows against the all-publications literature graph.

Three workflow sources, each a list of bio.tools pairs "tool A feeds tool B":
  galaxy     WorkflowHub Galaxy workflows   workflowhub/data/biotools/galaxy/biotools_edges.csv
  nextflow   WorkflowHub Nextflow workflows, without the nf-core space (counted below)
                                            workflowhub/data/biotools/nextflow/biotools_edges_without_nfcore.csv
  nf-core    nf-core pipelines, from their code   nf-core/final tables/biotools_edges.csv
Mapping confidence (1 best, 3 worst): Galaxy and Nextflow give it per edge; for nf-core
it comes from the nodes' matching_method (biotool_id 1, doi or url 2, name 3; an edge
has the worse of its two tools). nf-core pairs that only run side by side
(direction_type = parallel) count as co-use, not as a direction.

Both sides are mapped to the merged tools of step 12. Galaxy wrapper entries that step 12
dropped because they own no paper (samtools_sort, bedtools_bamtobed) go to the tool owning
their papers; other sub-tool IDs to a tool named in them (gatk_haplotype_caller -> gatk).

A literature step A -> B counts as in a workflow source if A and B are connected directly
(either way) or B is downstream of A through other steps in one workflow: papers often cite
STAR then DESeq2, while workflows put featureCounts in between.

Literature: the step-15 pairs (10+ articles), with the step weight (sum of 1 / tool gap,
~ articles citing the two as neighbours) and adjacency.

Outputs (patterns/data/):
  wf_vs_lit_pairs.parquet    one row per tool pair seen in a workflow source or the literature
  wf_vs_lit_tools.parquet    per tool: workflows per source, literature links
  wf_lit_gaps.parquet        literature steps (step weight >= GAP_MIN_STEP) in no workflow source
  wf_only.parquet            workflow links (WF_ONLY_MIN+ workflows, confidence <= 2) papers hardly cite
  wf_area_coverage.parquet   per workflow area and source: share of literature steps found
  wf_chain_coverage.parquet  step-15 chains and how many of their steps each source has
Needs steps 12, 13 and 15. Run from BH2026-InSoLiTo/:
    python directed-graph/patterns/16_workflowhub_vs_literature.py
"""

import re
from collections import defaultdict
from pathlib import Path

import duckdb
import networkx as nx
import numpy as np
import pandas as pd

DATA = Path('directed-graph/patterns/data')
WF = Path('directed-graph/workflowhub/data/biotools')
NFCORE = Path('nf-core/final tables')
TOOL_WORKS = 'directed-graph/citation-order/data/tool_works.parquet'
EXPORT = 'directed-graph/europepmc/data/export'
SOURCES = ['galaxy', 'nextflow', 'nf-core']
GAP_MIN_STEP = 50
WF_ONLY_MIN = 5
LIT_MIN = 5          # a workflow pair "is in the literature" if seen in this many articles
NF_CONF = {'biotool_id': 1, 'doi': 2, 'url': 2, 'name': 3}
AREAS = ['Enrichment, networks & pathways', 'RNA-seq & ChIP-seq', 'Structural biology & simulation',
         'Variant calling & population genetics', 'Phylogenetics', 'Genome assembly & annotation',
         'Metagenome binning & annotation', 'Microbial genomes', 'Phylogenomics', 'Transcript assembly',
         'Microbiome amplicons', 'Protein domains & localisation', 'Sequence search', 'Mass spectrometry & imaging']

con = duckdb.connect()
con.sql("SET preserve_insertion_order = false")

# ------------------------------------------------------------------ entry -> tool
tools = con.sql(f"SELECT tool, entries, strength FROM '{DATA}/allpub_tools.parquet'").df()
e2t = {e: t for t, es in zip(tools.tool, tools.entries) for e in es}
strength_of = tools.set_index('tool').strength.to_dict()
po = con.sql(f"SELECT work_id, owners FROM '{DATA}/paper_owner.parquet'").df()
owners = {w: list(o) for w, o in zip(po.work_id, po.owners)}
cites = con.sql(f"SELECT work_id, n_citations FROM '{EXPORT}/works.parquet' WHERE work_id IN "
                f"(SELECT work_id FROM '{TOOL_WORKS}')").df().set_index('work_id').n_citations.to_dict()
lists = con.sql(f"SELECT biotoolsID, list(DISTINCT work_id) AS w FROM '{TOOL_WORKS}' GROUP BY 1").df()
for e, ws in zip(lists.biotoolsID, lists.w):
    if e in e2t:
        continue
    for w in sorted(ws, key=lambda x: -cites.get(x, 0)):
        o = [e2t[x] for x in owners.get(w, []) if x in e2t]
        if o:
            e2t[e] = o[0]
            break


# Copies merged in step 12 that mostly carry another entry's citations: a workflow step with that
# entry is the copy (bio.tools' hisat2 lists Bowtie and SpliceMap; HISAT is cited via quasr/splicemap).
cp = con.sql(f"SELECT tool, citations_mostly_of FROM '{DATA}/allpub_tools.parquet' WHERE citations_mostly_of IS NOT NULL").df()
WF_OVERRIDE = {o: t for t, o in zip(cp.tool, cp.citations_mostly_of) if o and ',' not in o and o not in t.split('/')}


def to_tool(e):
    if e in WF_OVERRIDE:
        return WF_OVERRIDE[e]
    if e in e2t:
        return e2t[e]
    parts = [x for x in re.split(r'[^A-Za-z0-9]+', e) if len(x) >= 4]
    cands = [e2t[x] for x in parts if x in e2t] + [e2t[x.lower()] for x in parts if x.lower() in e2t]
    t = max(cands, key=lambda t: strength_of.get(t, 0)) if cands else e
    e2t[e] = t
    return t


# ------------------------------------------------------------------ workflow sources
def load(src):
    if src == 'galaxy':
        d = pd.read_csv(WF / 'galaxy' / 'biotools_edges.csv', dtype=str)
        d['parallel'] = False
    elif src == 'nextflow':
        d = pd.read_csv(WF / 'nextflow' / 'biotools_edges_without_nfcore.csv', dtype=str)
        d['parallel'] = False
    else:
        d = pd.read_csv(NFCORE / 'biotools_edges.csv', dtype=str)
        nodes = pd.read_csv(NFCORE / 'biotools_nodes.csv', dtype=str)
        nc = dict(zip(nodes.biotoolsID, nodes.matching_method.map(NF_CONF).fillna(3)))
        d['confidence'] = [max(nc.get(a, 3), nc.get(b, 3)) for a, b in zip(d.source, d.target)]
        d['parallel'] = d.direction_type.eq('parallel')
    d['confidence'] = d.confidence.astype(int)
    d['s'] = d.source.map(to_tool)
    d['t'] = d.target.map(to_tool)
    d = d[d.s != d.t]
    rows = defaultdict(lambda: [set(), 9, True])  # (s, t) -> workflows, best confidence, only parallel
    for r in d.itertuples(index=False):
        x = rows[(r.s, r.t)]
        x[0].update(str(r.workflow_ids).split(';'))
        x[1] = min(x[1], r.confidence)
        x[2] = x[2] and bool(r.parallel)
    out = pd.DataFrame([{'s': s, 't': t, 'wf': w, 'conf': c, 'parallel': p} for (s, t), (w, c, p) in rows.items()])
    out['source'] = src
    return out


wf = pd.concat([load(s) for s in SOURCES], ignore_index=True)


def reachable(src):
    """(a, b) -> workflows where b is downstream of a (any number of steps), from the per-workflow pairs."""
    if src == 'galaxy':
        d = pd.read_csv(WF / 'galaxy' / 'biotools_pairs.csv', dtype=str)
    elif src == 'nextflow':
        d = pd.read_csv(WF / 'nextflow' / 'biotools_pairs_without_nfcore.csv', dtype=str)
    else:
        d = pd.read_csv(NFCORE / 'biotools_pairs.csv', dtype=str)
        d = d[d.direction_type != 'parallel']
    d['s'], d['t'] = d.source.map(to_tool), d.target.map(to_tool)
    out = defaultdict(set)
    for w, g in d[d.s != d.t].groupby('workflow_id'):
        G = nx.DiGraph(list(zip(g.s, g.t)))
        for a in G:
            for b in nx.descendants(G, a):
                if b != a:
                    out[(a, b)].add(w)
    return out


reach = {src: reachable(src) for src in SOURCES}
tool_wfs = defaultdict(lambda: defaultdict(set))
for r in wf.itertuples(index=False):
    tool_wfs[r.s][r.source] |= r.wf
    tool_wfs[r.t][r.source] |= r.wf

# ------------------------------------------------------------------ join with the literature
lit = con.sql(f"""SELECT source, target, n_total, share_forward, adjacency, step_fwd, step_total
                  FROM '{DATA}/allpub_steps.parquet'""").df()
key = lambda a, b: (a, b) if a < b else (b, a)
lit['k'] = [key(a, b) for a, b in zip(lit.source, lit.target)]
L = {k: r for k, r in zip(lit.k, lit.itertuples(index=False))}  # a dict: tuple keys confuse .loc
pairs = {}
for r in wf.itertuples(index=False):
    k = key(r.s, r.t)
    p = pairs.setdefault(k, {})
    for src in SOURCES:
        p.setdefault(src, [set(), set(), set(), 9])  # same way as literature arrow, reverse, parallel, conf
    lit_src = L[k].source if k in L else min(k)
    slot = 2 if r.parallel else (0 if r.s == lit_src else 1)
    p[r.source][slot] |= r.wf
    p[r.source][3] = min(p[r.source][3], r.conf)
rows = []
for k, p in pairs.items():
    a, b = (L[k].source, L[k].target) if k in L else k
    row = {'source': a, 'target': b}
    for src in SOURCES:
        f, rv, par, c = p[src]
        row[f'{src}_fwd'], row[f'{src}_rev'], row[f'{src}_parallel'] = len(f), len(rv - f), len(par - f - rv)
        row[f'{src}_conf'] = c if c < 9 else np.nan
    if k in L:
        x = L[k]
        row.update(lit_articles=x.n_total, lit_share_forward=x.share_forward, adjacency=x.adjacency, lit_step=x.step_fwd)
    else:
        row.update(lit_articles=0, lit_share_forward=np.nan, adjacency=np.nan, lit_step=0)
    rows.append(row)
P = pd.DataFrame(rows)
for src in SOURCES:
    P[f'{src}_any'] = P[f'{src}_fwd'] + P[f'{src}_rev'] + P[f'{src}_parallel']
P['wf_total'] = P[[f'{s}_any' for s in SOURCES]].sum(1)
P['best_conf'] = P[[f'{s}_conf' for s in SOURCES]].min(1)

comm = con.sql(f"SELECT tool, community FROM '{DATA}/allpub_hubs.parquet'").df().set_index('tool').community
in_src = {src: {t for t, v in tool_wfs.items() if v.get(src)} for src in SOURCES}
in_any = set().union(*in_src.values())
wf_keys = {src: set(P.loc[P[f'{src}_any'] > 0, ['source', 'target']].itertuples(index=False, name=None)) for src in SOURCES}
steps = lit[lit.step_fwd >= GAP_MIN_STEP].copy()
for src in SOURCES:
    steps[f'direct_{src}'] = [(a, b) in wf_keys[src] for a, b in zip(steps.source, steps.target)]
    # in a workflow: connected directly (either way) or B downstream of A through other steps
    steps[f'in_{src}'] = steps[f'direct_{src}'] | pd.Series([(a, b) in reach[src] for a, b in zip(steps.source, steps.target)], index=steps.index)
steps['direct_any'] = steps[[f'direct_{s}' for s in SOURCES]].any(axis=1)
steps['in_any'] = steps[[f'in_{s}' for s in SOURCES]].any(axis=1)
steps['both_tools_in_workflows'] = steps.source.isin(in_any) & steps.target.isin(in_any)
steps['area'] = steps.source.map(comm)
gaps = steps[~steps.in_any].sort_values('step_fwd', ascending=False)
area = []
for c in range(len(AREAS)):
    s = steps[steps.area == c]
    row = {'area': AREAS[c], 'lit_steps': len(s)}
    for src in SOURCES:
        row[src] = float(s[f'in_{src}'].mean()) if len(s) else np.nan
    row['any'] = float(s.in_any.mean()) if len(s) else np.nan
    area.append(row)
area = pd.DataFrame(area)
wf_only = P[(P.wf_total >= WF_ONLY_MIN) & (P.lit_articles < LIT_MIN) & (P.best_conf <= 2)].sort_values('wf_total', ascending=False)

tl = pd.DataFrame([{'tool': t, **{f'{s}_workflows': len(v.get(s, ())) for s in SOURCES}} for t, v in tool_wfs.items()])
tl['workflows'] = tl[[f'{s}_workflows' for s in SOURCES]].sum(1)
tl['lit_links'] = tl.tool.map(strength_of).fillna(0)
absent = tools[~tools.tool.isin(in_any)].nlargest(25, 'strength')[['tool', 'strength']]

chains = con.sql(f"SELECT * FROM '{DATA}/allpub_chains_steps.parquet'").df()
ch_rows = []
for c in chains.chain:
    st = c.split(' > ')
    pr = list(zip(st, st[1:]))
    hit = lambda p, src: p in wf_keys[src] or p in reach[src]
    ch_rows.append({'chain': c, 'steps': len(pr), **{src: sum(hit(p, src) for p in pr) for src in SOURCES},
                    'any': sum(any(hit(p, s) for s in SOURCES) for p in pr)})
ch_cov = pd.DataFrame(ch_rows)

for name, df in (('wf_vs_lit_pairs', P), ('wf_vs_lit_tools', tl), ('wf_lit_gaps', gaps), ('wf_only', wf_only),
                 ('wf_area_coverage', area), ('wf_chain_coverage', ch_cov), ('wf_lit_steps', steps)):
    con.register('df', df)
    con.sql(f"COPY df TO '{DATA / name}.parquet' (FORMAT parquet)")

# ------------------------------------------------------------------ report
pd.set_option('display.width', 230)
pd.set_option('display.max_colwidth', 90)
print('Sources:')
for src in SOURCES:
    s = P[P[f'{src}_any'] > 0]
    n_wf = len(set().union(*(v.get(src, set()) for v in tool_wfs.values())))
    directed = s[(s[f'{src}_fwd'] + s[f'{src}_rev']) > 0]
    good = directed[directed[f'{src}_conf'] <= 2]
    seen = (good.lit_articles >= LIT_MIN).mean()
    clear = good[(good.lit_articles >= 30) & (np.abs(good.lit_share_forward - 0.5) >= 0.1) & (good[f'{src}_fwd'] != good[f'{src}_rev'])]
    agree = ((clear[f'{src}_fwd'] > clear[f'{src}_rev']) == (clear.lit_share_forward > 0.5)).mean()
    print(f'  {src:8s} {n_wf:4d} workflows, {len(in_src[src]):4d} tools, {len(s):5d} pairs ({len(directed)} directed); '
          f'confidence <= 2: {len(good)} directed, {seen:.0%} seen in {LIT_MIN}+ articles; '
          f'direction agrees with the literature for {agree:.0%} of {len(clear)} clear pairs')
print(f'\nLiterature steps (step weight >= {GAP_MIN_STEP}): {len(steps):,}; directly connected in a workflow: {steps.direct_any.mean():.0%}; '
      f'connected or downstream: {steps.in_any.mean():.0%} '
      f'(galaxy {steps.in_galaxy.mean():.0%}, nextflow {steps.in_nextflow.mean():.0%}, nf-core {steps["in_nf-core"].mean():.0%})')
print(f'  missing from all, with both tools in some workflow: {(gaps.both_tools_in_workflows).sum():,}; with a tool in no workflow: {(~gaps.both_tools_in_workflows).sum():,}')
top100 = steps.nlargest(100, 'step_fwd')
print(f'  of the 100 strongest literature steps: {top100.direct_any.mean():.0%} directly connected, {top100.in_any.mean():.0%} connected or downstream')
both_in = steps[steps.both_tools_in_workflows]
print(f'  steps whose two tools are both in workflows: {len(both_in):,}; of these in a workflow: {both_in.in_any.mean():.0%}')
print('\nCoverage by workflow area (share of literature steps found):')
print(area.round(2).to_string(index=False))
print('\nStrongest literature steps in no workflow, both tools used in workflows:')
print(gaps[gaps.both_tools_in_workflows].head(20)[['source', 'target', 'n_total', 'step_fwd', 'adjacency', 'share_forward']].round(2).to_string(index=False))
print('\nStrongest literature steps with a tool in no workflow:')
print(gaps[~gaps.both_tools_in_workflows].head(12)[['source', 'target', 'n_total', 'step_fwd']].round(0).to_string(index=False))
print('\nMost-cited tools in no workflow:')
print(absent.head(20).to_string(index=False))
print('\nWorkflow links papers hardly cite (confidence <= 2):')
print(wf_only.head(20)[['source', 'target', 'galaxy_any', 'nextflow_any', 'nf-core_any', 'lit_articles']].to_string(index=False))
print('\nTool use, workflows vs literature (Spearman):', round(tl[['workflows', 'lit_links']].corr(method='spearman').iloc[0, 1], 2))
print(tl.nlargest(12, 'workflows')[['tool', 'galaxy_workflows', 'nextflow_workflows', 'nf-core_workflows', 'lit_links']].to_string(index=False))
print('\nLiterature chains and their steps in workflows:')
print(ch_cov.sort_values('steps', ascending=False).head(20).to_string(index=False))
