"""Step 16: registered Galaxy workflows (WorkflowHub) against the all-publications graph.

Workflow edges come from directed-graph/workflowhub/ (download_workflowhub.py,
galaxy_tool_pairs.py, map_galaxy_tools.py): one row per bio.tools pair, tool A feeds
tool B, with the workflows it appears in. Both sides are mapped to the same tools:
  - bio.tools entries to the merged tools of step 12;
  - Galaxy wrapper entries that step 12 dropped because they own no paper
    (samtools_sort, bedtools_bamtobed, ...) to the tool that owns their papers, so a
    workflow step "samtools_sort" counts as SAMtools;
  - other sub-tool entries to a tool named in their ID (gatk_haplotype_caller -> gatk,
    ensembl_genomes_bacteria_blast -> blast).
Each workflow pair keeps the weakest mapping method of its two ends (codex > name >
prefix); results are given for all pairs and for codex-only pairs.

Outputs (patterns/data/):
  wf_vs_lit_pairs.parquet   every tool pair seen in either, with workflows, literature
                            articles, step weight, adjacency and both directions
  wf_vs_lit_tools.parquet   per tool: workflows, literature links, ranks
  wf_lit_gaps.parquet       common literature steps not in any workflow, split into
                            "both tools in Galaxy" (a workflow nobody registered) and
                            "a tool not in Galaxy"
  wf_only.parquet           workflow links the literature hardly has
  wf_chain_coverage.parquet literature chains of step 15 and how much of each is in workflows

Needs steps 12, 13, 15 and the workflowhub outputs. Run from BH2026-InSoLiTo/:
    python directed-graph/patterns/16_workflowhub_vs_literature.py
    python directed-graph/patterns/16_workflowhub_vs_literature.py --wf-dir other/data   # elsewhere
"""

import argparse
import re
from collections import defaultdict
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

DATA = Path('directed-graph/patterns/data')
TOOL_WORKS = 'directed-graph/citation-order/data/tool_works.parquet'
EXPORT = 'directed-graph/europepmc/data/export'
METHOD_RANK = {'codex': 0, 'name': 1, 'prefix': 2}
GAP_MIN_STEP = 50    # literature step weight (~ articles citing the two as neighbours)
WF_ONLY_MIN = 5      # workflows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--wf-dir', default='directed-graph/workflowhub/data')
    args = ap.parse_args(argv)
    wf_dir = Path(args.wf_dir)
    edges_csv, map_csv = wf_dir / 'biotools' / 'biotools_edges.csv', wf_dir / 'mapping' / 'galaxy_tool_map.csv'
    if not edges_csv.exists():
        raise SystemExit(f'{edges_csv} not found: run the three workflowhub scripts first (see directed-graph/workflowhub/README)')

    con = duckdb.connect()
    con.sql("SET preserve_insertion_order = false")
    tools = con.sql(f"SELECT tool, entries, strength FROM '{DATA}/allpub_tools.parquet'").df()
    e2t = {e: t for t, es in zip(tools.tool, tools.entries) for e in es}
    # Entries step 12 dropped: send them to the owner of their most-cited listed paper.
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

    wf = pd.read_csv(edges_csv, dtype=str)
    # Still unmapped (sub-tools such as gatk_haplotype_caller, ensembl_genomes_bacteria_blast):
    # a name part that is itself a tool's entry, the most-linked such tool.
    strength_of = tools.set_index('tool').strength.to_dict()
    for e in set(wf.source) | set(wf.target):
        if e in e2t:
            continue
        parts = [x for x in re.split(r'[^A-Za-z0-9]+', e) if len(x) >= 4]
        cands = [e2t[x] for x in parts if x in e2t] + [e2t[x.lower()] for x in parts if x.lower() in e2t]
        if cands:
            e2t[e] = max(cands, key=lambda t: strength_of.get(t, 0))
    wf['n_workflows'] = wf.n_workflows.astype(int)
    method = {}
    if map_csv.exists():
        m = pd.read_csv(map_csv, dtype=str).dropna(subset=['biotoolsID'])
        for b, meth in zip(m.biotoolsID, m.method):
            method[b] = min(method.get(b, 9), METHOD_RANK.get(meth, 9))
    wf['s'] = wf.source.map(lambda x: e2t.get(x, x))
    wf['t'] = wf.target.map(lambda x: e2t.get(x, x))
    wf['mapped_lit'] = wf.source.isin(e2t) & wf.target.isin(e2t)
    wf['method'] = [max(method.get(a, 9), method.get(b, 9)) for a, b in zip(wf.source, wf.target)]
    wf = wf[wf.s != wf.t]
    # Workflows per tool pair, counting each workflow once.
    ids = defaultdict(set)
    meth = {}
    for r in wf.itertuples(index=False):
        k = (r.s, r.t)
        ids[k].update(str(r.workflow_ids).split(';'))
        meth[k] = min(meth.get(k, 9), r.method)
    wfp = pd.DataFrame([{'s': s, 't': t, 'wf_fwd': len(v), 'method': meth[(s, t)]} for (s, t), v in ids.items()])
    rev = wfp.set_index(['s', 't']).wf_fwd
    wfp['wf_rev'] = [rev.get((t, s), 0) for s, t in zip(wfp.s, wfp.t)]
    wfp = wfp[(wfp.wf_fwd > wfp.wf_rev) | ((wfp.wf_fwd == wfp.wf_rev) & (wfp.s < wfp.t))]
    tool_wfs = defaultdict(set)  # workflows each tool appears in
    for r in wf.itertuples(index=False):
        for x in str(r.workflow_ids).split(';'):
            tool_wfs[r.s].add(x)
            tool_wfs[r.t].add(x)

    lit = con.sql(f"SELECT source, target, n_total, share_forward, adjacency, step_fwd, step_rev, step_total "
                  f"FROM '{DATA}/allpub_steps.parquet'").df()
    lit_k = {(a, b): i for i, (a, b) in enumerate(zip(lit.source, lit.target))}
    rows = []
    seen = set()
    for r in wfp.itertuples(index=False):
        i = lit_k.get((r.s, r.t))
        flip = False
        if i is None:
            i, flip = lit_k.get((r.t, r.s)), True
        L = lit.iloc[i] if i is not None else None
        a, b = (r.s, r.t)
        lit_fwd = None
        if L is not None:
            seen.add((L.source, L.target))
            lit_fwd = (1 - L.share_forward) if flip else L.share_forward
        rows.append({'source': a, 'target': b, 'wf_fwd': r.wf_fwd, 'wf_rev': r.wf_rev, 'wf_total': r.wf_fwd + r.wf_rev,
                     'method': {0: 'codex', 1: 'name', 2: 'prefix'}.get(r.method, 'unmapped'),
                     'lit_articles': L.n_total if L is not None else 0, 'lit_share_same_way': lit_fwd,
                     'adjacency': L.adjacency if L is not None else np.nan,
                     'lit_step': L.step_total if L is not None else 0})
    both = pd.DataFrame(rows)
    in_galaxy = set(tool_wfs)
    lit_only = lit[[k not in seen for k in zip(lit.source, lit.target)]].copy()
    lit_only['both_in_galaxy'] = lit_only.source.isin(in_galaxy) & lit_only.target.isin(in_galaxy)
    gaps = lit_only[lit_only.step_fwd >= GAP_MIN_STEP].sort_values('step_fwd', ascending=False)
    wf_only = both[(both.wf_total >= WF_ONLY_MIN) & (both.lit_articles < 5)].sort_values('wf_total', ascending=False)

    # Tool level
    lit_strength = tools.set_index('tool').strength
    tl = pd.DataFrame({'tool': list(tool_wfs), 'workflows': [len(v) for v in tool_wfs.values()]})
    tl['lit_links'] = tl.tool.map(lit_strength).fillna(0)
    tl['wf_rank'] = tl.workflows.rank(ascending=False)
    tl['lit_rank'] = tl.lit_links.rank(ascending=False)

    # Chains of step 15
    chains_f = DATA / 'allpub_chains_steps.parquet'
    cov = None
    if chains_f.exists():
        ch = con.sql(f"SELECT * FROM '{chains_f}'").df()
        wf_set = {(a, b) for a, b in zip(both.source, both.target)} | {(b, a) for a, b in zip(both.source, both.target)}
        cov = ch.assign(steps_in_workflows=[sum((x, y) in wf_set for x, y in zip(c.split(' > '), c.split(' > ')[1:]))
                                            for c in ch.chain])
        cov['steps'] = cov.length - 1

    out = {'wf_vs_lit_pairs': both, 'wf_vs_lit_tools': tl, 'wf_lit_gaps': gaps, 'wf_only': wf_only}
    if cov is not None:
        out['wf_chain_coverage'] = cov
    for name, df in out.items():
        con.register('df', df)
        con.sql(f"COPY df TO '{DATA / name}.parquet' (FORMAT parquet)")

    pd.set_option('display.width', 220)
    n_wf = len(set().union(*tool_wfs.values())) if tool_wfs else 0
    print(f'Workflows: {n_wf}; tools in workflows: {len(tool_wfs)}; workflow tool pairs: {len(both)} '
          f'(codex-mapped {(both.method == "codex").mean():.0%})')
    for label, sub in (('all pairs', both), ('codex only', both[both.method == 'codex'])):
        inlit = sub.lit_articles >= 5
        agree = sub[inlit & sub.lit_share_same_way.notna() & (sub.wf_fwd > sub.wf_rev)]
        print(f'  {label}: {inlit.mean():.0%} of workflow pairs seen in 5+ articles; '
              f'direction agrees with the literature for {(agree.lit_share_same_way > 0.5).mean():.0%} of {len(agree)}')
    top = lit[lit.step_fwd >= GAP_MIN_STEP]
    print(f'Literature steps (step weight >= {GAP_MIN_STEP}): {len(top)}; in some workflow: '
          f'{1 - len(gaps) / max(len(top), 1):.0%}; missing with both tools in Galaxy: {gaps.both_in_galaxy.sum()}')
    print(f'Spearman, tool use in workflows vs literature: {tl[["wf_rank", "lit_rank"]].corr(method="spearman").iloc[0, 1]:.2f}')
    print('\nCommon literature steps with both tools in Galaxy but in no workflow:')
    print(gaps[gaps.both_in_galaxy].head(20)[['source', 'target', 'n_total', 'step_fwd', 'adjacency', 'share_forward']].round(2).to_string(index=False))
    print('\nCommon literature steps with a tool not in any workflow:')
    print(gaps[~gaps.both_in_galaxy].head(12)[['source', 'target', 'n_total', 'step_fwd']].round(0).to_string(index=False))
    print('\nWorkflow links the literature hardly has:')
    print(wf_only.head(15)[['source', 'target', 'wf_total', 'lit_articles', 'method']].to_string(index=False))
    print('\nUsed in workflows far more than cited:')
    print(tl[tl.workflows >= 10].assign(gap=lambda d: d.lit_rank - d.wf_rank).nlargest(10, 'gap')[['tool', 'workflows', 'lit_links']].to_string(index=False))
    print('Cited far more than used in workflows (in Galaxy at all):')
    print(tl[tl.lit_links >= 5000].assign(gap=lambda d: d.wf_rank - d.lit_rank).nlargest(10, 'gap')[['tool', 'workflows', 'lit_links']].to_string(index=False))
    if cov is not None:
        print('\nLiterature chains and their steps found in workflows:')
        print(cov.head(20)[['chain', 'steps', 'steps_in_workflows']].to_string(index=False))


if __name__ == '__main__':
    main()
