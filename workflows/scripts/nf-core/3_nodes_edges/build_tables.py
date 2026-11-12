#!/usr/bin/env python3
"""build_tables.py - the final tables (pairs, edges, nodes) from the merged pair table.

  python build_tables.py [--results DIR] [--out DIR] [--dump FILE]

Reads  (from --results, default workflows/data/nf-core):
  tools_with_ids.csv        from `run.py idtable`
  pairs_merged_unique.csv   from `run.py merge` (must contain direction_type)
  the bio.tools dump        biotools_raw*.jsonl: found in --results, or give its place with --dump (it can stay where it is, it is only read)
Writes to --out (default: the same folder, workflows/data/nf-core):
  biotools_pairs.csv   one row per WORKFLOW and directed tool pair: workflow_id, source, target, parallel (`parallel` = the two tools run side by side
                       in that workflow, `not parallel` = one comes after the other). The flag is exact for every workflow
  biotools_edges.csv   source, target, n_workflows, workflow_ids; one row per direction, each with its own count (A->B and B->A are two rows)
  biotools_nodes.csv   one row per tool
and to --results:
  pairs_biotools.csv   one row per tool pair with all columns (the table the other three are made from; not a final table)

Only pairs where BOTH tools have a bio.tools ID and that are valid in at least one workflow are kept.
direction_type: one_direction (A feeds B, never B feeds A) | both_directions | parallel (side by side, no path between them).
Needs pandas.
"""
import argparse, collections, glob, os, sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '2_pairs'))
import merge_pairs as mp   # the same tool-name matching as the merge, so the pairs line up


def find_dump(results, dump):
    if dump: return os.path.expanduser(dump)
    found = sorted(glob.glob(os.path.join(results, 'biotools_raw*.jsonl')), key=os.path.getmtime, reverse=True)
    return found[0] if found else None


def pair_table(tools_file, pairs_file):
    """tool pairs where both tools have a bio.tools ID and the pair is valid in at least one workflow"""
    tools = pd.read_csv(tools_file)
    has_id = tools['biotools_id'].notna() & tools['biotools_id'].astype(str).str.strip().ne('')
    with_id = set(tools.loc[has_id, 'tool'].dropna().astype(str).str.strip())
    print(f'Tools in the table:        {len(tools):,}\n  with a bio.tools ID:     {has_id.sum():,}\n  without a bio.tools ID:  {(~has_id).sum():,}')

    pairs = pd.read_csv(pairs_file)
    if 'direction_type' not in pairs.columns: sys.exit('no direction_type column in ' + pairs_file + ': run `run.py merge` with the current scripts')
    both = pairs['source_software'].astype(str).str.strip().isin(with_id) & pairs['target_software'].astype(str).str.strip().isin(with_id)
    valid = pairs['n_pipelines_positive'] > 0
    out = pairs.loc[both & valid, [
        'source_key', 'target_key', 'source_software', 'target_software', 'source_biotools_id', 'target_biotools_id', 'source_match_method', 'target_match_method', 'direction_type',
        'n_pipelines_positive', 'n_pipelines_negative', 'n_process_edges_positive', 'reverse_direction_pipelines', 'pipelines_positive', 'pipelines_negative',
    ]].rename(columns={'source_software': 'source', 'target_software': 'target', 'source_biotools_id': 'source_biotool_id', 'target_biotools_id': 'target_biotool_id',
                       'source_match_method': 'source_match', 'target_match_method': 'target_match'}).reset_index(drop=True)
    out[['source_match', 'target_match']] = out[['source_match', 'target_match']].replace('meta_id', 'biotool_id')   # an ID written in the module's meta.yml
    print(f'All pairs: {len(pairs):,} | both tools have a bio.tools ID: {both.sum():,} | and valid in at least one workflow: {len(out):,}')
    print(out['direction_type'].value_counts().to_string())
    return out


def edge_table(pairs):
    """one row per DIRECTED tool pair. n_workflows / workflow_ids: the workflows where this direction (source -> target) appears.
    A -> B and B -> A are two rows with their own counts; a parallel pair appears once forward and once reverse in each workflow, so both rows have the same count"""
    return pairs[['source_biotool_id', 'target_biotool_id', 'n_pipelines_positive', 'pipelines_positive']].rename(columns={
        'source_biotool_id': 'source', 'target_biotool_id': 'target', 'n_pipelines_positive': 'n_workflows', 'pipelines_positive': 'workflow_ids'})


def workflow_types(results):
    """(workflow, source key, target key) -> the pair types found in that workflow, read from results/pipelines/<workflow>/software_pairs.csv"""
    mp.load_ids(results)
    types = collections.defaultdict(set)
    rows, _ = mp.read_pipelines(results)
    for pipe, src, dst, t, _n in rows:
        ks, kd = mp.key(src), mp.key(dst)
        if ks != kd: types[(pipe, ks, kd)].add(t)
    return types


def workflow_pair_table(pairs, types):
    """one row per workflow and directed pair. parallel = the pair is only side by side in that workflow; a sequence (data flow or gate / metadata) makes it `not parallel`"""
    rows, missing = [], 0
    for r in pairs.itertuples():
        for wf in r.pipelines_positive.split(';'):
            t = types.get((wf, r.source_key, r.target_key))
            if not t: missing += 1; continue
            rows.append({'workflow_id': wf, 'source': r.source_biotool_id, 'target': r.target_biotool_id, 'parallel': 'not parallel' if t & mp.SEQUENCE else 'parallel'})
    if missing: sys.exit(f'{missing} workflow-pair rows have no pair type in results/pipelines: run `run.py analyze` and `run.py merge` with the current scripts first')
    return pd.DataFrame(rows)


def node_table(pairs, dump):
    ends = pd.concat([pairs[['source_biotool_id', 'source_match', 'pipelines_positive']].set_axis(['biotoolsID', 'matching_method', 'pipelines_positive'], axis=1),
                      pairs[['target_biotool_id', 'target_match', 'pipelines_positive']].set_axis(['biotoolsID', 'matching_method', 'pipelines_positive'], axis=1)])
    ends['workflow_id'] = ends['pipelines_positive'].str.split(';')
    nodes = ends.explode('workflow_id').groupby('biotoolsID').agg(n_workflows=('workflow_id', 'nunique'), matching_method=('matching_method', 'first')).reset_index()
    raw = pd.read_json(dump, lines=True)
    names = (raw[['biotoolsID', 'name']].dropna().assign(key=lambda d: d['biotoolsID'].astype(str).str.strip().str.casefold()).drop_duplicates('key').set_index('key')['name'])
    nodes['biotools_name'] = nodes['biotoolsID'].astype(str).str.strip().str.casefold().map(names)
    print(f'bio.tools records: {len(raw):,} | nodes: {len(nodes):,} | names found: {nodes["biotools_name"].notna().sum():,} | missing: {nodes["biotools_name"].isna().sum():,}')
    if nodes['biotools_name'].isna().any(): print('  no name for: ' + ', '.join(nodes.loc[nodes['biotools_name'].isna(), 'biotoolsID']))
    return nodes[['biotoolsID', 'biotools_name', 'n_workflows', 'matching_method']]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    DATA = os.path.join(HERE, '..', '..', '..', 'data', 'nf-core')
    ap.add_argument('--results', default=DATA); ap.add_argument('--out', default=DATA); ap.add_argument('--dump')
    a = ap.parse_args()
    results, out = os.path.abspath(a.results), os.path.abspath(a.out)
    tools_file, pairs_file, dump = os.path.join(results, 'tools_with_ids.csv'), os.path.join(results, 'pairs_merged_unique.csv'), find_dump(results, a.dump)
    missing = [f for f in (tools_file, pairs_file) if not os.path.exists(f)] + ([] if dump and os.path.exists(dump) else ['bio.tools dump (biotools_raw*.jsonl: put it into the data folder (workflows/data/nf-core) or pass --dump FILE)'])
    for f in (tools_file, pairs_file, dump): print(('found    ' if f and os.path.exists(f) else 'MISSING  ') + str(f))
    if missing: sys.exit('missing: ' + '; '.join(missing))

    pairs = pair_table(tools_file, pairs_file)
    edges, workflow_pairs, nodes = edge_table(pairs), workflow_pair_table(pairs, workflow_types(results)), node_table(pairs, dump)
    assert len(edges) == len(pairs) and len(workflow_pairs) == pairs['n_pipelines_positive'].sum()
    assert set(nodes['biotoolsID']) == set(edges['source']) | set(edges['target'])
    assert set(workflow_pairs['parallel']) <= {'parallel', 'not parallel'}

    os.makedirs(out, exist_ok=True)
    pairs.to_csv(os.path.join(results, 'pairs_biotools.csv'), index=False)
    for name, df in (('biotools_edges.csv', edges), ('biotools_pairs.csv', workflow_pairs), ('biotools_nodes.csv', nodes)): df.to_csv(os.path.join(out, name), index=False)
    print(f'\nnodes {len(nodes):,} | edges {len(edges):,} | workflow pairs {len(workflow_pairs):,}\nwritten to {out} (pairs_biotools.csv to {results})')


if __name__ == '__main__':
    main()
