#!/usr/bin/env python3
"""build_tables.py - the final tables (pairs, edges, nodes) from the merged pair table.

  python build_tables.py [--results DIR] [--out DIR] [--dump FILE]

Reads  (from --results, default nf-core/data):
  tools_with_ids.csv        from `run.py idtable`
  pairs_merged_unique.csv   from `run.py merge` (must contain direction_type)
  the bio.tools dump        biotools_raw*.jsonl: found in --results, or give its place with --dump (it can stay where it is, it is only read)
Writes to --out (default nf-core/final tables):
  biotools_pairs.csv   one row per WORKFLOW and tool pair (the same pair appears once for every workflow where it is valid)
  biotools_edges.csv   one row per tool pair, with the list of workflows
  biotools_nodes.csv   one row per tool
and to --results:
  pairs_biotools.csv   one row per tool pair with all columns (the table the other three are made from; not a final table)

Only pairs where BOTH tools have a bio.tools ID and that are valid in at least one workflow are kept.
direction_type: one_direction (A feeds B, never B feeds A) | both_directions | parallel (side by side, no path between them).
Needs pandas.
"""
import argparse, glob, os, sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))


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
        'source_software', 'target_software', 'source_biotools_id', 'target_biotools_id', 'source_match_method', 'target_match_method', 'direction_type',
        'n_pipelines_positive', 'n_pipelines_negative', 'n_process_edges_positive', 'reverse_direction_pipelines', 'pipelines_positive', 'pipelines_negative',
    ]].rename(columns={'source_software': 'source', 'target_software': 'target', 'source_biotools_id': 'source_biotool_id', 'target_biotools_id': 'target_biotool_id',
                       'source_match_method': 'source_match', 'target_match_method': 'target_match'}).reset_index(drop=True)
    out[['source_match', 'target_match']] = out[['source_match', 'target_match']].replace('meta_id', 'biotool_id')   # an ID written in the module's meta.yml
    print(f'All pairs: {len(pairs):,} | both tools have a bio.tools ID: {both.sum():,} | and valid in at least one workflow: {len(out):,}')
    print(out['direction_type'].value_counts().to_string())
    return out


def edge_table(pairs):
    return pairs[['source_biotool_id', 'target_biotool_id', 'direction_type', 'n_pipelines_positive', 'pipelines_positive']].rename(columns={
        'source_biotool_id': 'source', 'target_biotool_id': 'target', 'n_pipelines_positive': 'n_workflows', 'pipelines_positive': 'workflow_ids'})


def workflow_pair_table(pairs):
    """one row per workflow and tool pair; source_biotools_name / target_biotools_name are the tool names used in the pipelines"""
    return (pairs.assign(workflow_id=pairs['pipelines_positive'].str.split(';')).explode('workflow_id')
            [['workflow_id', 'source_biotool_id', 'target_biotool_id', 'direction_type', 'source', 'target']]
            .rename(columns={'source_biotool_id': 'source', 'target_biotool_id': 'target', 'source': 'source_biotools_name', 'target': 'target_biotools_name'})
            .reset_index(drop=True))


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
    ap.add_argument('--results', default=os.path.join(HERE, '..', '..', 'data')); ap.add_argument('--out', default=os.path.join(HERE, '..', '..', 'final tables')); ap.add_argument('--dump')
    a = ap.parse_args()
    results, out = os.path.abspath(a.results), os.path.abspath(a.out)
    tools_file, pairs_file, dump = os.path.join(results, 'tools_with_ids.csv'), os.path.join(results, 'pairs_merged_unique.csv'), find_dump(results, a.dump)
    missing = [f for f in (tools_file, pairs_file) if not os.path.exists(f)] + ([] if dump and os.path.exists(dump) else ['bio.tools dump (biotools_raw*.jsonl: put it into the data folder or pass --dump FILE)'])
    for f in (tools_file, pairs_file, dump): print(('found    ' if f and os.path.exists(f) else 'MISSING  ') + str(f))
    if missing: sys.exit('missing: ' + '; '.join(missing))

    pairs = pair_table(tools_file, pairs_file)
    edges, workflow_pairs, nodes = edge_table(pairs), workflow_pair_table(pairs), node_table(pairs, dump)
    assert len(edges) == len(pairs) and len(workflow_pairs) == pairs['n_pipelines_positive'].sum()
    assert set(nodes['biotoolsID']) == set(edges['source']) | set(edges['target'])
    assert edges['direction_type'].notna().all() and workflow_pairs['direction_type'].notna().all()

    os.makedirs(out, exist_ok=True)
    pairs.to_csv(os.path.join(results, 'pairs_biotools.csv'), index=False)
    for name, df in (('biotools_edges.csv', edges), ('biotools_pairs.csv', workflow_pairs), ('biotools_nodes.csv', nodes)): df.to_csv(os.path.join(out, name), index=False)
    print(f'\nnodes {len(nodes):,} | edges {len(edges):,} | workflow pairs {len(workflow_pairs):,}\nwritten to {out} (pairs_biotools.csv to {results})')


if __name__ == '__main__':
    main()
