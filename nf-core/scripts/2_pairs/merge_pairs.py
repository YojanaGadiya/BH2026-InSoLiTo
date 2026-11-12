#!/usr/bin/env python3
"""merge_pairs.py - ONE table of directed software pairs from many pipelines, with the number of pipelines each pair occurs in.

  python merge_pairs.py --results RESULTS_DIR [--out OUT_DIR] [--pool-ids any|biotools]

RESULTS_DIR is the folder written by `run.py analyze` (RESULTS_DIR/pipelines/<pipeline>/software_pairs.csv). If tool_identifiers.csv (from
`run.py ids` / `match`) and tools_with_ids.csv (from `run.py idtable`) are there, tools are merged by bio.tools ID and the IDs are added.

pair_type of a pair in one pipeline (rank = strength of the evidence)
  1 data_flow        an output of A is a data input of B
  2 gate_metadata    B depends on A only through a key, a verdict or metadata (valid, weaker)
  3 parallel_co_use  A and B read the same upstream output, no path either way, same run; valid in BOTH directions (a row for each direction)
  4 rejected         the graph says A cannot feed B (never valid)
A pair is VALID when it is 1, 2 or 3 in at least one pipeline.

direction_type of a valid pair
  one_direction    A -> B is a sequence (type 1 or 2) and B -> A never is
  both_directions  A -> B and B -> A are both sequences (in the same or in different pipelines)
  parallel         A and B only run side by side (type 3 only)

Writes
  pairs_merged_unique.csv   one row per distinct directed pair: in how many pipelines it is valid / rejected, which, direction_type, bio.tools IDs
  control_pool.csv          the pool: valid pairs where BOTH tools have an identifier (see --pool-ids); a short list of columns
  pairs_merged_summary.md   counts and the most frequent pairs

Tool names are compared in a normalised form (lower case, text in brackets and non-letters removed: 'SAMtools' = 'samtools'); tools with the same
bio.tools ID are one tool, and `star_align` is `star` when `star` has an ID.
"""
import argparse, collections, csv, os, re, sys

PAIR_TYPES = [(1, 'data_flow'), (2, 'gate_metadata'), (3, 'parallel_co_use'), (4, 'rejected')]
RANK = {t: r for r, t in PAIR_TYPES}
SEQUENCE = {'data_flow', 'gate_metadata'}   # A is used before B (as opposed to side by side)
POSITIVE = SEQUENCE | {'parallel_co_use'}   # a VALID pair
UTILITY = {'multiqc', 'gunzip', 'gzip', 'cat', 'untar', 'tar', 'xz', 'unzip', 'bzip2', 'pigz', 'gawk', 'sed', 'python', 'rbase', 'coreutils', 'wget', 'curl', 'custom', 'gzipgunzip', 'untarunzip'}

BTID = {}      # name key -> bio.tools ID
PARENT = {}    # name key -> key of the parent tool, for names that are a tool plus an operation (star_align -> star)
METHOD = {}    # bio.tools ID (lower case) -> (rank, method, detail): how the ID was found, strongest evidence over all names that carry it
ANYID = set()  # name keys with ANY identifier (bio.tools ID, bioconda package, DOI, other ID)


def rawkey(name):
    n = re.sub(r'\(.*?\)', '', name or '').lower()
    k = re.sub(r'[^a-z0-9]', '', n) or re.sub(r'[^a-z0-9]', '', (name or '').lower())
    return PARENT.get(k, k)


def key(name):
    k = rawkey(name)
    return ('bt:' + BTID[k].lower()) if k in BTID else k


def method(k, i=1):
    """match method (i=1) or detail (i=2) of the bio.tools ID behind a pair key `bt:<id>`"""
    return METHOD.get(k[3:], (0, '', ''))[i] if k.startswith('bt:') else ''


def load_ids(results):
    f = os.path.join(results, 'tool_identifiers.csv')
    if not os.path.exists(f): return False
    for r in csv.DictReader(open(f, encoding='utf-8')):
        if r['biotools_id']:
            BTID[r['tool_key']] = r['biotools_id']
            src = r.get('id_source', '') or ''
            if src == 'meta.yml': rk, m, dt = 1, 'meta_id', 'identifier in the nf-core module meta.yml'
            elif src.startswith('biotools-jsonl-'):
                kinds = src[len('biotools-jsonl-'):].split('+'); dt = 'bio.tools dump: ' + '+'.join(kinds)
                rk, m = (2, 'doi') if 'doi' in kinds else ((3, 'other_id') if set(kinds) & {'pmid', 'pmcid', 'rrid'} else ((4, 'url') if 'url' in kinds else (5, 'name')))
            else: rk, m, dt = 5, 'name', 'bio.tools API exact ' + src.replace('bio.tools-exact-', '')
            bi = r['biotools_id'].lower()
            if bi not in METHOD or rk < METHOD[bi][0]: METHOD[bi] = (rk, m, dt)
        if r.get('parent_confidence') == 'operation-suffix' and r.get('parent_key'): PARENT[r['tool_key']] = r['parent_key']
        if any(r.get(c) for c in ('biotools_id', 'bioconda_package', 'doi', 'publications', 'other_ids')): ANYID.add(r['tool_key'])
    f = os.path.join(results, 'tools_with_ids.csv')
    if os.path.exists(f):   # the table of tools with at least one ID is the final word on which tools count as identified
        ANYID.clear(); ANYID.update(r['tool_key'] for r in csv.DictReader(open(f, encoding='utf-8')))
        print('tool list:', len(ANYID), 'tools with at least one identifier (tools_with_ids.csv)')
    print('tool IDs loaded:', len(BTID))
    return True


def write_csv(path, rows, cols=None):
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=cols or (list(rows[0].keys()) if rows else ['x']), extrasaction='ignore', lineterminator='\n'); w.writeheader(); w.writerows(rows)


def read_pipelines(results):
    """-> list of (pipeline, source, target, pair_type, n_process_edges); a parallel pair is added in both directions"""
    rows, n = [], 0
    pdir = os.path.join(results, 'pipelines')
    for p in sorted(os.listdir(pdir)) if os.path.isdir(pdir) else []:
        f = os.path.join(pdir, p, 'software_pairs.csv')
        if not os.path.exists(f): continue
        n += 1
        for r in csv.DictReader(open(f, encoding='utf-8')):
            t = r['pair_type']
            if t not in RANK: continue
            rows.append((p, r['source_software'], r['target_software'], t, int(r['n_process_edges'] or 0)))
            if t == 'parallel_co_use': rows.append((p, r['target_software'], r['source_software'], t, int(r['n_process_edges'] or 0)))
    return rows, n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--results', required=True); ap.add_argument('--out')
    ap.add_argument('--pool-ids', choices=['any', 'biotools'], default='any', help='control pool: both tools need at least one identifier of any kind (any) or a bio.tools ID (biotools)')
    a = ap.parse_args()
    out = a.out or a.results; os.makedirs(out, exist_ok=True)
    have_ids = load_ids(a.results)
    raw, n_pipes = read_pipelines(a.results)

    long_rows = []
    for pipe, src, dst, t, n_edges in raw:
        ks, kd = key(src), key(dst)
        if ks == kd: continue   # a tool feeding itself (also: two names of the same tool)
        long_rows.append({'pipeline': pipe, 'source_software': src, 'target_software': dst, 'source_key': ks, 'target_key': kd,
                          'source_biotools_id': ks[3:] if ks.startswith('bt:') else '', 'target_biotools_id': kd[3:] if kd.startswith('bt:') else '',
                          'source_match_method': method(ks), 'target_match_method': method(kd), 'source_match_detail': method(ks, 2), 'target_match_detail': method(kd, 2),
                          'source_has_any_id': ks.startswith('bt:') or rawkey(src) in ANYID, 'target_has_any_id': kd.startswith('bt:') or rawkey(dst) in ANYID,
                          'pair_type': t, 'pair_type_rank': RANK[t], 'n_process_edges': n_edges,
                          'involves_utility_or_reporting': ks.replace('bt:', '') in UTILITY or kd.replace('bt:', '') in UTILITY})

    groups = collections.defaultdict(list)
    for r in long_rows: groups[(r['source_key'], r['target_key'])].append(r)
    names = collections.defaultdict(collections.Counter)
    for r in long_rows: names[r['source_key']][r['source_software']] += 1; names[r['target_key']][r['target_software']] += 1
    shown = lambda k: names[k].most_common(1)[0][0]
    pipes_of = lambda rs, types: sorted({r['pipeline'] for r in rs if r['pair_type'] in types})
    uniq = []
    for (ks, kd), rs in groups.items():
        pos, neg = pipes_of(rs, POSITIVE), pipes_of(rs, {'rejected'})
        rev = pipes_of(groups.get((kd, ks), []), POSITIVE)
        seq, rev_seq = pipes_of(rs, SEQUENCE), pipes_of(groups.get((kd, ks), []), SEQUENCE)
        direction = ('both_directions' if rev_seq else 'one_direction') if seq else ('parallel' if pos else '')   # the reverse of a parallel pair is only its mirror, it does not make a sequence
        best = min(rs, key=lambda r: r['pair_type_rank'])['pair_type']
        uniq.append({'source_software': shown(ks), 'target_software': shown(kd), 'source_key': ks, 'target_key': kd,
                     'source_biotools_id': rs[0]['source_biotools_id'], 'target_biotools_id': rs[0]['target_biotools_id'],
                     'source_match_method': rs[0]['source_match_method'], 'target_match_method': rs[0]['target_match_method'],
                     'source_match_detail': rs[0]['source_match_detail'], 'target_match_detail': rs[0]['target_match_detail'],
                     'source_has_any_id': rs[0]['source_has_any_id'], 'target_has_any_id': rs[0]['target_has_any_id'],
                     'n_pipelines_positive': len(pos), 'n_pipelines_data_flow': len(pipes_of(rs, {'data_flow'})), 'n_pipelines_metadata_gate': len(pipes_of(rs, {'gate_metadata'})),
                     'n_pipelines_parallel': len(pipes_of(rs, {'parallel_co_use'})), 'n_pipelines_negative': len(neg),
                     'n_process_edges_positive': sum(r['n_process_edges'] for r in rs if r['pair_type'] in POSITIVE),
                     'best_pair_type': best, 'best_pair_type_rank': RANK[best], 'direction_type': direction,
                     'reverse_direction_pipelines': len(rev), 'reverse_sequence_pipelines': len(rev_seq),
                     'involves_utility_or_reporting': ks.replace('bt:', '') in UTILITY or kd.replace('bt:', '') in UTILITY,
                     'pipelines_positive': ';'.join(pos), 'pipelines_negative': ';'.join(neg)})
    uniq.sort(key=lambda r: (-r['n_pipelines_positive'], r['best_pair_type_rank'], r['source_key'], r['target_key']))
    write_csv(os.path.join(out, 'pairs_merged_unique.csv'), uniq)

    if have_ids:   # control pool: valid in >= 1 pipeline AND both tools have an identifier; A->B and B->A are separate rows
        ok_s = (lambda u: bool(u['source_biotools_id'])) if a.pool_ids == 'biotools' else (lambda u: bool(u['source_has_any_id']))
        ok_t = (lambda u: bool(u['target_biotools_id'])) if a.pool_ids == 'biotools' else (lambda u: bool(u['target_has_any_id']))
        valid = [u for u in uniq if u['n_pipelines_positive'] >= 1]
        pool = [u for u in valid if ok_s(u) and ok_t(u)]
        slim = ['source_software', 'target_software', 'source_key', 'target_key', 'source_biotools_id', 'target_biotools_id', 'source_match_method', 'target_match_method', 'n_pipelines_positive', 'direction_type', 'reverse_direction_pipelines']
        write_csv(os.path.join(out, 'control_pool.csv'), pool, slim)
        left_out = [u for u in valid if not (ok_s(u) and ok_t(u))]
        print('control pool:', len(pool), 'directed pairs | best_pair_type', dict(collections.Counter(u['best_pair_type'] for u in pool)), '| direction_type', dict(collections.Counter(u['direction_type'] for u in pool)))
        print('pool rule:', 'both tools have a bio.tools ID' if a.pool_ids == 'biotools' else 'both tools have at least one identifier', '| valid pairs left out because a tool has no ID:', len(left_out))
    else: print('NOTE: no tool_identifiers.csv in the results folder, so no control pool (run `run.py ids` first)')

    valid = [u for u in uniq if u['n_pipelines_positive']]
    core = [u for u in valid if not u['involves_utility_or_reporting']]
    multi = [u for u in valid if u['n_pipelines_positive'] >= 2]
    lines = ['# Merged software pairs', '', f'- pipelines merged: {n_pipes}', f'- distinct directed valid pairs: {len(valid)} ({len(core)} without utilities/MultiQC)', f'- valid pairs found in 2+ pipelines: {len(multi)}',
             f"- valid pairs by direction_type: {dict(collections.Counter(u['direction_type'] for u in valid))}", f"- pairs by best pair_type: {dict(collections.Counter(u['best_pair_type'] for u in uniq))}", '', '## pair_type', '']
    lines += [f'- `{t}` (rank {r})' for r, t in PAIR_TYPES]
    lines += ['', '## Most frequent valid pairs (analysis tools only)', '', '| source | target | pipelines | process edges | direction | best pair_type |', '|---|---|---|---|---|---|']
    for u in core[:40]: lines.append(f"| {u['source_software']} | {u['target_software']} | {u['n_pipelines_positive']} | {u['n_process_edges_positive']} | {u['direction_type']} | {u['best_pair_type']} |")
    open(os.path.join(out, 'pairs_merged_summary.md'), 'w', encoding='utf-8', newline='\n').write('\n'.join(lines) + '\n')
    print(f'{n_pipes} pipelines merged: {len(long_rows)} rows, {len(uniq)} distinct pairs, {len(valid)} valid, {len(multi)} valid in 2+ pipelines')
    print('most frequent (analysis tools only):')
    for u in core[:15]: print(f"  {u['n_pipelines_positive']:4d} pipelines  {u['source_software']} -> {u['target_software']}  [{u['direction_type']}, {u['best_pair_type']}]")
    print('written:', os.path.join(out, 'pairs_merged_unique.csv'))


if __name__ == '__main__':
    main()
