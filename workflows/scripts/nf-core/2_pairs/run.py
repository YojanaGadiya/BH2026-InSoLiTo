#!/usr/bin/env python3
"""run.py - software pairs from a folder of downloaded Nextflow pipelines, in a few steps.

Expected layout of the catalogue:  <catalogue>/<pipeline>/<version>/main.nf   (what scripts/1_download/downloader.py writes)

  python run.py setup                      # installs pyyaml and networkx (once)
  python run.py survey                     # what is in the catalogue: which pipelines can be analysed; changes nothing
  python run.py analyze [--limit 3]        # software pairs of every pipeline -> RESULTS/pipelines/<pipeline>/software_pairs.csv (a few minutes; --limit for a test)
  python run.py ids [--lookup]             # tool identifiers -> RESULTS/tool_identifiers.csv (--lookup also asks bio.tools, needs internet)
  python run.py match [--dump FILE]        # fill missing bio.tools IDs from a local bio.tools dump (biotools_raw*.jsonl in RESULTS, or --dump; the file can stay where it is)
  python run.py idtable                    # RESULTS/tools_with_ids.csv: every tool with at least one identifier
  python run.py merge [--pool-ids biotools]  # ONE table of all pairs over all pipelines -> RESULTS/pairs_merged_unique.csv, control_pool.csv
  python run.py tables [--dump FILE]       # the final tables (pairs, edges, nodes) -> nf-core/final tables (runs ../3_nodes_edges/build_tables.py)
  python run.py all [--lookup]             # analyze, ids, match, idtable, merge and tables one after the other

Options: --catalogue DIR (default: current folder)   --out DIR (default: nf-core/data, next to scripts/ and final tables/)
         --jobs N (pipelines analysed in parallel, default 4)   --timeout SECONDS (per pipeline, default 900)   --only TEXT (only pipelines whose path contains TEXT)
Nothing from the pipelines is executed and nothing is downloaded (except bio.tools answers with --lookup). Only text files are read.
"""
import argparse, collections, concurrent.futures as cf, csv, glob, json, os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..'))   # 2_pairs -> nf-core -> scripts -> workflows -> root
DATA = os.path.join(ROOT, 'workflows', 'data', 'nf-core')
BIOTOOLS_DATA = os.path.join(ROOT, 'directed-graph', 'biotools', 'data')
SKIP_PARENTS = ('modules', 'subworkflows', 'tests', 'test', '.git', 'work')   # a main.nf inside these is a module, not a pipeline
PRUNE = ('modules', 'subworkflows', 'workflows', 'conf', 'tests')

def libs_ok():
    try: import yaml, networkx; return True
    except ImportError: return False


def sh(*cmd):
    print('>', ' '.join(os.path.basename(c) if os.path.exists(c) else c for c in cmd), flush=True)
    return subprocess.run(cmd).returncode


def name_of(root, path):
    return os.path.relpath(path, root).replace('\\', '/').strip('./').replace('/', '__') or os.path.basename(path)


def find_pipelines(root):
    """-> (runnable, skipped): lists of (path, kind) and (path, reason). A pipeline is a folder with a main.nf (not inside modules/ ...)."""
    from pairs import pipeline_kind
    found, skipped = [], []
    for dp, dn, fn in os.walk(root):
        dn[:] = sorted(d for d in dn if d not in ('.git', 'work', 'node_modules'))
        rel = os.path.relpath(dp, root).replace('\\', '/').split('/')
        if rel == ['.']: continue
        if any(p in SKIP_PARENTS for p in rel): dn[:] = []; continue
        if 'main.nf' in fn:
            kind, why = pipeline_kind(dp)
            (found.append((dp, kind)) if kind else skipped.append((dp, why)))
            dn[:] = [d for d in dn if d not in PRUNE]
        elif len(rel) == 2 and not any(dp.startswith(f[0] + os.sep) for f in found): skipped.append((dp, 'no main.nf' + (' (only a download log)' if any(f.endswith('.log') for f in fn) else '')))
    return sorted(found), sorted(skipped)


def survey(a):
    found, skipped = find_pipelines(a.cat)
    c = collections.Counter(k for _, k in found)
    print(f'catalogue: {a.cat}\n{len(found)} pipelines can be analysed (modern DSL2: {c["dsl2"]}, old-style DSL1: {c["dsl1"]}); {len(skipped)} folders cannot:')
    for r, n in collections.Counter(w for _, w in skipped).most_common(): print(f'  {n:4d}  {r}')
    for p, w in skipped[:15]: print('   ', os.path.relpath(p, a.cat), '->', w)


def analyze_one(a, repo):
    name = name_of(a.cat, repo); out = os.path.join(a.out, 'pipelines', name); os.makedirs(out, exist_ok=True)
    t0 = time.time(); row = {'pipeline': name, 'path': repo, 'status': '', 'seconds': '', 'process_instances': '', 'process_edges': '', 'software_pairs': '', 'error': ''}
    try:
        r = subprocess.run([sys.executable, os.path.join(HERE, 'pairs.py'), '--repo', repo, '--out', out], capture_output=True, text=True, timeout=a.timeout)
        open(os.path.join(out, 'pairs.log'), 'w', encoding='utf-8').write(r.stdout + '\n' + r.stderr)
        if r.returncode != 0: row['status'] = 'FAILED'; row['error'] = (r.stderr.strip().splitlines() or r.stdout.strip().splitlines() or ['?'])[-1][:300]
        else:
            s = json.load(open(os.path.join(out, 'summary.json')))
            row.update(status=s['status'], process_instances=s['process_instances'], process_edges=sum(s['process_edges'].values()), software_pairs=sum(s['software_pairs'].values()),
                       error='; '.join(s['warnings'])[:300] or (('unresolved statements: %d' % s['unresolved_statements']) if s['unresolved_statements'] else ''))
    except subprocess.TimeoutExpired: row['status'] = 'TIMEOUT'
    except Exception as ex: row['status'] = 'ERROR'; row['error'] = str(ex)[:300]
    row['seconds'] = round(time.time() - t0, 1)
    return row


def analyze(a):
    if not libs_ok(): sys.exit('pyyaml / networkx missing - run: python run.py setup')
    pipes = [p for p, _ in find_pipelines(a.cat)[0]]
    if a.only: pipes = [p for p in pipes if a.only in p]
    if a.limit: pipes = pipes[:a.limit]
    print(len(pipes), 'pipelines under', a.cat); os.makedirs(a.out, exist_ok=True); rows = []
    with cf.ThreadPoolExecutor(max_workers=a.jobs) as ex:
        for i, row in enumerate(ex.map(lambda p: analyze_one(a, p), pipes), 1):
            rows.append(row); print(f"[{i}/{len(pipes)}] {row['status']:8} {row['pipeline']}  edges={row['process_edges']} {row['error'][:80]}", flush=True)
    with open(os.path.join(a.out, 'batch_summary.csv'), 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else ['pipeline'], lineterminator='\n'); w.writeheader(); w.writerows(rows)
    print('summary:', os.path.join(a.out, 'batch_summary.csv'), '|', dict(collections.Counter(r['status'] for r in rows)))


def ids(a):
    return sh(sys.executable, os.path.join(HERE, 'tool_ids.py'), 'collect', '--catalogue', a.cat, '--results', a.out, '--out', a.out, *(['--lookup'] if a.lookup else []))


def match(a):
    f = a.dump or next(iter(sorted(glob.glob(os.path.join(a.out, 'biotools_raw*.jsonl')), key=os.path.getmtime, reverse=True)), None)
    if not f or not os.path.exists(f): print('no bio.tools dump found - pass --dump FILE or put biotools_raw*.jsonl into ' + a.out + ' (step skipped)'); return 1
    return sh(sys.executable, os.path.join(HERE, 'tool_ids.py'), 'match', '--tools', os.path.join(a.out, 'tool_identifiers.csv'), '--jsonl', f, '--out', a.out)


def idtable(a):
    return sh(sys.executable, os.path.join(HERE, 'tool_ids.py'), 'table', '--dir', a.out)


def merge(a):
    return sh(sys.executable, os.path.join(HERE, 'merge_pairs.py'), '--results', a.out, '--pool-ids', a.pool_ids)


def tables(a):
    cmd = [sys.executable, os.path.join(HERE, '..', '3_nodes_edges', 'build_tables.py'), '--results', a.out] + (['--dump', a.dump] if a.dump else [])
    return sh(*cmd)


def setup(a):
    if not libs_ok(): sh(sys.executable, '-m', 'pip', 'install', '--user', '--quiet', 'pyyaml', 'networkx')
    print('libraries:', 'OK' if libs_ok() else 'STILL MISSING - run: pip install --user pyyaml networkx')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('step', choices=['setup', 'survey', 'analyze', 'ids', 'match', 'idtable', 'merge', 'tables', 'all'])
    ap.add_argument('--catalogue', default=os.path.join(DATA, 'download')); ap.add_argument('--out'); ap.add_argument('--jobs', type=int, default=4); ap.add_argument('--timeout', type=int, default=900)
    ap.add_argument('--limit', type=int, default=0); ap.add_argument('--only', default=''); ap.add_argument('--lookup', action='store_true'); ap.add_argument('--dump')
    ap.add_argument('--pool-ids', choices=['any', 'biotools'], default='any', help='merge: control pool needs both tools to have any identifier (any) or a bio.tools ID (biotools)')
    a = ap.parse_args()
    a.cat = os.path.abspath(a.catalogue)
    a.out = os.path.abspath(a.out or os.path.join(HERE, '..', '..', 'data'))   # nf-core/data, next to scripts/ and final tables/
    steps = {'setup': [setup], 'survey': [survey], 'analyze': [analyze], 'ids': [ids], 'match': [match], 'idtable': [idtable], 'merge': [merge], 'tables': [tables], 'all': [analyze, ids, match, idtable, merge, tables]}[a.step]
    for f in steps:
        rc = f(a)
        if rc and a.step != 'all': sys.exit(rc)


if __name__ == '__main__':
    main()
