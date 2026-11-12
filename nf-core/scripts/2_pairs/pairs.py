#!/usr/bin/env python3
"""pairs.py - the directed software pairs of ONE Nextflow pipeline, by static analysis (nothing is run, nothing is downloaded).

  python pairs.py --repo PIPELINE_DIR --out OUT_DIR

PIPELINE_DIR is a folder with a main.nf (modern DSL2 or old-style DSL1). The same rules are used for every pipeline.

Steps
  1. graph      which process output goes into which process input (graph_dsl2.py, or graph_dsl1.py for old-style pipelines)
  2. software   each process -> one software name (module meta.yml tool, else the conda package of environment.yml, else the module name;
                old-style pipelines: container / conda directives and script commands)
  3. edges      every process edge is classified (direct, merged channel, gate / metadata, ...) and checked: edges that cannot carry data
                (a dropped tuple element, a multiMap split) are rejected
  4. pairs      process edges are collapsed to software pairs; two tools that read the same upstream output and can run in the same run
                but have no path between them are a PARALLEL pair; mutually exclusive alternatives (--aligner star / hisat2) are dropped

pair_type of a software pair
  data_flow        an output of A is a data input of B
  gate_metadata    B depends on A only through a key, a verdict or metadata (valid, weaker)
  parallel_co_use  A and B read the same upstream output, no path either way, same run (valid in BOTH directions, written once here)
  rejected         every process edge between them was rejected (graph says A cannot feed B)

Writes OUT_DIR/software_pairs.csv, process_edges.csv, process_software_map.csv, summary.json
Needs Python >= 3.10, pyyaml, networkx.
"""
import argparse, collections, csv, json, os, re, sys

import networkx as nx
import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import graph_dsl1, graph_dsl2
from graph_dsl2 import cond_atoms, infeasible

UNNAMED_WORKFLOW = re.compile(r'^\s*workflow\s*\{', re.M)
DSL1_PROCESS = re.compile(r'^[ \t]*process[ \t]+\w+[ \t]*\{', re.M)
HELPERS = {'htslib', 'coreutils', 'grep', 'gzip', 'lbzip2', 'sed', 'tar', 'pigz', 'python', 'r-base'}   # never the "main" tool of a process


# ---------------------------------------------------------------------------------------------------------------- 0. kind of pipeline
def read_text(path):
    return open(path, encoding='utf-8', errors='replace').read()


def is_dsl1(main_text):
    """old-style pipeline: no unnamed `workflow {}` block, the processes are wired with `from` / `into`"""
    return not UNNAMED_WORKFLOW.search(main_text) and bool(DSL1_PROCESS.search(main_text))


def pipeline_kind(path):
    """'dsl2' | 'dsl1' | '' (empty: not analysable, the second value says why)"""
    main = os.path.join(path, 'main.nf')
    if not os.path.exists(main): return '', 'no main.nf'
    txt = read_text(main)
    if UNNAMED_WORKFLOW.search(txt):
        if any(os.path.isdir(os.path.join(path, d)) for d in ('workflows', 'subworkflows', 'modules')): return 'dsl2', ''
        return '', 'no workflows / subworkflows / modules folder'
    if DSL1_PROCESS.search(txt) and re.search(r'\binto\b', txt): return 'dsl1', ''
    return '', 'main.nf has no unnamed workflow {} block and no old-style processes (different entry point)'


# ---------------------------------------------------------------------------------------------------------------- 1. graph
def extract_graph(repo, entry='main.nf'):
    if is_dsl1(read_text(os.path.join(repo, entry))):
        inst, edges, unresolved, dropped = graph_dsl1.analyse(repo, entry)
        g = {'instances': inst, 'edges': edges, 'unresolved': unresolved, 'infeasible_edges_dropped': dropped, 'dsl1': True}
    else:
        a, edges = graph_dsl2.build(repo)
        g = {'instances': list(a.instances.values()), 'edges': edges, 'unresolved': a.unresolved, 'infeasible_edges_dropped': a.infeasible}
    return json.loads(json.dumps(g))


# ---------------------------------------------------------------------------------------------------------------- 2. software of a process
def module_meta(repo, module_file):
    """tool names from the module's meta.yml and conda packages from its environment.yml"""
    d = os.path.join(repo, os.path.dirname(module_file))
    tools, env = [], []
    try:
        m = yaml.safe_load(open(os.path.join(d, 'meta.yml'), encoding='utf-8')) or {}
        for t in m.get('tools') or []: tools += list(t.keys()) if isinstance(t, dict) else [str(t)]
    except Exception: pass
    try:
        y = yaml.safe_load(open(os.path.join(d, 'environment.yml'), encoding='utf-8')) or {}
        env = [re.sub(r'^.*::', '', x).split('=')[0] for x in y.get('dependencies', []) if isinstance(x, str)]
    except Exception: pass
    return tools, env


def software_map(repo, instances):
    rows = {}
    for path, i in instances.items():
        if i['kind'] != 'process': continue
        if i.get('dsl1'): tools, env, how = [], list(i.get('software_hint') or []), 'script and directives'   # old style: no module metadata
        else: tools, env = module_meta(repo, i['module_file']); how = 'environment.yml'
        sw = ([e for e in env if e not in HELPERS][:1] or tools[:1] or [i['module'].lower()])[0]   # last resort: the module name itself
        rows[path] = {'process_instance': path, 'process_module': i['module'], 'module_file': i['module_file'], 'module_line': i['module_line'], 'software': sw,
                      'meta_yml_tools': ';'.join(tools), 'environment_deps': ';'.join(env), 'method': how if env else ('meta.yml' if tools else 'module name')}
    return rows


# ---------------------------------------------------------------------------------------------------------------- 3. payload check of one edge
# Rules for nf-core module inputs. They reject edges that a union-semantics channel analysis invents:
#   PROJECTION-DROPPED  a map closure drops the element (`{ meta, _bai -> ... }`)
#   MULTIMAP-SPLIT      a multiMap separates the bam and the index branches
# INCOMPATIBLE is only a hint and never rejects an edge.
REJECT = {'PROJECTION-DROPPED', 'MULTIMAP-SPLIT'}
FAMILIES = [('fai', r'\bfai\b|\.fai'), ('sizes', r'sizes'), ('bai', r'\bbai\b|\bcsi\b|\bcrai\b|\bindex\b.*(bam|cram)|bam_index'),
            ('bam', r'\bbam|\bcram|\bsam\b|bam_transcript'), ('fastq', r'fastq|fastx|reads|\.fq'), ('gtf', r'\bgtf\b|\.gtf|annotation'),
            ('fasta', r'fasta|\.fa\b|genome|primary_ref|gunzip'), ('tx2gene', r'tx2gene'), ('counts', r'counts|matrix|tpm|lengths|quants|results'),
            ('bed', r'\bbed\b|bedgraph|genomecov'), ('index', r'\bindex\b|\bdb\b|database')]
DROPPABLE = {'bai': ('bai', 'index', 'csi'), 'fai': ('fai',), 'lint': ('lint',), 'sizes': ('sizes',)}
BAM_CONSUMERS = ('SALMON_QUANT', 'RSEM_CALCULATEEXPRESSION', 'SENTIEON_RSEMCALCULATEEXPRESSION', 'PICARD_MARKDUPLICATES')   # `reads` means an alignment here
FASTQ_CONSUMERS = ('FQ_LINT', 'BBMAP_BBSPLIT', 'SORTMERNA', 'RIBODETECTOR', 'BOWTIE2_ALIGN')
CLOSURE = re.compile(r'\{\s*([\w\s,_]+?)\s*->')


def families(text):
    t = text.lower()
    return {f for f, rx in FAMILIES if re.search(rx, t)}


class Payload:
    def __init__(self, repo): self.repo = repo; self._lines = {}

    def lines(self, f):
        if f not in self._lines: self._lines[f] = open(os.path.join(self.repo, f), encoding='utf-8').read().splitlines()
        return self._lines[f]

    def statement(self, at):
        f, ln = at.rsplit(':', 1); L = self.lines(f); i = int(ln) - 1; buf = []; depth = 0
        for k in range(i, min(i + 12, len(L))):
            buf.append(L[k]); depth += L[k].count('{') - L[k].count('}') + L[k].count('(') - L[k].count(')')
            if depth <= 0 and k > i - 1 and ('}' in L[k] or ')' in L[k] or k > i): break
        return ' '.join(buf)

    def dropped(self, trail, pfam, consumer_call):
        """a closure with an underscore-prefixed parameter named like the producer's output family throws that element away"""
        for st in trail:
            if 'map' not in st.get('operators', []) or ':' not in st.get('at', ''): continue
            try: txt = self.statement(st['at'])
            except Exception: continue
            for m in CLOSURE.finditer(txt):
                for prm in [x.strip() for x in m.group(1).split(',')]:
                    if prm.startswith('_') and prm != '_meta' and prm.lstrip('_') != '':
                        nm = prm.lstrip('_').lower()
                        for fam, keys in DROPPABLE.items():
                            if fam in pfam and any(k in nm for k in keys): return st['at'] + ' drops ' + prm
        try:   # closures written inline in the consumer call
            for m in CLOSURE.finditer(self.statement(consumer_call)):
                for prm in [x.strip() for x in m.group(1).split(',')]:
                    nm = prm.lstrip('_').lower()
                    if prm.startswith('_') and prm != '_meta' and any(fam in pfam and any(k in nm for k in keys) for fam, keys in DROPPABLE.items()): return consumer_call + ' (inline) drops ' + prm
        except Exception: pass
        return ''

    def verdict(self, e, src_module, dst_module, dst_instance):
        """-> (verdict, reason): compatible | unknown-producer-family | INCOMPATIBLE | PROJECTION-DROPPED | MULTIMAP-SPLIT"""
        emit, odecl, idecl = e['src_emit'] or '', e['src_output_decl'] or '', e['dst_input_decl'] or ''   # outputs used by position (.out[0]) have no emit name
        pf = families(emit + ' ' + odecl)
        slots = re.findall(r'(?:path|val)\s*\(?\s*([A-Za-z_][\w]*)', idecl)
        cf = families(' '.join(slots) + ' ' + idecl)
        if any(x in ('input', 'file', 'files') for x in slots): cf |= {'bam', 'fastq', 'fasta', 'gtf', 'bed', 'counts'}
        if any(x in ('input_index', 'index') for x in slots): cf |= {'bai'}
        if 'reads' in slots and dst_module in BAM_CONSUMERS: cf |= {'bam'}
        if any(x in ('intervals', 'srbam', 'lrbam') for x in slots): cf |= {'bam'}
        if dst_module in FASTQ_CONSUMERS and 'fastq' in pf: cf |= {'fastq'}
        why = self.dropped(e['channel_trail'], pf, e['consumer_call'])
        multimap = any(st.get('var', '').startswith('map_filtered') or 'multiMap' in st.get('operators', []) for st in e['channel_trail'])
        if why: v = 'PROJECTION-DROPPED'
        elif not pf: v = 'unknown-producer-family'
        elif pf & cf: v = 'compatible'
        else: v = 'INCOMPATIBLE'
        if v == 'INCOMPATIBLE' and multimap and 'bai' in pf: v = 'MULTIMAP-SPLIT'; why = 'a multiMap on the channel path separates the bam and index outputs; the analyzer unions the branches'
        return v, why


# ---------------------------------------------------------------------------------------------------------------- 4. process edges
GATE_CLASSES = ('gating-or-metadata', 'metadata-dependency', 'key-only-synchronization')
METADATA_EMITS = {'lib_format_counts', 'json_info', 'lint'}   # inferred strandedness / lint verdicts travel in meta or act as gates
PARSE_OPS = {'splitCsv', 'splitText', 'splitJson', 'splitYaml', 'splitFasta', 'splitFastq', 'readLines'}
MERGE_OPS = {'mix', 'join', 'combine', 'merge', 'cross', 'concat'}
GATING_EMIT = re.compile(r'^(log|lint|stats|flagstat|idxstats|json|html|zip|tsv|txt|summary|report|metrics)', re.I)


def classify(e, dst_module):
    ops, flags = set(e['operators']), list(e['flags'])
    emit = e['src_emit'] or ''
    if dst_module == 'MULTIQC': return 'reporting-aggregation'
    if any(f.startswith('key-only') for f in flags): return 'key-only-synchronization'
    if emit in METADATA_EMITS: return 'metadata-dependency'
    if ops & PARSE_OPS and not ops & MERGE_OPS: return 'parsed-content'   # file CONTENT is parsed into values that feed the consumer: a data dependency
    if GATING_EMIT.match(emit) and ops & {'join', 'mix', 'combine', 'filter'}: return 'gating-or-metadata'
    if ops & MERGE_OPS: return 'merged-channel'
    return 'direct'


def process_edges(g, insts, sw, repo):
    payload, rows = Payload(repo), []
    for e in g['edges']:
        s, d = e['src_instance'], e['dst_instance']
        if s not in insts or d not in insts: continue
        dst_mod, src_mod = insts[d]['module'], insts[s]['module']
        cls = classify(e, dst_mod)
        alt = 'alternative-branch-in-channel-path' in e['flags']
        v, why = ('not-checked', '')   # the payload check is defined for direct and merged-channel edges only
        if cls in ('direct', 'merged-channel'): v, why = payload.verdict(e, src_mod, dst_mod, d)
        if v in REJECT: pair_type, reason = 'rejected', 'payload check: ' + v + (' - ' + why if why else '')
        elif cls in GATE_CLASSES: pair_type, reason = 'gate_metadata', 'no payload consumed (synchronisation / gate / metadata)'
        elif cls == 'direct': pair_type, reason = 'data_flow', 'direct single-source path' + ('; exclusive alternative branches recorded as conditions' if alt else '')
        else: pair_type, reason = 'data_flow', cls + ' (payload check: ' + v + ')'
        atoms = sorted({f'{k} {o} {x}' for c in (e['conditions_consumer'] + e['path_branch_conditions'] + e['conditions_producer']) for k, o, x in cond_atoms(c)})
        rows.append({'edge_id': f'E{len(rows) + 1:05d}', 'pair_type': pair_type, 'reason': reason, 'edge_class': cls, 'payload_check': v,
                     'src_instance': s, 'src_module': src_mod, 'src_emit': e['src_emit'] or '', 'src_software': sw[s]['software'],
                     'dst_instance': d, 'dst_module': dst_mod, 'dst_input_index': e['dst_input_index'], 'dst_software': sw[d]['software'],
                     'required_params': ' & '.join(atoms), 'operators': ';'.join(e['operators']), 'producer_call': e['producer_call'], 'consumer_call': e['consumer_call'],
                     'src_output_decl': e['src_output_decl'], 'dst_input_decl': e['dst_input_decl']})
    return rows


# ---------------------------------------------------------------------------------------------------------------- 5. software pairs
def software_pairs(prow, sw, insts):
    by = collections.defaultdict(list)
    for r in prow:
        if r['src_software'] != r['dst_software']: by[(r['src_software'], r['dst_software'])].append(r)   # a tool feeding itself is not a pair
    out = {}

    def add(k, pair_type, rows, evidence=''):
        out[k] = {'source_software': k[0], 'target_software': k[1], 'pair_type': pair_type, 'n_process_edges': len(rows),
                  'n_data_flow_edges': sum(1 for r in rows if r['pair_type'] == 'data_flow'), 'example_edge_id': rows[0]['edge_id'] if rows else '', 'evidence': evidence}

    for k, rows in sorted(by.items()):
        if any(r['pair_type'] == 'data_flow' for r in rows): add(k, 'data_flow', rows)
        elif any(r['pair_type'] == 'gate_metadata' for r in rows): add(k, 'gate_metadata', rows)
        else: add(k, 'rejected', rows, 'every extracted edge is rejected')

    # parallel: both tools consume the same upstream output, no path in either direction, and the two can run in the same run
    G = nx.DiGraph([(r['src_instance'], r['dst_instance']) for r in prow])
    desc = {n: nx.descendants(G, n) for n in G}
    members = collections.defaultdict(list)
    for p, m in sw.items(): members[m['software']].append(p)
    consumers = collections.defaultdict(set)
    for r in prow:
        if r['pair_type'] == 'data_flow': consumers[(r['src_instance'], r['src_emit'])].add((r['dst_software'], r['dst_instance'], r['edge_id']))
    shared = {}
    for port, cs in sorted(consumers.items()):
        cs = sorted(cs)
        for x in cs:
            for y in cs:
                if x[0] != y[0] and (x[0], y[0]) not in shared: shared[(x[0], y[0])] = (port, x, y)
    conds = {p: i.get('conds', []) for p, i in insts.items()}
    for (s1, s2), (port, x, y) in sorted(shared.items()):
        if (s1, s2) in out: continue
        A, B = members.get(s1, []), members.get(s2, [])
        if any(q in desc.get(p, ()) or p in desc.get(q, ()) for p in A for q in B): continue
        if bool(A and B) and all(infeasible(conds.get(p, []) + conds.get(q, [])) for p in A for q in B): continue   # exclusive alternatives: never in one run, not a pair
        add((s1, s2), 'parallel_co_use', [], 'both consume ' + port[0].split(':')[-1] + '.' + str(port[1]) + ' (' + x[2] + ', ' + y[2] + '); no directed path between any process of the two tools')
    return list(out.values()), G


# ---------------------------------------------------------------------------------------------------------------- 6. checks and output
def checks(repo, prow, G):
    warn, nlines = [], {}
    cyc = [sorted(c)[:2] for c in nx.strongly_connected_components(G) if len(c) > 1]
    if cyc: warn.append('process graph has cycles: ' + str(cyc[:3]))
    ref = re.compile(r'([^\s:]+):(\d+)')
    bad = []
    for r in prow:
        for c in ('producer_call', 'consumer_call'):
            m = ref.fullmatch(r[c] or ''); pth = os.path.join(repo, m.group(1)) if m else None
            if not m or not os.path.exists(pth): bad.append((r['edge_id'], c, r[c])); continue
            if pth not in nlines: nlines[pth] = sum(1 for _ in open(pth, encoding='utf-8', errors='replace'))
            if int(m.group(2)) > nlines[pth]: bad.append((r['edge_id'], c, r[c]))
    if bad: warn.append('call sites that are not an existing file:line: ' + str(bad[:3]))
    return warn


def write_csv(path, rows, cols=None):
    with open(path, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=cols or (list(rows[0].keys()) if rows else ['x']), extrasaction='ignore', lineterminator='\n'); w.writeheader(); w.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--repo', required=True); ap.add_argument('--out', required=True); ap.add_argument('--entry', default='main.nf')
    a = ap.parse_args()
    repo, out = os.path.abspath(a.repo), os.path.abspath(a.out)
    if not os.path.exists(os.path.join(repo, a.entry)): sys.exit(a.entry + ' not found in ' + repo)
    os.makedirs(out, exist_ok=True)
    g = extract_graph(repo, a.entry)
    insts = {i['path']: i for i in g['instances']}
    n_proc = sum(1 for i in g['instances'] if i['kind'] == 'process')
    summary = {'pipeline_kind': 'dsl1' if g.get('dsl1') else 'dsl2', 'process_instances': n_proc, 'unresolved_statements': len(g['unresolved'])}
    if not g['edges']:   # e.g. a pipeline made of independent processes: a valid result, nothing to tabulate
        summary.update(status='NO_EDGES', process_edges={}, software_pairs={}, warnings=[])
        json.dump(summary, open(os.path.join(out, 'summary.json'), 'w'), indent=1)
        print('NO_EDGES: %d process instances but no data flows between them' % n_proc); return
    sw = software_map(repo, insts)
    prow = process_edges(g, insts, sw, repo)
    srows, G = software_pairs(prow, sw, insts)
    write_csv(os.path.join(out, 'process_software_map.csv'), list(sw.values()))
    write_csv(os.path.join(out, 'process_edges.csv'), prow)
    write_csv(os.path.join(out, 'software_pairs.csv'), srows, ['source_software', 'target_software', 'pair_type', 'n_process_edges', 'n_data_flow_edges', 'example_edge_id', 'evidence'])
    summary.update(status='OK', process_edges=dict(collections.Counter(r['pair_type'] for r in prow)), software_pairs=dict(collections.Counter(r['pair_type'] for r in srows)), warnings=checks(repo, prow, G))
    json.dump(summary, open(os.path.join(out, 'summary.json'), 'w'), indent=1)
    print('process edges:', summary['process_edges']); print('software pairs:', summary['software_pairs'])
    for w in summary['warnings']: print('WARNING:', w)


if __name__ == '__main__':
    main()
