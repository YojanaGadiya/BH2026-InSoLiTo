#!/usr/bin/env python3
"""Static channel-provenance analysis for OLD-STYLE (Nextflow DSL1) pipelines: one main.nf, `process x { input: ... from ch  output: ... into ch2 }`,
channels wired with `.set{}`, `.into{}`, `ch = a.mix(b)` and operator chains. There are no modules, so software is guessed from the process
directives (container / conda) and the commands of the script block.

Returns the same graph structure as graph_dsl2.py, so pairs.py treats both kinds of pipeline the same way (analyse() is called by pairs.py; the command line below is for debugging).
Semantics (same as the DSL2 analyzer): an edge P -> Q exists when a channel that Q reads in `from` can carry (through any chain of operators,
union semantics) an output declared by P in `into`. Source order, shared inputs and co-occurrence create no edge.
Limitations: only the entry file is read; for/while/try blocks are not interpreted; channels built inside functions are invisible.
"""
import argparse, json, os, re, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from graph_dsl2 import strip_comments, infeasible

IDENT = re.compile(r'[A-Za-z_]\w*')
SECTION = re.compile(r'^[ \t]*(input|output|when|script|shell|exec|stub)[ \t]*:', re.M)
PROC = re.compile(r'^[ \t]*process[ \t]+(\w+)[ \t\r\n]*\{', re.M)
OP_CALL = re.compile(r'\.\s*(\w+)\s*[\(\{]')

SHELL_SKIP = set('''echo cat mkdir cd ln mv cp rm touch set export true false exit test read sed awk gawk grep egrep fgrep sort cut head tail wc tr xargs gzip gunzip zcat
tee basename dirname date printf source trap local eval if then else elif fi for do done while case esac return ls find chmod pwd uniq paste join env which time bc expr
unset shift wait sleep tar unzip zip pigz mktemp nproc command declare let tput ulimit rev column diff cmp comm split def print println params final new task workflow'''.split())
INTERP = {'python', 'python2', 'python3', 'rscript', 'r', 'perl', 'bash', 'sh', 'java', 'ruby', 'node'}


def mask_strings(t):
    """Blank the CONTENT of string literals (keeps offsets and newlines)."""
    out = list(t); i = 0; n = len(t)
    while i < n:
        c = t[i]
        if t.startswith("'''", i) or t.startswith('"""', i):
            q = t[i:i + 3]; j = t.find(q, i + 3); j = n if j < 0 else j
            for k in range(i + 3, j):
                if out[k] != '\n': out[k] = ' '
            i = j + 3; continue
        if c in '\'"':
            j = i + 1
            while j < n and t[j] != c and t[j] != '\n': j += 2 if t[j] == '\\' else 1
            for k in range(i + 1, min(j, n)): out[k] = ' '
            i = j + 1; continue
        i += 1
    return ''.join(out)


def match_close(t, i, o, c):
    """index of the bracket that closes the one at t[i]; -1 if none"""
    d = 0
    for k in range(i, len(t)):
        if t[k] == o: d += 1
        elif t[k] == c:
            d -= 1
            if d == 0: return k
    return -1


def squeeze(s): return re.sub(r'\s+', ' ', s).strip()


def if_blocks(sk):
    """[(start, end, cond)] for every `if (c) { ... }` and its `else { ... }` (cond NOT(c))"""
    blocks = []
    for m in re.finditer(r'\bif\s*\(', sk):
        p = sk.index('(', m.start()); q = match_close(sk, p, '(', ')')
        if q < 0: continue
        j = q + 1
        while j < len(sk) and sk[j] in ' \t\r\n': j += 1
        if j >= len(sk) or sk[j] != '{': continue
        e = match_close(sk, j, '{', '}')
        if e < 0: continue
        cond = squeeze(sk[p + 1:q]); blocks.append((j, e, cond, 'if'))
        k = e + 1
        while k < len(sk) and sk[k] in ' \t\r\n': k += 1
        if sk.startswith('else', k):
            k2 = k + 4
            while k2 < len(sk) and sk[k2] in ' \t\r\n': k2 += 1
            if k2 < len(sk) and sk[k2] == '{':
                e2 = match_close(sk, k2, '{', '}')
                if e2 > 0: blocks.append((k2, e2, 'NOT(' + cond + ')', 'else'))
    return blocks


def split_decls(sk, orig):
    """split an input:/output: section into logical declarations -> [(text, skeleton_text, start_offset)]"""
    res = []; cur_o = []; cur_s = []; depth = 0; start = None; pos = 0
    lines = sk.split('\n'); olines = orig.split('\n')
    for idx, (ls, lo) in enumerate(zip(lines, olines)):
        s = ls.strip()
        off = pos; pos += len(ls) + 1
        if not s:
            continue
        if not cur_s: start = off + (len(ls) - len(ls.lstrip()))
        cur_s.append(s.rstrip('\\')); cur_o.append(lo.strip().rstrip('\\'))
        depth += sum(s.count(c) for c in '([{') - sum(s.count(c) for c in ')]}')
        nxt = ''
        for l2 in lines[idx + 1:]:
            if l2.strip(): nxt = l2.strip(); break
        if depth > 0 or s.endswith((',', '\\', '.', '=')) or re.match(r'(from|into|\.)\b|\.', nxt):
            continue
        res.append((' '.join(cur_o), ' '.join(cur_s), start)); cur_o = []; cur_s = []; depth = 0
    if cur_s: res.append((' '.join(cur_o), ' '.join(cur_s), start))
    return res


def top_level_find(sk, word, last=True):
    """index of `word` (as a whole word) at bracket depth 0, last or first occurrence; -1 if none"""
    depth = 0; found = -1
    for m in re.finditer(r'[(\[{)\]}]|\b' + word + r'\b', sk):
        t = m.group(0)
        if t in '([{': depth += 1
        elif t in ')]}': depth -= 1
        elif depth == 0:
            found = m.start()
            if not last: return found
    return found


def script_hint(script, directives, procname=''):
    tools = []
    for m in re.finditer(r"container\s+['\"]([^'\"]+)['\"]", directives):
        n = m.group(1).split('/')[-1].split(':')[0].split('@')[0]
        if n and '$' not in n: tools.append(n.lower())
    for m in re.finditer(r"conda\s+['\"]([^'\"]+)['\"]", directives):
        for pk in m.group(1).split():
            n = re.sub(r'^.*::', '', pk).split('=')[0]
            if n and '$' not in n and '/' not in n: tools.append(n.lower())
    for line in script.splitlines():
        for seg in re.split(r'\|\||\||&&|;', line):
            toks = seg.split()
            while toks and re.match(r'^\w+=', toks[0]): toks.pop(0)
            if not toks or (len(toks) > 1 and toks[1].startswith(('=', '+=', '?='))) or re.match(r'^\w+=', toks[0]): continue   # variable assignment, not a command
            t = toks[0]
            if '$' in t or t[0] in '-({[<>"\'\\#!%' or t.startswith(('/*', '*')): continue
            t = t.split('/')[-1]
            if t.lower() in INTERP:
                cand = None
                if t.lower() == 'java':
                    for k, x in enumerate(toks):
                        if x == '-jar' and k + 1 < len(toks): cand = toks[k + 1]
                else:
                    cand = next((x for x in toks[1:] if not x.startswith('-')), None)
                if cand and '$' not in cand:
                    cand = re.sub(r'\.(py|R|r|pl|sh|jar|rb|js)$', '', cand.split('/')[-1])
                    if cand: tools.append(cand.lower())
                continue
            if t.lower() in SHELL_SKIP or not re.match(r'^[A-Za-z][\w.+-]*$', t): continue
            tools.append(re.sub(r'\.(py|R|r|pl|sh)$', '', t).lower())
    tools = [x for x in tools if re.match(r'^[a-z][\w+-]*$', x)]
    n_dir = len(re.findall(r"container|conda", directives))
    first = tools[:n_dir] if n_dir else []              # container/conda names come first and keep their order
    rest = tools[n_dir:] if n_dir else tools
    pn = re.sub(r'[^a-z0-9]', '', procname.lower())
    rest = sorted(dict.fromkeys(rest), key=lambda x: (not (re.sub(r'[^a-z0-9]', '', x) in pn and len(x) > 2), -rest.count(x)))   # a command named like the process first, then most-used (stable)
    out = list(dict.fromkeys(first + rest))
    return out[:4]


def analyse(root, entry='main.nf'):
    path = os.path.join(root, entry)
    raw = open(path, encoding='utf-8', errors='replace').read()
    txt = strip_comments(raw); sk = mask_strings(txt)
    line = lambda off: txt.count('\n', 0, off) + 1
    blocks = if_blocks(sk)

    # ---------------------------------------------------------------- processes
    procs = []
    for m in PROC.finditer(sk):
        b0 = sk.index('{', m.end() - 1); b1 = match_close(sk, b0, '{', '}')
        if b1 < 0: continue
        procs.append({'name': m.group(1), 'start': m.start(), 'b0': b0, 'b1': b1, 'line': line(m.start() + len(m.group(0)) - len(m.group(0).lstrip()))})
    masked = list(sk)
    for p in procs:
        for k in range(p['start'], p['b1'] + 1):
            if masked[k] != '\n': masked[k] = ' '
    outside = ''.join(masked)

    instances = {}; unresolved = []; env = {}
    for p in procs:
        body_sk = sk[p['b0'] + 1:p['b1']]; body = txt[p['b0'] + 1:p['b1']]; base = p['b0'] + 1
        heads = list(SECTION.finditer(body_sk))
        secs = {}
        for i, h in enumerate(heads):
            end = heads[i + 1].start() if i + 1 < len(heads) else len(body_sk)
            secs.setdefault(h.group(1), (h.end(), end))
        directives = body[:heads[0].start()] if heads else body
        def sect(n): return secs.get(n)
        conds = [c for (a, b, c, _k) in sorted(blocks) if a < p['start'] < b]
        if sect('when'):
            a, b = sect('when'); w = squeeze(body[a:b])
            if w: conds.append(w)
        sc = sect('script') or sect('shell') or sect('exec')
        inputs, outputs = [], []
        if sect('input'):
            a, b = sect('input')
            for dtxt, dsk, off in split_decls(body_sk[a:b], body[a:b]):
                f = top_level_find(dsk, 'from')
                if f < 0: inputs.append({'decl': squeeze(dtxt), 'expr': '', 'line': line(base + a + off)}); continue
                inputs.append({'decl': squeeze(dtxt[:f]), 'expr': squeeze(dtxt[f + 4:]), 'line': line(base + a + off)})
        if sect('output'):
            a, b = sect('output')
            for k, (dtxt, dsk, off) in enumerate(split_decls(body_sk[a:b], body[a:b])):
                f = top_level_find(dsk, 'into')
                names = []; emit = None; decl = dtxt
                if f >= 0:
                    decl = dtxt[:f]; tail = dtxt[f + 4:]
                    mm = re.search(r'\bemit\s*:\s*(\w+)', tail)
                    if mm: emit = mm.group(1)
                    tail = re.sub(r'\b(emit|mode|optional|saveAs)\s*:\s*\S+', '', tail)
                    names = IDENT.findall(tail)
                outputs.append({'emit': emit or (names[0] if names else None), 'decl': squeeze(decl), 'topic': False, 'topic_name': None, 'channels': names, 'line': line(base + a + off)})
        hint = script_hint(txt[base + sc[0]:base + sc[1]] if sc else '', directives, p['name'])
        instances[p['name']] = {'path': p['name'], 'kind': 'process', 'name': p['name'], 'module': p['name'], 'module_file': entry, 'module_line': p['line'], 'call_file': entry, 'call_line': p['line'],
                                'conds': conds, 'inputs': inputs, 'outputs': outputs, 'dsl1': True, 'software_hint': hint}
        for k, o in enumerate(outputs):
            for ch in o['channels']:
                env.setdefault(ch, {})[(p['name'], k, o['emit'])] = ((ch, f'{entry}:{o["line"]}', ('output',)),)

    # ---------------------------------------------------------------- top-level statements
    work = list(outside)
    for a, b, _c, _k in blocks:
        if a < len(work) and not any(p['start'] <= a <= p['b1'] for p in procs):
            work[a] = ';'
            if b < len(work): work[b] = ';'
    work = ''.join(work)
    stmts = []; depth = 0; s0 = None; n = len(work)
    for i, c in enumerate(work):
        if s0 is None and not c.isspace() and c != ';': s0 = i
        if c in '([{': depth += 1
        elif c in ')]}': depth = max(0, depth - 1)
        end = False
        if depth == 0 and c == ';': end = True
        elif depth == 0 and c == '\n' and s0 is not None:
            prev = work[s0:i].rstrip(); k = i + 1
            while k < n and work[k] in ' \t\r\n': k += 1
            nxt = work[k:k + 2]
            if not (prev.endswith(('.', ',', '\\', '+', '?', ':', '=', '&&', '||', '|')) or nxt.startswith(('.', '?', ':', '|', '&&', '+')) or re.match(r'(into|from)\b', work[k:k + 5])): end = True
        if end and s0 is not None:
            stmts.append((s0, i)); s0 = None
    if s0 is not None: stmts.append((s0, n))

    def strip_lead(s):
        while True:
            s = s.lstrip()
            m = re.match(r'if\s*\(', s)
            if m:
                p = s.index('('); q = match_close(s, p, '(', ')')
                if q < 0: return s
                s = s[q + 1:]; continue
            if re.match(r'else\b', s): s = s[4:]; continue
            return s

    def eval_stmt(a, b):
        s = strip_lead(work[a:b])
        if not s.strip(): return 0
        orig = txt[a:b]; ln = line(a)
        names = []
        rhs = s
        m = re.match(r'(?:def\s+|final\s+)?\(\s*([\w\s,]+?)\s*\)\s*=(?!=)\s*(.*)', s, re.S)
        if m: names += IDENT.findall(m.group(1)); rhs = m.group(2)
        else:
            m = re.match(r'(?:def\s+|final\s+|Channel\s+)?([A-Za-z_]\w*)\s*=(?!=)\s*(.*)', s, re.S)
            if m: names.append(m.group(1)); rhs = m.group(2)
        spans = []
        for mm in re.finditer(r'\.\s*(set|tap)\s*[\{\(]\s*(\w+)\s*[\}\)]', rhs):
            names.append(mm.group(2)); spans.append(mm.span())
        for mm in re.finditer(r'\.\s*into\s*[\{\(]([^\}\)]*)[\}\)]', rhs):
            spans.append(mm.span())
            if not mm.group(1).strip().isdigit(): names += IDENT.findall(mm.group(1))
        if not names: return 0
        keep = list(rhs)
        for x, y in spans:
            for k in range(x, y): keep[k] = ' '
        rest = ''.join(keep)
        ops = tuple(sorted({o for o in OP_CALL.findall(rest) if o not in ('set', 'into', 'tap')}))
        src = {}
        for im in IDENT.finditer(rest):
            before = rest[:im.start()].rstrip(); after = rest[im.end():].lstrip()
            if before.endswith('.') or after.startswith('('): continue
            for key, tr in env.get(im.group(0), {}).items(): src.setdefault(key, tr)
        added = 0
        for nm in dict.fromkeys(names):
            tgt = env.setdefault(nm, {})
            for key, tr in src.items():
                if key not in tgt: tgt[key] = tr + ((nm, f'{entry}:{ln}', ops),); added += 1
        return added

    for _pass in range(8):
        if not sum(eval_stmt(a, b) for a, b in stmts): break

    # ---------------------------------------------------------------- edges
    edges = []; dropped = []; seen = set()
    for q, qi in instances.items():
        for idx, inp in enumerate(qi['inputs']):
            expr = inp['expr']
            if not expr: continue
            sk_e = mask_strings(expr); ops = tuple(sorted(set(OP_CALL.findall(sk_e))))
            srcs = {}; used = False
            for im in IDENT.finditer(sk_e):
                before = sk_e[:im.start()].rstrip(); after = sk_e[im.end():].lstrip()
                if before.endswith('.') or after.startswith('('): continue
                if im.group(0) in env: used = True
                for key, tr in env.get(im.group(0), {}).items(): srcs.setdefault(key, (im.group(0), tr))
            if not used and re.match(r'^[a-z]\w*$', expr.strip()) and expr.strip() not in ('params', 'null', 'true', 'false'):
                unresolved.append({'kind': 'dsl1-channel-without-producer', 'file': entry, 'line': inp['line'], 'process': q, 'channel': expr.strip()})
            for (p, k, emit), (var, tr) in srcs.items():
                if p == q: continue
                pi = instances[p]
                allc = list(qi['conds']) + list(pi['conds'])
                od = pi['outputs'][k]
                if infeasible(allc):
                    dropped.append({'src': p, 'emit': emit, 'dst': q, 'call': f'{entry}:{inp["line"]}', 'conditions': allc}); continue
                key = (p, emit, q, idx)
                if key in seen: continue
                seen.add(key)
                flags = []
                if qi['conds']: flags.append('conditional-consumer')
                if pi['conds']: flags.append('conditional-producer')
                trail = [{'var': s[0], 'at': s[1], 'operators': list(s[2])} for s in tr] + [{'var': var, 'at': f'{entry}:{inp["line"]}', 'operators': list(ops)}]
                alle = sorted({o for t in trail for o in t['operators']} - {'output'})
                edges.append({'src_instance': p, 'src_emit': emit, 'src_output_decl': od['decl'], 'dst_instance': q, 'dst_input_index': idx, 'dst_input_decl': inp['decl'],
                              'consumer_call': f'{entry}:{inp["line"]}', 'producer_call': f'{entry}:{pi["call_line"]}', 'conditions_consumer': qi['conds'], 'conditions_producer': pi['conds'],
                              'channel_trail': trail, 'path_branch_conditions': [], 'operators': alle, 'flags': flags})
    return list(instances.values()), edges, unresolved, dropped


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--root', required=True); ap.add_argument('--out', required=True); ap.add_argument('--entry', default='main.nf')
    a = ap.parse_args()
    inst, edges, unres, dropped = analyse(a.root, a.entry)
    os.makedirs(a.out, exist_ok=True)
    json.dump({'instances': inst, 'edges': edges, 'unresolved': unres, 'infeasible_edges_dropped': dropped, 'dsl1': True}, open(os.path.join(a.out, 'static_graph.json'), 'w'), indent=1)
    print(len(inst), 'process instances (DSL1);', len(edges), 'edges;', len(unres), 'unresolved;', len(dropped), 'infeasible edges dropped')
