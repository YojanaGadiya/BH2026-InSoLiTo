#!/usr/bin/env python3
"""Static channel-provenance analyzer for Nextflow DSL2 pipelines.

Goal: derive PROCESS-level directed data dependencies  P --(emit)--> Q(input i)
meaning "a channel carrying an output of process instance P is passed as the
i-th argument of process instance Q".  No transitive closure is performed: an
edge only exists where the producer's channel reaches the consumer's call
through channel variables, operators and workflow take/emit interfaces.

Method (all verified by reading the source; heuristics are flagged):
  * Instances are named by call path (e.g. NFCORE_PIPELINE:WORKFLOW:ALIGN_STAR:STAR_ALIGN),
    identical to Nextflow's own fully-qualified process names.
  * Each channel variable carries a *provenance set* of (instance, emit) ports.
  * Operators are treated conservatively as union of receiver + argument
    provenance (mix/join/combine/map/...); the operator names are kept as evidence.
  * Version/topic emits are excluded from data ports.
  * if/else branches are both explored; the branch conditions are recorded.
"""
import re, json, sys, os, argparse
from collections import OrderedDict

IDENT = re.compile(r'[A-Za-z_]\w*')
CASE_RX = re.compile(r'(case\b[^:\n]*|default\b)\s*:(?!:)')
TOPIC_RX = re.compile(r'[Cc]hannel\s*\.\s*topic\s*\(\s*[\x27\x22](\w+)[\x27\x22]\s*\)')   # channel.topic('name')

# ----------------------------------------------------------------- text utils
def strip_comments(t):
    """Blank out comments (keeping offsets/newlines). Strings are respected."""
    out = list(t); i = 0; n = len(t)
    while i < n:
        c = t[i]
        if t.startswith("'''", i) or t.startswith('"""', i):
            q = t[i:i+3]; j = t.find(q, i+3); i = (j + 3) if j >= 0 else n; continue
        if c in '\'"':
            j = i + 1
            while j < n and t[j] != c and t[j] != '\n':
                j += 2 if t[j] == '\\' else 1
            i = j + 1; continue
        if t.startswith('//', i):
            j = t.find('\n', i); j = n if j < 0 else j
            for k in range(i, j): out[k] = ' '
            i = j; continue
        if t.startswith('/*', i):
            j = t.find('*/', i + 2); j = n if j < 0 else j + 2
            for k in range(i, j):
                if out[k] != '\n': out[k] = ' '
            i = j; continue
        i += 1
    return ''.join(out)

def skip_string(t, i):
    c = t[i]
    if t.startswith(c * 3, i):
        j = t.find(c * 3, i + 3); return (j + 3) if j >= 0 else len(t)
    j = i + 1
    while j < len(t) and t[j] != c:
        j += 2 if t[j] == '\\' else 1
    return j + 1

def match(t, i):
    """index just after the bracket matching t[i]"""
    op = t[i]; cl = {'(': ')', '[': ']', '{': '}'}[op]
    d = 0; j = i
    while j < len(t):
        c = t[j]
        if c in '\'"': j = skip_string(t, j); continue
        if c == op: d += 1
        elif c == cl:
            d -= 1
            if d == 0: return j + 1
        j += 1
    return len(t)

def blank_strings(t):
    out = list(t); i = 0
    while i < len(t):
        if t[i] in '\'"':
            j = skip_string(t, i)
            for k in range(i + 1, j - 1):
                if out[k] != '\n': out[k] = ' '
            i = j
        else: i += 1
    return ''.join(out)

def split_top(t, s, e, sep=','):
    parts = []; d = 0; st = s; j = s
    while j < e:
        c = t[j]
        if c in '\'"': j = skip_string(t, j); continue
        if c in '([{': d += 1
        elif c in ')]}': d -= 1
        elif c == sep and d == 0:
            parts.append((st, j)); st = j + 1
        j += 1
    if t[st:e].strip(): parts.append((st, e))
    return parts

def split_statements(t, s, e):
    out = []; i = s
    while i < e:
        while i < e and t[i] in ' \t\r\n;': i += 1
        if i >= e: break
        st = i; d = 0; j = i
        while j < e:
            c = t[j]
            if c in '\'"': j = skip_string(t, j); continue
            if c in '([{': d += 1
            elif c in ')]}': d -= 1
            elif c == '\n' and d == 0:
                prev = t[st:j].rstrip()
                k = j
                while k < e and t[k] in ' \t\r\n': k += 1
                nxt = t[k:k + 6]
                if prev and prev[-1] in ',+=&|?:' : j += 1; continue
                if nxt[:1] in ('.', '?', ':', '|') or nxt.startswith(('&&', 'else')): j += 1; continue
                break
            j += 1
        out.append((st, min(j, e))); i = j + 1
    return out


def chain_start(txt, pos):
    """walk backwards from pos over a method chain (identifiers, dots, balanced brackets)"""
    j = pos - 1
    while j >= 0:
        while j >= 0 and txt[j].isspace(): j -= 1
        if j < 0: break
        c = txt[j]
        if c in ')]}':
            d = 0
            while j >= 0:
                if txt[j] in ')]}': d += 1
                elif txt[j] in '([{':
                    d -= 1
                    if d == 0: break
                j -= 1
            j -= 1
        elif c.isalnum() or c == '_':
            while j >= 0 and (txt[j].isalnum() or txt[j] == '_'): j -= 1
        elif c == '.':
            j -= 1
        else:
            break
    return j + 1

def key_only_spans(txt):
    """spans of `<chain>.map { a, _b -> a }` style expressions that keep only the sample key"""
    spans = []
    for m in re.finditer(r'\.map\s*\{\s*([\w\s,]*?)\s*->', txt):
        o = txt.index('{', m.start()); e = match(txt, o)
        body = txt[m.end():e - 1].strip()
        params = [x for x in re.split(r'[\s,]+', m.group(1)) if x]
        if params and re.fullmatch(r'\[?\s*(%s)(\.id)?\s*\]?' % re.escape(params[0]), body):
            spans.append((chain_start(txt, m.start()), e))
    return spans

def is_key(tr): return any(st[0] == 'keyonly' for st in tr)
def better(new, old):
    if is_key(old) != is_key(new): return is_key(old)
    return len(new) < len(old)


# ------------------------------------------------ light-weight condition feasibility
_ATOM_EQ = re.compile(r"^!?\(?\s*(?:params\.)?([A-Za-z_]\w*)\s*(==|!=)\s*['\"]([^'\"]+)['\"]\s*\)?$")
_ATOM_BOOL = re.compile(r"^(!?)\(?\s*(?:params\.)?([A-Za-z_]\w*)\s*\)?$")

def _atom(t):
    t = t.strip()
    m = _ATOM_EQ.match(t)
    if m:
        neg = t.startswith('!'); op = m.group(2)
        eq = (op == '==') != neg
        return (m.group(1), 'eq' if eq else 'ne', m.group(3))
    m = _ATOM_BOOL.match(t)
    if m: return (m.group(2), 'truthy', m.group(1) != '!')
    return None

def cond_atoms(cond):
    cond = cond.strip()
    if cond.startswith('NOT(') and cond.endswith(')'):
        inner = cond[4:-1]
        if '&&' in inner or '||' in inner: return []
        a = _atom(inner)
        if not a: return []
        k, op, v = a
        return [(k, {'eq': 'ne', 'ne': 'eq'}.get(op, op), (not v) if op == 'truthy' else v)]
    if '||' in cond: return []
    atoms = []
    for part in cond.split('&&'):
        a = _atom(part)
        if a: atoms.append(a)
    return atoms

def infeasible(conds):
    atoms = [a for c in conds for a in cond_atoms(c)]
    for i, (k1, o1, v1) in enumerate(atoms):
        for (k2, o2, v2) in atoms[i + 1:]:
            if k1 != k2: continue
            if o1 == 'eq' and o2 == 'eq' and v1 != v2: return True
            if {o1, o2} == {'eq', 'ne'} and v1 == v2: return True
            if o1 == o2 == 'truthy' and v1 != v2: return True
    return False

# --------------------------------------------------------------- file model
class NFFile:
    def __init__(self, path, root):
        self.path = os.path.abspath(path)
        self.rel = os.path.relpath(self.path, root).replace('\\', '/')
        raw = open(path, encoding='utf-8').read()
        self.raw = raw
        self.text = strip_comments(raw)
        self.includes = {}      # alias -> (orig name, resolved file path)
        self.processes = {}     # name -> dict
        self.workflows = {}     # name -> dict
        self._parse()

    def line(self, off): return self.text.count('\n', 0, off) + 1

    def _parse(self):
        t = self.text
        for m in re.finditer(r'include\s*\{(.*?)\}\s*from\s*[\'"]([^\'"]+)[\'"]', t, re.S):
            p = m.group(2)
            base = os.path.normpath(os.path.join(os.path.dirname(self.path), p))
            cands = [base + '.nf', os.path.join(base, 'main.nf'), base]
            res = next((c for c in cands if os.path.isfile(c)), None)
            for item in re.split(r'[;\n]', m.group(1)):
                item = item.strip()
                if not item: continue
                mm = re.match(r'(\w+)(?:\s+as\s+(\w+))?$', item)
                if mm: self.includes[mm.group(2) or mm.group(1)] = (mm.group(1), res, self.line(m.start()))
        for m in re.finditer(r'^\s*process\s+(\w+)\s*\{', t, re.M):
            o = t.index('{', m.start()); e = match(t, o)
            self.processes[m.group(1)] = self._process(m.group(1), o + 1, e - 1, m.start())
        for m in re.finditer(r'^\s*workflow\s*(\w*)\s*\{', t, re.M):
            o = t.index('{', m.start()); e = match(t, o)
            self.workflows[m.group(1) or '<entry>'] = self._workflow(o + 1, e - 1, m.start())

    def _sections(self, s, e, labels):
        t = self.text
        pos = []
        for m in re.finditer(r'^[ \t]*(%s)[ \t]*:' % '|'.join(labels), t[s:e], re.M):
            pos.append((m.group(1), s + m.end(), s + m.start()))
        secs = {}
        for k, (lab, st, ls) in enumerate(pos):
            en = pos[k + 1][2] if k + 1 < len(pos) else e
            secs[lab] = (st, en)
        return secs

    def _process(self, name, s, e, defpos):
        t = self.text
        secs = self._sections(s, e, ['input', 'output', 'when', 'script', 'shell', 'exec', 'stub', 'tag', 'label', 'conda', 'container'])
        ins = []
        if 'input' in secs:
            for a, b in split_statements(t, *secs['input']):
                ins.append({'decl': ' '.join(t[a:b].split()), 'line': self.line(a)})
        outs = []
        if 'output' in secs:
            for a, b in split_statements(t, *secs['output']):
                d = ' '.join(t[a:b].split())
                em = re.search(r'emit:\s*(\w+)', d)
                outs.append({'decl': d, 'emit': em.group(1) if em else None,
                             'topic': bool(re.search(r'topic:\s*versions', d)), 'topic_name': (re.search(r'topic:\s*(\w+)', d) or [None, None])[1], 'line': self.line(a)})
        return {'name': name, 'file': self.rel, 'line': self.line(defpos), 'inputs': ins, 'outputs': outs}

    def _workflow(self, s, e, defpos):
        t = self.text
        secs = self._sections(s, e, ['take', 'main', 'emit', 'publish'])
        takes = []
        if 'take' in secs:
            for a, b in split_statements(t, *secs['take']):
                nm = t[a:b].strip().split()
                if nm: takes.append(nm[0])
        if 'main' in secs: body = secs['main']
        else:
            body = (secs['take'][1] if 'take' in secs else s, secs['emit'][0] if 'emit' in secs else e)
            if 'take' in secs:  # no main label: body is after takes -> ambiguous, keep statements after take section
                body = (secs['take'][0], secs['emit'][0] - 5 if 'emit' in secs else e)
        emits = []
        if 'emit' in secs:
            ee = secs.get('publish', (0, e))[0] - 8 if 'publish' in secs and secs['publish'][0] > secs['emit'][0] else secs['emit'][1]
            for a, b in split_statements(t, secs['emit'][0], secs['emit'][1]):
                emits.append((a, b))
        return {'takes': takes, 'body': body, 'emits': emits, 'line': self.line(defpos)}

# ----------------------------------------------------------------- analyzer
class Analyzer:
    def __init__(self, root, entry='main.nf'):
        self.root = os.path.abspath(root)
        self.files = {}
        self.instances = OrderedDict()
        self.edges = []
        self.registry = {}      # wf instance -> {emit: prov}
        self.unresolved = []
        self.infeasible = []
        self.entry = self.load(os.path.join(self.root, entry))

    def load(self, p):
        p = os.path.abspath(p)
        if p not in self.files: self.files[p] = NFFile(p, self.root)
        return self.files[p]

    # ---- scope resolution
    def callee(self, f, name):
        """Return ('process'|'workflow', defining NFFile, def dict, orig_name) or None"""
        if name in f.workflows and name != '<entry>': return ('workflow', f, f.workflows[name], name)
        if name in f.processes: return ('process', f, f.processes[name], name)
        if name in f.includes:
            orig, res, _ = f.includes[name]
            if res:
                g = self.load(res)
                if orig in g.processes: return ('process', g, g.processes[orig], orig)
                if orig in g.workflows: return ('workflow', g, g.workflows[orig], orig)
        return None

    def callable_names(self, f):
        names = set(f.includes) | set(n for n in f.workflows if n != '<entry>') | set(f.processes)
        return names

    # ---- run
    def run(self):
        f = self.entry
        wf = f.workflows['<entry>']
        ctx = dict(f=f, path='', env={}, conds=[], names=self.callable_names(f))
        self.exec_range(ctx, *wf['body'])
        return self

    def exec_range(self, ctx, s, e):
        f = ctx['f']
        for a, b in split_statements(f.text, s, e):
            self.exec_stmt(ctx, a, b)

    def exec_stmt(self, ctx, a, b):
        f = ctx['f']; t = f.text
        s = t[a:b]
        if re.match(r'if\s*\(', s): return self.exec_if(ctx, a, b)
        if re.match(r'(def|final)\s+\w+\s*=\s*\{', s) or re.match(r'(def|final)\s+\w+\s*\(', s): return
        if re.match(r'switch\s*\(', s): return self.exec_switch(ctx, a, b)
        if re.match(r'(for|while|try)\b', s):
            self.unresolved.append({'kind': 'unsupported-statement', 'file': f.rel, 'line': f.line(a), 'text': ' '.join(s.split())[:120]})
        m = re.match(r'(?:def\s+|final\s+)?(\w+)\s*=(?!=)\s*', s)
        if m:
            var = m.group(1); rs = a + m.end()
            prov = self.eval_expr(ctx, rs, b, var)
            ctx['env'][var] = prov
            return
        m = re.match(r'(?:def\s+|final\s+)?\(([\w\s,]+)\)\s*=\s*', s)
        if m:
            prov = self.eval_expr(ctx, a + m.end(), b, None)
            for v in re.split(r'\s*,\s*', m.group(1).strip()): ctx['env'][v] = prov
            return
        # expression statement; handle .set { name }
        ms = re.search(r'\.set\s*\{\s*(\w+)\s*\}\s*$', s)
        if ms:
            prov = self.eval_expr(ctx, a, a + ms.start(), ms.group(1))
            ctx['env'][ms.group(1)] = prov
            return
        mi = re.search(r'\.into\s*\{([^}]*)\}\s*$', s)
        if mi:
            prov = self.eval_expr(ctx, a, a + mi.start(), None)
            for v in re.split(r'[;\s]+', mi.group(1).strip()):
                if v: ctx['env'][v] = prov
            return
        self.eval_expr(ctx, a, b, None)

    def exec_if(self, ctx, a, b):
        f = ctx['f']; t = f.text
        po = t.index('(', a); pe = match(t, po)
        cond = ' '.join(t[po + 1:pe - 1].split())
        i = pe
        while i < b and t[i] in ' \t\r\n': i += 1
        def block(i):
            if t[i] == '{':
                e = match(t, i); return (i + 1, e - 1), e
            j = t.find('\n', i); j = b if j < 0 else min(j, b)
            return (i, j), j
        (ts, te), i = block(i)
        k = i
        while k < b and t[k] in ' \t\r\n': k += 1
        els = None
        if t.startswith('else', k):
            k += 4
            while k < b and t[k] in ' \t\r\n': k += 1
            if t.startswith('if', k) and re.match(r'if\s*\(', t[k:k + 10]):
                els = ('chain', k, b)
            else:
                (es, ee), _ = block(k); els = ('block', es, ee)
        base = {k_: dict(v) for k_, v in ctx['env'].items()}
        def run(rng_s, rng_e, extra, chain=False):
            c2 = dict(ctx); c2['env'] = {k_: dict(v) for k_, v in base.items()}; c2['conds'] = ctx['conds'] + [extra]
            if chain: self.exec_if(c2, rng_s, rng_e)
            else: self.exec_range(c2, rng_s, rng_e)
            return c2['env']
        env_then = run(ts, te, cond)
        env_else = None
        if els:
            if els[0] == 'chain': env_else = run(els[1], els[2], 'NOT(' + cond + ')', chain=True)
            else: env_else = run(els[1], els[2], 'NOT(' + cond + ')')
        else: env_else = {k_: dict(v) for k_, v in base.items()}
        for var in set(env_then) | set(env_else):
            merged = {}
            for env, lab in ((env_then, cond), (env_else, 'NOT(' + cond + ')')):
                for it, tr in env.get(var, {}).items():
                    if it not in base.get(var, {}):   # value introduced inside this branch
                        tr = tr + (('branch:' + lab, f.rel, f.line(a), (), False),)
                    if it not in merged or better(tr, merged[it]): merged[it] = tr
            ctx['env'][var] = merged

    def exec_switch(self, ctx, a, b):
        """switch (expr) { case X: ... break  case Y: ... default: ... } explored like if / else-if chains: every case is a branch with condition `expr == X`."""
        f = ctx['f']; t = f.text
        po = t.index('(', a); pe = match(t, po)
        expr = ' '.join(t[po + 1:pe - 1].split())
        i = pe
        while i < b and t[i] in ' \t\r\n': i += 1
        if i >= b or t[i] != '{': return
        bs, bend = i + 1, match(t, i) - 1
        body = blank_strings(t[bs:bend]); marks = []; depth = 0; j = 0
        while j < len(body):
            ch = body[j]
            if ch in '({[': depth += 1
            elif ch in ')}]': depth -= 1
            elif depth == 0 and (j == 0 or body[j - 1] in ' \t\r\n;}'):
                m = CASE_RX.match(body, j)
                if m:
                    marks.append((bs + m.start(1), bs + m.end(1), bs + m.end())); j = m.end(); continue
            j += 1
        groups, pending = [], []
        for n_, (ls, le, cs) in enumerate(marks):
            ce = marks[n_ + 1][0] if n_ + 1 < len(marks) else bend
            lab = ' '.join(t[ls:le].split()); pending.append(lab)
            if t[cs:ce].strip():      # a non-empty body closes the group (empty body = fall-through label)
                groups.append((pending, cs, ce)); pending = []
        if pending: groups.append((pending, bend, bend))
        base = {k_: dict(v) for k_, v in ctx['env'].items()}
        envs = []
        for labs, cs, ce in groups:
            conds = [expr + ' == ' + l[4:].strip() if l.startswith('case') else 'default' for l in labs]
            cond = ' || '.join(conds) if len(conds) > 1 else conds[0]
            c2 = dict(ctx); c2['env'] = {k_: dict(v) for k_, v in base.items()}; c2['conds'] = ctx['conds'] + [cond]
            self.exec_range(c2, cs, ce)
            envs.append((cond, c2['env']))
        if not any(l == 'default' for labs, _, _ in groups for l in labs): envs.append(('NOT(' + expr + ' matches a case)', {k_: dict(v) for k_, v in base.items()}))
        for var in set().union(*[set(e) for _, e in envs]) if envs else []:
            merged = {}
            for lab, env in envs:
                for it, tr in env.get(var, {}).items():
                    if it not in base.get(var, {}): tr = tr + (('branch:' + lab, f.rel, f.line(a), (), False),)
                    if it not in merged or better(tr, merged[it]): merged[it] = tr
            ctx['env'][var] = merged

    # ---- expression evaluation
    @staticmethod
    def closure_masks(txt):
        masks = []
        for m in re.finditer(r'\{\s*([\w\s,]*?)\s*->', txt):
            e = match(txt, m.start())
            names = set(x for x in re.split(r'[\s,]+', m.group(1)) if x)
            if names: masks.append((m.start(), e, names))
        return masks

    def add_step(self, prov, step):
        out = {}
        for it, tr in prov.items():
            out[it] = (tr + (step,))[-14:]
        return out

    def eval_expr(self, ctx, a, b, assign_var):
        f = ctx['f']; raw = f.text[a:b]
        txt = blank_strings(raw)
        masks = self.closure_masks(txt)
        ops = re.findall(r'\.(\w+)\s*[\(\{]', txt)
        unused_param = bool(re.search(r'\{[^{}]*?\b_\w+[^{}]*?->', txt))
        prov = self.scan(ctx, txt, a, masks)
        if not prov: return prov
        step = (assign_var or '<expr>', f.rel, f.line(a), tuple(sorted(set(o for o in ops if o not in ('out',)))), unused_param)
        return self.add_step(prov, step)

    def scan(self, ctx, txt, base, masks):
        f = ctx['f']; env = ctx['env']; names = ctx['names']
        res = {}
        # topic channels used as DATA channels: channel.topic('x') = every process output declared `topic: x` (instantiated so far).
        # The 'versions' topic is bookkeeping, never data.
        raw = f.text[base:base + len(txt)]
        for mt in TOPIC_RX.finditer(raw):
            tn = mt.group(1)
            if tn == 'versions': continue
            found_ = {}
            for inst_, info_ in self.instances.items():
                if info_['kind'] != 'process': continue
                for k_, o_ in enumerate(info_['outputs']):
                    if o_.get('topic_name') == tn: found_[('port', inst_, o_['emit'] or f'#{k_}')] = inst_
            # a topic is global in Nextflow, but pipelines read it right after running the producing tool inside the same (sub)workflow instance:
            # prefer producers under the current call path (avoids linking the WGTS and TARGETED runs of the same tool); otherwise all producers
            local_ = {k: v for k, v in found_.items() if ctx['path'] and v.startswith(ctx['path'] + ':')}
            for k_ in (local_ or found_): res.setdefault(k_, ())
        kspans = key_only_spans(txt)
        i = 0
        def masked(pos, name):
            return any(s <= pos < e and name in nm for s, e, nm in masks)
        while i < len(txt):
            m = IDENT.search(txt, i)
            if not m: break
            name = m.group(0); st = m.start(); en = m.end()
            if st > 0 and txt[st - 1] in '.@$' or (st > 0 and txt[st - 1].isalnum()):
                i = en; continue
            if masked(st, name): i = en; continue
            # ----- call of process / workflow
            if name in names:
                k = en
                while k < len(txt) and txt[k] in ' \t': k += 1
                head = None
                if k < len(txt) and txt[k] == '(':
                    ce = match(txt, k)
                    head = self.do_call(ctx, name, base + k + 1, base + ce - 1, txt, k + 1, ce - 1, base)
                    pos = ce
                else:
                    inst = self.inst_name(ctx, name)
                    head = {('proc' if self.callee(f, name) and self.callee(f, name)[0] == 'process' else 'wf', inst): ()}
                    pos = en
                got, pos = self.consume_chain(ctx, head, txt, pos)
                if any(a_ <= st < b_ for a_, b_ in kspans): got = self.add_step(got, ('keyonly', f.rel, f.line(base + st), (), False))
                self.merge(res, got); i = pos; continue
            if name in env:
                head = env[name]
                got, pos = self.consume_chain(ctx, head, txt, en)
                if any(a_ <= st < b_ for a_, b_ in kspans): got = self.add_step(got, ('keyonly', f.rel, f.line(base + st), (), False))
                self.merge(res, got); i = pos; continue
            i = en
        return res

    @staticmethod
    def merge(dst, src):
        for it, tr in src.items():
            if it not in dst or better(tr, dst[it]): dst[it] = tr

    def inst_name(self, ctx, name):
        return (ctx['path'] + ':' + name) if ctx['path'] else name

    def data_ports(self, inst):
        info = self.instances[inst]
        return [(('port', inst, o['emit'] or f'#{k}'), ()) for k, o in enumerate(info['outputs']) if not o['topic'] and not (o['emit'] or '').startswith('versions')]

    def consume_chain(self, ctx, head, txt, pos):
        """Resolve .out / .emit on process / workflow references in `head`."""
        refs = {it: tr for it, tr in head.items() if it[0] in ('proc', 'wf')}
        plain = {it: tr for it, tr in head.items() if it[0] == 'port'}
        if not refs: return plain, pos
        m = re.compile(r'\s*\.\s*out(?!\w)').match(txt, pos)
        p2 = m.end() if m else pos
        m2 = re.compile(r'\s*\.\s*(\w+)(?!\w)(?!\s*[\(\{])').match(txt, p2)
        em = m2.group(1) if m2 else None
        m3 = re.compile(r'\s*\[\s*(\d+)\s*\]').match(txt, p2)
        out = dict(plain); newpos = pos
        for it, tr in refs.items():
            inst = it[1]; info = self.instances.get(inst)
            if not info: continue
            if it[0] == 'proc':
                dp = dict(self.data_ports(inst))
                if m3 and m:
                    idx = int(m3.group(1)); outs = info['outputs']
                    if idx < len(outs) and not outs[idx]['topic']:
                        out[('port', inst, outs[idx]['emit'] or f'#{idx}')] = tr
                    newpos = max(newpos, m3.end()); continue
                if em and any(k[2] == em for k in dp):          # PROC.out.emit  /  PROC(...).emit
                    out[('port', inst, em)] = tr; newpos = max(newpos, m2.end())
                elif em and m and any(o['emit'] == em for o in info['outputs']):
                    newpos = max(newpos, m2.end())               # versions/topic emit: no data port
                elif m:                                          # PROC.out  (all data ports)
                    for k in dp: out[k] = tr
                    newpos = max(newpos, m.end())
                else:
                    out[it] = tr                                 # bare reference kept as reference
            else:
                reg = self.registry.get(inst, {})
                if em and em in reg and (m or True):
                    for k, v in reg[em].items(): out[k] = v
                    newpos = max(newpos, m2.end())
                elif m:
                    for e_, pv in reg.items():
                        for k, v in pv.items(): out[k] = v
                    newpos = max(newpos, m.end())
                else:
                    out[it] = tr
        return out, newpos

    def do_call(self, ctx, name, a, b, txt, ls, le, base):
        f = ctx['f']
        info = self.callee(f, name)
        inst = self.inst_name(ctx, name)
        if not info:
            return {}
        kind, df, d, orig = info
        argspans = split_top(txt, ls, le)
        argprovs = []
        for (s_, e_) in argspans:
            seg = txt[s_:e_]
            argprovs.append(self.scan(ctx, seg, base + s_, self.closure_masks(seg)))
        line = f.line(base + ls)
        if kind == 'process':
            if inst not in self.instances:
                self.instances[inst] = {'path': inst, 'kind': 'process', 'name': name, 'module': orig,
                                        'module_file': df.rel, 'module_line': d['line'], 'call_file': f.rel, 'call_line': line,
                                        'conds': list(ctx['conds']), 'inputs': d['inputs'], 'outputs': d['outputs']}
            for idx, pv in enumerate(argprovs):
                decl = d['inputs'][idx]['decl'] if idx < len(d['inputs']) else None
                for it, tr in pv.items():
                    if it[0] == 'port': self.add_edge(inst, idx, decl, it, tr, ctx, f, line)
                    elif it[0] == 'proc':
                        for k, _ in self.data_ports(it[1]): self.add_edge(inst, idx, decl, k, tr, ctx, f, line)
                    elif it[0] == 'wf':
                        for _e, pv in self.registry.get(it[1], {}).items():
                            for k2, tr2 in pv.items():
                                if k2[0] == 'port': self.add_edge(inst, idx, decl, k2, tr2, ctx, f, line)
            if len(argprovs) != len(d['inputs']):
                self.unresolved.append({'kind': 'arg-count-mismatch', 'file': f.rel, 'line': line, 'process': inst, 'args': len(argprovs), 'declared': len(d['inputs'])})
            return {('proc', inst): ()}
        # subworkflow
        if inst not in self.instances:
            self.instances[inst] = {'path': inst, 'kind': 'workflow', 'name': name, 'module': orig, 'module_file': df.rel,
                                    'module_line': d['line'], 'call_file': f.rel, 'call_line': line, 'conds': list(ctx['conds']),
                                    'inputs': [], 'outputs': []}
        env = {}
        for k, tk in enumerate(d['takes']):
            env[tk] = argprovs[k] if k < len(argprovs) else {}
        c2 = dict(f=df, path=inst, env=env, conds=list(ctx['conds']), names=self.callable_names(df))
        self.exec_range(c2, *d['body'])
        reg = OrderedDict()
        for a_, b_ in d['emits']:
            s = df.text[a_:b_]
            m = re.match(r'(\w+)\s*=(?!=)\s*', s)
            if m:
                reg[m.group(1)] = self.eval_expr(c2, a_ + m.end(), b_, m.group(1))
            else:
                nm = s.strip().split()[0] if s.strip() else None
                if nm: reg[nm] = dict(c2['env'].get(nm, {}))
        # tag subworkflow crossing
        for nm, pv in reg.items():
            reg[nm] = self.add_step(pv, ('emit:' + inst.split(':')[-1] + '.' + nm, df.rel, d['line'], (), False))
        self.registry[inst] = reg
        return {('wf', inst): ()}

    def add_edge(self, dst, idx, decl, port, trail, ctx, f, line):
        allc = list(ctx['conds']) + [st[0][7:] for st in trail if st[0].startswith('branch:')] + list(self.instances.get(port[1], {}).get('conds', []))
        if infeasible(allc):
            self.infeasible.append({'src': port[1], 'emit': port[2], 'dst': dst, 'call': f'{f.rel}:{line}', 'conditions': allc})
            return
        self.edges.append({'src_instance': port[1], 'src_emit': port[2], 'dst_instance': dst, 'dst_input_index': idx,
                           'dst_input_decl': decl, 'call_file': f.rel, 'call_line': line, 'conds': list(ctx['conds']), 'trail': trail})

# ------------------------------------------------------------------- output
def build(root):
    a = Analyzer(root).run()
    # port declaration lookup
    out_decl = {}
    for p, info in a.instances.items():
        for o in info['outputs']: out_decl[(p, o['emit'])] = o
    ed = {}
    for e in a.edges:
        k = (e['src_instance'], e['src_emit'], e['dst_instance'], e['dst_input_index'])
        if k not in ed: ed[k] = e
    edges = []
    for k, e in ed.items():
        od = out_decl.get((e['src_instance'], e['src_emit']))
        ops = sorted({o for st in e['trail'] for o in st[3]})
        flags = []
        if any(st[4] for st in e['trail']): flags.append('closure-with-unused-param(weak)')
        if is_key(e['trail']): flags.append('key-only(sample key kept, payload dropped)')
        if e['conds']: flags.append('conditional-consumer')
        bconds = [st[0][7:] for st in e['trail'] if st[0].startswith('branch:')]
        if bconds: flags.append('alternative-branch-in-channel-path')
        si = a.instances.get(e['src_instance'])
        if si and si['conds']: flags.append('conditional-producer')
        edges.append({'src_instance': e['src_instance'], 'src_emit': e['src_emit'], 'src_output_decl': od['decl'] if od else None,
                      'dst_instance': e['dst_instance'], 'dst_input_index': e['dst_input_index'], 'dst_input_decl': e['dst_input_decl'],
                      'consumer_call': f"{e['call_file']}:{e['call_line']}",
                      'producer_call': f"{si['call_file']}:{si['call_line']}" if si else None,
                      'conditions_consumer': e['conds'], 'conditions_producer': si['conds'] if si else [],
                      'channel_trail': [{'var': s[0], 'at': f'{s[1]}:{s[2]}', 'operators': list(s[3])} for s in e['trail']],
                      'path_branch_conditions': bconds, 'operators': ops, 'flags': flags})
    return a, edges

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', required=True); ap.add_argument('--out', required=True)
    args = ap.parse_args()
    a, edges = build(args.root)
    os.makedirs(args.out, exist_ok=True)
    json.dump({'instances': list(a.instances.values()), 'edges': edges, 'unresolved': a.unresolved, 'infeasible_edges_dropped': a.infeasible},
              open(os.path.join(args.out, 'static_graph.json'), 'w'), indent=1)
    procs = [i for i in a.instances.values() if i['kind'] == 'process']
    print(len(procs), 'process instances;', len([i for i in a.instances.values() if i['kind'] == 'workflow']), 'workflow instances;', len(edges), 'edges;', len(a.unresolved), 'unresolved;', len(a.infeasible), 'infeasible emit-edges dropped')
