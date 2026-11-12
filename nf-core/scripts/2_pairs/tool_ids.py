#!/usr/bin/env python3
"""tool_ids.py - identifiers of the software tools: bio.tools IDs, DOIs, bioconda packages ... in three steps.

  python tool_ids.py collect --catalogue DIR --results RESULTS_DIR [--lookup]    -> tool_identifiers.csv
  python tool_ids.py match   --tools tool_identifiers.csv --jsonl biotools_raw.jsonl  -> fills missing bio.tools IDs from a local bio.tools dump
  python tool_ids.py table   --dir RESULTS_DIR                                   -> tools_with_ids.csv

collect  Sources, in this order: (1) every modules/**/meta.yml in the catalogue (tool name, `identifier: biotools:<id>`, homepage, DOI),
         (2) modules/**/environment.yml (bioconda package), (3) the software names found in RESULTS_DIR/<pipeline>/software_pairs.csv
         (names of old-style pipelines have no meta.yml), (4) with --lookup (needs internet) the bio.tools API: a tool WITHOUT an ID from meta.yml
         gets one only on an EXACT name match. Writes tool_identifiers.csv and biotools_cache.json (API answers, reused on re-runs).

match    Compares every identifier and link of a tool that still has no bio.tools ID with a local bio.tools dump (JSON lines, one record per line):
           doi / pmid / pmcid / rrid   publications and other IDs        (strongest)
           url                         homepage, documentation, links    (a bare hosting site such as github.com never counts)
           name                        tool name, name variants, bioconda package   (weakest)
         Among the records with strong evidence the one whose name also agrees wins; several records and no name tie-break = `ambiguous` (not accepted).
         Updates tool_identifiers.csv in place (backup tool_identifiers.before_jsonl.csv), writes biotools_match_report.csv and biotools_match_summary.csv.

table    tools_with_ids.csv: every tool with at least one identifier, its bio.tools ID or the reason there is none, and its other IDs.
         A homepage alone is a link, not an ID.

Words like `cat`, `gzip`, `find`, `filter`, `align` are shell utilities or operations of a tool, not tools: they never get an ID by name.
A name that is a tool plus an operation (`star_align`) is the tool itself (`star`) when that tool has an ID.
Nothing is executed; only text files are read.
"""
import argparse, collections, concurrent.futures as cf, csv, json, os, re, shutil, sys, time, urllib.parse, urllib.request

try: import yaml
except ImportError: sys.exit('pyyaml missing - run: pip install pyyaml')


# ================================================================== collect
HELPERS = {'htslib', 'coreutils', 'grep', 'gzip', 'lbzip2', 'sed', 'tar', 'pigz', 'python', 'r-base', 'bzip2', 'xz', 'pip', 'conda', 'bash', 'gawk', 'perl', 'findutils', 'zlib'}


def key(name):
    n = re.sub(r'\(.*?\)', '', name or '').lower()
    return re.sub(r'[^a-z0-9]', '', n) or re.sub(r'[^a-z0-9]', '', (name or '').lower())


# shell utilities and reporting wrappers: a name match on bio.tools would be a homonym (e.g. `cat` is not the CAT taxonomy tool), so they are never looked up by name
NO_NAME_LOOKUP = {'cat', 'gzip', 'gunzip', 'tar', 'untar', 'unzip', 'sed', 'awk', 'gawk', 'grep', 'sort', 'cut', 'head', 'tail', 'echo', 'python', 'rbase', 'r', 'perl', 'bash', 'java', 'custom', 'print', 'def', 'zcat', 'xz', 'bzip2', 'pigz', 'wget', 'curl', 'ln', 'cp', 'mv', 'rm', 'mkdir', 'touch', 'cd', 'wc', 'tee', 'join', 'paste'}


# words that name an OPERATION of a tool, not another tool: `star_align`, `samtools_sort`, `bwa_index` are STAR / SAMtools / BWA doing that step
OPS = set('align alignment aligner index indexing build genomegenerate generate sort view merge stats flagstat idxstats quant quantification count counts run call calling mem extract dedup markduplicates mark duplicates trim trimming filter convert fastq bam sam cram vcf prepare reference calculate expression profile pileup coverage depth report summary make fai faidx bgzip tabix stat metrics lint check collapse annotate detect predict classify assemble assembly map mapping split cat'.split())

# ordinary words that are also names of bio.tools records or of sub-commands (`bcftools filter`, `bcftools merge`): never looked up by name, never a parent tool
GENERIC_EXTRA = set('find search get list make generate select plot report check test parse read write load rename replace remove clean combine concat compare create design fix format import export label mask match normalize pack process rank refine sample scan score simulate subset summarize tag transform unpack update validate visualize'.split())
GENERIC = OPS | GENERIC_EXTRA | NO_NAME_LOOKUP


def strip_pkg_prefix(name):
    """bioconda package names carry an ecosystem prefix that bio.tools does not: bioconductor-deseq2 -> deseq2, r-ggplot2 -> ggplot2, perl-bioperl -> bioperl"""
    return re.sub(r'^(bioconductor|r|perl|python|py|ucsc)-', '', (name or '').strip().lower())


def pick(counter):
    return counter.most_common(1)[0][0] if counter else ''


def scan(cat):
    tools = {}
    def rec(k, name):
        return tools.setdefault(k, {'names': collections.Counter(), 'ids': collections.Counter(), 'homepage': collections.Counter(), 'doc': collections.Counter(), 'doi': collections.Counter(),
                                    'conda': collections.Counter(), 'modules': set(), 'pipelines': set()})
    for dp, dn, fn in os.walk(cat):
        dn[:] = [d for d in dn if d not in ('.git', 'tests', 'work')]
        rel = os.path.relpath(dp, cat).replace('\\', '/'); pipe = rel.split('/')[0]
        if 'meta.yml' in fn and '/modules/' in '/' + rel + '/':
            try: m = yaml.safe_load(open(os.path.join(dp, 'meta.yml'), encoding='utf-8')) or {}
            except Exception: continue
            tl = m.get('tools') or []
            if not isinstance(tl, list): continue
            envp = []
            if 'environment.yml' in fn:
                try: envp = [re.sub(r'^.*::', '', x).split('=')[0] for x in (yaml.safe_load(open(os.path.join(dp, 'environment.yml'), encoding='utf-8')) or {}).get('dependencies', []) if isinstance(x, str)]
                except Exception: envp = []
            for t in tl:
                if not isinstance(t, dict): continue
                for name, d in t.items():
                    d = d if isinstance(d, dict) else {}
                    k = key(name)
                    if not k: continue
                    r = rec(k, name); r['names'][str(name)] += 1; r['modules'].add(rel); r['pipelines'].add(pipe)
                    ident = d.get('identifier'); ident = ' '.join(map(str, ident)) if isinstance(ident, list) else str(ident or '')
                    for part in re.findall(r'biotools:([\w.\-]+)', ident): r['ids'][part] += 1
                    if d.get('homepage'): r['homepage'][str(d['homepage']).strip()] += 1
                    if d.get('documentation'): r['doc'][str(d['documentation']).strip()] += 1
                    doi = d.get('doi'); doi = ' '.join(map(str, doi)) if isinstance(doi, list) else str(doi or '')
                    for dd in re.findall(r'10\.\d{4,9}/[^\s,;]+', doi): r['doi'][dd.rstrip('.')] += 1
                    for pkg in envp:
                        if key(pkg) == k and pkg not in HELPERS: r['conda'][pkg] += 1
    return tools


def add_results_names(tools, results):
    n = 0
    for p in sorted(os.listdir(os.path.join(results, 'pipelines'))) if os.path.isdir(os.path.join(results, 'pipelines')) else []:
        f = os.path.join(results, 'pipelines', p, 'software_pairs.csv')
        if not os.path.exists(f): continue
        for row in csv.DictReader(open(f, encoding='utf-8')):
            for nm in (row['source_software'], row['target_software']):
                k = key(nm)
                if not k: continue
                if k not in tools:
                    tools[k] = {'names': collections.Counter(), 'ids': collections.Counter(), 'homepage': collections.Counter(), 'doc': collections.Counter(), 'doi': collections.Counter(),
                                'conda': collections.Counter(), 'modules': set(), 'pipelines': set()}; n += 1
                tools[k]['names'][nm] += 1; tools[k]['pipelines'].add(p.split('__')[0])
    return n


def api(url, tries=3):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={'Accept': 'application/json', 'User-Agent': 'nfrel-tool-ids'})
            with urllib.request.urlopen(req, timeout=25) as r: return json.loads(r.read().decode('utf-8'))
        except urllib.error.HTTPError as e:
            if e.code == 404: return None
        except Exception: pass
        time.sleep(1.5 * (i + 1))
    return 'ERROR'


def lookup(k, name, known_id, cache):
    ck = known_id or ('name:' + k)
    if not known_id and k in GENERIC: return {'found': False}   # `filter`, `merge`, `find`, `align` ... are operation words: a bio.tools record of that name is a homonym
    if ck in cache: return cache[ck]
    rec = None; how = ''
    if known_id:
        rec = api('https://bio.tools/api/tool/' + urllib.parse.quote(known_id) + '/?format=json'); how = 'meta.yml'
    else:
        if k in NO_NAME_LOOKUP: cache[ck] = {'found': False}; return cache[ck]
        base = strip_pkg_prefix(re.sub(r'\(.*?\)', '', name))
        ks = {k, key(base)}   # the key with and without the ecosystem prefix (bioconductor-deseq2 / deseq2)
        for cand in dict.fromkeys([base.replace(' ', '_'), base.replace(' ', '-'), base.replace('_', '-'), base.replace('-', '_'), base.replace(' ', ''), k]):
            if not cand: continue
            r = api('https://bio.tools/api/tool/' + urllib.parse.quote(cand) + '/?format=json')
            if isinstance(r, dict) and r.get('biotoolsID') and key(r['biotoolsID']) in ks: rec = r; how = 'bio.tools-exact-id'; break
        if rec is None:
            r = api('https://bio.tools/api/tool/?q=' + urllib.parse.quote(name) + '&format=json')
            if isinstance(r, dict):
                for it in r.get('list', []):
                    if key(it.get('name', '')) in ks or key(it.get('biotoolsID', '')) in ks: rec = it; how = 'bio.tools-exact-name'; break
    out = {'found': False}
    if isinstance(rec, dict) and rec.get('biotoolsID'):
        pubs = []
        for pb in rec.get('publication') or []:
            if pb.get('doi'): pubs.append('doi:' + str(pb['doi']))
            if pb.get('pmid'): pubs.append('pmid:' + str(pb['pmid']))
        others = ['%s:%s' % (o.get('type', ''), o.get('value', '')) for o in rec.get('otherID') or []]
        out = {'found': True, 'id': rec['biotoolsID'], 'name': rec.get('name', ''), 'how': how, 'homepage': rec.get('homepage', ''), 'other_ids': others, 'publications': list(dict.fromkeys(pubs))}
    cache[ck] = out
    return out


def resolve_parents(rows):
    """names that are a known tool plus an operation (star_align -> star) get parent_key / parent_biotools_id / operation_suffix / parent_confidence"""
    with_id = {r['tool_key']: r['biotools_id'] for r in rows if r['biotools_id'] and r['tool_key'] not in GENERIC and len(r['tool_key']) >= 3}
    for r in rows:
        r.update({'parent_key': '', 'parent_biotools_id': '', 'operation_suffix': '', 'parent_confidence': ''})
        if r['biotools_id']: continue
        toks = [t for t in re.split(r'[_\-\s]+', re.sub(r'\(.*?\)', '', r['tool_name']).lower()) if t]
        if len(toks) < 2: continue
        par = key(toks[0])
        if par in with_id and par != r['tool_key']:
            rest = toks[1:]
            r.update({'parent_key': par, 'parent_biotools_id': with_id[par], 'operation_suffix': ' '.join(rest), 'parent_confidence': 'operation-suffix' if all(t in OPS or t.isdigit() for t in rest) else 'first-word-review'})


def collect(a):
    os.makedirs(a.out, exist_ok=True)
    tools = scan(a.catalogue); n_meta = len(tools)
    extra = add_results_names(tools, a.results) if a.results and os.path.isdir(a.results) else 0
    print(f'{n_meta} tools from meta.yml; {extra} more software names from the results')
    cpath = os.path.join(a.out, 'biotools_cache.json')
    cache = json.load(open(cpath, encoding='utf-8')) if os.path.exists(cpath) else {}
    info = {}
    if a.lookup:
        todo = [(k, pick(t['names']), pick(t['ids'])) for k, t in sorted(tools.items())]
        print(f'bio.tools lookup for {len(todo)} tools (cached: {sum(1 for k, n, i in todo if (i or "name:" + k) in cache)}) ...', flush=True)
        with cf.ThreadPoolExecutor(max_workers=a.jobs) as ex:
            for i, ((k, nm, idd), res) in enumerate(zip(todo, ex.map(lambda x: lookup(x[0], x[1], x[2], cache), todo)), 1):
                info[k] = res
                if i % 100 == 0: print(f'  {i}/{len(todo)}', flush=True); json.dump(cache, open(cpath, 'w', encoding='utf-8'))
        json.dump(cache, open(cpath, 'w', encoding='utf-8'), indent=0)
    rows = []
    for k, t in sorted(tools.items()):
        res = info.get(k, {}); meta_id = pick(t['ids'])
        bid = meta_id or (res.get('id', '') if res.get('found') else '')
        src = 'meta.yml' if meta_id else (res.get('how', '') if res.get('found') else ('none' if a.lookup else 'none (run with --lookup)'))
        rows.append({'tool_key': k, 'tool_name': pick(t['names']), 'biotools_id': bid, 'biotools_url': ('https://bio.tools/' + bid) if bid else '', 'id_source': src,
                     'bioconda_package': pick(t['conda']), 'doi': ';'.join(d for d, _ in t['doi'].most_common(3)), 'homepage': pick(t['homepage']) or res.get('homepage', ''),
                     'documentation': pick(t['doc']), 'other_ids': ';'.join(res.get('other_ids', [])), 'publications': ';'.join(res.get('publications', [])),
                     'n_modules': len(t['modules']), 'n_pipelines': len(t['pipelines']), 'name_variants': ' | '.join(n for n, _ in t['names'].most_common(4))})
    resolve_parents(rows)
    rows.sort(key=lambda r: (-r['n_pipelines'], r['tool_key']))
    f = os.path.join(a.out, 'tool_identifiers.csv')
    with open(f, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else ['tool_key'], lineterminator='\n'); w.writeheader(); w.writerows(rows)
    c = collections.Counter(r['id_source'] for r in rows)
    print(f'{len(rows)} tools -> {f}\nwith a bio.tools ID: {sum(1 for r in rows if r["biotools_id"])} | by source: {dict(c)}')
    print('names resolved to a parent tool (e.g. star_align -> star):', sum(1 for r in rows if r['parent_confidence'] == 'operation-suffix'), '| first-word matches to review:', sum(1 for r in rows if r['parent_confidence'] == 'first-word-review'))
    print('with a bioconda package:', sum(1 for r in rows if r['bioconda_package']), '| with a DOI:', sum(1 for r in rows if r['doi'] or r['publications']))




# ================================================================== match
DOI_RX = re.compile(r'10\.\d{4,9}/[^\s"<>,;]+', re.I)
STRONG = ('doi', 'pmid', 'pmcid', 'rrid')
ORDER = ('doi', 'pmid', 'pmcid', 'rrid', 'url', 'name')


GENERIC_HOSTS = {'github.com', 'gitlab.com', 'bitbucket.org', 'sourceforge.net', 'bioconductor.org', 'cran.r-project.org', 'r-project.org', 'pypi.org', 'anaconda.org', 'bioconda.github.io', 'hub.docker.com', 'docker.com', 'zenodo.org', 'figshare.com', 'ncbi.nlm.nih.gov', 'ebi.ac.uk', 'galaxyproject.org', 'usegalaxy.org', 'nf-co.re', 'doi.org', 'google.com', 'sites.google.com', 'github.io', 'readthedocs.io', 'readthedocs.org'}


def norm_doi(x):
    x = (x or '').strip().lower()
    x = re.sub(r'^(https?://(dx\.)?doi\.org/|doi:\s*)', '', x)
    m = DOI_RX.search(x)
    return m.group(0).rstrip('.)]}') if m else ''


def norm_url(u):
    """host + path, lower case, without scheme / www / trailing slash / .git / index page / query; '' for a bare host"""
    u = (u or '').strip().lower()
    if not u: return ''
    u = re.sub(r'^[a-z]+://', '', u)
    u = re.sub(r'^www\.', '', u)
    u = re.split(r'[?#]', u)[0].rstrip('/')
    u = re.sub(r'(/index\.(html?|php))$', '', u)
    u = re.sub(r'\.git$', '', u)
    if '/' not in u and u in GENERIC_HOSTS: return ''   # a bare hosting site (github.com ...) identifies nothing; a bare tool site (htslib.org) does
    return u


def norm_rrid(x):
    m = re.search(r'scr_\d+', (x or '').lower())
    return m.group(0) if m else ''


def records(obj):
    """yield every dict that has a biotoolsID, wherever it is nested"""
    if isinstance(obj, dict):
        if obj.get('biotoolsID'): yield obj; return
        for v in obj.values(): yield from records(v)
    elif isinstance(obj, list):
        for v in obj: yield from records(v)


def record_ids(rec):
    """-> {kind: set(values)} for doi / pmid / pmcid / rrid / url"""
    d = collections.defaultdict(set)
    for p in rec.get('publication') or []:
        if not isinstance(p, dict): continue
        if norm_doi(p.get('doi')): d['doi'].add(norm_doi(p.get('doi')))
        if p.get('pmid'): d['pmid'].add(str(p['pmid']).strip())
        if p.get('pmcid'): d['pmcid'].add(str(p['pmcid']).strip().lower())
    for o in rec.get('otherID') or []:
        if not isinstance(o, dict): continue
        v = str(o.get('value', ''))
        if norm_doi(v): d['doi'].add(norm_doi(v))
        if norm_rrid(v): d['rrid'].add(norm_rrid(v))
    for k in ('link', 'documentation'):
        for x in rec.get(k) or []:
            if isinstance(x, dict):
                if norm_url(x.get('url')): d['url'].add(norm_url(x['url']))
                if norm_doi(x.get('url')): d['doi'].add(norm_doi(x['url']))
    for u in (rec.get('homepage'),):
        if norm_url(u): d['url'].add(norm_url(u))
        if norm_doi(u): d['doi'].add(norm_doi(u))
    return d


def tool_ids_of(row):
    d = collections.defaultdict(set)
    for col in ('doi', 'publications', 'other_ids'):
        for x in re.split(r'[;|]', row.get(col, '') or ''):
            x = x.strip()
            if not x: continue
            if norm_doi(x): d['doi'].add(norm_doi(x))
            if x.lower().startswith('pmid:'): d['pmid'].add(x[5:].strip())
            if x.lower().startswith('pmcid:'): d['pmcid'].add(x[6:].strip().lower())
            if norm_rrid(x): d['rrid'].add(norm_rrid(x))
    for col in ('homepage', 'documentation'):
        for x in re.split(r'[;|\s]+', row.get(col, '') or ''):
            if norm_url(x): d['url'].add(norm_url(x))
    return d


def match(a):
    out = a.out or os.path.dirname(os.path.abspath(a.tools))
    info = {}; idx = {k: collections.defaultdict(set) for k in STRONG + ('url',)}; by_name = collections.defaultdict(set); n_lines = 0; bad = 0
    with open(a.jsonl, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            line = line.strip()
            if not line: continue
            n_lines += 1
            try: obj = json.loads(line)
            except Exception: bad += 1; continue
            for rec in records(obj):
                bid = rec['biotoolsID']
                info[bid] = {'name': rec.get('name', ''), 'homepage': rec.get('homepage', '')}
                for kind, vals in record_ids(rec).items():
                    for v in vals: idx[kind][v].add(bid)
                for k in {key(rec.get('name', '')), key(bid)}:
                    if k: by_name[k].add(bid)
    print(f'bio.tools dump: {n_lines} lines ({bad} unreadable), {len(info)} tool records; index sizes: ' + ', '.join(f'{k} {len(v)}' for k, v in idx.items()))
    if not info: sys.exit('no records with a biotoolsID found - is this the right file? first characters: ' + open(a.jsonl, encoding='utf-8', errors='replace').read(200))
    rows = list(csv.DictReader(open(a.tools, encoding='utf-8'))); fields = list(rows[0].keys())
    bk = os.path.join(out, 'tool_identifiers.before_jsonl.csv')
    if not os.path.exists(bk): shutil.copy(a.tools, bk)   # kept from the FIRST run, so a re-run never overwrites the real 'before'
    for c in ('match_evidence', 'match_doi', 'match_candidates'):
        if c not in fields: fields.append(c)
    report = []; counts = collections.Counter()
    for r in rows:
        r.setdefault('match_evidence', ''); r.setdefault('match_doi', ''); r.setdefault('match_candidates', '')
        if r['biotools_id'] or r['tool_key'] in NO_NAME_LOOKUP: continue
        names = {r['tool_key']} | {key(n) for n in (r.get('name_variants') or '').split('|') if n.strip()} | {key(strip_pkg_prefix(n)) for n in (r.get('name_variants') or '').split('|') + [r.get('bioconda_package', '')] if n.strip()} | {key(r.get('bioconda_package', ''))}
        names.discard('')
        mine = tool_ids_of(r)
        cand = collections.defaultdict(lambda: collections.defaultdict(set))   # biotoolsID -> kind -> matched values
        for kind in STRONG + ('url',):
            for v in mine.get(kind, ()):
                for bid in idx[kind].get(v, ()): cand[bid][kind].add(v)
        for bid in set(cand) | set().union(*[by_name.get(k, set()) for k in names]):
            if {key(info[bid]['name']), key(bid)} & names: cand[bid]['name'].add('name')
        def pick(kinds):
            pool = {b: k for b, k in cand.items() if set(k) & set(kinds)}
            if not pool: return None, ''
            named = [b for b, k in pool.items() if 'name' in k]
            if len(named) == 1: return named[0], ''
            if len(named) > 1: return None, 'ambiguous'
            if len(pool) == 1: return next(iter(pool)), ''
            return None, 'ambiguous'
        chosen, amb = pick(STRONG)
        if not chosen and not amb: chosen, amb = pick(('url',))
        generic_only = all(k in GENERIC for k in names)   # a name like `filter` or `merge` alone is never enough
        if not chosen and not amb and not a.no_name_only and not generic_only:
            nm = {b: k for b, k in cand.items() if 'name' in k}
            if len(nm) == 1: chosen = next(iter(nm))
            elif len(nm) > 1: amb = 'ambiguous'
        evidence = ''
        if chosen: evidence = '+'.join(k for k in ORDER if k in cand[chosen])
        elif amb: evidence = 'ambiguous'
        counts[evidence or 'no-match'] += 1
        if chosen:
            r['biotools_id'] = chosen; r['biotools_url'] = 'https://bio.tools/' + chosen; r['id_source'] = 'biotools-jsonl-' + evidence
            r['match_evidence'] = evidence; r['match_doi'] = ';'.join(sorted(cand[chosen].get('doi', set())))
            if not r.get('homepage'): r['homepage'] = info[chosen]['homepage']
        if evidence:
            r['match_candidates'] = ';'.join(sorted(cand))[:300]
            report.append({'tool_key': r['tool_key'], 'tool_name': r['tool_name'], 'evidence': evidence, 'accepted': bool(chosen), 'chosen_biotools_id': chosen or '', 'chosen_name': info[chosen]['name'] if chosen else '',
                           'candidates': ';'.join(f"{b} ({info[b]['name']}: {'+'.join(k for k in ORDER if k in cand[b])})" for b in sorted(cand))[:500], 'tool_ids_and_links': ';'.join(f'{k}={sorted(v)[0]}' for k, v in mine.items() if v), 'name_variants': r.get('name_variants', '')})
    for r in rows:
        for c in ('parent_key', 'parent_biotools_id', 'operation_suffix', 'parent_confidence'): r[c] = ''
    resolve_parents(rows)   # a tool that just got an ID can be the parent of operation-suffixed names (star_align -> star)
    for c in ('parent_key', 'parent_biotools_id', 'operation_suffix', 'parent_confidence'):
        if c not in fields: fields.append(c)
    with open(os.path.join(out, 'tool_identifiers.csv'), 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction='ignore', lineterminator='\n'); w.writeheader(); w.writerows(rows)
    rf = os.path.join(out, 'biotools_match_report.csv')
    with open(rf, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=['tool_key', 'tool_name', 'evidence', 'accepted', 'chosen_biotools_id', 'chosen_name', 'candidates', 'tool_ids_and_links', 'name_variants'], lineterminator='\n'); w.writeheader(); w.writerows(report)
    was = lambda r: bool(r['biotools_id']) and not r['id_source'].startswith('biotools-jsonl')   # had an ID before any dump matching (meta.yml or API lookup)
    npipe = lambda r: int(r.get('n_pipelines') or 0)
    summ = [('tools in the table', len(rows), sum(npipe(r) for r in rows)), ('with a bio.tools ID BEFORE matching', sum(1 for r in rows if was(r)), sum(npipe(r) for r in rows if was(r))),
            ('with a bio.tools ID AFTER matching', sum(1 for r in rows if r['biotools_id']), sum(npipe(r) for r in rows if r['biotools_id'])),
            ('gained by matching', sum(1 for r in rows if r['biotools_id'] and not was(r)), sum(npipe(r) for r in rows if r['biotools_id'] and not was(r))),
            ('still without an ID', sum(1 for r in rows if not r['biotools_id']), sum(npipe(r) for r in rows if not r['biotools_id']))]
    for ev, n in collections.Counter(r['id_source'][len('biotools-jsonl-'):] for r in rows if r['id_source'].startswith('biotools-jsonl-')).most_common(): summ.append(('  gained via ' + ev, n, sum(npipe(r) for r in rows if r['id_source'] == 'biotools-jsonl-' + ev)))
    with open(os.path.join(out, 'biotools_match_summary.csv'), 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh, lineterminator='\n'); w.writerow(['metric', 'n_tools', 'sum_of_n_pipelines (a tool used in 10 pipelines counts 10)', 'share_of_tools']); [w.writerow([m, n, pw, f'{n / max(1, len(rows)):.1%}']) for m, n, pw in summ]
    print('\nBEFORE vs AFTER matching (tools | tool-uses over pipelines | share of tools)')
    for m, n, pw in summ: print(f'  {m:40s} {n:6d} | {pw:7d} | {n / max(1, len(rows)):6.1%}')
    print('\ntools without an ID before:', sum(counts.values()), '| by evidence:', dict(counts.most_common()))
    print('now with a bio.tools ID:', sum(1 for r in rows if r['biotools_id']), 'of', len(rows))
    print('updated:', os.path.join(out, 'tool_identifiers.csv'), '| report (review `ambiguous` and name-only rows):', rf)




# ================================================================== table
UTIL = {'cat', 'gzip', 'gunzip', 'tar', 'untar', 'unzip', 'sed', 'awk', 'gawk', 'grep', 'sort', 'cut', 'head', 'tail', 'echo', 'python', 'rbase', 'r', 'perl', 'bash', 'java', 'custom', 'zcat', 'xz', 'bzip2', 'pigz', 'wget', 'curl', 'ln', 'cp', 'mv', 'rm', 'mkdir', 'touch', 'wc', 'tee', 'join', 'paste'}


def others(r):
    out = []
    if r.get('bioconda_package'): out.append('bioconda:' + r['bioconda_package'])
    for x in re.split(r'[;|]', (r.get('doi', '') or '') + ';' + (r.get('publications', '') or '') + ';' + (r.get('other_ids', '') or '')):
        x = x.strip()
        if not x: continue
        m = re.search(r'10\.\d{4,9}/\S+', x)
        if m: out.append('doi:' + m.group(0).rstrip('.)]}'))
        elif x.lower().startswith(('pmid:', 'pmcid:')): out.append(x.lower())
        elif 'scr_' in x.lower(): out.append('rrid:' + re.search(r'SCR_\d+', x, re.I).group(0).upper())
        elif ':' in x and not x.lower().startswith('doi:'): out.append(x)
    return list(dict.fromkeys(out))


def origin(r):
    s = r.get('id_source', '')
    if s == 'meta.yml': return 'meta.yml (nf-core module metadata)'
    if s.startswith('bio.tools-exact'): return 'bio.tools API (' + s.replace('bio.tools-', '') + ')'
    if s.startswith('biotools-jsonl-'): return 'bio.tools dump (' + s[len('biotools-jsonl-'):] + ')'
    return s



def table(a):
    tf = os.path.join(a.dir, 'tool_identifiers.csv'); rf = os.path.join(a.dir, 'biotools_match_report.csv')
    if not os.path.exists(tf): sys.exit('tool_identifiers.csv not found in ' + os.path.abspath(a.dir) + ' - run: python run.py ids')
    rows = list(csv.DictReader(open(tf, encoding='utf-8')))
    report = {x['tool_key']: x for x in csv.DictReader(open(rf, encoding='utf-8'))} if os.path.exists(rf) else None
    out = []
    for r in rows:
        oth = others(r)
        if not (r['biotools_id'] or oth): continue
        reason = ''
        if not r['biotools_id']:
            k = r['tool_key']
            if k in UTIL: reason = 'shell utility, deliberately not looked up on bio.tools (a name match would be a homonym)'
            elif r.get('parent_confidence') == 'operation-suffix': reason = f"name is an operation of the tool '{r['parent_key']}' ({r['operation_suffix']}); the merge uses the parent's ID"
            elif r.get('parent_confidence') == 'first-word-review': reason = f"name starts with the tool '{r['parent_key']}' but '{r['operation_suffix']}' is not a known operation - needs review"
            elif report is None: reason = 'matching against the bio.tools dump has not been run yet (run: python run.py match)'
            elif k in report and report[k]['evidence'] == 'ambiguous': reason = 'ambiguous: several bio.tools records fit - ' + report[k]['candidates'][:200]
            elif k in report: reason = 'candidate found but not accepted (' + report[k]['evidence'] + ')'
            else: reason = 'no record in bio.tools matched this tool (searched by DOI, PubMed ID, RRID, link and name)'
        out.append({'tool': r['tool_name'], 'tool_key': r['tool_key'], 'biotools_id': r['biotools_id'], 'biotools_url': r['biotools_url'], 'biotools_id_origin': origin(r) if r['biotools_id'] else '',
                    'reason_no_biotools_id': reason, 'biotools_id_via_parent': r.get('parent_biotools_id', '') if (not r['biotools_id'] and r.get('parent_confidence') == 'operation-suffix') else '',
                    'other_ids_available': ' ; '.join(oth), 'n_other_ids': len(oth), 'n_pipelines': r.get('n_pipelines', '')})
    out.sort(key=lambda x: (x['biotools_id'] == '', -int(x['n_pipelines'] or 0), x['tool_key']))
    f = os.path.join(a.dir, 'tools_with_ids.csv')
    with open(f, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()) if out else ['tool'], lineterminator='\n'); w.writeheader(); w.writerows(out)
    nb = sum(1 for x in out if x['biotools_id'])
    print(f'{len(out)} tools with at least one identifier of {len(rows)} -> {f}')
    print(f'  with a bio.tools ID: {nb} | without one but with another ID: {len(out) - nb}')
    print('  reasons for no bio.tools ID:', dict(collections.Counter((x['reason_no_biotools_id'].split(' - ')[0].split(':')[0][:60]) for x in out if not x['biotools_id']).most_common()))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest='cmd', required=True)
    p = sp.add_parser('collect'); p.add_argument('--catalogue', required=True); p.add_argument('--results'); p.add_argument('--out', default='.'); p.add_argument('--lookup', action='store_true'); p.add_argument('--jobs', type=int, default=4)
    p = sp.add_parser('match'); p.add_argument('--tools', required=True); p.add_argument('--jsonl', required=True); p.add_argument('--out'); p.add_argument('--no-name-only', action='store_true', help='accept matches that rest on identifiers or links only, never on a name alone')
    p = sp.add_parser('table'); p.add_argument('--dir', default='.')
    a = ap.parse_args()
    {'collect': collect, 'match': match, 'table': table}[a.cmd](a)


if __name__ == '__main__':
    main()
