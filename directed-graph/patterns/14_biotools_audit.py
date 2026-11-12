"""Step 14: a bio.tools fix list, from how tools are cited.

Four kinds of problem, each ranked by how many citations or graph links it affects:

  missing_main_paper   a heavily cited paper whose title names a tool ("MZmine 2: ...",
                       "SIRIUS 4: ...", "The PRIDE database ... 2019") is listed by no
                       bio.tools entry. If an entry with that name exists, the paper should
                       be added to it; if not, the tool may have no entry at all.
  lists_other_paper    an entry lists a paper owned by another tool (step 1 rules), and that
                       paper is cited more than the entry's own papers, so the entry gets the
                       other tool's citations (quasr and splicemap list HISAT; hisat2 lists
                       Bowtie and SpliceMap but not HISAT).
  owns_no_paper        the entry lists papers but owns none: a wrapper, sub-tool or platform
                       listing (samtools_sort, the EBI BLAST entries). Its graph links are
                       copies of another tool's. Ranked by links in the all-publications graph.
  near_copy            two entries with near-identical neighbours and weight in the
                       all-publications graph (step 12), e.g. boltz / chai-1 / colabfold.

For missing_main_paper, severity 'high' means the missing paper is cited at least as often
as all the papers the entry lists put together (XCMS, MZmine, PRIDE); 'low' means the
entry already has a well-cited paper and this is an extra version or update.

Citations are counts of articles in the omics corpus citing the paper (any section);
graph links are from the all-publications edge list.

Writes directed-graph/biotools/data/biotools_fix_list.csv. Needs steps 1 and 12.
Run from BH2026-InSoLiTo/:
    python directed-graph/patterns/14_biotools_audit.py
"""

import re
from collections import defaultdict

import duckdb
import pandas as pd

EXPORT = 'directed-graph/europepmc/data/export'
TOOL_WORKS = 'directed-graph/citation-order/data/tool_works.parquet'
PUBS = 'directed-graph/biotools/data/biotools_publications.csv'
RAW_EDGES = 'directed-graph/citation-order/data/all_publications/directed_edges.parquet'
DATA = 'directed-graph/patterns/data'
OUT = 'directed-graph/biotools/data/biotools_fix_list.csv'
MIN_CITES = 50       # papers cited by at least this many corpus articles
AFFIXES = {'ebi', 'rest', 'ws', 'api', 'toolshed', 'wrapper', 'galaxy', 'hgv', 'bioconductor',
           'r', 'py', 'web', 'server', 'online', 'biolib', 'app', 'cli', 'plus', 'tool', 'tools'}
# Lead words that look tool-like but name a field, method or consortium, not a tool.
NOT_TOOLS = {'covid19', 'sarscov2', 'crispr', 'crisprcas9', 'rnaseq', 'chipseq', 'scrnaseq', 'mirna', 'lncrna',
             'circrna', 'mirnas', 'gwas', 'dna', 'rna', 'ngs', 'hiv1', 'mrsa', 'ms2', 'lcms', 'mrna', 'cas9', 'atacseq',
             'pcr', 'qpcr', 'nmr', 'mri', 'ct', 'hla', 'icd10', 'snps', 'cnvs', 'tcga', 'gtex', 'encode',
             'ampk', 'mtor', 'mafld', 'nafld', 'smote', 'evtrack', 'plsregression', 'ff14sb', 'charmm36m', 'scopems'}


def name_forms(entry):
    tokens = [t for t in re.split(r'[^a-z0-9]+', entry.lower()) if t and t not in AFFIXES]
    joined = ''.join(tokens)
    return {f for f in {joined, joined.rstrip('0123456789')} if len(f) >= 3}


def lead_name(title):
    """'MZmine 2: modular ...' -> 'mzmine2'; 'The PRIDE database ...' -> 'pride'."""
    t = title or ''
    m = re.match(r'\s*([A-Za-z][A-Za-z0-9+.-]{1,}(?:\s\d+(?:\.\d+)?)?)\s*[:—–-]\s', t)
    if m:
        word = m.group(1)
        if re.search(r'[A-Z]', word[1:]) or re.search(r'\d', word):
            form = re.sub(r'[^a-z0-9]', '', word.lower().replace('+', 'plus'))
            return form if len(form) >= 3 and form not in NOT_TOOLS else None
    m = re.match(r'\s*(?:\d{4}\s+update of\s+)?the\s+([A-Z][A-Za-z0-9-]{2,})\s+(?:database|resource|knowledgebase|'
                 r'consortium|server|web server|portal)\b', t, re.I)
    if m and re.search(r'[A-Z].*[A-Z]|\d', m.group(1)):
        form = re.sub(r'[^a-z0-9]', '', m.group(1).lower())
        return form if form not in NOT_TOOLS else None
    return None


con = duckdb.connect()
con.sql("SET preserve_insertion_order = false")
works = con.sql(f"""
    SELECT work_id, doi, pmid, year, title, n_citations FROM '{EXPORT}/works.parquet'
    WHERE n_citations >= {MIN_CITES} OR work_id IN (SELECT work_id FROM '{TOOL_WORKS}')
""").df()
cites = dict(zip(works.work_id, works.n_citations))
titles = dict(zip(works.work_id, works.title))
tw = con.sql(f"SELECT DISTINCT biotoolsID, work_id FROM '{TOOL_WORKS}'").df()
listed_by = tw.groupby('work_id').biotoolsID.apply(set).to_dict()
lists = tw.groupby('biotoolsID').work_id.apply(set).to_dict()
all_entries = con.sql(f"SELECT DISTINCT biotoolsID FROM read_csv('{PUBS}')").df().biotoolsID
by_form = defaultdict(set)
for e in all_entries:
    for f in name_forms(e):
        by_form[f].add(e)
po = con.sql(f"SELECT work_id, owners, rule FROM '{DATA}/paper_owner.parquet'").df()
owners = {w: list(o) for w, o in zip(po.work_id, po.owners)}
owned = defaultdict(set)
for w, os_ in owners.items():
    for o in os_:
        owned[o].add(w)
strength = con.sql(f"""SELECT t, sum(n) AS n FROM (SELECT source t, n_articles n FROM '{RAW_EDGES}' UNION ALL
                      SELECT target, n_articles FROM '{RAW_EDGES}') GROUP BY 1""").df().set_index('t').n.to_dict()
rows = []

# ---- 1. heavily cited papers naming a tool, listed by no entry
for w in works[(works.n_citations >= MIN_CITES) & ~works.work_id.isin(listed_by)].itertuples(index=False):
    lead = lead_name(w.title)
    if not lead:
        continue
    cands = by_form.get(lead, set()) | by_form.get(lead.rstrip('0123456789'), set())
    if len(cands) > 4:
        cands = by_form.get(lead, set())
    listed_c = sum(cites.get(x, 0) for c in cands for x in lists.get(c, ()))
    rows.append({'issue': 'missing_main_paper', 'entry': ', '.join(sorted(cands)) or '(no entry with this name)',
                 'severity': 'high' if w.n_citations >= listed_c else 'low',
                 'paper': w.title, 'doi': w.doi, 'pmid': w.pmid, 'year': w.year, 'citations': int(w.n_citations),
                 'graph_links': sum(strength.get(c, 0) for c in cands),
                 'detail': f'title names "{lead}"; the entry\'s listed papers have {listed_c} citations'
                           if cands else f'title names "{lead}"'})

# ---- 2. entries listing a more-cited paper they don't own
for e, ws in lists.items():
    own = owned.get(e, set()) & ws
    if not own:
        continue
    own_c = sum(cites.get(w, 0) for w in own)
    for w in ws - own:
        o = owners.get(w, [])
        if o and cites.get(w, 0) > max(own_c, 20):
            rows.append({'issue': 'lists_other_paper', 'entry': e, 'paper': titles.get(w), 'doi': None, 'pmid': None,
                         'year': None, 'citations': int(cites.get(w, 0)), 'graph_links': strength.get(e, 0),
                         'detail': f'paper belongs to {", ".join(o)}; the entry\'s own papers have {own_c} citations'})

# ---- 3. entries owning nothing, with links in the all-publications graph
for e, ws in lists.items():
    if owned.get(e) or strength.get(e, 0) < 100:
        continue
    o = sorted({x for w in ws for x in owners.get(w, [])})
    rows.append({'issue': 'owns_no_paper', 'entry': e, 'paper': '; '.join(filter(None, (titles.get(w) for w in ws)))[:300],
                 'doi': None, 'pmid': None, 'year': None, 'citations': int(sum(cites.get(w, 0) for w in ws)),
                 'graph_links': strength.get(e, 0),
                 'detail': f'its papers belong to {", ".join(o)}' if o else 'only shared platform papers (Galaxy, EBI, ...)'})

# ---- 4. near-identical copies (step 12)
cp = con.sql(f"SELECT tool, entries, strength, citations_mostly_of FROM '{DATA}/allpub_tools.parquet' WHERE citations_mostly_of IS NOT NULL").df()
for r in cp.itertuples(index=False):
    rows.append({'issue': 'near_copy', 'entry': ', '.join(r.entries), 'paper': None, 'doi': None, 'pmid': None, 'year': None,
                 'citations': None, 'graph_links': r.strength,
                 'detail': f'same neighbours and weight in the graph; shared citations mostly of {r.citations_mostly_of or "unknown"}'})

fix = pd.DataFrame(rows)
# high = the missing paper is cited at least as much as everything the entry lists
fix['severity'] = fix.severity.fillna('high')
key = {'missing_main_paper': 'citations', 'lists_other_paper': 'citations', 'owns_no_paper': 'graph_links', 'near_copy': 'graph_links'}
fix['_k'] = [r[key[r['issue']]] for _, r in fix.iterrows()]
fix['_s'] = fix.severity.map({'high': 0, 'low': 1})
fix = fix.sort_values(['issue', '_s', '_k'], ascending=[True, True, False]).drop(columns=['_k', '_s'])
fix = fix[['issue', 'severity', 'entry', 'paper', 'year', 'doi', 'pmid', 'citations', 'graph_links', 'detail']]
fix.to_csv(OUT, index=False)

pd.set_option('display.width', 230)
pd.set_option('display.max_colwidth', 80)
print(fix.groupby(['issue', 'severity']).size().to_string())
for issue in ('missing_main_paper', 'lists_other_paper', 'owns_no_paper', 'near_copy'):
    print(f'\n{issue}:')
    print(fix[fix.issue == issue].head(25)[['entry', 'paper', 'year', 'citations', 'graph_links', 'detail']].to_string(index=False))
print(f'\nWritten {OUT}')
