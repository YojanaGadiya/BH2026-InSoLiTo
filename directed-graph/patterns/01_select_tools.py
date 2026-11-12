"""Step 1: pick the most-used tools and build the tables the pattern analyses read.

"Use" = a research article cites one of the tool's papers in its methods section (the
same articles and section as the directed graph).

bio.tools entries share papers in ways that would count one citation for many tools:
Galaxy and EBI wrapper entries list the Galaxy or EBI platform paper, sub-entries list
their parent tool's paper (samtools_sort -> SAMtools), and some entries list papers of
tools they use. So every cited paper is first given an owner:

  0. a title that opens with a tool's name ("HISAT: a fast spliced aligner") names its
     tool: if that is a different bio.tools entry from the ones listing the paper, the
     paper goes to that entry (bio.tools' hisat2 record lists QuasR's papers, and QuasR
     and SpliceMap list the HISAT paper);
  1. a paper listed by one entry belongs to it;
  2. a paper listed by several entries belongs to those whose name is in its title
     ("STAR: ultrafast universal RNA-seq aligner" -> star);
  3. otherwise, if at most 3 entries list it, to those that have it as their main
     publication (at most 2);
  4. otherwise it is a shared platform paper (Galaxy, EBI tools) and is left out.

Entries that then own the same papers are merged, an entry whose papers are a subset
of a more-used entry's is merged into it (mafft_ebi -> MAFFT), and so is an entry whose
users are at least 80% the same as a more-used entry's that shares a paper with it
(hgv_david -> DAVID). A tool is named after its
most-used entry; `entries` lists all of them.

Writes to directed-graph/patterns/data/:
  tools.parquet        one row per tool: tool, entries, n_articles, rank
  tool_use.parquet     one row per (article, tool): pmcid, tool, year, topics
  tool_papers.parquet  one row per (tool, paper): work_id, doi, pmid, pmcid, year, title, ...
  paper_owner.parquet  every matched paper, its owner tools and the rule that decided

Needs tool_works.parquet and kept_articles.parquet from citation-order/ (steps 01-02),
the export tables, manifest.tsv and the bio.tools main publication files.
Run from BH2026-InSoLiTo/:
    python directed-graph/patterns/01_select_tools.py
"""

import re
from collections import defaultdict
from pathlib import Path

import duckdb
import pandas as pd

TOP_N = 300

EXPORT = 'directed-graph/europepmc/data/export'
MANIFEST = 'directed-graph/europepmc/data/manifest.tsv'
TOOL_WORKS = 'directed-graph/citation-order/data/tool_works.parquet'
KEPT = 'directed-graph/citation-order/data/kept_articles.parquet'
MAIN_UNIQUE = 'directed-graph/biotools/data/biotools_main_publications_unique.csv'
MAIN_SHARED = 'directed-graph/biotools/data/biotools_main_publications_shared.csv'
ALIASES = 'directed-graph/europepmc/data/export/doi_aliases.parquet'
OUT = 'directed-graph/patterns/data'

# Name parts that mark a wrapper or deployment rather than the tool itself.
AFFIXES = {'ebi', 'rest', 'ws', 'api', 'toolshed', 'wrapper', 'galaxy', 'hgv', 'bioconductor',
           'r', 'py', 'web', 'server', 'online', 'biolib', 'app', 'cli', 'plus', 'tool', 'tools'}


def name_forms(entry: str) -> set[str]:
    """How an entry's name could appear in a title: 'hgv_david' -> {'david'}, 'hisat2' -> {'hisat2', 'hisat'}."""
    tokens = [t for t in re.split(r'[^a-z0-9]+', entry.lower()) if t and t not in AFFIXES]
    joined = ''.join(tokens)
    forms = {joined, joined.rstrip('0123456789')}
    return {f for f in forms if len(f) >= 3}


def lead_name(title: str) -> str | None:
    """'HISAT: a fast...' -> 'hisat'. Only tool-like names: one word with an inner capital
    or a digit (HISAT, QuasR, DESeq2), so 'Metabolomics: ...' is not taken for a tool."""
    m = re.match(r'\s*([A-Za-z][A-Za-z0-9+.-]{2,})\s*:', title or '')
    if not m:
        return None
    word = m.group(1)
    if not (re.search(r'[A-Z]', word[1:]) or re.search(r'\d', word)):
        return None
    form = re.sub(r'[^a-z0-9]', '', word.lower().replace('+', 'plus'))  # SCENIC+ -> scenicplus
    return form if len(form) >= 4 else None


def related(entry: str, lead: str) -> bool:
    """Is this entry the titled tool or a version/wrapper of it? (pecan_dia for PECAN)"""
    tokens = set(re.split(r'[^a-z0-9]+', entry.lower()))
    return lead in tokens or any(f.startswith(lead) for f in name_forms(entry) | {re.sub(r'[^a-z0-9]', '', entry.lower())})


def title_forms(title: str) -> set[str]:
    """Title words, and runs of 2-3 words joined ('Bowtie 2' -> 'bowtie2')."""
    words = re.findall(r'[a-z0-9]+', (title or '').lower())
    out = set(words) | {w.rstrip('0123456789') for w in words}
    for n in (2, 3):
        out |= {''.join(words[i:i + n]) for i in range(len(words) - n + 1)}
    return out


Path(OUT).mkdir(parents=True, exist_ok=True)
con = duckdb.connect()

# --- Which papers does each entry list, and which are its main publications? ---
entries = con.sql(f"SELECT DISTINCT biotoolsID, work_id FROM '{TOOL_WORKS}'").df()
works = con.sql(f"""
    SELECT work_id, doi, pmid, pmcid, year, title, first_author, n_citations
    FROM '{EXPORT}/works.parquet'
    WHERE work_id IN (SELECT work_id FROM '{TOOL_WORKS}')
""").df()
alias = f"(SELECT preprint_doi AS doi, canonical_doi AS canonical FROM '{ALIASES}')" if Path(ALIASES).exists() \
    else "(SELECT NULL::VARCHAR AS doi, NULL::VARCHAR AS canonical WHERE false)"
main = con.sql(f"""
    WITH m AS (
        SELECT biotoolsID, lower(doi) AS doi, CAST(pmid AS VARCHAR) AS pmid FROM '{MAIN_UNIQUE}'
        UNION
        SELECT unnest(string_split(biotoolsIDs, ';')), lower(doi), CAST(pmid AS VARCHAR) FROM '{MAIN_SHARED}'
    ), w AS (
        SELECT w.work_id, coalesce(a.canonical, w.doi) AS doi, CAST(w.pmid AS VARCHAR) AS pmid
        FROM '{EXPORT}/works.parquet' w LEFT JOIN {alias} a ON a.doi = w.doi
        WHERE w.work_id IN (SELECT work_id FROM '{TOOL_WORKS}')
    )
    SELECT DISTINCT m.biotoolsID, w.work_id FROM m LEFT JOIN {alias} a ON a.doi = m.doi
    JOIN w ON w.doi = coalesce(a.canonical, m.doi) OR (m.pmid IS NOT NULL AND w.pmid = m.pmid)
""").df()
main_of = defaultdict(set)
for r in main.itertuples():
    main_of[r.work_id].add(r.biotoolsID)
title = dict(zip(works.work_id, works.title))

# --- Rules 0-4: an owner for every paper ---
listed = entries.groupby('work_id').biotoolsID.apply(set)
all_entries = con.sql("SELECT DISTINCT biotoolsID FROM read_csv('directed-graph/biotools/data/biotools_publications.csv')").df().biotoolsID
by_form = defaultdict(set)
for e in all_entries:
    for f in name_forms(e):
        by_form[f].add(e)
owners, rule = {}, {}
for work_id, cands in listed.items():
    lead = lead_name(title.get(work_id))
    exact = by_form.get(lead, set()) if lead else set()
    if exact and len(exact) <= 3 and not exact <= cands:
        # The title names an entry that doesn't list this paper: give it the paper, keeping
        # listing entries that are versions or wrappers of the same tool (merged below).
        owners[work_id] = set(exact) | {c for c in cands if related(c, lead)}
        rule[work_id] = 'title_names_entry'
        continue
    if len(cands) == 1:
        owners[work_id], rule[work_id] = set(cands), 'only_entry'
        continue
    tf = title_forms(title.get(work_id))
    named = {c for c in cands if name_forms(c) & tf}
    if named:
        owners[work_id], rule[work_id] = named, 'name_in_title'
        continue
    mains = cands & main_of.get(work_id, set())
    if len(cands) <= 3 and 1 <= len(mains) <= 2:
        owners[work_id], rule[work_id] = mains, 'main_publication'
    else:
        owners[work_id], rule[work_id] = set(), 'shared_platform_dropped'

owned = defaultdict(set)
for work_id, os_ in owners.items():
    for o in os_:
        owned[o].add(work_id)

# --- Use per entry (methods sections of research articles) ---
paper_use = con.sql(f"""
    SELECT DISTINCT m.pmcid, c.work_id
    FROM '{EXPORT}/mentions.parquet' m
    JOIN '{KEPT}' k USING (pmcid)
    JOIN '{EXPORT}/citations.parquet' c USING (pmcid, ref_id)
    WHERE m.sec_class = 'methods' AND m.in_table IS NOT TRUE AND m.in_caption IS NOT TRUE
      AND c.work_id IS NOT NULL
""").df()
users = paper_use.groupby('work_id').pmcid.apply(set).to_dict()


def n_users(works_: set) -> int:
    return len(set().union(*(users.get(w, set()) for w in works_))) if works_ else 0


entry_n = {e: n_users(ws) for e, ws in owned.items()}

# --- Merge identical sets, then subsets into the most-used superset ---
by_set = defaultdict(list)
for e, ws in owned.items():
    by_set[frozenset(ws)].append(e)
groups = {}  # representative -> (works, members)
for ws, members in by_set.items():
    rep = sorted(members, key=lambda e: (-entry_n[e], len(e), e))[0]
    groups[rep] = (set(ws), set(members))

work_groups = defaultdict(set)
for rep, (ws, _) in groups.items():
    for w in ws:
        work_groups[w].add(rep)
for rep in sorted(groups, key=lambda r: entry_n[r]):  # smallest first
    ws, members = groups[rep]
    supers = set.intersection(*(work_groups[w] for w in ws)) - {rep} if ws else set()
    supers = {s for s in supers if s in groups and groups[s][0] > ws}
    if supers:
        target = max(supers, key=lambda s: (entry_n[s], -len(s)))
        groups[target][1].update(members)
        for w in ws:
            work_groups[w].discard(rep)
        del groups[rep]

# Entries whose users are mostly another entry's (e.g. a duplicate record of the same paper).
group_users = {r: set().union(*(users.get(w, set()) for w in ws)) if ws else set() for r, (ws, _) in groups.items()}
for rep_ in sorted(groups, key=lambda r: len(group_users[r])):
    ws, members = groups[rep_]
    mine = group_users[rep_]
    if not mine:
        continue
    near = {g for w in ws for g in work_groups[w]} - {rep_}
    best = max((g for g in near if g in groups and len(group_users[g]) > len(mine)
                and len(mine & group_users[g]) >= 0.8 * len(mine)),
               key=lambda g: len(group_users[g]), default=None)
    if best:
        groups[best][0].update(ws)
        groups[best][1].update(members)
        group_users[best] |= mine
        for w in ws:
            work_groups[w].discard(rep_)
            work_groups[w].add(best)
        del groups[rep_]

tool_rows, paper_rows = [], []
def display_name(members: set) -> str:
    """The plainest entry: fewest wrapper affixes, then most used, then shortest."""
    def affix_count(e):
        return sum(t in AFFIXES for t in re.split(r'[^a-z0-9]+', e.lower()))
    return min(members, key=lambda e: (affix_count(e), -entry_n.get(e, 0), len(e), e))


for rep, (ws, members) in groups.items():
    rep = display_name(members)
    tool_rows.append({'tool': rep, 'entries': sorted(members), 'n_articles': n_users(ws)})
    paper_rows += [{'tool': rep, 'work_id': w} for w in ws]
tools = pd.DataFrame(tool_rows).sort_values(['n_articles', 'tool'], ascending=[False, True])
tools['rank'] = range(1, len(tools) + 1)
tools = tools[tools['rank'] <= TOP_N]
top = set(tools.tool)
tool_papers = pd.DataFrame(paper_rows)
tool_papers = tool_papers[tool_papers.tool.isin(top)].merge(works, on='work_id')

use = paper_use.merge(tool_papers[['tool', 'work_id']], on='work_id')[['pmcid', 'tool']].drop_duplicates()
con.register('use_df', use)
con.register('tools_df', tools)
con.register('papers_df', tool_papers)
con.register('owner_df', pd.DataFrame([{'work_id': w, 'owners': sorted(o), 'rule': rule[w],
                                        'listed_by': sorted(listed[w])} for w, o in owners.items()]))

con.sql(f"COPY (SELECT * FROM tools_df ORDER BY rank) TO '{OUT}/tools.parquet' (FORMAT parquet)")
con.sql(f"""
    COPY (
        SELECT u.pmcid, u.tool, TRY_CAST(a.pub_year AS INTEGER) AS year, a.topics
        FROM use_df u
        LEFT JOIN read_csv('{MANIFEST}', delim='\t', header=true, all_varchar=true) a USING (pmcid)
        ORDER BY u.tool, year
    ) TO '{OUT}/tool_use.parquet' (FORMAT parquet)
""")
con.sql(f"COPY (SELECT * FROM papers_df ORDER BY tool, n_citations DESC) TO '{OUT}/tool_papers.parquet' (FORMAT parquet)")
con.sql(f"COPY (SELECT * FROM owner_df ORDER BY work_id) TO '{OUT}/paper_owner.parquet' (FORMAT parquet)")

print('How each matched paper got its owner:')
con.sql("SELECT rule, count(*) AS papers FROM owner_df GROUP BY 1 ORDER BY 2 DESC").show()
print(f'Top {TOP_N} tools by research articles citing them in methods:')
con.sql("SELECT rank, tool, n_articles, len(entries) AS n_entries FROM tools_df ORDER BY rank LIMIT 20").show()
con.sql(f"""
    SELECT (SELECT count(*) FROM tools_df) AS tools, (SELECT min(n_articles) FROM tools_df) AS min_articles,
           (SELECT count(*) FROM use_df WHERE tool IN (SELECT tool FROM tools_df)) AS uses,
           (SELECT count(DISTINCT pmcid) FROM use_df WHERE tool IN (SELECT tool FROM tools_df)) AS articles,
           (SELECT count(*) FROM papers_df) AS tool_papers
""").show()
print(f'Written {OUT}/tools.parquet, tool_use.parquet, tool_papers.parquet, paper_owner.parquet')
