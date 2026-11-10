"""Step 1: find which bio.tools papers appear among the works cited in the corpus.

Matches each bio.tools publication (Vedran's table) against Yojana's `works`
by DOI, PMID or PMCID. Prints counts and examples and writes the matches to
tool_works.parquet (one row per tool and paper), read by the next steps.

Run from BH2026-InSoLiTo/:
    .venv/bin/python directed-graph/citation-order/01_match_tools.py
"""

from pathlib import Path

import duckdb

BIOTOOLS = 'directed-graph/biotools/data/biotools_publications.csv'
WORKS = 'directed-graph/europepmc/data/export/works.parquet'
# Optional, from europepmc/resolve_preprints.py: preprint DOI -> published DOI.
ALIASES = 'directed-graph/europepmc/data/export/doi_aliases.parquet'
OUTPUT = 'directed-graph/citation-order/data/tool_works.parquet'

con = duckdb.connect()

# DOIs are compared after mapping each preprint to its published version, on
# both sides, so a tool listed with its bioRxiv DOI also matches papers that
# cite the journal version, and the reverse. Without the aliases file, every
# DOI maps to itself.
if Path(ALIASES).exists():
    # canonical_doi is the published version, or for an unpublished bioRxiv
    # preprint its DOI without a "v2" suffix.
    con.sql(f"""
        CREATE TEMP TABLE alias AS
        SELECT preprint_doi AS doi, canonical_doi AS canonical
        FROM '{ALIASES}' WHERE canonical_doi <> preprint_doi
    """)
    print(f'Using preprint aliases from {ALIASES}')
else:
    con.sql("CREATE TEMP TABLE alias (doi VARCHAR, canonical VARCHAR)")
    print(f'No {ALIASES}: matching DOIs as they are')

# One join per identifier: a single join with OR conditions forces DuckDB to
# compare every row with every row. UNION also removes the duplicates left
# when a paper matches by more than one identifier.
con.sql(f"""
    CREATE TEMP TABLE tool_works AS
    WITH v_doi AS (
        SELECT v.biotoolsID, coalesce(a.canonical, lower(v.doi)) AS doi
        FROM '{BIOTOOLS}' v LEFT JOIN alias a ON a.doi = lower(v.doi)
        WHERE v.doi IS NOT NULL
    ), w_doi AS (
        SELECT w.work_id, w.title, w.n_citations, coalesce(a.canonical, w.doi) AS doi
        FROM '{WORKS}' w LEFT JOIN alias a ON a.doi = w.doi
        WHERE w.doi IS NOT NULL
    )
    SELECT v.biotoolsID, w.work_id, w.title, w.n_citations
    FROM v_doi v
    JOIN w_doi w ON v.doi = w.doi
    UNION
    SELECT v.biotoolsID, w.work_id, w.title, w.n_citations
    FROM '{BIOTOOLS}' v
    JOIN '{WORKS}' w ON CAST(v.pmid AS VARCHAR) = CAST(w.pmid AS VARCHAR)
    UNION
    SELECT v.biotoolsID, w.work_id, w.title, w.n_citations
    FROM '{BIOTOOLS}' v
    JOIN '{WORKS}' w ON CAST(v.pmcid AS VARCHAR) = CAST(w.pmcid AS VARCHAR)
""")

print('Tools and papers found in the corpus:')
con.sql("""
    SELECT count(DISTINCT biotoolsID) AS tools, count(DISTINCT work_id) AS papers
    FROM tool_works
""").show()

print('Most cited tool papers:')
con.sql("""
    SELECT biotoolsID, work_id, n_citations, left(title, 60) AS title
    FROM tool_works
    ORDER BY n_citations DESC
    LIMIT 10
""").show()

Path(OUTPUT).parent.mkdir(parents=True, exist_ok=True)  # data/ is git-ignored, so may not exist
con.sql(f"COPY tool_works TO '{OUTPUT}' (FORMAT parquet)")
print(f'Written {OUTPUT}')
