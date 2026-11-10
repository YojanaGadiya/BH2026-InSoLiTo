"""Report bio.tools entries whose publication is a preprint, and whether its journal version is listed.

One row per (tool, preprint DOI) listed in bio.tools, with the journal version if one is
known and what that means for the entry. Sorted so the most useful fixes come first.

Where the journal version comes from (``link_source``):
  crossref     Crossref's ``is-preprint-of`` relation (via europepmc/resolve_preprints.py)
  biorxiv_api  bioRxiv/medRxiv's own API, for links Crossref does not have
  title_match  no recorded link, but a journal paper with exactly the same title is cited
               in the corpus. Likely, but check before curating.

``finding``:
  missing_published_version         journal version exists, not listed in this entry
  published_version_on_other_entry  ...and it is listed under a different bio.tools entry
                                    (often a duplicate entry, e.g. "rctd" / "bioconductor-spacexr")
  both_listed_preprint_primary      both listed, but the preprint is typed Primary
  both_listed                       both listed: nothing to do
  no_published_version_found        no journal version known (may still exist)
  doi_not_registered                preprint DOI unknown to Crossref: probably a typo
  arxiv_not_checked                 arXiv DOIs are not in Crossref; not checked

Needs doi_aliases.parquet from europepmc/resolve_preprints.py.

Run from BH2026-InSoLiTo/:
    python directed-graph/biotools/report_preprint_links.py
"""

from pathlib import Path

import duckdb

BIOTOOLS = 'directed-graph/biotools/data/biotools_publications.csv'
MAIN_UNIQUE = 'directed-graph/biotools/data/biotools_main_publications_unique.csv'
MAIN_SHARED = 'directed-graph/biotools/data/biotools_main_publications_shared.csv'
ALIASES = 'directed-graph/europepmc/data/export/doi_aliases.parquet'
WORKS = 'directed-graph/europepmc/data/export/works.parquet'
OUTPUT = 'directed-graph/biotools/data/biotools_preprint_links.csv'

con = duckdb.connect()

con.sql(f"""
    CREATE TEMP TABLE bt AS
    SELECT biotoolsID, lower(trim(doi)) AS doi, CAST(pmid AS VARCHAR) AS pmid, type, is_primary
    FROM '{BIOTOOLS}' WHERE doi IS NOT NULL
""")
con.sql(f"""
    CREATE TEMP TABLE main_pub AS
    SELECT biotoolsID, lower(doi) AS doi FROM '{MAIN_UNIQUE}'
    UNION
    SELECT unnest(string_split(biotoolsIDs, ';')), lower(doi) FROM '{MAIN_SHARED}'
""")
con.sql(f"""
    CREATE TEMP TABLE works AS
    SELECT work_id, doi, CAST(pmid AS VARCHAR) AS pmid, n_citations, year, title,
           regexp_replace(lower(title), '[^a-z0-9]', '', 'g') AS norm_title
    FROM '{WORKS}'
""")

# Journal paper with the same title as an unlinked preprint (both cited in the corpus).
# Titles under 25 characters are too generic to trust. The most cited candidate wins.
con.sql(f"""
    CREATE TEMP TABLE title_link AS
    SELECT a.preprint_doi, arg_max(j.doi, j.n_citations) AS journal_doi
    FROM '{ALIASES}' a
    JOIN works p ON p.work_id = a.preprint_work_id
    JOIN works j ON j.norm_title = p.norm_title AND j.work_id <> p.work_id
    WHERE a.status = 'no_link' AND length(p.title) >= 25
      AND j.doi IS NOT NULL AND j.doi NOT IN (SELECT preprint_doi FROM '{ALIASES}')
    GROUP BY 1
""")

con.sql(f"""
    CREATE TEMP TABLE report AS
    WITH pre AS (
        SELECT b.biotoolsID, b.doi AS preprint_doi, b.type AS preprint_type, b.is_primary AS preprint_is_primary,
               a.server AS preprint_server, a.status, a.preprint_work_id,
               coalesce(a.published_doi, t.journal_doi) AS published_doi,
               CASE WHEN a.published_doi IS NOT NULL THEN a.source
                    WHEN t.journal_doi IS NOT NULL THEN 'title_match' END AS link_source
        FROM bt b
        JOIN '{ALIASES}' a ON a.preprint_doi = b.doi
        LEFT JOIN title_link t ON t.preprint_doi = a.preprint_doi
    ), enriched AS (
        SELECT p.*,
               EXISTS (SELECT 1 FROM main_pub m WHERE m.biotoolsID = p.biotoolsID AND m.doi = p.preprint_doi)
                   AS preprint_is_main_publication,
               wj.title AS published_title, wj.year AS published_year,
               wj.n_citations AS corpus_citations_published, wp.n_citations AS corpus_citations_preprint,
               (SELECT any_value(b2.is_primary) FROM bt b2 WHERE b2.biotoolsID = p.biotoolsID
                    AND (b2.doi = p.published_doi OR (wj.pmid IS NOT NULL AND b2.pmid = wj.pmid)))
                   AS published_is_primary,
               (SELECT string_agg(DISTINCT b3.biotoolsID, ';') FROM bt b3
                    WHERE b3.biotoolsID <> p.biotoolsID AND b3.doi = p.published_doi)
                   AS published_listed_under
        FROM pre p
        LEFT JOIN works wj ON wj.doi = p.published_doi
        LEFT JOIN works wp ON wp.work_id = p.preprint_work_id
    )
    SELECT biotoolsID, preprint_doi, preprint_server, published_doi, link_source,
           CASE
               WHEN status = 'skipped_arxiv' THEN 'arxiv_not_checked'
               WHEN published_doi IS NULL AND status = 'not_found' THEN 'doi_not_registered'
               WHEN published_doi IS NULL THEN 'no_published_version_found'
               WHEN published_is_primary IS NOT NULL AND preprint_is_primary AND NOT published_is_primary
                   THEN 'both_listed_preprint_primary'
               WHEN published_is_primary IS NOT NULL THEN 'both_listed'
               WHEN published_listed_under IS NOT NULL THEN 'published_version_on_other_entry'
               ELSE 'missing_published_version'
           END AS finding,
           preprint_is_main_publication, preprint_type, published_listed_under,
           coalesce(corpus_citations_published, 0) AS corpus_citations_published,
           coalesce(corpus_citations_preprint, 0) AS corpus_citations_preprint,
           published_year, published_title
    FROM enriched
""")

print('bio.tools entries citing a preprint:')
con.sql("""
    SELECT count(*) AS rows, count(DISTINCT biotoolsID) AS tools,
           count(DISTINCT biotoolsID) FILTER (WHERE preprint_is_main_publication) AS preprint_is_main_pub
    FROM report
""").show()

print('Findings (citations = citations in the corpus to the journal version):')
con.sql("""
    SELECT finding, count(*) AS rows, count(DISTINCT biotoolsID) AS tools,
           sum(corpus_citations_published) AS citations_published,
           sum(corpus_citations_preprint) AS citations_preprint
    FROM report GROUP BY 1 ORDER BY rows DESC
""").show()

print('Journal version found by:')
con.sql("""
    SELECT link_source, count(*) AS rows FROM report WHERE link_source IS NOT NULL
    GROUP BY 1 ORDER BY 2 DESC
""").show()

print('Most cited missing journal versions:')
con.sql("""
    SELECT biotoolsID, preprint_doi, published_doi, link_source, corpus_citations_published AS cites,
           left(published_title, 45) AS title
    FROM report WHERE finding IN ('missing_published_version', 'published_version_on_other_entry')
    ORDER BY cites DESC LIMIT 15
""").show(max_width=250)

Path(OUTPUT).parent.mkdir(parents=True, exist_ok=True)
con.sql(f"""
    COPY (
        SELECT * FROM report
        ORDER BY CASE finding
                     WHEN 'missing_published_version' THEN 1
                     WHEN 'published_version_on_other_entry' THEN 2
                     WHEN 'both_listed_preprint_primary' THEN 3
                     WHEN 'doi_not_registered' THEN 4
                     WHEN 'no_published_version_found' THEN 5
                     WHEN 'arxiv_not_checked' THEN 6
                     ELSE 7 END,
                 corpus_citations_published DESC, biotoolsID
    ) TO '{OUTPUT}' (HEADER)
""")
print(f'Written {OUTPUT}')
