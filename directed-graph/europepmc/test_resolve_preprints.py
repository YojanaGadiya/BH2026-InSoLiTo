"""Tests for resolve_preprints.py. Run: pytest test_resolve_preprints.py"""

import pytest

from resolve_preprints import (
    is_preprint_doi,
    lookup_doi,
    parse_biorxiv,
    parse_crossref,
    server_from_prefix,
)


@pytest.mark.parametrize("doi", [
    "10.1101/765453",                    # bioRxiv, old numbering (GAMer)
    "10.1101/003905v2",                  # bioRxiv with version suffix
    "10.1101/2021.10.04.463034",         # bioRxiv, date numbering (AlphaFold-Multimer)
    "10.64898/2026.05.12.26353027",      # openRxiv prefix
    "10.21203/rs.3.rs-123456/v1",        # Research Square
    "10.20944/preprints202001.0001.v1",  # Preprints.org
    "10.26434/chemrxiv-2022-abc12",      # ChemRxiv
    "10.2139/ssrn.3456789",              # SSRN
    "10.48550/arxiv.1802.03426",         # arXiv (UMAP)
])
def test_preprint_dois(doi):
    assert is_preprint_doi(doi)


@pytest.mark.parametrize("doi", [
    "10.1101/gr.1239303",                # Genome Research shares the 10.1101 prefix
    "10.1101/gad.1234.5",                # Genes & Development
    "10.1101/pdb.prot5439",              # Cold Spring Harbor Protocols
    "10.1016/j.bpj.2019.11.3400",        # GAMer's journal paper
    "10.1002/pmic.201900012",            # Wiley journal, not ESS Open Archive
    None, "",
])
def test_not_preprint_dois(doi):
    assert not is_preprint_doi(doi)


@pytest.mark.parametrize("doi,expected", [
    ("10.1101/003905v2", "10.1101/003905"),
    ("10.1101/2021.10.04.463034v1", "10.1101/2021.10.04.463034"),
    # Junk picked up from reference text (real examples from the corpus)
    ("10.1101/039230v1.abstract", "10.1101/039230"),
    ("10.1101/060012v060013.full.pdf", "10.1101/060012"),
    ("10.1101/0187622015", "10.1101/018762"),                       # year glued on
    ("10.1101/072470:072470", "10.1101/072470"),
    ("10.1101/2021.05.26.445185v1%0ahttps:/www.biorxiv.org/content/10.1101/2021.05.26.445185v1.abstract",
     "10.1101/2021.05.26.445185"),
    ("10.1101/2021.08.10.21261747v1%0ahttps://www.medrxiv.org/x", "10.1101/2021.08.10.21261747"),  # medRxiv: 8 digits
    ("10.1101/303727.%202020", "10.1101/303727"),
    # medRxiv's first DOIs (2019-20) are 8 digits: keep all of them
    ("10.1101/19012229", "10.1101/19012229"),
    ("10.1101/19012229v1.full", "10.1101/19012229"),
    ("10.1101/20018572v2", "10.1101/20018572"),
    ("10.1101/765453", "10.1101/765453"),
    ("10.21203/rs.3.rs-123456/v1", "10.21203/rs.3.rs-123456/v1"),  # "/v1" is part of this DOI
    ("10.1101/gr.1239303", "10.1101/gr.1239303"),
])
def test_lookup_doi(doi, expected):
    assert lookup_doi(doi) == expected


def test_server_from_prefix():
    assert server_from_prefix("10.1101/765453") == "bioRxiv/medRxiv"
    assert server_from_prefix("10.1016/x") is None


def test_parse_crossref_linked():
    # Trimmed from the real Crossref record for 10.1101/765453 (GAMer)
    msg = {
        "type": "posted-content", "subtype": "preprint", "publisher": "openRxiv",
        "institution": [{"name": "bioRxiv"}],
        "relation": {"is-preprint-of": [
            {"id-type": "doi", "id": "10.1016/J.BPJ.2019.11.3400", "asserted-by": "subject"}]},
    }
    assert parse_crossref(msg) == {
        "status": "linked", "published_doi": "10.1016/j.bpj.2019.11.3400", "server": "bioRxiv"}


def test_parse_crossref_no_link():
    msg = {"type": "posted-content", "publisher": "Research Square Platform LLC", "relation": {}}
    assert parse_crossref(msg) == {
        "status": "no_link", "published_doi": None, "server": "Research Square Platform LLC"}


def test_parse_crossref_not_preprint():
    # Prefix looked like a preprint, but Crossref says it is a journal article
    msg = {"type": "journal-article", "publisher": "Cold Spring Harbor Laboratory"}
    assert parse_crossref(msg)["status"] == "not_preprint"


def test_parse_crossref_ignores_non_doi_relations():
    msg = {"type": "posted-content", "relation": {"is-preprint-of": [{"id-type": "uri", "id": "https://x"}]}}
    assert parse_crossref(msg)["status"] == "no_link"


def test_parse_biorxiv():
    # Trimmed from api.biorxiv.org/details/biorxiv/10.1101/001230 (libRoadRunner)
    linked = {"collection": [{"doi": "10.1101/001230", "version": "1", "server": "bioRxiv",
                              "published": "10.1093/bioinformatics/btv363"}]}
    assert parse_biorxiv(linked) == {
        "status": "linked", "published_doi": "10.1093/bioinformatics/btv363", "server": "bioRxiv"}
    # 10.1101/000745: two versions, never published
    unpublished = {"collection": [{"version": "1", "server": "bioRxiv", "published": "NA"},
                                  {"version": "2", "server": "bioRxiv", "published": "NA"}]}
    assert parse_biorxiv(unpublished)["status"] == "no_link"
    assert parse_biorxiv({"messages": [{"status": "no posts found"}], "collection": []}) is None


# --------------------------------------------------------------------------- #
# End to end, with Crossref and the bioRxiv API replaced by canned answers
# --------------------------------------------------------------------------- #

CR = "https://api.crossref.org/works/"
BX = "https://api.biorxiv.org/details/"
CANNED = {
    # GAMer: Crossref has the link
    CR + "10.1101/765453": (200, {"message": {
        "type": "posted-content", "institution": [{"name": "bioRxiv"}],
        "relation": {"is-preprint-of": [{"id-type": "doi", "id": "10.1016/j.bpj.2019.11.3400"}]}}}),
    # libRoadRunner: Crossref has no link, the bioRxiv API does
    CR + "10.1101/001230": (200, {"message": {"type": "posted-content", "relation": {}}}),
    BX + "biorxiv/10.1101/001230": (200, {"collection": [
        {"server": "bioRxiv", "published": "10.1093/bioinformatics/btv363"}]}),
    # Never published anywhere; cited both bare and as "v2"
    CR + "10.1101/003905": (200, {"message": {"type": "posted-content", "relation": {}}}),
    BX + "biorxiv/10.1101/003905": (200, {"collection": [{"server": "bioRxiv", "published": "NA"}]}),
    # Research Square: no bioRxiv fallback
    CR + "10.21203/rs.3.rs-1/v1": (200, {"message": {"type": "posted-content", "publisher": "Research Square"}}),
    # Unknown everywhere
    CR + "10.1101/2020.01.01.000001": (404, {}),
    BX + "biorxiv/10.1101/2020.01.01.000001": (200, {"collection": []}),
    BX + "medrxiv/10.1101/2020.01.01.000001": (200, {"collection": []}),
}


class FakeResponse:
    def __init__(self, status, body):
        self.status_code, self._body = status, body

    def json(self):
        return self._body


class FakeSession:
    def __init__(self):
        self.calls = []

    def get(self, url, timeout=None):
        self.calls.append(url)
        return FakeResponse(*CANNED[url])


def test_end_to_end(tmp_path, monkeypatch):
    import duckdb

    import download_europepmc as dl
    import resolve_preprints as rp

    works = tmp_path / "works.parquet"
    duckdb.sql(f"""COPY (SELECT * FROM (VALUES
        (1, '10.1016/j.bpj.2019.11.3400'),   -- GAMer, journal version
        (2, '10.1101/765453'),               -- GAMer, preprint
        (3, '10.1101/gr.1239303'),           -- Genome Research: not a preprint
        (4, '10.21203/rs.3.rs-1/v1'),        -- preprint, never published
        (5, '10.48550/arxiv.1802.03426'),    -- arXiv: skipped
        (6, '10.1101/001230'),               -- libRoadRunner preprint
        (7, '10.1093/bioinformatics/btv363'),-- libRoadRunner journal version
        (8, '10.1101/003905'),               -- unpublished preprint...
        (9, '10.1101/003905v2')              -- ...also cited with a version suffix
    ) t(work_id, doi)) TO '{works}' (FORMAT parquet)""")
    biotools = tmp_path / "biotools.csv"
    biotools.write_text("biotoolsID,pmid,pmcid,doi\nGAMer,,,10.1101/765453\nX,,,10.1101/2020.01.01.000001\n")
    session = FakeSession()
    monkeypatch.setattr(dl, "get_session", lambda ua: session)

    out, cache = tmp_path / "doi_aliases.parquet", tmp_path / "cache.jsonl"
    argv = ["--works", str(works), "--biotools", str(biotools), "--out", str(out), "--cache", str(cache),
            "--max-rps", "1000"]
    assert rp.main(argv) == 0
    assert CR + "10.1101/003905v2" not in session.calls  # looked up without the suffix, once

    cols = "preprint_doi, canonical_doi, status, source, listed_in, preprint_work_id, published_work_id"
    rows = {r[0]: r[1:] for r in duckdb.sql(f"SELECT {cols} FROM '{out}'").fetchall()}
    assert rows["10.1101/765453"] == (
        "10.1016/j.bpj.2019.11.3400", "linked", "crossref", "biotools,works", 2, 1)
    assert rows["10.1101/001230"] == (
        "10.1093/bioinformatics/btv363", "linked", "biorxiv_api", "works", 6, 7)
    assert rows["10.1101/003905"][:2] == ("10.1101/003905", "no_link")
    assert rows["10.1101/003905v2"][:2] == ("10.1101/003905", "no_link")  # merges with the bare DOI
    assert rows["10.21203/rs.3.rs-1/v1"][1] == "no_link"
    assert rows["10.48550/arxiv.1802.03426"][1] == "skipped_arxiv"
    assert rows["10.1101/2020.01.01.000001"][1:4] == ("not_found", "crossref", "biotools")
    assert "10.1101/gr.1239303" not in rows

    # Second run asks nothing: everything is cached
    def no_network(ua):
        raise AssertionError("second run should not call any API")

    monkeypatch.setattr(dl, "get_session", no_network)
    assert rp.main(argv) == 0


def test_biorxiv_error_is_retried(tmp_path, monkeypatch):
    """If the bioRxiv API fails, the DOI is recorded as an error, not as 'no link'."""
    import download_europepmc as dl
    import resolve_preprints as rp

    class Flaky(FakeSession):
        def get(self, url, timeout=None):
            if url.startswith(BX):
                return FakeResponse(503, {})
            return super().get(url, timeout)

    monkeypatch.setattr(dl, "get_session", lambda ua: Flaky())
    rec = rp.resolve_one("10.1101/001230", "ua", dl.RateLimiter(0))
    assert rec["status"] == "error: biorxiv 503"
