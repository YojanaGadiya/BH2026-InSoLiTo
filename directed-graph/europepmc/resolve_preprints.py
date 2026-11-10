"""Find preprint DOIs and look up their published versions, so both match the same tool.

A tool can be cited through its preprint (bioRxiv, Research Square, ...) or through
the journal paper that came later, and bio.tools may list either one. GAMer, for
example, lists only its bioRxiv preprint (10.1101/765453), while some papers cite the
Biophysical Journal version (10.1016/j.bpj.2019.11.3400). Matched on DOI alone, those
citations are lost.

This script:

1. collects preprint DOIs from ``works.parquet`` and the bio.tools publication list,
   recognised by DOI prefix (see ``PREPRINT_DOI``). bioRxiv/medRxiv DOIs are looked up
   without a version suffix ("10.1101/003905v2" -> "10.1101/003905"): the suffix is not
   part of the registered DOI,
2. asks Crossref for each one. Preprint servers record the published version under
   ``relation.is-preprint-of``. Crossref misses many bioRxiv/medRxiv links, so when it
   has none, bioRxiv's own API (api.biorxiv.org) is asked too. arXiv DOIs are
   registered with DataCite, not Crossref, so they are skipped,
3. writes ``doi_aliases.parquet``: one row per preprint DOI as it appears in the data,
   with its ``canonical_doi`` (the published version if there is one, else the DOI
   without its version suffix) and the ``work_id`` of both versions in ``works``.

Answers are cached in ``<cache>``, so re-runs only ask about new DOIs and earlier errors.

Matching then maps every DOI to its ``canonical_doi`` on both sides
(see ``citation-order/01_match_tools.py``).

Example::

    python resolve_preprints.py --limit 50            # quick test
    python resolve_preprints.py --contact you@example.org
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.parse import quote

import duckdb
import pyarrow as pa
import requests

import download_europepmc as dl  # HTTP session with retries, rate limiter, TLS setup
from parse_jats import norm_doi

logger = logging.getLogger("resolve_preprints")

CROSSREF = "https://api.crossref.org/works/"
BIORXIV_API = "https://api.biorxiv.org/details/"  # + {biorxiv|medrxiv}/{doi}
CACHE_VERSION = 2  # bump when the lookup logic changes, so old cache entries are redone

# DOI prefixes of preprint servers. 10.1101 is shared with Cold Spring Harbor journals
# (Genome Research is 10.1101/gr.…), so only bioRxiv/medRxiv's own patterns count:
# plain digits (old) or a date (since late 2019). 10.64898 is openRxiv's newer prefix.
PREPRINT_DOI = re.compile(
    r"^(?:"
    r"10\.1101/(?:\d{6,}|\d{4}\.\d{2}\.\d{2}\.\d+)"   # bioRxiv / medRxiv
    r"|10\.64898/"                                     # bioRxiv / medRxiv (openRxiv)
    r"|10\.21203/rs\."                                 # Research Square
    r"|10\.20944/preprints"                            # Preprints.org
    r"|10\.26434/chemrxiv"                             # ChemRxiv
    r"|10\.2139/ssrn\."                                # SSRN
    r"|10\.31219/osf\.io"                              # OSF Preprints
    r"|10\.22541/au\."                                 # Authorea
    r"|10\.1002/essoar\."                              # ESS Open Archive
    r"|10\.7287/peerj\.preprints"                      # PeerJ Preprints
    r"|10\.48550/arxiv\."                              # arXiv (DataCite, not Crossref)
    r")"
)
RXIV_DOI = re.compile(r"^10\.(?:1101|64898)/(?:\d{6,}|\d{4}\.\d{2}\.\d{2}\.\d+)")
# The registered part of a bioRxiv/medRxiv DOI has a fixed shape: a date plus 6
# (bioRxiv) or 8 (medRxiv) digits since late 2019; before that 6 digits for bioRxiv
# and 8 for medRxiv's first DOIs (2019-20, "19012229" = year, month, number). Anything
# after it is a version ("v2") or junk picked up from the reference text (".abstract",
# ".full.pdf", a glued-on year, "%0ahttps://..."), so it is cut off for the lookup.
RXIV_CORE = re.compile(
    r"^10\.(?:1101|64898)/(?:\d{4}\.\d{2}\.\d{2}\.\d{6,8}|(?:19|20)\d{6}(?=\D|$)|\d{6})"
)
NOT_IN_CROSSREF = ("10.48550/",)

SERVER_BY_PREFIX = {
    "10.1101/": "bioRxiv/medRxiv", "10.64898/": "bioRxiv/medRxiv", "10.21203/": "Research Square",
    "10.20944/": "Preprints.org", "10.26434/": "ChemRxiv", "10.2139/": "SSRN",
    "10.31219/": "OSF Preprints", "10.22541/": "Authorea", "10.1002/": "ESS Open Archive",
    "10.7287/": "PeerJ Preprints", "10.48550/": "arXiv",
}


def is_preprint_doi(doi: str | None) -> bool:
    return bool(doi and PREPRINT_DOI.match(doi))


def is_rxiv_doi(doi: str) -> bool:
    return bool(RXIV_DOI.match(doi))


def lookup_doi(doi: str) -> str:
    """The registered DOI: for bioRxiv/medRxiv, only the fixed-shape ID (see RXIV_CORE)."""
    m = RXIV_CORE.match(doi)
    return m.group(0) if m else doi


def server_from_prefix(doi: str) -> str | None:
    return next((s for p, s in SERVER_BY_PREFIX.items() if doi.startswith(p)), None)


def parse_crossref(msg: dict) -> dict:
    """Interpret one Crossref ``message``: is it a preprint, and what was it published as?

    status:
      linked        preprint with a published version (``published_doi``)
      no_link       preprint, no published version recorded (yet)
      not_preprint  Crossref says it is not posted content: the prefix pattern was wrong
    """
    institution = (msg.get("institution") or [{}])[0].get("name")
    server = institution or msg.get("publisher")
    if msg.get("type") != "posted-content":
        return {"status": "not_preprint", "published_doi": None, "server": server}
    published = None
    for rel in (msg.get("relation") or {}).get("is-preprint-of", []):
        if (rel.get("id-type") or "").lower() == "doi" and rel.get("id"):
            published = norm_doi(rel["id"])
            break
    return {"status": "linked" if published else "no_link", "published_doi": published, "server": server}


def parse_biorxiv(body: dict) -> dict | None:
    """Interpret api.biorxiv.org/details: None if the DOI is not on that server.

    The answer lists every version of the preprint; any of them may carry the
    published DOI ("NA" when there is none).
    """
    versions = body.get("collection") or []
    if not versions:
        return None
    server = versions[-1].get("server")
    for v in reversed(versions):
        published = norm_doi(v.get("published")) if v.get("published") not in (None, "", "NA") else None
        if published:
            return {"status": "linked", "published_doi": published, "server": server}
    return {"status": "no_link", "published_doi": None, "server": server}


def _get_json(url: str, user_agent: str, limiter: dl.RateLimiter) -> tuple[int | str, dict | None]:
    """(HTTP status or error name, parsed JSON or None)."""
    session = dl.get_session(user_agent)
    limiter.wait()
    try:
        resp = session.get(url, timeout=60)
    except requests.RequestException as exc:
        return type(exc).__name__, None
    if resp.status_code != 200:
        return resp.status_code, None
    try:
        return 200, resp.json()
    except ValueError:
        return "bad JSON", None


def resolve_one(doi: str, user_agent: str, limiter: dl.RateLimiter) -> dict:
    """Look up one (version-free) preprint DOI. Crossref first, bioRxiv API as fallback."""
    rec = {"doi": doi, "v": CACHE_VERSION, "source": "crossref"}
    if doi.startswith(NOT_IN_CROSSREF):
        return {**rec, "status": "skipped_arxiv", "published_doi": None, "server": "arXiv", "source": None}

    code, body = _get_json(CROSSREF + quote(doi, safe="/"), user_agent, limiter)
    if code == 200 and body and "message" in body:
        rec.update(parse_crossref(body["message"]))
    elif code == 404:
        rec.update(status="not_found", published_doi=None, server=server_from_prefix(doi))
    else:
        return {**rec, "status": f"error: crossref {code}", "published_doi": None, "server": None}

    if is_rxiv_doi(doi) and rec["status"] in ("no_link", "not_found"):
        for server in ("biorxiv", "medrxiv"):
            code, body = _get_json(f"{BIORXIV_API}{server}/{doi}", user_agent, limiter)
            if code != 200:  # the API answers 200 with an empty list for unknown DOIs
                return {**rec, "status": f"error: biorxiv {code}"}  # retried next run
            found = parse_biorxiv(body or {})
            if found is None:
                continue  # not on this server: try medRxiv
            if found["status"] == "linked" or rec["status"] == "not_found":
                rec.update(found, source="biorxiv_api")
            break
    return rec


# --------------------------------------------------------------------------- #
# Inputs
# --------------------------------------------------------------------------- #


def preprints_in_works(works: Path) -> set[str]:
    rows = duckdb.sql(f"""
        SELECT DISTINCT doi FROM read_parquet('{works}')
        WHERE regexp_matches(doi, '{PREPRINT_DOI.pattern}')
    """).fetchall()
    return {d for (d,) in rows}


def preprints_in_biotools(path: Path) -> set[str]:
    if not path.exists():
        logger.warning("bio.tools publications not found at %s, skipping them", path)
        return set()
    with path.open(newline="") as fh:
        dois = {norm_doi(row.get("doi")) for row in csv.DictReader(fh)}
    return {d for d in dois if is_preprint_doi(d)}


def load_cache(cache: Path) -> dict[str, dict]:
    """Cached answers by lookup DOI. Entries from older versions of the logic are dropped."""
    done: dict[str, dict] = {}
    if cache.exists():
        with cache.open() as fh:
            for line in fh:
                rec = json.loads(line)
                if rec.get("v") == CACHE_VERSION:
                    done[rec["doi"]] = rec  # later lines win
    return done


# --------------------------------------------------------------------------- #
# Output
# --------------------------------------------------------------------------- #

ALIAS_SCHEMA = pa.schema([
    ("preprint_doi", pa.string()),    # as it appears in works / bio.tools
    ("lookup_doi", pa.string()),      # without bioRxiv/medRxiv version suffix
    ("published_doi", pa.string()),   # journal version, if any
    ("canonical_doi", pa.string()),   # published_doi, else lookup_doi: match on this
    ("server", pa.string()),
    ("status", pa.string()),
    ("source", pa.string()),          # crossref | biorxiv_api
    ("listed_in", pa.string()),       # works, biotools
])


def alias_rows(sources: dict[str, set[str]], answers: dict[str, dict]) -> list[dict]:
    rows = []
    for doi in sorted(sources):
        key = lookup_doi(doi)
        ans = answers.get(key)
        if ans is None:
            continue  # not looked up yet (e.g. --limit)
        rows.append({
            "preprint_doi": doi, "lookup_doi": key, "published_doi": ans["published_doi"],
            "canonical_doi": ans["published_doi"] or key, "server": ans["server"],
            "status": ans["status"], "source": ans.get("source"),
            "listed_in": ",".join(sorted(sources[doi])),
        })
    return rows


def write_aliases(rows: list[dict], works: Path, out: Path) -> None:
    con = duckdb.connect()
    con.register("aliases", pa.Table.from_pylist(rows, schema=ALIAS_SCHEMA))
    tmp = out.with_suffix(".parquet.part")
    con.execute(f"""
        COPY (
            SELECT a.*, wp.work_id AS preprint_work_id, wj.work_id AS published_work_id
            FROM aliases a
            LEFT JOIN (SELECT doi, work_id FROM read_parquet('{works}')) wp ON wp.doi = a.preprint_doi
            LEFT JOIN (SELECT doi, work_id FROM read_parquet('{works}')) wj ON wj.doi = a.published_doi
            ORDER BY a.preprint_doi
        ) TO '{tmp}' (FORMAT parquet, COMPRESSION zstd)
    """)
    tmp.replace(out)

    logger.info("Status of %d preprint DOIs:", len(rows))
    for status, source, n in con.execute(f"""SELECT status, coalesce(source, '-'), count(*) FROM '{out}'
            GROUP BY 1, 2 ORDER BY 3 DESC""").fetchall():
        logger.info("  %-22s %-12s %7d", status, source, n)
    (n_both,) = con.execute(f"""SELECT count(*) FROM '{out}'
        WHERE preprint_work_id IS NOT NULL AND published_work_id IS NOT NULL""").fetchone()
    (n_bt,) = con.execute(f"""SELECT count(*) FROM '{out}'
        WHERE listed_in LIKE '%biotools%' AND status = 'linked'""").fetchone()
    (n_ver,) = con.execute(f"SELECT count(*) FROM '{out}' WHERE preprint_doi <> lookup_doi").fetchone()
    logger.info("Both versions cited in the corpus (works that now merge): %d", n_both)
    logger.info("Versioned DOIs (…v2) mapped to their bare DOI: %d", n_ver)
    logger.info("bio.tools preprint DOIs that now also match the published paper: %d", n_bt)
    logger.info("Written %s", out)


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #


def main(argv: list[str] | None = None) -> int:
    here = Path(__file__).parent
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--works", type=Path, default=here / "data" / "export" / "works.parquet")
    p.add_argument("--biotools", type=Path,
                   default=here.parent / "biotools" / "data" / "biotools_publications.csv")
    p.add_argument("--out", type=Path, default=None,
                   help="Output file (default: doi_aliases.parquet next to --works).")
    p.add_argument("--cache", type=Path, default=here / "data" / "preprints" / "crossref_cache.jsonl")
    p.add_argument("--contact", default="",
                   help="Email for Crossref's 'polite' pool (faster, higher limits). "
                        "Without it the public pool is used: 1 request at a time.")
    p.add_argument("--workers", type=int, default=None, help="Default: 3 with --contact, else 1.")
    p.add_argument("--max-rps", type=float, default=None, help="Default: 8 with --contact, else 4.")
    p.add_argument("--limit", type=int, default=None, help="Only look up this many new DOIs (testing).")
    p.add_argument("--ca-bundle", default=None, help="PEM CA bundle, e.g. behind a TLS-inspecting proxy.")
    args = p.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    out = args.out or args.works.with_name("doi_aliases.parquet")
    workers = args.workers or (3 if args.contact else 1)
    max_rps = args.max_rps or (8.0 if args.contact else 4.0)
    dl._CA_BUNDLE = args.ca_bundle
    user_agent = "BH2026-InSoLiTo preprint resolver (https://github.com/inab/BH2026-InSoLiTo)"
    if args.contact:
        user_agent += f" mailto:{args.contact}"

    sources: dict[str, set[str]] = {}
    for doi in preprints_in_works(args.works):
        sources.setdefault(doi, set()).add("works")
    for doi in preprints_in_biotools(args.biotools):
        sources.setdefault(doi, set()).add("biotools")
    keys = {lookup_doi(d) for d in sources}
    logger.info("%d preprint DOIs (%d in works, %d in bio.tools), %d after removing version suffixes",
                len(sources), sum("works" in s for s in sources.values()),
                sum("biotools" in s for s in sources.values()), len(keys))

    args.cache.parent.mkdir(parents=True, exist_ok=True)
    cached = load_cache(args.cache)
    pending = sorted(k for k in keys if k not in cached or cached[k]["status"].startswith("error"))
    todo = pending[: args.limit] if args.limit else pending
    logger.info("%d already cached, %d to look up%s (%d workers, max %.0f/s)",
                len(keys) - len(pending), len(todo),
                f" (--limit; {len(pending) - len(todo)} left for later)" if len(todo) < len(pending) else "",
                workers, max_rps)

    limiter = dl.RateLimiter(max_rps)
    t0 = time.time()
    with args.cache.open("a") as cache_fh, ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(resolve_one, d, user_agent, limiter) for d in todo]
        for i, fut in enumerate(as_completed(futures), 1):
            rec = fut.result()
            cached[rec["doi"]] = rec
            cache_fh.write(json.dumps(rec) + "\n")  # written by this thread only
            if i % 500 == 0 or i == len(todo):
                cache_fh.flush()
                rate = i / max(time.time() - t0, 1e-6)
                logger.info("%d/%d (%.1f/s, ETA %.0f min)", i, len(todo), rate, (len(todo) - i) / rate / 60)

    write_aliases(alias_rows(sources, cached), args.works, out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
