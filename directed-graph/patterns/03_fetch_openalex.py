"""Step 3: fetch authors, institutions and yearly counts from OpenAlex.

Three steps, each cached and resumable (re-run the same command after a stop):

  tools     the top tools' papers (from step 1): OpenAlex id, year, citations per year,
            and every author with their OpenAlex author id and institutions
  articles  every research article that uses a top tool: its authors and institutions,
            to tell self-citation and same-organization use from cross-lab use
  authors   the first, last and corresponding authors of the tool papers: papers per
            year over their whole career, to see developer output around a release

OpenAlex ids are disambiguated (one id per person or institution), so matching does not
rely on names. About 1,400 requests in total. OpenAlex reports a small cost per request
(about $0.0001) against a free daily allowance; pass --api-key if you have one.

Writes to directed-graph/patterns/data/:
  oa_tool_papers.parquet   work_id, openalex_id, year, cited_by_count, counts_by_year, authors
  oa_articles.parquet      pmcid, openalex_id, year, author_ids, institution_ids, countries
  oa_authors.parquet       author_id, name, works_count, cited_by_count, h_index
  oa_author_years.parquet  author_id, year, n_works
(raw answers are cached in patterns/data/openalex/*.jsonl)

Run from BH2026-InSoLiTo/:
    python directed-graph/patterns/03_fetch_openalex.py --contact you@example.org
    python directed-graph/patterns/03_fetch_openalex.py --step tools --limit 20   # quick test
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from concurrent.futures import TimeoutError as FutureTimeout
from pathlib import Path
from urllib.parse import urlencode

import threading

import duckdb
import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'europepmc'))
import download_europepmc as dl  # noqa: E402  session with retries, rate limiter, TLS setup
from parse_jats import norm_doi  # noqa: E402

logger = logging.getLogger('fetch_openalex')

API = 'https://api.openalex.org'
DATA = Path('directed-graph/patterns/data')
CACHE = DATA / 'openalex'
BATCH = 100  # ids per OR-filter request (OpenAlex maximum)
REQUEST_DEADLINE = 120  # seconds, including retries
PAPER_FIELDS = 'id,doi,ids,title,publication_year,cited_by_count,counts_by_year,authorships'
ARTICLE_FIELDS = 'id,ids,publication_year,authorships'


_local = threading.local()


def openalex_session(user_agent: str) -> requests.Session:
    """Own retry policy: short capped backoff and Retry-After ignored, so a rate limit
    shows up as an error in the log within about a minute instead of a silent wait."""
    if not hasattr(_local, 'session'):
        retry = Retry(total=4, backoff_factor=3, status_forcelist=(429, 500, 502, 503, 504),
                      allowed_methods=('GET',), respect_retry_after_header=False, raise_on_status=False)
        s = requests.Session()
        s.mount('https://', HTTPAdapter(max_retries=retry, pool_maxsize=4))
        s.headers['User-Agent'] = user_agent
        if dl._CA_BUNDLE:
            s.verify = dl._CA_BUNDLE
        _local.session = s
    return _local.session


class Client:
    def __init__(self, contact: str, api_key: str, max_rps: float):
        self.params = {k: v for k, v in (('mailto', contact), ('api_key', api_key)) if v}
        self.limiter = dl.RateLimiter(max_rps)
        self.ua = 'BH2026-InSoLiTo pattern analysis (https://github.com/inab/BH2026-InSoLiTo)'

    def get(self, path: str, **params) -> dict:
        url = f'{API}/{path}?{urlencode({**params, **self.params}, safe=":|")}'
        self.limiter.wait()
        t0 = time.time()
        label = f"{path} {str(params.get('filter', ''))[:40]}"
        logger.debug('-> %s', label)
        one = ThreadPoolExecutor(max_workers=1)  # hard deadline, whatever the socket does
        fut = one.submit(get_session(self.ua).get, url, timeout=(10, 60))
        try:
            resp = fut.result(timeout=REQUEST_DEADLINE)
        except FutureTimeout:
            raise RuntimeError(f'no answer from OpenAlex within {REQUEST_DEADLINE}s for {label}') from None
        finally:
            one.shutdown(wait=False)  # don't wait for a hung request
        logger.info('   %s %s in %.1fs', resp.status_code, label, time.time() - t0)
        if resp.status_code == 429:
            raise SystemExit(f'OpenAlex rate limit (HTTP 429) after retries: {resp.text[:300]}\n'
                             'Progress is cached: wait (or add --api-key) and run the same command again.')
        if resp.status_code != 200:
            raise RuntimeError(f'HTTP {resp.status_code} for {path}: {resp.text[:200]}')
        return resp.json()


get_session = openalex_session  # tests replace this


# --------------------------------------------------------------------------- helpers

def short_id(url: str | None) -> str | None:
    """'https://openalex.org/A123' -> 'A123'; PMCID/PMID URLs -> the bare id."""
    if not url:
        return None
    tail = url.rstrip('/').rsplit('/', 1)[-1]
    if 'ncbi.nlm.nih.gov/pmc' in url:
        return 'PMC' + tail.removeprefix('PMC')
    return tail


def compact_authorships(work: dict) -> list[dict]:
    out = []
    for a in work.get('authorships') or []:
        insts = a.get('institutions') or []
        out.append({
            'author_id': short_id((a.get('author') or {}).get('id')),
            'name': (a.get('author') or {}).get('display_name'),
            'position': a.get('author_position'),
            'is_corresponding': bool(a.get('is_corresponding')),
            'institution_ids': [short_id(i.get('id')) for i in insts if i.get('id')],
            'institution_roots': sorted({short_id(x) for i in insts for x in (i.get('lineage') or [])[:1]}),
            'countries': sorted({i.get('country_code') for i in insts if i.get('country_code')}),
        })
    return out


def read_cache(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open() as fh:
        return [json.loads(line) for line in fh if line.strip()]


def run_batches(client: Client, filter_key: str, values: list[str], fields: str, workers: int):
    """Yield (batch values, results) for OR-filter batches of up to BATCH values."""
    batches = [values[i:i + BATCH] for i in range(0, len(values), BATCH)]

    def fetch(batch):
        """One OR-filter request; a rejected batch is split in half until the bad value is alone."""
        try:
            res = client.get('works', filter=f'{filter_key}:' + '|'.join(batch), select=fields, **{'per-page': 200})
            return res.get('results', [])
        except RuntimeError as exc:
            if 'no answer from OpenAlex' in str(exc):
                logger.warning('  %s; retrying this batch once', exc)
                res = client.get('works', filter=f'{filter_key}:' + '|'.join(batch), select=fields, **{'per-page': 200})
                return res.get('results', [])
            if 'HTTP 400' not in str(exc):
                raise
            if len(batch) == 1:
                logger.warning('  skipped %s %s: OpenAlex rejected it', filter_key, batch[0])
                return []
            mid = len(batch) // 2
            return fetch(batch[:mid]) + fetch(batch[mid:])

    def one(batch):
        return batch, fetch(batch)

    t0 = time.time()
    logger.info('  %s: %d requests of up to %d ids, %d at a time', filter_key, len(batches), BATCH, workers)
    stop = threading.Event()
    done_count = [0]

    def heartbeat():
        while not stop.wait(30):
            logger.info('  still working: %d/%d %s requests done after %.0fs', done_count[0], len(batches),
                        filter_key, time.time() - t0)

    threading.Thread(target=heartbeat, daemon=True).start()
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(one, b) for b in batches]
        for i, fut in enumerate(as_completed(futures), 1):
            done_count[0] = i
            try:
                yield fut.result()
            except BaseException:
                stop.set()
                for f in futures:
                    f.cancel()
                raise
            if i % 50 == 0 or i == len(batches):
                rate = i / max(time.time() - t0, 1e-6)
                logger.info('  %s: %d/%d requests (%.1f/s, ETA %.0f min)', filter_key, i, len(batches),
                            rate, (len(batches) - i) / rate / 60)
    stop.set()


# --------------------------------------------------------------------------- steps

def step_tools(client: Client, workers: int, limit: int | None) -> None:
    papers = duckdb.sql(f"SELECT DISTINCT work_id, doi, pmid, pmcid FROM '{DATA}/tool_papers.parquet'").df()
    cache = CACHE / 'tool_papers.jsonl'
    done = {r['work_id'] for r in read_cache(cache)}
    todo = papers[~papers.work_id.isin(done)]
    if limit:
        todo = todo.head(limit)
    logger.info('Tool papers: %d cached, %d to fetch', len(done), len(todo))

    found: dict[int, dict] = {}
    for key, col, fmt in (('doi', 'doi', lambda v: v), ('pmid', 'pmid', str), ('pmcid', 'pmcid', str)):
        rest = todo[~todo.work_id.isin(found) & todo[col].notna()]
        if key == 'doi':  # ':' ',' '|' break OpenAlex's filter syntax (old Wiley DOIs): use PMID/PMCID
            rest = rest[~rest.doi.str.contains(r'[,|:]', regex=True)]
        if rest.empty:
            continue
        lookup = {fmt(v).lower(): w for w, v in zip(rest.work_id, rest[col])}
        for _, results in run_batches(client, key, list(lookup), PAPER_FIELDS, workers):
            for r in results:
                ids = r.get('ids') or {}
                cand = {'doi': norm_doi(r.get('doi') or ''), 'pmid': short_id(ids.get('pmid')),
                        'pmcid': short_id(ids.get('pmcid'))}[key]
                w = lookup.get((cand or '').lower())
                if w is not None and w not in found:
                    found[w] = r
    with cache.open('a') as fh:
        for w in todo.work_id:
            r = found.get(w)
            fh.write(json.dumps({'work_id': int(w), 'found': r is not None, 'openalex_id': short_id(r['id']) if r else None,
                                 'year': r.get('publication_year') if r else None,
                                 'cited_by_count': r.get('cited_by_count') if r else None,
                                 'counts_by_year': r.get('counts_by_year') if r else None,
                                 'authors': compact_authorships(r) if r else []}) + '\n')
    rows = read_cache(cache)
    logger.info('Tool papers found in OpenAlex: %d of %d', sum(r['found'] for r in rows), len(rows))
    df = pd.DataFrame(rows)
    df['counts_by_year'] = df.counts_by_year.apply(json.dumps)
    duckdb.register('df', df)
    duckdb.sql(f"COPY df TO '{DATA}/oa_tool_papers.parquet' (FORMAT parquet)")


def step_articles(client: Client, workers: int, limit: int | None) -> None:
    pmcids = [r[0] for r in duckdb.sql(f"SELECT DISTINCT pmcid FROM '{DATA}/tool_use.parquet' ORDER BY 1").fetchall()]
    cache = CACHE / 'articles.jsonl'
    done = {r['pmcid'] for r in read_cache(cache)}
    todo = [p for p in pmcids if p not in done]
    if limit:
        todo = todo[:limit]
    logger.info('Articles: %d cached, %d to fetch', len(done), len(todo))
    with cache.open('a') as fh:
        for batch, results in run_batches(client, 'pmcid', todo, ARTICLE_FIELDS, workers):
            got = {}
            for r in results:
                p = short_id((r.get('ids') or {}).get('pmcid'))
                if p:
                    got[p] = r
            for p in batch:
                r = got.get(p)
                auth = compact_authorships(r) if r else []
                fh.write(json.dumps({
                    'pmcid': p, 'found': r is not None, 'openalex_id': short_id(r['id']) if r else None,
                    'year': r.get('publication_year') if r else None,
                    'author_ids': sorted({a['author_id'] for a in auth if a['author_id']}),
                    'institution_ids': sorted({i for a in auth for i in a['institution_ids']}),
                    'institution_roots': sorted({i for a in auth for i in a['institution_roots']}),
                    'countries': sorted({c for a in auth for c in a['countries']}),
                }) + '\n')
            fh.flush()
    rows = read_cache(cache)
    logger.info('Articles found in OpenAlex: %d of %d', sum(r['found'] for r in rows), len(rows))
    duckdb.register('df', pd.DataFrame(rows))
    duckdb.sql(f"COPY df TO '{DATA}/oa_articles.parquet' (FORMAT parquet)")


def step_authors(client: Client, workers: int, limit: int | None) -> None:
    papers = read_cache(CACHE / 'tool_papers.jsonl')
    ids = sorted({a['author_id'] for p in papers for a in p['authors']
                  if a['author_id'] and (a['position'] in ('first', 'last') or a['is_corresponding'])})
    cache = CACHE / 'authors.jsonl'
    done = {r['author_id'] for r in read_cache(cache)}
    todo = [a for a in ids if a not in done]
    if limit:
        todo = todo[:limit]
    logger.info('Authors: %d cached, %d to fetch', len(done), len(todo))

    def one(aid):
        rec = client.get(f'authors/{aid}', select='id,display_name,works_count,cited_by_count,summary_stats')
        years = client.get('works', filter=f'author.id:{aid}', group_by='publication_year')
        return {'author_id': aid, 'name': rec.get('display_name'), 'works_count': rec.get('works_count'),
                'cited_by_count': rec.get('cited_by_count'),
                'h_index': (rec.get('summary_stats') or {}).get('h_index'),
                'years': {g['key']: g['count'] for g in years.get('group_by', []) if g['key'].isdigit()}}

    with cache.open('a') as fh, ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(one, a) for a in todo]
        for i, fut in enumerate(as_completed(futures), 1):
            try:
                fh.write(json.dumps(fut.result()) + '\n')
            except Exception as exc:  # noqa: BLE001 - log and retry on the next run
                logger.warning('author failed: %s', exc)
            if i % 100 == 0:
                fh.flush()
                logger.info('  authors: %d/%d', i, len(todo))
    rows = read_cache(cache)
    duckdb.register('a', pd.DataFrame([{k: r[k] for k in ('author_id', 'name', 'works_count', 'cited_by_count', 'h_index')} for r in rows]))
    duckdb.sql(f"COPY a TO '{DATA}/oa_authors.parquet' (FORMAT parquet)")
    duckdb.register('y', pd.DataFrame([{'author_id': r['author_id'], 'year': int(y), 'n_works': n}
                                       for r in rows for y, n in r['years'].items()]))
    duckdb.sql(f"COPY y TO '{DATA}/oa_author_years.parquet' (FORMAT parquet)")
    logger.info('Authors fetched: %d', len(rows))


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--step', choices=['tools', 'articles', 'authors', 'all'], default='all')
    p.add_argument('--contact', default='', help="Email for OpenAlex's polite pool.")
    p.add_argument('--api-key', default=os.environ.get('OPENALEX_API_KEY', ''),
                   help='OpenAlex API key (default: the OPENALEX_API_KEY environment variable).')
    p.add_argument('--workers', type=int, default=3)
    p.add_argument('--max-rps', type=float, default=8.0)
    p.add_argument('--limit', type=int, default=None, help='Only fetch this many new items per step (testing).')
    p.add_argument('--ca-bundle', default=None, help='PEM CA bundle, e.g. behind a TLS-inspecting proxy.')
    args = p.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
    logging.getLogger('urllib3').setLevel(logging.ERROR)  # its retry warnings print the URL, api_key included
    dl._CA_BUNDLE = args.ca_bundle
    CACHE.mkdir(parents=True, exist_ok=True)
    client = Client(args.contact, args.api_key, args.max_rps)
    if args.step in ('tools', 'all'):
        step_tools(client, args.workers, args.limit)
    if args.step in ('authors', 'all'):
        step_authors(client, args.workers, args.limit)
    if args.step in ('articles', 'all'):
        step_articles(client, args.workers, args.limit)
    return 0


if __name__ == '__main__':
    sys.exit(main())
