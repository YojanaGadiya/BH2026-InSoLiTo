"""Step 8: fetch EDAM annotations (operations, input and output data and formats) from
bio.tools for the entries behind the top tools.

Only these ~420 entries are fetched, one request each, into a separate file, so the
team's bio.tools CSVs (biotools/pull_biotools.py) are not touched. Cached and resumable.

Writes patterns/data/edam_functions.parquet: one row per (tool, entry, function, kind,
EDAM term): kind is operation, input_data, input_format, output_data or output_format.

With --allpub-top N it fetches the entries of the N most-linked tools of the
all-publications graph (step 12) instead, and writes edam_functions_allpub.parquet
(same cache file, so entries already fetched are not fetched again).

Run from BH2026-InSoLiTo/:
    python directed-graph/patterns/08_fetch_edam.py
    python directed-graph/patterns/08_fetch_edam.py --allpub-top 2500
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from pathlib import Path

import duckdb
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'europepmc'))
import download_europepmc as dl  # noqa: E402  session with retries, rate limiter, TLS setup

logger = logging.getLogger('fetch_edam')
DATA = Path('directed-graph/patterns/data')
CACHE = DATA / 'biotools_records.jsonl'
API = 'https://bio.tools/api/tool/{}?format=json'


def rows_for(tool: str, entry: str, record: dict) -> list[dict]:
    out = []
    for f_idx, fn in enumerate(record.get('function') or []):
        for op in fn.get('operation') or []:
            out.append({'tool': tool, 'entry': entry, 'function': f_idx, 'kind': 'operation',
                        'uri': op.get('uri'), 'term': op.get('term')})
        for side in ('input', 'output'):
            for io in fn.get(side) or []:
                d = io.get('data') or {}
                out.append({'tool': tool, 'entry': entry, 'function': f_idx, 'kind': f'{side}_data',
                            'uri': d.get('uri'), 'term': d.get('term')})
                for fmt in io.get('format') or []:
                    out.append({'tool': tool, 'entry': entry, 'function': f_idx, 'kind': f'{side}_format',
                                'uri': fmt.get('uri'), 'term': fmt.get('term')})
    for tp in record.get('topic') or []:
        out.append({'tool': tool, 'entry': entry, 'function': -1, 'kind': 'topic', 'uri': tp.get('uri'), 'term': tp.get('term')})
    return out


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--max-rps', type=float, default=3.0)
    p.add_argument('--ca-bundle', default=None)
    p.add_argument('--allpub-top', type=int, default=0, help='use the N most-linked tools of the all-publications graph')
    args = p.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
    dl._CA_BUNDLE = args.ca_bundle

    if args.allpub_top:
        tools = duckdb.sql(f"SELECT tool, entries FROM '{DATA}/allpub_tools.parquet' "
                           f"ORDER BY strength DESC LIMIT {args.allpub_top}").df()
        out_file = DATA / 'edam_functions_allpub.parquet'
    else:
        tools = duckdb.sql(f"SELECT tool, entries FROM '{DATA}/tools.parquet'").df()
        out_file = DATA / 'edam_functions.parquet'
    pairs = [(t, e) for t, es in zip(tools.tool, tools.entries) for e in es]
    done = {}
    if CACHE.exists():
        for line in CACHE.open():
            r = json.loads(line)
            done[r['entry']] = r
    todo = [(t, e) for t, e in pairs if e not in done]
    logger.info('%d entries for %d tools: %d cached, %d to fetch', len(pairs), len(tools), len(done), len(todo))

    limiter = dl.RateLimiter(args.max_rps)
    session = dl.get_session('BH2026-InSoLiTo pattern analysis (https://github.com/inab/BH2026-InSoLiTo)')
    t0 = time.time()
    with CACHE.open('a') as fh:
        for i, (t, e) in enumerate(todo, 1):
            limiter.wait()
            try:
                resp = session.get(API.format(e), timeout=(10, 60))
                rec = {'entry': e, 'status': resp.status_code,
                       'record': resp.json() if resp.status_code == 200 else None}
            except Exception as exc:  # noqa: BLE001 - logged, retried next run
                logger.warning('%s: %s', e, type(exc).__name__)
                continue
            fh.write(json.dumps(rec) + '\n')
            done[e] = rec
            if i % 25 == 0 or i == len(todo):
                fh.flush()
                logger.info('  %d/%d (%.1f/s)', i, len(todo), i / max(time.time() - t0, 1e-6))

    rows = [row for t, e in pairs if done.get(e, {}).get('record') for row in rows_for(t, e, done[e]['record'])]
    df = pd.DataFrame(rows)
    duckdb.register('df', df)
    duckdb.sql(f"COPY df TO '{out_file}' (FORMAT parquet)")
    have = df[df.kind == 'operation'].tool.nunique() if len(df) else 0
    fmt = df[df.kind.str.endswith('format')].tool.nunique() if len(df) else 0
    logger.info('Tools with EDAM operations: %d of %d; with input/output formats: %d', have, len(tools), fmt)
    logger.info('Written %s', out_file)
    return 0


if __name__ == '__main__':
    sys.exit(main())
