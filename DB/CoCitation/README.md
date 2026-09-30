# Co-citation corpus: Europe PMC full-text download

`download_europepmc.py` downloads Open Access JATS full-text XML from the
[Europe PMC REST API](https://europepmc.org/RestfulWebService) for one or more topics.
By default it downloads proteomics and metabolomics.

```bash
pip install -r requirements.txt

# Smoke test: 50 articles per topic
python download_europepmc.py --limit 50 --contact you@example.org

# Full run (resumable: re-run the same command after an interruption)
python download_europepmc.py --topics proteomics metabolomics --contact you@example.org
```

Query per topic (override with `--query-template`):
`({topic}) AND OPEN_ACCESS:Y AND HAS_REFLIST:Y`

## Output (`data/`, git-ignored)

| Path | Contents |
|---|---|
| `search/<topic>.jsonl` | One search hit per line: pmcid, pmid, doi, title, journal, year, pubType |
| `xml/PMCxxx/PMCxxxxxxx.xml.gz` | Gzipped JATS full text, 10k files per folder |
| `manifest.tsv` | One row per unique article: IDs, topics it matched, pub type, download status |

`status` in the manifest is `ok`, `not_available` (Europe PMC returned 404) or `error: ...`.
Re-running retries everything that isn't `ok` unless you pass `--errors skip`.
Stopping a run with Ctrl+C or `kill` saves the manifest first.

## Options

- `--step search|fetch|all`: run only one step. `fetch` reuses the existing search results.
- `--errors retry|skip|only`: what to do with papers that failed in an earlier run (listed in `data/errors.tsv`).
  `skip` leaves them for later, and `only` retries just those.
- `--workers` (default 6) and `--max-rps` (default 8): parallel downloads, and a global cap on requests per second.
  Keep these modest. For the full 7.3M corpus, use the bulk FTP/S3 dumps instead of this API.

## TLS errors behind a corporate proxy

`CERTIFICATE_VERIFY_FAILED: self-signed certificate in certificate chain` means a
TLS-inspecting proxy sits between you and EBI. With `truststore` installed (it's in
`requirements.txt`, Python 3.10+) the script uses your OS certificate store, which
already trusts the proxy. Alternatively, pass `--ca-bundle /path/to/corp-root.pem`.
Don't disable verification.
