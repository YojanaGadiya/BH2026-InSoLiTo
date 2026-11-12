"""
Download the workflows of the Snakemake Workflow Catalog (https://snakemake.github.io/snakemake-workflow-catalog).

Step 1 (catalog):  the catalog's data.json, one entry per workflow           -> data/download/catalog.json
                   (full_name, standardized, stars, topics, wrappers, and for some a Snakemake rule graph)
Step 2 (clone):    each selected repository: a shallow clone of its latest release (else the default branch),
                   with only the files needed to read the workflow (Snakefiles, .smk, yaml, scripts), so
                   data stored in git is not downloaded; the clone runs no code
                                                                              -> data/download/repos/<owner>__<repo>/
                   the commit, or the error                                    -> data/download/repos/<...>/.cloned_from
Step 3 (wrappers): meta.yaml and environment.yaml of every snakemake-wrapper the workflows use
                                                                    -> data/download/wrappers/<version>/<path>/
Step 4 (summary):  the catalog's rule graphs (DOT), to check the static DAGs  -> data/download/rulegraphs/<...>.dot
                   one row per downloaded workflow                            -> data/download/snakemake_workflows.csv

"Standardized" workflows follow the catalog's layout rules (.snakemake-workflow-catalog.yml, config/README.md,
usually workflow/Snakefile, workflow/rules/, workflow/envs/). Runs can be stopped and restarted: repositories
already cloned are not cloned again, and every run adds to the same folders.

Usage (from the repository root):
    python3 workflows/scripts/snakemake/download_snakemake_catalog.py --sample 50 --standardized
    python3 workflows/scripts/snakemake/download_snakemake_catalog.py --sample 50 --non-standardized
    python3 workflows/scripts/snakemake/download_snakemake_catalog.py --standardized     # all 479 standardized
    python3 workflows/scripts/snakemake/download_snakemake_catalog.py                    # everything
    python3 workflows/scripts/snakemake/download_snakemake_catalog.py --summary-only     # rebuild the CSV
    python3 workflows/scripts/snakemake/download_snakemake_catalog.py --resparse --summary-only  # new file patterns
    python3 workflows/scripts/snakemake/download_snakemake_catalog.py --workflow snakemake-workflows/rna-seq-star-deseq2
"""

import argparse
import csv
import json
import os
import random
import subprocess
import sys
import time
import urllib.error
import urllib.request
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Final

from smkdag import DOWNLOAD_DIR
from smkdag.snakefile import main_snakefile, workflow_files
from smkdag.wrappers import FILES, RAW_URL, WRAPPER_LINE, WrapperRef, wrapper_dir, wrapper_ref

HERE = Path(__file__).resolve().parent
CATALOG_URL: Final = "https://raw.githubusercontent.com/snakemake/snakemake-workflow-catalog/main/data.json"
USER_AGENT: Final = "BH2026-InSoLiTo-snakemake-catalog/0.1"
SEED: Final = 42  # the samples are the same in every run
# files a clone keeps: what is needed to read the rules, their environments and scripts
SPARSE: Final = (
    "Snakefile*",
    "snakefile*",
    "*.smk",
    "*.snakefile",
    "*.snake",
    "*.rules",
    "*.rule",
    "*.yaml",
    "*.yml",
    "*.py",
    "*.R",
    "*.r",
    "*.sh",
    "*.pl",
)
CLONE_TIMEOUT: Final = 600
SUMMARY_FIELDS: Final = (
    "full_name",
    "folder",
    "standardized",
    "stars",
    "updated",
    "latest_release",
    "topics",
    "description",
    "n_catalog_wrappers",
    "has_rulegraph",
    "ref",
    "commit",
    "clone_error",
    "snakefile",
    "n_snakefiles",
)


def folder_name(full_name: str) -> str:
    """owner/repo -> owner__repo"""
    return full_name.replace("/", "__")


def fetch(url: str) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=120) as response:
        return response.read()


def load_catalog(path: Path, refresh: bool) -> list[dict[str, Any]]:
    """The catalog entries, downloaded once (again with --refresh-catalog)."""
    if refresh or not path.exists():
        print(f"Downloading {CATALOG_URL}")
        path.write_bytes(fetch(CATALOG_URL))
    entries = json.loads(path.read_text(encoding="utf-8"))
    print(f"Catalog: {len(entries)} workflows, {sum(1 for e in entries if e.get('standardized'))} standardized")
    return entries


def select(entries: list[dict[str, Any]], standardized: bool, non_standardized: bool, sample: int | None) -> list:
    """The workflows to download. A sample is the same in every run: for standardized workflows those with a
    catalog rule graph come first (to check the static DAGs against), the rest in a fixed random order."""
    groups = []
    if standardized or not non_standardized:
        chosen = sorted((e for e in entries if e.get("standardized")), key=lambda e: e["full_name"])
        with_graph = [e for e in chosen if e.get("rulegraph")]
        rest = [e for e in chosen if not e.get("rulegraph")]
        random.Random(SEED).shuffle(rest)
        groups.append(with_graph + rest)
    if non_standardized or not standardized:
        chosen = sorted((e for e in entries if not e.get("standardized")), key=lambda e: e["full_name"])
        random.Random(SEED).shuffle(chosen)
        groups.append(chosen)
    if sample is None:
        return [e for group in groups for e in group]
    if len(groups) == 1:
        return groups[0][:sample]
    merged = sorted((e for group in groups for e in group), key=lambda e: e["full_name"])  # both: one mixed sample
    random.Random(SEED).shuffle(merged)
    return merged[:sample]


def _git(args: list[str], cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    env = {**os.environ, "GIT_TERMINAL_PROMPT": "0", "GIT_LFS_SKIP_SMUDGE": "1"}
    return subprocess.run(
        ["git", *args], cwd=cwd, env=env, capture_output=True, text=True, timeout=CLONE_TIMEOUT, check=False
    )


def clone(full_name: str, release: str | None, dest: Path) -> dict[str, str]:
    """Shallow, sparse clone of the latest release (else the default branch); returns what was cloned."""
    url = f"https://github.com/{full_name}.git"
    record = {"url": url, "ref": "", "commit": "", "error": ""}
    for ref in [release, None] if release else [None]:
        if dest.exists():
            subprocess.run(["rm", "-rf", str(dest)], check=False)
        branch = ["--branch", ref] if ref else []
        try:
            result = _git(["clone", "--depth", "1", "--filter=blob:none", "--sparse", *branch, url, str(dest)])
            if result.returncode == 0:
                checkout = _git(["sparse-checkout", "set", "--no-cone", *SPARSE], cwd=dest)
                if checkout.returncode != 0:
                    record["error"] = checkout.stderr.strip()[-300:]
                    continue
                record["ref"] = ref or "default branch"
                record["commit"] = _git(["rev-parse", "HEAD"], cwd=dest).stdout.strip()
                record["error"] = ""
                break
            record["error"] = result.stderr.strip()[-300:]
        except subprocess.TimeoutExpired:
            record["error"] = f"timeout after {CLONE_TIMEOUT} s"
    dest.mkdir(parents=True, exist_ok=True)
    (dest / ".cloned_from").write_text(json.dumps(record, indent=1), encoding="utf-8")
    return record


def resparse(repos: Path) -> None:
    """Check out the files the current SPARSE patterns name in every successful clone (fetches only those)."""
    done = [d for d in sorted(repos.iterdir()) if (r := cloned(d)) is not None and not r.get("error")]
    failed = 0
    for n, dest in enumerate(done, 1):
        try:
            result = _git(["sparse-checkout", "set", "--no-cone", *SPARSE], cwd=dest)
            failed += result.returncode != 0
        except subprocess.TimeoutExpired:
            failed += 1
        if n % 200 == 0:
            print(f"  {n}/{len(done)}")
    print(f"Re-applied the file patterns to {len(done)} clones ({failed} failed)")


def cloned(dest: Path) -> dict[str, str] | None:
    marker = dest / ".cloned_from"
    return json.loads(marker.read_text(encoding="utf-8")) if marker.is_file() else None


def used_wrappers(repo: Path, entry: dict[str, Any]) -> set[WrapperRef]:
    """The wrappers a workflow uses: `wrapper:` directives in its files, and the catalog's list."""
    refs = set()
    for path in workflow_files(repo):
        for m in WRAPPER_LINE.finditer(path.read_text(encoding="utf-8", errors="replace")):
            ref = wrapper_ref(m.group(1))
            if ref:
                refs.add(ref)
    for w in (entry.get("wrappers") or {}).values():
        ref = wrapper_ref(f"{w.get('wrapper_version', 'master')}/{w.get('wrapper_name', '')}")
        if ref:
            refs.add(ref)
    return refs


def download_wrapper(wrappers_dir: Path, ref: WrapperRef, delay: float) -> bool:
    """meta.yaml and environment.yaml of a wrapper, once; a `.missing` marker if neither exists."""
    folder = wrapper_dir(wrappers_dir, ref)
    if any((folder / f).is_file() for f in FILES) or (folder / ".missing").is_file():
        return False
    folder.mkdir(parents=True, exist_ok=True)
    found = False
    for file in FILES:
        time.sleep(delay)
        try:
            (folder / file).write_bytes(fetch(RAW_URL.format(version=ref.version, path=ref.path, file=file)))
            found = True
        except urllib.error.HTTPError as error:
            if error.code != 404:
                print(f"  {ref.version}/{ref.path}/{file}: {error}", file=sys.stderr)
        except urllib.error.URLError as error:
            print(f"  {ref.version}/{ref.path}/{file}: {error}", file=sys.stderr)
    if not found:
        (folder / ".missing").write_text("", encoding="utf-8")
    return True


def download_all(
    chosen: list[dict[str, Any]], repos: Path, wrappers: Path, graphs: Path, delay: float
) -> tuple[int, int, int]:
    """Clone the chosen workflows and fetch their wrappers; returns (new clones, failed, new wrappers)."""
    n_new = n_failed = n_wrappers = 0
    for n, entry in enumerate(chosen, 1):
        dest = repos / folder_name(entry["full_name"])
        record = cloned(dest)
        if record is None:
            record = clone(entry["full_name"], entry.get("latest_release"), dest)
            n_new += 1
            status = f"error: {record['error'][:120]}" if record["error"] else record["ref"]
            print(f"  {n}/{len(chosen)} {entry['full_name']}: {status}")
        n_failed += bool(record.get("error"))
        if entry.get("rulegraph"):
            (graphs / f"{dest.name}.dot").write_text(entry["rulegraph"], encoding="utf-8")
        if not record.get("error"):
            for ref in sorted(used_wrappers(dest, entry), key=lambda r: (r.version, r.path)):
                n_wrappers += download_wrapper(wrappers, ref, delay)
    return n_new, n_failed, n_wrappers


def summary_row(entry: dict[str, Any], repo: Path, record: dict[str, str]) -> dict[str, Any]:
    main = main_snakefile(repo) if not record.get("error") else None
    updated = entry.get("updated_at")
    return {
        "full_name": entry["full_name"],
        "folder": repo.name,
        "standardized": str(bool(entry.get("standardized"))).lower(),
        "stars": entry.get("stargazers_count", ""),
        "updated": datetime.fromtimestamp(updated, UTC).date().isoformat() if updated else "",
        "latest_release": entry.get("latest_release") or "",
        "topics": ";".join(entry.get("topics") or []),
        "description": (entry.get("description") or "").replace("\n", " "),
        "n_catalog_wrappers": len(entry.get("wrappers") or {}),
        "has_rulegraph": str(bool(entry.get("rulegraph"))).lower(),
        "ref": record.get("ref", ""),
        "commit": record.get("commit", ""),
        "clone_error": record.get("error", ""),
        "snakefile": main.relative_to(repo).as_posix() if main else "",
        "n_snakefiles": len(workflow_files(repo)) if main else 0,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--standardized", action="store_true", help="only standardized workflows")
    parser.add_argument("--non-standardized", action="store_true", help="only non-standardized workflows")
    parser.add_argument("--sample", type=int, help="a fixed sample of N workflows (the same in every run)")
    parser.add_argument("--workflow", action="append", help="only this workflow, owner/repo (repeatable)")
    parser.add_argument("--summary-only", action="store_true", help="rebuild snakemake_workflows.csv, no downloads")
    parser.add_argument("--refresh-catalog", action="store_true", help="download data.json again")
    parser.add_argument(
        "--resparse", action="store_true", help="apply the current file patterns to the repositories already cloned"
    )
    parser.add_argument("--delay", type=float, default=0.2, help="seconds between wrapper downloads")
    parser.add_argument(
        "--download-dir", type=Path, default=DOWNLOAD_DIR, help="default: workflows/data/snakemake/download"
    )
    args = parser.parse_args()

    download: Path = args.download_dir
    repos, wrappers, graphs = download / "repos", download / "wrappers", download / "rulegraphs"
    for folder in (repos, wrappers, graphs):
        folder.mkdir(parents=True, exist_ok=True)
    entries = load_catalog(download / "catalog.json", args.refresh_catalog)
    by_name = {e["full_name"]: e for e in entries}

    if args.resparse:
        resparse(repos)
    if not args.summary_only:
        if args.workflow:
            unknown = [w for w in args.workflow if w not in by_name]
            if unknown:
                sys.exit(f"Not in the catalog: {', '.join(unknown)}")
            chosen = [by_name[w] for w in args.workflow]
        else:
            chosen = select(entries, args.standardized, args.non_standardized, args.sample)
        print(f"Selected {len(chosen)} workflows")
        try:
            n_new, n_failed, n_wrappers = download_all(chosen, repos, wrappers, graphs, args.delay)
            print(f"Cloned {n_new} new repositories ({n_failed} of {len(chosen)} failed), {n_wrappers} new wrappers")
        except KeyboardInterrupt:  # stopped by hand: still write the summary of what is there
            print("\nInterrupted: writing the summary of the repositories downloaded so far")

    rows = []
    for dest in sorted(p for p in repos.iterdir() if p.is_dir()):
        record = cloned(dest)
        entry = by_name.get(dest.name.replace("__", "/", 1))
        if record is not None and entry is not None:
            rows.append(summary_row(entry, dest, record))
    summary = download / "snakemake_workflows.csv"
    with summary.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    standardized = sum(r["standardized"] == "true" for r in rows)
    print(
        f"Wrote {summary}: {len(rows)} workflows ({standardized} standardized), "
        f"{sum(1 for r in rows if r['clone_error'])} failed, {sum(1 for r in rows if r['snakefile'])} with a Snakefile"
    )


if __name__ == "__main__":
    main()
