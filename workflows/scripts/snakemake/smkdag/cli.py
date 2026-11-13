"""Command line: build the DAG of every downloaded workflow (or of the --only ones)."""

import argparse
import csv
import shutil
import sys
import time
import traceback
from pathlib import Path
from typing import Final

from nfdag.steps import read_excluded

from . import DAG_DIR, DOWNLOAD_DIR
from . import __doc__ as package_doc
from .build import NoSnakefileError, build
from .check import Check, compare
from .output import OverviewRow, write_dag, write_error, write_overview

HERE: Final = Path(__file__).resolve().parent.parent
OVERVIEW_NAME: Final = "snakemake_dags.csv"
SLOW: Final = 20.0  # seconds: a workflow taking longer gets its own line


def _duration(seconds: float) -> str:
    minutes, sec = divmod(round(seconds), 60)
    hours, minutes = divmod(minutes, 60)
    return f"{hours}h{minutes:02d}m" if hours else f"{minutes}m{sec:02d}s" if minutes else f"{sec}s"


class Progress:
    """One progress line on stderr: rewritten in place in a terminal (every 5 s), a new line every 30 s in a log."""

    def __init__(self, total: int, enabled: bool) -> None:
        self.total, self.enabled = total, enabled
        self.tty = sys.stderr.isatty()
        self.interval = 5.0 if self.tty else 30.0
        self.start = self.last = time.monotonic()
        self.done = self.failed = 0

    def step(self, name: str, seconds: float, failed: bool) -> None:
        self.done += 1
        self.failed += failed
        if not self.enabled:
            return
        if seconds > SLOW:
            self.clear()
            print(f"  slow: {name} {_duration(seconds)}", file=sys.stderr)
        now = time.monotonic()
        if now - self.last >= self.interval or self.done == self.total:
            self.last = now
            elapsed = now - self.start
            left = elapsed / self.done * (self.total - self.done)
            line = (
                f"  {self.done}/{self.total} workflows, {self.failed} failed, {_duration(elapsed)}, "
                f"~{_duration(left)} left (now: {name})"
            )
            if self.tty:
                print(f"\r\033[K{line}", end="" if self.done < self.total else "\n", file=sys.stderr, flush=True)
            else:
                print(line, file=sys.stderr, flush=True)

    def clear(self) -> None:
        if self.tty:
            print("\r\033[K", end="", file=sys.stderr)


def _workflows(download: Path, only: list[str] | None) -> list[tuple[Path, str]]:
    """(repository folder, standardized) of the downloaded workflows that were cloned."""
    summary = download / "snakemake_workflows.csv"
    if not summary.is_file():
        sys.exit(f"{summary} not found: run 01_download_snakemake_catalog.py first")
    with summary.open(encoding="utf-8") as f:
        rows = [r for r in csv.DictReader(f) if not r["clone_error"]]
    chosen = [(download / "repos" / r["folder"], r["standardized"]) for r in rows]
    return [c for c in chosen if c[0].name in only] if only else chosen


def _summary(rows: list[OverviewRow], out_dir: Path) -> str:
    done = [r for r in rows if r.rules is not None]
    no_snakefile = sum(1 for r in rows if r.snakefile == "none")
    checked = [r for r in done if r.catalog_edges]
    found = sum(r.catalog_found or 0 for r in checked)
    expected = sum(r.catalog_edges or 0 for r in checked)
    extra = sum(r.catalog_extra or 0 for r in checked)
    inputs, unresolved = sum(r.inputs or 0 for r in done), sum(r.unresolved_inputs or 0 for r in done)
    text = (
        f"Wrote {len(done)} DAGs to {out_dir}: {sum(r.rules or 0 for r in done)} rules, "
        f"{sum(r.tool_nodes or 0 for r in done)} tool nodes, {sum(r.rule_edges or 0 for r in done)} rule edges; "
        f"{unresolved} of {inputs} rule inputs unresolved; {sum(1 for r in done if not r.rule_edges)} workflows "
        f"without edges; {no_snakefile} without Snakefile; {len(rows) - len(done) - no_snakefile} failed (error.log)"
    )
    if checked:
        text += (
            f"\nCatalog rule graphs ({len(checked)} workflows): {found} of {expected} edges found "
            f"(recall {found / max(expected, 1):.1%}), {extra} extra (precision {found / max(found + extra, 1):.1%})"
        )
    return text


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=package_doc, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--only", action="append", help="only this workflow folder, e.g. owner__repo (repeatable)")
    parser.add_argument("--verbose", action="store_true", help="print the warnings of each workflow")
    parser.add_argument("--quiet", action="store_true", help="no progress line")
    parser.add_argument(
        "--download-dir", type=Path, default=DOWNLOAD_DIR, help="default: workflows/data/snakemake/download"
    )
    parser.add_argument("--out-dir", type=Path, default=DAG_DIR, help="default: workflows/data/snakemake/dags")
    parser.add_argument(
        "--exclude-file",
        type=Path,
        default=HERE / "snakemake_excluded_tools.txt",
        help="executables that never become tool nodes (default: snakemake_excluded_tools.txt)",
    )
    args = parser.parse_args(argv)

    download, out_dir = args.download_dir, args.out_dir
    excluded = read_excluded(args.exclude_file)
    if not args.only and out_dir.exists():
        shutil.rmtree(out_dir)  # no DAGs left over from workflows no longer downloaded
    out_dir.mkdir(parents=True, exist_ok=True)
    rows: list[OverviewRow] = []
    workflows = _workflows(download, args.only)
    progress = Progress(len(workflows), enabled=not (args.quiet or args.verbose))
    for repo, standardized in workflows:
        started = time.monotonic()
        try:
            result = build(repo, download / "wrappers", excluded)
        except NoSnakefileError:
            rows.append(OverviewRow(repo.name, standardized, snakefile="none"))
            progress.step(repo.name, time.monotonic() - started, failed=False)
            continue
        except Exception:  # one unparsable workflow must not stop the run: keep the traceback
            path = write_error(out_dir / repo.name, traceback.format_exc())
            progress.clear()
            print(f"  failed {repo.name}: see {path}", file=sys.stderr)
            rows.append(OverviewRow(repo.name, standardized))
            progress.step(repo.name, time.monotonic() - started, failed=True)
            continue
        graph = download / "rulegraphs" / f"{repo.name}.dot"
        check: Check | None = None
        if graph.is_file():
            check = compare(graph, set(result.rules), set(result.rule_edges))
        write_dag(out_dir / repo.name, result, check)
        rows.append(OverviewRow.of(repo.name, standardized, result, check))
        rows[-1].seconds = round(time.monotonic() - started, 2)
        progress.step(repo.name, rows[-1].seconds, failed=False)
        if args.verbose:
            print(f"{repo.name}: {len(result.rules)} rules, {len(result.nodes)} tool nodes, {len(result.edges)} edges")
            for warning in result.warnings:
                print(f"  {warning}")
    write_overview(out_dir / OVERVIEW_NAME, rows)
    print(_summary(rows, out_dir))
