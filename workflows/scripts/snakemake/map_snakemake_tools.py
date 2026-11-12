"""
Translates the Snakemake DAGs (tool nodes) to bio.tools IDs, like map_nextflow_tools.py does for Nextflow.

Reads workflows/data/snakemake/dags/<owner>__<repo>/{nodes,edges}.csv (from extract_snakemake_dags.py), the
hand-curated snakemake_curated_map.tsv next to this script, the bio.tools dump
(directed-graph/biotools/data/biotools_raw*.jsonl) and the Galaxy Codex (workflows/data/mapping/galaxy/).
The methods are those of workflows/scripts/nextflow/nfdag/mapping.py, first match wins:
curated > wrapper_homepage (the homepage of a snakemake-wrapper is a bio.tools homepage) > biotools_id >
biotools_name > biotools_url > codex > prefix (rule names: bwa_mem -> bwa).

It outputs (in workflows/):
  - data/mapping/snakemake/snakemake_tool_map.tsv    one row per tool name, with the ID and the method used
  - data/mapping/snakemake/snakemake_to_curate.tsv   tools without ID or mapped only by prefix, most used first, with
                                           suggestions; the first four columns are those of the curated TSV
  - data/biotools/snakemake/biotools_pairs.csv   one row per workflow and pair     } all workflows; the columns of
  - data/biotools/snakemake/biotools_edges.csv   one row per pair, over workflows  } the Galaxy and Nextflow files,
  - data/biotools/snakemake/biotools_nodes.csv   one row per bio.tools ID          } confidence 1 (best) - 3
  - data/biotools/snakemake/biotools_{pairs,edges,nodes}_standardized.csv   the same, standardized workflows only
  - data/biotools/snakemake/workflows/<owner>__<repo>.csv   the bio.tools DAG of each workflow (for Cytoscape)

Usage (from the repository root, after extract_snakemake_dags.py):
    python3 workflows/scripts/snakemake/map_snakemake_tools.py
"""

import argparse
import csv
import sys
from pathlib import Path
from typing import Final

from smkdag import DAG_DIR, DOWNLOAD_DIR, MAP_DIR, OUT_DIR  # also puts workflows/scripts/nextflow on the path

# isort: split
from map_nextflow_tools import load_lookup
from nfdag.mapping import BiotoolsPair, Lookup, ToolMap, translate
from nfdag.report import (
    ToolStats,
    count_lost,
    print_summary,
    record,
    write_biotools,
    write_to_curate,
    write_tool_map,
    write_workflow_csvs,
)

HERE: Final = Path(__file__).resolve().parent
METHODS: Final = (
    "curated",
    "wrapper_homepage",
    "biotools_id",
    "biotools_name",
    "biotools_url",
    "codex",
    "prefix",
)
LABELS: Final = ("rules", "files")  # the context columns of the tool map and the to-curate list


def _workflow_id(folder: str) -> str:
    """owner__repo -> owner/repo"""
    return folder.replace("__", "/", 1)


def map_workflows(lookup: Lookup, dag_dir: Path) -> tuple[dict[str, ToolStats], dict[str, list[BiotoolsPair]]]:
    """Map all tool nodes; returns per tool name its statistics, and per workflow its bio.tools pairs."""
    tools: dict[str, ToolStats] = {}
    pairs: dict[str, list[BiotoolsPair]] = {}
    for folder in sorted(p for p in dag_dir.iterdir() if (p / "nodes.csv").is_file()):
        with (folder / "nodes.csv").open(encoding="utf-8") as f:
            nodes = list(csv.DictReader(f))
        with (folder / "edges.csv").open(encoding="utf-8") as f:
            edges = [(e["source"], e["target"], e["kind"]) for e in csv.DictReader(f)]
        mapped: dict[str, ToolMap] = {}
        for row in nodes:
            m = lookup.map(
                row["tool"], row["tool_source"], "", row["rule"], row["declared_tool"], homepage=row["homepage"]
            )
            mapped[row["node"]] = m
            record(
                tools,
                row["tool"],
                m,
                folder.name,
                unit=row["rule"],
                file=f"{_workflow_id(folder.name)}:{row['file']}",
                containers=row["container"],
                conda=row["conda_packages"],
            )
        pairs[folder.name], _ = translate(mapped, edges)
        count_lost(tools, mapped, {r["node"]: r["tool"].lower() for r in nodes}, edges)
    return tools, pairs


def standardized_folders(summary: Path) -> set[str]:
    with summary.open(encoding="utf-8") as f:
        return {r["folder"] for r in csv.DictReader(f) if r["standardized"] == "true"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dags", type=Path, default=DAG_DIR, help="the DAGs (default: workflows/data/snakemake/dags)")
    parser.add_argument("--map-dir", type=Path, default=MAP_DIR, help="default: workflows/data/mapping/snakemake")
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR, help="default: workflows/data/biotools/snakemake")
    parser.add_argument("--curated", type=Path, default=HERE / "snakemake_curated_map.tsv", help="the curated TSV")
    args = parser.parse_args()
    dag_dir, map_dir, out_dir = args.dags, args.map_dir, args.out_dir
    if not dag_dir.is_dir():
        sys.exit(f"{dag_dir} not found - run extract_snakemake_dags.py first")
    map_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)
    index, with_paper, lookup = load_lookup(args.curated)

    tools, pairs = map_workflows(lookup, dag_dir)
    write_tool_map(map_dir / "snakemake_tool_map.tsv", tools, index, with_paper, LABELS)
    to_curate = map_dir / "snakemake_to_curate.tsv"
    print_summary(tools, METHODS, to_curate, write_to_curate(to_curate, tools, index, LABELS))
    write_workflow_csvs(out_dir / "workflows", pairs)
    write_biotools(out_dir, pairs, index, with_paper, workflow_id=_workflow_id)
    standardized = standardized_folders(DOWNLOAD_DIR / "snakemake_workflows.csv")
    chosen = {name: p for name, p in pairs.items() if name in standardized}
    write_biotools(out_dir, chosen, index, with_paper, workflow_id=_workflow_id, suffix="_standardized")


if __name__ == "__main__":
    main()
