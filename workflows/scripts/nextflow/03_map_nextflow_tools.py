"""
Translates the Nextflow DAGs (tool nodes) to bio.tools IDs, like map_galaxy_tools.py does for Galaxy.

Reads workflows/data/nextflow/dags/<id>_v<version>/{nodes,edges}.csv (from 02_extract_nextflow_dags.py) and the
hand-curated nextflow_curated_map.tsv next to this script. The reference is the bio.tools dump
(directed-graph/biotools/data/biotools_raw*.jsonl). The mapping methods are described in nfdag/mapping.py:
curated > meta (nf-core meta.yml) > biotools_id > biotools_name > biotools_url > codex (Galaxy Codex) > prefix.

It outputs (in workflows/):
  - data/mapping/nextflow/nextflow_tool_map.tsv   one row per tool name, with the ID and the method used
  - data/mapping/nextflow/nextflow_to_curate.tsv  tools without ID or mapped only by prefix, most used first,
                                                  with suggestions. The first four columns are those of the
                                                  curated TSV: fill in action (map/drop) and biotoolsID, and
                                                  copy the row into nextflow_curated_map.tsv
  - data/biotools/nextflow/biotools_pairs.csv     one row per workflow and pair       } the same columns as the
  - data/biotools/nextflow/biotools_edges.csv     one row per pair, over workflows    } Galaxy files, with bio.tools
  - data/biotools/nextflow/biotools_nodes.csv     one row per bio.tools ID            } names for ToolShed names
  - data/biotools/nextflow/biotools_{pairs,edges,nodes}_without_nfcore.csv   the same without the workflows of
                                                  the WorkflowHub space nf-core (column nf_core of
                                                  nextflow_workflows.csv), which another tool gets from nf-core
  - data/biotools/nextflow/workflows/<id>_v<version>.csv   the bio.tools DAG of each workflow (for Cytoscape)

Usage (from the repository root, after 02_extract_nextflow_dags.py):
    python3 workflows/scripts/nextflow/03_map_nextflow_tools.py
"""

import argparse
import csv
import sys
from pathlib import Path

from nfdag.lookup import load_lookup
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

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DATA = ROOT / "workflows" / "data"
DAG_DIR = DATA / "nextflow" / "dags"
CURATED_FILE = HERE / "nextflow_curated_map.tsv"  # tracked, next to this script
MAP_DIR = DATA / "mapping" / "nextflow"
OUT_DIR = DATA / "biotools" / "nextflow"
SUMMARY_FILE = DATA / "nextflow" / "download" / "nextflow_workflows.csv"  # 01_download_workflowhub_nextflow.py
METHODS = ("curated", "meta", "biotools_id", "biotools_name", "biotools_url", "codex", "prefix")
LABELS = ("processes", "modules")  # the context columns of the tool map and the to-curate list


def _workflow_dirs(dag_dir: Path) -> list[Path]:
    dirs = [p for p in dag_dir.iterdir() if (p / "nodes.csv").is_file()]
    return sorted(dirs, key=lambda p: int(p.name.split("_")[0]))


def _workflow_id(folder_name: str) -> str:
    """462_v3 -> 462"""
    return folder_name.split("_v", maxsplit=1)[0]


def map_workflows(lookup: Lookup, dag_dir: Path) -> tuple[dict[str, ToolStats], dict[str, list[BiotoolsPair]]]:
    """Map all tool nodes; returns per tool name its statistics, and per workflow its bio.tools pairs."""
    tools: dict[str, ToolStats] = {}
    pairs: dict[str, list[BiotoolsPair]] = {}
    for folder in _workflow_dirs(dag_dir):
        with (folder / "nodes.csv").open(encoding="utf-8") as f:
            nodes = list(csv.DictReader(f))
        with (folder / "edges.csv").open(encoding="utf-8") as f:
            reader = csv.DictReader(f)
            if "kind" not in (reader.fieldnames or []) or (nodes and "tool" not in nodes[0]):
                sys.exit(
                    f"{folder} holds process-level DAGs from an older version: run 02_extract_nextflow_dags.py again"
                )
            edges = [(e["source"], e["target"], e["kind"]) for e in reader]
        mapped: dict[str, ToolMap] = {}
        for row in nodes:
            m = lookup.map(
                row["tool"], row["tool_source"], row["meta_biotools_id"], row["process"], row["declared_tool"]
            )
            mapped[row["node"]] = m
            record(
                tools,
                row["tool"],
                m,
                folder.name,
                unit=row["process"],
                file=row["module"],
                containers=row["containers"],
                conda=row["conda_packages"],
            )
        pairs[folder.name], _ = translate(mapped, edges)
        count_lost(tools, mapped, {r["node"]: r["tool"].lower() for r in nodes}, edges)
    return tools, pairs


def load_nfcore_ids() -> set[str] | None:
    """IDs of the workflows in the WorkflowHub space nf-core (column nf_core of nextflow_workflows.csv),
    None if the summary has no such column."""
    if not SUMMARY_FILE.exists():
        return None
    with SUMMARY_FILE.open(encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    if not rows or "nf_core" not in rows[0]:
        return None
    return {r["id"] for r in rows if r["nf_core"] == "true"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dags", type=Path, default=DAG_DIR, help="the DAGs (default: data/nextflow)")
    parser.add_argument("--curated", type=Path, default=CURATED_FILE, help="the curated TSV")
    args = parser.parse_args()
    if not args.dags.is_dir():
        sys.exit(f"{args.dags} not found - run 02_extract_nextflow_dags.py first")
    MAP_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    index, with_paper, lookup = load_lookup(args.curated)

    tools, pairs = map_workflows(lookup, args.dags)
    write_tool_map(MAP_DIR / "nextflow_tool_map.tsv", tools, index, with_paper, LABELS)
    to_curate = MAP_DIR / "nextflow_to_curate.tsv"
    print_summary(tools, METHODS, to_curate, write_to_curate(to_curate, tools, index, LABELS))
    write_workflow_csvs(OUT_DIR / "workflows", pairs)
    write_biotools(OUT_DIR, pairs, index, with_paper, workflow_id=_workflow_id)
    nfcore = load_nfcore_ids()
    if nfcore is None:
        print(f"No nf_core column in {SUMMARY_FILE}: run 01_download_workflowhub_nextflow.py --summary-only")
    else:  # the nf-core pipelines also come directly from nf-core, with another tool
        others = {name: p for name, p in pairs.items() if _workflow_id(name) not in nfcore}
        write_biotools(OUT_DIR, others, index, with_paper, workflow_id=_workflow_id, suffix="_without_nfcore")


if __name__ == "__main__":
    main()
