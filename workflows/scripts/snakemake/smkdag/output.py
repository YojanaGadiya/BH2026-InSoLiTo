"""Writing the DAGs: nodes.csv, edges.csv, rule_edges.csv and warnings.log per workflow, and the overview."""

import csv
from dataclasses import asdict, dataclass, fields
from pathlib import Path
from typing import Final, Self

from .build import Result
from .check import Check

NODE_FIELDS: Final = (
    "node",
    "rule",
    "step",
    "tool",
    "tool_source",  # wrapper | executable | conda | container | process (the rule name)
    "declared_tool",
    "homepage",
    "rule_kind",  # wrapper | shell | script | run | notebook | none
    "wrapper",
    "file",
    "line",
    "conda_env",
    "conda_packages",
    "container",
    "executables",
    "scripts",
    "inputs",
    "unresolved_inputs",
)
EDGE_FIELDS: Final = ("source", "target", "kind", "via")


@dataclass(slots=True)
class OverviewRow:
    """One row of the overview CSV; the counts are None for a workflow that failed."""

    workflow: str
    standardized: str
    snakefile: str = ""
    files: int | None = None
    rules: int | None = None  # DAG nodes
    target_rules: int | None = None
    tool_nodes: int | None = None
    rule_edges: int | None = None
    edges: int | None = None  # between tool nodes
    inputs: int | None = None
    unresolved_inputs: int | None = None
    warnings: int | None = None
    catalog_edges: int | None = None  # the check against the catalog's rule graph
    catalog_found: int | None = None
    catalog_extra: int | None = None
    recall: float | None = None
    precision: float | None = None

    @classmethod
    def of(cls, workflow: str, standardized: str, result: Result, check: Check | None) -> Self:
        row = cls(
            workflow,
            standardized,
            result.wf.main.relative_to(result.wf.root).as_posix(),
            len(result.wf.files),
            len(result.rules),
            len(result.targets),
            len(result.nodes),
            len(result.rule_edges),
            len(result.edges),
            result.inputs,
            result.unresolved,
            len(result.warnings),
        )
        if check is not None:
            row.catalog_edges, row.catalog_found, row.catalog_extra = check.catalog_edges, check.found, check.extra
            row.recall, row.precision = round(check.recall, 3), round(check.precision, 3)
        return row


def _rel(result: Result, path: Path) -> str:
    try:
        return path.resolve().relative_to(result.wf.root.resolve()).as_posix()
    except ValueError:
        return str(path)


def write_dag(out_dir: Path, result: Result, check: Check | None) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    with (out_dir / "nodes.csv").open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=NODE_FIELDS)
        writer.writeheader()
        for node, (rule_name, index, step) in result.nodes.items():
            rule = result.rules[rule_name]
            tools = rule.tools
            writer.writerow(
                {
                    "node": node,
                    "rule": rule_name,
                    "step": index,
                    "tool": step.tool,
                    "tool_source": step.source,
                    "declared_tool": step.declared,
                    "homepage": tools.homepage if step.source == "wrapper" else "",
                    "rule_kind": tools.kind,
                    "wrapper": tools.wrapper,
                    "file": _rel(result, rule.rule.file),
                    "line": rule.rule.line,
                    "conda_env": tools.conda_env,
                    "conda_packages": ";".join(tools.conda_packages),
                    "container": tools.container,
                    "executables": ";".join(tools.executables),
                    "scripts": ";".join(tools.scripts),
                    "inputs": rule.inputs.items,
                    "unresolved_inputs": rule.inputs.unresolved,
                }
            )
    with (out_dir / "edges.csv").open("w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out)
        writer.writerow(EDGE_FIELDS)
        for (a, b), (kind, via) in sorted(result.edges.items()):
            writer.writerow([a, b, kind, via])
    with (out_dir / "rule_edges.csv").open("w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out)
        writer.writerow(["source", "target", "via"])
        for (a, b), vias in sorted(result.rule_edges.items()):
            writer.writerow([a, b, ";".join(sorted(vias))])
    lines = list(result.warnings)
    if check is not None and check.missing:
        lines += [f"catalog rule graph edge not found: {a} -> {b}" for a, b in check.missing]
    (out_dir / "warnings.log").write_text("".join(" ".join(w.split()) + "\n" for w in lines), encoding="utf-8")


def write_error(out_dir: Path, text: str) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "error.log"
    path.write_text(text, encoding="utf-8")
    return path


def write_overview(path: Path, rows: list[OverviewRow]) -> None:
    with path.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=[f.name for f in fields(OverviewRow)])
        writer.writeheader()
        writer.writerows(asdict(row) for row in rows)
