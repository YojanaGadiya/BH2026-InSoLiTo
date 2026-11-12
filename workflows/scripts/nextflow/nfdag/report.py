"""Write the mapping results: the tool map, the to-curate list and the bio.tools pairs, edges and nodes.

Shared by 03_map_nextflow_tools.py and the Snakemake catalog (workflows/scripts/snakemake/03_map_snakemake_tools.py).
Each workflow is a set of tool nodes, mapped with mapping.Lookup and translated into bio.tools pairs with
mapping.translate; this module collects the statistics per tool name and writes the files.
"""

import csv
import re
from collections import Counter, defaultdict
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Final

from .mapping import CONFIDENCE, CURATED_FIELDS, BiotoolsIndex, BiotoolsPair, ToolMap, tool_key

EXAMPLES: Final = 3  # example processes / workflows per tool in the TSVs
LEVELS: Final = sorted(set(CONFIDENCE.values()))  # 1, 2, 3


@dataclass
class ToolStats:
    """What is known about one tool name over all workflows."""

    mapping: ToolMap
    workflows: set[str] = field(default_factory=set)
    nodes: int = 0
    lost_edges: int = 0
    units: dict[str, None] = field(default_factory=dict)  # ordered set: processes / rules
    files: dict[str, None] = field(default_factory=dict)  # modules / files
    containers: dict[str, None] = field(default_factory=dict)
    conda: dict[str, None] = field(default_factory=dict)


def record(  # noqa: PLR0913
    tools: dict[str, ToolStats],
    tool: str,
    mapping: ToolMap,
    workflow: str,
    *,
    unit: str,
    file: str,
    containers: str,
    conda: str,
) -> None:
    """Count one tool node of one workflow; a name mapped anywhere (e.g. by meta.yml) keeps that mapping."""
    stats = tools.setdefault(tool.lower(), ToolStats(mapping))
    if stats.mapping.action == "unmapped" and mapping.action != "unmapped":
        stats.mapping = mapping
    stats.workflows.add(workflow)
    stats.nodes += 1
    stats.units.setdefault(unit)
    stats.files.setdefault(file)
    for container in containers.split(";")[:1]:
        stats.containers.setdefault(container)
    stats.conda.setdefault(conda)


def count_lost(
    tools: dict[str, ToolStats], mapped: dict[str, ToolMap], tool_of: dict[str, str], edges: list[tuple[str, str, str]]
) -> None:
    """The edges each unmapped tool cuts, to rank the curation."""
    for a, b, _ in edges:
        for node in (a, b):
            if mapped[node].action == "unmapped":
                tools[tool_of[node]].lost_edges += 1


def _examples(values: dict[str, None]) -> str:
    return ";".join(v for v in list(values)[:EXAMPLES] if v)


def _names(index: BiotoolsIndex, ids: list[str]) -> str:
    return ";".join(index.name_of(i) for i in ids)


def _id_order(workflow_id: str) -> tuple[int, int, str]:
    return (0, int(workflow_id), "") if workflow_id.isdigit() else (1, 0, workflow_id)


def write_tool_map(
    path: Path, tools: dict[str, ToolStats], index: BiotoolsIndex, with_paper: set[str], labels: tuple[str, str]
) -> None:
    """One row per tool name; labels name the unit and file columns (processes, modules)."""
    columns = [
        "tool",
        "biotoolsID",
        "biotools_name",
        "method",
        "codex_biotoolsID",
        "match_biotoolsID",
        "n_workflows",
        "n_nodes",
        "lost_edges",
        *labels,
        "containers",
        "conda",
        "has_publication",
    ]
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(columns)
        for tool, s in sorted(tools.items(), key=lambda kv: (-len(kv[1].workflows), kv[0])):
            ids = ";".join(s.mapping.ids) if s.mapping.action == "map" else s.mapping.action
            writer.writerow(
                [
                    tool,
                    ids,
                    _names(index, s.mapping.ids),
                    s.mapping.method,
                    s.mapping.codex_id,
                    s.mapping.match_id,
                    len(s.workflows),
                    s.nodes,
                    s.lost_edges,
                    _examples(s.units),
                    _examples(s.files),
                    _examples(s.containers),
                    _examples(s.conda),
                    str(any(i in with_paper for i in s.mapping.ids)).lower(),
                ]
            )


def _context(s: ToolStats) -> list[str]:
    """Package names the tool comes with: conda packages and container tools (bioconda::bedtools=2.31 -> bedtools)."""
    names = []
    for value in [*s.conda, *s.containers]:
        for item in value.split(";"):
            name = item.rsplit("::", 1)[-1].rsplit("/", 1)[-1]
            name = re.split(r"[=:<>]", name)[0].strip()
            if name and name not in names and not name.startswith(("sha256", "data")):
                names.append(name)
    return names


def write_to_curate(path: Path, tools: dict[str, ToolStats], index: BiotoolsIndex, labels: tuple[str, str]) -> int:
    """Tools a person should look at (unmapped, or only a prefix match of the process / rule name);
    the first columns are those of the curated TSV."""
    rows = []
    for tool, s in tools.items():
        m = s.mapping
        if m.action != "unmapped" and m.method != "prefix":
            continue
        reason = "unmapped" if m.action == "unmapped" else "prefix match"
        suggestions = index.suggest(tool_key(tool), _context(s)) if m.action == "unmapped" else []
        row = [
            tool,
            "",
            ";".join(m.ids),
            reason,
            "; ".join(f"{i} ({index.name_of(i)})" for i in suggestions),
            _names(index, m.ids),
            m.codex_id,
            m.match_id,
            len(s.workflows),
            s.nodes,
            s.lost_edges,
            _examples(s.units),
            _examples(s.files),
            _examples(s.containers),
            _examples(s.conda),
        ]
        rows.append((s, row))
    rows.sort(key=lambda r: (-len(r[0].workflows), -r[0].lost_edges, r[1][0]))
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(
            [
                *CURATED_FIELDS,
                "suggestions",
                "biotools_name",
                "codex_biotoolsID",
                "match_biotoolsID",
                "n_workflows",
                "n_nodes",
                "lost_edges",
                *labels,
                "containers",
                "conda",
            ]
        )
        writer.writerows(r for _, r in rows)
    return len(rows)


def print_summary(tools: dict[str, ToolStats], methods: tuple[str, ...], to_curate: Path, n_curate: int) -> None:
    """Coverage by method, and the most used unmapped tools."""
    print(f"Tool names: {len(tools)}, tool nodes: {sum(s.nodes for s in tools.values())}")
    for method in methods:
        chosen = [s for s in tools.values() if s.mapping.method == method]
        print(f"  {method:16s} {len(chosen):5d} tools, {sum(s.nodes for s in chosen):6d} nodes")
    for action in ("drop", "unmapped"):
        chosen = [s for s in tools.values() if s.mapping.action == action]
        print(f"  {action:16s} {len(chosen):5d} tools, {sum(s.nodes for s in chosen):6d} nodes")
    kept = sum(s.nodes for s in tools.values() if s.mapping.action != "drop")
    mapped = sum(s.nodes for s in tools.values() if s.mapping.action == "map")
    print(f"Coverage: {mapped} of {kept} tool nodes not dropped are bio.tools entries ({mapped / max(kept, 1):.1%})")
    print(f"DAG edges lost at unmapped tools: {sum(s.lost_edges for s in tools.values())}")
    print(f"Wrote {to_curate} ({n_curate} tools). Most used unmapped:")
    unmapped = sorted(
        ((t, s) for t, s in tools.items() if s.mapping.action == "unmapped"), key=lambda ts: -len(ts[1].workflows)
    )
    for tool, s in unmapped[:15]:
        print(f"  {len(s.workflows):5d}  {tool}  ({_examples(s.units)})")


def write_workflow_csvs(folder: Path, pairs: dict[str, list[BiotoolsPair]]) -> None:
    """The bio.tools DAG of each workflow, for Cytoscape."""
    folder.mkdir(parents=True, exist_ok=True)
    for name, workflow_pairs in pairs.items():
        unique = sorted({(p.source, p.target, p.kind, p.source_node, p.target_node) for p in workflow_pairs})
        with (folder / f"{name}.csv").open("w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["source", "target", "kind", "source_node", "target_node"])
            writer.writerows(unique)


def _by_level(counts: Counter[int]) -> str:
    return ", ".join(f"{level}: {counts[level]}" for level in LEVELS)


def write_biotools(  # noqa: PLR0913
    folder: Path,
    pairs: dict[str, list[BiotoolsPair]],
    index: BiotoolsIndex,
    with_paper: set[str],
    *,
    workflow_id: Callable[[str], str],
    suffix: str = "",
) -> None:
    """The same three files as map_galaxy_tools.py, with the bio.tools names in place of the ToolShed names,
    over the given workflows; workflow_id gives the ID of a workflow folder, suffix is added to the file names."""
    edge_workflows: dict[tuple[str, str], set[str]] = defaultdict(set)
    edge_levels: dict[tuple[str, str], Counter[int]] = defaultdict(Counter)  # confidence -> number of workflows
    node_workflows: dict[str, set[str]] = defaultdict(set)
    node_levels: dict[str, Counter[int]] = defaultdict(Counter)  # confidence -> number of workflows
    n_pairs = 0
    with (folder / f"biotools_pairs{suffix}.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(
            ["workflow_id", "source", "target", "source_biotools_name", "target_biotools_name", "confidence"]
        )
        for name, workflow_pairs in pairs.items():
            wf_id = workflow_id(name)
            best: dict[tuple[str, str], int] = {}  # like Galaxy: once per workflow, the best of its tool pairs
            node_best: dict[str, int] = {}  # per node, the best level of its own mapping in this workflow
            for p in workflow_pairs:
                best[(p.source, p.target)] = min(p.confidence, best.get((p.source, p.target), p.confidence))
                for biotools_id, level in ((p.source, p.source_confidence), (p.target, p.target_confidence)):
                    node_best[biotools_id] = min(level, node_best.get(biotools_id, level))
            for biotools_id, level in node_best.items():
                node_levels[biotools_id][level] += 1
            for (source, target), confidence in sorted(best.items()):
                writer.writerow([wf_id, source, target, index.name_of(source), index.name_of(target), confidence])
                n_pairs += 1
                edge_workflows[(source, target)].add(wf_id)
                edge_levels[(source, target)][confidence] += 1
                node_workflows[source].add(wf_id)
                node_workflows[target].add(wf_id)
    with (folder / f"biotools_edges{suffix}.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "source",
                "target",
                "n_workflows",
                "workflow_ids",
                "source_biotools_name",
                "target_biotools_name",
                "confidence",  # the best over the workflows
                *(f"confidence{level}" for level in LEVELS),  # workflows per level
            ]
        )
        for (source, target), ids in sorted(edge_workflows.items(), key=lambda kv: (-len(kv[1]), kv[0])):
            levels = edge_levels[(source, target)]
            writer.writerow(
                [
                    source,
                    target,
                    len(ids),
                    ";".join(sorted(ids, key=_id_order)),
                    index.name_of(source),
                    index.name_of(target),
                    min(levels),
                    *(levels[level] for level in LEVELS),
                ]
            )
    with (folder / f"biotools_nodes{suffix}.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "biotoolsID",
                "biotools_name",
                "n_workflows",
                "has_publication",
                "confidence",  # the best over the workflows
                *(f"confidence{level}" for level in LEVELS),  # workflows per level
            ]
        )
        for biotools_id, ids in sorted(node_workflows.items(), key=lambda kv: (-len(kv[1]), kv[0])):
            levels = node_levels[biotools_id]
            writer.writerow(
                [
                    biotools_id,
                    index.name_of(biotools_id),
                    len(ids),
                    str(biotools_id in with_paper).lower(),
                    min(levels),
                    *(levels[level] for level in LEVELS),
                ]
            )
    edges_by = Counter(min(levels) for levels in edge_levels.values())
    nodes_by = Counter(min(node_levels[n]) for n in node_workflows)
    print(
        f"Wrote {folder}/biotools_*{suffix}.csv: {len(pairs)} workflows, {n_pairs} pairs, "
        f"{len(edge_workflows)} edges (confidence {_by_level(edges_by)}), "
        f"{len(node_workflows)} bio.tools nodes (confidence {_by_level(nodes_by)})"
    )
