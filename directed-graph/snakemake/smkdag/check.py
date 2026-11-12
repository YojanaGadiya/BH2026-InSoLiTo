"""Compare a static rule graph with the one Snakemake drew for the catalog (`snakemake --rulegraph` on the
workflow's test data, so it only holds the rules that test run needs)."""

import re
from dataclasses import dataclass
from pathlib import Path

_NODE = re.compile(r"^\s*(\d+)\s*\[\s*label\s*=\s*\"([^\"]+)\"", re.M)
_EDGE = re.compile(r"^\s*(\d+)\s*->\s*(\d+)", re.M)


def read_rulegraph(path: Path) -> tuple[set[str], set[tuple[str, str]]]:
    """(rules, edges) of a DOT rule graph."""
    text = path.read_text(encoding="utf-8")
    names = {m.group(1): m.group(2) for m in _NODE.finditer(text)}
    edges = {(names[a], names[b]) for a, b in _EDGE.findall(text) if a in names and b in names}
    return set(names.values()), edges


@dataclass(slots=True)
class Check:
    catalog_edges: int  # edges between rules that are DAG nodes on both sides
    found: int  # of those, found statically
    extra: int  # static edges between rules of the catalog graph that it does not have
    missing: list[tuple[str, str]]

    @property
    def recall(self) -> float:
        return self.found / self.catalog_edges if self.catalog_edges else 1.0

    @property
    def precision(self) -> float:
        return self.found / (self.found + self.extra) if self.found + self.extra else 1.0


def compare(path: Path, rules: set[str], edges: set[tuple[str, str]]) -> Check:
    """Only rules both graphs have count: the catalog's graph lacks rules its test run does not need,
    the static graph lacks target rules (`all`)."""
    catalog_rules, catalog_edges = read_rulegraph(path)
    common = catalog_rules & rules
    expected = {(a, b) for a, b in catalog_edges if a in common and b in common}
    ours = {(a, b) for a, b in edges if a in common and b in common}
    return Check(len(expected), len(expected & ours), len(ours - expected), sorted(expected - ours))
