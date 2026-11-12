"""Build the rule DAG of one workflow, and split each rule into the tools it runs.

Rules are linked like Snakemake does it, by file names: an input path of rule B that matches an output pattern
of rule A gives the edge A -> B. If several rules could produce an input, the most specific output pattern
wins (the longest literal text), as Snakemake picks one producer. Wildcards an input fills in at run time
(expand(..., tool=config[...])) can stand for several producers, which all count, like every if/else branch
of a Nextflow pipeline. `rules.A.output` gives the edge directly.
Rules without output that run no known tool (`all`, `help`) are targets, not DAG nodes.
"""

import itertools
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Final

from nfdag.graph import Graph
from nfdag.process import Step

from .paths import UNKNOWN, PathReader, Section
from .snakefile import Rule, Workflow, main_snakefile, read_workflow
from .tools import RuleTools, own_script_names, rule_tools

_WILDCARD: Final = re.compile(r"\{([^{}]*)\}")
PROBE: Final = "w"  # stands for a wildcard value of an input when it is matched against output patterns
RUNS: Final = ("shell", "script", "run", "wrapper", "notebook", "cwl")


class NoSnakefileError(Exception):
    """The repository holds no Snakefile (at the release that was cloned)."""


@dataclass(slots=True)
class RuleNode:
    rule: Rule
    inputs: Section
    outputs: Section
    tools: RuleTools


@dataclass(slots=True)
class Result:
    """The DAG of one workflow."""

    wf: Workflow
    rules: dict[str, RuleNode]  # the rules that are DAG nodes
    targets: list[str]  # rules without output that run nothing (all)
    rule_edges: dict[tuple[str, str], set[str]] = field(default_factory=dict)  # (A, B) -> via
    nodes: dict[str, tuple[str, int, Step]] = field(default_factory=dict)  # tool node -> (rule, index, step)
    edges: dict[tuple[str, str], tuple[str, str]] = field(default_factory=dict)  # -> (kind, via)
    warnings: list[str] = field(default_factory=list)

    @property
    def inputs(self) -> int:
        return sum(n.inputs.items for n in self.rules.values())

    @property
    def unresolved(self) -> int:
        return sum(n.inputs.unresolved for n in self.rules.values())


def _normalise(pattern: str) -> str:
    pattern = pattern.replace("{{", "{").replace("}}", "}")
    return re.sub(r"/{2,}", "/", pattern.strip()).removeprefix("./")


def _regex(pattern: str) -> re.Pattern[str]:
    """An output pattern as a regex: {sample} -> .+, {x,[^/]+} -> [^/]+."""
    out, last = [], 0
    for m in _WILDCARD.finditer(pattern):
        out.append(re.escape(pattern[last : m.start()]))
        name, _, constraint = m.group(1).partition(",")
        out.append(f"(?:{constraint})" if constraint and name != "*" else ".+")
        last = m.end()
    out.append(re.escape(pattern[last:]))
    try:
        return re.compile("".join(out))
    except re.error:
        return re.compile(_WILDCARD.sub(".+", re.escape(pattern).replace(r"\{", "{").replace(r"\}", "}")))


def _literal_length(pattern: str) -> int:
    return len(_WILDCARD.sub("", pattern))


@dataclass(slots=True)
class _Output:
    rule: str
    pattern: str
    regex: re.Pattern[str]
    specificity: int


def _producers(path: str, outputs: list[_Output], consumer: str) -> list[str]:
    """The rules whose output pattern matches an input path; only the most specific ones.
    First the input against the output patterns; if none matches, the other way round, because an input's
    wildcards may stand for values filled in at run time (expand("results/{tool}/x.bam", tool=config[...])
    is produced by a rule writing results/bwa/x.bam)."""
    probe = _WILDCARD.sub(PROBE, path)
    matched = [o for o in outputs if o.rule != consumer and o.regex.fullmatch(probe)]
    if matched:
        best = max(o.specificity for o in matched)
        return list(dict.fromkeys(o.rule for o in matched if o.specificity == best))
    if not (_WILDCARD.search(path) or UNKNOWN in path):
        return []
    # every rule that can write a file of this pattern: the value is a run-time choice (all branches count)
    input_regex = _regex(path)
    reverse = [o for o in outputs if o.rule != consumer and input_regex.fullmatch(_WILDCARD.sub(PROBE, o.pattern))]
    return list(dict.fromkeys(o.rule for o in reverse))


def build(repo: Path, wrappers_dir: Path, excluded: frozenset[str]) -> Result:
    main = main_snakefile(repo)
    if main is None:
        raise NoSnakefileError(repo.name)
    wf = read_workflow(repo, main)
    reader = PathReader(wf)
    bin_names = own_script_names(wf)
    rules: dict[str, RuleNode] = {}
    targets = []
    for name, rule in wf.rules.items():
        outputs = reader.section(rule.sections.get("output", ""), mine_functions=False)
        if not outputs.refs and not any(k in rule.sections for k in RUNS):
            targets.append(name)
            continue
        tools = rule_tools(rule, wf, wrappers_dir, excluded, bin_names)
        if not outputs.refs and all(step.source == "process" for step in tools.steps):
            targets.append(name)  # `all` / `help` / `aggregate` printing a message: no output, no tool
            continue
        inputs = reader.section(rule.sections.get("input", ""))
        rules[name] = RuleNode(rule, inputs, outputs, tools)
    result = Result(wf, rules, targets, warnings=list(wf.warnings))
    _link(result)
    _expand(result)
    return result


def _link(result: Result) -> None:
    outputs = [
        _Output(name, p, _regex(p), _literal_length(p))
        for name, node in result.rules.items()
        for p in dict.fromkeys(_normalise(r.pattern) for r in node.outputs.refs if r.pattern)
    ]
    graph = Graph()
    for name, node in result.rules.items():
        producers: dict[str, str] = {}  # producer -> via
        for ref in node.inputs.refs:
            if ref.rule:
                if ref.rule in result.rules:
                    producers.setdefault(ref.rule, "rules")
                elif ref.rule not in result.targets:
                    result.warnings.append(f"rule {name}: rules.{ref.rule} not found")
                continue
            via = "function" if ref.via == "function" else "path"
            for producer in _producers(_normalise(ref.pattern), outputs, name):
                if producers.get(producer) != "rules":
                    producers[producer] = "path" if producers.get(producer) == "path" else via
        for producer in graph.connect(producers.items(), name):
            result.warnings.append(f"left out {producer} -> {name}: it would close a cycle")
        if node.inputs.unresolved:
            result.warnings.append(f"rule {name}: {node.inputs.unresolved} of {node.inputs.items} inputs unresolved")
    result.rule_edges = graph.edges


def _expand(result: Result) -> None:
    """Each rule becomes the chain of its tools; an edge A -> B leaves A's last tool and reaches B's first."""
    chains: dict[str, list[str]] = {}
    for name, node in result.rules.items():
        chain = []
        for index, step in enumerate(node.tools.steps, 1):
            tool_node = f"{name}/{step.tool}"
            if tool_node not in result.nodes:
                result.nodes[tool_node] = (name, index, step)
                chain.append(tool_node)
        chains[name] = chain
        for a, b in itertools.pairwise(chain):
            result.edges[(a, b)] = ("chain", "")
    for (a, b), vias in result.rule_edges.items():
        result.edges[(chains[a][-1], chains[b][0])] = ("dag", ";".join(sorted(vias)))
