"""The file paths of a rule section, read with Python's ast (nothing is evaluated).

A path is a pattern with Snakemake wildcards (`results/{sample}.bam`); parts that are only known at run time
(config values, function arguments) become the wildcard `{*}`. A reference `rules.X.output` is kept as such.
Input functions (`input: get_bams`, `lambda wc: ...`) are mined for the paths they can return.
"""

import ast
import re
import textwrap
from dataclasses import dataclass, field
from typing import Final

from .snakefile import Workflow

UNKNOWN: Final = "{*}"
MAX_DEPTH: Final = 4  # variables / functions followed this deep
MAX_PATTERNS: Final = 64  # per expression, against combinatorial explosion of `a + b` with many alternatives
# Snakemake functions that wrap a path (temp("x") is the path x)
_PASS_THROUGH_WORDS: Final = """
    temp temporary protected ancient directory touch pipe report unpack ensure local service storage dynamic
    Path PurePath PurePosixPath str abspath normpath realpath
"""
PASS_THROUGH: Final = frozenset(_PASS_THROUGH_WORDS.split())
JOIN: Final = frozenset({"join"})  # os.path.join, path.join, join
EXPAND: Final = frozenset({"expand", "collect"})
# a string that looks like a file path: a directory part or a file extension
_PATHLIKE: Final = re.compile(r"/|\.[A-Za-z][\w]{0,7}(?:\.gz|\.bz2)?$|\{[\w,*]+\}\.")


@dataclass(frozen=True, slots=True)
class PathRef:
    """One input or output of a rule."""

    pattern: str = ""  # path with wildcards
    rule: str = ""  # rules.X.output -> X
    via: str = "literal"  # literal | variable | function | rules


@dataclass(slots=True)
class Section:
    """The paths of a section, and how much of it could not be read."""

    refs: list[PathRef] = field(default_factory=list)
    items: int = 0  # top-level items (inputs) in the section
    unresolved: int = 0  # items that gave no usable path


class PathReader:
    """Reads paths from section texts of one workflow, with its variables and functions."""

    def __init__(self, wf: Workflow) -> None:
        self.wf = wf
        self._functions: dict[str, ast.FunctionDef | None] = {}

    def section(self, text: str, mine_functions: bool = True) -> Section:
        """The paths of an `input:` / `output:` section (its text holds call arguments: `a, b=..., c`)."""
        result = Section()
        if not text.strip():
            return result
        try:
            call = ast.parse(f"f(\n{text}\n)", mode="eval").body
        except SyntaxError:
            literals = re.findall(r"""["']([^"'\n]+)["']""", text)
            result.items = max(len(literals), 1)
            result.refs = [PathRef(lit) for lit in literals if _PATHLIKE.search(lit)]
            result.unresolved = result.items - len(result.refs)
            return result
        assert isinstance(call, ast.Call)
        items = [*call.args, *(k.value for k in call.keywords)]
        result.items = len(items)
        for item in items:
            refs = self.refs(item, "literal", 0, mine_functions)
            usable = [r for r in refs if r.rule or useful(r.pattern)]
            result.refs += usable
            result.unresolved += not usable
        return result

    def refs(self, node: ast.expr, via: str, depth: int, mine: bool = True) -> list[PathRef]:
        """The paths an expression can stand for."""
        if depth > MAX_DEPTH:
            return []
        if (rule := _rule_reference(node)) is not None:
            return [PathRef(rule=rule, via="rules")]
        patterns = self._patterns(node, via, depth, mine)
        return patterns[:MAX_PATTERNS]

    def _patterns(self, node: ast.expr, via: str, depth: int, mine: bool) -> list[PathRef]:  # noqa: PLR0911
        if isinstance(node, ast.Constant):
            return [PathRef(node.value, via=via)] if isinstance(node.value, str) else []
        if isinstance(node, ast.JoinedStr):
            return [PathRef(self._fstring(node, depth), via=via)]
        if isinstance(node, (ast.List, ast.Tuple, ast.Set)):
            return [r for e in node.elts for r in self.refs(e, via, depth, mine)]
        if isinstance(node, ast.Dict):
            return [r for v in node.values if v is not None for r in self.refs(v, via, depth, mine)]
        if isinstance(node, ast.Starred):
            return self.refs(node.value, via, depth, mine)
        if isinstance(node, ast.IfExp):
            return self.refs(node.body, via, depth, mine) + self.refs(node.orelse, via, depth, mine)
        if isinstance(node, (ast.ListComp, ast.GeneratorExp, ast.SetComp)):
            return self.refs(node.elt, via, depth, mine)
        if isinstance(node, ast.BinOp):
            return self._binop(node, via, depth, mine)
        if isinstance(node, ast.Call):
            return self._call(node, via, depth, mine)
        if isinstance(node, ast.Name):
            return self._name(node.id, via, depth, mine)
        if isinstance(node, ast.Lambda):
            return [replace_via(r, "function") for r in self._mine(node.body, depth)] if mine else []
        if isinstance(node, (ast.Subscript, ast.Attribute)):
            return [PathRef(UNKNOWN, via=via)]  # config["x"], params.x: known at run time
        return []

    def _fstring(self, node: ast.JoinedStr, depth: int) -> str:
        parts = []
        for value in node.values:
            if isinstance(value, ast.Constant):
                parts.append(str(value.value))
            elif isinstance(value, ast.FormattedValue):
                parts.append(self._placeholder(value.value, depth))
        return "".join(parts)

    def _placeholder(self, node: ast.expr, depth: int) -> str:
        """The text an f-string placeholder stands for: wildcards.sample -> {sample}, a string variable -> its
        value, anything else -> {*}."""
        if (
            isinstance(node, ast.Attribute)
            and isinstance(node.value, ast.Name)
            and node.value.id in {"wildcards", "wc"}
        ):
            return "{" + node.attr + "}"
        if isinstance(node, ast.Name):
            values = {r.pattern for r in self._name(node.id, "variable", depth + 1, False) if not r.rule}
            if len(values) == 1:
                return values.pop()
        return UNKNOWN

    def _binop(self, node: ast.BinOp, via: str, depth: int, mine: bool) -> list[PathRef]:
        if isinstance(node.op, ast.Mod):  # "x/%s.bam" % sample
            return [PathRef(re.sub(r"%\w", UNKNOWN, r.pattern), via=r.via) for r in self.refs(node.left, via, depth)]
        if not isinstance(node.op, (ast.Add, ast.Div)):
            return []
        sep = "/" if isinstance(node.op, ast.Div) else ""
        left = self.refs(node.left, via, depth, mine) or [PathRef(UNKNOWN, via=via)]
        right = self.refs(node.right, via, depth, mine) or [PathRef(UNKNOWN, via=via)]
        out = []
        for a in left:
            for b in right:
                if not a.rule and not b.rule:
                    out.append(PathRef(f"{a.pattern}{sep}{b.pattern}", via=a.via if a.via != "literal" else b.via))
        return out[:MAX_PATTERNS]

    def _call(self, node: ast.Call, via: str, depth: int, mine: bool) -> list[PathRef]:
        name = _call_name(node.func)
        args = node.args
        if name in PASS_THROUGH or name in EXPAND:
            if args:
                return self.refs(args[0], via, depth, mine)
            return [r for k in node.keywords for r in self.refs(k.value, via, depth, mine)]
        if name == "multiext" and args:
            exts = [a.value for a in args[1:] if isinstance(a, ast.Constant) and isinstance(a.value, str)]
            exts += [
                k.value.value
                for k in node.keywords
                if isinstance(k.value, ast.Constant) and isinstance(k.value.value, str)
            ]
            return [PathRef(f"{p.pattern}{ext}", via=p.via) for p in self.refs(args[0], via, depth) for ext in exts]
        if name in JOIN and args:
            joined = [PathRef("", via=via)]
            for arg in args:
                parts = self.refs(arg, via, depth, mine) or [PathRef(UNKNOWN, via=via)]
                joined = [
                    PathRef(f"{a.pattern}/{b.pattern}" if a.pattern else b.pattern, via=via)
                    for a in joined
                    for b in parts
                    if not b.rule
                ][:MAX_PATTERNS]
            return joined
        if name == "format" and isinstance(node.func, ast.Attribute):  # "x/{}.bam".format(...)
            return [
                PathRef(re.sub(r"\{\d*\}", UNKNOWN, r.pattern), via=r.via)
                for r in self.refs(node.func.value, via, depth)
            ]
        if name in {"get", "output"} and isinstance(node.func, ast.Attribute):
            return self.refs(node.func.value, via, depth, mine)
        if name and mine and self._function(name) is not None:
            return self._call_function(name, depth)
        return []

    def _name(self, name: str, via: str, depth: int, mine: bool) -> list[PathRef]:
        if mine and self._function(name) is not None:  # input: get_bams
            return self._call_function(name, depth)
        exprs = self.wf.variables.get(name)
        if not exprs:
            return []
        out = []
        for expr in exprs[-3:]:  # the last assignments (if/else branches)
            try:
                node = ast.parse(expr.strip(), mode="eval").body
            except SyntaxError:
                continue
            out += [
                replace_via(r, "variable" if via == "literal" else via) for r in self.refs(node, via, depth + 1, mine)
            ]
        return out

    def _function(self, name: str) -> ast.FunctionDef | None:
        if name not in self._functions:
            source = self.wf.functions.get(name)
            parsed = None
            if source:
                try:
                    tree = ast.parse(textwrap.dedent(source))
                    parsed = next((n for n in tree.body if isinstance(n, ast.FunctionDef)), None)
                except SyntaxError:
                    parsed = None
            self._functions[name] = parsed
        return self._functions[name]

    def _call_function(self, name: str, depth: int) -> list[PathRef]:
        function = self._function(name)
        if function is None or depth > MAX_DEPTH:  # also stops functions calling each other in a circle
            return []
        return [replace_via(r, "function") for stmt in function.body for r in self._mine(stmt, depth + 1)]

    def _mine(self, node: ast.AST, depth: int) -> list[PathRef]:
        """Paths an input function can return: every path-like expression in its body, rules.X references,
        and what the functions it calls return. Sub-expressions of a path already taken are not taken again."""
        out: list[PathRef] = []
        taken: set[int] = set()
        for sub in ast.walk(node):
            if id(sub) in taken or not isinstance(sub, ast.expr):
                continue
            refs: list[PathRef] = []
            if _rule_reference(sub) is not None or _is_path_expr(sub):
                refs = self.refs(sub, "function", depth, mine=True)
            elif isinstance(sub, ast.Call) and (name := _call_name(sub.func)) and self._function(name) is not None:
                refs = self._call_function(name, depth)
            if refs:
                out += [r for r in refs if r.rule or useful(r.pattern)]
                taken.update(id(n) for n in ast.walk(sub))
        return out[:MAX_PATTERNS]


def replace_via(ref: PathRef, via: str) -> PathRef:
    return ref if ref.rule else PathRef(ref.pattern, ref.rule, via)


def _call_name(func: ast.expr) -> str:
    if isinstance(func, ast.Name):
        return func.id
    if isinstance(func, ast.Attribute):
        return func.attr
    return ""


def _is_path_expr(node: ast.expr) -> bool:
    """Whether an expression builds a path: a path-like string, f-string, expand(), multiext(), os.path.join()."""
    if isinstance(node, ast.Constant):
        return isinstance(node.value, str) and bool(_PATHLIKE.search(node.value)) and " " not in node.value
    if isinstance(node, ast.JoinedStr):
        text = "".join(v.value for v in node.values if isinstance(v, ast.Constant) and isinstance(v.value, str))
        return bool(_PATHLIKE.search(text)) and " " not in text
    if isinstance(node, ast.Call):
        return _call_name(node.func) in EXPAND | JOIN | {"multiext"} | PASS_THROUGH
    return False


def _rule_reference(node: ast.expr) -> str | None:
    """rules.X.output[...] / rules.X.output.name / checkpoints.X.get(...).output -> X"""
    while True:
        if isinstance(node, ast.Subscript):
            node = node.value
        elif isinstance(node, ast.Call):
            node = node.func
        elif isinstance(node, ast.Attribute):
            if isinstance(node.value, ast.Name) and node.value.id in {"rules", "checkpoints"}:
                return node.attr
            node = node.value
        else:
            return None


def useful(pattern: str) -> bool:
    """Whether a path pattern says anything: some literal text beyond wildcards, slashes and dots."""
    literal = re.sub(r"\{[^{}]*\}", "", pattern)
    return bool(re.search(r"[A-Za-z0-9]", literal)) and " " not in pattern.strip()
