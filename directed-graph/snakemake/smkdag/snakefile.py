"""Read a Snakemake workflow: the main Snakefile, the files it includes, local modules, and per rule its sections.

Snakefiles are Python with extra keywords. They are split into logical lines (strings and brackets kept
together), and grouped by indentation: `rule NAME:` / `checkpoint NAME:` blocks with their sections
(`input:`, `output:`, `shell:`, ...), `def` functions, module-level variables, `include:`, `module` and
`use rule`. Nothing is executed.
"""

import ast
import re
from collections.abc import Iterator
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Final

MAIN_CANDIDATES: Final = ("workflow/Snakefile", "Snakefile", "workflow/snakefile", "snakefile", "workflow/main.smk")
SNAKEFILE_SUFFIXES: Final = (".smk", ".snakefile", ".snake", ".rules", ".rule")
NOT_SNAKEFILES: Final = (".tmpl", ".old", ".bak", ".orig", ".backup", ".py", ".yaml", ".yml", ".md", ".txt", ".sh")
SKIP_DIRS: Final = frozenset({".git", ".test", ".tests", "test", "tests", ".snakemake", "node_modules"})
_RULE_SECTION_WORDS: Final = """
    input output params log benchmark threads resources priority message conda container singularity
    envmodules wrapper shell script run notebook cwl shadow group localrule retries default_target cache
    handover wildcard_constraints name template_engine
"""
RULE_SECTIONS: Final = frozenset(_RULE_SECTION_WORDS.split())
_RULE: Final = re.compile(r"(rule|checkpoint)(?:\s+(\w+))?\s*:\s*(.*)$", re.S)
_USE: Final = re.compile(
    r"use\s+rule\s+(?P<rules>\*|[\w\s,]+?)(?:\s+from\s+(?P<module>\w+))?(?:\s+as\s+(?P<alias>[\w*]+))?"
    r"(?:\s+(?P<with>with))?\s*:?\s*$",
    re.S,
)
_SECTION: Final = re.compile(r"(\w+)\s*:(?!=)\s*(.*)$", re.S)
_ASSIGN: Final = re.compile(r"([A-Za-z_]\w*)\s*(?::[^=]+)?=(?!=)\s*(.+)$", re.S)
_DEF: Final = re.compile(r"(?:async\s+)?def\s+(\w+)\s*\(")
_KEYWORD_LINE: Final = re.compile(r"(include|module|configfile|workdir|conda|container|singularity)\s*:\s*(.*)$", re.S)


@dataclass(slots=True)
class Line:
    """A logical line: a statement with its continuation lines, brackets and strings kept together."""

    indent: int
    text: str  # without the leading indentation, comments outside strings removed
    number: int  # first physical line, from 1


def logical_lines(source: str) -> list[Line]:  # noqa: PLR0915  (one character scanner)
    """Split Python-like source into logical lines, aware of strings, brackets and backslash continuations."""
    lines: list[Line] = []
    buffer: list[str] = []
    depth = 0
    quote = ""  # the open string delimiter (' " ''' or \"\"\")
    start = number = 1
    i, n = 0, len(source)
    while i < n:
        c = source[i]
        if quote:
            if c == "\\":
                buffer.append(source[i : i + 2])
                number += source[i : i + 2].count("\n")
                i += 2
                continue
            if source.startswith(quote, i):
                buffer.append(quote)
                i += len(quote)
                quote = ""
                continue
            if c == "\n":
                number += 1
                if len(quote) == 1:  # unterminated one-line string: give up on it
                    quote = ""
            buffer.append(c)
            i += 1
            continue
        if c == "#":
            while i < n and source[i] != "\n":
                i += 1
            continue
        if c in "\"'":
            quote = source[i : i + 3] if source[i : i + 3] in ('"""', "'''") else c
            buffer.append(quote)
            i += len(quote)
            continue
        if c in "([{":
            depth += 1
        elif c in ")]}":
            depth = max(depth - 1, 0)
        elif c == "\\" and source[i + 1 : i + 2] == "\n":
            buffer.append(" ")
            i += 2
            number += 1
            continue
        if c == "\n":
            number += 1
            if depth == 0:
                _flush(buffer, start, lines)
                buffer, start = [], number
                i += 1
                continue
        buffer.append(c)
        i += 1
    _flush(buffer, start, lines)
    return lines


def _flush(buffer: list[str], start: int, lines: list[Line]) -> None:
    text = "".join(buffer)
    if text.strip():
        stripped = text.lstrip(" \t")
        indent = len(text[: len(text) - len(stripped)].expandtabs(4))
        lines.append(Line(indent, stripped.rstrip(), start))


def _block(lines: list[Line], i: int) -> int:
    """The index after the block of lines indented deeper than lines[i]."""
    j = i + 1
    while j < len(lines) and lines[j].indent > lines[i].indent:
        j += 1
    return j


@dataclass(slots=True)
class Rule:
    """One rule (or checkpoint)."""

    name: str
    file: Path
    line: int
    sections: dict[str, str] = field(default_factory=dict)  # section -> its text (Python expressions)
    checkpoint: bool = False
    module: str = ""  # the module the rule comes from (`use rule * from M`)


@dataclass(slots=True)
class UseRule:
    """A `use rule` statement, applied when everything is read."""

    rules: list[str]  # names, or ["*"]
    module: str  # "" for a rule of this workflow
    alias: str  # new name, `*` replaced by the rule name; "" keeps the name
    sections: dict[str, str]  # overriding sections (`with:`)
    file: Path
    line: int


@dataclass(slots=True)
class Workflow:
    """Everything read from one workflow (or a module of it)."""

    root: Path  # the repository
    main: Path  # the main Snakefile
    files: list[Path] = field(default_factory=list)
    rules: dict[str, Rule] = field(default_factory=dict)
    functions: dict[str, str] = field(default_factory=dict)  # name -> source of the def
    variables: dict[str, list[str]] = field(default_factory=dict)  # name -> assigned expressions
    global_conda: str = ""  # a top-level `conda:` (expression text)
    global_container: str = ""
    modules: dict[str, dict[str, str]] = field(default_factory=dict)  # module name -> its sections
    uses: list[UseRule] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


def _named_snakefile(path: Path) -> bool:
    """Snakefile, also Snakefile.experiment, Snakefile-bin (not templates or backups: .tmpl, .old, .bak)."""
    name = path.name.lower()
    return name == "snakefile" or (
        name.startswith(("snakefile.", "snakefile-", "snakefile_")) and not name.endswith(NOT_SNAKEFILES)
    )


def _is_snakefile(path: Path) -> bool:
    return _named_snakefile(path) or path.suffix in SNAKEFILE_SUFFIXES


def workflow_files(repo: Path) -> list[Path]:
    """All Snakefiles and rule files of a repository (not in tests or hidden folders)."""
    found = []
    for path in repo.rglob("*"):
        rel = path.relative_to(repo).parts
        if path.is_file() and _is_snakefile(path) and not any(p in SKIP_DIRS or p.startswith(".") for p in rel[:-1]):
            found.append(path)
    return sorted(found)


def main_snakefile(repo: Path) -> Path | None:
    """The main Snakefile: workflow/Snakefile, else Snakefile, else the least deep file named Snakefile."""
    for candidate in MAIN_CANDIDATES:
        if (repo / candidate).is_file():
            return repo / candidate
    files = workflow_files(repo)
    for named in ([p for p in files if p.name.lower() == "snakefile"], [p for p in files if _named_snakefile(p)]):
        if named:
            return min(named, key=lambda p: (len(p.relative_to(repo).parts), len(p.name), str(p)))
    smk = [p for p in files if p.suffix in SNAKEFILE_SUFFIXES]
    return min(smk, key=lambda p: (len(p.relative_to(repo).parts), str(p))) if smk else None


def string_value(text: str) -> str | None:
    """The value of a string expression, if it is one literal (also an f-string without placeholders)."""
    try:
        node = ast.parse(text.strip(), mode="eval").body
    except SyntaxError:
        return None
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    return None


class Reader:
    """Reads the Snakefiles of one workflow, following includes."""

    def __init__(self, repo: Path, main: Path) -> None:
        self.wf = Workflow(repo, main)
        self._anonymous = 0

    def read(self) -> Workflow:
        self.read_file(self.wf.main, "")
        self._apply_uses()
        return self.wf

    def read_file(self, path: Path, module: str) -> None:
        if path in self.wf.files:
            return
        self.wf.files.append(path)
        try:
            source = path.read_text(encoding="utf-8", errors="replace")
        except OSError as error:
            self.wf.warnings.append(f"{self._rel(path)}: cannot read: {error}")
            return
        lines = logical_lines(source)
        i = 0
        while i < len(lines):
            i = self._statement(lines, i, path, module)

    def _rel(self, path: Path) -> str:
        try:
            return path.resolve().relative_to(self.wf.root.resolve()).as_posix()
        except ValueError:
            return str(path)

    def _statement(self, lines: list[Line], i: int, path: Path, module: str) -> int:
        """Read the statement at lines[i]; returns the index of the next one to read."""
        line = lines[i]
        end = _block(lines, i)
        if m := _RULE.match(line.text):
            self._rule(line, lines[i + 1 : end], m, path, module)
            return end
        if line.text.startswith("use ") and (m := _USE.match(line.text)):
            self._use(lines, i, end, m, path)
            return end
        if m := re.match(r"module\s+(\w+)\s*:\s*$", line.text):
            self.wf.modules[m.group(1)] = self._sections(lines[i + 1 : end])
            return end
        if m := _DEF.match(line.text):
            self.wf.functions[m.group(1)] = "\n".join(" " * (ln.indent - line.indent) + ln.text for ln in lines[i:end])
            return end
        if m := _KEYWORD_LINE.match(line.text):
            keyword, value = m.group(1), m.group(2)
            if not value.strip() and end > i + 1:  # value on the next line(s)
                value = "\n".join(ln.text for ln in lines[i + 1 : end])
            self._keyword(keyword, value, path, module, line.number)
            return end
        if m := _ASSIGN.match(line.text):
            self.wf.variables.setdefault(m.group(1), []).append(m.group(2))
        # if/for/with/try blocks: read their body (rules and includes inside `if config[...]:` count)
        return i + 1

    def _keyword(self, keyword: str, value: str, path: Path, module: str, number: int) -> None:
        if keyword == "include":
            target = self._resolve_path(value, path)
            if target is None:
                self.wf.warnings.append(f"{self._rel(path)}:{number}: include not found: {value.strip()[:80]}")
            elif _is_snakefile(target) or target.suffix in {".py", ""}:
                self.read_file(target, module)
        elif keyword == "conda":
            self.wf.global_conda = self.wf.global_conda or value
        elif keyword in {"container", "singularity"}:
            self.wf.global_container = self.wf.global_container or value

    def _sections(self, body: list[Line]) -> dict[str, str]:
        """The `keyword: value` sections of a rule body (or module / use-rule body)."""
        sections: dict[str, str] = {}
        if not body:
            return sections
        base = body[0].indent
        k = 0
        while k < len(body):
            line = body[k]
            m = _SECTION.match(line.text) if line.indent == base else None
            j = k + 1
            while j < len(body) and body[j].indent > base:
                j += 1
            if m and m.group(1) in RULE_SECTIONS | {
                "snakefile",
                "config",
                "prefix",
                "replace_prefix",
                "skip_validation",
            }:
                parts = [m.group(2)] if m.group(2).strip() else []
                if m.group(1) == "run":  # Python code: keep the indentation
                    parts += [" " * (ln.indent - base) + ln.text for ln in body[k + 1 : j]]
                else:
                    parts += [ln.text for ln in body[k + 1 : j]]
                sections[m.group(1)] = "\n".join(parts)
            k = j
        return sections

    def _rule(self, header: Line, body: list[Line], m: re.Match[str], path: Path, module: str) -> None:
        name = m.group(2)
        sections = self._sections(body)
        if m.group(3).strip():  # `rule x: input: ...` on one line is rare; keep what is there
            inline = _SECTION.match(m.group(3).strip())
            if inline:
                sections.setdefault(inline.group(1), inline.group(2))
        if not name:
            named = string_value(sections.get("name", ""))
            if named and re.fullmatch(r"\w+", named):
                name = named
            else:
                self._anonymous += 1
                name = f"_rule{self._anonymous}"
        if module:
            name = f"{module}:{name}"  # renamed by `use rule * from module as ...`
        if name in self.wf.rules:
            name = f"{name}_{header.number}"  # the same rule in two branches (if/else)
        self.wf.rules[name] = Rule(name, path, header.number, sections, m.group(1) == "checkpoint", module)

    def _use(self, lines: list[Line], i: int, end: int, m: re.Match[str], path: Path) -> None:
        names = [r.strip() for r in m.group("rules").split(",") if r.strip()]
        sections = self._sections(lines[i + 1 : end]) if m.group("with") else {}
        self.wf.uses.append(
            UseRule(names, m.group("module") or "", m.group("alias") or "", sections, path, lines[i].number)
        )

    def _resolve_path(self, value: str, base: Path) -> Path | None:
        """A file named by an include / snakefile / conda expression, relative to the including file."""
        text = path_text(value)
        if not text:
            return None
        text = text.removeprefix("./")
        candidates = [base.parent / text, self.wf.root / text, self.wf.root / "workflow" / text]
        found = next((c for c in candidates if c.is_file()), None)
        if found is None:  # e.g. a variable prefix: look for the file name in the repository
            name = Path(text).name
            same = [p for p in self.wf.root.rglob(name) if p.is_file()] if name and "{" not in name else []
            found = min(same, key=lambda p: len(p.parts)) if same else None
        return found.resolve() if found else None

    def _apply_uses(self) -> None:
        """`use rule a as b with:` copies (and changes) a rule; `use rule * from m` takes a module's rules."""
        for use in self.wf.uses:
            if use.module:
                self._read_module(use)
            for name in self._used_names(use):
                rule = self.wf.rules.get(f"{use.module}:{name}" if use.module else name)
                if rule is None:
                    self.wf.warnings.append(f"{self._rel(use.file)}:{use.line}: use rule {name}: rule not found")
                    continue
                alias = use.alias.replace("*", name) if use.alias else name
                sections = {**rule.sections, **use.sections}
                self.wf.rules[alias] = replace(rule, name=alias, sections=sections, module="")
        for name in [n for n, r in self.wf.rules.items() if r.module]:  # module rules not used
            del self.wf.rules[name]

    def _used_names(self, use: UseRule) -> Iterator[str]:
        if use.rules != ["*"]:
            yield from use.rules
            return
        prefix = f"{use.module}:"
        yield from [n.removeprefix(prefix) for n, r in self.wf.rules.items() if r.module == use.module]

    def _read_module(self, use: UseRule) -> None:
        module = self.wf.modules.get(use.module)
        if module is None:
            self.wf.warnings.append(f"{self._rel(use.file)}:{use.line}: module {use.module} not defined")
            return
        if any(r.module == use.module for r in self.wf.rules.values()):
            return  # read already
        snakefile = module.get("snakefile", "")
        if re.search(r"\b(github|gitlab|https?)\b", snakefile):
            self.wf.warnings.append(f"module {use.module}: remote Snakefile, not read: {snakefile.strip()[:100]}")
            return
        target = self._resolve_path(snakefile, use.file)
        if target is None:
            self.wf.warnings.append(f"module {use.module}: Snakefile not found: {snakefile.strip()[:100]}")
            return
        self.read_file(target, use.module)


def path_text(value: str) -> str:
    """The path an include / conda / script value names, for simple expressions: a string, an f-string or
    os.path.join / workflow.source_path / Path() of strings and variables (variables are dropped)."""
    try:
        node = ast.parse(value.strip(), mode="eval").body
    except SyntaxError:
        m = re.search(r"""["']([^"']+)["']""", value)
        return m.group(1) if m else ""
    return _path_node(node)


def _path_node(node: ast.expr) -> str:
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    if isinstance(node, ast.JoinedStr):
        parts = [v.value for v in node.values if isinstance(v, ast.Constant) and isinstance(v.value, str)]
        return "".join(parts).lstrip("/")
    if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Add, ast.Div)):
        left, right = _path_node(node.left), _path_node(node.right)
        sep = "/" if isinstance(node.op, ast.Div) and left and right else ""
        return f"{left}{sep}{right}".lstrip("/")
    if isinstance(node, ast.Call):
        parts = [p for a in node.args if (p := _path_node(a))]
        return "/".join(p.strip("/") for p in parts)
    return ""


def read_workflow(repo: Path, main: Path) -> Workflow:
    return Reader(repo, main).read()
