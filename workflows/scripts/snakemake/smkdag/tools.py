"""The tools a rule runs, from (best first): its snakemake-wrapper, its conda environment, the commands of its
shell section restricted to the declared software, its container, and else the rule name itself."""

import ast
import re
from dataclasses import dataclass, field
from pathlib import Path

from nfdag.process import Step, choose_steps, environment_file, image_tool, is_excluded, same_tool
from nfdag.shell import commands

from .snakefile import Rule, Workflow, path_text
from .wrappers import WrapperRef, read_wrapper, wrapper_ref

# more declared packages than this, none of them run in the rule: an exported environment, not the rule's tools
MAX_DECLARED = 5
# conda packages of a wrapper that are not its tool
WRAPPER_HELPERS = frozenset({"snakemake-wrapper-utils", "snakemake", "snakemake-minimal"})


@dataclass(slots=True)
class RuleTools:
    """Everything known about the software of one rule."""

    wrapper: str = ""  # version/path
    wrapper_tool: str = ""  # bio/bwa/mem -> bwa
    homepage: str = ""  # the wrapper's meta.yaml url
    conda_env: str = ""  # the environment file, relative to the repository
    conda_packages: list[str] = field(default_factory=list)
    conda_tools: list[str] = field(default_factory=list)
    container: str = ""
    container_tool: str = ""
    executables: list[str] = field(default_factory=list)
    scripts: list[str] = field(default_factory=list)  # own scripts (script:, or called in shell:)
    kind: str = ""  # wrapper | shell | script | run | notebook | none
    steps: list[Step] = field(default_factory=list)


def _package_name(package: str) -> str:
    """bioconda::samtools=1.19 -> samtools"""
    return re.split(r"[=<>!~ ]", package.rsplit("::", maxsplit=1)[-1].strip())[0].lower()


def _string_text(text: str) -> str:
    """The text of a string expression (a shell command): literals and f-string parts joined."""
    try:
        node = ast.parse(f"(\n{text}\n)", mode="eval").body
    except SyntaxError:
        return text
    parts = [n.value for n in ast.walk(node) if isinstance(n, ast.Constant) and isinstance(n.value, str)]
    return "\n".join(parts)


def own_script_names(wf: Workflow) -> set[str]:
    """File names in scripts/ and bin/ folders: commands with these names are the workflow's own scripts."""
    names = set()
    for folder in ("scripts", "bin", "workflow/scripts", "workflow/bin"):
        path = wf.root / folder
        if path.is_dir():
            names |= {p.name for p in path.iterdir() if p.is_file()}
    return names


def _conda(rule: Rule, wf: Workflow, info: RuleTools) -> None:
    value = rule.sections.get("conda") or wf.global_conda
    if not value:
        return
    text = path_text(value)
    if text.endswith((".yaml", ".yml")):
        text = text.removeprefix("./")
        candidates = [rule.file.parent / text, wf.main.parent / text, wf.root / text, wf.root / "workflow" / text]
        env = next((c for c in candidates if c.is_file()), None)
        if env is None and Path(text).name:
            same = [p for p in wf.root.rglob(Path(text).name) if p.is_file()]
            env = min(same, key=lambda p: len(p.parts)) if same else None
        if env is not None:
            info.conda_env = env.resolve().relative_to(wf.root.resolve()).as_posix()
            info.conda_packages = environment_file(env)
    elif text and "/" not in text and " " not in text:
        info.conda_env = text  # a named, existing environment: nothing to read
    info.conda_tools = list(dict.fromkeys(_package_name(p) for p in info.conda_packages if not p.startswith("-")))


def _container(rule: Rule, wf: Workflow, info: RuleTools) -> None:
    value = rule.sections.get("container") or rule.sections.get("singularity") or wf.global_container
    image = path_text(value) if value else ""
    if image:
        info.container = image
        if re.search(r"biocontainers|depot\.galaxyproject\.org", image):
            info.container_tool = image_tool(image)[0]


def _wrapper(ref: WrapperRef, wrappers_dir: Path, excluded: frozenset[str], info: RuleTools) -> Step:
    """The step of a wrapper rule: the wrapper's tool, declared by its main conda package."""
    wrapper = read_wrapper(wrappers_dir, ref)
    info.wrapper, info.wrapper_tool, info.homepage = f"{ref.version}/{ref.path}", ref.tool, wrapper.homepage
    if wrapper.packages:
        info.conda_packages = wrapper.packages
        info.conda_tools = list(dict.fromkeys(_package_name(p) for p in wrapper.packages))
    packages = [t for t in info.conda_tools if t not in WRAPPER_HELPERS and not is_excluded(t, excluded)]
    tool = ref.tool
    parts = ref.path.split("/")
    if len(parts) > 2 and not any(same_tool(tool, p) for p in packages):
        # a category, not a tool: bio/benchmark/chm-eval -> chm-eval (if that is the package)
        tool = next((parts[2] for p in packages if same_tool(parts[2], p)), tool)
    info.wrapper_tool = tool
    main = next((p for p in packages if same_tool(tool, p)), packages[0] if len(packages) == 1 else "")
    return Step(tool, "wrapper", "", main)


def rule_tools(
    rule: Rule, wf: Workflow, wrappers_dir: Path, excluded: frozenset[str], bin_names: set[str]
) -> RuleTools:
    """The tools of one rule, as steps in the order they run."""
    info = RuleTools()
    _container(rule, wf, info)
    sections = rule.sections
    if "wrapper" in sections:
        info.kind = "wrapper"
        ref = wrapper_ref(path_text(sections["wrapper"]) or sections["wrapper"].strip("\"' "))
        if ref is not None:
            info.steps = [_wrapper(ref, wrappers_dir, excluded, info)]
            return info
    _conda(rule, wf, info)
    if "shell" in sections:
        info.kind = "shell"
        found = commands(_string_text(sections["shell"]), bin_names)
        info.executables, info.scripts = found.executables, found.scripts
    else:
        info.kind = next((k for k in ("script", "run", "notebook", "cwl") if k in sections), "none")
        if info.kind in {"script", "notebook"}:
            info.scripts = [Path(path_text(sections[info.kind])).name]
    sources = {"conda": info.conda_tools, "container": [info.container_tool] if info.container_tool else []}
    steps = choose_steps(info.executables, sources, excluded, rule.name)
    if len(steps) > MAX_DECLARED and all(s.source != "executable" for s in steps):
        # nothing in the rule matched, and the environment is a full `conda env export` (hundreds of pinned
        # libraries): the packages say nothing about the rule's tool, the rule name stands for it
        steps = [Step(rule.name, "process")]
    info.steps = steps
    return info
