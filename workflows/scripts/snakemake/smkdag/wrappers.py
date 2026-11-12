"""snakemake-wrappers: references in `wrapper:` directives, and their cached meta.yaml / environment.yaml."""

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Final

import yaml

RAW_URL: Final = "https://raw.githubusercontent.com/snakemake/snakemake-wrappers/{version}/{path}/{file}"
FILES: Final = ("meta.yaml", "environment.yaml")
# `wrapper: "v3.0.0/bio/bwa/mem"`, also `"0.72.0/bio/..."`, `"master/bio/..."` and full GitHub URLs
_REF: Final = re.compile(
    r"^(?:https?://github\.com/snakemake/snakemake-wrappers/(?:raw|blob|tree)/)?"
    r"(?P<version>[^/{}\s]+)/(?P<path>(?:bio|utils|meta)/[\w./+-]+?)/?$"
)
# `wrapper:` followed by a string, on the same or the next line
WRAPPER_LINE: Final = re.compile(r"""\bwrapper\s*:\s*\n?\s*[rf]?["']([^"'\n]+)["']""")


@dataclass(frozen=True, slots=True)
class WrapperRef:
    version: str  # tag (v3.0.0), commit or branch
    path: str  # bio/bwa/mem

    @property
    def tool(self) -> str:
        """The tool of the wrapper: bio/bwa/mem -> bwa, meta/bio/bwa_mapping -> bwa_mapping."""
        parts = self.path.split("/")
        if parts[0] == "meta":
            parts = parts[1:]
        return parts[1] if len(parts) > 1 else parts[0]


def wrapper_ref(text: str) -> WrapperRef | None:
    """The wrapper a `wrapper:` value refers to; None for local wrappers (file:...) or unresolvable versions."""
    text = text.strip()
    if text.startswith(("file:", "git+file:")):
        return None
    m = _REF.match(text)
    if m is None:
        # the version from a variable or config: `f"{WRAPPER_VERSION}/bio/bwa/mem"` -> master
        m2 = re.search(r"((?:bio|utils|meta)/[\w./+-]+?)/?$", text)
        return WrapperRef("master", m2.group(1)) if m2 and "{" in text.split(m2.group(1))[0] else None
    return WrapperRef(m.group("version"), m.group("path"))


def wrapper_dir(wrappers_dir: Path, ref: WrapperRef) -> Path:
    return wrappers_dir / ref.version / ref.path


@dataclass(slots=True)
class WrapperInfo:
    """What a wrapper's meta.yaml and environment.yaml say."""

    ref: WrapperRef
    name: str = ""  # meta.yaml name, e.g. "bwa mem"
    homepage: str = ""  # meta.yaml url
    packages: list[str] = field(default_factory=list)  # conda dependencies, e.g. bwa =0.7.19
    found: bool = False  # whether the files were downloaded


def read_wrapper(wrappers_dir: Path, ref: WrapperRef) -> WrapperInfo:
    """The cached meta.yaml / environment.yaml of a wrapper (01_download_snakemake_catalog.py fetches them)."""
    info = WrapperInfo(ref)
    folder = wrapper_dir(wrappers_dir, ref)
    for file in FILES:
        path = folder / file
        if not path.is_file():
            continue
        info.found = True
        try:
            content = yaml.safe_load(path.read_text(encoding="utf-8"))
        except yaml.YAMLError:
            continue
        if not isinstance(content, dict):
            continue
        if file == "meta.yaml":
            info.name, info.homepage = str(content.get("name") or ""), str(content.get("url") or "")
            if not info.homepage:  # older wrappers: `For more information see `BWA <http://bio-bwa...>`_`
                m = re.search(r"https?://[^\s<>`'\")]+", str(content.get("description") or ""))
                info.homepage = m.group(0).rstrip(".,;") if m else ""
        else:
            info.packages = [d for d in content.get("dependencies") or [] if isinstance(d, str)]
    return info
