"""Load everything the bio.tools mapping needs: the bio.tools dump, the Galaxy Codex and a curated TSV.

Shared by 03_map_nextflow_tools.py and the Snakemake catalog (workflows/scripts/snakemake/03_map_snakemake_tools.py).
The paths of the bio.tools dump and the Codex, and two helpers, come from ../galaxy/map_galaxy_tools.py.
"""

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "galaxy"))  # workflows/scripts/galaxy
from map_galaxy_tools import (
    BIOTOOLS_PUBS,
    BIOTOOLS_RAW,
    CODEX_FILE,
    PREFIX_BLOCKLIST,
    find_column,
    variants,
)

from .mapping import BiotoolsIndex, Lookup, read_curated


def load_codex_suites() -> dict[str, str]:
    """{Galaxy suite ID -> bio.tools ID}. Unlike map_galaxy_tools.load_codex, without the tool IDs inside a
    suite: those are Galaxy wrapper names (the srst2 suite has a tool `samtools`), not executables."""
    if not CODEX_FILE.exists():
        return {}
    csv.field_size_limit(10**9)
    with CODEX_FILE.open(encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        suite_col = find_column(reader.fieldnames, "suite", "id")
        biotools_col = find_column(reader.fieldnames, "bio.tool", "id")
        if not suite_col or not biotools_col:
            return {}
        out: dict[str, str] = {}
        for row in reader:
            suite, biotools_id = (row.get(suite_col) or "").strip().lower(), (row.get(biotools_col) or "").strip()
            if suite and biotools_id:
                for name in variants(suite):
                    out.setdefault(name, biotools_id)
    return out


def load_with_paper() -> set[str]:
    """The bio.tools IDs with a publication."""
    with BIOTOOLS_PUBS.open(encoding="utf-8") as f:
        return {row["biotoolsID"] for row in csv.DictReader(f)}


def _name_parts(name: str) -> set[str]:
    return {p for p in re.split(r"[-_.]", name.lower()) if len(p) > 2}


def load_codex_conda() -> dict[str, str]:
    """{bioconda package -> bio.tools ID} from the Codex column 'Suite conda package'.
    That is the main requirement of a Galaxy suite, often only a dependency (a suite needing samtools is about
    srst2, one needing pandas about image analysis). So per package the most common bio.tools ID of its suites
    is taken, and only if it or one of those suites shares a name part with the package
    (samtools -> samtools, ncbi-amrfinderplus -> amrfinderplus, pandas -> nothing)."""
    if not CODEX_FILE.exists():
        return {}
    csv.field_size_limit(10**9)
    with CODEX_FILE.open(encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        conda_col = "Suite conda package" if "Suite conda package" in (reader.fieldnames or []) else None
        suite_col = find_column(reader.fieldnames, "suite", "id")
        biotools_col = find_column(reader.fieldnames, "bio.tool", "id")
        if not conda_col or not suite_col or not biotools_col:
            return {}
        suites: dict[str, dict[str, list[str]]] = defaultdict(lambda: defaultdict(list))  # package -> ID -> suites
        for row in reader:
            biotools_id = (row.get(biotools_col) or "").strip()
            suite = (row.get(suite_col) or "").strip()
            for package in (row.get(conda_col) or "").split(",") if biotools_id else []:
                if package.strip():
                    suites[package.strip().lower()][biotools_id].append(suite)
    out = {}
    for package, by_id in suites.items():
        biotools_id, its_suites = max(by_id.items(), key=lambda kv: len(kv[1]))
        names = _name_parts(biotools_id).union(*(_name_parts(s) for s in its_suites))
        if _name_parts(package) & names:
            out[package] = biotools_id
    return out


def load_lookup(curated_file: Path) -> tuple[BiotoolsIndex, set[str], Lookup]:
    """The bio.tools dump, the IDs with a paper, and the lookup with the Codex and the curated TSV."""
    if BIOTOOLS_RAW is None:
        sys.exit("No bio.tools dump (directed-graph/biotools/data/biotools_raw*.jsonl): it is the mapping reference")
    print(f"bio.tools: {BIOTOOLS_RAW}")
    index = BiotoolsIndex.load(BIOTOOLS_RAW)
    print(f"Codex: {CODEX_FILE}" if CODEX_FILE.exists() else f"No {CODEX_FILE.name}: skipping Codex matching")
    curated = read_curated(curated_file)
    unknown = sorted({i for e in curated.values() for i in e.ids if i.lower() not in index.ids})
    if unknown:
        print(f"Curated IDs not found in bio.tools (check spelling): {', '.join(unknown)}")
    return index, load_with_paper(), Lookup(index, load_codex_conda(), load_codex_suites(), curated, PREFIX_BLOCKLIST)
