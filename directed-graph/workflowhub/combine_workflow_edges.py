"""
Combines the tool-pair edges of Galaxy, Nextflow and nf-core workflows into one table.

Reads one edges CSV per source (columns: source, target, n_workflows, ...) and writes:
  - data/combined/combined_edges.csv  one row per directed tool pair:
        source, target, n_wf_galaxy, n_wf_nextflow, n_wf_nfcore
    A pair missing from a source gets 0 for that source.

Tool IDs are lowercased so the same bio.tools ID matches across sources.
If an input has a 'confidence' column, only edges with confidence <= MAX_CONFIDENCE are kept.

Usage (from the repository root):
    .venv/bin/python directed-graph/workflowhub/combine_workflow_edges.py
"""
from pathlib import Path

import pandas as pd

DATA = Path(__file__).resolve().parent / "data/biotools"

edge_filename = "biotools_edges.csv"
INPUTS = {
    "galaxy": DATA / "galaxy" / edge_filename,
    "nextflow": DATA / "nextflow" / edge_filename,
    "nfcore": DATA / "nf-core" / edge_filename,
}
OUT_FILE = DATA / "combined" / "combined_edges.csv"

MAX_CONFIDENCE = 3  # only for galaxy: 1 = codex mapping only, 2 = codex + name match, 3 = all

def read_edges(path: Path, name: str) -> pd.DataFrame:
    """One row per (source, target) with that source's workflow count."""
    df = pd.read_csv(path)
    if "confidence" in df.columns:
        df = df[df["confidence"] <= MAX_CONFIDENCE]
    df["source"] = df["source"].str.lower()
    df["target"] = df["target"].str.lower()
    # lowercasing can merge rows, so sum them
    return (df.groupby(["source", "target"], as_index=False)["n_workflows"].sum()
              .rename(columns={"n_workflows": f"n_wf_{name}"}))


def main() -> None:
    combined = None
    for name, path in INPUTS.items():
        edges = read_edges(path, name)
        print(f"{name}: {len(edges)} edges from {path.relative_to(DATA)}")
        combined = edges if combined is None else combined.merge(edges, on=["source", "target"], how="outer")

    counts = [f"n_wf_{name}" for name in INPUTS]
    combined[counts] = combined[counts].fillna(0).astype(int)
    combined = combined.sort_values(counts, ascending=False)

    OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(OUT_FILE, index=False)
    print(f"Wrote {OUT_FILE.relative_to(DATA)} ({len(combined)} edges)")
    print("edges in all sources:", (combined[counts] > 0).all(axis=1).sum())


if __name__ == "__main__":
    main()