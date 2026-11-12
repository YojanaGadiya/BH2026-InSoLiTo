"""
Combines the tool-pair edges of Galaxy, Nextflow and nf-core workflows into one table.

Reads one edges CSV per source (e.g. columns: source, target, n_workflows, ...) and writes:
  - data/combined/combined_edges.csv  one row per directed tool pair:
        source, target, n_wf_galaxy, n_wf_nextflow, n_wf_nfcore, n_articles, median_delta, median_tool_gap
    A pair missing from a source gets 0 for that source.

Tool IDs are matched case-insensitively (tool == TOoL), so the same bio.tools ID matches across sources.
Different spellings in different files are merged and the output uses upper case spelling if there is one (tool and TOoL).
A name spelt in more than one way within the same file (tool and TOoL) means two different tools, which means we cant rely on case
insensitive merging so it stops the script with an error.
If an input has a 'confidence' column, only edges with confidence <= MAX_CONFIDENCE are kept.

Usage (from the repository root):
    .venv/bin/python directed-graph/workflowhub/combine_workflow_edges.py
"""
from collections import Counter
from pathlib import Path

import pandas as pd

DATA = Path(__file__).resolve().parent / "data/biotools"

edge_filename = "biotools_edges.csv"
INPUTS = {
    "galaxy": DATA / "galaxy" / edge_filename,
    "nextflow": DATA / "nextflow" / "biotools_edges_without_nfcore.csv",
    "nfcore": DATA / "nf-core" / edge_filename,
}
CITATIONS = Path(__file__).resolve().parent.parent / "edges-processed" / "directed_edges.csv"
OUT_FILE = DATA / "combined" / "combined_edges.csv"

MAX_CONFIDENCE = 3  # only for galaxy: 1 = codex mapping only, 2 = codex + name match, 3 = all


def pick_spellings() -> dict[str, str]:
    """Returns {lowercase ID: spelling to use}, and prints the tools spelt in more than one way.

    Of several spellings (e.g. tool and TOoL), the one with upper case wins;
    if there are several with upper case, the most common one (in edges, over all files).
    """
    counts = {}    # lowercase -> Counter({spelling: edges}), over all files
    per_file = {}  # lowercase -> {file: {spellings}}
    for name, path in INPUTS.items():
        df = pd.read_csv(path, usecols=["source", "target"])
        for n in pd.concat([df["source"], df["target"]]).astype(str):
            counts.setdefault(n.lower(), Counter())[n] += 1
            per_file.setdefault(n.lower(), {}).setdefault(name, set()).add(n)

    chosen = {low: max(c, key=lambda s: (s != low, c[s], s)) for low, c in counts.items()}

    clashes = {low: c for low, c in counts.items() if len(c) > 1}
    print(f"tools spelt in more than one way: {len(clashes)} of {len(counts)}")
    for low, c in sorted(clashes.items()):
        print(f"  {chosen[low]:<25} from {dict(c)}")  # edge counts per spelling
    print()

    within = {low: f for low, f in per_file.items() if any(len(s) > 1 for s in f.values())}
    print(f"of those, spelt in more than one way within the same file: {len(within)}")
    for low, f in sorted(within.items()):
        print(f"  {chosen[low]:<25} " + "; ".join(f"{name}: {sorted(s)}" for name, s in f.items() if len(s) > 1))
    print()

    # Two spellings in one file are two different tools, so they cannot be matched case-insensitively
    if within:
        raise ValueError(
            f"{len(within)} names are spelt in more than one way within the same file, so they are "
            f"different tools that cannot be matched across files: {sorted(within)}. See the list above."
        )
    return chosen


def read_edges(path: Path, name: str, spelling: dict[str, str]) -> pd.DataFrame:
    """One row per (source, target) with that source's workflow count."""
    df = pd.read_csv(path)
    if "confidence" in df.columns:
        df = df[df["confidence"] <= MAX_CONFIDENCE]
    df["source"] = df["source"].astype(str).str.lower().map(spelling)
    df["target"] = df["target"].astype(str).str.lower().map(spelling)
    # merging spellings can merge rows, so sum them
    return (df.groupby(["source", "target"], as_index=False)["n_workflows"].sum()
              .rename(columns={"n_workflows": f"n_wf_{name}"}))


def read_citations(spelling: dict[str, str]) -> pd.DataFrame:
    """One row per direction: source cited before target in n_articles articles, with the medians.

    A direction with its own row in the file keeps that row. A direction found only in the _rev
    columns of the opposite row is added from there, with n_articles but no medians (the file has none).
    """
    cols = ["source", "target", "n_articles", "median_delta", "median_tool_gap"]
    df = pd.read_csv(CITATIONS)
    fwd = df[cols]
    rev = (df[df["n_articles_rev"] > 0]
           .rename(columns={"source": "target", "target": "source"})
           .assign(n_articles=lambda d: d["n_articles_rev"], median_delta=float("nan"), median_tool_gap=float("nan"))
           [cols])
    out = pd.concat([fwd, rev], ignore_index=True)
    for col in ("source", "target"):
        low = out[col].astype(str).str.lower()
        out[col] = low.map(spelling).fillna(low)
    # forward rows come first, so a direction with its own row keeps it (and its medians)
    out = out.drop_duplicates(["source", "target"], keep="first")
    print(f"citations: {len(fwd)} rows, {len(out) - len(fwd)} directions added from the _rev columns")
    return out


def main() -> None:
    spelling = pick_spellings()
    combined = None
    for name, path in INPUTS.items():
        edges = read_edges(path, name, spelling)
        print(f"{name}: {len(edges)} edges from {path.relative_to(DATA)}")
        combined = edges if combined is None else combined.merge(edges, on=["source", "target"], how="outer")

    counts = [f"n_wf_{name}" for name in INPUTS]
    combined[counts] = combined[counts].fillna(0).astype(int)
    combined = combined.sort_values(counts, ascending=False)

    combined = combined.merge(read_citations(spelling), on=["source", "target"], how="left")
    combined["n_articles"] = combined["n_articles"].fillna(0).astype(int)
    print("workflow edges also cited in the same direction:", (combined["n_articles"] > 0).sum())

    OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(OUT_FILE, index=False)
    print(f"Wrote {OUT_FILE.relative_to(DATA)} ({len(combined)} edges)")
    print("edges in all sources:", (combined[counts] > 0).all(axis=1).sum())


if __name__ == "__main__":
    main()