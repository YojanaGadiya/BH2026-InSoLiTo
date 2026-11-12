"""Find bio.tools tools that are equivalent, i.e. list exactly the same set of publications.

Reads data/biotools_publications.csv (made by pull_biotools.py) and writes
data/biotools_equivalent_tools.csv, one row per group of two or more tools with identical
publication sets:
    group_id, n_tools, n_papers, papers, biotoolsIDs

Rows are the same paper when they share any identifier (pmid, pmcid or doi), so a tool citing
only the DOI and another citing only the PMID of a paper still match. Each paper is named like
collection_id in pull_biotools.py: its first identifier in the order pmid, pmcid, doi, with its
kind as prefix and the lowest value when there are several (e.g. "pmid:10827456").
papers is "|"-joined (";" occurs inside SICI-style DOIs), biotoolsIDs is ";"-joined.
group_id is the first of the group's papers. The file is sorted by n_tools, largest first.

Citation order cannot tell the tools of one group apart, since any citation of their papers
counts for all of them.

Usage:
    python3 find_equivalent_tools.py

Only the standard library is needed.
"""
import csv
from collections import defaultdict
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent / "data"
IDS_FILE = DATA_DIR / "biotools_publications.csv"
OUT_FILE = DATA_DIR / "biotools_equivalent_tools.csv"
KINDS = ("pmid", "pmcid", "doi")


def main():
    parent = {}

    def find(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    # Join the identifiers of each row, so identifiers of one paper end up in one set
    tool_ids = []
    with open(IDS_FILE, newline="") as f:
        for row in csv.DictReader(f):
            ids = [(kind, row[kind].strip()) for kind in KINDS if row[kind].strip()]
            if not ids:
                continue
            find(ids[0])  # registers papers with a single identifier too
            for other in ids[1:]:
                parent[find(ids[0])] = find(other)
            tool_ids.append((row["biotoolsID"], ids[0]))

    # Name each paper by its first identifier (pmid, pmcid, doi order, lowest value)
    members = defaultdict(list)
    for ident in list(parent):
        members[find(ident)].append(ident)
    paper_name = {root: min(ids, key=lambda i: (KINDS.index(i[0]), i[1])) for root, ids in members.items()}

    tool_papers = defaultdict(set)
    for tool, ident in tool_ids:
        tool_papers[tool].add(paper_name[find(ident)])

    groups = defaultdict(list)
    for tool, papers in tool_papers.items():
        groups[frozenset(papers)].append(tool)

    sort_key = lambda i: (KINDS.index(i[0]), i[1])
    rows = []
    for papers, tools in groups.items():
        if len(tools) < 2:
            continue
        papers = sorted(papers, key=sort_key)
        names = [f"{kind}:{value}" for kind, value in papers]
        rows.append([names[0], len(tools), len(names), "|".join(names), ";".join(sorted(tools))])
    rows.sort(key=lambda r: (-r[1], r[0]))

    with open(OUT_FILE, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["group_id", "n_tools", "n_papers", "papers", "biotoolsIDs"])
        writer.writerows(rows)

    in_groups = sum(r[1] for r in rows)
    print(f"{len(tool_papers)} tools with publications, {len(groups)} distinct publication sets")
    print(f"{len(rows)} groups of equivalent tools, {in_groups} tools in them")
    print(f"Written {OUT_FILE}")


if __name__ == "__main__":
    main()
