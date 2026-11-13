# Snakemake Workflow Catalog

Downloads the workflows of the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog), reads the rule DAG of each one without running it (one node per tool), and maps the tools to bio.tools IDs. The result is a list of tool pairs (tool A feeds tool B) to compare with the citation graph.

From the repository root:
```
.venv/bin/python workflows/scripts/snakemake/01_download_snakemake_catalog.py --standardized   # catalog + git clones
.venv/bin/python workflows/scripts/snakemake/02_extract_snakemake_dags.py                      # rule DAG, one node per tool
.venv/bin/python workflows/scripts/snakemake/03_map_snakemake_tools.py                         # tool nodes -> bio.tools IDs
```
Requirements:
- Python 3.12+, PyYAML and git.
- For the mapping: the bio.tools dump `directed-graph/biotools/data/biotools_raw*.jsonl` and the Galaxy Codex table `workflows/data/mapping/galaxy/galaxy_codex_tools.tsv`.
- Shared code: the shell parser, the tool selection and the bio.tools mapping are imported from `../nextflow/nfdag/`.

## What the catalog provides

The catalog's `data.json` lists about 6,700 GitHub repositories, about 480 of them "standardized". For each one it gives:
- the repository name, description, topics, stars and latest release;
- the snakemake-wrappers it uses;
- for about 50 standardized workflows, the rule graph Snakemake drew on the workflow's test data.

"Standardized" workflows follow the catalog's layout: `.snakemake-workflow-catalog.yml`, `config/README.md`, usually `workflow/Snakefile`, `workflow/rules/*.smk` and `workflow/envs/*.yaml`.

Snakemake has no meta.yml with bio.tools IDs. The tools of a rule come from:
- its snakemake-wrapper (`wrapper: "v3.0.0/bio/bwa/mem"`), whose `meta.yaml` names the tool's homepage and whose `environment.yaml` lists its conda packages;
- its conda environment (`conda: "envs/x.yaml"`);
- the commands of its `shell:`;
- its container.

## Download

`01_download_snakemake_catalog.py` fetches the catalog and clones the selected repositories into `workflows/data/snakemake/download/`.
- **Selection:**
  - `--standardized` or `--non-standardized` picks one group; neither picks all workflows.
  - `--sample N` takes a fixed sample, the same in every run. In the standardized sample, the workflows with a catalog rule graph come first.
  - `--workflow owner/repo` (repeatable) names single workflows.
- **Clones:**
  - A clone is a shallow clone of the latest release tag, else of the default branch.
  - Only the files needed to read the workflow are checked out (Snakefiles in all their spellings, `.smk`, yaml, scripts), so data stored in git isn't downloaded. Cloning runs no code.
  - `.cloned_from` records the commit or the error, and a failed clone isn't retried (delete its folder to retry).
  - Runs can be stopped (Ctrl+C) and restarted: the summary is still written, and finished clones are skipped.
- **Wrappers:** the `meta.yaml` and `environment.yaml` of every wrapper used are fetched. A wrapper that doesn't exist at its version gets a `.missing` file.
- **Summary:**
  - `--summary-only` rebuilds `snakemake_workflows.csv` from the saved files, without downloading.
  - `--resparse` applies changed file patterns to the clones already there.

## DAG

`02_extract_snakemake_dags.py` (code in `smkdag/`) reads the main Snakefile (`workflow/Snakefile`, else `Snakefile`, else the least deep one) and everything it includes: local modules, `use rule ... as ... with:`, functions and variables. Nothing is executed. The Snakefiles are split into logical lines, and each rule section is read with Python's `ast`.
- **Edges (`smkdag/build.py`):**
  - **By path pattern:** an input path of rule B that matches an output pattern of rule A gives A -> B, as Snakemake links rules. If several output patterns match, the most specific one wins (the longest literal text).
  - **Run-time values:** an input whose wildcards are filled in at run time (`expand("results/{tool}/x.bam", tool=config[...])`) counts every rule that can write such a file. The DAG is the union over all configurations.
  - **By reference:** `rules.A.output` gives the edge directly.
  - **Input functions** (`input: get_bams`, lambdas) are mined for the paths they can return (`via=function`). An input that stays unresolved (config or sample-sheet values only) loses only its own edges.
  - **Targets:** rules without output that run no known tool (`all`, `help`) aren't DAG nodes.
  - **Cycles:** an edge that would close a cycle is left out, with a warning.
- **Tools per rule (`smkdag/tools.py`):**
  - **wrapper:** the wrapper's tool (`bio/bwa/mem` -> `bwa`; `bio/benchmark/chm-eval` -> `chm-eval`), declared by its conda package, with the homepage from its `meta.yaml` (older wrappers: the first link of the description).
  - **conda environment:**
    - The executables of `shell:` are restricted to the environment's packages, and several make a chain: `bwa mem | samtools sort` -> `…/bwa -> …/samtools` (`kind=chain`).
    - Without a matching executable, the packages themselves are used, unless the environment is a full export of more than 5 packages.
  - **container:** for a biocontainers image, its tool.
  - **Fallback:** the rule name, for rules that only run their own code.
  - **Not read as commands:** quoted programs and texts (awk, perl, R, echo, …) are not read as commands. Runners (`uv run`, `pixi run`, `conda run`) give the command they run.
  - **Excluded:** tools in `snakemake_excluded_tools.txt` never become nodes. The list covers compression, runtimes, Snakemake itself, environment and container tools, and reporting tools; globs like `r-*` are allowed.
- **Check against the catalog:** for workflows with a catalog rule graph, `snakemake_dags.csv` gives recall and precision of the static rule edges, and `warnings.log` lists the edges not found. That graph holds only the rules of one test configuration, so the precision measured against it is a lower bound.

Output per workflow, in `workflows/data/snakemake/dags/<owner>__<repo>/`:
- `nodes.csv`: one row per tool node: tool, source, declared package, homepage, wrapper, conda environment and packages, container, executables, scripts, and inputs / unresolved inputs.
- `edges.csv`: between tool nodes, with `kind` (dag/chain) and `via` (path/function/rules).
- `rule_edges.csv`: the rule graph.
- `warnings.log`: what could not be resolved.

Options:
- `--only owner__repo` (repeatable) extracts single workflows.
- `--out-dir` writes the DAGs elsewhere, e.g. for a test.

A full run first clears `workflows/data/snakemake/dags/`.

## Mapping

`03_map_snakemake_tools.py` gives each tool node a bio.tools ID. The goal is an ID for every tool that is not excluded or dropped. The reference is the full bio.tools dump: every ID used must exist in it. The methods are tried in this order, the first match wins, and `snakemake_tool_map.tsv` records the method:

| method | the tool gets the bio.tools ID ... | example | confidence |
|---|---|---|---|
| `curated` | from its row in `snakemake_curated_map.tsv` (`map`, or `drop`: the node is removed and its neighbours are connected) | | 1 |
| `wrapper_homepage` | whose bio.tools homepage is the wrapper's homepage. Trailing path parts of the wrapper URL may be removed if the entry then resembles the tool, but never down to a host serving many tools (github.com needs owner and repository). | `bwameth` -> bwa-meth | 2 |
| `biotools_id` | whose ID is the tool name | `samtools` | 2 |
| `biotools_name` | whose name is the tool name, if exactly one entry has that name | `picard` -> picard_tools | 2 |
| `biotools_url` | whose links name the tool as a GitHub/GitLab/Bitbucket repository or a CRAN/Bioconductor/PyPI/bioconda package, if exactly one entry has it and it resembles the tool name | `hmmer` -> hmmer3 | 2 |
| `codex` | that the Galaxy Codex gives for a suite or bioconda package of that name (only a fallback) | | 2 |
| `prefix` | whose ID is the rule name before the first `_`, only for rules without a tool of their own | `bwa_mem` -> bwa | 3 |

Details:
- **Spellings:** names are compared in lower case, without script suffixes, and with `-`, `_` and space treated as equal.
- **Group prefixes:** the `biotools_*` methods also try the name without a bioconda group prefix (`r-ggplot2` -> ggplot2).
- **Declared package:** every method tries the tool name first, then the package declared for it.
- **Unmapped tools** are removed with their edges. This leaves a gap: `A -> unmapped -> B` gives no `A -> B` (`lost_edges` counts the edges lost).
- **Same ID at both ends:** an edge whose two ends get the same ID is dropped.
- **Coverage:** the script prints the share of tool nodes (not dropped) that have a bio.tools ID.

**Confidence levels:**
- **Pairs:** a pair of tool nodes has the worse level of its two tools. A pair of one workflow has the best level among its tool pairs.
- **Edges and nodes:** an edge has the best level over its workflows, a node the best level of its own mapping.
- **Per level:** `confidence1`, `confidence2` and `confidence3` count the workflows at each level (they add up to `n_workflows`).

## Outputs

| file | content |
|---|---|
| `workflows/data/mapping/snakemake/snakemake_tool_map.tsv` | one row per tool name: `biotoolsID`, `biotools_name`, `method`, the candidates `codex_biotoolsID` and `match_biotoolsID`, and the rules and files it is used in |
| `workflows/data/mapping/snakemake/snakemake_to_curate.tsv` | the tools to curate (see below) |
| `workflows/data/biotools/snakemake/biotools_pairs.csv` | one row per workflow and pair: `workflow_id` (owner/repo), `source`, `target`, their bio.tools names, `confidence` |
| `workflows/data/biotools/snakemake/biotools_edges.csv` | one row per pair: `n_workflows`, `workflow_ids`, names, `confidence`, `confidence1-3` |
| `workflows/data/biotools/snakemake/biotools_nodes.csv` | one row per bio.tools ID: name, `n_workflows`, `has_publication`, `confidence`, `confidence1-3` |
| `workflows/data/biotools/snakemake/biotools_{pairs,edges,nodes}_standardized.csv` | the same, standardized workflows only |
| `workflows/data/biotools/snakemake/workflows/<owner>__<repo>.csv` | the bio.tools DAG of each workflow (`kind` = dag/chain/bridged), for Cytoscape; not tracked |

## Curation

`snakemake_to_curate.tsv` lists the tools without an ID, and those mapped only by `prefix`, most used first. Unmapped tools get up to five `suggestions` from the dump, which are only hints.

Its first columns are those of `snakemake_curated_map.tsv` (next to this README, tracked). Fill in `action` and copy the row over:
- `map`, with a `biotoolsID`; several separated by `;` make a chain;
- `drop`: not a bio.tools tool, so its neighbours are connected.

Then run `03_map_snakemake_tools.py` again. The DAGs don't need to be rebuilt.

`snakemake_excluded_tools.txt` is for helpers that run inside a rule next to the real tool. It is read by `02_extract_snakemake_dags.py`, so the DAGs have to be rebuilt after a change. A rule that runs nothing but excluded tools becomes a node named after the rule, and needs a `drop` row.

## Data

```
workflows/data/snakemake/                 not tracked
├── download/                             01_download_snakemake_catalog.py
│   ├── catalog.json                      the catalog's data.json
│   ├── snakemake_workflows.csv           one row per downloaded workflow (standardized, stars, ref, commit, Snakefile)
│   ├── repos/<owner>__<repo>/            sparse clones
│   ├── wrappers/<version>/<path>/        meta.yaml, environment.yaml
│   └── rulegraphs/<owner>__<repo>.dot    the catalog's rule graphs
└── dags/                                 02_extract_snakemake_dags.py: <owner>__<repo>/, snakemake_dags.csv
workflows/data/mapping/snakemake/         03_map_snakemake_tools.py: tool map, to-curate list
workflows/data/biotools/snakemake/        03_map_snakemake_tools.py: bio.tools pairs, edges, nodes
```

Lint and type check (settings in `pyproject.toml`), from this folder: `ruff check . && ruff format --check . && ty check`
