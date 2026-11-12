# Nextflow workflows from WorkflowHub

Downloads the Nextflow workflows registered on WorkflowHub, reads the DAG of each pipeline without running it (one node per tool), and maps the tools to bio.tools IDs. The output has the same columns as the Galaxy output, so both can be combined (`workflows/scripts/combine_workflow_edges.py`).

From the repository root:
```
.venv/bin/python workflows/scripts/nextflow/download_workflowhub_nextflow.py   # crates (--ignore-nfcore skips the nf-core space)
.venv/bin/python workflows/scripts/nextflow/extract_nextflow_dags.py           # DAG of each pipeline, one node per tool
.venv/bin/python workflows/scripts/nextflow/map_nextflow_tools.py              # tool nodes -> bio.tools IDs
```
They need PyYAML and git. The mapping also needs the bio.tools dump `directed-graph/biotools/data/biotools_raw*.jsonl` and the Galaxy Codex table `workflows/data/mapping/galaxy/galaxy_codex_tools.tsv`, and it reuses helpers of `../galaxy/map_galaxy_tools.py`.

## Download

`download_workflowhub_nextflow.py` downloads every Nextflow workflow on WorkflowHub. Metadata and crates are shared with the Galaxy download, in `workflows/data/galaxy/download/`. Each workflow is unpacked to `workflows/data/nextflow/download/<id>_v<version>/`.
- **Incomplete crates** (only a link to a git repository, empty or HTML `.nf` files, included files missing) are filled from a shallow `git clone` of the registered repository. `--no-clone` skips that. Each crate is tried once, and `.cloned_from` records the result.
- **`nextflow_workflows.csv`** has one row per workflow. Its column `nf_core` marks the workflows of the WorkflowHub space nf-core. `--summary-only` rebuilds the file from the saved files.

## DAG

`extract_nextflow_dags.py` (code in `nfdag/`) reads the `.nf` files without running them, for DSL1 and DSL2, with every if/else branch (the union over all parameter choices).
- **One node per tool:** each process is split into the tools it runs, in script order. `bwa-mem2 | samtools` becomes `…/bwa-mem2 → …/samtools` (`kind=chain`). Edges between processes (`kind=dag`) leave the last tool of one process and reach the first tool of the next.
- **Which tools:** the executables of a process, restricted to the software its module declares (meta.yml, conda, biocontainers). A process that only runs its own scripts keeps one node, named after the process.
- **Excluded tools:** tools in `nextflow_excluded_tools.txt` never become nodes and don't cut the flow. The list covers compression, language runtimes, containers and reporting tools, and accepts globs like `r-*`.
- **Output** in `workflows/data/nextflow/dags/<id>_v<version>/`:
  - `nodes.csv`: everything found per tool node (executables, own scripts, containers, conda packages, versions.yml tools, nf-core meta.yml with bio.tools IDs);
  - `edges.csv`;
  - `warnings.log`: what could not be resolved;
  - an overview, `nextflow_dags.csv`.

## Mapping the tools to bio.tools

`map_nextflow_tools.py` (logic in `nfdag/mapping.py`) gives each tool node a bio.tools ID. The goal is an ID for every tool that is not excluded or dropped. The reference is the full bio.tools dump: every ID used must exist in it. IDs from meta.yml, the Codex or the curated file that are not in bio.tools are not used (curated ones are reported).

Each tool node is mapped by its tool name: `bwa-mem2`, `bedGraphToBigWig`, or the process name for processes that only run their own scripts. The methods are tried in this order, and the first match wins:

| method | the tool gets the bio.tools ID ... | example |
|---|---|---|
| `curated` | from its row in `nextflow_curated_map.tsv` (`map`, or `drop`: the node is removed and its neighbours are connected) | `bwa-mem2` -> bwa |
| `meta` | that the module's nf-core `meta.yml` gives for this tool (`identifier: biotools:...`) | `stitch` -> stitch-snijderlab |
| `biotools_id` | whose ID is the tool name | `samtools`, `numpy` |
| `biotools_name` | whose name is the tool name, if exactly one entry has that name | `picard` -> picard_tools ("Picard") |
| `biotools_url` | whose homepage, link or download URL names the tool as a GitHub/GitLab/Bitbucket repository or a CRAN/Bioconductor/PyPI/bioconda package. Exactly one entry may have it, and its ID or name must resemble the tool name. | `hmmer` -> hmmer3 |
| `codex` | that the Galaxy Codex gives for a suite or bioconda package of that name (only a fallback) | `wigToBigWig` -> UCSC_Genome_Browser_Utilities |
| `prefix` | whose ID is the process name before the first `_`, only for processes without a tool of their own, and not for generic words | `SAMTOOLS_SORT` -> samtools |

Details:
- **Spellings:** names are compared in lower case, without `.py`, `.sh`, `.pl`, `.r`, `.exe` or `.jar`, and with `-`, `_` and space treated as equal.
- **Group prefixes:** the `biotools_*` methods also try the name without a bioconda group prefix (`r-`, `bioconductor-`, `python-`, `py-`, `perl-`), so `r-ggplot2` -> ggplot2.
- **Declared package:** every method tries the tool name first, then the package the module declares for it (`bowtie2-build` -> bowtie2).
- **Unmapped tools** are removed with their edges, as Galaxy tools without an ID are. This leaves a gap: `A -> unmapped -> B` gives no `A -> B` (`lost_edges` counts the edges lost).
- **Same ID at both ends:** an edge whose two ends get the same ID is dropped (`samtools view -> samtools sort`).
- **Coverage:** the script prints the share of tool nodes (not dropped) that have a bio.tools ID.

**Confidence** (1 = best, 3 = worst), from the mapping method: 1 = `curated`, `meta`; 2 = `biotools_id`, `biotools_name`, `biotools_url`, `codex`; 3 = `prefix`.
- **Pairs:** a pair of tool nodes has the worse level of its two tools. A pair of one workflow has the best level among its tool pairs.
- **Edges and nodes:** an edge has the best level over its workflows, a node the best level of its own mapping.
- **Per level:** edges and nodes also have `confidence1`, `confidence2`, `confidence3`, the number of workflows at each level (they add up to `n_workflows`).

## Outputs

| file | content |
|---|---|
| `workflows/data/mapping/nextflow/nextflow_tool_map.tsv` | one row per tool name: `biotoolsID`, `biotools_name`, `method`, the candidates `codex_biotoolsID` and `match_biotoolsID`, and where the tool is used |
| `workflows/data/mapping/nextflow/nextflow_to_curate.tsv` | the tools to curate (see below) |
| `workflows/data/biotools/nextflow/biotools_pairs.csv`, `biotools_edges.csv`, `biotools_nodes.csv` | the columns of the Galaxy files, with the bio.tools names in place of the ToolShed names, and `confidence` |
| `workflows/data/biotools/nextflow/biotools_{pairs,edges,nodes}_without_nfcore.csv` | the same without the workflows of the WorkflowHub space nf-core, because the nf-core pipelines come directly from nf-core (`workflows/scripts/nf-core/`) |
| `workflows/data/biotools/nextflow/workflows/<id>_v<version>.csv` | the bio.tools DAG of each workflow (`kind` = dag/chain/bridged), for Cytoscape; not tracked |

## Curation

`nextflow_to_curate.tsv` lists the tools without an ID, and those mapped only by `prefix`, most used first. Unmapped tools get up to five `suggestions` from the dump, which are only hints.

Its first columns are those of `nextflow_curated_map.tsv` (next to this README, tracked). Fill in `action` and copy the row over:
- `map`, with a `biotoolsID`; several separated by `;` make a chain;
- `drop`: not a bio.tools tool, so its neighbours are connected.

Then run `map_nextflow_tools.py` again. The DAGs don't need to be rebuilt.

`nextflow_excluded_tools.txt` is for helpers that run inside a process next to the real tool (`megahit | pigz`). It is read by `extract_nextflow_dags.py`, so the DAGs have to be rebuilt after a change. A process that runs nothing but an excluded tool (nf-core `GUNZIP`, `UNTAR`) becomes a node named after the process, and needs a `drop` row.

## Shared with Snakemake

The Snakemake catalog scripts (`../snakemake`) import these parts, so changes there affect both:
- the shell parser (`nfdag/shell.py`);
- the tool selection (`nfdag/process.py`, `choose_steps`);
- the mapping (`nfdag/mapping.py`);
- the output files (`nfdag/report.py`);
- `load_lookup` in `map_nextflow_tools.py`.

Lint and type check (settings in `pyproject.toml`), from this folder: `ruff check . && ruff format --check nfdag && ty check`
