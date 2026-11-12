# nf-core software pairs

Directed software pairs ("A → B") taken from the nf-core pipelines, with a bio.tools ID for every tool.

A pair means: a process that uses software A is followed by a process that uses software B in the same workflow, either because B reads an output of A, or because B depends on A (a key, a verdict or metadata), or because A and B run side by side on the same input.

## Source

- **Catalogue:** the nf-core pipeline catalogue (`nfcore_pipelines_catalogue.json`), downloaded with `scripts/1_download/downloader.py`. Archived pipelines are included, and old-style (DSL1) pipelines are read too.
- **Workflows processed: 158**, every workflow version in the catalogue. All workflows, including nf-core/rnaseq, go through exactly the same steps.
- **Tool identifiers:** the bio.tools dump `biotools_raw_14-10-2026.jsonl` (34,488 records) and the bio.tools API.

## Method (short)

1. **Download.** Only the text files that describe a pipeline are fetched from each repository (`main.nf`, `workflows/`, `subworkflows/`, `modules/**/*.nf`, config, `meta.yml` and `environment.yml` files). Nothing is run.
2. **Static analysis.** For every workflow, the Nextflow code is read to build a graph of which process output goes into which process input. Conditions are used to drop paths that can never happen together.
3. **Software per process.** Each process is mapped to one software: the tool named in the module `meta.yml`, otherwise the conda package of its `environment.yml`, otherwise the module name. Old-style pipelines have no modules, so the software is taken from the container and conda directives and the commands in the script.
4. **Software pairs per workflow.** Process edges are checked (an edge that cannot carry data, for example because a `map` drops the element or a `multiMap` splits it, is rejected) and collapsed to software pairs. Every pair gets a `pair_type`:
   - `data_flow`: an output of A is an input of B;
   - `gate_metadata`: B depends on A only through a key, a verdict or metadata;
   - `parallel_co_use`: A and B read the same upstream output, no path between them, same run (stored in both directions);
   - `rejected`: every edge between them was rejected (not a valid pair).

   Mutually exclusive alternatives (A or B, never both, for example `--aligner star` or `hisat2`) are dropped.
5. **Merge over all workflows.** Pairs are merged by tool, and the number of workflows where each pair is valid is counted. Tools with the same bio.tools ID are one tool (`star_align` counts as `star`).
6. **bio.tools IDs.** From the module `meta.yml`, then from the local dump by DOI, other IDs (PMID, PMCID, RRID), homepage URL and last by name. Stronger evidence wins over weaker. Generic words (`cat`, `find`, `filter`, `gzip` ...) never get an ID by name.
7. **Filter.** Only pairs where both tools have a bio.tools ID and that are valid in at least one workflow are kept.
8. **Tables.** `3_nodes_edges/build_tables.py` writes the final tables.

### `direction_type`

| value | meaning |
|---|---|
| `one_direction` | A → B happens in a sequence, B → A never does |
| `both_directions` | A → B and B → A both happen in sequences (in the same or in different workflows) |
| `parallel` | A and B only run side by side, no path between them |

If a pair is a sequence in some workflows and parallel in others, it is `one_direction` or `both_directions`: the sequence is the stronger information.

## Final tables (in `workflows/data/nf-core/`)

| file | one row per | main columns |
|---|---|---|
| `biotools_pairs.csv` | workflow and directed tool pair (the same pair once for every workflow where it is valid) | `workflow_id`, `source`, `target`, `parallel`. `parallel` is `parallel` when the two tools run side by side in that workflow (same input, no path between them), and `not parallel` when one comes after the other. The flag is exact for every workflow |
| `biotools_edges.csv` | directed tool pair | `source`, `target`, `n_workflows`, `workflow_ids`. `n_workflows` is the number of workflows where this direction appears. A to B and B to A are two rows with their own counts. A `parallel` pair appears once forward and once reverse in each workflow, so its two rows have the same count |
| `biotools_nodes.csv` | tool | `biotoolsID`, `biotools_name`, `n_workflows`, `matching_method` |

`source` and `target` in the edges and pairs are bio.tools IDs. `matching_method` says how the ID was found: `biotool_id` (written in the module `meta.yml`), `doi`, `url` or `name`.

## Limits

- The analysis is static: a pair means the code allows it, not that it was seen in a run.
- No pair was checked by hand. The software of a process is read from the module metadata or guessed from the environment, the container or the script (old-style pipelines), so some names are not the tool a researcher would name.
- `name` matches to bio.tools are the weakest and can be wrong for tools with a common name.
- Utilities and reporting tools (for example MultiQC) are in the tables. They can be filtered by hand.

## Intermediate files (`workflows/data/nf-core/intermediate/`)

What the scripts write on the way to the final tables, so the tables can be traced back.

| file | content |
|---|---|
| `pairs_merged_unique.csv` | one row per directed pair over all workflows: in how many workflows it is valid, which, `direction_type`, bio.tools IDs |
| `pairs_merged_summary.md` | counts and the most frequent pairs |
| `control_pool.csv` | valid pairs where both tools have at least one identifier |
| `tools_with_ids.csv`, `tool_identifiers.csv` | tools with at least one ID and why a bio.tools ID is missing; every identifier we hold for every tool |
| `pairs_biotools.csv` | the pairs with a bio.tools ID on both sides, all columns; the final tables are made from it |
| `batch_summary.csv` | status of every pipeline |
| `pipelines/<pipeline>/` | `software_pairs.csv` (the software pairs of one pipeline with their `pair_type`), `process_edges.csv` (the process edges behind them: class, payload check, conditions, call sites), `process_software_map.csv` (software of each process), `summary.json` |

Not in git (see `.gitignore`): the bio.tools dump `biotools_raw_14-10-2026.jsonl` (about 100 MB, a public input: put it into the intermediate folder, or pass its place with `--dump`; it is only read and can stay where it is), the files that only serve to check the ID matching (`biotools_match_report.csv`, `biotools_match_summary.csv`, `biotools_cache.json`) and the logs.

## Scripts (`workflows/scripts/nf-core/`)

| folder | content |
|---|---|
| `1_download/` | `downloader.py`: downloads the pipeline files from GitHub. Edit the paths and options at the top of the file, run it with `DRY_RUN = True` first. Its output (the catalogue) is not part of this folder. |
| `2_pairs/` | the analysis, in five files. `run.py` is the entry point; `pairs.py` builds the pairs of one pipeline (using `graph_dsl2.py` for modern and `graph_dsl1.py` for old-style pipelines); `tool_ids.py` finds the bio.tools IDs; `merge_pairs.py` merges all pipelines and writes `direction_type`. |
| `3_nodes_edges/` | `build_tables.py`: filters the pairs and writes the final tables. |

Python 3.10 or newer and `pip install -r requirements.txt` (pyyaml, networkx; `build_tables.py` also needs pandas).

### How to reproduce

```bash
# 1. download (edit the settings at the top first)
python 1_download/downloader.py

# 2. pairs. CATALOGUE is the folder with the downloaded pipelines (<pipeline>/<version>/main.nf); results go to workflows/data/nf-core/intermediate
python 2_pairs/run.py setup
python 2_pairs/run.py survey  --catalogue CATALOGUE      # what is in the catalogue
python 2_pairs/run.py analyze --catalogue CATALOGUE      # software pairs of every pipeline, a few minutes
python 2_pairs/run.py ids     --catalogue CATALOGUE --lookup   # tool identifiers (needs internet)
python 2_pairs/run.py match                              # put the bio.tools dump biotools_raw*.jsonl into workflows/data/nf-core/intermediate first
python 2_pairs/run.py idtable                            # tools_with_ids.csv
python 2_pairs/run.py merge                              # pairs_merged_unique.csv, control_pool.csv

# 3. the final tables: reads workflows/data/nf-core/intermediate, writes workflows/data/nf-core. The bio.tools dump is found in the intermediate folder, or give its place (it can stay where it is)
python 2_pairs/run.py tables --dump /path/to/biotools_raw_14-10-2026.jsonl
```

`run.py all --catalogue CATALOGUE --lookup` runs all the steps from `analyze` to `tables` one after the other. The intermediate files go to `workflows/data/nf-core/intermediate` and the final tables to `workflows/data/nf-core`. The commands above are run from this folder (`workflows/scripts/nf-core`).
