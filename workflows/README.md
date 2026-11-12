# workflows

Tool pairs (tool A feeds tool B) from Galaxy, Nextflow and nf-core workflows, mapped to bio.tools IDs to compare with the citation graph.

## Folders

- `scripts/`  one folder per workflow source, plus the top level `combine_workflow_edges.py` which combines the outputs from all
- `data/`     inputs and outputs, see below

## Scripts

Scripts to download Galaxy and nextflow workflows from WorkflowHub, and nf-core workflows from ... TODO Xenie add 

### Run order

From the repository root:

```
# Galaxy
.venv/bin/python workflows/scripts/galaxy/download_workflowhub_galaxy.py    # .ga files from WorkflowHub
.venv/bin/python workflows/scripts/galaxy/galaxy_tool_pairs.py              # tool pairs per workflow
.venv/bin/python workflows/scripts/galaxy/map_galaxy_tools.py               # ToolShed names tranlsated to bio.tools IDs

# Nextflow
.venv/bin/python workflows/scripts/nextflow/download_workflowhub_nextflow.py  # downloads nextflow files from WorkflowHub
.venv/bin/python workflows/scripts/nextflow/extract_nextflow_dags.py          # extracts graphs from each workflow
.venv/bin/python workflows/scripts/nextflow/map_nextflow_tools.py             # Mapping the tools to bio.tools and creating pairs and edges

# nf-core
.venv/bin/python workflows/scripts/nf-core/1_download/downloader.py.    # downloads the pipeline files from GitHub.
.venv/bin/python workflows/scripts/nf-core/2_pairs/run.py all           # the analysis, in five files

# Combine outputs from all sources
.venv/bin/python workflows/scripts/combine_workflow_edges.py      # combines the workflows from the three sources and the citation data
```

## Data

- `data/galaxy/`    Data related fo downloading Galaxy workflows and their tool pairs (ToolShed names)
- `data/nextflow/`  tsv? TODO: Julian add descritpion
- `data/nf-core/`   Data related fo downloading nf-core workflows and their tool pairs
- `data/mapping/`   tool name -> bio.tools ID maps, per source
- `data/biotools/`  tool pairs with bio.tools IDs, 
  -  per source:
        - `biotools_pairs.csv`  one row per workflow and pair
        - `biotools_edges.csv`  one row per pair, summed over workflows
        - `biotools_nodes.csv`  one row per tool
  - `combined/combined_edges.csv`  all sources in one table, with citation counts