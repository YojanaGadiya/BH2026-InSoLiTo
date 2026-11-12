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
.venv/bin/python workflows/scripts/nextflow/download_workflowhub_nextflow.py
.venv/bin/python workflows/scripts/nextflow/extract_nextflow_dags.py        # TODO: julian add description
.venv/bin/python workflows/scripts/nextflow/map_nextflow_tools.py           # TODO: julian add description

# nf-core
TODO: Xenie add scripts

# Combine outputs from all sources
.venv/bin/python workflows/scripts/combine_workflow_edges.py
```

## Data

- `data/galaxy/`    downloaded Galaxy workflows and their tool pairs (ToolShed names)
- `data/nextflow/`  tsv? TODO: Julian add descritpion
- `data/mapping/`   tool name -> bio.tools ID maps, per source
- `data/biotools/`  tool pairs with bio.tools IDs, 
  -  per source:
        - `biotools_pairs.csv`  one row per workflow and pair
        - `biotools_edges.csv`  one row per pair, summed over workflows
        - `biotools_nodes.csv`  one row per tool
  - `combined/combined_edges.csv`  all sources in one table, with citation counts