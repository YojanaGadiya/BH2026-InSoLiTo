"""
Extract the rule DAG of each downloaded Snakemake workflow, by reading its Snakefiles (nothing is run).

Input:  workflows/data/snakemake/download/repos/<owner>__<repo>/       (from 01_download_snakemake_catalog.py)
        workflows/data/snakemake/download/wrappers/<version>/<path>/  meta.yaml / environment.yaml of the wrappers
Output: workflows/data/snakemake/dags/<owner>__<repo>/nodes.csv       one row per tool node, with all tool information
        workflows/data/snakemake/dags/<owner>__<repo>/edges.csv       source, target, kind (dag / chain), via
        workflows/data/snakemake/dags/<owner>__<repo>/rule_edges.csv  the rule graph, as Snakemake's --rulegraph
        workflows/data/snakemake/dags/<owner>__<repo>/warnings.log    what could not be resolved
        workflows/data/snakemake/dags/snakemake_dags.csv              one row per workflow: counts, and the check
                                                                      against the catalog's rule graphs

Step 1 (read):   the main Snakefile and everything it includes, local modules; rules with their sections,
                 functions and variables                                                    (snakefile.py)
Step 2 (paths):  input / output paths of every rule, read with Python's ast: literals, f-strings, expand(),
                 multiext(), rules.X.output, and the paths input functions can return       (paths.py)
Step 3 (graph):  an input path matching another rule's output pattern links the two rules  (build.py)
Step 4 (tools):  per rule: snakemake-wrapper, conda environment, shell commands, container   (tools.py)

The shared parts (shell parser, tool selection, bio.tools mapping) come from nfdag in workflows/scripts/nextflow.
"""

import sys
import warnings
from pathlib import Path

# Snakefiles often hold "\d", "\." in normal strings (shell commands, regexes); Python warns about such escapes
# whenever ast parses them, although they are read correctly
warnings.filterwarnings("ignore", category=SyntaxWarning)

ROOT = Path(__file__).resolve().parents[4]  # smkdag -> snakemake -> scripts -> workflows -> repository root
SCRIPTS = ROOT / "workflows" / "scripts"
DATA = ROOT / "workflows" / "data"
DOWNLOAD_DIR = DATA / "snakemake" / "download"  # catalog, clones, wrappers (01_download_snakemake_catalog.py)
DAG_DIR = DATA / "snakemake" / "dags"  # 02_extract_snakemake_dags.py
MAP_DIR = DATA / "mapping" / "snakemake"  # 03_map_snakemake_tools.py: tool map, to-curate list
OUT_DIR = DATA / "biotools" / "snakemake"  # 03_map_snakemake_tools.py: bio.tools pairs, edges, nodes

# nfdag (shell parser, tool selection, bio.tools mapping and lookup) is shared with the Nextflow scripts
if str(SCRIPTS / "nextflow") not in sys.path:
    sys.path.insert(0, str(SCRIPTS / "nextflow"))
