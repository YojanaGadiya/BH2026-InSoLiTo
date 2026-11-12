"""
Extract the rule DAG of each downloaded Snakemake workflow, by reading its Snakefiles (nothing is run).

Input:  data/download/repos/<owner>__<repo>/            (from download_snakemake_catalog.py)
        data/download/wrappers/<version>/<path>/       meta.yaml / environment.yaml of the snakemake-wrappers
Output: data/dags/<owner>__<repo>/nodes.csv            one row per tool node, with all tool information
        data/dags/<owner>__<repo>/edges.csv            source, target, kind (dag / chain), via
        data/dags/<owner>__<repo>/rule_edges.csv       the rule graph (rule -> rule), as Snakemake's --rulegraph
        data/dags/<owner>__<repo>/warnings.log         what could not be resolved
        data/dags/snakemake_dags.csv                   one row per workflow: counts, and the check against the
                                                       rule graph the catalog has for some workflows

Step 1 (read):   the main Snakefile and everything it includes, local modules; rules with their sections,
                 functions and variables                                                    (snakefile.py)
Step 2 (paths):  input / output paths of every rule, read with Python's ast: literals, f-strings, expand(),
                 multiext(), rules.X.output, and the paths input functions can return       (paths.py)
Step 3 (graph):  an input path matching another rule's output pattern links the two rules  (build.py)
Step 4 (tools):  per rule: snakemake-wrapper, conda environment, shell commands, container   (tools.py)

The shared parts (shell parser, tool selection, bio.tools mapping) come from nfdag in directed-graph/workflowhub.
"""

import sys
import warnings
from pathlib import Path

# Snakefiles often hold "\d", "\." in normal strings (shell commands, regexes); Python warns about such escapes
# whenever ast parses them, although they are read correctly
warnings.filterwarnings("ignore", category=SyntaxWarning)

WORKFLOWHUB = Path(__file__).resolve().parents[2] / "workflowhub"
if str(WORKFLOWHUB) not in sys.path:
    sys.path.insert(0, str(WORKFLOWHUB))
