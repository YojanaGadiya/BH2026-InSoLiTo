"""
Extract the rule DAG of each downloaded Snakemake workflow (see smkdag/__init__.py), without running anything.

Usage (from the repository root, after download_snakemake_catalog.py):
    python3 directed-graph/snakemake/extract_snakemake_dags.py
    python3 directed-graph/snakemake/extract_snakemake_dags.py --only snakemake-workflows__rna-seq-star-deseq2 --verbose
"""

from smkdag.cli import main

if __name__ == "__main__":
    main()
