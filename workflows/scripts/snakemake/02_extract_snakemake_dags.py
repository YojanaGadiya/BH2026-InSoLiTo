"""
Extract the rule DAG of each downloaded Snakemake workflow (see smkdag/__init__.py), without running anything.

Usage (from the repository root, after 01_download_snakemake_catalog.py):
    python3 workflows/scripts/snakemake/02_extract_snakemake_dags.py
    python3 workflows/scripts/snakemake/02_extract_snakemake_dags.py --only owner__repo --verbose
"""

from smkdag.cli import main

if __name__ == "__main__":
    main()
