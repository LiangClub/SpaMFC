#!/usr/bin/env python
"""
SpaMFC Command Line Interface Entry Point

Usage:
    SpaMFC run --input data.h5ad --celltype "CAFs"
    SpaMFC run-multi --input data.h5ad --celltypes "Malignant cells,CAFs,ILC"
    SpaMFC info --input data.h5ad
    SpaMFC config --output ./my_config.yaml
    SpaMFC corr --input data.h5ad --target-genes EGFR,KRAS --de-genes GENE1,GENE2
    SpaMFC --help
"""

import sys
from pathlib import Path

def main():
    sys.path.insert(0, str(Path(__file__).parent))
    from src.cli import main as cli_main
    cli_main()

if __name__ == "__main__":
    main()