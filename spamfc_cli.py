#!/usr/bin/env python
"""
SpaMFC Command Line Interface Entry Point

Usage:
    spamfc_cli run --input data.h5ad --celltype "CAFs"
    spamfc_cli run-multi --input data.h5ad --celltypes "Malignant cells,CAFs,ILC"
    spamfc_cli info --input data.h5ad
    spamfc_cli config --output ./my_config.yaml
    spamfc_cli --help
"""

import sys
import os

def main():
    from src.cli import main as cli_main
    cli_main()

if __name__ == "__main__":
    main()