"""
SpaMFC CNV Inference Module

Provides CNV inference functionality using inferCNVpy:
- CNVInferencer: Main CNV inference class
- genomic_positions: Gene position handling
- cnv_scoring: CNV score calculation
"""

from .infercnv import CNVInferencer, run_cnv_inference
from .genomic_positions import (
    add_gene_location,
    add_genomic_positions,
    load_gtf_file,
    validate_genomic_positions,
    get_chromosome_order,
    filter_genes_by_chromosome,
)
from .cnv_scoring import (
    compute_cnv_score,
    get_cnv_summary,
    compute_cnv_correlation,
    classify_cells_by_cnv,
    get_chromosome_cnv_stats,
)

__all__ = [
    "CNVInferencer",
    "run_cnv_inference",
    "add_gene_location",
    "add_genomic_positions",
    "load_gtf_file",
    "validate_genomic_positions",
    "get_chromosome_order",
    "filter_genes_by_chromosome",
    "compute_cnv_score",
    "get_cnv_summary",
    "compute_cnv_correlation",
    "classify_cells_by_cnv",
    "get_chromosome_cnv_stats",
]