"""
SpaMFC CNV Inference Module

Provides CNV inference functionality using inferCNVpy:
- CNVInferencer: Main CNV inference class
- genomic_positions: Gene position handling
- cnv_scoring: CNV score calculation
"""

from .infercnv import CNVInferencer, run_cnv_inference
from .genomic_positions import add_genomic_positions, load_gtf_file
from .cnv_scoring import compute_cnv_score, get_cnv_summary

__all__ = [
    "CNVInferencer",
    "run_cnv_inference",
    "add_genomic_positions",
    "load_gtf_file",
    "compute_cnv_score",
    "get_cnv_summary",
]