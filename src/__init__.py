"""
SpaMFC: Spatial Multi-Feature Clustering for Spatial Transcriptomics Subtype Analysis

A modular pipeline for spatial transcriptomics subtype identification with:
- Flexible feature selection (spatial, CNV, expression, niche)
- Multi-modal feature fusion with adaptive weighting
- NMF consensus clustering
- Subtype annotation (marker genes, functional enrichment, niche analysis)
- Cross-sample subtype unification
- Complete command line interface
- Gene correlation analysis (Pearson, Spearman, Kendall)
- CNV inference and visualization (inferCNVpy)
"""

from .config import ConfigManager, SpaMFCConfig, FeatureConfig, CorrelationConfig, CNVInferenceConfig
from .pipeline import SpaMFCPipeline
from .cli import create_parser, main as cli_main
from .correlation import (
    GeneCorrelationAnalyzer,
    gene_correlation,
    CorrelationVisualizer,
    plot_correlation,
    plot_correlation_heatmap,
    load_matrices_from_npz,
    load_matrices_from_csv_gz,
)
from .cnv_inference import (
    CNVInferencer,
    run_cnv_inference,
    add_genomic_positions,
    load_gtf_file,
    compute_cnv_score,
    get_cnv_summary,
)
from .visualization import (
    CNVVisualizer,
    plot_chromosome_heatmap,
    plot_cnv_umap,
    plot_cnv_scores,
)

__version__ = "2.2.0"
__author__ = "SpaMFC Team"

__all__ = [
    "SpaMFCPipeline",
    "ConfigManager",
    "SpaMFCConfig",
    "FeatureConfig",
    "CorrelationConfig",
    "CNVInferenceConfig",
    "create_parser",
    "cli_main",
    "GeneCorrelationAnalyzer",
    "gene_correlation",
    "CorrelationVisualizer",
    "plot_correlation",
    "plot_correlation_heatmap",
    "load_matrices_from_npz",
    "load_matrices_from_csv_gz",
    "CNVInferencer",
    "run_cnv_inference",
    "add_genomic_positions",
    "load_gtf_file",
    "compute_cnv_score",
    "get_cnv_summary",
    "CNVVisualizer",
    "plot_chromosome_heatmap",
    "plot_cnv_umap",
    "plot_cnv_scores",
]