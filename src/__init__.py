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
"""

from .config import ConfigManager, SpaMFCConfig, FeatureConfig, CorrelationConfig
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

__version__ = "2.1.0"
__author__ = "SpaMFC Team"

__all__ = [
    "SpaMFCPipeline",
    "ConfigManager",
    "SpaMFCConfig",
    "FeatureConfig",
    "CorrelationConfig",
    "create_parser",
    "cli_main",
    "GeneCorrelationAnalyzer",
    "gene_correlation",
    "CorrelationVisualizer",
    "plot_correlation",
    "plot_correlation_heatmap",
    "load_matrices_from_npz",
    "load_matrices_from_csv_gz",
]