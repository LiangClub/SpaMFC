"""
SpaMFC Gene Correlation Analysis Module

Provides gene correlation analysis functionality:
- GeneCorrelationAnalyzer: Main analyzer class
- gene_correlation: Convenience function
- CorrelationVisualizer: Visualization class
- plot_correlation: Single pair scatter plot
- plot_correlation_heatmap: Correlation matrix heatmap
- load_matrices_from_npz: Load NPZ format matrices
- load_matrices_from_csv_gz: Load CSV.GZ format matrices
"""

from .gene_correlation import (
    GeneCorrelationAnalyzer,
    gene_correlation,
    load_matrices_from_npz,
    load_matrices_from_csv_gz
)
from .correlation_plot import (
    CorrelationVisualizer,
    plot_correlation,
    plot_correlation_heatmap
)

__all__ = [
    "GeneCorrelationAnalyzer",
    "gene_correlation",
    "CorrelationVisualizer",
    "plot_correlation",
    "plot_correlation_heatmap",
    "load_matrices_from_npz",
    "load_matrices_from_csv_gz",
]