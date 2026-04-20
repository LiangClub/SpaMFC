"""
Visualization modules for SpaMFC

This package contains modules for:
- spatial_plot: Spatial distribution visualization
- report: Report generation
"""

from .spatial_plot import SpatialVisualizer
from .report import ReportGenerator

__all__ = [
    "SpatialVisualizer",
    "ReportGenerator",
]