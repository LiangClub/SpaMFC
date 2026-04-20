"""
SpaMFC: Spatial Multi-Feature Clustering for Spatial Transcriptomics Subtype Analysis

A modular pipeline for spatial transcriptomics subtype identification with:
- Flexible feature selection (spatial, CNV, expression, niche)
- Multi-modal feature fusion with adaptive weighting
- NMF consensus clustering
- Subtype annotation (marker genes, functional enrichment, niche analysis)
- Cross-sample subtype unification
- Complete command line interface
"""

from .config import ConfigManager, SpaMFCConfig, FeatureConfig
from .pipeline import SpaMFCPipeline
from .cli import create_parser, main as cli_main

__version__ = "2.0.0"
__author__ = "SpaMFC Team"

__all__ = [
    "SpaMFCPipeline",
    "ConfigManager",
    "SpaMFCConfig",
    "FeatureConfig",
    "create_parser",
    "cli_main",
]