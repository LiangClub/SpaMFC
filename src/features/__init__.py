"""
Feature extraction modules for SpaMFC

This package contains modules for extracting different types of features:
- spatial: Spatial neighborhood features
- cnv: CNV (Copy Number Variation) features
- expression: Gene expression features
- niche: Niche features (scNiche)
"""

from .spatial import SpatialFeatureExtractor
from .cnv import CNVFeatureProcessor
from .expression import ExpressionFeatureProcessor
from .niche import NicheFeatureProcessor

__all__ = [
    "SpatialFeatureExtractor",
    "CNVFeatureProcessor",
    "ExpressionFeatureProcessor",
    "NicheFeatureProcessor",
]