"""
Feature fusion and clustering modules for SpaMFC

This package contains modules for:
- weighting: Adaptive weight calculation and feature fusion (FeatureFusion class)
- nmf: NMF consensus decomposition
- clustering: Subtype clustering and stability validation
"""

from .weighting import AdaptiveWeightCalculator, FeatureFusion
from .nmf import NMFConsensus
from .clustering import SubtypeClusterer

__all__ = [
    "FeatureFusion",
    "AdaptiveWeightCalculator",
    "NMFConsensus",
    "SubtypeClusterer",
]