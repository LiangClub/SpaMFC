"""
Cross-sample subtype unification modules for SpaMFC

This package contains modules for:
- similarity: Multi-dimensional similarity calculation
- fusion: Similarity matrix fusion
- mapping: Subtype mapping table generation
"""

from .similarity import SimilarityCalculator
from .fusion import SimilarityFusion
from .mapping import SubtypeMapper

__all__ = [
    "SimilarityCalculator",
    "SimilarityFusion",
    "SubtypeMapper",
]