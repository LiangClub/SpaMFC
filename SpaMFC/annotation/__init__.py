"""
Subtype annotation modules for SpaMFC

This package contains modules for:
- markers: Marker gene analysis (Wilcoxon test)
- enrichment: Functional enrichment analysis (GO/KEGG)
- niche_analysis: Niche analysis (scNiche)
"""

from .markers import MarkerGeneAnalyzer
from .enrichment import EnrichmentAnalyzer
from .niche_analysis import NicheAnalyzer

__all__ = [
    "MarkerGeneAnalyzer",
    "EnrichmentAnalyzer",
    "NicheAnalyzer",
]