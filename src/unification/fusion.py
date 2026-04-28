"""
Similarity Matrix Fusion Module

Fuses multiple similarity matrices into a single unified similarity matrix.
Uses weighted combination based on user-configured weights.

Users can configure:
- marker_weight: Weight for marker gene similarity
- pathway_weight: Weight for pathway similarity
- niche_weight: Weight for niche similarity
"""

import numpy as np
import pandas as pd
from typing import Dict, Optional, List
import warnings


class SimilarityFusion:
    """Similarity matrix fusion"""
    
    def __init__(
        self,
        marker_weight: float = 0.4,
        pathway_weight: float = 0.3,
        niche_weight: float = 0.3
    ):
        """
        Initialize similarity fusion
        
        Parameters:
            marker_weight: Weight for marker similarity
            pathway_weight: Weight for pathway similarity
            niche_weight: Weight for niche similarity
        """
        self.marker_weight = marker_weight
        self.pathway_weight = pathway_weight
        self.niche_weight = niche_weight
    
    def fuse(
        self,
        marker_sim: np.ndarray,
        pathway_sim: Optional[np.ndarray] = None,
        niche_sim: Optional[np.ndarray] = None,
        weights: Optional[Dict[str, float]] = None
    ) -> np.ndarray:
        """
        Fuse multiple similarity matrices
        
        Parameters:
            marker_sim: Marker gene similarity matrix
            pathway_sim: Pathway similarity matrix (optional)
            niche_sim: Niche similarity matrix (optional)
            weights: Override weights dictionary
        
        Returns:
            Fused similarity matrix
        """
        if marker_sim is None or len(marker_sim) == 0:
            warnings.warn("Marker similarity matrix is empty")
            return np.ndarray([])
        
        n = marker_sim.shape[0]
        
        if weights:
            marker_weight = weights.get("marker", self.marker_weight)
            pathway_weight = weights.get("pathway", self.pathway_weight)
            niche_weight = weights.get("niche", self.niche_weight)
        else:
            marker_weight = self.marker_weight
            pathway_weight = self.pathway_weight
            niche_weight = self.niche_weight
        
        fused_sim = marker_weight * marker_sim
        
        if pathway_sim is not None and np.size(pathway_sim) > 0:
            if pathway_sim.shape == (n, n):
                fused_sim += pathway_weight * pathway_sim
            else:
                warnings.warn("Pathway similarity matrix shape mismatch")
        
        if niche_sim is not None and np.size(niche_sim) > 0:
            if niche_sim.shape == (n, n):
                fused_sim += niche_weight * niche_sim
            else:
                warnings.warn("Niche similarity matrix shape mismatch")
        
        fused_sim = np.clip(fused_sim, 0, 1)
        
        return fused_sim
    
    def normalize_weights(self) -> Dict[str, float]:
        """Normalize weights to sum to 1"""
        total = self.marker_weight + self.pathway_weight + self.niche_weight
        
        if total == 0:
            return {"marker": 0.33, "pathway": 0.33, "niche": 0.33}
        
        return {
            "marker": self.marker_weight / total,
            "pathway": self.pathway_weight / total,
            "niche": self.niche_weight / total
        }
    
    def get_similarity_dataframe(
        self,
        similarity_matrix: np.ndarray,
        subtypes: List[str]
    ) -> pd.DataFrame:
        """
        Convert similarity matrix to DataFrame
        
        Parameters:
            similarity_matrix: Similarity matrix
            subtypes: List of subtype labels
        
        Returns:
            DataFrame of similarity values
        """
        if len(similarity_matrix) == 0 or len(subtypes) == 0:
            return pd.DataFrame()
        
        return pd.DataFrame(
            similarity_matrix,
            index=subtypes,
            columns=subtypes
        )
    
    def get_top_similar_pairs(
        self,
        similarity_matrix: np.ndarray,
        subtypes: List[str],
        top_n: int = 10,
        exclude_self: bool = True
    ) -> pd.DataFrame:
        """
        Get top similar subtype pairs
        
        Parameters:
            similarity_matrix: Similarity matrix
            subtypes: List of subtype labels
            top_n: Number of top pairs
            exclude_self: Whether to exclude self-similarity
        
        Returns:
            DataFrame of top similar pairs
        """
        if len(similarity_matrix) == 0 or len(subtypes) == 0:
            return pd.DataFrame()
        
        pairs = []
        
        for i, subtype_i in enumerate(subtypes):
            for j, subtype_j in enumerate(subtypes):
                if exclude_self and i == j:
                    continue
                
                pairs.append({
                    "subtype_1": subtype_i,
                    "subtype_2": subtype_j,
                    "similarity": similarity_matrix[i, j]
                })
        
        pairs_df = pd.DataFrame(pairs)
        pairs_df = pairs_df.sort_values("similarity", ascending=False)
        
        return pairs_df.head(top_n)
    
    def visualize_similarity(
        self,
        similarity_matrix: np.ndarray,
        subtypes: List[str],
        output_path: Optional[str] = None
    ):
        """
        Visualize similarity matrix as heatmap
        
        Parameters:
            similarity_matrix: Similarity matrix
            subtypes: List of subtype labels
            output_path: Path to save figure
        """
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
        except ImportError:
            warnings.warn("matplotlib or seaborn not installed")
            return
        
        if len(similarity_matrix) == 0 or len(subtypes) == 0:
            return
        
        sim_df = self.get_similarity_dataframe(similarity_matrix, subtypes)
        
        plt.figure(figsize=(10, 8))
        sns.heatmap(
            sim_df,
            annot=True,
            cmap="RdYlBu_r",
            vmin=0,
            vmax=1,
            square=True
        )
        plt.title("Subtype Similarity Matrix")
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches="tight")
        
        plt.close()