"""
Subtype Mapping Table Generation Module

Generates mapping from original subtypes (sample_cluster) to unified subtypes.
Uses hierarchical clustering on fused similarity matrix.

Users can configure:
- similarity_threshold: Threshold for grouping subtypes
- naming_scheme: Naming scheme for unified subtypes ("Subtype_A", "Subtype_B", ...)
"""

import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering
from typing import Dict, List, Optional
import warnings


class SubtypeMapper:
    """Subtype mapper for cross-sample unification"""
    
    def __init__(
        self,
        similarity_threshold: float = 0.6,
        naming_scheme: str = "letter"
    ):
        """
        Initialize subtype mapper
        
        Parameters:
            similarity_threshold: Threshold for grouping (0-1)
            naming_scheme: Naming scheme ("letter" or "number")
        """
        self.similarity_threshold = similarity_threshold
        self.naming_scheme = naming_scheme
    
    def generate_mapping(
        self,
        fused_similarity: np.ndarray,
        subtypes: List[str],
        threshold: Optional[float] = None
    ) -> Dict[str, str]:
        """
        Generate subtype mapping based on similarity
        
        Parameters:
            fused_similarity: Fused similarity matrix
            subtypes: List of original subtype labels
            threshold: Override similarity threshold
        
        Returns:
            Dictionary mapping original subtype to unified subtype
        """
        if len(fused_similarity) == 0 or len(subtypes) == 0:
            warnings.warn("Empty similarity matrix or subtype list")
            return {}
        
        threshold = threshold or self.similarity_threshold
        
        distance_matrix = 1 - fused_similarity
        
        n_clusters = self._estimate_n_clusters(distance_matrix, threshold)
        
        if n_clusters >= len(subtypes):
            return {s: self._generate_name(i) for i, s in enumerate(subtypes)}
        
        try:
            clustering = AgglomerativeClustering(
                n_clusters=n_clusters,
                metric="precomputed",
                linkage="average"
            )
            cluster_labels = clustering.fit_predict(distance_matrix)
            
        except Exception as e:
            warnings.warn(f"Clustering failed: {e}, using individual mapping")
            return {s: self._generate_name(i) for i, s in enumerate(subtypes)}
        
        unified_names = {}
        
        for cluster_id in range(n_clusters):
            cluster_subtypes = [
                subtypes[i] for i, label in enumerate(cluster_labels)
                if label == cluster_id
            ]
            
            unified_name = self._generate_name(cluster_id)
            
            for subtype in cluster_subtypes:
                unified_names[subtype] = unified_name
        
        return unified_names
    
    def _estimate_n_clusters(
        self,
        distance_matrix: np.ndarray,
        threshold: float
    ) -> int:
        """Estimate number of clusters based on threshold"""
        n = distance_matrix.shape[0]
        
        if n <= 1:
            return n
        
        distances = distance_matrix.flatten()
        distances = distances[distances > 0]
        
        if len(distances) == 0:
            return 1
        
        median_dist = np.median(distances)
        
        if median_dist > threshold:
            return n
        else:
            return max(1, int(n * threshold / median_dist))
    
    def _generate_name(self, cluster_id: int) -> str:
        """Generate unified subtype name"""
        if self.naming_scheme == "letter":
            if cluster_id < 26:
                return f"Subtype_{chr(65 + cluster_id)}"
            else:
                return f"Subtype_{cluster_id + 1}"
        else:
            return f"Subtype_{cluster_id + 1}"
    
    def create_mapping_table(
        self,
        unified_names: Dict[str, str]
    ) -> pd.DataFrame:
        """
        Create mapping table DataFrame
        
        Parameters:
            unified_names: Dictionary of subtype mappings
        
        Returns:
            DataFrame of mapping table
        """
        if len(unified_names) == 0:
            return pd.DataFrame()
        
        rows = []
        
        for original, unified in unified_names.items():
            parts = original.split("_")
            sample = parts[0] if len(parts) > 0 else original
            cluster_id = parts[1] if len(parts) > 1 else "0"
            
            rows.append({
                "original_subtype": original,
                "unified_subtype": unified,
                "sample": sample,
                "cluster_id": cluster_id
            })
        
        return pd.DataFrame(rows)
    
    def apply_mapping(
        self,
        adata,
        subtype_col: str,
        unified_names: Dict[str, str],
        unified_col: Optional[str] = None
    ):
        """
        Apply mapping to AnnData object
        
        Parameters:
            adata: AnnData object
            subtype_col: Column name for original subtype
            unified_names: Dictionary of subtype mappings
            unified_col: Column name for unified subtype (default: f"{subtype_col}_unified")
        
        Returns:
            Updated adata with unified subtype labels
        """
        unified_col = unified_col or f"{subtype_col}_unified"
        
        if subtype_col not in adata.obs:
            warnings.warn(f"Subtype column '{subtype_col}' not found")
            return adata
        
        adata.obs[unified_col] = adata.obs[subtype_col].map(unified_names)
        
        unmapped = adata.obs[unified_col].isna().sum()
        if unmapped > 0:
            warnings.warn(f"{unmapped} cells have unmapped subtypes")
        
        return adata
    
    def get_unified_subtype_counts(
        self,
        adata,
        unified_col: str
    ) -> pd.DataFrame:
        """
        Get counts of each unified subtype
        
        Parameters:
            adata: AnnData object
            unified_col: Column name for unified subtype
        
        Returns:
            DataFrame of unified subtype counts
        """
        if unified_col not in adata.obs:
            warnings.warn(f"Unified column '{unified_col}' not found")
            return pd.DataFrame()
        
        counts = adata.obs[unified_col].value_counts()
        
        return pd.DataFrame({
            "unified_subtype": counts.index,
            "count": counts.values
        })
    
    def validate_mapping(
        self,
        unified_names: Dict[str, str],
        subtypes: List[str]
    ) -> bool:
        """
        Validate mapping completeness
        
        Parameters:
            unified_names: Dictionary of subtype mappings
            subtypes: List of original subtypes
        
        Returns:
            True if mapping is complete
        """
        missing = [s for s in subtypes if s not in unified_names]
        
        if len(missing) > 0:
            warnings.warn(f"Missing mappings for: {missing}")
            return False
        
        return True
    
    def get_mapping_summary(
        self,
        unified_names: Dict[str, str]
    ) -> str:
        """Get mapping summary string"""
        unified_counts = {}
        
        for original, unified in unified_names.items():
            unified_counts[unified] = unified_counts.get(unified, 0) + 1
        
        summary = "Subtype unification summary:\n"
        for unified, count in sorted(unified_counts.items()):
            summary += f"  {unified}: {count} original subtypes\n"
        
        return summary