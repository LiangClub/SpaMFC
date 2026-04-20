"""
Multi-dimensional Similarity Calculation Module

Calculates similarity between subtypes across samples based on:
1. Marker gene similarity (Jaccard index)
2. Functional pathway similarity (pathway overlap)
3. Niche composition similarity (cosine similarity)

Users can configure weights for each dimension:
- marker_weight: Weight for marker gene similarity
- pathway_weight: Weight for pathway similarity
- niche_weight: Weight for niche similarity
"""

import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
from typing import Dict, List, Optional
import warnings


class SimilarityCalculator:
    """Multi-dimensional similarity calculator"""
    
    def __init__(
        self,
        marker_weight: float = 0.4,
        pathway_weight: float = 0.3,
        niche_weight: float = 0.3
    ):
        """
        Initialize similarity calculator
        
        Parameters:
            marker_weight: Weight for marker gene similarity
            pathway_weight: Weight for pathway similarity
            niche_weight: Weight for niche similarity
        """
        self.marker_weight = marker_weight
        self.pathway_weight = pathway_weight
        self.niche_weight = niche_weight
    
    def calculate_marker_similarity(
        self,
        markers_dict: Dict[str, List[str]]
    ) -> np.ndarray:
        """
        Calculate marker gene similarity using Jaccard index
        
        Parameters:
            markers_dict: Dictionary of marker genes for each subtype
        
        Returns:
            Similarity matrix
        """
        subtypes = list(markers_dict.keys())
        n = len(subtypes)
        
        if n == 0:
            return np.ndarray([])
        
        similarity = np.zeros((n, n))
        
        for i, subtype_i in enumerate(subtypes):
            for j, subtype_j in enumerate(subtypes):
                genes_i = set(markers_dict.get(subtype_i, []))
                genes_j = set(markers_dict.get(subtype_j, []))
                
                union = genes_i | genes_j
                intersection = genes_i & genes_j
                
                if len(union) > 0:
                    similarity[i, j] = len(intersection) / len(union)
                else:
                    similarity[i, j] = 0.0
        
        return similarity
    
    def calculate_pathway_similarity(
        self,
        enrichment_results: Dict[str, Dict],
        gene_set: str = "KEGG_2021_Human"
    ) -> np.ndarray:
        """
        Calculate functional pathway similarity
        
        Parameters:
            enrichment_results: Dictionary of enrichment results
            gene_set: Gene set to use for comparison
        
        Returns:
            Similarity matrix
        """
        subtypes = list(enrichment_results.keys())
        n = len(subtypes)
        
        if n == 0:
            return np.ndarray([])
        
        similarity = np.zeros((n, n))
        
        for i, subtype_i in enumerate(subtypes):
            for j, subtype_j in enumerate(subtypes):
                results_i = enrichment_results.get(subtype_i, {})
                results_j = enrichment_results.get(subtype_j, {})
                
                if gene_set in results_i and gene_set in results_j:
                    df_i = results_i[gene_set]
                    df_j = results_j[gene_set]
                    
                    if len(df_i) > 0 and len(df_j) > 0:
                        pathways_i = set(df_i["Term"].tolist())
                        pathways_j = set(df_j["Term"].tolist())
                        
                        union = pathways_i | pathways_j
                        intersection = pathways_i & pathways_j
                        
                        if len(union) > 0:
                            similarity[i, j] = len(intersection) / len(union)
                        else:
                            similarity[i, j] = 0.0
                    else:
                        similarity[i, j] = 0.0
                else:
                    similarity[i, j] = 0.0
        
        return similarity
    
    def calculate_niche_similarity(
        self,
        niche_profiles: Dict[str, Dict]
    ) -> np.ndarray:
        """
        Calculate niche composition similarity using cosine similarity
        
        Parameters:
            niche_profiles: Dictionary of niche profiles for each subtype
        
        Returns:
            Similarity matrix
        """
        subtypes = list(niche_profiles.keys())
        n = len(subtypes)
        
        if n == 0:
            return np.ndarray([])
        
        compositions = []
        
        for subtype in subtypes:
            profile = niche_profiles.get(subtype, {})
            composition = profile.get("composition", {})
            
            if isinstance(composition, dict):
                comp_values = list(composition.values())
            else:
                comp_values = [0]
            
            compositions.append(comp_values)
        
        max_len = max(len(c) for c in compositions) if compositions else 0
        
        compositions_padded = []
        for c in compositions:
            padded = c + [0] * (max_len - len(c))
            compositions_padded.append(padded)
        
        compositions_matrix = np.array(compositions_padded)
        
        if compositions_matrix.shape[1] == 0:
            return np.zeros((n, n))
        
        similarity = cosine_similarity(compositions_matrix)
        
        return similarity
    
    def calculate_expression_similarity(
        self,
        adata,
        subtype_col: str,
        markers_dict: Dict[str, List[str]]
    ) -> np.ndarray:
        """
        Calculate expression profile similarity for marker genes
        
        Parameters:
            adata: AnnData object
            subtype_col: Column name for subtype labels
            markers_dict: Dictionary of marker genes
        
        Returns:
            Similarity matrix
        """
        subtypes = list(markers_dict.keys())
        n = len(subtypes)
        
        if n == 0 or subtype_col not in adata.obs:
            return np.ndarray([])
        
        all_markers = set()
        for genes in markers_dict.values():
            all_markers.update(genes)
        
        available_markers = [g for g in all_markers if g in adata.var_names]
        
        if len(available_markers) == 0:
            return np.zeros((n, n))
        
        subtype_means = {}
        
        for subtype in subtypes:
            subtype_mask = adata.obs[subtype_col] == subtype
            
            if subtype_mask.sum() == 0:
                subtype_means[subtype] = np.zeros(len(available_markers))
                continue
            
            subtype_data = adata[subtype_mask, available_markers].X
            
            if hasattr(subtype_data, "toarray"):
                subtype_data = subtype_data.toarray()
            
            mean_expr = np.mean(subtype_data, axis=0)
            subtype_means[subtype] = mean_expr
        
        means_matrix = np.array([subtype_means[s] for s in subtypes])
        
        similarity = cosine_similarity(means_matrix)
        
        return similarity
    
    def get_subtype_labels(self, markers_dict: Dict) -> List[str]:
        """Get list of subtype labels"""
        return list(markers_dict.keys())