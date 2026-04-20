"""
Niche Analysis Module

Performs niche analysis using scNiche for multi-view clustering.
Analyzes niche composition and subtype-niche associations.

Users can configure:
- n_views: Number of views for scNiche
- target_k: Target number of niches
- resolution: Clustering resolution
"""

import numpy as np
import pandas as pd
from typing import Dict, Optional, List
from sklearn.metrics.pairwise import cosine_similarity
import warnings

try:
    import scniche as sn
except ImportError:
    warnings.warn("scNiche not installed, niche analysis will not work")


class NicheAnalyzer:
    """Niche analyzer using scNiche"""
    
    def __init__(
        self,
        n_views: int = 3,
        target_k: int = 15,
        resolution: float = 0.3
    ):
        """
        Initialize niche analyzer
        
        Parameters:
            n_views: Number of views for scNiche
            target_k: Target number of niches
            resolution: Clustering resolution
        """
        self.n_views = n_views
        self.target_k = target_k
        self.resolution = resolution
    
    def run_niche_clustering(
        self,
        adata,
        celltype_col: str,
        add_key: str = "scNiche"
    ):
        """
        Run scNiche clustering
        
        Parameters:
            adata: AnnData object
            celltype_col: Column name for cell type annotation
            add_key: Key to store niche labels
        
        Returns:
            Updated adata with niche labels
        """
        try:
            import scniche as sn
        except ImportError:
            warnings.warn("scNiche not installed, skipping niche clustering")
            return adata
        
        try:
            sn.pp.setup_anndata(adata, groupby=celltype_col)
            
            model = sn.tr.scNiche(
                adata,
                config={
                    "n_views": self.n_views,
                    "view_names": ["X_C2L", "X_data", "X_data_nbr"]
                }
            )
            
            model.train()
            
            adata = sn.tr.clustering(
                adata,
                target_k=self.target_k,
                clustering_method="kmeans",
                resolution=self.resolution,
                add_key=add_key
            )
            
            return adata
            
        except Exception as e:
            warnings.warn(f"scNiche clustering failed: {e}")
            return adata
    
    def get_niche_composition(
        self,
        adata,
        niche_col: str = "scNiche",
        celltype_col: str
    ) -> pd.DataFrame:
        """
        Calculate niche composition (cell type proportions in each niche)
        
        Parameters:
            adata: AnnData object
            niche_col: Column name for niche labels
            celltype_col: Column name for cell type annotation
        
        Returns:
            DataFrame of niche composition
        """
        if niche_col not in adata.obs:
            warnings.warn(f"Niche column '{niche_col}' not found")
            return pd.DataFrame()
        
        niche_labels = adata.obs[niche_col].unique()
        cell_types = adata.obs[celltype_col].unique()
        
        composition = {}
        
        for niche in niche_labels:
            niche_cells = adata[adata.obs[niche_col] == niche]
            
            if len(niche_cells) == 0:
                continue
            
            type_counts = niche_cells.obs[celltype_col].value_counts()
            type_props = type_counts / len(niche_cells)
            
            composition[niche] = type_props
        
        composition_df = pd.DataFrame(composition).T
        composition_df = composition_df.fillna(0)
        
        return composition_df
    
    def get_subtype_niche_association(
        self,
        adata,
        subtype_col: str,
        niche_col: str = "scNiche"
    ) -> pd.DataFrame:
        """
        Calculate subtype-niche association matrix
        
        Parameters:
            adata: AnnData object
            subtype_col: Column name for subtype labels
            niche_col: Column name for niche labels
        
        Returns:
            DataFrame of subtype-niche association
        """
        if subtype_col not in adata.obs:
            warnings.warn(f"Subtype column '{subtype_col}' not found")
            return pd.DataFrame()
        
        if niche_col not in adata.obs:
            warnings.warn(f"Niche column '{niche_col}' not found")
            return pd.DataFrame()
        
        subtype_mask = adata.obs[subtype_col].notna()
        adata_subtype = adata[subtype_mask]
        
        subtypes = adata_subtype.obs[subtype_col].unique()
        niches = adata_subtype.obs[niche_col].unique()
        
        association = {}
        
        for subtype in subtypes:
            subtype_cells = adata_subtype[adata_subtype.obs[subtype_col] == subtype]
            
            if len(subtype_cells) == 0:
                continue
            
            niche_counts = subtype_cells.obs[niche_col].value_counts()
            niche_props = niche_counts / len(subtype_cells)
            
            association[subtype] = niche_props
        
        association_df = pd.DataFrame(association).T
        association_df = association_df.fillna(0)
        
        return association_df
    
    def get_niche_similarity(
        self,
        niche_composition: pd.DataFrame
    ) -> np.ndarray:
        """
        Calculate niche similarity based on composition
        
        Parameters:
            niche_composition: DataFrame of niche composition
        
        Returns:
            Similarity matrix
        """
        if niche_composition.empty:
            return np.ndarray([])
        
        composition_matrix = niche_composition.values
        
        similarity = cosine_similarity(composition_matrix)
        
        return similarity
    
    def get_subtype_niche_profiles(
        self,
        adata,
        subtype_col: str,
        niche_col: str = "scNiche"
    ) -> Dict[str, Dict]:
        """
        Get niche profiles for each subtype
        
        Parameters:
            adata: AnnData object
            subtype_col: Column name for subtype labels
            niche_col: Column name for niche labels
        
        Returns:
            Dictionary of niche profiles for each subtype
        """
        association_df = self.get_subtype_niche_association(adata, subtype_col, niche_col)
        
        profiles = {}
        
        for subtype in association_df.index:
            niche_props = association_df.loc[subtype].to_dict()
            
            dominant_niches = [
                niche for niche, prop in niche_props.items()
                if prop > 0.1
            ]
            
            profiles[subtype] = {
                "composition": niche_props,
                "dominant_niches": dominant_niches,
                "primary_niche": max(niche_props, key=niche_props.get) if niche_props else None
            }
        
        return profiles
    
    def compare_niche_profiles(
        self,
        profiles: Dict[str, Dict]
    ) -> pd.DataFrame:
        """
        Compare niche profiles across subtypes
        
        Parameters:
            profiles: Dictionary of niche profiles
        
        Returns:
            DataFrame of niche profile comparison
        """
        comparison_rows = []
        
        for subtype, profile in profiles.items():
            comparison_rows.append({
                "subtype": subtype,
                "primary_niche": profile.get("primary_niche", None),
                "dominant_niches": ",".join(profile.get("dominant_niches", [])),
                "niche_diversity": len(profile.get("dominant_niches", []))
            })
        
        return pd.DataFrame(comparison_rows)