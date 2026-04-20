"""
Niche Feature Processing Module (scNiche)

Users can choose whether to use niche features via configuration:
- config.use_niche = True/False

Niche features use scNiche for multi-view clustering:
- X_C2L: cell2location cell type abundance
- X_data: gene expression features
- X_data_nbr: neighborhood expression features
"""

import numpy as np
import pandas as pd
from typing import Optional, List, Dict
import warnings

try:
    import scniche as sn
except ImportError:
    warnings.warn("scNiche not installed, niche features will be skipped")


class NicheFeatureProcessor:
    """Niche feature processor using scNiche"""
    
    def __init__(self, config):
        """
        Initialize niche feature processor
        
        Parameters:
            config: FeatureConfig object containing niche feature settings
        """
        self.config = config
        self.use = config.use
        self.n_views = config.n_views
        self.target_k = config.target_k
        self.resolution = config.resolution
    
    def process(
        self,
        adata,
        target_cells: List[str],
        celltype_col: str,
        sample_col: str = "sample_id"
    ) -> Optional[pd.DataFrame]:
        """
        Process niche features using scNiche
        
        Parameters:
            adata: AnnData object
            target_cells: List of target cell names
            celltype_col: Column name for cell type annotation
            sample_col: Column name for sample ID
        
        Returns:
            DataFrame of niche features, or None if use=False
        """
        if not self.use:
            return None
        
        try:
            import scniche as sn
        except ImportError:
            warnings.warn("scNiche not installed, skipping niche features")
            return None
        
        target_cells_valid = [c for c in target_cells if c in adata.obs_names]
        
        if len(target_cells_valid) == 0:
            warnings.warn("No valid target cells for niche processing")
            return None
        
        adata_target = adata[target_cells_valid, :].copy()
        
        try:
            sn.pp.setup_anndata(adata_target, groupby=celltype_col)
            
            model = sn.tr.scNiche(
                adata_target,
                config={
                    "n_views": self.n_views,
                    "view_names": ["X_C2L", "X_data", "X_data_nbr"]
                }
            )
            
            model.train()
            
            adata_target = sn.tr.clustering(
                adata_target,
                target_k=self.target_k,
                clustering_method="kmeans",
                resolution=self.resolution,
                add_key="scNiche"
            )
            
            niche_labels = adata_target.obs["scNiche"].values
            
            niche_df = pd.DataFrame(
                {"niche_label": niche_labels},
                index=target_cells_valid
            )
            
            return niche_df
            
        except Exception as e:
            warnings.warn(f"scNiche processing failed: {e}")
            return None
    
    def get_niche_composition(
        self,
        adata,
        niche_col: str = "scNiche",
        celltype_col: str
    ) -> Optional[pd.DataFrame]:
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
            warnings.warn(f"Niche labels not found in adata.obs['{niche_col}']")
            return None
        
        niche_labels = adata.obs[niche_col].unique()
        cell_types = adata.obs[celltype_col].unique()
        
        composition = {}
        for niche in niche_labels:
            niche_cells = adata[adata.obs[niche_col] == niche]
            type_counts = niche_cells.obs[celltype_col].value_counts()
            type_props = type_counts / len(niche_cells)
            composition[niche] = type_props
        
        composition_df = pd.DataFrame(composition).T
        composition_df = composition_df.fillna(0)
        
        return composition_df
    
    def get_feature_names(self) -> List[str]:
        """Get names of niche features"""
        return ["niche_label"]