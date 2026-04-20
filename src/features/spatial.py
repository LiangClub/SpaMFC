"""
Spatial Neighborhood Feature Extraction Module

Users can choose whether to use spatial features via configuration:
- config.use_spatial = True/False

Features extracted:
- Neighborhood cell type composition (weighted by distance)
- Edge detection (whether cell is at tumor edge)
- Distance to tumor center
"""

import numpy as np
import pandas as pd
from sklearn.neighbors import KDTree
from sklearn.preprocessing import MinMaxScaler
from typing import Optional, List
import warnings

try:
    import scanpy as sc
except ImportError:
    warnings.warn("scanpy not installed, some functions may not work")


class SpatialFeatureExtractor:
    """Spatial neighborhood feature extractor"""
    
    def __init__(self, config):
        """
        Initialize spatial feature extractor
        
        Parameters:
            config: FeatureConfig object containing spatial feature settings
        """
        self.config = config
        self.radius = config.radius
        self.k_neighbors = config.k_neighbors
        self.use = config.use
    
    def extract(
        self,
        adata,
        target_cells: List[str],
        celltype_col: str,
        global_cell_types: List[str],
        spatial_key: str = "spatial"
    ) -> Optional[pd.DataFrame]:
        """
        Extract spatial neighborhood features
        
        Parameters:
            adata: AnnData object
            target_cells: List of target cell names
            celltype_col: Column name for cell type annotation
            global_cell_types: List of all cell types in the dataset
            spatial_key: Key for spatial coordinates in adata.obsm
        
        Returns:
            DataFrame of spatial features, or None if use=False
        """
        if not self.use:
            return None
        
        if spatial_key not in adata.obsm:
            raise ValueError(f"Spatial coordinates not found in adata.obsm['{spatial_key}']")
        
        coords = adata.obsm[spatial_key]
        
        if len(target_cells) == 0:
            warnings.warn("No target cells provided")
            return None
        
        tree = KDTree(coords)
        sigma = self.radius / 2
        
        target_mask = adata.obs_names.isin(target_cells)
        target_coords = coords[target_mask]
        
        if len(target_coords) == 0:
            warnings.warn("No target cells found in adata")
            return None
        
        tumor_center = target_coords.mean(axis=0)
        
        cell_center_dist = np.linalg.norm(target_coords - tumor_center, axis=1)
        if cell_center_dist.max() > 0:
            cell_center_dist = MinMaxScaler().fit_transform(
                cell_center_dist.reshape(-1, 1)
            ).flatten()
        else:
            cell_center_dist = np.zeros(len(target_cells))
        
        env_features = []
        
        for idx, cell in enumerate(target_cells):
            try:
                cell_idx = adata.obs_names.get_loc(cell)
            except KeyError:
                continue
            
            distances, neighbor_indices = tree.query(
                coords[cell_idx].reshape(1, -1),
                k=min(self.k_neighbors, len(coords))
            )
            distances = distances[0]
            neighbor_indices = neighbor_indices[0]
            
            if distances.max() > 0:
                weights = np.exp(-(distances**2) / (2 * sigma**2))
                weights = weights / weights.sum()
            else:
                weights = np.ones(len(distances)) / len(distances)
            
            neighbor_types = adata.obs.iloc[neighbor_indices][celltype_col]
            type_weights = pd.Series(weights, index=neighbor_indices)
            weighted_counts = type_weights.groupby(neighbor_types).sum()
            
            env_vec = pd.Series(0.0, index=list(global_cell_types) + ["is_edge", "distance_to_center"])
            for ct in weighted_counts.index:
                if ct in env_vec.index:
                    env_vec[ct] = weighted_counts[ct]
            
            neighbor_cell_types = adata.obs.iloc[neighbor_indices][celltype_col]
            target_cell_type = adata.obs.loc[cell, celltype_col]
            is_edge = (neighbor_cell_types != target_cell_type).any()
            env_vec["is_edge"] = 1.0 if is_edge else 0.0
            env_vec["distance_to_center"] = cell_center_dist[idx]
            
            env_features.append(env_vec)
        
        if len(env_features) == 0:
            return None
        
        env_df = pd.DataFrame(env_features, index=target_cells)
        env_df = env_df.fillna(0.0)
        
        return env_df
    
    def extract_per_sample(
        self,
        adata,
        target_cells: List[str],
        celltype_col: str,
        sample_col: str,
        spatial_key: str = "spatial"
    ) -> Optional[pd.DataFrame]:
        """
        Extract spatial features per sample (to avoid cross-sample interference)
        
        Parameters:
            adata: AnnData object
            target_cells: List of target cell names
            celltype_col: Column name for cell type annotation
            sample_col: Column name for sample ID
            spatial_key: Key for spatial coordinates
        
        Returns:
            DataFrame of spatial features
        """
        if not self.use:
            return None
        
        global_cell_types = adata.obs[celltype_col].unique().tolist()
        
        all_features = []
        
        samples = adata.obs.loc[target_cells, sample_col].unique()
        
        for sample in samples:
            sample_cells = [
                c for c in target_cells
                if adata.obs.loc[c, sample_col] == sample
            ]
            
            if len(sample_cells) == 0:
                continue
            
            sample_features = self.extract(
                adata,
                sample_cells,
                celltype_col,
                global_cell_types,
                spatial_key
            )
            
            if sample_features is not None:
                all_features.append(sample_features)
        
        if len(all_features) == 0:
            return None
        
        return pd.concat(all_features, axis=0)
    
    def get_feature_names(self, global_cell_types: List[str]) -> List[str]:
        """Get names of spatial features"""
        return list(global_cell_types) + ["is_edge", "distance_to_center"]