"""
CNV (Copy Number Variation) Feature Processing Module

Users can choose whether to use CNV features via configuration:
- config.use_cnv = True/False

CNV features are only applicable to malignant cells.
Processing steps:
1. Extract CNV matrix from adata.obsm
2. Marker redundancy filtering (optional)
3. PCA dimensionality reduction
4. UMAP dimensionality reduction (optional)
"""

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from typing import Optional, List
from scipy import sparse
import warnings

try:
    import umap
except ImportError:
    warnings.warn("umap not installed, UMAP reduction will be skipped")


class CNVFeatureProcessor:
    """CNV feature processor"""
    
    def __init__(self, config):
        """
        Initialize CNV feature processor
        
        Parameters:
            config: FeatureConfig object containing CNV feature settings
        """
        self.config = config
        self.use = config.use
        self.pca_dim = config.pca_dim
        self.use_umap = config.use_umap
        self.umap_dim = config.umap_dim
        self.marker_filter = config.marker_filter
    
    def process(
        self,
        adata,
        target_cells: List[str],
        cnv_key: str = "X_cnv",
        markers: Optional[List[str]] = None
    ) -> Optional[pd.DataFrame]:
        """
        Process CNV features
        
        Parameters:
            adata: AnnData object
            target_cells: List of target cell names
            cnv_key: Key for CNV matrix in adata.obsm
            markers: List of tumor markers for redundancy filtering
        
        Returns:
            DataFrame of CNV features, or None if use=False
        """
        if not self.use:
            return None
        
        if cnv_key not in adata.obsm:
            warnings.warn(f"CNV matrix not found in adata.obsm['{cnv_key}'], skipping CNV features")
            return None
        
        cnv_data = adata.obsm[cnv_key]
        
        if sparse.issparse(cnv_data):
            cnv_data = cnv_data.toarray()
        cnv_data = np.asarray(cnv_data, dtype=np.float64)
        
        target_mask = adata.obs_names.isin(target_cells)
        cnv_target = cnv_data[target_mask]
        
        if len(cnv_target) == 0:
            warnings.warn("No target cells found for CNV processing")
            return None
        
        cnv_df = pd.DataFrame(
            cnv_target,
            index=target_cells,
            columns=[f"CNV_bin_{i}" for i in range(cnv_target.shape[1])]
        )
        
        cnv_df = cnv_df.fillna(0.0)
        cnv_clipped = np.clip(cnv_df.values, -2.0, 2.0)
        
        cnv_scaled = StandardScaler().fit_transform(cnv_clipped)
        
        if self.marker_filter and markers is not None:
            cnv_scaled = self._filter_marker_correlation(
                adata, target_cells, cnv_scaled, markers
            )
        
        n_pca = min(self.pca_dim, len(target_cells) - 1, cnv_scaled.shape[1])
        n_output = min(self.pca_dim, n_pca)
        if n_output < 2:
            warnings.warn("Not enough samples for PCA, returning raw features")
            return pd.DataFrame(
                cnv_scaled[:, :min(10, cnv_scaled.shape[1])],
                index=target_cells,
                columns=[f"CNV_raw_{i}" for i in range(min(10, cnv_scaled.shape[1]))]
            )
        
        pca = PCA(n_components=n_pca, random_state=0)
        cnv_pca = pca.fit_transform(cnv_scaled)
        
        if self.use_umap:
            try:
                umap_model = umap.UMAP(
                    n_components=min(self.umap_dim, n_pca),
                    random_state=0,
                    n_neighbors=min(5, len(target_cells) - 1)
                )
                cnv_umap = umap_model.fit_transform(cnv_pca)
                cnv_dim_df = pd.DataFrame(
                    cnv_umap,
                    index=target_cells,
                    columns=[f"CNV_UMAP{i+1}" for i in range(cnv_umap.shape[1])]
                )
            except (ValueError, RuntimeError, ImportError) as e:
                warnings.warn(f"UMAP failed: {e}, using PCA features instead")
                cnv_dim_df = pd.DataFrame(
                    cnv_pca[:, :n_output],
                    index=target_cells,
                    columns=[f"CNV_PC{i+1}" for i in range(n_output)]
                )
        else:
            cnv_dim_df = pd.DataFrame(
                cnv_pca[:, :n_output],
                index=target_cells,
                columns=[f"CNV_PC{i+1}" for i in range(n_output)]
            )
        
        cnv_dim_df = cnv_dim_df.fillna(0.0)
        
        return cnv_dim_df
    
    def _filter_marker_correlation(
        self,
        adata,
        target_cells: List[str],
        cnv_data: np.ndarray,
        markers: List[str]
    ) -> np.ndarray:
        """
        Filter CNV bins based on correlation with tumor markers
        
        Parameters:
            adata: AnnData object
            target_cells: List of target cell names
            cnv_data: CNV data matrix
            markers: List of tumor marker genes
        
        Returns:
            Filtered CNV data
        """
        target_mask = adata.obs_names.isin(target_cells)
        
        available_markers = [m for m in markers if m in adata.var_names]
        
        if len(available_markers) == 0:
            return cnv_data
        
        expr_data = adata[target_mask, available_markers].X
        if sparse.issparse(expr_data):
            expr_data = expr_data.toarray()
        expr_data = np.asarray(expr_data, dtype=np.float64)
        
        if expr_data.shape[1] == 0:
            return cnv_data
        
        marker_expr_mean = expr_data.mean(axis=1)
        
        correlations = []
        for i in range(cnv_data.shape[1]):
            cnv_bin = cnv_data[:, i]
            if np.std(cnv_bin) > 0 and np.std(marker_expr_mean) > 0:
                corr = np.corrcoef(cnv_bin, marker_expr_mean)[0, 1]
            else:
                corr = 0
            correlations.append(abs(corr))
        
        keep_indices = [i for i, corr in enumerate(correlations) if corr < 0.7]
        
        if len(keep_indices) == 0:
            return cnv_data
        
        return cnv_data[:, keep_indices]
    
    def process_per_sample(
        self,
        adata,
        target_cells: List[str],
        cnv_key: str = "X_cnv",
        sample_col: str = "sample_id",
        markers: Optional[List[str]] = None
    ) -> Optional[pd.DataFrame]:
        """
        Process CNV features per sample
        
        Parameters:
            adata: AnnData object
            target_cells: List of target cell names
            cnv_key: Key for CNV matrix
            sample_col: Column name for sample ID
            markers: List of tumor markers
        
        Returns:
            DataFrame of CNV features
        """
        if not self.use:
            return None
        
        all_features = []
        
        samples = adata.obs.loc[target_cells, sample_col].unique()
        
        for sample in samples:
            sample_cells = [
                c for c in target_cells
                if adata.obs.loc[c, sample_col] == sample
            ]
            
            if len(sample_cells) == 0:
                continue
            
            sample_features = self.process(
                adata,
                sample_cells,
                cnv_key,
                markers
            )
            
            if sample_features is not None:
                all_features.append(sample_features)
        
        if len(all_features) == 0:
            return None
        
        return pd.concat(all_features, axis=0)
    
    def get_feature_names(self, n_features: int) -> List[str]:
        """Get names of CNV features"""
        if self.use_umap:
            return [f"CNV_UMAP{i+1}" for i in range(self.umap_dim)]
        else:
            return [f"CNV_PC{i+1}" for i in range(min(10, n_features))]