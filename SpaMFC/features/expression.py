"""
Gene Expression Feature Processing Module

Users can choose whether to use expression features via configuration:
- config.use_expression = True/False

Expression features are applicable to all cell types.
Processing steps:
1. Extract expression matrix for target cells
2. Low expression gene filtering
3. PCA dimensionality reduction
4. UMAP dimensionality reduction (optional)
5. Keep selected PCA components
"""

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from typing import Optional, List
from scipy import sparse
import warnings

try:
    import scanpy as sc
except ImportError:
    warnings.warn("scanpy not installed, some functions may not work")

try:
    import umap
except ImportError:
    warnings.warn("umap not installed, UMAP reduction will be skipped")


class ExpressionFeatureProcessor:
    """Gene expression feature processor"""
    
    def __init__(self, config):
        """
        Initialize expression feature processor
        
        Parameters:
            config: FeatureConfig object containing expression feature settings
        """
        self.config = config
        self.use = config.use
        self.min_expr_pct = config.min_expr_pct
        self.pca_dim = config.pca_dim
        self.use_umap = config.use_umap
        self.umap_dim = config.umap_dim
        self.pca_keep_dim = config.pca_keep_dim
    
    def process(
        self,
        adata,
        target_cells: List[str]
    ) -> Optional[pd.DataFrame]:
        """
        Process gene expression features
        
        Parameters:
            adata: AnnData object
            target_cells: List of target cell names
        
        Returns:
            DataFrame of expression features, or None if use=False
        """
        if not self.use:
            return None
        
        if len(target_cells) == 0:
            warnings.warn("No target cells provided")
            return None
        
        target_cells_valid = [c for c in target_cells if c in adata.obs_names]
        
        if len(target_cells_valid) == 0:
            warnings.warn("No valid target cells found in adata")
            return None
        
        adata_target = adata[target_cells_valid, :].copy()
        
        min_cells = max(1, int(len(adata_target) * self.min_expr_pct))
        
        try:
            sc.pp.filter_genes(adata_target, min_cells=min_cells)
        except Exception:
            pass
        
        if adata_target.n_vars == 0:
            n_vars_to_keep = min(100, adata.n_vars)
            adata_target = adata[target_cells_valid, :n_vars_to_keep].copy()
        
        expr_data = adata_target.X
        if sparse.issparse(expr_data):
            expr_data = expr_data.toarray()
        expr_data = np.asarray(expr_data, dtype=np.float64)
        
        expr_scaled = StandardScaler().fit_transform(expr_data)
        
        n_pca = min(self.pca_dim, len(target_cells_valid) - 1, expr_scaled.shape[1])
        if n_pca < 2:
            warnings.warn("Not enough samples for PCA, returning raw features")
            return pd.DataFrame(
                expr_scaled[:, :min(10, expr_scaled.shape[1])],
                index=target_cells_valid,
                columns=[f"Expr_raw_{i}" for i in range(min(10, expr_scaled.shape[1]))]
            )
        
        pca = PCA(n_components=n_pca, random_state=0)
        expr_pca = pca.fit_transform(expr_scaled)
        
        expr_dim_df = pd.DataFrame(index=target_cells_valid)
        
        if self.use_umap:
            try:
                umap_model = umap.UMAP(
                    n_components=self.umap_dim,
                    random_state=0,
                    n_neighbors=min(5, len(target_cells_valid) - 1)
                )
                expr_umap = umap_model.fit_transform(expr_pca)
                expr_umap_df = pd.DataFrame(
                    expr_umap,
                    index=target_cells_valid,
                    columns=[f"Expr_UMAP{i+1}" for i in range(self.umap_dim)]
                )
                expr_dim_df = pd.concat([expr_dim_df, expr_umap_df], axis=1)
            except (ValueError, RuntimeError, ImportError) as e:
                warnings.warn(f"UMAP failed: {e}")
        
        pca_keep = min(self.pca_keep_dim, n_pca)
        expr_pca_df = pd.DataFrame(
            expr_pca[:, :pca_keep],
            index=target_cells_valid,
            columns=[f"Expr_PC{i+1}" for i in range(pca_keep)]
        )
        expr_dim_df = pd.concat([expr_dim_df, expr_pca_df], axis=1)
        
        expr_dim_df = expr_dim_df.fillna(0.0)
        
        return expr_dim_df
    
    def process_per_sample(
        self,
        adata,
        target_cells: List[str],
        sample_col: str = "sample_id"
    ) -> Optional[pd.DataFrame]:
        """
        Process expression features per sample
        
        Parameters:
            adata: AnnData object
            target_cells: List of target cell names
            sample_col: Column name for sample ID
        
        Returns:
            DataFrame of expression features
        """
        if not self.use:
            return None
        
        all_features = []
        
        samples = adata.obs.loc[target_cells, sample_col].unique()
        
        for sample in samples:
            sample_cells = [
                c for c in target_cells
                if c in adata.obs_names and adata.obs.loc[c, sample_col] == sample
            ]
            
            if len(sample_cells) == 0:
                continue
            
            sample_features = self.process(adata, sample_cells)
            
            if sample_features is not None:
                all_features.append(sample_features)
        
        if len(all_features) == 0:
            return None
        
        return pd.concat(all_features, axis=0)
    
    def get_feature_names(self) -> List[str]:
        """Get names of expression features"""
        names = []
        if self.use_umap:
            names.extend([f"Expr_UMAP{i+1}" for i in range(self.umap_dim)])
        names.extend([f"Expr_PC{i+1}" for i in range(self.pca_keep_dim)])
        return names