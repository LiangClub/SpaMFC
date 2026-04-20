"""
Subtype Clustering and Stability Validation Module

Supports multiple clustering methods:
- KMeans: Traditional K-means clustering
- Leiden: Graph-based Leiden clustering

Users can configure:
- method: Clustering method ("kmeans" or "leiden")
- n_clusters: Number of clusters for KMeans
- resolution: Resolution for Leiden
- per_sample: Whether to cluster per sample (recommended)
- n_bootstrap: Number of bootstrap runs for stability validation
"""

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.utils import resample
from typing import Optional, List, Tuple
import warnings

try:
    import scanpy as sc
except ImportError:
    warnings.warn("scanpy not installed, Leiden clustering will not work")


class SubtypeClusterer:
    """Subtype clusterer with stability validation"""
    
    def __init__(
        self,
        method: str = "kmeans",
        n_clusters: int = 5,
        resolution: float = 0.1,
        per_sample: bool = True,
        n_bootstrap: int = 30,
        stability_threshold: float = 0.7
    ):
        """
        Initialize subtype clusterer
        
        Parameters:
            method: Clustering method ("kmeans" or "leiden")
            n_clusters: Number of clusters for KMeans
            resolution: Resolution for Leiden
            per_sample: Whether to cluster per sample
            n_bootstrap: Number of bootstrap runs
            stability_threshold: Stability threshold
        """
        self.method = method
        self.n_clusters = n_clusters
        self.resolution = resolution
        self.per_sample = per_sample
        self.n_bootstrap = n_bootstrap
        self.stability_threshold = stability_threshold
    
    def cluster(
        self,
        adata,
        W_df: pd.DataFrame,
        sample_col: str,
        celltype: str
    ):
        """
        Perform subtype clustering
        
        Parameters:
            adata: AnnData object
            W_df: NMF embedding DataFrame
            sample_col: Column name for sample ID
            celltype: Cell type name
        
        Returns:
            Updated adata with subtype labels
        """
        if W_df is None or len(W_df) == 0:
            warnings.warn("Empty NMF embedding")
            return adata
        
        subtype_col = f"{celltype}_subtype"
        adata.obs[subtype_col] = np.nan
        
        if self.per_sample:
            adata = self._cluster_per_sample(adata, W_df, sample_col, subtype_col)
        else:
            adata = self._cluster_global(adata, W_df, sample_col, subtype_col)
        
        return adata
    
    def _cluster_per_sample(
        self,
        adata,
        W_df: pd.DataFrame,
        sample_col: str,
        subtype_col: str
    ):
        """Cluster per sample to avoid batch effects"""
        samples = adata.obs[sample_col].unique()
        
        for sample in samples:
            sample_cells = [
                c for c in W_df.index
                if c in adata.obs_names and adata.obs.loc[c, sample_col] == sample
            ]
            
            if len(sample_cells) == 0:
                continue
            
            W_sample = W_df.loc[sample_cells]
            
            n_clusters_sample = min(self.n_clusters, len(sample_cells) - 1)
            if n_clusters_sample < 2:
                for cell in sample_cells:
                    adata.obs.loc[cell, subtype_col] = f"{sample}_0"
                continue
            
            labels = self._cluster_sample(W_sample, n_clusters_sample)
            
            for cell, label in zip(sample_cells, labels):
                adata.obs.loc[cell, subtype_col] = f"{sample}_{label}"
        
        return adata
    
    def _cluster_global(
        self,
        adata,
        W_df: pd.DataFrame,
        sample_col: str,
        subtype_col: str
    ):
        """Cluster globally (may have batch effects)"""
        n_clusters = min(self.n_clusters, len(W_df) - 1)
        
        if n_clusters < 2:
            for cell in W_df.index:
                sample = adata.obs.loc[cell, sample_col]
                adata.obs.loc[cell, subtype_col] = f"{sample}_0"
            return adata
        
        labels = self._cluster_sample(W_df, n_clusters)
        
        for cell, label in zip(W_df.index, labels):
            sample = adata.obs.loc[cell, sample_col]
            adata.obs.loc[cell, subtype_col] = f"{sample}_{label}"
        
        return adata
    
    def _cluster_sample(
        self,
        W_sample: pd.DataFrame,
        n_clusters: int
    ) -> np.ndarray:
        """Cluster a single sample"""
        if self.method == "kmeans":
            kmeans = KMeans(n_clusters=n_clusters, random_state=0, n_init=10)
            labels = kmeans.fit_predict(W_sample.values)
        elif self.method == "leiden":
            try:
                adata_nmf = sc.AnnData(W_sample)
                sc.pp.neighbors(adata_nmf, n_neighbors=min(10, len(W_sample) - 1))
                sc.tl.leiden(adata_nmf, resolution=self.resolution)
                labels = adata_nmf.obs["leiden"].values.astype(int)
            except Exception as e:
                warnings.warn(f"Leiden clustering failed: {e}, using KMeans")
                kmeans = KMeans(n_clusters=n_clusters, random_state=0)
                labels = kmeans.fit_predict(W_sample.values)
        else:
            warnings.warn(f"Unknown method: {self.method}, using KMeans")
            kmeans = KMeans(n_clusters=n_clusters, random_state=0)
            labels = kmeans.fit_predict(W_sample.values)
        
        return labels
    
    def validate_stability(
        self,
        W_df: pd.DataFrame,
        n_bootstrap: Optional[int] = None
    ) -> Tuple[List[int], float]:
        """
        Validate clustering stability using bootstrap
        
        Parameters:
            W_df: NMF embedding DataFrame
            n_bootstrap: Override number of bootstrap runs
        
        Returns:
            Tuple of (stable cluster labels, stability score)
        """
        n_bootstrap = n_bootstrap or self.n_bootstrap
        
        if W_df is None or len(W_df) < 10:
            warnings.warn("Not enough samples for stability validation")
            return [], 0.0
        
        cluster_freq = {}
        
        for i in range(n_bootstrap):
            sample_indices = resample(
                range(len(W_df)),
                replace=True,
                n_samples=int(0.8 * len(W_df))
            )
            
            W_sample = W_df.iloc[sample_indices]
            
            n_clusters = min(self.n_clusters, len(W_sample) - 1)
            if n_clusters < 2:
                continue
            
            labels = self._cluster_sample(W_sample, n_clusters)
            
            for label in set(labels):
                cluster_freq[label] = cluster_freq.get(label, 0) + 1
        
        freq_threshold = n_bootstrap * self.stability_threshold
        stable_labels = [
            label for label, freq in cluster_freq.items()
            if freq >= freq_threshold
        ]
        
        if len(cluster_freq) > 0:
            stability_score = len(stable_labels) / len(cluster_freq)
        else:
            stability_score = 0.0
        
        return stable_labels, stability_score
    
    def get_subtype_counts(
        self,
        adata,
        subtype_col: str
    ) -> pd.Series:
        """Get counts of each subtype"""
        if subtype_col not in adata.obs:
            return pd.Series()
        
        return adata.obs[subtype_col].value_counts()
    
    def get_subtype_statistics(
        self,
        adata,
        subtype_col: str,
        sample_col: str
    ) -> pd.DataFrame:
        """Get statistics of subtypes per sample"""
        if subtype_col not in adata.obs:
            return pd.DataFrame()
        
        stats = []
        
        samples = adata.obs[sample_col].unique()
        
        for sample in samples:
            sample_mask = adata.obs[sample_col] == sample
            subtype_mask = adata.obs[subtype_col].notna()
            
            sample_subtypes = adata.obs[sample_mask & subtype_mask][subtype_col]
            
            for subtype in sample_subtypes.unique():
                count = (sample_subtypes == subtype).sum()
                stats.append({
                    "sample": sample,
                    "subtype": subtype,
                    "count": count,
                    "proportion": count / len(sample_subtypes) if len(sample_subtypes) > 0 else 0
                })
        
        return pd.DataFrame(stats)