"""
NMF Consensus Decomposition Module

Performs multiple NMF runs and takes consensus embedding to improve stability.
Users can configure:
- n_runs: Number of NMF runs (default: 10)
- n_components: Number of NMF factors (default: 5)
"""

import numpy as np
import pandas as pd
from sklearn.decomposition import NMF
from sklearn.preprocessing import MinMaxScaler
from typing import Optional
import warnings


class NMFConsensus:
    """NMF consensus decomposition"""
    
    def __init__(
        self,
        n_runs: int = 10,
        n_components: int = 5,
        init: str = "nndsvda",
        max_iter: int = 1000,
        tol: float = 1e-4
    ):
        """
        Initialize NMF consensus decomposer
        
        Parameters:
            n_runs: Number of NMF runs
            n_components: Number of NMF factors
            init: Initialization method
            max_iter: Maximum iterations
            tol: Convergence tolerance
        """
        self.n_runs = n_runs
        self.n_components = n_components
        self.init = init
        self.max_iter = max_iter
        self.tol = tol
    
    def decompose(
        self,
        fused_df: pd.DataFrame,
        n_components: Optional[int] = None
    ) -> pd.DataFrame:
        """
        Perform NMF consensus decomposition
        
        Parameters:
            fused_df: Fused feature DataFrame
            n_components: Override number of components
        
        Returns:
            NMF consensus embedding DataFrame
        """
        if fused_df is None or len(fused_df) == 0:
            warnings.warn("Empty input DataFrame")
            return None
        
        n_components = n_components or self.n_components
        n_components = min(n_components, len(fused_df) - 1, fused_df.shape[1])
        
        if n_components < 2:
            warnings.warn("Not enough samples for NMF, returning raw features")
            return fused_df.iloc[:, :min(5, fused_df.shape[1])]
        
        scaler = MinMaxScaler()
        fused_scaled = scaler.fit_transform(fused_df.values)
        
        fused_scaled = np.clip(fused_scaled, 0, 1)
        
        nmf_embeds = []
        
        for i in range(self.n_runs):
            try:
                model = NMF(
                    n_components=n_components,
                    init=self.init,
                    random_state=i,
                    max_iter=self.max_iter,
                    tol=self.tol
                )
                W = model.fit_transform(fused_scaled)
                nmf_embeds.append(W)
            except Exception as e:
                warnings.warn(f"NMF run {i} failed: {e}")
                continue
        
        if len(nmf_embeds) == 0:
            warnings.warn("All NMF runs failed, returning PCA features")
            from sklearn.decomposition import PCA
            pca = PCA(n_components=n_components)
            W = pca.fit_transform(fused_scaled)
            nmf_embeds.append(W)
        
        W_consensus = np.mean(np.stack(nmf_embeds, axis=2), axis=2)
        
        W_df = pd.DataFrame(
            W_consensus,
            index=fused_df.index,
            columns=[f"NMF_Factor_{i+1}" for i in range(n_components)]
        )
        
        W_df = W_df.fillna(0.0)
        
        return W_df
    
    def decompose_per_sample(
        self,
        fused_df: pd.DataFrame,
        sample_labels: pd.Series,
        n_components: Optional[int] = None
    ) -> pd.DataFrame:
        """
        Perform NMF decomposition per sample
        
        Parameters:
            fused_df: Fused feature DataFrame
            sample_labels: Series of sample labels for each cell
            n_components: Override number of components
        
        Returns:
            NMF embedding DataFrame
        """
        if fused_df is None or len(fused_df) == 0:
            return None
        
        all_embeddings = []
        
        samples = sample_labels.unique()
        
        for sample in samples:
            sample_cells = sample_labels[sample_labels == sample].index
            sample_cells = [c for c in sample_cells if c in fused_df.index]
            
            if len(sample_cells) == 0:
                continue
            
            sample_df = fused_df.loc[sample_cells]
            
            sample_embedding = self.decompose(sample_df, n_components)
            
            if sample_embedding is not None:
                all_embeddings.append(sample_embedding)
        
        if len(all_embeddings) == 0:
            return None
        
        return pd.concat(all_embeddings, axis=0)
    
    def get_factor_names(self, n_components: int) -> list:
        """Get names of NMF factors"""
        return [f"NMF_Factor_{i+1}" for i in range(n_components)]