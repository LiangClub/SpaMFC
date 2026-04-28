"""
SpaMFC CNV Inference Core Module

Provides CNV inference functionality using inferCNVpy.
Core features:
1. CNV inference from single-cell transcriptomics
2. Reference cell selection
3. CNV clustering and visualization preparation
"""

import numpy as np
import pandas as pd
from typing import Optional, List, Dict, Any
import warnings
from pathlib import Path

try:
    import infercnvpy as cnv
    INFERCNVPY_AVAILABLE = True
except ImportError:
    INFERCNVPY_AVAILABLE = False
    warnings.warn("infercnvpy not installed. CNV inference will be disabled.")

try:
    import scanpy as sc
except ImportError:
    warnings.warn("scanpy not installed")


class CNVInferencer:
    def __init__(
        self,
        window_size: int = 100,
        step_size: int = 10,
        dynamic_threshold: float = 0.5,
        exclude_genes: Optional[List[str]] = None,
        compute_scores: bool = True,
        clustering_resolution: float = 0.5,
        n_pcs: int = 30,
        n_neighbors: int = 15,
        verbose: bool = True
    ):
        self.window_size = window_size
        self.step_size = step_size
        self.dynamic_threshold = dynamic_threshold
        self.exclude_genes = exclude_genes or []
        self.compute_scores = compute_scores
        self.clustering_resolution = clustering_resolution
        self.n_pcs = n_pcs
        self.n_neighbors = n_neighbors
        self.verbose = verbose
        
        if not INFERCNVPY_AVAILABLE:
            raise ImportError("infercnvpy is required for CNV inference. Install with: pip install infercnvpy")
    
    def add_genomic_positions(
        self,
        adata,
        gtf_file: str
    ):
        if self.verbose:
            print(f"[CNV] Adding genomic positions from {gtf_file}")
        
        cnv.io.genomic_position_from_gtf(gtf_file, adata)
        
        if "chromosome" not in adata.var.columns:
            warnings.warn("Chromosome information not found in adata.var")
        
        if self.verbose:
            n_genes_with_pos = adata.var["chromosome"].notna().sum()
            print(f"[CNV] {n_genes_with_pos} genes with genomic positions")
        
        return adata
    
    def infercnv(
        self,
        adata,
        reference_key: str,
        reference_cat: List[str],
        layer: Optional[str] = None,
        inplace: bool = True
    ):
        if self.verbose:
            print(f"[CNV] Running infercnv with reference: {reference_cat}")
        
        if "chromosome" not in adata.var.columns:
            raise ValueError("Genomic positions not found. Run add_genomic_positions first.")
        
        cnv.tl.infercnv(
            adata,
            reference_key=reference_key,
            reference_cat=reference_cat,
            window_size=self.window_size,
            step=self.step_size,
            dynamic_threshold=self.dynamic_threshold,
            layer=layer,
            inplace=inplace
        )
        
        if self.verbose:
            cnv_shape = adata.obsm["X_cnv"].shape
            print(f"[CNV] CNV matrix shape: {cnv_shape}")
        
        return adata
    
    def compute_cnv_score(
        self,
        adata,
        inplace: bool = True
    ):
        if self.verbose:
            print("[CNV] Computing CNV scores")
        
        cnv.tl.cnv_score(adata, inplace=inplace)
        
        if self.verbose and inplace:
            if "cnv_score" in adata.obs.columns:
                print(f"[CNV] CNV score range: {adata.obs['cnv_score'].min():.2f} - {adata.obs['cnv_score'].max():.2f}")
        
        return adata
    
    def pca(
        self,
        adata,
        use_rep: str = "cnv",
        n_pcs: Optional[int] = None,
        inplace: bool = True
    ):
        n_pcs = n_pcs or self.n_pcs
        
        if self.verbose:
            print(f"[CNV] Running PCA with {n_pcs} components")
        
        cnv.tl.pca(adata, use_rep=use_rep, n_comps=n_pcs, inplace=inplace)
        
        return adata
    
    def neighbors(
        self,
        adata,
        use_rep: str = "cnv_pca",
        n_neighbors: Optional[int] = None,
        inplace: bool = True
    ):
        n_neighbors = n_neighbors or self.n_neighbors
        
        if self.verbose:
            print(f"[CNV] Computing neighbors with {n_neighbors} neighbors")
        
        cnv.pp.neighbors(adata, use_rep=use_rep, n_neighbors=n_neighbors, inplace=inplace)
        
        return adata
    
    def leiden(
        self,
        adata,
        resolution: Optional[float] = None,
        key_added: str = "cnv_leiden",
        inplace: bool = True
    ):
        resolution = resolution or self.clustering_resolution
        
        if self.verbose:
            print(f"[CNV] Running Leiden clustering with resolution {resolution}")
        
        cnv.tl.leiden(adata, resolution=resolution, key_added=key_added, inplace=inplace)
        
        if self.verbose and inplace:
            n_clusters = adata.obs[key_added].nunique()
            print(f"[CNV] Found {n_clusters} CNV clusters")
        
        return adata
    
    def umap(
        self,
        adata,
        inplace: bool = True
    ):
        if self.verbose:
            print("[CNV] Computing UMAP")
        
        cnv.tl.umap(adata, inplace=inplace)
        
        return adata
    
    def run_pipeline(
        self,
        adata,
        gtf_file: str,
        reference_key: str,
        reference_cat: List[str],
        layer: Optional[str] = None,
        run_clustering: bool = True,
        run_umap: bool = True
    ):
        if self.verbose:
            print("\n" + "="*60)
            print("SpaMFC CNV Inference Pipeline")
            print("="*60)
        
        self.add_genomic_positions(adata, gtf_file)
        
        self.infercnv(adata, reference_key, reference_cat, layer)
        
        if self.compute_scores:
            self.compute_cnv_score(adata)
        
        self.pca(adata)
        
        self.neighbors(adata)
        
        if run_clustering:
            self.leiden(adata)
        
        if run_umap:
            self.umap(adata)
        
        if self.verbose:
            print("="*60)
            print("CNV Inference Complete")
            print("="*60 + "\n")
        
        return adata
    
    def get_cnv_summary(self, adata) -> Dict[str, Any]:
        summary = {}
        
        if "X_cnv" in adata.obsm:
            cnv_matrix = adata.obsm["X_cnv"]
            summary["cnv_matrix_shape"] = cnv_matrix.shape
            if hasattr(cnv_matrix, "toarray"):
                cnv_dense = cnv_matrix.toarray()
                summary["cnv_mean"] = np.mean(cnv_dense)
                summary["cnv_std"] = np.std(cnv_dense)
            else:
                summary["cnv_mean"] = np.mean(cnv_matrix)
                summary["cnv_std"] = np.std(cnv_matrix)
        
        if "cnv_score" in adata.obs.columns:
            summary["cnv_score_mean"] = adata.obs["cnv_score"].mean()
            summary["cnv_score_std"] = adata.obs["cnv_score"].std()
            summary["cnv_score_range"] = (adata.obs["cnv_score"].min(), adata.obs["cnv_score"].max())
        
        if "cnv_leiden" in adata.obs.columns:
            summary["n_cnv_clusters"] = adata.obs["cnv_leiden"].nunique()
            summary["cluster_sizes"] = adata.obs["cnv_leiden"].value_counts().to_dict()
        
        return summary


def run_cnv_inference(
    adata,
    gtf_file: str,
    reference_key: str,
    reference_cat: List[str],
    window_size: int = 100,
    step_size: int = 10,
    resolution: float = 0.5,
    verbose: bool = True,
    **kwargs
):
    inferencer = CNVInferencer(
        window_size=window_size,
        step_size=step_size,
        clustering_resolution=resolution,
        verbose=verbose
    )
    
    return inferencer.run_pipeline(
        adata,
        gtf_file,
        reference_key,
        reference_cat,
        **kwargs
    )