"""
Marker Gene Analysis Module

Performs differential gene analysis to identify subtype-specific marker genes.
Uses Wilcoxon rank-sum test (scanpy's rank_genes_groups).

Users can configure:
- method: Differential test method ("wilcoxon", "t-test", "logreg")
- pval_threshold: P-value threshold for significant genes
- logfc_threshold: Log fold change threshold
- top_n: Number of top marker genes to extract
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional
import warnings
import sys

try:
    import scanpy as sc
except ImportError:
    warnings.warn("scanpy not installed, marker analysis will not work")

MEMORY_WARNING_THRESHOLD_MB = 500


def _estimate_memory_usage(n_cells: int, n_genes: int) -> float:
    estimated_mb = (n_cells * n_genes * 8) / (1024 * 1024)
    return estimated_mb


class MarkerGeneAnalyzer:
    """Marker gene analyzer for subtype identification"""
    
    def __init__(
        self,
        method: str = "wilcoxon",
        pval_threshold: float = 0.05,
        logfc_threshold: float = 0.25,
        top_n: int = 50
    ):
        """
        Initialize marker gene analyzer
        
        Parameters:
            method: Differential test method
            pval_threshold: P-value threshold
            logfc_threshold: Log fold change threshold
            top_n: Number of top markers to extract
        """
        self.method = method
        self.pval_threshold = pval_threshold
        self.logfc_threshold = logfc_threshold
        self.top_n = top_n
    
    def analyze(
        self,
        adata,
        celltype: str,
        sample_col: str,
        subtype_col: Optional[str] = None
    ) -> Dict[str, List[str]]:
        """
        Analyze marker genes for each subtype
        
        Parameters:
            adata: AnnData object
            celltype: Cell type name
            sample_col: Column name for sample ID
            subtype_col: Column name for subtype (default: f"{celltype}_subtype")
        
        Returns:
            Dictionary of marker genes for each subtype
        """
        subtype_col = subtype_col or f"{celltype}_subtype"
        
        if subtype_col not in adata.obs:
            warnings.warn(f"Subtype column '{subtype_col}' not found")
            return {}
        
        all_markers = {}
        
        samples = adata.obs[sample_col].unique()
        
        for sample in samples:
            sample_mask = adata.obs[sample_col] == sample
            subtype_mask = adata.obs[subtype_col].notna()
            
            adata_subtype = adata[sample_mask & subtype_mask].copy()
            
            if len(adata_subtype) == 0:
                continue
            
            subtypes = adata_subtype.obs[subtype_col].unique()
            
            if len(subtypes) < 2:
                continue
            
            try:
                if hasattr(adata_subtype.X, "toarray"):
                    estimated_mb = _estimate_memory_usage(adata_subtype.n_obs, adata_subtype.n_vars)
                    if estimated_mb > MEMORY_WARNING_THRESHOLD_MB:
                        warnings.warn(
                            f"Large matrix conversion: {adata_subtype.n_obs} cells x {adata_subtype.n_vars} genes "
                            f"(estimated {estimated_mb:.1f} MB). Consider filtering genes first."
                        )
                    adata_subtype.X = adata_subtype.X.toarray()
                
                sc.pp.normalize_total(adata_subtype, target_sum=1e4)
                sc.pp.log1p(adata_subtype)
                
                sc.tl.rank_genes_groups(
                    adata_subtype,
                    groupby=subtype_col,
                    method=self.method,
                    key_added=f"markers_{sample}"
                )
                
                for subtype in subtypes:
                    try:
                        markers_df = sc.get.rank_genes_groups_df(
                            adata_subtype,
                            group=subtype,
                            key=f"markers_{sample}"
                        )
                        
                        sig_markers = markers_df[
                            (markers_df["pvals_adj"] < self.pval_threshold) &
                            (markers_df["logfoldchanges"].abs() > self.logfc_threshold)
                        ]
                        
                        marker_genes = sig_markers["names"].tolist()[:self.top_n]
                        
                        all_markers[subtype] = marker_genes
                        
                    except Exception as e:
                        warnings.warn(f"Failed to get markers for subtype {subtype}: {e}")
                        continue
                
            except Exception as e:
                warnings.warn(f"Marker analysis failed for sample {sample}: {e}")
                continue
        
        return all_markers
    
    def analyze_global(
        self,
        adata,
        celltype: str,
        subtype_col: Optional[str] = None
    ) -> Dict[str, List[str]]:
        """
        Analyze marker genes globally (all samples combined)
        
        Parameters:
            adata: AnnData object
            celltype: Cell type name
            subtype_col: Column name for subtype
        
        Returns:
            Dictionary of marker genes for each subtype
        """
        subtype_col = subtype_col or f"{celltype}_subtype"
        
        if subtype_col not in adata.obs:
            warnings.warn(f"Subtype column '{subtype_col}' not found")
            return {}
        
        subtype_mask = adata.obs[subtype_col].notna()
        adata_subtype = adata[subtype_mask].copy()
        
        if len(adata_subtype) == 0:
            return {}
        
        subtypes = adata_subtype.obs[subtype_col].unique()
        
        if len(subtypes) < 2:
            return {}
        
        all_markers = {}
        
        try:
            if hasattr(adata_subtype.X, "toarray"):
                adata_subtype.X = adata_subtype.X.toarray()
            
            sc.pp.normalize_total(adata_subtype, target_sum=1e4)
            sc.pp.log1p(adata_subtype)
            
            sc.tl.rank_genes_groups(
                adata_subtype,
                groupby=subtype_col,
                method=self.method,
                key_added="markers_global"
            )
            
            for subtype in subtypes:
                try:
                    markers_df = sc.get.rank_genes_groups_df(
                        adata_subtype,
                        group=subtype,
                        key="markers_global"
                    )
                    
                    sig_markers = markers_df[
                        (markers_df["pvals_adj"] < self.pval_threshold) &
                        (markers_df["logfoldchanges"].abs() > self.logfc_threshold)
                    ]
                    
                    marker_genes = sig_markers["names"].tolist()[:self.top_n]
                    all_markers[subtype] = marker_genes
                    
                except Exception as e:
                    warnings.warn(f"Failed to get markers for subtype {subtype}: {e}")
                    continue
        
        except Exception as e:
            warnings.warn(f"Global marker analysis failed: {e}")
        
        return all_markers
    
    def get_marker_dataframe(
        self,
        markers_dict: Dict[str, List[str]]
    ) -> pd.DataFrame:
        """
        Convert marker genes dictionary to DataFrame
        
        Parameters:
            markers_dict: Dictionary of marker genes
        
        Returns:
            DataFrame of marker genes
        """
        rows = []
        
        for subtype, genes in markers_dict.items():
            for i, gene in enumerate(genes):
                rows.append({
                    "subtype": subtype,
                    "rank": i + 1,
                    "gene": gene
                })
        
        return pd.DataFrame(rows)
    
    def get_common_markers(
        self,
        markers_dict: Dict[str, List[str]],
        min_occurrence: int = 2
    ) -> List[str]:
        """
        Get marker genes common across multiple subtypes
        
        Parameters:
            markers_dict: Dictionary of marker genes
            min_occurrence: Minimum number of subtypes a gene must appear in
        
        Returns:
            List of common marker genes
        """
        gene_counts = {}
        
        for subtype, genes in markers_dict.items():
            for gene in genes:
                gene_counts[gene] = gene_counts.get(gene, 0) + 1
        
        common_genes = [
            gene for gene, count in gene_counts.items()
            if count >= min_occurrence
        ]
        
        return common_genes
    
    def get_unique_markers(
        self,
        markers_dict: Dict[str, List[str]]
    ) -> Dict[str, List[str]]:
        """
        Get subtype-specific unique marker genes
        
        Parameters:
            markers_dict: Dictionary of marker genes
        
        Returns:
            Dictionary of unique marker genes for each subtype
        """
        all_genes = set()
        for genes in markers_dict.values():
            all_genes.update(genes)
        
        unique_markers = {}
        
        for subtype, genes in markers_dict.items():
            other_genes = set()
            for other_subtype, other_genes_list in markers_dict.items():
                if other_subtype != subtype:
                    other_genes.update(other_genes_list)
            
            unique = [g for g in genes if g not in other_genes]
            unique_markers[subtype] = unique
        
        return unique_markers