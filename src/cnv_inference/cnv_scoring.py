"""
SpaMFC CNV Scoring Module

Provides CNV score calculation and analysis functions.
Features:
1. CNV score computation
2. CNV correlation analysis
3. CNV summary statistics
"""

import numpy as np
import pandas as pd
from typing import Optional, List, Dict, Any
import warnings
from scipy import stats

try:
    import infercnvpy as cnv
    INFERCNVPY_AVAILABLE = True
except ImportError:
    INFERCNVPY_AVAILABLE = False
    warnings.warn("infercnvpy not installed")


def compute_cnv_score(
    adata,
    method: str = "default",
    inplace: bool = True,
    key_added: str = "cnv_score"
):
    if not INFERCNVPY_AVAILABLE:
        raise ImportError("infercnvpy is required")
    
    cnv.tl.cnv_score(adata, inplace=inplace, key_added=key_added)
    
    return adata


def compute_cnv_correlation(
    adata,
    gene: str,
    cnv_key: str = "X_cnv"
) -> Dict[str, float]:
    if cnv_key not in adata.obsm:
        raise ValueError(f"CNV matrix not found in adata.obsm['{cnv_key}']")
    
    if gene not in adata.var_names:
        raise ValueError(f"Gene '{gene}' not found in adata.var_names")
    
    cnv_matrix = adata.obsm[cnv_key]
    gene_idx = adata.var_names.get_loc(gene)
    gene_expr = adata[:, gene].X
    
    if hasattr(gene_expr, 'toarray'):
        gene_expr = gene_expr.toarray().flatten()
    else:
        gene_expr = gene_expr.flatten()
    
    correlations = []
    for i in range(cnv_matrix.shape[1]):
        cnv_bin = cnv_matrix[:, i]
        
        valid_mask = ~(np.isnan(cnv_bin) | np.isnan(gene_expr))
        if valid_mask.sum() < 10:
            correlations.append(0)
            continue
        
        corr, _ = stats.spearmanr(cnv_bin[valid_mask], gene_expr[valid_mask])
        correlations.append(corr if not np.isnan(corr) else 0)
    
    return {
        'gene': gene,
        'max_correlation': max(correlations),
        'min_correlation': min(correlations),
        'mean_correlation': np.mean(correlations),
        'correlations': correlations
    }


def get_cnv_summary(adata) -> Dict[str, Any]:
    summary = {}
    
    if "X_cnv" in adata.obsm:
        cnv_matrix = adata.obsm["X_cnv"]
        
        summary["cnv_matrix_shape"] = cnv_matrix.shape
        summary["cnv_mean"] = float(np.mean(cnv_matrix))
        summary["cnv_std"] = float(np.std(cnv_matrix))
        summary["cnv_min"] = float(np.min(cnv_matrix))
        summary["cnv_max"] = float(np.max(cnv_matrix))
        
        amplification_threshold = 0.1
        deletion_threshold = -0.1
        
        amplification_pct = (cnv_matrix > amplification_threshold).sum() / cnv_matrix.size * 100
        deletion_pct = (cnv_matrix < deletion_threshold).sum() / cnv_matrix.size * 100
        
        summary["amplification_pct"] = float(amplification_pct)
        summary["deletion_pct"] = float(deletion_pct)
    
    if "cnv_score" in adata.obs.columns:
        cnv_scores = adata.obs["cnv_score"]
        
        summary["cnv_score_mean"] = float(cnv_scores.mean())
        summary["cnv_score_std"] = float(cnv_scores.std())
        summary["cnv_score_min"] = float(cnv_scores.min())
        summary["cnv_score_max"] = float(cnv_scores.max())
        
        high_cnv_threshold = cnv_scores.quantile(0.75)
        summary["high_cnv_cells"] = int((cnv_scores > high_cnv_threshold).sum())
    
    if "cnv_leiden" in adata.obs.columns:
        clusters = adata.obs["cnv_leiden"]
        
        summary["n_cnv_clusters"] = int(clusters.nunique())
        summary["cluster_sizes"] = clusters.value_counts().to_dict()
        
        if "cnv_score" in adata.obs.columns:
            cluster_scores = adata.obs.groupby("cnv_leiden")["cnv_score"].mean()
            summary["cluster_mean_scores"] = cluster_scores.to_dict()
    
    return summary


def classify_cells_by_cnv(
    adata,
    score_key: str = "cnv_score",
    threshold_high: Optional[float] = None,
    threshold_low: Optional[float] = None,
    key_added: str = "cnv_classification"
) -> pd.Series:
    if score_key not in adata.obs.columns:
        raise ValueError(f"CNV score not found in adata.obs['{score_key}']")
    
    scores = adata.obs[score_key]
    
    if threshold_high is None:
        threshold_high = scores.quantile(0.75)
    
    if threshold_low is None:
        threshold_low = scores.quantile(0.25)
    
    classification = pd.Series('normal', index=adata.obs_names)
    
    classification[scores > threshold_high] = 'high_cnv'
    classification[scores < threshold_low] = 'low_cnv'
    
    adata.obs[key_added] = classification
    
    return classification


def get_chromosome_cnv_stats(
    adata,
    cnv_key: str = "X_cnv"
) -> pd.DataFrame:
    if cnv_key not in adata.obsm:
        raise ValueError(f"CNV matrix not found in adata.obsm['{cnv_key}']")
    
    if "cnv" not in adata.uns or "chr_pos" not in adata.uns["cnv"]:
        warnings.warn("Chromosome positions not found in adata.uns")
        return pd.DataFrame()
    
    cnv_matrix = adata.obsm[cnv_key]
    chr_pos = adata.uns["cnv"]["chr_pos"]
    
    chr_stats = []
    
    for chrom, (start, end) in chr_pos.items():
        chr_cnv = cnv_matrix[:, start:end]
        
        chr_stats.append({
            'chromosome': chrom,
            'mean_cnv': float(np.mean(chr_cnv)),
            'std_cnv': float(np.std(chr_cnv)),
            'amplification_pct': float((chr_cnv > 0.1).sum() / chr_cnv.size * 100),
            'deletion_pct': float((chr_cnv < -0.1).sum() / chr_cnv.size * 100),
            'n_bins': end - start
        })
    
    return pd.DataFrame(chr_stats)