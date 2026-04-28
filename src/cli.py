"""
SpaMFC Command Line Interface Module (Optimized)

Provides complete command line interface for SpaMFC analysis.
Simplified feature selection with single --features parameter.

Usage:
    SpaMFC run --input data.h5ad --celltype-col "anno_cell2location_res" --celltype "CAFs" --features spatial,expression
    SpaMFC run-multi --input data.h5ad --celltype-col "anno_cell2location_res" --celltypes "Malignant cells,CAFs,ILC"
    SpaMFC info --input data.h5ad --celltype-col "anno_cell2location_res"
    SpaMFC config --output ./my_config.yaml
    SpaMFC corr --input data.h5ad --target-genes EGFR,KRAS --de-genes GENE1,GENE2 --method spearman
"""

import argparse
import sys
import os
import logging
from pathlib import Path
from typing import Optional, List, Dict
import warnings

try:
    import scanpy as sc
except ImportError:
    warnings.warn("scanpy not installed")

from src.utils.logger import setup_logger, log_section


VALID_FEATURES = ["spatial", "cnv", "expression", "niche"]
VALID_FUSION_METHODS = ["adaptive", "fixed", "concat"]

cli_logger = setup_logger("SpaMFC.CLI", save_log=False)


def create_parser() -> argparse.ArgumentParser:
    """Create main argument parser with all subcommands"""
    
    parser = argparse.ArgumentParser(
        prog="SpaMFC",
        description="SpaMFC: Spatial Multi-Feature Clustering for Spatial Transcriptomics Subtype Analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # CAFs analysis with spatial+expression features
  SpaMFC run --input data.h5ad --celltype-col "anno_cell2location_res" --celltype "CAFs" --features spatial,expression
  
  # Malignant cells with all features and fixed weights
  SpaMFC run --input data.h5ad --celltype-col "anno_cell2location_res" --celltype "Malignant cells" --features spatial,cnv,expression --fusion-method fixed --weights 0.3,0.3,0.4
  
  # Multi-cell type analysis
  SpaMFC run-multi --input data.h5ad --celltype-col "anno_cell2location_res" --celltypes "Malignant cells,CAFs,ILC" --features spatial,expression
  
  # Show data info
  SpaMFC info --input data.h5ad --celltype-col "anno_cell2location_res"
  
  # Generate config file
  SpaMFC config --output ./my_config.yaml --template cafs
  
  # Gene correlation analysis
  SpaMFC corr --input data.h5ad --target-genes EGFR,KRAS --de-genes GENE1,GENE2 --method spearman
"""
    )
    
    parser.add_argument(
        "--version", "-v",
        action="version",
        version="SpaMFC v2.2.0"
    )
    
    subparsers = parser.add_subparsers(
        dest="command",
        title="Available Commands",
        description="Use one of the following commands:"
    )
    
    _add_run_subparser(subparsers)
    _add_run_multi_subparser(subparsers)
    _add_info_subparser(subparsers)
    _add_config_subparser(subparsers)
    _add_corr_subparser(subparsers)
    _add_cnv_subparser(subparsers)
    
    return parser


def _add_run_subparser(subparsers):
    """Add 'run' subcommand for single cell type analysis"""
    
    run_parser = subparsers.add_parser(
        "run",
        help="Run subtype analysis for a single cell type",
        description="Execute SpaMFC analysis for a single cell type with customizable parameters."
    )
    
    _add_data_args(run_parser)
    _add_feature_args(run_parser)
    _add_fusion_args(run_parser)
    _add_analysis_args(run_parser)
    _add_annotation_args(run_parser)
    _add_visualization_args(run_parser)
    _add_runtime_args(run_parser)
    
    run_parser.add_argument(
        "--celltype", "-c",
        type=str,
        required=True,
        help="Target cell type name to analyze (e.g., 'Malignant cells', 'CAFs', 'ILC')"
    )
    
    run_parser.add_argument(
        "--markers",
        type=str,
        nargs="+",
        default=None,
        help="Tumor marker genes for CNV filtering (e.g., EPCAM KRT8 KRT18)"
    )


def _add_run_multi_subparser(subparsers):
    """Add 'run-multi' subcommand for multiple cell types analysis"""
    
    multi_parser = subparsers.add_parser(
        "run-multi",
        help="Run subtype analysis for multiple cell types",
        description="Execute SpaMFC analysis for multiple cell types in batch mode."
    )
    
    _add_data_args(multi_parser)
    _add_feature_args(multi_parser)
    _add_fusion_args(multi_parser)
    _add_runtime_args(multi_parser)
    
    multi_parser.add_argument(
        "--celltypes",
        type=str,
        required=True,
        help="Comma-separated list of cell types (e.g., 'Malignant cells,CAFs,ILC')"
    )


def _add_info_subparser(subparsers):
    """Add 'info' subcommand for displaying data information"""
    
    info_parser = subparsers.add_parser(
        "info",
        help="Display data information",
        description="Show basic information about the input h5ad file."
    )
    
    info_parser.add_argument(
        "--input", "-i",
        type=str,
        required=True,
        help="Input h5ad file path"
    )
    
    info_parser.add_argument(
        "--celltype-col",
        type=str,
        required=True,
        help="Column name for cell type annotation (required to show cell types)"
    )
    
    info_parser.add_argument(
        "--sample-col",
        type=str,
        default="sample_id",
        help="Column name for sample ID (default: sample_id)"
    )


def _add_config_subparser(subparsers):
    """Add 'config' subcommand for generating default config file"""
    
    config_parser = subparsers.add_parser(
        "config",
        help="Generate default configuration file",
        description="Generate a default YAML configuration file for SpaMFC."
    )
    
    config_parser.add_argument(
        "--output", "-o",
        type=str,
        default="./spamfc_config.yaml",
        help="Output config file path (default: ./spamfc_config.yaml)"
    )
    
    config_parser.add_argument(
        "--template",
        type=str,
        choices=["default", "malignant", "cafs", "immune"],
        default="default",
        help="Config template to use (default: default)"
    )


def _add_corr_subparser(subparsers):
    """Add 'corr' subcommand for gene correlation analysis"""
    
    corr_parser = subparsers.add_parser(
        "corr",
        help="Gene correlation analysis",
        description="Perform gene correlation analysis on spatial transcriptomics data."
    )
    
    corr_parser.add_argument(
        "--input", "-i",
        type=str,
        required=True,
        help="Input h5ad file path"
    )
    
    corr_parser.add_argument(
        "--target-genes",
        type=str,
        required=True,
        help="Target genes (comma-separated or file path)"
    )
    
    corr_parser.add_argument(
        "--de-genes",
        type=str,
        required=True,
        help="DE genes (comma-separated or file path)"
    )
    
    corr_parser.add_argument(
        "--method",
        type=str,
        choices=["pearson", "spearman", "kendall"],
        default="spearman",
        help="Correlation method (default: spearman)"
    )
    
    corr_parser.add_argument(
        "--p-adjust",
        type=str,
        choices=["fdr_bh", "bonferroni", "holm", "none"],
        default="fdr_bh",
        help="P-value adjustment method (default: fdr_bh)"
    )
    
    corr_parser.add_argument(
        "--threshold-p",
        type=float,
        default=0.05,
        help="P-value threshold (default: 0.05)"
    )
    
    corr_parser.add_argument(
        "--min-corr",
        type=float,
        default=0.0,
        help="Minimum correlation threshold (default: 0.0)"
    )
    
    corr_parser.add_argument(
        "--output", "-o",
        type=str,
        default="./correlation_results/",
        help="Output directory (default: ./correlation_results/)"
    )
    
    corr_parser.add_argument(
        "--save-matrices",
        action="store_true",
        default=False,
        help="Save full correlation matrices"
    )
    
    corr_parser.add_argument(
        "--matrix-format",
        type=str,
        choices=["npz", "csv", "csv.gz"],
        default="npz",
        help="Matrix format (default: npz)"
    )
    
    corr_parser.add_argument(
        "--n-workers",
        type=int,
        default=None,
        help="Number of parallel workers (default: auto)"
    )
    
    corr_parser.add_argument(
        "--max-memory",
        type=int,
        default=512,
        help="Maximum memory limit in MB (default: 512)"
    )
    
    corr_parser.add_argument(
        "--batch-size",
        type=int,
        default=500,
        help="Batch size for processing (default: 500)"
    )
    
    corr_parser.add_argument(
        "--sample-spots",
        type=int,
        default=None,
        help="Sample spots for acceleration (default: None)"
    )
    
    corr_parser.add_argument(
        "--verbose",
        action="store_true",
        default=True,
        help="Verbose output (default: True)"
    )
    
    corr_parser.add_argument(
        "--quiet",
        action="store_true",
        default=False,
        help="Quiet mode"
    )


def _add_cnv_subparser(subparsers):
    """Add 'cnv' subcommand for CNV inference"""
    
    cnv_parser = subparsers.add_parser(
        "cnv",
        help="CNV inference and visualization",
        description="Perform CNV inference using inferCNVpy for malignant cell analysis."
    )
    
    cnv_parser.add_argument(
        "--input", "-i",
        type=str,
        required=True,
        help="Input h5ad file path"
    )
    
    cnv_parser.add_argument(
        "--gtf",
        type=str,
        default=None,
        help="Gene annotation file (GTF or TSV format: gene_symbol, chromosome, start, end). Required for infercnv method."
    )
    
    cnv_parser.add_argument(
        "--method",
        type=str,
        default="infercnv",
        choices=["infercnv", "copykat"],
        help="CNV inference method (default: infercnv). copykat auto-infers gene positions."
    )
    
    cnv_parser.add_argument(
        "--species",
        type=str,
        default="human",
        choices=["human", "mouse"],
        help="Species for gene annotation (default: human)"
    )
    
    cnv_parser.add_argument(
        "--reference-key",
        type=str,
        required=True,
        help="Reference key in adata.obs (e.g., 'cell_type')"
    )
    
    cnv_parser.add_argument(
        "--reference-cat",
        type=str,
        required=True,
        help="Reference categories (comma-separated, e.g., 'Normal cells,Fibroblasts')"
    )
    
    cnv_parser.add_argument(
        "--output", "-o",
        type=str,
        default="./cnv_results/",
        help="Output directory (default: ./cnv_results/)"
    )
    
    cnv_parser.add_argument(
        "--window-size",
        type=int,
        default=100,
        help="Window size for CNV smoothing (default: 100)"
    )
    
    cnv_parser.add_argument(
        "--step-size",
        type=int,
        default=10,
        help="Step size for CNV smoothing (default: 10)"
    )
    
    cnv_parser.add_argument(
        "--resolution",
        type=float,
        default=0.5,
        help="Leiden clustering resolution (default: 0.5)"
    )
    
    cnv_parser.add_argument(
        "--n-pcs",
        type=int,
        default=30,
        help="Number of PCA components (default: 30)"
    )
    
    cnv_parser.add_argument(
        "--plot-heatmap",
        action="store_true",
        default=True,
        help="Generate chromosome heatmap (default: True)"
    )
    
    cnv_parser.add_argument(
        "--plot-umap",
        action="store_true",
        default=True,
        help="Generate CNV UMAP (default: True)"
    )
    
    cnv_parser.add_argument(
        "--plot-spatial",
        action="store_true",
        default=False,
        help="Generate CNV spatial distribution plots"
    )
    
    cnv_parser.add_argument(
        "--sample-col",
        type=str,
        default="sample_id",
        help="Sample column name for spatial plots (default: sample_id)"
    )
    
    cnv_parser.add_argument(
        "--spatial-key",
        type=str,
        default="spatial",
        help="Spatial coordinates key in adata.obsm (default: spatial)"
    )
    
    cnv_parser.add_argument(
        "--cmap",
        type=str,
        default="bwr",
        help="Colormap for CNV heatmap (default: bwr)"
    )
    
    cnv_parser.add_argument(
        "--groupby",
        type=str,
        default=None,
        help="Column to group cells for heatmap (default: use reference-key column)"
    )
    
    cnv_parser.add_argument(
        "--verbose",
        action="store_true",
        default=True,
        help="Verbose output (default: True)"
    )
    
    cnv_parser.add_argument(
        "--quiet",
        action="store_true",
        default=False,
        help="Quiet mode"
    )


def _add_data_args(parser):
    """Add data input/output arguments"""
    
    data_group = parser.add_argument_group("Data Parameters")
    
    data_group.add_argument(
        "--input", "-i",
        type=str,
        required=True,
        help="Input h5ad file path"
    )
    
    data_group.add_argument(
        "--output", "-o",
        type=str,
        default="./results/",
        help="Output directory (default: ./results/)"
    )
    
    data_group.add_argument(
        "--celltype-col",
        type=str,
        required=True,
        help="Column name for cell type annotation (required, e.g., 'anno_cell2location_res')"
    )
    
    data_group.add_argument(
        "--sample-col",
        type=str,
        default="sample_id",
        help="Column name for sample ID (default: sample_id)"
    )
    
    data_group.add_argument(
        "--spatial-key",
        type=str,
        default="spatial",
        help="Key for spatial coordinates in adata.obsm (default: spatial)"
    )
    
    data_group.add_argument(
        "--cnv-key",
        type=str,
        default="X_cnv",
        help="Key for CNV matrix in adata.obsm (default: X_cnv)"
    )


def _add_feature_args(parser):
    """Add feature selection arguments - simplified single parameter"""
    
    feature_group = parser.add_argument_group("Feature Selection")
    
    feature_group.add_argument(
        "--features",
        type=str,
        default="spatial",
        help="Features to use, comma-separated (default: spatial). "
             "Valid options: spatial, cnv, expression, niche. "
             "Examples: 'spatial', 'spatial,expression', 'spatial,cnv,expression'"
    )
    
    feature_group.add_argument(
        "--radius",
        type=int,
        default=100,
        help="Spatial neighborhood radius (default: 100)"
    )
    
    feature_group.add_argument(
        "--k-neighbors",
        type=int,
        default=15,
        help="Number of neighbors for spatial features (default: 15)"
    )


def _add_fusion_args(parser):
    """Add feature fusion method arguments"""
    
    fusion_group = parser.add_argument_group("Feature Fusion")
    
    fusion_group.add_argument(
        "--fusion-method",
        type=str,
        choices=VALID_FUSION_METHODS,
        default="adaptive",
        help="Feature fusion method (default: adaptive). "
             "adaptive: adaptive weighting based on feature importance; "
             "fixed: fixed weights from --weights parameter; "
             "concat: direct concatenation without weighting"
    )
    
    fusion_group.add_argument(
        "--weights",
        type=str,
        default=None,
        help="Fixed weights for features, comma-separated, must match --features order (used when fusion-method=fixed). "
             "Example: --features spatial,cnv,expression --weights 0.3,0.3,0.4"
    )


def _add_analysis_args(parser):
    """Add analysis parameters"""
    
    analysis_group = parser.add_argument_group("Analysis Parameters")
    
    analysis_group.add_argument(
        "--n-clusters",
        type=int,
        default=5,
        help="Number of clusters for KMeans (default: 5)"
    )
    
    analysis_group.add_argument(
        "--nmf-components",
        type=int,
        default=5,
        help="Number of NMF factors (default: 5)"
    )
    
    analysis_group.add_argument(
        "--nmf-runs",
        type=int,
        default=10,
        help="Number of NMF runs for consensus (default: 10)"
    )
    
    analysis_group.add_argument(
        "--clustering-method",
        type=str,
        choices=["kmeans", "leiden"],
        default="kmeans",
        help="Clustering method (default: kmeans)"
    )
    
    analysis_group.add_argument(
        "--resolution",
        type=float,
        default=0.1,
        help="Resolution for Leiden clustering (default: 0.1)"
    )
    
    analysis_group.add_argument(
        "--per-sample",
        action="store_true",
        default=True,
        help="Cluster per sample to avoid batch effects (default: True)"
    )
    
    analysis_group.add_argument(
        "--global-cluster",
        action="store_true",
        default=False,
        help="Cluster globally (may have batch effects)"
    )


def _add_annotation_args(parser):
    """Add annotation parameters"""
    
    annotation_group = parser.add_argument_group("Annotation Parameters")
    
    annotation_group.add_argument(
        "--enable-marker-analysis",
        action="store_true",
        default=True,
        help="Enable marker gene analysis (default: True)"
    )
    
    annotation_group.add_argument(
        "--disable-marker-analysis",
        action="store_true",
        default=False,
        help="Disable marker gene analysis"
    )
    
    annotation_group.add_argument(
        "--enable-enrichment",
        action="store_true",
        default=True,
        help="Enable functional enrichment analysis (default: True)"
    )
    
    annotation_group.add_argument(
        "--disable-enrichment",
        action="store_true",
        default=False,
        help="Disable functional enrichment analysis"
    )
    
    annotation_group.add_argument(
        "--enable-unification",
        action="store_true",
        default=True,
        help="Enable cross-sample subtype unification (default: True)"
    )
    
    annotation_group.add_argument(
        "--disable-unification",
        action="store_true",
        default=False,
        help="Disable cross-sample subtype unification"
    )
    
    annotation_group.add_argument(
        "--marker-top-n",
        type=int,
        default=50,
        help="Number of top marker genes to extract (default: 50)"
    )
    
    annotation_group.add_argument(
        "--pval-threshold",
        type=float,
        default=0.05,
        help="P-value threshold for significance (default: 0.05)"
    )


def _add_visualization_args(parser):
    """Add visualization parameters"""
    
    viz_group = parser.add_argument_group("Visualization Parameters")
    
    viz_group.add_argument(
        "--save-plots",
        action="store_true",
        default=True,
        help="Save visualization plots (default: True)"
    )
    
    viz_group.add_argument(
        "--no-plots",
        action="store_true",
        default=False,
        help="Do not save visualization plots"
    )
    
    viz_group.add_argument(
        "--plot-format",
        type=str,
        choices=["pdf", "png", "svg"],
        default="pdf",
        help="Plot format (default: pdf)"
    )
    
    viz_group.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Plot resolution DPI (default: 300)"
    )


def _add_runtime_args(parser):
    """Add runtime parameters"""
    
    runtime_group = parser.add_argument_group("Runtime Parameters")
    
    runtime_group.add_argument(
        "--config",
        type=str,
        default=None,
        help="Path to YAML configuration file"
    )
    
    runtime_group.add_argument(
        "--verbose",
        action="store_true",
        default=True,
        help="Print detailed progress information (default: True)"
    )
    
    runtime_group.add_argument(
        "--quiet",
        action="store_true",
        default=False,
        help="Suppress detailed output"
    )
    
    runtime_group.add_argument(
        "--n-jobs",
        type=int,
        default=1,
        help="Number of parallel jobs (default: 1)"
    )


def parse_features(features_str: str) -> Dict[str, bool]:
    """Parse features string into feature selection dict"""
    
    features_list = [f.strip().lower() for f in features_str.split(",")]
    
    features = {
        "use_spatial": False,
        "use_cnv": False,
        "use_expression": False,
        "use_niche": False
    }
    
    for feature in features_list:
        if feature in VALID_FEATURES:
            features[f"use_{feature}"] = True
        else:
            cli_logger.warning(f"Unknown feature '{feature}', ignoring. Valid features: {VALID_FEATURES}")
    
    if not any(features.values()):
        cli_logger.warning(f"No valid features specified, using default 'spatial'")
        features["use_spatial"] = True
    
    return features


def parse_weights(weights_str: str, features_str: str) -> Dict[str, float]:
    """Parse weights string into weight dict, matching features order"""
    
    features_list = [f.strip().lower() for f in features_str.split(",")]
    features_list = [f for f in features_list if f in VALID_FEATURES]
    
    if weights_str is None:
        default_weight = 1.0 / len(features_list) if features_list else 0.25
        weights = {}
        for feature in features_list:
            weights[feature] = default_weight
        return weights
    
    try:
        weights_list = [float(w.strip()) for w in weights_str.split(",")]
    except ValueError as e:
        raise ValueError(f"Invalid weight value in '{weights_str}'. Weight values must be numbers: {e}")
    
    for w in weights_list:
        if w < 0:
            raise ValueError(f"Weight values must be non-negative: found {w}")
    
    if len(weights_list) != len(features_list):
        cli_logger.warning(f"Number of weights ({len(weights_list)}) does not match number of features ({len(features_list)})")
        if len(weights_list) < len(features_list):
            remaining = len(features_list) - len(weights_list)
            default_weight = 1.0 / len(features_list)
            weights_list.extend([default_weight] * remaining)
        else:
            weights_list = weights_list[:len(features_list)]
    
    weights = {}
    for i, feature in enumerate(features_list):
        weights[feature] = weights_list[i]
    
    total = sum(weights.values())
    if total > 0:
        weights = {k: v / total for k, v in weights.items()}
    else:
        raise ValueError("Total weight cannot be zero")
    
    return weights


def validate_features_for_celltype(celltype: str, features: Dict[str, bool]) -> Dict[str, bool]:
    """Validate and adjust features based on cell type"""
    
    malignant_keywords = ["malignant", "tumor", "cancer", "恶性", "肿瘤"]
    is_malignant = any(kw in celltype.lower() for kw in malignant_keywords)
    
    if features["use_cnv"] and not is_malignant:
        cli_logger.warning(f"CNV features are only recommended for malignant cells. "
              f"Cell type '{celltype}' may not be malignant. CNV will still be used if data exists.")
    
    return features


def parse_gene_list(genes_str: str) -> List[str]:
    """Parse gene list from comma-separated string or file path"""
    
    if os.path.isfile(genes_str):
        file_path = os.path.abspath(genes_str)
        file_ext = os.path.splitext(file_path)[1].lower()
        valid_extensions = ['.txt', '.csv', '.tsv', '.list', '.genes']
        if file_ext not in valid_extensions and file_ext != '':
            warnings.warn(f"Unexpected file extension '{file_ext}' for gene list file")
        try:
            with open(genes_str, 'r', encoding='utf-8') as f:
                genes = [line.strip() for line in f if line.strip() and not line.startswith('#')]
        except (IOError, OSError, PermissionError) as e:
            raise ValueError(f"Failed to read gene file '{genes_str}': {e}")
    else:
        genes = [g.strip() for g in genes_str.split(",")]
    
    if not genes:
        raise ValueError("No genes found in input")
    
    return genes


def corr_command(args):
    """Execute 'corr' subcommand"""
    
    try:
        import scanpy as sc
    except ImportError:
        cli_logger.error("scanpy is not installed. Please install it first.")
        sys.exit(1)
    
    from src.correlation import gene_correlation
    
    target_genes = parse_gene_list(args.target_genes)
    de_genes = parse_gene_list(args.de_genes)
    
    log_section(cli_logger, "SpaMFC Gene Correlation Analysis")
    cli_logger.info(f"Input: {args.input}")
    cli_logger.info(f"Output: {args.output}")
    cli_logger.info(f"Target genes: {len(target_genes)} genes")
    cli_logger.info(f"DE genes: {len(de_genes)} genes")
    cli_logger.info(f"Method: {args.method}")
    cli_logger.info(f"P-value adjustment: {args.p_adjust}")
    cli_logger.info(f"P-value threshold: {args.threshold_p}")
    cli_logger.info(f"Min correlation: {args.min_corr}")
    
    adata = sc.read_h5ad(args.input)
    
    p_adjust = None if args.p_adjust == "none" else args.p_adjust
    
    corr_df, pval_df, sig_pairs = gene_correlation(
        adata,
        target_genes=target_genes,
        de_genes=de_genes,
        method=args.method,
        p_adjust=p_adjust,
        threshold_p=args.threshold_p,
        min_corr_threshold=args.min_corr,
        output_dir=args.output,
        max_memory_mb=args.max_memory,
        n_workers=args.n_workers,
        batch_size=args.batch_size,
        sample_spots=args.sample_spots,
        save_full_matrices=args.save_matrices,
        matrix_format=args.matrix_format,
        verbose=args.verbose and not args.quiet
    )
    
    log_section(cli_logger, f"Found {len(sig_pairs)} significant correlation pairs")
    cli_logger.info(f"Results saved to: {args.output}")


def run_command(args):
    """Execute 'run' subcommand"""
    
    try:
        import scanpy as sc
    except ImportError:
        cli_logger.error("scanpy is not installed. Please install it first.")
        sys.exit(1)
    
    from src.pipeline import SpaMFCPipeline
    from src.config import ConfigManager
    
    if args.config:
        pipeline = SpaMFCPipeline(args.config)
    else:
        pipeline = SpaMFCPipeline()
    
    pipeline.config.output_dir = args.output
    pipeline.config.celltype_col = args.celltype_col
    pipeline.config.sample_col = args.sample_col
    pipeline.config.spatial_key = args.spatial_key
    pipeline.config.cnv_key = args.cnv_key
    
    features = parse_features(args.features)
    features = validate_features_for_celltype(args.celltype, features)
    pipeline.set_feature_usage(**features)
    
    pipeline.config.spatial_config.radius = args.radius
    pipeline.config.spatial_config.k_neighbors = args.k_neighbors
    
    pipeline.config.clustering_config.n_clusters = args.n_clusters
    pipeline.config.nmf_config.n_components = args.nmf_components
    pipeline.config.nmf_config.n_runs = args.nmf_runs
    pipeline.config.clustering_config.method = args.clustering_method
    pipeline.config.clustering_config.resolution = args.resolution
    pipeline.config.clustering_config.per_sample = args.per_sample and not args.global_cluster
    
    pipeline.config.marker_config.top_n = args.marker_top_n
    pipeline.config.marker_config.pval_threshold = args.pval_threshold
    
    pipeline.config.unification_config.enable = args.enable_unification and not args.disable_unification
    
    pipeline.config.visualization_config.save_plots = args.save_plots and not args.no_plots
    pipeline.config.visualization_config.format = args.plot_format
    pipeline.config.visualization_config.dpi = args.dpi
    
    pipeline.config.runtime_config.verbose = args.verbose and not args.quiet
    pipeline.config.runtime_config.n_jobs = args.n_jobs
    
    pipeline.config.fusion_method = args.fusion_method
    
    weights = parse_weights(args.weights, args.features)
    if "spatial" in weights:
        pipeline.config.spatial_config.fixed_weight = weights["spatial"]
    if "cnv" in weights:
        pipeline.config.cnv_config.fixed_weight = weights["cnv"]
    if "expression" in weights:
        pipeline.config.expression_config.fixed_weight = weights["expression"]
    if "niche" in weights:
        pipeline.config.niche_config.fixed_weight = weights["niche"]
    
    log_section(cli_logger, "SpaMFC Command Line Analysis")
    cli_logger.info(f"Input: {args.input}")
    cli_logger.info(f"Output: {args.output}")
    cli_logger.info(f"Cell type column: {args.celltype_col}")
    cli_logger.info(f"Cell type: {args.celltype}")
    cli_logger.info(f"Features: {args.features}")
    cli_logger.info(f"Fusion method: {args.fusion_method}")
    if args.weights:
        cli_logger.info(f"Fixed weights: {args.weights}")
    
    adata = sc.read_h5ad(args.input)
    
    if args.celltype_col not in adata.obs:
        cli_logger.error(f"Cell type column '{args.celltype_col}' not found in adata.obs")
        cli_logger.error(f"Available columns: {list(adata.obs.columns)}")
        sys.exit(1)
    
    if args.celltype not in adata.obs[args.celltype_col].values:
        cli_logger.error(f"Cell type '{args.celltype}' not found in column '{args.celltype_col}'")
        cli_logger.error(f"Available cell types: {list(adata.obs[args.celltype_col].unique())}")
        sys.exit(1)
    
    markers = args.markers if args.markers else None
    
    adata = pipeline.run(adata, args.celltype, markers)
    
    output_path = Path(args.output) / args.celltype.replace(" ", "_")
    output_path.mkdir(parents=True, exist_ok=True)
    
    output_file = output_path / f"{args.celltype.replace(' ', '_')}_subtype_annotated.h5ad"
    adata.write_h5ad(output_file)
    
    cli_logger.info(f"Results saved to: {output_file}")


def run_multi_command(args):
    """Execute 'run-multi' subcommand"""
    
    try:
        import scanpy as sc
    except ImportError:
        cli_logger.error("scanpy is not installed. Please install it first.")
        sys.exit(1)
    
    from src.pipeline import SpaMFCPipeline
    
    celltypes = [c.strip() for c in args.celltypes.split(",")]
    
    if args.config:
        pipeline = SpaMFCPipeline(args.config)
    else:
        pipeline = SpaMFCPipeline()
    
    pipeline.config.output_dir = args.output
    pipeline.config.celltype_col = args.celltype_col
    pipeline.config.sample_col = args.sample_col
    pipeline.config.runtime_config.verbose = args.verbose and not args.quiet
    
    features = parse_features(args.features)
    pipeline.config.fusion_method = args.fusion_method
    
    weights = parse_weights(args.weights, args.features)
    if "spatial" in weights:
        pipeline.config.spatial_config.fixed_weight = weights["spatial"]
    if "cnv" in weights:
        pipeline.config.cnv_config.fixed_weight = weights["cnv"]
    if "expression" in weights:
        pipeline.config.expression_config.fixed_weight = weights["expression"]
    if "niche" in weights:
        pipeline.config.niche_config.fixed_weight = weights["niche"]
    
    log_section(cli_logger, "SpaMFC Multi-Cell Type Analysis")
    cli_logger.info(f"Input: {args.input}")
    cli_logger.info(f"Output: {args.output}")
    cli_logger.info(f"Cell type column: {args.celltype_col}")
    cli_logger.info(f"Cell types: {celltypes}")
    cli_logger.info(f"Features: {args.features}")
    cli_logger.info(f"Fusion method: {args.fusion_method}")
    if args.weights:
        cli_logger.info(f"Fixed weights: {args.weights}")
    
    adata = sc.read_h5ad(args.input)
    
    if args.celltype_col not in adata.obs:
        cli_logger.error(f"Cell type column '{args.celltype_col}' not found in adata.obs")
        cli_logger.error(f"Available columns: {list(adata.obs.columns)}")
        sys.exit(1)
    
    available_celltypes = list(adata.obs[args.celltype_col].unique())
    missing_celltypes = [ct for ct in celltypes if ct not in available_celltypes]
    
    if missing_celltypes:
        cli_logger.warning(f"Cell types not found: {missing_celltypes}")
        cli_logger.warning(f"Available cell types: {available_celltypes}")
        celltypes = [ct for ct in celltypes if ct in available_celltypes]
        
        if not celltypes:
            cli_logger.error("No valid cell types to analyze")
            sys.exit(1)
    
    adata = pipeline.run_multi_celltypes(adata, celltypes)
    
    output_file = Path(args.output) / "all_celltypes_subtype_annotated.h5ad"
    adata.write_h5ad(output_file)
    
    cli_logger.info(f"Results saved to: {output_file}")


def info_command(args):
    """Execute 'info' subcommand"""
    
    try:
        import scanpy as sc
    except ImportError:
        cli_logger.error("scanpy is not installed. Please install it first.")
        sys.exit(1)
    
    log_section(cli_logger, f"Data Information: {args.input}")
    
    adata = sc.read_h5ad(args.input)
    
    cli_logger.info(f"Cells: {adata.n_obs}")
    cli_logger.info(f"Genes: {adata.n_vars}")
    
    if args.celltype_col in adata.obs:
        celltypes = adata.obs[args.celltype_col].value_counts()
        cli_logger.info(f"Cell Type Distribution (from '{args.celltype_col}'):")
        for ct, count in celltypes.items():
            cli_logger.info(f"  {ct}: {count} ({count/adata.n_obs*100:.1f}%)")
    else:
        cli_logger.warning(f"Cell type column '{args.celltype_col}' not found")
        cli_logger.warning(f"Available columns: {list(adata.obs.columns)}")
    
    if args.sample_col in adata.obs:
        samples = adata.obs[args.sample_col].unique()
        cli_logger.info(f"Samples: {len(samples)} ({', '.join(samples[:5])}...)")
    
    if "spatial" in adata.obsm:
        cli_logger.info(f"Spatial coordinates: Available")
    
    if "X_cnv" in adata.obsm:
        cli_logger.info(f"CNV matrix: Available")


def config_command(args):
    """Execute 'config' subcommand"""
    
    from src.config import ConfigManager
    
    template_files = {
        "default": "configs/default_config.yaml",
        "malignant": "configs/malignant_config.yaml",
        "cafs": "configs/cafs_config.yaml",
        "immune": "configs/immune_config.yaml"
    }
    
    template_path = template_files.get(args.template)
    
    if template_path and Path(template_path).exists():
        import shutil
        shutil.copy(template_path, args.output)
        cli_logger.info(f"Config file generated: {args.output}")
        cli_logger.info(f"Template: {args.template}")
    else:
        config_manager = ConfigManager()
        config_manager.save_to_yaml(args.output)
        cli_logger.info(f"Default config file generated: {args.output}")
    
    cli_logger.info(f"Edit the config file to customize your analysis parameters.")
    cli_logger.info(f"Usage: SpaMFC run --config {args.output} --input data.h5ad --celltype-col 'anno_cell2location_res' --celltype 'CAFs'")


def cnv_command(args):
    """Execute 'cnv' subcommand"""
    
    try:
        import scanpy as sc
    except ImportError:
        cli_logger.error("scanpy is not installed. Please install it first.")
        sys.exit(1)
    
    if args.method == "infercnv":
        try:
            import infercnvpy as cnv
        except ImportError:
            cli_logger.error("infercnvpy is not installed. Please install it first.")
            cli_logger.error("Install with: pip install infercnvpy")
            sys.exit(1)
    
    from src.cnv_inference import CNVInferencer, run_cnv_inference, add_gene_location
    from src.visualization import CNVVisualizer
    
    reference_cat = [c.strip() for c in args.reference_cat.split(",")]
    
    log_section(cli_logger, "SpaMFC CNV Inference Analysis")
    cli_logger.info(f"Input: {args.input}")
    cli_logger.info(f"Method: {args.method}")
    cli_logger.info(f"Species: {args.species}")
    cli_logger.info(f"Gene file: {args.gtf}")
    cli_logger.info(f"Reference key: {args.reference_key}")
    cli_logger.info(f"Reference categories: {reference_cat}")
    cli_logger.info(f"Output: {args.output}")
    cli_logger.info(f"Window size: {args.window_size}")
    cli_logger.info(f"Resolution: {args.resolution}")
    
    adata = sc.read_h5ad(args.input)
    
    if args.reference_key not in adata.obs:
        cli_logger.error(f"Reference key '{args.reference_key}' not found in adata.obs")
        cli_logger.error(f"Available columns: {list(adata.obs.columns)}")
        sys.exit(1)
    
    adata = add_gene_location(
        adata,
        gene_file=args.gtf,
        method=args.method,
        species=args.species,
        inplace=True,
        verbose=args.verbose and not args.quiet
    )
    
    if args.method == "infercnv":
        inferencer = CNVInferencer(
            window_size=args.window_size,
            step_size=args.step_size,
            clustering_resolution=args.resolution,
            n_pcs=args.n_pcs,
            verbose=args.verbose and not args.quiet
        )
        
        adata = inferencer.infercnv(
            adata,
            reference_key=args.reference_key,
            reference_cat=reference_cat
        )
        
        if args.plot_umap:
            inferencer.pca(adata)
            inferencer.neighbors(adata)
            inferencer.leiden(adata)
            inferencer.umap(adata)
        
        inferencer.compute_cnv_score(adata)
    
    output_path = Path(args.output)
    output_path.mkdir(parents=True, exist_ok=True)
    
    output_file = output_path / "cnv_annotated.h5ad"
    adata.write_h5ad(output_file)
    cli_logger.info(f"[CNV] Results saved to: {output_file}")
    
    groupby_col = args.groupby if args.groupby else args.reference_key
    
    if args.plot_heatmap and args.method == "infercnv":
        visualizer = CNVVisualizer(
            save_plots=True,
            cmap=args.cmap
        )
        
        visualizer.plot_cluster_cnv_heatmap(
            adata,
            groupby=groupby_col,
            output_dir=str(output_path),
            celltype="CNV Analysis"
        )
    
    if args.plot_spatial and args.spatial_key in adata.obsm:
        visualizer = CNVVisualizer(save_plots=True)
        
        visualizer.plot_cnv_spatial(
            adata,
            output_dir=str(output_path / "spatial"),
            sample_col=args.sample_col,
            spatial_key=args.spatial_key,
            title="CNV Score"
        )
    
    if args.method == "infercnv":
        summary = inferencer.get_cnv_summary(adata)
        
        log_section(cli_logger, "CNV Analysis Summary")
        if "n_cnv_clusters" in summary:
            cli_logger.info(f"CNV clusters: {summary['n_cnv_clusters']}")
        if "cnv_score_mean" in summary:
            cli_logger.info(f"Mean CNV score: {summary['cnv_score_mean']:.3f}")
        if "cnv_score_range" in summary:
            cli_logger.info(f"CNV score range: [{summary['cnv_score_range'][0]:.3f}, {summary['cnv_score_range'][1]:.3f}]")


def main():
    """Main entry point for command line interface"""
    
    parser = create_parser()
    args = parser.parse_args()
    
    if args.command is None:
        parser.print_help()
        sys.exit(0)
    
    if args.command == "run":
        run_command(args)
    elif args.command == "run-multi":
        run_multi_command(args)
    elif args.command == "info":
        info_command(args)
    elif args.command == "config":
        config_command(args)
    elif args.command == "corr":
        corr_command(args)
    elif args.command == "cnv":
        cnv_command(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()