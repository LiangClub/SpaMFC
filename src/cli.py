"""
SpaMFC Command Line Interface Module (Optimized)

Provides complete command line interface for SpaMFC analysis.
Simplified feature selection with single --features parameter.

Usage:
    spamfc_cli run --input data.h5ad --celltype-col "anno_cell2location_res" --celltype "CAFs" --features spatial,expression
    spamfc_cli run-multi --input data.h5ad --celltype-col "anno_cell2location_res" --celltypes "Malignant cells,CAFs,ILC"
    spamfc_cli info --input data.h5ad --celltype-col "anno_cell2location_res"
    spamfc_cli config --output ./my_config.yaml
"""

import argparse
import sys
import os
from pathlib import Path
from typing import Optional, List, Dict
import warnings

try:
    import scanpy as sc
except ImportError:
    warnings.warn("scanpy not installed")


VALID_FEATURES = ["spatial", "cnv", "expression", "niche"]
VALID_FUSION_METHODS = ["adaptive", "fixed", "concat"]


def create_parser() -> argparse.ArgumentParser:
    """Create main argument parser with all subcommands"""
    
    parser = argparse.ArgumentParser(
        prog="spamfc_cli",
        description="SpaMFC: Spatial Multi-Feature Clustering for Spatial Transcriptomics Subtype Analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # CAFs analysis with spatial+expression features
  spamfc_cli run --input data.h5ad --celltype-col "anno_cell2location_res" --celltype "CAFs" --features spatial,expression
  
  # Malignant cells with all features and fixed weights
  spamfc_cli run --input data.h5ad --celltype-col "anno_cell2location_res" --celltype "Malignant cells" --features spatial,cnv,expression --fusion-method fixed --weights 0.3,0.3,0.4
  
  # Multi-cell type analysis
  spamfc_cli run-multi --input data.h5ad --celltype-col "anno_cell2location_res" --celltypes "Malignant cells,CAFs,ILC" --features spatial,expression
  
  # Show data info
  spamfc_cli info --input data.h5ad --celltype-col "anno_cell2location_res"
  
  # Generate config file
  spamfc_cli config --output ./my_config.yaml --template cafs
"""
    )
    
    parser.add_argument(
        "--version", "-v",
        action="version",
        version="SpaMFC v2.0.0"
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
            print(f"Warning: Unknown feature '{feature}', ignoring. Valid features: {VALID_FEATURES}")
    
    if not any(features.values()):
        print(f"Warning: No valid features specified, using default 'spatial'")
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
    
    weights_list = [float(w.strip()) for w in weights_str.split(",")]
    
    if len(weights_list) != len(features_list):
        print(f"Warning: Number of weights ({len(weights_list)}) does not match number of features ({len(features_list)})")
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
    
    return weights


def validate_features_for_celltype(celltype: str, features: Dict[str, bool]) -> Dict[str, bool]:
    """Validate and adjust features based on cell type"""
    
    malignant_keywords = ["malignant", "tumor", "cancer", "恶性", "肿瘤"]
    is_malignant = any(kw in celltype.lower() for kw in malignant_keywords)
    
    if features["use_cnv"] and not is_malignant:
        print(f"Warning: CNV features are only recommended for malignant cells. "
              f"Cell type '{celltype}' may not be malignant. CNV will still be used if data exists.")
    
    return features


def run_command(args):
    """Execute 'run' subcommand"""
    
    try:
        import scanpy as sc
    except ImportError:
        print("Error: scanpy is not installed. Please install it first.")
        sys.exit(1)
    
    from ..pipeline import SpaMFCPipeline
    from ..config import ConfigManager
    
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
    
    print(f"\n{'='*60}")
    print(f"SpaMFC Command Line Analysis")
    print(f"{'='*60}")
    print(f"Input: {args.input}")
    print(f"Output: {args.output}")
    print(f"Cell type column: {args.celltype_col}")
    print(f"Cell type: {args.celltype}")
    print(f"Features: {args.features}")
    print(f"Fusion method: {args.fusion_method}")
    if args.weights:
        print(f"Fixed weights: {args.weights}")
    print(f"{'='*60}\n")
    
    adata = sc.read_h5ad(args.input)
    
    if args.celltype_col not in adata.obs:
        print(f"Error: Cell type column '{args.celltype_col}' not found in adata.obs")
        print(f"Available columns: {list(adata.obs.columns)}")
        sys.exit(1)
    
    if args.celltype not in adata.obs[args.celltype_col].values:
        print(f"Error: Cell type '{args.celltype}' not found in column '{args.celltype_col}'")
        print(f"Available cell types: {list(adata.obs[args.celltype_col].unique())}")
        sys.exit(1)
    
    markers = args.markers if args.markers else None
    
    adata = pipeline.run(adata, args.celltype, markers)
    
    output_path = Path(args.output) / args.celltype.replace(" ", "_")
    output_path.mkdir(parents=True, exist_ok=True)
    
    output_file = output_path / f"{args.celltype.replace(' ', '_')}_subtype_annotated.h5ad"
    adata.write_h5ad(output_file)
    
    print(f"\nResults saved to: {output_file}")
    print(f"{'='*60}\n")


def run_multi_command(args):
    """Execute 'run-multi' subcommand"""
    
    try:
        import scanpy as sc
    except ImportError:
        print("Error: scanpy is not installed. Please install it first.")
        sys.exit(1)
    
    from ..pipeline import SpaMFCPipeline
    
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
    
    print(f"\n{'='*60}")
    print(f"SpaMFC Multi-Cell Type Analysis")
    print(f"{'='*60}")
    print(f"Input: {args.input}")
    print(f"Output: {args.output}")
    print(f"Cell type column: {args.celltype_col}")
    print(f"Cell types: {celltypes}")
    print(f"Features: {args.features}")
    print(f"Fusion method: {args.fusion_method}")
    if args.weights:
        print(f"Fixed weights: {args.weights}")
    print(f"{'='*60}\n")
    
    adata = sc.read_h5ad(args.input)
    
    if args.celltype_col not in adata.obs:
        print(f"Error: Cell type column '{args.celltype_col}' not found in adata.obs")
        print(f"Available columns: {list(adata.obs.columns)}")
        sys.exit(1)
    
    available_celltypes = list(adata.obs[args.celltype_col].unique())
    missing_celltypes = [ct for ct in celltypes if ct not in available_celltypes]
    
    if missing_celltypes:
        print(f"Warning: Cell types not found: {missing_celltypes}")
        print(f"Available cell types: {available_celltypes}")
        celltypes = [ct for ct in celltypes if ct in available_celltypes]
        
        if not celltypes:
            print("Error: No valid cell types to analyze")
            sys.exit(1)
    
    adata = pipeline.run_multi_celltypes(adata, celltypes)
    
    output_file = Path(args.output) / "all_celltypes_subtype_annotated.h5ad"
    adata.write_h5ad(output_file)
    
    print(f"\nResults saved to: {output_file}")
    print(f"{'='*60}\n")


def info_command(args):
    """Execute 'info' subcommand"""
    
    try:
        import scanpy as sc
    except ImportError:
        print("Error: scanpy is not installed. Please install it first.")
        sys.exit(1)
    
    print(f"\n{'='*60}")
    print(f"Data Information: {args.input}")
    print(f"{'='*60}\n")
    
    adata = sc.read_h5ad(args.input)
    
    print(f"Cells: {adata.n_obs}")
    print(f"Genes: {adata.n_vars}")
    
    if args.celltype_col in adata.obs:
        celltypes = adata.obs[args.celltype_col].value_counts()
        print(f"\nCell Type Distribution (from '{args.celltype_col}'):")
        for ct, count in celltypes.items():
            print(f"  {ct}: {count} ({count/adata.n_obs*100:.1f}%)")
    else:
        print(f"\nWarning: Cell type column '{args.celltype_col}' not found")
        print(f"Available columns: {list(adata.obs.columns)}")
    
    if args.sample_col in adata.obs:
        samples = adata.obs[args.sample_col].unique()
        print(f"\nSamples: {len(samples)} ({', '.join(samples[:5])}...)")
    
    if "spatial" in adata.obsm:
        print(f"\nSpatial coordinates: Available")
    
    if "X_cnv" in adata.obsm:
        print(f"CNV matrix: Available")
    
    print(f"\n{'='*60}\n")


def config_command(args):
    """Execute 'config' subcommand"""
    
    from ..config import ConfigManager
    
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
        print(f"\nConfig file generated: {args.output}")
        print(f"Template: {args.template}")
    else:
        config_manager = ConfigManager()
        config_manager.save_to_yaml(args.output)
        print(f"\nDefault config file generated: {args.output}")
    
    print(f"\nEdit the config file to customize your analysis parameters.")
    print(f"Usage: spamfc_cli run --config {args.output} --input data.h5ad --celltype-col 'anno_cell2location_res' --celltype 'CAFs'\n")


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
    else:
        parser.print_help()


if __name__ == "__main__":
    main()