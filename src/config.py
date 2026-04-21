"""
Configuration Management Module for SpaMFC

Supports YAML config file loading, parameter validation, and user feature selection.
Users can freely choose which features to use:
- use_spatial: Spatial neighborhood features
- use_cnv: CNV features (only for malignant cells)
- use_expression: Gene expression features
- use_niche: Niche features (scNiche)
"""

import yaml
from pathlib import Path
from typing import Dict, Any, Optional, List
from dataclasses import dataclass, field
import os


@dataclass
class FeatureConfig:
    """Feature configuration class for individual feature settings"""
    use: bool = True
    radius: int = 100
    k_neighbors: int = 15
    weight_method: str = "adaptive"
    fixed_weight: float = 0.3
    pca_dim: int = 50
    use_umap: bool = True
    umap_dim: int = 2
    marker_filter: bool = True
    min_expr_pct: float = 0.1
    pca_keep_dim: int = 30
    n_views: int = 3
    target_k: int = 15
    resolution: float = 0.3


@dataclass
class NMFConfig:
    """NMF configuration"""
    n_runs: int = 10
    n_components: int = 5
    init: str = "nndsvda"
    max_iter: int = 1000
    tol: float = 1e-4


@dataclass
class ClusteringConfig:
    """Clustering configuration"""
    method: str = "kmeans"
    n_clusters: int = 5
    resolution: float = 0.1
    per_sample: bool = True
    n_bootstrap: int = 30
    stability_threshold: float = 0.7


@dataclass
class MarkerAnalysisConfig:
    """Marker gene analysis configuration"""
    method: str = "wilcoxon"
    pval_threshold: float = 0.05
    logfc_threshold: float = 0.25
    top_n: int = 50


@dataclass
class EnrichmentConfig:
    """Functional enrichment configuration"""
    gene_sets: List[str] = field(default_factory=lambda: [
        "GO_Biological_Process_2023",
        "KEGG_2021_Human"
    ])
    pval_threshold: float = 0.05
    top_n: int = 20


@dataclass
class UnificationConfig:
    """Cross-sample unification configuration"""
    enable: bool = True
    marker_weight: float = 0.4
    pathway_weight: float = 0.3
    niche_weight: float = 0.3
    similarity_threshold: float = 0.6


@dataclass
class VisualizationConfig:
    """Visualization configuration"""
    save_plots: bool = True
    dpi: int = 300
    format: str = "pdf"
    show_plots: bool = False


@dataclass
class RuntimeConfig:
    """Runtime configuration"""
    n_jobs: int = 2
    verbose: bool = True
    chunksize: int = 500


@dataclass
class CorrelationConfig:
    """Gene correlation analysis configuration"""
    method: str = "spearman"
    p_adjust: str = "fdr_bh"
    threshold_p: float = 0.05
    min_corr_threshold: float = 0.0
    max_memory_mb: int = 512
    n_workers: Optional[int] = None
    batch_size: int = 500
    enable_numba: bool = True
    sample_spots: Optional[int] = None
    save_full_matrices: bool = False
    matrix_format: str = "npz"


@dataclass
class CNVInferenceConfig:
    """CNV inference configuration using inferCNVpy"""
    enable: bool = False
    gtf_file: str = ""
    reference_key: str = "cell_type"
    reference_cat: List[str] = field(default_factory=list)
    window_size: int = 100
    step_size: int = 10
    dynamic_threshold: float = 0.5
    exclude_genes: List[str] = field(default_factory=list)
    compute_scores: bool = True
    clustering_resolution: float = 0.5
    n_pcs: int = 30
    n_neighbors: int = 15
    run_umap: bool = True
    cmap: str = "bwr"


@dataclass
class SpaMFCConfig:
    """SpaMFC main configuration class"""
    input_path: str = ""
    output_dir: str = "./results/"
    sample_col: str = "sample_id"
    celltype_col: str = "anno_cell2location_res"
    spatial_key: str = "spatial"
    cnv_key: str = "X_cnv"
    
    use_spatial: bool = True
    use_cnv: bool = False
    use_expression: bool = True
    use_niche: bool = False
    
    fusion_method: str = "adaptive"
    
    spatial_config: FeatureConfig = field(default_factory=FeatureConfig)
    cnv_config: FeatureConfig = field(default_factory=FeatureConfig)
    expression_config: FeatureConfig = field(default_factory=FeatureConfig)
    niche_config: FeatureConfig = field(default_factory=FeatureConfig)
    
    nmf_config: NMFConfig = field(default_factory=NMFConfig)
    clustering_config: ClusteringConfig = field(default_factory=ClusteringConfig)
    marker_config: MarkerAnalysisConfig = field(default_factory=MarkerAnalysisConfig)
    enrichment_config: EnrichmentConfig = field(default_factory=EnrichmentConfig)
    unification_config: UnificationConfig = field(default_factory=UnificationConfig)
    visualization_config: VisualizationConfig = field(default_factory=VisualizationConfig)
    runtime_config: RuntimeConfig = field(default_factory=RuntimeConfig)
    correlation_config: CorrelationConfig = field(default_factory=CorrelationConfig)
    cnv_inference_config: CNVInferenceConfig = field(default_factory=CNVInferenceConfig)


class ConfigManager:
    """Configuration manager for loading and managing SpaMFC settings"""
    
    def __init__(self, config_path: Optional[str] = None):
        self.config = SpaMFCConfig()
        if config_path:
            self.load_from_yaml(config_path)
    
    def load_from_yaml(self, config_path: str) -> None:
        """Load configuration from YAML file"""
        with open(config_path, "r", encoding="utf-8") as f:
            yaml_config = yaml.safe_load(f)
        self._parse_yaml_config(yaml_config)
    
    def _parse_yaml_config(self, yaml_config: Dict) -> None:
        """Parse YAML configuration into SpaMFCConfig"""
        data_cfg = yaml_config.get("data", {})
        self.config.input_path = data_cfg.get("input_path", "")
        self.config.output_dir = data_cfg.get("output_dir", "./results/")
        self.config.sample_col = data_cfg.get("sample_col", "sample_id")
        self.config.celltype_col = data_cfg.get("celltype_col", "anno_cell2location_res")
        self.config.spatial_key = data_cfg.get("spatial_key", "spatial")
        self.config.cnv_key = data_cfg.get("cnv_key", "X_cnv")
        
        self.config.fusion_method = yaml_config.get("fusion_method", "adaptive")
        
        features_cfg = yaml_config.get("features", {})
        
        spatial_cfg = features_cfg.get("spatial", {})
        self.config.use_spatial = spatial_cfg.get("use", True)
        self.config.spatial_config = FeatureConfig(
            use=spatial_cfg.get("use", True),
            radius=spatial_cfg.get("radius", 100),
            k_neighbors=spatial_cfg.get("k_neighbors", 15),
            weight_method=spatial_cfg.get("weight_method", "adaptive"),
            fixed_weight=spatial_cfg.get("fixed_weight", 0.3)
        )
        
        cnv_cfg = features_cfg.get("cnv", {})
        self.config.use_cnv = cnv_cfg.get("use", False)
        self.config.cnv_config = FeatureConfig(
            use=cnv_cfg.get("use", False),
            pca_dim=cnv_cfg.get("pca_dim", 50),
            use_umap=cnv_cfg.get("use_umap", True),
            umap_dim=cnv_cfg.get("umap_dim", 2),
            marker_filter=cnv_cfg.get("marker_filter", True),
            fixed_weight=cnv_cfg.get("fixed_weight", 0.3)
        )
        
        expr_cfg = features_cfg.get("expression", {})
        self.config.use_expression = expr_cfg.get("use", True)
        self.config.expression_config = FeatureConfig(
            use=expr_cfg.get("use", True),
            min_expr_pct=expr_cfg.get("min_expr_pct", 0.1),
            pca_dim=expr_cfg.get("pca_dim", 50),
            use_umap=expr_cfg.get("use_umap", True),
            umap_dim=expr_cfg.get("umap_dim", 2),
            pca_keep_dim=expr_cfg.get("pca_keep_dim", 30),
            fixed_weight=expr_cfg.get("fixed_weight", 0.4)
        )
        
        niche_cfg = features_cfg.get("niche", {})
        self.config.use_niche = niche_cfg.get("use", False)
        self.config.niche_config = FeatureConfig(
            use=niche_cfg.get("use", False),
            n_views=niche_cfg.get("n_views", 3),
            target_k=niche_cfg.get("target_k", 15),
            resolution=niche_cfg.get("resolution", 0.3)
        )
        
        nmf_cfg = yaml_config.get("nmf", {})
        self.config.nmf_config = NMFConfig(
            n_runs=nmf_cfg.get("n_runs", 10),
            n_components=nmf_cfg.get("n_components", 5),
            init=nmf_cfg.get("init", "nndsvda"),
            max_iter=nmf_cfg.get("max_iter", 1000),
            tol=nmf_cfg.get("tol", 1e-4)
        )
        
        cluster_cfg = yaml_config.get("clustering", {})
        self.config.clustering_config = ClusteringConfig(
            method=cluster_cfg.get("method", "kmeans"),
            n_clusters=cluster_cfg.get("n_clusters", 5),
            resolution=cluster_cfg.get("resolution", 0.1),
            per_sample=cluster_cfg.get("per_sample", True),
            n_bootstrap=cluster_cfg.get("n_bootstrap", 30),
            stability_threshold=cluster_cfg.get("stability_threshold", 0.7)
        )
        
        marker_cfg = yaml_config.get("marker_analysis", {})
        self.config.marker_config = MarkerAnalysisConfig(
            method=marker_cfg.get("method", "wilcoxon"),
            pval_threshold=marker_cfg.get("pval_threshold", 0.05),
            logfc_threshold=marker_cfg.get("logfc_threshold", 0.25),
            top_n=marker_cfg.get("top_n", 50)
        )
        
        enrich_cfg = yaml_config.get("enrichment", {})
        self.config.enrichment_config = EnrichmentConfig(
            gene_sets=enrich_cfg.get("gene_sets", ["GO_Biological_Process_2023", "KEGG_2021_Human"]),
            pval_threshold=enrich_cfg.get("pval_threshold", 0.05),
            top_n=enrich_cfg.get("top_n", 20)
        )
        
        unif_cfg = yaml_config.get("unification", {})
        self.config.unification_config = UnificationConfig(
            enable=unif_cfg.get("enable", True),
            marker_weight=unif_cfg.get("marker_weight", 0.4),
            pathway_weight=unif_cfg.get("pathway_weight", 0.3),
            niche_weight=unif_cfg.get("niche_weight", 0.3),
            similarity_threshold=unif_cfg.get("similarity_threshold", 0.6)
        )
        
        vis_cfg = yaml_config.get("visualization", {})
        self.config.visualization_config = VisualizationConfig(
            save_plots=vis_cfg.get("save_plots", True),
            dpi=vis_cfg.get("dpi", 300),
            format=vis_cfg.get("format", "pdf"),
            show_plots=vis_cfg.get("show_plots", False)
        )
        
        runtime_cfg = yaml_config.get("runtime", {})
        self.config.runtime_config = RuntimeConfig(
            n_jobs=runtime_cfg.get("n_jobs", 2),
            verbose=runtime_cfg.get("verbose", True),
            chunksize=runtime_cfg.get("chunksize", 500)
        )
    
    def set_feature_usage(
        self,
        use_spatial: Optional[bool] = None,
        use_cnv: Optional[bool] = None,
        use_expression: Optional[bool] = None,
        use_niche: Optional[bool] = None
    ) -> None:
        """Dynamically set feature usage by user"""
        if use_spatial is not None:
            self.config.use_spatial = use_spatial
            self.config.spatial_config.use = use_spatial
        if use_cnv is not None:
            self.config.use_cnv = use_cnv
            self.config.cnv_config.use = use_cnv
        if use_expression is not None:
            self.config.use_expression = use_expression
            self.config.expression_config.use = use_expression
        if use_niche is not None:
            self.config.use_niche = use_niche
            self.config.niche_config.use = use_niche
    
    def validate_config(self) -> bool:
        """Validate configuration validity"""
        active_features = self.get_active_features()
        if len(active_features) == 0:
            raise ValueError("At least one feature type must be enabled")
        
        if self.config.use_cnv:
            if self.config.cnv_key not in ["X_cnv", "cnv"]:
                print(f"Warning: CNV key '{self.config.cnv_key}' may not exist in adata.obsm")
        
        output_dir = Path(self.config.output_dir)
        if not output_dir.exists():
            output_dir.mkdir(parents=True, exist_ok=True)
        
        return True
    
    def get_active_features(self) -> List[str]:
        """Get list of currently active features"""
        active = []
        if self.config.use_spatial:
            active.append("spatial")
        if self.config.use_cnv:
            active.append("cnv")
        if self.config.use_expression:
            active.append("expression")
        if self.config.use_niche:
            active.append("niche")
        return active
    
    def get_feature_weights(self) -> Dict[str, float]:
        """Get feature weights based on configuration"""
        weights = {}
        if self.config.use_spatial:
            if self.config.spatial_config.weight_method == "fixed":
                weights["spatial"] = self.config.spatial_config.fixed_weight
        if self.config.use_cnv:
            if self.config.cnv_config.weight_method == "fixed" or not hasattr(self.config.cnv_config, "weight_method"):
                weights["cnv"] = self.config.cnv_config.fixed_weight
        if self.config.use_expression:
            weights["expression"] = self.config.expression_config.fixed_weight
        if self.config.use_niche:
            weights["niche"] = self.config.niche_config.fixed_weight
        return weights
    
    def save_to_yaml(self, output_path: str) -> None:
        """Save current configuration to YAML file"""
        yaml_config = {
            "data": {
                "input_path": self.config.input_path,
                "output_dir": self.config.output_dir,
                "sample_col": self.config.sample_col,
                "celltype_col": self.config.celltype_col,
                "spatial_key": self.config.spatial_key,
                "cnv_key": self.config.cnv_key
            },
            "fusion_method": self.config.fusion_method,
            "features": {
                "spatial": {
                    "use": self.config.use_spatial,
                    "radius": self.config.spatial_config.radius,
                    "k_neighbors": self.config.spatial_config.k_neighbors,
                    "weight_method": self.config.spatial_config.weight_method,
                    "fixed_weight": self.config.spatial_config.fixed_weight
                },
                "cnv": {
                    "use": self.config.use_cnv,
                    "pca_dim": self.config.cnv_config.pca_dim,
                    "use_umap": self.config.cnv_config.use_umap,
                    "umap_dim": self.config.cnv_config.umap_dim,
                    "marker_filter": self.config.cnv_config.marker_filter,
                    "fixed_weight": self.config.cnv_config.fixed_weight
                },
                "expression": {
                    "use": self.config.use_expression,
                    "min_expr_pct": self.config.expression_config.min_expr_pct,
                    "pca_dim": self.config.expression_config.pca_dim,
                    "use_umap": self.config.expression_config.use_umap,
                    "umap_dim": self.config.expression_config.umap_dim,
                    "pca_keep_dim": self.config.expression_config.pca_keep_dim,
                    "fixed_weight": self.config.expression_config.fixed_weight
                },
                "niche": {
                    "use": self.config.use_niche,
                    "n_views": self.config.niche_config.n_views,
                    "target_k": self.config.niche_config.target_k,
                    "resolution": self.config.niche_config.resolution
                }
            },
            "nmf": {
                "n_runs": self.config.nmf_config.n_runs,
                "n_components": self.config.nmf_config.n_components,
                "init": self.config.nmf_config.init,
                "max_iter": self.config.nmf_config.max_iter,
                "tol": self.config.nmf_config.tol
            },
            "clustering": {
                "method": self.config.clustering_config.method,
                "n_clusters": self.config.clustering_config.n_clusters,
                "resolution": self.config.clustering_config.resolution,
                "per_sample": self.config.clustering_config.per_sample,
                "n_bootstrap": self.config.clustering_config.n_bootstrap,
                "stability_threshold": self.config.clustering_config.stability_threshold
            },
            "marker_analysis": {
                "method": self.config.marker_config.method,
                "pval_threshold": self.config.marker_config.pval_threshold,
                "logfc_threshold": self.config.marker_config.logfc_threshold,
                "top_n": self.config.marker_config.top_n
            },
            "enrichment": {
                "gene_sets": self.config.enrichment_config.gene_sets,
                "pval_threshold": self.config.enrichment_config.pval_threshold,
                "top_n": self.config.enrichment_config.top_n
            },
            "unification": {
                "enable": self.config.unification_config.enable,
                "marker_weight": self.config.unification_config.marker_weight,
                "pathway_weight": self.config.unification_config.pathway_weight,
                "niche_weight": self.config.unification_config.niche_weight,
                "similarity_threshold": self.config.unification_config.similarity_threshold
            },
            "visualization": {
                "save_plots": self.config.visualization_config.save_plots,
                "dpi": self.config.visualization_config.dpi,
                "format": self.config.visualization_config.format,
                "show_plots": self.config.visualization_config.show_plots
            },
            "runtime": {
                "n_jobs": self.config.runtime_config.n_jobs,
                "verbose": self.config.runtime_config.verbose,
                "chunksize": self.config.runtime_config.chunksize
            },
            "correlation": {
                "method": self.config.correlation_config.method,
                "p_adjust": self.config.correlation_config.p_adjust,
                "threshold_p": self.config.correlation_config.threshold_p,
                "min_corr_threshold": self.config.correlation_config.min_corr_threshold,
                "max_memory_mb": self.config.correlation_config.max_memory_mb,
                "n_workers": self.config.correlation_config.n_workers,
                "batch_size": self.config.correlation_config.batch_size,
                "enable_numba": self.config.correlation_config.enable_numba,
                "sample_spots": self.config.correlation_config.sample_spots,
                "save_full_matrices": self.config.correlation_config.save_full_matrices,
                "matrix_format": self.config.correlation_config.matrix_format
            }
        }
        
        with open(output_path, "w", encoding="utf-8") as f:
            yaml.dump(yaml_config, f, default_flow_style=False, allow_unicode=True)
    
    def print_config(self) -> None:
        """Print current configuration summary"""
        print("=" * 50)
        print("SpaMFC Configuration Summary")
        print("=" * 50)
        print(f"Input path: {self.config.input_path}")
        print(f"Output directory: {self.config.output_dir}")
        print(f"Sample column: {self.config.sample_col}")
        print(f"Cell type column: {self.config.celltype_col}")
        print("-" * 50)
        print(f"Feature fusion method: {self.config.fusion_method}")
        print("-" * 50)
        print("Feature Selection (User Configurable):")
        print(f"  use_spatial: {self.config.use_spatial}")
        print(f"  use_cnv: {self.config.use_cnv}")
        print(f"  use_expression: {self.config.use_expression}")
        print(f"  use_niche: {self.config.use_niche}")
        print(f"  Active features: {self.get_active_features()}")
        print("-" * 50)
        print("NMF Configuration:")
        print(f"  n_runs: {self.config.nmf_config.n_runs}")
        print(f"  n_components: {self.config.nmf_config.n_components}")
        print("-" * 50)
        print("Clustering Configuration:")
        print(f"  method: {self.config.clustering_config.method}")
        print(f"  n_clusters: {self.config.clustering_config.n_clusters}")
        print(f"  per_sample: {self.config.clustering_config.per_sample}")
        print("=" * 50)