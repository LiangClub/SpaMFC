"""
SpaMFC Main Pipeline Entry Point

Integrates all modules to provide a complete analysis workflow:
1. Feature extraction (spatial, CNV, expression, niche) - user configurable
2. Feature fusion with adaptive weighting
3. NMF consensus decomposition
4. Subtype clustering (per-sample)
5. Marker gene analysis
6. Functional enrichment analysis
7. Niche analysis
8. Cross-sample subtype unification
9. Visualization and report generation

Users can freely choose which features to use via configuration.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional, Dict, List
import warnings
import time

try:
    import scanpy as sc
except ImportError:
    warnings.warn("scanpy not installed")

from .config import ConfigManager, SpaMFCConfig
from .features.spatial import SpatialFeatureExtractor
from .features.cnv import CNVFeatureProcessor
from .features.expression import ExpressionFeatureProcessor
from .features.niche import NicheFeatureProcessor
from .fusion.weighting import FeatureFusion, AdaptiveWeightCalculator
from .fusion.nmf import NMFConsensus
from .fusion.clustering import SubtypeClusterer
from .annotation.markers import MarkerGeneAnalyzer
from .annotation.enrichment import EnrichmentAnalyzer
from .annotation.niche_analysis import NicheAnalyzer
from .unification.similarity import SimilarityCalculator
from .unification.fusion import SimilarityFusion
from .unification.mapping import SubtypeMapper
from .visualization.spatial_plot import SpatialVisualizer
from .visualization.report import ReportGenerator
from .utils.logger import setup_logger, log_section, log_step


class SpaMFCPipeline:
    """SpaMFC analysis pipeline"""
    
    def __init__(self, config_path: Optional[str] = None):
        """
        Initialize SpaMFC pipeline
        
        Parameters:
            config_path: Path to YAML configuration file
        """
        self.config_manager = ConfigManager(config_path)
        self.config = self.config_manager.config
        
        self.logger = setup_logger(
            "SpaMFC",
            log_dir=Path(self.config.output_dir) / "logs",
            verbose=self.config.runtime_config.verbose
        )
        
        self._init_modules()
        
        self.results = {}
    
    def _init_modules(self):
        """Initialize all functional modules"""
        self.spatial_extractor = SpatialFeatureExtractor(self.config.spatial_config)
        self.cnv_processor = CNVFeatureProcessor(self.config.cnv_config)
        self.expr_processor = ExpressionFeatureProcessor(self.config.expression_config)
        self.niche_processor = NicheFeatureProcessor(self.config.niche_config)
        
        self.feature_fusion = FeatureFusion(
            fusion_method=self.config.fusion_method
        )
        self.weight_calculator = AdaptiveWeightCalculator(
            weight_method=self.config.fusion_method
        )
        self.nmf_consensus = NMFConsensus(
            n_runs=self.config.nmf_config.n_runs,
            n_components=self.config.nmf_config.n_components,
            init=self.config.nmf_config.init,
            max_iter=self.config.nmf_config.max_iter,
            tol=self.config.nmf_config.tol
        )
        self.clusterer = SubtypeClusterer(
            method=self.config.clustering_config.method,
            n_clusters=self.config.clustering_config.n_clusters,
            resolution=self.config.clustering_config.resolution,
            per_sample=self.config.clustering_config.per_sample,
            n_bootstrap=self.config.clustering_config.n_bootstrap,
            stability_threshold=self.config.clustering_config.stability_threshold
        )
        
        self.marker_analyzer = MarkerGeneAnalyzer(
            method=self.config.marker_config.method,
            pval_threshold=self.config.marker_config.pval_threshold,
            logfc_threshold=self.config.marker_config.logfc_threshold,
            top_n=self.config.marker_config.top_n
        )
        self.enrichment_analyzer = EnrichmentAnalyzer(
            gene_sets=self.config.enrichment_config.gene_sets,
            pval_threshold=self.config.enrichment_config.pval_threshold,
            top_n=self.config.enrichment_config.top_n
        )
        self.niche_analyzer = NicheAnalyzer(
            n_views=self.config.niche_config.n_views,
            target_k=self.config.niche_config.target_k,
            resolution=self.config.niche_config.resolution
        )
        
        self.similarity_calculator = SimilarityCalculator(
            marker_weight=self.config.unification_config.marker_weight,
            pathway_weight=self.config.unification_config.pathway_weight,
            niche_weight=self.config.unification_config.niche_weight
        )
        self.similarity_fusion = SimilarityFusion(
            marker_weight=self.config.unification_config.marker_weight,
            pathway_weight=self.config.unification_config.pathway_weight,
            niche_weight=self.config.unification_config.niche_weight
        )
        self.subtype_mapper = SubtypeMapper(
            similarity_threshold=self.config.unification_config.similarity_threshold
        )
        
        self.visualizer = SpatialVisualizer(
            save_plots=self.config.visualization_config.save_plots,
            dpi=self.config.visualization_config.dpi,
            format=self.config.visualization_config.format,
            show_plots=self.config.visualization_config.show_plots
        )
        self.report_generator = ReportGenerator()
    
    def set_feature_usage(
        self,
        use_spatial: Optional[bool] = None,
        use_cnv: Optional[bool] = None,
        use_expression: Optional[bool] = None,
        use_niche: Optional[bool] = None
    ):
        """
        Dynamically set feature usage (user configurable)
        
        Parameters:
            use_spatial: Whether to use spatial features
            use_cnv: Whether to use CNV features
            use_expression: Whether to use expression features
            use_niche: Whether to use niche features
        """
        self.config_manager.set_feature_usage(
            use_spatial, use_cnv, use_expression, use_niche
        )
        
        self.config = self.config_manager.config
        
        self.spatial_extractor.use = self.config.use_spatial
        self.cnv_processor.use = self.config.use_cnv
        self.expr_processor.use = self.config.use_expression
        self.niche_processor.use = self.config.use_niche
        
        self.feature_fusion = FeatureFusion(
            fusion_method=self.config.fusion_method
        )
        self.weight_calculator = AdaptiveWeightCalculator(
            weight_method=self.config.fusion_method
        )
    
    def run(
        self,
        adata,
        celltype: str,
        markers: Optional[List[str]] = None
    ):
        """
        Execute complete analysis workflow for a single cell type
        
        Parameters:
            adata: AnnData object
            celltype: Target cell type name
            markers: Tumor markers for CNV filtering (optional)
        
        Returns:
            Updated adata with subtype labels and analysis results
        """
        verbose = self.config.runtime_config.verbose
        
        if verbose:
            log_section(self.logger, f"SpaMFC Analysis for: {celltype}")
            self.logger.info(f"Active features: {self.config_manager.get_active_features()}")
        
        start_time = time.time()
        
        self.config_manager.validate_config()
        
        celltype_col = self.config.celltype_col
        sample_col = self.config.sample_col
        
        if not adata.obs_names.is_unique:
            adata.obs_names_make_unique()
        
        target_mask = adata.obs[celltype_col] == celltype
        target_cells = adata.obs[target_mask].index.tolist()
        target_indices = np.where(target_mask.values)[0].tolist()
        
        if len(target_cells) == 0:
            raise ValueError(f"Cell type '{celltype}' not found in adata.obs['{celltype_col}']")
        
        if verbose:
            self.logger.info(f"Found {len(target_cells)} cells for {celltype}")
        
        feature_dfs = {}
        
        if self.config.use_spatial:
            if verbose:
                log_step(self.logger, 1, 14, "Extracting spatial neighborhood features")
            
            global_cell_types = adata.obs[celltype_col].unique().tolist()
            
            spatial_df = self.spatial_extractor.extract_per_sample(
                adata, target_cells, target_indices, celltype_col, sample_col,
                self.config.spatial_key
            )
            
            if spatial_df is not None:
                feature_dfs["spatial"] = spatial_df
                if verbose:
                    self.logger.info(f"  Spatial features: {spatial_df.shape}")
        
        if self.config.use_cnv:
            if verbose:
                log_step(self.logger, 2, 14, "Processing CNV features")
            
            cnv_df = self.cnv_processor.process_per_sample(
                adata, target_cells, target_indices,
                self.config.cnv_key,
                sample_col,
                markers
            )
            
            if cnv_df is not None:
                feature_dfs["cnv"] = cnv_df
                if verbose:
                    self.logger.info(f"  CNV features: {cnv_df.shape}")
        
        if self.config.use_expression:
            if verbose:
                log_step(self.logger, 3, 14, "Processing expression features")
            
            expr_df = self.expr_processor.process_per_sample(
                adata, target_cells, target_indices, sample_col
            )
            
            if expr_df is not None:
                feature_dfs["expression"] = expr_df
                if verbose:
                    self.logger.info(f"  Expression features: {expr_df.shape}")
        
        if self.config.use_niche:
            if verbose:
                log_step(self.logger, 4, 14, "Processing niche features")
            
            niche_df = self.niche_processor.process(
                adata, target_cells, celltype_col, sample_col
            )
            
            if niche_df is not None:
                feature_dfs["niche"] = niche_df
                if verbose:
                    self.logger.info(f"  Niche features: {niche_df.shape}")
        
        if len(feature_dfs) == 0:
            raise ValueError("No features extracted. Please enable at least one feature type.")
        
        if verbose:
            log_step(self.logger, 5, 14, "Calculating feature weights")
        
        fixed_weights = {
            "spatial": self.config.spatial_config.fixed_weight,
            "cnv": self.config.cnv_config.fixed_weight,
            "expression": self.config.expression_config.fixed_weight,
            "niche": self.config.niche_config.fixed_weight
        }
        
        weights = self.feature_fusion.get_weights(feature_dfs, fixed_weights)
        
        if verbose:
            self.logger.info(f"  Fusion method: {self.config.fusion_method}")
            self.logger.info(f"  Weights: {weights}")
        
        if verbose:
            log_step(self.logger, 6, 14, "Fusing features")
        
        fused_df = self.feature_fusion.fuse(feature_dfs, fixed_weights)
        
        if verbose:
            self.logger.info(f"  Fused features: {fused_df.shape}")
        
        if verbose:
            log_step(self.logger, 7, 14, "NMF consensus decomposition")
        
        W_df = self.nmf_consensus.decompose(fused_df)
        
        if verbose:
            self.logger.info(f"  NMF embedding: {W_df.shape}")
        
        if verbose:
            log_step(self.logger, 8, 14, "Subtype clustering")
        
        adata = self.clusterer.cluster(adata, W_df, sample_col, celltype)
        
        subtype_col = f"{celltype}_subtype"
        
        if verbose:
            subtype_counts = adata.obs[subtype_col].value_counts()
            self.logger.info(f"  Subtypes identified: {len(subtype_counts)}")
            self.logger.info(f"  Subtype distribution:")
            for subtype, count in subtype_counts.items():
                self.logger.info(f"    {subtype}: {count}")
        
        if verbose:
            log_step(self.logger, 9, 14, "Marker gene analysis")
        
        markers_dict = self.marker_analyzer.analyze(adata, celltype, sample_col)
        
        self.results[celltype] = {
            "markers": markers_dict,
            "weights": weights,
            "nmf_embedding": W_df
        }
        
        if verbose:
            self.logger.info(f"  Marker genes identified for {len(markers_dict)} subtypes")
        
        if verbose:
            log_step(self.logger, 10, 14, "Functional enrichment analysis")
        
        enrichment_results = self.enrichment_analyzer.analyze(markers_dict)
        
        self.results[celltype]["enrichment"] = enrichment_results
        
        if verbose:
            n_enriched = sum(1 for r in enrichment_results.values() if r is not None)
            self.logger.info(f"  Enrichment results for {n_enriched} subtypes")
        
        niche_profiles = None
        if self.config.use_niche:
            if verbose:
                log_step(self.logger, 11, 14, "Niche analysis")
            
            niche_profiles = self.niche_analyzer.get_subtype_niche_profiles(
                adata, subtype_col, "scNiche"
            )
            
            self.results[celltype]["niche_profiles"] = niche_profiles
        
        if self.config.unification_config.enable:
            if verbose:
                log_step(self.logger, 12, 14, "Cross-sample subtype unification")
            
            subtypes = list(markers_dict.keys())
            
            if len(subtypes) >= 2:
                marker_sim = self.similarity_calculator.calculate_marker_similarity(markers_dict)
                pathway_sim = self.similarity_calculator.calculate_pathway_similarity(enrichment_results)
                
                if niche_profiles is not None:
                    niche_sim = self.similarity_calculator.calculate_niche_similarity(niche_profiles)
                else:
                    niche_sim = None
                
                fused_sim = self.similarity_fusion.fuse(marker_sim, pathway_sim, niche_sim)
                
                unified_names = self.subtype_mapper.generate_mapping(fused_sim, subtypes)
                
                unified_col = f"{celltype}_subtype_unified"
                adata = self.subtype_mapper.apply_mapping(adata, subtype_col, unified_names, unified_col)
                
                mapping_table = self.subtype_mapper.create_mapping_table(unified_names)
                
                self.results[celltype]["mapping_table"] = mapping_table
                
                if verbose:
                    self.logger.info(f"  Unified subtypes: {len(set(unified_names.values()))}")
                    self.logger.info(f"  Mapping summary:")
                    for unified, count in pd.Series(unified_names.values()).value_counts().items():
                        self.logger.info(f"    {unified}: {count} original subtypes")
        
        if self.config.visualization_config.save_plots:
            if verbose:
                log_step(self.logger, 13, 14, "Generating visualizations")
            
            output_dir = self.config.output_dir
            
            self.visualizer.plot_spatial_subtype(adata, celltype, output_dir, subtype_col, sample_col=sample_col)
            
            self.visualizer.plot_marker_genes(adata, celltype, output_dir, markers_dict, sample_col=sample_col)
            
            self.visualizer.plot_subtype_proportion(adata, celltype, output_dir, subtype_col, sample_col=sample_col)
        
        if verbose:
            log_step(self.logger, 14, 14, "Generating report")
        
        mapping_table = self.results[celltype].get("mapping_table")
        
        self.report_generator.generate(
            adata, celltype, markers_dict, enrichment_results,
            self.config.output_dir, mapping_table, niche_profiles
        )
        
        elapsed_time = time.time() - start_time
        
        if verbose:
            log_section(self.logger, f"Analysis completed in {elapsed_time:.2f} seconds")
        
        return adata
    
    def run_multi_celltypes(
        self,
        adata,
        celltypes: List[str],
        markers_dict: Optional[Dict[str, List[str]]] = None
    ):
        """
        Execute analysis for multiple cell types
        
        Parameters:
            adata: AnnData object
            celltypes: List of cell type names
            markers_dict: Dictionary of tumor markers for each cell type (optional)
        
        Returns:
            Updated adata with subtype labels for all cell types
        """
        verbose = self.config.runtime_config.verbose
        
        if verbose:
            log_section(self.logger, "SpaMFC Multi-Cell Type Analysis")
            self.logger.info(f"Cell types: {celltypes}")
        
        for celltype in celltypes:
            markers = markers_dict.get(celltype) if markers_dict else None
            
            original_config = {
                "use_spatial": self.config.use_spatial,
                "use_cnv": self.config.use_cnv,
                "use_expression": self.config.use_expression,
                "use_niche": self.config.use_niche
            }
            
            if celltype in ["Malignant cells", "Malignant", "Tumor cells"]:
                self.set_feature_usage(use_spatial=True, use_cnv=True, use_expression=True)
            elif celltype in ["CAFs", "Fibroblasts", "Cancer Associated Fibroblasts"]:
                self.set_feature_usage(use_spatial=True, use_cnv=False, use_expression=True)
            elif celltype in ["ILC", "T cells", "B cells", "Immune cells", "Lymphocytes"]:
                self.set_feature_usage(use_spatial=False, use_cnv=False, use_expression=True)
            elif celltype in ["ECs", "Endothelial cells", "Endothelial"]:
                self.set_feature_usage(use_spatial=False, use_cnv=False, use_expression=True)
            else:
                self.set_feature_usage(use_spatial=True, use_cnv=False, use_expression=True)
            
            adata = self.run(adata, celltype, markers)
            
            self.set_feature_usage(**original_config)
        
        if verbose:
            log_section(self.logger, "All cell types analyzed successfully")
        
        return adata
    
    def get_results(self, celltype: str) -> Dict:
        """Get analysis results for a specific cell type"""
        return self.results.get(celltype, {})
    
    def get_all_results(self) -> Dict:
        """Get all analysis results"""
        return self.results
    
    def save_results(self, output_path: str):
        """Save all results to files"""
        output_dir = Path(output_path)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        for celltype, result in self.results.items():
            celltype_safe = celltype.replace(" ", "_")
            celltype_dir = output_dir / celltype_safe
            celltype_dir.mkdir(parents=True, exist_ok=True)
            
            if "markers" in result:
                markers_df = self.marker_analyzer.get_marker_dataframe(result["markers"])
                markers_df.to_csv(celltype_dir / "markers.csv", index=False)
            
            if "weights" in result:
                weights_df = pd.DataFrame({
                    "feature": list(result["weights"].keys()),
                    "weight": list(result["weights"].values())
                })
                weights_df.to_csv(celltype_dir / "weights.csv", index=False)
            
            if "mapping_table" in result:
                result["mapping_table"].to_csv(celltype_dir / "mapping.csv", index=False)
    
    def print_config(self):
        """Print current configuration"""
        self.config_manager.print_config()


def run_spamfc(
    adata,
    celltypes: List[str],
    config_path: Optional[str] = None,
    output_dir: str = "./results/",
    use_spatial: bool = True,
    use_cnv: bool = False,
    use_expression: bool = True,
    use_niche: bool = False
):
    """
    Convenience function to run SpaMFC analysis
    
    Parameters:
        adata: AnnData object
        celltypes: List of cell type names to analyze
        config_path: Path to configuration file (optional)
        output_dir: Output directory
        use_spatial: Whether to use spatial features
        use_cnv: Whether to use CNV features
        use_expression: Whether to use expression features
        use_niche: Whether to use niche features
    
    Returns:
        Updated adata with subtype labels
    """
    pipeline = SpaMFCPipeline(config_path)
    
    pipeline.config.output_dir = output_dir
    
    pipeline.set_feature_usage(use_spatial, use_cnv, use_expression, use_niche)
    
    adata = pipeline.run_multi_celltypes(adata, celltypes)
    
    return adata