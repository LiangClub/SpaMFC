"""
SpaMFC CNV Visualization Module

Provides CNV visualization functions using inferCNVpy plotting.
Features:
1. Chromosome heatmap visualization
2. CNV UMAP/t-SNE visualization
3. CNV score distribution plots
4. CNV spatial distribution plots
"""

import os
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional, List, Dict, Tuple, Any
import warnings

try:
    import matplotlib.pyplot as plt
    import seaborn as sns
except ImportError:
    warnings.warn("matplotlib or seaborn not installed")

try:
    import scanpy as sc
except ImportError:
    warnings.warn("scanpy not installed")

try:
    import infercnvpy as cnv
    INFERCNVPY_AVAILABLE = True
except ImportError:
    INFERCNVPY_AVAILABLE = False
    warnings.warn("infercnvpy not installed")


class CNVVisualizer:
    def __init__(
        self,
        save_plots: bool = True,
        dpi: int = 300,
        format: str = "pdf",
        show_plots: bool = False,
        cmap: str = "bwr"
    ):
        self.save_plots = save_plots
        self.dpi = dpi
        self.format = format
        self.show_plots = show_plots
        self.cmap = cmap
        
        if not INFERCNVPY_AVAILABLE:
            warnings.warn("infercnvpy not installed. CNV visualization will be limited.")
    
    def plot_chromosome_heatmap(
        self,
        adata,
        groupby: str = "cnv_leiden",
        use_rep: str = "cnv",
        figsize: Tuple[int, int] = (16, 10),
        output_path: Optional[str] = None,
        title: Optional[str] = None,
        **kwargs
    ):
        if not INFERCNVPY_AVAILABLE:
            warnings.warn("infercnvpy required for chromosome heatmap")
            return None
        
        if f"X_{use_rep}" not in adata.obsm:
            warnings.warn(f"CNV data not found in adata.obsm['X_{use_rep}']")
            return None
        
        if groupby not in adata.obs.columns:
            warnings.warn(f"Groupby column '{groupby}' not found in adata.obs")
            return None
        
        fig_axes = cnv.pl.chromosome_heatmap(
            adata,
            groupby=groupby,
            use_rep=use_rep,
            cmap=self.cmap,
            figsize=figsize,
            show=False,
            **kwargs
        )
        
        if title:
            if isinstance(fig_axes, dict):
                for ax in fig_axes.values():
                    ax.set_title(title)
            elif hasattr(fig_axes, 'figure'):
                fig_axes.figure.suptitle(title)
        
        if self.save_plots and output_path:
            output_dir = Path(output_path).parent
            output_dir.mkdir(parents=True, exist_ok=True)
            
            if isinstance(fig_axes, dict):
                fig = fig_axes.get('heatmap_ax', list(fig_axes.values())[0]).figure
            else:
                fig = fig_axes.figure if hasattr(fig_axes, 'figure') else fig_axes
            
            fig.savefig(output_path, dpi=self.dpi, bbox_inches="tight")
            print(f"[CNV Plot] Saved: {output_path}")
        
        if self.show_plots:
            plt.show()
        else:
            plt.close('all')
        
        return fig_axes
    
    def plot_chromosome_heatmap_summary(
        self,
        adata,
        groupby: str = "cnv_leiden",
        use_rep: str = "cnv",
        figsize: Tuple[int, int] = (16, 5),
        output_path: Optional[str] = None,
        title: Optional[str] = None,
        **kwargs
    ):
        if not INFERCNVPY_AVAILABLE:
            warnings.warn("infercnvpy required for chromosome heatmap summary")
            return None
        
        if f"X_{use_rep}" not in adata.obsm:
            warnings.warn(f"CNV data not found in adata.obsm['X_{use_rep}']")
            return None
        
        fig_axes = cnv.pl.chromosome_heatmap_summary(
            adata,
            groupby=groupby,
            use_rep=use_rep,
            cmap=self.cmap,
            figsize=figsize,
            show=False,
            **kwargs
        )
        
        if title:
            if isinstance(fig_axes, dict):
                for ax in fig_axes.values():
                    ax.set_title(title)
            elif hasattr(fig_axes, 'figure'):
                fig_axes.figure.suptitle(title)
        
        if self.save_plots and output_path:
            output_dir = Path(output_path).parent
            output_dir.mkdir(parents=True, exist_ok=True)
            
            if isinstance(fig_axes, dict):
                fig = fig_axes.get('heatmap_ax', list(fig_axes.values())[0]).figure
            else:
                fig = fig_axes.figure if hasattr(fig_axes, 'figure') else fig_axes
            
            fig.savefig(output_path, dpi=self.dpi, bbox_inches="tight")
            print(f"[CNV Plot] Saved: {output_path}")
        
        if self.show_plots:
            plt.show()
        else:
            plt.close('all')
        
        return fig_axes
    
    def plot_cnv_umap(
        self,
        adata,
        color: str = "cnv_leiden",
        figsize: Tuple[int, int] = (8, 6),
        output_path: Optional[str] = None,
        title: Optional[str] = None,
        **kwargs
    ):
        if not INFERCNVPY_AVAILABLE:
            warnings.warn("infercnvpy required for CNV UMAP")
            return None
        
        if "X_umap_cnv" not in adata.obsm:
            warnings.warn("CNV UMAP not found. Run cnv.tl.umap() first.")
            return None
        
        fig, ax = plt.subplots(figsize=figsize)
        
        cnv.pl.umap(
            adata,
            color=color,
            ax=ax,
            show=False,
            **kwargs
        )
        
        if title:
            ax.set_title(title)
        
        if self.save_plots and output_path:
            output_dir = Path(output_path).parent
            output_dir.mkdir(parents=True, exist_ok=True)
            fig.savefig(output_path, dpi=self.dpi, bbox_inches="tight")
            print(f"[CNV Plot] Saved: {output_path}")
        
        if self.show_plots:
            plt.show()
        else:
            plt.close(fig)
        
        return fig
    
    def plot_cnv_scores(
        self,
        adata,
        score_key: str = "cnv_score",
        groupby: str = "cnv_leiden",
        figsize: Tuple[int, int] = (10, 6),
        output_path: Optional[str] = None,
        title: Optional[str] = None,
        **kwargs
    ):
        if score_key not in adata.obs.columns:
            warnings.warn(f"CNV score not found in adata.obs['{score_key}']")
            return None
        
        fig, ax = plt.subplots(figsize=figsize)
        
        sc.pl.violin(
            adata,
            score_key,
            groupby=groupby,
            ax=ax,
            show=False,
            **kwargs
        )
        
        if title:
            ax.set_title(title)
        
        if self.save_plots and output_path:
            output_dir = Path(output_path).parent
            output_dir.mkdir(parents=True, exist_ok=True)
            fig.savefig(output_path, dpi=self.dpi, bbox_inches="tight")
            print(f"[CNV Plot] Saved: {output_path}")
        
        if self.show_plots:
            plt.show()
        else:
            plt.close(fig)
        
        return fig
    
    def plot_cnv_score_distribution(
        self,
        adata,
        score_key: str = "cnv_score",
        figsize: Tuple[int, int] = (8, 6),
        output_path: Optional[str] = None,
        bins: int = 50,
        title: Optional[str] = None
    ):
        if score_key not in adata.obs.columns:
            warnings.warn(f"CNV score not found in adata.obs['{score_key}']")
            return None
        
        scores = adata.obs[score_key].dropna()
        
        fig, ax = plt.subplots(figsize=figsize)
        
        ax.hist(scores, bins=bins, color='steelblue', edgecolor='white', alpha=0.7)
        
        ax.axvline(x=0, color='red', linestyle='--', linewidth=1.5, label='Normal (0)')
        
        ax.set_xlabel('CNV Score')
        ax.set_ylabel('Number of Cells')
        
        if title:
            ax.set_title(title)
        else:
            ax.set_title('CNV Score Distribution')
        
        ax.legend()
        
        if self.save_plots and output_path:
            output_dir = Path(output_path).parent
            output_dir.mkdir(parents=True, exist_ok=True)
            fig.savefig(output_path, dpi=self.dpi, bbox_inches="tight")
            print(f"[CNV Plot] Saved: {output_path}")
        
        if self.show_plots:
            plt.show()
        else:
            plt.close(fig)
        
        return fig
    
    def plot_cnv_spatial(
        self,
        adata,
        score_key: str = "cnv_score",
        sample_col: str = "sample_id",
        spatial_key: str = "spatial",
        figsize: Tuple[int, int] = (10, 8),
        output_dir: Optional[str] = None,
        cmap: str = "RdBu_r",
        title: Optional[str] = None
    ):
        if score_key not in adata.obs.columns:
            warnings.warn(f"CNV score not found in adata.obs['{score_key}']")
            return None
        
        if spatial_key not in adata.obsm:
            warnings.warn(f"Spatial coordinates not found in adata.obsm['{spatial_key}']")
            return None
        
        if output_dir:
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)
        
        samples = adata.obs[sample_col].unique()
        figs = []
        
        for sample in samples:
            sample_mask = adata.obs[sample_col] == sample
            adata_sample = adata[sample_mask].copy()
            
            fig, ax = plt.subplots(figsize=figsize)
            
            coords = adata_sample.obsm[spatial_key]
            scores = adata_sample.obs[score_key]
            
            scatter = ax.scatter(
                coords[:, 0],
                coords[:, 1],
                c=scores,
                cmap=cmap,
                s=10,
                alpha=0.8
            )
            
            plt.colorbar(scatter, ax=ax, label='CNV Score')
            
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            
            if title:
                ax.set_title(f"{sample} - {title}")
            else:
                ax.set_title(f"{sample} - CNV Score Spatial Distribution")
            
            ax.set_aspect('equal')
            
            if self.save_plots and output_dir:
                fig_path = output_path / f"{sample}_cnv_spatial.{self.format}"
                fig.savefig(fig_path, dpi=self.dpi, bbox_inches="tight")
                print(f"[CNV Plot] Saved: {fig_path}")
            
            figs.append(fig)
            
            if self.show_plots:
                plt.show()
            else:
                plt.close(fig)
        
        return figs
    
    def plot_cluster_cnv_heatmap(
        self,
        adata,
        groupby: str = "cnv_leiden",
        use_rep: str = "cnv",
        output_dir: Optional[str] = None,
        celltype: Optional[str] = None
    ):
        if output_dir:
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)
        
        title_prefix = f"{celltype} - " if celltype else ""
        
        heatmap_path = output_path / f"chromosome_heatmap.{self.format}" if output_dir else None
        self.plot_chromosome_heatmap(
            adata,
            groupby=groupby,
            use_rep=use_rep,
            output_path=heatmap_path,
            title=f"{title_prefix}CNV Chromosome Heatmap"
        )
        
        summary_path = output_path / f"chromosome_heatmap_summary.{self.format}" if output_dir else None
        self.plot_chromosome_heatmap_summary(
            adata,
            groupby=groupby,
            use_rep=use_rep,
            output_path=summary_path,
            title=f"{title_prefix}CNV Summary by Cluster"
        )
    
    def generate_all_cnv_plots(
        self,
        adata,
        output_dir: str,
        celltype: Optional[str] = None,
        sample_col: str = "sample_id",
        spatial_key: str = "spatial"
    ):
        output_path = Path(output_dir) / "cnv"
        output_path.mkdir(parents=True, exist_ok=True)
        
        self.plot_cluster_cnv_heatmap(
            adata,
            output_dir=str(output_path),
            celltype=celltype
        )
        
        if "cnv_score" in adata.obs.columns:
            score_dist_path = output_path / f"cnv_score_distribution.{self.format}"
            self.plot_cnv_score_distribution(
                adata,
                output_path=str(score_dist_path),
                title=f"{celltype} - CNV Score Distribution" if celltype else "CNV Score Distribution"
            )
            
            if "cnv_leiden" in adata.obs.columns:
                violin_path = output_path / f"cnv_score_by_cluster.{self.format}"
                self.plot_cnv_scores(
                    adata,
                    output_path=str(violin_path),
                    title=f"{celltype} - CNV Score by Cluster" if celltype else "CNV Score by Cluster"
                )
        
        if "X_umap_cnv" in adata.obsm and "cnv_leiden" in adata.obs.columns:
            umap_path = output_path / f"cnv_umap.{self.format}"
            self.plot_cnv_umap(
                adata,
                output_path=str(umap_path),
                title=f"{celltype} - CNV UMAP" if celltype else "CNV UMAP"
            )
        
        if spatial_key in adata.obsm and "cnv_score" in adata.obs.columns:
            spatial_dir = output_path / "spatial"
            self.plot_cnv_spatial(
                adata,
                output_dir=str(spatial_dir),
                sample_col=sample_col,
                spatial_key=spatial_key,
                title=f"{celltype} CNV Score" if celltype else "CNV Score"
            )
        
        print(f"[CNV Plot] All plots saved to: {output_path}")


def plot_chromosome_heatmap(adata, **kwargs):
    visualizer = CNVVisualizer()
    return visualizer.plot_chromosome_heatmap(adata, **kwargs)


def plot_cnv_umap(adata, **kwargs):
    visualizer = CNVVisualizer()
    return visualizer.plot_cnv_umap(adata, **kwargs)


def plot_cnv_scores(adata, **kwargs):
    visualizer = CNVVisualizer()
    return visualizer.plot_cnv_scores(adata, **kwargs)