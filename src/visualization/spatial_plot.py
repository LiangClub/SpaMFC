"""
Spatial Distribution Visualization Module

Visualizes subtype spatial distribution using scanpy's spatial plotting.
Supports multiple visualization types:
- Subtype spatial distribution
- Marker gene expression
- Niche distribution
- CNV score distribution

Users can configure:
- save_plots: Whether to save plots
- dpi: Plot resolution
- format: Plot format (pdf, png, svg)
- show_plots: Whether to display plots
"""

import os
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional, List, Dict
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


class SpatialVisualizer:
    """Spatial distribution visualizer"""
    
    def __init__(
        self,
        save_plots: bool = True,
        dpi: int = 300,
        format: str = "pdf",
        show_plots: bool = False
    ):
        """
        Initialize spatial visualizer
        
        Parameters:
            save_plots: Whether to save plots
            dpi: Plot resolution
            format: Plot format
            show_plots: Whether to display plots
        """
        self.save_plots = save_plots
        self.dpi = dpi
        self.format = format
        self.show_plots = show_plots
    
    def plot_spatial_subtype(
        self,
        adata,
        celltype: str,
        output_dir: str,
        subtype_col: Optional[str] = None,
        unified_col: Optional[str] = None,
        sample_col: str = "sample_id"
    ):
        """
        Plot subtype spatial distribution
        
        Parameters:
            adata: AnnData object
            celltype: Cell type name
            output_dir: Output directory
            subtype_col: Column name for original subtype
            unified_col: Column name for unified subtype
            sample_col: Column name for sample ID
        """
        try:
            import matplotlib.pyplot as plt
            import scanpy as sc
        except ImportError:
            warnings.warn("matplotlib or scanpy not installed")
            return
        
        subtype_col = subtype_col or f"{celltype}_subtype"
        unified_col = unified_col or f"{celltype}_subtype_unified"
        
        celltype_safe = celltype.replace(" ", "_")
        output_path = Path(output_dir) / celltype_safe / "spatial"
        output_path.mkdir(parents=True, exist_ok=True)
        
        samples = adata.obs[sample_col].unique()
        
        for sample in samples:
            sample_mask = adata.obs[sample_col] == sample
            adata_sample = adata[sample_mask].copy()
            
            if "spatial" not in adata_sample.obsm:
                warnings.warn(f"Spatial coordinates not found for sample {sample}")
                continue
            
            if subtype_col in adata_sample.obs:
                subtype_mask = adata_sample.obs[subtype_col].notna()
                adata_subtype = adata_sample[subtype_mask].copy()
                
                if len(adata_subtype) > 0:
                    fig, ax = plt.subplots(figsize=(10, 8))
                    
                    try:
                        sc.pl.spatial(
                            adata_subtype,
                            color=subtype_col,
                            ax=ax,
                            show=False,
                            title=f"{sample} - {celltype} Subtypes"
                        )
                        
                        if self.save_plots:
                            fig_path = output_path / f"{sample}_subtype.{self.format}"
                            fig.savefig(fig_path, dpi=self.dpi, bbox_inches="tight")
                        
                        if self.show_plots:
                            plt.show()
                        else:
                            plt.close()
                            
                    except Exception as e:
                        warnings.warn(f"Failed to plot spatial for {sample}: {e}")
                        plt.close()
            
            if unified_col in adata_sample.obs:
                unified_mask = adata_sample.obs[unified_col].notna()
                adata_unified = adata_sample[unified_mask].copy()
                
                if len(adata_unified) > 0:
                    fig, ax = plt.subplots(figsize=(10, 8))
                    
                    try:
                        sc.pl.spatial(
                            adata_unified,
                            color=unified_col,
                            ax=ax,
                            show=False,
                            title=f"{sample} - {celltype} Unified Subtypes"
                        )
                        
                        if self.save_plots:
                            fig_path = output_path / f"{sample}_unified.{self.format}"
                            fig.savefig(fig_path, dpi=self.dpi, bbox_inches="tight")
                        
                        if self.show_plots:
                            plt.show()
                        else:
                            plt.close()
                            
                    except Exception as e:
                        warnings.warn(f"Failed to plot unified spatial for {sample}: {e}")
                        plt.close()
    
    def plot_marker_genes(
        self,
        adata,
        celltype: str,
        output_dir: str,
        markers_dict: Dict[str, List[str]],
        sample_col: str = "sample_id",
        top_n: int = 5
    ):
        """
        Plot marker gene expression
        
        Parameters:
            adata: AnnData object
            celltype: Cell type name
            output_dir: Output directory
            markers_dict: Dictionary of marker genes
            sample_col: Column name for sample ID
            top_n: Number of top markers to plot
        """
        try:
            import matplotlib.pyplot as plt
            import scanpy as sc
        except ImportError:
            warnings.warn("matplotlib or scanpy not installed")
            return
        
        celltype_safe = celltype.replace(" ", "_")
        output_path = Path(output_dir) / celltype_safe / "markers"
        output_path.mkdir(parents=True, exist_ok=True)
        
        all_markers = set()
        for genes in markers_dict.values():
            all_markers.update(genes[:top_n])
        
        available_markers = [g for g in all_markers if g in adata.var_names]
        
        if len(available_markers) == 0:
            warnings.warn("No marker genes found in adata")
            return
        
        samples = adata.obs[sample_col].unique()
        
        for sample in samples:
            sample_mask = adata.obs[sample_col] == sample
            adata_sample = adata[sample_mask].copy()
            
            if "spatial" not in adata_sample.obsm:
                continue
            
            for marker in available_markers[:min(10, len(available_markers))]:
                fig, ax = plt.subplots(figsize=(8, 6))
                
                try:
                    sc.pl.spatial(
                        adata_sample,
                        color=marker,
                        ax=ax,
                        show=False,
                        cmap="Reds",
                        title=f"{sample} - {marker}"
                    )
                    
                    if self.save_plots:
                        fig_path = output_path / f"{sample}_{marker}.{self.format}"
                        fig.savefig(fig_path, dpi=self.dpi, bbox_inches="tight")
                    
                    if self.show_plots:
                        plt.show()
                    else:
                        plt.close()
                        
                except Exception as e:
                    warnings.warn(f"Failed to plot marker {marker}: {e}")
                    plt.close()
    
    def plot_niche_distribution(
        self,
        adata,
        celltype: str,
        output_dir: str,
        niche_col: str = "scNiche",
        sample_col: str = "sample_id"
    ):
        """
        Plot niche distribution
        
        Parameters:
            adata: AnnData object
            celltype: Cell type name
            output_dir: Output directory
            niche_col: Column name for niche labels
            sample_col: Column name for sample ID
        """
        try:
            import matplotlib.pyplot as plt
            import scanpy as sc
        except ImportError:
            warnings.warn("matplotlib or scanpy not installed")
            return
        
        if niche_col not in adata.obs:
            warnings.warn(f"Niche column '{niche_col}' not found")
            return
        
        celltype_safe = celltype.replace(" ", "_")
        output_path = Path(output_dir) / celltype_safe / "niche"
        output_path.mkdir(parents=True, exist_ok=True)
        
        samples = adata.obs[sample_col].unique()
        
        for sample in samples:
            sample_mask = adata.obs[sample_col] == sample
            adata_sample = adata[sample_mask].copy()
            
            if "spatial" not in adata_sample.obsm:
                continue
            
            fig, ax = plt.subplots(figsize=(10, 8))
            
            try:
                sc.pl.spatial(
                    adata_sample,
                    color=niche_col,
                    ax=ax,
                    show=False,
                    title=f"{sample} - Niche Distribution"
                )
                
                if self.save_plots:
                    fig_path = output_path / f"{sample}_niche.{self.format}"
                    fig.savefig(fig_path, dpi=self.dpi, bbox_inches="tight")
                
                if self.show_plots:
                    plt.show()
                else:
                    plt.close()
                    
            except (ValueError, RuntimeError, ImportError) as e:
                warnings.warn(f"Failed to plot niche for {sample}: {e}")
                plt.close()
    
    def plot_subtype_proportion(
        self,
        adata,
        celltype: str,
        output_dir: str,
        subtype_col: Optional[str] = None,
        sample_col: str = "sample_id"
    ):
        """
        Plot subtype proportion bar chart
        
        Parameters:
            adata: AnnData object
            celltype: Cell type name
            output_dir: Output directory
            subtype_col: Column name for subtype
            sample_col: Column name for sample ID
        """
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
        except ImportError:
            warnings.warn("matplotlib or seaborn not installed")
            return
        
        subtype_col = subtype_col or f"{celltype}_subtype"
        
        if subtype_col not in adata.obs:
            warnings.warn(f"Subtype column '{subtype_col}' not found")
            return
        
        celltype_safe = celltype.replace(" ", "_")
        output_path = Path(output_dir) / celltype_safe / "statistics"
        output_path.mkdir(parents=True, exist_ok=True)
        
        subtype_mask = adata.obs[subtype_col].notna()
        adata_subtype = adata[subtype_mask].copy()
        
        samples = adata_subtype.obs[sample_col].unique()
        
        proportion_data = []
        
        for sample in samples:
            sample_mask = adata_subtype.obs[sample_col] == sample
            sample_subtypes = adata_subtype.obs[sample_mask][subtype_col]
            
            counts = sample_subtypes.value_counts()
            proportions = counts / counts.sum()
            
            for subtype, prop in proportions.items():
                proportion_data.append({
                    "sample": sample,
                    "subtype": subtype,
                    "proportion": prop
                })
        
        proportion_df = pd.DataFrame(proportion_data)
        
        fig, ax = plt.subplots(figsize=(12, 6))
        
        try:
            sns.barplot(
                data=proportion_df,
                x="sample",
                y="proportion",
                hue="subtype",
                ax=ax
            )
            
            ax.set_title(f"{celltype} Subtype Proportion by Sample")
            ax.set_xlabel("Sample")
            ax.set_ylabel("Proportion")
            ax.legend(title="Subtype", bbox_to_anchor=(1.05, 1), loc="upper left")
            
            plt.tight_layout()
            
            if self.save_plots:
                fig_path = output_path / f"{celltype}_proportion.{self.format}"
                fig.savefig(fig_path, dpi=self.dpi, bbox_inches="tight")
            
            if self.show_plots:
                plt.show()
            else:
                plt.close()
                
        except (ValueError, RuntimeError, ImportError) as e:
            warnings.warn(f"Failed to plot proportion: {e}")
            plt.close()