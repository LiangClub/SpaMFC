"""
SpaMFC Gene Correlation Visualization Module

Provides visualization functions for gene correlation analysis results:
- Single gene pair scatter plot
- Top N gene pairs grid plot
- Multiple gene pairs heatmap
- Correlation matrix heatmap
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.sparse as sp
from scipy.stats import pearsonr, spearmanr, kendalltau
from typing import Union, Tuple, Optional, List
import logging
import os


def plot_correlation(
    x: Union[list, np.ndarray, pd.Series, int, float],
    y: Union[list, np.ndarray, pd.Series, int, float],
    xlabel: str = 'X',
    ylabel: str = 'Y',
    title: str = 'Correlation Plot',
    method: str = 'pearson',
    trendline_color: str = 'red',
    trendline_style: str = '-',
    trendline_width: float = 2,
    scatter_color: str = 'black',
    scatter_size: float = 30,
    scatter_alpha: float = 1.0,
    ci_alpha: float = 0.2,
    ci_color: Optional[str] = None,
    font_scale: float = 1.0,
    title_fontsize: float = 12,
    label_fontsize: float = 10,
    tick_fontsize: float = 8,
    annot_fontsize: float = 10,
    figure_size: Tuple[float, float] = (3, 3),
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None,
    show_trendline: bool = True,
    show_ci: bool = True,
    save_path: Optional[str] = None,
    dpi: int = 300,
    tight_layout: bool = True
) -> Tuple[plt.Figure, plt.Axes]:
    def to_1d_array(data):
        arr = np.asarray(data)
        if arr.ndim == 0:
            arr = arr.reshape(1)
        elif arr.ndim > 1:
            raise ValueError(f"Input dimension error! Only 1D supported, got {arr.ndim}")
        return arr
    
    x_arr = to_1d_array(x)
    y_arr = to_1d_array(y)
    
    if len(x_arr) != len(y_arr):
        raise ValueError(f"X and Y length mismatch! X={len(x_arr)}, Y={len(y_arr)}")
    if len(x_arr) == 0:
        raise ValueError("Input data cannot be empty!")
    if not np.issubdtype(x_arr.dtype, np.number) or not np.issubdtype(y_arr.dtype, np.number):
        raise ValueError("X and Y must be numeric!")
    
    valid_methods = ['pearson', 'spearman', 'kendall']
    if method not in valid_methods:
        raise ValueError(f"method must be one of {valid_methods}, got {method}")
    
    valid_mask = ~(np.isnan(x_arr) | np.isnan(y_arr))
    valid_count = np.sum(valid_mask)
    
    if valid_count < 2:
        show_trendline = False
        show_ci = False
        if valid_count == 0:
            raise ValueError("No valid data points (all NaN)!")
        elif valid_count == 1:
            print("Warning: Only 1 valid point, cannot compute correlation")
    
    def get_times_font():
        import matplotlib.font_manager as fm
        font_names = ['Times New Roman', 'Times', 'DejaVu Serif']
        for name in font_names:
            if name in [f.name for f in fm.fontManager.ttflist]:
                return name
        return 'serif'
    times_font = get_times_font()
    
    sns.set_style("white", {"font_scale": font_scale})
    plt.rcParams.update({
        'font.family': times_font,
        'font.size': annot_fontsize,
        'axes.labelsize': label_fontsize,
        'axes.titlesize': title_fontsize,
        'xtick.labelsize': tick_fontsize,
        'ytick.labelsize': tick_fontsize,
        'axes.unicode_minus': False
    })
    
    fig, ax = plt.subplots(figsize=figure_size)
    
    if show_trendline and valid_count >= 2:
        ci = 95 if show_ci else None
        ci_color = ci_color or trendline_color
        
        sns.regplot(
            x=x_arr,
            y=y_arr,
            ax=ax,
            scatter=False,
            line_kws={
                'color': trendline_color,
                'linestyle': trendline_style,
                'linewidth': trendline_width
            },
            ci=ci,
            color=ci_color,
            scatter_kws={'alpha': 0}
        )
        
        for child in ax.get_children():
            if isinstance(child, plt.Polygon):
                child.set_alpha(ci_alpha)
                child.set_color(ci_color)
    
    ax.scatter(
        x_arr,
        y_arr,
        color=scatter_color,
        s=scatter_size,
        alpha=scatter_alpha,
        edgecolors='none'
    )
    
    def get_best_annot_pos(ax, x_data, y_data):
        candidate_pos = [
            (0.1, 0.9),
            (0.1, 0.1),
            (0.9, 0.9),
            (0.9, 0.1)
        ]
        
        def trans_rel_to_data(rel_pos):
            return ax.transAxes.transform(rel_pos)
        
        data_coords = np.column_stack((x_data, y_data))
        avg_distances = []
        for pos in candidate_pos:
            pos_data = trans_rel_to_data(pos)
            distances = np.sqrt(np.sum((data_coords - pos_data)**2, axis=1))
            avg_distances.append(np.mean(distances))
        
        best_idx = np.argmax(avg_distances)
        return candidate_pos[best_idx]
    
    if valid_count >= 2:
        x_valid = x_arr[valid_mask]
        y_valid = y_arr[valid_mask]
        
        if method == 'pearson':
            corr, p_value = pearsonr(x_valid, y_valid)
        elif method == 'spearman':
            corr, p_value = spearmanr(x_valid, y_valid)
        elif method == 'kendall':
            corr, p_value = kendalltau(x_valid, y_valid)
        
        if p_value < 1e-10:
            p_text = 'P < 1e-10'
        elif p_value < 0.001:
            p_text = 'P < 0.001'
        else:
            p_text = f'P = {p_value:.3f}'
        
        best_pos = get_best_annot_pos(ax, x_valid, y_valid)
        
        annot_text = f'{method.capitalize()} $r$ = {corr:.4f}\n{p_text}'
        ax.text(
            best_pos[0],
            best_pos[1],
            annot_text,
            transform=ax.transAxes,
            fontsize=annot_fontsize,
            fontstyle='italic',
            fontfamily=times_font,
            verticalalignment='center',
            horizontalalignment='center'
        )
    
    ax.set_xlabel(xlabel, fontfamily=times_font, fontsize=label_fontsize)
    ax.set_ylabel(ylabel, fontfamily=times_font, fontsize=label_fontsize)
    ax.set_title(title, fontfamily=times_font, fontsize=title_fontsize)
    
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    ax.grid(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if tight_layout:
        plt.tight_layout()
    
    if save_path is not None:
        fig.savefig(
            save_path,
            dpi=dpi,
            bbox_inches='tight' if tight_layout else None,
            facecolor='white',
            edgecolor='none'
        )
        print(f"Plot saved: {save_path}")
    
    return fig, ax


def sparse_to_1d_float(sparse_view: sp.spmatrix) -> np.ndarray:
    if not sp.issparse(sparse_view):
        raise TypeError(f"Input is not sparse! Type: {type(sparse_view)}")
    dense_arr = sparse_view.toarray()
    arr_1d = dense_arr.ravel()
    arr_float = arr_1d.astype(np.float64)
    arr_float[np.isinf(arr_float)] = np.nan
    return arr_float


def plot_correlation_heatmap(
    corr_matrix: pd.DataFrame,
    pval_matrix: pd.DataFrame,
    threshold_p: float = 0.05,
    figsize: Tuple[int, int] = (10, 8),
    cmap: str = "RdBu_r",
    title: str = "Target Genes vs DE Genes Correlation Heatmap",
    save_path: Optional[str] = None,
    dpi: int = 300
) -> plt.Figure:
    sig_mask = pval_matrix < threshold_p
    annot_text = np.empty_like(corr_matrix.values, dtype=object)
    for i in range(corr_matrix.shape[0]):
        for j in range(corr_matrix.shape[1]):
            corr = f"{corr_matrix.iloc[i, j]:.2f}"
            sig = "*" if sig_mask.iloc[i, j] else ""
            annot_text[i, j] = f"{corr}{sig}"

    plt.figure(figsize=figsize)
    sns.heatmap(
        corr_matrix,
        annot=annot_text,
        fmt="",
        cmap=cmap,
        vmin=-1,
        vmax=1,
        center=0,
        linewidths=0.5,
        cbar_kws={"label": "Correlation Coefficient"},
        square=True
    )
    plt.title(title, fontsize=14, pad=20)
    plt.xlabel("Differentially Expressed Genes", fontsize=12)
    plt.ylabel("Target Genes", fontsize=12)
    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0)
    plt.tight_layout()

    plt.text(
        0.02, 0.02, 
        "* P < 0.05", 
        transform=plt.gca().transAxes,
        fontsize=10,
        bbox=dict(facecolor="white", alpha=0.8)
    )

    if save_path is not None:
        plt.savefig(save_path, dpi=dpi, bbox_inches="tight")
        print(f"Heatmap saved: {save_path}")
    
    return plt.gcf()


class CorrelationVisualizer:
    def __init__(
        self,
        st_expr_matrix: Union[pd.DataFrame, np.ndarray],
        output_dir: str = "correlation_plots",
        logger: Optional[logging.Logger] = None
    ):
        self.st_expr_matrix = st_expr_matrix
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        
        if logger is None:
            logger = logging.getLogger(__name__)
            logger.setLevel(logging.INFO)
            handler = logging.StreamHandler()
            handler.setFormatter(logging.Formatter('%(message)s'))
            logger.addHandler(handler)
        self.logger = logger
        
        if hasattr(st_expr_matrix, 'obs_names') and hasattr(st_expr_matrix, 'var_names'):
            self.logger.info("[Data Prep] AnnData object detected")
            if sp.issparse(st_expr_matrix.X):
                expr = st_expr_matrix.X.toarray()
            else:
                expr = st_expr_matrix.X
            self.expr_df = pd.DataFrame(
                expr.T,
                index=[str(g).upper() for g in st_expr_matrix.var_names],
                columns=st_expr_matrix.obs_names
            )
        elif isinstance(st_expr_matrix, np.ndarray):
            self.expr_df = pd.DataFrame(st_expr_matrix)
        else:
            self.expr_df = st_expr_matrix.copy()
    
    def _get_gene_expression(self, gene_name: str) -> np.ndarray:
        if gene_name in self.expr_df.columns:
            expr = self.expr_df[gene_name].values
        elif gene_name in self.expr_df.index:
            expr = self.expr_df.loc[gene_name].values
        else:
            raise ValueError(f"Gene '{gene_name}' not in expression matrix")
        
        if sp.issparse(expr):
            expr = sparse_to_1d_float(expr)
        
        return expr
    
    def plot_single_pair_scatter(
        self,
        gene1: str,
        gene2: str,
        method: str = 'pearson',
        sample_size: Optional[int] = None,
        show_plot: bool = False,
        dpi: int = 300
    ) -> Tuple[plt.Figure, plt.Axes]:
        self.logger.info(f"[Scatter] Plotting: {gene1} vs {gene2}")
        
        expr1 = self._get_gene_expression(gene1)
        expr2 = self._get_gene_expression(gene2)
        
        if sample_size is not None and len(expr1) > sample_size:
            idx = np.random.choice(len(expr1), sample_size, replace=False)
            expr1 = expr1[idx]
            expr2 = expr2[idx]
            self.logger.info(f"[Scatter] Sampled: {sample_size} / {len(expr1)}")
        
        save_path = os.path.join(self.output_dir, f"scatter_{gene1}_{gene2}.png")
        fig, ax = plot_correlation(
            x=expr1,
            y=expr2,
            xlabel=gene1,
            ylabel=gene2,
            title=f"{gene1} vs {gene2}",
            method=method,
            save_path=save_path,
            dpi=dpi
        )
        
        if not show_plot:
            plt.close(fig)
        
        return fig, ax
    
    def plot_top_pairs_scatter_grid(
        self,
        sig_pairs_df: pd.DataFrame,
        top_n: int = 9,
        n_cols: int = 3,
        method: str = 'pearson',
        sample_size: Optional[int] = None,
        dpi: int = 300
    ) -> plt.Figure:
        self.logger.info(f"[Grid] Plotting Top {top_n} pairs")
        
        top_pairs = sig_pairs_df.head(top_n)
        
        n_rows = (top_n + n_cols - 1) // n_cols
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols*4, n_rows*4))
        if n_rows == 1 and n_cols == 1:
            axes = np.array([axes])
        axes = axes.flatten()
        
        for idx, (_, row) in enumerate(top_pairs.iterrows()):
            gene1, gene2 = row['target_gene'], row['de_gene']
            ax = axes[idx]
            
            expr1 = self._get_gene_expression(gene1)
            expr2 = self._get_gene_expression(gene2)
            
            if sample_size is not None and len(expr1) > sample_size:
                idx_rand = np.random.choice(len(expr1), sample_size, replace=False)
                expr1 = expr1[idx_rand]
                expr2 = expr2[idx_rand]
            
            valid_mask = ~(np.isnan(expr1) | np.isnan(expr2))
            if method == 'pearson':
                corr, pval = pearsonr(expr1[valid_mask], expr2[valid_mask])
            elif method == 'spearman':
                corr, pval = spearmanr(expr1[valid_mask], expr2[valid_mask])
            else:
                corr, pval = kendalltau(expr1[valid_mask], expr2[valid_mask])
            
            ax.scatter(expr1, expr2, s=10, alpha=0.5, color='black', edgecolors='none')
            
            from sklearn.linear_model import LinearRegression
            valid_idx = np.where(valid_mask)[0]
            if len(valid_idx) > 1:
                lr = LinearRegression()
                lr.fit(expr1[valid_idx].reshape(-1, 1), expr2[valid_idx])
                x_line = np.linspace(expr1.min(), expr1.max(), 100)
                y_line = lr.predict(x_line.reshape(-1, 1))
                ax.plot(x_line, y_line, color='red', linewidth=1.5)
            
            p_text = 'P<0.001' if pval < 0.001 else f'P={pval:.3f}'
            ax.text(0.05, 0.95, f'{method.capitalize()} r={corr:.2f}\n{p_text}',
                   transform=ax.transAxes, fontsize=9, verticalalignment='top')
            ax.set_xlabel(gene1, fontsize=10)
            ax.set_ylabel(gene2, fontsize=10)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
        
        for idx in range(top_n, len(axes)):
            axes[idx].axis('off')
        
        plt.tight_layout()
        
        save_path = os.path.join(self.output_dir, f"scatter_grid_top{top_n}.png")
        fig.savefig(save_path, dpi=dpi, bbox_inches='tight')
        self.logger.info(f"[Grid] Saved: {save_path}")
        plt.close(fig)
        
        return fig
    
    def plot_multiple_pairs_heatmap(
        self,
        gene_pairs: List[Tuple[str, str]],
        sig_pairs_df: Optional[pd.DataFrame] = None,
        corr_df: Optional[pd.DataFrame] = None,
        method: str = 'pearson',
        figsize: Tuple[int, int] = (10, 8),
        dpi: int = 300
    ) -> plt.Figure:
        self.logger.info(f"[Heatmap] Plotting {len(gene_pairs)} pairs")
        
        corrs = []
        for gene1, gene2 in gene_pairs:
            if sig_pairs_df is not None:
                mask = (sig_pairs_df['target_gene'] == gene1) & (sig_pairs_df['de_gene'] == gene2)
                if mask.any():
                    corr = sig_pairs_df.loc[mask, 'correlation'].values[0]
                    corrs.append(corr)
                    continue
            
            if corr_df is not None:
                if gene1 in corr_df.index and gene2 in corr_df.columns:
                    corr = corr_df.loc[gene1, gene2]
                    corrs.append(corr)
                    continue
                elif gene2 in corr_df.index and gene1 in corr_df.columns:
                    corr = corr_df.loc[gene2, gene1]
                    corrs.append(corr)
                    continue
            
            expr1 = self._get_gene_expression(gene1)
            expr2 = self._get_gene_expression(gene2)
            valid_mask = ~(np.isnan(expr1) | np.isnan(expr2))
            
            if method == 'pearson':
                corr, _ = pearsonr(expr1[valid_mask], expr2[valid_mask])
            elif method == 'spearman':
                corr, _ = spearmanr(expr1[valid_mask], expr2[valid_mask])
            else:
                corr, _ = kendalltau(expr1[valid_mask], expr2[valid_mask])
            
            corrs.append(corr)
        
        labels = [f"{g1}\nvs\n{g2}" for g1, g2 in gene_pairs]
        corr_matrix = pd.DataFrame([corrs], columns=labels)
        
        fig, ax = plt.subplots(figsize=figsize)
        sns.heatmap(
            corr_matrix,
            annot=True,
            fmt='.2f',
            cmap='RdBu_r',
            vmin=-1,
            vmax=1,
            center=0,
            ax=ax,
            cbar_kws={'label': 'Correlation'}
        )
        ax.set_title('Gene Pairs Correlation Heatmap', fontsize=14)
        ax.set_ylabel('Gene Pairs')
        plt.xticks(rotation=0)
        plt.tight_layout()
        
        save_path = os.path.join(self.output_dir, "heatmap_gene_pairs.png")
        fig.savefig(save_path, dpi=dpi, bbox_inches='tight')
        self.logger.info(f"[Heatmap] Saved: {save_path}")
        plt.close(fig)
        
        return fig
    
    def plot_correlation_matrix_heatmap(
        self,
        genes: List[str],
        corr_df: Optional[pd.DataFrame] = None,
        method: str = 'pearson',
        figsize: Tuple[int, int] = (10, 8),
        dpi: int = 300
    ) -> plt.Figure:
        self.logger.info(f"[Matrix Heatmap] Plotting {len(genes)} genes")
        
        if corr_df is not None:
            available_genes = [g for g in genes if g in corr_df.index and g in corr_df.columns]
            if len(available_genes) < 2:
                raise ValueError(f"Insufficient genes: {len(available_genes)}/{len(genes)}")
            corr_matrix = corr_df.loc[available_genes, available_genes]
        else:
            expr_matrix = np.column_stack([self._get_gene_expression(g) for g in genes])
            corr_matrix = pd.DataFrame(
                np.corrcoef(expr_matrix.T),
                index=genes,
                columns=genes
            )
        
        fig = plot_correlation_heatmap(
            corr_matrix=corr_matrix,
            pval_matrix=pd.DataFrame(np.ones_like(corr_matrix), index=corr_matrix.index, columns=corr_matrix.columns),
            figsize=figsize,
            title='Gene Correlation Matrix',
            save_path=os.path.join(self.output_dir, "heatmap_correlation_matrix.png"),
            dpi=dpi
        )
        plt.close(fig)
        
        return fig