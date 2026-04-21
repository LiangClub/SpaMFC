# inferCNVpy 功能集成规范

## 1. 项目概述

将 inferCNVpy 的 CNV 推断和可视化功能集成到 SpaMFC 项目中，增强对恶性细胞的 CNV 分析能力。

## 2. inferCNVpy 功能分析

### 2.1 核心功能

| 功能模块 | 函数 | 说明 |
|----------|------|------|
| **CNV推断** | `cnv.tl.infercnv()` | 从单细胞转录组数据推断CNV |
| **基因位置** | `cnv.io.genomic_position_from_gtf()` | 从GTF文件获取基因染色体位置 |
| **CNV评分** | `cnv.tl.cnv_score()` | 计算CNV评分 |
| **PCA降维** | `cnv.tl.pca()` | CNV数据PCA降维 |
| **聚类** | `cnv.tl.leiden()` | 基于CNV的Leiden聚类 |

### 2.2 可视化功能

| 函数 | 说明 |
|------|------|
| `cnv.pl.chromosome_heatmap()` | 染色体热图（显示所有细胞CNV模式） |
| `cnv.pl.chromosome_heatmap_summary()` | 染色体热图摘要（按组平均） |
| `cnv.pl.umap()` | CNV UMAP可视化 |
| `cnv.pl.tsne()` | CNV t-SNE可视化 |

### 2.3 数据结构

inferCNVpy 使用以下 AnnData 结构存储 CNV 数据：

- `adata.obsm['X_cnv']` - CNV矩阵（cells × genomic_bins）
- `adata.uns['cnv']['chr_pos']` - 染色体位置信息
- `adata.var['chromosome']` - 基因染色体注释
- `adata.var['start']`, `adata.var['end']` - 基因起始/结束位置

## 3. 集成方案

### 3.1 新增模块结构

```
src/
├── cnv_inference/              ← 新增CNV推断模块
│   ├── __init__.py
│   ├── infercnv.py             ← CNV推断核心功能
│   ├── genomic_positions.py    ← 基因位置处理
│   └── cnv_scoring.py          ← CNV评分计算
├── visualization/
│   ├── cnv_plot.py             ← 新增CNV可视化模块
│   └── __init__.py             ← 更新导出
```

### 3.2 功能设计

#### 3.2.1 CNV推断模块 (infercnv.py)

```python
class CNVInferencer:
    def infercnv(
        self,
        adata,
        reference_key: str,
        reference_cat: List[str],
        window_size: int = 100,
        step_size: int = 10,
        dynamic_threshold: float = 0.5,
        exclude_genes: Optional[List[str]] = None,
        ...
    ) -> AnnData
    
    def run_infercnv_pipeline(
        self,
        adata,
        gtf_file: str,
        reference_key: str,
        reference_cat: List[str],
        ...
    ) -> AnnData
```

#### 3.2.2 CNV可视化模块 (cnv_plot.py)

```python
class CNVVisualizer:
    def plot_chromosome_heatmap(
        self,
        adata,
        groupby: str = "cnv_leiden",
        use_rep: str = "cnv",
        cmap: str = "bwr",
        figsize: Tuple[int, int] = (16, 10),
        ...
    )
    
    def plot_chromosome_heatmap_summary(
        self,
        adata,
        groupby: str,
        ...
    )
    
    def plot_cnv_umap(
        self,
        adata,
        color: str,
        ...
    )
    
    def plot_cnv_scores(
        self,
        adata,
        score_key: str,
        ...
    )
```

### 3.3 CLI 子命令设计

新增 `cnv` 子命令：

```bash
SpaMFC cnv --input data.h5ad \
    --gtf genes.gtf \
    --reference-key "cell_type" \
    --reference-cat "Normal cells,Fibroblasts" \
    --output ./cnv_results/
```

### 3.4 配置更新

在 `SpaMFCConfig` 中添加 `CNVInferenceConfig`：

```python
@dataclass
class CNVInferenceConfig:
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
```

## 4. 任务分解

| 任务ID | 任务名称 | 优先级 | 说明 |
|--------|----------|--------|------|
| T-001 | 创建 cnv_inference 模块目录 | 高 | 新建 src/cnv_inference/ |
| T-002 | 实现 infercnv.py | 高 | CNV推断核心功能 |
| T-003 | 实现 genomic_positions.py | 高 | 基因位置处理 |
| T-004 | 实现 cnv_scoring.py | 中 | CNV评分计算 |
| T-005 | 实现 cnv_plot.py | 高 | CNV可视化模块 |
| T-006 | 更新 visualization/__init__.py | 中 | 导出CNVVisualizer |
| T-007 | 更新 config.py | 中 | 添加CNVInferenceConfig |
| T-008 | 更新 cli.py | 高 | 添加cnv子命令 |
| T-009 | 更新 src/__init__.py | 中 | 导出新模块 |
| T-010 | 更新 README.md | 低 | 文档更新 |
| T-011 | 语法验证 | 高 | py_compile验证 |

## 5. 技术规范

### 5.1 infercnv.py 核心实现

```python
import infercnvpy as cnv
import scanpy as sc
from typing import Optional, List

class CNVInferencer:
    def __init__(self, config: CNVInferenceConfig):
        self.config = config
    
    def add_genomic_positions(self, adata, gtf_file: str):
        cnv.io.genomic_position_from_gtf(gtf_file, adata)
        return adata
    
    def infercnv(
        self,
        adata,
        reference_key: str,
        reference_cat: List[str],
        **kwargs
    ):
        cnv.tl.infercnv(
            adata,
            reference_key=reference_key,
            reference_cat=reference_cat,
            window_size=self.config.window_size,
            step_size=self.config.step_size,
            **kwargs
        )
        return adata
    
    def compute_cnv_scores(self, adata):
        cnv.tl.cnv_score(adata)
        return adata
    
    def run_pipeline(self, adata, gtf_file: str, reference_key: str, reference_cat: List[str]):
        self.add_genomic_positions(adata, gtf_file)
        self.infercnv(adata, reference_key, reference_cat)
        if self.config.compute_scores:
            self.compute_cnv_scores(adata)
        cnv.tl.pca(adata, use_rep="cnv")
        cnv.pp.neighbors(adata, use_rep="cnv")
        cnv.tl.leiden(adata, resolution=self.config.clustering_resolution)
        return adata
```

### 5.2 cnv_plot.py 可视化实现

```python
import infercnvpy as cnv
import matplotlib.pyplot as plt
from typing import Optional, Tuple

class CNVVisualizer:
    def __init__(self, save_plots: bool = True, dpi: int = 300, format: str = "pdf"):
        self.save_plots = save_plots
        self.dpi = dpi
        self.format = format
    
    def plot_chromosome_heatmap(
        self,
        adata,
        groupby: str = "cnv_leiden",
        use_rep: str = "cnv",
        output_path: Optional[str] = None,
        **kwargs
    ):
        fig = cnv.pl.chromosome_heatmap(
            adata,
            groupby=groupby,
            use_rep=use_rep,
            show=False,
            **kwargs
        )
        if self.save_plots and output_path:
            fig.savefig(output_path, dpi=self.dpi, bbox_inches="tight")
        return fig
    
    def plot_chromosome_heatmap_summary(
        self,
        adata,
        groupby: str = "cnv_leiden",
        use_rep: str = "cnv",
        output_path: Optional[str] = None,
        **kwargs
    ):
        fig = cnv.pl.chromosome_heatmap_summary(
            adata,
            groupby=groupby,
            use_rep=use_rep,
            show=False,
            **kwargs
        )
        if self.save_plots and output_path:
            fig.savefig(output_path, dpi=self.dpi, bbox_inches="tight")
        return fig
    
    def plot_cnv_umap(
        self,
        adata,
        color: str = "cnv_leiden",
        output_path: Optional[str] = None,
        **kwargs
    ):
        fig = cnv.pl.umap(adata, color=color, show=False, **kwargs)
        if self.save_plots and output_path:
            fig.savefig(output_path, dpi=self.dpi, bbox_inches="tight")
        return fig
    
    def plot_cnv_scores(
        self,
        adata,
        score_key: str = "cnv_score",
        output_path: Optional[str] = None,
        **kwargs
    ):
        fig, ax = plt.subplots(figsize=(8, 6))
        sc.pl.violin(adata, score_key, groupby="cnv_leiden", ax=ax, show=False)
        if self.save_plots and output_path:
            fig.savefig(output_path, dpi=self.dpi, bbox_inches="tight")
        return fig
```

### 5.3 CLI 子命令参数

```python
def _add_cnv_subparser(subparsers):
    cnv_parser = subparsers.add_parser(
        "cnv",
        help="CNV inference and visualization",
        description="Perform CNV inference using inferCNVpy."
    )
    
    cnv_parser.add_argument("--input", "-i", required=True, help="Input h5ad file")
    cnv_parser.add_argument("--gtf", required=True, help="GTF file with gene positions")
    cnv_parser.add_argument("--reference-key", required=True, help="Reference key in adata.obs")
    cnv_parser.add_argument("--reference-cat", required=True, help="Reference categories (comma-separated)")
    cnv_parser.add_argument("--window-size", type=int, default=100)
    cnv_parser.add_argument("--step-size", type=int, default=10)
    cnv_parser.add_argument("--output", "-o", default="./cnv_results/")
    cnv_parser.add_argument("--plot-heatmap", action="store_true", default=True)
    cnv_parser.add_argument("--plot-umap", action="store_true", default=True)
    cnv_parser.add_argument("--resolution", type=float, default=0.5)
```

## 6. 验证标准

- 所有 Python 文件语法正确
- infercnvpy 正确导入和使用
- CLI 命令 `SpaMFC cnv --help` 正常工作
- CNV热图可视化正常生成
- 与现有SpaMFC流程兼容