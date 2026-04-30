# SpaMFC: Spatial Multi-Feature Clustering

**空间转录组亚群精细分群分析工具**

[![Python](https://img.shields.io/badge/Python-3.10+-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Version](https://img.shields.io/badge/Version-2.2.0-orange.svg)](https://github.com/spamfc/SpaMFC)

SpaMFC (Spatial Multi-Feature Clustering) 是一个用于空间转录组数据亚群精细分群分析的Python工具包。通过整合多种特征（空间邻域特征、CNV特征、基因表达特征、生态位特征），实现准确且稳健的亚群识别与注释。同时提供基因相关性分析和CNV推断功能。

## 功能特点

- **多模态特征融合** - 整合空间位置、基因表达、CNV拷贝数变异等多种特征
- **灵活的特征选择** - 用户可自由选择使用哪些特征进行分析
- **三种融合方法** - 支持自适应权重、固定权重、直接拼接三种融合策略
- **分样本聚类** - 避免批次效应，每个样本独立处理
- **亚群功能注释** - 特征基因分析、功能富集分析、生态位分析
- **跨样本统一** - 多样本亚群识别后统一注释
- **基因相关性分析** - 支持Pearson、Spearman、Kendall三种方法，Numba加速
- **CNV推断分析** - 集成inferCNVpy，支持染色体热图、CNV聚类、CNV评分可视化
- **完整CLI接口** - 提供命令行工具，方便批量分析

## 安装

### 创建 Conda 环境（推荐）

```bash
# 创建新的conda环境
conda create -n spamfc python=3.10 -y
conda activate spamfc
```

### 从源码安装

```bash
git clone https://github.com/spamfc/SpaMFC.git
cd SpaMFC

# 安装
pip install -e .

# 或安装所有可选依赖
pip install -e ".[all]"
```

### 安装可选依赖

```bash
# 功能富集分析
pip install gseapy>=1.0

# CNV分析（恶性细胞）
pip install infercnvpy>=0.4

# Excel报告生成
pip install openpyxl>=3.0

# 安装所有可选依赖
pip install SpaMFC[all]
```

### 安装 scNiche（生态位分析）

scNiche 需要从 GitHub 安装，并依赖 PyTorch 和 DGL：

```bash
# 1. 安装 PyTorch（根据您的CUDA版本选择）
# CUDA 11.3 示例：
pip install torch==1.12.1+cu113 --extra-index-url https://download.pytorch.org/whl/cu113

# CUDA 12.1 示例：
pip install torch==2.1.0+cu121 --extra-index-url https://download.pytorch.org/whl/cu121

# CPU版本：
pip install torch

# 2. 安装 DGL（根据您的CUDA版本选择）
# CUDA 11.3 示例：
pip install dgl==1.1.0+cu113 -f https://data.dgl.ai/wheels/cu113/repo.html

# CUDA 12.1 示例：
pip install dgl==1.1.0+cu121 -f https://data.dgl.ai/wheels/cu121/repo.html

# CPU版本：
pip install dgl

# 3. 安装 scNiche
pip install scNiche @ git+https://github.com/ZJUFanLab/scNiche.git

# 或手动安装：
git clone https://github.com/ZJUFanLab/scNiche.git
cd scNiche
python setup.py build
python setup.py install
```

> **注意**: scNiche 需要 Python >= 3.9，建议使用 CUDA 环境以获得更好的性能。详细安装说明请参考 [scNiche GitHub](https://github.com/ZJUFanLab/scNiche)。

## 快速开始

### 命令行使用

```bash
# 查看帮助信息
SpaMFC --help

# 查看数据信息
SpaMFC info --input data.h5ad --celltype-col "anno_cell2location_res"

# CAFs分析（空间+表达特征）
SpaMFC run --input data.h5ad \
    --celltype-col "anno_cell2location_res" \
    --celltype "CAFs" \
    --features spatial,expression \
    --fusion-method adaptive

# 恶性肿瘤细胞分析（所有特征+固定权重）
SpaMFC run --input data.h5ad \
    --celltype-col "anno_cell2location_res" \
    --celltype "Malignant cells" \
    --features spatial,cnv,expression \
    --fusion-method fixed \
    --weights 0.3,0.3,0.4 \
    --markers EPCAM KRT8 KRT18

# 多细胞类型批量分析
SpaMFC run-multi --input data.h5ad \
    --celltype-col "anno_cell2location_res" \
    --celltypes "Malignant cells,CAFs,ILC,ECs" \
    --features spatial,expression

# 生成配置文件
SpaMFC config --output ./my_config.yaml --template cafs

# 基因相关性分析
SpaMFC corr --input data.h5ad \
    --target-genes EGFR,KRAS,TP53 \
    --de-genes GENE1,GENE2,GENE3 \
    --method spearman \
    --output ./correlation_results/

# CNV推断分析（恶性细胞）
SpaMFC cnv --input data.h5ad \
    --gtf genes.gtf \
    --reference-key "cell_type" \
    --reference-cat "Normal cells,Fibroblasts" \
    --output ./cnv_results/ \
    --plot-heatmap \
    --plot-umap
```

### Python API使用

```python
from SpaMFC import (
    SpaMFCPipeline, 
    gene_correlation, 
    CorrelationVisualizer,
    CNVInferencer,
    CNVVisualizer
)

# 创建分析流程
pipeline = SpaMFCPipeline("./configs/default_config.yaml")

# 设置特征（用户自由选择）
pipeline.set_feature_usage(
    use_spatial=True,
    use_cnv=False,
    use_expression=True,
    use_niche=False
)

# 分析单个细胞类型
adata = pipeline.run(adata, "CAFs")

# 批量分析多种细胞类型
adata = pipeline.run_multi_celltypes(adata, ["Malignant cells", "CAFs", "ILC"])

# 获取分析结果
results = pipeline.get_results("CAFs")
print(results["markers"])

# 基因相关性分析
corr_df, pval_df, sig_pairs = gene_correlation(
    adata,
    target_genes=["EGFR", "KRAS", "TP53"],
    de_genes=["GENE1", "GENE2", "GENE3"],
    method="spearman",
    threshold_p=0.05,
    output_dir="./correlation_results/"
)

# 可视化
visualizer = CorrelationVisualizer(adata, output_dir="./plots/")
visualizer.plot_single_pair_scatter("EGFR", "KRAS")
visualizer.plot_top_pairs_scatter_grid(sig_pairs, top_n=12)

# CNV推断分析
cnv_inferencer = CNVInferencer(
    window_size=100,
    step_size=10,
    clustering_resolution=0.5
)

adata = cnv_inferencer.run_pipeline(
    adata,
    gtf_file="genes.gtf",
    reference_key="cell_type",
    reference_cat=["Normal cells", "Fibroblasts"]
)

# CNV可视化
cnv_visualizer = CNVVisualizer(save_plots=True, dpi=300)
cnv_visualizer.plot_chromosome_heatmap(adata, groupby="cnv_leiden")
cnv_visualizer.plot_cnv_umap(adata, color="cnv_leiden")
cnv_visualizer.plot_cnv_scores(adata, score_key="cnv_score")
```

## 命令行参数说明

### run 子命令 - 单细胞类型分析

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--input, -i` | 输入h5ad文件路径 | 必填 |
| `--celltype-col` | 细胞类型列名 | 必填 |
| `--celltype, -c` | 目标细胞类型名称 | 必填 |
| `--output, -o` | 输出目录 | `./results/` |
| `--features` | 使用的特征，逗号分隔 | `spatial` |
| `--fusion-method` | 融合方法: adaptive/fixed/concat | `adaptive` |
| `--weights` | 固定权重，与features一一对应 | 自动均等 |
| `--n-clusters` | 聚类数量 | 5 |
| `--nmf-components` | NMF因子数 | 5 |
| `--clustering-method` | 聚类方法: kmeans/leiden | `kmeans` |
| `--markers` | 肿瘤标志物基因列表 | 无 |

### run-multi 子命令 - 多细胞类型分析

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--input, -i` | 输入h5ad文件路径 | 必填 |
| `--celltype-col` | 细胞类型列名 | 必填 |
| `--celltypes` | 细胞类型列表，逗号分隔 | 必填 |
| `--features` | 使用的特征 | `spatial` |
| `--fusion-method` | 融合方法 | `adaptive` |

### corr 子命令 - 基因相关性分析

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--input, -i` | 输入h5ad文件路径 | 必填 |
| `--target-genes` | 目标基因列表（逗号分隔或文件） | 必填 |
| `--de-genes` | 差异基因列表（逗号分隔或文件） | 必填 |
| `--method` | 相关性方法: pearson/spearman/kendall | `spearman` |
| `--p-adjust` | P值校正: fdr_bh/bonferroni/holm/none | `fdr_bh` |
| `--threshold-p` | P值阈值 | 0.05 |
| `--min-corr` | 最小相关性阈值 | 0.0 |
| `--output, -o` | 输出目录 | `./correlation_results/` |
| `--save-matrices` | 是否保存完整矩阵 | False |
| `--matrix-format` | 矩阵格式: npz/csv/csv.gz | `npz` |
| `--n-workers` | 并行进程数 | 自动 |
| `--max-memory` | 最大内存限制(MB) | 512 |

### cnv 子命令 - CNV推断分析

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--input, -i` | 输入h5ad文件路径 | 必填 |
| `--gtf` | GTF基因位置文件 | 必填 |
| `--reference-key` | 参考细胞类型列名 | 必填 |
| `--reference-cat` | 参考细胞类型（逗号分隔） | 必填 |
| `--output, -o` | 输出目录 | `./cnv_results/` |
| `--window-size` | CNV平滑窗口大小 | 100 |
| `--step-size` | CNV平滑步长 | 10 |
| `--resolution` | Leiden聚类分辨率 | 0.5 |
| `--n-pcs` | PCA组件数 | 30 |
| `--plot-heatmap` | 生成染色体热图 | True |
| `--plot-umap` | 生成CNV UMAP | True |
| `--plot-spatial` | 生成CNV空间分布图 | False |
| `--cmap` | 热图颜色映射 | `bwr` |

### info 子命令 - 显示数据信息

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--input, -i` | 输入h5ad文件路径 | 必填 |
| `--celltype-col` | 细胞类型列名 | 必填 |
| `--sample-col` | 样本ID列名 | `sample_id` |

### config 子命令 - 生成配置文件

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--output, -o` | 输出配置文件路径 | `./spamfc_config.yaml` |
| `--template` | 配置模板: default/malignant/cafs/immune | `default` |

## 特征说明

| 特征 | 说明 | 适用细胞类型 |
|------|------|--------------|
| `spatial` | 空间邻域特征 | 所有细胞类型 |
| `cnv` | CNV拷贝数变异特征 | 恶性肿瘤细胞 |
| `expression` | 基因表达特征 | 所有细胞类型 |
| `niche` | 生态位特征 | 所有细胞类型（可选） |

## 融合方法说明

| 方法 | 说明 |
|------|------|
| `adaptive` | 使用随机森林计算特征重要性，动态调整权重 |
| `fixed` | 使用用户指定的固定权重 |
| `concat` | 直接拼接特征，不做权重加权 |

## 基因相关性分析说明

### 相关性方法

| 方法 | 说明 | 性能 |
|------|------|------|
| `pearson` | 线性相关性，适用于正态分布数据 | Numba加速 ~20x |
| `spearman` | 秩相关性，适用于非线性关系 | Numba加速 ~10-15x |
| `kendall` | 秩相关性，适用于小样本 | 无加速 |

### 输出文件

| 文件 | 说明 |
|------|------|
| `significant_pairs.csv` | 显著相关对表格 |
| `statistics.json` | 分析统计信息 |
| `matrices.npz` | 相关性矩阵（可选） |
| `matrices_meta.json` | 矩阵元数据 |

## CNV推断分析说明

### 分析流程

1. **添加基因位置** - 从GTF文件获取基因染色体位置
2. **CNV推断** - 使用inferCNVpy推断CNV模式
3. **CNV评分** - 计算每个细胞的CNV评分
4. **降维聚类** - PCA降维 + Leiden聚类
5. **可视化** - 染色体热图、UMAP、空间分布

### 输出文件

| 文件 | 说明 |
|------|------|
| `cnv_annotated.h5ad` | 带CNV注释的AnnData文件 |
| `chromosome_heatmap.pdf` | 染色体热图 |
| `chromosome_heatmap_summary.pdf` | 染色体热图摘要 |
| `cnv_umap.pdf` | CNV UMAP可视化 |
| `cnv_score_distribution.pdf` | CNV评分分布图 |

## 配置文件说明

配置文件采用YAML格式，项目提供多个预设模板：

- `default_config.yaml` - 默认配置
- `malignant_config.yaml` - 恶性细胞配置（启用CNV）
- `cafs_config.yaml` - CAFs配置（空间+表达）
- `immune_config.yaml` - 免疫细胞配置（仅表达）

### 配置文件结构

```yaml
data:
  input_path: "./data/input.h5ad"
  output_dir: "./results/"
  sample_col: "sample_id"
  celltype_col: "anno_cell2location_res"

fusion_method: "adaptive"

features:
  spatial:
    use: true
    radius: 100
    k_neighbors: 15
    fixed_weight: 0.3
  cnv:
    use: false
    fixed_weight: 0.3
  expression:
    use: true
    fixed_weight: 0.4
  niche:
    use: false
    fixed_weight: 0.0

nmf:
  n_runs: 10
  n_components: 5

clustering:
  method: "kmeans"
  n_clusters: 5
  per_sample: true

unification:
  enable: true
  marker_weight: 0.4
  pathway_weight: 0.3
  niche_weight: 0.3
```

## 项目结构

```
SpaMFC/
├── SpaMFC/                      # 核心代码目录（Python包）
│   ├── __init__.py              # 包初始化
│   ├── config.py                # 配置管理模块
│   ├── pipeline.py              # 主流程入口
│   ├── cli.py                   # 命令行接口模块
│   ├── features/                # 特征提取模块
│   │   ├── spatial.py           # 空间邻域特征
│   │   ├── cnv.py               # CNV特征
│   │   ├── expression.py        # 表达特征
│   │   └── niche.py             # 生态位特征
│   ├── fusion/                  # 特征融合模块
│   │   ├── weighting.py         # 权重计算与融合
│   │   ├── nmf.py               # NMF共识分解
│   │   └── clustering.py        # 聚类
│   ├── annotation/              # 功能注释模块
│   │   ├── markers.py           # 特征基因分析
│   │   ├── enrichment.py        # 功能富集
│   │   └── niche_analysis.py    # 生态位分析
│   ├── correlation/             # 基因相关性分析模块
│   │   ├── gene_correlation.py  # 核心计算模块
│   │   └── correlation_plot.py  # 可视化模块
│   ├── cnv_inference/           # CNV推断模块
│   │   ├── infercnv.py          # CNV推断核心
│   │   ├── genomic_positions.py # 基因位置处理
│   │   └── cnv_scoring.py       # CNV评分计算
│   ├── unification/             # 跨样本统一模块
│   │   ├── similarity.py        # 相似性计算
│   │   ├── fusion.py            # 相似度融合
│   │   └── mapping.py           # 亚型映射
│   └── visualization/           # 可视化模块
│       ├── spatial_plot.py      # 空间可视化
│       ├── cnv_plot.py          # CNV可视化
│       └── report.py            # 报告生成
├── configs/                     # 配置文件目录
│   ├── default_config.yaml      # 默认配置
│   ├── malignant_config.yaml    # 恶性细胞配置
│   ├── cafs_config.yaml         # CAFs配置
│   └── immune_config.yaml       # 免疫细胞配置
├── examples/                    # 示例代码
│   ├── example_malignant.py     # 恶性细胞分析示例
│   ├── example_cafs.py          # CAFs分析示例
│   └── example_multi.py         # 多细胞类型分析示例
├── docker/                      # Docker配置
│   ├── Dockerfile               # GPU版本
│   ├── Dockerfile.cpu           # CPU版本
│   └── Dockerfile.cuda12.4      # CUDA 12.4版本
├── Dockerfile                   # 主Docker文件
├── setup.py                     # 安装脚本
├── requirements.txt             # 依赖列表
├── README.md                    # 项目文档
└── LICENSE                      # 许可证
```

## 依赖说明

### 核心依赖
- Python >= 3.10
- scanpy >= 1.10
- numpy >= 1.20
- pandas >= 2.0
- scikit-learn >= 1.0
- scipy >= 1.7
- matplotlib >= 3.5
- seaborn >= 0.12
- pyyaml >= 6.0
- umap-learn >= 0.5

### 基因相关性分析依赖
- numba >= 0.56 - JIT加速
- tqdm >= 4.60 - 进度条
- psutil >= 5.8 - 内存监控
- statsmodels >= 0.13 - P值校正

### CNV推断分析依赖
- infercnvpy >= 0.4 - CNV推断

### 可选依赖
- gseapy >= 1.0 - 功能富集分析
- scNiche >= 1.1 - 生态位分析
- openpyxl >= 3.0 - Excel报告生成

## 输出结果

### 亚群分析输出

分析完成后，输出目录包含以下内容：

- `*_subtype_annotated.h5ad` - 带有亚群注释的AnnData文件
- `spatial_subtype.pdf` - 空间分布可视化图
- `marker_genes.pdf` - 特征基因热图
- `subtype_proportion.pdf` - 亚群比例图
- `markers.csv` - 特征基因列表
- `weights.csv` - 特征权重
- `mapping.csv` - 亚群映射表

### 基因相关性分析输出

- `significant_pairs.csv` - 显著相关对表格
- `statistics.json` - 分析统计信息
- `matrices.npz` - 相关性矩阵（可选）
- `correlation_matrix.csv.gz` - 相关性矩阵（可选）
- `pvalue_matrix.csv.gz` - P值矩阵（可选）

### CNV推断分析输出

- `cnv_annotated.h5ad` - 带CNV注释的AnnData文件
- `chromosome_heatmap.pdf` - 染色体热图
- `chromosome_heatmap_summary.pdf` - 染色体热图摘要
- `cnv_score_distribution.pdf` - CNV评分分布
- `cnv_umap.pdf` - CNV UMAP可视化
- `spatial/` - CNV空间分布图（可选）

## 许可证

本项目采用 MIT 许可证。详见 [LICENSE](LICENSE) 文件。

## 贡献指南

欢迎提交Issue和Pull Request来改进项目。

## 联系方式

- 项目主页: https://github.com/spamfc/SpaMFC
- 问题反馈: https://github.com/spamfc/SpaMFC/issues