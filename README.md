# SpaMFC: Spatial Multi-Feature Clustering

**空间转录组亚群精细分群分析工具**

[![Python](https://img.shields.io/badge/Python-3.10+-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Version](https://img.shields.io/badge/Version-0.0.1-orange.svg)](https://github.com/spamfc/SpaMFC)

SpaMFC (Spatial Multi-Feature Clustering) 是一个用于空间转录组数据亚群精细分群分析的Python工具包。通过整合多种特征（空间邻域特征、CNV特征、基因表达特征、生态位特征），实现准确且稳健的亚群识别与注释。

## 功能特点

- **多模态特征融合** - 整合空间位置、基因表达、CNV拷贝数变异等多种特征
- **灵活的特征选择** - 用户可自由选择使用哪些特征进行分析
- **三种融合方法** - 支持自适应权重、固定权重、直接拼接三种融合策略
- **分样本聚类** - 避免批次效应，每个样本独立处理
- **亚群功能注释** - 特征基因分析、功能富集分析、生态位分析
- **跨样本统一** - 多样本亚群识别后统一注释
- **完整CLI接口** - 提供命令行工具，方便批量分析

## 安装

### 从源码安装

```bash
git clone https://github.com/spamfc/SpaMFC.git
cd SpaMFC
pip install -e .
```

### 使用pip安装
```bash
pip install SpaMFC
```

### 安装可选依赖
```bash
# 功能富集分析
pip install gseapy>=1.0

# CNV分析（恶性细胞）
pip install infercnvpy>=0.4

# 生态位分析
pip install scNiche>=1.1

# Excel报告生成
pip install openpyxl>=3.0

# 安装所有可选依赖
pip install SpaMFC[all]
```

## 快速开始

### 命令行使用

```bash
# 查看帮助信息
spamfc_cli --help

# 查看数据信息
spamfc_cli info --input data.h5ad --celltype-col "anno_cell2location_res"

# CAFs分析（空间+表达特征）
spamfc_cli run --input data.h5ad \
    --celltype-col "anno_cell2location_res" \
    --celltype "CAFs" \
    --features spatial,expression \
    --fusion-method adaptive

# 恶性肿瘤细胞分析（所有特征+固定权重）
spamfc_cli run --input data.h5ad \
    --celltype-col "anno_cell2location_res" \
    --celltype "Malignant cells" \
    --features spatial,cnv,expression \
    --fusion-method fixed \
    --weights 0.3,0.3,0.4 \
    --markers EPCAM KRT8 KRT18

# 多细胞类型批量分析
spamfc_cli run-multi --input data.h5ad \
    --celltype-col "anno_cell2location_res" \
    --celltypes "Malignant cells,CAFs,ILC,ECs" \
    --features spatial,expression

# 生成配置文件
spamfc_cli config --output ./my_config.yaml --template cafs
```

### Python API使用

```python
from SpaMFC import SpaMFCPipeline

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
├── src/                          # 核心代码目录
│   ├── __init__.py               # 包初始化
│   ├── config.py                 # 配置管理模块
│   ├── pipeline.py               # 主流程入口
│   ├── cli.py                    # 命令行接口模块
│   ├── features/                 # 特征提取模块
│   │   ├── spatial.py            # 空间邻域特征
│   │   ├── cnv.py                # CNV特征
│   │   ├── expression.py         # 表达特征
│   │   └── niche.py              # 生态位特征
│   ├── fusion/                   # 特征融合模块
│   │   ├── weighting.py          # 权重计算与融合
│   │   ├── nmf.py                # NMF共识分解
│   │   └── clustering.py         # 聚类
│   ├── annotation/               # 功能注释模块
│   │   ├── markers.py            # 特征基因分析
│   │   ├── enrichment.py         # 功能富集
│   │   └── niche_analysis.py     # 生态位分析
│   ├── unification/              # 跨样本统一模块
│   │   ├── similarity.py         # 相似性计算
│   │   ├── fusion.py             # 相似度融合
│   │   └── mapping.py            # 亚型映射
│   └── visualization/            # 可视化模块
│       ├── spatial_plot.py       # 空间可视化
│       └── report.py             # 报告生成
├── configs/                      # 配置文件目录
│   ├── default_config.yaml       # 默认配置
│   ├── malignant_config.yaml     # 恶性细胞配置
│   ├── cafs_config.yaml          # CAFs配置
│   └── immune_config.yaml        # 免疫细胞配置
├── examples/                     # 示例代码
│   ├── example_malignant.py      # 恶性细胞分析示例
│   ├── example_cafs.py           # CAFs分析示例
│   └── example_multi.py          # 多细胞类型分析示例
├── spamfc_cli.py                 # 命令行入口
├── setup.py                      # 安装脚本
├── requirements.txt              # 依赖列表
├── README.md                     # 项目文档
└── LICENSE                        # 许可证
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

### 可选依赖
- gseapy >= 1.0 - 功能富集分析
- infercnvpy >= 0.4 - CNV分析（恶性细胞）
- scNiche >= 1.1 - 生态位分析
- openpyxl >= 3.0 - Excel报告生成

## 输出结果

分析完成后，输出目录包含以下内容：

- `*_subtype_annotated.h5ad` - 带有亚群注释的AnnData文件
- `spatial_subtype.pdf` - 空间分布可视化图
- `marker_genes.pdf` - 特征基因热图
- `subtype_proportion.pdf` - 亚群比例图
- `markers.csv` - 特征基因列表
- `weights.csv` - 特征权重
- `mapping.csv` - 亚群映射表

## 许可证

本项目采用 MIT 许可证。详见 [LICENSE](LICENSE) 文件。

## 贡献指南

欢迎提交Issue和Pull Request来改进项目。

## 联系方式

- 项目主页: https://github.com/spamfc/SpaMFC
- 问题反馈: https://github.com/spamfc/SpaMFC/issues