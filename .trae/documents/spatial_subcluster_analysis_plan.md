# 空间转录组亚群精细分群分析项目计划（SpaMFC v2.0）

## 项目概述

本项目旨在开发一套**通用化、模块化**的空间转录组亚群精细分群分析流程，支持所有细胞类型的亚群识别，并提供完整的亚群功能注释和跨样本亚群统一注释功能。

---

## 一、核心设计理念

### 1.1 通用性设计
- **兼容所有细胞类型**：肿瘤细胞、CAFs、免疫细胞、内皮细胞、上皮细胞等
- **灵活特征选择**：用户可自由选择是否使用每种特征
- **模块化架构**：各分析步骤独立封装，可按需组合

### 1.2 用户特征选择设计

用户可以通过配置文件或命令行参数选择是否使用以下特征：

| 特征类型 | 参数名 | 默认值 | 适用细胞类型 | 说明 |
|---------|-------|-------|-------------|------|
| **空间邻域特征** | `use_spatial` | True | 所有类型 | 邻域细胞类型组成、边缘检测、中心距离 |
| **CNV特征** | `use_cnv` | False | 仅恶性细胞 | CNV拷贝数变异特征（需infercnvpy） |
| **基因表达特征** | `use_expression` | True | 所有类型 | PCA/UMAP降维后的表达特征 |
| **生态位特征** | `use_niche` | False | 所有类型 | scNiche多视图特征（需cell2location） |

**特征组合推荐**：
```
恶性肿瘤细胞: use_spatial=True, use_cnv=True, use_expression=True
CAFs/基质细胞: use_spatial=True, use_expression=True
免疫细胞: use_expression=True, use_spatial=False（可选）
内皮细胞: use_expression=True, use_spatial=False（可选）
上皮细胞: use_expression=True, use_spatial=False（可选）
```

### 1.3 三大核心功能
1. **亚群识别**：多模态特征融合 + NMF共识聚类
2. **亚群注释**：特征基因分析 + 功能富集 + 生态位分析
3. **跨样本统一**：基于特征基因、功能和生态位的亚群统一注释

---

## 二、项目代码结构

```
SpaMFC/
├── src/
│   ├── __init__.py
│   ├── config.py              # 配置管理模块
│   ├── features/
│   │   ├── __init__.py
│   │   ├── spatial.py         # 空间特征提取
│   │   ├── cnv.py             # CNV特征处理
│   │   ├── expression.py      # 表达特征处理
│   │   └── niche.py           # 生态位特征
│   ├── fusion/
│   │   ├── __init__.py
│   │   ├── weighting.py       # 自适应权重计算
│   │   ├── nmf.py             # NMF共识分解
│   │   └── clustering.py      # 聚类与稳定性验证
│   ├── annotation/
│   │   ├── __init__.py
│   │   ├── markers.py         # 特征基因分析
│   │   ├── enrichment.py      # 功能富集分析
│   │   └── niche_analysis.py  # 生态位分析
│   ├── unification/
│   │   ├── __init__.py
│   │   ├── similarity.py      # 多维度相似性计算
│   │   ├── fusion.py          # 相似度融合
│   │   └── mapping.py         # 亚型映射表生成
│   ├── visualization/
│   │   ├── __init__.py
│   │   ├── spatial_plot.py    # 空间可视化
│   │   └── report.py          # 报告生成
│   └── pipeline.py            # 主流程入口
├── configs/
│   ├── default_config.yaml    # 默认配置文件
│   ├── malignant_config.yaml  # 恶性细胞配置
│   ├── cafs_config.yaml       # CAFs配置
│   └── immune_config.yaml     # 免疫细胞配置
├── examples/
│   ├── example_malignant.py   # 恶性细胞分析示例
│   ├── example_cafs.py        # CAFs分析示例
│   └── example_multi.py       # 多细胞类型分析示例
├── tests/
│   ├── test_features.py
│   ├── test_fusion.py
│   └── test_annotation.py
├── requirements.txt
├── setup.py
└── README.md
```

---

## 三、配置系统设计

### 3.1 YAML配置文件格式

```yaml
# configs/default_config.yaml
# SpaMFC 默认配置文件

# ==================== 数据输入配置 ====================
data:
  input_path: "./data/input.h5ad"
  output_dir: "./results/"
  sample_col: "sample_id"
  celltype_col: "anno_cell2location_res"
  spatial_key: "spatial"
  cnv_key: "X_cnv"

# ==================== 特征选择配置（用户可自由选择） ====================
features:
  # 空间邻域特征
  spatial:
    use: true                    # 是否使用空间特征
    radius: 100                  # 邻域半径
    k_neighbors: 15              # 邻域细胞数
    weight_method: "adaptive"     # 权重方法: "adaptive" 或 "fixed"
    fixed_weight: 0.3            # 固定权重（仅当weight_method="fixed"时生效）
    
  # CNV特征（仅恶性细胞）
  cnv:
    use: false                   # 是否使用CNV特征
    pca_dim: 50                  # PCA降维维度
    use_umap: true               # 是否使用UMAP
    umap_dim: 2                  # UMAP维度
    marker_filter: true          # 是否进行标志物冗余过滤
    fixed_weight: 0.3            # 固定权重
    
  # 基因表达特征
  expression:
    use: true                    # 是否使用表达特征
    min_expr_pct: 0.1            # 低表达基因过滤阈值
    pca_dim: 50                  # PCA降维维度
    use_umap: true               # 是否使用UMAP
    umap_dim: 2                  # UMAP维度
    pca_keep_dim: 30             # 保留的PCA维度
    fixed_weight: 0.4            # 固定权重
    
  # 生态位特征（scNiche）
  niche:
    use: false                   # 是否使用生态位特征
    n_views: 3                   # 视图数量
    target_k: 15                 # 目标聚类数
    resolution: 0.3              # 聚类分辨率

# ==================== NMF配置 ====================
nmf:
  n_runs: 10                     # NMF运行次数
  n_components: 5                # NMF因子数
  init: "nndsvda"                # 初始化方法
  max_iter: 1000                 # 最大迭代次数
  tol: 1e-4                      # 收敛阈值

# ==================== 聚类配置 ====================
clustering:
  method: "kmeans"               # 聚类方法: "kmeans" 或 "leiden"
  n_clusters: 5                  # KMeans聚类数
  resolution: 0.1                # Leiden分辨率
  per_sample: true               # 是否分样本聚类
  n_bootstrap: 30                # Bootstrap验证次数
  stability_threshold: 0.7       # 稳定性阈值

# ==================== 特征基因分析配置 ====================
marker_analysis:
  method: "wilcoxon"             # 差异基因检验方法
  pval_threshold: 0.05           # p值阈值
  logfc_threshold: 0.25          # logFC阈值
  top_n: 50                      # 提取Top N标志基因

# ==================== 功能富集配置 ====================
enrichment:
  gene_sets:
    - "GO_Biological_Process_2023"
    - "KEGG_2021_Human"
  pval_threshold: 0.05           # 富集p值阈值
  top_n: 20                      # Top N富集结果

# ==================== 跨样本统一配置 ====================
unification:
  enable: true                   # 是否启用跨样本统一
  marker_weight: 0.4             # 特征基因相似度权重
  pathway_weight: 0.3            # 功能通路相似度权重
  niche_weight: 0.3              # 生态位相似度权重
  similarity_threshold: 0.6      # 统一阈值

# ==================== 可视化配置 ====================
visualization:
  save_plots: true               # 是否保存图像
  dpi: 300                       # 图像分辨率
  format: "pdf"                  # 图像格式
  show_plots: false              # 是否显示图像

# ==================== 运行配置 ====================
runtime:
  n_jobs: 2                      # 并行任务数
  verbose: true                  # 是否输出详细日志
  chunksize: 500                 # 数据分块大小
```

### 3.2 细胞类型特定配置

```yaml
# configs/malignant_config.yaml - 恶性肿瘤细胞配置
celltype: "Malignant cells"

features:
  spatial:
    use: true
    radius: 100
    k_neighbors: 15
  cnv:
    use: true                    # 恶性细胞启用CNV
    pca_dim: 50
    use_umap: true
  expression:
    use: true
    pca_dim: 50
    pca_keep_dim: 30

clustering:
  n_clusters: 5
  per_sample: true

# configs/cafs_config.yaml - CAFs配置
celltype: "CAFs"

features:
  spatial:
    use: true                    # CAFs启用空间特征
    radius: 100
    k_neighbors: 15
  cnv:
    use: false                   # CAFs不使用CNV
  expression:
    use: true
    pca_dim: 50

clustering:
  n_clusters: 4

# configs/immune_config.yaml - 免疫细胞配置
celltype: "ILC"

features:
  spatial:
    use: false                   # 免疫细胞可选空间特征
  cnv:
    use: false                   # 免疫细胞不使用CNV
  expression:
    use: true                    # 表达特征主导
    pca_dim: 50
    pca_keep_dim: 30

clustering:
  n_clusters: 3
```

---

## 四、核心模块代码设计

### 4.1 配置管理模块 (config.py)

```python
"""
配置管理模块
支持YAML配置文件加载、参数验证、用户特征选择
"""

import yaml
from pathlib import Path
from typing import Dict, Any, Optional
from dataclasses import dataclass, field

@dataclass
class FeatureConfig:
    """特征配置类"""
    use: bool = True
    # 空间特征参数
    radius: int = 100
    k_neighbors: int = 15
    weight_method: str = "adaptive"
    fixed_weight: float = 0.3
    # CNV参数
    pca_dim: int = 50
    use_umap: bool = True
    umap_dim: int = 2
    marker_filter: bool = True
    # 表达参数
    min_expr_pct: float = 0.1
    pca_keep_dim: int = 30

@dataclass
class SpaMFCConfig:
    """SpaMFC总配置类"""
    # 数据配置
    input_path: str = ""
    output_dir: str = "./results/"
    sample_col: str = "sample_id"
    celltype_col: str = "anno_cell2location_res"
    
    # 特征选择（用户可配置）
    use_spatial: bool = True
    use_cnv: bool = False
    use_expression: bool = True
    use_niche: bool = False
    
    # 各特征详细配置
    spatial_config: FeatureConfig = field(default_factory=FeatureConfig)
    cnv_config: FeatureConfig = field(default_factory=FeatureConfig)
    expression_config: FeatureConfig = field(default_factory=FeatureConfig)
    
    # NMF配置
    nmf_runs: int = 10
    nmf_components: int = 5
    
    # 聚类配置
    clustering_method: str = "kmeans"
    n_clusters: int = 5
    per_sample: bool = True
    
    # 功能注释配置
    enable_marker_analysis: bool = True
    enable_enrichment: bool = True
    enable_niche_analysis: bool = False
    
    # 跨样本统一配置
    enable_unification: bool = True
    
    # 运行配置
    verbose: bool = True
    n_jobs: int = 2

class ConfigManager:
    """配置管理器"""
    
    def __init__(self, config_path: Optional[str] = None):
        self.config = SpaMFCConfig()
        if config_path:
            self.load_from_yaml(config_path)
    
    def load_from_yaml(self, config_path: str):
        """从YAML文件加载配置"""
        with open(config_path, 'r') as f:
            yaml_config = yaml.safe_load(f)
        self._parse_yaml_config(yaml_config)
    
    def _parse_yaml_config(self, yaml_config: Dict):
        """解析YAML配置"""
        # 数据配置
        data_cfg = yaml_config.get('data', {})
        self.config.input_path = data_cfg.get('input_path', '')
        self.config.output_dir = data_cfg.get('output_dir', './results/')
        
        # 特征选择（核心：用户可自由选择）
        features_cfg = yaml_config.get('features', {})
        
        # 空间特征
        spatial_cfg = features_cfg.get('spatial', {})
        self.config.use_spatial = spatial_cfg.get('use', True)
        self.config.spatial_config = FeatureConfig(
            use=spatial_cfg.get('use', True),
            radius=spatial_cfg.get('radius', 100),
            k_neighbors=spatial_cfg.get('k_neighbors', 15),
            weight_method=spatial_cfg.get('weight_method', 'adaptive'),
            fixed_weight=spatial_cfg.get('fixed_weight', 0.3)
        )
        
        # CNV特征
        cnv_cfg = features_cfg.get('cnv', {})
        self.config.use_cnv = cnv_cfg.get('use', False)
        self.config.cnv_config = FeatureConfig(
            use=cnv_cfg.get('use', False),
            pca_dim=cnv_cfg.get('pca_dim', 50),
            use_umap=cnv_cfg.get('use_umap', True),
            umap_dim=cnv_cfg.get('umap_dim', 2),
            marker_filter=cnv_cfg.get('marker_filter', True),
            fixed_weight=cnv_cfg.get('fixed_weight', 0.3)
        )
        
        # 表达特征
        expr_cfg = features_cfg.get('expression', {})
        self.config.use_expression = expr_cfg.get('use', True)
        self.config.expression_config = FeatureConfig(
            use=expr_cfg.get('use', True),
            min_expr_pct=expr_cfg.get('min_expr_pct', 0.1),
            pca_dim=expr_cfg.get('pca_dim', 50),
            use_umap=expr_cfg.get('use_umap', True),
            umap_dim=expr_cfg.get('umap_dim', 2),
            pca_keep_dim=expr_cfg.get('pca_keep_dim', 30),
            fixed_weight=expr_cfg.get('fixed_weight', 0.4)
        )
        
        # 生态位特征
        niche_cfg = features_cfg.get('niche', {})
        self.config.use_niche = niche_cfg.get('use', False)
        
        # NMF配置
        nmf_cfg = yaml_config.get('nmf', {})
        self.config.nmf_runs = nmf_cfg.get('n_runs', 10)
        self.config.nmf_components = nmf_cfg.get('n_components', 5)
        
        # 聚类配置
        cluster_cfg = yaml_config.get('clustering', {})
        self.config.clustering_method = cluster_cfg.get('method', 'kmeans')
        self.config.n_clusters = cluster_cfg.get('n_clusters', 5)
        self.config.per_sample = cluster_cfg.get('per_sample', True)
        
        # 功能注释配置
        marker_cfg = yaml_config.get('marker_analysis', {})
        self.config.enable_marker_analysis = marker_cfg.get('enable', True)
        
        enrich_cfg = yaml_config.get('enrichment', {})
        self.config.enable_enrichment = enrich_cfg.get('enable', True)
        
        # 跨样本统一配置
        unif_cfg = yaml_config.get('unification', {})
        self.config.enable_unification = unif_cfg.get('enable', True)
    
    def set_feature_usage(self, use_spatial: bool, use_cnv: bool, 
                          use_expression: bool, use_niche: bool):
        """用户动态设置特征使用"""
        self.config.use_spatial = use_spatial
        self.config.use_cnv = use_cnv
        self.config.use_expression = use_expression
        self.config.use_niche = use_niche
    
    def validate_config(self) -> bool:
        """验证配置有效性"""
        # 至少需要一种特征
        if not any([self.config.use_spatial, self.config.use_cnv, 
                    self.config.use_expression, self.config.use_niche]):
            raise ValueError("至少需要启用一种特征类型")
        
        # CNV特征需要CNV数据
        if self.config.use_cnv:
            if self.config.cnv_key not in self.config.required_keys:
                raise ValueError("使用CNV特征需要提供CNV数据")
        
        return True
    
    def get_active_features(self) -> list:
        """获取当前启用的特征列表"""
        active = []
        if self.config.use_spatial:
            active.append('spatial')
        if self.config.use_cnv:
            active.append('cnv')
        if self.config.use_expression:
            active.append('expression')
        if self.config.use_niche:
            active.append('niche')
        return active
```

### 4.2 特征提取模块设计

#### spatial.py - 空间特征提取
```python
"""
空间邻域特征提取模块
用户可通过配置选择是否使用
"""

import numpy as np
import pandas as pd
from sklearn.neighbors import KDTree
from sklearn.preprocessing import MinMaxScaler
from typing import Optional
import scanpy as sc

class SpatialFeatureExtractor:
    """空间特征提取器"""
    
    def __init__(self, config: FeatureConfig):
        self.config = config
        self.radius = config.radius
        self.k_neighbors = config.k_neighbors
    
    def extract(self, adata: sc.AnnData, target_cells: list,
                celltype_col: str, global_cell_types: list) -> pd.DataFrame:
        """
        提取空间邻域特征
        
        参数:
            adata: AnnData对象
            target_cells: 目标细胞列表
            celltype_col: 细胞类型列名
            global_cell_types: 全局细胞类型列表
        
        返回:
            空间特征DataFrame
        """
        if not self.config.use:
            return None
        
        coords = adata.obsm['spatial']
        tree = KDTree(coords)
        sigma = self.radius / 2
        
        env_features = []
        
        # 计算肿瘤拓扑特征（如果目标细胞是肿瘤细胞）
        target_mask = adata.obs_names.isin(target_cells)
        target_coords = coords[target_mask]
        tumor_center = target_coords.mean(axis=0) if len(target_coords) > 0 else np.zeros(2)
        
        cell_center_dist = np.linalg.norm(
            target_coords - tumor_center, axis=1
        )
        cell_center_dist = MinMaxScaler().fit_transform(
            cell_center_dist.reshape(-1, 1)
        ).flatten()
        
        for idx, cell in enumerate(target_cells):
            cell_idx = adata.obs_names.get_loc(cell)
            
            # KDTree邻域查询
            distances, neighbor_indices = tree.query(
                coords[cell_idx].reshape(1, -1), 
                k=self.k_neighbors
            )
            distances = distances[0]
            neighbor_indices = neighbor_indices[0]
            
            # 距离加权
            weights = np.exp(-(distances**2) / (2*sigma**2))
            weights = weights / weights.sum()
            
            # 邻域细胞类型组成
            neighbor_types = adata.obs.iloc[neighbor_indices][celltype_col]
            type_weights = pd.Series(weights, index=neighbor_indices)
            weighted_counts = type_weights.groupby(neighbor_types).sum()
            
            # 构建全局特征向量
            env_vec = pd.Series(0.0, index=global_cell_types)
            env_vec[weighted_counts.index] = weighted_counts.values
            
            # 拓扑特征
            neighbor_cell_types = adata.obs.iloc[neighbor_indices][celltype_col]
            is_edge = (neighbor_cell_types != celltype_col).any()
            env_vec['is_edge'] = 1.0 if is_edge else 0.0
            env_vec['distance_to_center'] = cell_center_dist[idx]
            
            env_features.append(env_vec)
        
        return pd.DataFrame(env_features, index=target_cells)
```

#### cnv.py - CNV特征处理
```python
"""
CNV特征处理模块
仅适用于恶性细胞，用户可选择是否使用
"""

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import MinMaxScaler, StandardScaler
import umap
from scipy import sparse
import scanpy as sc

class CNVFeatureProcessor:
    """CNV特征处理器"""
    
    def __init__(self, config: FeatureConfig):
        self.config = config
        self.pca_dim = config.pca_dim
        self.use_umap = config.use_umap
        self.umap_dim = config.umap_dim
        self.marker_filter = config.marker_filter
    
    def process(self, adata: sc.AnnData, target_cells: list,
                cnv_key: str = 'X_cnv',
                markers: list = None) -> pd.DataFrame:
        """
        处理CNV特征
        
        参数:
            adata: AnnData对象
            target_cells: 目标细胞列表
            cnv_key: CNV矩阵键名
            markers: 肿瘤标志物列表（用于冗余过滤）
        
        返回:
            CNV特征DataFrame
        """
        if not self.config.use:
            return None
        
        # 提取CNV矩阵
        cnv_data = adata.obsm[cnv_key]
        if sparse.issparse(cnv_data):
            cnv_data = cnv_data.toarray()
        cnv_data = np.asarray(cnv_data, dtype=np.float64)
        
        # 筛选目标细胞
        target_mask = adata.obs_names.isin(target_cells)
        cnv_target = cnv_data[target_mask]
        
        # 构建DataFrame
        cnv_df = pd.DataFrame(
            cnv_target,
            index=target_cells,
            columns=[f'CNV_bin_{i}' for i in range(cnv_target.shape[1])]
        )
        
        # 填充与裁剪异常值
        cnv_df = cnv_df.fillna(0.0)
        cnv_clipped = np.clip(cnv_df.values, -2.0, 2.0)
        
        # 标志物冗余过滤（可选）
        if self.marker_filter and markers:
            cnv_df_scaled = StandardScaler().fit_transform(cnv_clipped)
            # 计算CNV与标志物表达的相关性
            # ...（实现标志物过滤逻辑）
        
        # PCA降维
        n_pca = min(self.pca_dim, len(target_cells)-1, cnv_df.shape[1])
        pca = PCA(n_components=n_pca, random_state=0)
        cnv_pca = pca.fit_transform(cnv_df_scaled if self.marker_filter else cnv_clipped)
        
        # UMAP降维（可选）
        if self.use_umap:
            umap_model = umap.UMAP(
                n_components=self.umap_dim,
                random_state=0,
                n_neighbors=5
            )
            cnv_umap = umap_model.fit_transform(cnv_pca)
            cnv_dim_df = pd.DataFrame(
                cnv_umap,
                index=target_cells,
                columns=[f'CNV_UMAP{i+1}' for i in range(self.umap_dim)]
            )
        else:
            cnv_dim_df = pd.DataFrame(
                cnv_pca[:, :min(10, n_pca)],
                index=target_cells,
                columns=[f'CNV_PC{i+1}' for i in range(min(10, n_pca))]
            )
        
        return cnv_dim_df
```

#### expression.py - 表达特征处理
```python
"""
基因表达特征处理模块
所有细胞类型都可使用，用户可选择是否启用
"""

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import MinMaxScaler
import umap
import scanpy as sc

class ExpressionFeatureProcessor:
    """表达特征处理器"""
    
    def __init__(self, config: FeatureConfig):
        self.config = config
        self.min_expr_pct = config.min_expr_pct
        self.pca_dim = config.pca_dim
        self.use_umap = config.use_umap
        self.umap_dim = config.umap_dim
        self.pca_keep_dim = config.pca_keep_dim
    
    def process(self, adata: sc.AnnData, target_cells: list) -> pd.DataFrame:
        """
        处理基因表达特征
        
        参数:
            adata: AnnData对象
            target_cells: 目标细胞列表
        
        返回:
            表达特征DataFrame
        """
        if not self.config.use:
            return None
        
        # 提取目标细胞子集
        adata_target = adata[target_cells, :].copy()
        
        # 低表达基因过滤
        min_cells = int(len(adata_target) * self.min_expr_pct)
        sc.pp.filter_genes(adata_target, min_cells=min_cells)
        
        if adata_target.n_vars == 0:
            # 兜底：保留前100个基因
            adata_target = adata[target_cells, :min(100, adata.n_vars)].copy()
        
        # 提取表达矩阵
        expr_data = adata_target.X
        if hasattr(expr_data, 'toarray'):
            expr_data = expr_data.toarray()
        expr_data = np.asarray(expr_data, dtype=np.float64)
        
        expr_df = pd.DataFrame(
            expr_data,
            index=target_cells,
            columns=adata_target.var_names
        )
        
        # PCA降维
        n_pca = min(self.pca_dim, len(target_cells)-1, expr_df.shape[1])
        pca = PCA(n_components=n_pca, random_state=0)
        expr_pca = pca.fit_transform(expr_df.values)
        
        # 构建特征DataFrame
        expr_dim_df = pd.DataFrame()
        
        # UMAP降维（可选）
        if self.use_umap:
            umap_model = umap.UMAP(
                n_components=self.umap_dim,
                random_state=0,
                n_neighbors=5
            )
            expr_umap = umap_model.fit_transform(expr_pca)
            expr_umap_df = pd.DataFrame(
                expr_umap,
                index=target_cells,
                columns=[f'Expr_UMAP{i+1}' for i in range(self.umap_dim)]
            )
            expr_dim_df = pd.concat([expr_dim_df, expr_umap_df], axis=1)
        
        # 保留PCA特征
        pca_keep = expr_pca[:, :min(self.pca_keep_dim, n_pca)]
        expr_pca_df = pd.DataFrame(
            pca_keep,
            index=target_cells,
            columns=[f'Expr_PC{i+1}' for i in range(pca_keep.shape[1])]
        )
        expr_dim_df = pd.concat([expr_dim_df, expr_pca_df], axis=1)
        
        return expr_dim_df
```

### 4.3 特征融合与聚类模块设计

#### weighting.py - 自适应权重计算
```python
"""
自适应权重计算模块
根据特征重要性动态调整各特征权重
"""

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import MinMaxScaler

class AdaptiveWeightCalculator:
    """自适应权重计算器"""
    
    def __init__(self, weight_method: str = 'adaptive'):
        self.weight_method = weight_method
    
    def calculate(self, feature_dfs: dict, fixed_weights: dict = None) -> dict:
        """
        计算特征权重
        
        参数:
            feature_dfs: 各特征DataFrame字典 {'spatial': df, 'cnv': df, 'expression': df}
            fixed_weights: 固定权重字典
        
        返回:
            权重字典 {'spatial': weight, 'cnv': weight, 'expression': weight}
        """
        if self.weight_method == 'fixed':
            return fixed_weights or {'spatial': 0.3, 'cnv': 0.3, 'expression': 0.4}
        
        # 自适应权重计算
        active_features = [k for k, v in feature_dfs.items() if v is not None]
        
        if len(active_features) == 1:
            return {active_features[0]: 1.0}
        
        # 拼接特征
        fused_df = pd.concat([feature_dfs[k] for k in active_features], axis=1)
        
        # 生成伪标签（基于空间距离）
        if 'spatial' in feature_dfs:
            spatial_df = feature_dfs['spatial']
            if 'distance_to_center' in spatial_df.columns:
                distance = spatial_df['distance_to_center']
                pseudo_labels = (distance > distance.median()).astype(int).values
            else:
                pseudo_labels = np.random.randint(0, 2, len(fused_df))
        else:
            pseudo_labels = np.random.randint(0, 2, len(fused_df))
        
        # 训练随机森林
        rf = RandomForestClassifier(n_estimators=200, random_state=0)
        rf.fit(fused_df.values, pseudo_labels)
        
        # 计算各特征重要性
        feature_importances = rf.feature_importances_
        
        # 划分各特征的重要性区间
        weights = {}
        start_idx = 0
        for feature_name in active_features:
            feature_dim = feature_dfs[feature_name].shape[1]
            end_idx = start_idx + feature_dim
            importance = np.mean(feature_importances[start_idx:end_idx])
            weights[feature_name] = importance
            start_idx = end_idx
        
        # 归一化权重
        total_importance = sum(weights.values())
        if total_importance > 0:
            weights = {k: v/total_importance for k, v in weights.items()}
        else:
            weights = {k: 1.0/len(active_features) for k in active_features}
        
        return weights
```

#### nmf.py - NMF共识分解
```python
"""
NMF共识分解模块
多次运行NMF取共识嵌入，提高稳定性
"""

import numpy as np
import pandas as pd
from sklearn.decomposition import NMF
from sklearn.preprocessing import MinMaxScaler

class NMFConsensus:
    """NMF共识分解器"""
    
    def __init__(self, n_runs: int = 10, n_components: int = 5):
        self.n_runs = n_runs
        self.n_components = n_components
    
    def decompose(self, fused_df: pd.DataFrame) -> pd.DataFrame:
        """
        NMF共识分解
        
        参数:
            fused_df: 融合特征DataFrame
        
        返回:
            NMF共识嵌入DataFrame
        """
        # 标准化特征
        scaler = MinMaxScaler()
        fused_scaled = scaler.fit_transform(fused_df.values)
        
        # 多次NMF运行
        nmf_embeds = []
        for i in range(self.n_runs):
            model = NMF(
                n_components=self.n_components,
                init='nndsvda',
                random_state=i,
                max_iter=1000,
                tol=1e-4
            )
            W = model.fit_transform(fused_scaled)
            nmf_embeds.append(W)
        
        # 取共识嵌入
        W_consensus = np.mean(np.stack(nmf_embeds, axis=2), axis=2)
        
        # 构建DataFrame
        W_df = pd.DataFrame(
            W_consensus,
            index=fused_df.index,
            columns=[f'NMF_Factor_{i+1}' for i in range(self.n_components)]
        )
        
        return W_df
```

#### clustering.py - 聚类与稳定性验证
```python
"""
聚类与稳定性验证模块
支持KMeans和Leiden聚类，分样本处理
"""

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.utils import resample
import scanpy as sc

class SubtypeClusterer:
    """亚型聚类器"""
    
    def __init__(self, method: str = 'kmeans', n_clusters: int = 5,
                 resolution: float = 0.1, per_sample: bool = True,
                 n_bootstrap: int = 30, stability_threshold: float = 0.7):
        self.method = method
        self.n_clusters = n_clusters
        self.resolution = resolution
        self.per_sample = per_sample
        self.n_bootstrap = n_bootstrap
        self.stability_threshold = stability_threshold
    
    def cluster(self, adata: sc.AnnData, W_df: pd.DataFrame,
               sample_col: str, celltype: str) -> sc.AnnData:
        """
        执行聚类
        
        参数:
            adata: AnnData对象
            W_df: NMF嵌入DataFrame
            sample_col: 样本列名
            celltype: 细胞类型名称
        
        返回:
            更新后的adata
        """
        subtype_col = f'{celltype}_subtype'
        adata.obs[subtype_col] = np.nan
        
        samples = adata.obs[sample_col].unique()
        
        if self.per_sample:
            # 分样本聚类
            for sample in samples:
                sample_cells = [c for c in W_df.index 
                               if adata.obs.loc[c, sample_col] == sample]
                if not sample_cells:
                    continue
                
                W_sample = W_df.loc[sample_cells]
                
                # 聚类
                if self.method == 'kmeans':
                    kmeans = KMeans(n_clusters=self.n_clusters, random_state=0)
                    labels = kmeans.fit_predict(W_sample.values)
                else:
                    adata_nmf = sc.AnnData(W_sample)
                    sc.pp.neighbors(adata_nmf, n_neighbors=10)
                    sc.tl.leiden(adata_nmf, resolution=self.resolution)
                    labels = adata_nmf.obs['leiden'].values
                
                # 写回adata
                for cell, label in zip(sample_cells, labels):
                    adata.obs.loc[cell, subtype_col] = f'{sample}_{label}'
        else:
            # 全局聚类
            if self.method == 'kmeans':
                kmeans = KMeans(n_clusters=self.n_clusters, random_state=0)
                labels = kmeans.fit_predict(W_df.values)
            else:
                adata_nmf = sc.AnnData(W_df)
                sc.pp.neighbors(adata_nmf, n_neighbors=10)
                sc.tl.leiden(adata_nmf, resolution=self.resolution)
                labels = adata_nmf.obs['leiden'].values
            
            # 写回adata
            for cell, label in zip(W_df.index, labels):
                sample = adata.obs.loc[cell, sample_col]
                adata.obs.loc[cell, subtype_col] = f'{sample}_{label}'
        
        return adata
    
    def validate_stability(self, W_df: pd.DataFrame) -> list:
        """
        验证聚类稳定性
        
        返回:
            稳定亚型列表
        """
        cluster_freq = {}
        
        for i in range(self.n_bootstrap):
            # Bootstrap重采样
            sample_indices = resample(
                range(len(W_df)), 
                replace=True,
                n_samples=int(0.8*len(W_df))
            )
            W_sample = W_df.iloc[sample_indices]
            
            # 聚类
            if self.method == 'kmeans':
                kmeans = KMeans(n_clusters=self.n_clusters, random_state=i)
                labels = kmeans.fit_predict(W_sample.values)
            else:
                adata_nmf = sc.AnnData(W_sample)
                sc.pp.neighbors(adata_nmf, n_neighbors=10)
                sc.tl.leiden(adata_nmf, resolution=self.resolution, random_state=i)
                labels = adata_nmf.obs['leiden'].values
            
            # 统计频率
            for label in set(labels):
                cluster_freq[label] = cluster_freq.get(label, 0) + 1
        
        # 稳定亚型筛选
        freq_threshold = self.n_bootstrap * self.stability_threshold
        stable_labels = [label for label, freq in cluster_freq.items() 
                        if freq >= freq_threshold]
        
        return stable_labels
```

### 4.4 亚群功能注释模块设计

#### markers.py - 特征基因分析
```python
"""
特征基因分析模块
Wilcoxon差异基因检验
"""

import pandas as pd
import scanpy as sc
from typing import Dict, List

class MarkerGeneAnalyzer:
    """特征基因分析器"""
    
    def __init__(self, method: str = 'wilcoxon',
                 pval_threshold: float = 0.05,
                 logfc_threshold: float = 0.25,
                 top_n: int = 50):
        self.method = method
        self.pval_threshold = pval_threshold
        self.logfc_threshold = logfc_threshold
        self.top_n = top_n
    
    def analyze(self, adata: sc.AnnData, celltype: str,
               sample_col: str) -> Dict[str, List[str]]:
        """
        分析特征基因
        
        参数:
            adata: AnnData对象
            celltype: 细胞类型名称
            sample_col: 样本列名
        
        返回:
            各亚型的标志基因字典
        """
        subtype_col = f'{celltype}_subtype'
        all_markers = {}
        
        samples = adata.obs[sample_col].unique()
        
        for sample in samples:
            # 提取该样本的亚型细胞
            sample_mask = adata.obs[sample_col] == sample
            subtype_mask = adata.obs[subtype_col].notna()
            adata_subtype = adata[sample_mask & subtype_mask].copy()
            
            if len(adata_subtype) == 0:
                continue
            
            # 表达预处理
            sc.pp.normalize_total(adata_subtype, target_sum=1e4)
            sc.pp.log1p(adata_subtype)
            
            # 差异基因分析
            sc.tl.rank_genes_groups(
                adata_subtype,
                groupby=subtype_col,
                method=self.method,
                key_added=f'markers_{sample}'
            )
            
            # 提取标志基因
            subtypes = adata_subtype.obs[subtype_col].unique()
            for subtype in subtypes:
                markers_df = sc.get.rank_genes_groups_df(
                    adata_subtype,
                    group=subtype,
                    key=f'markers_{sample}'
                )
                
                # 筛选显著基因
                sig_markers = markers_df[
                    (markers_df['pvals_adj'] < self.pval_threshold) &
                    (markers_df['logfoldchanges'].abs() > self.logfc_threshold)
                ]
                
                all_markers[subtype] = sig_markers['names'].tolist()[:self.top_n]
        
        return all_markers
```

#### enrichment.py - 功能富集分析
```python
"""
功能富集分析模块
GO/KEGG通路富集
"""

import gseapy as gp
from typing import Dict, List

class EnrichmentAnalyzer:
    """功能富集分析器"""
    
    def __init__(self, gene_sets: List[str] = None,
                 pval_threshold: float = 0.05,
                 top_n: int = 20):
        self.gene_sets = gene_sets or [
            'GO_Biological_Process_2023',
            'KEGG_2021_Human'
        ]
        self.pval_threshold = pval_threshold
        self.top_n = top_n
    
    def analyze(self, markers_dict: Dict[str, List[str]]) -> Dict:
        """
        功能富集分析
        
        参数:
            markers_dict: 各亚型标志基因字典
        
        返回:
            各亚型的富集结果字典
        """
        enrichment_results = {}
        
        for subtype, genes in markers_dict.items():
            if len(genes) < 5:
                continue
            
            try:
                # GO富集
                go_result = gp.enrichr(
                    gene_list=genes,
                    gene_sets='GO_Biological_Process_2023',
                    cutoff=self.pval_threshold
                )
                
                # KEGG富集
                kegg_result = gp.enrichr(
                    gene_list=genes,
                    gene_sets='KEGG_2021_Human',
                    cutoff=self.pval_threshold
                )
                
                enrichment_results[subtype] = {
                    'GO': go_result.results.head(self.top_n),
                    'KEGG': kegg_result.results.head(self.top_n)
                }
            except Exception as e:
                print(f"亚型 {subtype} 富集分析失败: {e}")
                enrichment_results[subtype] = None
        
        return enrichment_results
```

### 4.5 跨样本统一模块设计

#### similarity.py - 多维度相似性计算
```python
"""
多维度相似性计算模块
计算特征基因、功能通路、生态位相似性
"""

import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity

class SimilarityCalculator:
    """相似性计算器"""
    
    def calculate_marker_similarity(self, markers_dict: Dict) -> np.ndarray:
        """
        计算特征基因相似性（Jaccard）
        """
        subtypes = list(markers_dict.keys())
        n = len(subtypes)
        similarity = np.zeros((n, n))
        
        for i, subtype_i in enumerate(subtypes):
            for j, subtype_j in enumerate(subtypes):
                genes_i = set(markers_dict[subtype_i])
                genes_j = set(markers_dict[subtype_j])
                
                if len(genes_i | genes_j) > 0:
                    similarity[i, j] = len(genes_i & genes_j) / len(genes_i | genes_j)
        
        return similarity
    
    def calculate_pathway_similarity(self, enrichment_dict: Dict) -> np.ndarray:
        """
        计算功能通路相似性
        """
        subtypes = list(enrichment_dict.keys())
        n = len(subtypes)
        similarity = np.zeros((n, n))
        
        for i, subtype_i in enumerate(subtypes):
            for j, subtype_j in enumerate(subtypes):
                if enrichment_dict[subtype_i] is None or enrichment_dict[subtype_j] is None:
                    continue
                
                pathways_i = set(enrichment_dict[subtype_i]['KEGG']['Term'].tolist())
                pathways_j = set(enrichment_dict[subtype_j]['KEGG']['Term'].tolist())
                
                if len(pathways_i | pathways_j) > 0:
                    similarity[i, j] = len(pathways_i & pathways_j) / len(pathways_i | pathways_j)
        
        return similarity
    
    def calculate_niche_similarity(self, niche_dict: Dict) -> np.ndarray:
        """
        计算生态位组成相似性（余弦相似度）
        """
        subtypes = list(niche_dict.keys())
        n = len(subtypes)
        similarity = np.zeros((n, n))
        
        compositions = []
        for subtype in subtypes:
            compositions.append(niche_dict[subtype]['composition'])
        
        compositions = np.array(compositions)
        similarity = cosine_similarity(compositions)
        
        return similarity
```

#### mapping.py - 亚型映射表生成
```python
"""
亚型映射表生成模块
生成原始亚型到统一亚型的映射关系
"""

import pandas as pd
from sklearn.cluster import AgglomerativeClustering

class SubtypeMapper:
    """亚型映射器"""
    
    def __init__(self, similarity_threshold: float = 0.6):
        self.similarity_threshold = similarity_threshold
    
    def generate_mapping(self, fused_similarity: np.ndarray,
                        subtypes: list) -> Dict[str, str]:
        """
        生成亚型映射
        
        参数:
            fused_similarity: 融合相似度矩阵
            subtypes: 原始亚型列表
        
        返回:
            映射字典 {原始亚型: 统一亚型}
        """
        # 相似度转距离
        distance_matrix = 1 - fused_similarity
        
        # 层次聚类
        clustering = AgglomerativeClustering(
            n_clusters=None,
            distance_threshold=1 - self.similarity_threshold,
            affinity='precomputed',
            linkage='average'
        )
        cluster_labels = clustering.fit_predict(distance_matrix)
        
        # 生成统一命名
        unified_names = {}
        for cluster_id in set(cluster_labels):
            cluster_subtypes = [subtypes[i] for i, label in enumerate(cluster_labels) 
                               if label == cluster_id]
            
            unified_name = f'Subtype_{chr(65 + cluster_id)}'
            
            for subtype in cluster_subtypes:
                unified_names[subtype] = unified_name
        
        return unified_names
    
    def create_mapping_table(self, unified_names: Dict) -> pd.DataFrame:
        """
        创建映射表DataFrame
        """
        mapping_table = pd.DataFrame({
            'original_subtype': list(unified_names.keys()),
            'unified_subtype': list(unified_names.values()),
            'sample': [s.split('_')[0] for s in unified_names.keys()],
            'cluster_id': [s.split('_')[1] for s in unified_names.keys()]
        })
        
        return mapping_table
```

### 4.6 主流程入口 (pipeline.py)

```python
"""
SpaMFC主流程入口
整合所有模块，提供完整的分析流程
"""

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Optional, Dict

from .config import ConfigManager, SpaMFCConfig
from .features.spatial import SpatialFeatureExtractor
from .features.cnv import CNVFeatureProcessor
from .features.expression import ExpressionFeatureProcessor
from .fusion.weighting import AdaptiveWeightCalculator
from .fusion.nmf import NMFConsensus
from .fusion.clustering import SubtypeClusterer
from .annotation.markers import MarkerGeneAnalyzer
from .annotation.enrichment import EnrichmentAnalyzer
from .unification.similarity import SimilarityCalculator
from .unification.mapping import SubtypeMapper
from .visualization.spatial_plot import SpatialVisualizer
from .visualization.report import ReportGenerator

class SpaMFCPipeline:
    """SpaMFC分析流程"""
    
    def __init__(self, config_path: Optional[str] = None):
        self.config_manager = ConfigManager(config_path)
        self.config = self.config_manager.config
        
        # 初始化各模块
        self._init_modules()
    
    def _init_modules(self):
        """初始化各功能模块"""
        # 特征提取模块
        self.spatial_extractor = SpatialFeatureExtractor(self.config.spatial_config)
        self.cnv_processor = CNVFeatureProcessor(self.config.cnv_config)
        self.expr_processor = ExpressionFeatureProcessor(self.config.expression_config)
        
        # 特征融合模块
        self.weight_calculator = AdaptiveWeightCalculator()
        self.nmf_consensus = NMFConsensus(
            self.config.nmf_runs, 
            self.config.nmf_components
        )
        self.clusterer = SubtypeClusterer(
            self.config.clustering_method,
            self.config.n_clusters,
            per_sample=self.config.per_sample
        )
        
        # 功能注释模块
        self.marker_analyzer = MarkerGeneAnalyzer()
        self.enrichment_analyzer = EnrichmentAnalyzer()
        
        # 跨样本统一模块
        self.similarity_calculator = SimilarityCalculator()
        self.subtype_mapper = SubtypeMapper()
        
        # 可视化模块
        self.visualizer = SpatialVisualizer()
        self.report_generator = ReportGenerator()
    
    def run(self, adata: sc.AnnData, celltype: str) -> sc.AnnData:
        """
        执行完整分析流程
        
        参数:
            adata: AnnData对象
            celltype: 目标细胞类型
        
        返回:
            更新后的adata
        """
        if self.config.verbose:
            print(f"开始分析细胞类型: {celltype}")
            print(f"启用的特征: {self.config_manager.get_active_features()}")
        
        # Step 1: 提取目标细胞
        target_cells = adata.obs[
            adata.obs[self.config.celltype_col] == celltype
        ].index.tolist()
        
        if len(target_cells) == 0:
            raise ValueError(f"未找到细胞类型: {celltype}")
        
        # Step 2: 特征提取（根据用户配置）
        feature_dfs = {}
        
        # 空间特征（用户可选择）
        if self.config.use_spatial:
            if self.config.verbose:
                print("提取空间邻域特征...")
            global_cell_types = adata.obs[self.config.celltype_col].unique()
            feature_dfs['spatial'] = self.spatial_extractor.extract(
                adata, target_cells,
                self.config.celltype_col,
                global_cell_types
            )
        
        # CNV特征（用户可选择，仅恶性细胞）
        if self.config.use_cnv:
            if self.config.verbose:
                print("处理CNV特征...")
            feature_dfs['cnv'] = self.cnv_processor.process(
                adata, target_cells,
                cnv_key='X_cnv'
            )
        
        # 表达特征（用户可选择）
        if self.config.use_expression:
            if self.config.verbose:
                print("处理表达特征...")
            feature_dfs['expression'] = self.expr_processor.process(
                adata, target_cells
            )
        
        # Step 3: 自适应权重计算
        weights = self.weight_calculator.calculate(feature_dfs)
        if self.config.verbose:
            print(f"特征权重: {weights}")
        
        # Step 4: 特征融合
        fused_df = self._fuse_features(feature_dfs, weights)
        
        # Step 5: NMF共识分解
        W_df = self.nmf_consensus.decompose(fused_df)
        
        # Step 6: 聚类
        adata = self.clusterer.cluster(
            adata, W_df,
            self.config.sample_col,
            celltype
        )
        
        # Step 7: 特征基因分析（可选）
        if self.config.enable_marker_analysis:
            markers_dict = self.marker_analyzer.analyze(
                adata, celltype, self.config.sample_col
            )
        
        # Step 8: 功能富集分析（可选）
        if self.config.enable_enrichment:
            enrichment_dict = self.enrichment_analyzer.analyze(markers_dict)
        
        # Step 9: 跨样本统一（可选）
        if self.config.enable_unification:
            # 计算多维度相似性
            marker_sim = self.similarity_calculator.calculate_marker_similarity(markers_dict)
            pathway_sim = self.similarity_calculator.calculate_pathway_similarity(enrichment_dict)
            
            # 融合相似度
            fused_sim = (
                0.4 * marker_sim +
                0.3 * pathway_sim
            )
            
            # 生成映射
            subtypes = list(markers_dict.keys())
            unified_names = self.subtype_mapper.generate_mapping(fused_sim, subtypes)
            
            # 应用映射
            subtype_col = f'{celltype}_subtype'
            unified_col = f'{celltype}_subtype_unified'
            adata.obs[unified_col] = adata.obs[subtype_col].map(unified_names)
        
        # Step 10: 可视化（可选）
        if self.config.save_plots:
            self.visualizer.plot_spatial(
                adata, celltype,
                self.config.output_dir
            )
        
        # Step 11: 生成报告
        self.report_generator.generate(
            adata, celltype,
            markers_dict, enrichment_dict,
            self.config.output_dir
        )
        
        if self.config.verbose:
            print(f"细胞类型 {celltype} 分析完成")
        
        return adata
    
    def _fuse_features(self, feature_dfs: Dict, weights: Dict) -> pd.DataFrame:
        """融合多模态特征"""
        fused_parts = []
        
        for feature_name, df in feature_dfs.items():
            if df is not None:
                # 标准化
                scaled = MinMaxScaler().fit_transform(df.values)
                scaled_df = pd.DataFrame(
                    scaled,
                    index=df.index,
                    columns=df.columns
                )
                # 加权
                weighted = scaled_df * weights.get(feature_name, 0.3)
                fused_parts.append(weighted)
        
        fused_df = pd.concat(fused_parts, axis=1)
        fused_df = fused_df.fillna(0.0)
        
        return fused_df
    
    def run_multi_celltypes(self, adata: sc.AnnData,
                            celltypes: list) -> sc.AnnData:
        """
        执行多细胞类型分析
        
        参数:
            adata: AnnData对象
            celltypes: 细胞类型列表
        
        返回:
            更新后的adata
        """
        for celltype in celltypes:
            adata = self.run(adata, celltype)
        
        return adata
```

---

## 五、使用示例

### 5.1 基本使用

```python
import scanpy as sc
from SpaMFC import SpaMFCPipeline

# 加载数据
adata = sc.read_h5ad("./data/input.h5ad")

# 创建分析流程（使用默认配置）
pipeline = SpaMFCPipeline()

# 分析恶性肿瘤细胞（启用所有特征）
pipeline.config.use_spatial = True
pipeline.config.use_cnv = True
pipeline.config.use_expression = True
adata = pipeline.run(adata, 'Malignant cells')

# 分析CAFs（仅使用空间+表达特征）
pipeline.config.use_spatial = True
pipeline.config.use_cnv = False
pipeline.config.use_expression = True
adata = pipeline.run(adata, 'CAFs')

# 分析免疫细胞（仅使用表达特征）
pipeline.config.use_spatial = False
pipeline.config.use_cnv = False
pipeline.config.use_expression = True
adata = pipeline.run(adata, 'ILC')
```

### 5.2 使用配置文件

```python
from SpaMFC import SpaMFCPipeline

# 使用恶性细胞配置文件
pipeline = SpaMFCPipeline('configs/malignant_config.yaml')
adata = pipeline.run(adata, 'Malignant cells')

# 使用CAFs配置文件
pipeline = SpaMFCPipeline('configs/cafs_config.yaml')
adata = pipeline.run(adata, 'CAFs')
```

### 5.3 多细胞类型批量分析

```python
from SpaMFC import SpaMFCPipeline

pipeline = SpaMFCPipeline('configs/default_config.yaml')

# 批量分析多种细胞类型
celltypes = ['Malignant cells', 'CAFs', 'ILC', 'ECs', 'Epithelial Cell']
adata = pipeline.run_multi_celltypes(adata, celltypes)
```

---

## 六、项目执行计划

### Phase 1: 项目结构搭建（第1周）
- [ ] 创建项目目录结构
- [ ] 编写setup.py和requirements.txt
- [ ] 实现配置管理模块(config.py)
- [ ] 编写默认配置文件

### Phase 2: 特征提取模块（第2周）
- [ ] 实现空间特征提取(spatial.py)
- [ ] 实现CNV特征处理(cnv.py)
- [ ] 实现表达特征处理(expression.py)
- [ ] 实现生态位特征(niche.py)

### Phase 3: 特征融合与聚类模块（第3周）
- [ ] 实现自适应权重计算(weighting.py)
- [ ] 实现NMF共识分解(nmf.py)
- [ ] 实现聚类与稳定性验证(clustering.py)

### Phase 4: 功能注释模块（第4周）
- [ ] 实现特征基因分析(markers.py)
- [ ] 实现功能富集分析(enrichment.py)
- [ ] 实现生态位分析(niche_analysis.py)

### Phase 5: 跨样本统一模块（第5周）
- [ ] 实现相似性计算(similarity.py)
- [ ] 实现相似度融合(fusion.py)
- [ ] 实现亚型映射(mapping.py)

### Phase 6: 可视化与主流程（第6周）
- [ ] 实现空间可视化(spatial_plot.py)
- [ ] 实现报告生成(report.py)
- [ ] 实现主流程入口(pipeline.py)
- [ ] 编写使用示例和测试

---

## 七、质量控制要点

### 7.1 配置验证
- 至少启用一种特征类型
- CNV特征需要CNV数据存在
- 生态位分析需要cell2location数据

### 7.2 特征质量检查
- 空间特征：邻域细胞数 >= 5
- CNV特征：PCA解释方差 >= 10%
- 表达特征：高变基因数 >= 1000

### 7.3 聚类质量检查
- NMF收敛性（迭代次数 < 1000）
- 聚类稳定性（稳定亚型数 >= 3）
- 各亚型细胞数分布均衡

---

## 八、注意事项

1. **特征选择灵活性**：用户可自由组合特征，至少启用一种
2. **CNV特征限制**：仅适用于恶性细胞，其他类型自动禁用
3. **分样本聚类**：建议启用以避免批次效应
4. **跨样本统一**：需要足够的样本数（>= 3）
5. **内存管理**：大样本数据建议分块处理
6. **配置文件优先级**：命令行参数 > 配置文件 > 默认值

---

*计划完善日期：2026-04-20*
*项目路径：D:\01.工作\02.研发\SpaMFC*