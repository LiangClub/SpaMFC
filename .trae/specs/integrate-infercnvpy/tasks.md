# inferCNVpy 功能集成任务列表

## 任务概览

| 任务ID | 任务名称 | 优先级 | 状态 |
|--------|----------|--------|------|
| T-001 | 创建 cnv_inference 模块目录 | 高 | pending |
| T-002 | 实现 infercnv.py | 高 | pending |
| T-003 | 实现 genomic_positions.py | 高 | pending |
| T-004 | 实现 cnv_scoring.py | 中 | pending |
| T-005 | 实现 cnv_plot.py | 高 | pending |
| T-006 | 更新 visualization/__init__.py | 中 | pending |
| T-007 | 更新 config.py | 中 | pending |
| T-008 | 更新 cli.py | 高 | pending |
| T-009 | 更新 src/__init__.py | 中 | pending |
| T-010 | 更新 README.md | 低 | pending |
| T-011 | 语法验证 | 高 | pending |

## 任务详情

### T-001: 创建 cnv_inference 模块目录

**描述**: 创建新的 CNV 推断模块目录结构

**输出**: `src/cnv_inference/__init__.py`

---

### T-002: 实现 infercnv.py

**描述**: 实现 CNV 推断核心功能类 CNVInferencer

**功能**:
- `add_genomic_positions()` - 添加基因位置信息
- `infercnv()` - 执行CNV推断
- `run_pipeline()` - 完整CNV分析流程

**输出**: `src/cnv_inference/infercnv.py`

---

### T-003: 实现 genomic_positions.py

**描述**: 实现基因位置处理功能

**功能**:
- `load_gtf()` - 加载GTF文件
- `add_chromosome_info()` - 添加染色体信息到adata.var
- `validate_genomic_positions()` - 验证基因位置完整性

**输出**: `src/cnv_inference/genomic_positions.py`

---

### T-004: 实现 cnv_scoring.py

**描述**: 实现 CNV 评分计算功能

**功能**:
- `compute_cnv_score()` - 计算CNV评分
- `compute_cnv_correlation()` - CNV相关性分析
- `get_cnv_summary()` - CNV统计摘要

**输出**: `src/cnv_inference/cnv_scoring.py`

---

### T-005: 实现 cnv_plot.py

**描述**: 实现 CNV 可视化模块

**功能**:
- `plot_chromosome_heatmap()` - 染色体热图
- `plot_chromosome_heatmap_summary()` - 染色体热图摘要
- `plot_cnv_umap()` - CNV UMAP可视化
- `plot_cnv_scores()` - CNV评分分布图
- `plot_cnv_spatial()` - CNV空间分布图

**输出**: `src/visualization/cnv_plot.py`

---

### T-006: 更新 visualization/__init__.py

**描述**: 导出 CNVVisualizer 类

**修改**: 添加 `from .cnv_plot import CNVVisualizer`

---

### T-007: 更新 config.py

**描述**: 添加 CNVInferenceConfig 配置类

**新增内容**:
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

---

### T-008: 更新 cli.py

**描述**: 添加 cnv 子命令

**新增内容**:
- `_add_cnv_subparser()` 函数
- `cnv_command()` 函数
- 更新 `main()` 函数处理 cnv 命令

---

### T-009: 更新 src/__init__.py

**描述**: 导出 CNV 推断模块

**新增导出**:
- `CNVInferencer`
- `CNVVisualizer`
- `CNVInferenceConfig`

---

### T-010: 更新 README.md

**描述**: 添加 CNV 分析文档

**新增内容**:
- CNV分析使用说明
- CLI命令示例
- Python API示例

---

### T-011: 语法验证

**描述**: 验证所有新增文件语法正确

**验证命令**: `python -m py_compile src/cnv_inference/*.py src/visualization/cnv_plot.py`

## 执行顺序

1. T-001 → T-002 → T-003 → T-004 (CNV推断模块)
2. T-005 → T-006 (CNV可视化模块)
3. T-007 → T-008 → T-009 (配置和CLI)
4. T-010 (文档)
5. T-011 (验证)