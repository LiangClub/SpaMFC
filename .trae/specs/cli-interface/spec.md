# SpaMFC 命令行接口封装 Spec

## Why
用户需要将SpaMFC项目封装成一个可通过命令行直接调用的软件，方便在不同环境下使用，无需编写Python代码即可完成分析。

## What Changes
- 创建命令行入口脚本 `spamfc_cli.py`
- 支持完整的参数配置（数据输入、特征选择、分析参数等）
- 支持多种分析模式（单细胞类型、多细胞类型、批量分析）
- 支持配置文件和命令行参数两种方式
- 提供帮助信息和参数验证

## Impact
- Affected specs: 项目使用方式
- Affected code: 新增 `src/cli.py` 和 `spamfc_cli.py`

## ADDED Requirements
### Requirement: 命令行接口
系统 SHALL 提供完整的命令行接口，支持以下功能：

#### Scenario: 基本分析
- **WHEN** 用户执行 `spamfc_cli run --input data.h5ad --celltype "CAFs"`
- **THEN** 系统执行CAFs亚群分析并输出结果

#### Scenario: 多细胞类型分析
- **WHEN** 用户执行 `spamfc_cli run-multi --input data.h5ad --celltypes "Malignant cells,CAFs,ILC"`
- **THEN** 系统批量执行多种细胞类型的亚群分析

#### Scenario: 特征选择
- **WHEN** 用户执行 `spamfc_cli run --input data.h5ad --celltype "CAFs" --use-spatial --use-expression --no-cnv`
- **THEN** 系统使用指定的特征组合进行分析

#### Scenario: 配置文件模式
- **WHEN** 用户执行 `spamfc_cli run --config configs/cafs_config.yaml --input data.h5ad`
- **THEN** 系统使用配置文件中的参数进行分析

### Requirement: 参数完整性
命令行接口 SHALL 支持以下参数类别：

1. **数据参数**
   - `--input`: 输入h5ad文件路径
   - `--output`: 输出目录
   - `--sample-col`: 样本列名
   - `--celltype-col`: 细胞类型列名

2. **特征选择参数**
   - `--use-spatial`: 启用空间特征
   - `--use-cnv`: 启用CNV特征
   - `--use-expression`: 启用表达特征
   - `--use-niche`: 启用生态位特征
   - `--no-spatial`: 禁用空间特征
   - `--no-cnv`: 禁用CNV特征
   - `--no-expression`: 禁用表达特征
   - `--no-niche`: 禁用生态位特征

3. **分析参数**
   - `--n-clusters`: 聚类数量
   - `--nmf-components`: NMF因子数
   - `--nmf-runs`: NMF运行次数
   - `--clustering-method`: 聚类方法（kmeans/leiden）
   - `--per-sample`: 分样本聚类

4. **注释参数**
   - `--enable-marker-analysis`: 启用特征基因分析
   - `--enable-enrichment`: 启用功能富集
   - `--enable-unification`: 启用跨样本统一

5. **可视化参数**
   - `--save-plots`: 保存图像
   - `--plot-format`: 图像格式（pdf/png/svg）
   - `--dpi`: 图像分辨率

6. **运行参数**
   - `--verbose`: 详细输出
   - `--n-jobs`: 并行任务数
   - `--config`: 配置文件路径

### Requirement: 子命令
命令行接口 SHALL 提供以下子命令：

- `run`: 单细胞类型分析
- `run-multi`: 多细胞类型批量分析
- `info`: 显示数据信息
- `config`: 生成默认配置文件
- `version`: 显示版本信息
- `help`: 显示帮助信息