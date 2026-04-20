# CLI参数优化与特征融合改进 Spec

## Why
当前CLI参数设计不够简洁，特征选择参数过于分散，需要优化为更直观的单参数模式。同时特征融合方法需要支持多种选择，让用户可以灵活配置。

## What Changes
- 简化特征选择参数：从多个`--use-xxx/--no-xxx`改为单个`--features`参数
- 添加细胞类型列参数：用户必须指定细胞类型列名才能选择目标细胞
- 添加特征融合方法选择参数
- 优化特征提取和融合的相关代码

## Impact
- Affected specs: CLI接口设计、特征融合模块
- Affected code: src/cli.py, src/fusion/weighting.py, src/pipeline.py

## ADDED Requirements
### Requirement: 简化的特征选择参数
系统 SHALL 提供单个`--features`参数用于特征选择，支持以下格式：
- 单特征：`--features spatial`
- 双特征组合：`--features spatial,expression`
- 三特征组合：`--features spatial,cnv,expression`
- 默认值：`spatial`

#### Scenario: 特征选择
- **WHEN** 用户执行 `spamfc_cli run --input data.h5ad --celltype "CAFs" --features spatial,expression`
- **THEN** 系统使用空间特征和表达特征进行分析

### Requirement: 细胞类型列参数
系统 SHALL 要求用户指定细胞类型列名参数`--celltype-col`，用于选择目标细胞。

#### Scenario: 细胞类型列指定
- **WHEN** 用户执行 `spamfc_cli run --input data.h5ad --celltype-col "anno_cell2location_res" --celltype "CAFs"`
- **THEN** 系统从指定列中筛选CAFs细胞进行分析

### Requirement: 特征融合方法选择
系统 SHALL 提供`--fusion-method`参数，支持以下融合方法：
- `adaptive`: 自适应权重融合（默认）
- `fixed`: 固定权重融合
- `concat`: 直接拼接（无权重）

#### Scenario: 融合方法选择
- **WHEN** 用户执行 `spamfc_cli run --input data.h5ad --celltype "CAFs" --fusion-method adaptive`
- **THEN** 系统使用自适应权重方法融合特征

## MODIFIED Requirements
### Requirement: 特征选择参数
原有多个特征开关参数（`--use-spatial`, `--no-spatial`等）改为单个`--features`参数。

### Requirement: 特征融合代码
特征融合模块需要支持多种融合方法，包括自适应权重、固定权重和直接拼接。

## REMOVED Requirements
### Requirement: 多个特征开关参数
**Reason**: 参数过于分散，用户体验不佳
**Migration**: 使用单个`--features`参数替代