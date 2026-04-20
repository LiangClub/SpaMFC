# 项目代码优化与多余文件清理 Spec

## Why
项目中有多个参考用的notebook文件和旧的计划文件，这些文件只是开发参考，不属于最终项目代码，需要清理以保持项目整洁。

## What Changes
- 删除参考用的notebook文件（CAFs.NMF.ipynb等）
- 删除旧的plan.md文件
- 保留核心项目代码结构（src/, configs/, examples/）

## Impact
- Affected specs: 项目文件结构
- Affected code: 无代码变更，仅文件删除

## ADDED Requirements
### Requirement: 项目文件清理
系统 SHALL 保持整洁的项目结构，仅包含必要的项目代码文件。

#### Scenario: 清理参考文件
- **WHEN** 项目开发完成
- **THEN** 删除所有仅用于参考的notebook文件和旧计划文件

## REMOVED Requirements
### Requirement: 参考notebook文件
**Reason**: 这些文件仅用于开发参考，不属于最终项目代码
**Migration**: 已将核心逻辑整合到src/目录下的Python模块中

待删除文件列表：
1. CAFs.NMF.ipynb - CAFs分析参考
2. Cafs.nich.ipynb - CAFs生态位分析参考
3. Malignan.nich.ipynb - 恶性细胞生态位分析参考
4. Tumor.NMF.CNV_1.ipynb - 肿瘤CNV分析参考1
5. Tumor.NMF.CNV_2.ipynb - 肿瘤CNV分析参考2
6. tumor.NMF.CNV_3.ipynb - 肿瘤CNV分析参考3
7. tumor.NMF.ipynb - 肿瘤NMF分析参考
8. tumor.cnv.NMF.ipynb - 肿瘤CNV NMF分析参考
9. plan.md - 旧计划文件（已有新版本在.trae/documents/）