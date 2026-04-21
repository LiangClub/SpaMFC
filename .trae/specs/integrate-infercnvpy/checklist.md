# inferCNVpy 功能集成检查清单

## 预检查

- [x] 了解 inferCNVpy 功能
- [x] 分析当前 CNV 模块
- [x] 设计集成方案
- [x] 创建 spec.md
- [x] 创建 tasks.md
- [x] 创建 checklist.md

## CNV推断模块

- [ ] 创建 src/cnv_inference/__init__.py
- [ ] 实现 infercnv.py (CNVInferencer类)
- [ ] 实现 genomic_positions.py
- [ ] 实现 cnv_scoring.py

## CNV可视化模块

- [ ] 实现 cnv_plot.py (CNVVisualizer类)
- [ ] 更新 visualization/__init__.py

## 配置和CLI

- [ ] 添加 CNVInferenceConfig 到 config.py
- [ ] 添加 cnv 子命令到 cli.py
- [ ] 更新 src/__init__.py 导出

## 文档更新

- [ ] 更新 README.md CNV分析说明

## 最终验证

- [ ] py_compile 验证所有文件
- [ ] CLI cnv 子命令测试
- [ ] 确认所有任务完成