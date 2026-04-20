# Tasks

- [ ] Task 1: 优化CLI参数设计
  - [ ] SubTask 1.1: 将多个特征开关参数改为单个`--features`参数
  - [ ] SubTask 1.2: 添加`--celltype-col`必选参数
  - [ ] SubTask 1.3: 添加`--fusion-method`参数
  - [ ] SubTask 1.4: 更新参数解析逻辑

- [ ] Task 2: 优化特征融合模块
  - [ ] SubTask 2.1: 添加直接拼接融合方法
  - [ ] SubTask 2.2: 优化自适应权重计算逻辑
  - [ ] SubTask 2.3: 支持固定权重配置

- [ ] Task 3: 更新pipeline代码
  - [ ] SubTask 3.1: 更新特征选择逻辑
  - [ ] SubTask 3.2: 更新特征融合调用逻辑
  - [ ] SubTask 3.3: 添加融合方法参数传递

- [ ] Task 4: 更新帮助文档和示例
  - [ ] SubTask 4.1: 更新CLI帮助信息
  - [ ] SubTask 4.2: 更新使用示例

# Task Dependencies
- Task 2 depends on Task 1
- Task 3 depends on Task 2
- Task 4 depends on Task 3