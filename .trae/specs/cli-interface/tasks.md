# Tasks

- [x] Task 1: 创建命令行接口模块
  - [x] SubTask 1.1: 创建src/cli.py模块
  - [x] SubTask 1.2: 实现参数解析器（argparse）
  - [x] SubTask 1.3: 实现参数验证逻辑
  - [x] SubTask 1.4: 实现配置文件与命令行参数合并

- [x] Task 2: 实现子命令
  - [x] SubTask 2.1: 实现`run`子命令（单细胞类型分析）
  - [x] SubTask 2.2: 实现`run-multi`子命令（多细胞类型分析）
  - [x] SubTask 2.3: 实现`info`子命令（显示数据信息）
  - [x] SubTask 2.4: 实现`config`子命令（生成默认配置）
  - [x] SubTask 2.5: 实现`version`子命令

- [x] Task 3: 创建命令行入口脚本
  - [x] SubTask 3.1: 创建spamfc_cli.py入口脚本
  - [x] SubTask 3.2: 配置setup.py入口点

- [x] Task 4: 完善帮助文档
  - [x] SubTask 4.1: 添加完整参数帮助信息
  - [x] SubTask 4.2: 添加使用示例说明

- [x] Task 5: 测试命令行功能
  - [x] SubTask 5.1: 测试基本分析命令
  - [x] SubTask 5.2: 测试多细胞类型分析
  - [x] SubTask 5.3: 测试配置文件模式

# Task Dependencies
- Task 2 depends on Task 1
- Task 3 depends on Task 2
- Task 4 depends on Task 3
- Task 5 depends on Task 4