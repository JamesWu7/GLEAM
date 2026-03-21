# GLEAM 子代理分工（中文）

本目录用于定义自治开发流水线中的子代理角色，避免职责重叠。

- 机器可读清单：`subagents_zh.yaml`
- 目标流程：`分析 -> 修复 -> 验证 -> PR -> GitHub人审`

## 一览表

| 代理ID | 中文名称 | 阶段 | 是否可写代码 | 核心职责 |
|---|---|---|---|---|
| `supervisor_agent` | 总控代理 | 全流程 | 否 | 编排流程、创建分支/PR（不直接合并） |
| `explorer_agent` | 代码探索代理 | 分析 | 否 | 扫描代码与潜在风险 |
| `reviewer_agent` | 规范评审代理 | 分析 | 否 | 检查 API/文档一致性 |
| `ci_explorer_agent` | CI日志探索代理 | 分析 | 否 | 解析 build/check 失败点 |
| `pkgdown_explorer_agent` | 文档构建探索代理 | 分析 | 否 | 解析 pkgdown/教程页问题 |
| `bug_analyzer_agent` | 根因分析代理 | 根因定位 | 否 | 输出 `error_type/file/root_cause/fix_strategy` |
| `fix_engineer_agent` | 修复工程代理 | 修复 | 是 | 最小补丁修复（唯一可写代码） |
| `test_runner_agent` | 测试验证代理 | 验证 | 否 | 执行 check/tests/pkgdown |
| `pr_writer_agent` | PR文案代理 | PR生成 | 否 | 生成 PR 标题、提交信息、说明 |
| `final_reviewer_agent` | 最终复核代理 | 复核 | 否 | 回归评估，输出是否可人审 |

## 分工原则

1. 只有 `fix_engineer_agent` 可以改代码。
2. 分析阶段 4 个代理并行，先找全量问题。
3. 根因先统一收敛，再进入修复。
4. 验证结果必须上传为 GitHub artifact。
5. 合并前必须 GitHub 人工审核。
