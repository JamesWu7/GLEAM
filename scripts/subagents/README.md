# GLEAM 子代理分工（精简常驻模式）

本目录定义自治开发流水线中的子代理角色，目标是减少重复分析角色并保持自动修复闭环。

- 机器可读清单：`subagents_zh.yaml`
- 目标流程：`分析 -> 根因 -> 修复 -> 验证 -> 复核 -> 自动推送/PR`

## 一览表

| 代理ID | 中文名称 | 阶段 | 是否可写代码 | 核心职责 |
|---|---|---|---|---|
| `supervisor_agent` | 总控代理 | 全流程 | 否 | 编排流程、自动推送分支、创建PR |
| `ci_explorer_agent` | CI日志探索代理 | 分析 | 否 | 解析 build/check 失败点 |
| `pkgdown_explorer_agent` | 文档构建探索代理 | 分析 | 否 | 解析 pkgdown/教程页问题 |
| `bug_analyzer_agent` | 根因分析代理 | 根因定位 | 否 | 输出 `error_type/file/root_cause/fix_strategy` |
| `fix_engineer_agent` | 修复工程代理 | 修复 | 是 | 最小补丁修复（唯一可写代码） |
| `test_runner_agent` | 测试验证代理 | 验证 | 否 | 执行 check/tests/pkgdown |
| `final_reviewer_agent` | 最终复核代理 | 复核 | 否 | 回归评估，输出 readiness |
| `pr_writer_agent` | PR文案代理 | PR生成 | 否 | 生成 PR 标题、提交信息、说明 |

## 分工原则

1. 只有 `fix_engineer_agent` 可以改代码。
2. 分析阶段仅保留 CI 与 pkgdown 两条主链，去掉重复角色。
3. 验证结果必须上传为 GitHub artifact。
4. 审核通过后自动推送到 GitHub 分支并创建 PR。
5. 若 `R-CMD-check`/`pkgdown` 在 PR 分支失败，自动抓取日志并回推修复。
