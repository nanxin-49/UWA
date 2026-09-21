# PR #1 main baseline cleanup report

日期：2026-09-21

## 清理内容

- 删除 `results/validation/**/build/`、`bin/`、`toolchain_tmp/` 中误提交或可重新生成的 Bellhop/Fortran 编译器、Acoustics Toolbox 副本和工具链残留。
- 删除 `scripts/**/__pycache__/` 以及根目录空的临时 `test.m`。
- `.gitignore` 现在在全局 `.m/.md` allow-rule 之后再次排除 `build/`、`toolchain_tmp/`、`bin/`、`__pycache__/` 和 `.cache/`，因此这些目录中的源码/文档也不会被重新纳入 Git。
- 按用户已执行的移动保留 `old/`→`docs/old/` 的历史文件变化；忽略的旧 MAT/PNG/diary 文件未纳入 Git。

## 保留资产与边界

- 保留 `build_manifest.json`、轻量 CSV/JSON、authoritative reports，以及
  `scripts/validation/support` 中人为维护的 Bellhop overlay 源码。
- `src/`、根目录兼容入口和 `examples/` 未发现 Model-1 kz-aware、Model-2
  angle+slope 或 Helmholtz BIE 接入生产 PE 的情况；公共入口未隐式调用
  `scripts/validation/`。
- 生产 PE、PM 海面模型、Bellhop 验证逻辑和通信算法未做物理逻辑改动。
  Bellhop 活动脚本改为通过 `BELLHOP_EXE` 或
  `BELLHOP_TOOLBOX_ROOT` 解析外部依赖。

## 回归检查

- `git diff --check`：通过。
- 生成物目录扫描：通过；活动树无 `build/`、`toolchain_tmp/`、`bin/` 或
  `__pycache__/` 残留。
- `.gitignore` 规则探针：通过，build 内 `.m/.md` 也被忽略。
- 生产代码、入口和活动脚本的机器绝对路径扫描：通过。历史报告和保留的
  结果 CSV 中仍可能有运行时 provenance 路径，它们不被运行时读取，未改动其
  验证结论。
- `build_manifest.json` JSON 结构：已检查；路径字段已改为可移植占位符/相对路径，哈希和验证身份字段保留。
- MATLAB `setup_vertical_project`、静态检查、`main_vertical`/demo、PE/PM、通信 smoke：未执行成功。共享 MATLAB MCP 会话两次均返回 `failed to attach to MATLAB session`；按仓库规则未通过 PowerShell 启动 MATLAB 作为替代。
- Bellhop：未伪造运行成功；仅完成配置解析/依赖入口审查。当前环境未要求运行外部 Bellhop。

## 结论

由于 MATLAB 共享会话不可用，当前不能声明 `MAIN_BASELINE_READY`。除该外部验证环境阻塞外，未发现需要扩大范围或改变既有验证结论的问题。

阻塞项：

1. MATLAB MCP 共享会话不可 attach，因此无法完成本次要求的 MATLAB 静态检查和基础 smoke。
