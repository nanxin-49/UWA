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
- MATLAB shared session：通过，版本为 MATLAB R2025b Update 6；仓库根目录为
  `E:\MISC\CARPE3D_matlab\Explain`。`setup_vertical_project` 成功且未改变当前目录；
  `scripts/`、`examples/`、`scripts/validation/` 及各 `src` 模块路径可解析。
- MATLAB Code Analyzer：本次涉及的正式入口、核心 `src` 文件和 PE validation
  入口均无 error；`examples/comm_main_vertical_psk_demo.m` 仅有一个既有的
  unreachable-statement warning，未发现 syntax、undefined local function、路径解析
  或 function/file name mismatch blocker。
- 主入口 smoke：`main_vertical` 成功调用 `examples/main_vertical_demo.m`；默认
  4 kHz uniform、1024x1024 demo 在临时目录完成，无仓库生成物。
- PE smoke：128x128、4 kHz、uniform、`save_mode=rx_only`、无海面反射的轻量 case
  通过，`pass_1_over_R=true`、`Nf=1`；同一轻量参数启用 PM 海面后通过，输出
  `surface_elevation`/`h_reflect` 有限且有效（最大海面幅度约 0.44528 m）。
- Communication smoke：兼容入口在 `COMM_NX=256`、`COMM_NY=256`、`Nf=8`、
  `n_sym=32`、`EbN0=10 dB` 的轻量配置下完成 `direct_only` 与
  `direct_plus_reflect` 两个场景；两者 BER/SER 均为 0，入口逻辑完成。demo 内部
  的 `clear` 会清除外层临时变量，但不影响入口执行结果。
- Validation isolation：production/demo 与正式入口扫描未发现对
  `scripts/validation`、Helmholtz BIE、Model-1 kz-aware 或 Model-2 angle+slope
  的隐式调用。
- Bellhop：`BELLHOP_EXE`/`BELLHOP_TOOLBOX_ROOT` 未配置；解析结果为空，轻量缺失
  依赖 guard 清晰返回 `Set overrides.bellhop_exe or BELLHOP_EXE.`，未启动外部
  Bellhop，也未伪造运行成功。

## 结论

MAIN_BASELINE_READY

本次 MATLAB baseline regression 未发现阻塞 PR #1 合入 `main` 的问题；未修改物理模型、验证算法或仓库结构。
