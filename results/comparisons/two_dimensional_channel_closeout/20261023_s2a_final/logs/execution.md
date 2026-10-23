# S2a 执行记录

- 会话：用户恢复的既有 MATLAB R2025b Update 6，经 MATLAB MCP 运行；运行前后 cwd 均为 `E:\MISC\CARPE3D_matlab\Explain`。初次连接报 `failed to attach to MATLAB session`，恢复后没有从 shell 启动 MATLAB。
- 运行：`setup_vertical_project; r=run_two_dimensional_closeout_s2a('20261023_s2a_final',1718734848); a=audit_two_dimensional_closeout_s2a(r.run_dir);`。`1718734848` 字节取自运行前 Windows `\Memory\Available Bytes` 计数器；保守增量估计见 `../data/cost.json`，满足 30% 可用物理 RAM 保留门。运行名是唯一目录标签，不宣称精确测量时刻。
- 物理工作：正入射常数源 3 个实际 PE 入口控制；保存 Gaussian 源的平面和三档 PE 实际入口核对；完整 FullFK 平面与弱档线复算及中档新线；强档复用已核验 E12 保存线。没有 BIE 边界求解或保存密度接收场评价，没有 Bellhop 作业、PM、整套 x-z 场、宽带、FastFK 或生产 API 运行。
- 诊断运行 `20261023_s2a_a` 在 `memory` 返回结构缺 `PhysicalMemory` 字段时停止，未生成正式输出。改为 Windows 物理可用计数器显式传入。`20261023_s2a_b` 在此基础上通过审计；最终运行增加 E01–E09 继承哈希与源/版本元数据并独立复核，旧诊断目录均保留。
- 六个改动/新增 `.m` 的 MATLAB Code Analyzer 均为零消息。旧验证 PE helper 的默认 Gaussian 与显式传入同数组的 128 点/32 m/Tx10/Rx3/步长0.5 m 回归完全相同，入射和反射最大差都为零。`git diff --check` 在最终文档更新后再次执行。
- 最终文档更新后再次运行六文件 Code Analyzer，仍全部零消息。对已存在 `audit/` 的正式目录直接重复调用审计入口，被入口第8行的非覆盖断言预期拒绝；首次正式审计已通过并保存在 `../audit/independent_audit.json`，未删除或覆盖审计证据。`git diff --check` 通过（仅Git提示两份已有CRLF文档将按配置处理行尾）。
