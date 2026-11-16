# Bellhop远程拷贝完成（2026-11-16）

本地 `C:\Users\ASUS\.codex\worktrees\bellhop_rotated_fig12b\Explain` 统一映射到 `/home/misc/cwx/bellhop_rotated_C_transfer_20261116_01/Explain`，子目录结构完全保留。逐文件对应表在本地20261116_transfer_01/file_mapping.csv，以及远端同根transfer_metadata/file_mapping.csv。

共2736文件、9532044433字节（8.8774GiB），直接流式传输1156.228秒，无本地大压缩包。tar/SSH解包均退出0；远端50秒逐项SHA校验2736/2736通过，sha256sum退出0，完整日志已取回。准备阶段另核原A/B/C2148项SHA及8个overlay/实际构建源码一致。gfortran已复查为13.3.0，make4.3、ar2.42可用。

范围包括实际修复Fortran构建源码、相关MATLAB脚本、冻结计划、A31引用的full03分量和A/B/C现有证据、原25项保持证明、101束3列Windows基准、源码快照和SHA；按允许清单拷贝，无cash、Git或凭据。少量Windows二进制因历史SHA绑定作为证据保留，不能在Linux运行。

校验完成后的壳脚本遇到CR空行，外层退出1，附加environment_after_transfer.log未生成；不影响已保存的2736项校验PASS。取回三个实际存在的校验记录，第四可选日志下载失败已记录，不重传或重复全量核验。不得据此宣称Linux构建或运行时准备通过。

本轮只拷贝，C仍30/39，未重算、未编译/最小物理验证/求解/ABC图；Windows已有结果和CP均未覆盖，不D/FFP、不提交/推送。原声学预算余额799.165599秒继续保留，文件传输及完整性校验I/O独立记时，不增加声学预算。Windows绝对路径记录未改写，后续Linux调度和路径适配仍需完成。

[完成CP](../results/validation/bellhop_rotated_fig12b/20261116_transfer_01/checkpoint.json)；[完整对应表](../results/validation/bellhop_rotated_fig12b/20261116_transfer_01/file_mapping.csv)；[远端SHA回执](../results/validation/bellhop_rotated_fig12b/20261116_transfer_01/remote_verification/checkpoint_verify.json)。生产接口/物理公式不变，无需修改vertical_comm_guide；进度看板只作导航。
