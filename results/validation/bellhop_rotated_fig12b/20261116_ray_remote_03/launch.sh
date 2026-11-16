#!/bin/bash
set -u
ROOT=/home/misc/cwx/bellhop_rotated_C_transfer_20261116_01/Explain
OUT=$ROOT/results/validation/bellhop_rotated_fig12b/20261116_ray_remote_03
: "${BH_RAY_REMAINING_S:?requires approved remaining budget}"
: "${BH_RAY_TOTAL_LIMIT_S:?requires approved total limit}"
export BH_RAY_REMAINING_S BH_RAY_TOTAL_LIMIT_S BH_RAY_EXE BH_RAY_EXE_SHA
export BH_RAY_START_EPOCH=$(date +%s)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
printf '%s\n' "$BH_RAY_START_EPOCH" > "$OUT/start_epoch.txt"
/home/misc/MATLAB/R2025b/bin/matlab -singleCompThread -batch "addpath('$OUT/code'); run_bellhop_rays_remote('$ROOT','$OUT');" > "$OUT/stdout.log" 2> "$OUT/stderr.log"
code=$?
printf '%s\n' "$code" > "$OUT/exit_code.txt"
exit "$code"
