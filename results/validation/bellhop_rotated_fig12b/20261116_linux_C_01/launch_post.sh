#!/bin/bash
set -u
ROOT=/home/misc/cwx/bellhop_rotated_C_transfer_20261116_01/Explain
OUT=$ROOT/results/validation/bellhop_rotated_fig12b/20261116_linux_C_post01
RUN=$ROOT/results/validation/bellhop_rotated_fig12b/20261116_linux_C_01
export BL_JOB_START_EPOCH=$(date +%s)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
/home/misc/MATLAB/R2025b/bin/matlab -singleCompThread -batch "addpath('$OUT/code'); summarize_bellhop_rotated_C_saved('$ROOT','$RUN','$OUT/audit');" > "$OUT/stdout.log" 2> "$OUT/stderr.log"
code=$?
printf '%s\n' "$code" > "$OUT/exit_code.txt"
exit "$code"
