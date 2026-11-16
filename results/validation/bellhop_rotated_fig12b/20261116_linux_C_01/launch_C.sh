#!/bin/bash
set -u
ROOT=/home/misc/cwx/bellhop_rotated_C_transfer_20261116_01/Explain
OUT=$ROOT/results/validation/bellhop_rotated_fig12b/20261116_linux_C_01
export BL_JOB_START_EPOCH=$(date +%s)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
printf '%s\n' "$BL_JOB_START_EPOCH" > "$OUT/C_start_epoch.txt"
/home/misc/MATLAB/R2025b/bin/matlab -singleCompThread -batch "addpath('$OUT/runtime03'); linux_run_C_remaining('$ROOT','$OUT');" > "$OUT/C_stdout.log" 2> "$OUT/C_stderr.log"
code=$?
printf '%s\n' "$code" > "$OUT/C_exit_code.txt"
exit "$code"
