#!/bin/bash
set -euo pipefail
OUT=/home/misc/cwx/bellhop_1500m_DS4062p5_C801_20261119_01
ROOT=/home/misc/cwx/bellhop_rotated_C_transfer_20261116_01/Explain
printf '%s\n' "$$" > "$OUT/MATLAB_pid.txt"
cat "/proc/$$/stat" > "$OUT/MATLAB_identity.txt"
exec /home/misc/MATLAB/R2025b/bin/matlab -singleCompThread -batch "addpath('$OUT/code','-begin'); execute_p1_ds4062('$ROOT','$OUT');"