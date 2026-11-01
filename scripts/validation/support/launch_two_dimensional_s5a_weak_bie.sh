#!/usr/bin/env bash
set -euo pipefail
TASK=/home/misc/cwx/uwa2d_s5a_20261027
MATLAB=/home/misc/MATLAB2023b/bin/matlab
JOB=${1:?job id required}
ids=(smoke model main96 self64 space14 window60 window70)
upper=(80 120 600 450 800 450 600)
nodes=(0 0 5760 5760 6720 4800 5600)
ix=-1
for i in "${!ids[@]}"; do
  if [[ ${ids[$i]} == "$JOB" ]]; then ix=$i; break; fi
done
(( ix >= 0 )) || { echo 'Unknown job'; exit 2; }
[[ $(realpath "$TASK") == "$TASK" && -x "$MATLAB" ]] || exit 3
sha256sum -c "$TASK/input/SHA256SUMS" || exit 4
spent=0
for ((i=0;i<ix;i++)); do
  previous="$TASK/jobs/${ids[$i]}"
  [[ -f "$previous/exit_code.txt" && -f "$previous/result.mat" ]] || exit 5
  [[ $(cat "$previous/exit_code.txt") == 0 ]] || exit 6
  elapsed=$(cat "$previous/process_elapsed_seconds.txt")
  spent=$((spent+elapsed))
done
projected=$((spent+300))
for ((i=ix;i<${#ids[@]};i++)); do projected=$((projected+upper[$i])); done
(( projected <= 3600 )) || { echo "Cumulative 3600 s gate: $projected"; exit 7; }
available=$(awk '/MemAvailable:/ {printf "%.0f", $2*1024}' /proc/meminfo)
if (( ix == 0 )); then
  increment=$((128*1024*1024))
elif (( ix == 1 )); then
  increment=$((256*1024*1024))
else
  n=${nodes[$ix]}
  increment=$((64*n*n+128*1024*1024))
fi
(( increment*10 <= available*7 )) || { echo "Physical RAM gate: $increment $available"; exit 8; }
jobdir="$TASK/jobs/$JOB"
mkdir -- "$jobdir" # atomic duplicate prevention: no -p
printf 'job=%s\navailable_physical_bytes=%s\nestimated_increment_bytes=%s\nprojected_cumulative_seconds=%s\n' \
  "$JOB" "$available" "$increment" "$projected" > "$jobdir/shell_preflight.txt"
date -u +'%Y-%m-%dT%H:%M:%SZ' > "$jobdir/start_utc.txt"
start=$(date +%s)
set +e
timeout --signal=TERM --kill-after=15s 900s "$MATLAB" -batch \
  "addpath('$TASK/code'); run_two_dimensional_s5a_weak_bie_remote('$TASK','$JOB',$available)" \
  > "$jobdir/stdout.log" 2> "$jobdir/stderr.log"
status=$?
set -e
end=$(date +%s)
printf '%s\n' "$status" > "$jobdir/exit_code.txt"
printf '%s\n' "$((end-start))" > "$jobdir/process_elapsed_seconds.txt"
date -u +'%Y-%m-%dT%H:%M:%SZ' > "$jobdir/end_utc.txt"
if (( status != 0 )); then
  echo "JOB_FAILED $JOB exit=$status elapsed=$((end-start))" >&2
  exit "$status"
fi
[[ -f "$jobdir/result.mat" && -f "$jobdir/receipt.json" ]] || exit 9
sha256sum "$jobdir/result.mat" "$jobdir/receipt.json" > "$jobdir/output_sha256.txt"
echo "JOB_COMPLETE $JOB elapsed=$((end-start))"
