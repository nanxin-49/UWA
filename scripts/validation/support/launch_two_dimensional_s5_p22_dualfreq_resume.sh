#!/usr/bin/env bash
set -euo pipefail

# Separate continuation after the retained jobs/freeze exit=124. Never edits it.
# Invoke only after the user confirms the remote noninteractive license is fixed.
TASK=/home/misc/cwx/uwa2d_s5_p22_dualfreq_20261029_01
MATLAB_BIN=/home/misc/MATLAB2023b/bin/matlab
JOB=${1:?One continuation job ID is required}
[[ ${LICENSE_RESTORED_ACK:-} == YES ]] || { echo 'License restoration not acknowledged' >&2; exit 2; }
ids=(freeze_retry1 f4031p25_model f4031p25_main96 f4031p25_self64 \
  f4031p25_space14 f4031p25_window60 f4031p25_window70)
upper=(180 180 300 300 350 300 350)
ix=-1
for i in "${!ids[@]}"; do
  if [[ ${ids[$i]} == "$JOB" ]]; then ix=$i; break; fi
done
(( ix >= 0 )) || { echo 'Unknown continuation job' >&2; exit 3; }
[[ -d "$TASK" && $(realpath "$TASK") == "$TASK" && -x "$MATLAB_BIN" ]] || exit 4
cd "$TASK"
sha256sum -c input/SHA256SUMS >/dev/null || exit 5

prior_ids=(smoke f4062p5_model f4062p5_main96 f4062p5_self64 \
  f4062p5_space14 f4062p5_window60 f4062p5_window70 freeze)
spent=0
for previous in "${prior_ids[@]}"; do
  if [[ $previous == smoke || $previous == freeze ]]; then
    prior_dir=$TASK/jobs/$previous
  else
    frequency=${previous%%_*}
    layer=${previous#*_}
    prior_dir=$TASK/jobs/$frequency/$layer
  fi
  [[ -f "$prior_dir/exit_code.txt" && -f "$prior_dir/process_elapsed_seconds.txt" ]] || exit 6
  status=$(cat "$prior_dir/exit_code.txt")
  if [[ $previous == freeze ]]; then
    [[ $status == 124 && ! -e "$TASK/prediction" ]] || {
      # A successful retry creates prediction; retain failed freeze status.
      [[ $status == 124 && -f "$TASK/prediction/prediction.json" ]] || exit 7
    }
  else
    [[ $status == 0 && -f "$prior_dir/result.mat" ]] || exit 8
  fi
  elapsed=$(cat "$prior_dir/process_elapsed_seconds.txt")
  [[ $elapsed =~ ^[0-9]+$ ]] || exit 9
  spent=$((spent+elapsed))
done
for ((i=0;i<ix;i++)); do
  previous=${ids[$i]}
  if [[ $previous == freeze_retry1 ]]; then
    prior_dir=$TASK/jobs/freeze_retry1
    prior_result=$TASK/prediction/prediction.json
  else
    layer=${previous#*_}
    prior_dir=$TASK/jobs/f4031p25/$layer
    prior_result=$prior_dir/result.mat
  fi
  [[ -f "$prior_dir/exit_code.txt" && -f "$prior_dir/process_elapsed_seconds.txt" \
    && -f "$prior_result" ]] || exit 10
  [[ $(cat "$prior_dir/exit_code.txt") == 0 ]] || exit 11
  elapsed=$(cat "$prior_dir/process_elapsed_seconds.txt")
  [[ $elapsed =~ ^[0-9]+$ ]] || exit 12
  spent=$((spent+elapsed))
done
(( spent + upper[ix] <= 1800 )) || {
  echo "1800 s process budget gate: spent=$spent next=${upper[$ix]}" >&2; exit 13;
}
disk_free=$(df -B1 --output=avail "$TASK" | tail -1 | tr -d '[:space:]')
tmp_free=$(df -B1 --output=avail /tmp | tail -1 | tr -d '[:space:]')
[[ $disk_free =~ ^[0-9]+$ && $tmp_free =~ ^[0-9]+$ ]] || exit 14
(( disk_free >= 10*1024*1024*1024 && tmp_free >= 10*1024*1024*1024 )) || exit 15
available=$(awk '/MemAvailable:/ {printf "%.0f", $2*1024}' /proc/meminfo)
[[ $available =~ ^[0-9]+$ ]] || exit 16
active_matlab=$(pgrep -c -x MATLAB || true)
(( active_matlab == 0 )) || { echo "Other MATLAB process detected: $active_matlab" >&2; exit 17; }

if [[ $JOB == freeze_retry1 ]]; then
  [[ ! -e "$TASK/jobs/freeze_retry1" && ! -e "$TASK/prediction" \
    && ! -e "$TASK/jobs/f4031p25" ]] || exit 18
  jobdir=$TASK/jobs/freeze_retry1
  estimated_increment=$((128*1024*1024))
  call="freeze_two_dimensional_s5_p22_halfstep_prediction('$TASK')"
else
  layer=${JOB#*_}
  jobdir=$TASK/jobs/f4031p25/$layer
  if [[ $layer == model ]]; then
    estimated_increment=$((256*1024*1024))
  else
    case $layer in
      main96) L=60; ppw=12;;
      self64) L=60; ppw=12;;
      space14) L=60; ppw=14;;
      window60) L=60; ppw=10;;
      window70) L=70; ppw=10;;
      *) exit 19;;
    esac
    nodes=$(awk -v L="$L" -v p="$ppw" 'BEGIN {v=2*L*p*4031.25/(1500*8); n=int(v); if (v>n) n++; print 8*n}')
    estimated_increment=$((64*nodes*nodes+128*1024*1024))
  fi
  [[ -f "$TASK/prediction/prediction.json" && -f "$TASK/prediction/prediction_sha256.txt" ]] || exit 20
  predicted=$(cat "$TASK/prediction/prediction_sha256.txt")
  actual=$(sha256sum "$TASK/prediction/prediction.json" | cut -d' ' -f1)
  [[ $predicted == "$actual" ]] || exit 21
  call="run_two_dimensional_s5_p22_dualfreq_remote('$TASK',4031.25,'$layer',$available)"
fi
(( estimated_increment*10 <= available*7 )) || exit 22

mkdir -p -- "$(dirname "$jobdir")"
mkdir -- "$jobdir" # Atomic duplicate prevention; failed original job is untouched.
printf 'job=%s\nprior_process_seconds=%s\nprocess_upper_seconds=%s\nram_available_bytes=%s\nestimated_increment_bytes=%s\ndisk_available_bytes=%s\ntmp_available_bytes=%s\n' \
  "$JOB" "$spent" "${upper[$ix]}" "$available" "$estimated_increment" \
  "$disk_free" "$tmp_free" > "$jobdir/shell_preflight.txt"
if [[ $JOB == f4031p25_* ]]; then
  sha256sum "$TASK/prediction/prediction.json" > "$jobdir/prediction_sha256_at_start.txt"
fi
date -u +'%Y-%m-%dT%H:%M:%SZ' > "$jobdir/start_utc.txt"
start=$(date +%s)
limit=$((upper[ix]-15))
set +e
timeout --signal=TERM --kill-after=10s "${limit}s" "$MATLAB_BIN" -singleCompThread -batch \
  "addpath('$TASK/code'); $call" > "$jobdir/stdout.log" 2> "$jobdir/stderr.log"
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
if [[ $JOB == freeze_retry1 ]]; then
  [[ -f "$TASK/prediction/prediction.json" && -f "$TASK/prediction/prediction.mat" \
    && -f "$TASK/prediction/prediction_sha256.txt" ]] || exit 23
  sha256sum "$TASK/prediction/prediction.json" "$TASK/prediction/prediction.mat" \
    > "$jobdir/output_sha256.txt"
else
  [[ -f "$jobdir/result.mat" && -f "$jobdir/receipt.json" ]] || exit 24
  sha256sum "$jobdir/result.mat" "$jobdir/receipt.json" > "$jobdir/output_sha256.txt"
fi
echo "JOB_COMPLETE $JOB elapsed=$((end-start)) spent_before=$spent"
