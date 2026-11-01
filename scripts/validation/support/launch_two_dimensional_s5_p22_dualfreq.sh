#!/usr/bin/env bash
set -euo pipefail

TASK=/home/misc/cwx/uwa2d_s5_p22_dualfreq_20261029_01
MATLAB_BIN=/home/misc/MATLAB2023b/bin/matlab
JOB=${1:?One frozen job ID is required}
ids=(smoke f4062p5_model f4062p5_main96 f4062p5_self64 \
  f4062p5_space14 f4062p5_window60 f4062p5_window70 freeze \
  f4031p25_model f4031p25_main96 f4031p25_self64 \
  f4031p25_space14 f4031p25_window60 f4031p25_window70)
upper=(120 180 300 300 350 300 350 90 180 300 300 350 300 350)
ix=-1
for i in "${!ids[@]}"; do
  if [[ ${ids[$i]} == "$JOB" ]]; then ix=$i; break; fi
done
(( ix >= 0 )) || { echo 'Unknown job ID' >&2; exit 2; }
[[ -d "$TASK" && $(realpath "$TASK") == "$TASK" && -x "$MATLAB_BIN" ]] || exit 3
[[ -f "$TASK/input/SHA256SUMS" ]] || exit 4
cd "$TASK"
sha256sum -c input/SHA256SUMS >/dev/null || exit 5

spent=0
for ((i=0;i<ix;i++)); do
  previous=${ids[$i]}
  if [[ $previous == freeze ]]; then
    prior_dir=$TASK/jobs/freeze
    prior_result=$TASK/prediction/prediction.json
  elif [[ $previous == smoke ]]; then
    prior_dir=$TASK/jobs/smoke
    prior_result=$prior_dir/result.mat
  else
    frequency=${previous%%_*}
    layer=${previous#*_}
    prior_dir=$TASK/jobs/$frequency/$layer
    prior_result=$prior_dir/result.mat
  fi
  [[ -f "$prior_dir/exit_code.txt" && -f "$prior_dir/process_elapsed_seconds.txt" \
    && -f "$prior_result" ]] || { echo "Incomplete prior job: $previous" >&2; exit 6; }
  [[ $(cat "$prior_dir/exit_code.txt") == 0 ]] || { echo "Prior job failed: $previous" >&2; exit 7; }
  elapsed=$(cat "$prior_dir/process_elapsed_seconds.txt")
  [[ $elapsed =~ ^[0-9]+$ ]] || exit 8
  spent=$((spent+elapsed))
done
(( spent + upper[ix] <= 1800 )) || { echo "1800 s process budget gate: spent=$spent next=${upper[$ix]}" >&2; exit 9; }

disk_free=$(df -B1 --output=avail "$TASK" | tail -1 | tr -d '[:space:]')
tmp_free=$(df -B1 --output=avail /tmp | tail -1 | tr -d '[:space:]')
[[ $disk_free =~ ^[0-9]+$ && $tmp_free =~ ^[0-9]+$ ]] || exit 10
(( disk_free >= 10*1024*1024*1024 && tmp_free >= 10*1024*1024*1024 )) || {
  echo "Disk gate: task=$disk_free tmp=$tmp_free" >&2; exit 11;
}
available=$(awk '/MemAvailable:/ {printf "%.0f", $2*1024}' /proc/meminfo)
[[ $available =~ ^[0-9]+$ ]] || exit 12
active_matlab=$(pgrep -c -x MATLAB || true)
(( active_matlab == 0 )) || { echo "Other MATLAB process detected: $active_matlab" >&2; exit 13; }
if [[ $JOB == smoke ]]; then
  jobdir=$TASK/jobs/smoke
  estimated_increment=$((128*1024*1024))
  call="run_two_dimensional_s5a_weak_bie_remote('$TASK','smoke',$available)"
elif [[ $JOB == freeze ]]; then
  jobdir=$TASK/jobs/freeze
  estimated_increment=$((128*1024*1024))
  call="freeze_two_dimensional_s5_p22_halfstep_prediction('$TASK')"
else
  frequency_id=${JOB%%_*}
  layer=${JOB#*_}
  jobdir=$TASK/jobs/$frequency_id/$layer
  if [[ $frequency_id == f4062p5 ]]; then frequency_hz=4062.5
  else frequency_hz=4031.25; fi
  if [[ $layer == model ]]; then
    estimated_increment=$((256*1024*1024))
  else
    case $layer in
      main96) L=60; ppw=12;;
      self64) L=60; ppw=12;;
      space14) L=60; ppw=14;;
      window60) L=60; ppw=10;;
      window70) L=70; ppw=10;;
      *) exit 14;;
    esac
    node_estimate=$(awk -v L="$L" -v p="$ppw" -v f="$frequency_hz" \
      'BEGIN {v=2*L*p*f/(1500*8); n=int(v); if (v>n) n++; print 8*n}')
    estimated_increment=$((64*node_estimate*node_estimate+128*1024*1024))
  fi
  call="run_two_dimensional_s5_p22_dualfreq_remote('$TASK',$frequency_hz,'$layer',$available)"
fi
(( estimated_increment*10 <= available*7 )) || {
  echo "30 percent physical RAM reserve gate: increment=$estimated_increment available=$available" >&2; exit 15;
}

mkdir -p -- "$(dirname "$jobdir")"
mkdir -- "$jobdir" # no reuse or overwrite
printf 'job=%s\ncode_manifest=%s\nmatlab=%s\nram_available_bytes=%s\nestimated_increment_bytes=%s\ndisk_available_bytes=%s\ntmp_available_bytes=%s\nprior_process_seconds=%s\nprocess_upper_seconds=%s\n' \
  "$JOB" "$TASK/input/SHA256SUMS" "$MATLAB_BIN" "$available" \
  "$estimated_increment" "$disk_free" "$tmp_free" "$spent" "${upper[$ix]}" \
  > "$jobdir/shell_preflight.txt"
if [[ $JOB == f4031p25_* ]]; then
  [[ -f "$TASK/prediction/prediction.json" && -f "$TASK/prediction/prediction_sha256.txt" ]] || exit 16
  sha256sum "$TASK/prediction/prediction.json" > "$jobdir/prediction_sha256_at_start.txt"
fi
date -u +'%Y-%m-%dT%H:%M:%SZ' > "$jobdir/start_utc.txt"
start=$(date +%s)
limit=$((upper[ix]-15))
set +e
timeout --signal=TERM --kill-after=10s "${limit}s" "$MATLAB_BIN" -singleCompThread -batch \
  "addpath('$TASK/code'); $call" \
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
if [[ $JOB == freeze ]]; then
  [[ -f "$TASK/prediction/prediction.json" && -f "$TASK/prediction/prediction.mat" \
    && -f "$TASK/prediction/prediction_sha256.txt" ]] || exit 17
  sha256sum "$TASK/prediction/prediction.json" "$TASK/prediction/prediction.mat" \
    > "$jobdir/output_sha256.txt"
else
  [[ -f "$jobdir/result.mat" && -f "$jobdir/receipt.json" ]] || exit 18
  sha256sum "$jobdir/result.mat" "$jobdir/receipt.json" > "$jobdir/output_sha256.txt"
fi
echo "JOB_COMPLETE $JOB elapsed=$((end-start)) spent_before=$spent"
