#!/usr/bin/env bash
set -u
root=/home/misc/cwx/deepsea3d_65grid_four_20261027
result="$root/results/validation/deepsea3d_65grid_four_20261027"
out="$result/run_render_01"
test -d "$root" && test -f "$result/four_point_blind_metrics.csv" && test ! -e "$out" || exit 66
cd "$root" || exit 66
sha256sum -c remote_input.sha256 > /dev/null || exit 67
test "$(sha256sum scripts/reporting/render_deepsea3d_65grid_four.m | cut -d' ' -f1)" = 5c47593b051073455f8f2509f9fb952222be6fb3018c0ae5ec9a2b43da57343b || exit 67
test "$(sha256sum results/validation/deepsea3d_65grid_four_20261027/four_point_blind_metrics.csv | cut -d' ' -f1)" = 58c8947faa25340d7201fadfcf44c456b1e39cd34710dc676c270da13c0479a1 || exit 67
used=95
for item in "$result"/run_f*/process_wall_s.txt "$result"/run_audit_01/process_wall_s.txt; do
  if test -f "$item"; then n=$(cat "$item"); used=$((used+n)); fi
done
test $((used+120)) -le 720 || exit 69
mem_total=$(awk '/^MemTotal:/ {print $2}' /proc/meminfo)
mem_available=$(awk '/^MemAvailable:/ {print $2}' /proc/meminfo)
awk -v a="$mem_available" -v t="$mem_total" 'BEGIN {exit !(a-3145728 >= 0.30*t)}' || exit 68
disk_available=$(df -Pk "$root" | tail -1 | awk '{print $4}')
test "$disk_available" -ge 2097152 || exit 68
mkdir "$out" || exit 66
date -u '+%Y-%m-%dT%H:%M:%SZ' > "$out/start_utc.txt"
cat /proc/meminfo > "$out/mem_before.txt"
cat /proc/loadavg > "$out/load_before.txt"
df -Pk "$root" > "$out/disk_before.txt"
start=$(date +%s)
timeout --signal=TERM --kill-after=5s 120s \
  /usr/bin/time -v /home/misc/MATLAB2023b/bin/matlab -batch \
  "addpath('$root'); setup_vertical_project(); render_deepsea3d_65grid_four('$root');" \
  > "$out/stdout.log" 2> "$out/stderr.log"
status=$?
end=$(date +%s)
printf '%s\n' "$status" > "$out/exit_status.txt"
printf '%s\n' "$((end-start))" > "$out/process_wall_s.txt"
date -u '+%Y-%m-%dT%H:%M:%SZ' > "$out/end_utc.txt"
cat /proc/meminfo > "$out/mem_after.txt"
cat /proc/loadavg > "$out/load_after.txt"
df -Pk "$root" > "$out/disk_after.txt"
exit "$status"
