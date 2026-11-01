#!/usr/bin/env bash
set -u
root=/home/misc/cwx/deepsea3d_65grid_four_20261027
case "${1:-}" in
  4750) freq=4750; label=4750; remaining=4; run_id=02 ;;
  7250) freq=7250; label=7250; remaining=3; run_id=01 ;;
  4125) freq=4125; label=4125; remaining=2; run_id=01 ;;
  6562.5) freq=6562.5; label=6562.5; remaining=1; run_id=01 ;;
  *) exit 64 ;;
esac
out="$root/results/validation/deepsea3d_65grid_four_20261027/run_f${label}_${run_id}"
test -d "$root" || exit 66
test ! -e "$out" || exit 65
cd "$root" || exit 66
sha256sum -c remote_input.sha256 > /dev/null || exit 67
mem_total=$(awk '/^MemTotal:/ {print $2}' /proc/meminfo)
mem_available=$(awk '/^MemAvailable:/ {print $2}' /proc/meminfo)
awk -v a="$mem_available" -v t="$mem_total" 'BEGIN {exit !(a-3145728 >= 0.30*t)}' || exit 68
disk_available=$(df -Pk "$root" | tail -1 | awk '{print $4}')
test "$disk_available" -ge 2097152 || exit 68
used=95 # Conservative charge for this Goal's MATLAB R2023b version probe.
for item in "$root"/results/validation/deepsea3d_65grid_four_20261027/run_f*/process_wall_s.txt; do
  if test -f "$item"; then n=$(cat "$item"); used=$((used+n)); fi
done
test $((used+remaining*120+25)) -le 720 || exit 69
mkdir "$out" || exit 66
date -u '+%Y-%m-%dT%H:%M:%SZ' > "$out/start_utc.txt"
cat /proc/meminfo > "$out/mem_before.txt"
cat /proc/loadavg > "$out/load_before.txt"
df -Pk "$root" > "$out/disk_before.txt"
start=$(date +%s)
timeout --signal=TERM --kill-after=5s 150s \
  /usr/bin/time -v /home/misc/MATLAB2023b/bin/matlab -batch \
  "addpath('$root'); setup_vertical_project(); run_deepsea3d_65grid_four_frequency($freq);" \
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
