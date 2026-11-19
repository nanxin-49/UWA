#!/bin/bash
set -euo pipefail
OUT=/home/misc/cwx/bellhop_1500m_DS4062p5_C801_20261119_01
ROOT=/home/misc/cwx/bellhop_rotated_C_transfer_20261116_01/Explain
cd "$OUT"
export P1_START_EPOCH=$(date +%s.%N)
printf '%s\n' "$P1_START_EPOCH" > remote_start_epoch.txt
free -b > resources_before.txt
df -B1 "$OUT" >> resources_before.txt
cat /proc/loadavg > load_before.txt
command -v setsid timeout sha256sum tar > required_commands.txt
sha256sum -c input_SHA256SUMS > input_SHA_checks.log
printf '%s  %s\n' 5d87293a27e7fafc3aedefb569f3f35065294f42426503c4881d6eb55e6769e9 "$ROOT/results/validation/bellhop_rotated_fig12b/20261116_linux_C_01/bin/bellhop_rotated_receiver64_2020" | sha256sum -c - > ELF_SHA_check.log
launch_available=$(awk '/MemAvailable:/ {printf "%.0f",$2*1024}' /proc/meminfo)
awk -v a="$launch_available" 'BEGIN {exit !(a>=1610612736)}'
freebytes=$(df -B1 --output=avail "$OUT" | tail -1)
awk -v a="$freebytes" 'BEGIN {exit !(a>=21474836480)}'
reserve=$(awk -v a="$launch_available" 'BEGIN {v=.3*a;if(v<268435456)v=268435456;printf "%.0f",v}')
setsid bash "$OUT/code/matlab_exec.sh" > MATLAB_stdout.log 2> MATLAB_stderr.log &
owned=$!
printf '%s\n' "$owned" > MATLAB_owned_session_pid.txt
for t in $(seq 1 100); do test -f MATLAB_identity.txt && break; sleep .01; done
sid=$(ps -o sid= -p "$owned" | tr -d ' ')
test "$sid" = "$owned"
setsid bash "$OUT/code/remote_deadline_watchdog.sh" "$OUT" "$owned" &
watchdog=$!
printf '%s\n' "$watchdog" > remote_watchdog_pid.txt
rc=0; reason=''; min_available=$launch_available; SECONDS=0; lastdisk=-1
while kill -0 "$owned" 2>/dev/null; do
 available=$(awk '/MemAvailable:/ {printf "%.0f",$2*1024}' /proc/meminfo)
 if (( available < min_available )); then min_available=$available; fi
 if (( available < reserve )); then reason=SHARED_RAM_RESERVE; fi
 if (( SECONDS >= 170 )); then reason=REMOTE_170S_DEADLINE; fi
 if (( SECONDS != lastdisk )); then
  lastdisk=$SECONDS; bytes=$(du -sb "$OUT" | cut -f1); freebytes=$(df -B1 --output=avail "$OUT" | tail -1)
  if (( bytes > 2147483648 )); then reason=OUTPUT_CAP; fi
  if (( freebytes < 21474836480 )); then reason=FREE_DISK_GATE; fi
 fi
 if test -n "$reason"; then
  printf '%s\n' "$reason" > remote_abort_reason.txt
  kill -TERM -- "-$owned" 2>/dev/null || true
  sleep .2
  kill -KILL -- "-$owned" 2>/dev/null || true
  break
 fi
 sleep .1
done
wait "$owned" || rc=$?
kill -TERM -- "-$watchdog" 2>/dev/null || true
wait "$watchdog" 2>/dev/null || true
printf '%s\n' "$rc" > MATLAB_exit_code.txt
printf '{"launch_available_bytes":%s,"minimum_observed_available_bytes":%s,"reserve_bytes":%s,"abort_reason":"%s","owned_session_pid":%s}\n' "$launch_available" "$min_available" "$reserve" "$reason" "$owned" > aggregate_resource_receipt.json
free -b > resources_after.txt
du -sb "$OUT" > output_bytes.txt
read -r finalbytes rest < output_bytes.txt
if (( finalbytes > 1073741824 )); then printf "OUTPUT_CAP_INCLUDING_RETURN_TAR\n" > remote_abort_reason.txt; exit 124; fi
date +%s.%N > remote_end_epoch.txt
find . -type f ! -name return_SHA256SUMS ! -name returned.tar -print0 | sort -z | xargs -0 sha256sum > return_SHA256SUMS
tar --exclude=./returned.tar -cf returned.tar .
if test -n "$reason"; then exit 124; fi
exit "$rc"