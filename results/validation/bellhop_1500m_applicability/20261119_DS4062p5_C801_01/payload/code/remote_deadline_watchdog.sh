#!/bin/bash
set -euo pipefail
OUT=$1
owner=$2
has_owned_live_pid() {
 local p=$1 idfile=$2 saved live savedtail livetail
 test -f "$idfile" && test -r "/proc/$p/stat" || return 1
 IFS= read -r saved < "$idfile" || return 1
 IFS= read -r live < "/proc/$p/stat" || return 1
 savedtail=${saved##*) }; livetail=${live##*) }
 set -- $savedtail; local oldstart=${20}
 set -- $livetail; local newstart=${20} session=${4}
 test "$oldstart" = "$newstart" && test "$session" = "$owner"
}
sleep 170
valid=0
if has_owned_live_pid "$owner" "$OUT/MATLAB_identity.txt"; then valid=1; fi
for slot in D S; do
 if test -f "$OUT/$slot/case_data/pid.txt"; then
  IFS= read -r native < "$OUT/$slot/case_data/pid.txt" || true
  if has_owned_live_pid "$native" "$OUT/$slot/case_data/process_identity.txt"; then valid=1; fi
 fi
done
if test "$valid" = 1; then
 kill -TERM -- "-$owner" 2>/dev/null || true
 sleep .2
 kill -KILL -- "-$owner" 2>/dev/null || true
 printf 'INDEPENDENT_REMOTE_170S_WATCHDOG\n' > "$OUT/remote_watchdog_fired.txt"
fi