#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "$0")" && pwd)
rc=0
timeout --foreground --signal=TERM --kill-after=1 120s bash "$HERE/native_leaf_exec.sh" "$1" "$2" || rc=$?
printf '%s\n' "$rc" > native_guard_exit.txt
exit "$rc"