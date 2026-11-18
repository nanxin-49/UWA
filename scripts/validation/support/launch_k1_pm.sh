#!/usr/bin/env bash
set -euo pipefail
action=${1:---dry-run}; root=${2:?absolute package root}
export PYTHONDONTWRITEBYTECODE=1
if [[ "$action" == --dry-run ]]; then
  python3 "$root/code/k1_pm_preflight.py" "$root"
  printf '%s\n' 'K1_DRY_RUN_ONLY: one MATLAB process, one computational thread; no MATLAB/network/solve invoked.'
  exit 0
fi
[[ "$root" =~ ^/home/misc/cwx/k1_pm_[A-Za-z0-9_-]+$ && -d "$root" && ! -L "$root" ]] || exit 64
root=$(realpath -e "$root")
case "$action" in
  start)
    python3 "$root/code/k1_pm_preflight.py" "$root"
    [[ ! -L /home/misc/cwx/k1_pm_serial.lock ]] || exit 65
    exec 9>/home/misc/cwx/k1_pm_serial.lock
    flock -n 9 || { echo 'Another task-owned K1 session holds the single process slot'; exit 75; }
    session="$root/sessions/$(date -u +%Y%m%dT%H%M%S)_$$"
    mkdir -p "$session"
    nohup bash "$root/code/launch_k1_pm.sh" run "$root" "$session" 9>&9 </dev/null >"$session/launcher.log" 2>&1 &
    printf '%s\n' "$session" >"$root/latest_session.txt"
    printf 'SESSION=%s\nCONTROLLER_PID=%s\n' "$session" "$!"
    ;;
  run)
    session=${3:?session}; [[ "$session" == "$root"/sessions/* && -d "$session" ]] || exit 64
    flock -n 9 || exit 75
    exec python3 -B "$root/code/k1_pm_supervisor.py" "$root" "$session"
    ;;
  status)
    if [[ ! -f "$root/latest_session.txt" ]]; then echo 'NOT_STARTED'; exit 0; fi
    session=$(cat "$root/latest_session.txt"); [[ "$session" == "$root"/sessions/* ]] || exit 65
    printf 'SESSION=%s\n' "$session"
    if [[ -f "$session/exit.json" ]]; then cat "$session/exit.json"; else echo 'NO_EXIT_RECEIPT_YET: starting/running/interrupted; never infer PASS'; fi
    [[ ! -f "$root/checkpoint.json" ]] || cat "$root/checkpoint.json"
    ;;
  fetch)
    [[ ! -L /home/misc/cwx/k1_pm_serial.lock ]] || exit 65
    exec 9>/home/misc/cwx/k1_pm_serial.lock
    flock -n 9 || { echo 'K1 slot active; Fetch is only for stopped tasks'; exit 75; }
    python3 -B "$root/code/k1_pm_preflight.py" "$root" archive
    ;;
  *) exit 64;;
esac
