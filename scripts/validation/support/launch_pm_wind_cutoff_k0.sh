#!/usr/bin/env bash
set -euo pipefail
action=${1:?action}; root=${2:?root}
[[ "$root" =~ ^/home/misc/cwx/pm_wind_cutoff_k0_[A-Za-z0-9_-]+$ ]] || exit 64
[[ -d "$root" && ! -L "$root" ]] || exit 66
root=$(realpath "$root")
case "$action" in
  start)
    python3 "$root/code/pm_k0_preflight.py" "$root"
    # One shared K0 slot across task-owned uploads; no lock-file deletion.
    exec 9>/home/misc/cwx/pm_wind_cutoff_k0_serial.lock
    flock -n 9 || { echo 'K0 slot busy; no new MATLAB started'; exit 75; }
    session="$root/sessions/$(date -u +%Y%m%dT%H%M%S)_$$"
    mkdir -p "$session"
    nohup bash "$root/code/launch_pm_wind_cutoff_k0.sh" run "$root" "$session" 9>&9 </dev/null >"$session/launcher.log" 2>&1 &
    pid=$!
    printf '%s\n' "$session" >"$root/latest_session.txt"
    printf 'SESSION=%s\nLAUNCHER_PID=%s\n' "$session" "$pid"
    ;;
  run)
    session=${3:?session}
    [[ "$session" == "$root"/sessions/* && -d "$session" ]] || exit 64
    # FD9 is inherited from start and retained throughout the supervisor.
    flock -n 9 || exit 75
    exec python3 "$root/code/pm_k0_supervisor.py" "$root" "$session"
    ;;
  status)
    if [[ ! -f "$root/latest_session.txt" ]]; then echo 'NOT_STARTED'; exit 0; fi
    session=$(cat "$root/latest_session.txt")
    [[ "$session" == "$root"/sessions/* && -d "$session" ]] || exit 65
    printf 'SESSION=%s\n' "$session"
    if [[ -f "$session/exit.json" ]]; then cat "$session/exit.json"; else
      echo 'NO_EXIT_RECEIPT_YET (starting/running or interrupted; do not infer PASS)'
      [[ ! -f "$session/identity.json" ]] || cat "$session/identity.json"
    fi
    ;;
  fetch)
    exec 9>/home/misc/cwx/pm_wind_cutoff_k0_serial.lock
    flock -n 9 || { echo 'K0 active; stopped-task Fetch only'; exit 75; }
    python3 "$root/code/pm_k0_preflight.py" "$root" archive
    ;;
  *) exit 64;;
esac
