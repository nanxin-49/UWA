#!/usr/bin/env bash
set -euo pipefail
action=${1:?action}; root=${2:?root}
[[ "$root" =~ ^/home/misc/cwx/d1_depth_[A-Za-z0-9_-]+$ ]] || exit 64
[[ -d "$root" && ! -L "$root" && "$(realpath "$root")" == "$root" ]] || exit 66
case "$action" in
 start)
  python3 "$root/code/d1_package.py" verify "$root"
  exec 9>/home/misc/cwx/d1_depth_serial.lock
  flock -n 9 || { echo 'D1 serial slot busy'; exit 75; }
  [[ ! -e "$root/launch_started.json" && ! -e "$root/output" ]] || { echo 'Already attempted; inspect Status/Fetch, no repeat'; exit 73; }
  session="$root/sessions/$(date -u +%Y%m%dT%H%M%S)_$$"
  mkdir -p "$session"
  printf '{"status":"LAUNCH_ATTEMPTED","session":"%s"}\n' "$session" >"$root/launch_started.json"
  printf '%s\n' "$session" >"$root/latest_session.txt"
  nohup bash "$root/code/launch_d1_depth.sh" run "$root" "$session" 9>&9 </dev/null >"$session/launcher.log" 2>&1 &
  printf 'SESSION=%s\nLAUNCHER_PID=%s\n' "$session" "$!"
  ;;
 run)
  session=${3:?session}
  [[ "$session" == "$root"/sessions/* && -d "$session" ]] || exit 64
  flock -n 9 || exit 75
  exec python3 "$root/code/d1_supervisor.py" "$root" "$session"
  ;;
 status)
  if [[ ! -f "$root/latest_session.txt" ]]; then echo NOT_STARTED; exit 0; fi
  session=$(cat "$root/latest_session.txt")
  [[ "$session" == "$root"/sessions/* && -d "$session" && ! -L "$session" ]] || exit 65
  printf 'SESSION=%s\n' "$session"
  if [[ -f "$session/exit.json" ]]; then cat "$session/exit.json"; else
   echo 'NO_EXIT_RECEIPT_YET: starting/running/interrupted; no PASS inferred'
   [[ ! -f "$session/identity.json" ]] || cat "$session/identity.json"
  fi
  [[ ! -f "$root/output/checkpoint.json" ]] || python3 -c 'import json,sys; c=json.load(open(sys.argv[1])); print(json.dumps({k:c.get(k) for k in ["status","completed_steps","elapsed_s","next_allowed_step"]}))' "$root/output/checkpoint.json"
  ;;
 fetch)
  exec 9>/home/misc/cwx/d1_depth_serial.lock
  flock -n 9 || { echo 'D1 active: Fetch refused'; exit 75; }
  python3 "$root/code/d1_package.py" archive "$root"
  ;;
 *) exit 64;;
esac
