#!/usr/bin/env bash
# Synthetic local status fixtures only. No SSH or MATLAB.
set -eu
repo=$(realpath "${1:?repository required}")
out=${2:?new output directory required}
[[ ! -e "$out" ]] || exit 2
mkdir -p "$out"; out=$(realpath "$out")
source "$repo/scripts/validation/support/status_r2_pm_session.sh"
if ! command -v python3 >/dev/null; then python3() { python "$@"; }; fi
root="$out/task_fixture"; sid=fixture_session_1
mkdir -p "$root/control" "$root/sessions/$sid" "$root/jobs/bie_window"
session="$root/sessions/$sid"; ps_mode=none
ps() {
    case "$ps_mode" in
        matlab) printf '123 00:40 Sl 99.0 1450000 MATLAB /opt/MATLAB -singleCompThread -batch run_r2_pm_seed %s\n' "$root" ;;
        wrapper) printf '123 00:01 Ss 0.0 1400 time /usr/bin/time run_r2_pm_seed %s\n' "$root" ;;
        reused) printf '123 00:40 Sl 99.0 1450000 MATLAB /opt/MATLAB run_r2_pm_seed /unrelated/task\n' ;;
    esac
}
passed=0
expect() {
    local label=$1 state=$2
    r1_c4_status_render "$root" > "$out/$label.log"
    grep -q "^STATE=$state$" "$out/$label.log"
    passed=$((passed+1))
}
expect not_started NOT_STARTED
printf '%s\n' "$sid" > "$root/control/active_session.txt"
printf '123\n' > "$session/process_group_leader.pid"
printf 'PREFLIGHT\n' > "$root/control/stage_id.txt"
ps_mode=matlab
expect startup MATLAB_INITIALIZING
grep -q 'This is not an error. Do not restart.' "$out/startup.log"
printf '{"session_id":"%s","stage_state":"C4_SOLVING","status":"RUNNING"}\n' "$sid" > "$root/checkpoint.json"
printf '{"invocations":1}\n' > "$root/jobs/bie_window/solve_invocation.json"
printf 'C4\n' > "$root/control/stage_id.txt"
expect solving RUNNING
grep -q 'TOTAL_C4_SOLVE_INVOCATIONS=1' "$out/solving.log"
printf '{"session_id":"old_session","stage_state":"COMPLETE"}\n' > "$root/checkpoint.json"
expect stale_checkpoint MATLAB_INITIALIZING
grep -q 'Previous-session' "$out/stale_checkpoint.log"
printf '{bad json}\n' > "$root/checkpoint.json"
expect invalid_checkpoint RUNNING_CHECKPOINT_UNREADABLE
ps_mode=wrapper
expect wrapper LAUNCH_WRAPPER_ACTIVE
ps_mode=reused
expect reused_pid NO_LIVE_PROCESS_NO_EXIT_RECEIPT
if grep -q '^MATLAB_PID=' "$out/reused_pid.log"; then exit 3; fi
ps_mode=none
expect no_process NO_LIVE_PROCESS_NO_EXIT_RECEIPT
printf '0\n' > "$session/matlab_exit_code.txt"
printf '304\n' > "$session/process_elapsed_s.txt"
expect exit_zero EXITED_0_VERIFY_FINAL_ARTIFACTS
grep -q 'Exit 0 alone does not establish' "$out/exit_zero.log"
grep -q '^ACTUAL_PROCESS_ELAPSED_SECONDS=304$' "$out/exit_zero.log"
printf '1\n' > "$session/matlab_exit_code.txt"
printf 'fixture failure retained\n' > "$session/stderr.log"
expect exit_failure FAILED_STOPPED
grep -q 'fixture failure retained' "$out/exit_failure.log"
printf 'not-an-exit-code\n' > "$session/matlab_exit_code.txt"
expect invalid_exit INVALID_EXIT_RECEIPT
printf '../invalid\n' > "$root/control/active_session.txt"
if r1_c4_status_render "$root" > "$out/invalid_session.log" 2>&1; then exit 3; fi
passed=$((passed+1))
printf 'PASS tests=%d MATLAB_starts=0 SSH_connections=0\n' "$passed"
