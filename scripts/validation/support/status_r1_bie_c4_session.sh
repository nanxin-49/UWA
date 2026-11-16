#!/usr/bin/env bash
# Read-only status presentation. No locks, writes, signals, or MATLAB launch.
# Can be sent by the local manager through SSH without changing a sealed task.
r1_status_json_field() {
    python3 -c 'import json,sys; v=json.load(open(sys.argv[1],encoding="utf-8-sig")); x=v.get(sys.argv[2],""); print(x if isinstance(x,(str,int,float)) else "")' "$1" "$2" 2>/dev/null
}
r1_c4_status_render() {
    local root=$1 sid session leader='' row pid elapsed procstate cpu rss name args
    local live=0 wrapper=0 checkpoint_session='' checkpoint_state='' stage='UNKNOWN' rc=''
    local phase='UNKNOWN' epoch='' now age='unknown' baseline='' sample=''
    printf 'TASK_ROOT=%s\n' "$root"
    if [[ ! -f "$root/control/active_session.txt" ]]; then
        printf 'STATE=NOT_STARTED\nNo active session has been recorded; no computation was started by this query.\n'
        return 0
    fi
    sid=$(cat "$root/control/active_session.txt") || return 2
    [[ "$sid" =~ ^[A-Za-z0-9_-]+$ ]] || { printf 'STATE=INVALID_SESSION_ID\n' >&2; return 2; }
    session="$root/sessions/$sid"
    [[ -d "$session" ]] || { printf 'STATE=SESSION_DIRECTORY_MISSING\n' >&2; return 2; }
    printf 'SESSION=%s\n' "$sid"
    if [[ -f "$session/start_epoch.txt" ]]; then
        epoch=$(cat "$session/start_epoch.txt")
        if [[ "$epoch" =~ ^[0-9]+$ ]]; then now=$(date +%s); age=$((now-epoch)); fi
    fi
    printf 'SESSION_AGE_SECONDS=%s\n' "$age"
    if [[ -f "$session/process_elapsed_s.txt" ]]; then
        local actual_elapsed
        actual_elapsed=$(cat "$session/process_elapsed_s.txt")
        [[ ! "$actual_elapsed" =~ ^[0-9]+$ ]] || printf 'ACTUAL_PROCESS_ELAPSED_SECONDS=%s\n' "$actual_elapsed"
    fi
    if [[ -f "$root/control/stage_id.txt" ]]; then
        stage=$(cat "$root/control/stage_id.txt")
        [[ "$stage" =~ ^[A-Z0-9_]+$ ]] || stage=UNKNOWN
    fi
    if [[ -f "$root/checkpoint.json" ]]; then
        checkpoint_session=$(r1_status_json_field "$root/checkpoint.json" session_id) || checkpoint_session='INVALID_JSON'
        if [[ "$checkpoint_session" == "$sid" ]]; then
            checkpoint_state=$(r1_status_json_field "$root/checkpoint.json" stage_state) || checkpoint_state='INVALID_JSON'
            printf 'CHECKPOINT_STATE=%s\n' "$checkpoint_state"
            printf 'CHECKPOINT_STATUS=%s\n' "$(r1_status_json_field "$root/checkpoint.json" status)"
        else
            printf 'CHECKPOINT_NOTE=Previous-session or invalid checkpoint; not evidence of progress in this session.\n'
        fi
    fi
    if [[ -f "$session/process_group_leader.pid" ]]; then
        leader=$(cat "$session/process_group_leader.pid")
        if [[ "$leader" =~ ^[0-9]+$ ]]; then
            # Only inspect this recorded SID, and display processes whose argv
            # identifies this exact task and entrypoint. Ignore unrelated/reused PIDs.
            while read -r pid elapsed procstate cpu rss name args; do
                [[ "$args" == *"$root"* && "$args" == *run_r1_bie_c4_continuation* ]] || continue
                if [[ "$name" == MATLAB || "$name" == matlab ]]; then
                    live=1
                    printf 'MATLAB_PID=%s ELAPSED=%s PROCESS_STATE=%s CPU_PERCENT=%s RSS_KIB=%s\n' "$pid" "$elapsed" "$procstate" "$cpu" "$rss"
                else
                    wrapper=1
                fi
            done < <(ps --sid "$leader" -o pid=,etime=,stat=,pcpu=,rss=,comm=,args= 2>/dev/null)
        else
            printf 'PROCESS_NOTE=Invalid recorded PID; process state is unknown.\n'
        fi
    fi
    if [[ -f "$session/matlab_exit_code.txt" ]]; then
        rc=$(cat "$session/matlab_exit_code.txt")
        if [[ ! "$rc" =~ ^[0-9]+$ ]]; then
            phase='INVALID_EXIT_RECEIPT'
        elif [[ "$rc" != 0 ]]; then
            phase='FAILED_STOPPED'
        elif ((live)); then
            phase='EXIT_RECEIPT_PROCESS_CONFLICT'
        else
            phase='EXITED_0_VERIFY_FINAL_ARTIFACTS'
        fi
        printf 'STATE=%s\nREAL_EXIT_CODE=%s\n' "$phase" "$rc"
        [[ "$rc" != 0 ]] || printf 'Exit 0 alone does not establish numerical/reference PASS; inspect final audits and seals.\n'
    elif ((live)); then
        if [[ -n "$checkpoint_state" && "$checkpoint_session" == "$sid" ]]; then
            phase='RUNNING'
            printf 'STATE=RUNNING\nSTAGE=%s\nMATLAB is active; no exit receipt is expected while it runs.\n' "$stage"
        elif [[ "$checkpoint_session" == INVALID_JSON ]]; then
            phase='RUNNING_CHECKPOINT_UNREADABLE'
            printf 'STATE=RUNNING_CHECKPOINT_UNREADABLE\nMATLAB is active but checkpoint JSON could not be parsed; inspect the recorded file.\n'
        else
            phase='MATLAB_INITIALIZING'
            printf 'STATE=MATLAB_INITIALIZING\nSTAGE=%s\nMATLAB is active before its first checkpoint (startup/input/SHA validation). This is not an error. Do not restart.\n' "$stage"
        fi
    elif ((wrapper)); then
        phase='LAUNCH_WRAPPER_ACTIVE'
        printf 'STATE=LAUNCH_WRAPPER_ACTIVE\nThe task-owned launch wrapper is present; MATLAB startup has not yet been observed. Do not start a duplicate.\n'
    else
        phase='NO_LIVE_PROCESS_NO_EXIT_RECEIPT'
        printf 'STATE=NO_LIVE_PROCESS_NO_EXIT_RECEIPT\nNo matching live MATLAB or real exit receipt was found. This is UNKNOWN, not PASS; inspect launcher logs before any restart.\n'
    fi
    if [[ -f "$session/resources.csv" ]]; then
        sample=$(tail -n 1 "$session/resources.csv")
        printf 'LATEST_RESOURCE_SAMPLE_EPOCH_AVAILABLE_BYTES_RSS_KIB_STAGE=%s\n' "$sample"
        if [[ -f "$session/runtime.json" ]]; then
            baseline=$(r1_status_json_field "$session/runtime.json" RAM_baseline_bytes) || baseline='unknown'
            printf 'RAM_BASELINE_BYTES=%s RESERVE_FRACTION=0.30\n' "$baseline"
        fi
    fi
    if [[ -f "$root/jobs/bie_window/solve_invocation.json" ]]; then
        printf 'TOTAL_C4_SOLVE_INVOCATIONS=%s\n' "$(r1_status_json_field "$root/jobs/bie_window/solve_invocation.json" invocations)"
    fi
    [[ ! -f "$session/external_stop.txt" ]] || { printf 'EXTERNAL_STOP='; cat "$session/external_stop.txt"; }
    printf 'STDOUT=%s/stdout.log\nSTDERR=%s/stderr.log\n' "$session" "$session"
    printf 'LAUNCHER_STDOUT=%s/launcher_stdout.log\nLAUNCHER_STDERR=%s/launcher_stderr.log\n' "$session" "$session"
    if [[ "$phase" == FAILED_STOPPED || "$phase" == NO_LIVE_PROCESS_NO_EXIT_RECEIPT || "$phase" == INVALID_EXIT_RECEIPT ]]; then
        for name in launcher_stderr.log stderr.log; do
            if [[ -s "$session/$name" ]]; then printf 'LAST_ERROR_LINES=%s\n' "$name"; tail -n 16 "$session/$name"; fi
        done
    fi
}
r1_c4_status_main() {
    local root
    root=$(realpath -e "${1:?task root required}") || return 2
    [[ "$root" =~ ^/home/misc/cwx/r1_bie_manual_c4_policy_v2_[A-Za-z0-9_-]+$ ]] || { printf 'Invalid task-owned status root\n' >&2; return 2; }
    r1_c4_status_render "$root"
}
if [[ "${BASH_SOURCE[0]-}" == "$0" ]]; then r1_c4_status_main "$@"; fi
