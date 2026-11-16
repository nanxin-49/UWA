#!/usr/bin/env bash
# One MATLAB process, no cumulative/forecast/runtime timeout. RAM and
# numerical/audit failures stop the task. Only the user invokes start/resume.
set -u
mode=${1:---dry-run}
root=$(realpath -e "${2:?package root required}") || exit 2
script="$root/code/launch_r1_bie_c4_continuation.sh"
source "$root/code/r1_bie_runtime.sh"
if [[ "$mode" == --dry-run ]]; then
    cd "$root" || exit 2
    python3 "$root/code/r1_bie_c4_preflight.py" inputs "$root" || exit 3
    if [[ -f package.sha256 ]]; then sha256sum -c package.sha256 || exit 3; fi
    printf 'DRY_RUN_INPUTS_ONLY: no MATLAB, no solve, no network; readiness requires separate ready check.\n'
    printf 'Would use: /home/misc/MATLAB/R2025b/bin/matlab -singleCompThread -batch run_r1_bie_c4_continuation\n'
    printf 'Stages: verify migrated C0-C3; ONLY C4 70/12/96; final audits; reference=C3.\n'
    printf 'No C0-C3 solve branch. C4 invocation marker forbids automatic re-solve.\n'
    exit 0
fi
[[ "$root" =~ ^/home/misc/cwx/r1_bie_manual_c4_policy_v2_[A-Za-z0-9_-]+$ ]] || { echo 'Invalid task-owned root' >&2; exit 2; }
cd "$root" || exit 2
if [[ "$mode" == status ]]; then
    if [[ -f checkpoint.json ]]; then cat checkpoint.json; else echo 'No MATLAB checkpoint yet.'; fi
    if [[ -f control/active_session.txt ]]; then
        sid=$(cat control/active_session.txt)
        [[ "$sid" =~ ^[A-Za-z0-9_-]+$ ]] || exit 2
        printf '\nLatest session: %s\n' "$sid"
        if [[ -f "sessions/$sid/matlab_exit_code.txt" ]]; then
            printf 'REAL MATLAB EXIT CODE: '; cat "sessions/$sid/matlab_exit_code.txt"
            [[ ! -f "sessions/$sid/external_stop.txt" ]] || cat "sessions/$sid/external_stop.txt"
        else
            echo 'No exit receipt yet; inspect session logs/lock (not evidence of numerical PASS).'
        fi
        printf 'stdout: %s/sessions/%s/stdout.log\nstderr: %s/sessions/%s/stderr.log\n' "$root" "$sid" "$root" "$sid"
    fi
    exit 0
fi
for command in sha256sum flock setsid ps awk date; do command -v "$command" >/dev/null || exit 3; done
if [[ "$mode" == start ]]; then
    python3 "$root/code/r1_bie_c4_preflight.py" ready "$root" || exit 78
    [[ -f package.sha256 ]] || exit 78
    mkdir -p sessions control
    sid="$(date -u +%Y%m%dT%H%M%S_%N)_$$"
    mkdir "sessions/$sid" || exit 3
    nohup bash "$script" run "$root" "$sid" > "sessions/$sid/launcher_stdout.log" \
        2> "sessions/$sid/launcher_stderr.log" < /dev/null &
    echo "START_REQUEST session=$sid controller_pid=$!"
    echo "Check: bash '$script' status '$root'"
    exit 0
fi
[[ "$mode" == run ]] || { echo 'Use --dry-run, start or status' >&2; exit 2; }
sid=${3:?session ID required}
[[ "$sid" =~ ^[A-Za-z0-9_-]+$ && -d "sessions/$sid" ]] || exit 2
session="sessions/$sid"
exec 9> workflow.lock
if ! flock -n 9; then echo 75 > "$session/matlab_exit_code.txt"; echo 'Another task session holds the lock' >&2; exit 75; fi
first_run=1
if [[ -f control/active_session.txt ]]; then
    first_run=0
    previous=$(cat control/active_session.txt)
    [[ "$previous" =~ ^[A-Za-z0-9_-]+$ ]] || exit 2
    if [[ -f "sessions/$previous/process_group_leader.pid" ]]; then
        previous_pid=$(cat "sessions/$previous/process_group_leader.pid")
        [[ "$previous_pid" =~ ^[0-9]+$ ]] || exit 2
        if [[ -n "$(ps -o pid= --sid "$previous_pid" 2>/dev/null)" ]]; then
            echo 'Previous task-owned process session is still alive; refusing duplicate MATLAB' >&2
            echo 75 > "$session/matlab_exit_code.txt"; exit 75
        fi
    fi
fi
printf '%s\n' "$sid" > control/active_session.txt
mkdir -p control/resources jobs finalization
printf 'PREFLIGHT\n' > control/stage_id.txt
sha256sum -c package.sha256 > "$session/input_verification.log" || { echo 76 > "$session/matlab_exit_code.txt"; exit 76; }
python3 "$root/code/r1_bie_c4_preflight.py" ready "$root" >> "$session/input_verification.log" || { echo 78 > "$session/matlab_exit_code.txt"; exit 78; }
[[ -x /home/misc/MATLAB/R2025b/bin/matlab && -x /usr/bin/time ]] || { echo 77 > "$session/matlab_exit_code.txt"; exit 77; }
config_stop() {
    echo "RUNTIME_CONFIGURATION_STOP: $*" > "$session/external_stop.txt"
    echo 78 > "$session/matlab_exit_code.txt"
    exit 78
}
available=$(r1_read_available /proc/meminfo) || config_stop 'system memory read failed'
[[ -f control/initial_available_bytes.txt ]] || config_stop 'inherited initial memory file missing'
r1_prepare_runtime "$available" control/initial_available_bytes.txt "$session/runtime.json" \
    "$root/code/r1_bie_runtime_json.py" || config_stop 'memory validation or JSON serialization/parse failed'
# Python validates the independently named fields; Bash computes the same max.
initial=$(cat control/initial_available_bytes.txt) || config_stop 'initial memory read failed'
initial=$(r1_uint "$initial") || config_stop 'invalid initial memory'
if ((available>initial)); then baseline=$available; else baseline=$initial; fi
if ((available*7/10 < 12*1024*1024*1024)); then config_stop '12 GiB RAM admission estimate does not fit'; fi
start=$(date +%s); echo "$start" > "$session/start_epoch.txt"
code="addpath('$root/code'); run_r1_bie_c4_continuation('$root','run','$sid');"
printf '%s\n' "$code" > "$session/command.txt"
export R1_C4_ACTIVE_SESSION="$sid"
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 setsid \
    /usr/bin/time -v -o "$session/time.txt" \
    /home/misc/MATLAB/R2025b/bin/matlab -singleCompThread -batch "$code" \
    > "$session/stdout.log" 2> "$session/stderr.log" &
taskpid=$!; echo "$taskpid" > "$session/process_group_leader.pid"
stop_owned_group() {
    local pgid
    pgid=$(ps -o pgid= -p "$taskpid" | tr -d ' ')
    if [[ "$pgid" == "$taskpid" ]]; then
        kill -TERM -- "-$taskpid" 2>/dev/null || true
        sleep 3
        kill -KILL -- "-$taskpid" 2>/dev/null || true
    fi
}
trap 'echo USER_SIGNAL_STOP > "$session/external_stop.txt"; stop_owned_group' INT TERM
while kill -0 "$taskpid" 2>/dev/null; do
    mem=$(r1_read_available /proc/meminfo) || { echo SYSTEM_MEMORY_READ_STOP > "$session/external_stop.txt"; stop_owned_group; break; }
    stage=$(cat control/stage_id.txt 2>/dev/null || echo UNKNOWN)
    [[ "$stage" =~ ^[A-Z0-9_]+$ ]] || stage=UNKNOWN
    rss=$(ps -o rss= --sid "$taskpid" 2>/dev/null | awk 'BEGIN{s=0}{s+=$1}END{print s}')
    printf '%s,%s,%s,%s\n' "$(date +%s)" "$mem" "$rss" "$stage" >> "$session/resources.csv"
    printf '%s,%s,%s\n' "$(date +%s)" "$mem" "$rss" >> "control/resources/$stage.csv"
    if ! awk -v a="$mem" -v b="$baseline" 'BEGIN{exit !(a>=.3*b)}'; then
        echo RAM_RESERVE_STOP > "$session/external_stop.txt"; stop_owned_group; break
    fi
    sleep 1
done
wait "$taskpid"; rc=$?
trap - INT TERM
echo "$rc" > "$session/matlab_exit_code.txt"
echo "$(($(date +%s)-start))" > "$session/process_elapsed_s.txt"
date -Is > "$session/finished.txt"
find jobs finalization "$session" control -type f ! -name output_sha256.txt -print0 \
    | sort -z | xargs -0 sha256sum > "$session/output_sha256.txt"
exit "$rc"
