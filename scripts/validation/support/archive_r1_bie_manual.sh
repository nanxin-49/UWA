#!/usr/bin/env bash
# User-triggered retrieval snapshot. Refuses to snapshot a running task.
set -eu
root=$(realpath -e "${1:?task root required}")
[[ "$root" =~ ^/home/misc/cwx/r1_bie_manual_[A-Za-z0-9_-]+$ ]] || exit 2
cd "$root"
exec 9> workflow.lock
flock -n 9 || { echo 'Task is still running; retrieval waits for a stopped task.' >&2; exit 75; }
if [[ -f control/active_session.txt ]]; then
    sid=$(cat control/active_session.txt); [[ "$sid" =~ ^[A-Za-z0-9_-]+$ ]] || exit 2
    if [[ -f "sessions/$sid/process_group_leader.pid" ]]; then
        leader=$(cat "sessions/$sid/process_group_leader.pid"); [[ "$leader" =~ ^[0-9]+$ ]] || exit 2
        [[ -z "$(ps -o pid= --sid "$leader" 2>/dev/null)" ]] || { echo 'Previous MATLAB session is still alive' >&2; exit 75; }
    fi
fi
mkdir -p retrieval
stamp="$(date -u +%Y%m%dT%H%M%S_%N)_$$"
list="retrieval/files_$stamp.list"
find . -type f ! -path './retrieval/*' ! -name workflow.lock ! -name input_package.tar -print0 | sort -z > "$list"
manifest="retrieval/SHA256SUMS_$stamp.txt"
xargs -0 sha256sum < "$list" > "$manifest"
archive="retrieval/result_$stamp.tar"
tar --null -T "$list" -cf "$archive"
tar -rf "$archive" "$manifest"
printf 'ARCHIVE_PATH=%s/%s\n' "$root" "$archive"
