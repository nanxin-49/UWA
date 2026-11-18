#!/usr/bin/env bash
# Local synthetic RAM/lock tests only. No SSH, MATLAB or process termination.
set -eu
repo=$(realpath "${1:?repository}"); out=${2:?new output directory}
[[ ! -e "$out" ]] || exit 2
mkdir -p "$out"; out=$(realpath "$out")
code="$repo/scripts/validation/support"; if [[ -d "$repo/code" ]]; then code="$repo/code"; fi
source "$code/r1_bie_runtime.sh"
source "$code/r2_pm_parallel.sh"
passed=0
for tuple in '64 64 yes' '40 40 yes' '32 32 no' '30 40 no' '40 64 no' '48 64 yes'; do
    read -r current baseline expected <<< "$tuple"
    if r2_parallel_admit "$((current*1024*1024*1024))" "$((baseline*1024*1024*1024))"; then actual=yes; else actual=no; fi
    [[ "$actual" == "$expected" ]]; passed=$((passed+1))
done
if r2_parallel_admit x 1; then exit 3; fi; passed=$((passed+1))
if r2_parallel_admit 2 1; then exit 3; fi; passed=$((passed+1))
base=$((40*1024*1024*1024)); boundary=$((base*3/10+24*1024*1024*1024))
r2_parallel_admit "$boundary" "$base"; passed=$((passed+1))
if r2_parallel_admit "$((boundary-1))" "$base"; then exit 3; fi; passed=$((passed+1))
# Two independent per-seed slots coexist; a third copy of either seed fails.
if ! command -v flock >/dev/null; then
    printf '{"status":"PASS_AVAILABLE_CHECKS","tests":%s,"failed_checks":0,"new_acoustic_solves":0,"SSH_connections":0,"native_flock_checks":"UNRUN: flock unavailable locally; mandatory before remote MATLAB"}\n' "$passed" > "$out/receipt.json"
    cat "$out/receipt.json"; exit 0
fi
exec 8> "$out/seed260002.lock"; flock -n 8
exec 9> "$out/seed260003.lock"; flock -n 9; passed=$((passed+1))
if (exec 8>&-; exec 9>&-; flock -n "$out/seed260002.lock" true); then exit 3; fi; passed=$((passed+1))
if (exec 8>&-; exec 9>&-; flock -n "$out/seed260003.lock" true); then exit 3; fi; passed=$((passed+1))
exec 8>&-; exec 9>&-
flock -n "$out/seed260002.lock" true; passed=$((passed+1))
printf '{"status":"PASS","tests":%s,"failed_checks":0,"new_acoustic_solves":0,"SSH_connections":0}\n' "$passed" > "$out/receipt.json"
cat "$out/receipt.json"
