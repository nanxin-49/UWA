#!/usr/bin/env bash
set -eu
repo=$(realpath "${1:?repository required}")
out=${2:?new output directory required}
[[ ! -e "$out" ]] || exit 2
mkdir -p "$out"; out=$(realpath "$out")
source "$repo/scripts/validation/support/r1_bie_runtime.sh"
# Git Bash on Windows has python.exe; emulate only the executable name in
# this local test. Production still requires an installed python3.
if ! command -v python3 >/dev/null; then python3() { python "$@"; }; fi
serializer="$repo/scripts/validation/support/r1_bie_runtime_json.py"
passed=0
for pair in '20 10 20' '10 10 10' '5 10 10' '0 0 0' '0005 0010 10'; do
    read -r current initial expected <<< "$pair"
    printf '%s\n' "$initial" > "$out/initial.txt"
    json="$out/valid_$passed.json"
    r1_prepare_runtime "$current" "$out/initial.txt" "$json" "$serializer"
    python3 -c 'import json,sys; r=json.load(open(sys.argv[1])); assert r["RAM_baseline_bytes"]==int(sys.argv[2]); assert r["starting_available_bytes"]==int(sys.argv[3]); assert r["initial_available_bytes"]==int(sys.argv[4])' "$json" "$expected" "$current" "$initial"
    passed=$((passed+1))
done
for bad in '' x -1 1.5 1e3 '1 2' 9007199254740992; do
    printf '10\n' > "$out/initial.txt"
    if r1_prepare_runtime "$bad" "$out/initial.txt" "$out/bad_current.json" "$serializer"; then exit 3; fi
    printf '%s\n' "$bad" > "$out/initial.txt"
    if r1_prepare_runtime 10 "$out/initial.txt" "$out/bad_initial.json" "$serializer"; then exit 3; fi
    passed=$((passed+2))
done
if r1_prepare_runtime 10 "$out/missing.txt" "$out/missing.json" "$serializer"; then exit 3; fi
passed=$((passed+1))
for bad in '' 'MemTotal: 10 kB' 'MemAvailable: x kB' 'MemAvailable: 5 bytes' 'MemAvailable: 5 kB extra'; do
    printf '%s\n' "$bad" > "$out/meminfo"
    if r1_read_available "$out/meminfo"; then exit 3; fi
    passed=$((passed+1))
done
if r1_read_available "$out/no_meminfo"; then exit 3; fi
passed=$((passed+1))
printf 'MemAvailable: 12 kB\n' > "$out/meminfo"
[[ $(r1_read_available "$out/meminfo") == 12288 ]]; passed=$((passed+1))
printf '{"RAM_baseline_bytes":}\n' > "$out/invalid.json"
if python3 "$serializer" verify "$out/invalid.json"; then exit 3; fi
passed=$((passed+1))
printf '10\n' > "$out/initial.txt"
# Existing output is immutable, never replace an old session's runtime file.
if r1_prepare_runtime 10 "$out/initial.txt" "$out/valid_0.json" "$serializer"; then exit 3; fi
passed=$((passed+1))
printf 'PASS tests=%s MATLAB_starts=0 SSH_connections=0\n' "$passed"
