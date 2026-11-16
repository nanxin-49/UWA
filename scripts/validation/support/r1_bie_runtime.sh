#!/usr/bin/env bash
# Configuration only; never starts MATLAB. All byte counts must be exactly
# representable in MATLAB double and safe for Bash signed arithmetic.
r1_uint() {
    local value=${1-}
    [[ "$value" =~ ^[0-9]+$ ]] || { echo 'Empty/non-integer memory value' >&2; return 1; }
    while [[ ${#value} -gt 1 && "$value" == 0* ]]; do value=${value:1}; done
    [[ ${#value} -lt 16 || ( ${#value} -eq 16 && "$value" < 9007199254740992 ) ]] || {
        echo 'Memory value exceeds exact integer range' >&2; return 1;
    }
    printf '%s\n' "$value"
}
r1_read_available() {
    local file=${1:-/proc/meminfo} key value unit rest found='' count=0 kb
    [[ -r "$file" ]] || { echo 'Cannot read system memory file' >&2; return 1; }
    while read -r key value unit rest; do
        if [[ "$key" == MemAvailable: ]]; then
            [[ "$unit" == kB && -z "$rest" ]] || { echo 'Invalid MemAvailable units/record' >&2; return 1; }
            found=$value; count=$((count+1))
        fi
    done < "$file"
    [[ $count -eq 1 ]] || { echo 'Missing/duplicate MemAvailable record' >&2; return 1; }
    kb=$(r1_uint "$found") || return 1
    ((kb<=8796093022207)) || { echo 'MemAvailable overflows exact byte range' >&2; return 1; }
    printf '%s\n' "$((kb*1024))"
}
r1_prepare_runtime() {
    local available initial baseline
    available=$(r1_uint "${1-}") || return 1
    [[ -r "${2-}" && -f "${2-}" ]] || { echo 'Initial memory file missing/unreadable' >&2; return 1; }
    initial=$(cat -- "$2") || return 1
    initial=$(r1_uint "$initial") || return 1
    # Parenthesized shell arithmetic, never AWK printf/redirection syntax.
    if ((available>initial)); then baseline=$available; else baseline=$initial; fi
    command -v python3 >/dev/null || { echo 'python3 required for JSON serialization' >&2; return 1; }
    python3 "$4" write "$3" "$baseline" "$available" "$initial" || return 1
    python3 "$4" verify "$3" || return 1
}
