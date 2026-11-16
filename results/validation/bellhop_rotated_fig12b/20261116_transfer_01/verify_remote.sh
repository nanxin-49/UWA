set -eu
root=/home/misc/cwx/bellhop_rotated_C_transfer_20261116_01
meta="$root/transfer_metadata"
test ! -e "$meta/checkpoint_verify.json"
start=$(date +%s)
set +e
(cd "$root/Explain" && LC_ALL=C sha256sum -c ../transfer_metadata/SHA256SUMS.lf) > "$meta/sha256_check.log" 2> "$meta/sha256_check.stderr.log"
code=$?
set -e
end=$(date +%s)
count=$(grep -c ': OK$' "$meta/sha256_check.log" || true)
expected=$(wc -l < "$meta/SHA256SUMS.lf")
if test "$code" -eq 0 && test "$count" -eq "$expected"; then status=TRANSFER_SHA_PASS; else status=STOPPED_SHA_FAILED; fi
printf '{"status":"%s","verified_files":%s,"expected_files":%s,"sha256_exit_code":%s,"elapsed_s":%s,"new_Bellhop_processes":0,"C_completed":30,"C_recomputed":0,"remote_project_root":"%s/Explain"}\n' "$status" "$count" "$expected" "$code" "$((end-start))" "$root" > "$meta/checkpoint_verify.json"
cat "$meta/checkpoint_verify.json"
if test "$status" != TRANSFER_SHA_PASS; then tail -n 12 "$meta/sha256_check.stderr.log"; exit 1; fi
for rel in src/channel src/propagation src/surface src/receiver src/statistics src/communication src/bubble examples scripts/comparisons scripts/experiments; do mkdir -p "$root/Explain/$rel"; done
{
 date -u +%FT%TZ
 gfortran --version | head -n 1
 make --version | head -n 1
 ar --version | head -n 1
 awk '/MemAvailable/ {print}' /proc/meminfo
 df -Pk "$root"
} > "$meta/environment_after_transfer.log"
