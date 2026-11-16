set -eu
run=/home/misc/cwx/bellhop_rotated_C_20261116_remote01
umask 077
mkdir "$run"
{
 date -u +%FT%TZ
 printf 'TARGET_USER '; id -un
 cat /etc/os-release
 uname -m
 getconf _NPROCESSORS_ONLN
 uptime
 printf 'MAKE '; make --version | head -n 1
 printf 'GFORTRAN_PATH\n'
 for cmd in gfortran gfortran-15 gfortran-14 gfortran-13 gfortran-12 gfortran-11; do
  if command -v "$cmd" >/dev/null 2>&1; then command -v "$cmd"; "$cmd" --version | head -n 1; fi
 done
 printf 'INSTALLED_GFORTRAN_PACKAGES\n'
 dpkg-query -W -f='${Package} ${db:Status-Abbrev} ${Version}\n' 'gfortran*' 2>/dev/null || true
 printf 'RAM_KIB\n'; head -n 5 /proc/meminfo
 printf 'DISK_KIB\n'; df -Pk "$run"
 printf 'MATLAB_FILE\n'; test -x /home/misc/MATLAB/R2025b/bin/matlab; readlink -f /home/misc/MATLAB/R2025b/bin/matlab
} > "$run/environment.log" 2> "$run/environment.stderr.log"
cat > "$run/matlab_probe.m" <<'MATLAB'
tic; diary('/home/misc/cwx/bellhop_rotated_C_20261116_remote01/matlab_diary.log');
r=struct('status','MATLAB_HEALTH_PASS','release',version,'architecture',computer,'license_test',license('test','MATLAB'),'pid',feature('getpid'),'new_physics_solves',0);
assert(r.license_test==1); assert(isequal(1+1,2)); r.elapsed_s=toc;
fid=fopen('/home/misc/cwx/bellhop_rotated_C_20261116_remote01/matlab_receipt.json','w'); assert(fid>=0); fprintf(fid,'%s',jsonencode(r,PrettyPrint=true)); fclose(fid); disp(jsonencode(r)); diary off
MATLAB
set +e
/usr/bin/time -v -o "$run/matlab_resources.txt" timeout 60s /home/misc/MATLAB/R2025b/bin/matlab -singleCompThread -batch "run('$run/matlab_probe.m')" > "$run/matlab_stdout.log" 2> "$run/matlab_stderr.log"
code=$?
set -e
printf '%s\n' "$code" > "$run/matlab_exit_code.txt"
(cd "$run" && sha256sum environment.log environment.stderr.log matlab_probe.m matlab_resources.txt matlab_stdout.log matlab_stderr.log matlab_exit_code.txt matlab_diary.log matlab_receipt.json > preflight_sha256.txt)
cat "$run/environment.log"
printf 'MATLAB_EXIT %s\n' "$code"
if test -f "$run/matlab_receipt.json"; then cat "$run/matlab_receipt.json"; fi
exit "$code"
