#!/bin/bash
set -euo pipefail
ROOT=/home/misc/cwx/bellhop_rotated_C_transfer_20261116_01/Explain
OUT=$ROOT/results/validation/bellhop_rotated_fig12b/20261116_ray_remote_03
SRC=$ROOT/results/validation/bellhop_rotated_fig12b/20261116_linux_C_01/build
start=$(date +%s)
test ! -e "$OUT/build"
cp -r "$SRC" "$OUT/build"
mkdir -p "$OUT/bin"
python3 - "$OUT/build/Bellhop/bellhop.f90" <<'PY'
import sys
p=sys.argv[1]
s=open(p).read()
a='  ray2D( 1 : NPost ) = ray2D( WallBranchStart : Beam%Nsteps )\n  Beam%Nsteps = NPost'
b="  IF ( Beam%RunType(1:1) /= 'R' ) THEN\n"+a+"\n  END IF"
assert s.count(a)==1
open(p,'w').write(s.replace(a,b))
PY
diff -u "$SRC/Bellhop/bellhop.f90" "$OUT/build/Bellhop/bellhop.f90" > "$OUT/R_only_export.patch" || test $? -eq 1
cd "$OUT/build/Bellhop"
flags=(-O2 -std=legacy -fno-range-check -ffree-line-length-none -static -I../misc -I../tslib)
gfortran -c "${flags[@]}" bellhop.f90
gfortran -o "$OUT/bin/bellhop_ray_export_2020" "${flags[@]}" angleMod.o ArrMod.o bdryMod.o bellhopMod.o sspMod.o RotatedFieldState.o ReflectMod.o WriteRay.o influence.o Step.o ReadEnvironmentBell.o bellhop.o ../misc/libmisc.a
sha256sum "$OUT/bin/bellhop_ray_export_2020" > "$OUT/binary_SHA256SUMS"
find "$OUT/build" -name '*.f90' -print0 | sort -z | xargs -0 sha256sum > "$OUT/source_SHA256SUMS"
echo $(($(date +%s)-start)) > "$OUT/build_elapsed_s.txt"
