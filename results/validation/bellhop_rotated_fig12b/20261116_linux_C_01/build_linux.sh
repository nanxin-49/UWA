#!/bin/bash
set -euo pipefail
ROOT=/home/misc/cwx/bellhop_rotated_C_transfer_20261116_01/Explain
OUT=$ROOT/results/validation/bellhop_rotated_fig12b/20261116_linux_C_01
SRC=$ROOT/results/validation/bellhop_rotated_fig12b/20261115_performance_review_01/optimized_plain_build01/build/AcousticsToolbox_2020
OVER=$ROOT/results/validation/bellhop_rotated_fig12b/20261115_performance_review_01/optimized_plain_overlay01
mkdir -p "$OUT"
test ! -e "$OUT/build"
start=$(date +%s)
mkdir -p "$OUT/build/misc" "$OUT/build/Bellhop" "$OUT/build/tslib" "$OUT/bin"
cp "$SRC"/misc/*.f90 "$OUT/build/misc/"
cp "$SRC"/Bellhop/*.f90 "$OUT/build/Bellhop/"
if compgen -G "$SRC/tslib/*.f90" > /dev/null; then cp "$SRC"/tslib/*.f90 "$OUT/build/tslib/"; fi
for file in SourceReceiverPositions ReadEnvironmentBell Step bellhop influence bellhopMod sspMod RotatedFieldState; do
 if [ "$file" = SourceReceiverPositions ]; then dest=misc; else dest=Bellhop; fi
 cmp "$OVER/$file.f90" "$OUT/build/$dest/$file.f90"
done
find "$OUT/build" -name '*.f90' -type f -print0 | sort -z | xargs -0 sha256sum > "$OUT/source_SHA256SUMS"
gfortran --version > "$OUT/compiler_version.txt"
flags=(-O2 -std=legacy -fno-range-check -ffree-line-length-none -static -I../misc -I../tslib)
cd "$OUT/build/misc"
for file in FatalError MathConstants pchipMod AttenMod PolyMod RefCoef beampattern subtabulate calculateweights monotonicMod SortMod SourceReceiverPositions sspMod RWSHDFile interpolation MergeVectorsMod munk splinec norms cross_products PekRoot RootFinderSecantMod ReadEnvironmentMod; do
 gfortran -c "${flags[@]}" "$file.f90"
done
ar rcs libmisc.a FatalError.o beampattern.o MathConstants.o RefCoef.o SourceReceiverPositions.o pchipMod.o AttenMod.o sspMod.o RWSHDFile.o interpolation.o MergeVectorsMod.o munk.o ReadEnvironmentMod.o SortMod.o splinec.o subtabulate.o calculateweights.o norms.o cross_products.o monotonicMod.o PolyMod.o PekRoot.o RootFinderSecantMod.o
cd "$OUT/build/Bellhop"
for file in bellhopMod angleMod ArrMod bdryMod sspMod RotatedFieldState Cone ReflectMod WriteRay Step influence ReadEnvironmentBell RayNormals bellhop; do
 gfortran -c "${flags[@]}" "$file.f90"
done
gfortran -o "$OUT/bin/bellhop_rotated_receiver64_2020" "${flags[@]}" angleMod.o ArrMod.o bdryMod.o bellhopMod.o sspMod.o RotatedFieldState.o ReflectMod.o WriteRay.o influence.o Step.o ReadEnvironmentBell.o bellhop.o ../misc/libmisc.a
sha256sum "$OUT/bin/bellhop_rotated_receiver64_2020" > "$OUT/binary_SHA256SUMS"
end=$(date +%s)
echo $((end-start)) > "$OUT/build_elapsed_s.txt"
printf 'BUILD_PASS elapsed=%ss\n' "$((end-start))"
