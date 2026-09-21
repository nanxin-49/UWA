# Bellhop 2020 tilted straight internal-wall visual-enhancement overlay

This directory contains the complete modified `bellhop.f90` and `Step.f90`
files used by the validation-only tilted-wall visual-enhancement experiment. Their unmodified bases are
from `AcousticsToolbox_2020/Bellhop` (`2020_11_4`).

The build script copies the official source into
`results/validation/bellhop_internal_tilted_wall_poc/build`, applies these two
overlays, and compiles a separate Windows executable with a portable
MinGW-w64/GFortran toolchain. It never replaces the official Bellhop source or
executable.

The validation executable requires a three-line `.iw2` sidecar:

1. line-intercept `R0` in metres (must be `100`);
2. straight-wall slope `a` (finite and `|a|<=1`), for `r=R0+a*z`;
3. mapped receiver range in metres (must be `103`).

The wall frame is `t=(a,1)/sqrt(1+a^2)`, `n=(1,-a)/sqrt(1+a^2)` with
`kappa=0`, and the physical reflection calls the native 2-D `Reflect2D` with
TOP semantics and pressure-release conditions. After that one reflection the
overlay applies the already validated proper pi rotation and passes only the
transformed, monotonically increasing branch to the unchanged Cartesian
geometric-hat influence routine.

It is deliberately limited to a uniform, lossless 1500 m/s environment, one
source at transverse depth zero, coherent TL, one wall reflection, and the
existing `.sbp` Gaussian source pattern.
