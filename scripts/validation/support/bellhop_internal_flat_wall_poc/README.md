# Bellhop 2020 flat internal-wall validation overlay

This directory contains the complete modified `bellhop.f90` and `Step.f90`
files used by the validation-only flat-wall POC. Their unmodified bases are
from `AcousticsToolbox_2020/Bellhop` (`2020_11_4`).

The build script copies the official source into
`results/validation/bellhop_internal_flat_wall_poc/build`, applies these two
overlays, and compiles a separate Windows executable with a portable
MinGW-w64/GFortran toolchain. It never replaces the official Bellhop source or
executable.

The validation executable requires a two-line `.iw2` sidecar:

1. flat wall range in metres;
2. mapped post-wall receiver range in metres, strictly greater than the wall.

The historical flat regression uses `100` and `103` m.  The same
validation-only binary also accepts translated flat walls for the Stage 0D
constant-height sign audit; no official Bellhop executable is modified.

It is deliberately limited to a uniform, lossless 1500 m/s environment, one
source at transverse depth zero, coherent TL, one receiver range at 103 m, and
Cartesian geometric-hat influence with the existing `.sbp` source pattern.
