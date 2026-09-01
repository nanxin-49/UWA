# Bellhop 2020 curved-wall beam/frame audit overlay

This directory contains the validation-only Bellhop 2020 source overlays and
isolated Windows build script used by the curved-wall covariance audit. The
build copies the official 2020 source tree and produces a separate audit
executable; it does not replace the formal Bellhop executable.

`bellhop.f90` and `Step.f90` retain the sinusoidal internal-wall geometry,
native `Reflect2D` call, and proper post-reflection pi rotation. The
`influence.f90` overlay adds diagnostic logging around the existing
`InfluenceGeoHatCart` calculation; its field-contribution formula is not
changed.

The audit sidecars record reflection state, post-rotation frame state, and
per-ray receiver contributions. The permanent MATLAB regression
`validate_bellhop_shd_receiver_range_pairing_vertical` separately enforces
the 103 m rotated / 97 m native receiver-column pairing and preserves 102 m as
a negative control. This overlay does not implement a PM random wall.
