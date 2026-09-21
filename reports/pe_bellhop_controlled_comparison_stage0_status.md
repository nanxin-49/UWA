# PE--Bellhop controlled comparison execution status

日期：2026-09-10  
执行方式：MATLAB MCP；Bellhop 2020 official executable plus rebuilt
validation-only flat/parametric binaries。

## Stage 0

**PASS**. Authoritative artifacts are under
`results/validation/pe_bellhop_controlled_comparison/stage0/`.

| check/result | value |
|---|---:|
| PE--exact-AS M99 complex L2 | `4.2785525e-13` |
| flat internal-wall M99 complex L2 | `0.0058993922` |
| flat internal-wall M99 phase RMS | `0.0058992479 rad` |
| 5001→10001 beam M99 L2 | `5.854857e-06` |
| generic-flat ↔ official free-field M99 L2 | `0` |
| half-dx shared-node M99 L2 | PASS |
| wall residual / kappa / phase / positive range | PASS |
| NaN/Inf | `0` |

The Stage-0 floor is frozen for subsequent stages; no PE or Bellhop core was
modified.

## Stage 1 status

The requested weak sinusoid (`A=0.01 m`, `K=0.10 rad/m`) was started, but the
first high-cost generic internal-wall run did not finish in the bounded
execution window and was stopped. No physical PASS/FAIL conclusion is assigned
to the incomplete primary run, and Stage 2--7 remain locked.

The partial run reached valid input/diagnostic output for the lower-cost smoke
configuration (`1001→2001` beams); its model comparison is not authoritative
because the primary beam gate was not used. The incomplete primary status is
recorded in `results/validation/pe_bellhop_controlled_comparison/stage1/`.

This is a runtime-completion limitation of the current generic validation
binary, not evidence for a PE--Bellhop model discrepancy. A future continuation
should resume Stage 1 with a bounded, separately monitored Bellhop invocation or
an explicitly approved cost reduction, while preserving the frozen Goal gates.
