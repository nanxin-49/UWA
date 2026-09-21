# Stage 7 existing M=50 interpretation

状态：**PASS**；solver calls/seeds rerun: `0/0`。

**Provenance: historical Bellhop point-source R.** These values are not renamed X and receive no empirical R-to-X correction.

| M | seed range | mean Delta TL (dB) | std (dB) | PE/BH mean power | mean power difference | circular mean phase (rad) | circular std | R |
|---:|:---:|---:|---:|---:|---:|---:|---:|---:|
| 50 | 260001--260050 | -0.0074785801 | 0.962232465 | 1.00450991 / 1.00700722 | -0.0024973135 | -0.375344713 | 0.594700972 | 0.837918346 |

- Convention: Stage 1Y closes the active comparison convention; M=50 remains historical R and is not relabeled or corrected.
- Power: Near-equal ensemble mean powers are compatible with signed realization-to-realization amplitude differences averaging out.
- Phase: Nonzero circular mean phase is not removed by the closed convention and is consistent with spatial reflection-model discrepancy.
- Validity: The canonical fixed-PM case lies in Region III, beyond the sampled smooth-sinusoid Region-I/II transition.
- Source limitation: No empirical R-to-X correction is applied; M=50 cannot replace the new X-source controlled result.

The M=50 ensemble therefore supports a model-discrepancy interpretation with historical-source provenance limits: near-equal mean power does not imply per-realization field closure, and the nonzero circular phase is consistent with the spatial distortion identified in Stage 5 and the Region-III fixed-PM result in Stage 6.
