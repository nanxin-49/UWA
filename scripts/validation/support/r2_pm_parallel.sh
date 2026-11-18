#!/usr/bin/env bash
# R2 aggregate reservation, separate from immutable R1 runtime serializer.
r2_parallel_admit() {
    local available baseline
    available=$(r1_uint "${1-}") || return 2
    baseline=$(r1_uint "${2-}") || return 2
    ((baseline>=available)) || return 2
    ((available-baseline*3/10 >= 24*1024*1024*1024))
}
