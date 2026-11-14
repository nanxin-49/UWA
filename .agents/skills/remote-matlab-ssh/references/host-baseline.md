# Prior verified host facts (not a guarantee of current readiness)

- As of 2026-11-09, approved Windows SSH alias `matlab-ubuntu-new` reached a shared Ubuntu account `misc`, with MATLAB R2025b Update 1 at `/home/misc/MATLAB/R2025b/bin/matlab`. The machine and license must be checked again for each run.
- Historical root `/home/misc/cwx` was not a full Git checkout. Do not replace or `git init` its shared contents. Prior baseline lived in a unique snapshot under `/home/misc/cwx/project_baselines/20261110_121507_new_ubuntu_minimal/project/`.
- A small 23-file cross-platform project baseline passed limited CPU double checks on 2026-11-10; it did not qualify full project dependencies, arbitrary Linux MEX, large windows, or broadband work. See `reports/new_ubuntu_minimal_project_baseline_2026-11-10.md`.
- The prior piston P0/P1 remote run `20261112_remote_p0p1_01` used a unique path under `/home/misc/cwx/piston_p0p1_20261112_01`, up to two independent one-thread MATLAB processes **because that specific scope was authorized**. It does not grant future two-process approval.
- A prior new MATLAB process stopped at MathWorks account/email authentication despite other processes having succeeded. Treat this as an unresolved per-process license preflight when it recurs; do not enter email/credentials, assume active compute, or restart duplicate jobs.
- The old `matlab-ubuntu` host and R2023b results are historical; preserve their provenance and never substitute them for the approved new host.
