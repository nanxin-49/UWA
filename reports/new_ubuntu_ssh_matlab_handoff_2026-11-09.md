# New Ubuntu SSH/MATLAB handoff (2026-11-09)

## Current target and historical boundary

- All future **remote SSH** MATLAB work uses the Windows alias `matlab-ubuntu-new`, login `misc`, and executable `/home/misc/MATLAB/R2025b/bin/matlab`. The designated project directory is `/home/misc/cwx`. Use the alias rather than recording its IP in project files.
- The old alias `matlab-ubuntu` and its R2023b machine are historical. Keep the old Windows SSH entry, key, and dated results intact, but do not use that machine for new work without explicit user reauthorization. Historical reports retain their original machine and release attribution.
- The new alias uses its own encrypted Ed25519 key at Windows `~/.ssh/id_ed25519_matlab_ubuntu_new`. The user verified the new host's ED25519 fingerprint on the Ubuntu console before it was trusted. Windows `ssh-agent` was set to `Manual`, started, and loaded with the new key for the verified session. This proves non-interactive access in that session, not persistence across a reboot or agent restart. Never put private keys, passphrases, or passwords in the repository or run logs.

## What was actually verified

| Check | Observation |
| --- | --- |
| New SSH alias | `ssh -o BatchMode=yes matlab-ubuntu-new 'whoami && hostname'` returned `misc` and hostname `misc`; a terminal allocation also succeeded. |
| Old SSH alias | `ssh -o BatchMode=yes matlab-ubuntu 'whoami && hostname'` returned `misc` and hostname `misc-MS-7E27` after the new alias was configured. Its configuration was not changed. |
| New MATLAB release | `/home/misc/MATLAB/R2025b/bin/matlab -batch 'disp(version)'` returned `25.2.0.3042426 (R2025b) Update 1`, exit code 0. |
| Standalone calculation | A deterministic `fft`/`ifft` round trip with an error assertion passed, reconstruction error 0, exit code 0, elapsed 33.167 s. The full receipt is retained locally at `results/validation/new_ubuntu_matlab_smoke/20261109_232540_new_ubuntu/checkpoint.json` and is excluded from Git. |
| Project directory | `/home/misc/cwx` exists on the new machine, but `git -C /home/misc/cwx rev-parse --show-toplevel` reported that it is not a Git repository. No project code was migrated or run there. |

The standalone calculation checks MATLAB startup and basic numeric execution only. It does not qualify project code, toolboxes, Linux dependencies, model numerics, or earlier R2023b reference results on the new host.

## Procedure for the next authorized remote MATLAB task

1. Read the current-state index in `PROJECT_CONTEXT.md`, the applicable active plan and stop/resource gates, and the required entrypoint. Keep the old alias out of new execution commands.
2. Check `ssh -o BatchMode=yes matlab-ubuntu-new 'whoami && hostname'`, the new MATLAB path and release, and the Windows `ssh-agent` state. If the agent or key is unavailable, restore access through a user-authorized local method; do not switch to the old host or change keys silently.
3. Before any project sync or project execution, verify the canonical `/home/misc/cwx` path, its contents and ownership within the authorized scope, intended code revision, and required inputs. It was not a Git checkout on 2026-11-09. Stop and ask before initializing, replacing, or overwriting unknown remote work; transfer only task-required files after that gate is resolved, excluding `cash/`, credentials, caches, and unrelated results.
4. Check Linux/R2025b syntax and toolbox compatibility, MEX/external binaries, case-sensitive paths, and output formats. Start with a small deterministic **real-project** case and a fixed-seed, explicit-tolerance comparison against an accepted baseline. The earlier standalone test cannot replace this gate.
5. For a substantial run, inspect shared RAM, CPU/load, disk, and task-relevant MATLAB jobs; honor the active plan's budget and stop rules. Use a unique result directory, explicit MATLAB `-batch` executable and input/output paths, full stdout/stderr and real exit status, SHA-256 evidence, bounded receipt, and checkpoint. Keep short tests foreground; use a detached task-specific launcher only for an approved long run.
6. On timeout, disconnect, failed assertion, missing artifact, or compatibility uncertainty, inspect the existing job/checkpoint before retrying. Verify numerical assertions and retained artifact hashes; never infer success from SSH or MATLAB exit code alone.

This handoff changes the remote execution target and workflow only. No physics, public MATLAB interface, or communication method changed.
