---
name: remote-matlab-ssh
description: >-
  Execute or manage explicitly authorized MATLAB jobs on the laboratory Ubuntu
  host from Windows Codex over SSH, including connection checks, snapshots,
  batch execution, approved parallel jobs, recovery, and verified retrieval.
  Use only for explicit remote SSH/MATLAB operations; do not trigger for local
  MATLAB work, theory discussion, report reading, or mere mention of SSH.
---

# Remote MATLAB via SSH: Windows Codex -> shared Ubuntu

## Scope and execution choice

- Apply this workflow only for a task explicitly authorized to use SSH or the remote MATLAB host. The existing local MATLAB MCP session remains the default for lightweight local code checks. A failed MCP connection is **not** approval to switch to SSH or to launch another local MATLAB instance.
- This skill is an execution procedure, **not authorization** to run a solve, copy unrelated files, advance a research stage, spend additional compute, or run multiple workers. Use the user's most recent specific scope, time/RAM ceiling, and stop gates. Read `AGENTS.md` for repository-level physics and artifact conventions.

## Pinned target and identity

- Use Windows SSH alias `matlab-ubuntu-new` with non-interactive key authentication: `ssh -o BatchMode=yes matlab-ubuntu-new ...`. The login on this host is shared account `misc`; approved working root `/home/misc/cwx`.
- Invoke `/home/misc/MATLAB/R2025b/bin/matlab` explicitly, normally with `-batch` and `-singleCompThread`. Verify the actual release and license before relying on a prior run's configuration.
- Never use historical alias `matlab-ubuntu` (old R2023b host) as fallback without new user authorization. Do not hard-code an IP address, passwords, account-email credentials, private keys or tokens in source, reports, logs or Git.
- Host location and key material come from the user's local SSH configuration; never attempt to discover credentials or modify SSH agent/server settings on the user's behalf.

## Safety on a shared machine

- Treat `/home/misc/cwx` as a shared-account **convention, not a security boundary**. Work only in the exact user-approved project subtree and task-owned unique run directories; do not enumerate or modify other users' directories, session data or jobs.
- Do not change `/home/misc/.ssh`, `/home/misc/.codex`, global dotfiles, the MATLAB installation, licenses, firewall, or the remote SSH server. No `sudo`, unattended package installs, destructive sync, shared-root cleanup, or implicit remote Codex login.
- Do not kill by process name (`pkill matlab`, `killall matlab`) or stop unidentified PIDs. Only terminate a positively identified task-owned job after specific approval.

## Before any run

1. Confirm the user actually approved the remote task and whether it is inspection, small verification, heavy solve, parallel work, transfer, or result retrieval. Limit commands to that scope; if approvals are absent, perform static/local planning only.
2. Test non-interactive SSH (`BatchMode=yes`), confirm the canonical root exists and inspect only the intended run subtree. Do **not** assume `/home/misc/cwx` is a complete Git checkout or overwrite it.
3. Verify the explicit MATLAB executable, license for the intended process count, required dependencies/platform behavior (Linux MEX, external programs, toolboxes), and task inputs. Record code revision plus any uncommitted local changes.
4. Check aggregate available RAM, CPU/load and free disk; measure only task-owned MATLAB jobs rather than inspecting unrelated user activity. Reserve >=30% of starting available RAM unless the user approved a stricter limit. Estimate memory and runtime of each planned job. Heavy, parallel, or multi-process work needs explicit budget/worker approval, not historical precedent.

## Transfer and execution

1. Build a **non-destructive allowlisted code/data snapshot** for the authorized run. Exclude repository `cash/`, credentials, other experiments, local caches and unrelated outputs. Preserve local edits. Send/check SHA-256 manifests and create a unique remote run ID; no `rsync --delete` or overwrite of preexisting remote files.
2. Resolve project root explicitly; call `setup_vertical_project` where appropriate. Use absolute output paths and never depend on interactive shell state, implicit MATLAB path, or `cd` inside production code.
3. For short probes, foreground `ssh ... matlab -batch ...` is acceptable. For approved longer jobs, use a task-owned detached launcher with durable stdout/stderr, MATLAB diary, PID, process identity/start time, command, receipt and real exit code. Use `-singleCompThread` where appropriate.
4. If the user authorized parallel work, schedule **bounded independent batch processes**, each with its own job ID and outputs; cap concurrency per explicit approval and actual peak-RAM estimates. Do not create a MATLAB `parpool`, modify thread counts, or silently increase workers.
5. For transient failures/SSH disconnects, first inspect **the original task-owned** job's PID, stdout/stderr, receipt and checkpoint. A MATLAB license login/email prompt is **not** evidence that the calculation started. Do not duplicate an existing job or fill credentials; request user action if license authentication blocks non-interactive startup.

## Closeout

- Verify numerical assertions and convergence/physical gates; check real exit code, full logs, expected MAT/CSV/figures, input/code/output SHA and remote-to-local transfer integrity. `exit code 0` alone is not numerical PASS.
- Store code snapshot, config, provenance, logs, receipts, checksums, peak memory/elapsed time and audit under an isolated `results/<category>/<case>/<run_id>/` path; keep local project reports in `reports/` as needed. Do not overwrite prior runs or reinterpret old-machine/old-source results.
- Return a compact receipt: task/run ID, server/MATLAB release, code SHA, configuration and source width, stage, status, key metrics/failed gates, RAM/time, local result/checkpoint links, and explicit limitations. Update `reports/active_task_tracker.md` and (if project status materially changed) project progress records under `AGENTS.md`.
- Do not push/merge/commit Git unless the user explicitly asks. On safety, access, license, RAM, or scientific-reference failure, stop the affected job and preserve evidence instead of widening the task.

## Historical reference

For the new-host baseline and pitfalls that are relevant only when planning compatibility or recovery, read `references/host-baseline.md`; do not use prior smoke-test results as automatic qualification for a new full-workload run.
