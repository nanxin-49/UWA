"""Bounded Linux-only K1 infrastructure diagnostic; no MATLAB or acoustics.

Only launches this file's tiny, finite diagnostic children. Signal tests use
the frozen controller's identity-checked stop_owned; never inspect job args.
The lock is inside the unique diagnostic directory, not the live K1 lease.
"""
import datetime
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time


def put(path, value):
    Path(path).write_text(json.dumps(value, indent=2), encoding="utf-8")


def ticks(pid):
    try:
        return Path(f"/proc/{pid}/stat").read_text().rsplit(") ", 1)[1].split()[19]
    except FileNotFoundError:
        return None


def worker(mode, root):
    root = Path(root)
    if mode == "--lock-worker":
        import fcntl
        assert Path(os.readlink("/proc/self/fd/9")) == root / "diagnostic.lock"
        fcntl.flock(9, fcntl.LOCK_EX | fcntl.LOCK_NB)
    if mode == "--ignore-term-worker":
        signal.signal(signal.SIGTERM, signal.SIG_IGN)
    put(root / "ready.json", dict(pid=os.getpid(), start_ticks=ticks(os.getpid())))
    end = time.monotonic() + 24
    while time.monotonic() < end and not (root / "release").exists():
        time.sleep(.05)


def main(root):
    root = Path(root).resolve()
    assert str(root).startswith("/home/misc/cwx/k1_native_probe_") and root.is_dir()
    sys.path.insert(0, str(root / "code"))
    import k1_pm_supervisor as sup
    import k1_pm_preflight as pre
    started = time.monotonic()
    tests, details, owned = [], [], []
    receipt = dict(status="K1_NATIVE_DIAGNOSTIC_FAILED", failed_checks=1,
                   new_acoustic_solves=0, MATLAB_processes_started=0,
                   shared_K1_lease_accessed=False, other_processes_signalled=0)

    def ok(name):
        tests.append(name)

    def spawn(label, mode, lock=False):
        folder = root / label
        folder.mkdir(exist_ok=False)
        script = root / "code/test_k1_pm_native_linux.py"
        if lock:
            command = ["bash", "-c", 'exec 9>"$1/diagnostic.lock"; flock -n 9 || exit 75; exec nohup python3 -B "$2" --lock-worker "$1"',
                       "k1-native-lock", str(folder), str(script)]
        else:
            command = [sys.executable, "-B", str(script), mode, str(folder)]
        out = (folder / "stdout.log").open("wb")
        err = (folder / "stderr.log").open("wb")
        try:
            child = subprocess.Popen(command, stdout=out, stderr=err, start_new_session=True,
                                     env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1"))
        finally:
            out.close()
            err.close()
        identity = dict(pid=child.pid, process_group=child.pid, start_ticks=ticks(child.pid), root=str(root))
        owned.append((child, identity))
        put(folder / "identity.json", identity)
        until = time.monotonic() + 5
        while not (folder / "ready.json").exists() and child.poll() is None and time.monotonic() < until:
            time.sleep(.05)
        assert (folder / "ready.json").is_file() and child.poll() is None, label + " failed readiness"
        assert pre.read(folder / "ready.json")["pid"] == child.pid
        assert sup.owned_alive(identity)
        return folder, child, identity

    try:
        assert sys.platform == "linux"
        expected = pre.read(root / "source_hashes.json")
        for item in expected:
            assert pre.sha(root / item["path"]) == item["sha256"], item["path"]
        ok("diagnostic_frozen_code_SHA")
        mem = sup.available()
        mp, ip = sup.pressure("memory"), sup.pressure("io")
        cfg = pre.read(root / "config.json")
        cfg["_root"] = str(root)
        try:
            admission = dict(status="PASS", measurements=sup.admission(cfg, mem))
        except AssertionError as exc:
            admission = dict(status="REFUSED_NO_DISPATCH", reason=str(exc))
        details.append(dict(test="actual_admission_read_only", result=admission,
                            MemAvailable_bytes=mem, memory_full_avg10=mp, io_full_avg10=ip))
        ok("native_proc_PSI_statvfs_admission_read_only")

        folder, child, identity = spawn("lock", "--lock-worker", lock=True)
        contended = subprocess.run(["flock", "-n", str(folder / "diagnostic.lock"), "true"], timeout=5)
        assert contended.returncode == 1
        ok("fd9_flock_survives_bash_nohup_exec")
        ok("independent_native_flock_contender_refused")
        (folder / "release").touch()
        assert child.wait(timeout=5) == 0
        released = subprocess.run(["flock", "-n", str(folder / "diagnostic.lock"), "true"], timeout=5)
        assert released.returncode == 0
        ok("native_lock_released_after_owner_exit")

        folder, child, identity = spawn("term", "--term-worker")
        assert sup.group_rss(child.pid) > 0
        ok("native_owned_session_RSS_sample")
        before = time.monotonic()
        sup.stop_owned(child, identity)
        assert child.returncode == -signal.SIGTERM
        details.append(dict(test="real_verified_TERM", elapsed_s=time.monotonic()-before, exit_code=child.returncode))
        ok("real_TERM_only_verified_diagnostic_session")

        folder, child, identity = spawn("identity_refusal", "--term-worker")
        bad = dict(identity, start_ticks=str(int(identity["start_ticks"]) + 1))
        try:
            sup.stop_owned(child, bad)
        except RuntimeError as exc:
            assert "identity mismatch" in str(exc)
        else:
            raise AssertionError("Incorrect identity was not refused")
        assert child.poll() is None
        ok("wrong_start_ticks_refused_live_child_unaffected")
        sup.stop_owned(child, identity)
        assert child.returncode == -signal.SIGTERM

        folder, child, identity = spawn("kill", "--ignore-term-worker")
        before = time.monotonic()
        sup.stop_owned(child, identity)
        elapsed = time.monotonic() - before
        assert child.returncode == -signal.SIGKILL and elapsed >= 14
        details.append(dict(test="real_verified_KILL_after_TERM_timeout", elapsed_s=elapsed, exit_code=child.returncode))
        ok("real_KILL_after_15s_TERM_timeout_identity_rechecked")

        # CLI must reject a missing live K1 lease before verify()/MATLAB dispatch.
        isolated = root / "lease_refusal"
        session = isolated / "sessions/diagnostic"
        session.mkdir(parents=True)
        cli = subprocess.run([sys.executable, "-B", str(root / "code/k1_pm_supervisor.py"), str(isolated), str(session)],
                             capture_output=True, text=True, timeout=5)
        (isolated / "stdout.log").write_text(cli.stdout, encoding="utf-8")
        (isolated / "stderr.log").write_text(cli.stderr, encoding="utf-8")
        value = pre.read(session / "exit.json")
        assert cli.returncode == 78 and value["case_worker_dispatches"] == 0
        assert value["status"] == "K1_CONTROLLER_START_REFUSED"
        ok("actual_controller_CLI_no_lease_refused_before_dispatch")
        assert all(c.poll() is not None for c, _ in owned)
        ok("all_diagnostic_children_reaped_no_surviving_jobs")
        receipt.update(status="K1_NATIVE_DIAGNOSTIC_PASS", failed_checks=0)
    except Exception as exc:
        receipt["error"] = repr(exc)
    finally:
        for child, identity in owned:
            if child.poll() is None:
                try:
                    sup.stop_owned(child, identity)
                except Exception as exc:
                    receipt.setdefault("cleanup_errors", []).append(repr(exc))
                    receipt.update(status="K1_NATIVE_DIAGNOSTIC_FAILED", failed_checks=1)
        receipt.update(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), elapsed_s=time.monotonic()-started,
                       tests=len(tests), test_names=tests, diagnostic_child_processes=len(owned),
                       all_diagnostic_children_reaped=all(c.poll() is not None for c, _ in owned), details=details,
                       limitations=["Uses isolated diagnostic lease, not shared live K1 lease or a formal launcher Start",
                                    "No additional parallel MATLAB license seats tested", "Formal K1 acoustics unrun"])
        put(root / "receipt.json", receipt)
        files = [p for p in root.rglob("*") if p.is_file() and p.name != "result_sha256.json"]
        put(root / "result_sha256.json", [dict(path=p.relative_to(root).as_posix(), sha256=pre.sha(p)) for p in sorted(files)])
    print(json.dumps({k: v for k, v in receipt.items() if k not in ["details", "test_names"]}))
    return 0 if receipt["failed_checks"] == 0 else 1


if __name__ == "__main__":
    if sys.argv[1].startswith("--"):
        worker(sys.argv[1], sys.argv[2])
    else:
        sys.exit(main(sys.argv[1]))
