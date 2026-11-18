"""One owned R2025b MATLAB process, no retries, no network, RAM/time stops."""
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time
from d0_package import preflight, sha, verify


def available():
    for line in Path('/proc/meminfo').read_text().splitlines():
        if line.startswith('MemAvailable:'):
            return int(line.split()[1]) * 1024
    raise RuntimeError('MemAvailable missing')


def write(path, value):
    temp = path.with_suffix('.pending')
    temp.write_text(json.dumps(value, indent=2), encoding='utf-8')
    temp.replace(path)


def stop_reason(now, baseline, elapsed, depth_elapsed, cfg):
    if now < baseline * cfg['RAM_reserve_fraction']:
        return 'RAM_RESERVE_STOP'
    if elapsed > cfg['max_wall_s']:
        return 'WALL_TIME_STOP'
    if depth_elapsed > cfg['max_depth_wall_s']:
        return 'DEPTH_WALL_TIME_STOP'
    return None


def validate_output(root, binding):
    output = root / 'output'
    cp = json.loads((output / 'checkpoint.json').read_text())
    assert cp['status'] == 'D0_COMPLETE_INCIDENT_ONLY' and cp['manifest_sha256'] == binding
    assert cp['failed_or_unrun_checks'] == []
    verify(output, cp['artifacts_sha256'])
    for depth in [100, 300, 600, 1000, 1490]:
        folder = output / f'depth_{depth:04d}m'
        depth_cp = json.loads((folder / 'checkpoint.json').read_text())
        assert depth_cp['status'] == 'AUDITED_PASS' and depth_cp['manifest_sha256'] == binding
        verify(folder, depth_cp['artifacts_sha256'])
        audit = json.loads((folder / 'audit.json').read_text())
        assert audit['failed_checks'] == 0 and all(c['passed'] for c in audit['checks'])
    assert json.loads((output / 'figures/figure_audit.json').read_text())['failed_checks'] == 0


def main(root, session):
    root = Path(root).resolve()
    session = Path(session).resolve()
    assert root / 'sessions' in session.parents and session.is_dir()
    started = time.monotonic()
    child = None
    result = dict(status='FAILED_BEFORE_MATLAB', matlab_exit_code=None, supervisor_exit_code=78)
    interrupted = [False]
    signal.signal(signal.SIGINT, lambda *_: interrupted.__setitem__(0, True))
    signal.signal(signal.SIGTERM, lambda *_: interrupted.__setitem__(0, True))
    try:
        cfg, _ = preflight(root)
        binding = sha(root / 'manifest.json')
        baseline = available()
        assert baseline >= cfg['minimum_available_RAM_GiB'] * 1024**3
        assert baseline * .7 >= cfg['estimated_working_RAM_GiB'] * 1024**3
        assert shutil.disk_usage(root).free >= 2*1024**3
        assert not (root / 'output').exists(), 'Existing output requires manual inspection, no retry'
        executable = '/home/misc/MATLAB/R2025b/bin/matlab'
        assert Path(executable).is_file()
        expression = f"addpath('{root}/code'); run_d0_incident_depth('{root}','execute');"
        command = [executable, '-singleCompThread', '-batch', expression]
        env = os.environ.copy()
        env.update(OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
                   D0_SUPERVISED_MANIFEST_SHA256=binding)
        write(session / 'runtime.json', dict(command=command, initial_available_bytes=baseline,
              reserve_fraction=.30, max_wall_s=cfg['max_wall_s'], max_depth_wall_s=cfg['max_depth_wall_s'],
              manifest_sha256=binding, load_average=os.getloadavg(), workers=1))
        minimum = baseline
        peak_rss = 0
        reason = None
        with (session / 'stdout.log').open('wb') as stdout, (session / 'stderr.log').open('wb') as stderr:
            child = subprocess.Popen(command, cwd=root, env=env, stdout=stdout, stderr=stderr, start_new_session=True)
            write(session / 'identity.json', dict(matlab_pid=child.pid, process_group=child.pid,
                  supervisor_pid=os.getpid(), start_monotonic=started, root=str(root)))
            with (session / 'resource.tsv').open('w', encoding='ascii', buffering=1) as log:
                log.write('elapsed_s\tMemAvailable_bytes\towned_MATLAB_PID_RSS_bytes\n')
                while child.poll() is None:
                    now = available()
                    minimum = min(minimum, now)
                    elapsed = time.monotonic() - started
                    rss = 0
                    try:
                        for line in Path(f'/proc/{child.pid}/status').read_text().splitlines():
                            if line.startswith('VmRSS:'):
                                rss = int(line.split()[1])*1024
                    except FileNotFoundError:
                        pass
                    peak_rss = max(peak_rss, rss)
                    active = [p for p in (root/'output').glob('depth_*/started.json')
                              if not (p.parent/'checkpoint.json').exists()]
                    depth_elapsed = max([time.time()-p.stat().st_mtime for p in active] + [0])
                    reason = stop_reason(now, baseline, elapsed, depth_elapsed, cfg)
                    if interrupted[0]:
                        reason = 'OWNED_SESSION_INTERRUPTED'
                    log.write(f'{elapsed:.3f}\t{now}\t{rss}\n')
                    if reason:
                        os.killpg(child.pid, signal.SIGTERM)
                        try:
                            child.wait(timeout=15)
                        except subprocess.TimeoutExpired:
                            os.killpg(child.pid, signal.SIGKILL)
                        break
                    time.sleep(1)
            code = child.wait()
        valid = False
        if code == 0 and reason is None:
            validate_output(root, binding)
            valid = True
        result.update(status='D0_COMPLETE_INCIDENT_ONLY' if valid else 'FAILED_STOPPED',
                      matlab_exit_code=code, supervisor_exit_code=0 if valid else 1,
                      minimum_available_bytes=minimum, owned_PID_peak_RSS_bytes=peak_rss,
                      stop_reason=reason)
    except Exception as exc:
        result['error'] = repr(exc)
        if child is not None:
            if child.poll() is None:
                os.killpg(child.pid, signal.SIGTERM)
                try:
                    child.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    os.killpg(child.pid, signal.SIGKILL)
                    child.wait()
            result['matlab_exit_code'] = child.returncode
    finally:
        result['elapsed_s'] = time.monotonic()-started
        write(session / 'exit.json', result)
    return result['supervisor_exit_code']


if __name__ == '__main__':
    sys.exit(main(*sys.argv[1:3]))
