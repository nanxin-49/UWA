"""Single owned MATLAB process; durable exit status, RAM/time guard, no SSH."""
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time
from pm_k0_preflight import verify


def available():
    with open('/proc/meminfo', encoding='ascii') as f:
        for line in f:
            if line.startswith('MemAvailable:'):
                return int(line.split()[1]) * 1024
    raise RuntimeError('MemAvailable missing')


def write(path, value):
    temp = path.with_suffix('.pending')
    temp.write_text(json.dumps(value, indent=2), encoding='utf-8')
    temp.replace(path)


def main():
    root = Path(sys.argv[1]).resolve()
    session = Path(sys.argv[2]).resolve()
    assert root / 'sessions' in session.parents
    started = time.monotonic()
    result = dict(status='FAILED_BEFORE_MATLAB', matlab_exit_code=None, supervisor_exit_code=78)
    child = None
    interrupted = [False]
    signal.signal(signal.SIGTERM, lambda *_: interrupted.__setitem__(0, True))
    signal.signal(signal.SIGINT, lambda *_: interrupted.__setitem__(0, True))
    try:
        manifest = json.loads((root / 'manifest.json').read_text())
        verify(root, manifest['files'])
        cfg = json.loads((root / 'config.json').read_text())
        initial_available = available()
        assert initial_available >= 2 * 1024**3 and shutil.disk_usage(root).free >= 1024**3
        baseline_path = root / 'ram_baseline.json'
        from pm_k0_preflight import sha
        binding = sha(root / 'manifest.json')
        baseline = initial_available
        if baseline_path.is_file():
            prior = json.loads(baseline_path.read_text())
            assert prior['manifest_sha256'] == binding
            baseline = max(baseline, prior['available_bytes'])
        assert initial_available >= baseline * cfg['RAM_reserve_fraction'], 'RAM reserve before launch'
        write(baseline_path, dict(manifest_sha256=binding, available_bytes=baseline))
        matlab = '/home/misc/MATLAB/R2025b/bin/matlab'
        assert Path(matlab).is_file()
        # Root was restricted to safe ASCII in the shell launcher.
        expression = f"addpath('{root}/code'); run_pm_wind_cutoff_k0('{root}','execute');"
        env = os.environ.copy()
        for name in ['OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS']:
            env[name] = '1'
        command = [matlab, '-singleCompThread', '-batch', expression]
        write(session / 'runtime.json', dict(command=command, initial_available_bytes=initial_available,
              baseline_available_bytes=baseline,
              reserve_fraction=cfg['RAM_reserve_fraction'], max_wall_s=cfg['max_wall_s'],
              load_average=os.getloadavg(), platform=sys.platform, concurrency=1))
        with (session / 'stdout.log').open('wb') as stdout, (session / 'stderr.log').open('wb') as stderr:
            child = subprocess.Popen(command, cwd=root, env=env, stdout=stdout, stderr=stderr, start_new_session=True)
            write(session / 'identity.json', dict(matlab_pid=child.pid, process_group=child.pid,
                  supervisor_pid=os.getpid(), session=str(session), root=str(root)))
            minimum = baseline
            peak_rss = 0
            stop_reason = None
            with (session / 'resource.tsv').open('w', encoding='ascii', buffering=1) as resource:
                resource.write('elapsed_s\tMemAvailable_bytes\towned_matlab_RSS_bytes\n')
                while child.poll() is None:
                    current = available()
                    minimum = min(minimum, current)
                    rss = 0
                    try:
                        for line in Path(f'/proc/{child.pid}/status').read_text().splitlines():
                            if line.startswith('VmRSS:'):
                                rss = int(line.split()[1])*1024
                    except FileNotFoundError:
                        pass
                    peak_rss = max(peak_rss, rss)
                    elapsed = time.monotonic()-started
                    resource.write(f'{elapsed:.3f}\t{current}\t{rss}\n')
                    if current < baseline*cfg['RAM_reserve_fraction']:
                        stop_reason = 'RAM_RESERVE_STOP'
                    if elapsed > cfg['max_wall_s']:
                        stop_reason = 'WALL_TIME_STOP'
                    if interrupted[0]:
                        stop_reason = 'OWNED_SESSION_INTERRUPTED'
                    if stop_reason:
                        os.killpg(child.pid, signal.SIGTERM)
                        try:
                            child.wait(timeout=15)
                        except subprocess.TimeoutExpired:
                            os.killpg(child.pid, signal.SIGKILL)
                        break
                    time.sleep(1)
            code = child.wait()
        cp = root / 'output/checkpoint.json'
        valid = False
        if code == 0 and cp.is_file() and not stop_reason:
            data = json.loads(cp.read_text())
            verify(root / 'output', data['artifacts_sha256'])
            for c in cfg['cases']:
                case = root / 'output' / f"seed{c['seed']}_U{c['U_PM_mps']}_N{c['N']:02d}"
                commit = json.loads((case / 'checkpoint.json').read_text())
                assert commit['status'] == 'AUDITED_PASS'
                verify(case, commit['artifacts_sha256'])
            valid = data['status'] == 'K0_COMPLETE_SURFACE_ONLY' and data['case_count'] == 45
        result.update(status='K0_COMPLETE_SURFACE_ONLY' if valid else 'FAILED_STOPPED',
                      matlab_exit_code=code, supervisor_exit_code=0 if valid else 1,
                      minimum_available_bytes=minimum, owned_pid_peak_rss_bytes=peak_rss,
                      stop_reason=stop_reason, audit_checkpoint=str(cp))
    except Exception as exc:
        result['error'] = repr(exc)
        if child is not None and child.poll() is None:
            os.killpg(child.pid, signal.SIGTERM)
            child.wait(timeout=15)
    finally:
        result['elapsed_s'] = time.monotonic()-started
        write(session / 'exit.json', result)
    return result['supervisor_exit_code']


if __name__ == '__main__':
    sys.exit(main())
