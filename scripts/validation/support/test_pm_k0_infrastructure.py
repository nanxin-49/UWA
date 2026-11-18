"""Offline tests; all process/network actions mocked. No MATLAB launch."""
import contextlib
import importlib.util
import io
import json
from pathlib import Path
import shutil
import sys
from unittest.mock import patch, Mock


def main():
    package = Path(sys.argv[1]).resolve()
    output = Path(sys.argv[2]).resolve()
    output.mkdir(parents=True, exist_ok=False)
    sys.path.insert(0, str(package / 'code'))
    import pm_k0_preflight as pre
    import pm_k0_supervisor as sup
    tests = 0
    log = []
    root = output / 'fixture'
    shutil.copytree(package / 'code', root / 'code')
    shutil.copytree(package / 'inputs', root / 'inputs')
    shutil.copytree(package / 'provenance', root / 'provenance')
    for name in ['config.json', 'manifest.json']:
        shutil.copy2(package / name, root / name)
    entries = json.loads((root / 'manifest.json').read_text())['files']
    pre.verify(root, entries); tests += 1
    for bad in ['../outside', '/outside', 'C:/outside', 'cash/file', 'a\\b']:
        try:
            pre.verify(root, [dict(path=bad, sha256='0'*64)])
        except (AssertionError, FileNotFoundError):
            tests += 1
        else:
            raise AssertionError('Unsafe path accepted: '+bad)
    victim = root / 'inputs/seed260001.csv'
    original = victim.read_bytes()
    victim.write_bytes(original+b' ')
    try:
        pre.verify(root, entries)
    except AssertionError:
        tests += 1
    else:
        raise AssertionError('Input drift accepted')
    victim.write_bytes(original)
    # No real Popen: lifecycle/resource tests use a mock child and own temp tree.
    for scenario in ['nonzero_exit', 'missing_checkpoint', 'ram_stop', 'wall_stop', 'startup_failure', 'audited_complete']:
        session = root / 'sessions' / scenario
        session.mkdir(parents=True)
        child = Mock()
        child.pid = 123456789
        child.poll.return_value = None if scenario in ['ram_stop','wall_stop'] else 0
        child.wait.return_value = 2 if scenario == 'nonzero_exit' else 0
        if scenario == 'audited_complete':
            # Pure metadata fixture, with no surfaces or physics computation.
            cfg = json.loads((root/'config.json').read_text())
            aggregate = []
            for c in cfg['cases']:
                folder = root/'output'/f"seed{c['seed']}_U{c['U_PM_mps']}_N{c['N']:02d}"
                folder.mkdir(parents=True)
                payload = folder/'mock.txt'
                payload.write_text('MOCK_ONLY_NO_NUMERICAL_EVIDENCE')
                checkpoint = folder/'checkpoint.json'
                checkpoint.write_text(json.dumps(dict(status='AUDITED_PASS',artifacts_sha256=[dict(path='mock.txt',sha256=pre.sha(payload))])))
                aggregate.append(dict(path=checkpoint.relative_to(root/'output').as_posix(),sha256=pre.sha(checkpoint)))
            (root/'output/checkpoint.json').write_text(json.dumps(dict(status='K0_COMPLETE_SURFACE_ONLY',case_count=45,artifacts_sha256=aggregate)))
        baseline = 4*1024**3
        memory = [baseline, baseline//10] if scenario == 'ram_stop' else [baseline, baseline]
        clock = [0, 1801, 1802] if scenario == 'wall_stop' else [0, 1, 2]
        import os
        # Patch platform probes without touching another process or installation.
        with contextlib.ExitStack() as stack:
            stack.enter_context(patch.object(sys, 'argv', ['supervisor',str(root),str(session)]))
            stack.enter_context(patch.object(sup, 'available', side_effect=memory))
            stack.enter_context(patch.object(sup.os, 'getloadavg', return_value=(0,0,0), create=True))
            stack.enter_context(patch.object(sup.os, 'killpg', create=True))
            stack.enter_context(patch.object(sup.signal, 'signal'))
            stack.enter_context(patch.object(sup.time, 'monotonic', side_effect=clock))
            stack.enter_context(patch.object(sup.shutil, 'disk_usage', return_value=Mock(free=2*1024**3)))
            real_file = Path.is_file
            stack.enter_context(patch.object(Path,'is_file',lambda p: True if p.as_posix().endswith('/bin/matlab') else real_file(p)))
            popen = stack.enter_context(patch.object(sup.subprocess,'Popen',return_value=child))
            if scenario == 'startup_failure':
                popen.side_effect = OSError('synthetic launch failure')
            result = sup.main()
        receipt = json.loads((session/'exit.json').read_text())
        if scenario == 'audited_complete':
            assert result == 0 and receipt['status'] == 'K0_COMPLETE_SURFACE_ONLY'
        else:
            assert result != 0 and receipt['status'] != 'K0_COMPLETE_SURFACE_ONLY'
        if scenario == 'ram_stop':
            assert receipt['stop_reason'] == 'RAM_RESERVE_STOP'
        if scenario == 'wall_stop':
            assert receipt['stop_reason'] == 'WALL_TIME_STOP'
        log.append(dict(scenario=scenario, receipt=receipt))
        tests += 1
    receipt = dict(status='OFFLINE_INFRASTRUCTURE_PASS',tests=tests,SSH_connections=0,MATLAB_launches=0,
                   limitations=['Native Linux flock/nohup not executed','No real process-group termination exercised'])
    (output/'details.json').write_text(json.dumps(log,indent=2),encoding='utf-8')
    (output/'receipt.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
    print(json.dumps(receipt))


if __name__ == '__main__':
    main()
