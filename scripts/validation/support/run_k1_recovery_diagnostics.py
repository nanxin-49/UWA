"""A1/A2: bounded, serial, nonacoustic MATLAB diagnostics only."""
import datetime
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time


def run(root):
    root=Path(root).resolve()
    assert root.parent==Path('/home/misc/cwx') and root.name.startswith('k1_recovery_diag_')
    sys.path.insert(0,str(root/'code'))
    import k1_pm_supervisor as sup
    import k1_pm_preflight as pre
    started=time.monotonic(); cfg=pre.read(root/'config.json'); cfg['_root']=str(root)
    assert cfg['max_processes']==1 and cfg['threads']==1
    pre.verify_rows(root,pre.read(root/'diagnostic_manifest.json'))
    baseline=max(119021486080,sup.available()); results=[]; child=None; identity=None
    receipt=dict(status='A_DIAGNOSTICS_FAILED_STOPPED',failed_checks=1,new_acoustic_solves=0,old_K1_used_wall_s=636,
                 RAM_baseline_bytes=baseline,completed_modes=[],additional_parallel_license_tested=False)
    def write(path,value): sup.write(path,value)
    try:
        for mode in ['A1','A2']:
            assert time.monotonic()-started<240
            snap=sup.admission(cfg,baseline)
            folder=root/mode; folder.mkdir(exist_ok=False); (folder/'prefs').mkdir()
            command=['/home/misc/MATLAB/R2025b/bin/matlab','-singleCompThread','-batch',f"addpath('{root}/code'); probe_k1_recovery_environment('{root}','{mode}');"]
            env=dict(os.environ,OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MATLAB_PREFDIR=str(folder/'prefs'),PYTHONDONTWRITEBYTECODE='1')
            write(folder/'command.json',dict(argv=command,admission=snap))
            begin=time.monotonic(); stop=None; bad=0; peak=0; low=snap['MemAvailable_bytes']
            with (folder/'stdout.log').open('wb') as stdout,(folder/'stderr.log').open('wb') as stderr:
                child=subprocess.Popen(command,cwd=root,env=env,stdout=stdout,stderr=stderr,start_new_session=True)
                identity=dict(pid=child.pid,process_group=child.pid,start_ticks=sup.start_ticks(child.pid),root=str(root),mode=mode)
                assert sup.owned_alive(identity); write(folder/'identity.json',identity)
                with (folder/'resources.csv').open('w',buffering=1) as f:
                    f.write('elapsed_s,MemAvailable_bytes,owned_group_RSS_bytes,memory_full_avg10,io_full_avg10\n')
                    while child.poll() is None:
                        mem=sup.available(); rss=sup.group_rss(child.pid); mp=sup.pressure('memory'); ip=sup.pressure('io')
                        low=min(low,mem); peak=max(peak,rss)
                        f.write(f'{time.monotonic()-begin:.3f},{mem},{rss},{mp},{ip}\n')
                        bad=bad+1 if mp>10 or ip>30 else 0
                        if mem<.3*baseline: stop='RAM_RESERVE_STOP'
                        elif rss>cfg['task_RSS_cap_bytes']: stop='RSS_STOP'
                        elif bad>=10: stop='PRESSURE_STOP'
                        elif time.monotonic()-begin>=120: stop='MODE_120S_STOP'
                        elif time.monotonic()-started>=240: stop='DIAGNOSTIC_240S_STOP'
                        if stop: sup.stop_owned(child,identity); break
                        time.sleep(1)
                code=child.wait()
            entry=dict(mode=mode,exit_code=code,stop_reason=stop,elapsed_s=time.monotonic()-begin,peak_RSS_bytes=peak,minimum_available_bytes=low,identity=identity)
            write(folder/'exit.json',entry); results.append(entry); child=None; identity=None
            assert code==0 and stop is None, mode+' failed; no retry'
            checked=pre.read(folder/'receipt.json'); assert checked['status']=='NONACOUSTIC_DIAGNOSTIC_PASS' and checked['failed_checks']==0 and checked['new_acoustic_solves']==0
            receipt['completed_modes'].append(mode)
        receipt.update(status='A1_A2_NONACOUSTIC_PASS',failed_checks=0)
    except Exception as exc:
        receipt['error']=repr(exc)
        if child is not None and child.poll() is None:
            try: sup.stop_owned(child,identity)
            except Exception as stop_error: receipt['owned_stop_refused_or_failed']=repr(stop_error)
    finally:
        elapsed=time.monotonic()-started
        receipt.update(elapsed_s=elapsed,charged_wall_s=math.ceil(elapsed),K1_used_wall_s_after=636+math.ceil(elapsed),results=results,
                       UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),all_owned_children_reaped=child is None or child.poll() is not None)
        write(root/'receipt.json',receipt)
    print(json.dumps(receipt)); return 0 if receipt['failed_checks']==0 else 1


if __name__=='__main__': sys.exit(run(sys.argv[1]))
