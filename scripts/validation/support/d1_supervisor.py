"""Single owned MATLAB process, fixed RAM/time/pressure gates, no network/retries."""
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time
from d1_package import preflight, sha, verify

def write(path,value):
    tmp=path.with_suffix('.pending'); tmp.write_text(json.dumps(value,indent=2),encoding='utf-8'); tmp.replace(path)
def available():
    return int(next(s.split()[1] for s in Path('/proc/meminfo').read_text().splitlines() if s.startswith('MemAvailable:')))*1024
def pressure(kind):
    text=(Path('/proc/pressure')/kind).read_text()
    full=next(s for s in text.splitlines() if s.startswith('full '))
    return float(next(w.split('=')[1] for w in full.split() if w.startswith('avg10=')))
def identity(pid):
    t=Path(f'/proc/{pid}/stat').read_text(); tail=t[t.rfind(')')+2:].split()
    return dict(pid=pid,group=int(tail[2]),session=int(tail[3]),start_ticks=int(tail[19]))
def owned_group_rss(group):
    rss=0
    for p in Path('/proc').iterdir():
        if p.name.isdigit():
            try:
                info=identity(int(p.name))
                if info['group']==group and info['session']==group:
                    for line in (p/'status').read_text().splitlines():
                        if line.startswith('VmRSS:'): rss+=int(line.split()[1])*1024
            except (OSError,ValueError,StopIteration): pass
    return rss
def stop_reason(now,baseline,elapsed,depth_elapsed,rss,cfg,active_depth=0):
    prior=cfg.get('recovery',{})
    elapsed+=prior.get('prior_supervisor_elapsed_s',0)
    if active_depth==100: depth_elapsed+=prior.get('prior_supervisor_elapsed_s',0)
    if now<baseline*.3: return 'RAM_RESERVE_STOP'
    if rss>cfg['task_RSS_cap_GiB']*1024**3: return 'OWNED_GROUP_RSS_STOP'
    if elapsed>cfg['max_wall_s']: return 'WALL_TIME_STOP'
    if depth_elapsed>cfg['max_depth_wall_s']: return 'DEPTH_WALL_TIME_STOP'
    return None
def terminate_owned(child,original):
    if child.poll() is not None: return
    assert identity(child.pid)==original and original['group']==original['session']==child.pid
    os.killpg(child.pid,signal.SIGTERM)
    try: child.wait(timeout=15)
    except subprocess.TimeoutExpired:
        assert identity(child.pid)==original
        os.killpg(child.pid,signal.SIGKILL); child.wait()
def validate_output(root,binding):
    cfg,_=preflight(root)
    output=root/'output'; cp=json.loads((output/'checkpoint.json').read_text())
    local_figures='final_figures_mode' in cfg
    expected_status='D1A_NUMERICS_COMPLETE_LOCAL_FIGURES_REQUIRED' if local_figures else 'D1A_COMPLETE_WAITING_REVIEW'
    assert cp['status']==expected_status and cp['manifest_sha256']==binding
    imported=10-cfg['recovery']['remaining_BIE_solves'] if 'recovery' in cfg else 0
    expected_unrun=['Local PNG/FIG reconstruction and zero-pixel audit required'] if local_figures else []
    assert cp['failed_or_unrun_checks']==expected_unrun and cp['new_BIE_solves']==10-imported
    assert cp.get('imported_BIE_solves',0)==imported
    verify(output,cp['artifacts_sha256'])
    for depth in [100,1490]:
        folder=output/f'depth_{depth:04d}m'; dc=json.loads((folder/'checkpoint.json').read_text())
        assert dc['status']=='AUDITED_PASS' and dc['stage']==('DEPTH_NUMERICAL_AUDITED' if local_figures else 'DEPTH_AUDITED') and dc['manifest_sha256']==binding
        verify(folder,dc['artifacts_sha256'])
        for name in ['audit/receipt.json','independent_final_metrics_audit.json']+([] if local_figures else ['figures/figure_audit.json']):
            a=json.loads((folder/name).read_text()); assert a['failed_checks']==0
        for stage in ['bie_coarse','bie_main','bie_self','bie_space','bie_window']:
            p=folder/'jobs'/stage; sc=json.loads((p/'checkpoint.json').read_text())
            assert sc['status']=='AUDITED_PASS' and sc['manifest_sha256']==binding
            verify(p,sc['artifacts_sha256']); assert json.loads((p/'audit.json').read_text())['failed_checks']==0
    preflight(root)

def main(root,session):
    root=Path(root).resolve(); session=Path(session).resolve()
    assert root.parent==Path('/home/misc/cwx') and root.name.startswith('d1_depth_')
    assert session.parent==root/'sessions' and session.is_dir()
    # fd9 must refer to the actual inherited flock lease, not any open descriptor.
    assert os.readlink('/proc/self/fd/9')=='/home/misc/cwx/d1_depth_serial.lock'
    started=time.monotonic(); child=None; original=None
    result=dict(status='FAILED_BEFORE_MATLAB',matlab_exit_code=None,supervisor_exit_code=78)
    interrupted=[False]
    signal.signal(signal.SIGINT,lambda *_:interrupted.__setitem__(0,True))
    signal.signal(signal.SIGTERM,lambda *_:interrupted.__setitem__(0,True))
    try:
        cfg,_=preflight(root); binding=sha(root/'manifest.json'); start_available=available()
        baseline=max(start_available,cfg.get('recovery',{}).get('prior_RAM_baseline_bytes',0))
        assert start_available>=18*1024**3 and start_available-baseline*.3>=12*1024**3, 'RAM admission failed'
        assert pressure('memory')<=5 and pressure('io')<=10 and os.getloadavg()[0]<=21, 'Pressure/load admission failed'
        assert shutil.disk_usage(root).free>=5*1024**3 and os.statvfs(root).f_favail>=5000
        assert not (root/'output').exists(), 'One shot: inspect Status/Fetch; no retry'
        executable='/home/misc/MATLAB/R2025b/bin/matlab'; assert Path(executable).is_file()
        expression=f"addpath('{root}/code'); run_d1_depth_comparison('{root}','dry-run'); run_d1_depth_comparison('{root}','execute');"
        command=[executable,'-singleCompThread','-batch',expression]
        env=os.environ.copy(); env.update(OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',
            D1_SUPERVISED_MANIFEST_SHA256=binding,D1_RAM_BASELINE_BYTES=str(baseline))
        write(session/'runtime.json',dict(command=command,initial_available_bytes=baseline,measured_start_available_bytes=start_available,config=cfg,manifest_sha256=binding))
        minimum=start_available; peak=0; high_pressure=0; reason=None
        with (session/'stdout.log').open('wb') as stdout,(session/'stderr.log').open('wb') as stderr:
            child=subprocess.Popen(command,cwd=root,env=env,stdout=stdout,stderr=stderr,start_new_session=True)
            original=identity(child.pid); assert original['group']==original['session']==child.pid
            write(session/'identity.json',dict(matlab=original,supervisor_pid=os.getpid(),root=str(root)))
            with (session/'resource.tsv').open('w',encoding='ascii',buffering=1) as log:
                log.write('elapsed_s\tMemAvailable_bytes\towned_group_RSS_bytes\tmemory_full_avg10\tio_full_avg10\n')
                while child.poll() is None:
                    now=available(); minimum=min(minimum,now); elapsed=time.monotonic()-started
                    rss=owned_group_rss(child.pid); peak=max(peak,rss); mem=pressure('memory'); io=pressure('io')
                    active=[p for p in (root/'output').glob('depth_*/started.json') if not (p.parent/'checkpoint.json').exists() or json.loads((p.parent/'checkpoint.json').read_text()).get('stage')!='DEPTH_AUDITED']
                    depth_elapsed=max([time.time()-p.stat().st_mtime for p in active]+[0])
                    active_depth=int(active[0].parent.name.split('_')[1][:-1]) if active else 0
                    reason=stop_reason(now,baseline,elapsed,depth_elapsed,rss,cfg,active_depth)
                    high_pressure=high_pressure+1 if mem>10 or io>30 else 0
                    if high_pressure>=10: reason='SUSTAINED_PRESSURE_STOP'
                    if interrupted[0]: reason='OWNED_SESSION_INTERRUPTED'
                    log.write(f'{elapsed:.3f}\t{now}\t{rss}\t{mem}\t{io}\n')
                    if reason: terminate_owned(child,original); break
                    time.sleep(1)
            code=child.wait()
        if code==0 and reason is None: validate_output(root,binding)
        valid=code==0 and reason is None
        success='D1A_NUMERICS_COMPLETE_LOCAL_FIGURES_REQUIRED' if 'final_figures_mode' in cfg else 'D1A_COMPLETE_WAITING_REVIEW'
        result.update(status=success if valid else 'FAILED_STOPPED',matlab_exit_code=code,
            supervisor_exit_code=0 if valid else 1,stop_reason=reason,minimum_available_bytes=minimum,owned_group_peak_RSS_bytes=peak)
    except Exception as exc:
        result['error']=repr(exc)
        if child is not None:
            if original is not None: terminate_owned(child,original)
            result['matlab_exit_code']=child.poll()
    finally:
        result['elapsed_s']=time.monotonic()-started; write(session/'exit.json',result)
        if 'cfg' in locals() and 'recovery' in cfg:
            result['cumulative_elapsed_s']=result['elapsed_s']+cfg['recovery']['prior_supervisor_elapsed_s']
            write(session/'exit.json',result)
    return result['supervisor_exit_code']
if __name__=='__main__': sys.exit(main(*sys.argv[1:3]))
