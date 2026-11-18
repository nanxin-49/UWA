"""Bounded serial K1 controller. Only explicit user Start reaches this file."""
import datetime
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
from k1_pm_preflight import read,sha,verify,verify_rows,completed_case


def write(path,value):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_name(path.name+'.pending')
    temp.write_text(json.dumps(value,indent=2),encoding='utf-8'); temp.replace(path)


def available():
    for line in Path('/proc/meminfo').read_text().splitlines():
        if line.startswith('MemAvailable:'): return int(line.split()[1])*1024
    raise RuntimeError('MemAvailable unavailable')


def pressure(kind):
    lines=Path('/proc/pressure/'+kind).read_text().splitlines()
    for line in lines:
        if line.startswith('full '): return float(dict(x.split('=') for x in line.split()[1:])['avg10'])
    raise RuntimeError('PSI full unavailable: '+kind)


def start_ticks(pid):
    try: return Path(f'/proc/{pid}/stat').read_text().rsplit(') ',1)[1].split()[19]
    except FileNotFoundError: return None


def owned_alive(identity):
    pid=identity['pid']
    return start_ticks(pid)==identity['start_ticks'] and start_ticks(pid) is not None and os.getpgid(pid)==identity['process_group']==pid


def stop_owned(child,identity):
    # Popen object plus Linux start time + unique session/group prevent PID reuse.
    if child.poll() is not None: return
    if not owned_alive(identity): raise RuntimeError('Refuse stop: owned process identity mismatch')
    os.killpg(identity['pid'],signal.SIGTERM)
    try: child.wait(timeout=15)
    except subprocess.TimeoutExpired:
        if not owned_alive(identity): raise RuntimeError('Refuse SIGKILL: identity changed')
        os.killpg(identity['pid'],signal.SIGKILL); child.wait(timeout=10)


def group_rss(pid):
    p=subprocess.run(['ps','-o','rss=','--sid',str(pid)],capture_output=True,text=True,timeout=5)
    if p.returncode!=0: return 0
    return sum(int(v) for v in p.stdout.split())*1024


def admission(cfg,baseline):
    mem=available(); memory_psi=pressure('memory'); io_psi=pressure('io')
    assert mem-.3*baseline>=cfg['case_working_bytes'], 'RAM_ADMISSION_STOP'
    assert memory_psi<=cfg['admission_memory_PSI_full_avg10_max'], 'MEMORY_PRESSURE_ADMISSION_STOP'
    assert io_psi<=cfg['admission_io_PSI_full_avg10_max'], 'IO_PRESSURE_ADMISSION_STOP'
    assert os.getloadavg()[0]<=cfg['load1_max'], 'CPU_LOAD_ADMISSION_STOP'
    disk=os.statvfs(cfg['_root'])
    assert disk.f_bavail*disk.f_frsize>=cfg['minimum_disk_bytes'] and disk.f_favail>=cfg['minimum_free_inodes'], 'DISK_ADMISSION_STOP'
    return dict(MemAvailable_bytes=mem,memory_full_avg10=memory_psi,io_full_avg10=io_psi,load1=os.getloadavg()[0])


def main(root,session):
    root=Path(root).resolve(); session=Path(session).resolve()
    assert root/'sessions' in session.parents
    started=time.monotonic(); child=None; identity=None; cfg=None; base_used=0; binding=None
    result=dict(status='K1_FAILED_BEFORE_DISPATCH',exit_code=78,case_worker_dispatches=0,metadata_repair_dispatches=0,case_exit_codes={},completed_cases=[])
    interrupted=[False]
    signal.signal(signal.SIGTERM,lambda *_:interrupted.__setitem__(0,True))
    signal.signal(signal.SIGINT,lambda *_:interrupted.__setitem__(0,True))
    control=root/'control'; control.mkdir(exist_ok=True)
    try:
        verify(root); cfg=read(root/'config.json'); cfg['_root']=str(root)
        assert root.name==cfg['remote_task_name'] and cfg['max_processes']==1 and cfg['threads']==1
        binding=sha(root/'manifest.json')
        # Reject a surviving previous task-owned child or uncertain fork window.
        pending=control/'pending_spawn.json'
        assert not pending.exists(), 'UNCERTAIN_PREVIOUS_SPAWN_STOP: inspect evidence; no automatic relaunch'
        active=control/'active_child.json'
        if active.exists():
            old=read(active); assert old['root']==str(root)
            assert not owned_alive(old), 'PREVIOUS_OWNED_PROCESS_STILL_ALIVE'
        accounting=control/'accounting.json'
        if accounting.exists():
            prior=read(accounting); assert prior['manifest_sha256']==binding; base_used=prior['used_wall_s']
        baseline=available(); baseline_path=control/'resource_baseline.json'
        if baseline_path.exists():
            prior=read(baseline_path); assert prior['manifest_sha256']==binding; baseline=max(baseline,prior['RAM_baseline_bytes'])
        write(baseline_path,dict(manifest_sha256=binding,RAM_baseline_bytes=baseline))
        (control/'RAM_baseline_bytes.txt').write_text(str(baseline)+'\n')
        write(session/'policy_and_baseline.json',dict(config=cfg,manifest_sha256=binding,baseline_bytes=baseline,prior_used_s=base_used))
        all_complete=[]
        def checkpoint_complete(case):
            return (case/'checkpoint.json').is_file() and read(case/'checkpoint.json').get('status')=='K1_CASE_COMPLETE'
        for c in cfg['cases']:
            case=root/'cases'/c['case_id']
            if (case/'finalization/committed.json').is_file():
                completed_case(case,require_checkpoint=False)
                if checkpoint_complete(case): all_complete.append(completed_case(case))
        if (root/'checkpoint.json').is_file():
            cp=read(root/'checkpoint.json')
            if cp['status']=='K1_COMPLETE_SEVEN_CONDITIONS':
                assert cp['manifest_sha256']==binding and len(all_complete)==6
                verify_rows(root,cp['artifacts_sha256'])
                result.update(status='K1_ALREADY_COMPLETE_VERIFIED_NO_REPEAT',exit_code=0,completed_cases=[c['case_id'] for c in all_complete]); return result
        min_available=baseline; peak_rss=0

        def account():
            used=base_used+math.ceil(time.monotonic()-started)
            write(accounting,dict(manifest_sha256=binding,used_wall_s=used,session=session.name))
            return used

        def run_matlab(label,expression,limit,case=None):
            nonlocal child,identity,min_available,peak_rss
            assert account()<cfg['max_total_wall_s'], 'TOTAL_BUDGET_STOP'
            snap=admission(cfg,baseline)
            folder=session/label; folder.mkdir(exist_ok=False)
            if case is not None:
                case_session=case/'sessions'/session.name; case_session.mkdir(parents=True,exist_ok=False)
                (case/'control').mkdir(exist_ok=True)
                (case/'control/active_session.txt').write_text(session.name+'\n')
                (case_session/'start_epoch.txt').write_text(str(int(time.time()))+'\n')
                write(case_session/'shared_runtime.json',dict(RAM_baseline_bytes=baseline,starting_available_bytes=snap['MemAvailable_bytes'],reserve_fraction=.3))
                folder=case_session
            matlab='/home/misc/MATLAB/R2025b/bin/matlab'
            command=[matlab,'-singleCompThread','-batch',expression]
            env=dict(os.environ,OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',
                     PYTHONDONTWRITEBYTECODE='1',K1_ACTIVE_SESSION=session.name,K1_MASTER_ROOT=str(root))
            # Separate process-local preferences; never modify the shared profile.
            prefs=folder/'prefs'; prefs.mkdir(); env['MATLAB_PREFDIR']=str(prefs)
            write(folder/'command.json',dict(argv=command,case_id=label,admission=snap,thread_environment={k:env[k] for k in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS']}))
            write(pending,dict(root=str(root),session=session.name,label=label,command=command))
            begin=time.monotonic(); stop_reason=None; pressure_bad=0
            with (folder/'stdout.log').open('wb') as out,(folder/'stderr.log').open('wb') as err:
                child=subprocess.Popen(command,cwd=root,env=env,stdout=out,stderr=err,start_new_session=True)
                identity=dict(root=str(root),session=session.name,label=label,pid=child.pid,process_group=child.pid,start_ticks=start_ticks(child.pid),command=command)
                assert identity['start_ticks'] is not None and os.getpgid(child.pid)==child.pid
                write(folder/'identity.json',identity); write(active,identity); pending.unlink()
                with (folder/'resource.csv').open('w',encoding='ascii',buffering=1) as resources:
                    resources.write('elapsed_s,MemAvailable_bytes,owned_group_RSS_bytes,memory_full_avg10,io_full_avg10,stage\n')
                    while child.poll() is None:
                        mem=available(); rss=group_rss(child.pid); mp=pressure('memory'); ip=pressure('io')
                        min_available=min(min_available,mem); peak_rss=max(peak_rss,rss)
                        stage='PREFLIGHT' if case is None else ((case/'control/stage_id.txt').read_text().strip() if (case/'control/stage_id.txt').exists() else 'STARTING')
                        resources.write(f'{time.monotonic()-begin:.3f},{mem},{rss},{mp},{ip},{stage}\n')
                        pressure_bad=pressure_bad+1 if mp>cfg['runtime_memory_PSI_full_avg10_max'] or ip>cfg['runtime_io_PSI_full_avg10_max'] else 0
                        if mem<.3*baseline: stop_reason='RAM_RESERVE_STOP'
                        elif rss>cfg['task_RSS_cap_bytes']: stop_reason='OWNED_RSS_BUDGET_STOP'
                        elif pressure_bad>=cfg['pressure_stop_consecutive_samples']: stop_reason='SHARED_PRESSURE_STOP'
                        elif time.monotonic()-begin>limit: stop_reason='CASE_OR_PHASE_WALL_STOP'
                        elif account()>=cfg['max_total_wall_s']: stop_reason='TOTAL_WALL_STOP'
                        elif interrupted[0]: stop_reason='OWNED_CONTROLLER_SIGNAL_STOP'
                        if stop_reason:
                            stop_owned(child,identity); break
                        time.sleep(1)
                code=child.wait()
            write(folder/'exit.json',dict(matlab_exit_code=code,stop_reason=stop_reason,elapsed_s=time.monotonic()-begin,identity=identity))
            result['case_exit_codes'][label]=code
            child=None; identity=None
            assert code==0 and stop_reason is None, 'FAILED_STOPPED '+label+' '+str(stop_reason)
            verify(root)

        # Pre-acoustic tests exercise immutable geometry, anchor, state and gates.
        acceptance=session/'pre_acoustic_acceptance'
        run_matlab('preflight',f"addpath('{root}/code'); validate_k1_pm_preparation('{root}','{acceptance}');",180)
        accept=read(acceptance/'receipt.json'); assert accept['status']=='K1_PREPARATION_TESTS_PASS' and accept['failed_checks']==0 and accept['new_acoustic_solves']==0
        for c in cfg['cases']:
            case=root/'cases'/c['case_id']
            has_final=(case/'finalization/committed.json').is_file()
            if has_final:
                completed_case(case,require_checkpoint=False)
            if has_final and checkpoint_complete(case):
                passed=completed_case(case)
            else:
                # One process only; any failure aborts before later dispatch.
                result['case_worker_dispatches']+=1
                if has_final: result['metadata_repair_dispatches']+=1
                run_matlab(c['case_id'],f"addpath('{case}/code'); run_k1_pm_case('{case}','run','{session.name}');",cfg['max_case_wall_s'],case)
                passed=completed_case(case)
            result['completed_cases'].append(passed['case_id'])
            write(root/'checkpoint.json',dict(status='K1_PARTIAL_CASES_AUDITED',manifest_sha256=binding,completed_cases=result['completed_cases'],failed_checks=0,next_allowed_step='Remaining frozen cases only; no retries'))
        folder=root/'summary'; folder.mkdir(exist_ok=True)
        attempts=list(folder.glob('attempt_*')); target=folder/f'attempt_{len(attempts)+1:04d}'
        run_matlab('summary',f"addpath('{root}/code'); summarize_k1_pm_cases('{root}','{target}');",cfg['max_summary_wall_s'])
        summary=read(target/'receipt.json'); assert summary['status']=='K1_SUMMARY_AUDITED_PASS' and summary['failed_checks']==0 and summary['new_cases']==6 and summary['conditions_including_anchor']==7
        artifacts=[dict(path=p.relative_to(root).as_posix(),sha256=sha(p)) for p in target.rglob('*') if p.is_file()]
        for c in cfg['cases']:
            for name in ['checkpoint.json','finalization/committed.json']:
                p=root/'cases'/c['case_id']/name; artifacts.append(dict(path=p.relative_to(root).as_posix(),sha256=sha(p)))
        write(root/'checkpoint.json',dict(status='K1_COMPLETE_SEVEN_CONDITIONS',manifest_sha256=binding,completed_cases=result['completed_cases'],
              artifacts_sha256=artifacts,failed_checks=0,summary_path=target.relative_to(root).as_posix(),new_BIE_invocations=30,next_allowed_step='STOP; user Fetch and review; no expansion'))
        result.update(status='K1_COMPLETE_SEVEN_CONDITIONS',exit_code=0,minimum_available_bytes=min_available,owned_group_peak_RSS_bytes=peak_rss)
    except Exception as exc:
        result['error']=repr(exc)
        if child is not None and child.poll() is None:
            try: stop_owned(child,identity)
            except Exception as stop_error: result['identity_stop_refused_or_failed']=repr(stop_error)
        result['status']='K1_FAILED_STOPPED_NO_FURTHER_DISPATCH'
    finally:
        result['elapsed_s']=time.monotonic()-started
        if binding is not None:
            write(control/'accounting.json',dict(manifest_sha256=binding,used_wall_s=base_used+math.ceil(result['elapsed_s']),session=session.name))
        write(session/'exit.json',result)
    return result


if __name__=='__main__':
    try:
        # The CLI is only valid with the launcher's shared K1 lease inherited.
        import fcntl
        assert os.readlink('/proc/self/fd/9')=='/home/misc/cwx/k1_pm_serial.lock', 'Use launch_k1_pm.sh start; no inherited task lease'
        fcntl.flock(9,fcntl.LOCK_EX|fcntl.LOCK_NB)
        value=main(sys.argv[1],sys.argv[2])
    except Exception as exc:
        value=dict(status='K1_CONTROLLER_START_REFUSED',exit_code=78,error=repr(exc),case_worker_dispatches=0)
        task=Path(sys.argv[1]).resolve(); session=Path(sys.argv[2]).resolve()
        if task/'sessions' in session.parents: write(session/'exit.json',value)
    print(json.dumps(value)); sys.exit(value['exit_code'])
