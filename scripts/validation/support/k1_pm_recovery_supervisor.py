"""One-shot serial recovery: preserve old N7, inherited time, original gates."""
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
import k1_pm_supervisor as safety
from k1_pm_recovery_preflight import read,sha,verify,verify_rows,completed_case

write=safety.write


def capture_owned_crash(folder,root,identity):
    pid=str(identity['pid']); rows=[]
    for parent in [Path(root),Path('/home/misc'),Path('/tmp')]:
        for name in ['matlab_crash_dump.'+pid,'matlab_crash_dump.'+pid+'-1','hs_err_pid'+pid+'.log']:
            path=parent/name
            if path.is_file() and not path.is_symlink() and path.stat().st_size<=2_000_000:
                raw=path.read_text(errors='replace')
                clean='\n'.join(line for line in raw.splitlines() if 'Session Key' not in line and 'Entitlement ID' not in line)
                rows.append(dict(path=str(path),original_sha256=sha(path),bytes=path.stat().st_size,redacted_text=clean))
    write(folder/'owned_crash_evidence.json',dict(identity=identity,matches=rows,search_scope='exact owned PID only; no directory enumeration'))


def main(root,session):
    root=Path(root).resolve(); session=Path(session).resolve(); assert root/'sessions' in session.parents
    started=time.monotonic(); cfg=None; binding=None; baseline=None; child=None; identity=None
    used_before=724; active_case=None; case_begin=None; case_prior=0; case_used={}; interrupted=[False]; accounting_ready=False
    result=dict(status='K1_RECOVERY_FAILED_STOPPED',exit_code=78,case_worker_dispatches=0,metadata_repair_dispatches=0,
                completed_cases=[],case_exit_codes={},imported_BIE_invocations=5,new_BIE_invocations=0)
    signal.signal(signal.SIGTERM,lambda *_:interrupted.__setitem__(0,True)); signal.signal(signal.SIGINT,lambda *_:interrupted.__setitem__(0,True))
    control=root/'control'; control.mkdir(exist_ok=True)
    def account():
        used=used_before+math.ceil(time.monotonic()-started)
        if accounting_ready:
            write(control/'accounting.json',dict(manifest_sha256=binding,used_wall_s=used,session=session.name,inherited_wall_s=724))
            if active_case is not None:
                case_used[active_case]=case_prior+math.ceil(time.monotonic()-case_begin)
            write(control/'case_accounting.json',dict(manifest_sha256=binding,used_wall_s=case_used))
        return used
    try:
        verify(root); cfg=read(root/'config.json'); cfg['_root']=str(root)
        assert root.name==cfg['remote_task_name']
        binding=sha(root/'manifest.json')
        pending=control/'pending_spawn.json'; active=control/'active_child.json'
        assert not pending.exists(), 'UNCERTAIN_PREVIOUS_SPAWN_STOP'
        if active.exists():
            previous=read(active); assert previous['root']==str(root)
            assert not safety.owned_alive(previous), 'PREVIOUS_OWNED_CHILD_ALIVE'
        if (control/'accounting.json').exists():
            prior=read(control/'accounting.json'); assert prior['manifest_sha256']==binding
            used_before=prior['used_wall_s']; assert used_before>=724
        case_used=dict(cfg['prior_case_wall_s'])
        if (control/'case_accounting.json').exists():
            prior=read(control/'case_accounting.json'); assert prior['manifest_sha256']==binding
            assert set(prior['used_wall_s'])==set(case_used)
            assert all(prior['used_wall_s'][k]>=case_used[k] for k in case_used)
            case_used=prior['used_wall_s']
        accounting_ready=True
        baseline=max(cfg['inherited_RAM_baseline_bytes'],safety.available())
        if (control/'resource_baseline.json').exists():
            prior=read(control/'resource_baseline.json'); assert prior['manifest_sha256']==binding
            baseline=max(baseline,prior['RAM_baseline_bytes'])
        write(control/'resource_baseline.json',dict(manifest_sha256=binding,RAM_baseline_bytes=baseline))
        (control/'RAM_baseline_bytes.txt').write_text(str(baseline)+'\n')
        write(session/'policy_and_baseline.json',dict(config=cfg,manifest_sha256=binding,baseline_bytes=baseline,prior_used_s=used_before,case_prior_s=case_used))
        if (root/'checkpoint.json').exists():
            cp=read(root/'checkpoint.json')
            if cp['status']=='K1_COMPLETE_SEVEN_CONDITIONS':
                assert cp['manifest_sha256']==binding and len(cp['completed_cases'])==6
                for c in cfg['cases']: completed_case(root/'cases'/c['case_id'])
                verify_rows(root,cp['artifacts_sha256'])
                result.update(status='K1_RECOVERY_ALREADY_COMPLETE_VERIFIED_NO_REPEAT',exit_code=0,completed_cases=cp['completed_cases']); return result
        minimum_available=baseline; peak_rss=0

        def run_matlab(label,expression,limit,case=None):
            nonlocal child,identity,active_case,case_begin,case_prior,minimum_available,peak_rss
            assert account()<cfg['max_total_wall_s'], 'TOTAL_BUDGET_STOP'
            snap=safety.admission(cfg,baseline); folder=session/label; folder.mkdir(exist_ok=False)
            if case is not None:
                assert label!='U6_N07', 'IMPORTED_N7_MUST_NEVER_RUN'
                active_case=label; case_begin=time.monotonic(); case_prior=case_used[label]
                assert case_prior<cfg['max_case_wall_s'], 'CASE_CUMULATIVE_BUDGET_STOP'
                folder=case/'sessions'/session.name; folder.mkdir(parents=True,exist_ok=False)
                (case/'control').mkdir(exist_ok=True); (case/'control/active_session.txt').write_text(session.name+'\n')
                (folder/'start_epoch.txt').write_text(str(int(time.time()))+'\n')
                write(folder/'shared_runtime.json',dict(RAM_baseline_bytes=baseline,starting_available_bytes=snap['MemAvailable_bytes'],reserve_fraction=.3,prior_case_wall_s=case_prior))
            prefs=folder/'prefs'; prefs.mkdir()
            command=['/home/misc/MATLAB/R2025b/bin/matlab','-singleCompThread','-batch',expression]
            env=dict(os.environ,MATLAB_PREFDIR=str(prefs),OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',
                     K1_ACTIVE_SESSION=session.name,K1_MASTER_ROOT=str(root),PYTHONDONTWRITEBYTECODE='1')
            write(folder/'command.json',dict(argv=command,admission=snap,case_prior_wall_s=case_prior if case else None))
            write(pending,dict(root=str(root),session=session.name,label=label,command=command))
            begin=time.monotonic(); model_begin=None; bad=0; stop=None
            with (folder/'stdout.log').open('wb') as out,(folder/'stderr.log').open('wb') as err:
                child=subprocess.Popen(command,cwd=root,env=env,stdout=out,stderr=err,start_new_session=True)
                identity=dict(root=str(root),session=session.name,label=label,pid=child.pid,process_group=child.pid,start_ticks=safety.start_ticks(child.pid),command=command)
                assert safety.owned_alive(identity); write(folder/'identity.json',identity); write(active,identity); pending.unlink()
                with (folder/'resource.csv').open('w',buffering=1) as f:
                    f.write('elapsed_s,MemAvailable_bytes,owned_group_RSS_bytes,memory_full_avg10,io_full_avg10,stage,model_event\n')
                    while child.poll() is None:
                        mem=safety.available(); rss=safety.group_rss(child.pid); mp=safety.pressure('memory'); ip=safety.pressure('io')
                        minimum_available=min(minimum_available,mem); peak_rss=max(peak_rss,rss)
                        stage='PREFLIGHT' if case is None else ((case/'control/stage_id.txt').read_text().strip() if (case/'control/stage_id.txt').exists() else 'STARTING')
                        event='' if case is None or not (case/'control/model_event.txt').exists() else (case/'control/model_event.txt').read_text().strip()
                        f.write(f'{time.monotonic()-begin:.3f},{mem},{rss},{mp},{ip},{stage},{event}\n')
                        bad=bad+1 if mp>10 or ip>30 else 0
                        if label=='U6_N40' and stage=='MODELS' and model_begin is None: model_begin=time.monotonic()
                        used=account()
                        if mem<.3*baseline: stop='RAM_RESERVE_STOP'
                        elif rss>cfg['task_RSS_cap_bytes']: stop='OWNED_RSS_BUDGET_STOP'
                        elif bad>=10: stop='SHARED_PRESSURE_STOP'
                        elif time.monotonic()-begin>limit: stop='CASE_OR_PHASE_WALL_STOP'
                        elif case is not None and case_used[label]>=cfg['max_case_wall_s']: stop='CASE_CUMULATIVE_WALL_STOP'
                        elif label=='U6_N40' and stage=='MODELS' and model_begin is not None and time.monotonic()-model_begin>=180: stop='N40_MODELS_180S_STOP'
                        elif used>=cfg['max_total_wall_s']: stop='TOTAL_WALL_STOP'
                        elif interrupted[0]: stop='OWNED_CONTROLLER_SIGNAL_STOP'
                        if stop: safety.stop_owned(child,identity); break
                        time.sleep(1)
                code=child.wait()
            used=account()
            if stop is None and used>=cfg['max_total_wall_s']: stop='TOTAL_WALL_STOP_ON_EXIT'
            if stop is None and case is not None and case_used[label]>=cfg['max_case_wall_s']: stop='CASE_WALL_STOP_ON_EXIT'
            active_case=None; case_begin=None
            write(folder/'exit.json',dict(matlab_exit_code=code,stop_reason=stop,elapsed_s=time.monotonic()-begin,identity=identity))
            if code!=0:
                try: capture_owned_crash(folder,root,identity)
                except Exception as error: write(folder/'crash_capture_failure.json',dict(error=repr(error),identity=identity,matlab_exit_code=code))
            result['case_exit_codes'][label]=code; child=None; identity=None
            assert code==0 and stop is None, 'FAILED_STOPPED '+label+' '+str(stop)
            verify(root)

        acceptance=session/'pre_acoustic_acceptance'
        run_matlab('preflight',f"addpath('{root}/code'); validate_k1_pm_recovery_preparation('{root}','{acceptance}');",180)
        a=read(acceptance/'receipt.json'); assert a['failed_checks']==0 and a['new_acoustic_solves']==0
        graphics=session/'saved_graphics_acceptance'
        run_matlab('saved_graphics',f"addpath('{root}/code'); validate_k1_recovery_runtime_graphics('{root}','{graphics}');",180)
        a=read(graphics/'receipt.json'); assert a['status']=='RECOVERY_RUNTIME_SAVED_GRAPHICS_PASS' and a['failed_checks']==0
        for c in cfg['cases']:
            case=root/'cases'/c['case_id']; final=case/'finalization/committed.json'
            has_final=final.is_file()
            if has_final:
                completed_case(case,require_checkpoint=False)
            cp=read(case/'checkpoint.json') if (case/'checkpoint.json').exists() else {}
            if has_final and cp.get('status')=='K1_CASE_COMPLETE':
                completed_case(case)
            else:
                assert c['case_id']!='U6_N07'
                result['case_worker_dispatches']+=1
                if has_final: result['metadata_repair_dispatches']+=1
                run_matlab(c['case_id'],f"addpath('{case}/code'); run_k1_pm_recovery_case('{case}','run','{session.name}');",cfg['max_case_wall_s'],case)
                completed_case(case)
            result['completed_cases'].append(c['case_id'])
            write(root/'checkpoint.json',dict(status='K1_PARTIAL_CASES_AUDITED',manifest_sha256=binding,completed_cases=result['completed_cases'],failed_checks=0))
        folder=root/'summary'; folder.mkdir(exist_ok=True); target=folder/f"attempt_{len(list(folder.glob('attempt_*')))+1:04d}"
        run_matlab('summary',f"addpath('{root}/code'); summarize_k1_recovery_cases('{root}','{target}');",600)
        s=read(target/'receipt.json'); assert s['status']=='K1_SUMMARY_AUDITED_PASS' and s['failed_checks']==0 and s['new_BIE_invocations']==25
        artifacts=[dict(path=p.relative_to(root).as_posix(),sha256=sha(p)) for p in target.rglob('*') if p.is_file()]
        for c in cfg['cases']:
            for name in ['checkpoint.json','finalization/committed.json']:
                p=root/'cases'/c['case_id']/name; artifacts.append(dict(path=p.relative_to(root).as_posix(),sha256=sha(p)))
        # Count durable new BIE markers; the old imported case has its own 5.
        count=sum((root/'cases'/c['case_id']/'jobs'/j/'solve_invocation.json').is_file() for c in cfg['cases'][1:] for j in ['bie_coarse','bie_main','bie_self','bie_space','bie_window'])
        assert count==25
        write(root/'checkpoint.json',dict(status='K1_COMPLETE_SEVEN_CONDITIONS',manifest_sha256=binding,completed_cases=result['completed_cases'],artifacts_sha256=artifacts,failed_checks=0,new_BIE_invocations=25,imported_BIE_invocations=5,total_BIE_invocations=30,summary_path=target.relative_to(root).as_posix(),next_allowed_step='STOP; user Fetch and review'))
        result.update(status='K1_COMPLETE_SEVEN_CONDITIONS',exit_code=0,new_BIE_invocations=25,minimum_available_bytes=minimum_available,owned_group_peak_RSS_bytes=peak_rss)
    except Exception as exc:
        result['error']=repr(exc)
        if child is not None and child.poll() is None:
            try: safety.stop_owned(child,identity)
            except Exception as error: result['owned_stop_refused_or_failed']=repr(error)
    finally:
        if accounting_ready: account()
        result['elapsed_s']=time.monotonic()-started; result['inherited_used_wall_s']=724
        write(session/'exit.json',result)
    return result


def goal_main(root,session):
    """Codex-authorized bounded case parallelism; stages remain sequential."""
    root=Path(root).resolve(); session=Path(session).resolve(); assert root/'sessions' in session.parents
    started=time.monotonic(); cfg=read(root/'config.json'); cfg['_root']=str(root)
    verify(root); assert cfg.get('codex_target_mode') and cfg['max_processes']==4 and root.name==cfg['remote_task_name']
    binding=sha(root/'manifest.json'); control=root/'control'; control.mkdir(exist_ok=True)
    assert not (control/'pending_spawn.json').exists(), 'UNCERTAIN_PREVIOUS_SPAWN_STOP'
    if (control/'active_children.json').exists():
        for old in read(control/'active_children.json').values():
            assert old['root']==str(root) and not safety.owned_alive(old), 'PREVIOUS_OWNED_CHILD_ALIVE'
    inherited=cfg['inherited_used_wall_s']; used_before=inherited; case_used=dict(cfg['prior_case_wall_s'])
    if (control/'accounting.json').exists():
        prior=read(control/'accounting.json'); assert prior['manifest_sha256']==binding and prior['used_wall_s']>=inherited
        used_before=prior['used_wall_s']
    if (control/'case_accounting.json').exists():
        prior=read(control/'case_accounting.json'); assert prior['manifest_sha256']==binding
        assert set(prior['used_wall_s'])==set(case_used) and all(prior['used_wall_s'][k]>=case_used[k] for k in case_used)
        case_used=prior['used_wall_s']
    baseline=max(cfg['inherited_RAM_baseline_bytes'],safety.available())
    if (control/'resource_baseline.json').exists():
        prior=read(control/'resource_baseline.json'); assert prior['manifest_sha256']==binding; baseline=max(baseline,prior['RAM_baseline_bytes'])
    write(control/'resource_baseline.json',dict(manifest_sha256=binding,RAM_baseline_bytes=baseline))
    (control/'RAM_baseline_bytes.txt').write_text(str(baseline)+'\n')
    write(session/'policy_and_baseline.json',dict(config=cfg,manifest_sha256=binding,baseline_bytes=baseline,prior_used_s=used_before,case_prior_s=case_used))
    interrupted=[False]; signal.signal(signal.SIGTERM,lambda *_:interrupted.__setitem__(0,True)); signal.signal(signal.SIGINT,lambda *_:interrupted.__setitem__(0,True))
    active={}; failure=None; completed=['U6_N07']; peak_total=0; minimum_available=baseline
    result=dict(status='K1_TARGET_FAILED_STOPPED',exit_code=78,case_worker_dispatches=0,completed_cases=completed,case_exit_codes={},imported_BIE_invocations=5,new_BIE_invocations=0,max_active_case_processes=0)
    def account():
        used=used_before+math.ceil(time.monotonic()-started)
        for label,item in active.items():
            if item['case'] is not None: case_used[label]=item['prior_s']+math.ceil(time.monotonic()-item['begin'])
        write(control/'accounting.json',dict(manifest_sha256=binding,used_wall_s=used,inherited_wall_s=inherited,session=session.name))
        write(control/'case_accounting.json',dict(manifest_sha256=binding,used_wall_s=case_used))
        write(control/'active_children.json',{label:item['identity'] for label,item in active.items()})
        return used
    def start(label,expression,limit,case=None):
        assert account()<cfg['max_total_wall_s'] and not failure and label not in active
        snap=safety.admission(cfg,baseline)
        required=cfg['case_working_bytes']*(1+sum(item['case'] is not None for item in active.values()))
        assert safety.available()-.3*baseline>=required, 'CONCURRENT_RAM_ADMISSION_STOP'
        folder=session/label
        prior_s=0
        if case is not None:
            assert label!='U6_N07' and len(active)<cfg['max_processes']; prior_s=case_used[label]
            assert prior_s<1800, 'CASE_CUMULATIVE_BUDGET_STOP'
            folder=case/'sessions'/session.name; folder.mkdir(parents=True,exist_ok=False)
            (case/'control').mkdir(exist_ok=True); (case/'control/active_session.txt').write_text(session.name+'\n')
            (folder/'start_epoch.txt').write_text(str(int(time.time()))+'\n')
            write(folder/'shared_runtime.json',dict(RAM_baseline_bytes=baseline,starting_available_bytes=snap['MemAvailable_bytes'],reserve_fraction=.3,prior_case_wall_s=prior_s))
        else: folder.mkdir(exist_ok=False)
        prefs=folder/'prefs'; prefs.mkdir()
        command=['/home/misc/MATLAB/R2025b/bin/matlab','-singleCompThread','-batch',expression]
        env=dict(os.environ,MATLAB_PREFDIR=str(prefs),OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',K1_ACTIVE_SESSION=session.name,K1_MASTER_ROOT=str(root),PYTHONDONTWRITEBYTECODE='1')
        write(folder/'command.json',dict(argv=command,admission=snap,prior_case_wall_s=prior_s))
        write(control/'pending_spawn.json',dict(root=str(root),label=label,session=session.name,command=command))
        with (folder/'stdout.log').open('wb') as out,(folder/'stderr.log').open('wb') as err:
            child=subprocess.Popen(command,cwd=root,env=env,stdout=out,stderr=err,start_new_session=True)
        identity=dict(root=str(root),session=session.name,label=label,pid=child.pid,process_group=child.pid,start_ticks=safety.start_ticks(child.pid),command=command)
        assert safety.owned_alive(identity); write(folder/'identity.json',identity)
        resources=(folder/'resource.csv').open('w',buffering=1); resources.write('elapsed_s,MemAvailable_bytes,owned_group_RSS_bytes,memory_full_avg10,io_full_avg10,stage,model_event\n')
        active[label]=dict(child=child,identity=identity,folder=folder,case=case,begin=time.monotonic(),prior_s=prior_s,limit=limit,model_begin=None,resources=resources,stop=None)
        write(control/'active_children.json',{key:item['identity'] for key,item in active.items()}); (control/'pending_spawn.json').unlink()
        if case is not None:
            result['case_worker_dispatches']+=1
            result['max_active_case_processes']=max(result['max_active_case_processes'],sum(item['case'] is not None for item in active.values()))
    def stop_all(reason):
        for item in active.values():
            if item['child'].poll() is None:
                item['stop']=reason; safety.stop_owned(item['child'],item['identity'])
    def tick():
        nonlocal failure,peak_total,minimum_available
        mem=safety.available(); mp=safety.pressure('memory'); ip=safety.pressure('io'); minimum_available=min(minimum_available,mem)
        total_rss=0; used=account(); tick.bad=tick.bad+1 if mp>10 or ip>30 else 0
        reason=None
        if mem<.3*baseline: reason='RAM_RESERVE_STOP'
        elif tick.bad>=10: reason='SHARED_PRESSURE_STOP'
        elif used>=10800: reason='TOTAL_WALL_STOP'
        elif interrupted[0]: reason='OWNED_CONTROLLER_SIGNAL_STOP'
        for label,item in list(active.items()):
            p=item['child']; case=item['case']; elapsed=time.monotonic()-item['begin']
            rss=safety.group_rss(p.pid) if p.poll() is None else 0; total_rss+=rss
            stage='PREFLIGHT' if case is None else ((case/'control/stage_id.txt').read_text().strip() if (case/'control/stage_id.txt').exists() else 'STARTING')
            event='' if case is None or not (case/'control/model_event.txt').exists() else (case/'control/model_event.txt').read_text().strip()
            item['resources'].write(f'{elapsed:.3f},{mem},{rss},{mp},{ip},{stage},{event}\n')
            if label=='U6_N40' and stage=='MODELS' and item['model_begin'] is None: item['model_begin']=time.monotonic()
            own=None
            if rss>12884901888: own='OWNED_RSS_BUDGET_STOP'
            elif elapsed>=item['limit']: own='CASE_OR_PHASE_WALL_STOP'
            elif case is not None and case_used[label]>=1800: own='CASE_CUMULATIVE_WALL_STOP'
            elif label=='U6_N40' and stage=='MODELS' and item['model_begin'] is not None and time.monotonic()-item['model_begin']>=180: own='N40_MODELS_180S_STOP'
            if own and p.poll() is None:
                item['stop']=own; safety.stop_owned(p,item['identity']); failure=failure or own
        peak_total=max(peak_total,total_rss)
        if total_rss>cfg['max_processes']*12884901888: reason='TOTAL_OWNED_RSS_STOP'
        if reason: failure=failure or reason; stop_all(reason)
        for label,item in list(active.items()):
            p=item['child']
            if p.poll() is None: continue
            code=p.wait(); elapsed=time.monotonic()-item['begin']; item['resources'].close()
            if item['case'] is not None: case_used[label]=item['prior_s']+math.ceil(elapsed)
            write(item['folder']/'exit.json',dict(matlab_exit_code=code,stop_reason=item['stop'],elapsed_s=elapsed,identity=item['identity']))
            result['case_exit_codes'][label]=code
            if code!=0 or item['stop']:
                failure=failure or ('FAILED_STOPPED '+label+' '+str(item['stop']))
                try: capture_owned_crash(item['folder'],root,item['identity'])
                except Exception as exc: write(item['folder']/'crash_capture_failure.json',dict(error=repr(exc)))
            else:
                try:
                    verify(root)
                    if item['case'] is not None:
                        completed_case(item['case']); completed.append(label)
                        write(root/'checkpoint.json',dict(status='K1_TARGET_PARTIAL_CASES_AUDITED',manifest_sha256=binding,completed_cases=completed,failed_checks=0))
                except Exception as exc: failure=failure or repr(exc)
            del active[label]
        account()
    tick.bad=0
    def standalone(label,expression,limit):
        start(label,expression,limit)
        while active:
            tick()
            if active: time.sleep(1)
        assert not failure, failure
    try:
        completed_case(root/'cases/U6_N07')
        if (root/'checkpoint.json').exists():
            prior=read(root/'checkpoint.json')
            if prior['status']=='K1_COMPLETE_SEVEN_CONDITIONS':
                assert prior['manifest_sha256']==binding; verify_rows(root,prior['artifacts_sha256'])
                for c in cfg['cases']: completed_case(root/'cases'/c['case_id'])
                result.update(status='K1_TARGET_ALREADY_COMPLETE_VERIFIED',exit_code=0,completed_cases=prior['completed_cases']); return result
        accept=session/'pre_acoustic_acceptance'
        standalone('preflight',f"addpath('{root}/code'); validate_k1_pm_recovery_preparation('{root}','{accept}');",180)
        acceptance=read(accept/'receipt.json')
        assert acceptance['status']=='K1_RECOVERY_PREPARATION_TESTS_PASS' and acceptance['failed_checks']==0 and acceptance['new_acoustic_solves']==0
        graph=session/'saved_graphics_acceptance'
        standalone('saved_graphics',f"addpath('{root}/code'); validate_k1_recovery_runtime_graphics('{root}','{graph}');",180)
        acceptance=read(graph/'receipt.json')
        assert acceptance['status']=='RECOVERY_RUNTIME_SAVED_GRAPHICS_PASS' and acceptance['failed_checks']==0 and acceptance['new_acoustic_solves']==0
        queue=[]
        for c in cfg['cases'][1:]:
            case=root/'cases'/c['case_id']
            if (case/'finalization/committed.json').exists() and (case/'checkpoint.json').exists() and read(case/'checkpoint.json').get('status')=='K1_CASE_COMPLETE':
                completed_case(case); completed.append(c['case_id'])
            else: queue.append(c)
        while queue or active:
            if not failure:
                while queue and len(active)<cfg['max_processes']:
                    n40=root/'cases/U6_N40/model_stage/committed.json'
                    if 'U6_N40' in active and not n40.exists(): break
                    if n40.exists():
                        fp=sha(root/'cases/U6_N40/manifest.json')
                        from k1_pm_recovery_preflight import legacy
                        legacy.seal(n40.parent,'committed.json',fp+':MODELS:'+fp)
                    required=cfg['case_working_bytes']*(len(active)+1)
                    if safety.available()-.3*baseline<required:
                        if active: break
                        raise RuntimeError('CASE_RAM_ADMISSION_STOP')
                    c=queue.pop(0); case=root/'cases'/c['case_id']; label=c['case_id']
                    start(label,f"addpath('{case}/code'); run_k1_pm_recovery_case('{case}','run','{session.name}');",1800,case)
                    if label=='U6_N40' and not n40.exists(): break
            if active: tick()
            if failure and not active: break
            if active: time.sleep(1)
        assert not failure and not queue and len(completed)==6, failure or 'Incomplete case set'
        folder=root/'summary'; folder.mkdir(exist_ok=True); target=folder/f"attempt_{len(list(folder.glob('attempt_*')))+1:04d}"
        standalone('summary',f"addpath('{root}/code'); summarize_k1_recovery_cases('{root}','{target}');",600)
        s=read(target/'receipt.json'); assert s['status']=='K1_SUMMARY_AUDITED_PASS' and s['failed_checks']==0 and s['new_BIE_invocations']==25
        count=sum((root/'cases'/c['case_id']/'jobs'/j/'solve_invocation.json').exists() for c in cfg['cases'][1:] for j in ['bie_coarse','bie_main','bie_self','bie_space','bie_window']); assert count==25
        artifacts=[dict(path=p.relative_to(root).as_posix(),sha256=sha(p)) for p in target.rglob('*') if p.is_file()]
        for c in cfg['cases']:
            for name in ['checkpoint.json','finalization/committed.json']:
                p=root/'cases'/c['case_id']/name; artifacts.append(dict(path=p.relative_to(root).as_posix(),sha256=sha(p)))
        write(root/'checkpoint.json',dict(status='K1_COMPLETE_SEVEN_CONDITIONS',manifest_sha256=binding,completed_cases=completed,artifacts_sha256=artifacts,failed_checks=0,new_BIE_invocations=25,imported_BIE_invocations=5,total_BIE_invocations=30,summary_path=target.relative_to(root).as_posix(),next_allowed_step='Codex Fetch and independently review saved results; no expansion'))
        result.update(status='K1_COMPLETE_SEVEN_CONDITIONS',exit_code=0,new_BIE_invocations=25)
    except Exception as exc:
        result['error']=repr(exc)
        # Unexpected controller exceptions stop only verified task-owned jobs.
        if active:
            try: stop_all('CONTROLLER_EXCEPTION_STOP')
            except Exception as error: result['owned_stop_refused_or_failed']=repr(error)
            for label,item in active.items():
                item['resources'].close(); write(item['folder']/'exit.json',dict(matlab_exit_code=item['child'].poll(),stop_reason=item['stop'],identity=item['identity']))
    finally:
        account(); result.update(elapsed_s=time.monotonic()-started,inherited_used_wall_s=inherited,owned_total_peak_RSS_bytes=peak_total,minimum_available_bytes=minimum_available)
        write(session/'exit.json',result)
    return result


if __name__=='__main__':
    try:
        import fcntl
        assert os.readlink('/proc/self/fd/9')=='/home/misc/cwx/k1_pm_serial.lock', 'No inherited K1 task lease'
        fcntl.flock(9,fcntl.LOCK_EX|fcntl.LOCK_NB)
        entry=goal_main if read(Path(sys.argv[1])/'config.json').get('codex_target_mode') else main
        result=entry(sys.argv[1],sys.argv[2])
    except Exception as exc:
        result=dict(status='K1_RECOVERY_START_REFUSED',exit_code=78,error=repr(exc),case_worker_dispatches=0)
        root=Path(sys.argv[1]).resolve(); session=Path(sys.argv[2]).resolve()
        if root/'sessions' in session.parents: write(session/'exit.json',result)
    print(json.dumps(result)); sys.exit(result['exit_code'])
