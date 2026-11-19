"""Offline recovery orchestration tests; Popen/signals entirely mocked."""
import contextlib
import copy
import io
import json
from pathlib import Path
import re
import shutil
import sys
from unittest.mock import patch,Mock


def main(package,output):
    package=Path(package).resolve(); output=Path(output).resolve(); assert not output.exists(); output.mkdir(parents=True)
    sys.path.insert(0,str(package/'code'))
    import k1_pm_recovery_preflight as pre
    import k1_pm_recovery_supervisor as sup
    tests=[]; details=[]
    def ok(name): tests.append(name)
    def expect_fail(fn,name):
        try: fn()
        except Exception: ok(name)
        else: raise AssertionError('Missing rejection '+name)
    pre.verify(package); ok('real_1459_file_dual_chain_and_575_legacy_SHA')
    cfg=pre.read(package/'config.json'); assert cfg['inherited_used_wall_s']==724; ok('636_plus_30_plus_58_budget_inherited')
    altered=output/'tamper_snapshot'; shutil.copytree(package,altered)
    for key,bad in [('max_processes',2),('threads',2),('inherited_used_wall_s',636),('inherited_RAM_baseline_bytes',1),('N40_models_limit_s',181),('new_BIE_invocations_max',30)]:
        changed=copy.deepcopy(cfg); changed[key]=bad
        (altered/'config.json').write_text(json.dumps(changed))
        manifest=pre.read(package/'manifest.json')
        for row in manifest['files']:
            if row['path']=='config.json': row['sha256']=pre.sha(altered/'config.json')
        (altered/'manifest.json').write_text(json.dumps(manifest))
        expect_fail(lambda:pre.verify(altered),'semantic_tamper_'+key)
    shutil.copy2(package/'config.json',altered/'config.json'); shutil.copy2(package/'manifest.json',altered/'manifest.json')
    p=altered/'cases/U6_N07/finalization/committed.json'; p.write_bytes(p.read_bytes()+b' ')
    expect_fail(lambda:pre.verify(altered),'imported_N7_seal_drift_rejected')
    jobnames=['bie_coarse','bie_main','bie_self','bie_space','bie_window']
    def put(path,value):
        path.parent.mkdir(parents=True,exist_ok=True); path.write_text(json.dumps(value))
    scenarios=['success','case_failure','graphics_failure','preflight_failure','models_180s',
               'RAM_stop','RSS_stop','PSI_stop','pending_fork','prior_child_alive','bad_accounting',
               'case_budget','total_budget','skip_sealed_N40','metadata_repair_N40']
    for scenario in scenarios:
        task=output/scenario/cfg['remote_task_name']; task.mkdir(parents=True); (task/'sessions/test').mkdir(parents=True)
        shutil.copy2(package/'config.json',task/'config.json'); shutil.copy2(package/'manifest.json',task/'manifest.json')
        control=task/'control'; control.mkdir(); launches=[]; stops=[]; passed={'U6_N07'}; sealed={'U6_N07'}
        for c in cfg['cases']: (task/'cases'/c['case_id']).mkdir(parents=True)
        def finish_case(name):
            passed.add(name); sealed.add(name); case=task/'cases'/name
            put(case/'finalization/committed.json',dict(TEST_ONLY=True)); put(case/'checkpoint.json',dict(status='K1_CASE_COMPLETE'))
            if name!='U6_N07':
                for job in jobnames: put(case/'jobs'/job/'solve_invocation.json',dict(TEST_ONLY=True,invocations=1))
        finish_case('U6_N07')
        if scenario in ['skip_sealed_N40','metadata_repair_N40']:
            finish_case('U6_N40')
            if scenario=='metadata_repair_N40': passed.remove('U6_N40'); put(task/'cases/U6_N40/checkpoint.json',dict(status='INCOMPLETE_METADATA'))
        binding=pre.sha(task/'manifest.json')
        if scenario=='pending_fork': put(control/'pending_spawn.json',dict(TEST_ONLY=True))
        if scenario=='prior_child_alive': put(control/'active_child.json',dict(root=str(task),pid=99,start_ticks='prior'))
        if scenario=='bad_accounting': put(control/'accounting.json',dict(manifest_sha256='WRONG',used_wall_s=999))
        if scenario=='total_budget': put(control/'accounting.json',dict(manifest_sha256=binding,used_wall_s=10800))
        if scenario=='case_budget':
            used=dict(cfg['prior_case_wall_s']); used['U6_N40']=1800
            put(control/'case_accounting.json',dict(manifest_sha256=binding,used_wall_s=used))
        class Clock:
            now=0
            def monotonic(self): self.now+=1; return self.now
        clock=Clock(); current=[None]; counter=[0]; children={}
        class FakeChild:
            def __init__(self,command,**kwargs):
                assert kwargs['start_new_session'] and kwargs['env']['OMP_NUM_THREADS']==kwargs['env']['MKL_NUM_THREADS']=='1'
                assert kwargs['env']['K1_MASTER_ROOT']==str(task)
                expr=command[-1]
                if 'validate_k1_pm_recovery_preparation' in expr: self.label='preflight'
                elif 'validate_k1_recovery_runtime_graphics' in expr: self.label='saved_graphics'
                elif 'summarize_k1_recovery_cases' in expr: self.label='summary'
                else: self.label=next(c['case_id'] for c in cfg['cases'][1:] if c['case_id'] in expr)
                assert self.label!='U6_N07'; launches.append(self.label); current[0]=self.label
                counter[0]+=1; self.pid=1000+counter[0]; self.returncode=None; self.polls=0; children[self.pid]=self
                self.code=1 if (scenario=='case_failure' and self.label=='U6_N40') or (scenario=='graphics_failure' and self.label=='saved_graphics') or (scenario=='preflight_failure' and self.label=='preflight') else 0
                if self.label.startswith('U'):
                    p=task/'cases'/self.label/'control/stage_id.txt'; p.parent.mkdir(parents=True,exist_ok=True); p.write_text('MODELS')
            def poll(self):
                self.polls+=1
                if self.returncode is not None: return self.returncode
                if self.polls==1 or (scenario in ['models_180s','PSI_stop'] and (self.label=='U6_N40' if scenario=='models_180s' else self.label=='preflight')): return None
                return self.code
            def wait(self,timeout=None):
                if self.returncode is None: self.returncode=self.code
                if self.returncode==0:
                    if self.label=='preflight': put(task/'sessions/test/pre_acoustic_acceptance/receipt.json',dict(failed_checks=0,new_acoustic_solves=0))
                    elif self.label=='saved_graphics': put(task/'sessions/test/saved_graphics_acceptance/receipt.json',dict(status='RECOVERY_RUNTIME_SAVED_GRAPHICS_PASS',failed_checks=0))
                    elif self.label=='summary':
                        folder=task/'summary/attempt_0001'; put(folder/'receipt.json',dict(status='K1_SUMMARY_AUDITED_PASS',failed_checks=0,new_BIE_invocations=25)); (folder/'TEST_ONLY.txt').write_text('No scientific arrays computed')
                    else: finish_case(self.label)
                return self.returncode
        def stop(child,identity):
            assert identity['pid']==child.pid; stops.append(child.label); child.returncode=-15
        def completed(case,require_checkpoint=True):
            name=Path(case).name
            assert name in (passed if require_checkpoint else sealed)
            return dict(status='AUDITED_CASE_PASS',case_id=name,BIE_invocations=5)
        def available(): return 1 if scenario=='RAM_stop' and current[0] else 150*1024**3
        def pressure(kind): return 99 if scenario=='PSI_stop' and current[0] else 0
        stat=Mock(f_bavail=10**9,f_frsize=4096,f_favail=100000)
        with contextlib.ExitStack() as stack:
            for target,name,value in [(sup,'verify',lambda *_:None),(sup,'verify_rows',lambda *_:None),(sup,'completed_case',completed),(sup,'capture_owned_crash',lambda *_:None),
                                      (sup.subprocess,'Popen',FakeChild),(sup.time,'monotonic',clock.monotonic),(sup.time,'sleep',lambda *_:None),
                                      (sup.safety,'available',available),(sup.safety,'pressure',pressure),(sup.safety,'group_rss',lambda *_:13*1024**3 if scenario=='RSS_stop' else 1024**3),
                                      (sup.safety,'owned_alive',lambda ident: ident['pid']==99 or (ident['pid'] in children and children[ident['pid']].returncode is None)),(sup.safety,'start_ticks',lambda *_:'TEST_ONLY'),(sup.safety,'stop_owned',stop),
                                      (sup.safety.os,'getloadavg',lambda:(1,1,1)),(sup.safety.os,'statvfs',lambda *_:stat)]:
                stack.enter_context(patch.object(target,name,value,create=True))
            with contextlib.redirect_stdout(io.StringIO()): result=sup.main(task,task/'sessions/test')
            if scenario=='success':
                assert result['exit_code']==0 and result['case_worker_dispatches']==5 and result['new_BIE_invocations']==25
                assert launches==['preflight','saved_graphics','U6_N40','U5_N12','U8_N12','U5_N40','U8_N40','summary']
                account=pre.read(control/'accounting.json'); assert account['used_wall_s']>=724
                used=pre.read(control/'case_accounting.json')['used_wall_s']; assert used['U6_N40']>=18.833152669016272
                prior_launches=len(launches); second=sup.main(task,task/'sessions/test')
                assert second['exit_code']==0 and len(launches)==prior_launches
                ok('complete_restart_verifies_and_launches_nothing')
            elif scenario=='skip_sealed_N40': assert result['exit_code']==0 and result['case_worker_dispatches']==4 and 'U6_N40' not in launches
            elif scenario=='metadata_repair_N40': assert result['exit_code']==0 and result['metadata_repair_dispatches']==1
            elif scenario in ['case_failure','graphics_failure','preflight_failure']:
                assert result['exit_code']!=0 and 'U5_N12' not in launches
            elif scenario in ['models_180s','RAM_stop','RSS_stop','PSI_stop']:
                assert result['exit_code']!=0 and len(stops)==1
                if scenario=='models_180s': assert 'N40_MODELS_180S_STOP' in result['error']
            else:
                assert result['exit_code']!=0
                if scenario in ['pending_fork','prior_child_alive','bad_accounting','total_budget']: assert not launches
                if scenario=='bad_accounting': assert pre.read(control/'accounting.json')['used_wall_s']==999
            ok('scenario_'+scenario)
        details.append(dict(scenario=scenario,result=result,mocked_launches=launches,mocked_stops=stops))
    receipt=dict(status='K1_RECOVERY_OFFLINE_CONTROLLER_PASS',failed_checks=0,tests=len(tests),test_names=tests,real_MATLAB_processes=0,real_SSH=0,real_signals=0)
    (output/'receipt.json').write_text(json.dumps(receipt,indent=2)); (output/'scenarios.json').write_text(json.dumps(details,indent=2))
    print(json.dumps({k:v for k,v in receipt.items() if k!='test_names'}))


if __name__=='__main__': main(*sys.argv[1:])
