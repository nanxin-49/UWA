"""Offline K1 controller tests. All Popen/kill actions mocked; no MATLAB/SSH."""
import contextlib
import copy
import io
import json
from pathlib import Path
import shutil
import sys
import tarfile
from unittest.mock import Mock,patch


def main():
    package=Path(sys.argv[1]).resolve(); output=Path(sys.argv[2]).resolve(); output.mkdir(parents=True,exist_ok=False)
    sys.path.insert(0,str(package/'code'))
    import k1_pm_preflight as pre
    import k1_pm_supervisor as sup
    tests=[]

    def ok(name): tests.append(name)
    def reject(function,name):
        try: function()
        except Exception: ok(name)
        else: raise AssertionError('Expected rejection: '+name)

    pre.verify(package); ok('real_snapshot_preflight')
    for bad in ['../outside','/outside','C:/outside','a\\b','cash/a']:
        reject(lambda:pre.bound(package,bad),'bound_'+bad)
    cfg=pre.read(package/'config.json'); cfg['_root']=str(package)
    stat=Mock(f_bavail=10**9,f_frsize=4096,f_favail=100000)
    with patch.object(sup,'available',return_value=45*1024**3),patch.object(sup,'pressure',return_value=0),patch.object(sup.os,'getloadavg',return_value=(3,3,3),create=True),patch.object(sup.os,'statvfs',return_value=stat,create=True):
        sup.admission(cfg,45*1024**3); ok('single_case_RAM_admission')
    for name,mem,mp,ip in [('RAM_admission_reject',20,0,0),('measured_memory_PSI_reject',45,27.49,0),('measured_io_PSI_reject',45,0,38.24)]:
        with patch.object(sup,'available',return_value=mem*1024**3),patch.object(sup,'pressure',side_effect=lambda k:mp if k=='memory' else ip),patch.object(sup.os,'getloadavg',return_value=(3,3,3),create=True),patch.object(sup.os,'statvfs',return_value=stat,create=True):
            reject(lambda:sup.admission(cfg,45*1024**3),name)
    identity=dict(pid=123456789,process_group=123456789,start_ticks='old')
    child=Mock(); child.poll.return_value=None
    with patch.object(sup,'start_ticks',return_value='new'),patch.object(sup.os,'getpgid',return_value=123456789,create=True),patch.object(sup.os,'killpg',create=True) as kill:
        reject(lambda:sup.stop_owned(child,identity),'PID_reuse_stop_refused'); kill.assert_not_called()
    with patch.object(sup,'owned_alive',return_value=True),patch.object(sup.os,'killpg',create=True) as kill:
        child.wait.return_value=0; sup.stop_owned(child,identity); assert kill.call_count==1; ok('verified_owned_TERM_only')
    child.wait.side_effect=[sup.subprocess.TimeoutExpired('mock',15),0]
    with patch.object(sup,'owned_alive',side_effect=[True,False]),patch.object(sup.os,'killpg',create=True) as kill:
        reject(lambda:sup.stop_owned(child,identity),'identity_change_before_KILL_refused'); assert kill.call_count==1
    child.wait.side_effect=None

    # Immutable-input mutation test; it never changes the real package.
    altered=output/'mutated_snapshot'; shutil.copytree(package,altered)
    p=altered/'cases/U6_N07/input/receivers.csv'; p.write_bytes(p.read_bytes()+b' ')
    reject(lambda:pre.verify(altered),'receiver_input_SHA_drift_reject')
    altered_config=output/'parallel_mutation'; shutil.copytree(package,altered_config)
    c=pre.read(altered_config/'config.json'); c['max_processes']=2
    (altered_config/'config.json').write_text(json.dumps(c),encoding='utf-8')
    m=pre.read(altered_config/'manifest.json')
    for r in m['files']:
        if r['path']=='config.json': r['sha256']=pre.sha(altered_config/'config.json')
    (altered_config/'manifest.json').write_text(json.dumps(m),encoding='utf-8')
    reject(lambda:pre.verify(altered_config),'parallel_policy_cannot_be_silently_raised')

    # Exercise the real file/checkpoint verifier with a metadata-only fixture.
    fixture=output/'metadata_case'; shutil.copytree(package/'cases/U6_N07',fixture)
    fp=pre.sha(fixture/'manifest.json')
    def json_file(p,value):
        p.parent.mkdir(parents=True,exist_ok=True); p.write_text(json.dumps(value),encoding='utf-8')
    def make_seal(folder,name,signature,files,status):
        json_file(folder/name,dict(schema_version=1,status=status,signature=signature,files=[dict(path=x,sha256=pre.sha(folder/x)) for x in files]))
    model=fixture/'model_stage'; json_file(model/'jobs/models/receipt.json',dict(status='COMPLETE',failed_checks=0))
    json_file(model/'solve_invocation.json',dict(invocations=1,stage='MODELS',manifest_sha256=fp))
    (model/'common.mat').write_text('TEST_ONLY_METADATA_NOT_MATLAB_DATA')
    signature=fp+':MODELS:'+fp
    make_seal(model,'solved.json',signature,['common.mat','jobs/models/receipt.json','solve_invocation.json'],'SOLVED_UNAUDITED')
    make_seal(model,'committed.json',signature,['common.mat','jobs/models/receipt.json','solved.json','solve_invocation.json'],'AUDITED_PASS')
    pred=pre.sha(model/'committed.json')
    names=['bie_coarse','bie_main','bie_self','bie_space','bie_window']
    for j,name in enumerate(names):
        dest=fixture/'jobs'/name; dest.mkdir(parents=True)
        for file in ['data.mat','points.csv','audit_reconstruction.mat']:
            (dest/file).write_text('TEST_ONLY_METADATA_NOT_NUMERICAL_EVIDENCE')
        json_file(dest/'config.json',dict(test_only=True)); json_file(dest/'timings.json',dict(test_only=True))
        json_file(dest/'solve_invocation.json',dict(invocations=1,stage='C'+str(j),manifest_sha256=fp))
        json_file(dest/'audit.json',dict(status='PASS',failed_checks=0,platform_midpoints=[340,408,408,476,476][j],independent_boundary_points=19,
                  checks={k:True for k in ['density_receiver','MAT_CSV','component_closure','incident_rhs','all_midpoints','independent_boundary']}))
        json_file(dest/'gates.json',dict(status='PASS'))
        signature=pred+':C'+str(j)+':'+fp
        make_seal(dest,'solved.json',signature,['data.mat','points.csv','config.json','timings.json','solve_invocation.json'],'SOLVED_UNAUDITED')
        make_seal(dest,'committed.json',signature,['data.mat','points.csv','config.json','timings.json','solve_invocation.json','solved.json','audit.json','audit_reconstruction.mat','gates.json'],'AUDITED_PASS')
        pred=pre.sha(dest/'committed.json')
    final=fixture/'finalization'; work=final/'attempts/attempt_0001/work'
    json_file(final/'pointer.json',dict(relative_work_dir='attempts/attempt_0001/work'))
    json_file(work/'audit/receipt.json',dict(failed_checks=0))
    json_file(work/'independent_final_metrics_audit.json',dict(failed_checks=0,reference_sensitivities=[1e-7]*4,policy_version='r1-bie-reference-v2'))
    json_file(work/'figure_audit/receipt.json',dict(failed_checks=0,checked_error_curves=12,pixel_tolerance=0,records={'F02':dict(pixel_differences=[0,0]),'F04':dict(pixel_differences=[0,0])}))
    paths=['pointer.json']+[p.relative_to(final).as_posix() for p in work.rglob('*') if p.is_file()]
    make_seal(final,'committed.json',pred+':FINAL:'+fp,paths,'AUDITED_PASS')
    json_file(fixture/'checkpoint.json',dict(status='K1_CASE_COMPLETE',failed_checks=0,manifest_sha256=fp,case_id='U6_N07',final_work_dir='finalization/attempts/attempt_0001/work',completed_steps=['MODELS','C0','C1','C2','C3','C4']))
    pre.completed_case(fixture); ok('real_checkpoint_verifier_metadata_fixture')
    p=fixture/'jobs/bie_space/audit.json'; p.write_text(p.read_text()+' ')
    reject(lambda:pre.completed_case(fixture),'sealed_independent_audit_drift_reject')
    uncertain=output/'uncertain_fetch'; (uncertain/'control').mkdir(parents=True)
    (uncertain/'control/pending_spawn.json').write_text('{}')
    reject(lambda:pre.ensure_stopped(uncertain),'Fetch_uncertain_fork_rejected')

    results=[]
    for scenario in ['success_resume','metadata_resume','case_failure','missing_audit','RAM_runtime','RSS_runtime','PSI_runtime','budget','pending_spawn']:
        root=output/scenario/cfg['remote_task_name']; shutil.copytree(package,root)
        session=root/'sessions'/'TEST_ONLY_01'; session.mkdir(parents=True)
        passed=set(); launch_labels=[]; live=[None]; stops=[]
        control=root/'control'; control.mkdir(exist_ok=True)
        if scenario=='budget': sup.write(control/'accounting.json',dict(manifest_sha256=pre.sha(root/'manifest.json'),used_wall_s=10800))
        if scenario=='pending_spawn': sup.write(control/'pending_spawn.json',dict(test_only=True))
        if scenario=='metadata_resume':
            passed.add('U6_N07'); old=root/'cases/U6_N07'; (old/'finalization').mkdir()
            sup.write(old/'finalization/committed.json',dict(test_only=True)); sup.write(old/'checkpoint.json',dict(status='RUNNING',test_only=True))
        base_memory_calls=[0]; pressure_calls=[0]

        class FakeChild:
            def __init__(self,label,code,running):
                self.pid=123456789; self.label=label; self.code=code; self.running=running; self.finished=False
            def poll(self):
                if self.running and not self.finished: return None
                self.finished=True; return self.code
            def wait(self,timeout=None):
                self.finished=True; return self.code

        def fake_completed(case,require_checkpoint=True):
            name=Path(case).name
            if name not in passed: raise AssertionError('No independently audited case fixture')
            if require_checkpoint: assert sup.read(Path(case)/'checkpoint.json').get('status')=='K1_CASE_COMPLETE'
            return dict(case_id=name,status='AUDITED_CASE_PASS',BIE_invocations=5,work='TEST_ONLY')

        def popen(argv,**kwargs):
            expression=argv[-1]
            if 'validate_k1_pm_preparation' in expression:
                label='preflight'; dest=session/'pre_acoustic_acceptance'; dest.mkdir()
                sup.write(dest/'receipt.json',dict(status='K1_PREPARATION_TESTS_PASS',failed_checks=0,new_acoustic_solves=0))
            elif 'summarize_k1_pm_cases' in expression:
                label='summary'; dest=root/'summary/attempt_0001'; dest.mkdir(parents=True)
                sup.write(dest/'receipt.json',dict(status='K1_SUMMARY_AUDITED_PASS',failed_checks=0,new_cases=6,conditions_including_anchor=7))
            else:
                label=next(c['case_id'] for c in cfg['cases'] if c['case_id'] in expression)
                if not (scenario in ['case_failure','missing_audit'] and label=='U6_N40'):
                    passed.add(label)
                    case=root/'cases'/label; (case/'finalization').mkdir(exist_ok=True)
                    sup.write(case/'finalization/committed.json',dict(test_only=True))
                    sup.write(case/'checkpoint.json',dict(test_only=True,status='K1_CASE_COMPLETE'))
            launch_labels.append(label)
            code=1 if scenario=='case_failure' and label=='U6_N40' else 0
            fake=FakeChild(label,code,scenario in ['RAM_runtime','RSS_runtime','PSI_runtime'])
            live[0]=fake; return fake

        def ticks(pid): return 'owned_start' if live[0] is not None and not live[0].finished else None
        def mem():
            base_memory_calls[0]+=1
            return 1024**3 if scenario=='RAM_runtime' and base_memory_calls[0]>=3 else 45*1024**3
        def psi(kind):
            if kind=='memory': pressure_calls[0]+=1
            return 11 if scenario=='PSI_runtime' and kind=='memory' and pressure_calls[0]>1 else 0
        def kill(pid,sig):
            assert pid==123456789; stops.append([pid,int(sig)])
        with contextlib.ExitStack() as stack:
            stack.enter_context(patch.object(sup.subprocess,'Popen',side_effect=popen))
            stack.enter_context(patch.object(sup,'completed_case',side_effect=fake_completed))
            stack.enter_context(patch.object(sup,'available',side_effect=mem))
            stack.enter_context(patch.object(sup,'pressure',side_effect=psi))
            stack.enter_context(patch.object(sup,'start_ticks',side_effect=ticks))
            stack.enter_context(patch.object(sup,'group_rss',return_value=13*1024**3 if scenario=='RSS_runtime' else 1024**3))
            stack.enter_context(patch.object(sup.os,'getpgid',return_value=123456789,create=True))
            stack.enter_context(patch.object(sup.os,'killpg',side_effect=kill,create=True))
            stack.enter_context(patch.object(sup.os,'getloadavg',return_value=(3,3,3),create=True))
            stack.enter_context(patch.object(sup.os,'statvfs',return_value=stat,create=True))
            stack.enter_context(patch.object(sup.signal,'signal'))
            stack.enter_context(patch.object(sup.time,'sleep'))
            r=sup.main(root,session)
            if scenario in ['success_resume','metadata_resume']:
                assert r['exit_code']==0 and len(passed)==6 and launch_labels==['preflight']+[c['case_id'] for c in cfg['cases']]+['summary']
                before=len(launch_labels); resume=root/'sessions/TEST_ONLY_02'; resume.mkdir()
                rr=sup.main(root,resume); assert rr['status']=='K1_ALREADY_COMPLETE_VERIFIED_NO_REPEAT' and len(launch_labels)==before
                ok('full_dispatch_six_cases_then_summary_no_overlap'); ok('completed_resume_no_MATLAB_repeat')
                if scenario=='metadata_resume':
                    assert r['metadata_repair_dispatches']==1
                    ok('audited_final_seal_stale_top_checkpoint_metadata_only_resume')
            elif scenario in ['case_failure','missing_audit']:
                assert r['exit_code']!=0 and launch_labels==['preflight','U6_N07','U6_N40'] and len(passed)==1
                ok(scenario+'_stops_later_dispatch')
            elif scenario in ['RAM_runtime','RSS_runtime','PSI_runtime']:
                assert r['exit_code']!=0 and launch_labels==['preflight'] and len(stops)==1
                ok(scenario+'_stops_only_verified_owned_group')
            else:
                assert r['exit_code']!=0 and not launch_labels
                ok(scenario+'_rejects_before_MATLAB')
        results.append(dict(scenario=scenario,exit=r,launch_labels=launch_labels,mocked_stop_calls=stops))

    # Safe extraction rejects traversal and links before materializing payloads.
    for name,kind in [('traversal','../outside'),('link','symlink')]:
        path=output/(name+'.tar')
        with tarfile.open(path,'w') as tar:
            info=tarfile.TarInfo(kind)
            if name=='link': info.type=tarfile.SYMTYPE; info.linkname='/outside'
            tar.addfile(info)
        reject(lambda:pre.safe_extract(path,output/(name+'_dest')),'unsafe_archive_'+name)
    receipt=dict(status='K1_OFFLINE_INFRASTRUCTURE_PASS',failed_checks=0,tests=len(tests),test_names=tests,
                 real_MATLAB_processes=0,SSH_connections=0,real_process_termination=0,
                 limitations=['Linux native lock inheritance/contention and real process-group signals not exercised'])
    (output/'receipt.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
    (output/'scenario_details.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in receipt.items() if k!='test_names'}))


if __name__=='__main__': main()
