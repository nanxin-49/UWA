"""Build a new manual K1 recovery snapshot; no network/MATLAB/physics calls."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import tarfile

REPO=Path(__file__).resolve().parents[3]
OLD=REPO/'results/comparisons/pm_wind_cutoff_k1/preparation_20261118_05'
SOURCE=REPO/'results/comparisons/pm_wind_cutoff_k1/k1_pm_20261118_05_retrieved_20261118_175259'
DIAG=REPO/'results/validation/k1_recovery_preparation/20261118_01'


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p): return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def write(p,v): Path(p).write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8')
def replace(text,old,new):
    assert text.count(old)==1, 'Expected one code fragment: '+old[:100]
    return text.replace(old,new,1)


def generate():
    """Copy frozen orchestration; change only logging/recovery provenance."""
    original=(OLD/'code/run_k1_piston_pm_job.m').read_text(encoding='utf-8-sig')
    text=original.replace('run_k1_piston_pm_job','run_k1_recovery_model_job')
    for old,new in [
        ("if strcmp(job_id,'models')", "if strcmp(job_id,'models')\n    k1_recovery_event(run_dir,'MODELS','BEGIN');\n    k1_recovery_event(run_dir,'INPUT','BEGIN');"),
        ("    direct=baffled_strip_piston_field_vertical(x,z,source);", "    k1_recovery_event(run_dir,'INPUT','END');\n    k1_recovery_event(run_dir,'DIRECT','BEGIN');\n    direct=baffled_strip_piston_field_vertical(x,z,source);\n    k1_recovery_event(run_dir,'DIRECT','END');"),
        ("    flat=-baffled_strip_piston_field_vertical(x,-z,source);", "    k1_recovery_event(run_dir,'FLAT','BEGIN');\n    flat=-baffled_strip_piston_field_vertical(x,-z,source);\n    k1_recovery_event(run_dir,'FLAT','END');"),
        ("    for j=1:numel(PS), PS{j}=local_ps(x,flat,source,profile,recipe(j,:)); end", "    for j=1:numel(PS)\n        k1_recovery_event(run_dir,sprintf('PS%d',j),'BEGIN');\n        PS{j}=local_ps(x,flat,source,profile,recipe(j,:));\n        k1_recovery_event(run_dir,sprintf('PS%d',j),'END');\n    end"),
        ("    for j=1:4, FK{j}=closeout_s3a_pm_fullfk(src,profile,x,z,recipes(j,1),recipes(j,2),flat); end", "    for j=1:4\n        k1_recovery_event(run_dir,sprintf('FullFK%d',j),'BEGIN');\n        FK{j}=closeout_s3a_pm_fullfk(src,profile,x,z,recipes(j,1),recipes(j,2),flat);\n        k1_recovery_event(run_dir,sprintf('FullFK%d',j),'END');\n    end"),
        ("    save(fullfile(run_dir,'common.mat')", "    k1_recovery_event(run_dir,'COMMON_SAVE','BEGIN');\n    save(fullfile(run_dir,'common.mat')"),
        ("    receipt.metrics=struct('ps_window'", "    k1_recovery_event(run_dir,'COMMON_SAVE','END');\n    k1_recovery_event(run_dir,'MODEL_GATES','BEGIN');\n    receipt.metrics=struct('ps_window'"),
        ("    assert(all(fk_values<cfg.fk_gate),'FullFK convergence gate');", "    assert(all(fk_values<cfg.fk_gate),'FullFK convergence gate');\n    k1_recovery_event(run_dir,'MODEL_GATES','END');\n    k1_recovery_event(run_dir,'MODELS','END');")]: text=replace(text,old,new)
    (REPO/'scripts/comparisons/run_k1_recovery_model_job.m').write_text(text,encoding='utf-8')
    worker=(OLD/'code/run_k1_pm_case.m').read_text(encoding='utf-8-sig').replace('run_k1_pm_case','run_k1_pm_recovery_case').replace('run_k1_piston_pm_job','run_k1_recovery_model_job')
    old="        local_models_check(dest,m.seed);\n        k1_pm_io('seal',dest,'committed.json'"
    assert worker.count(old)==2
    worker=worker.replace(old,"        k1_recovery_event(dest,'SAVED_MODEL_AUDIT','BEGIN');\n        local_models_check(dest,m.seed);\n        k1_recovery_event(dest,'SAVED_MODEL_AUDIT','END');\n        k1_pm_io('seal',dest,'committed.json'")
    old="{'common.mat','jobs/models/receipt.json','solved.json','solve_invocation.json'}"
    assert worker.count(old)==2
    worker=worker.replace(old,"{'common.mat','jobs/models/receipt.json','solved.json','solve_invocation.json','recovery_events.jsonl'}")
    (REPO/'scripts/validation/run_k1_pm_recovery_case.m').write_text(worker,encoding='utf-8')
    summary=(OLD/'code/summarize_k1_pm_cases.m').read_text(encoding='utf-8-sig').replace('summarize_k1_pm_cases','summarize_k1_recovery_cases')
    summary=summary.replace("strcmp(cfg.kind,'K1_PM_SIX_NEW_SURFACES')","strcmp(cfg.kind,'K1_PM_RECOVERY_FIVE_REMAINING')")
    summary=replace(summary,"writetable(metrics,fullfile(out,'seven_condition_metrics.csv'));", "metrics.BIE_origin=repmat(\"this_recovery_independent\",height(metrics),1);\nmetrics.BIE_origin(metrics.case_id==\"U6_N07\")=\"imported_audited_N7\";\nmetrics.BIE_origin(metrics.case_id==\"U6_N12_ANCHOR\")=\"qualified_R1_anchor\";\nwritetable(metrics,fullfile(out,'seven_condition_metrics.csv'));")
    summary=summary.replace("'new_BIE_invocations',30", "'new_BIE_invocations',25,'total_BIE_invocations_including_import',30,'imported_BIE_invocations',5")
    summary=summary.replace("receipt.new_BIE_invocations=0; receipt.new_physical_cases=0;", "receipt.new_BIE_invocations=0; receipt.total_BIE_invocations_including_import=0; receipt.imported_BIE_invocations=0; receipt.new_physical_cases=0;")
    summary=summary.replace('六个新海面各完成独立PS六层、BIE C0–C4和完整FullFK四层及原独立审计；U6/N12是严格匹配的历史R1锚点。','N7从已审查旧K1任务导入；其余五面在本恢复任务独立完成原PS/BIE/FullFK层级。U6/N12为严格匹配的历史R1锚点。')
    (REPO/'scripts/reporting/summarize_k1_recovery_cases.m').write_text(summary,encoding='utf-8')
    manage=(OLD/'code/manage_k1_pm.ps1').read_text(encoding='utf-8-sig').replace('code/k1_pm_preflight.py','code/k1_pm_recovery_preflight.py').replace('code/launch_k1_pm.sh','code/launch_k1_pm_recovery.sh')
    (REPO/'scripts/validation/support/manage_k1_pm_recovery.ps1').write_text(manage,encoding='utf-8')
    launch=(OLD/'code/launch_k1_pm.sh').read_text(encoding='utf-8-sig').replace('k1_pm_preflight.py','k1_pm_recovery_preflight.py').replace('k1_pm_supervisor.py','k1_pm_recovery_supervisor.py').replace('launch_k1_pm.sh','launch_k1_pm_recovery.sh')
    (REPO/'scripts/validation/support/launch_k1_pm_recovery.sh').write_text(launch,encoding='utf-8',newline='\n')


def prepare(pkg):
    assert not pkg.exists(); pkg.mkdir(parents=True)
    # Input snapshot is not a result cache: exact lineage evidence for this recovery.
    lineage=pkg/'lineage'; lineage.mkdir()
    oldcopy=lineage/'original_inputs'; oldcopy.mkdir()
    for name in ['code','cases','anchor','provenance','resource_evidence']: shutil.copytree(OLD/name,oldcopy/name)
    for name in ['config.json','manifest.json']: shutil.copy2(OLD/name,oldcopy/name)
    shutil.copytree(SOURCE,lineage/'failed_task',ignore=shutil.ignore_patterns('result.tar','__pycache__','prefs'))
    for item in read(SOURCE/'retrieval/manifest.json'): assert sha(SOURCE/item['path'])==item['sha256']
    shutil.copytree(OLD/'anchor',pkg/'anchor')
    (pkg/'cases').mkdir(); shutil.copytree(SOURCE/'cases/U6_N07',pkg/'cases/U6_N07')
    write(lineage/'N7_import.json',dict(status='EXACT_IMMUTABLE_N7_IMPORT',source_master_manifest_sha256=sha(SOURCE/'manifest.json'),files=[dict(path=p.relative_to(SOURCE/'cases/U6_N07').as_posix(),sha256=sha(p)) for p in sorted((SOURCE/'cases/U6_N07').rglob('*')) if p.is_file()]))
    shutil.copy2(REPO/'results/validation/k1_retrieval_review/20261118_01/numeric_01/numeric_receipt.json',lineage/'N7_local_numeric_receipt.json')
    shutil.copy2(REPO/'results/validation/k1_retrieval_review/20261118_01/owned_crash_summary.json',lineage/'N40_crash_summary.json')
    (lineage/'diagnostics').mkdir()
    shutil.copytree(DIAG/'A1_A2_retrieved',lineage/'diagnostics/A1_A2')
    shutil.copytree(DIAG/'A0_Linux_retrieved',lineage/'diagnostics/Linux_graphics')
    shutil.copytree(OLD/'code',pkg/'code')
    shutil.copy2(OLD/'code/k1_pm_preflight.py',pkg/'code/k1_legacy_preflight.py')
    additions=[
        'scripts/comparisons/run_k1_recovery_model_job.m','scripts/validation/run_k1_pm_recovery_case.m',
        'scripts/validation/validate_k1_pm_recovery_preparation.m','scripts/validation/validate_k1_recovery_runtime_graphics.m',
        'scripts/reporting/summarize_k1_recovery_cases.m','scripts/validation/support/k1_recovery_event.m',
        'scripts/validation/support/k1_pm_recovery_preflight.py','scripts/validation/support/k1_pm_recovery_supervisor.py',
        'scripts/validation/support/manage_k1_pm_recovery.ps1','scripts/validation/support/launch_k1_pm_recovery.sh']
    for name in additions: shutil.copy2(REPO/name,pkg/'code'/Path(name).name)
    cfg=read(OLD/'config.json'); cfg.update(kind='K1_PM_RECOVERY_FIVE_REMAINING',remote_task_name=pkg.name.replace('recovery_preparation_','k1_pm_recovery_'),
        inherited_used_wall_s=724,inherited_RAM_baseline_bytes=119021486080,N40_models_limit_s=180,new_BIE_invocations_max=25,imported_BIE_invocations=5,
        prior_case_wall_s=dict(U6_N40=18.833152669016272,U5_N12=0,U8_N12=0,U5_N40=0,U8_N40=0),graphics_qualification='R2025b_Update1_GLNXA64_original_gates_PASS')
    write(pkg/'config.json',cfg)
    for c in cfg['cases'][1:]:
        case=pkg/'cases'/c['case_id']; case.mkdir(); shutil.copytree(OLD/'cases'/c['case_id']/'input',case/'input')
        shutil.copytree(pkg/'code',case/'code'); shutil.copy2(OLD/'cases'/c['case_id']/'policy.json',case/'policy.json')
        manifest=read(OLD/'cases'/c['case_id']/'manifest.json')
        manifest.update(parent_case_manifest_sha256=sha(OLD/'cases'/c['case_id']/'manifest.json'),parent_master_manifest_sha256=sha(OLD/'manifest.json'),recovery_attempt_family='explicit_after_preserved_failure_v1')
        manifest['files']=[dict(path=p.relative_to(case).as_posix(),sha256=sha(p)) for p in sorted(case.rglob('*')) if p.is_file()]
        write(case/'manifest.json',manifest)
    (pkg/'provenance').mkdir(); shutil.copy2(OLD/'code/run_k1_piston_pm_job.m',pkg/'provenance/original_model_job.m')
    write(pkg/'provenance/recovery.json',dict(original_failed_source=str(SOURCE),old_used_wall_s=636,nonacoustic_charged_s=88,inherited_used_wall_s=724,new_acoustic_solves_started=0,Windows_layout_fix='UNRESOLVED_RETAINED_FAIL',execution_platform_graphics='Ubuntu original gates verified; runtime pre-acoustic replay required'))
    files=[dict(path=p.relative_to(pkg).as_posix(),sha256=sha(p)) for p in sorted(pkg.rglob('*')) if p.is_file()]
    write(pkg/'manifest.json',dict(schema_version=1,scope='K1_PM_RECOVERY_MASTER',files=files,anchor_proof_sha256=sha(pkg/'anchor/anchor_proof.json')))
    print(json.dumps(dict(status='RECOVERY_SNAPSHOT_READY_UNEXECUTED',files=len(files),path=str(pkg))))


def archive(pkg):
    assert not (pkg/'input_package.tar').exists()
    target=pkg/'input_package.tar'
    with tarfile.open(target,'x') as tar:
        for name in ['code','cases','anchor','lineage','provenance','config.json','manifest.json']: tar.add(pkg/name,arcname=name)
    (pkg/'input_package.tar.sha256').write_text(sha(target)+'  input_package.tar\n')
    print(json.dumps(dict(archive_sha256=sha(target),bytes=target.stat().st_size)))


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('action',choices=['generate','prepare','archive']); parser.add_argument('package',nargs='?')
    args=parser.parse_args()
    if args.action=='generate': generate()
    elif args.action=='prepare': prepare(Path(args.package).resolve())
    else: archive(Path(args.package).resolve())
