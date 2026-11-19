"""Dual-chain recovery input verification. No MATLAB or acoustic execution."""
import json
from pathlib import Path
import sys
import tarfile
import time
import k1_legacy_preflight as legacy

read=legacy.read; sha=legacy.sha; bound=legacy.bound; verify_rows=legacy.verify_rows
completed_case=legacy.completed_case


def verify(root):
    root=Path(root).resolve(); m=read(root/'manifest.json'); cfg=read(root/'config.json')
    assert m['scope']=='K1_PM_RECOVERY_MASTER' and cfg['kind']=='K1_PM_RECOVERY_FIVE_REMAINING'
    count=verify_rows(root,m['files'])
    target=bool(cfg.get('codex_target_mode',False))
    if target:
        seats=read(root/'lineage/goal_evidence/parallel_license_probe.json')
        graphics_goal=read(root/'lineage/goal_evidence/canonical_graphics.json')
        assert seats['status']=='FOUR_INDEPENDENT_BASE_MATLAB_PROCESSES_OVERLAP_PASS' and seats['failed_checks']==0
        assert seats['max_usable_extra_processes_demonstrated']==4 and seats['all_owned_processes_exited']
        assert graphics_goal['status']=='CANONICAL_SAVED_FIG_TWO_PROCESS_PASS' and graphics_goal['failed_checks']==0
        assert len(graphics_goal['passes'])==2 and len(graphics_goal['PNG_cross_process_SHA'])==2
        assert cfg['max_processes']==4 and cfg['inherited_used_wall_s']==graphics_goal['inherited_after_s']>=818
        prior=root/'lineage/recovery04_failed'
        assert read(prior/'control/accounting.json')['used_wall_s']==790
        verify_rows(prior,read(prior/'retrieval/manifest.json'))
    else:
        assert cfg['max_processes']==1 and cfg['inherited_used_wall_s']==724
    assert cfg['threads']==1 and not cfg['parpool'] and not cfg['automatic_retry']
    assert cfg['RAM_reserve_fraction']==.3 and cfg['case_working_bytes']==cfg['task_RSS_cap_bytes']==12884901888
    assert cfg['max_total_wall_s']==10800 and cfg['max_case_wall_s']==1800 and cfg['max_summary_wall_s']==600
    assert cfg['admission_memory_PSI_full_avg10_max']==5 and cfg['admission_io_PSI_full_avg10_max']==10 and cfg['load1_max']==21
    assert cfg['runtime_memory_PSI_full_avg10_max']==10 and cfg['runtime_io_PSI_full_avg10_max']==30 and cfg['pressure_stop_consecutive_samples']==10
    assert cfg['inherited_RAM_baseline_bytes']>=119021486080
    assert cfg['N40_models_limit_s']==180 and cfg['new_BIE_invocations_max']==25 and cfg['imported_BIE_invocations']==5
    original=root/'lineage/original_inputs'; frozen=legacy.verify(original)
    assert sha(original/'manifest.json')=='815e44af687af278bcd45b5a05ee7be60b638aaf0c97901e2ac117dd082b11ec'
    failed=root/'lineage/failed_task'; oldrows=read(failed/'retrieval/manifest.json'); oldcount=verify_rows(failed,oldrows)
    assert oldcount==575 and read(failed/'control/accounting.json')['used_wall_s']==636
    assert read(root/'lineage/diagnostics/A1_A2/receipt.json')['K1_used_wall_s_after']==666
    graphics=read(root/'lineage/diagnostics/Linux_graphics/receipt.json')
    assert graphics['status']=='LINUX_SAVED_GRAPHICS_ORIGINAL_GATES_PASS' and graphics['failed_checks']==0
    assert graphics['K1_used_wall_s_after']==724
    assert read(root/'lineage/diagnostics/A1_A2/receipt.json')['failed_checks']==0
    local=read(root/'lineage/N7_local_numeric_receipt.json'); assert local['failed_checks']==0
    oldcase=failed/'cases/U6_N07'; imported=root/'cases/U6_N07'
    inherited=completed_case(imported)
    rows=read(root/'lineage/N7_import.json')['files']; verify_rows(imported,rows)
    assert sha(imported/'manifest.json')==sha(oldcase/'manifest.json')
    assert sha(imported/'finalization/committed.json')==sha(oldcase/'finalization/committed.json')
    assert sha(root/'anchor/anchor_proof.json')==sha(original/'anchor/anchor_proof.json')
    assert [(c['case_id'],c['U_PM_mps'],c['modes']) for c in cfg['cases']]==legacy.MATRIX
    for c in cfg['cases'][1:]:
        case=root/'cases'/c['case_id']; old=original/'cases'/c['case_id']; cm=read(case/'manifest.json')
        assert cm['scope']=='K1_PM_CASE' and cm['seed']==260001 and cm['case_id']==c['case_id']
        assert cm['threads']==1 and cm['max_case_processes']==1 and not cm['automatic_retry']
        verify_rows(case,cm['files'])
        assert cm['parent_case_manifest_sha256']==sha(old/'manifest.json')
        assert cm['parent_master_manifest_sha256']==sha(original/'manifest.json')
        for p in (old/'input').iterdir():
            if p.is_file(): assert sha(case/'input'/p.name)==sha(p), 'Input drift '+c['case_id']+'/'+p.name
        assert sha(case/'policy.json')==sha(old/'policy.json')
    mark=read(failed/'cases/U6_N40/model_stage/solve_invocation.json')
    assert mark['invocations']==1 and not (failed/'cases/U6_N40/model_stage/solved.json').exists()
    assert cfg['prior_case_wall_s']['U6_N40']==read(failed/'cases/U6_N40/sessions/20261118T093555_1124145/exit.json')['elapsed_s']
    return dict(status='K1_RECOVERY_DUAL_CHAIN_PREFLIGHT_PASS',files=count,old_files_verified=oldcount,original_inputs=frozen['files'],imported_case=inherited['case_id'],new_cases=5,acoustic_calls=0)


def archive(root):
    root=Path(root).resolve(); legacy.ensure_stopped(root); verify(root)
    active=root/'control/active_children.json'
    if active.exists():
        import os
        for identity in read(active).values():
            path=Path('/proc')/str(identity['pid'])/'stat'
            if path.is_file():
                ticks=path.read_text().rsplit(') ',1)[1].split()[19]
                assert not (ticks==identity['start_ticks'] and os.getpgid(identity['pid'])==identity['process_group']==identity['pid']), 'Owned parallel child still alive; Fetch refused'
    folder=root/'retrieval'; folder.mkdir(exist_ok=True)
    stem='result_'+time.strftime('%Y%m%d_%H%M%S'); target=folder/(stem+'.tar'); assert not target.exists()
    files=[]
    for name in ['code','cases','anchor','lineage','provenance','resource_evidence','sessions','control','summary']:
        if (root/name).exists(): files += [p for p in (root/name).rglob('*') if p.is_file() and not any(x in ['prefs','__pycache__'] for x in p.relative_to(root).parts)]
    files += [root/'config.json',root/'manifest.json']
    for name in ['checkpoint.json','latest_session.txt']:
        if (root/name).is_file(): files.append(root/name)
    rows=[]
    for p in sorted(files):
        assert not p.is_symlink(); rows.append(dict(path=p.relative_to(root).as_posix(),sha256=sha(p)))
    mf=folder/(stem+'.json'); mf.write_text(json.dumps(rows,indent=2))
    with tarfile.open(target,'x') as tar:
        for p in files: tar.add(p,arcname=p.relative_to(root).as_posix(),recursive=False)
        tar.add(mf,arcname='retrieval/manifest.json',recursive=False)
    print('ARCHIVE_PATH='+str(target)); print('ARCHIVE_SHA256='+sha(target))


if __name__=='__main__':
    root=Path(sys.argv[1]).resolve(); action=sys.argv[2] if len(sys.argv)>2 else 'verify'
    if action=='verify': print(json.dumps(verify(root)))
    elif action=='archive': archive(root)
    elif action=='extract': print(json.dumps(legacy.safe_extract(sys.argv[3],root)))
    else: raise ValueError('Unknown action')
