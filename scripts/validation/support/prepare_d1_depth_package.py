"""Allowlisted D1 snapshot builder; no MATLAB, subprocesses, network or archive discovery."""
from pathlib import Path
import hashlib
import json
import shutil
import sys
import tarfile

ROOT = Path(__file__).resolve().parents[3]
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p, v): p.write_text(json.dumps(v, indent=2), encoding='utf-8', newline='\n')

def build(destination,origin_dir=None,audit_dir=None,*attempt_dirs):
    dest=Path(destination).resolve()
    assert ROOT/'results/comparisons/d1_depth' in dest.parents and not dest.exists()
    dest.mkdir(parents=True); (dest/'code').mkdir(); (dest/'inputs').mkdir(); (dest/'provenance').mkdir()
    files=[
        'scripts/comparisons/run_d1_depth_comparison.m','scripts/comparisons/run_d1_model_job.m',
        'scripts/validation/audit_d1_models.m','scripts/reporting/render_d1_depth_results.m',
        'scripts/validation/support/d1_bie_gates.m','scripts/validation/support/audit_d1_final.m',
        'scripts/validation/support/write_d1_csv.m','scripts/validation/audit_d1_recovery_origin.m',
        'scripts/validation/support/verify_d1_saved_import.m',
        'scripts/validation/validate_d1_recovery_preparation.m',
        'scripts/validation/audit_d1_saved_tx100.m',
        'scripts/validation/refresh_d1_saved_figure_audit.m',
        'scripts/validation/support/d1_package.py','scripts/validation/support/d1_supervisor.py',
        'scripts/validation/support/launch_d1_depth.sh','scripts/validation/support/manage_d1_depth.ps1',
        'scripts/validation/support/invoke_r1_bie_status_stdin.ps1',
        'scripts/validation/support/solve_helmholtz_bie_halfplane_vertical.m',
        'scripts/validation/support/solve_helmholtz_bie_halfplane_cached_vertical.m',
        'scripts/validation/support/audit_r1_bie_manual_stage.m','scripts/validation/support/r1_bie_manual_io.m',
        'scripts/validation/support/r1_bie_continuation_io.m','scripts/validation/support/r2_pm_io.m',
        'scripts/validation/support/closeout_s3a_pm_surface.m','scripts/validation/support/closeout_s3a_pm_fullfk.m',
        'scripts/validation/support/fk_single_receiver_incident.m',
        'src/propagation/baffled_strip_piston_field_vertical.m',
        'scripts/validation/support/piston_independent_angular_spectrum_vertical.m',
        'scripts/validation/support/d0_incident_spectrum_reference.m','scripts/experiments/pm_sensitivity_io.m']
    sources=[]
    for rel in files:
        p=ROOT/rel; shutil.copyfile(p,dest/'code'/p.name); sources.append(dict(path=rel,sha256=sha(p)))
    k0=ROOT/'results/experiments/pm_wind_cutoff/pm_wind_cutoff_k0_20261117_05_retrieved_20261117_232303/output/seed260001_U6_N12'
    k0cp=json.loads((k0/'checkpoint.json').read_text(encoding='utf-8-sig'))
    assert k0cp['status']=='AUDITED_PASS'
    # Copy only exact frozen inputs, never density or a precomputed acoustic field.
    for name in ['coefficients.csv','checkpoint.json']:
        shutil.copyfile(k0/name,dest/'provenance'/('k0_'+name))
        sources.append(dict(path=(k0/name).relative_to(ROOT).as_posix(),sha256=sha(k0/name)))
    origin=k0/'coefficients.csv'; lines=origin.read_text().splitlines()
    mapped=['k_rad_per_m,cos_coeff_m,sin_coeff_m,Sk_m3']
    for line in lines[1:]:
        row=line.split(','); mapped.append(','.join([row[1],row[3],row[4],row[2]]))
    (dest/'inputs/coefficients.csv').write_text('\n'.join(mapped)+'\n',encoding='utf-8')
    # Each depth evaluates the same 411 coordinates directly, not an interpolated R1 line.
    (dest/'inputs/receivers.csv').write_text('x_m,z_m\n'+''.join(f'{(-36+72*i/410):.17g},6\n' for i in range(411)),encoding='utf-8')
    policy=ROOT/'results/comparisons/piston_2d_pm/r1_bie_manual_c4_policy_v2_20261116_04_retrieved_20261116_173609/policy.json'
    shutil.copyfile(policy,dest/'inputs/policy.json'); sources.append(dict(path=policy.relative_to(ROOT).as_posix(),sha256=sha(policy)))
    for rel in ['reports/d0_incident_depth_review_2026-11-18.md','reports/k1_resources_strategy_2026-11-18.md',
                'reports/codex_token_execution_protocol.md','AGENTS.md']:
        p=ROOT/rel; sources.append(dict(path=rel,sha256=sha(p)))
    write(dest/'provenance/sources.json',sources)
    cfg=dict(kind='D1A_FIXED_PM_DEPTH_COMPARISON',depths_m=[100,1490],receiver_z_m=6,receiver_halfwidth_m=36,receivers=411,
        source=dict(frequency_hz=4000,c_mps=1500,rho_kgpm3=1025,v0_mps=.001,total_width_m=.1),
        PM=dict(U_PM_mps=6,N=12,seed=260001,period_m=160,plateau_m=42,support_m=50),
        time_convention='exp(-i omega t)',phase_contract_version='2.0.0',representation='absolute_physical',
        workers=1,threads=1,RAM_reserve_fraction=.3,estimated_working_RAM_GiB=12,minimum_available_RAM_GiB=18,
        max_wall_s=5400,max_depth_wall_s=2400,task_RSS_cap_GiB=12,
        bottom_enabled=False,downward_radiation=False,multiple_reflections=False,structure_scattering=False,
        D1b_enabled=False,automatic_retry=False,user_manual_start_only=True,
        BIE_recipes=[[60,10,96],[60,12,96],[60,12,64],[60,14,96],[70,12,96]],
        PS_recipes=[[128,.05,.1875],[256,.05,.1875],[512,.05,.1875],[1024,.05,.1875],[1024,.025,.1875],[1024,.025,.09375]],
        FullFK_recipes=[[400,50],[800,50],[960,60],[1440,90]],
        gates=dict(sensitivity=1e-4,linear=1e-10,boundary=1e-4,source=1e-7,PS=1e-4,PS_step=1e-10,FullFK=1e-5,
            reconstruction=1e-11,resolution_factor=5,phase_mask_fraction=.01))
    if origin_dir is not None:
        from d1_package import preflight,verify
        origin=Path(origin_dir).resolve(); audit=Path(audit_dir).resolve()
        prior,_=preflight(origin)
        rows=json.loads((origin/'retrieval/manifest.json').read_text(encoding='utf-8-sig')); verify(origin,rows)
        proof=json.loads((audit/'receipt.json').read_text(encoding='utf-8-sig'))
        checkpoint=json.loads((audit/'checkpoint.json').read_text(encoding='utf-8-sig'))
        verify(audit,checkpoint['artifacts_sha256'])
        for row in checkpoint['code_versions_sha256']: assert sha(ROOT/row['path'])==row['sha256']
        full=proof['status']=='SAVED_TX100_DEPTH_AUDIT_PASS'
        assert proof['status'] in ['SAVED_C0_INDEPENDENT_AUDIT_PASS','SAVED_TX100_DEPTH_AUDIT_PASS'] and proof['failed_checks']==0
        assert proof['new_BIE_solves']==0
        if not full: assert proof['fixed_coordinate_max_difference_m']==0
        assert cfg=={k:v for k,v in prior.items() if k!='recovery'}, 'Recovery physics/budgets must match original configuration'
        for row in rows:
            target=dest/'origin'/row['path']; target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(origin/row['path'],target)
        (dest/'origin/retrieval').mkdir(exist_ok=True)
        shutil.copyfile(origin/'retrieval/manifest.json',dest/'origin/retrieval/manifest.json')
        shutil.copytree(audit,dest/'recovery_audit')
        for file in (audit/'audited_work').rglob('*'):
            rel=file.relative_to(audit/'audited_work')
            generated={'independent_reference_resolution.csv','independent_final_metrics_audit.json','monotonicity_diagnostic.json'}
            if file.is_file() and rel.parts[0] not in ['audit','figures'] and rel.as_posix() not in generated:
                target=dest/'imported/depth_0100m'/rel; target.parent.mkdir(parents=True,exist_ok=True); shutil.copyfile(file,target)
        imported=dest/'imported/depth_0100m'
        for name in ['source_geometry_audit.json','source_geometry_audit.mat']:
            shutil.copyfile(origin/'output/depth_0100m'/name,imported/name)
        (imported/'jobs/models').mkdir(exist_ok=True)
        shutil.copyfile(origin/'output/depth_0100m/jobs/models/receipt.json',imported/'jobs/models/receipt.json')
        for name in ['receivers.csv','coefficients.csv','policy.json']:
            shutil.copyfile(origin/'inputs'/name,dest/'inputs'/name)
        session='20261118T133851_1508804' if full else '20261118T090204_1101719'
        old=json.loads((origin/'sessions'/session/'exit.json').read_text())
        runtime=json.loads((origin/'sessions'/session/'runtime.json').read_text())
        cfg['recovery']=dict(kind='IMPORT_AUDITED_TX100_NO_RESOLVE' if full else 'IMPORT_AUDITED_D1_C0_NO_RESOLVE',
            remote_task_name='d1_depth_recovery_20261118_02' if full else 'd1_depth_recovery_20261118_01',
            origin_manifest_sha256=sha(origin/'manifest.json'),remaining_BIE_solves=5 if full else 9,
            common_sha256=sha(imported/'common.mat'),C0_data_sha256=sha(imported/'jobs/bie_coarse/data.mat'),
            prior_supervisor_elapsed_s=old.get('cumulative_elapsed_s',old['elapsed_s']),prior_RAM_baseline_bytes=runtime['initial_available_bytes'])
        if full:
            cfg['recovery']['layer_data_sha256']=[sha(imported/'jobs'/stage/'data.mat') for stage in ['bie_coarse','bie_main','bie_self','bie_space','bie_window']]
            cfg['recovery']['comparison_sha256']=sha(imported/'comparison.mat')
        if attempt_dirs:
            assert full
            records=[]
            for index,attempt_dir in enumerate(attempt_dirs,1):
                attempt=Path(attempt_dir).resolve(); packet=json.loads((attempt/'remote_stdout.json').read_text(encoding='utf-8-sig'))
                assert packet['owned_matlab_still_alive'] is False
                sid=Path(packet['session']).name; key='sessions/'+sid+'/'
                exit_info=json.loads(packet['files'][key+'exit.json']); runtime_info=json.loads(packet['files'][key+'runtime.json'])
                stopped=json.loads(packet['files']['output/checkpoint.json'])
                assert stopped['new_BIE_solves']==0 and len(stopped['completed_steps'])==5
                assert abs(exit_info['cumulative_elapsed_s']-cfg['recovery']['prior_supervisor_elapsed_s']-exit_info['elapsed_s'])<1e-9
                rel='attempt_evidence/attempt_'+str(index); shutil.copytree(attempt,dest/rel); records.append(rel+'/remote_stdout.json')
                cfg['recovery']['prior_supervisor_elapsed_s']+=exit_info['elapsed_s']
                cfg['recovery']['prior_RAM_baseline_bytes']=max(cfg['recovery']['prior_RAM_baseline_bytes'],runtime_info['initial_available_bytes'])
            cfg['recovery']['additional_failed_attempts']=records
            cfg['recovery']['remote_task_name']='d1_depth_recovery_20261118_04'
            cfg['final_figures_mode']='LOCAL_MCP_REQUIRED_AFTER_REMOTE_NUMERICAL'
    write(dest/'config.json',cfg)
    entries=[dict(path=p.relative_to(dest).as_posix(),sha256=sha(p)) for p in sorted(dest.rglob('*')) if p.is_file()]
    scope='D1A_RECOVERY_WITH_AUDITED_TX100' if origin_dir is not None and full else 'D1A_RECOVERY_WITH_AUDITED_C0' if origin_dir is not None else 'D1A_INPUT_CODE_ONLY_NO_DENSITY'
    write(dest/'manifest.json',dict(schema_version=1,scope=scope,files=entries))
    archive=dest/'input_package.tar'
    with tarfile.open(archive,'x') as tf:
        for rel in ['code','inputs','provenance','origin','imported','recovery_audit','attempt_evidence','config.json','manifest.json']:
            if (dest/rel).exists(): tf.add(dest/rel,arcname=rel)
    (dest/'input_package.tar.sha256').write_text(sha(archive)+'  input_package.tar\n',encoding='ascii')
    print(json.dumps(dict(status='D1_RECOVERY_PREPARED_NOT_STARTED' if origin_dir is not None else 'D1_PREPARED_UNRUN',files=len(entries),archive_bytes=archive.stat().st_size,archive_sha256=sha(archive))))
    return dest

if __name__=='__main__': build(*sys.argv[1:])
