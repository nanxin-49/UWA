"""K1 file/CSV/checkpoint integrity checks; no MATLAB or acoustic calls."""
import csv
import hashlib
import json
import os
from pathlib import Path
import sys
import tarfile
import time

MATRIX=[('U6_N07',6,7),('U6_N40',6,40),('U5_N12',5,12),('U8_N12',8,12),('U5_N40',5,40),('U8_N40',8,40)]
SOURCE=dict(frequency_hz=4000,c_mps=1500,rho_kgpm3=1025,v0_mps=.001,total_width_m=.1,z_source_m=100)


def read(p):
    return json.loads(Path(p).read_text(encoding='utf-8-sig'))


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):
            h.update(block)
    return h.hexdigest()


def bound(root, relative):
    rel=Path(relative)
    assert not rel.is_absolute() and '..' not in rel.parts and ':' not in relative and '\\' not in relative
    assert rel.parts and rel.parts[0]!='cash'
    p=(root/rel).resolve()
    assert root.resolve() in p.parents and not (root/rel).is_symlink()
    return p


def verify_rows(root, rows):
    assert rows
    seen=set()
    for row in rows:
        rel=row['path']; assert rel not in seen; seen.add(rel)
        p=bound(root,rel); assert p.is_file() and sha(p)==row['sha256'], 'SHA drift: '+rel
    return len(rows)


def csv_rows(p):
    with Path(p).open(encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f))


def verify(root):
    root=Path(root).resolve(); m=read(root/'manifest.json')
    assert m['scope']=='K1_PM_MASTER'
    count=verify_rows(root,m['files']); cfg=read(root/'config.json')
    assert cfg['kind']=='K1_PM_SIX_NEW_SURFACES'
    assert [(c['case_id'],c['U_PM_mps'],c['modes']) for c in cfg['cases']]==MATRIX
    assert cfg['max_processes']==1 and cfg['threads']==1 and not cfg['parpool'] and not cfg['automatic_retry']
    assert cfg['RAM_reserve_fraction']==.3 and cfg['case_working_bytes']==12*1024**3 and cfg['task_RSS_cap_bytes']==12*1024**3
    assert cfg['max_case_wall_s']==1800 and cfg['max_total_wall_s']==10800
    assert cfg['admission_memory_PSI_full_avg10_max']==5 and cfg['admission_io_PSI_full_avg10_max']==10
    proof=read(root/'anchor/anchor_proof.json')
    assert proof['status']=='STRICT_R1_K0_ANCHOR_MATCH_PASS' and proof['failed_checks']==0
    assert sha(root/'anchor/anchor_proof.json')==m['anchor_proof_sha256']
    assert proof['R1_formal_qualification_sha256']=='eb6d86a8e845ee2785efab4b93470f91d26026ed0d415bdbd808023152b43265'
    assert sha(root/'anchor/R1_formal_review.json')==proof['R1_formal_qualification_sha256']
    assert sha(root/'anchor/R1_final_seal.json')==proof['R1_final_seal_sha256']
    for name,u,n in MATRIX:
        case=root/'cases'/name; cm=read(case/'manifest.json'); verify_rows(case,cm['files'])
        assert cm['scope']=='K1_PM_CASE' and cm['case_id']==name and cm['N']==n and cm['U_PM_mps']==u
        assert cm['seed']==260001 and cm['threads']==1 and cm['max_case_processes']==1
        assert sha(case/'input/r1_anchor_proof.json')==m['anchor_proof_sha256']==cm['anchor_proof_sha256']
        f=read(case/'input/frozen_config.json')
        assert f['source']==SOURCE and f['seed']==260001 and f['case_id']==name and f['U_PM_mps']==u and f['modes']==n
        assert f['wind_label']=='U_PM' and f['family']=='legacy12_gaussian_extension_v1'
        assert f['period_m']==160 and f['C2_plateau_halfwidth_m']==42 and f['C2_support_halfwidth_m']==50
        assert f['phase_contract_version']=='2.0.0' and f['phase_mask_fraction_of_BIE_RMS']==.01
        assert f['BIE_recipes']==[[60,10,96],[60,12,96],[60,12,64],[60,14,96],[70,12,96]]
        assert f['PS_recipes']==[[128,.05,.1875],[256,.05,.1875],[512,.05,.1875],[1024,.05,.1875],[1024,.025,.1875],[1024,.025,.09375]]
        assert f['FullFK_recipes']==[[400,50],[800,50],[960,60],[1440,90]]
        coeff=csv_rows(case/'input/coefficients.csv'); original=csv_rows(case/'input/k0_coefficients.csv')
        assert len(coeff)==len(original)==n
        for a,b in zip(coeff,original):
            for x,y in [('k_rad_per_m','k_rad_per_m'),('Sk_m3','Sk_m3'),('cos_coeff_m','cos_depth_coeff_m'),('sin_coeff_m','sin_depth_coeff_m')]:
                assert float(a[x])==float(b[y]), 'K0 coefficient drift'
        rx=csv_rows(case/'input/receivers.csv'); assert len(rx)==411
        assert all(float(r['x_m'])==(j-205)*.1953125 and float(r['z_m'])==3 for j,r in enumerate(rx))
        k0=read(case/'input/k0_checkpoint.json')
        assert k0['status']=='AUDITED_PASS' and k0['manifest_sha256']=='cf48069437260281d34a560ba6c3bcee1047d61acdc18154e8c1d68aeb6fd6c6'
        names={r['path']:r['sha256'] for r in k0['artifacts_sha256']}
        for old,new in [('coefficients.csv','k0_coefficients.csv'),('spectra.csv','k0_spectra.csv'),('surface.csv','k0_surface.csv'),('metrics.json','k0_metrics.json')]:
            assert sha(case/'input'/new)==names[old]
    return dict(status='K1_INPUT_PREFLIGHT_PASS',files=count,cases=6,anchor='STRICT_MATCH',acoustic_calls=0,SSH_calls=0)


def seal(root, name, signature):
    s=read(root/name); assert s['schema_version']==1 and s['signature']==signature and s['status']=='AUDITED_PASS'
    verify_rows(root,s['files']); return s


def completed_case(case,require_checkpoint=True):
    case=Path(case).resolve(); m=read(case/'manifest.json'); fp=sha(case/'manifest.json')
    verify_rows(case,m['files'])
    model=seal(case/'model_stage','committed.json',fp+':MODELS:'+fp)
    assert {'common.mat','jobs/models/receipt.json','solved.json','solve_invocation.json'} <= {x['path'] for x in model['files']}
    model_solved=read(case/'model_stage/solved.json')
    assert model_solved['status']=='SOLVED_UNAUDITED' and model_solved['signature']==fp+':MODELS:'+fp
    verify_rows(case/'model_stage',model_solved['files'])
    model_receipt=read(case/'model_stage/jobs/models/receipt.json')
    assert model_receipt['status']=='COMPLETE' and model_receipt['failed_checks']==0
    marker=read(case/'model_stage/solve_invocation.json'); assert marker['invocations']==1 and marker['stage']=='MODELS' and marker['manifest_sha256']==fp
    pred=sha(case/'model_stage/committed.json')
    for j,name in enumerate(['bie_coarse','bie_main','bie_self','bie_space','bie_window']):
        dest=case/'jobs'/name; signed=seal(dest,'committed.json',pred+':C'+str(j)+':'+fp)
        assert {'data.mat','points.csv','config.json','timings.json','solved.json','solve_invocation.json','audit.json','audit_reconstruction.mat','gates.json'} <= {x['path'] for x in signed['files']}
        audit=read(dest/'audit.json'); gates=read(dest/'gates.json'); marker=read(dest/'solve_invocation.json')
        solved=read(dest/'solved.json')
        assert solved['status']=='SOLVED_UNAUDITED' and solved['signature']==pred+':C'+str(j)+':'+fp
        verify_rows(dest,solved['files'])
        assert audit['status']=='PASS' and audit['failed_checks']==0 and gates['status']=='PASS'
        assert all(audit['checks'].values()) and set(audit['checks'])=={'density_receiver','MAT_CSV','component_closure','incident_rhs','all_midpoints','independent_boundary'}
        assert audit['platform_midpoints']==[340,408,408,476,476][j] and audit['independent_boundary_points']==19
        assert marker['invocations']==1 and marker['stage']=='C'+str(j) and marker['manifest_sha256']==fp
        pred=sha(dest/'committed.json')
    seal(case/'finalization','committed.json',pred+':FINAL:'+fp)
    pointer=read(case/'finalization/pointer.json'); work=bound(case/'finalization',pointer['relative_work_dir'])
    if (case/'checkpoint.json').is_file():
        cp=read(case/'checkpoint.json'); assert cp['manifest_sha256']==fp and cp['case_id']==m['case_id']
    else:
        cp=None
    if require_checkpoint:
        assert cp is not None and cp['status']=='K1_CASE_COMPLETE' and cp['failed_checks']==0
        assert (case/cp['final_work_dir']).resolve()==work and cp['completed_steps']==['MODELS','C0','C1','C2','C3','C4']
    for path in ['audit/receipt.json','independent_final_metrics_audit.json','figure_audit/receipt.json']:
        r=read(work/path); assert r['failed_checks']==0
    final=read(work/'independent_final_metrics_audit.json')
    assert len(final['reference_sensitivities'])==4 and all(x<1e-4 for x in final['reference_sensitivities']) and final['policy_version']=='r1-bie-reference-v2'
    figures=read(work/'figure_audit/receipt.json')
    assert figures['checked_error_curves']==12 and figures['pixel_tolerance']==0
    assert all(r['pixel_differences']==[0,0] for r in figures['records'].values())
    return dict(case_id=m['case_id'],work=str(work),BIE_invocations=5,status='AUDITED_CASE_PASS')


def ensure_stopped(root):
    assert not (root/'control/pending_spawn.json').exists(), 'Uncertain spawn; inspect before Fetch'
    path=root/'control/active_child.json'
    if path.is_file():
        identity=read(path); assert identity['root']==str(root.resolve())
        stat=Path(f"/proc/{identity['pid']}/stat")
        if stat.is_file():
            ticks=stat.read_text().rsplit(') ',1)[1].split()[19]
            assert not (ticks==identity['start_ticks'] and os.getpgid(identity['pid'])==identity['process_group']==identity['pid']), 'Owned child still active; Fetch refused'


def archive(root):
    ensure_stopped(root)
    folder=root/'retrieval'; folder.mkdir(exist_ok=True)
    stem='result_'+time.strftime('%Y%m%d_%H%M%S')
    target=folder/(stem+'.tar'); assert not target.exists()
    files=[]
    for name in ['code','cases','anchor','provenance','resource_evidence','sessions','control','summary']:
        if (root/name).exists():
            files += [p for p in (root/name).rglob('*') if p.is_file() and not any(x in ['prefs','__pycache__'] for x in p.relative_to(root).parts)]
    files += [root/'config.json',root/'manifest.json']
    for name in ['checkpoint.json','latest_session.txt']:
        if (root/name).is_file(): files.append(root/name)
    rows=[]
    for p in sorted(files):
        assert not p.is_symlink(); rows.append(dict(path=p.relative_to(root).as_posix(),sha256=sha(p)))
    mf=folder/(stem+'.json'); mf.write_text(json.dumps(rows,indent=2),encoding='utf-8')
    with tarfile.open(target,'x') as tar:
        for p in files: tar.add(p,arcname=p.relative_to(root).as_posix(),recursive=False)
        tar.add(mf,arcname='retrieval/manifest.json',recursive=False)
    print('ARCHIVE_PATH='+str(target)); print('ARCHIVE_SHA256='+sha(target))


def safe_extract(archive_path,destination):
    with tarfile.open(archive_path,'r') as tar:
        for member in tar.getmembers():
            p=Path(member.name)
            assert not p.is_absolute() and '..' not in p.parts and ':' not in member.name and '\\' not in member.name
            assert p.parts and p.parts[0]!='cash'
            assert member.isfile() or member.isdir(), 'Nonregular archive member'
        tar.extractall(destination)
    rows=read(Path(destination)/'retrieval/manifest.json'); verify_rows(Path(destination).resolve(),rows)
    return dict(status='FETCH_INTERNAL_SHA_PASS',files=len(rows))


if __name__=='__main__':
    root=Path(sys.argv[1]).resolve(); command=sys.argv[2] if len(sys.argv)>2 else 'verify'
    if command=='verify': print(json.dumps(verify(root)))
    elif command=='archive': verify(root); archive(root)
    elif command=='case': print(json.dumps(completed_case(root)))
    elif command=='extract': print(json.dumps(safe_extract(sys.argv[3],root)))
    else: raise ValueError('Unknown preflight action')
