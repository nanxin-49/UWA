"""D1 infrastructure only: SHA, frozen configuration, safe stopped-task archives."""
import hashlib
import json
from pathlib import Path, PurePosixPath
import sys
import tarfile
import time


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def verify(root, entries):
    root = Path(root).resolve()
    assert entries, 'Empty manifest'
    seen = set()
    for row in entries:
        rel = row['path']
        p = PurePosixPath(rel)
        assert rel and not p.is_absolute() and '..' not in p.parts
        assert ':' not in rel and '\\' not in rel and p.parts[0] != 'cash' and rel not in seen
        seen.add(rel)
        dest = root.joinpath(*p.parts)
        assert root in dest.resolve().parents and dest.is_file()
        assert not any(q.is_symlink() for q in [dest, *dest.parents] if q != root.parent)
        assert sha(dest) == row['sha256'], f'SHA drift: {rel}'


def preflight(root):
    import csv, math
    root=Path(root).resolve()
    manifest=json.loads((root/'manifest.json').read_text(encoding='utf-8-sig'))
    assert manifest['scope'] in ['D1A_INPUT_CODE_ONLY_NO_DENSITY','D1A_RECOVERY_WITH_AUDITED_C0','D1A_RECOVERY_WITH_AUDITED_TX100']
    verify(root,manifest['files'])
    cfg=json.loads((root/'config.json').read_text(encoding='utf-8-sig'))
    assert cfg['kind']=='D1A_FIXED_PM_DEPTH_COMPARISON'
    assert cfg['depths_m']==[100,1490] and cfg['receiver_z_m']==6 and cfg['receiver_halfwidth_m']==36 and cfg['receivers']==411
    assert cfg['source']==dict(frequency_hz=4000,c_mps=1500,rho_kgpm3=1025,v0_mps=.001,total_width_m=.1)
    assert cfg['PM']==dict(U_PM_mps=6,N=12,seed=260001,period_m=160,plateau_m=42,support_m=50)
    assert cfg['workers']==cfg['threads']==1 and cfg['RAM_reserve_fraction']==.3
    assert cfg['max_wall_s']==5400 and cfg['max_depth_wall_s']==2400
    assert cfg['estimated_working_RAM_GiB']==12 and cfg['minimum_available_RAM_GiB']==18 and cfg['task_RSS_cap_GiB']==12
    assert cfg['time_convention']=='exp(-i omega t)' and cfg['phase_contract_version']=='2.0.0' and cfg['representation']=='absolute_physical'
    assert not any(cfg[k] for k in ['bottom_enabled','downward_radiation','multiple_reflections','structure_scattering','D1b_enabled','automatic_retry'])
    assert cfg['BIE_recipes']==[[60,10,96],[60,12,96],[60,12,64],[60,14,96],[70,12,96]]
    assert cfg['PS_recipes']==[[128,.05,.1875],[256,.05,.1875],[512,.05,.1875],[1024,.05,.1875],[1024,.025,.1875],[1024,.025,.09375]]
    assert cfg['FullFK_recipes']==[[400,50],[800,50],[960,60],[1440,90]]
    assert cfg['gates']==dict(sensitivity=1e-4,linear=1e-10,boundary=1e-4,source=1e-7,PS=1e-4,PS_step=1e-10,FullFK=1e-5,reconstruction=1e-11,resolution_factor=5,phase_mask_fraction=.01)
    def rows(p):
        with p.open(encoding='utf-8-sig',newline='') as f: return list(csv.DictReader(f))
    orig=rows(root/'provenance/k0_coefficients.csv'); coeff=rows(root/'inputs/coefficients.csv')
    cp=json.loads((root/'provenance/k0_checkpoint.json').read_text(encoding='utf-8-sig'))
    assert cp['status']=='AUDITED_PASS' and cp['case']==dict(U_PM_mps=6,seed=260001,N=12)
    declared={r['path']:r['sha256'] for r in cp['artifacts_sha256']}
    assert declared['coefficients.csv']==sha(root/'provenance/k0_coefficients.csv')
    assert len(orig)==len(coeff)==12
    for j,(x,y) in enumerate(zip(orig,coeff)):
        for old,new in [('k_rad_per_m','k_rad_per_m'),('Sk_m3','Sk_m3'),('cos_depth_coeff_m','cos_coeff_m'),('sin_depth_coeff_m','sin_coeff_m')]:
            assert x[old]==y[new], 'Original decimal coefficient token drift'
        assert abs(float(y['k_rad_per_m'])-(j+1)*2*math.pi/160)<1e-14
        assert all(math.isfinite(float(v)) for v in y.values())
    rx=rows(root/'inputs/receivers.csv'); assert len(rx)==411
    assert all(abs(float(r['x_m'])-(-36+72*j/410))<1e-13 and float(r['z_m'])==6 for j,r in enumerate(rx))
    policy=json.loads((root/'inputs/policy.json').read_text(encoding='utf-8-sig'))
    assert policy['version']=='r1-bie-reference-v2' and policy['grid_monotonicity_mode']=='diagnostic'
    assert policy['gates']['sensitivity']==1e-4 and policy['gates']['resolution_factor']==5
    if 'recovery' in cfg:
        r=cfg['recovery']; origin=root/'origin'
        full=r['kind']=='IMPORT_AUDITED_TX100_NO_RESOLVE'
        assert manifest['scope']==('D1A_RECOVERY_WITH_AUDITED_TX100' if full else 'D1A_RECOVERY_WITH_AUDITED_C0')
        if str(root.parent)=='/home/misc/cwx': assert root.name==r['remote_task_name']
        assert r['kind'] in ['IMPORT_AUDITED_D1_C0_NO_RESOLVE','IMPORT_AUDITED_TX100_NO_RESOLVE']
        assert r['remaining_BIE_solves']==(5 if full else 9)
        expected='6221c13050770c38392be06d9675328673585a96000b5ba00e399235441ac727' if full else '74b1a47e6bd59ece8589510b7287a3b3de00dc6bf3774142f99f2081be19d81a'
        assert r['origin_manifest_sha256']==expected
        assert sha(origin/'manifest.json')==r['origin_manifest_sha256']
        origin_cfg,_=preflight(origin)
        assert ('recovery' in origin_cfg)==full
        retrieved=json.loads((origin/'retrieval/manifest.json').read_text(encoding='utf-8-sig'))
        verify(origin,retrieved)
        session='20261118T133851_1508804' if full else '20261118T090204_1101719'
        old=json.loads((origin/'sessions'/session/'exit.json').read_text())
        assert old['matlab_exit_code']==old['supervisor_exit_code']==1 and old['stop_reason'] is None
        previous=old.get('cumulative_elapsed_s',old['elapsed_s'])
        assert abs(previous-(old['elapsed_s']+origin_cfg.get('recovery',{}).get('prior_supervisor_elapsed_s',0)))<1e-9
        runtime=json.loads((origin/'sessions'/session/'runtime.json').read_text())
        baseline=runtime['initial_available_bytes']
        attempts=r.get('additional_failed_attempts',[r['additional_failed_attempt']] if 'additional_failed_attempt' in r else [])
        seen_attempts=set()
        for attempt in attempts:
            assert full and attempt.startswith('attempt_evidence/') and '..' not in Path(attempt).parts
            packet=json.loads((root/attempt).read_text(encoding='utf-8-sig'))
            assert packet['owned_matlab_still_alive'] is False
            sid=Path(packet['session']).name; assert sid in ['20261118T141216_1539908','20261118T141851_1543007'] and sid not in seen_attempts
            seen_attempts.add(sid); prefix='sessions/'+sid+'/'
            for key in ['manifest.json','output/checkpoint.json',prefix+'exit.json',prefix+'runtime.json']:
                assert hashlib.sha256(packet['files'][key].encode('utf-8')).hexdigest()==packet['sha256'][key]
            known={'20261118T141216_1539908':'d2026cf1e66aa995337ed3bbf83f2e98688c8361c1ea78160e01553d51384a5b','20261118T141851_1543007':'89d4a449e9b6096f6d1c2cd53ea25bfbc3da7f3c69274904016bcf06e315a5e5'}
            assert packet['sha256']['manifest.json']==known[sid]
            failed=json.loads(packet['files']['output/checkpoint.json'])
            assert failed['new_BIE_solves']==0 and len(failed['completed_steps'])==5
            assert 'd1:CSVOverwrite' in failed['failed_or_unrun_checks'] or 'Reopened PNG pixel mismatch' in failed['failed_or_unrun_checks']
            extra=json.loads(packet['files'][prefix+'exit.json']); extra_runtime=json.loads(packet['files'][prefix+'runtime.json'])
            assert extra['matlab_exit_code']==extra['supervisor_exit_code']==1 and extra['stop_reason'] is None
            assert abs(extra['cumulative_elapsed_s']-previous-extra['elapsed_s'])<1e-9
            previous+=extra['elapsed_s']; baseline=max(baseline,extra_runtime['initial_available_bytes'])
        assert r['prior_supervisor_elapsed_s']==previous and r['prior_RAM_baseline_bytes']==baseline
        if 'final_figures_mode' in cfg: assert cfg['final_figures_mode']=='LOCAL_MCP_REQUIRED_AFTER_REMOTE_NUMERICAL' and full
        assert sha(root/'imported/depth_0100m/common.mat')==r['common_sha256']==sha(origin/'output/depth_0100m/common.mat')
        assert sha(root/'imported/depth_0100m/jobs/bie_coarse/data.mat')==r['C0_data_sha256']==sha(origin/'output/depth_0100m/jobs/bie_coarse/data.mat')
        proof=json.loads((root/'recovery_audit/receipt.json').read_text())
        assert proof['status']==('SAVED_TX100_DEPTH_AUDIT_PASS' if full else 'SAVED_C0_INDEPENDENT_AUDIT_PASS')
        assert proof['failed_checks']==0 and proof['new_BIE_solves']==0
        assert proof['origin_manifest_sha256']==r['origin_manifest_sha256']
        if full:
            assert proof['completed_BIE_layers']==5 and proof['figure_audit']['differing_pixel_channels']==0
            assert proof['figure_audit']['reopened_curves']==24
            ids=['bie_coarse','bie_main','bie_self','bie_space','bie_window']
            for j,stage in enumerate(ids):
                rel='output/depth_0100m/jobs/'+stage
                assert sha(root/'imported/depth_0100m/jobs'/stage/'data.mat')==r['layer_data_sha256'][j]==sha(origin/rel/'data.mat')
                seal=json.loads((origin/rel/'checkpoint.json').read_text()); assert seal['status']=='AUDITED_PASS'; verify(origin/rel,seal['artifacts_sha256'])
            assert sha(root/'imported/depth_0100m/comparison.mat')==r['comparison_sha256']==sha(origin/'output/depth_0100m/comparison.mat')
        else:
            assert proof['fixed_coordinate_max_difference_m']==0
        assert (root/'inputs/receivers.csv').read_bytes()==(origin/'inputs/receivers.csv').read_bytes()
    else:
        assert manifest['scope']=='D1A_INPUT_CODE_ONLY_NO_DENSITY'
    return cfg,manifest


def archive(root):
    root = Path(root).resolve()
    preflight(root)
    retrieval = root / 'retrieval'
    retrieval.mkdir(exist_ok=True)
    stem = 'result_' + time.strftime('%Y%m%d_%H%M%S')
    destination = retrieval / (stem + '.tar')
    files = []
    for prefix in ['code', 'inputs', 'provenance', 'sessions', 'output', 'origin', 'imported', 'recovery_audit','attempt_evidence']:
        folder = root / prefix
        if folder.exists():
            assert not folder.is_symlink()
            for p in folder.rglob('*'):
                assert not p.is_symlink()
                if p.is_file():
                    files.append(p)
    files += [root / 'config.json', root / 'manifest.json']
    for name in ['launch_started.json', 'latest_session.txt']:
        if (root / name).is_file():
            files.append(root / name)
    rows = [dict(path=p.relative_to(root).as_posix(), sha256=sha(p)) for p in sorted(files)]
    index = retrieval / (stem + '.json')
    with index.open('x', encoding='utf-8') as f:
        json.dump(rows, f, indent=2)
    with tarfile.open(destination, 'x') as tf:
        for p in files:
            tf.add(p, arcname=p.relative_to(root).as_posix(), recursive=False)
        tf.add(index, arcname='retrieval/manifest.json')
    return destination, sha(destination)


def extract_verified(archive_path, destination):
    destination = Path(destination).resolve()
    assert not destination.exists(), 'New extraction directory required'
    with tarfile.open(archive_path) as tf:
        members = tf.getmembers()
        seen = set()
        for m in members:
            p = PurePosixPath(m.name)
            assert not p.is_absolute() and '..' not in p.parts and ':' not in m.name and '\\' not in m.name
            assert p.parts[0] != 'cash' and m.name not in seen
            assert m.isfile() or m.isdir(), 'Links/devices prohibited'
            seen.add(m.name)
        destination.mkdir(parents=True)
        tf.extractall(destination, members=members, filter='data')
    result_index = destination / 'retrieval/manifest.json'
    if result_index.exists():
        verify(destination, json.loads(result_index.read_text(encoding='utf-8-sig')))
    preflight(destination)


if __name__ == '__main__':
    command, root = sys.argv[1:3]
    if command == 'verify':
        cfg, manifest = preflight(root)
        print(json.dumps(dict(status='PREFLIGHT_PASS', files=len(manifest['files']), depths=cfg['depths_m'], field_calls=0)))
    elif command == 'archive':
        path, digest = archive(root)
        print('ARCHIVE_PATH=' + str(path))
        print('ARCHIVE_SHA256=' + digest)
    elif command == 'extract':
        extract_verified(root, sys.argv[3])
        print('EXTRACTION_SHA_VERIFIED=' + sys.argv[3])
    else:
        raise ValueError(command)
