"""D0 infrastructure only: SHA, frozen configuration, safe stopped-task archives."""
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
    root = Path(root).resolve()
    manifest = json.loads((root / 'manifest.json').read_text(encoding='utf-8-sig'))
    verify(root, manifest['files'])
    cfg = json.loads((root / 'config.json').read_text(encoding='utf-8-sig'))
    assert cfg['kind'] == 'D0_INCIDENT_ONLY'
    assert cfg['depths_m'] == [100, 300, 600, 1000, 1490] and cfg['fine_depths_m'] == [100, 1490]
    assert cfg['source'] == dict(frequency_hz=4000, c_mps=1500, rho_kgpm3=1025, v0_mps=.001, total_width_m=.1)
    assert cfg['dx_m'] == .05 and cfg['fine_dx_m'] == .025 and cfg['x_halfwidth_m'] == 220
    assert cfg['water_depth_context_m'] == 1500 and cfg['time_convention'] == 'exp(-i omega t)'
    assert cfg['x_source_m'] == 0
    assert cfg['source_model'] == '2D infinite uniform in-phase strip, ideal rigid rear baffle, upward halfspace only'
    assert cfg['medium_model'] == 'uniform lossless water; no SSP, absorption or current'
    assert cfg['phase_contract_version'] == '2.0.0' and cfg['representation'] == 'absolute_physical'
    assert cfg['workers'] == cfg['threads'] == 1 and cfg['RAM_reserve_fraction'] == .30
    assert cfg['max_wall_s'] == 3600 and cfg['max_depth_wall_s'] == 900
    assert cfg['minimum_available_RAM_GiB'] == 4 and cfg['estimated_working_RAM_GiB'] == 2
    assert not cfg['reflection_enabled'] and not cfg['bottom_enabled']
    assert cfg['PM'] == dict(U_PM_mps=6, N=12, seed=260001, L_m=160, C2_plateau_m=42, C2_zero_m=50,
                            geometry='original saved embedded samples abs(x)<=50; z_actual=-eta; no extension')
    assert cfg['gates'] == dict(source_relative=1e-7, gradient_relative=1e-6, phase_absolute_Pa=1e-10,
                               csv_relative=1e-12, fft_relative=1e-12, dft_relative=1e-10,
                               common_point_relative=1e-12, sampling_relative=1e-4,
                               geometry_absolute_m=1e-11, max_phase_step_rad=3.141592653589793/4,
                               phase_mask_fraction=1e-4, spectrum_phase_mask_fraction=1e-6)
    assert cfg['audit_x_m'] == [-220, -200, -150, -100, -50, -20, 0, 20, 50, 100, 150, 200, 220]
    case = root / 'inputs/k0'
    cp = json.loads((case / 'checkpoint.json').read_text(encoding='utf-8-sig'))
    assert cp['status'] == 'AUDITED_PASS' and cp['case'] == dict(U_PM_mps=6, seed=260001, N=12)
    verify(case, cp['artifacts_sha256'])
    return cfg, manifest


def archive(root):
    root = Path(root).resolve()
    preflight(root)
    retrieval = root / 'retrieval'
    retrieval.mkdir(exist_ok=True)
    stem = 'result_' + time.strftime('%Y%m%d_%H%M%S')
    destination = retrieval / (stem + '.tar')
    files = []
    for prefix in ['code', 'inputs', 'provenance', 'sessions', 'output']:
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
