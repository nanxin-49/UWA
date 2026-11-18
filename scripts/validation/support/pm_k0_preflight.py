"""K0 allowlist verifier / stopped-task archive; Python standard library only."""
import hashlib
import json
from pathlib import Path
import sys
import tarfile
import time


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def verify(root, entries):
    assert entries, 'Empty manifest'
    seen = set()
    for entry in entries:
        rel = entry['path']
        p = Path(rel)
        assert not p.is_absolute() and '..' not in p.parts and ':' not in rel and '\\' not in rel
        assert p.parts[0] != 'cash' and rel not in seen
        seen.add(rel)
        target = (root / p).resolve()
        assert root in target.parents and target.is_file() and not (root / p).is_symlink()
        assert sha(target) == entry['sha256'], f'SHA drift: {rel}'


def main():
    root = Path(sys.argv[1]).resolve()
    command = sys.argv[2] if len(sys.argv) > 2 else 'verify'
    manifest = json.loads((root / 'manifest.json').read_text(encoding='utf-8-sig'))
    verify(root, manifest['files'])
    cfg = json.loads((root / 'config.json').read_text(encoding='utf-8-sig'))
    assert cfg['kind'] == 'K0_SURFACE_ONLY' and cfg['wind_label'] == 'U_PM'
    assert cfg['N'] == [7, 12, 19, 32, 40] and cfg['U_PM_mps'] == [5, 6, 8]
    assert cfg['seeds'] == [260001, 260002, 260003] and cfg['L_m'] == 160
    expected = [dict(seed=s, U_PM_mps=u, N=n) for s in cfg['seeds'] for u in cfg['U_PM_mps'] for n in cfg['N']]
    assert cfg['cases'] == expected and cfg['sample_count'] == 8192
    assert cfg['max_wall_s'] == 1800 and cfg['RAM_reserve_fraction'] == .30
    if command == 'archive':
        # Caller holds the same nonblocking workflow lock as the supervisor.
        folder = root / 'retrieval'
        folder.mkdir(exist_ok=True)
        name = 'result_' + time.strftime('%Y%m%d_%H%M%S')
        archive = folder / (name + '.tar')
        assert not archive.exists()
        files = []
        for prefix in ['code', 'inputs', 'provenance', 'sessions', 'output']:
            if (root / prefix).exists():
                files += [p for p in (root / prefix).rglob('*') if p.is_file()]
        files += [root / 'config.json', root / 'manifest.json']
        for name in ['ram_baseline.json', 'latest_session.txt']:
            if (root / name).is_file():
                files.append(root / name)
        entries = []
        for p in sorted(files):
            assert not p.is_symlink()
            entries.append(dict(path=p.relative_to(root).as_posix(), sha256=sha(p)))
        (folder / (name + '.json')).write_text(json.dumps(entries, indent=2), encoding='utf-8')
        with tarfile.open(archive, 'x') as tar:
            for p in files:
                tar.add(p, arcname=p.relative_to(root).as_posix(), recursive=False)
            tar.add(folder / (name + '.json'), arcname='retrieval/manifest.json')
        print('ARCHIVE_PATH=' + str(archive))
        print('ARCHIVE_SHA256=' + sha(archive))
    else:
        assert command == 'verify'
        print(json.dumps(dict(status='PREFLIGHT_PASS', files=len(manifest['files']), cases=45, network_calls=0)))


if __name__ == '__main__':
    main()
