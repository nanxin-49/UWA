"""Read-only v2 input/readiness checks for the launcher and Windows packager.

No MAT loading, MATLAB execution, subprocesses, or network access.
"""
import hashlib
import json
import re
import sys
from pathlib import Path


def sha(path):
    with path.open('rb') as stream:
        digest = hashlib.sha256()
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
        return digest.hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def artifact(root, relative):
    if not isinstance(relative, str) or re.search(r'(^/|:|(^|[\\/])\.\.([\\/]|$))', relative):
        raise ValueError('unsafe artifact path')
    path = (root / relative).resolve()
    if not path.is_relative_to(root.resolve()) or path == root.resolve():
        raise ValueError('artifact escapes root')
    return path


def files(root, entries, required=()):
    names = [e['path'] for e in entries]
    assert entries and len(names) == len(set(names)) and set(required) <= set(names)
    for entry in entries:
        assert sha(artifact(root, entry['path'])) == entry['sha256'], entry['path']


def seal(root, name, signature, required, schema=2):
    record = read(root / name)
    assert record['schema_version'] == schema and record['signature'] == signature
    files(root, record['files'], required)
    return record


def policy_check(p):
    assert p['schema_version'] == 2 and p['version'] == 'r1-bie-reference-v2'
    assert p['grid_monotonicity_mode'] == 'diagnostic' and p['reference_stage'] == 3
    assert p['pair_names'] == ['grid_10_12', 'grid_12_14', 'self_64_96', 'window_60_70']
    assert p['gates'] == dict(sensitivity=1e-4, linear=1e-10, boundary=1e-4,
        reconstruction=1e-11, geometry=1e-12, source=1e-7, PS=1e-4,
        PS_step=1e-10, FullFK=1e-5, resolution_factor=5, phase_mask_fraction=.01)


def inputs(root):
    m = read(root / 'manifest.json')
    assert m['schema_version'] == 2 and m['scope'] == 'R1_BIE_C4_CONTINUATION'
    assert m['allowed_new_solve_stage'] == 4 and m['single_process'] and m['single_thread']
    assert not m['parpool'] and not m['automatic_retry'] and m['RAM_reserve_fraction'] == .3
    files(root, m['files'], ['policy.json', 'approval.json', 'migration_sources.json', 'common.mat'])
    policy_check(read(root / 'policy.json'))
    assert sha(root / 'policy.json') == m['policy_sha256']
    approval = read(root / 'approval.json')
    assert approval['implementation_approved'] is True and approval['C4_started_by_agent'] is False
    assert approval['policy_version'] == 'r1-bie-reference-v2'
    origin = root / 'origin_evidence'
    sources = read(root / 'migration_sources.json')
    assert sha(origin / 'manifest.json') == sources['legacy_manifest_sha256']
    listing = artifact(origin, sources['retrieval_manifest_path'])
    assert sha(listing) == sources['retrieval_manifest_sha256']
    entries = []
    for line in listing.read_text().splitlines():
        match = re.fullmatch(r'([0-9a-f]{64})  (.+)', line)
        assert match
        entries.append(dict(sha256=match[1], path=match[2]))
    assert len(entries) == 110
    files(origin, entries)
    bound = {str(artifact(origin, e['path'])) for e in entries}
    legacy = read(origin / 'manifest.json')
    assert legacy['schema_version'] == 1 and legacy['scope'] == 'R1_BIE_C1_C4_MANUAL'
    assert len(legacy['files']) == 36
    files(origin, legacy['files'])
    pred = sources['legacy_manifest_sha256']
    ids = ['bie_coarse', 'bie_main', 'bie_self', 'bie_space']
    payload = ['data.mat', 'points.csv', 'config.json', 'timings.json']
    for i, name in enumerate(ids):
        dest = origin / 'jobs' / name
        signature = f'{pred}:C{i}:{sources["legacy_manifest_sha256"]}'
        s = seal(dest, 'solved.json', signature, payload, 1)
        assert s['status'] == 'SOLVED_UNAUDITED' and len(s['files']) == 4
        audit = read(dest / 'audit.json')
        assert audit['status'] == 'PASS' and audit['failed_checks'] == 0
        assert len(audit['checks']) == 6 and all(v is True for v in audit['checks'].values())
        required = payload + ['audit.json', 'audit_reconstruction.mat', 'gates.json']
        assert all(str((dest / f).resolve()) in bound for f in required + ['solved.json'])
        gates = read(dest / 'gates.json')
        if i < 3:
            c = seal(dest, 'committed.json', signature, required, 1)
            assert c['status'] == 'AUDITED_PASS' and gates['status'] == 'PASS'
            pred = sha(dest / 'committed.json')
        else:
            assert not (dest / 'committed.json').exists()
            assert gates['status'] == 'FAIL' and gates['grid_nonmonotonic'] is True
            assert len(gates['pairs']) == 3 and all(p['relative'] < 1e-4 for p in gates['pairs'])
    cp = read(origin / 'checkpoint.json')
    assert cp['error_identifier'] == 'r1manual:Nonmonotonic' and cp['new_acoustic_solves'] == 3
    assert not (origin / 'jobs/bie_window/data.mat').exists()
    return m


def ready(root):
    inputs(root)
    fp, ph = sha(root / 'manifest.json'), sha(root / 'policy.json')
    pred = fp
    payload = ['data.mat', 'points.csv', 'config.json', 'timings.json']
    for i, name in enumerate(['bie_coarse', 'bie_main', 'bie_self', 'bie_space']):
        dest = root / 'jobs' / name
        signature = f'{pred}:C{i}:{fp}:{ph}'
        record = seal(dest, 'migration_commit.json', signature,
            payload + ['audit.json', 'audit_reconstruction.mat', 'gates_v2.json', 'effective_config.json'])
        expected = 'IMPORTED_AUDITED_PASS' if i < 3 else 'IMPORTED_AUDITED_PASS_WITH_NONMONOTONIC_DIAGNOSTIC'
        assert record['status'] == expected
        for f in payload:
            assert sha(dest / f) == sha(root / 'origin_evidence/jobs' / name / f)
        a, g = read(dest / 'audit.json'), read(dest / 'gates_v2.json')
        assert a['status'] == 'PASS' and a['failed_checks'] == 0
        assert g['status'] == 'PASS' and g['policy_sha256'] == ph
        pred = sha(dest / 'migration_commit.json')
    complete = seal(root, 'migration_complete.json', f'{pred}:MIGRATED:{fp}:{ph}',
        ['migration_checkpoint.json', 'migration_audits/C3/audit.json', 'migration_audits/C3/audit_reconstruction.mat'])
    assert complete['status'] == 'MIGRATION_AUDITED_PASS'
    v = read(root / 'local_validation.json')
    assert v['status'] == 'PASS' and v['failed_checks'] == 0 and v['new_acoustic_solves'] == 0
    assert v['manifest_sha256'] == fp and v['policy_sha256'] == ph
    validation = seal(root, 'local_validation_commit.json', f'{pred}:VALIDATED:{fp}:{ph}',
        ['local_validation.json', 'preparation_validation/receipt.json', 'preparation_validation/code_analyzer.json'])
    assert validation['status'] == 'LOCAL_ACCEPTANCE_PASS'
    initial = root / 'control/initial_available_bytes.txt'
    assert sha(initial) == sha(root / 'origin_evidence/control/initial_available_bytes.txt')
    return v


if __name__ == '__main__':
    mode, directory = sys.argv[1:3]
    root = Path(directory).resolve()
    if mode == 'inputs':
        inputs(root)
    elif mode == 'ready':
        ready(root)
    else:
        raise ValueError('use inputs or ready')
    print(json.dumps(dict(status='PASS', mode=mode, MATLAB_starts=0, new_acoustic_solves=0)))
