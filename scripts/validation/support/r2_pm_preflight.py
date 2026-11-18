"""Read-only R2 packaging guard, standard library only; never starts MATLAB/SSH."""
import csv
import hashlib
import json
import math
import re
import sys
from pathlib import Path

sys.dont_write_bytecode = True


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def safe(root, rel):
    assert rel and not re.search(r'(^/|:|\\|(^|/)\.\.(/|$))', rel)
    p = (root / rel).resolve()
    assert p.is_relative_to(root.resolve())
    return p


def check(root):
    root = Path(root).resolve()
    m = read(root / 'manifest.json')
    assert m['schema_version'] == 1 and m['scope'] == 'R2_PM_TWO_SEED_REVIEW'
    assert m['seed'] in (260002, 260003)
    assert re.fullmatch(r'r2_pm_seed' + str(m['seed']) + r'_[A-Za-z0-9_-]+', m['remote_task_name'])
    assert m['single_process'] and m['single_thread'] and not m['parpool'] and not m['automatic_retry']
    assert m['RAM_reserve_fraction'] == .3 and not m['cumulative_time_cap'] and not m['cost_prediction_gate']
    assert m['max_parallel_seed_processes'] == 2 and m['estimated_working_bytes_per_seed'] == 12*1024**3
    assert m['aggregate_admission_bytes'] == 24*1024**3
    paths = [f['path'] for f in m['files']]
    assert len(paths) == len(set(paths)) and len(paths) >= 30
    for f in m['files']:
        assert re.fullmatch('[a-f0-9]{64}', f['sha256'])
        assert sha(safe(root, f['path'])) == f['sha256'], f['path']
    assert m['policy_sha256'] == sha(root / 'policy.json')
    policy = read(root / 'policy.json')
    assert policy['version'] == 'r1-bie-reference-v2' and policy['grid_monotonicity_mode'] == 'diagnostic'
    assert policy['reference_stage'] == 3
    assert policy['pair_names'] == ['grid_10_12', 'grid_12_14', 'self_64_96', 'window_60_70']
    assert policy['gates'] == dict(sensitivity=1e-4, linear=1e-10, boundary=1e-4, reconstruction=1e-11,
                                 geometry=1e-12, source=1e-7, PS=1e-4, PS_step=1e-10, FullFK=1e-5,
                                 resolution_factor=5, phase_mask_fraction=.01)
    cfg = read(root / 'input/frozen_config.json')
    assert cfg['seed'] == m['seed']
    assert cfg['source'] == dict(frequency_hz=4000, c_mps=1500, rho_kgpm3=1025, v0_mps=.001, total_width_m=.1, z_source_m=100)
    for key, value in dict(period_m=160, modes=12, C2_plateau_halfwidth_m=42, C2_support_halfwidth_m=50,
                           receivers=411, receiver_z_m=3, phase_mask_fraction_of_BIE_RMS=.01).items():
        assert cfg[key] == value
    assert cfg['phase_contract_version'] == '2.0.0'
    assert cfg['BIE_recipes'] == [[60, 10, 96], [60, 12, 96], [60, 12, 64], [60, 14, 96], [70, 12, 96]]
    assert cfg['PS_recipes'] == [[128,.05,.1875],[256,.05,.1875],[512,.05,.1875],[1024,.05,.1875],[1024,.025,.1875],[1024,.025,.09375]]
    assert cfg['FullFK_recipes'] == [[400,50],[800,50],[960,60],[1440,90]]
    def rows(name):
        with (root / 'input' / name).open(encoding='utf-8-sig', newline='') as stream:
            return list(csv.DictReader(stream))
    p, anchor, receivers = rows('coefficients.csv'), rows('r1_coefficients.csv'), rows('receivers.csv')
    assert len(p) == len(anchor) == 12 and len(receivers) == 411
    for i, (a, b) in enumerate(zip(p, anchor), 1):
        assert a['k_rad_per_m'] == b['k_rad_per_m'] and a['Sk_m3'] == b['Sk_m3']
        assert abs(float(a['k_rad_per_m']) - i*2*math.pi/160) < 1e-14
        assert all(math.isfinite(float(v)) for v in a.values())
    assert p != anchor
    for i, r in enumerate(receivers):
        assert float(r['x_m']) == (i-205)*.1953125 and float(r['z_m']) == 3
    src = read(root / 'provenance/sources.json')
    assert src['coefficient_sha256'] == src['coefficient_parent_sha256'] == sha(root / 'input/coefficients.csv')
    assert src['policy_sha256'] == sha(root / 'policy.json') and not src['density_reuse'] and not src['R1_results_modified']
    if (root / 'package.sha256').exists():
        for row in (root / 'package.sha256').read_text().splitlines():
            digest, rel = row.split('  ', 1)
            assert sha(safe(root, rel)) == digest, rel
    return dict(status='PASS', seed=m['seed'], SHA_bound_files=len(paths), failed_checks=0,
                new_acoustic_solves=0, SSH_connections=0, mode='static_input_verification_only')


if __name__ == '__main__':
    print(json.dumps(check(sys.argv[1]), ensure_ascii=False))
