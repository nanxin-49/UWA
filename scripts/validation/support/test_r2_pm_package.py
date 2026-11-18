"""Local R2 static/rejection tests. No acoustic solve, MATLAB or network."""
import csv
import importlib.util
import json
import math
import shutil
import sys
from pathlib import Path

sys.dont_write_bytecode = True


def main():
    root, out = (Path(s).resolve() for s in sys.argv[1:3])
    out.mkdir(parents=True, exist_ok=False)
    spec = importlib.util.spec_from_file_location('r2checker', root / 'code/r2_pm_preflight.py')
    checker = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(checker)
    m = checker.read(root / 'manifest.json')
    passed = ['original input preflight']
    initial = checker.check(root)
    fixture = out / 'NONPHYSICAL_TAMPER_FIXTURE'
    shutil.copytree(root, fixture, ignore=shutil.ignore_patterns('input_package.tar', 'input_package.tar.sha256', 'preparation_validation'))
    checker.check(fixture)
    passed.append('relocated input preflight')

    def reject(name, rel, replacement):
        path = fixture / rel
        old = path.read_bytes()
        try:
            if replacement is None:
                path.unlink()
            else:
                path.write_bytes(replacement)
            try:
                checker.check(fixture)
            except (AssertionError, ValueError, FileNotFoundError, KeyError):
                passed.append(name)
            else:
                raise AssertionError('invalid input accepted: ' + name)
        finally:
            path.write_bytes(old)

    for rel in ['input/coefficients.csv', 'input/r1_coefficients.csv', 'input/receivers.csv',
                'input/frozen_config.json', 'policy.json', 'code/run_r2_pm_seed.m',
                'code/solve_helmholtz_bie_halfplane_vertical.m', 'code/closeout_s3a_pm_fullfk.m',
                'code/r1_bie_continuation_gates.m', 'code/r2_pm_parallel.sh']:
        reject('changed SHA: ' + rel, rel, b'tampered')
    reject('missing own seed coefficients', 'input/coefficients.csv', None)
    for key, value in [('seed', 260004), ('scope', 'R1_BIE_C4_CONTINUATION'), ('single_thread', False),
                       ('single_process', False), ('automatic_retry', True), ('RAM_reserve_fraction', .1),
                       ('max_parallel_seed_processes', 3), ('aggregate_admission_bytes', 12*1024**3)]:
        changed = dict(m, **{key: value})
        reject('manifest policy refusal: ' + key, 'manifest.json', json.dumps(changed).encode())
    changed = dict(m, files=m['files']+[m['files'][0]])
    reject('duplicate manifest path', 'manifest.json', json.dumps(changed).encode())
    changed = dict(m, files=[dict(m['files'][0], path='../escape')]+m['files'][1:])
    reject('unsafe manifest path', 'manifest.json', json.dumps(changed).encode())
    # Even a rehashed policy/config cannot silently relax approved constants.
    for rel, key in [('policy.json', 'gates'), ('input/frozen_config.json', 'BIE_recipes')]:
        path = fixture / rel; old = path.read_bytes(); old_manifest = (fixture / 'manifest.json').read_bytes()
        try:
            obj = checker.read(path)
            if key == 'gates':
                obj[key]['sensitivity'] = 1e-3
            else:
                obj[key][4][0] = 80
            path.write_text(json.dumps(obj))
            changed = checker.read(fixture / 'manifest.json')
            for item in changed['files']:
                if item['path'] == rel:
                    item['sha256'] = checker.sha(path)
            if rel == 'policy.json':
                changed['policy_sha256'] = checker.sha(path)
            (fixture / 'manifest.json').write_text(json.dumps(changed))
            try:
                checker.check(fixture)
            except AssertionError:
                passed.append('rehashed numerical policy/config refusal: ' + rel)
            else:
                raise AssertionError('relaxed frozen setting accepted')
        finally:
            path.write_bytes(old); (fixture / 'manifest.json').write_bytes(old_manifest)
    # PS formula remains byte-for-byte the R1 function, apart from its separator.
    r1 = (root / 'provenance/run_piston_pm_r1_job.m').read_text()
    r2 = (root / 'code/run_r2_piston_pm_job.m').read_text()
    extract = lambda s: s[s.index('function r=local_ps('):s.index('function p=local_boundary_input(') if 'function p=local_boundary_input(' in s else s.index('function v=local_worst(')].strip()
    assert extract(r1) == extract(r2)
    passed.append('R1 PS formula exact text equivalence')
    assert 'reference_data.mat' not in r2 and 'local_gaussian' not in r2 and 'seed_260001_coefficients.csv' not in r2
    assert not list(root.glob('**/data.mat')) and not (root / 'common.mat').exists()
    passed.append('no density/other-seed field imports')
    # Non-acoustic geometry arithmetic only: no normalization or random draws.
    with (root / 'input/coefficients.csv').open() as f:
        profile = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(f)]
    depths, slopes = [], []
    for i in range(14001):
        x = -70 + i*.01; r = abs(x)
        raw = sum(v['cos_coeff_m']*math.cos(v['k_rad_per_m']*x)+v['sin_coeff_m']*math.sin(v['k_rad_per_m']*x) for v in profile)
        derivative = sum(v['k_rad_per_m']*(-v['cos_coeff_m']*math.sin(v['k_rad_per_m']*x)+v['sin_coeff_m']*math.cos(v['k_rad_per_m']*x)) for v in profile)
        chi, cp = (1., 0.) if r <= 42 else (0., 0.)
        if 42 < r < 50:
            t = (r-42)/8; chi = 1-10*t**3+15*t**4-6*t**5
            cp = (-30*t*t+60*t**3-30*t**4)*math.copysign(1, x)/8
        depths.append(raw*chi); slopes.append(derivative*chi+raw*cp)
    assert min(depths)>-2 and max(depths)<3 and all(math.isfinite(v) for v in depths+slopes)
    passed.append('geometry-only clearance/finite grid check')
    report = dict(status='PASS', failed_checks=0, seed=m['seed'], tests=len(passed), checks=passed,
                  new_acoustic_solves=0, SSH_connections=0, MATLAB_executed=False,
                  Code_Analyzer='UNRUN: MATLAB MCP failed to attach',
                  MATLAB_synthetic_resume_tests='UNRUN; delivered and mandatory before acoustic entry',
                  manifest_sha256=checker.sha(root / 'manifest.json'), SHA_bound_files=initial['SHA_bound_files'],
                  geometry_only=dict(min_depth_m=min(depths), max_depth_m=max(depths), max_abs_slope=max(map(abs, slopes))))
    (out / 'receipt.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ('checks',)}, ensure_ascii=False))


if __name__ == '__main__':
    main()
