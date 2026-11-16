"""Local provenance/readiness rejection tests. No MATLAB, solve, SSH, or network."""
import importlib.util
import json
import shutil
import sys
from pathlib import Path

sys.dont_write_bytecode = True


def main():
    root, out = map(lambda p: Path(p).resolve(), sys.argv[1:3])
    out.mkdir(parents=True, exist_ok=False)
    spec = importlib.util.spec_from_file_location('r1preflight', root / 'code/r1_bie_c4_preflight.py')
    checker = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(checker)
    checker.ready(root)
    fixture = out / 'NONPHYSICAL_COPY_FOR_TAMPER_TESTS'
    shutil.copytree(root, fixture, ignore=shutil.ignore_patterns('input_package.tar', 'input_package.tar.sha256'))
    checker.ready(fixture)
    passed = ['original readiness', 'byte-identical fixture readiness']

    def reject(name, relative, replacement):
        path = fixture / relative
        previous = path.read_bytes()
        try:
            if replacement is None:
                path.unlink()
            else:
                path.write_bytes(replacement)
            try:
                checker.ready(fixture)
            except (AssertionError, ValueError, FileNotFoundError, KeyError):
                passed.append(name)
            else:
                raise RuntimeError('invalid evidence accepted: ' + name)
        finally:
            path.write_bytes(previous)

    for name, relative in [
        ('changed original density MAT', 'origin_evidence/jobs/bie_space/data.mat'),
        ('changed original CSV', 'origin_evidence/jobs/bie_main/points.csv'),
        ('missing C3 original audit', 'origin_evidence/jobs/bie_space/audit.json'),
        ('changed original predecessor commit', 'origin_evidence/jobs/bie_self/committed.json'),
        ('changed migrated data', 'jobs/bie_space/data.mat'),
        ('changed migration predecessor', 'jobs/bie_self/migration_commit.json'),
        ('missing migration completion', 'migration_complete.json'),
        ('missing new C3 audit', 'migration_audits/C3/audit.json'),
        ('missing validation', 'local_validation.json'),
        ('missing validation seal', 'local_validation_commit.json'),
        ('missing inherited initial RAM', 'control/initial_available_bytes.txt')]:
        reject(name, relative, None)

    for name, relative, key, value in [
        ('old schema rejected', 'manifest.json', 'schema_version', 1),
        ('old scope rejected', 'manifest.json', 'scope', 'R1_BIE_C1_C4_MANUAL'),
        ('unauthorized solve stage rejected', 'manifest.json', 'allowed_new_solve_stage', 3),
        ('changed policy mode', 'policy.json', 'grid_monotonicity_mode', 'ignore'),
        ('missing implementation approval', 'approval.json', 'implementation_approved', False),
        ('false acceptance rejected', 'local_validation.json', 'status', 'FAIL')]:
        record = checker.read(fixture / relative)
        record[key] = value
        reject(name, relative, json.dumps(record).encode())

    for unsafe in ['../escape', '/absolute', 'C:/escape', 'a/../../escape']:
        try:
            checker.artifact(fixture, unsafe)
        except ValueError:
            passed.append('unsafe path rejected: ' + unsafe)
        else:
            raise RuntimeError('unsafe path accepted')
    for relative, expected in [('policy.json', dict(sensitivity=1.1e-4)), ('policy.json', dict(resolution_factor=4))]:
        record = checker.read(fixture / relative)
        record['gates'].update(expected)
        try:
            checker.policy_check(record)
        except AssertionError:
            passed.append('numerical policy change rejected: ' + str(expected))
        else:
            raise RuntimeError('altered numerical gates accepted')
    checker.ready(root)
    receipt = dict(status='PASS', failed_checks=0, passed_checks=len(passed), checks=passed,
                   execution_mode='Python JSON/SHA/readiness only', MATLAB_starts=0,
                   new_acoustic_solves=0, SSH_connections=0)
    (out / 'receipt.json').write_text(json.dumps(receipt, indent=2), encoding='utf-8')
    print(json.dumps({k: receipt[k] for k in ['status', 'passed_checks', 'new_acoustic_solves', 'SSH_connections']}))


if __name__ == '__main__':
    main()
