"""Offline infrastructure tests; never invokes SSH or MATLAB."""
import importlib.util
import io
import json
from pathlib import Path
import shutil
import sys
import tarfile
import time


def main(package, out):
    package, out = Path(package).resolve(), Path(out).resolve()
    assert not out.exists()
    out.mkdir(parents=True)
    sys.path.insert(0, str(package/'code'))
    import d0_package as p
    import d0_supervisor as s
    tests = []

    def check(name, fn):
        fn()
        tests.append(dict(name=name, passed=True))

    def rejects(fn):
        try:
            fn()
        except (AssertionError, ValueError, FileNotFoundError):
            return
        raise AssertionError('Expected rejection')

    cfg, manifest = p.preflight(package)
    check('manifest_28_and_five_depths', lambda: require(len(manifest['files']) == 28 and len(cfg['depths_m']) == 5))
    relocated = out/'relocated'
    p.extract_verified(package/'input_package.tar', relocated)
    check('archive_relocation', lambda: p.preflight(relocated))
    target = relocated/'code/d0_supervisor.py'
    before = target.read_bytes()
    target.write_bytes(before+b'\n#tampered\n')
    check('tampered_code_rejected', lambda: rejects(lambda: p.preflight(relocated)))
    target.write_bytes(before)
    check('path_traversal_rejected', lambda: rejects(lambda: p.verify(package, [dict(path='../x',sha256='0'*64)])))
    check('cash_reference_rejected', lambda: rejects(lambda: p.verify(package, [dict(path='cash/x',sha256='0'*64)])))
    for name, now, elapsed, depth_elapsed, expected in [
        ('RAM_below',29,0,0,'RAM_RESERVE_STOP'),('RAM_equal',30,0,0,None),
        ('wall_over',100,3601,0,'WALL_TIME_STOP'),('wall_equal',100,3600,0,None),
        ('depth_over',100,1,901,'DEPTH_WALL_TIME_STOP'),('normal',100,1,1,None)]:
        check(name, lambda now=now,elapsed=elapsed,depth_elapsed=depth_elapsed,expected=expected:
              require(s.stop_reason(now,100,elapsed,depth_elapsed,cfg)==expected))
    # Malicious link is rejected before creation of the extraction directory.
    evil=out/'bad_link.tar'
    with tarfile.open(evil,'w') as tf:
        row=tarfile.TarInfo('code/link'); row.type=tarfile.SYMTYPE; row.linkname='../outside'; tf.addfile(row)
    check('tar_symlink_rejected',lambda: rejects(lambda:p.extract_verified(evil,out/'bad_extracted')))
    # Stopped partial output can be archived; archive includes the real exit receipt.
    session=relocated/'sessions/test'; session.mkdir(parents=True)
    (session/'exit.json').write_text(json.dumps(dict(status='FAILED_STOPPED',matlab_exit_code=1)))
    archive,digest=p.archive(relocated)
    fetched=out/'fetched_partial'; p.extract_verified(archive,fetched)
    check('partial_Fetch_SHA_and_exit',lambda: require(p.sha(archive)==digest and json.loads((fetched/'sessions/test/exit.json').read_text())['matlab_exit_code']==1))
    check('partial_output_not_PASS',lambda: rejects(lambda:s.validate_output(fetched,p.sha(fetched/'manifest.json'))))
    # Frozen gates remain rejected even if the outer SHA list were regenerated.
    c=relocated/'config.json'; original=c.read_bytes(); bad=json.loads(original); bad['gates']['source_relative']=1e-3
    c.write_text(json.dumps(bad)); m=relocated/'manifest.json'; original_manifest=m.read_bytes(); rows=json.loads(original_manifest)
    for row in rows['files']:
        if row['path']=='config.json': row['sha256']=p.sha(c)
    m.write_text(json.dumps(rows))
    check('relaxed_gate_rejected',lambda:rejects(lambda:p.preflight(relocated)))
    c.write_bytes(original); m.write_bytes(original_manifest)
    receipt=dict(status='OFFLINE_INFRASTRUCTURE_PASS',checks=tests,failed_checks=0,network_calls=0,MATLAB_processes=0)
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2))
    print(json.dumps(dict(status=receipt['status'],checks=len(tests),network_calls=0,MATLAB_processes=0)))


def require(value):
    assert value


if __name__=='__main__':
    main(*sys.argv[1:3])
