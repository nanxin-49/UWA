"""Offline D1 integrity/stop/relocation tests; never invokes SSH or MATLAB."""
import ast
import io
import json
from pathlib import Path
import sys
import tarfile

def require(x): assert x
def main(package,out):
    package,out=Path(package).resolve(),Path(out).resolve(); assert not out.exists(); out.mkdir(parents=True)
    sys.path.insert(0,str(package/'code'))
    import d1_package as p
    import d1_supervisor as s
    checks=[]
    def test(name,fn): fn(); checks.append(dict(name=name,passed=True))
    def rejects(fn):
        try: fn()
        except (AssertionError,ValueError,FileNotFoundError,KeyError): return
        raise AssertionError('Expected rejection')
    cfg,m=p.preflight(package)
    test('31 manifest inputs and two frozen depths',lambda:require(len(m['files'])==31 and cfg['depths_m']==[100,1490]))
    relocated=out/'relocated'; p.extract_verified(package/'input_package.tar',relocated)
    test('relocated package passes',lambda:p.preflight(relocated))
    def mutate(rel,edit,rehash=False):
        f=relocated/rel; before=f.read_bytes(); mf=relocated/'manifest.json'; mb=mf.read_bytes()
        try:
            f.write_text(edit(before.decode()),encoding='utf-8')
            if rehash:
                mm=json.loads(mb)
                for row in mm['files']:
                    if row['path']==rel: row['sha256']=p.sha(f)
                mf.write_text(json.dumps(mm))
            rejects(lambda:p.preflight(relocated))
        finally: f.write_bytes(before); mf.write_bytes(mb)
    test('code tamper rejected',lambda:mutate('code/d1_supervisor.py',lambda t:t+'\n# tamper'))
    test('relaxed sensitivity rejected even if rehashed',lambda:mutate('config.json',lambda t:t.replace('"sensitivity": 0.0001','"sensitivity": 0.01'),True))
    test('third depth rejected even if rehashed',lambda:mutate('config.json',lambda t:t.replace('1490','1491'),True))
    test('Rx3 rejected even if rehashed',lambda:mutate('config.json',lambda t:t.replace('"receiver_z_m": 6','"receiver_z_m": 3'),True))
    test('coefficient drift rejected even if rehashed',lambda:mutate('inputs/coefficients.csv',lambda t:t.replace(t.splitlines()[1],t.splitlines()[1].replace(',','0,',1)),True))
    test('traversal rejected',lambda:rejects(lambda:p.verify(package,[dict(path='../escape',sha256='0'*64)])))
    test('forbidden archive reference rejected',lambda:rejects(lambda:p.verify(package,[dict(path='cash/x',sha256='0'*64)])))
    for label,mem,elapsed,depth,rss,expected in [
        ('RAM below',29,0,0,0,'RAM_RESERVE_STOP'),('RAM equal',30,0,0,0,None),
        ('wall over',100,5401,0,0,'WALL_TIME_STOP'),('wall equal',100,5400,0,0,None),
        ('depth over',100,1,2401,0,'DEPTH_WALL_TIME_STOP'),('depth equal',100,1,2400,0,None),
        ('RSS over',100,1,1,12*1024**3+1,'OWNED_GROUP_RSS_STOP'),('RSS equal',100,1,1,12*1024**3,None)]:
        test(label,lambda mem=mem,elapsed=elapsed,depth=depth,rss=rss,expected=expected:require(s.stop_reason(mem,100,elapsed,depth,rss,cfg)==expected))
    for label,type_,name in [('symlink',tarfile.SYMTYPE,'code/link'),('traversal',tarfile.REGTYPE,'../escape'),('duplicate',tarfile.REGTYPE,'same')]:
        evil=out/(label+'.tar')
        with tarfile.open(evil,'w') as tf:
            ti=tarfile.TarInfo(name); ti.type=type_; ti.linkname='../escape'; tf.addfile(ti,io.BytesIO(b''))
            if label=='duplicate': tf.addfile(ti,io.BytesIO(b''))
        test('tar '+label+' rejected',lambda evil=evil,label=label:rejects(lambda:p.extract_verified(evil,out/(label+'_extract'))))
    session=relocated/'sessions/test'; session.mkdir(parents=True)
    (session/'exit.json').write_text(json.dumps(dict(status='FAILED_STOPPED',matlab_exit_code=1)))
    archive,digest=p.archive(relocated); fetched=out/'fetched_partial'; p.extract_verified(archive,fetched)
    test('failed partial archive SHA retained',lambda:require(p.sha(archive)==digest and json.loads((fetched/'sessions/test/exit.json').read_text())['matlab_exit_code']==1))
    test('partial outputs cannot pass',lambda:rejects(lambda:s.validate_output(fetched,p.sha(fetched/'manifest.json'))))
    for f in (package/'code').glob('*.py'): ast.parse(f.read_text())
    test('Python AST syntax valid',lambda:None)
    # Audit exact-depth substitutions, preventing the known Rx3 propagation bug.
    model=(package/'code/run_d1_model_job.m').read_text(); audit=(package/'code/audit_d1_models.m').read_text()
    test('PS propagates actual Rx6 and preserves carrier',lambda:require('depth=rz(1)' in model and 'exp(1i*k*depth)' in model and 'ceil(3/' not in model))
    test('independent Green uses actual Rx depth',lambda:require("hypot(q.x(jj)-xp.',q.z(jj))" in audit and '.*q.z(jj)./rr' in audit))
    test('independent deep source segmented quadrature',lambda:require('d0_incident_spectrum_reference(xx,zz,s,32,4)' in audit))
    receipt=dict(status='OFFLINE_INFRASTRUCTURE_PASS',checks=checks,failed_checks=0,network_calls=0,MATLAB_processes=0,
        unrun_checks=['MATLAB parser/Analyzer/dry-run','MATLAB numerical audits and FIG reconstruction','Linux flock/signals/resource sampling and real supervisor workload'])
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)); (out/'fixture_archive.txt').write_text(str(archive))
    print(json.dumps(dict(status=receipt['status'],checks=len(checks),network_calls=0,MATLAB_processes=0)))
if __name__=='__main__': main(*sys.argv[1:3])
