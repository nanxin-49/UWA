"""Offline recovery provenance, no-resolve admission, carried budgets, safe extraction."""
import ast,json,sys,tarfile,io
from pathlib import Path

def main(package,out):
    package,out=Path(package).resolve(),Path(out).resolve(); assert not out.exists(); out.mkdir(parents=True)
    sys.path.insert(0,str(package/'code'))
    import d1_package as p
    import d1_supervisor as s
    checks=[]
    def check(name,f): f(); checks.append(dict(name=name,passed=True))
    def require(v): assert v
    def rejects(f):
        try:f()
        except (AssertionError,ValueError,FileNotFoundError,KeyError):return
        raise AssertionError('Expected rejection')
    cfg,m=p.preflight(package); recovered=out/'relocated'; p.extract_verified(package/'input_package.tar',recovered)
    check('full manifest and bounded remaining solves',lambda:require(len(m['files'])>=125 and cfg['recovery']['remaining_BIE_solves'] in [5,9]))
    check('relocated full recovery validates original chain',lambda:p.preflight(recovered))
    def mutate(rel,fn):
        target=recovered/rel; original=target.read_bytes(); mf=recovered/'manifest.json'; mb=mf.read_bytes()
        try:
            target.write_bytes(fn(original)); mm=json.loads(mb)
            for r in mm['files']:
                if r['path']==rel:r['sha256']=p.sha(target)
            mf.write_text(json.dumps(mm)); rejects(lambda:p.preflight(recovered))
        finally:target.write_bytes(original);mf.write_bytes(mb)
    for key,value in [('remaining_BIE_solves',10),('prior_supervisor_elapsed_s',0),('prior_RAM_baseline_bytes',1)]:
        def edit(raw,key=key,value=value):
            c=json.loads(raw);c['recovery'][key]=value;return json.dumps(c).encode()
        check('reject reset '+key,lambda edit=edit:mutate('config.json',edit))
    check('reject changed imported C0 even with new outer SHA',lambda:mutate('imported/depth_0100m/jobs/bie_coarse/data.mat',lambda b:b+b'corrupt'))
    check('reject changed imported common',lambda:mutate('imported/depth_0100m/common.mat',lambda b:b+b'corrupt'))
    def fail_proof(raw):
        rec=json.loads(raw);rec['status']='FAILED_STOPPED';return json.dumps(rec).encode()
    check('reject failed numerical proof',lambda:mutate('recovery_audit/receipt.json',fail_proof))
    check('reject changed original coordinate input',lambda:mutate('origin/inputs/receivers.csv',lambda b:b.replace(b'36',b'35',1)))
    prior=cfg['recovery']['prior_supervisor_elapsed_s']
    check('total budget carries original elapsed',lambda:require(s.stop_reason(100,100,5400-prior+.01,0,0,cfg)=='WALL_TIME_STOP'))
    check('remaining budget below total accepted',lambda:require(s.stop_reason(100,100,5400-prior-.01,0,0,cfg) is None))
    check('Tx100 depth budget carries original elapsed',lambda:require(s.stop_reason(100,100,1,2400-prior+.01,0,cfg,100)=='DEPTH_WALL_TIME_STOP'))
    check('Tx1490 does not inherit Tx100 depth time',lambda:require(s.stop_reason(100,100,1,2399,0,cfg,1490) is None))
    check('RAM reserve unchanged',lambda:require(s.stop_reason(29,100,0,0,0,cfg)=='RAM_RESERVE_STOP'))
    check('RSS cap unchanged',lambda:require(s.stop_reason(100,100,0,0,12*1024**3+1,cfg)=='OWNED_GROUP_RSS_STOP'))
    check('partial output never qualifies',lambda:rejects(lambda:s.validate_output(recovered,p.sha(recovered/'manifest.json'))))
    code=(package/'code/run_d1_depth_comparison.m').read_text(encoding='utf-8')
    check('C0 reuse branch has no solver invocation',lambda:require('solve_helmholtz' not in code.split('if reuse\n',1)[1].split('else\n',1)[0]))
    check('C0 independent auditor retained',lambda:require('audit_r1_bie_manual_stage(dest' in code))
    if cfg['recovery']['remaining_BIE_solves']==5:
        check('reject changed C3 data even if outer SHA updated',lambda:mutate('imported/depth_0100m/jobs/bie_space/data.mat',lambda b:b+b'corrupt'))
        check('complete Tx100 analysis is reused',lambda:require("if ~imported_depth, run_d1_model_job(work,'analyze')" in code))
    check('no default-precision output in D1 numerical entry',lambda:require('writetable(' not in code and 'writetable(' not in (package/'code/run_d1_model_job.m').read_text()))
    for f in (package/'code').glob('*.py'):ast.parse(f.read_text(encoding='utf-8'))
    check('Python syntax valid',lambda:None)
    for label,kind,name in [('link',tarfile.SYMTYPE,'link'),('escape',tarfile.REGTYPE,'../bad')]:
        archive=out/(label+'.tar')
        with tarfile.open(archive,'w') as tf:
            t=tarfile.TarInfo(name);t.type=kind;t.linkname='../bad';tf.addfile(t,io.BytesIO(b''))
        check('reject tar '+label,lambda archive=archive,label=label:rejects(lambda:p.extract_verified(archive,out/(label+'_extracted'))))
    session=recovered/'sessions/test';session.mkdir(parents=True)
    (session/'exit.json').write_text(json.dumps(dict(status='FAILED_STOPPED',matlab_exit_code=1)))
    archive,digest=p.archive(recovered);dest=out/'fetched';p.extract_verified(archive,dest)
    check('recovery Fetch includes origin and audit chain',lambda:require((dest/'origin/retrieval/manifest.json').exists() and (dest/'recovery_audit/receipt.json').exists() and p.sha(archive)==digest))
    receipt=dict(status='RECOVERY_OFFLINE_PASS',failed_checks=0,checks=checks,network_calls=0,new_MATLAB_processes=0,new_BIE_solves=0,
        unrun_checks=['Linux native continuation and remaining numerical solves'])
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
    (out/'fixture_archive.txt').write_text(str(archive),encoding='utf-8')
    print(json.dumps(dict(status=receipt['status'],checks=len(checks),network_calls=0,new_BIE_solves=0)))
if __name__=='__main__': main(*sys.argv[1:3])
