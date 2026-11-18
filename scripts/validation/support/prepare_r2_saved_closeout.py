"""Freeze provenance for local saved-data R2 closeout; never MATLAB/SSH/solve."""
import hashlib,json,re,shutil,sys,time
from pathlib import Path

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p): return json.loads(p.read_text(encoding='utf-8-sig'))
def safe(root,rel):
    assert not re.search(r'(^/|:|\\|(^|/)\.\.(/|$))',rel)
    p=(root/rel).resolve(); assert p.is_relative_to(root.resolve()); return p
def files(root,rows):
    assert len(rows)==len({f['path'] for f in rows})
    for f in rows: assert sha(safe(root,f['path']))==f['sha256'],f['path']
def seal(root,name,sig=None):
    s=read(root/name); assert s['files']
    if sig is not None: assert s['signature']==sig
    files(root,s['files']); return s

def main():
    project,out=map(lambda s:Path(s).resolve(),sys.argv[1:3]); mode=sys.argv[3]
    assert out.is_relative_to(project/'results/comparisons/piston_2d_pm')
    assert mode in ('prepare','postcheck'); out.mkdir(parents=True,exist_ok=True)
    before=time.perf_counter(); records=[]
    for seed in (260001,260002,260003):
        if seed==260001:
            source=project/'results/comparisons/piston_2d_pm/r1_bie_manual_c4_policy_v2_20261116_04_retrieved_20261116_173609'
        else:
            source=project/f'results/comparisons/piston_2d_pm/r2_seed{seed}_retrieved_20261117_01'
        lists=list((source/'retrieval').glob('SHA256SUMS_*.txt')); assert len(lists)==1
        rows=lists[0].read_text().splitlines(); rels=[]
        for row in rows:
            digest,rel=row.split('  ',1); assert sha(safe(source,rel))==digest,rel; rels.append(rel)
        assert len(rows)=={260001:287,260002:229,260003:239}[seed] and len(set(rels))==len(rows)
        manifest=read(source/'manifest.json'); files(source,manifest['files']); fp=sha(source/'manifest.json')
        cp=read(source/'checkpoint.json'); final=source/'finalization'; work=final/'attempts/attempt_0001/work'
        if seed!=260001:
            assert fp==sha(project/f'results/comparisons/piston_2d_pm/r2_preparation_20261117_03/seed{seed}/manifest.json')
            d=source/'model_stage'; seal(d,'solved.json',fp+':MODELS:'+fp);seal(d,'committed.json',fp+':MODELS:'+fp)
            pred=sha(d/'committed.json')
            for i,name in enumerate(['bie_coarse','bie_main','bie_self','bie_space','bie_window']):
                d=source/'jobs'/name; sig=pred+f':C{i}:'+fp
                assert seal(d,'committed.json',sig)['status']=='AUDITED_PASS'
                seal(d,'solved.json',sig)
                assert read(d/'audit.json')['failed_checks']==0 and read(d/'gates.json')['status']=='PASS'
                assert sha(d/'data.mat')==sha(work/'jobs'/name/'data.mat')
                pred=sha(d/'committed.json')
            final_sig=pred+':FINAL:'+fp
        else:
            final_sig=sha(source/'jobs/bie_window/committed.json')+':FINAL:'+fp+':'+sha(source/'policy.json')
            prior=project/'results/comparisons/piston_2d_pm/20261116_r1_v2_review_01/review_receipt.json'
            r=read(prior); assert r['status']=='SAVED_RESULT_REVIEW_PASS' and r['failed_checks']==0
        has_final=(final/'committed.json').exists()
        if has_final: assert seal(final,'committed.json',final_sig)['status']=='AUDITED_PASS'
        assert has_final==(seed!=260002)
        session=source/'sessions'/cp['session_id'] if 'session_id' in cp else None
        rec=dict(seed=seed,source=str(source),work=str(work),retrieval_sha_files=len(rows),
            retrieval_manifest_sha256=sha(lists[0]),manifest_sha256=fp,policy_sha256=sha(source/'policy.json'),
            original_status=cp['status'],original_final_seal=has_final,
            original_real_exit_code=int((session/'matlab_exit_code.txt').read_text()) if session else None,prior_formal_review=None)
        if seed==260001: rec['prior_formal_review']=dict(path=str(prior),sha256=sha(prior),status=r['status'])
        records.append(rec)
        if mode=='prepare' and seed!=260001:
            dest=out/f'seed{seed}'; dest.mkdir(exist_ok=False)
            (dest/'provenance').mkdir()
            for name in ['manifest.json','policy.json','checkpoint.json']:
                shutil.copyfile(source/name,dest/'provenance'/('original_'+name))
            if seed==260002:
                failure=dest/'provenance/original_failure';failure.mkdir()
                for name in ['F02_pointwise_errors.png','F02_audit_reexport.png','F02_pointwise_errors.fig']:
                    shutil.copyfile(work/'figures_final'/name,failure/name)
                shutil.copyfile(session/'stderr.log',failure/'original_stderr.log')
            (dest/'provenance/source.json').write_text(json.dumps(rec,ensure_ascii=False,indent=2),encoding='utf-8')
    report=dict(status='SOURCE_SHA_AND_CHAIN_PASS',failed_checks=0,mode=mode,elapsed_s=time.perf_counter()-before,
        sources=records,new_acoustic_solves=0,SSH_connections=0,source_results_unchanged=True)
    if mode=='postcheck':
        pre=read(out/'source_preflight.json')
        assert pre['sources']==report['sources']
    target=out/('source_preflight.json' if mode=='prepare' else 'source_postcheck.json')
    assert not target.exists();target.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:report[k] for k in ['status','mode','elapsed_s','failed_checks','new_acoustic_solves','SSH_connections']}))

if __name__=='__main__': main()
