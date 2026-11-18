"""Bind completed local artifacts after MATLAB diaries close; no physics."""
import hashlib,json,sys
from pathlib import Path

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def verify(root):
    s=read(root/'committed.json');assert s['status']=='AUDITED_PASS' and s['files']
    for f in s['files']:
        p=(root/f['path']).resolve();assert p.is_relative_to(root.resolve()) and sha(p)==f['sha256']
    return s
def bind(root,signature):
    assert not (root/'committed.json').exists()
    paths=sorted(p for p in root.rglob('*') if p.is_file())
    records=[dict(path=p.relative_to(root).as_posix(),sha256=sha(p)) for p in paths]
    record=dict(schema_version=1,status='AUDITED_PASS',signature=signature,files=records)
    (root/'committed.json').write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf-8')
    verify(root);return len(records)
def main():
    case=Path(sys.argv[1]).resolve();mode=sys.argv[2]
    if mode=='seeds':
        reports=[]
        for seed in (260002,260003):
            root=case/f'seed{seed}'; receipt=read(root/'closeout_receipt.json')
            assert receipt['status']=='FORMAL_LOCAL_SAVED_CLOSEOUT_PASS' and receipt['failed_checks']==0
            n=read(root/'numeric_receipt.json'); f=read(root/'figures_receipt.json')
            assert n['failed_checks']==f['failed_checks']==0 and len(n['stages'])==5
            assert all(s['failed_checks']==0 for s in n['stages'])
            assert f['new_figure_audit']['pixel_tolerance']==0
            for fig in f['new_figure_audit']['records'].values(): assert fig['pixel_differences']==[0,0] and fig['layout_exact']
            plan=read(root/'seal_plan.json');count=bind(root,plan['signature'])
            reports.append(dict(seed=seed,status=receipt['status'],SHA_bound_files=count,committed_sha256=sha(root/'committed.json')))
        (case/'seed_seal_receipt.json').write_text(json.dumps(dict(status='BOTH_SEED_SEALS_VERIFIED',failed_checks=0,seeds=reports),indent=2),encoding='utf-8')
        print(json.dumps(reports))
    elif mode=='verify':
        for seed in (260002,260003):verify(case/f'seed{seed}')
        if (case/'committed.json').exists():verify(case)
        print('ALL_PRESENT_LOCAL_SEALS_SHA_PASS')
    elif mode=='master':
        for seed in (260002,260003):verify(case/f'seed{seed}')
        assert read(case/'source_postcheck.json')['failed_checks']==0
        summary=read(case/'three_seed_summary/qualified_receipt.json');assert summary['failed_checks']==0
        assert read(case/'final_checkpoint.json')['status']=='R2_SAVED_CLOSEOUT_COMPLETE'
        count=bind(case,'R2:SAVED_CLOSEOUT:20261117:260002:260003')
        print(json.dumps(dict(status='R2_MASTER_FINAL_SEAL_PASS',SHA_bound_files=count,sha256=sha(case/'committed.json'))))
    else:raise ValueError(mode)
if __name__=='__main__':main()
