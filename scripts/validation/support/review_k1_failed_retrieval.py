"""Verify a stopped K1 retrieval and inventory evidence, without executing it."""
import csv
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys
import time


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''): h.update(block)
    return h.hexdigest()


def read(path): return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def main(source, package, output):
    start=time.monotonic()
    source=Path(source).resolve(); package=Path(package).resolve(); output=Path(output).resolve()
    assert not output.exists(); output.mkdir(parents=True)
    repo=Path(__file__).resolve().parents[3]
    # Execute only the already-reviewed local file-chain checker, not fetched code.
    file=repo/'scripts/validation/support/k1_pm_preflight.py'
    assert sha(file)==sha(source/'code/k1_pm_preflight.py')==sha(package/'code/k1_pm_preflight.py')
    spec=importlib.util.spec_from_file_location('k1_review_preflight',file)
    pre=importlib.util.module_from_spec(spec); spec.loader.exec_module(pre)
    transport=read(source/'transport_receipt.json')
    archive_sha=sha(source/'result.tar')
    assert transport['transport_pass'] and archive_sha==transport['expected_remote_archive_sha256']==transport['local_archive_sha256']
    rows=read(source/'retrieval/manifest.json'); count=pre.verify_rows(source,rows)
    assert pre.sha(source/'manifest.json')==pre.sha(package/'manifest.json')=='815e44af687af278bcd45b5a05ee7be60b638aaf0c97901e2ac117dd082b11ec'
    frozen=pre.verify(source); completed=pre.completed_case(source/'cases/U6_N07')
    checkpoint=read(source/'checkpoint.json'); assert checkpoint['completed_cases']==['U6_N07']
    sid=Path((source/'latest_session.txt').read_text().strip()).name
    session=read(source/'sessions'/sid/'exit.json'); assert session['case_exit_codes']==dict(preflight=0,U6_N07=0,U6_N40=-9)
    case=source/'cases/U6_N40'; exit_record=read(case/'sessions'/sid/'exit.json')
    err=(case/'sessions'/sid/'stderr.log').read_text(); assert 'double free or corruption (out)' in err
    assert exit_record['stop_reason'] is None and not (case/'model_stage/solved.json').exists() and not (case/'model_stage/committed.json').exists()
    marker=read(case/'model_stage/solve_invocation.json'); assert marker['invocations']==1
    assert not (case/'jobs').exists()
    case_status=[]
    for c in read(source/'config.json')['cases']:
        p=source/'cases'/c['case_id']; cp=read(p/'checkpoint.json') if (p/'checkpoint.json').exists() else None
        case_status.append(dict(case_id=c['case_id'],checkpoint=cp,has_MODELS_invocation=(p/'model_stage/solve_invocation.json').exists(),has_final_seal=(p/'finalization/committed.json').exists()))
    inventory=[]
    for p in source.rglob('*'):
        if p.is_file() and p.name!='result.tar':
            inventory.append(dict(path=p.relative_to(source).as_posix(),bytes=p.stat().st_size,sha256=sha(p)))
    crashes=[r for r in inventory if any(x in Path(r['path']).name.lower() for x in ['crash','core.','hs_err','stacktrace'])]
    resources=[]
    for c in ['U6_N07','U6_N40']:
        path=source/'cases'/c/'sessions'/sid/'resource.csv'
        with path.open(newline='') as f: records=list(csv.DictReader(f))
        resources.append(dict(case_id=c,samples=len(records),min_available_bytes=min(int(r['MemAvailable_bytes']) for r in records),peak_RSS_bytes=max(int(r['owned_group_RSS_bytes']) for r in records),last=records[-1]))
    # Header inventory only here; native whos('-file') handles MATLAB tables.
    mats=[]
    for p in [source/'cases/U6_N07/finalization/attempts/attempt_0001/work/common.mat',source/'cases/U6_N07/finalization/attempts/attempt_0001/work/comparison.mat'] + [source/'cases/U6_N07/jobs'/n/'data.mat' for n in ['bie_coarse','bie_main','bie_self','bie_space','bie_window']]:
        with p.open('rb') as f: header=f.read(128)
        assert header.startswith(b'MATLAB 5.0 MAT-file')
        mats.append(dict(path=p.relative_to(source).as_posix(),header=header[:116].decode('ascii',errors='replace').rstrip(),variables='Native MATLAB whos-file required; scipy whosmat cannot read existing table serialization'))
    (output/'retrieved_inventory_sha256.json').write_text(json.dumps(inventory,indent=2))
    (output/'MAT_inventory.json').write_text(json.dumps(mats,indent=2))
    (output/'case_states.json').write_text(json.dumps(case_status,indent=2))
    receipt=dict(status='K1_FAILED_RETRIEVAL_SHA_AND_STATE_PASS',failed_checks=0,source=str(source),session=sid,transport_archive_sha256=archive_sha,verified_retrieval_files=count,frozen_input_files=frozen['files'],U6_N07=completed,U6_N40=dict(exit=exit_record,stderr=err,stage='MODELS',new_BIE_stages=0,unsealed_invocation=True),crash_artifacts_in_retrieval=crashes,resources=resources,elapsed_s=time.monotonic()-start,new_acoustic_solves=0)
    (output/'receipt.json').write_text(json.dumps(receipt,indent=2))
    print(json.dumps({k:v for k,v in receipt.items() if k not in ['U6_N07','U6_N40','resources']}))


if __name__=='__main__': main(*sys.argv[1:])
