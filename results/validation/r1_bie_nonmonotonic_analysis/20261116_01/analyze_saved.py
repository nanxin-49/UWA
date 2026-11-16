"""Read-only CSV/JSON/SHA analysis. No MAT loading, MATLAB, SSH or solver."""
import csv, hashlib, json, math, sys
from pathlib import Path

root, out = map(Path, sys.argv[1:3])
def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
    return json.loads(p.read_text(encoding='utf-8-sig'))
def verify_seal(folder, name, signature):
    seal = read(folder/name)
    assert seal['signature'] == signature
    for f in seal['files']:
        assert sha(folder/f['path']) == f['sha256'], f['path']
    return seal

checks = 0
for line in next((root/'retrieval').glob('SHA256SUMS_*.txt')).read_text().splitlines():
    expected, rel = line.split('  ', 1)
    assert sha(root/rel) == expected, rel
    checks += 1
manifest = read(root/'manifest.json')
for f in manifest['files']:
    assert sha(root/f['path']) == f['sha256']
fingerprint = sha(root/'manifest.json')
names = ['bie_coarse','bie_main','bie_self','bie_space']
fields = []
stages = []
pred = fingerprint
for i, name in enumerate(names):
    folder=root/'jobs'/name
    signature=f'{pred}:C{i}:{fingerprint}'
    solved=verify_seal(folder,'solved.json',signature)
    committed=(folder/'committed.json').exists()
    if committed:
        seal=verify_seal(folder,'committed.json',signature)
        assert seal['status']=='AUDITED_PASS'
        pred=sha(folder/'committed.json')
    rows=list(csv.DictReader((folder/'points.csv').open(encoding='utf-8-sig')))
    coords=[float(r['x_m']) for r in rows]
    if fields: assert coords==x
    x=coords
    assert len(x)==411 and all(float(r['z_m'])==3 for r in rows)
    fields.append([complex(float(r['real_Pa']),float(r['imag_Pa'])) for r in rows])
    stages.append(dict(stage=f'C{i}',committed=committed,solved_SHA_valid=True,
                       config=read(folder/'config.json'),gates=read(folder/'gates.json'),
                       audit=read(folder/'audit.json'),timings=read(folder/'timings.json')))
edges=[math.floor(1+410*j/4+0.5)-1 for j in range(5)]
supports=[('line',0,410)]+[(f'Q{j+1}',edges[j],edges[j+1]) for j in range(4)]
norm=lambda v: math.sqrt(math.fsum(abs(z)**2 for z in v))
summary=[]; pointrows=[]
for label, ia, ib in [('C0_C1',0,1),('C2_C1',2,1),('C1_C3',1,3),('C0_C3',0,3)]:
    a,b=fields[ia],fields[ib]; residual=[u-v for u,v in zip(a,b)]
    regions=[]
    for name,lo,hi in supports:
        regions.append(dict(support=name,left_x_m=x[lo],right_x_m=x[hi],
                            relative_L2=norm(residual[lo:hi+1])/norm(b[lo:hi+1])))
    worst=max(regions,key=lambda r:r['relative_L2'])
    peak=max(range(411),key=lambda j:abs(residual[j]))
    energy=math.fsum(abs(z)**2 for z in residual)
    top=sorted(range(411),key=lambda j:abs(residual[j]),reverse=True)
    summary.append(dict(pair=label,regions=regions,worst_support=worst,
        peak_x_m=x[peak],peak_abs_Pa=abs(residual[peak]),
        largest_point_energy_fraction=abs(residual[peak])**2/energy,
        top10_points_energy_fraction=sum(abs(residual[j])**2 for j in top[:10])/energy,
        reference_min_abs_Pa=min(map(abs,b))))
    for j in range(411):
        pointrows.append(dict(pair=label,x_m=x[j],residual_abs_Pa=abs(residual[j])))
gates0=stages[1]['gates']['pairs']
if isinstance(gates0,dict): gates0=[gates0]
for index,pair in enumerate(summary[:3]):
    expected = stages[3]['gates']['pairs'][index]['relative']
    assert abs(pair['worst_support']['relative_L2']-expected)<1e-14
old=summary[0]['worst_support']['relative_L2']; new=summary[2]['worst_support']['relative_L2']
outcome=dict(status='STATIC_SAVED_EVIDENCE_ANALYSIS_COMPLETE',retrieval_SHA_files=checks,
    package_SHA_files=len(manifest['files']),stage_evidence=stages,pairs=summary,
    trigger=dict(first=old,second=new,increase=new-old,ratio=new/old,percent_increase=100*(new/old-1)),
    points_with_increasing_absolute_difference=sum(abs(fields[1][j]-fields[3][j])>abs(fields[0][j]-fields[1][j]) for j in range(411)),
    checkpoint=read(root/'checkpoint.json'),analysis_mode='CSV JSON and SHA only; no MAT load or acoustic calculation')
(out/'analysis.json').write_text(json.dumps(outcome,ensure_ascii=False,indent=2),encoding='utf-8')
with (out/'point_residuals.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(pointrows[0])); w.writeheader(); w.writerows(pointrows)
print(json.dumps({k:outcome[k] for k in ['retrieval_SHA_files','trigger','points_with_increasing_absolute_difference']}))
for p in summary: print(json.dumps(p))
