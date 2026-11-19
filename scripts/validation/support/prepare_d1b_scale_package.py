"""Freeze the 4096/8192 manual pilot. Infrastructure only; no SSH/MATLAB."""
import json
from pathlib import Path
import shutil
import sys
import tarfile
from d1_package import sha, verify

ROOT = Path(__file__).resolve().parents[3]
SUP = ROOT/'scripts/validation/support'
OLD = ROOT/'results/validation/d1b_v2/preparation_20261119_02'
SAVED = ROOT/'results/validation/d1b_v2/d1b_v2_20261119_02_retrieved_20261119_01'
REVIEW = ROOT/'results/validation/d1b_v2/review_20261119_02/audit.json'


def once(s, a, b):
    assert s.count(a) == 1, f'Review frozen upstream before adapting: {a[:90]}'
    return s.replace(a, b)


def instrument(s):
    s = once(s, 'd1b_bie_prototype(cfg,workers)', 'd1b_scale_bie(cfg,workers)')
    s = once(s, 'overall_clock=tic;', "d1b_scale_phase(cfg,'geometry');\noverall_clock=tic;")
    s = once(s, 'cache=local_prepare_cache(edges,gx,cfg);', "d1b_scale_phase(cfg,'cache');\ncache=local_prepare_cache(edges,gx,cfg);")
    s = once(s, 'A = complex(zeros(node_count,node_count));', "d1b_scale_phase(cfg,'base_assembly');\nA = complex(zeros(node_count,node_count));")
    anchor = 'timings.base_assembly_s=toc(phase_clock); phase_clock=tic;'
    # Recompute every uncorrected base row on the serial client. This is an
    # algorithm-equivalence audit, separate from independent GL112 physics.
    s = once(s, anchor, anchor + '''
d1b_scale_phase(cfg,'serial_base_row_audit');
ids=cfg.worker_ids(:);
assert(isequal(sort(unique(assembly_owner)),sort(ids)));
stats=zeros(numel(ids),5);
for first=1:32:node_count
    rows=first:min(first+31,node_count);
    baseline=local_assembly(rows,x,z,nx,nz,k,cfg.halfplane_h_m,weighted_arc);
    delta=A(rows,:)-baseline;
    for w=1:numel(ids)
        take=assembly_owner(rows)==ids(w);
        if ~any(take), continue, end
        err=abs(delta(take,:)); ref=abs(baseline(take,:));
        stats(w,1:3)=stats(w,1:3)+[nnz(take),sum(err.^2,'all'),sum(ref.^2,'all')];
        stats(w,4:5)=max(stats(w,4:5),[max(err,[],'all'),max(ref,[],'all')]);
    end
end
checks=repmat(struct('pid',0,'rows',0,'relative_L2',0,'max_scaled_abs',0),numel(ids),1);
for w=1:numel(ids)
    checks(w)=struct('pid',ids(w),'rows',stats(w,1), ...
        'relative_L2',sqrt(stats(w,2))/max(sqrt(stats(w,3)),realmin), ...
        'max_scaled_abs',stats(w,4)/max(stats(w,5),1));
end
serial_audit=struct('relative_L2',max([checks.relative_L2]), ...
    'max_scaled_abs',max([checks.max_scaled_abs]),'rows',sum([checks.rows]), ...
    'source','serial same-formula row dispatch','workers',checks);
assert(serial_audit.relative_L2<1e-11 && serial_audit.max_scaled_abs<1e-11);
timings.serial_base_row_audit_s=toc(phase_clock);
d1b_scale_phase(cfg,'near_replacement'); phase_clock=tic;
''')
    for marker, stage in [
        ('u_inc_nodes = cfg.incident_fn(x,z);','rhs'),
        ('psi = A\\rhs;','dense_solve'),
        ('[u_receiver,receiver_owner] = d1b_tiles','receivers'),
        ('x_check = 0.5*(edges(1:end-1)+edges(2:end)).\';','boundary_residual')]:
        s = once(s, marker, f"d1b_scale_phase(cfg,'{stage}');\n"+marker)
    s = once(s, 'out.worker_rows = struct(', '''out.serial_assembly_audit = serial_audit;
d1b_scale_phase(cfg,'fixed_density_serial_check');
fixed_checks=struct('pid',{},'receiver',{},'midpoint',{});
if workers==2
    assert(isequal(sort(unique(receiver_owner)),sort(ids)) && isequal(sort(unique(boundary_owner)),sort(ids)));
    for w=1:numel(ids)
        ri=local_pick(receiver_owner,ids(w)); mi=local_pick(boundary_owner,ids(w));
        sr=local_evaluate_off_surface(receiver_x(ri),receiver_z(ri),x,z,nx,nz, ...
            weighted_arc,psi,k,cfg.halfplane_h_m);
        sb=u_check(mi,:)+local_boundary(x_check(mi),z_check(mi),edges,cache,psi,k,cfg);
        fixed_checks(w)=struct('pid',ids(w),'receiver',local_match(u_receiver(ri,:),sr,1e-11), ...
            'midpoint',local_match(boundary_total(mi,:),sb,1e-10));
    end
end
out.fixed_density_serial_checks=fixed_checks;
out.worker_rows = struct(''')
    return s + '''
function ii=local_pick(owner,pid)
rows=find(owner==pid);
ii=rows(unique(round(linspace(1,numel(rows),min(12,numel(rows))))));
end

function checks=local_match(actual,ref,tol)
err=actual-ref;
checks=struct('rows',size(actual,1),'relative_L2',norm(err(:))/max(norm(ref(:)),realmin), ...
    'max_scaled_abs',max(abs(err),[],'all')/max(max(abs(ref),[],'all'),1));
assert(checks.relative_L2<tol && checks.max_scaled_abs<tol);
end
'''


def build(dest):
    dest = Path(dest).resolve()
    assert ROOT/'results/validation/d1b_scale' in dest.parents and not dest.exists()
    review = json.loads(REVIEW.read_text(encoding='utf-8-sig'))
    assert review['status']=='PASS_SAVED_RESOURCE_PROTOTYPE_REVIEW' and review['failed_checks']==0
    oldmanifest = json.loads((OLD/'manifest.json').read_text(encoding='utf-8-sig'))
    assert sha(OLD/'manifest.json') == review['source_manifest_sha256']
    verify(OLD, oldmanifest['files'])
    exit_info=json.loads((SAVED/'output/exit.json').read_text(encoding='utf-8-sig'))
    assert exit_info['matlab_exit_code']==exit_info['supervisor_exit_code']==0
    # Verify the exact selected saved surface against the original exit seal.
    declared={r['path']:r['sha256'] for r in exit_info['stage_files_sha256']}
    code=dest/'code'; inputs=dest/'inputs'; provenance=dest/'provenance'
    for p in [code,inputs,provenance]: p.mkdir(parents=True)
    origins={}
    def cp(src,target):
        shutil.copyfile(src,target); origins[src.relative_to(ROOT).as_posix()]=sha(src)
    for name in ['d1b_tiles.m','d1b_pm_surface.m','d1b_boundary_audit.m',
                 'baffled_strip_piston_field_vertical.m','d1b_remote.py','d1_package.py']:
        cp(OLD/'code'/name,code/name)
    for name in ['d1b_scale_phase.m','d1b_scale_json.m','d1b_scale_hash.m',
                 'audit_d1b_scale_saved.m','d1b_scale_remote.py']:
        cp(SUP/name,code/name)
    cp(ROOT/'scripts/validation/run_d1b_scale_pilot.m',code/'run_d1b_scale_pilot.m')
    cp(ROOT/'scripts/reporting/render_d1b_scale_saved.m',code/'render_d1b_scale_saved.m')
    cp(OLD/'code/d1b_bie_prototype.m',provenance/'d1b_bie_prototype.m')
    (code/'d1b_scale_bie.m').write_text(instrument((OLD/'code/d1b_bie_prototype.m').read_text(encoding='utf-8-sig')),encoding='utf-8',newline='\n')
    for rel,target in [('surfaces/N160_low0/profile.mat','profile.mat'),
                       ('surfaces/N160_low0/coefficients.csv','coefficients.csv'),
                       ('surfaces/surface_receipt.json','surface_receipt.json')]:
        src=SAVED/'output'/rel
        assert declared['output/'+rel]==sha(src)
        cp(src,inputs/target)
    for src in [REVIEW,ROOT/'reports/d1b_v2_resource_review_2026-11-19.md']:
        cp(src,provenance/src.name)
    cp(Path(__file__),provenance/Path(__file__).name)
    cfg=dict(scope='D1B_SCALE_PILOT_ONLY',full_d1b_enabled=False,nodes=[4096,8192],
             control_nodes=1536,half_widths_m=[64,128],workers=2,threads=1,
             max_wall_s=1800,action_cutoff_s=1770,min_available_GiB=18,
             rss_cap_GiB=12,reserve_fraction=.3,min_disk_GiB=20,output_cap_GiB=8,
             source=dict(frequency_hz=4000,c_mps=1500,rho_kgpm3=1025,v0_mps=.001,total_width_m=.1,z_source_m=1490),
             PM=dict(period_m=640,U_PM_mps=6,seed=260001,N=160,low_frequency_enabled=False,plateau_m=220,support_m=240),
             receiver_count=2279,receiver_depth_m=6,receiver_half_width_m=200,
             BIE=dict(ppw=12,panel_order=8,self_order=96,h_m=-2,window_plateau_ratio=.85),
             gates=dict(matrix=1e-11,rhs=1e-11,density=1e-10,receiver=1e-11,midpoint_equivalence=1e-10,
                        independent_equivalence=1e-10,linear=1e-10,boundary=1e-4),
             automatic_retry=False,increase_workers=False,full_large_serial_solve=False,
             large_serial_comparison='all base rows by worker; 12 fixed-density receiver and midpoint rows per worker; all GL112 audit rows; full corrected matrix and fields at 1536 control')
    (dest/'config.json').write_text(json.dumps(cfg,indent=2)+'\n',encoding='utf-8')
    paths=sorted(p for p in dest.rglob('*') if p.is_file())
    manifest=dict(scope=cfg['scope'],origins=origins,source_02_manifest_sha256=sha(OLD/'manifest.json'),
                  files=[dict(path=p.relative_to(dest).as_posix(),sha256=sha(p)) for p in paths])
    mf=dest/'manifest.json'; mf.write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    archive=dest/'input_package.tar'
    with tarfile.open(archive,'x') as tar:
        for p in paths+[mf]: tar.add(p,arcname=p.relative_to(dest).as_posix(),recursive=False)
    (dest/'input_package.tar.sha256').write_text(f'{sha(archive)}  {archive.name}\n',encoding='utf-8')
    print(json.dumps(dict(status='PACKAGE_ONLY_NO_SSH',archive=str(archive),sha256=sha(archive),files=len(paths)),indent=2))

if __name__=='__main__': build(sys.argv[1])
