"""Create an allowlisted D1b prototype snapshot from the reviewed research kernels."""
import json
from pathlib import Path
import shutil
import sys
import tarfile

from d1_package import sha


ROOT = Path(__file__).resolve().parents[3]
SUPPORT = ROOT / 'scripts/validation/support'
K0 = ROOT / ('results/experiments/pm_wind_cutoff/'
             'pm_wind_cutoff_k0_20261117_05_retrieved_20261117_232303/output')
BASE = 'ab04a1d544d1a3be05d1e8a67d326ca4e6c2a3bd'


def replace(source, old, new):
    if source.count(old) != 1:
        raise ValueError(f'Review upstream kernel before adapting: {old[:110]}')
    return source.replace(old, new)


def bie_variant(source):
    source = replace(source, 'solve_helmholtz_bie_halfplane_cached_vertical(cfg)',
                     'd1b_bie_prototype(cfg,workers)')
    source = replace(source, '    cfg (1,1) struct',
                     '    cfg (1,1) struct\n    workers (1,1) double {mustBeMember(workers,[0 2 4])}')
    start = source.index('for first = 1:block_size:node_count')
    end = source.index('timings.base_assembly_s', start)
    source = source[:start] + '''assembly_owner = zeros(node_count,1);
for first = 1:block_size:node_count
    local_guard(cfg);
    rows = first:min(first+block_size-1,node_count);
    [A(rows,:),assembly_owner(rows)] = d1b_tiles(numel(rows),32,workers, ...
        @(ii) local_assembly(rows(ii),x,z,nx,nz,k,cfg.halfplane_h_m,weighted_arc));
end

''' + source[end:]
    source = replace(source, '''u_receiver = local_evaluate_off_surface(receiver_x,receiver_z,x,z,nx,nz, ...
    weighted_arc,psi,k,cfg.halfplane_h_m);''', '''[u_receiver,receiver_owner] = d1b_tiles(numel(receiver_x),32,workers, ...
    @(ii) local_evaluate_off_surface(receiver_x(ii),receiver_z(ii),x,z,nx,nz, ...
        weighted_arc,psi,k,cfg.halfplane_h_m));''')
    start = source.index('psi_check = complex(zeros(size(u_check)));')
    end = source.index('boundary_residual = vecnorm',start)
    source = source[:start] + '''[correction,boundary_owner] = d1b_tiles(numel(x_check),8,workers, ...
    @(ii) local_boundary(x_check(ii),z_check(ii),edges,cache,psi,k,cfg));
boundary_total = u_check + correction;
''' + source[end:]
    source = replace(source, 'if isfield(cfg,\'return_system_for_validation\')', '''out.worker_rows = struct('assembly',assembly_owner,'receiver',receiver_owner,'boundary',boundary_owner);
if isfield(cfg,'return_system_for_validation')''')
    return source + '''
function A = local_assembly(rows,x,z,nx,nz,k,h,arc)
[S,K] = local_halfplane_kernels(x(rows),z(rows),x,z,nx,nz,k,h);
A = 2*(K-1i*k*S).*arc.';
end

function correction = local_boundary(x,z,edges,cache,psi,k,cfg)
n = cfg.panel_order;
correction = complex(zeros(numel(x),size(psi,2)));
for tt = 1:numel(x)
    pt = find(x(tt)>=edges(1:end-1) & x(tt)<=edges(2:end),1,'first');
    layer = complex(zeros(1,size(psi,2)));
    for ps = 1:numel(edges)-1
        cols = (ps-1)*n+(1:n);
        template = 1;
        if ps==pt, template = n+2; end
        weights = local_cached_weights(x(tt),z(tt),cache,ps,template,k,cfg.halfplane_h_m);
        layer = layer+weights*psi(cols,:);
    end
    correction(tt,:) = .5*(cache.basis0*psi((pt-1)*n+(1:n),:))+layer;
end
end
'''


def fk_variant(source):
    source = replace(source, 'function r = closeout_s3a_pm_fullfk(src,profile,x,z,panels,extent_m,flat_exact)',
                     'function [field,owner,setup_s,receiver_s,terms] = d1b_fullfk_prototype(src,profile,x,z,panels,extent_m,flat_exact,workers)')
    source = replace(source, '    profile table','    profile (1,1) struct')
    source = replace(source, '    flat_exact (:,1) double\nend','    flat_exact (:,1) double\n    workers (1,1) double\nend\nt0 = tic;')
    source = replace(source, 'closeout_s3a_pm_surface(profile,xs)','d1b_pm_surface(profile,xs)')
    start = source.index('t1=complex(zeros(size(x)))')
    kernel_start = source.index('    g=1i/4',start)
    kernel_end = source.index('\nend\nr=struct',kernel_start)
    kernel = source[kernel_start:kernel_end].replace('z(ii)','z').replace('t1(ii)=','t1=').replace('t2(ii)=','t2=')
    return source[:start] + '''setup_s = toc(t0); t0 = tic;
[terms,owner] = d1b_tiles(numel(x),32,workers, ...
    @(ii) local_terms(x(ii),z(ii),xs,eta,nx,nz,jac,ui,dni,u0,uz0,qw,k));
field = flat_exact+terms(:,1)+terms(:,2);
receiver_s = toc(t0);
end

function terms = local_terms(x,z,xs,eta,nx,nz,jac,ui,dni,u0,uz0,qw,k)
    dx=x-xs.'; dz=z-eta.'; rr=hypot(dx,dz);
''' + kernel + '\n    terms = [t1 t2];\nend\n'


def audit_variant(source):
    start = source.index('[g,w]=local_gauss(112);')
    end = source.index('values.independent_19_boundary_relative', start)
    loop = source[start:end].replace('q.profile','profile').replace('bn.density','density')
    loop = loop.replace('closeout_s3a_pm_surface','d1b_pm_surface')
    return '''function [residual,owner,relative] = d1b_boundary_audit(bn,profile,src,workers)
% GL112 independent product-basis reconstruction at 19 platform midpoints.
L = bn.half_width_m; P = bn.panel_count; n = bn.panel_order; k = bn.k_radpm;
[gx,~] = local_gauss(n);
edges = linspace(-L,L,P+1)'; centers = (edges(1:end-1)+edges(2:end))/2;
half = diff(edges)/2;
mid = centers(abs(centers)<=.85*L);
assert(isequal(mid,bn.boundary_check_x_m));
xt = mid(unique(round(linspace(1,numel(mid),19)))); zt = d1b_pm_surface(profile,xt);
drive = baffled_strip_piston_field_vertical(xt,zt,src,128)-baffled_strip_piston_field_vertical(xt,-zt,src,128);
density = bn.density;
[term,owner] = d1b_tiles(numel(xt),4,workers, ...
    @(ii) local_terms(xt(ii),zt(ii),profile,density,gx,edges,centers,half,n,P,L,k));
residual = drive+term;
relative = norm(residual)/max(norm(drive),realmin);
end

function term = local_terms(xt,zt,profile,density,gx,edges,centers,half,n,P,L,k)
''' + loop + 'term = .5*density0+layer;\nend\n' + source[source.index('function [g,w]=local_gauss'):]


def build(dest):
    dest = Path(dest).resolve()
    assert ROOT/'results/validation/d1b_v2' in dest.parents and not dest.exists()
    code = dest/'code'; inputs = dest/'inputs'
    code.mkdir(parents=True); inputs.mkdir()
    files = [SUPPORT/name for name in (
        'd1b_tiles.m','d1b_ps_prototype.m','d1b_pm_surface.m','d1b_remote.py',
        'd1_package.py','write_d1_csv.m','fk_single_receiver_incident.m',
        'solve_helmholtz_bie_halfplane_cached_vertical.m')]
    files += [ROOT/rel for rel in ('scripts/validation/run_d1b_resource_prototype.m',
              'scripts/experiments/build_d1b_pm_design.m','src/propagation/baffled_strip_piston_field_vertical.m')]
    origins = {}
    for file in files:
        shutil.copyfile(file,code/file.name)
        origins[file.relative_to(ROOT).as_posix()] = sha(file)
    for name,base,modify in (
        ('d1b_bie_prototype.m','solve_helmholtz_bie_halfplane_cached_vertical.m',bie_variant),
        ('d1b_fullfk_prototype.m','closeout_s3a_pm_fullfk.m',fk_variant),
        ('d1b_boundary_audit.m','audit_r1_bie_manual_stage.m',audit_variant)):
        source = SUPPORT/base
        (code/name).write_text(modify(source.read_text(encoding='utf-8-sig')),encoding='utf-8',newline='\n')
        origins[source.relative_to(ROOT).as_posix()] = sha(source)
    origins[Path(__file__).relative_to(ROOT).as_posix()] = sha(__file__)
    for n in (12,40):
        original = K0/f'seed260001_U6_N{n:02d}'
        cp = json.loads((original/'checkpoint.json').read_text(encoding='utf-8-sig'))
        assert cp['status']=='AUDITED_PASS' and cp['case']==dict(U_PM_mps=6,seed=260001,N=n)
        expected = {r['path']: r['sha256'] for r in cp['artifacts_sha256']}
        assert expected['coefficients.csv']==sha(original/'coefficients.csv')
        for src,suffix in [('coefficients.csv','.csv'),('checkpoint.json','_checkpoint.json')]:
            shutil.copyfile(original/src,inputs/f'k0_n{n}{suffix}')
            origins[(original/src).relative_to(ROOT).as_posix()] = sha(original/src)
    cfg = dict(scope='D1B_RESOURCE_PROTOTYPE_ONLY',full_d1b_enabled=False,tx_depth_m=1490,rx_depth_m=6,
               frequency_hz=4000,total_width_m=.1,receiver_count=2279,period_m=640,
               pm_U_mps=6,pm_seed=260001,prototype_case='N160_low0',pool_workers=[0,2,4],
               repeats=2,prototype_half_widths_m=[6,12,24],max_nodes=1536,max_wall_s=1800,
               min_available_GiB=18,rss_cap_GiB=12,reserve_fraction=.3)
    (dest/'config.json').write_text(json.dumps(cfg,indent=2)+'\n',encoding='utf-8')
    paths = sorted(p for p in dest.rglob('*') if p.is_file())
    manifest = dict(scope=cfg['scope'],github_compare_commit=BASE,origins=origins,
                    files=[dict(path=p.relative_to(dest).as_posix(),sha256=sha(p)) for p in paths])
    mf = dest/'manifest.json'; mf.write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    archive = dest/'input_package.tar'
    with tarfile.open(archive,'x') as tar:
        for p in [*paths,mf]:
            tar.add(p,arcname=p.relative_to(dest).as_posix(),recursive=False)
    (dest/'input_package.tar.sha256').write_text(f'{sha(archive)}  {archive.name}\n',encoding='utf-8')
    print(archive,sha(archive))


if __name__=='__main__':
    build(sys.argv[1])
