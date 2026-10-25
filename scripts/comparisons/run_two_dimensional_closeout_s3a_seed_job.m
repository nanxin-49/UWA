function receipt = run_two_dimensional_closeout_s3a_seed_job(root,run_dir,job_id,available_physical_bytes)
%RUN_TWO_DIMENSIONAL_CLOSEOUT_S3A_SEED_JOB Guarded seed260001 native-field job.
arguments
    root (1,:) char
    run_dir (1,:) char
    job_id (1,:) char
    available_physical_bytes (1,1) double {mustBePositive}
end
assert(isfolder(root) && isfolder(run_dir) && contains(run_dir,'s3a_seed260001'));
ids={'main96','self64','space14','window60','window70','pe','fullfk'};
ix=find(strcmp(ids,job_id)); assert(isscalar(ix),'Unknown frozen S3A job');
out=fullfile(run_dir,'jobs',[job_id '.mat']); assert(~isfile(out),'No job overwrite');
for j=1:ix-1
    assert(isfile(fullfile(run_dir,'jobs',[ids{j} '.mat'])),'Earlier job is missing');
end
budget=jsondecode(fileread(fullfile(run_dir,'protocol','budget.json')));
assert(budget.seed==260001 && budget.subbudget_s==1700 && budget.shared_prior_s==3401.424104);
assert(budget.shared_cap_s==7506.601 && budget.S3A_total_s==3000);
upper=budget.upper_job_s(:); actual=zeros(7,1);
for j=1:ix-1
    a=load(fullfile(run_dir,'jobs',[ids{j} '.mat']),'receipt');
    assert(strcmp(a.receipt.status,'COMPLETE'));
    actual(j)=a.receipt.wall_seconds;
end
upper(1:ix-1)=actual(1:ix-1);
assert(sum(upper)+budget.failure_reserve_s<=budget.subbudget_s,'S3A seed cost gate');
assert(max(upper(ix:end))<900,'Single-job time gate');
assert(budget.shared_prior_s+sum(upper)+budget.failure_reserve_s <= budget.shared_cap_s);
coeff_file=fullfile(run_dir,'input','seed_260001_coefficients.csv');
source_file=fullfile(run_dir,'input','frozen_source.mat');
reference_file=fullfile(run_dir,'input','reference_data.mat');
profile=readtable(coeff_file); s=load(source_file,'src'); src=s.src;
t=load(reference_file,'x','z','flat_exact','U','ix9');
x=t.x(:); z=t.z(:); assert(numel(x)==417 && numel(z)==417 && numel(t.ix9)==9);
assert(numel(src.coeff)==984 && max(profile.k_rad_per_m)<.4712389+1e-9);
assert(all(abs(x(1:411))<42) && all(z(1:411)==3));
N=[3840 3840 4480 3200 3736];
peak=[4*16*N.^2+128*1024^2,128*1024^2,256*1024^2];
assert(peak(ix)<=.7*available_physical_bytes,'30 percent physical RAM headroom');
assert(available_physical_bytes>0);
pre=struct('job_id',job_id,'available_physical_bytes',available_physical_bytes, ...
    'estimated_increment_bytes',peak(ix),'reserve_fraction',.3, ...
    'projected_seed_core_seconds',sum(upper)+budget.failure_reserve_s, ...
    'seed_subbudget_s',budget.subbudget_s,'same_source',true,'same_receivers',true);
local_json(fullfile(run_dir,'jobs',[job_id '_preflight.json']),pre);
clock=tic;
try
    if ix<=5
        L=[60 60 60 60 70]; ppw=[12 12 14 10 10]; qa=[96 64 96 96 96];
        cfg=struct('frequency_hz',4000,'c0_mps',1500,'half_width_m',L(ix), ...
            'points_per_wavelength',ppw(ix),'panel_order',8, ...
            'self_quadrature_order',qa(ix),'halfplane_h_m',-2, ...
            'window_plateau_ratio',.85, ...
            'eta_fn',@(xx) closeout_s3a_pm_surface(profile,xx), ...
            'eta_prime_fn',@(xx) local_slope(profile,xx), ...
            'incident_fn',@(xx,zz) fk_single_receiver_incident(src,xx,zz), ...
            'receiver_x_m',x,'receiver_z_m',z);
        boundary=solve_helmholtz_bie_halfplane_vertical(cfg);
        assert(boundary.finite && boundary.node_count==N(ix));
        payload=struct('u_native',boundary.receiver_field(:),'x_m',x,'z_m',z, ...
            'boundary',boundary,'half_width_m',L(ix),'ppw',ppw(ix),'qa',qa(ix), ...
            'version','current frozen BIE combined-layer R0');
        key=struct('node_count',boundary.node_count, ...
            'linear_residual',boundary.linear_residual, ...
            'offgrid_boundary_residual',boundary.offgrid_boundary_residual);
    elseif ix==6
        W=192.1875; x0=(-W/2)+(0:983)'*(W/984);
        source=ifft(src.coeff(:));
        assert(max(abs(source-exp(-.5*(x0/.3).^2)))<1e-12);
        [eta,~,~,~]=closeout_s3a_pm_surface(profile,x0);
        cfg=struct('frequency_hz',4000,'c0_mps',1500,'xw_m',W,'nx',984, ...
            'z_tx_m',100,'z_rx_m',3,'sigma_src_m',.3, ...
            'surface_reflect_coeff',-1,'step_m',.05,'x_rx_m',0, ...
            'reflection_model','model0_normal','source_field_override',source, ...
            'eta_z_m',eta);
        flat_cfg=cfg; flat_cfg.eta_z_m=zeros(size(eta));
        a=run_pe_physical_surface_model0_validation(flat_cfg);
        [ok,ii]=ismember(x(1:411),x0); assert(all(ok));
        flat_native=a.reflected_native(:);
        flat_relative=norm(flat_native(ii)-t.U(1:411,1,2))/norm(t.U(1:411,1,2));
        assert(flat_relative<1e-10,'Same-source flat PE/BIE control');
        u=complex(zeros(417,1)); half=u;
        for zz=unique(z).'
            cfg.z_rx_m=zz; a=run_pe_physical_surface_model0_validation(cfg);
            rows=z==zz; [ok,jj]=ismember(x(rows),x0); assert(all(ok));
            v=a.reflected_native(:); u(rows)=v(jj);
            cfg2=cfg; cfg2.step_m=.025;
            b=run_pe_physical_surface_model0_validation(cfg2);
            v=b.reflected_native(:); half(rows)=v(jj);
        end
        payload=struct('u_native',u,'half_step_native',half,'x_m',x,'z_m',z, ...
            'flat_relative_to_BIE',flat_relative, ...
            'version','PE-mode0 S2a-depth-positive-native-v1');
        key=struct('flat_relative_to_BIE',flat_relative, ...
            'half_step_relative',norm(u-half)/norm(u));
    else
        flat=t.flat_exact(:);
        a=closeout_s3a_pm_fullfk(src,profile,x,z,400,50,flat);
        b=closeout_s3a_pm_fullfk(src,profile,x,z,800,50,flat);
        c=closeout_s3a_pm_fullfk(src,profile,x,z,960,60,flat);
        d=closeout_s3a_pm_fullfk(src,profile,x,z,1440,90,flat);
        payload=struct('u_native',b.receiver_field(:),'x_m',x,'z_m',z, ...
            'GL8_400',a,'GL8_800',b,'extent60',c,'extent90',d, ...
            'version','FullFK_complete_analytic_flat_exterior_GL8_v1');
        key=struct('quadrature_relative',norm(b.receiver_field-a.receiver_field)/norm(b.receiver_field), ...
            'extent60_90_relative',norm(c.receiver_field-d.receiver_field)/norm(b.receiver_field), ...
            'component_max_abs',max(abs(b.receiver_field-b.term_incident_dn-b.term_green_dn)));
    end
    seconds=toc(clock); assert(seconds<900,'Single job exceeded 900 s');
    assert(all(isfinite(payload.u_native)) && numel(payload.u_native)==417);
    receipt=struct('job_id',job_id,'status','COMPLETE','wall_seconds',seconds, ...
        'available_physical_bytes',available_physical_bytes,'estimated_increment_bytes',peak(ix), ...
        'primary_points',417,'nine_points',numel(t.ix9),'key_metrics',key);
    save(out,'payload','receipt','-v7.3');
    local_json(fullfile(run_dir,'jobs',[job_id '.json']),receipt);
catch ME
    seconds=toc(clock);
    failure=struct('job_id',job_id,'status','FAILED','wall_seconds',seconds, ...
        'message',ME.message,'identifier',ME.identifier);
    local_json(fullfile(run_dir,'jobs',[job_id '_failure.json']),failure);
    rethrow(ME);
end
disp(receipt);
end

function y=local_slope(profile,x)
[~,y]=closeout_s3a_pm_surface(profile,x);
end

function local_json(path,value)
fid=fopen(path,'w'); assert(fid>0); c=onCleanup(@() fclose(fid));
fprintf(fid,'%s',jsonencode(value));
end
