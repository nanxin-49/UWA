function receipt = run_two_dimensional_s5a_weak_bie_remote(task_root,job_id,available_physical_bytes)
%RUN_TWO_DIMENSIONAL_S5A_WEAK_BIE_REMOTE One guarded 6 kHz P22 reference layer.
% Task-root code and source are transferred only into the approved isolated dir.
arguments
    task_root (1,:) char
    job_id (1,:) char
    available_physical_bytes (1,1) double {mustBeFinite,mustBePositive}
end
assert(isfolder(task_root) && isfolder(fullfile(task_root,'code')));
addpath(fullfile(task_root,'code'));
if strcmp(job_id,'smoke')
    receipt=local_smoke(task_root,available_physical_bytes);
    return
end
if strcmp(job_id,'model')
    receipt=local_model(task_root,available_physical_bytes);
    return
end
ids={'main96','self64','space14','window60','window70'};
ix=find(strcmp(ids,job_id)); assert(isscalar(ix),'Unknown job id');
out=fullfile(task_root,'jobs',job_id);
assert(isfolder(out) && ~isfile(fullfile(out,'result.mat')));
assert(isfile(fullfile(task_root,'jobs','model','result.mat')), ...
    'Remote paired PE/FullFK model control must pass first');
for jj=1:ix-1
    assert(isfile(fullfile(task_root,'jobs',ids{jj},'result.mat')), ...
        'Earlier reference layer is missing');
end
source_path=fullfile(task_root,'input','frozen_source.mat');
assert(isfile(source_path));
s=load(source_path,'src'); src=closeout2d_source_at_frequency(s.src,6000,1500);
assert(src.n_original==984 && src.z_tx_m==100);
tx_x=src.x0_m+(0:983)'*192.1875/984;
tx_expected=ifft(src.coeff(:));
tx_actual=fk_single_receiver_incident(src,tx_x,100+zeros(984,1));
assert(norm(tx_actual-tx_expected)/norm(tx_expected)<1e-10);
L=[60 60 60 60 70]; ppw=[12 12 14 10 10];
qa=[96 64 96 96 96];
lambda=1500/6000;
panels=ceil(2*L(ix)/(lambda/ppw(ix)*8));
node_estimate=8*panels;
estimated_increment_bytes=64*node_estimate^2+128*1024^2;
assert(estimated_increment_bytes<=0.7*available_physical_bytes, ...
    '30 percent physical RAM headroom gate');
cfg=struct('frequency_hz',6000,'c0_mps',1500, ...
    'half_width_m',L(ix),'points_per_wavelength',ppw(ix), ...
    'panel_order',8,'self_quadrature_order',qa(ix), ...
    'halfplane_h_m',-2,'window_plateau_ratio',.85, ...
    'eta_fn',@(x) closeout_s2b_surface(x,.01), ...
    'eta_prime_fn',@(x) local_slope(x), ...
    'incident_fn',@(x,z) fk_single_receiver_incident(src,x,z), ...
    'receiver_x_m',0,'receiver_z_m',3);
pre=struct('job_id',job_id,'frequency_hz',6000,'amplitude_m',.01, ...
    'receiver_x_m',0,'receiver_z_m',3,'node_estimate',node_estimate, ...
    'estimated_increment_bytes',estimated_increment_bytes, ...
    'available_physical_bytes',available_physical_bytes, ...
    'ram_reserve_fraction',.3,'solver','frozen S3A halfplane BIE', ...
    'source','frozen 984-term Gaussian; kz updated at 6 kHz');
local_write_json(fullfile(out,'preflight.json'),pre);
tic_job=tic;
boundary=solve_helmholtz_bie_halfplane_vertical(cfg);
seconds=toc(tic_job);
assert(boundary.node_count==node_estimate && boundary.finite);
assert(all(isfinite(boundary.receiver_field)));
receipt=struct('status','COMPLETE','job_id',job_id, ...
    'wall_seconds',seconds,'node_count',boundary.node_count, ...
    'linear_residual',boundary.linear_residual, ...
    'offgrid_boundary_residual',boundary.offgrid_boundary_residual, ...
    'receiver_native_real',real(boundary.receiver_field), ...
    'receiver_native_imag',imag(boundary.receiver_field), ...
    'available_physical_bytes',available_physical_bytes, ...
    'estimated_increment_bytes',estimated_increment_bytes, ...
    'matlab_version',version);
save(fullfile(out,'result.mat'),'boundary','receipt','cfg','pre','src','-v7.3');
local_write_json(fullfile(out,'receipt.json'),receipt);
disp(receipt);
end

function receipt=local_model(task_root,available_physical_bytes)
out=fullfile(task_root,'jobs','model');
assert(isfolder(out) && ~isfile(fullfile(out,'result.mat')));
assert(isfile(fullfile(task_root,'jobs','smoke','result.mat')));
estimate=256*1024^2;
assert(estimate<=.7*available_physical_bytes,'Model RAM gate');
s=load(fullfile(task_root,'input','frozen_source.mat'),'src');
src=closeout2d_source_at_frequency(s.src,6000,1500);
N=src.n_original; W=192.1875;
x=(-W/2)+(0:N-1)'*W/N;
[eta,~]=closeout_s2b_surface(x,.01);
cfg=struct('frequency_hz',6000,'c0_mps',1500,'xw_m',W,'nx',N, ...
    'z_tx_m',100,'z_rx_m',3,'sigma_src_m',.3, ...
    'surface_reflect_coeff',-1,'step_m',.05,'x_rx_m',0, ...
    'reflection_model','model0_normal','source_field_override',ifft(src.coeff(:)), ...
    'eta_z_m',eta);
tic_model=tic;
pe=run_pe_physical_surface_model0_validation(cfg);
u_flat=-fk_single_receiver_incident(src,0,-3);
fk400=closeout_s2b_fullfk_analytic_exterior(src,.01,0,3,400,u_flat);
fk800=closeout_s2b_fullfk_analytic_exterior(src,.01,0,3,800,u_flat);
phase=exp(1i*2*pi*6000*97/1500);
H_pe=conj(pe.reflected_native(pe.ix_rx))*phase;
H_fk=conj(fk800.receiver_field)*phase;
local_pe=complex(-.105917460930372,-.103617351605008);
local_fk=complex(-.105917469807197,-.10361735723127);
pe_crossplatform=abs(H_pe-local_pe)/abs(local_pe);
fk_crossplatform=abs(H_fk-local_fk)/abs(local_fk);
fk_quadrature=abs(fk400.receiver_field-fk800.receiver_field)/ ...
    abs(fk800.receiver_field);
assert(pe_crossplatform<1e-10 && fk_crossplatform<1e-10, ...
    'Remote/local paired model compatibility gate');
receipt=struct('status','PASS','job_id','model', ...
    'matlab_version',version,'frequency_hz',6000,'amplitude_m',.01, ...
    'pe_crossplatform_relative',pe_crossplatform, ...
    'fullfk_crossplatform_relative',fk_crossplatform, ...
    'fullfk_400_800_relative',fk_quadrature, ...
    'pe_reflect_native_real',real(pe.reflected_native(pe.ix_rx)), ...
    'pe_reflect_native_imag',imag(pe.reflected_native(pe.ix_rx)), ...
    'fk_reflect_native_real',real(fk800.receiver_field), ...
    'fk_reflect_native_imag',imag(fk800.receiver_field), ...
    'wall_seconds',toc(tic_model), ...
    'available_physical_bytes',available_physical_bytes, ...
    'estimated_increment_bytes',estimate,'bie_boundary_solves',0);
save(fullfile(out,'result.mat'),'receipt','pe','fk400','fk800','src');
local_write_json(fullfile(out,'receipt.json'),receipt);
disp(receipt);
end

function receipt=local_smoke(task_root,available_physical_bytes)
out=fullfile(task_root,'jobs','smoke');
assert(isfolder(out) && ~isfile(fullfile(out,'result.mat')));
assert(contains(version,'R2023b'),'Remote MATLAB version gate');
estimate=128*1024^2;
assert(estimate<=.7*available_physical_bytes,'Smoke RAM gate');
s=load(fullfile(task_root,'input','frozen_source.mat'),'src');
src=closeout2d_source_at_frequency(s.src,6000,1500);
N=src.n_original; source=ifft(src.coeff(:));
cfg=struct('frequency_hz',6000,'c0_mps',1500,'xw_m',192.1875, ...
    'nx',N,'z_tx_m',100,'z_rx_m',3,'sigma_src_m',.3, ...
    'surface_reflect_coeff',-1,'step_m',.05,'x_rx_m',0, ...
    'reflection_model','model0_normal','source_field_override',source, ...
    'eta_z_m',zeros(N,1));
tic_smoke=tic; p=run_pe_physical_surface_model0_validation(cfg);
k=2*pi*6000/1500; phase=exp(1i*2*pi*6000*97/1500);
hd=conj(p.direct_receiver*exp(1i*k*97))*phase;
hr=conj(p.reflected_native(p.ix_rx))*phase;
local_hd=complex(.109231394104757,.106681137941899);
local_hr=complex(-.105932667223712,-.103601879791101);
direct_relative=abs(hd-local_hd)/abs(local_hd);
reflect_relative=abs(hr-local_hr)/abs(local_hr);
assert(direct_relative<1e-10 && reflect_relative<1e-10, ...
    'R2025b-to-R2023b flat source compatibility gate');
receipt=struct('status','PASS','job_id','smoke','matlab_version',version, ...
    'direct_relative_to_local',direct_relative, ...
    'reflect_relative_to_local',reflect_relative, ...
    'wall_seconds',toc(tic_smoke),'available_physical_bytes', ...
    available_physical_bytes,'estimated_increment_bytes',estimate, ...
    'bie_boundary_solves',0);
save(fullfile(out,'result.mat'),'receipt','hd','hr','src');
local_write_json(fullfile(out,'receipt.json'),receipt);
disp(receipt);
end

function slope=local_slope(x)
[~,slope]=closeout_s2b_surface(x,.01);
end

function local_write_json(path,value)
fid=fopen(path,'w'); assert(fid>0,'Cannot open JSON output');
guard=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(value));
end
