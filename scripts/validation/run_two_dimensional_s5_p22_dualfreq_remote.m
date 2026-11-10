function receipt = run_two_dimensional_s5_p22_dualfreq_remote(task_root,frequency_hz,job_id,available_physical_bytes)
%RUN_TWO_DIMENSIONAL_S5_P22_DUALFREQ_REMOTE One isolated P22 pilot job.
% Validation-only code. The caller owns remote authorization, job creation,
% process wall-clock/disk checks, full logs, SHA manifests and exit status.
arguments
    task_root (1,:) char
    frequency_hz (1,1) double {mustBeFinite,mustBePositive}
    job_id (1,:) char
    available_physical_bytes (1,1) double {mustBeFinite,mustBePositive}
end
assert(any(frequency_hz==[4062.5,4031.25]),'Unapproved pilot frequency');
assert(isfolder(task_root) && isfolder(fullfile(task_root,'code')));
addpath(fullfile(task_root,'code'));
if frequency_hz==4062.5
    frequency_id='f4062p5';
else
    frequency_id='f4031p25';
end
source_path=fullfile(task_root,'input','frozen_source.mat');
assert(isfile(source_path),'Frozen source missing');
s=load(source_path,'src');
src=closeout2d_source_at_frequency(s.src,frequency_hz,1500);
assert(src.n_original==984 && src.z_tx_m==100 && src.x0_m==-96.09375);
assert(max(abs(src.coeff(:)-s.src.coeff(:)))==0 && ...
    max(abs(src.kx(:)-s.src.kx(:)))==0,'Frozen source changed');
tau_ref=96.875/1500;
if strcmp(job_id,'model')
    receipt=local_model(task_root,frequency_id,frequency_hz,src,tau_ref,available_physical_bytes);
    return
end
ids={'main96','self64','space14','window60','window70'};
ix=find(strcmp(ids,job_id));
assert(isscalar(ix),'Unknown reference layer');
out=fullfile(task_root,'jobs',frequency_id,job_id);
assert(isfolder(out) && ~isfile(fullfile(out,'result.mat')),'Job directory missing or already used');
assert(isfile(fullfile(task_root,'jobs',frequency_id,'model','result.mat')), ...
    'Same-frequency PE/FullFK model result required');
for jj=1:ix-1
    assert(isfile(fullfile(task_root,'jobs',frequency_id,ids{jj},'result.mat')), ...
        'Earlier BIE reference layer missing');
end
if frequency_hz==4031.25
    prediction_path=fullfile(task_root,'prediction','prediction.json');
    assert(isfile(prediction_path),'Half-step prediction must be frozen first');
    prediction=jsondecode(fileread(prediction_path));
    prediction_sha_path=fullfile(task_root,'prediction','prediction_sha256.txt');
    assert(isfile(prediction_sha_path) && ...
        strcmp(strtrim(fileread(prediction_sha_path)),local_sha256(prediction_path)), ...
        'Half-step prediction SHA mismatch');
    assert(strcmp(prediction.status,'FROZEN_BEFORE_HALFSTEP_BIE') && ...
        prediction.frequency_hz==4031.25 && ...
        isfinite(prediction.native_real) && isfinite(prediction.native_imag), ...
        'Half-step prediction contract invalid');
end
L=[60 60 60 60 70]; ppw=[12 12 14 10 10];
qa=[96 64 96 96 96];
lambda=1500/frequency_hz;
panels=ceil(2*L(ix)/(lambda/ppw(ix)*8));
node_estimate=8*panels;
estimated_increment_bytes=64*node_estimate^2+128*1024^2;
assert(estimated_increment_bytes<=.7*available_physical_bytes, ...
    '30 percent physical RAM reserve failed');
cfg=struct('frequency_hz',frequency_hz,'c0_mps',1500, ...
    'half_width_m',L(ix),'points_per_wavelength',ppw(ix), ...
    'panel_order',8,'self_quadrature_order',qa(ix), ...
    'halfplane_h_m',-2,'window_plateau_ratio',.85, ...
    'eta_fn',@(x) closeout_s2b_surface(x,.01), ...
    'eta_prime_fn',@(x) local_slope(x), ...
    'incident_fn',@(x,z) fk_single_receiver_incident(src,x,z), ...
    'receiver_x_m',0,'receiver_z_m',3);
pre=struct('job_id',job_id,'frequency_hz',frequency_hz, ...
    'receiver_x_m',0,'receiver_z_m',3,'surface_A_m',.01, ...
    'surface_K_radpm',.1,'surface_phi_rad',0,'tau_ref_s',tau_ref, ...
    'node_estimate',node_estimate, ...
    'estimated_increment_bytes',estimated_increment_bytes, ...
    'available_physical_bytes',available_physical_bytes, ...
    'ram_reserve_fraction',.3,'source_modes',984);
local_write_json(fullfile(out,'preflight.json'),pre);
tx_x=src.x0_m+(0:983)'*192.1875/984;
tx_expected=ifft(src.coeff(:));
tx_actual=fk_single_receiver_incident(src,tx_x,100+zeros(984,1));
assert(norm(tx_actual-tx_expected)/norm(tx_expected)<1e-10,'Source grid mismatch');
tic_job=tic;
boundary=solve_helmholtz_bie_halfplane_vertical(cfg);
seconds=toc(tic_job);
assert(boundary.node_count==node_estimate && boundary.finite);
assert(all(isfinite(boundary.receiver_field)));
H_reflect_native=boundary.receiver_field;
H_direct_native=fk_single_receiver_incident(src,0,3);
H_total_native=H_direct_native+H_reflect_native;
phase=exp(1i*2*pi*frequency_hz*tau_ref);
H_direct_dsp=conj(H_direct_native)*phase;
H_reflect_dsp=conj(H_reflect_native)*phase;
H_total_dsp=conj(H_total_native)*phase;
assert(abs(H_total_dsp-H_direct_dsp-H_reflect_dsp) ...
    <=1e-12*max(1,abs(H_total_dsp)),'Complex component closure failed');
receipt=struct('status','COMPLETE_UNAUDITED','job_id',job_id, ...
    'frequency_hz',frequency_hz,'wall_seconds',seconds, ...
    'node_count',boundary.node_count,'linear_residual',boundary.linear_residual, ...
    'offgrid_boundary_residual',boundary.offgrid_boundary_residual, ...
    'receiver_native_real',real(H_reflect_native), ...
    'receiver_native_imag',imag(H_reflect_native), ...
    'available_physical_bytes',available_physical_bytes, ...
    'estimated_increment_bytes',estimated_increment_bytes, ...
    'matlab_version',version);
save(fullfile(out,'result.mat'),'boundary','receipt','cfg','pre','src', ...
    'H_direct_native','H_reflect_native','H_total_native', ...
    'H_direct_dsp','H_reflect_dsp','H_total_dsp','tau_ref','-v7.3');
local_write_json(fullfile(out,'receipt.json'),receipt);
disp(receipt);
end

function receipt=local_model(task_root,frequency_id,frequency_hz,src,tau_ref,available_physical_bytes)
out=fullfile(task_root,'jobs',frequency_id,'model');
assert(isfolder(out) && ~isfile(fullfile(out,'result.mat')),'Job directory missing or already used');
estimated_increment_bytes=256*1024^2;
assert(estimated_increment_bytes<=.7*available_physical_bytes,'Model RAM reserve failed');
N=984; W=192.1875;
x=(-W/2)+(0:N-1)'*W/N;
[eta,~]=closeout_s2b_surface(x,.01);
cfg=struct('frequency_hz',frequency_hz,'c0_mps',1500,'xw_m',W,'nx',N, ...
    'z_tx_m',100,'z_rx_m',3,'sigma_src_m',.3, ...
    'surface_reflect_coeff',-1,'step_m',.05,'x_rx_m',0, ...
    'reflection_model','model0_normal','source_field_override',ifft(src.coeff(:)), ...
    'eta_z_m',eta);
tic_job=tic;
pe=run_pe_physical_surface_model0_validation(cfg);
H_flat_native=-fk_single_receiver_incident(src,0,-3);
fk400=closeout_s2b_fullfk_analytic_exterior(src,.01,0,3,400,H_flat_native);
fk800=closeout_s2b_fullfk_analytic_exterior(src,.01,0,3,800,H_flat_native);
seconds=toc(tic_job);
H_direct_native=fk_single_receiver_incident(src,0,3);
H_pe_reflect_native=pe.reflected_native(pe.ix_rx);
H_fk_reflect_native=fk800.receiver_field;
H_pe_total_native=H_direct_native+H_pe_reflect_native;
H_fk_total_native=H_direct_native+H_fk_reflect_native;
phase=exp(1i*2*pi*frequency_hz*tau_ref);
H_pe_dsp=conj(H_pe_total_native)*phase;
H_fk_dsp=conj(H_fk_total_native)*phase;
fk_quadrature_relative=abs(fk400.receiver_field-fk800.receiver_field)/abs(fk800.receiver_field);
assert(isfinite(fk_quadrature_relative));
receipt=struct('status','COMPLETE_UNAUDITED','job_id','model', ...
    'frequency_hz',frequency_hz,'wall_seconds',seconds, ...
    'pe_reflect_native_real',real(H_pe_reflect_native), ...
    'pe_reflect_native_imag',imag(H_pe_reflect_native), ...
    'fk_reflect_native_real',real(H_fk_reflect_native), ...
    'fk_reflect_native_imag',imag(H_fk_reflect_native), ...
    'fk_400_800_relative',fk_quadrature_relative, ...
    'available_physical_bytes',available_physical_bytes, ...
    'estimated_increment_bytes',estimated_increment_bytes, ...
    'matlab_version',version);
save(fullfile(out,'result.mat'),'receipt','cfg','pe','fk400','fk800','src', ...
    'H_flat_native','H_direct_native','H_pe_reflect_native','H_fk_reflect_native', ...
    'H_pe_total_native','H_fk_total_native','H_pe_dsp','H_fk_dsp','tau_ref','-v7.3');
local_write_json(fullfile(out,'receipt.json'),receipt);
disp(receipt);
end

function slope=local_slope(x)
[~,slope]=closeout_s2b_surface(x,.01);
end

function local_write_json(path,value)
fid=fopen(path,'w');
assert(fid>0,'Cannot open JSON output');
guard=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(value));
end

function h=local_sha256(path)
md=java.security.MessageDigest.getInstance('SHA-256');
fid=fopen(path,'rb'); assert(fid>=0);
guard=onCleanup(@() fclose(fid)); %#ok<NASGU>
while ~feof(fid)
    md.update(fread(fid,1048576,'*uint8'));
end
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
