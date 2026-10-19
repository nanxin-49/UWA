function state = prepare_pe_bellhop_bie_xz_data(mode, run_dir)
%PREPARE_PE_BELLHOP_BIE_XZ_DATA One strong-surface 411x45 true field grid.
% init snapshots numeric inputs. stage writes ONE Bellhop grid job only.
% pebie evaluates saved PE starting field and saved BIE density, no BIE solve.
% finish imports a completed external Bellhop process; never launches it.
if nargin<1, mode='init'; end
if nargin<2, run_dir=''; end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
if strcmp(mode,'init')
    if isempty(run_dir)
        run_dir=fullfile(root,'results','visualization','pe_bellhop_bie_xz',char(datetime('now','Format','yyyyMMdd_HHmmss_SSS')));
    end
    assert(~isfolder(run_dir),'Refusing to overwrite existing run.');
    state=initialize(root,run_dir); checkpoint(state,run_dir); return
end
t=load(fullfile(run_dir,'state.mat'),'state'); state=t.state;
switch mode
    case 'stage'
        stage(state,run_dir);
    case 'pebie'
        pebie(run_dir);
        state.status='PE_BIE_EVALUATED';
    case 'finish'
        finish(run_dir); state.status='FIELD_DATA_COMPLETE';
    case 'status'
        disp(state);
    otherwise
        error('Use init, stage, pebie, finish, or status.');
end
checkpoint(state,run_dir);
end

function s=initialize(root,out)
for folder={'source_data','plot_data','metrics','figures','source_data/bellhop','source_data/workflow_code'}
    mkdir(fullfile(out,folder{1}));
end
base=fullfile(root,'results','visualization','pe_bellhop_bie_comparison','20261018_194331_523');
archive=fullfile(root,'results','archive','cl-verification-worktree-unmanaged-20260930','validation','pe_bellhop_helmholtz_bie_reference');
paths={fullfile(base,'source_data','stage0_validation.mat'), ...
    fullfile(archive,'R5_stronger_height','R5_stronger_height_validation.mat'), ...
    fullfile(archive,'R2_weak_rough','R2_weak_rough_convergence.mat'), ...
    fullfile(base,'source_data','job_08.mat'),fullfile(base,'source_data','job_02.mat'), ...
    fullfile(base,'plot_data','comparison_data.mat'),fullfile(base,'source_data','bellhop_internal_wall.exe')};
names={'stage0.mat','strong_reference.mat','weak_reference.mat','previous_strong_job.mat', ...
    'previous_flat_job.mat','previous_comparison.mat','bellhop_internal_wall.exe'};
hashes=strings(numel(paths),1);
for j=1:numel(paths)
    assert(isfile(paths{j})); hashes(j)=string(file_hash(paths{j}));
    copyfile(paths{j},fullfile(out,'source_data',names{j}));
    assert(strcmp(file_hash(fullfile(out,'source_data',names{j})),hashes(j)));
end
assert(strcmpi(hashes(end),'81de6a512037d63f06ede0537d5a591fde5e41b4ab9a95ff57f53f66917a5360'));
writetable(table(string(paths(:)),string(names(:)),hashes, ...
    'VariableNames',{'source','snapshot','sha256'}),fullfile(out,'source_data','input_manifest.csv'));
t=load(fullfile(out,'source_data','stage0.mat'),'validation'); s0=t.validation;
t=load(fullfile(out,'source_data','previous_comparison.mat'),'data'); old=t.data;
ref=fullfile(out,'source_data','strong_reference.mat');
pe_names={'surface_reflected_field','reflected_field','surface_incident_field','source_field','kx_rad_per_m'};
pe=struct;
for j=1:numel(pe_names)
    value=read_numeric(ref,['/validation/PE/' pe_names{j}]); pe.(pe_names{j})=value(:);
end
fine=read_boundary(ref,4); coarse=read_boundary(ref,3);
mk=logical(s0.footprint.m99.mask(:)); xfull=s0.x_m(:); x=xfull(mk); z=(.5:.125:6)';
assert(numel(x)==411 && numel(z)==45 && all(abs(diff(x)-.1953125)<1e-12));
assert(s0.config.frequency_hz==4000 && s0.config.z_tx_m==100 && s0.config.z_rx_m==3);
r5x=read_numeric(ref,'/validation/x_m'); r5x=r5x(:);
assert(max(abs(x-r5x(mk)))<1e-6);
flat=read_numeric(fullfile(out,'source_data','weak_reference.mat'),'/validation/flat_reference');
flat=flat(:);
k=2*pi*4000/1500; phase_restore=exp(1i*k*103);
flat_pe=s0.pe.reflected_field(:)*phase_restore;
flat_difference=norm(flat_pe(mk)-flat(mk))/norm(flat(mk)); assert(flat_difference<1e-11);
[~,axis_index]=min(abs(xfull)); bhflat=s0.bellhop_internal_high.field(:);
source_diagnostic=struct('PE_restored_flat_vs_BIE_l2',flat_difference, ...
    'BH_over_BIE_flat_center_abs',abs(bhflat(axis_index)/flat(axis_index)), ...
    'BH_over_BIE_flat_center_phase_rad',angle(bhflat(axis_index)/flat(axis_index)), ...
    'fitted_scalar_applied',false,'phase_alignment_applied',false);
inputs=struct('x_m',x,'z_m',z,'xfull_m',xfull,'original_M99',mk,'pe',pe, ...
    'boundary_fine',fine,'boundary_coarse',coarse,'source_pattern',s0.source_pattern, ...
    'old_BIE_ratio',old.cases(3).G_BIE,'BIE_flat_native',flat, ...
    'old_PE_ratio',old.cases(3).G_PE,'PE_flat_envelope',s0.pe.reflected_field(:), ...
    'source_diagnostic',source_diagnostic,'config',s0.config);
assert(norm(conj(fine.receiver_field(mk)./flat(mk))-old.cases(3).G_BIE(mk))<1e-12);
save(fullfile(out,'source_data','numeric_inputs.mat'),'inputs','-v7.3');
fid=fopen(fullfile(out,'source_data','source_conventions.json'),'w'); cl=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(source_diagnostic,PrettyPrint=true)); clear cl
s=struct('status','INITIALIZED','run_dir',out,'grid_size',[411 45], ...
    'new_Bellhop_jobs_allowed',1,'new_BIE_boundary_solves',0,'amplitude_m',.2, ...
    'frequency_hz',4000,'created',char(datetime('now')));
fprintf('XZ run: %s\n',out); disp(source_diagnostic);
end

function stage(s,out)
assert(s.new_Bellhop_jobs_allowed==1);
folder=fullfile(out,'source_data','bellhop'); base=fullfile(folder,'field');
assert(~isfile([base '.env']) && ~isfile(fullfile(folder,'console.log')),'Do not restage or duplicate a job.');
t=load(fullfile(out,'source_data','numeric_inputs.mat'),'inputs'); d=t.inputs;
xx=linspace(-80,80,4097)'; a=abs(xx); r=(a-42)/8; w=double(a<=42); mid=a>42&a<50;
w(mid)=1-10*r(mid).^3+15*r(mid).^4-6*r(mid).^5;
eta=.2*sin(.1*xx).*w;
cfg=struct('bellhop_exe',fullfile(out,'source_data','bellhop_internal_wall.exe'), ...
    'case_root',base,'run_type','C','source_geometry','X','frequency_hz',4000,'c0_mps',1500, ...
    'source_depth_m',0,'receiver_depths_m',sort(-d.x_m).', ...
    'receiver_ranges_m',(100+d.z_m).','mapped_receiver_range_m',103, ...
    'beam_count',10001,'angle_limits_deg',[-30 30],'step_m',.05,'domain_half_depth_m',1000, ...
    'source_pattern_angles_deg',d.source_pattern.angles_deg,'source_pattern_level_db',d.source_pattern.level_db, ...
    'wall_r0_m',100,'wall_seed',0,'wall_profile_r_m',100-eta,'wall_profile_z_m',xx,'reuse_existing',true);
write_bellhop_unfolded_gaussian_env_vertical(base,cfg);
fid=fopen([base '.iwpm'],'w'); cl=onCleanup(@()fclose(fid));
fprintf(fid,'%.17g\n%.17g\n%d\n%d\n',100,103,0,numel(xx)); fprintf(fid,'%.17g %.17g\n',[100-eta xx].'); clear cl
save(fullfile(out,'source_data','bellhop_config.mat'),'cfg','-v7.3');
disp('Staged one 10001-beam 411x45 receiver-grid job; no process launched.');
end

function pebie(out)
done=fullfile(out,'plot_data','pe_bie_fields.mat');
if isfile(done), disp('PE/BIE field evaluations already checkpointed.'); return; end
t=load(fullfile(out,'source_data','numeric_inputs.mat'),'inputs'); d=t.inputs;
k=2*pi*4000/1500; kz=sqrt(complex(k^2-d.pe.kx_rad_per_m(:).^2,0));
spec=fft(d.pe.surface_reflected_field(:)); pe=complex(zeros(411,45));
for j=1:45
    full=ifft(spec.*exp(1i*d.z_m(j)*(kz-k))); pe(:,j)=full(d.original_M99);
end
[X,Z]=ndgrid(d.x_m,d.z_m); clk=tic;
line=evaluate_saved_halfplane_bie_field(d.boundary_fine,d.x_m,3+zeros(size(d.x_m)));
line_error=norm(line-d.boundary_fine.receiver_field(d.original_M99))/norm(line);
assert(line_error<1e-12,'Reconstructed BIE quadrature differs from saved receiver line.');
pe_error=norm(pe(:,d.z_m==3)-d.pe.reflected_field(d.original_M99))/norm(d.pe.reflected_field(d.original_M99));
assert(pe_error<1e-11,'Saved PE starting field does not reproduce old receiver line.');
[PX,PZ]=ndgrid(d.x_m(1:41:end),[.5 1 3 6]);
fine=evaluate_saved_halfplane_bie_field(d.boundary_fine,PX,PZ);
coarse=evaluate_saved_halfplane_bie_field(d.boundary_coarse,PX,PZ);
probe_l2=norm(fine(:)-coarse(:))/norm(fine(:));
assert(probe_l2<1e-4,'Saved BIE density levels differ too much on new-grid probes.');
bie=evaluate_saved_halfplane_bie_field(d.boundary_fine,X,Z); seconds=toc(clk);
fields=struct('x_m',d.x_m,'z_m',d.z_m,'PE_envelope',pe, ...
    'PE_carrier_restored',pe.*exp(1i*k*(100+d.z_m.')), ...
    'BIE_native',bie,'source_diagnostic',d.source_diagnostic, ...
    'checks',struct('PE_z3_closure',pe_error,'BIE_z3_closure',line_error, ...
    'BIE_12_vs_14_ppw_probe_l2',probe_l2,'probe_count',numel(PX),'BIE_field_evaluation_seconds',seconds), ...
    'new_BIE_boundary_solves',0);
save(done,'fields','-v7.3'); disp(fields.checks);
writetable(table(PX(:),PZ(:),real(fine(:)),imag(fine(:)),real(coarse(:)),imag(coarse(:)), ...
    'VariableNames',{'x_m','z_m','fine_real','fine_imag','coarse_real','coarse_imag'}),fullfile(out,'metrics','BIE_density_level_probes.csv'));
end

function finish(out)
receipt=fullfile(out,'source_data','bellhop','process_result.json');
assert(isfile(receipt),'External Bellhop process not yet completed; do not duplicate it.');
p=jsondecode(fileread(receipt)); assert(p.exit_code==0);
t=load(fullfile(out,'source_data','bellhop_config.mat'),'cfg'); cfg=t.cfg;
t=load(fullfile(out,'plot_data','pe_bie_fields.mat'),'fields'); data=t.fields;
cwd=pwd; run=run_bellhop_internal_pm_wall_poc_vertical(cfg); assert(strcmp(cwd,pwd));
data.BH_raw_SHD=complex(zeros(411,45)); max_error=0;
for j=1:45
    for i=1:411
        [value,m]=select_bellhop_shd_pressure_at_range_vertical(run.data,100+data.z_m(j),-data.x_m(i),1e-6);
        data.BH_raw_SHD(i,j)=value; max_error=max(max_error,max(abs([m.actual_range_m-100-data.z_m(j),m.actual_depth_m+data.x_m(i)])));
    end
end
data.BH_native=conj(data.BH_raw_SHD); data.BH_solver_wall_seconds=p.wall_seconds;
d=run.diagnostics;
geometry=struct('wall_residual_m',max(abs(d.wall_residual)), ...
    'phase_jump_error_rad',max(abs(d.phase_delta-pi)),'p_rotation_error',max(abs(d.p_rot_error)), ...
    'q_rotation_error',max(abs(d.q_rot_error)),'min_post_dr_m',min(d.min_post_dr),'coordinate_error_m',max_error);
assert(geometry.wall_residual_m<1e-9 && geometry.phase_jump_error_rad<1e-10 ...
    && geometry.p_rotation_error<1e-12 && geometry.q_rotation_error<1e-12 ...
    && geometry.min_post_dr_m>0 && max_error<=1e-6 && all(isfinite(d{:,:}),'all'));
old=load(fullfile(out,'source_data','previous_strong_job.mat'),'job_result');
t=load(fullfile(out,'source_data','numeric_inputs.mat'),'inputs'); mk=t.inputs.original_M99;
data.checks.BH_z3_closure=norm(data.BH_native(:,data.z_m==3)-old.job_result.native_field(mk))/norm(old.job_result.native_field(mk));
assert(data.checks.BH_z3_closure<2e-5,'New Bellhop grid differs from old receiver line.');
data.geometry=geometry; data.schema_version='1.0';
data.methods={'PE_carrier_restored','BH_native','BIE_native'};
data.conventions=struct('PE','saved reduced envelope times exp(+i*k0*(100+z)); no fitted phase', ...
    'BH','conj(raw SHD), existing native exp(-i*omega*t) convention; native source units', ...
    'BIE','saved density Green evaluation, native exp(-i*omega*t); source peak=1', ...
    'cross_method_calibration','No fitted amplitude/phase. Bellhop source scale and phase reference differ; no absolute pressure or TL claims.');
data.phase_rad=struct; data.magnitude=struct;
[X,Z]=ndgrid(data.x_m,data.z_m); T=table(X(:),Z(:),'VariableNames',{'x_m','z_m'});
T.PE_envelope_real=real(data.PE_envelope(:)); T.PE_envelope_imag=imag(data.PE_envelope(:));
T.BH_raw_SHD_real=real(data.BH_raw_SHD(:)); T.BH_raw_SHD_imag=imag(data.BH_raw_SHD(:));
for j=1:3
    name=data.methods{j}; u=data.(name); assert(all(isfinite(u),'all'));
    data.magnitude.(name)=abs(u); data.phase_rad.(name)=angle(u);
    T.([name '_real'])=real(u(:)); T.([name '_imag'])=imag(u(:));
    T.([name '_magnitude'])=abs(u(:)); T.([name '_phase_rad'])=angle(u(:));
end
data.new_Bellhop_jobs=1; data.new_BIE_boundary_solves=0;
data.new_PE_full_runs=0; data.PE_saved_start_depth_evaluations=45;
save(fullfile(out,'plot_data','xz_data.mat'),'data','-v7.3');
writetable(T,fullfile(out,'plot_data','xz_fields.csv'));
writetable(d,fullfile(out,'metrics','wall_diagnostics.csv'));
fid=fopen(fullfile(out,'metrics','field_checks.json'),'w'); cl=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(struct('checks',data.checks,'geometry',geometry),PrettyPrint=true));
disp(data.checks);
end

function b=read_boundary(file,index)
fid=H5F.open(file,'H5F_ACC_RDONLY','H5P_DEFAULT'); cl=onCleanup(@()H5F.close(fid));
did=H5D.open(fid,'/validation/spatial_runs'); cdid=onCleanup(@()H5D.close(did)); refs=H5D.read(did);
obj=H5R.dereference(did,'H5R_OBJECT',refs(:,index)); group=H5I.get_name(obj); H5O.close(obj);
names={'node_x_m','node_z_m','node_normal','window','density','half_width_m','halfplane_h_m', ...
    'panel_order','panel_count','k_radpm','receiver_field','receiver_x_m'}; b=struct;
for j=1:numel(names)
    value=read_numeric(file,[group '/' names{j}]);
    if strcmp(names{j},'node_normal'), b.(names{j})=value; else, b.(names{j})=value(:); end
end
b.hdf5_group=group; b.spatial_run_index=index;
end

function value=read_numeric(file,path)
value=h5read(file,path);
if isstruct(value), assert(all(isfield(value,{'real','imag'}))); value=complex(value.real,value.imag); end
assert(isnumeric(value) && all(isfinite(value),'all'));
end

function checkpoint(state,out)
save(fullfile(out,'state.mat'),'state','-v7.3');
fid=fopen(fullfile(out,'status.json'),'w'); cl=onCleanup(@()fclose(fid)); fprintf(fid,'%s',jsonencode(state,PrettyPrint=true));
end

function value=file_hash(path)
fid=fopen(path,'rb'); assert(fid>=0); cl=onCleanup(@()fclose(fid)); md=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid), md.update(fread(fid,1048576,'*uint8')); end
value=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
