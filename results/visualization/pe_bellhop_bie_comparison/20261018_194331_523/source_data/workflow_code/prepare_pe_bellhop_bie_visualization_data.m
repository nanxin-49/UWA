function state = prepare_pe_bellhop_bie_visualization_data(mode, run_dir, bellhop_exe)
%PREPARE_PE_BELLHOP_BIE_VISUALIZATION_DATA Reuse PE/BIE; checkpoint eight BH jobs.
% init creates a new independent run. next executes one missing Bellhop job.
% finalize computes plot data without invoking PE or BIE; status only reads.
if nargin<1, mode='init'; end
if nargin<2, run_dir=''; end
if nargin<3, bellhop_exe=''; end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
if strcmp(mode,'init')
    if isempty(run_dir)
        run_dir=fullfile(root,'results','visualization','pe_bellhop_bie_comparison', ...
            char(datetime('now','Format','yyyyMMdd_HHmmss_SSS')));
    end
    assert(~isfolder(run_dir),'Refusing to overwrite an existing run directory.');
    state=initialize(root,run_dir);
    checkpoint(state,run_dir);
    return
end
assert(isfolder(run_dir),'Supply the exact run directory returned by init.');
q=load(fullfile(run_dir,'state.mat'),'state'); state=q.state;
switch mode
    case 'stage'
        stage_inputs(state,run_dir,bellhop_exe);
    case 'next'
        state=next_job(state,run_dir,bellhop_exe); checkpoint(state,run_dir);
    case 'finalize'
        assert(all([state.jobs.complete]),'Finish all eight Bellhop jobs first.');
        state=finalize(state,run_dir); checkpoint(state,run_dir);
    case 'status'
        disp(struct2table(state.jobs)); disp(state.status);
    otherwise
        error('Unknown mode. Use init, next, finalize, or status.');
end
end

function s=initialize(root,out)
for folder={'source_data','plot_data','metrics','figures','figures/bellhop_bie','figures/three_way'}
    mkdir(fullfile(out,folder{1}));
end
paths={'results/validation/pe_bellhop_controlled_comparison/stage0/stage0_validation.mat', ...
    'results/validation/pe_bie_error_decomposition/pe_bie_error_decomposition_validation.mat', ...
    'reports/pe_bellhop_helmholtz_bie_final_report.md', ...
    'reports/pe_bellhop_helmholtz_bie_R3_weak_three_way_closure_report.md', ...
    'reports/pe_bellhop_helmholtz_bie_R4_region_II_adjudication_report.md', ...
    'reports/pe_bellhop_helmholtz_bie_R5_stronger_height_report.md', ...
    'scripts/validation/support/bellhop_internal_pm_wall_poc/bellhop.f90', ...
    'scripts/validation/support/bellhop_internal_pm_wall_poc/Step.f90', ...
    'scripts/validation/support/bellhop_internal_pm_wall_poc/build_bellhop_internal_pm_wall_poc.ps1', ...
    'results/validation/bellhop_internal_pm_wall_poc/build_manifest.json'};
hashes=strings(numel(paths),1); names=hashes;
for j=1:numel(paths)
    src=fullfile(root,paths{j}); assert(isfile(src),'Missing active source: %s',src);
    [~,nm,ex]=fileparts(src); names(j)=string([nm ex]);
    dest=fullfile(out,'source_data',char(names(j))); copyfile(src,dest);
    hashes(j)=string(file_hash(src)); assert(strcmp(file_hash(dest),hashes(j)));
end
writetable(table(string(paths(:)),names,hashes,'VariableNames',{'source','snapshot','sha256'}), ...
    fullfile(out,'source_data','input_manifest.csv'));
a=load(fullfile(out,'source_data','stage0_validation.mat'),'validation'); s0=a.validation;
b=load(fullfile(out,'source_data','pe_bie_error_decomposition_validation.mat'),'validation'); old=b.validation;
cfg=s0.config;
assert(cfg.frequency_hz==4000 && cfg.c0_mps==1500 && cfg.z_tx_m==100 && cfg.z_rx_m==3);
assert(cfg.sigma_src_m==.3 && cfg.nx==984 && cfg.xw_m==192.1875);
assert(strcmp(cfg.source_geometry,'X') && cfg.pe_step_m==.05 && cfg.bellhop_step_m==.05);
assert(cfg.mapped_receiver_range_m==103 && cfg.wall_r0_m==100 && cfg.phase_sign==-1);
x=s0.x_m(:); mask=logical(s0.footprint.m99.mask(:)); weights=s0.footprint.energy_weights(:);
assert(numel(x)==984 && numel(mask)==984 && all(isfinite(weights)) && all(weights>=0));
labels={'weak_low_K','region_II_low_K','strong_height_low_K'}; amplitudes=[.01 .05 .2];
cases=struct([]);
for j=1:3
    ix=find(strcmp(string({old.details.case_name}),labels{j})); assert(isscalar(ix)); d=old.details(ix);
    assert(max(abs(x-d.x_m(:)))<=1e-6,'Derived receiver grid differs from Stage-0.');
    eta=surface(x,amplitudes(j)); assert(max(abs(eta-d.surface_elevation_m(:)))<1e-12);
    pe=d.G_Model0(:); bie=d.G_BIE(:); assert(numel(pe)==984 && numel(bie)==984);
    assert(all(isfinite(pe)) && all(isfinite(bie)),'Invalid saved source curve.');
    cases(j).name=labels{j}; cases(j).amplitude_m=amplitudes(j);
    cases(j).eta_m=eta; cases(j).G_PE=pe; cases(j).G_BIE=bie;
    cases(j).G_BIE_native=conj(bie); cases(j).BIE_source=char(d.source_file);
end
source_coeff=fft(s0.pe.source_field(:)); kx=s0.pe.kx_rad_per_m(:);
jobs=repmat(struct('case_index',0,'beams',0,'complete',false,'file',''),1,8); n=0;
for ci=0:3
    for beams=[5001 10001]
        n=n+1; jobs(n)=struct('case_index',ci,'beams',beams,'complete',false, ...
            'file',sprintf('job_%02d.mat',n));
    end
end
s=struct('schema_version','1.0','status','READY_FOR_BELLHOP','created',char(datetime('now')), ...
    'matlab_version',version,'config',cfg,'x_m',x,'mask',mask,'weights',weights, ...
    'cases',cases,'source_pattern',s0.source_pattern,'source_coefficients',source_coeff, ...
    'source_kx',kx,'source_initial',s0.pe.source_field(:), ...
    'saved_flat_BH_native',s0.bellhop_internal_high.field(:),'jobs',jobs, ...
    'source_note','PE/BIE reused from active numeric derivative MAT; original R3-R5 directories absent', ...
    'convention','G_PE unchanged; G_BIE stored comparison; G_BH_comparison=conj(G_BH_native)');
numeric_inputs=rmfield(s,{'jobs','saved_flat_BH_native'});
save(fullfile(out,'source_data','numeric_inputs.mat'),'numeric_inputs','-v7.3');
fprintf('Created independent run: %s\n',out);
end

function s=next_job(s,out,exe)
ix=find(~[s.jobs.complete],1);
if isempty(ix), disp('All eight jobs already complete.'); return; end
% Import any finished staged job before considering a new invocation.
for candidate=find(~[s.jobs.complete])
    jj=s.jobs(candidate); folder=fullfile(out,'source_data','bellhop',sprintf('case_%d_B%d',jj.case_index,jj.beams));
    receipt=fullfile(folder,'process_result.json');
    if isfile(receipt)
        p=jsondecode(fileread(receipt)); assert(p.exit_code==0,'Staged Bellhop process failed.');
        ix=candidate; break
    end
end
assert(isfile(exe),'Supply a provenance-checked internal-wall executable.');
stamped=fullfile(out,'source_data','bellhop_internal_wall.exe');
if ~isfile(stamped), copyfile(exe,stamped); end
assert(strcmp(file_hash(exe),file_hash(stamped)),'Executable changed between checkpoint jobs.');
if ~isfield(s,'executable_sha256'), s.executable_sha256=file_hash(stamped); end
job=s.jobs(ix); done=fullfile(out,'source_data',job.file);
if isfile(done)
    q=load(done,'job_result'); assert(q.job_result.case_index==job.case_index && q.job_result.beams==job.beams);
    assert(strcmp(q.job_result.executable_sha256,s.executable_sha256));
    s.jobs(ix).complete=true; return
end
cfg=s.config; ci=job.case_index; amp=0; if ci>0, amp=s.cases(ci).amplitude_m; end
xx=linspace(-80,80,4097)'; wall=100-surface(xx,amp);
[zsort,ord]=sort(-s.x_m); invord=zeros(size(ord)); invord(ord)=(1:numel(ord))';
case_dir=fullfile(out,'source_data','bellhop',sprintf('case_%d_B%d',ci,job.beams));
if ~isfolder(case_dir), mkdir(case_dir); end
if isfile(fullfile(case_dir,'console.log'))
    assert(isfile(fullfile(case_dir,'process_result.json')), ...
        'Staged process may still be running: do not launch a duplicate.');
end
c=struct('bellhop_exe',stamped,'case_root',fullfile(case_dir,'field'), ...
    'run_type',cfg.run_type,'source_geometry','X','frequency_hz',4000,'c0_mps',1500, ...
    'source_depth_m',0,'receiver_depths_m',zsort.', ...
    'receiver_ranges_m',[cfg.guard_receiver_range_m 103],'mapped_receiver_range_m',103, ...
    'beam_count',job.beams,'angle_limits_deg',cfg.angle_limits_deg,'step_m',.05, ...
    'domain_half_depth_m',cfg.domain_half_depth_m, ...
    'source_pattern_angles_deg',s.source_pattern.angles_deg, ...
    'source_pattern_level_db',s.source_pattern.level_db, ...
    'wall_r0_m',100,'wall_seed',cfg.wall_seed,'wall_profile_r_m',wall, ...
    'wall_profile_z_m',xx,'reuse_existing',true);
fprintf('Bellhop job %d/8: A=%.3g, beams=%d\n',ix,amp,job.beams);
clk=tic; original_cwd=pwd; run=run_bellhop_internal_pm_wall_poc_vertical(c); wall_seconds=toc(clk);
assert(strcmp(pwd,original_cwd),'Bellhop runner changed the caller cwd.');
raw=complex(zeros(numel(s.x_m),1)); coordinate_error=0;
for j=1:numel(raw)
    [raw(j),meta]=select_bellhop_shd_pressure_at_range_vertical(run.data,103,zsort(j),1e-6);
    coordinate_error=max(coordinate_error,max(abs([meta.actual_range_m-103,meta.actual_depth_m-zsort(j)])));
end
raw=raw(invord); native=conj(raw); d=run.diagnostics;
geometry=struct('wall_residual_m',max(abs(d.wall_residual)), ...
    'phase_jump_error_rad',max(abs(d.phase_delta-pi)), ...
    'p_rotation_error',max(abs(d.p_rot_error)),'q_rotation_error',max(abs(d.q_rot_error)), ...
    'min_post_dr_m',min(d.min_post_dr),'coordinate_error_m',coordinate_error);
geometry.pass=all(isfinite(d{:,:}),'all') && geometry.wall_residual_m<=1e-9 ...
    && geometry.phase_jump_error_rad<=1e-10 && geometry.p_rotation_error<=1e-12 ...
    && geometry.q_rotation_error<=1e-12 && geometry.min_post_dr_m>0 && coordinate_error<=1e-6;
assert(geometry.pass,'Bellhop geometry or pressure-release diagnostics failed.');
assert(contains(fileread(run.files.prt),'Line source','IgnoreCase',true),'PRT does not identify X line source.');
assert(all(isfinite(raw)),'Bellhop returned invalid complex samples.');
job_result=struct('case_index',ci,'beams',job.beams,'raw_SHD',raw,'native_field',native, ...
    'geometry',geometry,'wall_seconds',wall_seconds,'config',c, ...
    'executable_sha256',s.executable_sha256);
save(done,'job_result','-v7.3'); writetable(d,fullfile(case_dir,'wall_diagnostics.csv'));
s.jobs(ix).complete=true; s.status='BELLHOP_PARTIAL';
if all([s.jobs.complete]), s.status='READY_TO_FINALIZE'; end
end

function s=finalize(s,out)
data=struct('schema_version','1.0','config',s.config,'x_m',s.x_m,'mask',s.mask, ...
    'weights',s.weights,'cases',s.cases,'source_initial',s.source_initial, ...
    'source_coefficients',s.source_coefficients,'source_kx',s.source_kx, ...
    'source_pattern',s.source_pattern,'source_note',s.source_note,'convention',s.convention);
allmetrics=struct([]); sensitivity=struct([]); geometry=struct([]); timings=struct([]);
for j=1:8
    q=load(fullfile(out,'source_data',s.jobs(j).file),'job_result'); runs(j)=q.job_result; %#ok<AGROW>
    gg=q.job_result.geometry; gg.case_index=q.job_result.case_index; gg.beams=q.job_result.beams;
    geometry=[geometry;gg]; %#ok<AGROW>
    jr=q.job_result;
    receipt=fullfile(out,'source_data','bellhop',sprintf('case_%d_B%d',jr.case_index,jr.beams),'process_result.json');
    seconds=jr.wall_seconds; timing_source='MATLAB runner wall clock';
    if isfile(receipt)
        pr=jsondecode(fileread(receipt)); seconds=pr.wall_seconds; timing_source='external process receipt';
    end
    timings=[timings;struct('case_index',jr.case_index,'beams',jr.beams, ...
        'solver_wall_seconds',seconds,'source',timing_source)]; %#ok<AGROW>
end
flat=runs([runs.case_index]==0 & [runs.beams]==10001).native_field;
data.flat_control=metric(flat,s.saved_flat_BH_native,data.weights,data.mask);
assert(data.flat_control.complex_l2<1e-12,'New flat control differs from saved Stage-0.');
data.flat_control.reference='saved Stage-0 Bellhop native flat, 10001 beams';
writetable(struct2table(data.flat_control),fullfile(out,'metrics','flat_control.csv'));
writetable(struct2table(timings),fullfile(out,'metrics','solver_wall_times.csv'));
for ci=1:3
    ratios=cell(1,2);
    for bi=1:2
        beams=[5001 10001]; br=beams(bi);
        flat=runs([runs.case_index]==0 & [runs.beams]==br).native_field;
        rough=runs([runs.case_index]==ci & [runs.beams]==br).native_field;
        valid_denominator=isfinite(flat)&isfinite(rough)&abs(flat)>0;
        assert(all(valid_denominator(data.mask)),'Invalid flat denominator inside fixed M99.');
        value=complex(nan(size(flat)),nan(size(flat)));
        value(valid_denominator)=rough(valid_denominator)./flat(valid_denominator);
        ratios{bi}=value;
    end
    c=data.cases(ci); c.G_BH_native=ratios{2}; c.G_BH=conj(ratios{2});
    c.G_BH_5001=conj(ratios{1});
    valid=isfinite(c.G_BH)&isfinite(c.G_PE)&isfinite(c.G_BIE);
    assert(all(valid(data.mask)),'Invalid curves inside fixed M99.');
    c.valid_mask=valid; c.BH5001_valid_mask=isfinite(c.G_BH_5001);
    c.plot_mask=data.mask & valid;
    assert(all(abs(c.G_BH(c.plot_mask))>0 & abs(c.G_BIE(c.plot_mask))>0 & abs(c.G_PE(c.plot_mask))>0));
    a={c.G_PE,c.G_BH,c.G_PE}; b={c.G_BIE,c.G_BIE,c.G_BH}; tags={'PE_BIE','BH_BIE','PE_BH'};
    c.residuals=repmat(struct('pair','','complex_abs',[],'phase_rad',[], ...
        'magnitude_ratio_db',[]),1,3);
    c.transforms=struct;
    for method={'PE','BH','BIE'}
        label=method{1}; value=c.(['G_' label]);
        phase=angle(value); gaps=[false;abs(diff(phase))>pi];
        plotted_phase=phase; plotted_phase(gaps)=NaN;
        c.transforms.(label)=struct('magnitude',abs(value),'phase_rad',phase, ...
            'phase_mrad',phase*1000,'phase_wrap_gap',gaps,'plotted_phase_mrad',plotted_phase*1000);
    end
    for p=1:3
        c.residuals(p)=struct('pair',tags{p},'complex_abs',abs(a{p}-b{p}), ...
            'phase_rad',angle(a{p}.*conj(b{p})), ...
            'magnitude_ratio_db',20*log10(abs(a{p})./abs(b{p})));
        row=metric(a{p},b{p},data.weights,c.plot_mask); row.case_name=c.name; row.pair=tags{p};
        allmetrics=[allmetrics;row]; %#ok<AGROW>
    end
    row=metric(c.G_BH_5001,c.G_BH,data.weights,c.plot_mask); row.case_name=c.name; row.pair='BH_5001_vs_10001';
    sensitivity=[sensitivity;row]; %#ok<AGROW>
    data.cases(ci).G_BH=c.G_BH; data.cases(ci).G_BH_native=c.G_BH_native;
    data.cases(ci).G_BH_5001=c.G_BH_5001; data.cases(ci).plot_mask=c.plot_mask;
    data.cases(ci).residuals=c.residuals;
    data.cases(ci).transforms=c.transforms;
    data.cases(ci).valid_mask=c.valid_mask; data.cases(ci).BH5001_valid_mask=c.BH5001_valid_mask;
    T=table(data.x_m,c.eta_m,data.weights,c.plot_mask,real(c.G_PE),imag(c.G_PE), ...
        real(c.G_BH),imag(c.G_BH),real(c.G_BIE),imag(c.G_BIE), ...
        real(c.G_BH_native),imag(c.G_BH_native),real(c.G_BH_5001),imag(c.G_BH_5001), ...
        'VariableNames',{'x_m','eta_m','incident_weight','M99_mask','PE_real','PE_imag', ...
        'BH_real','BH_imag','BIE_real','BIE_imag','BH_native_real','BH_native_imag','BH5001_real','BH5001_imag'});
    for p=1:3
        tag=tags{p}; T.([tag '_complex_abs'])=c.residuals(p).complex_abs;
        T.([tag '_phase_rad'])=c.residuals(p).phase_rad;
        T.([tag '_magnitude_ratio_db'])=c.residuals(p).magnitude_ratio_db;
    end
    for method={'PE','BH','BIE'}
        tag=method{1}; value=c.transforms.(tag);
        T.([tag '_magnitude'])=value.magnitude;
        T.([tag '_principal_phase_rad'])=value.phase_rad;
        T.([tag '_phase_wrap_gap'])=value.phase_wrap_gap;
        T.([tag '_plotted_phase_mrad'])=value.plotted_phase_mrad;
    end
    T.ratio_valid=c.valid_mask; T.BH5001_ratio_valid=c.BH5001_valid_mask;
    writetable(T,fullfile(out,'plot_data',[c.name '.csv']));
end
data.metrics=struct2table(allmetrics); data.beam_sensitivity=struct2table(sensitivity);
data.geometry=struct2table(geometry);
historical=table(repelem(string({s.cases.name}).',3),repmat(["PE_BIE";"BH_BIE";"PE_BH"],3,1), ...
    [.00315757;8.09198e-6;.00315855;.01608996;3.88955e-5;.01609469;.08090811;9.97602e-5;.08092279], ...
    repelem([6.00191e-8;6.02438e-8;6.17857e-8],3), ...
    'VariableNames',{'case_name','pair','historical_complex_l2','historical_BIE_U'});
data.historical=historical;
data.definitions=struct('complex_l2','sqrt(sum(w*abs(A-B)^2)/sum(w*abs(B)^2))', ...
    'phase_rad','angle(A*conj(B))','magnitude_ratio_db','20*log10(abs(A)/abs(B))', ...
    'weight_rule','Stage-0 incident energy, normalized over fixed M99 mask', ...
    'invalid_rule','Zero Bellhop flat denominator gives explicit complex NaN outside fixed M99; valid masks saved; no interpolation; all M99 samples finite', ...
    'phase_plot_rule','Principal phase in mrad; omit sample after a consecutive principal-phase jump exceeding pi');
writetable(data.metrics,fullfile(out,'metrics','current_metrics.csv'));
writetable(data.beam_sensitivity,fullfile(out,'metrics','beam_sensitivity.csv'));
writetable(data.geometry,fullfile(out,'metrics','geometry_checks.csv'));
writetable(historical,fullfile(out,'metrics','historical_report_values.csv'));
save(fullfile(out,'plot_data','comparison_data.mat'),'data','-v7.3');
s.status='PLOT_DATA_COMPLETE';
end

function m=metric(a,b,weights,mask)
w=weights(mask); w=w/sum(w); a=a(mask); b=b(mask);
m=struct('complex_l2',sqrt(sum(w.*abs(a-b).^2)/sum(w.*abs(b).^2)), ...
    'phase_rms_rad',sqrt(sum(w.*angle(a.*conj(b)).^2)), ...
    'magnitude_ratio_rms_db',sqrt(sum(w.*(20*log10(abs(a)./abs(b))).^2)), ...
    'sample_count',nnz(mask));
end

function stage_inputs(s,out,exe)
% Stage independent executable inputs; never launch MATLAB or a solver.
stamped=fullfile(out,'source_data','bellhop_internal_wall.exe');
assert(isfile(stamped) && strcmp(file_hash(stamped),file_hash(exe)));
cfg=s.config; xx=linspace(-80,80,4097)'; zsort=sort(-s.x_m);
for ix=find(~[s.jobs.complete])
    job=s.jobs(ix); ci=job.case_index; amp=0; if ci>0, amp=s.cases(ci).amplitude_m; end
    folder=fullfile(out,'source_data','bellhop',sprintf('case_%d_B%d',ci,job.beams));
    if ~isfolder(folder), mkdir(folder); end
    c=struct('bellhop_exe',stamped,'case_root',fullfile(folder,'field'), ...
        'run_type',cfg.run_type,'source_geometry','X','frequency_hz',4000,'c0_mps',1500, ...
        'source_depth_m',0,'receiver_depths_m',zsort.', ...
        'receiver_ranges_m',[cfg.guard_receiver_range_m 103],'mapped_receiver_range_m',103, ...
        'beam_count',job.beams,'angle_limits_deg',cfg.angle_limits_deg,'step_m',.05, ...
        'domain_half_depth_m',cfg.domain_half_depth_m, ...
        'source_pattern_angles_deg',s.source_pattern.angles_deg, ...
        'source_pattern_level_db',s.source_pattern.level_db, ...
        'wall_r0_m',100,'wall_seed',cfg.wall_seed,'wall_profile_r_m',100-surface(xx,amp), ...
        'wall_profile_z_m',xx);
    assert(~isfile([c.case_root '.shd']),'Do not restage a started job.');
    write_bellhop_unfolded_gaussian_env_vertical(c.case_root,c);
    fid=fopen([c.case_root '.iwpm'],'w'); clean=onCleanup(@()fclose(fid));
    fprintf(fid,'%.17g\n%.17g\n%d\n%d\n',100,103,cfg.wall_seed,numel(xx));
    fprintf(fid,'%.17g %.17g\n',[c.wall_profile_r_m(:) xx].'); clear clean
end
disp('Remaining Bellhop inputs staged; no solver launched.');
end

function eta=surface(x,A)
t=(abs(x)-42)/8; w=double(abs(x)<=42); ix=abs(x)>42 & abs(x)<50;
w(ix)=1-10*t(ix).^3+15*t(ix).^4-6*t(ix).^5; eta=A*sin(.1*x).*w;
end

function checkpoint(state,out)
save(fullfile(out,'state.mat'),'state','-v7.3');
fid=fopen(fullfile(out,'status.json'),'w'); clean=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(struct('status',state.status,'completed_jobs',nnz([state.jobs.complete]), ...
    'total_jobs',8,'new_PE_solves',0,'new_BIE_solves',0),PrettyPrint=true));
end

function value=file_hash(path)
fid=fopen(path,'rb'); assert(fid>=0); clean=onCleanup(@()fclose(fid));
md=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid), bytes=fread(fid,1048576,'*uint8'); md.update(bytes); end
value=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
