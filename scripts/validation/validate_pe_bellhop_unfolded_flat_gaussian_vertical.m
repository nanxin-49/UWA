function validation = validate_pe_bellhop_unfolded_flat_gaussian_vertical(overrides)
%VALIDATE_PE_BELLHOP_UNFOLDED_FLAT_GAUSSIAN_VERTICAL
% Validation-only unfolded-coordinate PE--Bellhop comparison.

if nargin<1 || isempty(overrides), overrides=struct(); end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
run(fullfile(root,'scripts','bootstrap_project.m'));
cfg=local_config(root,overrides);
if ~exist(cfg.output_dir,'dir'), mkdir(cfg.output_dir); end

if cfg.use_saved_normalization_audit && exist(cfg.normalization_audit_fallback_file,'file')==2
    loaded=load(cfg.normalization_audit_fallback_file,'validation');
    if ~isfield(loaded,'validation') || ~isfield(loaded.validation,'passed') || ~loaded.validation.passed
        error('Saved Bellhop normalization audit is missing or failed.');
    end
    bh_norm=loaded.validation;
    bh_norm.audit_execution='validated_saved_fallback';
elseif cfg.run_normalization_audit
    try
        bh_norm=validate_bellhop_freefield_normalization_vertical(struct( ...
            'output_dir',fullfile(cfg.output_dir,'normalization_audit'), ...
            'bellhop_exe',cfg.bellhop_exe,'frequencies_hz',cfg.audit_frequencies_hz, ...
            'beam_counts',cfg.normalization_beam_counts,'step_values_m',0.05,'ranges_m',[20 40 70 100]));
        bh_norm.audit_execution='fresh';
    catch audit_error
        if isempty(cfg.normalization_audit_fallback_file) || exist(cfg.normalization_audit_fallback_file,'file')~=2
            rethrow(audit_error);
        end
        loaded=load(cfg.normalization_audit_fallback_file,'validation');
        if ~isfield(loaded,'validation') || ~isfield(loaded.validation,'passed') || ~loaded.validation.passed
            rethrow(audit_error);
        end
        bh_norm=loaded.validation;
        bh_norm.audit_execution='validated_saved_fallback';
        bh_norm.audit_fallback_reason=audit_error.message;
    end
else
    bh_norm=struct('passed',true,'selected_spatial_sign',1,'bellhop_source_constant',1, ...
        'normalization_constant',1,'green_conversion_factor',1/(4*pi));
end
if ~bh_norm.passed, error('Bellhop normalization audit failed.'); end

freq=cfg.frequencies_hz(:).';
pe=local_run_pe(cfg,freq);
bh=local_run_bellhop(cfg,freq,bh_norm);
as_d=zeros(size(freq)); as_r=zeros(size(freq));
for ii=1:numel(freq)
    as_d(ii)=local_gaussian_weyl_axis(freq(ii),cfg.c0_mps,cfg.sigma_src_m,cfg.direct_range_m);
    as_r(ii)=-local_gaussian_weyl_axis(freq(ii),cfg.c0_mps,cfg.sigma_src_m,cfg.reflect_range_m);
end
relative=local_relative(freq,pe.H_direct,pe.H_reflect,bh.H_direct,bh.H_reflect,as_d,as_r);
spatial=local_spatial(cfg,bh_norm);
beam=local_beam(cfg,bh_norm);
identity=local_identity(cfg,pe);
native=local_native(cfg);
checks=local_checks(relative,spatial,beam,identity,native,pe.closure,cfg);
validation=struct('schema_version','1.0.0','config',cfg,'bellhop_normalization',bh_norm, ...
    'pe',pe,'bellhop',bh,'continuous_gaussian_reference',struct('H_direct',as_d, ...
    'H_reflect',as_r,'source_scale_to_bellhop',as_d./bh.H_direct), ...
    'relative_metrics',relative,'spatial_profiles',spatial,'beam_convergence',beam, ...
    'unfolded_identity',identity,'native_auxiliary',native,'checks',checks, ...
    'passed',all(checks.passed));
validation.files=local_outputs(validation,root);
if cfg.fail_on_check && ~validation.passed
    error('Unfolded PE--Bellhop validation failed; see %s.',validation.files.report);
end
end

function cfg=local_config(root,o)
exe=getenv('BELLHOP_EXE');
if isempty(exe)
    candidates={'E:\MISC\BELLHOP\AcousticsToolbox_2017\Bellhop\bellhop.exe', ...
        'E:\MISC\fxx\Bellhop相关\Bellhop例程包\atWin10_2020_11_4\atWin10_2020_11_4\windows-bin-20201102\bellhop.exe'};
    for ii=1:numel(candidates), if exist(candidates{ii},'file')==2, exe=candidates{ii}; break; end, end
end
cfg=struct('output_dir',fullfile(root,'results','validation','pe_bellhop_unfolded_flat_gaussian'), ...
    'report_path',fullfile(root,'reports','pe_bellhop_unfolded_flat_gaussian_report.md'), ...
    'bellhop_exe',exe,'c0_mps',1500,'z_tx_m',100,'z_rx_m',3,'sigma_src_m',0.3, ...
    'frequencies_hz',linspace(4000,8000,65),'audit_frequencies_hz',[4000 6000 8000], ...
    'direct_range_m',97,'reflect_range_m',103,'pe_width_m',192.1875,'pe_n',984, ...
    'pe_stepz_lamb',0.5,'sponge_ratio',0,'alpha_max_np_per_m',0, ...
    'receiver_offsets_m',[0 1.953125 4.8828125 9.765625 14.6484375 19.53125], ...
    'profile_frequencies_hz',[4000 6000 8000],'bellhop_beam_count',10001, ...
    'beam_audit_count',5001,'bellhop_step_m',0.05,'bellhop_domain_half_depth_m',1000, ...
    'normalization_beam_counts',[501 2001], ...
    'bellhop_angle_limits_deg',[-30 30],'source_pattern_clip_db',-120, ...
    'run_normalization_audit',true,'use_saved_normalization_audit',true, ...
    'run_native_auxiliary',true,'fail_on_check',true, ...
    'tl_limit_db',0.25,'phase_limit_rad',0.05,'spatial_complex_l2_limit',0.02, ...
    'group_delay_limit_s',20e-6,'pdp_delay_limit_s',0.125e-3,'beam_tl_limit_db',0.1, ...
    'beam_phase_limit_rad',0.02,'unfolded_l2_limit',1e-10,'native_delay_limit_s',0.1e-3, ...
    'normalization_audit_fallback_file',fullfile(root,'results','validation','pe_bellhop_freefield', ...
    'bellhop_normalization','bellhop_freefield_normalization.mat'));
names=fieldnames(o);
for ii=1:numel(names)
    if ~isfield(cfg,names{ii})
        error('Unknown override: %s',names{ii});
    end
    cfg.(names{ii})=o.(names{ii});
end
if isfield(o,'output_dir') && ~isfield(o,'report_path'), cfg.report_path=fullfile(cfg.output_dir,'unfolded_validation_report.md'); end
if any(diff(cfg.frequencies_hz)<=0), error('frequencies_hz must be strictly increasing.'); end
if any(~ismember(cfg.profile_frequencies_hz,cfg.frequencies_hz)), error('Profile frequencies must be in frequencies_hz.'); end
end

function out=local_run_pe(cfg,freq)
flat=zeros(cfg.pe_n,cfg.pe_n);
p=local_pe_params(cfg,cfg.z_tx_m,cfg.z_rx_m,true,freq,flat,'rx_only'); chain=vertical_channel_model(p);
% The reflection-enabled chain already exposes its deterministic direct
% branch. Avoid a second complete PE march solely to reproduce H_direct.
out=struct('chain',chain,'direct_only',[],'H_direct',chain.H_direct_physical_f(:).', ...
    'H_reflect',chain.H_reflect_physical_f(:).','H_total',chain.H_physical_f(:).', ...
    'closure',max(abs(chain.H_physical_f-chain.H_direct_physical_f-chain.H_reflect_physical_f)));
end

function p=local_pe_params(cfg,ztx,zrx,reflection,freq,eta,save_mode)
p=struct('f0',freq,'c0',cfg.c0_mps,'z_max',ztx,'z_tx',ztx,'z_rx',zrx, ...
    'xw',cfg.pe_width_m,'yw',cfg.pe_width_m,'nx',cfg.pe_n,'ny',cfg.pe_n, ...
    'x_tx',0,'y_tx',0,'x_rx',0,'y_rx',0,'sigma_src_m',cfg.sigma_src_m, ...
    'source_mode','gaussian','stepz_lamb',cfg.pe_stepz_lamb, ...
    'sponge_ratio',cfg.sponge_ratio,'alpha_max_np_per_m',cfg.alpha_max_np_per_m, ...
    'validation_allow_extended_window',true,'validation_allow_extended_sponge_ratio',true, ...
    'env_mode','uniform','show_figures',false,'save_mode',save_mode,'enforce_1_over_R',false, ...
    'enable_surface_reflection',reflection,'surface_boundary_model','kirchhoff_spatial', ...
    'surface_reflect_coeff',-1,'surface_phase_mode','normal','surface_elevation_override_xy',eta, ...
    'surface_ssa_random_scatter',false,'surface_kstat_random_scatter',false, ...
    'surface_wavefield_diagnostics',false,'enable_bubbles',false,'channel_phase_reference','direct_dsp', ...
    'use_gpu',false);
end

function bh=local_run_bellhop(cfg,freq,norm)
bh=struct('frequency_hz',freq,'H_direct',complex(zeros(size(freq))), ...
    'H_reflect',complex(zeros(size(freq))),'runs',{cell(size(freq))});
for ii=1:numel(freq)
    pat=local_pattern(freq(ii),cfg);
    c=local_bh_cfg(cfg,freq(ii),cfg.receiver_offsets_m,[cfg.direct_range_m cfg.reflect_range_m], ...
        'C',cfg.bellhop_beam_count,pat,fullfile(cfg.output_dir,'bellhop','formal',local_frequency_tag(freq(ii))));
    r=run_bellhop_unfolded_gaussian_vertical(c); bh.runs{ii}=r;
    p=local_convert(r.data.pressure,norm); [~,i0]=min(abs(r.data.receiver_depth_m));
    [~,id]=min(abs(r.data.receiver_range_m-cfg.direct_range_m)); [~,ir]=min(abs(r.data.receiver_range_m-cfg.reflect_range_m));
    bh.H_direct(ii)=p(i0,id); bh.H_reflect(ii)=-p(i0,ir);
end
end

function c=local_bh_cfg(cfg,f,depths,ranges,run_type,beams,pat,case_root)
c=struct('bellhop_exe',cfg.bellhop_exe,'case_root',case_root,'frequency_hz',f, ...
    'c0_mps',cfg.c0_mps,'source_depth_m',0,'receiver_depths_m',depths, ...
    'receiver_ranges_m',ranges,'run_type',run_type,'beam_count',beams, ...
    'angle_limits_deg',cfg.bellhop_angle_limits_deg,'step_m',cfg.bellhop_step_m, ...
    'domain_half_depth_m',cfg.bellhop_domain_half_depth_m, ...
    'source_pattern_angles_deg',pat.angles_deg,'source_pattern_level_db',pat.level_db);
end

function p=local_convert(raw,norm)
if norm.selected_spatial_sign<0, p=conj(raw)./conj(norm.bellhop_source_constant); else, p=raw./norm.bellhop_source_constant; end
end

function pat=local_pattern(f,cfg)
angles=linspace(cfg.bellhop_angle_limits_deg(1),cfg.bellhop_angle_limits_deg(2),2401).';
k=2*pi*f/cfg.c0_mps; th=deg2rad(angles); D=cos(th).*exp(-0.5*(k*cfg.sigma_src_m*sin(th)).^2);
D(abs(angles)>89.9)=0; D=abs(D)/max(abs(D));
pat=struct('angles_deg',angles,'level_db',20*log10(max(D,10^(cfg.source_pattern_clip_db/20))), ...
    'formula','cos(theta)*exp(-0.5*(k*sigma*sin(theta))^2)','frequency_hz',f);
end

function m=local_relative(f,pe_d,pe_r,bh_d,bh_r,as_d,as_r)
qp=pe_r./pe_d; qb=bh_r./bh_d; qa=as_r./as_d;
m=struct('frequency_hz',f,'q_pe',qp,'q_bellhop',qb,'q_as',qa, ...
    'tl_difference_db',20*log10(max(abs(qp),realmin)./max(abs(qb),realmin)), ...
    'phase_difference_rad',angle(qp.*conj(qb)), ...
    'as_tl_difference_db',20*log10(max(abs(qp),realmin)./max(abs(qa),realmin)), ...
    'as_phase_difference_rad',angle(qp.*conj(qa)), ...
    'relative_total_pe',1+qp,'relative_total_bellhop',1+qb);
m.total_tl_difference_db=20*log10(max(abs(1+qp),realmin)./max(abs(1+qb),realmin));
m.total_phase_difference_rad=angle((1+qp).*conj(1+qb));
if numel(f)>=3
    m.group_delay_pe=local_group_delay(f,qp); m.group_delay_bellhop=local_group_delay(f,qb);
    m.group_delay_difference_s=m.group_delay_pe-m.group_delay_bellhop;
    m.internal_group_delay_max_abs_s=max(abs(m.group_delay_difference_s(2:end-1)));
    m.pdp_pe=local_pdp(qp,f); m.pdp_bellhop=local_pdp(qb,f);
    m.pdp_delay_difference_s=m.pdp_pe.delay_s-m.pdp_bellhop.delay_s;
else
    m.group_delay_pe=NaN(size(f)); m.group_delay_bellhop=NaN(size(f)); m.group_delay_difference_s=NaN(size(f));
    m.internal_group_delay_max_abs_s=NaN; m.pdp_pe=struct('delay_s',NaN); m.pdp_bellhop=struct('delay_s',NaN); m.pdp_delay_difference_s=NaN;
end
end

function gd=local_group_delay(f,h)
gd=gradient(unwrap(angle(h)),f)/(2*pi);
end

function p=local_pdp(h,f)
n=numel(h); w=0.5-0.5*cos(2*pi*(0:n-1)/(n-1)); nfft=8*n;
z=ifft(h(:).'.*w,nfft); [~,ix]=max(abs(z).^2); dt=1/(nfft*(f(2)-f(1)));
delay=(ix-1)*dt; if delay>0.5*nfft*dt, delay=delay-nfft*dt; end
p=struct('delay_s',delay,'time_s',((0:nfft-1)-double(delay<0)*nfft)*dt,'power',abs(z).^2);
end

function s=local_spatial(cfg,norm)
freq=cfg.profile_frequencies_hz(:).'; off=cfg.receiver_offsets_m(:).';
tmpl=struct('frequency_hz',NaN,'offset_m',NaN,'pe_direct',complex(NaN),'pe_reflect',complex(NaN), ...
    'bh_direct',complex(NaN),'bh_reflect',complex(NaN),'direct_tl_error_db',NaN, ...
    'reflect_tl_error_db',NaN,'direct_phase_error_rad',NaN,'reflect_phase_error_rad',NaN, ...
    'direct_complex_error',NaN,'reflect_complex_error',NaN);
rows=repmat(tmpl,1,numel(off)*numel(freq));
nrow=0;
for ff=1:numel(freq)
    od=vertical_channel_model(local_pe_params(cfg,cfg.z_tx_m,cfg.z_rx_m,false,freq(ff),[],'slice'));
    pr=local_pe_params(cfg,cfg.z_tx_m,cfg.z_rx_m,true,freq(ff),zeros(cfg.pe_n,cfg.pe_n),'rx_only');
    pr.surface_wavefield_diagnostics=true; pr.surface_wavefield_max_z_samples=3;
    or=vertical_channel_model(pr); meta=or.surface_wavefield_meta;
    if ~isfield(meta,'receiver_reflected_xy'), error('Missing receiver_reflected_xy.'); end
    [~,iy]=min(abs(od.y)); [~,ix0]=min(abs(od.x)); pd=od.psifinal_xy(iy,:); prf=meta.receiver_reflected_xy(iy,:);
    pat=local_pattern(freq(ff),cfg); c=local_bh_cfg(cfg,freq(ff),off,[cfg.direct_range_m cfg.reflect_range_m], ...
        'C',cfg.bellhop_beam_count,pat,fullfile(cfg.output_dir,'bellhop','profiles',local_frequency_tag(freq(ff))));
    rb=run_bellhop_unfolded_gaussian_vertical(c); pb=local_convert(rb.data.pressure,norm); [~,ib0]=min(abs(rb.data.receiver_depth_m));
    for jj=1:numel(off)
        [~,ix]=min(abs(od.x-off(jj))); [~,ib]=min(abs(rb.data.receiver_depth_m-off(jj)));
        a=pd(ix)/pd(ix0); b=pb(ib,1)/pb(ib0,1); ar=prf(ix)/prf(ix0); br=pb(ib,2)/pb(ib0,2);
        q=tmpl; q.frequency_hz=freq(ff); q.offset_m=off(jj); q.pe_direct=pd(ix); q.pe_reflect=prf(ix); q.bh_direct=pb(ib,1); q.bh_reflect=-pb(ib,2);
        q.direct_tl_error_db=20*log10(max(abs(a),realmin)/max(abs(b),realmin)); q.reflect_tl_error_db=20*log10(max(abs(ar),realmin)/max(abs(br),realmin));
        q.direct_phase_error_rad=angle(a*conj(b)); q.reflect_phase_error_rad=angle(ar*conj(br)); q.direct_complex_error=abs(a-b)/max(abs(b),realmin); q.reflect_complex_error=abs(ar-br)/max(abs(br),realmin);
        nrow=nrow+1;
        rows(nrow)=q;
    end
end
rows=rows(1:nrow);
s=struct('rows',rows,'table',struct2table(rows),'max_direct_tl_error_db',max(abs([rows.direct_tl_error_db])), ...
    'max_reflect_tl_error_db',max(abs([rows.reflect_tl_error_db])),'max_direct_phase_error_rad',max(abs([rows.direct_phase_error_rad])), ...
    'max_reflect_phase_error_rad',max(abs([rows.reflect_phase_error_rad])),'max_complex_l2',max([[rows.direct_complex_error] [rows.reflect_complex_error]]));
end

function b=local_beam(cfg,norm)
freq=cfg.audit_frequencies_hz(:).'; rows=repmat(struct('frequency_hz',NaN,'tl_difference_db',NaN,'phase_difference_rad',NaN),numel(freq),1);
for ii=1:numel(freq)
    pat=local_pattern(freq(ii),cfg); depths=[0 cfg.receiver_offsets_m]; ranges=[cfg.direct_range_m cfg.reflect_range_m];
    c1=local_bh_cfg(cfg,freq(ii),depths,ranges,'C',cfg.beam_audit_count,pat,fullfile(cfg.output_dir,'bellhop','beam_audit',sprintf('%s_b%d',local_frequency_tag(freq(ii)),cfg.beam_audit_count)));
    c2=c1; c2.beam_count=cfg.bellhop_beam_count; c2.case_root=fullfile(cfg.output_dir,'bellhop','beam_audit',sprintf('%s_b%d',local_frequency_tag(freq(ii)),cfg.bellhop_beam_count));
    a=run_bellhop_unfolded_gaussian_vertical(c1); z=run_bellhop_unfolded_gaussian_vertical(c2); pa=local_convert(a.data.pressure,norm); pz=local_convert(z.data.pressure,norm);
    rows(ii).frequency_hz=freq(ii); rows(ii).tl_difference_db=max(abs(20*log10(max(abs(pa),realmin)./max(abs(pz),realmin)))); rows(ii).phase_difference_rad=sqrt(mean(angle(pa(:).*conj(pz(:))).^2));
end
b=struct('rows',rows,'table',struct2table(rows),'max_tl_difference_db',max(abs([rows.tl_difference_db])),'max_phase_difference_rad',max([rows.phase_difference_rad]));
end

function tag=local_frequency_tag(f)
% Bellhop case roots must not contain a decimal point: Bellhop treats the
% suffix after the first dot as a file extension and truncates the basename.
tag=strrep(sprintf('f%.10g',f),'.','p');
tag=strrep(tag,'-','m');
end

function u=local_identity(cfg,pe)
freq=cfg.frequencies_hz(:).'; rows=repmat(struct('frequency_hz',NaN,'relative_l2',NaN),numel(freq),1);
[X,Y]=meshgrid(pe.chain.x,pe.chain.y);
[psi0,source_meta]=gaussian_source_initial_field_vertical(X,Y,pe.chain.config);
[~,ix0]=min(abs(pe.chain.x)); [~,iy0]=min(abs(pe.chain.y));
for ii=1:numel(freq)
    [psi_as,as_meta]=exact_angular_spectrum_one_step_vertical(psi0,cfg.pe_width_m, ...
        cfg.pe_width_m,freq(ii),cfg.c0_mps,cfg.reflect_range_m);
    a=pe.chain.H_reflect_reduced_f(ii); b=-psi_as(iy0,ix0);
    rows(ii).frequency_hz=freq(ii); rows(ii).relative_l2=abs(a-b)/max(abs(b),realmin);
end
u=struct('rows',rows,'table',struct2table(rows),'max_relative_l2',max([rows.relative_l2]), ...
    'source_meta',source_meta,'as_meta',as_meta,'reference','-one-step exact discrete angular spectrum at 103 m');
end

function n=local_native(cfg)
n=struct('executed',false,'passed',true,'paths',table(),'geometry',struct(),'files',struct()); if ~cfg.run_native_auxiliary, return; end
c=struct('bellhop_exe',cfg.bellhop_exe,'case_root',fullfile(cfg.output_dir,'native_auxiliary','tx80_rx10_r6'),'frequency_hz',4000, ...
    'water_depth_m',100,'c0_mps',cfg.c0_mps,'z_tx_m',80,'z_rx_m',10,'offset_m',6,'beam_count',10001, ...
    'angle_min_deg',-89.5,'angle_max_deg',-60,'arrival_cluster_tolerance_ms',0.2);
r=run_bellhop_flat_surface_arrivals_vertical(c); n=struct('executed',true,'passed',height(r.paths)==2 && all(r.paths.bottom_bounce_count==0), ...
    'paths',r.paths,'geometry',r.geometry,'files',r.files);
end

function checks=local_checks(r,s,b,u,n,closure,cfg)
names=["bellhop_beam_tl";"bellhop_beam_phase";"spatial_tl";"spatial_phase";"spatial_complex_l2"; ...
    "wideband_tl";"wideband_phase";"group_delay";"pdp_delay";"unfolded_pe_identity";"native_auxiliary";"H_closure"];
values=[b.max_tl_difference_db;b.max_phase_difference_rad;max(s.max_direct_tl_error_db,s.max_reflect_tl_error_db); ...
    max(s.max_direct_phase_error_rad,s.max_reflect_phase_error_rad);s.max_complex_l2;max(abs(r.tl_difference_db)); ...
    max(abs(r.phase_difference_rad));r.internal_group_delay_max_abs_s;abs(r.pdp_delay_difference_s);u.max_relative_l2;double(~n.passed);closure];
limits=[cfg.beam_tl_limit_db;cfg.beam_phase_limit_rad;cfg.tl_limit_db;cfg.phase_limit_rad;cfg.spatial_complex_l2_limit; ...
    cfg.tl_limit_db;cfg.phase_limit_rad;cfg.group_delay_limit_s;cfg.pdp_delay_limit_s;cfg.unfolded_l2_limit;0;1e-12];
checks=table(names,values,limits,values<=limits,'VariableNames',{'check_name','value','limit','passed'});
end

function files=local_outputs(v,root)
out=v.config.output_dir; if ~exist(out,'dir'), mkdir(out); end
rt=table(v.relative_metrics.frequency_hz(:),v.relative_metrics.tl_difference_db(:), ...
    v.relative_metrics.phase_difference_rad(:),v.relative_metrics.group_delay_difference_s(:), ...
    'VariableNames',{'frequency_hz','tl_difference_db','phase_difference_rad','group_delay_difference_s'});
writetable(rt,fullfile(out,'relative_wideband.csv')); writetable(v.spatial_profiles.table,fullfile(out,'spatial_profiles.csv'));
writetable(v.beam_convergence.table,fullfile(out,'beam_convergence.csv')); writetable(v.unfolded_identity.table,fullfile(out,'unfolded_identity.csv')); writetable(v.checks,fullfile(out,'checks.csv'));
mat_file=fullfile(out,'pe_bellhop_unfolded_flat_gaussian_validation.mat'); validation=v; save(mat_file,'validation','-v7.3');
fig_file=fullfile(out,'pe_bellhop_unfolded_flat_gaussian.png'); local_plot(v,fig_file); local_report(v,root);
files=struct('mat',mat_file,'figure',fig_file,'report',v.config.report_path,'relative_csv',fullfile(out,'relative_wideband.csv'),'checks_csv',fullfile(out,'checks.csv'));
end

function local_plot(v,file)
fig=figure('Visible','off','Color','w','Position',[100 100 1200 800]); cleanup=onCleanup(@()close(fig));
subplot(2,2,1); plot(v.relative_metrics.frequency_hz/1000,20*log10(abs(v.relative_metrics.q_pe)),'o-',v.relative_metrics.frequency_hz/1000,20*log10(abs(v.relative_metrics.q_bellhop)),'s--'); grid on; xlabel('f (kHz)'); ylabel('20log10|H_r/H_d| (dB)'); legend('PE','Bellhop');
subplot(2,2,2); plot(v.relative_metrics.frequency_hz/1000,unwrap(angle(v.relative_metrics.q_pe)),'o-',v.relative_metrics.frequency_hz/1000,unwrap(angle(v.relative_metrics.q_bellhop)),'s--'); grid on; xlabel('f (kHz)'); ylabel('phase (rad)'); legend('PE','Bellhop');
subplot(2,2,3); plot(v.relative_metrics.frequency_hz/1000,1e6*v.relative_metrics.group_delay_pe,'o-',v.relative_metrics.frequency_hz/1000,1e6*v.relative_metrics.group_delay_bellhop,'s--'); grid on; xlabel('f (kHz)'); ylabel('relative group delay (us)'); legend('PE','Bellhop');
subplot(2,2,4); bar(categorical(v.checks.check_name),v.checks.value); grid on; ylabel('check value'); title(sprintf('passed=%d',v.passed));
exportgraphics(fig,file,'Resolution',180); clear cleanup
end

function local_report(v,~)
fid=fopen(v.config.report_path,'w','n','UTF-8'); if fid<0, error('Cannot create report.'); end
cleanup=onCleanup(@()fclose(fid)); c=v.config;
fprintf(fid,'# 展开坐标 PE--Bellhop 平面海面 Gaussian 验证\n\n');
fprintf(fid,'状态：`passed=%s`。主方法将原竖直传播展开为 Bellhop 水平自由场，并在展开接收端施加压力释放系数 `-1`。\n\n',string(v.passed));
fprintf(fid,'## 环境\n\n- `c=%.0f m/s`，Tx 深度 `%.3g m`，Rx 深度 `%.3g m`，Gaussian `sigma=%.3g m`。\n',c.c0_mps,c.z_tx_m,c.z_rx_m,c.sigma_src_m);
fprintf(fid,'- 直达展开距离 `%.3g m`，反射展开距离 `%.3g m`；PE 窗口 `%.6g m`、`N=%d`、sponge off。\n',c.direct_range_m,c.reflect_range_m,c.pe_width_m,c.pe_n);
fprintf(fid,'- 频率轴 `%d` 点：`%.3g--%.3g kHz`；Bellhop 使用 `%d` 条波束，5001 条用于收敛审计。\n\n',numel(c.frequencies_hz),c.frequencies_hz(1)/1000,c.frequencies_hz(end)/1000,c.bellhop_beam_count);
fprintf(fid,'## 结果摘要\n\n| check | value | limit | pass |\n|---|---:|---:|:---:|\n');
for ii=1:height(v.checks), fprintf(fid,'| %s | %.8g | %.8g | %d |\n',v.checks.check_name(ii),v.checks.value(ii),v.checks.limit(ii),v.checks.passed(ii)); end
fprintf(fid,'\n');
if v.passed, fprintf(fid,'PE marching、平面反射展开和 Bellhop 在当前 Gaussian、均匀介质、4--8 kHz 条件下的相对复传播结果通过预设门槛。\n\n'); else, fprintf(fid,'至少一个门槛未通过；失败项只作为诊断，不修改阈值，也不反向改变 PE 主线。\n\n'); end
fprintf(fid,'本结论不覆盖粗糙海面、深度相关声速、海底、气泡、Doppler 或实际换能器。更换换能器后必须用测量/建模的幅相指向性重新生成 `.sbp`，并重新确认 PE 横向窗口和边缘能量。\n\n');
fprintf(fid,'![summary](../results/validation/pe_bellhop_unfolded_flat_gaussian/pe_bellhop_unfolded_flat_gaussian.png)\n'); clear cleanup
end

function value=local_gaussian_weyl_axis(f,c0,sigma,L)
k=2*pi*f/c0; qmax=k+max(40/sigma,20*k);
fn=@(q) q.*exp(-0.5*sigma^2*q.^2).*exp(1i*L*sqrt(complex(k^2-q.^2,0)));
value=sigma^2*integral(fn,0,qmax,'ArrayValued',true,'RelTol',2e-7,'AbsTol',1e-10);
end
