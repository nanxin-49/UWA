function validation = validate_reflected_chain_wideband_vertical(overrides)
%VALIDATE_REFLECTED_CHAIN_WIDEBAND_VERTICAL Wideband and post-gate sponge audit.
if nargin<1, overrides=struct(); end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
run(fullfile(root,'scripts','bootstrap_project.m'));
cfg=local_config(root,overrides);
gate_path=fullfile(cfg.output_dir,'reflected_chain_window_validation.mat');
if exist(gate_path,'file')~=2, error('Single-frequency reflected-chain gate is missing.'); end
d=load(gate_path,'validation'); gate=d.validation;
if ~gate.converged, error('Wideband validation is blocked until the no-sponge chain converges.'); end
master=load(fullfile(cfg.output_dir,'fixed_master_surface.mat'),'master_eta','master_meta');

required_nominal=gate.cases(gate.required_index).width_nominal_m;
reference_nominal=gate.cases(gate.reference_index).width_nominal_m;
specs=struct('name',{'reference_no_sponge','nominal160_no_sponge','production_default'}, ...
    'nominal_width_m',{reference_nominal,160,50}, ...
    'ratio',{0,0,0.12},'alpha',{0,0,0.15});
if abs(required_nominal-reference_nominal)>eps
    specs(end+1)=struct('name','required_no_sponge','nominal_width_m',required_nominal,'ratio',0,'alpha',0);
end
specs=local_unique_specs(specs);
cases=repmat(local_empty_case(),numel(specs),1);
for ii=1:numel(specs)
    checkpoint=fullfile(cfg.output_dir,sprintf('wideband_case_%s.mat',specs(ii).name));
    if cfg.resume_cases && exist(checkpoint,'file')==2
        saved=load(checkpoint,'wideband_case');
        cases(ii)=saved.wideband_case;
    else
        cases(ii)=local_run_case(specs(ii),cfg,master.master_eta);
        wideband_case=cases(ii);
        save(checkpoint,'wideband_case','-v7.3');
    end
end
iref=find(strcmp({cases.name},'reference_no_sponge'),1);
if isempty(iref), [~,iref]=max([cases.width_actual_m]); end
[wideband_table,summary_table]=local_compare(cases,iref,cfg);
writetable(wideband_table,fullfile(cfg.output_dir,'wideband_reflected_total_comparison.csv'));
writetable(summary_table,fullfile(cfg.output_dir,'wideband_summary.csv'));

sponge_table=table(); sponge_recommendation=struct('use_sponge',false,'reason','No scan required.');
candidate_index=find(strcmp({cases.name},'required_no_sponge'),1);
if isempty(candidate_index), candidate_index=iref; end
candidate=cases(candidate_index);
if ~gate.cases(gate.required_index).reference_edge_valid
    [sponge_table,sponge_recommendation]=local_sponge_scan(candidate,cases(iref),cfg,master.master_eta);
    writetable(sponge_table,fullfile(cfg.output_dir,'weak_sponge_scan.csv'));
end

validation=struct('schema_version','1.0.0','config',cfg,'window_gate',gate, ...
    'master_surface_meta',master.master_meta,'cases',cases,'reference_index',iref, ...
    'wideband_table',wideband_table,'summary_table',summary_table, ...
    'sponge_table',sponge_table,'sponge_recommendation',sponge_recommendation, ...
    'recommendation',local_recommendation(gate,sponge_recommendation));
save(fullfile(cfg.output_dir,'reflected_chain_wideband_validation.mat'),'validation','-v7.3');
local_plot(validation);
local_append_report(validation,root);
end

function cfg=local_config(root,o)
cfg=struct('output_dir',fullfile(root,'results','validation','pe_reflected_chain_window'), ...
    'f_axis_hz',linspace(3000,5000,33),'f_ref_hz',4000,'c0_mps',1500, ...
    'z_tx_m',100,'z_rx_m',3,'sigma_src_m',0.3,'stepz_lamb',0.5,'dx_m',50/256, ...
    'wind_speed_mps',5,'target_hs_m',0.5,'sea_seed',12345,'reflect_coeff',-1, ...
    'tl_limit_db',0.1,'phase_limit_rad',0.02,'group_delay_limit_s',20e-6, ...
    'center_l2_limit',0.01,'use_gpu',true,'resume_cases',true, ...
    'weak_ratios',[0.05 0.08 0.10 0.12], ...
    'weak_alphas',[0.005 0.01 0.025 0.05],'edge_reduction_target_db',10);
n=fieldnames(o); for ii=1:numel(n), cfg.(n{ii})=o.(n{ii}); end
end

function c=local_run_case(s,cfg,master)
n=local_n(s.nominal_width_m,cfg.dx_m); w=n*cfg.dx_m; eta=local_crop(master,n);
p=local_params(cfg,w,n,eta,s.ratio,s.alpha);
out=vertical_channel_model(p);
c=local_empty_case(); c.name=s.name; c.nominal_width_m=s.nominal_width_m; c.width_actual_m=w;
c.ratio=s.ratio; c.alpha=s.alpha; c.f_axis_hz=out.f_axis(:);
c.H_reflect_physical_f=out.H_reflect_physical_f(:); c.H_total_physical_f=out.H_physical_f(:);
c.H_direct_f=out.H_direct_f(:); c.H_reflect_f=out.H_reflect_f(:); c.H_f=out.H_f(:);
c.closure_max_abs=max(abs(c.H_f-c.H_direct_f-c.H_reflect_f));
end

function p=local_params(c,w,n,eta,ratio,alpha)
p=struct('f0',c.f_axis_hz,'f_ref_hz',c.f_ref_hz,'c0',c.c0_mps, ...
    'z_max',c.z_tx_m,'z_tx',c.z_tx_m,'z_rx',c.z_rx_m,'x_tx',0,'y_tx',0,'x_rx',0,'y_rx',0, ...
    'xw',w,'yw',w,'nx',n,'ny',n,'source_mode','gaussian','sigma_src_m',c.sigma_src_m, ...
    'stepz_lamb',c.stepz_lamb,'env_mode','uniform','show_figures',false, ...
    'save_mode','rx_only','enforce_1_over_R',false,'enable_surface_reflection',true, ...
    'surface_boundary_model','kirchhoff_spatial','surface_reflect_coeff',c.reflect_coeff, ...
    'surface_phase_mode','normal','sea_wind_speed',c.wind_speed_mps,'sea_hs_target',c.target_hs_m, ...
    'surface_elevation_override_xy',eta,'surface_wavefield_diagnostics',false, ...
    'validation_allow_extended_window',w>100,'validation_allow_extended_sponge_ratio',true, ...
    'sponge_ratio',ratio,'alpha_max_np_per_m',alpha,'enable_bubbles',false, ...
    'use_gpu',c.use_gpu);
end

function [rows,summary]=local_compare(cases,iref,cfg)
ref=cases(iref); n=numel(cases)*numel(ref.f_axis_hz);
r=repmat(struct('case_name',"",'frequency_hz',NaN,'reflect_magnitude_db',NaN, ...
    'reflect_phase_rad',NaN,'reflect_group_delay_s',NaN,'reflect_tl_error_db',NaN, ...
    'reflect_phase_error_rad',NaN,'reflect_group_delay_error_s',NaN, ...
    'total_magnitude_db',NaN,'total_tl_error_db',NaN,'total_phase_error_rad',NaN, ...
    'dsp_total_magnitude_db',NaN,'dsp_total_tl_error_db',NaN, ...
    'dsp_total_phase_error_rad',NaN),n,1); q=0;
s=repmat(struct('case_name',"",'width_actual_m',NaN,'max_reflect_tl_error_db',NaN, ...
    'max_reflect_phase_error_rad',NaN,'max_internal_group_delay_error_s',NaN, ...
    'max_total_tl_error_db',NaN,'max_total_phase_error_rad',NaN, ...
    'max_dsp_total_tl_error_db',NaN,'max_dsp_total_phase_error_rad',NaN, ...
    'closure_max_abs',NaN,'passed',false),numel(cases),1);
phase_ref=unwrap(angle(ref.H_reflect_physical_f)); gd_ref=gradient(phase_ref,ref.f_axis_hz)/(2*pi);
for ii=1:numel(cases)
    c=cases(ii); ph=unwrap(angle(c.H_reflect_physical_f)); gd=gradient(ph,c.f_axis_hz)/(2*pi);
    tl=-20*log10(abs(c.H_reflect_physical_f./ref.H_reflect_physical_f));
    pe=unwrap(angle(c.H_reflect_physical_f.*conj(ref.H_reflect_physical_f)));
    ge=gd-gd_ref; ttl=-20*log10(abs(c.H_total_physical_f./ref.H_total_physical_f));
    tpe=unwrap(angle(c.H_total_physical_f.*conj(ref.H_total_physical_f)));
    dsp_tl=-20*log10(abs(c.H_f./ref.H_f));
    dsp_pe=unwrap(angle(c.H_f.*conj(ref.H_f)));
    for ff=1:numel(c.f_axis_hz)
        q=q+1; r(q).case_name=string(c.name); r(q).frequency_hz=c.f_axis_hz(ff);
        r(q).reflect_magnitude_db=20*log10(abs(c.H_reflect_physical_f(ff)));
        r(q).reflect_phase_rad=ph(ff); r(q).reflect_group_delay_s=gd(ff);
        r(q).reflect_tl_error_db=tl(ff); r(q).reflect_phase_error_rad=pe(ff);
        r(q).reflect_group_delay_error_s=ge(ff); r(q).total_magnitude_db=20*log10(abs(c.H_total_physical_f(ff)));
        r(q).total_tl_error_db=ttl(ff); r(q).total_phase_error_rad=tpe(ff);
        r(q).dsp_total_magnitude_db=20*log10(abs(c.H_f(ff)));
        r(q).dsp_total_tl_error_db=dsp_tl(ff); r(q).dsp_total_phase_error_rad=dsp_pe(ff);
    end
    internal=2:numel(ge)-1;
    s(ii).case_name=string(c.name); s(ii).width_actual_m=c.width_actual_m;
    s(ii).max_reflect_tl_error_db=max(abs(tl)); s(ii).max_reflect_phase_error_rad=max(abs(pe));
    s(ii).max_internal_group_delay_error_s=max(abs(ge(internal))); s(ii).closure_max_abs=c.closure_max_abs;
    s(ii).max_total_tl_error_db=max(abs(ttl)); s(ii).max_total_phase_error_rad=max(abs(tpe));
    s(ii).max_dsp_total_tl_error_db=max(abs(dsp_tl));
    s(ii).max_dsp_total_phase_error_rad=max(abs(dsp_pe));
    s(ii).passed=s(ii).max_reflect_tl_error_db<=cfg.tl_limit_db && ...
        s(ii).max_reflect_phase_error_rad<=cfg.phase_limit_rad && ...
        s(ii).max_internal_group_delay_error_s<=cfg.group_delay_limit_s;
end
rows=struct2table(r); summary=struct2table(s);
end

function [tbl,recommendation]=local_sponge_scan(candidate,ref,cfg,master)
r=repmat(struct('ratio',NaN,'alpha_max_np_per_m',NaN,'max_tl_error_db',NaN, ...
    'max_phase_error_rad',NaN,'max_internal_group_delay_error_s',NaN,'passed',false), ...
    numel(cfg.weak_ratios)*numel(cfg.weak_alphas),1); q=0;
for ratio=cfg.weak_ratios
    for alpha=cfg.weak_alphas
        q=q+1; s=struct('name','weak_sponge','nominal_width_m',candidate.nominal_width_m,'ratio',ratio,'alpha',alpha);
        c=local_run_case(s,cfg,master); ph=unwrap(angle(c.H_reflect_physical_f.*conj(ref.H_reflect_physical_f)));
        tl=-20*log10(abs(c.H_reflect_physical_f./ref.H_reflect_physical_f));
        gd=gradient(unwrap(angle(c.H_reflect_physical_f)),c.f_axis_hz)/(2*pi);
        gd0=gradient(unwrap(angle(ref.H_reflect_physical_f)),ref.f_axis_hz)/(2*pi); internal=2:numel(gd)-1;
        r(q).ratio=ratio; r(q).alpha_max_np_per_m=alpha; r(q).max_tl_error_db=max(abs(tl));
        r(q).max_phase_error_rad=max(abs(ph)); r(q).max_internal_group_delay_error_s=max(abs(gd(internal)-gd0(internal)));
        r(q).passed=r(q).max_tl_error_db<=cfg.tl_limit_db && r(q).max_phase_error_rad<=cfg.phase_limit_rad && r(q).max_internal_group_delay_error_s<=cfg.group_delay_limit_s;
    end
end
tbl=struct2table(r); ok=find(tbl.passed,1);
if isempty(ok), recommendation=struct('use_sponge',false,'reason','No weak sponge met the fixed center-channel gates.');
else, recommendation=struct('use_sponge',true,'ratio',tbl.ratio(ok),'alpha_max_np_per_m',tbl.alpha_max_np_per_m(ok),'reason','Weakest scanned sponge meeting wideband gates; edge-reduction must also be confirmed from single-frequency diagnostics.'); end
end

function r=local_recommendation(gate,sp)
if isfield(sp,'use_sponge') && sp.use_sponge
    r=struct('window_m',gate.required_width_m,'use_sponge',true,'ratio',sp.ratio,'alpha_max_np_per_m',sp.alpha_max_np_per_m);
else
    r=struct('window_m',gate.required_width_m,'use_sponge',false,'ratio',0,'alpha_max_np_per_m',0);
end
end

function local_plot(v)
out=fullfile(v.config.output_dir,'figures'); if ~exist(out,'dir'),mkdir(out);end
f=figure('Visible','off','Color','w','Position',[100 100 1150 760]); t=tiledlayout(3,1,'TileSpacing','compact');
nexttile; hold on; grid on; for ii=1:numel(v.cases), c=v.cases(ii); plot(c.f_axis_hz/1000,20*log10(abs(c.H_reflect_physical_f)),'DisplayName',c.name); end; ylabel('|H reflect| (dB)'); legend('Location','best');
nexttile; hold on; grid on; for ii=1:numel(v.cases), c=v.cases(ii); plot(c.f_axis_hz/1000,unwrap(angle(c.H_reflect_physical_f)),'DisplayName',c.name); end; ylabel('phase (rad)');
nexttile; hold on; grid on; for ii=1:numel(v.cases), c=v.cases(ii); gd=gradient(unwrap(angle(c.H_reflect_physical_f)),c.f_axis_hz)/(2*pi); plot(c.f_axis_hz/1000,gd*1e6,'DisplayName',c.name); end; ylabel('group delay (us)'); xlabel('frequency (kHz)'); title(t,'Wideband reflected-path convergence'); exportgraphics(f,fullfile(out,'05_wideband_reflected_response.png'),'Resolution',180); close(f);
f=figure('Visible','off','Color','w'); hold on; grid on; for ii=1:numel(v.cases), c=v.cases(ii); plot(c.f_axis_hz/1000,20*log10(abs(c.H_total_physical_f)),'DisplayName',c.name); end; xlabel('frequency (kHz)'); ylabel('|H total physical| (dB)'); legend('Location','best'); title('Total-channel window stability'); exportgraphics(f,fullfile(out,'06_wideband_total_response.png'),'Resolution',180); close(f);
end

function local_append_report(v,root)
path=fullfile(root,'reports','pe_reflected_chain_window_validation_report.md'); fid=fopen(path,'a'); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'\n## Wideband qualification\n\n');
fprintf(fid,'| case | width (m) | max reflect TL (dB) | max reflect phase (rad) | max internal GD (us) | max total TL (dB) | max total phase (rad) | pass |\n');
fprintf(fid,'|---|---:|---:|---:|---:|---:|---:|:---:|\n');
for ii=1:height(v.summary_table)
    row=v.summary_table(ii,:);
    fprintf(fid,'| %s | %.6f | %.6g | %.6g | %.6g | %.6g | %.6g | %s |\n', ...
        string(row.case_name),row.width_actual_m,row.max_reflect_tl_error_db, ...
        row.max_reflect_phase_error_rad,1e6*row.max_internal_group_delay_error_s, ...
        row.max_total_tl_error_db,row.max_total_phase_error_rad, ...
        local_word(row.passed,'PASS','FAIL'));
end
fprintf(fid,'\nThe physical-total and communication `H_f` comparison columns are both saved in `wideband_reflected_total_comparison.csv`; their window errors agree under the common direct-DSP phase reference. ');
fprintf(fid,'Recommended configuration: W=%.6g m, sponge=%s.\n',v.recommendation.window_m,local_word(v.recommendation.use_sponge,'on','off'));
end
function s=local_word(tf,a,b),if tf,s=a;else,s=b;end,end
function s=local_unique_specs(s), key=strcat(string([s.nominal_width_m]),"_",string([s.ratio]),"_",string([s.alpha])); [~,i]=unique(key,'stable'); s=s(sort(i)); end
function n=local_n(w,dx),n=2*round(w/dx/2);end
function c=local_crop(a,n),s=(size(a,1)-n)/2+1;c=a(s:s+n-1,s:s+n-1);end
function c=local_empty_case(),c=struct('name','','nominal_width_m',NaN,'width_actual_m',NaN,'ratio',NaN,'alpha',NaN,'f_axis_hz',[],'H_reflect_physical_f',[],'H_total_physical_f',[],'H_direct_f',[],'H_reflect_f',[],'H_f',[],'closure_max_abs',NaN);end
