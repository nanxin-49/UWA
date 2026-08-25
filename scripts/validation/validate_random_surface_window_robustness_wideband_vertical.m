function validation = validate_random_surface_window_robustness_wideband_vertical(overrides)
%VALIDATE_RANDOM_SURFACE_WINDOW_ROBUSTNESS_WIDEBAND_VERTICAL Selected-case audit.
if nargin<1, overrides=struct(); end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
run(fullfile(root,'scripts','bootstrap_project.m'));
cfg=local_config(root,overrides);
gate_path=fullfile(cfg.output_dir,'random_surface_window_robustness_4khz.mat');
if exist(gate_path,'file')~=2, error('The formal 4 kHz random-surface gate is missing.'); end
d=load(gate_path,'validation'); gate=d.validation;
if ~gate.all_references_valid, error('Wideband is blocked because a 256 m reference is edge-invalid.'); end
if ~exist(cfg.checkpoint_dir,'dir'),mkdir(cfg.checkpoint_dir);end

specs=local_specs(gate.selection_table);
cases=repmat(local_empty_case(),numel(specs),1);
for ii=1:numel(specs)
    checkpoint=local_checkpoint(cfg,specs(ii));
    fprintf('\nWideband random-surface case: %s, seed=%d, W=%g m.\n', ...
        specs(ii).sea_state,specs(ii).seed,specs(ii).window_nominal_m);
    if cfg.resume_cases && exist(checkpoint,'file')==2
        s=load(checkpoint,'wideband_case'); cases(ii)=s.wideband_case;
    else
        cases(ii)=local_run_case(specs(ii),gate,cfg,checkpoint);
        wideband_case=cases(ii); save(checkpoint,'wideband_case','-v7.3');
    end
end

[frequency_table,summary_table]=local_compare(cases,cfg);
writetable(frequency_table,fullfile(cfg.output_dir,'wideband_selected_frequency_comparison.csv'));
writetable(summary_table,fullfile(cfg.output_dir,'wideband_selected_summary.csv'));
validation=struct('schema_version','1.0.0','config',cfg,'gate',gate, ...
    'specs',specs,'cases',cases,'frequency_table',frequency_table, ...
    'summary_table',summary_table,'recommendation',local_recommendation(gate,summary_table,cfg));
save(fullfile(cfg.output_dir,'random_surface_window_robustness_wideband.mat'),'validation','-v7.3');
local_plot(validation); local_append_report(validation);
end

function cfg=local_config(root,o)
cfg=struct('output_dir',fullfile(root,'results','validation','pe_random_surface_window_robustness'), ...
    'checkpoint_dir',fullfile(root,'results','validation','pe_random_surface_window_robustness','checkpoints_wideband'), ...
    'report_path',fullfile(root,'reports','pe_random_surface_window_robustness_report.md'), ...
    'f_axis_hz',linspace(3000,5000,33),'edge_probe_hz',[3000 5000], ...
    'tl_limit_db',0.1,'phase_limit_rad',0.02,'group_delay_limit_s',20e-6, ...
    'center_l2_limit',0.01,'edge5_fraction_limit',1e-4,'boundary_limit_db',-40, ...
    'use_gpu',true,'resume_cases',true,'max_z_samples',17,'frequency_chunk_size',4);
n=fieldnames(o); for ii=1:numel(n),cfg.(n{ii})=o.(n{ii});end
if isfield(o,'output_dir') && ~isfield(o,'checkpoint_dir'),cfg.checkpoint_dir=fullfile(cfg.output_dir,'checkpoints_wideband');end
if isfield(o,'output_dir') && ~isfield(o,'report_path'),cfg.report_path=fullfile(cfg.output_dir,'reduced_random_surface_report.md');end
end

function specs=local_specs(selection)
specs=repmat(struct('sea_state','','target_hs_m',NaN,'wind_speed_mps',NaN, ...
    'seed',NaN,'window_nominal_m',NaN,'role','','reason',''),0,1);
for ii=1:height(selection)
    base=struct('sea_state',char(selection.sea_state(ii)),'target_hs_m',selection.target_hs_m(ii), ...
        'wind_speed_mps',selection.wind_speed_mps(ii),'seed',selection.seed(ii), ...
        'window_nominal_m',cfg_reference_width(),'role','reference','reason',char(selection.reason(ii)));
    specs=local_add_spec(specs,base);
    if contains(selection.reason(ii),'160 m')
        target=160;
    else
        target=192;
    end
    base.window_nominal_m=target; base.role='candidate'; specs=local_add_spec(specs,base);
end
end

function s=local_add_spec(s,r)
if isempty(s),s=r;return,end
key=strcmp({s.sea_state},r.sea_state)&[s.seed]==r.seed&[s.window_nominal_m]==r.window_nominal_m;
if ~any(key),s(end+1,1)=r;end
end

function c=local_run_case(spec,gate,cfg,checkpoint)
master_path=fullfile(cfg.output_dir,sprintf('master_surface_%s_seed%d.mat', ...
    matlab.lang.makeValidName(spec.sea_state),spec.seed));
d=load(master_path,'master_eta'); n=local_n(spec.window_nominal_m,gate.config.dx_m);
eta=local_crop(d.master_eta,n); w=n*gate.config.dx_m;
c=local_empty_case(); copy={'sea_state','target_hs_m','wind_speed_mps','seed','window_nominal_m','role','reason'};
for ii=1:numel(copy),c.(copy{ii})=spec.(copy{ii});end
c.window_actual_m=w;
partial_path=replace(checkpoint,'.mat','_partial.mat');
if cfg.resume_cases && exist(partial_path,'file')==2
    d=load(partial_path,'wideband_case_partial'); c=d.wideband_case_partial;
end
done=c.f_axis_hz(:).'; requested=cfg.f_axis_hz(:).'; pending=requested(~ismember(requested,done));
while ~isempty(pending)
    take=pending(1:min(cfg.frequency_chunk_size,numel(pending)));
    p=local_params(gate.config,spec,w,n,eta,take,false,cfg); out=vertical_channel_model(p);
    c.f_axis_hz=[c.f_axis_hz(:);out.f_axis(:)];
    c.H_reflect_physical_f=[c.H_reflect_physical_f(:);out.H_reflect_physical_f(:)];
    c.H_total_physical_f=[c.H_total_physical_f(:);out.H_physical_f(:)];
    c.H_direct_f=[c.H_direct_f(:);out.H_direct_f(:)];
    c.H_reflect_f=[c.H_reflect_f(:);out.H_reflect_f(:)]; c.H_f=[c.H_f(:);out.H_f(:)];
    [c.f_axis_hz,order]=sort(c.f_axis_hz); fields={'H_reflect_physical_f','H_total_physical_f','H_direct_f','H_reflect_f','H_f'};
    for jj=1:numel(fields),c.(fields{jj})=c.(fields{jj})(order);end
    wideband_case_partial=c; save(partial_path,'wideband_case_partial','-v7.3');
    done=c.f_axis_hz(:).'; pending=requested(~ismember(requested,done)); clear out
end
c.closure_max_abs=max(abs(c.H_f-c.H_direct_f-c.H_reflect_f));
if isempty(c.edge_probe)
    c.edge_probe=local_edge_probe(gate.config,spec,w,n,eta,cfg);
    wideband_case_partial=c; save(partial_path,'wideband_case_partial','-v7.3');
end
end

function probe=local_edge_probe(base,spec,w,n,eta,cfg)
probe=repmat(struct('frequency_hz',NaN,'max_edge5_fraction',NaN,'max_boundary_db',NaN, ...
    'upward_max_edge5_fraction',NaN,'downward_max_edge5_fraction',NaN),numel(cfg.edge_probe_hz)+1,1);
freq=[cfg.edge_probe_hz(:);base.frequency_hz];
for ii=1:numel(freq)
    if freq(ii)==base.frequency_hz
        gate_path=fullfile(cfg.output_dir,'random_surface_window_robustness_4khz.mat');
        d=load(gate_path,'validation'); q=d.validation.comparison_table;
        row=q(strcmp(q.sea_state,spec.sea_state)&q.seed==spec.seed&q.window_nominal_m==spec.window_nominal_m,:);
        probe(ii)=struct('frequency_hz',freq(ii),'max_edge5_fraction',row.max_edge5_fraction, ...
            'max_boundary_db',row.max_boundary_db,'upward_max_edge5_fraction',row.upward_max_edge5_fraction, ...
            'downward_max_edge5_fraction',row.downward_max_edge5_fraction);
    else
        p=local_params(base,spec,w,n,eta,freq(ii),true,cfg); o=vertical_channel_model(p);
        a=o.surface_wavefield_meta.incident_energy_trace; b=o.surface_wavefield_meta.reflected_energy_trace;
        probe(ii)=struct('frequency_hz',freq(ii),'max_edge5_fraction',max([a.edge5_fraction,b.edge5_fraction]), ...
            'max_boundary_db',max([a.boundary_max_relative_db,b.boundary_max_relative_db]), ...
            'upward_max_edge5_fraction',max(a.edge5_fraction), ...
            'downward_max_edge5_fraction',max(b.edge5_fraction));
    end
end
end

function p=local_params(base,spec,w,n,eta,freq,diagnostics,cfg)
p=struct('f0',freq,'f_ref_hz',4000,'c0',base.c0_mps,'z_max',base.z_tx_m, ...
    'z_tx',base.z_tx_m,'z_rx',base.z_rx_m,'x_tx',0,'y_tx',0,'x_rx',0,'y_rx',0, ...
    'xw',w,'yw',w,'nx',n,'ny',n,'source_mode','gaussian','sigma_src_m',base.sigma_src_m, ...
    'stepz_lamb',base.stepz_lamb,'env_mode','uniform','show_figures',false, ...
    'save_mode','rx_only','enforce_1_over_R',false,'enable_surface_reflection',true, ...
    'surface_boundary_model','kirchhoff_spatial','surface_reflect_coeff',base.reflect_coeff, ...
    'surface_phase_mode','normal','sea_wind_speed',spec.wind_speed_mps, ...
    'sea_hs_target',spec.target_hs_m,'surface_roughness_scale_mode','target_hs', ...
    'sea_seed',spec.seed,'surface_elevation_override_xy',eta, ...
    'surface_wavefield_diagnostics',diagnostics,'surface_wavefield_max_z_samples',cfg.max_z_samples, ...
    'validation_allow_extended_window',true,'validation_allow_extended_sponge_ratio',true, ...
    'sponge_ratio',0,'alpha_max_np_per_m',0,'enable_bubbles',false,'use_gpu',cfg.use_gpu);
if isscalar(freq),p.f_ref_hz=freq;end
end

function [rows,summary]=local_compare(cases,cfg)
rows=repmat(struct('sea_state',"",'seed',NaN,'window_nominal_m',NaN,'frequency_hz',NaN, ...
    'reflect_tl_error_db',NaN,'reflect_phase_error_rad',NaN,'group_delay_error_s',NaN, ...
    'total_tl_error_db',NaN,'total_phase_error_rad',NaN,'dsp_total_tl_error_db',NaN, ...
    'dsp_total_phase_error_rad',NaN),0,1);
summary=repmat(struct('sea_state',"",'target_hs_m',NaN,'seed',NaN,'window_nominal_m',NaN, ...
    'max_reflect_tl_error_db',NaN,'max_reflect_phase_error_rad',NaN, ...
    'max_internal_group_delay_error_s',NaN,'max_total_tl_error_db',NaN, ...
    'max_total_phase_error_rad',NaN,'max_dsp_total_tl_error_db',NaN, ...
    'max_dsp_total_phase_error_rad',NaN,'max_edge5_fraction',NaN,'max_boundary_db',NaN, ...
    'closure_max_abs',NaN,'response_passed',false,'edge_passed',false,'passed',false),0,1);
candidates=find(strcmp({cases.role},'candidate'));
for jj=1:numel(candidates)
    c=cases(candidates(jj)); ir=find(strcmp({cases.role},'reference')&strcmp({cases.sea_state},c.sea_state)&[cases.seed]==c.seed,1);
    ref=cases(ir); phase=unwrap(angle(c.H_reflect_physical_f.*conj(ref.H_reflect_physical_f)));
    tl=-20*log10(abs(c.H_reflect_physical_f./ref.H_reflect_physical_f));
    gd=gradient(unwrap(angle(c.H_reflect_physical_f)),c.f_axis_hz)/(2*pi);
    gd0=gradient(unwrap(angle(ref.H_reflect_physical_f)),ref.f_axis_hz)/(2*pi); ge=gd-gd0;
    ttl=-20*log10(abs(c.H_total_physical_f./ref.H_total_physical_f));
    tph=unwrap(angle(c.H_total_physical_f.*conj(ref.H_total_physical_f)));
    dtl=-20*log10(abs(c.H_f./ref.H_f)); dph=unwrap(angle(c.H_f.*conj(ref.H_f)));
    for ff=1:numel(c.f_axis_hz)
        r=struct('sea_state',string(c.sea_state),'seed',c.seed,'window_nominal_m',c.window_nominal_m, ...
            'frequency_hz',c.f_axis_hz(ff),'reflect_tl_error_db',tl(ff), ...
            'reflect_phase_error_rad',phase(ff),'group_delay_error_s',ge(ff), ...
            'total_tl_error_db',ttl(ff),'total_phase_error_rad',tph(ff), ...
            'dsp_total_tl_error_db',dtl(ff),'dsp_total_phase_error_rad',dph(ff));
        rows(end+1,1)=r; %#ok<AGROW>
    end
    internal=2:numel(ge)-1; edge=max([c.edge_probe.max_edge5_fraction]); boundary=max([c.edge_probe.max_boundary_db]);
    s=summary_template(); s.sea_state=string(c.sea_state);s.target_hs_m=c.target_hs_m;s.seed=c.seed;s.window_nominal_m=c.window_nominal_m;
    s.max_reflect_tl_error_db=max(abs(tl));s.max_reflect_phase_error_rad=max(abs(phase));
    s.max_internal_group_delay_error_s=max(abs(ge(internal)));s.max_total_tl_error_db=max(abs(ttl));
    s.max_total_phase_error_rad=max(abs(tph));s.max_dsp_total_tl_error_db=max(abs(dtl));
    s.max_dsp_total_phase_error_rad=max(abs(dph));s.max_edge5_fraction=edge;s.max_boundary_db=boundary;
    s.closure_max_abs=max(c.closure_max_abs,ref.closure_max_abs);
    s.response_passed=s.max_reflect_tl_error_db<=cfg.tl_limit_db&&s.max_reflect_phase_error_rad<=cfg.phase_limit_rad&&s.max_internal_group_delay_error_s<=cfg.group_delay_limit_s;
    s.edge_passed=edge<=cfg.edge5_fraction_limit&&boundary<=cfg.boundary_limit_db;s.passed=s.response_passed&&s.edge_passed;
    summary(end+1,1)=s; %#ok<AGROW>
end
rows=struct2table(rows);summary=struct2table(summary);
end

function r=local_recommendation(gate,summary,cfg)
q192=gate.comparison_table(gate.comparison_table.window_nominal_m==192,:);
q160=gate.comparison_table(gate.comparison_table.window_nominal_m==160,:);
wb192=summary(summary.window_nominal_m==192,:); wb160=summary(summary.window_nominal_m==160,:);
if all(q192.passed)&&all(wb192.passed)
    choice='A';text='192 m + no sponge passed every 4 kHz realization and all selected wideband stress cases.';
elseif any(~q192.passed)||any(~wb192.passed)
    choice='D';text='192 m failed at least one strict random-surface or selected-wideband gate.';
elseif all(q160.passed)&&(~isempty(wb160)&&all(wb160.passed))
    choice='B';text='160 m passed all available gates.';
else
    choice='C';text='Use the automated outer-energy and boundary-amplitude gates to select the window.';
end
r=struct('choice',choice,'text',text,'edge5_limit',cfg.edge5_fraction_limit,'boundary_limit_db',cfg.boundary_limit_db);
end

function local_plot(v)
out=fullfile(v.config.output_dir,'figures');if ~exist(out,'dir'),mkdir(out);end
f=figure('Visible','off','Color','w','Position',[100 100 1200 760]);t=tiledlayout(3,1,'TileSpacing','compact');
for metric={'reflect_tl_error_db','reflect_phase_error_rad','group_delay_error_s'}
    nexttile;hold on;grid on
    for ii=1:height(v.summary_table)
        q=v.frequency_table(strcmp(v.frequency_table.sea_state,v.summary_table.sea_state(ii))&v.frequency_table.seed==v.summary_table.seed(ii)&v.frequency_table.window_nominal_m==v.summary_table.window_nominal_m(ii),:);
        y=q.(metric{1});if strcmp(metric{1},'group_delay_error_s'),y=y*1e6;end
        plot(q.frequency_hz/1000,y,'DisplayName',sprintf('%s s%d W%d',v.summary_table.sea_state(ii),v.summary_table.seed(ii),v.summary_table.window_nominal_m(ii)));
    end
    ylabel(strrep(metric{1},'_',' '));
end
xlabel('frequency (kHz)');legend('Location','best');title(t,'Selected random-surface wideband errors versus 256 m');
exportgraphics(f,fullfile(out,'05_selected_wideband_errors.png'),'Resolution',180);close(f);
end

function local_append_report(v)
fid=fopen(v.config.report_path,'a');cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'\n## Selected 3--5 kHz qualification\n\n');
fprintf(fid,'| sea | seed | W (m) | max TL (dB) | max phase (rad) | max internal GD (us) | max edge5 | boundary (dB) | pass |\n');
fprintf(fid,'|---|---:|---:|---:|---:|---:|---:|---:|:---:|\n');
for ii=1:height(v.summary_table)
    r=v.summary_table(ii,:);fprintf(fid,'| %s | %d | %.0f | %.6g | %.6g | %.6g | %.6g | %.6g | %s |\n', ...
        string(r.sea_state),r.seed,r.window_nominal_m,r.max_reflect_tl_error_db, ...
        r.max_reflect_phase_error_rad,1e6*r.max_internal_group_delay_error_s, ...
        r.max_edge5_fraction,r.max_boundary_db,local_word(r.passed,'PASS','FAIL'));
end
fprintf(fid,'\n## Final Q1--Q6 decision\n\nRecommendation: **%s**. %s\n',v.recommendation.choice,v.recommendation.text);
fprintf(fid,'\nAutomated dynamic fallback gate: use the candidate only if outer-5%% energy is below %.1e and boundary amplitude is no higher than %.0f dB, in addition to the fixed response/center-field gates.\n',v.recommendation.edge5_limit,v.recommendation.boundary_limit_db);
end

function p=local_checkpoint(cfg,s)
p=fullfile(cfg.checkpoint_dir,sprintf('%s_seed%d_W%d.mat',matlab.lang.makeValidName(s.sea_state),s.seed,round(s.window_nominal_m)));
end
function n=local_n(w,dx),n=2*round(w/dx/2);end
function c=local_crop(a,n),s=(size(a,1)-n)/2+1;c=a(s:s+n-1,s:s+n-1);end
function w=cfg_reference_width(),w=256;end
function c=local_empty_case()
c=struct('sea_state','','target_hs_m',NaN,'wind_speed_mps',NaN,'seed',NaN, ...
    'window_nominal_m',NaN,'window_actual_m',NaN,'role','','reason','', ...
    'f_axis_hz',[],'H_reflect_physical_f',[],'H_total_physical_f',[], ...
    'H_direct_f',[],'H_reflect_f',[],'H_f',[],'closure_max_abs',NaN,'edge_probe',struct([]));
end
function s=summary_template()
s=struct('sea_state',"",'target_hs_m',NaN,'seed',NaN,'window_nominal_m',NaN, ...
    'max_reflect_tl_error_db',NaN,'max_reflect_phase_error_rad',NaN, ...
    'max_internal_group_delay_error_s',NaN,'max_total_tl_error_db',NaN, ...
    'max_total_phase_error_rad',NaN,'max_dsp_total_tl_error_db',NaN, ...
    'max_dsp_total_phase_error_rad',NaN,'max_edge5_fraction',NaN,'max_boundary_db',NaN, ...
    'closure_max_abs',NaN,'response_passed',false,'edge_passed',false,'passed',false);
end
function s=local_word(tf,a,b),if tf,s=a;else,s=b;end,end
