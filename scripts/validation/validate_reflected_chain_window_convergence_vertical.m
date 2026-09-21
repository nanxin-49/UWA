function validation = validate_reflected_chain_window_convergence_vertical(overrides)
%VALIDATE_REFLECTED_CHAIN_WINDOW_CONVERGENCE_VERTICAL Fixed-surface window audit.
% The default run is the preregistered production-scale 4 kHz validation.

if nargin < 1, overrides = struct(); end
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
run(fullfile(root,'scripts','bootstrap_project.m'));
cfg = local_config(root,overrides);
if ~exist(cfg.output_dir,'dir'), mkdir(cfg.output_dir); end

[master_eta,master_meta] = local_master_surface(cfg);
master_file = fullfile(cfg.output_dir,'fixed_master_surface.mat');
save(master_file,'master_eta','master_meta','-v7.3');

nominal = cfg.initial_windows_m(:).';
initial_count = numel(nominal);
escalation = cfg.escalation_windows_m(:).';
cases = repmat(local_empty_case(),0,1);
pair_rows = repmat(local_empty_pair(),0,1);
next_windows = nominal;
converged = false;

while ~isempty(next_windows)
    width_nominal = next_windows(1);
    next_windows(1) = [];
    n = local_n(width_nominal,cfg.dx_m);
    width_actual = n*cfg.dx_m;
    eta_xy = local_center_crop(master_eta,n);
    fprintf('\nReflected-chain window: nominal %.3f m, actual %.6f m, N=%d.\n', ...
        width_nominal,width_actual,n);
    out = vertical_channel_model(local_params(cfg,width_actual,n,eta_xy,0,0));
    cases(end+1,1) = local_extract_case(out,width_nominal,eta_xy,cfg); %#ok<AGROW>
    clear out eta_xy

    if numel(cases) >= 2
        pair_rows(end+1,1) = local_compare_pair(cases(end-1),cases(end),cfg); %#ok<AGROW>
        if numel(cases)>=initial_count && local_pair_passes(pair_rows(end),cfg) && ...
                cases(end).reference_edge_valid
            converged = true;
            break
        end
    end
    if isempty(next_windows) && ~isempty(escalation)
        next_windows = escalation(1);
        escalation(1) = [];
    end
end

if converged
    reference_index = numel(cases);
    reference_width_m = cases(reference_index).width_actual_m;
    metric_rows = local_reference_metrics(cases,reference_index,cfg);
    [required_index,required_width_m] = local_required_window(cases,metric_rows,cfg);
else
    reference_index = NaN;
    reference_width_m = NaN;
    required_index = NaN;
    required_width_m = NaN;
    metric_rows = repmat(local_empty_metric(),0,1);
end

case_table = struct2table(local_case_rows(cases));
pair_table = struct2table(pair_rows);
metric_table = struct2table(metric_rows);
energy_table = local_energy_table(cases);
writetable(case_table,fullfile(cfg.output_dir,'single_frequency_cases.csv'));
writetable(pair_table,fullfile(cfg.output_dir,'adjacent_window_convergence.csv'));
writetable(metric_table,fullfile(cfg.output_dir,'reference_window_metrics.csv'));
writetable(energy_table,fullfile(cfg.output_dir,'stage_energy_traces.csv'));

validation = struct('schema_version','1.0.0','config',cfg, ...
    'master_surface_meta',master_meta,'cases',cases,'case_table',case_table, ...
    'pair_table',pair_table,'metric_table',metric_table,'energy_table',energy_table, ...
    'converged',converged,'reference_index',reference_index, ...
    'reference_width_m',reference_width_m,'required_index',required_index, ...
    'required_width_m',required_width_m,'width_160_sufficient',false, ...
    'recommendation',local_recommendation(converged,required_width_m));
if converged
    validation.width_160_sufficient = local_width_passes(160,cases,metric_rows,cfg);
end
save(fullfile(cfg.output_dir,'reflected_chain_window_validation.mat'),'validation','-v7.3');
generate_reflected_chain_window_validation_figures(validation);
local_write_report(validation);
if ~converged && cfg.fail_on_nonconvergence
    error('Reflected chain did not converge within the tested window limit.');
end
end

function cfg = local_config(root,overrides)
cfg = struct('output_dir',fullfile(root,'results','validation','pe_reflected_chain_window'), ...
    'report_path',fullfile(root,'reports','pe_reflected_chain_window_validation_report.md'), ...
    'frequency_hz',4000,'c0_mps',1500,'z_tx_m',100,'z_rx_m',3, ...
    'sigma_src_m',0.3,'stepz_lamb',0.5,'dx_m',50/256, ...
    'wind_speed_mps',5,'target_hs_m',0.5,'sea_seed',12345, ...
    'reflect_coeff',-1,'master_window_nominal_m',384, ...
    'initial_windows_m',[50 128 160 192], ...
    'escalation_windows_m',[256 320 384], ...
    'center_radii_m',[0.5 1 2],'max_z_samples',65, ...
    'tl_limit_db',0.1,'phase_limit_rad',0.02,'center_l2_limit',0.01, ...
    'boundary_limit_db',-40,'edge5_fraction_limit',1e-4, ...
    'fail_on_nonconvergence',true);
names = fieldnames(overrides);
for ii=1:numel(names), cfg.(names{ii})=overrides.(names{ii}); end
if isfield(overrides,'output_dir') && ~isfield(overrides,'report_path')
    cfg.report_path=fullfile(cfg.output_dir,'reduced_validation_report.md');
end
cfg.initial_windows_m = unique(cfg.initial_windows_m,'stable');
cfg.escalation_windows_m = setdiff(cfg.escalation_windows_m,cfg.initial_windows_m,'stable');
end

function [eta,meta] = local_master_surface(cfg)
n = local_n(cfg.master_window_nominal_m,cfg.dx_m);
w = n*cfg.dx_m;
spec = raw_pm_spectrum_grid_vertical(cfg.wind_speed_mps,n,n,w,w);
[eta_raw,sample_meta] = sample_raw_pm_surface_vertical(spec,cfg.sea_seed);
raw_hs = 4*std(eta_raw(:));
if raw_hs>0, eta=eta_raw*(cfg.target_hs_m/raw_hs); else, eta=zeros(size(eta_raw)); end
meta=struct('nominal_width_m',cfg.master_window_nominal_m,'actual_width_m',w, ...
    'nx',n,'ny',n,'dx_m',cfg.dx_m,'seed',cfg.sea_seed, ...
    'wind_speed_mps',cfg.wind_speed_mps,'target_hs_m',cfg.target_hs_m, ...
    'raw_hs_m',raw_hs,'scaled_hs_m',4*std(eta(:)), ...
    'mean_m',mean(eta(:)),'min_m',min(eta(:)),'max_m',max(eta(:)), ...
    'sum_m',sum(eta(:)),'sum_squares_m2',sum(eta(:).^2), ...
    'generation_meta',sample_meta,'normalization_count',1, ...
    'crop_rule','exact central crop on common dx grid; crops are never renormalized');
end

function p=local_params(cfg,w,n,eta,ratio,alpha)
p=struct('f0',cfg.frequency_hz,'c0',cfg.c0_mps,'z_max',cfg.z_tx_m, ...
    'z_tx',cfg.z_tx_m,'z_rx',cfg.z_rx_m,'x_tx',0,'y_tx',0,'x_rx',0,'y_rx',0, ...
    'xw',w,'yw',w,'nx',n,'ny',n,'source_mode','gaussian', ...
    'sigma_src_m',cfg.sigma_src_m,'stepz_lamb',cfg.stepz_lamb, ...
    'env_mode','uniform','show_figures',false,'save_mode','rx_only', ...
    'enforce_1_over_R',false,'enable_surface_reflection',true, ...
    'surface_boundary_model','kirchhoff_spatial','surface_reflect_coeff',cfg.reflect_coeff, ...
    'surface_phase_mode','normal','sea_wind_speed',cfg.wind_speed_mps, ...
    'sea_hs_target',cfg.target_hs_m,'surface_roughness_scale_mode','target_hs', ...
    'sea_seed',cfg.sea_seed,'surface_elevation_override_xy',eta, ...
    'surface_wavefield_diagnostics',true,'surface_wavefield_max_z_samples',cfg.max_z_samples, ...
    'surface_wavefield_slice_axis','x','validation_allow_extended_window',w>100, ...
    'validation_allow_extended_sponge_ratio',true,'sponge_ratio',ratio, ...
    'alpha_max_np_per_m',alpha,'enable_bubbles',false);
end

function c=local_extract_case(out,wn,eta,cfg)
m=out.surface_wavefield_meta;
phase_screen_error=norm(abs(m.surface_reflected_xy(:))- ...
    abs(cfg.reflect_coeff)*abs(m.surface_incident_xy(:)))/max(norm(abs(m.surface_incident_xy(:))),eps);
trace_boundary=max([m.incident_energy_trace.boundary_max_relative_db, ...
    m.reflected_energy_trace.boundary_max_relative_db]);
trace_edge=max([m.incident_energy_trace.edge5_fraction, ...
    m.reflected_energy_trace.edge5_fraction]);
c=local_empty_case();
c.width_nominal_m=wn; c.width_actual_m=out.config.xw; c.nx=out.config.nx;
c.dx_m=out.config.dx; c.crop_hs_m=4*std(eta(:)); c.crop_mean_m=mean(eta(:));
c.H_reflect_physical=out.H_reflect_physical_f(out.idx_f_ref);
c.H_reflect_reduced=out.H_reflect_reduced_f(out.idx_f_ref);
c.H_direct=out.H_direct_f(out.idx_f_ref); c.H_reflect=out.H_reflect_f(out.idx_f_ref);
c.H_total=out.H_f(out.idx_f_ref); c.H_total_physical=out.H_physical_f(out.idx_f_ref);
c.closure_error=abs(c.H_total-c.H_direct-c.H_reflect);
c.phase_screen_magnitude_rel_l2=phase_screen_error;
c.endpoint_error=max(struct2array(m.endpoint_consistency));
c.reference_max_boundary_db=trace_boundary; c.reference_max_edge5_fraction=trace_edge;
c.reference_edge_valid=trace_boundary<=cfg.boundary_limit_db && trace_edge<=cfg.edge5_fraction_limit;
c.x_m=m.x_m; c.y_m=m.y_m; c.surface_elevation_xy=eta;
c.surface_incident_xy=m.surface_incident_xy; c.surface_reflected_xy=m.surface_reflected_xy;
c.receiver_reflected_xy=m.receiver_reflected_xy;
c.incident_energy_trace=m.incident_energy_trace; c.reflected_energy_trace=m.reflected_energy_trace;
end

function row=local_compare_pair(a,b,cfg)
row=local_empty_pair(); row.smaller_width_m=a.width_actual_m; row.larger_width_m=b.width_actual_m;
row.reflect_tl_difference_db=-20*log10(abs(a.H_reflect_physical/b.H_reflect_physical));
row.reflect_phase_difference_rad=angle(a.H_reflect_physical*conj(b.H_reflect_physical));
for name={'surface_incident_xy','surface_reflected_xy','receiver_reflected_xy'}
    key=name{1}; suffix=erase(key,'_xy');
    row.([suffix '_center_l2'])=local_center_l2(a.(key),a.x_m,a.y_m,b.(key),b.x_m,b.y_m,2);
end
row.larger_boundary_max_db=b.reference_max_boundary_db;
row.larger_edge5_fraction=b.reference_max_edge5_fraction;
row.passed=local_pair_passes(row,cfg) && b.reference_edge_valid;
end

function tf=local_pair_passes(row,cfg)
tf=abs(row.reflect_tl_difference_db)<=cfg.tl_limit_db && ...
    abs(row.reflect_phase_difference_rad)<=cfg.phase_limit_rad && ...
    row.surface_incident_center_l2<=cfg.center_l2_limit && ...
    row.surface_reflected_center_l2<=cfg.center_l2_limit && ...
    row.receiver_reflected_center_l2<=cfg.center_l2_limit;
end

function rows=local_reference_metrics(cases,iref,cfg)
rows=repmat(local_empty_metric(),numel(cases)*numel(cfg.center_radii_m),1); q=0;
ref=cases(iref);
for ii=1:numel(cases)
    for rr=1:numel(cfg.center_radii_m)
        q=q+1; r=cfg.center_radii_m(rr); cur=cases(ii);
        rows(q).width_nominal_m=cur.width_nominal_m; rows(q).width_actual_m=cur.width_actual_m;
        rows(q).radius_m=r;
        rows(q).reflect_tl_error_db=-20*log10(abs(cur.H_reflect_physical/ref.H_reflect_physical));
        rows(q).reflect_phase_error_rad=angle(cur.H_reflect_physical*conj(ref.H_reflect_physical));
        for name={'surface_incident_xy','surface_reflected_xy','receiver_reflected_xy'}
            key=name{1}; suffix=erase(key,'_xy');
            [l2,mean_amp,max_amp,phase_rms,energy_db]=local_field_metrics( ...
                cur.(key),cur.x_m,cur.y_m,ref.(key),ref.x_m,ref.y_m,r);
            rows(q).([suffix '_complex_l2'])=l2;
            rows(q).([suffix '_mean_amplitude_error_db'])=mean_amp;
            rows(q).([suffix '_max_amplitude_error_db'])=max_amp;
            rows(q).([suffix '_phase_rms_rad'])=phase_rms;
            rows(q).([suffix '_energy_error_db'])=energy_db;
        end
    end
end
end

function [idx,w]=local_required_window(cases,rows,cfg)
idx=NaN; w=NaN;
for ii=1:numel(cases)
    mask=[rows.width_actual_m]==cases(ii).width_actual_m & [rows.radius_m]==2;
    r=rows(mask);
    if isempty(r), continue, end
    field_ok=r.surface_incident_complex_l2<=cfg.center_l2_limit && ...
        r.surface_reflected_complex_l2<=cfg.center_l2_limit && ...
        r.receiver_reflected_complex_l2<=cfg.center_l2_limit;
    scalar_ok=abs(r.reflect_tl_error_db)<=cfg.tl_limit_db && abs(r.reflect_phase_error_rad)<=cfg.phase_limit_rad;
    later_ok=true;
    for jj=ii:numel(cases)
        mask_j=[rows.width_actual_m]==cases(jj).width_actual_m & [rows.radius_m]==2;
        rr=rows(mask_j);
        later_ok=later_ok && ~isempty(rr) && abs(rr.reflect_tl_error_db)<=cfg.tl_limit_db && ...
            abs(rr.reflect_phase_error_rad)<=cfg.phase_limit_rad && ...
            rr.surface_incident_complex_l2<=cfg.center_l2_limit && ...
            rr.surface_reflected_complex_l2<=cfg.center_l2_limit && ...
            rr.receiver_reflected_complex_l2<=cfg.center_l2_limit;
    end
    if scalar_ok && field_ok && later_ok && cases(ii).reference_edge_valid
        idx=ii; w=cases(ii).width_actual_m; return
    end
end
end

function tf=local_width_passes(nominal,cases,rows,cfg)
[~,ii]=min(abs([cases.width_nominal_m]-nominal));
if isempty(ii) || abs(cases(ii).width_nominal_m-nominal)>1, tf=false; return, end
r=rows([rows.width_actual_m]==cases(ii).width_actual_m & [rows.radius_m]==2);
tf=~isempty(r) && abs(r.reflect_tl_error_db)<=cfg.tl_limit_db && ...
    abs(r.reflect_phase_error_rad)<=cfg.phase_limit_rad && ...
    r.surface_incident_complex_l2<=cfg.center_l2_limit && ...
    r.surface_reflected_complex_l2<=cfg.center_l2_limit && ...
    r.receiver_reflected_complex_l2<=cfg.center_l2_limit && cases(ii).reference_edge_valid;
end

function e=local_center_l2(a,xa,ya,b,xb,yb,r)
[va,vb]=local_center_vectors(a,xa,ya,b,xb,yb,r);
e=norm(va-vb)/max(norm(vb),eps);
end

function [l2,mean_amp,max_amp,phase_rms,energy_db]=local_field_metrics(a,xa,ya,b,xb,yb,r)
[va,vb]=local_center_vectors(a,xa,ya,b,xb,yb,r);
l2=norm(va-vb)/max(norm(vb),eps);
amp_db=20*log10(max(abs(va),eps)./max(abs(vb),eps));
mean_amp=mean(abs(amp_db)); max_amp=max(abs(amp_db));
phase_rms=sqrt(mean(angle(va.*conj(vb)).^2));
energy_db=10*log10(sum(abs(va).^2)/max(sum(abs(vb).^2),eps));
end

function [va,vb]=local_center_vectors(a,xa,ya,b,xb,yb,r)
[Xa,Ya]=meshgrid(xa,ya); [Xb,Yb]=meshgrid(xb,yb);
va=a(hypot(Xa,Ya)<=r); vb=b(hypot(Xb,Yb)<=r);
if numel(va)~=numel(vb), error('Common-dx center regions are not aligned.'); end
end

function t=local_energy_table(cases)
rows=repmat(struct('width_actual_m',NaN,'stage',"",'z_m',NaN, ...
    'total_energy',NaN,'center_energy_rho_le_2m',NaN,'edge5_energy',NaN, ...
    'edge10_energy',NaN,'edge5_fraction',NaN,'edge10_fraction',NaN, ...
    'boundary_max_relative_db',NaN),0,1);
for ii=1:numel(cases)
    for stage={'upward_incident','downward_reflected'}
        if strcmp(stage{1},'upward_incident'), tr=cases(ii).incident_energy_trace; else, tr=cases(ii).reflected_energy_trace; end
        for jj=1:numel(tr.z_m)
            q=numel(rows)+1; rows(q).width_actual_m=cases(ii).width_actual_m; rows(q).stage=string(stage{1}); rows(q).z_m=tr.z_m(jj);
            for f={'total_energy','center_energy_rho_le_2m','edge5_energy','edge10_energy','edge5_fraction','edge10_fraction','boundary_max_relative_db'}
                rows(q).(f{1})=tr.(f{1})(jj);
            end
        end
    end
end
t=struct2table(rows);
end

function rows=local_case_rows(cases)
rows=repmat(struct('width_nominal_m',NaN,'width_actual_m',NaN,'nx',NaN,'dx_m',NaN, ...
    'crop_hs_m',NaN,'crop_mean_m',NaN,'reflect_magnitude',NaN,'reflect_tl_db',NaN, ...
    'reflect_phase_rad',NaN,'closure_error',NaN,'phase_screen_magnitude_rel_l2',NaN, ...
    'endpoint_error',NaN,'max_boundary_db',NaN,'max_edge5_fraction',NaN,'edge_valid',false),numel(cases),1);
for ii=1:numel(cases)
    c=cases(ii); rows(ii).width_nominal_m=c.width_nominal_m; rows(ii).width_actual_m=c.width_actual_m;
    rows(ii).nx=c.nx; rows(ii).dx_m=c.dx_m; rows(ii).crop_hs_m=c.crop_hs_m; rows(ii).crop_mean_m=c.crop_mean_m;
    rows(ii).reflect_magnitude=abs(c.H_reflect_physical); rows(ii).reflect_tl_db=-20*log10(abs(c.H_reflect_physical));
    rows(ii).reflect_phase_rad=angle(c.H_reflect_physical); rows(ii).closure_error=c.closure_error;
    rows(ii).phase_screen_magnitude_rel_l2=c.phase_screen_magnitude_rel_l2; rows(ii).endpoint_error=c.endpoint_error;
    rows(ii).max_boundary_db=c.reference_max_boundary_db; rows(ii).max_edge5_fraction=c.reference_max_edge5_fraction;
    rows(ii).edge_valid=c.reference_edge_valid;
end
end

function r=local_recommendation(ok,w)
if ok
    r=struct('window_m',w,'use_sponge',false,'sponge_ratio',0,'alpha_max_np_per_m',0, ...
        'reason','Smallest tested no-sponge window satisfying receiver, center-field, and edge gates.');
else
    r=struct('window_m',NaN,'use_sponge',false,'sponge_ratio',NaN,'alpha_max_np_per_m',NaN, ...
        'reason','No converged no-sponge reference; production change is blocked.');
end
end

function local_write_report(v)
path=v.config.report_path;
fid=fopen(path,'w'); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'# Full PE -> surface -> PE window convergence validation\n\n');
fprintf(fid,'Status: **%s**. Production defaults were not changed.\n\n',local_word(v.converged,'CONVERGED','NOT CONVERGED'));
fprintf(fid,'## Environment\n\n- Uniform c=1500 m/s; Tx=(0,0,100 m); Rx=(0,0,3 m).\n');
fprintf(fid,'- 4 kHz, production Gaussian sigma=0.3 m, Kirchhoff spatial, U=5 m/s, master Hs=0.5 m.\n');
fprintf(fid,'- One seed-12345 master surface, normalized once and centrally cropped; all single-frequency cases use sponge off.\n\n');
fprintf(fid,'## Single-frequency results\n\n');
fprintf(fid,'| nominal W (m) | actual W (m) | crop Hs (m) | reflect TL (dB) | phase (rad) | max edge-5 fraction | boundary (dB) | edge gate |\n');
fprintf(fid,'|---:|---:|---:|---:|---:|---:|---:|:---:|\n');
for ii=1:height(v.case_table)
    row=v.case_table(ii,:);
    fprintf(fid,'| %.3f | %.6f | %.6g | %.6g | %.6g | %.6g | %.6g | %s |\n', ...
        row.width_nominal_m,row.width_actual_m,row.crop_hs_m,row.reflect_tl_db, ...
        row.reflect_phase_rad,row.max_edge5_fraction,row.max_boundary_db, ...
        local_word(row.edge_valid,'PASS','FAIL'));
end
fprintf(fid,'\n| window pair (m) | reflect TL difference (dB) | phase difference (rad) | incident center L2 | receiver center L2 | pass |\n');
fprintf(fid,'|---|---:|---:|---:|---:|:---:|\n');
for ii=1:height(v.pair_table)
    row=v.pair_table(ii,:);
    fprintf(fid,'| %.3f -> %.3f | %.6g | %.6g | %.6g | %.6g | %s |\n', ...
        row.smaller_width_m,row.larger_width_m,row.reflect_tl_difference_db, ...
        row.reflect_phase_difference_rad,row.surface_incident_center_l2, ...
        row.receiver_reflected_center_l2,local_word(row.passed,'PASS','FAIL'));
end
fprintf(fid,'## Direct answers\n\n');
if v.converged
    fprintf(fid,'1. Required reflected-path window: **%.6g m actual**.\n',v.required_width_m);
    fprintf(fid,'2. Is nominal 160 m sufficient under every preregistered gate: **%s**. Its receiver response may still satisfy the response-only tolerances; inspect the edge columns above.\n',local_word(v.width_160_sufficient,'yes','no'));
    fprintf(fid,'3. Sponge decision: **larger window with sponge off** under the preregistered gates.\n');
else
    fprintf(fid,'1. Required reflected-path window: **not established** within tested limits.\n');
    fprintf(fid,'2. Is nominal 160 m sufficient: **not demonstrated**.\n');
    fprintf(fid,'3. Sponge decision: **deferred**; strong sponge is not used to mask nonconvergence.\n');
end
fprintf(fid,'\n## Artifacts\n\nCSV/MAT and figures are in `results/validation/pe_reflected_chain_window/`.\n');
end

function s=local_word(tf,a,b), if tf, s=a; else, s=b; end, end
function n=local_n(w,dx), n=2*round(w/dx/2); end
function crop=local_center_crop(a,n), s=(size(a,1)-n)/2+1; crop=a(s:s+n-1,s:s+n-1); end
function c=local_empty_case()
c=struct('width_nominal_m',NaN,'width_actual_m',NaN,'nx',NaN,'dx_m',NaN, ...
    'crop_hs_m',NaN,'crop_mean_m',NaN,'H_reflect_physical',complex(NaN), ...
    'H_reflect_reduced',complex(NaN),'H_direct',complex(NaN),'H_reflect',complex(NaN), ...
    'H_total',complex(NaN),'H_total_physical',complex(NaN),'closure_error',NaN, ...
    'phase_screen_magnitude_rel_l2',NaN,'endpoint_error',NaN, ...
    'reference_max_boundary_db',NaN,'reference_max_edge5_fraction',NaN, ...
    'reference_edge_valid',false,'x_m',[],'y_m',[],'surface_elevation_xy',[], ...
    'surface_incident_xy',[],'surface_reflected_xy',[],'receiver_reflected_xy',[], ...
    'incident_energy_trace',struct(),'reflected_energy_trace',struct());
end
function r=local_empty_pair()
r=struct('smaller_width_m',NaN,'larger_width_m',NaN,'reflect_tl_difference_db',NaN, ...
    'reflect_phase_difference_rad',NaN,'surface_incident_center_l2',NaN, ...
    'surface_reflected_center_l2',NaN,'receiver_reflected_center_l2',NaN, ...
    'larger_boundary_max_db',NaN,'larger_edge5_fraction',NaN,'passed',false);
end
function r=local_empty_metric()
r=struct('width_nominal_m',NaN,'width_actual_m',NaN,'radius_m',NaN, ...
    'reflect_tl_error_db',NaN,'reflect_phase_error_rad',NaN, ...
    'surface_incident_complex_l2',NaN,'surface_incident_mean_amplitude_error_db',NaN, ...
    'surface_incident_max_amplitude_error_db',NaN,'surface_incident_phase_rms_rad',NaN, ...
    'surface_incident_energy_error_db',NaN,'surface_reflected_complex_l2',NaN, ...
    'surface_reflected_mean_amplitude_error_db',NaN,'surface_reflected_max_amplitude_error_db',NaN, ...
    'surface_reflected_phase_rms_rad',NaN,'surface_reflected_energy_error_db',NaN, ...
    'receiver_reflected_complex_l2',NaN,'receiver_reflected_mean_amplitude_error_db',NaN, ...
    'receiver_reflected_max_amplitude_error_db',NaN,'receiver_reflected_phase_rms_rad',NaN, ...
    'receiver_reflected_energy_error_db',NaN);
end
