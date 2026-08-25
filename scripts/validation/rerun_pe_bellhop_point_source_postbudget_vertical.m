function validation=rerun_pe_bellhop_point_source_postbudget_vertical(overrides)
%RERUN_PE_BELLHOP_POINT_SOURCE_POSTBUDGET_VERTICAL Final gated comparison.
% Runs only after the continuous Weyl reference gate is available.

if nargin<1 || isempty(overrides), overrides=struct(); end
root=fileparts(fileparts(fileparts(mfilename('fullpath')))); addpath(root); setup_vertical_project();
cfg=local_defaults(root); cfg=local_overrides(cfg,overrides);
budget_file=fullfile(cfg.error_budget_dir,'point_source_error_budget.mat');
if exist(budget_file,'file')~=2, error('Run validate_pe_point_source_error_budget_vertical first.'); end
d=load(budget_file,'validation'); budget=d.validation;
if ~budget.continuous_weyl_reliable
    error('Continuous Weyl free-space reference gate failed; Bellhop rerun is blocked.');
end
if ~exist(cfg.output_dir,'dir'), mkdir(cfg.output_dir); end

bh_audit=validate_bellhop_freefield_normalization_vertical(struct( ...
    'output_dir',fullfile(cfg.output_dir,'bellhop_normalization'), ...
    'bellhop_exe',cfg.bellhop_exe,'frequencies_hz',cfg.frequency_hz, ...
    'beam_counts',cfg.normalization_beam_counts, ...
    'step_values_m',cfg.normalization_step_values_m));

Raxial=cfg.s0_m+cfg.L_m;
bh_cfg=struct('bellhop_exe',cfg.bellhop_exe, ...
    'case_root',fullfile(cfg.output_dir,'bellhop_rerun'), ...
    'frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps,'source_depth_m',0, ...
    'receiver_depths_m',cfg.offsets_m,'receiver_ranges_m',Raxial, ...
    'run_type','C','beam_count',cfg.bellhop_beam_count,'angle_limits_deg',[-180 180], ...
    'step_m',cfg.bellhop_step_m,'domain_half_depth_m',1000);
bh=run_bellhop_freefield_vertical(bh_cfg); p_bh=bh.data.pressure(:,1);
if bh_audit.selected_spatial_sign<0
    H_bh=conj(p_bh)./conj(bh_audit.bellhop_source_constant)/(4*pi);
else
    H_bh=p_bh./bh_audit.bellhop_source_constant/(4*pi);
end

alpha_values=[0 cfg.sponge_alpha_np_per_m]; H_pe=complex(zeros(numel(cfg.offsets_m),2));
for aa=1:2
    p=local_pe_params(cfg,alpha_values(aa)); out=vertical_channel_model(p);
    [~,iy]=min(abs(out.y)); k=2*pi*cfg.frequency_hz/cfg.c0_mps;
    for xx=1:numel(cfg.offsets_m)
        [~,ix]=min(abs(out.x-cfg.offsets_m(xx)));
        H_pe(xx,aa)=out.psifinal_xy(iy,ix)*exp(1i*k*Raxial);
    end
end
k=2*pi*cfg.frequency_hz/cfg.c0_mps; ranges=hypot(Raxial,cfg.offsets_m(:));
H_exact=exp(1i*k*ranges)./(4*pi*ranges);
rows=table(cfg.offsets_m(:),ranges,H_pe(:,1),H_pe(:,2),H_bh,H_exact, ...
    'VariableNames',{'offset_m','range_m','H_pe_no_sponge','H_pe_sponge','H_bellhop','H_exact'});
rows.finite_window_amplitude_error_db=20*log10(abs(rows.H_pe_no_sponge./rows.H_exact));
rows.finite_window_tl_error_db=-rows.finite_window_amplitude_error_db;
rows.sponge_amplitude_change_db=20*log10(abs(rows.H_pe_sponge./rows.H_pe_no_sponge));
rows.sponge_tl_change_db=-rows.sponge_amplitude_change_db;
rows.total_pe_bellhop_amplitude_error_db=20*log10(abs(rows.H_pe_sponge./rows.H_bellhop));
rows.total_pe_bellhop_tl_error_db=-rows.total_pe_bellhop_amplitude_error_db;
rows.no_sponge_phase_error_rad=angle(rows.H_pe_no_sponge.*conj(rows.H_exact));
rows.sponge_phase_increment_rad=angle(rows.H_pe_sponge.*conj(rows.H_pe_no_sponge));
rows.total_pe_bellhop_phase_rad=angle(rows.H_pe_sponge.*conj(rows.H_bellhop));
rows.bellhop_analytic_amplitude_error_db=20*log10(abs(rows.H_bellhop./rows.H_exact));
rows.bellhop_analytic_tl_error_db=-rows.bellhop_analytic_amplitude_error_db;

energy_rows=budget.sponge_table(budget.sponge_table.width_m==cfg.width_m & ...
    abs(budget.sponge_table.sponge_ratio-cfg.sponge_ratio)<1e-12 & ...
    ismember(budget.sponge_table.alpha_max_np_per_m,[0 cfg.sponge_alpha_np_per_m]),:);

validation=struct('schema_version','1.0.0','config',cfg,'error_budget',budget, ...
    'bellhop_normalization',bh_audit,'comparison_table',rows, ...
    'default_sponge_energy_table',energy_rows, ...
    'passed',max(abs(rows.total_pe_bellhop_tl_error_db))<=cfg.tl_tolerance_db && ...
    sqrt(mean(rows.total_pe_bellhop_phase_rad.^2))<=cfg.phase_tolerance_rad);
validation.files=local_outputs(validation,root);
end

function cfg=local_defaults(root)
cfg=struct('bellhop_exe',getenv('BELLHOP_EXE'), ...
    'error_budget_dir',fullfile(root,'results','validation','pe_point_source_error_budget'), ...
    'output_dir',fullfile(root,'results','validation','pe_bellhop_freefield','post_error_budget'), ...
    'frequency_hz',4000,'c0_mps',1500,'s0_m',10,'L_m',70,'offsets_m',[0 0.5 1 2], ...
    'width_m',100,'nx',400,'stepz_lamb',0.5,'sponge_ratio',0.12, ...
    'sponge_alpha_np_per_m',0.15,'bellhop_beam_count',10001,'bellhop_step_m',0.05, ...
    'normalization_beam_counts',[5001 10001],'normalization_step_values_m',[0.1 0.05], ...
    'tl_tolerance_db',0.5,'phase_tolerance_rad',0.1);
end

function cfg=local_overrides(cfg,o)
names=fieldnames(o); for ii=1:numel(names), if ~isfield(cfg,names{ii}), error('Unknown override: %s',names{ii}); end, cfg.(names{ii})=o.(names{ii}); end
end

function p=local_pe_params(cfg,alpha)
p=struct('f0',cfg.frequency_hz,'c0',cfg.c0_mps,'z_max',cfg.L_m,'z_tx',cfg.L_m, ...
    'z_rx',0,'xw',cfg.width_m,'yw',cfg.width_m,'nx',cfg.nx,'ny',cfg.nx, ...
    'x_tx',0,'y_tx',0,'x_rx',0,'y_rx',0,'stepz_lamb',cfg.stepz_lamb, ...
    'source_mode','custom_field_fn','source_field_fn',@virtual_point_source_initial_field_vertical, ...
    'virtual_source_distance_m',cfg.s0_m,'source_green_amplitude',1/(4*pi), ...
    'sponge_ratio',cfg.sponge_ratio,'alpha_max_np_per_m',alpha,'env_mode','uniform', ...
    'enable_surface_reflection',false,'enable_bubbles',false,'enforce_1_over_R',false, ...
    'show_figures',false,'save_mode','slice','use_gpu',false);
end

function files=local_outputs(v,root)
out=v.config.output_dir; csv=fullfile(out,'postbudget_pe_bellhop.csv');
mat=fullfile(out,'postbudget_pe_bellhop.mat'); report=fullfile(root,'reports','pe_bellhop_post_error_budget_report.md');
writetable(v.comparison_table,csv); validation=v; schema_version=v.schema_version;
save(mat,'validation','schema_version','-v7.3');
fid=fopen(report,'w','n','UTF-8'); c=onCleanup(@()fclose(fid));
fprintf(fid,'# Post-error-budget PE--Bellhop rerun\n\nStatus: `passed=%s`.\n\n',string(v.passed));
fprintf(fid,'The continuous Weyl reference passed independently. The PE input remains the finite spatial-plane source so its error can be separated from sponge error.\n\n');
fprintf(fid,'Definitions: `dA = 20 log10(|H_sponge|/|H_no_sponge|)` and `dTL = TL_sponge - TL_no_sponge = -dA`. Positive dTL means additional loss.\n\n');
fprintf(fid,'- Maximum finite-window TL error: `%.6g dB`.\n',max(abs(v.comparison_table.finite_window_tl_error_db)));
fprintf(fid,'- Maximum sponge-only TL change: `%.6g dB`.\n',max(abs(v.comparison_table.sponge_tl_change_db)));
fprintf(fid,'- Maximum final PE--Bellhop TL error: `%.6g dB`.\n',max(abs(v.comparison_table.total_pe_bellhop_tl_error_db)));
fprintf(fid,'- Bellhop--analytic maximum TL error: `%.6g dB`.\n',max(abs(v.comparison_table.bellhop_analytic_tl_error_db)));
fprintf(fid,'\n## Receiver-by-receiver decomposition\n\n');
fprintf(fid,'| Offset (m) | Finite-window dA (dB) | Finite-window dTL (dB) | Sponge dA (dB) | Sponge dTL (dB) | Sponge dphase (rad) | Final PE-BH dTL (dB) |\n');
fprintf(fid,'|---:|---:|---:|---:|---:|---:|---:|\n');
for ii=1:height(v.comparison_table)
    r=v.comparison_table(ii,:);
    fprintf(fid,'| %.3f | %.6f | %.6f | %.6f | %.6f | %.6f | %.6f |\n', ...
        r.offset_m,r.finite_window_amplitude_error_db,r.finite_window_tl_error_db, ...
        r.sponge_amplitude_change_db,r.sponge_tl_change_db, ...
        r.sponge_phase_increment_rad,r.total_pe_bellhop_tl_error_db);
end
fprintf(fid,'\n## Default sponge terminal-plane energy\n\n');
fprintf(fid,'Center energy is integrated over `rho <= %.3g m`; edge energy is integrated over the ratio-defined sponge band.\n\n',v.error_budget.config.center_energy_radius_m);
fprintf(fid,'| alpha (Np/m) | |H axis| | phase axis (rad) | center energy | edge energy | total energy | dE center (dB) | dE edge (dB) |\n');
fprintf(fid,'|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for ii=1:height(v.default_sponge_energy_table)
    r=v.default_sponge_energy_table(ii,:);
    fprintf(fid,'| %.3f | %.9g | %.9f | %.9g | %.9g | %.9g | %.6f | %.6f |\n', ...
        r.alpha_max_np_per_m,abs(r.H_axis),angle(r.H_axis),r.center_energy, ...
        r.edge_energy,r.total_energy,r.center_energy_change_db,r.edge_energy_change_db);
end
clear c
files=struct('mat',mat,'csv',csv,'report',report);
end
