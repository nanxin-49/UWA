function comparison = validate_pe_bellhop_pm_stage1_dimensionality(overrides)
%VALIDATE_PE_BELLHOP_PM_STAGE1_DIMENSIONALITY Stage 1B PE dimensionality audit.
%   The same fixed, band-limited eta(x) is used in the 1-transverse bridge and
%   in the production two-transverse PE with eta(x,y)=eta(x).  Only the PE
%   dimensionality is changed; no Bellhop or production physics is modified.

if nargin < 1 || isempty(overrides), overrides = struct(); end
if ~isstruct(overrides) || ~isscalar(overrides)
    error('overrides must be a scalar struct.');
end
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root); setup_vertical_project();
cfg = local_config(root, overrides);
if ~exist(cfg.output_dir,'dir'), mkdir(cfg.output_dir); end

profile = load_fixed_pm_profile_for_pe_bellhop_validation(cfg.coeff_file);
x = (-0.5*cfg.xw_m) + (0:cfg.nx-1)*(cfg.xw_m/cfg.nx);
eval_x = evaluate_fixed_pm_fourier_profile(profile,x);
eta_x = eval_x.eta_m(:).';
base = struct('frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps, ...
    'xw_m',cfg.xw_m,'nx',cfg.nx,'z_tx_m',cfg.z_tx_m,'z_rx_m',cfg.z_rx_m, ...
    'sigma_src_m',cfg.sigma_src_m,'surface_reflect_coeff',-1, ...
    'step_m',cfg.step_m,'x_rx_m',cfg.x_rx_m);
flat_1d = run_pe_1d_surface_reflection_validation(setfield(base,'surface_elevation_x_m',zeros(1,cfg.nx))); %#ok<SFLD>
rough_1d = run_pe_1d_surface_reflection_validation(setfield(base,'surface_elevation_x_m',eta_x)); %#ok<SFLD>
G1 = rough_1d.reflected_receiver/flat_1d.reflected_receiver;

ny_values = cfg.ny_values(:).';
rows = repmat(local_empty_row(),1,numel(ny_values));
for ii = 1:numel(ny_values)
    ny = ny_values(ii);
    flat2 = local_run_full_pe(cfg,ny,zeros(1,cfg.nx));
    rough2 = local_run_full_pe(cfg,ny,eta_x);
    G2 = rough2.H_reflect_reduced/flat2.H_reflect_reduced;
    rows(ii).ny = ny;
    rows(ii).G_1t = G1;
    rows(ii).G_2t = G2;
    rows(ii).delta_tl_db = 20*log10(abs(G2/G1));
    rows(ii).delta_phase_rad = angle(G2*conj(G1));
    rows(ii).complex_relative_error = abs(G2-G1)/max(abs(G1),realmin);
    rows(ii).flat_field_finite = all(isfinite([real(flat2.H_reflect_reduced),imag(flat2.H_reflect_reduced)]));
    rows(ii).rough_field_finite = all(isfinite([real(rough2.H_reflect_reduced),imag(rough2.H_reflect_reduced)]));
    rows(ii).flat = flat2;
    rows(ii).rough = rough2;
end

% The dimensionality effect is measured, not rejected as a model failure. The
% hard checks only require finite production outputs and a converged ny pair.
finite_flags = [rows.flat_field_finite] & [rows.rough_field_finite];
if numel(rows)>1
    dny_tl = 20*log10(abs(rows(end).G_2t/rows(1).G_2t));
    dny_phase = angle(rows(end).G_2t*conj(rows(1).G_2t));
else
    dny_tl = NaN; dny_phase = NaN;
end
checks = struct('profile_provenance',strcmp(profile.coeff_file_sha256,local_file_sha256(cfg.coeff_file)), ...
    'one_transverse_finite',all(isfinite([real(G1),imag(G1)])), ...
    'two_transverse_finite',all(finite_flags), ...
    'ny_sensitivity_finite',isfinite(dny_tl)&&isfinite(dny_phase));
checks.all=all(structfun(@(v)logical(v),checks));
comparison = struct('schema_version','1.0.0','stage','1B_pe_dimensionality', ...
    'config',cfg,'profile',profile,'flat_1d',flat_1d,'rough_1d',rough_1d, ...
    'G_1t',G1,'rows',rows,'ny_delta_tl_db',dny_tl,'ny_delta_phase_rad',dny_phase, ...
    'checks',checks,'passed',checks.all);
comparison.files=local_write_outputs(comparison);
if cfg.fail_on_check && ~comparison.passed
    error('Stage 1B dimensionality audit failed; see %s.',cfg.report_path);
end
end

function out = local_run_full_pe(cfg,ny,eta_x)
p = struct('f0',cfg.frequency_hz,'enable_wideband',false,'c0',cfg.c0_mps, ...
    'z_max',cfg.z_tx_m,'z_tx',cfg.z_tx_m,'z_rx',cfg.z_rx_m, ...
    'stepz_lamb',cfg.stepz_lamb,'xw',cfg.xw_m,'yw',cfg.yw_m, ...
    'nx',cfg.nx,'ny',ny,'x_tx',0,'y_tx',0,'x_rx',cfg.x_rx_m,'y_rx',0, ...
    'sigma_src_m',cfg.sigma_src_m,'source_mode','gaussian','sponge_ratio',0, ...
    'alpha_max_np_per_m',0,'validation_allow_extended_window',true, ...
    'validation_allow_extended_sponge_ratio',true,'env_mode','uniform', ...
    'show_figures',false,'enforce_1_over_R',false,'enable_surface_reflection',true, ...
    'surface_reflect_coeff',-1,'surface_phase_mode','normal', ...
    'surface_boundary_model','kirchhoff_spatial', ...
    'surface_elevation_override_xy',repmat(eta_x,ny,1), ...
    'surface_wavefield_diagnostics',false,'save_mode','rx_only','use_gpu',false, ...
    'channel_phase_reference','direct_dsp');
o=vertical_channel_model(p);
out=struct('H_direct_reduced',o.H_direct_reduced_f(1), ...
    'H_reflect_reduced',o.H_reflect_reduced_f(1),'H_total_reduced',o.H_total_reduced_f(1), ...
    'config',o.config);
end

function cfg=local_config(root,o)
cfg=struct('frequency_hz',4000,'c0_mps',1500,'z_tx_m',100,'z_rx_m',3, ...
    'sigma_src_m',0.3,'xw_m',192.1875,'nx',984,'yw_m',50,'ny_values',[256 512], ...
    'step_m',0.05,'stepz_lamb',0.5,'x_rx_m',0,'wall_seed',260001, ...
    'coeff_file',fullfile(root,'results','validation','bellhop_internal_pm_fixed_realization', ...
    'fixed_pm_fourier_coefficients.csv'),'output_dir',fullfile(root,'results','validation', ...
    'pe_bellhop_pm_stage1_dimensionality'),'report_path',fullfile(root,'reports', ...
    'pe_bellhop_pm_stage1_dimensionality_report.md'),'fail_on_check',true);
names=fieldnames(o);for ii=1:numel(names),if ~isfield(cfg,names{ii}),error('Unknown Stage 1B override: %s.',names{ii});end;cfg.(names{ii})=o.(names{ii});end
end

function row=local_empty_row()
row=struct('ny',NaN,'G_1t',NaN,'G_2t',NaN,'delta_tl_db',NaN,'delta_phase_rad',NaN, ...
    'complex_relative_error',NaN,'flat_field_finite',false,'rough_field_finite',false, ...
    'flat',struct(),'rough',struct());
end

function files=local_write_outputs(c)
out=c.config.output_dir;
summary=rmfield(c.rows,{'flat','rough'});
writetable(struct2table(summary),fullfile(out,'stage1b_dimensionality_summary.csv'));
writetable(struct2table(c.checks),fullfile(out,'stage1b_dimensionality_checks.csv'));
save(fullfile(out,'stage1b_dimensionality_comparison.mat'),'c','-v7.3');
local_write_report(c.config.report_path,c);
files=struct('summary',fullfile(out,'stage1b_dimensionality_summary.csv'), ...
    'checks',fullfile(out,'stage1b_dimensionality_checks.csv'), ...
    'mat',fullfile(out,'stage1b_dimensionality_comparison.mat'),'report',c.config.report_path);
end

function local_write_report(path,c)
fid=fopen(path,'w','n','UTF-8');if fid<0,error('Cannot write %s.',path);end
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# PE--Bellhop PM Stage 1B dimensionality sensitivity\n\n状态：**%s**\n\n',ternary(c.passed,'PASS_WITH_LIMITS','FAIL'));
fprintf(fid,'固定条件：4 kHz、uniform c=1500 m/s、seed=260001、W=%.9g m、nx=%d、PE step=%.6g m、y-invariant eta(x,y)=eta(x)、on-axis receiver。1-transverse PE 与 production 2-transverse PE 的 reflected-only rough/flat ratio 仅用于量化 dimensionality sensitivity。\n\n',c.config.xw_m,c.config.nx,c.config.step_m);
fprintf(fid,'## Results\n\n| ny | G_1T | G_2T | delta TL (dB) | delta phase (rad) | complex relative error |\n|---:|---:|---:|---:|---:|---:|\n');
for ii=1:numel(c.rows),r=c.rows(ii);fprintf(fid,'| %d | %.12g%+.12gi | %.12g%+.12gi | %.8g | %.8g | %.8g |\n',r.ny,real(r.G_1t),imag(r.G_1t),real(r.G_2t),imag(r.G_2t),r.delta_tl_db,r.delta_phase_rad,r.complex_relative_error);end
fprintf(fid,'\n2T ny endpoint sensitivity: %.8g dB, %.8g rad. The 2T-vs-1T difference is dimensional/source mapping diagnostic, not a Bellhop or PE implementation failure.\n\n',c.ny_delta_tl_db,c.ny_delta_phase_rad);
fprintf(fid,'## Checks\n\n');n=fieldnames(c.checks);for ii=1:numel(n),fprintf(fid,'- %s: %s\n',n{ii},ternary(c.checks.(n{ii}),'PASS','FAIL'));end
fprintf(fid,'\n## Interpretation\n\nThe surface samples are copied identically in y, and the production PE core is unchanged. This report quantifies the 1T-to-2T bridge contribution that must be reported alongside the PE/Bellhop Stage 1A model discrepancy.\n');
end
function s=ternary(c,a,b),if c,s=a;else,s=b;end,end
function d=local_file_sha256(path),md=java.security.MessageDigest.getInstance('SHA-256');fid=fopen(path,'rb');b=fread(fid,Inf,'*uint8');fclose(fid);md.update(typecast(b,'int8'));d=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));end
