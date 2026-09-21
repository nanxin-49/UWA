function bridge = validate_pe_1d_validation_bridge(overrides)
%VALIDATE_PE_1D_VALIDATION_BRIDGE Stage 0B 1-D PE/ky=0 audit.
%   Runs an independent one-transverse-dimensional split-step bridge and
%   compares it with a one-step exact 1-D angular spectrum and the ky=0
%   projection of the production two-transverse-dimensional PE.  No Bellhop
%   propagation is performed in this stage.

if nargin < 1 || isempty(overrides), overrides = struct(); end
if ~isstruct(overrides) || ~isscalar(overrides)
    error('overrides must be a scalar struct.');
end
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup_vertical_project();
cfg = local_config(root, overrides);
if ~exist(cfg.output_dir, 'dir'), mkdir(cfg.output_dir); end

x = (-0.5*cfg.xw_m) + (0:cfg.nx-1) * (cfg.xw_m/cfg.nx);
eta_cases = {zeros(1,cfg.nx), cfg.weak_amplitude_m*cos(2*pi*x/cfg.xw_m)};
case_names = {'flat','weak_phase_screen'};
case_rows = repmat(local_empty_row(), 1, numel(case_names));
for ic = 1:numel(case_names)
    one_cfg = struct('frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps, ...
        'xw_m',cfg.xw_m,'nx',cfg.nx,'z_tx_m',cfg.z_tx_m,'z_rx_m',cfg.z_rx_m, ...
        'sigma_src_m',cfg.sigma_src_m,'surface_elevation_x_m',eta_cases{ic}, ...
        'surface_reflect_coeff',cfg.surface_reflect_coeff,'step_m', ...
        cfg.stepz_lamb*(cfg.c0_mps/cfg.frequency_hz)/(2*pi), ...
        'x_rx_m',cfg.x_rx_m);
    one = run_pe_1d_surface_reflection_validation(one_cfg);
    exact = local_exact_1d(one_cfg);
    full = local_run_full_pe(cfg, eta_cases{ic});
    ky0 = local_extract_ky0(full);
    ay = local_discrete_y_source_factor(full.y, cfg.sigma_src_m);
    row = local_compare_case(case_names{ic}, one, exact, ky0, ay, cfg);
    row.full_pe = full;
    row.one_d = one;
    row.exact_1d = exact;
    row.eta_x_m = eta_cases{ic};
    case_rows(ic) = row;
end

checks = local_build_checks(case_rows, cfg);
bridge = struct('schema_version','1.0.0','stage','0B_1_transverse_PE_bridge', ...
    'config',cfg,'cases',case_rows,'checks',checks,'passed',all(checks.passed));
bridge.files = local_write_outputs(root, cfg, bridge);
if cfg.fail_on_check && ~bridge.passed
    error('Stage 0B 1-D PE bridge validation failed; see %s.', cfg.report_path);
end
end

function cfg = local_config(root, overrides)
cfg = struct('frequency_hz',4000,'c0_mps',1500,'z_tx_m',100,'z_rx_m',3, ...
    'sigma_src_m',0.3,'surface_reflect_coeff',-1,'xw_m',192.1875,'nx',984, ...
    'yw_m',50,'ny',256,'x_rx_m',0,'stepz_lamb',0.5, ...
    'weak_amplitude_m',1e-3,'output_dir',fullfile(root,'results','validation', ...
    'pe_1d_validation_bridge'),'report_path',fullfile(root,'reports', ...
    'pe_1d_validation_bridge_report.md'),'fail_on_check',true, ...
    'one_d_as_limit',1e-11,'ky0_direct_limit',1e-10,'ky0_reflect_limit',1e-10, ...
    'flat_sign_limit',1e-12,'run_full_pe',true);
names = fieldnames(overrides);
for ii = 1:numel(names)
    if ~isfield(cfg,names{ii}), error('Unknown bridge override: %s.',names{ii}); end
    cfg.(names{ii}) = overrides.(names{ii});
end
if ~cfg.run_full_pe
    error('Stage 0B requires run_full_pe=true for the ky=0 cross-check.');
end
end

function one = local_exact_1d(cfg)
x = (-0.5*cfg.xw_m) + (0:cfg.nx-1)*(cfg.xw_m/cfg.nx);
kx = (2*pi/cfg.xw_m)*[0:(cfg.nx/2-1),-cfg.nx/2:-1];
k0 = 2*pi*cfg.frequency_hz/cfg.c0_mps;
kz = sqrt(complex(k0^2-kx.^2,0));
source = exp(-0.5*(x/cfg.sigma_src_m).^2);
prop = @(distance) exp(1i*distance*(kz-k0));
inc = ifft(fft(source).*prop(cfg.z_tx_m));
screen = cfg.surface_reflect_coeff .* exp(1i*2*k0*cfg.surface_elevation_x_m(:).');
ref_surface = screen .* inc;
direct = ifft(fft(source).*prop(cfg.z_tx_m-cfg.z_rx_m));
reflected = ifft(fft(ref_surface).*prop(cfg.z_rx_m));
[~,ix_rx] = min(abs(x-cfg.x_rx_m));
one = struct('operator','independent one-step exact 1-D angular spectrum', ...
    'x_m',x,'direct_field',direct,'reflected_field',reflected, ...
    'surface_incident_field',inc,'surface_reflected_field',ref_surface, ...
    'direct_receiver',direct(ix_rx),'reflected_receiver',reflected(ix_rx), ...
    'k0_rad_per_m',k0,'ix_rx',ix_rx);
end

function full = local_run_full_pe(cfg, eta_x)
if ~cfg.run_full_pe, full = struct('skipped',true); return; end
p = struct('f0',cfg.frequency_hz,'enable_wideband',false,'c0',cfg.c0_mps, ...
    'z_max',cfg.z_tx_m,'z_tx',cfg.z_tx_m,'z_rx',cfg.z_rx_m, ...
    'stepz_lamb',cfg.stepz_lamb,'xw',cfg.xw_m,'yw',cfg.yw_m, ...
    'nx',cfg.nx,'ny',cfg.ny,'x_tx',0,'y_tx',0,'x_rx',cfg.x_rx_m, ...
    'y_rx',0,'sigma_src_m',cfg.sigma_src_m,'source_mode','gaussian', ...
    'sponge_ratio',0,'alpha_max_np_per_m',0,'validation_allow_extended_window',true, ...
    'validation_allow_extended_sponge_ratio',true, ...
    'env_mode','uniform','show_figures',false,'enforce_1_over_R',false, ...
    'enable_surface_reflection',true,'surface_reflect_coeff',cfg.surface_reflect_coeff, ...
    'surface_phase_mode','normal','surface_boundary_model','kirchhoff_spatial', ...
    'surface_elevation_override_xy',repmat(eta_x,cfg.ny,1), ...
    'surface_wavefield_diagnostics',true,'surface_wavefield_slice_axis','x', ...
    'surface_wavefield_max_z_samples',4,'save_mode','slice','use_gpu',false, ...
    'channel_phase_reference','direct_dsp');
out = vertical_channel_model(p);
if ~isfield(out,'psifinal_xy') || isempty(out.psifinal_xy)
    error('Full PE did not return psifinal_xy; save_mode=slice is required.');
end
if ~isfield(out.surface_wavefield_meta,'receiver_reflected_xy') || ...
        isempty(out.surface_wavefield_meta.receiver_reflected_xy)
    error('Full PE did not return receiver_reflected_xy diagnostics.');
end
full = struct('x',out.x(:).','y',out.y(:),'direct_xy',out.psifinal_xy, ...
    'reflected_xy',out.surface_wavefield_meta.receiver_reflected_xy, ...
    'h_direct_reduced',out.H_direct_reduced_f(1), ...
    'h_reflect_reduced',out.H_reflect_reduced_f(1), ...
    'surface_meta',out.surface_wavefield_meta,'config',out.config);
end

function ky0 = local_extract_ky0(full)
dy = full.y(2)-full.y(1);
ky0 = struct('direct_field',sum(full.direct_xy,1)*dy, ...
    'reflected_field',sum(full.reflected_xy,1)*dy, ...
    'dy_m',dy,'direct_receiver',[], 'reflected_receiver',[]);
end

function ay = local_discrete_y_source_factor(y, sigma)
y = y(:);
dy = y(2)-y(1);
ay = sum(exp(-0.5*(y/sigma).^2))*dy;
end

function row = local_compare_case(name, one, exact, ky0, ay, cfg)
scale = ay;
row = local_empty_row();
row.case_name = name;
row.one_d_as_direct_complex = local_rel(one.direct_field, exact.direct_field);
row.one_d_as_reflect_complex = local_rel(one.reflected_field, exact.reflected_field);
row.one_d_as_direct_phase_rad = angle(one.direct_receiver*conj(exact.direct_receiver));
row.one_d_as_reflect_phase_rad = angle(one.reflected_receiver*conj(exact.reflected_receiver));
row.ky0_direct_complex = local_rel(ky0.direct_field, scale*one.direct_field);
row.ky0_reflect_complex = local_rel(ky0.reflected_field, scale*one.reflected_field);
row.ky0_direct_receiver_complex = local_rel(ky0.direct_field(one.ix_rx), scale*one.direct_field(one.ix_rx));
row.ky0_reflect_receiver_complex = local_rel(ky0.reflected_field(one.ix_rx), scale*one.reflected_field(one.ix_rx));
row.ky0_direct_phase_rad = angle(ky0.direct_field(one.ix_rx)*conj(scale*one.direct_field(one.ix_rx)));
row.ky0_reflect_phase_rad = angle(ky0.reflected_field(one.ix_rx)*conj(scale*one.reflected_field(one.ix_rx)));
row.weak_height_max_m = max(abs(one.surface_elevation_x_m));
row.flat_reflection_sign = NaN;
if strcmp(name,'flat')
    row.flat_reflection_sign = one.reflected_receiver / ...
        local_nonreflecting_receiver(one, cfg) / cfg.surface_reflect_coeff;
end
end

function value = local_nonreflecting_receiver(one, cfg)
% Reconstruct the no-boundary direct value at the reflected image distance.
k0 = one.k0_rad_per_m;
kz = sqrt(complex(k0^2-one.kx_rad_per_m.^2,0));
prop = exp(1i*cfg.z_tx_m*(kz-k0));
inc = ifft(fft(one.source_field).*prop);
direct_image = ifft(fft(inc).*exp(1i*cfg.z_rx_m*(kz-k0)));
[~,ix] = min(abs(one.x_m-cfg.x_rx_m));
value = direct_image(ix);
end

function value = local_rel(a,b)
value = max(abs(a(:)-b(:))) / max(max(abs(b(:))),realmin);
end

function row = local_empty_row()
row = struct('case_name','','one_d_as_direct_complex',NaN, ...
    'one_d_as_reflect_complex',NaN,'one_d_as_direct_phase_rad',NaN, ...
    'one_d_as_reflect_phase_rad',NaN,'ky0_direct_complex',NaN, ...
    'ky0_reflect_complex',NaN,'ky0_direct_receiver_complex',NaN, ...
    'ky0_reflect_receiver_complex',NaN,'ky0_direct_phase_rad',NaN, ...
    'ky0_reflect_phase_rad',NaN,'flat_reflection_sign',NaN, ...
    'weak_height_max_m',NaN,'one_d',struct(),'exact_1d',struct(), ...
    'full_pe',struct(),'eta_x_m',[]);
end

function checks = local_build_checks(rows, cfg)
checks = table( ...
    ["flat_one_d_as";"weak_one_d_as";"flat_ky0_direct";"flat_ky0_reflect"; ...
     "weak_ky0_direct";"weak_ky0_reflect";"flat_reflection_sign"], ...
    [rows(1).one_d_as_direct_complex; rows(2).one_d_as_direct_complex; ...
     rows(1).ky0_direct_complex; rows(1).ky0_reflect_complex; ...
     rows(2).ky0_direct_complex; rows(2).ky0_reflect_complex; ...
     abs(rows(1).flat_reflection_sign-1)], ...
    [cfg.one_d_as_limit;cfg.one_d_as_limit;cfg.ky0_direct_limit;cfg.ky0_reflect_limit; ...
     cfg.ky0_direct_limit;cfg.ky0_reflect_limit;cfg.flat_sign_limit], ...
    'VariableNames',{'check_name','value','limit'});
checks.passed = checks.value <= checks.limit;
end

function files = local_write_outputs(root, cfg, bridge)
out = cfg.output_dir;
writetable(cfg_to_table(cfg),fullfile(out,'bridge_config.csv'));
writetable(bridge.checks,fullfile(out,'bridge_checks.csv'));
rows = bridge.cases;
summary = rmfield(rows,{'one_d','exact_1d','full_pe','eta_x_m'});
writetable(struct2table(summary),fullfile(out,'bridge_case_summary.csv'));
mat_file = fullfile(out,'pe_1d_validation_bridge.mat');
save(mat_file,'bridge','-v7.3');
local_write_report(cfg.report_path, bridge);
files = struct('mat',mat_file,'checks',fullfile(out,'bridge_checks.csv'), ...
    'summary',fullfile(out,'bridge_case_summary.csv'),'report',cfg.report_path);
end

function t = cfg_to_table(cfg)
names = {'frequency_hz','c0_mps','z_tx_m','z_rx_m','sigma_src_m','xw_m','nx','yw_m','ny','stepz_lamb','weak_amplitude_m'};
values = [cfg.frequency_hz,cfg.c0_mps,cfg.z_tx_m,cfg.z_rx_m,cfg.sigma_src_m,cfg.xw_m,cfg.nx,cfg.yw_m,cfg.ny,cfg.stepz_lamb,cfg.weak_amplitude_m];
t = table(names(:),values(:),'VariableNames',{'name','value'});
end

function local_write_report(path, bridge)
fid = fopen(path,'w','n','UTF-8');
if fid<0, error('Cannot create bridge report: %s',path); end
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'# PE 一横向维 validation bridge Stage 0B 报告\n\n');
fprintf(fid,'状态：**%s**\n\n',ternary(bridge.passed,'PASS','FAIL'));
fprintf(fid,'本阶段只验证一横向维 split-step PE、独立一步 exact angular spectrum，以及生产二维横向 PE 的 y-invariant ky=0 投影；没有运行 Bellhop，也没有修改生产 PE。\n\n');
fprintf(fid,'## 配置\n\n');
fprintf(fid,'- f/c：%.0f Hz / %.0f m/s；Tx/Rx：%.3g / %.3g m；sigma：%.3g m\n',bridge.config.frequency_hz,bridge.config.c0_mps,bridge.config.z_tx_m,bridge.config.z_rx_m,bridge.config.sigma_src_m);
fprintf(fid,'- PE bridge window/grid：%.9g m / %d；y cross-check window/grid：%.9g m / %d；sponge off\n',bridge.config.xw_m,bridge.config.nx,bridge.config.yw_m,bridge.config.ny);
fprintf(fid,'- stepz_lamb：%.6g；weak test eta amplitude：%.6g m\n\n',bridge.config.stepz_lamb,bridge.config.weak_amplitude_m);
fprintf(fid,'## Case results\n\n| case | 1D-vs-AS direct | 1D-vs-AS reflected | full-PE ky0 direct | full-PE ky0 reflected | flat sign |\n|---|---:|---:|---:|---:|---:|\n');
for ii=1:numel(bridge.cases)
    r=bridge.cases(ii);
    fprintf(fid,'| %s | %.6g | %.6g | %.6g | %.6g | %.12g |\n',r.case_name,r.one_d_as_direct_complex,r.one_d_as_reflect_complex,r.ky0_direct_complex,r.ky0_reflect_complex,r.flat_reflection_sign);
end
fprintf(fid,'\n## Checks\n\n');
for ii=1:height(bridge.checks)
    fprintf(fid,'- %s: value %.6g, limit %.6g, %s\n',bridge.checks.check_name(ii),bridge.checks.value(ii),bridge.checks.limit(ii),ternary(bridge.checks.passed(ii),'PASS','FAIL'));
end
fprintf(fid,'\n## Interpretation\n\n');
fprintf(fid,'The 1-D helper uses the same reduced square-root propagator and pressure-release phase-screen convention as the production chain. The full PE comparison is made after integrating the receiver plane over y, which selects the conserved ky=0 mode; the expected source factor is the discrete y Gaussian integral and is not fitted. This validates the dimensional bridge, not PE--Bellhop rough-surface agreement.\n');
end

function out = ternary(condition,yes_value,no_value)
if condition,out=yes_value;else,out=no_value;end
end
