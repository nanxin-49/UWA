function result = validate_bellhop_internal_large_sinusoidal_wall_visual_case(overrides)
%VALIDATE_BELLHOP_INTERNAL_LARGE_SINUSOIDAL_WALL_VISUAL_CASE
% Run one full-fan, two-period sinusoidal internal-wall case for visualization.
% This validation-only entrypoint runs Bellhop only; it does not run PE and
% does not modify Bellhop core physics.  The existing validated tilted case
% is reused only as a second panel in the trajectory renderer.
if nargin < 1 || isempty(overrides), overrides = struct(); end
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
run(fullfile(root,'setup_vertical_project.m'));
cfg = struct( ...
    'output_dir',fullfile(root,'results','validation','bellhop_internal_wall_two_period_full_fan'), ...
    'visual_output_dir',fullfile(root,'results','visualization','bellhop_internal_wall_visuals'), ...
    'internal_exe',fullfile(root,'results','validation','bellhop_internal_sinusoidal_wall_poc', ...
        'bin','bellhop_iwall_sinusoidal_2020.exe'), ...
    'official_exe','E:\MISC\BELLHOP\AcousticsToolbox_2020\windows-bin-20201102\bellhop.exe', ...
    'frequency_hz',4000,'c0_mps',1500,'wall_r0_m',100,'wall_amplitude_m',2.4, ...
    'wall_wavenumber_per_m',2*pi/30,'wall_z_support_m',[-75 75], ...
    'plot_z_limits_m',[-65 65],'wall_zoom_z_limits_m',[-30 30], ...
    'wall_profile_count',3001,'receiver_depths_m',-65:0.5:65, ...
    'mapped_receiver_ranges_m',sort([100.05:0.1:129.95 103]), ...
    'incident_receiver_ranges_m',0.05:0.1:104.05, ...
    'beam_count',5001,'angle_limits_deg',[-30 30],'step_m',0.05, ...
    'domain_half_depth_m',1000,'sigma_src_m',0.3,'source_pattern_clip_db',-120, ...
    'ray_angles_deg',[-28 -20 -12 -4 0 4 12 20 28], ...
    'target_eigenray_angle_deg',1.23855637169814, ...
    'physical_receiver_range_m',97,'mapped_receiver_range_m',103, ...
    'report_file',fullfile(root,'reports','bellhop_internal_wall_visualization_report.md'));
names = fieldnames(overrides);
for ii = 1:numel(names)
    if ~isfield(cfg,names{ii}), error('Unknown override: %s.',names{ii}); end
    cfg.(names{ii}) = overrides.(names{ii});
end
assert(isfile(cfg.internal_exe),'Missing internal-wall Bellhop executable: %s',cfg.internal_exe);
assert(isfile(cfg.official_exe),'Missing official Bellhop executable: %s',cfg.official_exe);
if ~isfolder(cfg.output_dir), mkdir(cfg.output_dir); end
if ~isfolder(cfg.visual_output_dir), mkdir(cfg.visual_output_dir); end

z_wall = linspace(cfg.wall_z_support_m(1),cfg.wall_z_support_m(2),cfg.wall_profile_count).';
r_wall = cfg.wall_r0_m-cfg.wall_amplitude_m*sin(cfg.wall_wavenumber_per_m*z_wall);
writetable(table(z_wall,r_wall,'VariableNames',{'z_m','r_wall_m'}), ...
    fullfile(cfg.output_dir,'large_sinusoidal_wall_profile.csv'));
preflight = local_geometry_preflight(cfg,z_wall,r_wall);
preflight_file = fullfile(cfg.output_dir,'full_fan_geometry_preflight.csv');
writetable(preflight.table,preflight_file);
fprintf('Geometry preflight: hits=%d/%d, min receiver margin=%.6g m, ', ...
    nnz(preflight.table.found),height(preflight.table),preflight.min_receiver_margin_m);
fprintf('min rotated u_r=%.6g, min |u dot n|=%.6g.\n', ...
    preflight.min_rotated_ur,preflight.min_abs_incident_normal);
if ~preflight.pass
    error('Full-fan geometry preflight failed; Bellhop was not started. Inspect %s.',preflight_file);
end
pat = local_pattern(cfg.frequency_hz,cfg);
base = struct('run_type','C','frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps, ...
    'source_depth_m',0,'receiver_depths_m',cfg.receiver_depths_m, ...
    'beam_count',cfg.beam_count,'angle_limits_deg',cfg.angle_limits_deg, ...
    'step_m',cfg.step_m,'domain_half_depth_m',cfg.domain_half_depth_m, ...
    'sigma_src_m',cfg.sigma_src_m,'source_pattern_clip_db',cfg.source_pattern_clip_db, ...
    'source_pattern_angles_deg',pat.angles_deg,'source_pattern_level_db',pat.level_db, ...
    'source_geometry','R');

wall = base;
wall.bellhop_exe = cfg.internal_exe;
wall.case_root = fullfile(cfg.output_dir,'cases','internal_large_sinusoid');
wall.wall_r0_m = cfg.wall_r0_m;
wall.wall_amplitude_m = cfg.wall_amplitude_m;
wall.wall_wavenumber_per_m = cfg.wall_wavenumber_per_m;
wall.wall_profile_r_m = r_wall;
wall.wall_profile_z_m = z_wall;
wall.mapped_receiver_range_m = cfg.mapped_receiver_range_m;
wall.receiver_ranges_m = cfg.mapped_receiver_ranges_m;
rs = run_bellhop_internal_sinusoidal_wall_poc_vertical(wall);

incident = base;
incident.bellhop_exe = cfg.official_exe;
incident.case_root = fullfile(cfg.output_dir,'cases','incident_dense');
incident.receiver_ranges_m = cfg.incident_receiver_ranges_m;
ri = run_bellhop_unfolded_gaussian_vertical(incident);

% Reuse the already generated tilted dense field only as the comparison panel
% in the existing renderer; no tilted propagation is rerun here.
old_mat = fullfile(root,'results','validation','bellhop_internal_wall_visuals','dense_wall_fields.mat');
if ~isfile(old_mat), error('Missing retained tilted visualization data: %s',old_mat); end
old = load(old_mat,'rt');
dense_mat = fullfile(cfg.output_dir,'large_sinusoidal_wall_dense_fields.mat');
rr = cfg.mapped_receiver_ranges_m;
z = cfg.receiver_depths_m;
save(dense_mat,'ri','rs','old','z','rr','-v7.3');

% Match the renderer's expected field names without modifying the old MAT.
S = load(dense_mat,'ri','rs','old','z','rr');
rt = S.old.rt; %#ok<NASGU>
save(dense_mat,'ri','rs','rt','z','rr','-v7.3');
visuals = generate_bellhop_internal_wall_visuals_vertical(struct( ...
    'input_mat',dense_mat,'output_dir',cfg.visual_output_dir,'r0_m',cfg.wall_r0_m, ...
    'physical_rx_range_m',cfg.physical_receiver_range_m, ...
    'z_plot_m',cfg.plot_z_limits_m,'wall_zoom_z_m',cfg.wall_zoom_z_limits_m, ...
    'wall_support_z_m',cfg.wall_z_support_m, ...
    'ray_angles_deg',cfg.ray_angles_deg, ...
    'target_eigenray_angle_deg',cfg.target_eigenray_angle_deg, ...
    'sinusoidal_amplitude_m',cfg.wall_amplitude_m, ...
    'sinusoidal_wavenumber_per_m',cfg.wall_wavenumber_per_m, ...
    'report_file',cfg.report_file));

result = struct('config',cfg,'preflight',preflight,'wall',rs,'incident',ri,'visuals',visuals, ...
    'files',struct('dense_mat',dense_mat,'profile_csv',fullfile(cfg.output_dir, ...
    'large_sinusoidal_wall_profile.csv'),'preflight_csv',preflight_file));
save(fullfile(cfg.output_dir,'large_sinusoidal_wall_visual_case.mat'),'result','-v7.3');
d = rs.diagnostics;
fprintf('Large sinusoidal Bellhop case complete: A=%.6g m, K=%.6g 1/m, N=%d.\n', ...
    cfg.wall_amplitude_m,cfg.wall_wavenumber_per_m,cfg.wall_profile_count);
fprintf('max wall residual=%.3e m, min post-wall dr=%.6g m, finite diagnostics=%d/%d.\n', ...
    max(abs(d.wall_residual)),min(d.min_post_dr),nnz(all(isfinite(d{:,1:end}),2)),height(d));
fprintf('No PE run was performed.\n');
end

function pat = local_pattern(f,cfg)
angles = linspace(cfg.angle_limits_deg(1),cfg.angle_limits_deg(2),2401).';
k = 2*pi*f/cfg.c0_mps; th = deg2rad(angles);
D = cos(th).*exp(-0.5*(k*cfg.sigma_src_m*sin(th)).^2); D = abs(D)/max(abs(D));
pat = struct('angles_deg',angles,'level_db',20*log10(max(D,10^(cfg.source_pattern_clip_db/20))));
end

function out = local_geometry_preflight(cfg,z_wall,r_wall)
% Match the validation binary's first-positive segment intersection, then
% verify the fixed receiver plane and proper-rotation range direction.
alpha_deg = linspace(cfg.angle_limits_deg(1),cfg.angle_limits_deg(2),cfg.beam_count).';
n_beams = numel(alpha_deg);
found = false(n_beams,1); n_intersections = zeros(n_beams,1);
hit_r_m = NaN(n_beams,1); hit_z_m = NaN(n_beams,1);
analytic_residual_m = NaN(n_beams,1); receiver_margin_m = NaN(n_beams,1);
support_margin_m = NaN(n_beams,1); incident_dot_normal = NaN(n_beams,1);
reflected_ur = NaN(n_beams,1); reflected_uz = NaN(n_beams,1);
rotated_ur = NaN(n_beams,1); rotated_uz = NaN(n_beams,1);
d_r = diff(r_wall); d_z = diff(z_wall);
r_seg = r_wall(1:end-1); z_seg = z_wall(1:end-1);
tol = 1e-10;
for ii = 1:n_beams
    alpha = deg2rad(alpha_deg(ii)); u = [cos(alpha);sin(alpha)];
    den = u(1)*d_z-u(2)*d_r;
    h = (r_seg.*d_z-z_seg.*d_r)./den;
    lambda = (r_seg*u(2)-z_seg*u(1))./den;
    valid = abs(den)>1e-14 & h>tol & lambda>=-tol & lambda<=1+tol;
    n_intersections(ii) = nnz(valid);
    if ~any(valid), continue; end
    candidates = find(valid); [h_hit,jj] = min(h(candidates));
    seg = candidates(jj); lam = min(1,max(0,lambda(seg)));
    found(ii) = true;
    hit_r_m(ii) = r_seg(seg)+lam*d_r(seg);
    hit_z_m(ii) = z_seg(seg)+lam*d_z(seg);
    analytic_residual_m(ii) = hit_r_m(ii)-(cfg.wall_r0_m- ...
        cfg.wall_amplitude_m*sin(cfg.wall_wavenumber_per_m*hit_z_m(ii)));
    slope = -cfg.wall_amplitude_m*cfg.wall_wavenumber_per_m* ...
        cos(cfg.wall_wavenumber_per_m*hit_z_m(ii));
    tangent = [slope;1]; tangent = tangent/norm(tangent);
    normal = [tangent(2);-tangent(1)];
    incident_dot_normal(ii) = dot(u,normal);
    reflected = u-2*incident_dot_normal(ii)*normal;
    rotated = -reflected;
    reflected_ur(ii) = reflected(1); reflected_uz(ii) = reflected(2);
    rotated_ur(ii) = rotated(1); rotated_uz(ii) = rotated(2);
    receiver_margin_m(ii) = cfg.mapped_receiver_range_m-(2*cfg.wall_r0_m-hit_r_m(ii));
    support_margin_m(ii) = min(hit_z_m(ii)-cfg.wall_z_support_m(1), ...
        cfg.wall_z_support_m(2)-hit_z_m(ii));
    if abs(h_hit-hypot(hit_r_m(ii),hit_z_m(ii)))>1e-6
        error('Internal preflight path-length inconsistency at alpha %.6g deg.',alpha_deg(ii));
    end
end
T = table(alpha_deg,found,n_intersections,hit_r_m,hit_z_m,analytic_residual_m, ...
    receiver_margin_m,support_margin_m,incident_dot_normal,reflected_ur,reflected_uz, ...
    rotated_ur,rotated_uz);
valid_rows = found & all(isfinite(T{:,4:end}),2);
out = struct('table',T, ...
    'min_receiver_margin_m',min(receiver_margin_m(valid_rows)), ...
    'min_support_margin_m',min(support_margin_m(valid_rows)), ...
    'min_rotated_ur',min(rotated_ur(valid_rows)), ...
    'min_abs_incident_normal',min(abs(incident_dot_normal(valid_rows))), ...
    'max_analytic_residual_m',max(abs(analytic_residual_m(valid_rows))));
out.pass = all(valid_rows) && out.min_receiver_margin_m>0.1 && ...
    out.min_support_margin_m>5 && out.min_rotated_ur>0.05 && ...
    out.min_abs_incident_normal>0.2 && out.max_analytic_residual_m<1e-3;
end
