function validation = validate_bellhop_internal_vertical_tangent_poc(options)
%VALIDATE_BELLHOP_INTERNAL_VERTICAL_TANGENT_POC
% Validate a smooth parametric wall at a point with dr/ds = 0.
% The PM validation binary is used only as a generic sampled-wall driver;
% this case is deterministic and does not generate a PM realization.
arguments
    options.output_dir (1,:) char = ''
    options.validation_exe (1,:) char = ''
    options.frequency_hz (1,1) double {mustBePositive} = 4000
    options.profile_counts (1,:) double {mustBeInteger,mustBePositive} = [65 129 257]
    options.step_m (1,1) double {mustBePositive} = 0.1
    options.beam_count (1,1) double {mustBeInteger,mustBePositive} = 2001
    options.radius_m (1,1) double {mustBePositive} = 100
    options.profile_half_span_m (1,1) double {mustBePositive} = 40
    options.wall_seed (1,1) double {mustBeFinite} = 424242
    options.fail_on_check (1,1) logical = true
end

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup_vertical_project;
if isempty(options.output_dir)
    options.output_dir = fullfile(root,'results','validation', ...
        'bellhop_internal_vertical_tangent_poc');
end
if isempty(options.validation_exe)
    options.validation_exe = fullfile(root,'results','validation', ...
        'bellhop_internal_pm_wall_poc','bin','bellhop_iwall_pm_2020.exe');
end
assert(isfile(options.validation_exe), ...
    'Missing generic sampled-wall validation binary: %s',options.validation_exe);
if ~exist(options.output_dir,'dir'), mkdir(options.output_dir); end

R0 = 100;
zc = 0;
cfg = struct('c0_mps',1500,'wall_r0_m',R0, ...
    'mapped_receiver_range_m',103,'source_depth_m',zc, ...
    'receiver_depths_m',zc,'receiver_ranges_m',[102 103], ...
    'run_type','C','angle_limits_deg',[-8 8], ...
    'domain_half_depth_m',1000,'sigma_src_m',0.3, ...
    'source_pattern_clip_db',-120,'frequency_hz',options.frequency_hz);
pat = local_pattern(cfg.frequency_hz,cfg);
cfg.source_pattern_angles_deg = pat.angles_deg;
cfg.source_pattern_level_db = pat.level_db;

template = struct('profile_count',NaN,'step_m',NaN,'beam_count',NaN, ...
    'hit_r_m',NaN,'hit_z_m',NaN,'intersection_residual_m',NaN, ...
    'tangent_error',NaN,'normal_error',NaN,'tangent_norm_error',NaN, ...
    'normal_norm_error',NaN,'tangent_normal_dot',NaN, ...
    'kappa_numeric_per_m',NaN,'kappa_analytic_per_m',NaN, ...
    'kappa_error_per_m',NaN,'u_dot_n',NaN,'specular_direction_error',NaN, ...
    'rn_model',NaN,'rn_estimated',NaN,'rn_error',NaN, ...
    'p_reflect_error',NaN,'q_reflect_error',NaN,'p_rotation_error',NaN, ...
    'q_rotation_error',NaN,'phase_jump_error_rad',NaN,'amp_jump_error',NaN, ...
    'travel_time_s',NaN,'travel_time_error_s',NaN, ...
    'min_post_range_increment_m',NaN,'all_post_range_positive',false, ...
    'finite_diagnostics',false,'case_root','');
rows = repmat(template,numel(options.profile_counts),1);
runs = cell(numel(options.profile_counts),1);

for ii = 1:numel(options.profile_counts)
    N = options.profile_counts(ii);
    z = linspace(-options.profile_half_span_m, ...
        options.profile_half_span_m,N).';
    a = options.radius_m;
    r = R0 + a - sqrt(max(a^2-z.^2,0));
    tag = sprintf('N_%d_step_%s_beams_%d',N,local_tag(options.step_m),options.beam_count);
    common = cfg;
    common.wall_seed = options.wall_seed;
    common.wall_profile_r_m = r;
    common.wall_profile_z_m = z;
    common.step_m = options.step_m;
    common.beam_count = options.beam_count;
    common.case_root = fullfile(options.output_dir,'cases',tag);
    common.bellhop_exe = options.validation_exe;
    runs{ii} = run_bellhop_internal_pm_wall_poc_vertical(common);
    d = runs{ii}.diagnostics;
    [~,i0] = min(abs(d.alpha_deg));

    hit = [d.hit_r(i0); d.hit_z(i0)];
    zhit = hit(2);
    root_term = sqrt(max(a^2-zhit^2,realmin));
    t_analytic = [zhit/root_term; 1];
    t_analytic = t_analytic / norm(t_analytic);
    n_analytic = [t_analytic(2); -t_analytic(1)];
    kappa_analytic = -1/a;

    t_num = [d.wall_t_r(i0); d.wall_t_z(i0)];
    n_num = [d.wall_n_r(i0); d.wall_n_z(i0)];
    % iwdiag stores c*t for the incident direction (a unit vector), while
    % Reflect2D's Th is the slowness component t*n = (u*n)/c.
    u_num = [d.inc_ur(i0); d.inc_uz(i0)];
    un = dot(u_num,n_num);
    ref_num = [d.ref_ur(i0); d.ref_uz(i0)];
    ref_expected = u_num - 2*dot(u_num,n_num)*n_num;
    qv = [d.q1_in(i0); d.q2_in(i0)];
    dp = [d.p1_ref(i0)-d.p1_in(i0); d.p2_ref(i0)-d.p2_in(i0)];
    if dot(qv,qv) > 1e-24
        rn_est = dot(dp,qv)/dot(qv,qv);
    else
        rn_est = NaN;
    end
    th = un / cfg.c0_mps;
    rn_model = -2*kappa_analytic/(cfg.c0_mps^2*th);

    q = template;
    q.profile_count = N; q.step_m = options.step_m; q.beam_count = options.beam_count;
    q.hit_r_m = hit(1); q.hit_z_m = hit(2);
    q.intersection_residual_m = abs(d.wall_residual(i0));
    q.tangent_error = norm(t_num-t_analytic);
    q.normal_error = norm(n_num-n_analytic);
    q.tangent_norm_error = abs(norm(t_num)-1);
    q.normal_norm_error = abs(norm(n_num)-1);
    q.tangent_normal_dot = dot(t_num,n_num);
    q.kappa_numeric_per_m = d.kappa(i0);
    q.kappa_analytic_per_m = kappa_analytic;
    q.kappa_error_per_m = d.kappa(i0)-kappa_analytic;
    q.u_dot_n = un;
    q.specular_direction_error = norm(ref_num-ref_expected);
    q.rn_model = rn_model; q.rn_estimated = rn_est; q.rn_error = rn_est-rn_model;
    q.p_reflect_error = d.p_ref_error(i0); q.q_reflect_error = d.q_ref_error(i0);
    q.p_rotation_error = d.p_rot_error(i0); q.q_rotation_error = d.q_rot_error(i0);
    q.phase_jump_error_rad = abs(angle(exp(1i*(d.phase_delta(i0)-pi))));
    q.amp_jump_error = abs(d.amp_delta(i0));
    q.travel_time_s = d.tau_receiver_real(i0);
    q.travel_time_error_s = q.travel_time_s - 103/cfg.c0_mps;
    q.min_post_range_increment_m = min(d.min_post_dr);
    q.all_post_range_positive = all(d.min_post_dr>0);
    numeric = d{:,setdiff(d.Properties.VariableNames,{'mode'})};
    q.finite_diagnostics = all(isfinite(numeric),'all');
    q.case_root = common.case_root;
    rows(ii) = q;
end

table_rows = struct2table(rows);
checks = struct;
checks.finite = all(table_rows.finite_diagnostics);
checks.intersection = max(table_rows.intersection_residual_m) <= 1e-8;
checks.frame = max([table_rows.tangent_error;table_rows.normal_error; ...
    table_rows.tangent_norm_error;table_rows.normal_norm_error; ...
    abs(table_rows.tangent_normal_dot)]) <= 5e-10;
checks.curvature = max(abs(table_rows.kappa_error_per_m)) <= 2e-4;
checks.non_grazing = min(table_rows.u_dot_n) >= 0.9;
checks.direction = max(table_rows.specular_direction_error) <= 1e-10;
checks.rn = max(abs(table_rows.rn_error)) <= 2e-8;
checks.dynamic_state = max([table_rows.q_reflect_error;table_rows.p_rotation_error; ...
    table_rows.q_rotation_error]) <= 1e-10;
checks.phase = max(table_rows.phase_jump_error_rad) <= 1e-10;
checks.amplitude = max(table_rows.amp_jump_error) <= 1e-10;
checks.travel_time = max(abs(table_rows.travel_time_error_s)) <= 2e-9;
checks.positive_range = all(table_rows.all_post_range_positive);
curv_abs = abs(table_rows.kappa_error_per_m);
checks.curvature_convergence = curv_abs(end) <= curv_abs(1) && ...
    abs(table_rows.rn_error(end)) <= abs(table_rows.rn_error(1)) + 1e-10;
checks.all = all(structfun(@(x)logical(x),checks));

csv_file = fullfile(options.output_dir,'vertical_tangent_convergence.csv');
writetable(table_rows,csv_file);
validation = struct('config',cfg,'options',options,'table',table_rows, ...
    'checks',checks,'runs',{runs},'files',struct('convergence_csv',csv_file, ...
    'mat',fullfile(options.output_dir,'vertical_tangent_poc_results.mat')));
save(validation.files.mat,'validation','-v7.3');
disp(table_rows(:,{'profile_count','intersection_residual_m','tangent_error', ...
    'normal_error','kappa_numeric_per_m','kappa_error_per_m','u_dot_n', ...
    'specular_direction_error','rn_error','phase_jump_error_rad', ...
    'travel_time_error_s','min_post_range_increment_m'}));
disp(checks);
if options.fail_on_check && ~checks.all
    error('Vertical-tangent internal-wall POC failed; inspect %s.',options.output_dir);
end
end

function pat = local_pattern(f,cfg)
angles = linspace(cfg.angle_limits_deg(1),cfg.angle_limits_deg(2),2401).';
k = 2*pi*f/cfg.c0_mps; th = deg2rad(angles);
D = cos(th).*exp(-0.5*(k*cfg.sigma_src_m*sin(th)).^2);
D = abs(D)/max(abs(D));
pat = struct('angles_deg',angles, ...
    'level_db',20*log10(max(D,10^(cfg.source_pattern_clip_db/20))));
end

function tag = local_tag(value)
tag = strrep(sprintf('%.10g',value),'.','p');
tag = strrep(tag,'-','m');
end
