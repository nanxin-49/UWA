run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'bootstrap_project.m'));
% Compare reduced-grid vertical channel responses across bubble models.

clear
format compact

result_file = 'compare_bubble_models_vertical_result.mat';
figure_prefix = 'compare_bubble_models_vertical_';
if exist(result_file, 'file')
    warning('Result file %s already exists and will be overwritten.', result_file);
end

paramsV = struct();
paramsV.f0 = 4000;
paramsV.enable_wideband = true;
paramsV.f_band_hz = [4000, 8000];
paramsV.Nf_min = 16;
paramsV.Nf_max = 16;
paramsV.f_ref_hz = 6000;
paramsV.c0 = 1500;
paramsV.z_max = 100;
paramsV.stepz_lamb = 0.5;
paramsV.xw = 50;
paramsV.yw = 50;
paramsV.nx = 128;
paramsV.ny = 128;
paramsV.x_tx = 0;
paramsV.y_tx = 0;
paramsV.z_tx = 100;
paramsV.x_rx = 0;
paramsV.y_rx = 0;
paramsV.z_rx = 3;
paramsV.nout = 2;
paramsV.sigma_src_m = 0.4;
paramsV.sponge_ratio = 0.12;
paramsV.alpha_max_np_per_m = 0.15;
paramsV.env_mode = 'uniform';
paramsV.enforce_1_over_R = false;
paramsV.show_figures = true;
paramsV.save_mode = 'rx_only';
paramsV.use_gpu = false;
paramsV.enable_surface_reflection = true;
paramsV.sea_wind_speed = 5.0;
paramsV.sea_hs_target = 0.5;
paramsV.sea_seed = 12345;

scenarios = struct([]);
scenarios(1).name = 'no_bubble';
scenarios(1).enable_bubbles = false;
scenarios(1).bubble_model = 'off';
scenarios(1).bubble_spatial_mode = 'none';
scenarios(1).bubble_strength_scale = 1;
scenarios(1).bubble_alpha0_np_per_m = 0;
scenarios(1).bubble_layer_decay_m = 0.4;
scenarios(1).bubble_delta_c0_mps = 0;

scenarios(2).name = 'level0_empirical';
scenarios(2).enable_bubbles = true;
scenarios(2).bubble_model = 'level0_empirical';
scenarios(2).bubble_spatial_mode = '1d';
scenarios(2).bubble_strength_scale = 1;
scenarios(2).bubble_alpha0_np_per_m = 0.02;
scenarios(2).bubble_layer_decay_m = 20;
scenarios(2).bubble_delta_c0_mps = 0;

scenarios(3).name = 'hall1d_default';
scenarios(3).enable_bubbles = true;
scenarios(3).bubble_model = 'hall1d';
scenarios(3).bubble_spatial_mode = '1d';
scenarios(3).bubble_strength_scale = 1;
scenarios(3).bubble_alpha0_np_per_m = 0;
scenarios(3).bubble_layer_decay_m = 0.4;
scenarios(3).bubble_delta_c0_mps = 0;

scenarios(4).name = 'hall1d_strong';
scenarios(4).enable_bubbles = true;
scenarios(4).bubble_model = 'hall1d';
scenarios(4).bubble_spatial_mode = '1d';
scenarios(4).bubble_strength_scale = 1e6;
scenarios(4).bubble_alpha0_np_per_m = 0;
scenarios(4).bubble_layer_decay_m = 0.4;
scenarios(4).bubble_delta_c0_mps = 0;

results = struct([]);
summary = struct([]);
for ss = 1:numel(scenarios)
    p = paramsV;
    p.enable_bubbles = scenarios(ss).enable_bubbles;
    p.bubble_model = scenarios(ss).bubble_model;
    p.bubble_spatial_mode = scenarios(ss).bubble_spatial_mode;
    p.bubble_strength_scale = scenarios(ss).bubble_strength_scale;
    p.bubble_alpha0_np_per_m = scenarios(ss).bubble_alpha0_np_per_m;
    p.bubble_layer_decay_m = scenarios(ss).bubble_layer_decay_m;
    p.bubble_delta_c0_mps = scenarios(ss).bubble_delta_c0_mps;

    channel = vertical_channel_model(p);
    invariant_error = max(abs(channel.H_f - (channel.H_direct_f + channel.H_reflect_f)));
    if invariant_error > 1e-10
        error('compare_bubble_models_vertical:Invariant', ...
            'H_f invariant failed for scenario %s: %.3e', ...
            scenarios(ss).name, invariant_error);
    end

    results(ss).name = scenarios(ss).name;
    results(ss).paramsV = p;
    results(ss).channel = channel;
    results(ss).f_axis = channel.f_axis;
    results(ss).H_f = channel.H_f;
    results(ss).H_direct_f = channel.H_direct_f;
    results(ss).H_reflect_f = channel.H_reflect_f;
    results(ss).bubble_meta = channel.bubble_meta;
    results(ss).invariant_error = invariant_error;

    summary(ss).name = scenarios(ss).name;
    summary(ss).h_total_abs = abs(channel.h_total);
    summary(ss).h_total_phase_rad = angle(channel.h_total);
    summary(ss).max_delta_TL_dB = 0;
    summary(ss).max_alpha_bub = local_get_stat(channel.bubble_meta, 'alpha_bub_stats', 'max');
    summary(ss).max_beta = local_get_stat(channel.bubble_meta, 'beta_stats', 'max');
    summary(ss).invariant_error = invariant_error;
end

f_axis = results(1).f_axis(:);
H0_f = results(1).H_f(:);
H0_direct_f = results(1).H_direct_f(:);
H0_reflect_f = results(1).H_reflect_f(:);
H0_safe_f = H0_f;
H0_safe_f(abs(H0_safe_f) <= eps) = eps;
H0_direct_safe_f = H0_direct_f;
H0_direct_safe_f(abs(H0_direct_safe_f) <= eps) = eps;

for ss = 1:numel(results)
    H_f = results(ss).H_f(:);
    H_direct_f = results(ss).H_direct_f(:);
    H_reflect_f = results(ss).H_reflect_f(:);

    results(ss).H_ratio = H_f ./ H0_safe_f;
    results(ss).direct_ratio = H_direct_f ./ H0_direct_safe_f;
    reflect_den = H0_reflect_f;
    reflect_mask = abs(reflect_den) > eps;
    reflect_ratio = complex(NaN(size(H_reflect_f)), NaN(size(H_reflect_f)));
    reflect_ratio(reflect_mask) = H_reflect_f(reflect_mask) ./ reflect_den(reflect_mask);
    results(ss).reflect_ratio = reflect_ratio;
    results(ss).delta_TL_total_dB = -20*log10(abs(H_f) ./ max(abs(H0_f), eps));
    results(ss).delta_phase_total_rad = unwrap(angle(H_f)) - unwrap(angle(H0_f));
    summary(ss).max_delta_TL_dB = max(results(ss).delta_TL_total_dB);
end

summary_table = struct2table(summary);
disp(summary_table)

local_plot_magnitude(f_axis, results, figure_prefix);
local_plot_phase(f_axis, results, figure_prefix);
local_plot_delta_tl(f_axis, results, figure_prefix);
local_plot_components(f_axis, results, figure_prefix);
local_plot_bubble_meta(results, figure_prefix);

save(result_file, 'paramsV', 'scenarios', 'results', 'summary_table');

function local_plot_magnitude(f_axis, results, figure_prefix)
figure(41); clf
hold on
for ss = 1:numel(results)
    plot(f_axis, 20*log10(max(abs(results(ss).H_f), eps)), 'LineWidth', 1.2)
end
grid on
xlabel('Frequency (Hz)')
ylabel('|H(f)| (dB)')
title('Total channel magnitude')
legend({results.name}, 'Interpreter', 'none', 'Location', 'best')
print(41, '-dpng', '-r200', [figure_prefix 'Figure41_H_magnitude.png'])
end

function local_plot_phase(f_axis, results, figure_prefix)
figure(42); clf
hold on
for ss = 1:numel(results)
    plot(f_axis, unwrap(angle(results(ss).H_f)), 'LineWidth', 1.2)
end
grid on
xlabel('Frequency (Hz)')
ylabel('Unwrapped phase (rad)')
title('Total channel phase')
legend({results.name}, 'Interpreter', 'none', 'Location', 'best')
print(42, '-dpng', '-r200', [figure_prefix 'Figure42_H_phase.png'])
end

function local_plot_delta_tl(f_axis, results, figure_prefix)
figure(43); clf
hold on
for ss = 1:numel(results)
    plot(f_axis, results(ss).delta_TL_total_dB, 'LineWidth', 1.2)
end
grid on
xlabel('Frequency (Hz)')
ylabel('\DeltaTL relative to no bubble (dB)')
title('Bubble excess transmission loss')
legend({results.name}, 'Interpreter', 'none', 'Location', 'best')
print(43, '-dpng', '-r200', [figure_prefix 'Figure43_delta_TL.png'])
end

function local_plot_components(f_axis, results, figure_prefix)
figure(44); clf
tiledlayout(2, 1)
nexttile
hold on
for ss = 1:numel(results)
    plot(f_axis, 20*log10(max(abs(results(ss).H_direct_f), eps)), 'LineWidth', 1.2)
end
grid on
xlabel('Frequency (Hz)')
ylabel('|H_{direct}| (dB)')
title('Direct component')
legend({results.name}, 'Interpreter', 'none', 'Location', 'best')

nexttile
hold on
for ss = 1:numel(results)
    plot(f_axis, 20*log10(max(abs(results(ss).H_reflect_f), eps)), 'LineWidth', 1.2)
end
grid on
xlabel('Frequency (Hz)')
ylabel('|H_{reflect}| (dB)')
title('Reflected component')
legend({results.name}, 'Interpreter', 'none', 'Location', 'best')
print(44, '-dpng', '-r200', [figure_prefix 'Figure44_components.png'])
end

function local_plot_bubble_meta(results, figure_prefix)
max_alpha = zeros(numel(results), 1);
max_beta = zeros(numel(results), 1);
for ss = 1:numel(results)
    max_alpha(ss) = local_get_stat(results(ss).bubble_meta, 'alpha_bub_stats', 'max');
    max_beta(ss) = local_get_stat(results(ss).bubble_meta, 'beta_stats', 'max');
end

figure(45); clf
tiledlayout(2, 1)
nexttile
bar(max_alpha)
grid on
set(gca, 'XTick', 1:numel(results), 'XTickLabel', {results.name})
xtickangle(20)
ylabel('max \alpha_{bub} (Np/m)')
title('Bubble attenuation summary')

nexttile
bar(max_beta)
grid on
set(gca, 'XTick', 1:numel(results), 'XTickLabel', {results.name})
xtickangle(20)
ylabel('max \beta')
title('Void fraction summary')
print(45, '-dpng', '-r200', [figure_prefix 'Figure45_bubble_meta.png'])
end

function value = local_get_stat(meta, stats_field, stat_name)
value = NaN;
if isstruct(meta) && isfield(meta, stats_field)
    stats = meta.(stats_field);
    if isstruct(stats) && isfield(stats, stat_name)
        value = stats.(stat_name);
    end
end
end

