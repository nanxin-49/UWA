function plot_adjoint_pe_receiver_projection_validation_vertical(v, out_dir, mode)
%PLOT_ADJOINT_PE_RECEIVER_PROJECTION_VALIDATION Write validation report figures.

arguments
    v (1,1) struct
    out_dir (1,:) char
    mode (1,:) char
end
if ~isfolder(out_dir)
    mkdir(out_dir);
end

fig = figure('Visible','off');
semilogy(v.adjoint.f_axis_hz,max(v.adjoint.relative_error,[],2),'o-'); grid on
xlabel('Frequency (Hz)'); ylabel('Max relative error');
title('Exact discrete adjoint');
local_export(fig,out_dir,['adjoint_error_' mode '.png']);

fig = figure('Visible','off'); tiledlayout(2,2);
nexttile; imagesc(abs(v.f9.analytic_stats.C_H)); axis image; colorbar; title('|C_H| analytic');
nexttile; imagesc(abs(v.f9.comparison.sample_stats.C)); axis image; colorbar; title('|C_H| sample');
nexttile; imagesc(abs(v.f9.analytic_stats.P_H)); axis image; colorbar; title('|P_H| analytic');
nexttile; imagesc(abs(v.f9.comparison.sample_stats.P)); axis image; colorbar; title('|P_H| sample');
local_export(fig,out_dir,['receiver_covariance_' mode '.png']);

fig = figure('Visible','off');
plot(v.f9.comparison.pdp_analytic,'LineWidth',1.4); hold on
plot(v.f9.comparison.pdp_sample,'--','LineWidth',1.2); grid on
legend('analytic','sample','Location','best'); title('Reflected-scatter PDP');
xlabel('Delay sample'); ylabel('Power');
local_export(fig,out_dir,['receiver_pdp_' mode '.png']);

fig = figure('Visible','off'); tiledlayout(1,2);
nexttile; plot(v.f9.comparison.lfm_analytic,'LineWidth',1.4); hold on
plot(v.f9.comparison.lfm_sample,'--','LineWidth',1.2); grid on
legend('analytic','sample','Location','best'); title('Reflected-scatter LFM power');
xlabel('Sample'); ylabel('Power');
nexttile; plot(v.f9.comparison.matched_filter_analytic,'LineWidth',1.4); hold on
plot(v.f9.comparison.matched_filter_sample,'--','LineWidth',1.2); grid on
legend('analytic','sample','Location','best'); title('Matched-filter power');
xlabel('Sample'); ylabel('Power');
local_export(fig,out_dir,['receiver_lfm_matched_filter_' mode '.png']);

fig = figure('Visible','off');
eigenvalues = sort(real(v.f9.analytic_stats.augmented_covariance_eigenvalues),'descend');
scale = max(max(abs(eigenvalues)),1);
semilogy(max(eigenvalues,eps(scale)),'o-'); grid on
xlabel('Augmented eigenvalue index'); ylabel('Eigenvalue');
title(sprintf('Augmented covariance spectrum, min = %.3g',min(eigenvalues)));
local_export(fig,out_dir,['augmented_covariance_eigenvalues_' mode '.png']);

[labels,timing_values,memory_values] = local_performance_arrays(v);
fig = figure('Visible','off'); tiledlayout(1,2);
nexttile; bar(timing_values); set(gca,'YScale','log','XTickLabel',labels); grid on
ylabel('Wall time (s, log scale)');
legend('kernel build','analytic C/P','forward / realization', ...
    'projection / realization','Location','best'); title('Timing');
nexttile; bar(memory_values/(1024^2)); set(gca,'YScale','log','XTickLabel',labels); grid on
ylabel('Memory (MiB, log scale)');
legend('q + a_{PE}','a_{PM}','FFT weights','streamed lag pair', ...
    'Location','best'); title('Feasible working arrays');
local_export(fig,out_dir,['timing_memory_' mode '.png']);
end

function [labels,timing_values,memory_values] = local_performance_arrays(v)
labels = {'F=9'};
cases = {v.f9};
if isfield(v,'f64') && isfield(v.f64,'comparison')
    labels{end+1} = 'F=64'; %#ok<AGROW>
    cases{end+1} = v.f64; %#ok<AGROW>
end
timing_values = zeros(numel(cases),4);
memory_values = zeros(numel(cases),4);
for ii = 1:numel(cases)
    item = cases{ii};
    timing_values(ii,:) = [item.projection_build.total_s, item.analytic_meta.total_s, ...
        item.ensemble_timing.cached_forward_s/item.sample_count, ...
        item.ensemble_timing.projection_s/item.sample_count];
    memory_values(ii,:) = [item.projection_build.projection_array_bytes, ...
        item.analytic_meta.a_pm_spatial_bytes, item.analytic_meta.a_pm_fft_bytes, ...
        item.analytic_meta.streamed_lag_pair_bytes];
end
end

function local_export(fig,out_dir,name)
exportgraphics(fig,fullfile(out_dir,name),'Resolution',180);
close(fig);
end
