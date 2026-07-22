%% Independent Li et al. (2009) explicit rough-surface acoustic validation
% No communication modulation, noise, circuit, kstat, or SSA path is used.
clear; clc;
root_dir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root_dir);
output_dir = fullfile(root_dir,'results','validation','li2009_explicit_surface');
if ~exist(output_dir,'dir'), mkdir(output_dir); end

%% Cost-controlled baseline: raw PM grid is wider than the PE crop.
baseline = struct( ...
    'run_label','baseline_128', ...
    'output_dir',output_dir, ...
    'H_m',256, ...
    'pe_nx',128,'pe_ny',128,'pe_xw_m',100,'pe_yw_m',100, ...
    'pm_nx',256,'pm_ny',256,'pm_xw_m',200,'pm_yw_m',200, ...
    'dz_m',4,'Nf',65,'sea_seed_list',41001+(0:127), ...
    'sample_convergence_counts',[8,16,32,64,128], ...
    'U10_list_mps',[5,10],'make_figures',true,'save_artifacts',true);
baseline_path = fullfile(output_dir,'baseline_128_result.mat');
if exist(baseline_path,'file')
    loaded = load(baseline_path,'result');
    baseline_result = loaded.result;
else
    baseline_result = li2009_explicit_surface_validation(baseline);
end

%% Numerical convergence at U10=10 m/s.
% Variants use 16 samples to distinguish numerical shifts from the full
% baseline's Monte Carlo convergence table. Fine-grid and paper-depth runs
% are intentionally separate because their costs scale differently.
common = baseline;
common.U10_list_mps = 10;
common.sea_seed_list = 42001+(0:15);
common.make_figures = false;
variants = { ...
    struct('run_label','matched_reference_16'), ...
    struct('run_label','pm_aperture_100m','pm_nx',128,'pm_ny',128, ...
        'pm_xw_m',100,'pm_yw_m',100), ...
    struct('run_label','dz_8m','dz_m',8), ...
    struct('run_label','sponge_10pct','sponge_ratio',0.10), ...
    struct('run_label','sponge_15pct','sponge_ratio',0.15), ...
    struct('run_label','frequency_33','Nf',33), ...
    struct('run_label','window_150m','pe_nx',192,'pe_ny',192, ...
        'pe_xw_m',150,'pe_yw_m',150), ...
    struct('run_label','grid_coarse_matched','Nf',33, ...
        'sea_seed_list',43001+(0:7)), ...
    struct('run_label','grid_dx_0p520833m','pe_nx',192,'pe_ny',192, ...
        'pe_xw_m',100,'pe_yw_m',100,'pm_nx',384,'pm_ny',384, ...
        'pm_xw_m',200,'pm_yw_m',200,'Nf',33,'sea_seed_list',43001+(0:7)), ...
    struct('run_label','grid_dx_0p390625m','pe_nx',256,'pe_ny',256, ...
        'pe_xw_m',100,'pe_yw_m',100,'pm_nx',512,'pm_ny',512, ...
        'pm_xw_m',200,'pm_yw_m',200,'Nf',33,'sea_seed_list',43001+(0:7))};

convergence_rows = cell(numel(variants)+1,1);
convergence_rows{1} = local_convergence_row('baseline',baseline_result,10);
variant_results = cell(numel(variants),1);
for iv = 1:numel(variants)
    p = local_merge(common,variants{iv});
    result_path = fullfile(output_dir,[p.run_label,'_result.mat']);
    if exist(result_path,'file')
        loaded = load(result_path,'result');
        variant_results{iv} = loaded.result;
    else
        variant_results{iv} = li2009_explicit_surface_validation(p);
    end
    convergence_rows{iv+1} = local_convergence_row(p.run_label,variant_results{iv},10);
end
convergence_table = struct2table(vertcat(convergence_rows{:}));
writetable(convergence_table,fullfile(output_dir,'numerical_convergence_summary.csv'));

%% Gradual paper-geometry approach: correct 1024-m depth/domain, coarse grid.
paper_depth = baseline;
paper_depth.run_label = 'paper_depth_coarse_grid';
paper_depth.H_m = 1024;
paper_depth.pe_nx = 256; paper_depth.pe_ny = 256;
paper_depth.pe_xw_m = 200; paper_depth.pe_yw_m = 200;
paper_depth.pm_nx = 256; paper_depth.pm_ny = 256;
paper_depth.pm_xw_m = 200; paper_depth.pm_yw_m = 200;
paper_depth.Nf = 33;
paper_depth.sea_seed_list = 44001+(0:7);
paper_depth.make_figures = false;
paper_depth_path = fullfile(output_dir,'paper_depth_coarse_grid_result.mat');
if exist(paper_depth_path,'file')
    loaded = load(paper_depth_path,'result');
    paper_depth_result = loaded.result;
else
    paper_depth_result = li2009_explicit_surface_validation(paper_depth);
end

save(fullfile(output_dir,'li2009_validation_suite.mat'), ...
    'baseline_result','variant_results','convergence_table','paper_depth_result','-v7.3');
local_plot_convergence(convergence_table,output_dir);

fprintf('\nLi2009 validation acceptance (baseline):\n');
disp(baseline_result.acceptance);
fprintf('Paper-depth coarse-grid flat peak error: %.6g ms\n', ...
    1e3*paper_depth_result.flat.tau0_error_peak_s);

function out=local_merge(base,override)
out=base; names=fieldnames(override);
for ii=1:numel(names), out.(names{ii})=override.(names{ii}); end
end

function row=local_convergence_row(label,r,U10)
[~,iu]=min(abs(r.U10_list_mps-U10));
row=struct('case_label',label,'H_m',r.config.H_m, ...
    'pe_nx',r.config.pe_nx,'pe_xw_m',r.config.pe_xw_m, ...
    'pm_nx',r.config.pm_nx,'pm_xw_m',r.config.pm_xw_m, ...
    'dx_m',r.config.dx_m,'dz_m',r.config.dz_m, ...
    'sponge_ratio',r.config.sponge_ratio,'Nf',r.config.Nf, ...
    'sample_count',numel(r.sea_seed_list), ...
    'pm_capture_ratio',r.surface_audit(iu).spectrum_capture_ratio, ...
    'flat_peak_error_ms',1e3*r.flat.tau0_error_peak_s, ...
    'peak_mean_offset_ms',1e3*(mean(r.peak_time_s(:,iu))-r.geometry.tau0_s), ...
    'peak_std_ms',1e3*std(r.peak_time_s(:,iu),1), ...
    'first_mean_offset_ms',1e3*(mean(r.first_threshold_time_s(:,iu))-r.geometry.tau0_s), ...
    'first_std_ms',1e3*std(r.first_threshold_time_s(:,iu),1));
end

function local_plot_convergence(tbl,out)
f=figure('Visible','off','Color','w');
yyaxis left; plot(1:height(tbl),tbl.peak_std_ms,'o-','LineWidth',1.2); hold on
plot(1:height(tbl),tbl.first_std_ms,'s-','LineWidth',1.2);
ylabel('Travel-time standard deviation (ms)');
yyaxis right; plot(1:height(tbl),tbl.pm_capture_ratio,'d--','LineWidth',1.2);
ylabel('PM discrete/infinite variance ratio');
xticks(1:height(tbl)); xticklabels(tbl.case_label); xtickangle(30); grid on
legend('peak width','first-threshold width','PM capture','Location','best');
title('Li2009 independent numerical convergence comparisons');
exportgraphics(f,fullfile(out,'numerical_convergence.png'),'Resolution',180); close(f)
end
