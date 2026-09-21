run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'bootstrap_project.m'));
%AUDIT_PE_CACHING_VERTICAL Count and benchmark cacheable PE work.

clear
format compact

cfg = struct();
cfg.f_axis_hz = linspace(4000,8000,32).';
cfg.c0_mps = 1500;
cfg.z_tx_m = 100;
cfg.z_rx_m = 3;
cfg.stepz_lamb = 0.5;
cfg.grid_n = 128;
cfg.aperture_m = 50;
cfg.L_list = [1,16,64,128];
cfg.batch_size = 8;

f = cfg.f_axis_hz;
step_m = cfg.stepz_lamb .* cfg.c0_mps ./ f;
n_tx_surface = ceil(cfg.z_tx_m ./ step_m);
n_tx_rx = ceil((cfg.z_tx_m-cfg.z_rx_m) ./ step_m);
n_surface_rx = ceil(cfg.z_rx_m ./ step_m);

% Current implementation independently computes direct, tx->surface, and
% surface->rx for every realization. Each march uses 2*n_step+2 transforms.
current_march_steps_per_realization = sum(n_tx_rx+n_tx_surface+n_surface_rx);
current_march_fft_per_realization = sum(2*n_tx_rx+2 + ...
    2*n_tx_surface+2 + 2*n_surface_rx+2);
current_kstat_boundary_fft_per_realization = 7*numel(f);

% Optimized design: one tx->surface march captures the direct receiver value
% at z_rx, and one coherent surface->rx march is cached. Each realization
% only marches its random scatter from surface to receiver.
cached_build_steps = sum(n_tx_surface+n_surface_rx);
cached_build_march_fft = sum(2*n_tx_surface+2 + 2*n_surface_rx+2);
joint_model_fft = numel(f)*(numel(f)+1); % C and P spectra for i<=j.
cached_sample_steps = sum(n_surface_rx);
cached_sample_march_fft = sum(2*n_surface_rx+2);
cached_sample_generator_ifft = numel(f);

rows = struct([]);
for ii=1:numel(cfg.L_list)
    L=cfg.L_list(ii);
    row = struct();
    row.L = L;
    row.current_march_steps = L*current_march_steps_per_realization;
    row.cached_march_steps = cached_build_steps + L*cached_sample_steps;
    row.step_reduction_factor = row.current_march_steps/row.cached_march_steps;
    row.current_fft_count = L*(current_march_fft_per_realization + ...
        current_kstat_boundary_fft_per_realization);
    row.cached_fft_count = cached_build_march_fft + joint_model_fft + ...
        L*(cached_sample_march_fft+cached_sample_generator_ifft);
    row.fft_reduction_factor = row.current_fft_count/row.cached_fft_count;
    if isempty(rows), rows=row; else, rows(end+1)=row; end %#ok<SAGROW>
end
cost_table = struct2table(rows);

% Measured primitive transform time on the target audit grid.
rng(60101,'twister');
field = randn(cfg.grid_n)+1i*randn(cfg.grid_n);
n_repeat = 100;
tic_fft = tic;
for ii=1:n_repeat, fft2(field); end
fft2_time_s = toc(tic_fft)/n_repeat;
tic_ifft = tic;
for ii=1:n_repeat, ifft2(field); end
ifft2_time_s = toc(tic_ifft)/n_repeat;
transform_time_s = 0.5*(fft2_time_s+ifft2_time_s);
cost_table.current_fft_time_projection_s = cost_table.current_fft_count*transform_time_s;
cost_table.cached_fft_time_projection_s = cost_table.cached_fft_count*transform_time_s;

% Small end-to-end baseline timing through the unmodified public channel.
benchmark = struct();
benchmark.grid_n = 128;
benchmark.f_axis_hz = linspace(4000,8000,4);
benchmark.repeat_count = 2;
benchmark.elapsed_s = zeros(benchmark.repeat_count,1);
paramsV = local_benchmark_params(benchmark);
for ii=1:benchmark.repeat_count
    paramsV.sea_seed = 62000+ii;
    t=tic; vertical_channel_model(paramsV); benchmark.elapsed_s(ii)=toc(t);
end
benchmark.mean_elapsed_s = mean(benchmark.elapsed_s);
benchmark.mean_per_frequency_s = benchmark.mean_elapsed_s/numel(benchmark.f_axis_hz);
benchmark.note = ['Measured current public path only. Cached runtime is a transparent ', ...
    'FFT-count projection because the cache execution path is not added to the public propagator in this stage.'];

nxy = cfg.grid_n^2;
memory = struct();
memory.incident_fields_bytes = nxy*numel(f)*16;
memory.propagation_kernels_bytes = nxy*numel(f)*16;
memory.static_medium_screen_bytes = nxy*numel(f)*16;
memory.direct_and_coherent_response_bytes = 2*numel(f)*16;
memory.full_joint_batch_L64_bytes = nxy*numel(f)*64*16;
memory.streamed_joint_batch_bytes = nxy*numel(f)*cfg.batch_size*16;
memory.total_recommended_cache_bytes = memory.incident_fields_bytes + ...
    memory.propagation_kernels_bytes + memory.static_medium_screen_bytes + ...
    memory.direct_and_coherent_response_bytes + memory.streamed_joint_batch_bytes;
memory.note = ['Use realization batches; do not retain all L boundary fields. ', ...
    'Counts are complex double and exclude MATLAB array headers.'];

cacheability = table( ...
    ["tx_to_surface_incident";"direct_channel";"coherent_reflection"; ...
     "propagation_kernel";"grid_medium_terms";"random_surface_to_rx"], ...
    [true;true;true;true;true;false], ...
    ["fixed env/f/source";"extract during cached tx-to-surface march"; ...
     "fixed env/f/U";"fixed env/f/grid/dz";"fixed env/f/grid"; ...
     "changes for every realization"], ...
    'VariableNames',{'quantity','cacheable','condition'});

result_file=project_result_file('validation','audit_pe_caching_vertical_result.mat');
csv_file=project_result_file('validation','audit_pe_caching_vertical_cost_table.csv');
save(result_file,'cfg','n_tx_surface','n_tx_rx','n_surface_rx','cost_table', ...
    'benchmark','memory','cacheability');
writetable(cost_table,csv_file);
disp(cost_table); disp(benchmark); disp(memory); disp(cacheability);
fprintf('Saved %s\n',result_file);

function paramsV=local_benchmark_params(benchmark)
paramsV=struct();
paramsV.f0=benchmark.f_axis_hz;
paramsV.enable_wideband=false;
paramsV.c0=1500; paramsV.z_max=100; paramsV.stepz_lamb=0.5;
paramsV.xw=50; paramsV.yw=50; paramsV.nx=benchmark.grid_n; paramsV.ny=benchmark.grid_n;
paramsV.x_tx=0; paramsV.y_tx=0; paramsV.z_tx=100;
paramsV.x_rx=0; paramsV.y_rx=0; paramsV.z_rx=3;
paramsV.nout=2; paramsV.sigma_src_m=0.4; paramsV.sponge_ratio=0.12;
paramsV.alpha_max_np_per_m=0.15; paramsV.env_mode='uniform';
paramsV.enforce_1_over_R=false; paramsV.show_figures=false; paramsV.save_mode='rx_only';
paramsV.use_gpu=false; paramsV.enable_surface_reflection=true;
paramsV.sea_wind_speed=5; paramsV.sea_hs_target=0.5;
paramsV.surface_roughness_scale_mode='raw_pm'; paramsV.sea_seed=62001;
paramsV.surface_reflect_coeff=-1; paramsV.surface_phase_mode='normal';
paramsV.surface_boundary_model='kirchhoff_kstat'; paramsV.surface_kstat_random_scatter=true;
end
