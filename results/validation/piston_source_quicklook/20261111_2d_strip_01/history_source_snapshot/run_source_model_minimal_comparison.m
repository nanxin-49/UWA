function results = run_source_model_minimal_comparison(run_id)
%RUN_SOURCE_MODEL_MINIMAL_COMPARISON Normalized, solver-free source comparison.
% Run after setup_vertical_project, from any working directory.
if nargin < 1 || isempty(run_id)
    run_id = datestr(now,'yyyymmdd_HHMMSS');
end
assert(ischar(run_id) || isstring(run_id));
run_id = char(run_id);
assert(~isempty(regexp(run_id,'^[A-Za-z0-9_]+$','once')),'Use an alphanumeric run_id.');
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
out = fullfile(root,'results','comparisons','source_model_minimal',run_id);
assert(~exist(out,'dir'),'Run directory already exists.');
mkdir(out);
diary_path = fullfile(out,'matlab_diary.txt');
diary(diary_path); diary_cleanup = onCleanup(@()diary('off')); %#ok<NASGU>
t_start = tic;
try
    cfg = struct('frequency_hz',4000,'c_mps',1500,'water_depth_m',1500, ...
        'z_tx_m',1490,'z_rx_m',6,'sigma_m',0.3,'piston_radius_m',0.05);
    cfg.delta_z_m = cfg.z_tx_m-cfg.z_rx_m;
    cfg.k_radpm = 2*pi*cfg.frequency_hz/cfg.c_mps;
    x_m = (0:1:200).';
    theta_deg = (0:0.1:90).';
    theta_rad = theta_deg*pi/180;
    theta_rx_rad = atan2(x_m,cfg.delta_z_m);
    theta_max_deg = theta_rx_rad(end)*180/pi;
    D_gaussian = exp(-0.5*(cfg.k_radpm*cfg.sigma_m*sin(theta_rad)).^2);
    D_point = ones(size(theta_rad));
    D_piston = local_piston(cfg.k_radpm*cfg.piston_radius_m*sin(theta_rad));
    D_gaussian_rx = exp(-0.5*(cfg.k_radpm*cfg.sigma_m*sin(theta_rx_rad)).^2);
    D_piston_rx = local_piston(cfg.k_radpm*cfg.piston_radius_m*sin(theta_rx_rad));
    R_m = hypot(x_m,cfg.delta_z_m);
    geometric_factor = cfg.delta_z_m./R_m;
    receiver_gaussian = D_gaussian_rx.*geometric_factor;
    receiver_point = geometric_factor;
    receiver_piston = D_piston_rx.*geometric_factor;
    relative_gaussian_db = 20*log10(D_gaussian_rx);
    relative_piston_db = 20*log10(D_piston_rx);
    checks = struct;
    checks.point_constant = all(D_point == 1);
    checks.piston_axis_finite_one = isfinite(D_piston(1)) && D_piston(1) == 1;
    checks.common_geometry = all(R_m >= cfg.delta_z_m) && cfg.delta_z_m == 1484;
    checks.common_axis_reference = all(abs([receiver_gaussian(1),receiver_point(1),receiver_piston(1)]-1)<1e-14);
    checks.physical_order = all(receiver_gaussian <= receiver_piston+1e-14) && ...
        all(receiver_piston <= receiver_point+1e-14);
    % Numerically transform the production 1-D slice. The separable y integral
    % cancels in the axis normalization of its 2-D Fourier transform.
    x_source = (-2:0.002:2).';
    psi = gaussian_source_initial_field_vertical(x_source,zeros(size(x_source)), ...
        struct('x_tx',0,'y_tx',0,'sigma_src_m',cfg.sigma_m));
    theta_check = [0;5;theta_max_deg;15;30]*pi/180;
    q = cfg.k_radpm*sin(theta_check);
    numerical_spectrum = zeros(size(q));
    for ii = 1:numel(q)
        numerical_spectrum(ii) = trapz(x_source,psi.*exp(-1i*q(ii)*x_source));
    end
    numerical_spectrum = numerical_spectrum/numerical_spectrum(1);
    analytic_spectrum = exp(-0.5*(cfg.sigma_m*q).^2);
    spectrum_error = max(abs(numerical_spectrum-analytic_spectrum));
    checks.gaussian_transform = spectrum_error < 1e-8;
    check_values = struct2cell(checks);
    failed_count = sum(~cell2mat(check_values));
    assert(failed_count == 0,'One or more model checks failed.');
    sample_ix = [1,51,101,151,201];
    sample_table = table(x_m(sample_ix),theta_rx_rad(sample_ix)*180/pi, ...
        receiver_gaussian(sample_ix),receiver_point(sample_ix),receiver_piston(sample_ix), ...
        relative_gaussian_db(sample_ix),relative_piston_db(sample_ix), ...
        'VariableNames',{'x_m','theta_deg','gaussian_rel_axis','point_rel_axis', ...
        'piston_rel_axis','gaussian_vs_point_db','piston_vs_point_db'});
    full_table = table(x_m,theta_rx_rad*180/pi,receiver_gaussian,receiver_point, ...
        receiver_piston,relative_gaussian_db,relative_piston_db, ...
        'VariableNames',sample_table.Properties.VariableNames);
    results = struct('run_id',run_id,'config',cfg,'theta_max_deg',theta_max_deg, ...
        'spectrum_error_max',spectrum_error,'checks',checks,'failed_count',failed_count, ...
        'sample_table',sample_table,'full_table',full_table,'out_dir',out);
    writetable(sample_table,fullfile(out,'sample_positions.csv'));
    writetable(full_table,fullfile(out,'receiver_curves.csv'));
    save(fullfile(out,'source_comparison.mat'),'results','theta_deg', ...
        'D_gaussian','D_point','D_piston','-v7');
    f1 = figure('Visible','off','Color','w','Position',[100 100 920 540]);
    ax1 = axes(f1); hold(ax1,'on');
    plot(ax1,theta_deg,20*log10(max(D_gaussian,realmin)),'-','LineWidth',2,'DisplayName','Gaussian initial-plane spectrum');
    plot(ax1,theta_deg,20*log10(D_point),'--','LineWidth',1.8,'DisplayName','Omnidirectional point');
    plot(ax1,theta_deg,20*log10(D_piston),'-.','LineWidth',1.8,'DisplayName','Uniform circular piston');
    geometry_line = xline(ax1,theta_max_deg,':','Geometry max','LineWidth',1.5);
    geometry_line.HandleVisibility = 'off';
    xlabel(ax1,'Angle from upward vertical (deg)'); ylabel(ax1,'Axis-normalized pressure (dB)');
    title(ax1,'Emission angle response: Gaussian spectrum vs physical piston pattern');
    xlim(ax1,[0 90]); ylim(ax1,[-30 1]); grid(ax1,'on'); legend(ax1,'Location','southwest');
    local_save_figure(f1,out,'emission_pattern');
    f2 = figure('Visible','off','Color','w','Position',[100 100 920 680]);
    layout = tiledlayout(f2,2,1,'TileSpacing','compact');
    ax2 = nexttile(layout); hold(ax2,'on');
    plot(ax2,x_m,20*log10(receiver_gaussian),'-','LineWidth',2,'DisplayName','Gaussian');
    plot(ax2,x_m,20*log10(receiver_point),'--','LineWidth',1.8,'DisplayName','Point');
    plot(ax2,x_m,20*log10(receiver_piston),'-.','LineWidth',1.8,'DisplayName','Piston');
    ylabel(ax2,'Receiver amplitude vs axis (dB)'); title(ax2,'Common 1/R geometry, each axis set to 0 dB');
    grid(ax2,'on'); legend(ax2,'Location','southwest');
    ax3 = nexttile(layout); hold(ax3,'on');
    plot(ax3,x_m,relative_gaussian_db,'-','LineWidth',2,'DisplayName','Gaussian / point');
    plot(ax3,x_m,relative_piston_db,'-.','LineWidth',1.8,'DisplayName','Piston / point');
    xlabel(ax3,'Horizontal receiver offset x (m)'); ylabel(ax3,'Amplitude difference (dB)');
    xlim(ax2,[0 200]); xlim(ax3,[0 200]); grid(ax3,'on'); legend(ax3,'Location','southwest');
    local_save_figure(f2,out,'receiver_amplitude');
    results.elapsed_s = toc(t_start);
    save(fullfile(out,'source_comparison.mat'),'results','theta_deg', ...
        'D_gaussian','D_point','D_piston','-v7');
    fprintf('RUN_ID=%s STATUS=PASS ELAPSED_S=%.3f FAILED_CHECKS=%d\n',run_id,results.elapsed_s,failed_count);
    fprintf('THETA_MAX_DEG=%.9f GAUSSIAN_MINUS_POINT_DB_200=%.9f PISTON_MINUS_POINT_DB_200=%.9f\n', ...
        theta_max_deg,relative_gaussian_db(end),relative_piston_db(end));
    disp(sample_table);
catch ME
    fprintf(2,'STATUS=FAIL ERROR=%s\n',ME.message);
    rethrow(ME);
end
end

function D = local_piston(u)
D = ones(size(u));
nonzero = u ~= 0;
D(nonzero) = 2*besselj(1,u(nonzero))./u(nonzero);
end

function local_save_figure(fig,out,stem)
drawnow;
fig_path = fullfile(out,[stem '.fig']);
png_path = fullfile(out,[stem '.png']);
assert(~exist(fig_path,'file') && ~exist(png_path,'file'));
savefig(fig,fig_path);
exportgraphics(fig,png_path,'Resolution',180);
close(fig);
end
