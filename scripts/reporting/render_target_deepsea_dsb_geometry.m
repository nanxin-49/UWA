function result = render_target_deepsea_dsb_geometry(run_id)
%RENDER_TARGET_DEEPSEA_DSB_GEOMETRY Analytic D/S/B mirror paths only.
% No pressure amplitudes, boundary solver, or PE/WAPE/BIE/FullFK calculations.
if nargin < 1 || isempty(run_id)
    run_id = char(datetime('now','Format','yyyyMMdd_HHmmss'));
end
run_id = char(run_id);
assert(~isempty(regexp(run_id,'^[A-Za-z0-9_]+$','once')), ...
    'run_id must contain only letters, digits, and underscores.');
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
out_dir = fullfile(root,'results','visualization', ...
    'rc1_target_deepsea_dsb_geometry',run_id);
assert(~isfolder(out_dir),'Result directory already exists.');

cfg = struct('water_depth_m',1500,'z_tx_m',1490, ...
    'z_rx_m',[3 6 10],'c_mps',1500,'x_m',(0:200).', ...
    'surface','flat_pressure_release', ...
    'bottom','flat_single_specular_path_geometry_only');
x_m = cfg.x_m;
z_rx_m = cfg.z_rx_m;
n_x = numel(x_m);
n_z = numel(z_rx_m);
L_D_m = zeros(n_x,n_z);
L_S_m = zeros(n_x,n_z);
L_B_m = zeros(n_x,n_z);
for j = 1:n_z
    z_rx = z_rx_m(j);
    L_D_m(:,j) = hypot(x_m,cfg.z_tx_m-z_rx);
    L_S_m(:,j) = hypot(x_m,cfg.z_tx_m+z_rx);
    % Mirror the source across z=water_depth_m: z_image=2H-z_tx.
    L_B_m(:,j) = hypot(x_m,2*cfg.water_depth_m-cfg.z_tx_m-z_rx);
end
tau_D_s = L_D_m/cfg.c_mps;
tau_S_s = L_S_m/cfg.c_mps;
tau_B_s = L_B_m/cfg.c_mps;
delta_S_D_s = tau_S_s-tau_D_s;
delta_B_D_s = tau_B_s-tau_D_s;
delta_B_S_s = tau_B_s-tau_S_s;
assert(all(delta_S_D_s(:)>0) && all(delta_B_D_s(:)>0));
assert(all(delta_B_S_s(:,1:2)>0,'all'));
assert(max(abs(delta_B_S_s(:,3)))<1e-14, ...
    'The 10 m receiver S/B geometry must coincide.');
assert(max(abs(delta_B_D_s(1,:)-20/cfg.c_mps))<1e-14);

rows = n_x*n_z;
data = table(repmat(x_m,n_z,1),repelem(z_rx_m(:),n_x), ...
    repmat(cfg.z_tx_m,rows,1),repmat(cfg.water_depth_m,rows,1), ...
    reshape(L_D_m,[],1),reshape(L_S_m,[],1),reshape(L_B_m,[],1), ...
    reshape(tau_D_s,[],1),reshape(tau_S_s,[],1),reshape(tau_B_s,[],1), ...
    reshape(delta_S_D_s,[],1),reshape(delta_B_D_s,[],1), ...
    reshape(delta_B_S_s,[],1), ...
    'VariableNames',{'x_m','z_rx_m','z_tx_m','water_depth_m', ...
    'L_D_m','L_S_m','L_B_m','tau_D_s','tau_S_s','tau_B_s', ...
    'delta_S_D_s','delta_B_D_s','delta_B_S_s'});
key_data = data(ismember(data.x_m,[0 100 150 200]),:);
mkdir(out_dir);
writetable(data,fullfile(out_dir,'dsb_geometry_samples.csv'));
writetable(key_data,fullfile(out_dir,'dsb_geometry_key_points.csv'));
fid = fopen(fullfile(out_dir,'config.json'),'w');
assert(fid>0,'Could not open config.json for writing.');
file_guard = onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(rmfield(cfg,'x_m')));
clear file_guard

fig = figure('Visible','off','Color','w','Position',[70 70 1420 530]);
t = tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
for j = 1:n_z
    ax = nexttile(t,j); hold(ax,'on');
    plot(ax,x_m,zeros(n_x,1),'-','Color',[0.35 0.35 0.35], ...
        'LineWidth',1.8,'DisplayName','D: direct');
    plot(ax,x_m,1e3*delta_B_D_s(:,j),'-','Color',[0.82 0.35 0.07], ...
        'LineWidth',2.0,'DisplayName','B: one flat-bottom image');
    plot(ax,x_m,1e3*delta_S_D_s(:,j),'--','Color',[0.00 0.38 0.70], ...
        'LineWidth',2.2,'DisplayName','S: one flat-surface image');
    if j == 3
        plot(ax,x_m(1:20:end),1e3*delta_B_D_s(1:20:end,j),'o', ...
            'Color',[0.82 0.35 0.07],'MarkerFaceColor','w', ...
            'MarkerSize',5,'HandleVisibility','off');
        text(ax,7,10.7,'S and B coincide geometrically', ...
            'Color',[0.20 0.20 0.20],'FontSize',10);
    end
    xlim(ax,[0 200]); ylim(ax,[0 14.5]); yticks(ax,0:2:14);
    xlabel(ax,'Horizontal offset x (m)');
    ylabel(ax,'Arrival time relative to D (ms)');
    title(ax,sprintf('Receiver depth %g m',z_rx_m(j)));
    grid(ax,'on');
    if j == 1
        legend(ax,'Location','east','FontSize',9);
    end
end
title(t,['Uniform-c, flat surface and bottom: analytic D/S/B path times ', ...
    '(not PE or measured pressure)']);
fig_path = fullfile(out_dir,'F01_DSB_relative_delays.fig');
png_path = fullfile(out_dir,'F01_DSB_relative_delays.png');
assert(~isfile(fig_path) && ~isfile(png_path), ...
    'Refusing to overwrite the figure.');
drawnow;
savefig(fig,fig_path);
exportgraphics(fig,png_path,'Resolution',220);
close(fig);

result = struct('run_id',run_id,'output_dir',out_dir, ...
    'sample_count',height(data),'figure_count',1,'failed_checks',0, ...
    'max_abs_S_B_delay_at_10m_s',max(abs(delta_B_S_s(:,3))), ...
    'config',cfg);
fprintf(['RC1_DSB_ANALYTIC_PASS run=%s rows=%d ', ...
    'max_abs_S_B_delay_at_10m_s=%.3g output=%s\n'], ...
    run_id,height(data),result.max_abs_S_B_delay_at_10m_s,out_dir);
end
