function result = render_target_deepsea_flat_geometry(run_id)
%RENDER_TARGET_DEEPSEA_FLAT_GEOMETRY Analytic image-source geometry only.
% No PE/WAPE, boundary solver, or measured field is used by this workflow.
if nargin < 1 || isempty(run_id)
    run_id = char(datetime('now','Format','yyyyMMdd_HHmmss'));
end
run_id = char(run_id);
assert(~isempty(regexp(run_id,'^[A-Za-z0-9_]+$','once')), ...
    'run_id must contain only letters, digits, and underscores.');
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
out_dir = fullfile(root,'results','visualization', ...
    'rc1_target_deepsea_flat_geometry',run_id);
assert(~isfolder(out_dir),'Result directory already exists.');

cfg = struct('water_depth_m',1500,'z_tx_m',1490, ...
    'z_rx_m',[3 6 10],'c_mps',1500,'x_m',(0:200).', ...
    'surface','flat_pressure_release', ...
    'model','uniform_sound_speed_image_source_geometry');
x_m = cfg.x_m;
z_rx_m = cfg.z_rx_m;
n_x = numel(x_m);
n_z = numel(z_rx_m);
L_D_m = zeros(n_x,n_z);
L_S_m = zeros(n_x,n_z);
tau_D_s = zeros(n_x,n_z);
tau_S_s = zeros(n_x,n_z);
delta_tau_s = zeros(n_x,n_z);
delta_f_hz = zeros(n_x,n_z);
theta_v_deg = zeros(n_x,n_z);
for j = 1:n_z
    z_rx = z_rx_m(j);
    L_D_m(:,j) = hypot(x_m,cfg.z_tx_m-z_rx);
    L_S_m(:,j) = hypot(x_m,cfg.z_tx_m+z_rx);
    tau_D_s(:,j) = L_D_m(:,j)/cfg.c_mps;
    tau_S_s(:,j) = L_S_m(:,j)/cfg.c_mps;
    delta_tau_s(:,j) = tau_S_s(:,j)-tau_D_s(:,j);
    delta_f_hz(:,j) = 1./delta_tau_s(:,j);
    theta_v_deg(:,j) = atan2d(x_m,cfg.z_tx_m-z_rx);
end
assert(all(isfinite(delta_f_hz(:))) && all(delta_tau_s(:)>0));
assert(max(abs(delta_tau_s(1,:)-2*z_rx_m/cfg.c_mps))<1e-12);
assert(all(diff(theta_v_deg,1,1)>=0,'all'));

rows = n_x*n_z;
data = table(repmat(x_m,n_z,1), ...
    repelem(z_rx_m(:),n_x),repmat(cfg.z_tx_m,rows,1), ...
    reshape(L_D_m,[],1),reshape(L_S_m,[],1), ...
    reshape(tau_D_s,[],1),reshape(tau_S_s,[],1), ...
    reshape(delta_tau_s,[],1),reshape(delta_f_hz,[],1), ...
    reshape(theta_v_deg,[],1), ...
    'VariableNames',{'x_m','z_rx_m','z_tx_m','L_D_m','L_S_m', ...
    'tau_D_s','tau_S_s','delta_tau_s','delta_f_hz','theta_v_deg'});
key_mask = ismember(data.x_m,[0 100 150 200]);
key_data = data(key_mask,:);
mkdir(out_dir);
writetable(data,fullfile(out_dir,'geometry_samples.csv'));
writetable(key_data,fullfile(out_dir,'geometry_key_points.csv'));
fid = fopen(fullfile(out_dir,'config.json'),'w');
assert(fid>0,'Could not open config.json for writing.');
file_guard = onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(rmfield(cfg,'x_m')));
clear file_guard

colors = [0.00 0.40 0.68; 0.84 0.36 0.08; 0.30 0.55 0.28];
labels = {'Rx 3 m','Rx 6 m','Rx 10 m'};
styles = {'-','--','-.'};
caption = 'Uniform-c, flat pressure-release surface: analytic geometry (not PE or measured field)';

% Figure 1: representative 150 m offset, 6 m receiver, with a surface inset.
x_rx = 150;
z_rx = 6;
x_surface = x_rx*cfg.z_tx_m/(cfg.z_tx_m+z_rx);
fig = local_figure();
t = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(t,1); hold(ax,'on');
plot(ax,[0 205],[0 0],'k-','LineWidth',1.4);
plot(ax,[0 205],[1500 1500],'k:','LineWidth',1.2);
plot(ax,[0 x_rx],[cfg.z_tx_m z_rx],'-','Color',colors(1,:),'LineWidth',2);
plot(ax,[0 x_surface x_rx],[cfg.z_tx_m 0 z_rx],'--', ...
    'Color',colors(2,:),'LineWidth',2);
plot(ax,0,cfg.z_tx_m,'ko','MarkerFaceColor','k');
plot(ax,x_rx,z_rx,'ks','MarkerFaceColor','k');
plot(ax,x_surface,0,'kd','MarkerFaceColor','k');
text(ax,3,1450,'Tx: (0,1490 m)');
text(ax,155,75,'Rx: (150,6 m)');
xlim(ax,[-5 205]); ylim(ax,[0 1500]); set(ax,'YDir','reverse');
xlabel(ax,'Horizontal offset x (m)'); ylabel(ax,'Depth z (m)');
title(ax,'Full water column'); grid(ax,'on');
ax = nexttile(t,2); hold(ax,'on');
plot(ax,[140 155],[0 0],'k-','LineWidth',1.4);
plot(ax,[140 x_rx],[cfg.z_tx_m-(cfg.z_tx_m-z_rx)*140/x_rx z_rx], ...
    '-','Color',colors(1,:),'LineWidth',2);
plot(ax,[140 x_surface x_rx], ...
    [cfg.z_tx_m*(1-140/x_surface) 0 z_rx],'--', ...
    'Color',colors(2,:),'LineWidth',2);
plot(ax,x_surface,0,'kd','MarkerFaceColor','k');
plot(ax,x_rx,z_rx,'ks','MarkerFaceColor','k');
xlim(ax,[140 155]); ylim(ax,[0 100]); set(ax,'YDir','reverse');
xlabel(ax,'Horizontal offset x (m)'); ylabel(ax,'Depth z (m)');
title(ax,'Near-surface detail'); grid(ax,'on');
title(t,caption);
local_save(fig,out_dir,'F01_geometry_paths');

% Figure 2: absolute geometric travel times; depth and path use redundant cues.
fig = local_figure();
t = tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
for j = 1:n_z
    ax = nexttile(t,j); hold(ax,'on');
    plot(ax,x_m,1e3*tau_D_s(:,j),'-','Color',colors(j,:), ...
        'LineWidth',1.8,'DisplayName','Direct');
    plot(ax,x_m,1e3*tau_S_s(:,j),'--','Color',colors(j,:), ...
        'LineWidth',1.8,'DisplayName','Surface image');
    xlabel(ax,'Horizontal offset x (m)'); ylabel(ax,'Travel time (ms)');
    title(ax,labels{j}); xlim(ax,[0 200]); grid(ax,'on');
    legend(ax,'Location','northwest');
end
title(t,caption);
local_save(fig,out_dir,'F02_path_travel_times');

local_curves(out_dir,'F03_delay_difference',x_m,1e3*delta_tau_s, ...
    z_rx_m,labels,styles,colors,'Direct-to-surface delay difference (ms)',caption);
local_curves(out_dir,'F04_interference_spacing',x_m,delta_f_hz, ...
    z_rx_m,labels,styles,colors,'1 / delay difference (Hz)',caption);
local_curves(out_dir,'F05_vertical_deviation_angle',x_m,theta_v_deg, ...
    z_rx_m,labels,styles,colors,'Angle from vertical upward axis (deg)',caption);

result = struct('run_id',run_id,'output_dir',out_dir, ...
    'sample_count',height(data),'figure_count',5, ...
    'failed_checks',0,'config',cfg);
fprintf('RC1_ANALYTIC_GEOMETRY_PASS run=%s rows=%d figures=%d output=%s\n', ...
    run_id,height(data),result.figure_count,out_dir);
end

function fig = local_figure()
fig = figure('Visible','off','Color','w','Position',[80 80 1200 520]);
end

function local_curves(out_dir,stem,x_m,y,z_rx_m,labels,styles,colors,y_label,caption)
fig = local_figure();
ax = axes(fig); hold(ax,'on');
for j = 1:numel(z_rx_m)
    plot(ax,x_m,y(:,j),styles{j},'Color',colors(j,:), ...
        'LineWidth',2,'DisplayName',labels{j});
end
xlabel(ax,'Horizontal offset x (m)'); ylabel(ax,y_label);
title(ax,caption); xlim(ax,[0 200]); grid(ax,'on');
legend(ax,'Location','best');
local_save(fig,out_dir,stem);
end

function local_save(fig,out_dir,stem)
fig_path = fullfile(out_dir,[stem '.fig']);
png_path = fullfile(out_dir,[stem '.png']);
assert(~isfile(fig_path) && ~isfile(png_path),'Refusing to overwrite a figure.');
drawnow;
savefig(fig,fig_path);
exportgraphics(fig,png_path,'Resolution',180);
close(fig);
end
