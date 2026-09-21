function outputs = generate_bellhop_pe_equivalent_ray_paths_vertical(overrides)
% Plot only the Bellhop paths matching the PE direct + one-surface topology.
% Uses the current 100 m / Tx 99 m / Rx 3 m representative environment.
if nargin < 1, overrides = struct(); end
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
source_dir = fullfile(root,'results','validation','bellhop_100m_rough_surface','source_1m');
cfg = struct('source_mat',fullfile(source_dir,'source1m_plot_data.mat'), ...
    'output_dir',fullfile(source_dir,'pe_equivalent_paths'), ...
    'show_figures',false,'step_m',0.02,'range_box_m',12);
names = fieldnames(overrides);
for ii = 1:numel(names)
    assert(isfield(cfg,names{ii}),'Unknown override: %s',names{ii});
    cfg.(names{ii}) = overrides.(names{ii});
end
assert(isfile(cfg.source_mat),'Run generate_bellhop_100m_source1m_visuals_vertical first.');
saved = load(cfg.source_mat,'data','cfg');
base_cfg = saved.cfg;
data = saved.data;
assert(base_cfg.water_depth_m == 100 && base_cfg.tx_depth_m == 99 && ...
    base_cfg.rx_depth_m == 3,'Unexpected cached physical geometry.');

direct = data.arrivals(data.arrivals.top_bounces == 0 & ...
    data.arrivals.bottom_bounces == 0,:);
surface = data.arrivals(data.arrivals.top_bounces == 1 & ...
    data.arrivals.bottom_bounces == 0,:);
assert(height(direct) == 1 && height(surface) == 1, ...
    'Expected exactly one direct and one once-surface arrival.');
selected = [direct; surface];
launch_angles = selected.launch_angle_deg.';

if ~isfolder(cfg.output_dir), mkdir(cfg.output_dir); end
base = fullfile(cfg.output_dir,'pe_equivalent_R');
local_run_r(base_cfg,data.surface,base,launch_angles,cfg);
[raw_rays,header] = local_read_rays([base '.ray']);
assert(numel(raw_rays) == 2,'Expected two R trajectories.');
% Bellhop writes R records in increasing-angle order; restore the selected
% arrival order explicitly so path labels cannot be inferred from file order.
rays = raw_rays;
used = false(1,numel(raw_rays));
for ii = 1:numel(launch_angles)
    delta = abs([raw_rays.angle_deg]-launch_angles(ii));
    delta(used) = inf;
    [best,index] = min(delta);
    assert(best < 2e-5,'R launch angle does not match selected A arrival.');
    rays(ii) = raw_rays(index);
    used(index) = true;
end

rx = [base_cfg.rx_range_m;base_cfg.rx_depth_m];
tx = [0;base_cfg.tx_depth_m];
rows = cell(2,1);
for ii = 1:2
    xy = rays(ii).xy - [0;base_cfg.datum_shift_m];
    assert(norm(xy(:,1)-tx) < 1e-6,'Bellhop source changed.');
    [xy,miss] = local_prefix_to_point(xy,rx);
    rays(ii).xy = xy;
    rays(ii).top = selected.top_bounces(ii);
    rays(ii).bottom = selected.bottom_bounces(ii);
    rows{ii} = table(string(local_label(ii)),launch_angles(ii), ...
        selected.receive_angle_deg(ii),selected.delay_s(ii), ...
        selected.amplitude(ii),selected.phase_deg(ii),miss,size(xy,2), ...
        'VariableNames',{'path','launch_angle_deg','receive_angle_deg', ...
        'arrival_delay_s','arrival_amplitude','arrival_phase_deg', ...
        'ray_closest_to_rx_m','plotted_vertices'});
end
summary = vertcat(rows{:});
writetable(summary,fullfile(cfg.output_dir,'pe_equivalent_paths.csv'));
outputs = struct('png',fullfile(cfg.output_dir,'bellhop_pe_direct_once_surface.png'), ...
    'csv',fullfile(cfg.output_dir,'pe_equivalent_paths.csv'), ...
    'mat',fullfile(cfg.output_dir,'pe_equivalent_paths.mat'), ...
    'ray',string([base '.ray']));
save(outputs.mat,'rays','summary','header','base_cfg','cfg','outputs');
local_plot(rays,data.surface,base_cfg,cfg,outputs.png,summary);
fprintf('Saved PE-equivalent Bellhop paths: direct=%d, once-surface=%d, bottom=0.\n', ...
    sum([rays.top] == 0),sum([rays.top] == 1));
disp(summary);
end

function local_run_r(c,surface,base,angles,plot_cfg)
x = surface.x_m;
eta = surface.surface_depth_m;
inside = abs(x) <= plot_cfg.range_box_m + 1e-8;
x = x(inside); eta = eta(inside);
bottom = c.water_depth_m + c.datum_shift_m;
sspmax = bottom + 0.1;
env = sprintf(['''100m water, PE-equivalent direct and once-surface paths''\n%.17g\n1\n''CVW *''\n' ...
    '2 0 %.17g\n0 %.17g 0 1 0 0 /\n%.17g %.17g 0 1 0 0 /\n' ...
    '''A*'' 0\n%.17g %.17g 0 %.17g %.17g 0 /\n' ...
    '1\n%.17g /\n1\n%.17g /\n1\n%.17g /\n''R  X''\n%d\n%s/\n%.17g %.17g %.17g\n'], ...
    c.frequency_hz,sspmax,c.c0_mps,sspmax,c.c0_mps,sspmax,c.bottom_c_mps, ...
    c.bottom_rho_g_cm3,c.bottom_alpha_db_lambda,c.tx_depth_m+c.datum_shift_m, ...
    c.rx_depth_m+c.datum_shift_m,c.rx_range_m/1000,numel(angles), ...
    sprintf('%.17g ',angles),plot_cfg.step_m,sspmax,plot_cfg.range_box_m/1000);
ati = sprintf('''C''\n%d\n%s',numel(x),sprintf('%.17g %.17g\n',[x/1000 eta+c.datum_shift_m].'));
bty = sprintf('''L''\n2\n%.17g %.17g\n%.17g %.17g\n', ...
    -plot_cfg.range_box_m/1000,bottom,plot_cfg.range_box_m/1000,bottom);
local_write([base '.env'],env);
local_write([base '.ati'],ati);
local_write([base '.bty'],bty);
[folder,name] = fileparts(base);
previous = pwd; cleanup = onCleanup(@()cd(previous));
cd(folder);
[status,message] = system(sprintf('"%s" "%s"',c.bellhop_exe,name));
local_write([base '.run.log'],message);
assert(status == 0 && isfile([base '.ray']),'Bellhop R run failed: %s',message);
prt = fileread([base '.prt']);
assert(contains(prt,'CPU Time') && ~contains(prt,'FATAL ERROR'),'Incomplete Bellhop R run.');
end

function [rays,header] = local_read_rays(file)
fid = fopen(file,'r'); assert(fid >= 0,'Cannot open %s.',file);
cleanup = onCleanup(@()fclose(fid));
fgetl(fid); frequency = fscanf(fid,'%f',1);
sources = fscanf(fid,'%d',3); beams = fscanf(fid,'%d',2);
top = fscanf(fid,'%f',1); bottom = fscanf(fid,'%f',1);
fgetl(fid); kind = strtrim(fgetl(fid));
assert(strcmp(kind,'''rz''') && all(sources == 1),'Expected one-source 2-D rays.');
header = struct('frequency_hz',frequency,'beam_counts',beams, ...
    'top_depth_m',top,'bottom_depth_m',bottom);
rays = struct('angle_deg',{},'top',{},'bottom',{},'xy',{});
while true
    angle = fscanf(fid,'%f',1);
    if isempty(angle), break; end
    h = fscanf(fid,'%d',3);
    xy = fscanf(fid,'%f',[2 h(1)]);
    assert(numel(h) == 3 && numel(xy) == 2*h(1),'Incomplete ray record.');
    rays(end+1) = struct('angle_deg',angle,'top',h(2),'bottom',h(3),'xy',xy); %#ok<AGROW>
end
end

function [prefix,miss] = local_prefix_to_point(xy,point)
d = diff(xy,1,2); d2 = sum(d.^2,1);
t = sum((point-xy(:,1:end-1)).*d,1)./max(d2,realmin);
t = max(0,min(1,t)); q = xy(:,1:end-1)+d.*t;
[miss2,ix] = min(sum((q-point).^2,1));
miss = sqrt(miss2);
prefix = [xy(:,1:ix) q(:,ix)];
end

function local_plot(rays,surface,c,plot_cfg,file,summary)
visible = 'off'; if plot_cfg.show_figures, visible = 'on'; end
fig = figure('Visible',visible,'Color','w','Position',[40 40 1500 850]);
cleanup = onCleanup(@()close(fig));
layout = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
colors = [0.05 0.40 0.72;0.90 0.36 0.10];
for pp = 1:2
    ax = nexttile(layout); hold(ax,'on');
    for ii = 1:2
        plot(ax,rays(ii).xy(1,:),rays(ii).xy(2,:),'-','Color',colors(ii,:), ...
            'LineWidth',2.2,'DisplayName',local_label(ii));
        plot(ax,rays(ii).xy(1,end),rays(ii).xy(2,end),'o','Color',colors(ii,:), ...
            'MarkerFaceColor','w','MarkerSize',5,'HandleVisibility','off');
    end
    plot(ax,surface.x_m,surface.surface_depth_m,'k-','LineWidth',1.2, ...
        'DisplayName','Rough pressure-release surface');
    yline(ax,c.water_depth_m,'Color',[0.35 0.24 0.10],'LineWidth',2, ...
        'DisplayName','Seabed');
    plot(ax,0,c.tx_depth_m,'kd','MarkerFaceColor',[1 .75 .1],'MarkerSize',9, ...
        'DisplayName','Tx');
    plot(ax,c.rx_range_m,c.rx_depth_m,'kp','MarkerFaceColor',[.2 .85 .3], ...
        'MarkerSize',13,'DisplayName','Rx');
    set(ax,'YDir','reverse','FontSize',11); grid(ax,'on'); box(ax,'on');
    xlabel(ax,'Horizontal x (m)'); ylabel(ax,'Physical depth z (m)');
    if pp == 1
        xlim(ax,[-0.5 max(1.5,c.rx_range_m+0.5)]); ylim(ax,[-0.2 100.5]);
        title(ax,'Full PE-equivalent geometry'); legend(ax,'Location','eastoutside');
    else
        xlim(ax,[-0.25 max(1.25,c.rx_range_m+0.25)]); ylim(ax,[-0.15 5]);
        title(ax,'Surface reflection and receiver detail');
    end
end
title(layout,{sprintf('Bellhop paths matching PE topology at %.0f Hz',c.frequency_hz), ...
    sprintf('Direct: %.6f ms | once-surface: %.6f ms | no bottom-reflected paths', ...
    1e3*summary.arrival_delay_s(1),1e3*summary.arrival_delay_s(2)), ...
    'A-mode launch angles retraced in R mode; circles are closest traced points to Rx'}, ...
    'FontSize',14,'FontWeight','bold');
exportgraphics(fig,file,'Resolution',170);
end

function label = local_label(ii)
if ii == 1, label = 'Direct (0 surface, 0 bottom)';
else, label = 'Once-surface (1 surface, 0 bottom)'; end
end

function local_write(file,content)
fid = fopen(file,'w','n','UTF-8'); assert(fid >= 0,'Cannot write %s.',file);
cleanup = onCleanup(@()fclose(fid)); fprintf(fid,'%s',content);
end
