function output = render_source_model_2d_beam_maps(run_id)
%RENDER_SOURCE_MODEL_2D_BEAM_MAPS Project analytic angle weights into 2-D.
% Colors show axis-normalized directional pressure, with spreading removed.
if nargin < 1 || isempty(run_id)
    run_id = char(datetime('now','Format','yyyyMMdd_HHmmss'));
end
run_id = char(run_id);
assert(~isempty(regexp(run_id,'^[A-Za-z0-9_]+$','once')),'Invalid run_id.');
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
source_file = fullfile(root,'results','comparisons','source_model_minimal', ...
    'run_20261030_02','source_comparison.mat');
piston_file = fullfile(root,'results','comparisons','source_model_minimal', ...
    'piston_frequency_20261030_01','piston_frequency_data.mat');
assert(isfile(source_file) && isfile(piston_file),'Required comparison data missing.');
s = load(source_file,'results');
p = load(piston_file,'result');
cfg = s.results.config;
assert(cfg.frequency_hz==4000 && cfg.c_mps==1500 && cfg.sigma_m==0.3);
assert(cfg.piston_radius_m==0.05 && cfg.delta_z_m==1484);
assert(p.result.config.radius_m==cfg.piston_radius_m && ...
    p.result.config.c_mps==cfg.c_mps && ...
    isequal(p.result.config.frequencies_hz,[4000 8000 12000]));
out_dir = fullfile(root,'results','visualization','source_model_2d_beams',run_id);
assert(~isfolder(out_dir),'Refusing to overwrite an existing result directory.');
mkdir(out_dir);
L = cfg.delta_z_m;
theta_max_deg = atan2(200,L)*180/pi;
assert(abs(theta_max_deg-7.6756)<1e-4);
% Full fan shows beam shape; lower row uses the actual +/-200 m link window.
x_full = linspace(-L,L,401);
x_link = linspace(-200,200,401);
s_up = linspace(0,L,301);
source_titles = {'Gaussian initial-plane spectrum','Omnidirectional point', ...
    'Uniform circular piston, 4 kHz'};
freq_titles = {'Piston 4 kHz','Piston 8 kHz','Piston 12 kHz'};
source_fn = {@(th)local_gaussian(th,cfg.frequency_hz,cfg.c_mps,cfg.sigma_m), ...
    @(th)ones(size(th)), ...
    @(th)local_piston(th,cfg.frequency_hz,cfg.c_mps,cfg.piston_radius_m)};
freqs = p.result.config.frequencies_hz;
freq_fn = cell(1,3);
for ii = 1:3
    f_hz = freqs(ii);
    freq_fn{ii} = @(th)local_piston(th,f_hz,cfg.c_mps,cfg.piston_radius_m);
end
source_paths = local_render_set(out_dir,'three_sources_2d',source_titles,source_fn, ...
    x_full,x_link,s_up,[-30 0],[-2.2 0],L,theta_max_deg);
frequency_paths = local_render_set(out_dir,'piston_frequencies_2d',freq_titles,freq_fn, ...
    x_full,x_link,s_up,[-9 0],[-0.16 0],L,theta_max_deg);
source_at_max = [source_fn{1}(theta_max_deg*pi/180),1, ...
    source_fn{3}(theta_max_deg*pi/180)];
piston_at_max = zeros(1,3);
for ii = 1:3
    piston_at_max(ii) = freq_fn{ii}(theta_max_deg*pi/180);
end
assert(max(abs(20*log10(piston_at_max(:))-p.result.summary.level_at_7p6756_deg_db))<1e-5);
output = struct('run_id',run_id,'source_data',{{source_file,piston_file}}, ...
    'output_dir',out_dir,'theta_max_deg',theta_max_deg, ...
    'source_at_max_db',20*log10(source_at_max), ...
    'piston_at_max_db',20*log10(piston_at_max), ...
    'full_fan_x_m',[x_full(1) x_full(end)],'link_x_m',[x_link(1) x_link(end)], ...
    'upward_distance_m',[s_up(1) s_up(end)], ...
    'source_figure',{source_paths},'frequency_figure',{frequency_paths}, ...
    'meaning','Axis-normalized angular pressure weight mapped to (x,upward range); no geometric spreading or solved field.');
save(fullfile(out_dir,'visualization_meta.mat'),'output','-v7');
fprintf('RUN_ID=%s STATUS=PASS THETA_MAX_DEG=%.6f\n',run_id,theta_max_deg);
fprintf('SOURCE_DB_AT_MAX=%.6f,%.6f,%.6f\n',output.source_at_max_db);
fprintf('PISTON_DB_AT_MAX=%.6f,%.6f,%.6f\n',output.piston_at_max_db);
end

function paths = local_render_set(out_dir,stem,titles,functions,x_full,x_link,s_up,full_clim,link_clim,L,theta_max_deg)
fig = figure('Visible','off','Color','w','Position',[40 40 1500 900]);
t = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
for row = 1:2
    if row==1, x_m=x_full; limits_db=full_clim; else, x_m=x_link; limits_db=link_clim; end
    [X,S] = meshgrid(x_m,s_up);
    theta = atan2(abs(X),S);
    for col = 1:3
        D = functions{col}(theta);
        assert(all(isfinite(D(:))) && all(D(:)>0) && all(D(:)<=1+1e-12));
        gain_db = 20*log10(D);
        ax = nexttile(t,(row-1)*3+col);
        img = imagesc(ax,x_m,s_up,gain_db);
        set(ax,'YDir','normal'); axis(ax,'tight');
        clim(ax,limits_db); colormap(ax,parula(256)); hold(ax,'on');
        if row==1
            plot(ax,[-200 0 200],[L 0 L],'w:','LineWidth',1.2,'HandleVisibility','off');
            plot(ax,[-200 200],[L L],'wo','MarkerSize',4,'HandleVisibility','off');
            subtitle = 'wide-angle view';
        else
            link_mask = abs(X) <= (200/L).*S;
            set(img,'AlphaData',double(link_mask));
            set(ax,'Color',[0.93 0.93 0.93]);
            plot(ax,[-200 0 200],[L 0 L],'k:','LineWidth',1,'HandleVisibility','off');
            plot(ax,[-200 200],[L L],'ko','MarkerFaceColor','w','MarkerSize',4,'HandleVisibility','off');
            subtitle = 'current link-angle cone';
        end
        plot(ax,0,0,'wp','MarkerFaceColor','w','MarkerSize',9,'HandleVisibility','off');
        title(ax,sprintf('%s | %s',titles{col},subtitle),'FontSize',10);
        if row==2, xlabel(ax,'Horizontal offset x (m)'); end
        if col==1, ylabel(ax,'Upward distance from Tx (m)'); end
        if col==3
            cb=colorbar(ax); cb.Label.String='Axis-normalized directional pressure (dB)';
        end
    end
end
title(t,sprintf('Analytic directional gain (spreading removed); Tx 1490 m, Rx 6 m; link half-angle %.4f deg',theta_max_deg), ...
    'FontSize',14);
drawnow;
paths = {fullfile(out_dir,[stem '.fig']),fullfile(out_dir,[stem '.png'])};
savefig(fig,paths{1});
exportgraphics(fig,paths{2},'Resolution',180);
close(fig);
end

function D = local_gaussian(theta,f_hz,c_mps,sigma_m)
k=2*pi*f_hz/c_mps;
D=exp(-0.5*(k*sigma_m*sin(theta)).^2);
end

function D = local_piston(theta,f_hz,c_mps,radius_m)
u=(2*pi*f_hz/c_mps)*radius_m*sin(theta);
D=ones(size(u));
nonzero=u~=0;
D(nonzero)=2*besselj(1,u(nonzero))./u(nonzero);
end
