function outputs = generate_bellhop_internal_wall_visuals_vertical(overrides)
%GENERATE_BELLHOP_INTERNAL_WALL_VISUALS_VERTICAL
% Render stored validation-only internal-wall rays and dense coherent fields.
% The stored Bellhop fields are in the mapped chart.  Plots are returned to
% physical coordinates with T^{-1}(r',z')=(200-r',-z'); no propagation is run.
if nargin < 1 || isempty(overrides), overrides = struct(); end
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
cfg = struct('input_mat',fullfile(root,'results','validation', ...
    'bellhop_internal_wall_visuals','dense_wall_fields.mat'), ...
    'output_dir',fullfile(root,'results','visualization', ...
    'bellhop_internal_wall_visuals'),'r0_m',100,'physical_rx_range_m',97, ...
    'z_plot_m',[-32 32],'ray_angles_deg',[-14 -10 -6 -2 0 2 6 10], ...
    'show_figures',false,'tl_floor_db',120,'valid_floor_db',80);
names = fieldnames(overrides);
for ii = 1:numel(names)
    if ~isfield(cfg,names{ii}), error('Unknown override: %s.',names{ii}); end
    cfg.(names{ii}) = overrides.(names{ii});
end
if ~isfile(cfg.input_mat), error('Missing dense field MAT: %s',cfg.input_mat); end
if ~isfolder(cfg.output_dir), mkdir(cfg.output_dir); end
S = load(cfg.input_mat,'rt','rs','ri','z','rr');
cases = {S.rt,'tilted','r=100+0.05z'; S.rs,'sinusoidal','r=100-2sin(0.04z)'};
ray_tables = cell(2,1); rays = cell(2,1);
for cc = 1:2
    rays{cc} = local_extract_rays(cases{cc,1},cfg);
    ray_tables{cc} = local_ray_table(rays{cc},cases{cc,2});
end
ray_table = vertcat(ray_tables{:});
writetable(ray_table,fullfile(cfg.output_dir,'ray_trajectory_audit.csv'));
outputs = struct('ray_png',fullfile(cfg.output_dir,'internal_wall_ray_trajectories.png'), ...
    'tl_png',fullfile(cfg.output_dir,'internal_wall_tl_2d.png'), ...
    'incident_tl_png',fullfile(cfg.output_dir,'internal_wall_incident_tl_2d.png'), ...
    'ray_csv',fullfile(cfg.output_dir,'ray_trajectory_audit.csv'), ...
    'manifest',fullfile(cfg.output_dir,'visualization_manifest.mat'));
local_plot_rays(cases,rays,cfg,outputs.ray_png);
local_plot_tl(cases,rays,cfg,outputs.tl_png);
local_plot_incident_tl(S.ri,cases,rays,cfg,outputs.incident_tl_png);
manifest = struct('cfg',cfg,'cases',{cases},'outputs',outputs, ...
    'description','Bellhop coherent 2-D TL fields mapped from the positive-range internal-wall chart'); %#ok<NASGU>
save(outputs.manifest,'manifest','-v7.3');
    local_write_report(root,cfg,ray_table,cases);
fprintf('Saved internal-wall trajectory and 2-D TL figures under %s.\n',cfg.output_dir);
end

function rays = local_extract_rays(result,cfg)
d = result.diagnostics;
residual_name = 'wall_residual';
if ~ismember(residual_name,d.Properties.VariableNames), residual_name = 'line_residual'; end
angles = cfg.ray_angles_deg(:);
rays = repmat(struct('alpha_deg',NaN,'hit',[],'incident',[],'reflected',[], ...
    'rotated',[],'incident_xy',[],'native_xy',[],'mapped_xy',[],'residual_m',NaN),numel(angles),1);
for ii = 1:numel(angles)
    [~,ix] = min(abs(d.alpha_deg-angles(ii)));
    hit = [d.hit_r(ix); d.hit_z(ix)];
    inc = [d.inc_ur(ix); d.inc_uz(ix)];
    ref = [d.ref_ur(ix); d.ref_uz(ix)];
    rot = [d.rot_ur(ix); d.rot_uz(ix)];
    if abs(ref(1)) < 1e-10 || abs(rot(1)) < 1e-10
        error('Selected ray has near-zero range direction at %.6g deg.',angles(ii));
    end
    % Incident segment in the mapped chart.
    % Native backward branch evaluated up to physical r=97 m.
    s_native = (cfg.physical_rx_range_m-hit(1))/ref(1);
    native_end = hit+s_native*ref;
    % Properly rotated positive-range branch, then inverse mapped to physics.
    mapped_hit = [2*cfg.r0_m-hit(1); -hit(2)];
    s_mapped = (2*cfg.r0_m-cfg.physical_rx_range_m-mapped_hit(1))/rot(1);
    mapped = mapped_hit + (0:1).*s_mapped.*rot;
    mapped_xy = [2*cfg.r0_m-mapped(1,:); -mapped(2,:)];
    % Densify only for display; all endpoints remain the stored geometry.
    tt = linspace(0,1,80);
    native_xy = hit + (native_end-hit).*tt;
    mapped_hit_phys = mapped_xy(:,1); mapped_end_phys = mapped_xy(:,end);
    mapped_xy = mapped_hit_phys + (mapped_end_phys-mapped_hit_phys).*tt;
    rays(ii).alpha_deg = d.alpha_deg(ix); rays(ii).hit = hit;
    rays(ii).incident = inc; rays(ii).reflected = ref; rays(ii).rotated = rot;
    rays(ii).incident_xy = [hit(1).*tt; hit(2).*tt];
    rays(ii).native_xy = native_xy; rays(ii).mapped_xy = mapped_xy;
    rays(ii).residual_m = d.(residual_name)(ix);
end
end

function t = local_ray_table(rays,kind)
n = numel(rays); t = table(strings(n,1),zeros(n,1),zeros(n,1),zeros(n,1), ...
    zeros(n,1),zeros(n,1),zeros(n,1),zeros(n,1), ...
    'VariableNames',{'wall','alpha_deg','hit_r_m','hit_z_m','intersection_residual_m', ...
    'inc_ur','inc_uz','ref_ur'});
for ii = 1:n
    t.wall(ii)=string(kind); t.alpha_deg(ii)=rays(ii).alpha_deg;
    t.hit_r_m(ii)=rays(ii).hit(1); t.hit_z_m(ii)=rays(ii).hit(2);
    t.intersection_residual_m(ii)=rays(ii).residual_m;
    t.inc_ur(ii)=rays(ii).incident(1); t.inc_uz(ii)=rays(ii).incident(2);
    t.ref_ur(ii)=rays(ii).reflected(1);
end
end

function local_plot_rays(cases,rays,cfg,file)
visible='off'; if cfg.show_figures, visible='on'; end
fig=figure('Visible',visible,'Color','w','Position',[40 40 1500 980]);
cleanup=onCleanup(@()close(fig)); tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
colors=lines(numel(rays{1}));
for cc=1:2
    wall_kind=cases{cc,2}; ax=nexttile(tl); hold(ax,'on');
    z=linspace(cfg.z_plot_m(1),cfg.z_plot_m(2),401); rw=local_wall(wall_kind,z);
    plot(ax,rw,z,'k-','LineWidth',2.0,'DisplayName','internal wall');
    for ii=1:numel(rays{cc})
        q=rays{cc}(ii);
        plot(ax,q.incident_xy(1,:),q.incident_xy(2,:),'-','Color',[0.05 0.35 0.80], ...
            'LineWidth',1.35,'HandleVisibility','off');
        plot(ax,q.native_xy(1,:),q.native_xy(2,:),'--','Color',[0.55 0.55 0.55], ...
            'LineWidth',2.2,'HandleVisibility','off');
        plot(ax,q.mapped_xy(1,:),q.mapped_xy(2,:),'-','Color',colors(ii,:), ...
            'LineWidth',0.9,'HandleVisibility','off');
        mid=round(size(q.incident_xy,2)/2);
        local_arrow(ax,q.incident_xy(:,mid),q.incident_xy(:,mid+1),[0.05 0.35 0.80]);
        local_arrow(ax,q.mapped_xy(:,mid),q.mapped_xy(:,mid+1),colors(ii,:));
    end
    plot(ax,0,0,'kd','MarkerFaceColor',[1 .75 .1],'MarkerSize',9,'DisplayName','Tx');
    plot(ax,cfg.physical_rx_range_m,0,'kp','MarkerFaceColor',[.2 .85 .3], ...
        'MarkerSize',12,'DisplayName','physical Rx (97 m)');
    set(ax,'YDir','reverse'); grid(ax,'on'); box(ax,'on'); axis(ax,'equal');
    xlim(ax,[0 103]); ylim(ax,cfg.z_plot_m); xlabel(ax,'physical range r (m)'); ylabel(ax,'physical depth z (m)');
    title(ax,sprintf('%s wall: physical (r,z)',wall_kind));
    if cc==1, legend(ax,'Location','southwest'); end
    axz=nexttile(tl); hold(axz,'on');
    plot(axz,z,rw,'k-','LineWidth',2.0);
    for ii=1:numel(rays{cc})
        q=rays{cc}(ii); plot(axz,q.incident_xy(2,:),q.incident_xy(1,:),'-','Color',[0.05 0.35 0.80],'LineWidth',1.0);
        plot(axz,q.native_xy(2,:),q.native_xy(1,:),'--','Color',[0.55 0.55 0.55],'LineWidth',1.7);
        plot(axz,q.mapped_xy(2,:),q.mapped_xy(1,:),'-','Color',colors(ii,:),'LineWidth',0.8);
    end
    set(axz,'YDir','normal'); grid(axz,'on'); box(axz,'on');
    xlim(axz,cfg.z_plot_m); ylim(axz,[96 104]); xlabel(axz,'z (m) — display-axis reorder'); ylabel(axz,'range r (m)');
    title(axz,'wall neighborhood (axes reordered for display)');
end
title(tl,{'Bellhop internal-wall ray trajectories','Blue: incident | gray dashed: native backward branch | colored: inverse-mapped proper-rotation branch'},'FontWeight','bold');
exportgraphics(fig,file,'Resolution',180); clear cleanup
end

function local_plot_tl(cases,rays,cfg,file)
visible='off'; if cfg.show_figures, visible='on'; end
fig=figure('Visible',visible,'Color','w','Position',[40 40 1500 720]); cleanup=onCleanup(@()close(fig));
tl=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
for cc=1:2
    ax=nexttile(tl); data=cases{cc,1}.data; rp=double(data.receiver_range_m(:)); zp=-double(data.receiver_depth_m(:));
    [rp,ix]=sort(2*cfg.r0_m-rp); [zp,iz]=sort(zp);
    p=double(data.pressure(iz,ix)); [tl_db,valid]=local_tl_map(p,cfg.valid_floor_db); %#ok<ASGLU>
    [R,Z]=meshgrid(rp,zp);
    surf(ax,R,Z,tl_db,'EdgeColor','none'); view(ax,2); shading(ax,'interp'); hold(ax,'on');
    zdraw=local_overlay_height(tl_db);
    ax.Color=[0.88 0.88 0.88]; colormap(ax,turbo); clim(ax,[40 105]); colorbar(ax);
    z=linspace(cfg.z_plot_m(1),cfg.z_plot_m(2),401); rw=local_wall(cases{cc,2},z);
    plot3(ax,rw,z,zdraw*ones(size(z)),'k-','LineWidth',2.0);
    r_end=min(rp);
    for ii=1:numel(rays{cc})
        q=rays{cc}(ii); xy=local_reflected_path_to_range(q,r_end);
        plot3(ax,xy(1,:),xy(2,:),zdraw*ones(1,size(xy,2)),'-','Color',[1 1 1], ...
            'LineWidth',0.85,'HandleVisibility','off');
    end
    plot3(ax,cfg.physical_rx_range_m,0,zdraw,'wp','MarkerFaceColor',[.1 .8 .2],'MarkerSize',10);
    xlabel(ax,'physical range r (m)'); ylabel(ax,'physical depth z (m)'); title(ax,sprintf('%s wall: coherent reflected-only TL=-20log10|p|',cases{cc,2}));
    xlim(ax,[min(rp) max(rp)]); ylim(ax,cfg.z_plot_m); set(ax,'YDir','reverse'); grid(ax,'on'); box(ax,'on');
end
title(tl,{'Bellhop-style 2-D coherent reflected-only TL','Smooth TL background; invalid/no-beam samples are masked; rays are inverse-mapped to physical (r,z)'},'FontWeight','bold');
exportgraphics(fig,file,'Resolution',180); clear cleanup
end

function local_plot_incident_tl(incident,cases,rays,cfg,file)
visible='off'; if cfg.show_figures, visible='on'; end
fig=figure('Visible',visible,'Color','w','Position',[40 40 1500 720]); cleanup=onCleanup(@()close(fig));
tl=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
data=incident.data; rp=double(data.receiver_range_m(:)); zp=double(data.receiver_depth_m(:));
[rp,ix]=sort(rp); [zp,iz]=sort(zp); p=double(data.pressure(iz,ix)); [tl_db,valid]=local_tl_map(p,cfg.valid_floor_db); %#ok<ASGLU>
[R,Z]=meshgrid(rp,zp);
for cc=1:2
    ax=nexttile(tl); surf(ax,R,Z,tl_db,'EdgeColor','none'); view(ax,2); shading(ax,'interp'); hold(ax,'on');
    zdraw=local_overlay_height(tl_db);
    ax.Color=[0.88 0.88 0.88]; colormap(ax,turbo); clim(ax,[-35 65]); colorbar(ax);
    z=linspace(cfg.z_plot_m(1),cfg.z_plot_m(2),401);
    plot3(ax,local_wall(cases{cc,2},z),z,zdraw*ones(size(z)),'k-','LineWidth',2.0);
    for ii=1:numel(rays{cc})
        q=rays{cc}(ii); plot3(ax,q.incident_xy(1,:),q.incident_xy(2,:), ...
            zdraw*ones(1,size(q.incident_xy,2)),'-','Color',[1 1 1], ...
            'LineWidth',0.85,'HandleVisibility','off');
    end
    plot3(ax,0,0,zdraw,'wd','MarkerFaceColor',[1 .75 .1],'MarkerSize',9);
    plot3(ax,cfg.physical_rx_range_m,0,zdraw,'wp','MarkerFaceColor',[.2 .85 .3],'MarkerSize',11);
    xlabel(ax,'physical range r (m)'); ylabel(ax,'physical depth z (m)');
    title(ax,sprintf('%s wall: coherent incident-only TL=-20log10|p|',cases{cc,2}));
    xlim(ax,[0 max(rp)]); ylim(ax,cfg.z_plot_m); set(ax,'YDir','reverse'); grid(ax,'on'); box(ax,'on');
end
title(tl,{'Bellhop-style 2-D coherent incident-only TL','Smooth TL background; no reflected branch; white curves are complete incident rays'},'FontWeight','bold');
exportgraphics(fig,file,'Resolution',180); clear cleanup
end

function [tl_db,valid]=local_tl_map(p,valid_floor_db)
amp=abs(p); scale=max(amp(:));
if isempty(scale) || ~isfinite(scale) || scale<=0, scale=1; end
valid=amp>=scale*10^(-valid_floor_db/20) & amp>0 & isfinite(amp);
tl_db=-20*log10(max(amp,realmin('double')));
tl_db(~valid)=NaN;
end

function zdraw=local_overlay_height(tl_db)
v=tl_db(isfinite(tl_db));
if isempty(v), zdraw=1; else, zdraw=max(v)+1; end
end

function xy=local_reflected_path_to_range(q,r_end)
if abs(q.reflected(1))<eps, xy=q.native_xy; return; end
s=(r_end-q.hit(1))/q.reflected(1); endpoint=q.hit+s*q.reflected;
tt=linspace(0,1,120); xy=q.hit+(endpoint-q.hit).*tt;
end

function r=local_wall(kind,z)
if strcmp(kind,'tilted'), r=100+0.05*z; else, r=100-2*sin(0.04*z); end
end

function local_arrow(ax,a,b,color)
v=b-a; if norm(v)==0, return; end
quiver(ax,a(1),a(2),0.18*v(1),0.18*v(2),0,'Color',color,'LineWidth',1.0,'MaxHeadSize',2,'HandleVisibility','off');
end

function local_write_report(root,cfg,ray_table,cases)
file=fullfile(root,'reports','bellhop_internal_wall_visualization_report.md'); fid=fopen(file,'w','n','UTF-8');
if fid<0, error('Cannot write %s.',file); end; cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'# Bellhop internal-wall visualization\n\n');
fprintf(fid,'本报告由已保存的 Bellhop validation-only SHD/diagnostic 数据生成；入射场 SHD 来自官方 matched-halfspace `C*X` run，绘图脚本本身不重新运行 PE 或 Bellhop。\n\n');
fprintf(fid,'- walls: `%s`, `%s`; z support: [%.0f, %.0f] m; physical Rx: `(%.0f,0)` m.\n',cases{1,3},cases{2,3},cfg.z_plot_m(1),cfg.z_plot_m(2),cfg.physical_rx_range_m);
fprintf(fid,'- inverse map: `r=200-r''`, `z=-z''`; dense SHD contains only the post-wall reflected branch and is plotted as `TL=-20 log10(abs(p))`.\n');
fprintf(fid,'- both TL figures use smooth interpolated Bellhop-style backgrounds; samples below %.0f dB relative to the field maximum (or with no finite beam support) are masked as NaN and shown with the axes background.\n',cfg.valid_floor_db);
fprintf(fid,'- incident-only field uses the same dense receiver depth grid and extends past the wall support for context; the two wall curves are overlays, not reflecting boundaries in that run.\n');
fprintf(fid,'- selected representative rays: %d per wall; gray dashed paths are native backward branches and colored paths are inverse-mapped proper-rotation branches.\n\n',numel(cfg.ray_angles_deg));
fprintf(fid,'## Ray diagnostics\n\n');
fprintf(fid,'| wall | alpha (deg) | hit r (m) | hit z (m) | residual (m) | inc u_r | inc u_z | ref u_r |\n|---|---:|---:|---:|---:|---:|---:|---:|\n');
for ii=1:height(ray_table)
 fprintf(fid,'| %s | %.6g | %.9g | %.9g | %.3e | %.6g | %.6g | %.6g |\n',ray_table.wall(ii),ray_table.alpha_deg(ii),ray_table.hit_r_m(ii),ray_table.hit_z_m(ii),ray_table.intersection_residual_m(ii),ray_table.inc_ur(ii),ray_table.inc_uz(ii),ray_table.ref_ur(ii));
end
fprintf(fid,'\n![ray trajectories](../results/visualization/bellhop_internal_wall_visuals/internal_wall_ray_trajectories.png)\n\n![reflected-only 2-D TL](../results/visualization/bellhop_internal_wall_visuals/internal_wall_tl_2d.png)\n\n![incident-only 2-D TL](../results/visualization/bellhop_internal_wall_visuals/internal_wall_incident_tl_2d.png)\n');
end
