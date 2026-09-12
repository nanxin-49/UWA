function outputs = generate_bellhop_internal_wall_visuals_vertical(overrides)
%GENERATE_BELLHOP_INTERNAL_WALL_VISUALS_VERTICAL
% Render stored validation-only internal-wall rays and dense coherent fields.
% The stored Bellhop fields are in the mapped chart.  Plots are returned to
% physical coordinates with T^{-1}(r',z')=(200-r',-z'); no propagation is run.
if nargin < 1 || isempty(overrides), overrides = struct(); end
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
cfg = struct('input_mat',fullfile(root,'results','validation', ...
    'bellhop_internal_wall_two_period_full_fan','large_sinusoidal_wall_dense_fields.mat'), ...
    'output_dir',fullfile(root,'results','visualization', ...
    'bellhop_internal_wall_visuals'),'r0_m',100,'physical_rx_range_m',97, ...
    'z_plot_m',[-65 65],'wall_zoom_z_m',[-30 30],'wall_support_z_m',[-75 75], ...
    'ray_angles_deg',[-28 -20 -12 -4 0 4 12 20 28], ...
    'target_eigenray_angle_deg',1.23855637169814, ...
    'show_figures',false,'tl_floor_db',120,'valid_floor_db',80, ...
    'reference_range_m',1,'common_clim_db',[-20 100], ...
    'tilted_slope_m_per_m',0.05,'sinusoidal_amplitude_m',2.4, ...
    'sinusoidal_wavenumber_per_m',2*pi/30, ...
    'report_file',fullfile(root,'reports','bellhop_internal_wall_visualization_report.md'));
names = fieldnames(overrides);
for ii = 1:numel(names)
    if ~isfield(cfg,names{ii}), error('Unknown override: %s.',names{ii}); end
    cfg.(names{ii}) = overrides.(names{ii});
end
if isempty(cfg.wall_zoom_z_m), cfg.wall_zoom_z_m=cfg.z_plot_m; end
if isempty(cfg.wall_support_z_m), cfg.wall_support_z_m=cfg.z_plot_m; end
if ~isfile(cfg.input_mat), error('Missing dense field MAT: %s',cfg.input_mat); end
if ~isfolder(cfg.output_dir), mkdir(cfg.output_dir); end
S = load(cfg.input_mat,'rt','rs','ri','z','rr');
cfg.reference_amp = local_reference_amplitude(S.ri,cfg.reference_range_m);
cases = {S.rt,'tilted',sprintf('r=%.6g%+.6gz',cfg.r0_m,cfg.tilted_slope_m_per_m); ...
    S.rs,'sinusoidal',sprintf('r=%.6g-%.6gsin(%.6gz)',cfg.r0_m, ...
    cfg.sinusoidal_amplitude_m,cfg.sinusoidal_wavenumber_per_m)};
ray_tables = cell(2,1); rays = cell(2,1);
for cc = 1:2
    cfg_case=cfg;
    if cc==2 && isfinite(cfg.target_eigenray_angle_deg)
        cfg_case.ray_angles_deg=sort(unique([cfg.ray_angles_deg(:);cfg.target_eigenray_angle_deg]));
    end
    rays{cc} = local_extract_rays(cases{cc,1},cfg_case);
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
    local_write_report(root,cfg,ray_table,cases,S.ri);
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
    for cc=1:2
        colors=lines(numel(rays{cc}));
        target_idx=[];
        if cc==2 && isfinite(cfg.target_eigenray_angle_deg)
            [~,target_idx]=min(abs([rays{cc}.alpha_deg]-cfg.target_eigenray_angle_deg));
        end
    wall_kind=cases{cc,2}; ax=nexttile(tl); hold(ax,'on');
    z=linspace(cfg.z_plot_m(1),cfg.z_plot_m(2),401); rw=local_wall(wall_kind,z,cfg);
    plot(ax,rw,z,'k-','LineWidth',2.0,'DisplayName','internal wall');
    for ii=1:numel(rays{cc})
        q=rays{cc}(ii);
        is_target=isequal(ii,target_idx);
        inc_color=[0.05 0.35 0.80]; mapped_color=colors(ii,:);
        inc_width=1.35; mapped_width=0.9; inc_visibility='off'; mapped_visibility='off';
        inc_name=''; mapped_name='';
        if is_target
            inc_color=[0.90 0.15 0.05]; mapped_color=[0.90 0.15 0.05];
            inc_width=2.6; mapped_width=2.6; inc_visibility='on'; mapped_visibility='on';
            inc_name=sprintf('Rx eigenray incident (%.3f deg)',q.alpha_deg);
            mapped_name='Rx eigenray reflected';
        end
        plot(ax,q.incident_xy(1,:),q.incident_xy(2,:),'-','Color',inc_color, ...
            'LineWidth',inc_width,'HandleVisibility',inc_visibility,'DisplayName',inc_name);
        plot(ax,q.native_xy(1,:),q.native_xy(2,:),'--','Color',[0.55 0.55 0.55], ...
            'LineWidth',2.2,'HandleVisibility','off');
        plot(ax,q.mapped_xy(1,:),q.mapped_xy(2,:),'-','Color',mapped_color, ...
            'LineWidth',mapped_width,'HandleVisibility',mapped_visibility,'DisplayName',mapped_name);
        mid=round(size(q.incident_xy,2)/2);
        local_arrow(ax,q.incident_xy(:,mid),q.incident_xy(:,mid+1),inc_color);
        local_arrow(ax,q.mapped_xy(:,mid),q.mapped_xy(:,mid+1),mapped_color);
    end
    plot(ax,0,0,'kd','MarkerFaceColor',[1 .75 .1],'MarkerSize',9,'DisplayName','Tx');
    plot(ax,cfg.physical_rx_range_m,0,'kp','MarkerFaceColor',[.2 .85 .3], ...
        'MarkerSize',12,'DisplayName','physical Rx (97 m)');
    set(ax,'YDir','reverse'); grid(ax,'on'); box(ax,'on'); axis(ax,'equal');
    xlim(ax,[0 103]); ylim(ax,cfg.z_plot_m); xlabel(ax,'physical range r (m)'); ylabel(ax,'physical depth z (m)');
    title(ax,sprintf('%s wall: physical (r,z)',wall_kind));
    legend(ax,'Location','southwest');
    axz=nexttile(tl); hold(axz,'on');
    z_zoom=linspace(cfg.wall_zoom_z_m(1),cfg.wall_zoom_z_m(2),401);
    rw_zoom=local_wall(wall_kind,z_zoom,cfg);
    plot(axz,z_zoom,rw_zoom,'k-','LineWidth',2.0);
    for ii=1:numel(rays{cc})
        q=rays{cc}(ii); is_target=isequal(ii,target_idx);
        inc_color=[0.05 0.35 0.80]; mapped_color=colors(ii,:); inc_width=1.0; mapped_width=0.8;
        if is_target, inc_color=[0.90 0.15 0.05]; mapped_color=inc_color; inc_width=2.4; mapped_width=2.4; end
        plot(axz,q.incident_xy(2,:),q.incident_xy(1,:),'-','Color',inc_color,'LineWidth',inc_width);
        plot(axz,q.native_xy(2,:),q.native_xy(1,:),'--','Color',[0.55 0.55 0.55],'LineWidth',1.7);
        plot(axz,q.mapped_xy(2,:),q.mapped_xy(1,:),'-','Color',mapped_color,'LineWidth',mapped_width);
        % In reordered display coordinates (x=z,y=r), draw the local wall
        % normal explicitly.  The equal-aspect setting below is essential:
        % otherwise the 64 m x 8 m limits visibly distort the physical angle.
        [~,nw]=local_wall_frame(wall_kind,q.hit(2),cfg); nd=[nw(2);nw(1)];
        hh=[q.hit(2);q.hit(1)]; ln=1.1;
        plot(axz,hh(1)+[-ln ln]*nd(1),hh(2)+[-ln ln]*nd(2),'m--','LineWidth',0.8,'HandleVisibility','off');
    end
    set(axz,'YDir','normal'); grid(axz,'on'); box(axz,'on'); axis(axz,'equal');
    % Use a wider local range window so equal physical units remain readable
    % while the wall and all reflected branches stay in view.
    wall_pad=max(1,0.35*cfg.sinusoidal_amplitude_m);
    xlim(axz,cfg.wall_zoom_z_m);
    ylim(axz,[cfg.r0_m-cfg.sinusoidal_amplitude_m-wall_pad, ...
        cfg.r0_m+cfg.sinusoidal_amplitude_m+wall_pad]);
    xlabel(axz,'z (m) — display-axis reorder'); ylabel(axz,'range r (m)');
    title(axz,'wall neighborhood (equal aspect; magenta = wall normal)');
end
title(tl,{'Bellhop internal-wall ray trajectories','Blue: incident | gray dashed: native backward branch | colored: inverse-mapped branch | red: Rx eigenray'},'FontWeight','bold');
exportgraphics(fig,file,'Resolution',180); clear cleanup
end

function local_plot_tl(cases,rays,cfg,file)
visible='off'; if cfg.show_figures, visible='on'; end
fig=figure('Visible',visible,'Color','w','Position',[40 40 1500 720]); cleanup=onCleanup(@()close(fig));
tl=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
for cc=1:2
    ax=nexttile(tl); data=cases{cc,1}.data; rp=double(data.receiver_range_m(:)); zp=-double(data.receiver_depth_m(:));
    [rp,ix]=sort(2*cfg.r0_m-rp); [zp,iz]=sort(zp);
    p=double(data.pressure(iz,ix)); [tl_db,valid]=local_tl_map(p,cfg.reference_amp,cfg.valid_floor_db); %#ok<ASGLU>
    [R,Z]=meshgrid(rp,zp);
    surf(ax,R,Z,tl_db,'EdgeColor','none'); view(ax,2); shading(ax,'interp'); hold(ax,'on');
    zdraw=local_overlay_height(tl_db);
    ax.Color=[0.88 0.88 0.88]; colormap(ax,turbo); clim(ax,cfg.common_clim_db); colorbar(ax);
    z=linspace(cfg.wall_support_z_m(1),cfg.wall_support_z_m(2),601); rw=local_wall(cases{cc,2},z,cfg);
    plot3(ax,rw,z,zdraw*ones(size(z)),'k-','LineWidth',2.0);
    r_end=min(rp);
    for ii=1:numel(rays{cc})
        q=rays{cc}(ii); xy=local_reflected_path_to_range(q,r_end);
        plot3(ax,xy(1,:),xy(2,:),zdraw*ones(1,size(xy,2)),'-','Color',[1 1 1], ...
            'LineWidth',0.85,'HandleVisibility','off');
    end
    plot3(ax,cfg.physical_rx_range_m,0,zdraw,'wp','MarkerFaceColor',[.1 .8 .2],'MarkerSize',10);
    xlabel(ax,'physical range r (m)'); ylabel(ax,'physical depth z (m)'); title(ax,sprintf('%s wall: reflected-only TL=-20log10(|p|/p_{1m})',cases{cc,2}));
    wall_pad=max(0.5,0.02*(max(rp)-min(rp)));
    xlim(ax,[min([rp;rw(:)])-wall_pad,max([rp;rw(:)])+wall_pad]);
    ylim(ax,cfg.wall_support_z_m); set(ax,'YDir','reverse'); grid(ax,'on'); box(ax,'on');
end
title(tl,{'Bellhop-style 2-D coherent reflected-only TL','Smooth TL background; invalid/no-beam samples are masked; rays are inverse-mapped to physical (r,z)'},'FontWeight','bold');
exportgraphics(fig,file,'Resolution',180); clear cleanup
end

function local_plot_incident_tl(incident,cases,rays,cfg,file)
visible='off'; if cfg.show_figures, visible='on'; end
fig=figure('Visible',visible,'Color','w','Position',[40 40 1500 720]); cleanup=onCleanup(@()close(fig));
tl=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
data=incident.data; rp=double(data.receiver_range_m(:)); zp=double(data.receiver_depth_m(:));
[rp,ix]=sort(rp); [zp,iz]=sort(zp); p=double(data.pressure(iz,ix)); [tl_db,valid]=local_tl_map(p,cfg.reference_amp,cfg.valid_floor_db); %#ok<ASGLU>
[R,Z]=meshgrid(rp,zp);
for cc=1:2
    ax=nexttile(tl); surf(ax,R,Z,tl_db,'EdgeColor','none'); view(ax,2); shading(ax,'interp'); hold(ax,'on');
    zdraw=local_overlay_height(tl_db);
    ax.Color=[0.88 0.88 0.88]; colormap(ax,turbo); clim(ax,cfg.common_clim_db); colorbar(ax);
    z=linspace(cfg.wall_support_z_m(1),cfg.wall_support_z_m(2),601);
    rw=local_wall(cases{cc,2},z,cfg);
    plot3(ax,rw,z,zdraw*ones(size(z)),'k-','LineWidth',2.0);
    for ii=1:numel(rays{cc})
        q=rays{cc}(ii); plot3(ax,q.incident_xy(1,:),q.incident_xy(2,:), ...
            zdraw*ones(1,size(q.incident_xy,2)),'-','Color',[1 1 1], ...
            'LineWidth',0.85,'HandleVisibility','off');
    end
    plot3(ax,0,0,zdraw,'wd','MarkerFaceColor',[1 .75 .1],'MarkerSize',9);
    plot3(ax,cfg.physical_rx_range_m,0,zdraw,'wp','MarkerFaceColor',[.2 .85 .3],'MarkerSize',11);
    xlabel(ax,'physical range r (m)'); ylabel(ax,'physical depth z (m)');
    title(ax,sprintf('%s wall: incident-only TL=-20log10(|p|/p_{1m})',cases{cc,2}));
    wall_pad=max(0.5,0.01*max(rp));
    xlim(ax,[0,max([rp;rw(:)])+wall_pad]); ylim(ax,cfg.wall_support_z_m);
    set(ax,'YDir','reverse'); grid(ax,'on'); box(ax,'on');
end
title(tl,{'Bellhop-style 2-D coherent incident-only TL','Smooth TL background; no reflected branch; white curves are complete incident rays'},'FontWeight','bold');
exportgraphics(fig,file,'Resolution',180); clear cleanup
end

function [tl_db,valid]=local_tl_map(p,p_ref,valid_floor_db)
amp=abs(p); scale=max(amp(:));
if isempty(scale) || ~isfinite(scale) || scale<=0, scale=1; end
valid=amp>=scale*10^(-valid_floor_db/20) & amp>0 & isfinite(amp);
tl_db=-20*log10(max(amp,realmin('double'))/p_ref);
tl_db(~valid)=NaN;
end

function p_ref=local_reference_amplitude(incident,reference_range_m)
% Shared axial source reference for both incident and reflected maps.
% The value is read/interpolated from the stored incident-only Bellhop run;
% no propagation or source re-fitting is performed by this renderer.
data=incident.data;
rp=double(data.receiver_range_m(:)); zp=double(data.receiver_depth_m(:));
p=double(data.pressure);
[rp,ix]=sort(rp); [~,iz]=min(abs(zp)); axis_p=p(iz,ix);
p_ref=abs(interp1(rp,axis_p,reference_range_m,'linear','extrap'));
if ~isfinite(p_ref) || p_ref<=0
    error('Invalid common source reference at %.6g m.',reference_range_m);
end
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

function r=local_wall(kind,z,cfg)
if strcmp(kind,'tilted')
    r=cfg.r0_m+cfg.tilted_slope_m_per_m*z;
else
    r=cfg.r0_m-cfg.sinusoidal_amplitude_m*sin(cfg.sinusoidal_wavenumber_per_m*z);
end
end

function [tw,nw]=local_wall_frame(kind,z,cfg)
if strcmp(kind,'tilted')
    drdz=cfg.tilted_slope_m_per_m;
else
    drdz=-cfg.sinusoidal_amplitude_m*cfg.sinusoidal_wavenumber_per_m* ...
        cos(cfg.sinusoidal_wavenumber_per_m*z);
end
tw=[drdz;1]; tw=tw/norm(tw);
nw=[tw(2);-tw(1)];
end

function local_arrow(ax,a,b,color)
v=b-a; if norm(v)==0, return; end
quiver(ax,a(1),a(2),0.18*v(1),0.18*v(2),0,'Color',color,'LineWidth',1.0,'MaxHeadSize',2,'HandleVisibility','off');
end

function local_write_report(root,cfg,ray_table,cases,incident)
file=cfg.report_file; fid=fopen(file,'w','n','UTF-8');
if fid<0, error('Cannot write %s.',file); end; cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'# Bellhop internal-wall visualization\n\n');
fprintf(fid,'本报告由已保存的 Bellhop validation-only SHD/diagnostic 数据生成；入射场 SHD 来自官方 matched-halfspace coherent run，绘图脚本本身不重新运行 PE 或 Bellhop。\n\n');
fprintf(fid,'- walls: `%s`, `%s`; sampled wall support: [%.0f, %.0f] m; plotted z window: [%.0f, %.0f] m; physical Rx: `(%.0f,0)` m.\n', ...
    cases{1,3},cases{2,3},cfg.wall_support_z_m(1),cfg.wall_support_z_m(2), ...
    cfg.z_plot_m(1),cfg.z_plot_m(2),cfg.physical_rx_range_m);
fprintf(fid,'- inverse map: `r=200-r''`, `z=-z''`; dense SHD contains only the post-wall reflected branch. Both TL maps use the same axial source reference `p_ref=|p_inc(r=%.6g m,z=0)|=%.9g` and are plotted as `TL=-20 log10(abs(p)/p_ref)`.\n',cfg.reference_range_m,cfg.reference_amp);
fprintf(fid,'- This common reference only harmonizes the visualization; it does not renormalize the reflected branch or alter Bellhop physics. The reflected-only field remains a separate post-wall branch and must not be interpreted as a direct-plus-reflected total field.\n');
fprintf(fid,'- both TL figures use smooth interpolated Bellhop-style backgrounds; samples below %.0f dB relative to the field maximum (or with no finite beam support) are masked as NaN and shown with the axes background.\n',cfg.valid_floor_db);
fprintf(fid,'- incident-only field uses the same dense receiver depth grid and extends past the wall support for context; the two wall curves are overlays, not reflecting boundaries in that run.\n');
fprintf(fid,'- TL axes include the complete sampled wall support. Regions outside the computed SHD receiver grid remain background/NaN; no field values are extrapolated to fill them.\n');
n_tilted=nnz(ray_table.wall=="tilted"); n_sinusoidal=nnz(ray_table.wall=="sinusoidal");
fprintf(fid,'- selected representative rays: tilted %d, sinusoidal %d; the latter includes the Rx eigenray near %.6g deg. Gray dashed paths are native backward branches and colored paths are inverse-mapped proper-rotation branches.\n\n', ...
    n_tilted,n_sinusoidal,cfg.target_eigenray_angle_deg);
fprintf(fid,'## 当前图件介绍与分析\n\n');
fprintf(fid,'两张二维 TL 图的背景都是 Bellhop 密集接收网格上的相干复声压；TL 仅表示压力幅度，不是单条声线能量或总场功率。入射图只包含 Tx 到墙前的场，反射图只包含一次 wall reflection 后、经 proper rotation 并逆映射回物理坐标的 reflected-only 分支。\n\n');
fprintf(fid,'### 统一参考与绝对幅度\n\n');
fprintf(fid,'两张图现在共同使用 `p_ref=|p_inc(r=%.6g m,z=0)|=%.9g`，色标统一为 `[%.6g, %.6g] dB`。这一统一只改变显示参考，不对任何分支做幅度拟合、重标定或物理修正。\n\n',cfg.reference_range_m,cfg.reference_amp,cfg.common_clim_db(1),cfg.common_clim_db(2));
st_inc=local_report_field_stats(incident,cfg);
st_t=local_report_field_stats(cases{1,1},cfg);
st_s=local_report_field_stats(cases{2,1},cfg);
fprintf(fid,'| field | max `|p|` (Bellhop units) | common-reference TL range (dB) | median TL (dB) | valid fraction |\n|---|---:|---:|---:|---:|\n');
fprintf(fid,'| incident-only | %.6g | [%.3f, %.3f] | %.3f | %.3f |\n',st_inc.max_amp,st_inc.min_tl,st_inc.max_tl,st_inc.median_tl,st_inc.valid_fraction);
fprintf(fid,'| tilted reflected-only | %.6g | [%.3f, %.3f] | %.3f | %.3f |\n',st_t.max_amp,st_t.min_tl,st_t.max_tl,st_t.median_tl,st_t.valid_fraction);
fprintf(fid,'| sinusoidal reflected-only | %.6g | [%.3f, %.3f] | %.3f | %.3f |\n\n',st_s.max_amp,st_s.min_tl,st_s.max_tl,st_s.median_tl,st_s.valid_fraction);
ratio_t_db=20*log10(st_t.max_amp/st_inc.max_amp);
ratio_s_db=20*log10(st_s.max_amp/st_inc.max_amp);
fprintf(fid,['统一参考后，入射、tilted reflected-only 和 sinusoidal reflected-only 的最大原始幅度分别为 ' ...
    '`%.6g`、`%.6g` 和 `%.6g`；两个反射图最大值相对入射图最大值分别为 `%.3f dB` 和 `%.3f dB`。' ...
    '这些是不同空间区域内的场最大值，只用于说明共同色标下的显示动态范围，不能直接解释为墙面的能量反射率。\n\n'], ...
    st_inc.max_amp,st_t.max_amp,st_s.max_amp,ratio_t_db,ratio_s_db);
fprintf(fid,'### 图形形状、灰色区域与平滑性\n\n');
sin_formula=cases{2,3};
fprintf(fid,'入射图使用有限 Gaussian 发射扇区，灰色三角区域对应没有有效 beam support 的位置（当前入射网格有效比例约为 `%.2f%%`）。反射图的物理 range 来自 `r=200-r''`，因此正 range mapped branch 映射回物理坐标后，反射场自然显示在墙的左侧。tilted case 的有效比例约为 `%.2f%%`；sinusoidal case（`%s`）约为 `%.2f%%`。\n\n',100*st_inc.valid_fraction,100*st_t.valid_fraction,sin_formula,100*st_s.valid_fraction);
fprintf(fid,['背景的平滑性来自相干 Bellhop beam influence 在密集二维接收网格上的累积以及绘图时的面内插值；当前没有海底、多次反射或直达/反射总场干涉，因此不应期待标准多途 Bellhop 图中常见的强烈干涉条纹。' ...
    '图上的代表性声线（tilted %d 条、sinusoidal %d 条）只是几何叠加，不是 TL 背景中的全部 beam 或能量脊线。掩膜边界的少量台阶来自接收网格和 NaN support 边界，而不是传播算子台阶。\n\n'],n_tilted,n_sinusoidal);
fprintf(fid,'### 声线几何解释\n\n');
g_t=local_report_geometry_stats(cases{1,1});
g_s=local_report_geometry_stats(cases{2,1});
fprintf(fid,['声线不需要垂直撞击墙面；正确条件是入射、反射方向关于局部墙面法线镜像对称。保存的全射线几何审计给出 tilted/sinusoidal 最大镜面方向误差分别为 ' ...
    '`%.3e` 和 `%.3e`，对应 `max|t dot n|` 分别为 `%.3e` 和 `%.3e`。因此轨迹图中看似斜入射的线是正确的 specular reflection，而不是法线或比例错误。\n\n'], ...
    g_t.max_specular_error,g_s.max_specular_error,g_t.max_t_dot_n,g_s.max_t_dot_n);
fprintf(fid,'## Ray diagnostics\n\n');
fprintf(fid,'| wall | alpha (deg) | hit r (m) | hit z (m) | residual (m) | inc u_r | inc u_z | ref u_r |\n|---|---:|---:|---:|---:|---:|---:|---:|\n');
for ii=1:height(ray_table)
 fprintf(fid,'| %s | %.6g | %.9g | %.9g | %.3e | %.6g | %.6g | %.6g |\n',ray_table.wall(ii),ray_table.alpha_deg(ii),ray_table.hit_r_m(ii),ray_table.hit_z_m(ii),ray_table.intersection_residual_m(ii),ray_table.inc_ur(ii),ray_table.inc_uz(ii),ray_table.ref_ur(ii));
end
viz_rel=strrep(erase(cfg.output_dir,[root filesep]),filesep,'/');
fprintf(fid,'\n![ray trajectories](../%s/internal_wall_ray_trajectories.png)\n\n![reflected-only 2-D TL](../%s/internal_wall_tl_2d.png)\n\n![incident-only 2-D TL](../%s/internal_wall_incident_tl_2d.png)\n', ...
    viz_rel,viz_rel,viz_rel);
end

function st=local_report_field_stats(result,cfg)
p=double(result.data.pressure); amp=abs(p);
finite_amp=amp(isfinite(amp));
if isempty(finite_amp), scale=1; else, scale=max(finite_amp); end
valid=amp>=scale*10^(-cfg.valid_floor_db/20) & amp>0 & isfinite(amp);
tl=-20*log10(max(amp,realmin('double'))/cfg.reference_amp); v=tl(valid);
if isempty(v), v=NaN; end
st=struct('max_amp',scale,'min_tl',min(v),'max_tl',max(v), ...
    'median_tl',median(v),'valid_fraction',nnz(valid)/numel(valid));
end

function st=local_report_geometry_stats(result)
d=result.diagnostics;
st=struct('max_specular_error',max(abs(d.specular_error)), ...
    'max_t_dot_n',max(abs(d.wall_t_r.*d.wall_n_r+d.wall_t_z.*d.wall_n_z)));
end
