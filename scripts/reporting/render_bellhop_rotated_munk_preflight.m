function output = render_bellhop_rotated_munk_preflight(input_root,output_dir)
%RENDER_BELLHOP_ROTATED_MUNK_PREFLIGHT Saved-evidence figures, no acoustic solve.
arguments
    input_root (1,:) char
    output_dir (1,:) char
end
assert(~isfolder(output_dir)); mkdir(output_dir); render_start=tic;
[~,mem]=memory; assert(mem.PhysicalMemory.Available>=512*2^20,'RAM preflight gate.');
ray_data=load(fullfile(input_root,'20261114_s1_rays_01','result.mat'),'runs','receipt');
axis_data=load(fullfile(input_root,'20261114_s1_direct_02','result.mat'),'runs','configs','receipt');
assert(strcmp(ray_data.receipt.status,'COMPLETE') && strcmp(axis_data.receipt.status,'COMPLETE'));
colors=[0 .35 .65;.75 .25 0;.2 .55 .3]; styles={'-','--',':'};
fig=figure('Visible','off','Color','w','Position',[80 80 1500 680]);
tiles=tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
a1=nexttile(tiles); hold(a1,'on');
up=ray_data.runs{1}; down=ray_data.runs{2};
for j=1:3
    rq=up.rays(j).rq;
    keep=rq(:,1)<=4950;
    plot(a1,rq(keep,2)-1000,4950-rq(keep,1),'LineStyle',styles{j}, ...
        'Color',colors(j,:),'DisplayName',sprintf('up %.2g deg',up.rays(j).angle_deg));
end
plot(a1,0,4950,'ko','MarkerFaceColor','k','DisplayName','source');
yline(a1,0,'k--','surface (not applied)'); yline(a1,5000,'k:','bottom (not applied)');
set(a1,'YDir','reverse'); xlim(a1,[-5 5]); ylim(a1,[0 5000]); grid(a1,'on');
xlabel(a1,'physical x (m)'); ylabel(a1,'physical depth z (m)');
title(a1,{'Physical coordinates','direct rays only; unequal axis scales'});
a2=nexttile(tiles); hold(a2,'on');
for j=1:3
    rq=up.rays(j).rq; keep=rq(:,1)<=4950;
    plot(a2,rq(keep,1),rq(keep,2)-1000,'LineStyle',styles{j},'Color',colors(j,:), ...
        'DisplayName',sprintf('%.2g deg',up.rays(j).angle_deg));
end
xline(a2,4950,'k--','future surface wall','HandleVisibility','off'); xlim(a2,[0 5000]); ylim(a2,[-5 5]); grid(a2,'on');
xlabel(a2,'upward chart r = 4950 - z (m)'); ylabel(a2,'q - q_0 = x (m)');
title(a2,{'Upward proper rotation','Q SSP: c(r) = c_{Munk}(4950-r)'});
legend(a2,'Location','southwest');
a3=nexttile(tiles); hold(a3,'on');
for j=1:3
    rq=down.rays(j).rq; keep=rq(:,1)<=50;
    plot(a3,rq(keep,1),rq(keep,2)-1000,'LineStyle',styles{j},'Color',colors(j,:), ...
        'DisplayName',sprintf('%.2g deg',down.rays(j).angle_deg));
end
xline(a3,50,'k:','future fluid-bottom wall','HandleVisibility','off'); xlim(a3,[0 50]); ylim(a3,[-5 5]); grid(a3,'on');
xlabel(a3,'downward chart r = z - 4950 (m)'); ylabel(a3,'q - q_0 = -x (m)');
title(a3,{'Downward proper rotation','Q SSP: c(r) = c_{Munk}(4950+r)'});
legend(a3,'Location','southwest');
title(tiles,'S1 coordinate/path validation; surface/bottom reflection and total field NOT computed');
save_pair(fig,'F01_rotation_direct_rays');
fig2=figure('Visible','off','Color','w','Position',[80 80 1250 600]);
tiles=tiledlayout(fig2,1,2,'TileSpacing','compact','Padding','compact');
for id=1:2
    ax=nexttile(tiles); hold(ax,'on');
    for j=1:3
        row=axis_data.runs{(id-1)*3+j};
        cfg=axis_data.configs{(id-1)*3+j};
        s=abs(row.axis_delay_s-row.reference_delay_s)*1e6;
        plot(ax,cfg.receiver_range_m,s,'LineStyle',styles{j},'Marker','o','Color',colors(j,:), ...
            'DisplayName',sprintf('Q %.1f m / step %.1f m',cfg.ssp_spacing_m,cfg.step_m));
    end
    yline(ax,.05/(2*pi*10000)*1e6,'k--','0.05 rad at 10 kHz','HandleVisibility','off');
    xlabel(ax,'rotated receiver range (m)'); ylabel(ax,'absolute time error (microseconds)');
    grid(ax,'on'); ylim(ax,[0 .85]); legend(ax,'Location','northwest');
    if id==1, title(ax,'Upward direct axis'); else, title(ax,'Downward direct axis'); end
end
title(tiles,{'Direct-axis arrival time vs continuous Munk integral','Joint refinement pilot; separate-factor results and independent replay are saved'});
save_pair(fig2,'F02_axis_travel_time');
save(fullfile(output_dir,'source_data.mat'),'ray_data','axis_data');
html=['<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>旋转 Munk S1 先导</title>' ...
    '<style>body{font:17px/1.6 sans-serif;max-width:1500px;margin:auto;padding:24px}img{width:100%}strong{color:#9c2200}</style>' ...
    '<h1>Bellhop 旋转 Munk：S1 直接传播先导</h1><p><strong>尚无反射场、完整复声场或 Fig.12(b) 复现结论。</strong>' ...
    '图中线条来自已保存 Bellhop 射线，不是插值补场。两种旋转均保持正定向；人工墙尚未应用。坐标轴比例不相等。</p>' ...
    '<h2>坐标与路径</h2><img src="F01_rotation_direct_rays.png" alt="物理坐标与向上、向下旋转坐标的三条直接传播射线">' ...
    '<p><a href="F01_rotation_direct_rays.fig">可编辑 FIG</a></p><h2>轴向时延</h2>' ...
    '<img src="F02_axis_travel_time.png" alt="Q采样和步长加密的上下行轴向时延误差，含0.05rad等效门">' ...
    '<p><a href="F02_axis_travel_time.fig">可编辑 FIG</a> · <a href="source_data.mat">原始绘图源</a></p></html>'];
fid=fopen(fullfile(output_dir,'index.html'),'w','n','UTF-8'); assert(fid>=0); fprintf(fid,'%s',html); fclose(fid);
output=struct('status','FIGURES_CREATED_NOT_YET_INDEPENDENTLY_AUDITED','output_dir',output_dir,'elapsed_s',toc(render_start), ...
    'scope','Direct geometry/time only; no surface, fluid-bottom, total field, or FFP field comparison');
fid=fopen(fullfile(output_dir,'checkpoint.json'),'w'); fprintf(fid,'%s',jsonencode(output,PrettyPrint=true)); fclose(fid);
disp(jsonencode(output));
    function save_pair(handle,stem)
        drawnow; savefig(handle,fullfile(output_dir,[stem '.fig']));
        exportgraphics(handle,fullfile(output_dir,[stem '.png']),'Resolution',120);
        close(handle);
    end
end