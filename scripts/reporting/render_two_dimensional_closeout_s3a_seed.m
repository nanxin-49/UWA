function files = render_two_dimensional_closeout_s3a_seed(root,run_dir)
%RENDER_TWO_DIMENSIONAL_CLOSEOUT_S3A_SEED Saved-data PM figures; no solve.
arguments
    root (1,:) char
    run_dir (1,:) char
end
assert(isfolder(root) && startsWith(run_dir,root));
dest=fullfile(run_dir,'figures'); assert(isfolder(dest));
q=load(fullfile(run_dir,'analysis','analysis_data.mat'),'x','z','ix9','flat','rough','bres');
p=readtable(fullfile(run_dir,'input','seed_260001_coefficients.csv'));
xx=linspace(-60,60,2401)';
[eta,slope,~,raw]=closeout_s3a_pm_surface(p,xx);
names={'F01_PM_geometry_v2','F02_nine_point_errors_v2','F03_method_line_v2'};
for j=1:numel(names)
    assert(~isfile(fullfile(dest,[names{j} '.fig'])) && ...
        ~isfile(fullfile(dest,[names{j} '.png'])));
end
files=cell(3,2);
fig=figure('Visible','off','Color','w','Position',[100 100 1160 720]);
t=tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
ax=nexttile(t); plot(ax,xx,raw,'Color',[.65 .65 .65], ...
    'DisplayName','原始周期PM'); hold(ax,'on');
plot(ax,xx,eta,'b','LineWidth',1.4,'DisplayName','C2嵌入PM');
xline(ax,-50,':k','HandleVisibility','off');
xline(ax,50,':k','HandleVisibility','off');
xlim(ax,[-60 60]); ylabel(ax,'海面深度 η (m)'); grid(ax,'on'); legend(ax,'Location','best');
title(ax,'seed260001：同一保存Fourier实现；|x|≥50 m严格平面');
ax=nexttile(t); plot(ax,xx,slope,'Color',[.2 .45 .2]);
xlim(ax,[-60 60]); xlabel(ax,'水平位置 x (m)'); ylabel(ax,'坡度 dη/dx'); grid(ax,'on');
local_save(fig,dest,names{1}); files(1,:)={fullfile(dest,[names{1} '.fig']),fullfile(dest,[names{1} '.png'])};
fig=figure('Visible','off','Color','w','Position',[100 100 1190 760]);
t=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
ii=q.ix9(:); x=q.x(ii); z=q.z(ii);
labels=compose('P%d%d',ceil((1:9)'/3),mod((1:9)'-1,3)+1);
% Use explicit tile positions to keep shared color meaning fixed.
for k=1:4
    ax=nexttile(t,k);
    m=1; if mod(k,2)==0, m=3; end
    e=q.rough(ii,m)-q.rough(ii,2);
    if k<=2
        val=abs(e); maxv=max(abs(q.rough(ii,1)-q.rough(ii,2)));
        scatter(ax,x,z,175,val,'filled','MarkerEdgeColor','k');
        clim(ax,[0 maxv]); colormap(ax,parula); cb=colorbar(ax); cb.Label.String='|u_m-u_{BIE}|';
        if m==1, title(ax,'PE-mode0：复场残差'); else, title(ax,'完整FullFK：复场残差（同色标）'); end
    else
        val=angle(q.rough(ii,m).*conj(q.rough(ii,2)));
        scatter(ax,x,z,175,val,'filled','MarkerEdgeColor','k');
        clim(ax,[-.14 .14]); colormap(ax,turbo); cb=colorbar(ax); cb.Label.String='相位差 (rad)';
        if m==1, title(ax,'PE-mode0：相位差'); else, title(ax,'完整FullFK：相位差（同色标）'); end
    end
    set(ax,'YDir','reverse'); xlim(ax,[-.28 .28]); ylim(ax,[2.82 3.18]);
    xlabel(ax,'x (m)'); ylabel(ax,'z (m，向下)'); grid(ax,'on');
    for h=1:9, text(ax,x(h)+.008,z(h),labels(h),'FontSize',8); end
end
local_save(fig,dest,names{2}); files(2,:)={fullfile(dest,[names{2} '.fig']),fullfile(dest,[names{2} '.png'])};
fig=figure('Visible','off','Color','w','Position',[100 100 1180 740]);
t=tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
line=1:411; x=q.x(line); b=q.rough(line,2);
ax=nexttile(t); plot(ax,x,abs(b),'k','LineWidth',1.5,'DisplayName','BIE PM'); hold(ax,'on');
plot(ax,x,abs(q.rough(line,1)),'--','Color',[.8 .25 .2],'DisplayName','PE-mode0 PM');
plot(ax,x,abs(q.rough(line,3)),':','Color',[.1 .45 .8],'LineWidth',1.5,'DisplayName','完整FullFK PM');
plot(ax,x,abs(q.flat(line,2)),'-.','Color',[.4 .4 .4],'DisplayName','BIE 平面');
xlabel(ax,'x (m)'); ylabel(ax,'|u|（同源原生复场）'); title(ax,'z=3 m 原411点线：PM与同源平面');
legend(ax,'Location','best'); grid(ax,'on');
ax=nexttile(t); plot(ax,x,abs(q.rough(line,1)-b),'Color',[.8 .25 .2], ...
    'DisplayName','|PE−BIE|'); hold(ax,'on');
plot(ax,x,abs(q.rough(line,3)-b),'Color',[.1 .45 .8], ...
    'DisplayName','|FullFK−BIE|');
plot(ax,x,q.bres(line),'k:','DisplayName','BIE五层最大逐点变化');
xlabel(ax,'x (m)'); ylabel(ax,'复场绝对差'); legend(ax,'Location','best'); grid(ax,'on');
local_save(fig,dest,names{3}); files(3,:)={fullfile(dest,[names{3} '.fig']),fullfile(dest,[names{3} '.png'])};
end

function local_save(fig,dest,name)
drawnow;
savefig(fig,fullfile(dest,[name '.fig']));
exportgraphics(fig,fullfile(dest,[name '.png']),'Resolution',180);
close(fig);
end
