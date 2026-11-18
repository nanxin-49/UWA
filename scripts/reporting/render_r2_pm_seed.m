function receipt=render_r2_pm_seed(root,output_name)
%RENDER_PISTON_PM_R1 Saved 411-point data, truthful per-source normalization.
arguments
    root (1,:) char
    output_name (1,:) char = 'figures'
end
assert(~contains(output_name,{'/','\','..'}));
t0=tic; dest=fullfile(root,output_name); assert(~isfolder(dest)); mkdir(dest);
q=load(fullfile(root,'common.mat')); a=load(fullfile(root,'comparison.mat'));
x=q.x; colors=[0 .447 .698;.1 .1 .1;.835 .369 0]; styles={'--','-','-.'}; names={'PS','BIE','FullFK'};
suffix=''; if ~a.reference_qualified, suffix=' | PROVISIONAL: BIE convergence incomplete'; end
f=figure('Visible','off','Position',[80 80 1150 850],'Color','w'); tl=tiledlayout(f,3,1,'TileSpacing','compact');
ax=nexttile(tl); hold(ax,'on');
for m=1:3, plot(ax,x,abs(a.U(:,m)),styles{m},'Color',colors(m,:),'DisplayName',names{m}); end
ylabel(ax,'Reflected |p| (Pa)'); legend(ax,'Location','eastoutside'); grid(ax,'on');
ax=nexttile(tl); hold(ax,'on');
for m=1:3, plot(ax,x,real(a.U(:,m)),styles{m},'Color',colors(m,:)); end
ylabel(ax,'Re(p) (Pa)'); grid(ax,'on');
ax=nexttile(tl); hold(ax,'on'); plot(ax,q.sx,q.eta,'k-'); yline(ax,0,':');
ylabel(ax,'Surface depth z_s (m)'); xlabel(ax,'Horizontal receiver / surface coordinate x (m)'); grid(ax,'on');
title(tl,[sprintf('R2 | 10 cm strip | 4 kHz | seed%d | Rx z = 3 m',q.cfg.seed) suffix]);
local_save(f,dest,'F01_fields_and_surface');
f=figure('Visible','off','Position',[80 80 1150 1000],'Color','w'); tl=tiledlayout(f,4,1,'TileSpacing','compact');
for panel=1:4
    ax=nexttile(tl); hold(ax,'on');
    for m=[1 3]
        e=a.U(:,m)-a.U(:,2);
        if panel==1, v=abs(e); lab='|p_m - p_{BIE}| (Pa)';
        elseif panel==2, v=abs(e)/a.scale; lab='|residual| / BIE line RMS';
        elseif panel==3, v=abs(a.U(:,m))-abs(a.U(:,2)); lab='|p_m| - |p_{BIE}| (Pa)';
        else, v=angle(a.U(:,m).*conj(a.U(:,2))); v(~a.mask)=NaN; lab='Wrapped phase difference (rad)'; end
        ln=plot(ax,x,v,styles{m},'Color',colors(m,:),'DisplayName',names{m});
        ln.UserData=struct('method',names{m},'panel',panel);
    end
    ylabel(ax,lab); grid(ax,'on'); if panel==1, legend(ax,'Location','eastoutside'); end
    if panel==4, xlabel(ax,'Receiver x (m)'); end
end
title(tl,{['R2 errors relative to BIE' suffix],sprintf('Phase mask: |BIE| >= 0.01 line RMS (%d/411)',nnz(a.mask))});
local_save(f,dest,'F02_pointwise_errors');
f=figure('Visible','off','Position',[80 80 1150 1000],'Color','w'); tl=tiledlayout(f,4,1,'TileSpacing','compact');
for panel=1:4
    ax=nexttile(tl); e=a.U(:,3)-a.U(:,2);
    switch panel
        case 1, v=abs(e); lab='|FullFK - BIE| (Pa)';
        case 2, v=abs(e)/a.scale; lab='|residual| / BIE line RMS';
        case 3, v=abs(a.U(:,3))-abs(a.U(:,2)); lab='|FullFK| - |BIE| (Pa)';
        case 4, v=angle(a.U(:,3).*conj(a.U(:,2))); v(~a.mask)=NaN; lab='Wrapped phase difference (rad)';
    end
    ln=plot(ax,x,v,'-.','Color',colors(3,:),'LineWidth',1.3);
    ln.UserData=struct('method','FullFK','panel',panel); ylabel(ax,lab); grid(ax,'on'); xlim(ax,[min(x) max(x)]);
end
xlabel(ax,'Receiver x (m)'); title(tl,{['FullFK errors: separate vertical scales' suffix],'Same 411 saved receiver points; no smoothing'});
local_save(f,dest,'F04_FullFK_pointwise_errors');
% Reopen and independently reconstruct all F02 plotted data from saved U.
fig=openfig(fullfile(dest,'F02_pointwise_errors.fig'),'new','invisible'); lines=findall(fig,'Type','line'); checked=0;
for j=1:numel(lines)
    info=lines(j).UserData; if ~isstruct(info) || ~isfield(info,'panel'), continue; end
    if strcmp(info.method,'PS'), m=1; else, m=3; end
    switch info.panel
        case 1, expected=abs(a.U(:,m)-a.U(:,2));
        case 2, expected=abs(a.U(:,m)-a.U(:,2))/a.scale;
        case 3, expected=abs(a.U(:,m))-abs(a.U(:,2));
        case 4, expected=angle(a.U(:,m).*conj(a.U(:,2))); expected(~a.mask)=NaN;
    end
    assert(isequaln(lines(j).XData(:),x) && isequaln(lines(j).YData(:),expected)); checked=checked+1;
end
assert(checked==8); fig.Position=[80 80 1150 1000]; drawnow;
local_export(fig,fullfile(dest,'F02_audit_reexport.png')); close(fig);
im1=imread(fullfile(dest,'F02_pointwise_errors.png')); im2=imread(fullfile(dest,'F02_audit_reexport.png'));
assert(isequal(size(im1),size(im2))); pixel_diff=nnz(im1~=im2); assert(pixel_diff==0);
fig=openfig(fullfile(dest,'F04_FullFK_pointwise_errors.fig'),'new','invisible'); lines=findall(fig,'Type','line'); count=0;
for j=1:numel(lines)
    info=lines(j).UserData; if ~isstruct(info) || ~isfield(info,'panel'), continue; end
    switch info.panel
        case 1, expected=abs(a.U(:,3)-a.U(:,2));
        case 2, expected=abs(a.U(:,3)-a.U(:,2))/a.scale;
        case 3, expected=abs(a.U(:,3))-abs(a.U(:,2));
        case 4, expected=angle(a.U(:,3).*conj(a.U(:,2))); expected(~a.mask)=NaN;
    end
    assert(isequaln(lines(j).XData(:),x) && isequaln(lines(j).YData(:),expected)); count=count+1;
end
assert(count==4); fig.Position=[80 80 1150 1000]; drawnow;
local_export(fig,fullfile(dest,'F04_audit_reexport.png')); close(fig);
im1=imread(fullfile(dest,'F04_FullFK_pointwise_errors.png')); im2=imread(fullfile(dest,'F04_audit_reexport.png'));
assert(isequal(size(im1),size(im2))); pixel_diff4=nnz(im1~=im2); assert(pixel_diff4==0);
receipt=struct('status','FIGURES_SAVED_AND_RECONSTRUCTED','failed_checks',0,'checked_error_curves',checked, ...
    'checked_FullFK_detail_curves',count,'reexport_pixel_differences',[pixel_diff pixel_diff4],'elapsed_s',toc(t0));
fid=fopen(fullfile(dest,'receipt.json'),'w'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid); disp(jsonencode(receipt));
end
function local_save(f,dest,name)
drawnow; savefig(f,fullfile(dest,[name '.fig'])); local_export(f,fullfile(dest,[name '.png'])); close(f);
end
function local_export(f,path)
exportgraphics(f,path,'Resolution',160,'Width',1800,'Height',1400,'Units','pixels');
end
