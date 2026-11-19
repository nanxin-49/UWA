function receipt=render_d1_depth_summary_saved(review,out)
%RENDER_D1_DEPTH_SUMMARY_SAVED Separate PS and FullFK spatial error comparisons.
assert(~isfolder(out)); mkdir(out); t0=tic;
r=jsondecode(fileread(fullfile(review,'receipt.json'))); assert(r.failed_checks==0);
Q=cell(2,1); A=Q;
for j=1:2
    depths=[100 1490]; root=fullfile(review,sprintf('depth_%04dm',depths(j)));
    Q{j}=load(fullfile(root,'common.mat')); A{j}=load(fullfile(root,'comparison.mat'));
end
assert(isequal(Q{1}.x,Q{2}.x)); colors=[0 .447 .698;.835 .369 0]; styles={'-','--'};
methods=[1 3]; names={'PS','FullFK'}; checked=0;
labels={'Absolute complex residual (Pa)','Complex residual / own BIE line RMS','Amplitude difference (Pa)','Wrapped phase difference (rad)'};
positions=[100 600 470 285;790 600 470 285;100 120 470 285;790 120 470 285];
for mm=1:2
    method=methods(mm); canvas=[40 40 1400 1050]; f=figure('Visible','off','Color','w','Units','pixels','Position',canvas,'Resize','off','MenuBar','none','ToolBar','none'); close_live=onCleanup(@()close(f));
    annotation(f,'textbox',[.04 .95 .92 .04],'String',sprintf('D1a %s/BIE | Tx 100 vs 1490 m | Rx 6 m | 0.10 m strip | 4 kHz',names{mm}), ...
        'EdgeColor','none','HorizontalAlignment','center','FontName','Arial','FontSize',15,'Interpreter','none');
    annotation(f,'textbox',[.04 .91 .92 .03],'String',sprintf('Same original PM/C2 realization; 411 direct points; phase valid: %d / %d',nnz(A{1}.mask),nnz(A{2}.mask)), ...
        'EdgeColor','none','HorizontalAlignment','center','FontName','Arial','FontSize',11,'Interpreter','none');
    for panel=1:4
        ax=axes('Parent',f,'Units','normalized','Position',positions(panel,:)./[1400 1050 1400 1050],'PositionConstraint','innerposition','FontName','Arial','FontSize',11); hold(ax,'on'); ys=[];
        for depth=1:2
            y=local_value(A{depth},method,panel); h=plot(ax,Q{depth}.x,y,styles{depth},'Color',colors(depth,:),'LineWidth',1.2);
            h.UserData=struct('depth_index',depth,'panel',panel,'method',method); ys=[ys;y(isfinite(y))]; %#ok<AGROW>
        end
        span=max(max(ys)-min(ys),1e-12); limits=[min(ys)-.06*span max(ys)+.06*span];
        if panel<=2, limits(1)=0; end
        ax.XLim=[-36 36]; ax.XTick=-36:12:36; ax.YLim=limits; ax.YTick=linspace(limits(1),limits(2),5);
        ax.XAxis.Exponent=0; ax.YAxis.Exponent=0;
        ax.XTickLabel=arrayfun(@(x)sprintf('%g',x),ax.XTick,'UniformOutput',false);
        ax.YTickLabel=arrayfun(@(x)sprintf('%.4g',x),ax.YTick,'UniformOutput',false);
        grid(ax,'on'); xlabel(ax,'Horizontal offset x (m)'); ylabel(ax,labels{panel});
        ax.UserData=panel;
        lg=legend(ax,{'Tx 100 m','Tx 1490 m'},'AutoUpdate','off','FontName','Arial','FontSize',10,'Orientation','horizontal');
        lg.Units='normalized'; lg.Position=[positions(panel,1)+50 positions(panel,2)+positions(panel,4)+20 370 45]./[1400 1050 1400 1050];
        lg.UserData=panel;
    end
    drawnow; stem=fullfile(out,sprintf('F%02d_%s_depth_errors',mm+4,names{mm})); savefig(f,[stem '.fig']); clear close_live;
    for pass=1:2
        g=openfig([stem '.fig'],'new','invisible'); close_reopened=onCleanup(@()close(g)); g.Units='pixels'; g.Position=canvas; drawnow;
        axs=findall(g,'Type','axes'); legends=findall(g,'Type','legend');
        for j=1:numel(axs), panel=axs(j).UserData; axs(j).Units='normalized'; axs(j).Position=positions(panel,:)./[1400 1050 1400 1050]; end
        for j=1:numel(legends), panel=legends(j).UserData; legends(j).Units='normalized'; legends(j).Position=[positions(panel,1)+50 positions(panel,2)+positions(panel,4)+20 370 45]./[1400 1050 1400 1050]; end
        drawnow;
        lines=findall(g,'Type','line'); assert(numel(lines)==8); keys=zeros(2,4);
        for j=1:numel(lines)
            u=lines(j).UserData; assert(u.method==method); keys(u.depth_index,u.panel)=keys(u.depth_index,u.panel)+1;
            expected=local_value(A{u.depth_index},method,u.panel);
            assert(isequaln(lines(j).XData(:),Q{u.depth_index}.x) && isequaln(lines(j).YData(:),expected)); checked=checked+1;
        end
        assert(all(keys==1,'all'));
        png=[stem '.png']; if pass==2, png=[stem '_reopened.png']; end
        exportgraphics(g,png,'Resolution',150); clear close_reopened;
    end
    assert(isequal(imread([stem '.png']),imread([stem '_reopened.png'])),'Summary PNG reopen pixel mismatch');
end
receipt=struct('status','PASS','failed_checks',0,'elapsed_s',toc(t0),'checked_curves',checked,'different_pixel_channels',0,'new_BIE_solves',0);
pm_sensitivity_io('write',fullfile(out,'receipt.json'),receipt); fprintf('%s\n',jsonencode(receipt));
end
function y=local_value(a,method,panel)
b=a.U(:,2); model=a.U(:,method);
switch panel
    case 1, y=abs(model-b);
    case 2, y=abs(model-b)/a.scale;
    case 3, y=abs(model)-abs(b);
    case 4, y=angle(model.*conj(b)); y(~a.mask)=NaN;
end
end
