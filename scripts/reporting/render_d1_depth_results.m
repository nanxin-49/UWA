function receipt=render_d1_depth_results(root)
%RENDER_D1_DEPTH_RESULTS Saved-field PNG/FIG, reopened curve reconstruction.
q=load(fullfile(root,'common.mat')); a=load(fullfile(root,'comparison.mat'));
dest=fullfile(root,'figures'); assert(~isfolder(dest)); mkdir(dest);
names={'PS','BIE','FullFK'}; valid=0; styles={'-','--','-.'}; colors=[0 .447 .698;.25 .25 .25;.835 .369 0];
for j=1:4
    if j==1
        Y=abs(a.U); label='Reflected amplitude (peak Pa)';
    elseif j==2
        Y=abs(a.U-a.U(:,2))/a.scale; label='Complex residual / BIE line RMS';
    elseif j==3
        Y=abs(a.U)-abs(a.U(:,2)); label='Amplitude difference (Pa)';
    else
        Y=angle(a.U.*conj(a.U(:,2))); Y(~a.mask,:)=NaN; label='Wrapped phase difference (rad)';
    end
    canvas=[40 40 1000 550];
    f=figure('Visible','off','Color','w','Units','pixels','Position',canvas,'Resize','off','MenuBar','none','ToolBar','none'); c=onCleanup(@()close(f));
    ax=axes('Parent',f,'Units','pixels','Position',[100 75 700 400],'PositionConstraint','innerposition','FontName','Arial','FontSize',11);
    h=plot(ax,q.x,Y,'LineWidth',1.1);
    for method=1:3, h(method).LineStyle=styles{method}; h(method).Color=colors(method,:); h(method).UserData=method; end
    grid(ax,'on'); xlabel(ax,'Horizontal offset x (m)'); ylabel(ax,label);
    yy=Y(isfinite(Y)); assert(~isempty(yy)); lo=min(yy); hi=max(yy); margin=.08*max(hi-lo,1e-12*max(1,max(abs(yy))));
    ax.XLim=[-36 36]; ax.YLim=[lo-margin hi+margin]; ax.XTick=-36:12:36; ax.YTick=linspace(ax.YLim(1),ax.YLim(2),5);
    ax.XAxis.Exponent=0; ax.YAxis.Exponent=0;
    ax.XTickLabel=arrayfun(@(v)sprintf('%g',v),ax.XTick,'UniformOutput',false);
    ax.YTickLabel=arrayfun(@(v)sprintf('%.4g',v),ax.YTick,'UniformOutput',false);
    lg=legend(ax,h,names,'AutoUpdate','off','FontName','Arial','FontSize',11); lg.Units='pixels'; lg.Position=[820 380 150 75];
    title(ax,sprintf('D1a Tx %g m, Rx 6 m, 0.10 m strip, 4 kHz',q.cfg.source.z_source_m)); drawnow;
    base=fullfile(dest,sprintf('F%02d',j)); savefig(f,[base '.fig']); clear c;
    % PNGs are rendered from two independent reads of the editable FIG.
    % This avoids comparing a live auto-layout legend with its serialized form.
    expected_layout=[];
    for pass=1:2
        g=openfig([base '.fig'],'new','invisible'); guard=onCleanup(@()close(g));
        g.Units='pixels'; g.Position=canvas; drawnow; lines=findall(g,'Type','line');
        assert(numel(lines)==3); seen=zeros(1,3);
        for k=1:3
            n=lines(k).UserData; assert(isscalar(n) && ismember(n,1:3)); seen(n)=seen(n)+1;
            assert(isequaln(lines(k).XData(:),q.x) && isequaln(lines(k).YData(:),Y(:,n)),'Saved FIG curve differs from MAT');
            assert(strcmp(lines(k).LineStyle,styles{n}) && isequal(lines(k).Color,colors(n,:)));
            valid=valid+1;
        end
        assert(all(seen==1));
        reopened_axes=findall(g,'Type','axes'); reopened_legend=findall(g,'Type','legend');
        assert(isscalar(reopened_axes) && isscalar(reopened_legend));
        assert(isequal(reopened_axes.XLim,[-36 36]) && isequal(reopened_axes.YLim,[lo-margin hi+margin]));
        assert(strcmp(reopened_axes.YLabel.String,label) && isequal(string(reopened_legend.String(:)),string(names(:))));
        layout=[reopened_axes.Position reopened_legend.Position reopened_legend.FontSize];
        if pass==1, expected_layout=layout; else, assert(isequal(layout,expected_layout)); end
        target=[base '.png']; if pass==2, target=[base '_reopened.png']; end
        exportgraphics(g,target,'Resolution',150,'Width',1000,'Height',550,'Units','pixels');
        clear guard;
    end
    A=imread([base '.png']); B=imread([base '_reopened.png']); assert(isequal(A,B),'Reopened PNG pixel mismatch');
end
receipt=struct('status','PASS','failed_checks',0,'reopened_curves',valid,'phase_valid_points',nnz(a.mask), ...
    'independent_reopens_per_FIG',2,'differing_pixel_channels',0,'export_source','saved FIG in both independent passes');
pm_sensitivity_io('write',fullfile(dest,'figure_audit.json'),receipt);
end
