function receipt=render_k1_pm_summary_fields(records,out)
%RENDER_K1_PM_SUMMARY_FIELDS Saved field plots; per-surface BIE RMS, no solves.
assert(numel(records)==7 && ~isfolder(out)); mkdir(out);
assert(all([records.test_only]) || ~any([records.test_only])); total=0;
for method=[1 3]
    fig=figure('Visible','off','Units','pixels','Position',[40 40 1800 1400], ...
        'Color','w','Resize','off','MenuBar','none','ToolBar','none'); guard=onCleanup(@()close(fig));
    name='PS'; if method==3, name='FullFK'; end
    tag='K1'; if records(1).test_only, tag='TEST_ONLY'; end
    annotation(fig,'textbox',[.07 .958 .88 .025],'String', ...
        sprintf('%s | %s vs independent BIE | seed260001 | 0.10 m strip, 4 kHz',tag,name), ...
        'EdgeColor','none','HorizontalAlignment','center','FontSize',15,'Interpreter','none');
    annotation(fig,'textbox',[.07 .934 .88 .022],'String', ...
        '411 actual receivers; each surface uses its own BIE RMS and fixed 0.01 RMS phase mask', ...
        'EdgeColor','none','HorizontalAlignment','center','FontSize',11,'Interpreter','none');
    bottoms=[1030 710 390 70]; kinds={'absolute','normalized','amplitude','phase'};
    expected=struct('x',{},'y',{},'case_index',{},'panel',{}); count=0;
    for panel=1:4
        ax=axes(fig,'Units','pixels','Position',[150 bottoms(panel) 1330 255],'PositionConstraint','innerposition', ...
            'FontSize',11,'Box','on','XGrid','on','YGrid','on','LooseInset',[0 0 0 0]); hold(ax,'on'); all_y=[];
        for j=1:7
            r=records(j); b=r.U(:,2); u=r.U(:,method);
            assert(numel(r.x)==411 && abs(r.scale-norm(b)/sqrt(411))<1e-12 && isequal(r.mask,abs(b)>=.01*r.scale));
            switch kinds{panel}
                case 'absolute', y=abs(u-b);
                case 'normalized', y=abs(u-b)/r.scale;
                case 'amplitude', y=abs(u)-abs(b);
                case 'phase', y=angle(u.*conj(b)); y(~r.mask)=NaN;
            end
            color=[0 .447 .698]; if r.U_PM_mps==5, color=[.835 .369 0]; elseif r.U_PM_mps==8, color=[.009 .55 .35]; end
            style='-'; if r.N==7, style='--'; elseif r.N==40, style=':'; end
            line=plot(ax,r.x,y,'Color',color,'LineStyle',style,'LineWidth',1.1,'DisplayName',r.case_id);
            line.UserData=struct('case_index',j,'panel',panel);
            count=count+1; expected(count)=struct('x',r.x(:),'y',y(:),'case_index',j,'panel',panel);
            all_y=[all_y;y(isfinite(y))]; %#ok<AGROW>
        end
        top=max(abs(all_y))*1.10; if isempty(top) || top==0, top=1; end
        limits=[-top top]; if panel<=2, limits=[0 top]; end
        ax.XLim=[-40.0390625 40.0390625]; ax.YLim=limits; ax.XTick=-40:10:40; ax.YTick=linspace(limits(1),limits(2),5);
        ax.XTickLabel=arrayfun(@(v)sprintf('%g',v),ax.XTick,'UniformOutput',false);
        ax.YTickLabel=arrayfun(@(v)sprintf('%.4g',v),ax.YTick,'UniformOutput',false);
        ax.XAxis.Exponent=0; ax.YAxis.Exponent=0;
        labels={'|p_m-p_{BIE}| (Pa)','|residual| / own BIE line RMS','|p_m|-|p_{BIE}| (Pa)','Wrapped phase difference (rad)'};
        ylabel(ax,labels{panel}); xlabel(ax,'Receiver x (m)');
        if panel==1
            lg=legend(ax,'Location','northeast','Interpreter','none','FontSize',10,'AutoUpdate','off');
            lg.Units='pixels'; lg.Position=[1500 1050 280 220];
        end
    end
    drawnow; stem=['K1_' name '_seven_conditions']; path=fullfile(out,[stem '.fig']); savefig(fig,path);
    exportgraphics(fig,fullfile(out,[stem '.png']),'Resolution',160,'Width',1800,'Height',1400,'Units','pixels');
    clear guard
    f=openfig(path,'new','invisible'); fg=onCleanup(@()close(f)); f.Units='pixels'; f.Position=[40 40 1800 1400]; f.Resize='off'; drawnow;
    lines=findall(f,'Type','line'); assert(numel(lines)==28); seen=false(7,4);
    for j=1:numel(lines)
        key=lines(j).UserData; ix=find([expected.case_index]==key.case_index & [expected.panel]==key.panel);
        assert(isscalar(ix) && ~seen(key.case_index,key.panel)); seen(key.case_index,key.panel)=true;
        assert(isequaln(lines(j).XData(:),expected(ix).x) && isequaln(lines(j).YData(:),expected(ix).y));
    end
    assert(all(seen,'all')); pixels=zeros(1,2); original=imread(fullfile(out,[stem '.png']));
    for pass=1:2
        target=fullfile(out,sprintf('%s_reopen_%d.png',stem,pass));
        exportgraphics(f,target,'Resolution',160,'Width',1800,'Height',1400,'Units','pixels');
        reconstructed=imread(target); assert(isequal(size(original),size(reconstructed)));
        pixels(pass)=nnz(original~=reconstructed); assert(pixels(pass)==0,'k1:SummaryPixels','Summary PNG pixel mismatch');
    end
    clear fg
    total=total+count;
end
receipt=struct('status','K1_SUMMARY_FIGURES_PASS','failed_checks',0,'curves_checked',total, ...
    'zero_pixel_reexports',4,'test_only',records(1).test_only,'acoustic_calls',0);
k1_pm_io('write',fullfile(out,'receipt.json'),receipt);
end
