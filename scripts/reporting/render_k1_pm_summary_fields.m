function receipt=render_k1_pm_summary_fields(records,out)
%RENDER_K1_PM_SUMMARY_FIELDS Seven saved conditions, fixed full-paper panels.
assert(numel(records)==7 && ~isfolder(out)); mkdir(out);
assert(all([records.test_only]) || ~any([records.test_only])); total=0; differences=struct;
kinds={'absolute','normalized','amplitude','phase'};
labels={'Absolute complex residual (Pa)','Complex residual / own BIE line RMS','Amplitude difference (Pa)','Wrapped phase difference (rad)'};
save(fullfile(out,'saved_field_records.mat'),'records','-v7');
for method=[1 3]
    name='PS'; if method==3, name='FullFK'; end
    for panel=1:4
        stem=['K1_' name '_' kinds{panel}];
        f=figure('Visible','off','WindowStyle','normal','Units','pixels','Position',[10 10 1000 550], ...
            'Color','w','Resize','off','MenuBar','none','ToolBar','none'); guard=onCleanup(@()close(f));
        ax=axes('Parent',f,'Units','pixels','Position',[105 75 680 400],'PositionConstraint','innerposition', ...
            'FontName','Arial','FontSize',11,'Box','on','XGrid','on','YGrid','on'); hold(ax,'on'); all_y=[];
        for j=1:7
            r=records(j); assert(numel(r.x)==411 && abs(r.scale-norm(r.U(:,2))/sqrt(411))<1e-12 && isequal(r.mask,abs(r.U(:,2))>=.01*r.scale));
            y=local_value(r,method,panel);
            color=[0 .447 .698]; if r.U_PM_mps==5, color=[.835 .369 0]; elseif r.U_PM_mps==8, color=[.009 .55 .35]; end
            style='-'; if r.N==7, style='--'; elseif r.N==40, style=':'; end
            h=plot(ax,r.x,y,'Color',color,'LineStyle',style,'LineWidth',1.1,'DisplayName',sprintf('U=%g, N=%d',r.U_PM_mps,r.N));
            if r.U_PM_mps==5, h.Marker='o'; elseif r.U_PM_mps==8, h.Marker='s'; end
            h.MarkerIndices=1:24:411; h.MarkerSize=3;
            h.UserData=struct('case_index',j,'panel',panel,'method',method); all_y=[all_y;y(isfinite(y))]; %#ok<AGROW>
        end
        bound=max(abs(all_y))*1.1; if bound==0, bound=1; end
        limits=[-bound bound]; if panel<=2, limits=[0 bound]; end
        ax.XLim=[-40.0390625 40.0390625]; ax.YLim=limits; ax.XTick=-40:10:40; ax.YTick=linspace(limits(1),limits(2),5);
        ax.XTickLabel=arrayfun(@(v)sprintf('%g',v),ax.XTick,'UniformOutput',false);
        ax.YTickLabel=arrayfun(@(v)sprintf('%.4g',v),ax.YTick,'UniformOutput',false);
        ax.XAxis.Exponent=0; ax.YAxis.Exponent=0; xlabel(ax,'Receiver x (m)'); ylabel(ax,labels{panel});
        tag='K1'; if records(1).test_only, tag='TEST ONLY'; end
        title(ax,sprintf('%s %s / own independent BIE | seed 260001 | 0.10 m strip | 4 kHz',tag,name),'FontSize',10,'Interpreter','none');
        lg=legend(ax,'Location','northeast','FontName','Arial','FontSize',10,'AutoUpdate','off','Interpreter','none');
        lg.Units='pixels'; lg.Position=[780 245 175 220]; drawnow;
        path=fullfile(out,[stem '.fig']); savefig(f,path); clear guard;
        layouts=cell(1,3); images=cell(1,3);
        for pass=1:3
            f=openfig(path,'new','invisible'); guard=onCleanup(@()close(f)); f.Units='pixels'; f.Position=[10 10 1000 550]; drawnow;
            lines=findall(f,'Type','line'); assert(numel(lines)==7); seen=false(1,7);
            for k=1:7
                u=lines(k).UserData; assert(u.method==method && u.panel==panel && ~seen(u.case_index)); seen(u.case_index)=true;
                r=records(u.case_index); expected=local_value(r,method,panel);
                assert(isequaln(lines(k).XData(:),r.x(:)) && isequaln(lines(k).YData(:),expected(:)));
            end
            assert(all(seen)); axs=findall(f,'Type','axes'); legends=findall(f,'Type','legend'); assert(isscalar(axs) && isscalar(legends));
            assert(isequal(axs.Position,[105 75 680 400]) && isequal(legends.Position,[780 245 175 220]));
            assert(isequal(axs.XLim,[-40.0390625 40.0390625]) && isequal(axs.YLim,limits));
            layout=[f.Position axs.Position legends.Position axs.XLim axs.YLim axs.FontSize legends.FontSize]; layouts{pass}=layout;
            f.PaperUnits='inches'; f.PaperPosition=[0 0 10 5.5]; f.PaperSize=[10 5.5]; f.PaperPositionMode='manual';
            target=fullfile(out,[stem '.png']); if pass>1, target=fullfile(out,sprintf('%s_reopen%d.png',stem,pass-1)); end
            print(f,target,'-dpng','-r100');
            assert(isequal(layout,[f.Position axs.Position legends.Position axs.XLim axs.YLim axs.FontSize legends.FontSize]));
            images{pass}=imread(target); assert(isequal(size(images{pass}),[550 1000 3])); clear guard;
        end
        assert(isequal(layouts{1},layouts{2}) && isequal(layouts{1},layouts{3}));
        pixels=[nnz(images{1}~=images{2}) nnz(images{1}~=images{3})]; differences.(stem)=pixels;
        assert(all(pixels==0),'k1:SummaryPixels','Summary PNG pixel mismatch'); total=total+7;
    end
end
receipt=struct('status','K1_SUMMARY_FIGURES_PASS','failed_checks',0,'curves_checked',total,'zero_pixel_reexports',16, ...
    'pixel_differences',differences,'pixel_tolerance',0,'layout_version','K1-small-panels-v3-fixed-paper','test_only',records(1).test_only,'acoustic_calls',0);
k1_pm_io('write',fullfile(out,'receipt.json'),receipt);
end
function y=local_value(r,method,panel)
u=r.U(:,method); b=r.U(:,2);
switch panel
    case 1, y=abs(u-b);
    case 2, y=abs(u-b)/r.scale;
    case 3, y=abs(u)-abs(b);
    case 4, y=angle(u.*conj(b)); y(~r.mask)=NaN;
end
end
