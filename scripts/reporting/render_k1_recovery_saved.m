function receipt=render_k1_recovery_saved(work_dirs,out)
%RENDER_K1_RECOVERY_SAVED Saved K1 arrays, fixed small single-panel layout.
arguments
    work_dirs (1,1) cell
    out (1,:) char
end
assert(~isfolder(out)); mkdir(out);
q=load(fullfile(work_dirs{1},'common.mat'),'cfg','x'); a=load(fullfile(work_dirs{1},'comparison.mat'));
assert(a.reference_qualified && numel(q.x)==411);
kinds={'absolute','normalized','amplitude','phase'};
labels={'|p_m - p_{BIE}| (Pa)','|residual| / BIE line RMS','|p_m| - |p_{BIE}| (Pa)','Wrapped phase difference (rad)'};
for family=1:2
    methods=[1 3]; prefix='F02'; if family==2, methods=3; prefix='F04'; end
    for panel=1:4
        stem=[prefix '_' kinds{panel}];
        f=figure('Visible','off','WindowStyle','normal','Units','pixels','Position',[10 10 1000 550], ...
            'Color','w','Resize','off','MenuBar','none','ToolBar','none'); guard=onCleanup(@()local_close(f));
        ax=axes('Parent',f,'Units','pixels','Position',[100 75 700 400],'PositionConstraint','innerposition', ...
            'FontName','Arial','FontSize',11,'Box','on','XGrid','on','YGrid','on','LooseInset',[0 0 0 0]); hold(ax,'on');
        values=[];
        for m=methods
            name='PS'; color=[0 .447 .698]; style='--';
            if m==3, name='FullFK'; color=[.835 .369 0]; style='-.'; end
            switch kinds{panel}
                case 'absolute', y=abs(a.U(:,m)-a.U(:,2));
                case 'normalized', y=abs(a.U(:,m)-a.U(:,2))/a.scale;
                case 'amplitude', y=abs(a.U(:,m))-abs(a.U(:,2));
                case 'phase', y=angle(a.U(:,m).*conj(a.U(:,2))); y(~a.mask)=NaN;
            end
            ln=plot(ax,q.x,y,style,'Color',color,'LineWidth',1.15,'DisplayName',name);
            ln.UserData=struct('seed',q.cfg.seed,'method',name,'kind',kinds{panel},'panel',panel);
            values=[values;y(isfinite(y))]; %#ok<AGROW>
        end
        bound=max(abs(values))*1.1; if bound==0, bound=1; end
        limits=[-bound bound]; if panel<=2, limits=[0 bound]; end
        ax.XLim=[min(q.x) max(q.x)]; ax.YLim=limits; ax.XTick=-40:10:40; ax.YTick=linspace(limits(1),limits(2),5);
        ax.XTickLabel=arrayfun(@(v)sprintf('%g',v),ax.XTick,'UniformOutput',false);
        ax.YTickLabel=arrayfun(@(v)sprintf('%.4g',v),ax.YTick,'UniformOutput',false); ax.XAxis.Exponent=0; ax.YAxis.Exponent=0;
        xlabel(ax,'Receiver x (m)'); ylabel(ax,labels{panel});
        title(ax,sprintf('U_PM=%g m/s | N=%d | seed %d | 0.10 m strip | 4 kHz',q.cfg.U_PM_mps,q.cfg.modes,q.cfg.seed),'Interpreter','none','FontSize',11);
        lg=legend(ax,'Location','northeast','FontName','Arial','FontSize',10,'AutoUpdate','off'); lg.Units='pixels'; lg.Position=[805 365 140 75];
        drawnow; savefig(f,fullfile(out,[stem '.fig'])); clear guard;
        f=openfig(fullfile(out,[stem '.fig']),'new','invisible'); guard=onCleanup(@()local_close(f));
        f.Units='pixels'; f.Position=[10 10 1000 550]; drawnow;
        f.PaperUnits='inches'; f.PaperPosition=[0 0 10 5.5]; f.PaperSize=[10 5.5]; f.PaperPositionMode='manual';
        print(f,fullfile(out,[stem '.png']),'-dpng','-r100'); clear guard;
    end
end
receipt=struct('status','SAVED_SINGLE_PANEL_GENERATION_COMPLETE_UNAUDITED','layout_version','K1-small-panels-v3-fixed-paper','error_curves',12,'new_acoustic_solves',0,'pixel_tolerance',0);
r1_bie_manual_io('write',fullfile(out,'generation_receipt.json'),receipt);
end
function local_close(f)
if isgraphics(f), close(f); end
end
