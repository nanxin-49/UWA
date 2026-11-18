function receipt=render_k1_pm_stable_saved(work_dirs,out)
%RENDER_K1_PM_STABLE_SAVED Fixed-position graphics, saved arrays only.
arguments
    work_dirs (1,:) cell
    out (1,:) char
end
assert(ismember(numel(work_dirs),[1 3]) && ~isfolder(out)); mkdir(out);
Q=cell(size(work_dirs)); A=Q;
for j=1:numel(Q)
    Q{j}=load(fullfile(work_dirs{j},'common.mat'),'cfg','x');
    A{j}=load(fullfile(work_dirs{j},'comparison.mat'));
    assert(A{j}.reference_qualified && numel(Q{j}.x)==411);
end
if isscalar(Q)
    for file=1:2
        specs=cell(4,1);
        for panel=1:4
            specs{panel}=struct('seed',{},'method',{},'kind',{},'x',{},'y',{},'label',{},'color',{},'style',{});
            methods=[1 3]; if file==2, methods=3; end
            kinds={'absolute','normalized','amplitude','phase'};
            for m=methods
                specs{panel}(end+1).seed=Q{1}.cfg.seed;
                specs{panel}(end).method=local_method(m);
                specs{panel}(end).kind=kinds{panel};
                specs{panel}(end).x=Q{1}.x; specs{panel}(end).y=local_value(A{1},m,kinds{panel});
                specs{panel}(end).label=local_method(m);
                specs{panel}(end).color=[0 .447 .698]; specs{panel}(end).style='--';
                if m==3, specs{panel}(end).color=[.835 .369 0]; specs{panel}(end).style='-.'; end
            end
        end
        stem='F02_pointwise_errors'; if file==2, stem='F04_FullFK_pointwise_errors'; end
        title_text=sprintf('K1 | U_PM=%g m/s, N=%d, seed %d | 0.10 m strip | 4 kHz',Q{1}.cfg.U_PM_mps,Q{1}.cfg.modes,Q{1}.cfg.seed);
        subtitle=sprintf('411 actual receivers; phase mask |BIE| >= 0.01 line RMS (%d/411)',nnz(A{1}.mask));
        local_figure(specs,out,stem,title_text,subtitle);
    end
else
    assert(isequal(cellfun(@(q)q.cfg.seed,Q),[260001 260002 260003]));
    specs=cell(4,1); colors=[0 .447 .698;.009 .55 .35;.835 .369 0]; styles={'-','--','-.'};
    for panel=1:4
        specs{panel}=struct('seed',{},'method',{},'kind',{},'x',{},'y',{},'label',{},'color',{},'style',{});
        m=1; if ismember(panel,[2 4]), m=3; end
        kind='normalized'; if panel>=3, kind='phase'; end
        for j=1:3
            specs{panel}(j)=struct('seed',Q{j}.cfg.seed,'method',local_method(m),'kind',kind, ...
                'x',Q{j}.x,'y',local_value(A{j},m,kind),'label',sprintf('%s / seed %d',local_method(m),Q{j}.cfg.seed), ...
                'color',colors(j,:),'style',styles{j});
        end
    end
    local_figure(specs,out,'F05_three_seed_errors','R1 + R2 | 0.10 m 2-D strip | 4 kHz', ...
        'Each seed uses its own C3 BIE line RMS; phase mask = 0.01 RMS; 411 actual receivers');
end
receipt=struct('status','FIXED_LAYOUT_FIGURES_GENERATED_PENDING_INDEPENDENT_AUDIT','new_acoustic_solves',0, ...
    'layout','independent pixel-position axes; fixed limits/ticks/labels; no tiledlayout auto geometry', ...
    'pixel_acceptance','unchanged: zero differing channel values after each FIG reopen', ...
    'data_transform','only original abs, whole-line RMS normalization, amplitude difference and masked wrapped phase');
r1_bie_manual_io('write',fullfile(out,'generation_receipt.json'),receipt);
end
function local_figure(specs,out,stem,title_text,subtitle)
f=figure('Visible','off','Units','pixels','Position',[40 40 1800 1400], ...
    'Color','w','Resize','off','MenuBar','none','ToolBar','none'); guard=onCleanup(@()close(f));
annotation(f,'textbox',[.07 .958 .88 .025],'String',title_text,'EdgeColor','none', ...
    'FontName','Arial','FontSize',15,'HorizontalAlignment','center','Interpreter','none');
annotation(f,'textbox',[.07 .934 .88 .022],'String',subtitle,'EdgeColor','none', ...
    'FontName','Arial','FontSize',11,'HorizontalAlignment','center','Interpreter','none');
bottoms=[1030 710 390 70];
for panel=1:4
    ax=axes('Parent',f,'Units','pixels','Position',[150 bottoms(panel) 1330 255], ...
        'PositionConstraint','innerposition','FontName','Arial','FontSize',11, ...
        'Box','on','XGrid','on','YGrid','on','LooseInset',[0 0 0 0]); hold(ax,'on');
    ax.UserData=struct('role','fixed_error_axis','panel',panel);
    all_y=[];
    for j=1:numel(specs{panel})
        s=specs{panel}(j); ln=plot(ax,s.x,s.y,s.style,'Color',s.color,'LineWidth',1.15,'DisplayName',s.label);
        ln.UserData=struct('seed',s.seed,'method',s.method,'kind',s.kind,'panel',panel);
        all_y=[all_y;s.y(isfinite(s.y))]; %#ok<AGROW>
    end
    kind=specs{panel}(1).kind;
    if ismember(kind,{'absolute','normalized'})
        upper=max(all_y)*1.10; if upper==0, upper=1; end
        limits=[0 upper];
    else
        upper=max(abs(all_y))*1.10; if upper==0, upper=1; end
        limits=[-upper upper];
    end
    ax.XLim=[min(specs{panel}(1).x) max(specs{panel}(1).x)]; ax.YLim=limits;
    ax.XTick=-40:10:40; ax.YTick=linspace(limits(1),limits(2),5);
    ax.XTickLabel=arrayfun(@(v)sprintf('%g',v),ax.XTick,'UniformOutput',false);
    ax.YTickLabel=arrayfun(@(v)sprintf('%.4g',v),ax.YTick,'UniformOutput',false);
    ax.XAxis.Exponent=0; ax.YAxis.Exponent=0;
    switch kind
        case 'absolute', label='|p_m - p_{BIE}| (Pa)';
        case 'normalized', label='|residual| / own BIE line RMS';
        case 'amplitude', label='|p_m| - |p_{BIE}| (Pa)';
        case 'phase', label='Wrapped phase difference (rad)';
    end
    ylabel(ax,label,'FontSize',12); xlabel(ax,'Receiver x (m)','FontSize',11);
    if panel==1 || numel(specs{panel})==3
        lg=legend(ax,'Location','northeast','FontSize',10,'AutoUpdate','off'); lg.Units='pixels';
        lg.Position=[1500 bottoms(panel)+175 280 72];
    end
end
drawnow;
savefig(f,fullfile(out,[stem '.fig']));
exportgraphics(f,fullfile(out,[stem '.png']),'Resolution',160,'Width',1800,'Height',1400,'Units','pixels');
end
function name=local_method(m)
name='PS'; if m==3, name='FullFK'; end
end
function v=local_value(a,m,kind)
switch kind
    case 'absolute', v=abs(a.U(:,m)-a.U(:,2));
    case 'normalized', v=abs(a.U(:,m)-a.U(:,2))/a.scale;
    case 'amplitude', v=abs(a.U(:,m))-abs(a.U(:,2));
    case 'phase', v=angle(a.U(:,m).*conj(a.U(:,2))); v(~a.mask)=NaN;
end
end
