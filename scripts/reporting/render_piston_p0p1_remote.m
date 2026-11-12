function receipt = render_piston_p0p1_remote(task_root,figure_directory)
%RENDER_PISTON_P0P1_REMOTE Plot saved remote evidence, never launch a solver.
% Run in the licensed remote MATLAB only. Preserve PNG and editable FIG.
arguments
    task_root (1,:) char
    figure_directory (1,:) char = 'figures'
end
run='20261112_remote_p0p1_01';
vr=fullfile(task_root,'results','validation','piston_source',run);
cr=fullfile(task_root,'results','comparisons','piston_2d_single_frequency',run);
if isfolder(cr), out=fullfile(cr,figure_directory); else, out=fullfile(vr,figure_directory); end
assert(~isfolder(out),'Use a fresh figure output directory.'); mkdir(out);
t0=tic; stems={}; stages={}; colors=[0 .35 .65;.8 .32 .05;.1 .55 .35]; styles={'-','--','-.'}; labels={'PE Model-0','BIE','FullFK'};
source_table=readtable(fullfile(vr,'jobs','p0_pe','metrics.csv'));
fig=figure('Visible','off','Color','w','Position',[40 40 1100 530]);
tl=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
title(tl,'P0: actual PE launch representation | each aperture checked separately','FontSize',14);
refined_job='p0_pe_refined';
selection_file=fullfile(vr,'p0_recovery_selection.json');
if isfile(fullfile(vr,'p0_active_selection_v5.json'))
    selection_file=fullfile(vr,'p0_active_selection_v5.json');
end
if isfile(selection_file)
    selection=jsondecode(fileread(selection_file)); refined_job=selection.active_refined_job_id;
end
refined_file=fullfile(vr,'jobs',refined_job,'metrics.csv');
for wi=1:2
    ax=nexttile(tl); hold(ax,'on');
    set(ax,'YScale','log');
    W=[.2 .1]; mask=abs(source_table.total_width_m-W(wi))<1e-12;
    semilogy(ax,source_table.frequency_hz(mask)/1000,source_table.actual_gradient_relative(mask),'o-','LineWidth',1.8,'DisplayName','512 m launch window');
    recovery_file=fullfile(vr,'jobs','p0_pe_refined_recovery01','metrics.csv');
    if isfile(recovery_file)
        recovery=readtable(recovery_file);
        recovery_mask=abs(recovery.total_width_m-W(wi))<1e-12;
        semilogy(ax,recovery.frequency_hz(recovery_mask)/1000,recovery.actual_gradient_relative(recovery_mask), ...
            'd-.','LineWidth',1.8,'DisplayName','8192 m launch window');
    end
    if isfile(refined_file)
        refined=readtable(refined_file); mask=abs(refined.total_width_m-W(wi))<1e-12;
        semilogy(ax,refined.frequency_hz(mask)/1000,refined.actual_gradient_relative(mask),'s--','LineWidth',1.8,'DisplayName','32768 m launch window');
    end
    h=yline(ax,1e-4,'k:','Frozen gradient gate'); h.HandleVisibility='off';
    title(ax,sprintf('Total aperture width %g cm',100*W(wi))); grid(ax,'on');
    xlabel(ax,'Frequency (kHz)'); ylabel(ax,'Relative pressure-gradient error');
    xticks(ax,[4 6 8]); ylim(ax,[1e-6 3e-3]); local_legend(ax);
end
local_save(fig,out,'F01_P0_launch_gradient'); stems{end+1}='F01_P0_launch_gradient';
q=load(fullfile(vr,'common.mat'));
gate_names=struct('free','gate_cont_free_r2','flat','gate_cont_flat_r2','weak','gate_cont_weak_r2');
if isfile(fullfile(cr,'stage_gate_selection.json'))
    gate_names=jsondecode(fileread(fullfile(cr,'stage_gate_selection.json')));
end
for si=1:3
    all_stages={'free','flat','weak'}; stage=all_stages{si};
    summary=fullfile(cr,'jobs',gate_names.(stage),'stage_summary.mat');
    if ~isfile(summary), continue; end
    data=load(summary,'selected','checks'); selected=data.selected;
    common_limits=[56 56 110]; component_names={'direct','reflect','total'};
    for mi=1:3
        for ci=1:3
            common_limits(ci)=max(common_limits(ci),1.05*max(abs(selected{mi}.(component_names{ci})),[],'all'));
        end
    end
    passed=all(cell2mat(struct2cell(data.checks)));
    if passed, status='stage gate passed'; else, status='INCONCLUSIVE: gate not passed'; end
    stages{end+1}=stage; %#ok<AGROW>
    if strcmp(stage,'free')
        fig=figure('Visible','off','Color','w','Position',[40 40 1420 760]);
        tl=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
        title(tl,sprintf('4 kHz free field | z=3 m | peak Pa | %s',status),'FontSize',14);
        for wi=1:2
            reference=q.p_direct(:,wi);
            for col=1:3
                ax=nexttile(tl,(wi-1)*3+col); hold(ax,'on');
                for mi=1:3
                    value=selected{mi}.direct(:,wi);
                    if col==1, yy=abs(value);
                    elseif col==2, yy=abs(value-reference);
                    else, yy=local_phase(value,reference); end
                    plot(ax,q.cfg.rx_x_m,yy,styles{mi},'Color',colors(mi,:),'LineWidth',1.6,'DisplayName',labels{mi});
                end
                if col==1
                    plot(ax,q.cfg.rx_x_m,abs(reference),'k:','LineWidth',1.5,'DisplayName','Aperture reference');
                    ylabel(ax,'Peak |p_D| (Pa)'); ylim(ax,[0 common_limits(1)]);
                elseif col==2, ylabel(ax,'Complex absolute error (Pa)');
                else, ylabel(ax,'Phase error (deg)'); end
                title(ax,sprintf('W=%g cm',100*q.cfg.widths_m(wi)));
                local_axes(ax); if wi==1 && col==1, local_legend(ax); end
            end
        end
        local_save(fig,out,'F02_free_field'); stems{end+1}='F02_free_field'; %#ok<AGROW>
    elseif strcmp(stage,'flat')
        components={'direct','reflect','total'};
        for error_plot=[false true]
            fig=figure('Visible','off','Color','w','Position',[40 40 1420 760]);
            tl=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
            if error_plot, description='complex absolute error'; else, description='absolute peak pressure'; end
            title(tl,sprintf('4 kHz pressure-release plane | %s | %s',description,status),'FontSize',14);
            for wi=1:2
                for ci=1:3
                    name=components{ci}; ax=nexttile(tl,(wi-1)*3+ci); hold(ax,'on');
                    if ci==1, reference=q.p_direct(:,wi);
                    elseif ci==2, reference=q.p_mirror(:,wi);
                    else, reference=q.p_direct(:,wi)+q.p_mirror(:,wi); end
                    for mi=1:3
                        matrix=selected{mi}.(name); value=matrix(:,wi);
                        if error_plot, yy=abs(value-reference); else, yy=abs(value); end
                        plot(ax,q.cfg.rx_x_m,yy,styles{mi},'Color',colors(mi,:),'LineWidth',1.6,'DisplayName',labels{mi});
                    end
                    if ~error_plot
                        plot(ax,q.cfg.rx_x_m,abs(reference),'k:','LineWidth',1.4,'DisplayName','Independent mirror');
                        ylim(ax,[0 common_limits(ci)]);
                    end
                    title(ax,sprintf('W=%g cm | %s',100*q.cfg.widths_m(wi),name));
                    if error_plot, ylabel(ax,'Complex absolute error (Pa)'); else, ylabel(ax,'Peak magnitude (Pa)'); end
                    local_axes(ax); if wi==1 && ci==1, local_legend(ax); end
                end
            end
            if error_plot, stem='F04_flat_field_errors'; else, stem='F03_flat_direct_reflect_total'; end
            local_save(fig,out,stem); stems{end+1}=stem; %#ok<AGROW>
        end
    else
        fig=figure('Visible','off','Color','w','Position',[40 40 1420 760]);
        tl=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
        title(tl,sprintf('4 kHz weak C2 sinusoid | reflected field, peak Pa | %s',status),'FontSize',14);
        for wi=1:2
            reference=selected{2}.reflect(:,wi);
            for ci=1:3
                ax=nexttile(tl,(wi-1)*3+ci); hold(ax,'on');
                for mi=1:3
                    value=selected{mi}.reflect(:,wi);
                    if ci==1, yy=abs(value);
                    elseif ci==2, yy=local_phase(value,reference);
                    else, yy=abs(value-reference); end
                    plot(ax,q.cfg.rx_x_m,yy,styles{mi},'Color',colors(mi,:),'LineWidth',1.7,'DisplayName',labels{mi});
                end
                title(ax,sprintf('W=%g cm',100*q.cfg.widths_m(wi))); local_axes(ax);
                if ci==1, ylabel(ax,'Peak |p_R| (Pa)'); ylim(ax,[0 common_limits(2)]);
                elseif ci==2, ylabel(ax,'Phase difference from BIE (deg)');
                else, ylabel(ax,'Complex difference from BIE (Pa)'); end
                if wi==1 && ci==1, local_legend(ax); end
            end
        end
        local_save(fig,out,'F05_weak_reflected_field'); stems{end+1}='F05_weak_reflected_field'; %#ok<AGROW>
    end
end
receipt=struct('status','SAVED_DATA_FIGURES_EXPORTED_NOT_PHYSICAL_QUALIFICATION', ...
    'elapsed_s',toc(t0),'stages_present',{stages},'stems',{stems},'png_count',numel(stems), ...
    'fig_count',numel(stems),'raw_normalization','none; peak Pa at fixed velocity, phase difference explicitly referenced', ...
    'weak_field_policy','phase only where reference magnitude >=1e-4 of its own line maximum; gaps retained', ...
    'unrun_checks',{{'independent FIG reconstruction','actual delivered PNG visual review'}});
fid=fopen(fullfile(out,'render_receipt.json'),'w'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid);
fprintf('%s\n',jsonencode(receipt));
end
function phase=local_phase(value,reference)
phase=NaN(size(value)); valid=abs(reference)>=1e-4*max(abs(reference));
phase(valid)=angle(value(valid).*conj(reference(valid)))*180/pi;
end
function local_axes(ax)
grid(ax,'on'); xlim(ax,[-40 40]); xlabel(ax,'Receiver x (m), z=3 m'); set(ax,'FontSize',10);
end
function local_save(fig,out,stem)
assert(~isfile(fullfile(out,[stem '.fig'])) && ~isfile(fullfile(out,[stem '.png'])));
drawnow; savefig(fig,fullfile(out,[stem '.fig']));
exportgraphics(fig,fullfile(out,[stem '.png']),'Resolution',150); close(fig);
end
function local_legend(ax)
curves=findobj(ax,'Type','line'); curves=flipud(curves);
names={curves.DisplayName};
lgd=legend(ax,curves,names,'Location','southwest','AutoUpdate','off');
set(lgd,'Visible','on','FontSize',9,'Color','white');
end
