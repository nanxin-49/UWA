function receipt = audit_piston_p0p1_figures_remote(task_root)
%AUDIT_PISTON_P0P1_FIGURES_REMOTE Independently match saved FIG curves to MAT.
arguments
    task_root (1,:) char
end
run='20261112_remote_p0p1_01';
vr=fullfile(task_root,'results','validation','piston_source',run);
cr=fullfile(task_root,'results','comparisons','piston_2d_single_frequency',run);
dir_fig=fullfile(cr,'figures'); q=load(fullfile(cr,'common.mat'));
out=fullfile(cr,'figure_audit'); assert(~isfolder(out)); mkdir(out);
rows=[]; names={};
gate_names=struct('free','gate_cont_free_r2','flat','gate_cont_flat_r2','weak','gate_cont_weak_r2');
if isfile(fullfile(cr,'stage_gate_selection.json'))
    gate_names=jsondecode(fileread(fullfile(cr,'stage_gate_selection.json')));
end
for stage={'free','flat','weak'}
    name=stage{1}; summary=fullfile(cr,'jobs',gate_names.(name),'stage_summary.mat');
    if ~isfile(summary), continue; end
    d=load(summary,'selected'); selected=d.selected;
    if strcmp(name,'free')
        stems={'F02_free_field'};
    elseif strcmp(name,'flat')
        stems={'F03_flat_direct_reflect_total','F04_flat_field_errors'};
    else
        stems={'F05_weak_reflected_field'};
    end
    for si=1:numel(stems)
        stem=stems{si}; expected={};
        for wi=1:2
            for ci=1:3
                if strcmp(name,'free')
                    reference=q.p_direct(:,wi);
                elseif strcmp(name,'flat')
                    if ci==1
                        reference=q.p_direct(:,wi);
                    elseif ci==2
                        reference=q.p_mirror(:,wi);
                    else
                        reference=q.p_direct(:,wi)+q.p_mirror(:,wi);
                    end
                else
                    reference=selected{2}.reflect(:,wi);
                end
                for mi=1:3
                    if strcmp(name,'free')
                        value=selected{mi}.direct(:,wi);
                        if ci==1
                            y=abs(value);
                        elseif ci==2
                            y=abs(value-reference);
                        else
                            y=local_phase(value,reference);
                        end
                    elseif strcmp(name,'flat')
                        fields={'direct','reflect','total'}; matrix=selected{mi}.(fields{ci}); value=matrix(:,wi);
                        if si==1, y=abs(value); else, y=abs(value-reference); end
                    else
                        value=selected{mi}.reflect(:,wi);
                        if ci==1
                            y=abs(value);
                        elseif ci==2
                            y=local_phase(value,reference);
                        else
                            y=abs(value-reference);
                        end
                    end
                    expected{end+1}=y; %#ok<AGROW>
                end
                if (strcmp(name,'free') && ci==1) || (strcmp(name,'flat') && si==1)
                    expected{end+1}=abs(reference); %#ok<AGROW>
                end
            end
        end
        fig=openfig(fullfile(dir_fig,[stem '.fig']),'new','invisible');
        lines=findobj(fig,'Type','line'); used=false(size(lines)); missing=0;
        for j=1:numel(expected)
            match=0;
            for k=1:numel(lines)
                if ~used(k) && isequal(lines(k).XData(:),q.cfg.rx_x_m) && isequaln(lines(k).YData(:),expected{j})
                    match=k; break;
                end
            end
            if match==0, missing=missing+1; else, used(match)=true; end
        end
        info=imfinfo(fullfile(dir_fig,[stem '.png']));
        rows=[rows;numel(expected) numel(lines) missing info.Width info.Height]; %#ok<AGROW>
        names{end+1}=stem; %#ok<AGROW>
        if strcmp(name,'weak')
            drawnow; savefig(fig,fullfile(out,'weak_reopened.fig'));
            exportgraphics(fig,fullfile(out,'weak_reopened.png'),'Resolution',150);
        end
        close(fig);
    end
end
fig=openfig(fullfile(dir_fig,'F01_P0_launch_gradient.fig'),'new','invisible');
ax=findobj(fig,'Type','axes'); source_axes_log=all(strcmp({ax.YScale},'log'));
close(fig);
T=array2table(rows,'VariableNames',{'expected_curves','saved_curves','missing_curves','PNG_width_px','PNG_height_px'});
T.figure_stem=names.'; writetable(T,fullfile(out,'curve_reconstruction.csv'));
receipt=struct('status','FIG_CURVE_RECONSTRUCTION_COMPLETE','figures',height(T), ...
    'curves_verified',sum(T.expected_curves),'missing_curves',sum(T.missing_curves), ...
    'source_gradient_axes_log',source_axes_log,'failed_checks',sum(T.missing_curves)+~source_axes_log, ...
    'PNG_visual_review','pending local actual image inspection');
fid=fopen(fullfile(out,'receipt.json'),'w'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid);
fprintf('%s\n',jsonencode(receipt)); assert(receipt.failed_checks==0);
end
function y=local_phase(p,r)
y=NaN(size(p)); valid=abs(r)>=1e-4*max(abs(r)); y(valid)=angle(p(valid).*conj(r(valid)))*180/pi;
end
