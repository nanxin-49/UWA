function receipt=render_pm_wind_cutoff_review(source_dir, numeric_review_dir, output_dir)
%RENDER_PM_WIND_CUTOFF_REVIEW Audited saved K0 data, descriptive contrasts.
arguments
    source_dir (1,:) char
    numeric_review_dir (1,:) char
    output_dir (1,:) char
end
pm_sensitivity_io('absolute',output_dir); assert(~isfolder(output_dir)); mkdir(output_dir);
timer=tic; local=jsondecode(fileread(fullfile(numeric_review_dir,'receipt.json')));
assert(strcmp(local.status,'K0_LOCAL_SAVED_REVIEW_PASS') && local.case_count==45 && local.failed_checks==0);
csv=fullfile(source_dir,'output','summary_attempt_001','statistics.csv'); rows=readtable(csv);
geometry=readtable(fullfile(numeric_review_dir,'continuous_geometry_bounds.csv'));
fields={'raw_rms_height_m','raw_rms_slope','raw_rms_second_per_m', ...
    'embedded_rms_height_m','embedded_rms_slope','embedded_rms_second_per_m'};
% Fixed-U cutoff ratios, fixed-N wind ratios, and paired difference-in-differences.
cutoff=table; wind=table; interaction=table;
for seed=260001:260003
    for u=[5 6 8]
        for pair=[12 32]
            a=local_row(rows,seed,u,pair); b=local_row(rows,seed,u,40);
            for j=1:numel(fields)
                name=fields{j}; cutoff=[cutoff;table(seed,u,pair,40,string(name),a.(name),b.(name), ...
                    b.(name)/a.(name),b.(name)-a.(name),'VariableNames', ...
                    {'seed','U_PM_mps','N_from','N_to','response','value_from','value_to','ratio','difference'})]; %#ok<AGROW>
            end
        end
    end
    for N=[12 40]
        a=local_row(rows,seed,5,N); b=local_row(rows,seed,8,N);
        for j=1:numel(fields)
            name=fields{j}; wind=[wind;table(seed,N,string(name),a.(name),b.(name),b.(name)/a.(name), ...
                'VariableNames',{'seed','N','response','U5_value','U8_value','ratio_U8_U5'})]; %#ok<AGROW>
        end
    end
    a=local_row(rows,seed,5,12); b=local_row(rows,seed,5,40);
    c=local_row(rows,seed,8,12); d=local_row(rows,seed,8,40);
    for j=1:numel(fields)
        name=fields{j}; interaction=[interaction;table(seed,string(name), ...
            (d.(name)-c.(name))-(b.(name)-a.(name)), ...
            'VariableNames',{'seed','response','I_U8minusU5_N40minusN12'})]; %#ok<AGROW>
    end
end
writetable(cutoff,fullfile(output_dir,'fixed_U_cutoff_contrasts.csv'));
writetable(wind,fullfile(output_dir,'fixed_N_wind_contrasts.csv'));
writetable(interaction,fullfile(output_dir,'paired_interactions.csv'));
candidate=rows(rows.seed==260001 & ((rows.U_PM_mps==6 & ismember(rows.N,[7 12 40])) | ...
    (ismember(rows.U_PM_mps,[5 8]) & ismember(rows.N,[12 40]))),:);
candidate=sortrows(candidate,{'U_PM_mps','N'}); assert(height(candidate)==7);
selected=table;
for j=1:height(candidate)
    g=geometry(geometry.seed==candidate.seed(j) & geometry.U_PM_mps==candidate.U_PM_mps(j) & geometry.N==candidate.N(j),:);
    assert(height(g)==1 && g.certified_auxiliary_plane_gap_m>0 && g.certified_receiver_clearance_m>0);
    selected=[selected;[candidate(j,{'seed','U_PM_mps','N','Kmax_rad_per_m','raw_rms_height_m', ...
        'raw_rms_slope','raw_rms_second_per_m','embedded_max_abs_slope','continuous_coverage_0_Kmax'}), ...
        g(:,{'certified_auxiliary_plane_gap_m','certified_receiver_clearance_m'})]]; %#ok<AGROW>
end
writetable(selected,fullfile(output_dir,'k1_selected_geometry.csv'));
fig=figure('Visible','off','Color','w','Position',[60 60 1550 1000]); guard=onCleanup(@()close(fig));
colors=[0 .35 .65;.8 .35 0;0 .55 .35]; markers={'o','s','^'}; styles={'-','--',':'};
winds=[5 6 8];
titles={'Raw RMS elevation (m)','Raw RMS slope','Raw RMS second derivative (1/m)', ...
    'C2 RMS elevation (m)','C2 RMS slope','C2 RMS second derivative (1/m)'};
axes_list=gobjects(6,1); handles=gobjects(9,1); expected=struct('x',{},'y',{}); n=0;
for j=1:6
    axes_list(j)=axes(fig,'Position',[.07+mod(j-1,3)*.32 .64-floor((j-1)/3)*.36 .26 .25]);
    ax=axes_list(j); hold(ax,'on');
    for u=1:3
        for s=1:3
            r=rows(rows.U_PM_mps==winds(u) & rows.seed==260000+s,:);
            r=sortrows(r,'N');
            h=plot(ax,r.Kmax_rad_per_m,r.(fields{j}),'Color',colors(u,:),'Marker',markers{u}, ...
                'LineStyle',styles{s},'LineWidth',1.1,'DisplayName',sprintf('U_{PM}=%d m/s, seed%d',winds(u),260000+s));
            if j==1, handles((u-1)*3+s)=h; end
            n=n+1; expected(n)=struct('x',r.Kmax_rad_per_m(:),'y',r.(fields{j})(:));
        end
    end
    title(ax,titles{j}); xlabel(ax,'Physical Kmax (rad/m)'); grid(ax,'on');
    xlim(ax,[.24 1.61]); set(ax,'XTick',[.275 .471 .746 1.257 1.571],'FontSize',10);
end
for j=1:3
    top=max([rows.(fields{j});rows.(fields{j+3})]);
    ylim(axes_list(j),[0 1.08*top]); ylim(axes_list(j+3),[0 1.08*top]);
end
legend(axes_list(1),handles,'NumColumns',3,'Position',[.07 .035 .9 .14],'FontSize',10,'AutoUpdate','off');
annotation(fig,'textbox',[.05 .93 .9 .045],'String', ...
    'K0: 45 audited surfaces | L=160 m, 8192 samples | paired seeds; no RMS rescaling', ...
    'Interpreter','none','EdgeColor','none','HorizontalAlignment','center','FontSize',14);
drawnow; savefig(fig,fullfile(output_dir,'K0_raw_C2_statistics.fig'));
exportgraphics(fig,fullfile(output_dir,'K0_raw_C2_statistics.png'),'Resolution',130);
reopened=openfig(fullfile(output_dir,'K0_raw_C2_statistics.fig'),'new','invisible'); rg=onCleanup(@()close(reopened));
lines=findall(reopened,'Type','line'); assert(numel(lines)==numel(expected)); used=false(numel(expected),1);
for j=1:numel(lines)
    found=false;
    for k=1:numel(expected)
        if ~used(k) && isequal(lines(j).XData(:),expected(k).x) && isequal(lines(j).YData(:),expected(k).y)
            used(k)=true; found=true; break
        end
    end
    assert(found,'Saved review FIG differs from CSV');
end
receipt=struct('status','K0_SAVED_SUMMARY_PASS','failed_checks',0,'rows',45,'curves_reopened',numel(lines), ...
    'cutoff_contrast_rows',height(cutoff),'wind_contrast_rows',height(wind),'interaction_rows',height(interaction), ...
    'K1_selected_cases',height(selected),'K1_execution_enabled',false,'elapsed_s',toc(timer), ...
    'source_statistics_sha256',pm_sensitivity_io('sha',csv), ...
    'independent_review_receipt_sha256',pm_sensitivity_io('sha',fullfile(numeric_review_dir,'receipt.json')), ...
    'interpretation','Surface statistics only; no acoustic approximation qualification');
pm_sensitivity_io('write',fullfile(output_dir,'receipt.json'),receipt); disp(jsonencode(receipt));
end

function r=local_row(rows,seed,u,n)
r=rows(rows.seed==seed & rows.U_PM_mps==u & rows.N==n,:); assert(height(r)==1);
end
