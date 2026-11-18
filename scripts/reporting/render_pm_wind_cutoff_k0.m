function receipt = render_pm_wind_cutoff_k0(output_root, summary_dir, cfg)
%RENDER_PM_WIND_CUTOFF_K0 Saved K0 statistics; all seeds shown, no fitted CI.
pm_sensitivity_io('absolute',summary_dir);
rows=table;
for j=1:numel(cfg.cases)
    c=cfg.cases(j); id=sprintf('seed%d_U%d_N%02d',c.seed,c.U_PM_mps,c.N);
    cp=jsondecode(fileread(fullfile(output_root,id,'checkpoint.json')));
    assert(strcmp(cp.status,'AUDITED_PASS'));
    pm_sensitivity_io('verify',fullfile(output_root,id),cp.artifacts_sha256);
    rows=[rows;struct2table(jsondecode(fileread(fullfile(output_root,id,'metrics.json'))))]; %#ok<AGROW>
end
writetable(rows,fullfile(summary_dir,'statistics.csv'));
fig=figure('Visible','off','Color','w','Position',[60 60 1500 1000]); guard=onCleanup(@()close(fig));
fields={'raw_rms_height_m','raw_rms_slope','raw_rms_second_per_m', ...
    'embedded_rms_height_m','embedded_rms_slope','embedded_rms_second_per_m', ...
    'continuous_coverage_0_Kmax','variance_above_Kmax_m2','embedded_power_above_Kmax_m2'};
titles={'Raw RMS height (m)','Raw RMS slope','Raw RMS second derivative (1/m)', ...
    'C2 RMS height (m)','C2 RMS slope','C2 RMS second derivative (1/m)', ...
    'Continuous variance coverage [0,Kmax]','PM variance above cutoff (m^2)','C2 power above raw cutoff (m^2)'};
colors=[0 .35 .65;.8 .35 0;0 .55 .35]; markers={'o','s','^'}; styles={'-','--',':'};
for j=1:9
    ax=axes(fig,'Position',[.065+mod(j-1,3)*.325 .73-floor((j-1)/3)*.31 .265 .20]); hold(ax,'on');
    for u=1:3
        for s=1:3
            r=rows(rows.U_PM_mps==cfg.U_PM_mps(u) & rows.seed==cfg.seeds(s),:);
            if isempty(r), continue; end
            [~,ix]=sort(r.N); r=r(ix,:);
            plot(ax,r.Kmax_rad_per_m,r.(fields{j}),'Color',colors(u,:),'Marker',markers{u}, ...
                'LineStyle',styles{s},'DisplayName',sprintf('U_{PM}=%d, seed%d',cfg.U_PM_mps(u),cfg.seeds(s)));
        end
    end
    title(ax,titles{j}); xlabel(ax,'Physical Kmax (rad/m)'); grid(ax,'on');
    if j==1, legend(ax,'Location','northwest','FontSize',7); end
end
annotation(fig,'textbox',[.05 .94 .9 .035],'String', ...
    sprintf('%s | surface statistics only | full 160 m period; no RMS rescaling',cfg.kind), ...
    'EdgeColor','none','HorizontalAlignment','center','FontSize',13,'Interpreter','none');
drawnow; savefig(fig,fullfile(summary_dir,'statistics.fig')); exportgraphics(fig,fullfile(summary_dir,'statistics.png'),'Resolution',130);
% FIG data audit against each plotted saved-data line, without re-solving.
other=openfig(fullfile(summary_dir,'statistics.fig'),'invisible'); other_guard=onCleanup(@()close(other));
ax=findall(fig,'Type','axes'); bx=findall(other,'Type','axes'); assert(numel(ax)==numel(bx)); curves=0;
for j=1:numel(ax)
    a=findall(ax(j),'Type','line'); b=findall(bx(j),'Type','line'); assert(numel(a)==numel(b));
    for k=1:numel(a)
        assert(isequal(a(k).XData,b(k).XData) && isequal(a(k).YData,b(k).YData)); curves=curves+1;
    end
end
clear other_guard guard
fig2=figure('Visible','off','Color','w','Position',[60 60 1500 1000]); guard2=onCleanup(@()close(fig2));
all_curves=struct('x',{},'y',{}); line_index=0;
for u=1:3
    a1=axes(fig2,'Position',[.07 .72-(u-1)*.30 .40 .20]); hold(a1,'on');
    a2=axes(fig2,'Position',[.57 .72-(u-1)*.30 .39 .20]); hold(a2,'on');
    for jj=1:3
        nn=[7 12 40]; N=nn(jj);
        ix=find([cfg.cases.seed]==260001 & [cfg.cases.U_PM_mps]==cfg.U_PM_mps(u) & [cfg.cases.N]==N);
        if isempty(ix), continue; end
        id=sprintf('seed260001_U%d_N%02d',cfg.U_PM_mps(u),N);
        s=readtable(fullfile(output_root,id,'surface.csv'));
        p=readtable(fullfile(output_root,id,'spectra.csv'));
        % Solid is raw; dashed is C2. Colors/markers identify physical N.
        for embed=0:1
            if embed==0
                yy=s.raw_elevation_m; power=p.raw_bin_power_m2; style='-';
            else
                yy=s.embedded_elevation_m; power=p.embedded_bin_power_m2; style='--';
            end
            plot(a1,s.x_m,yy,'Color',colors(jj,:),'LineStyle',style, ...
                'DisplayName',sprintf('N%d %s',N,local_label(embed)));
            line_index=line_index+1; all_curves(line_index)=struct('x',s.x_m.','y',yy.');
            visible=p.k_rad_per_m<=3 & power>1e-28;
            xx=p.k_rad_per_m(visible); yy=power(visible)/ (2*pi/160);
            plot(a2,xx,yy,'Color',colors(jj,:),'LineStyle',style, ...
                'Marker',markers{jj},'MarkerIndices',1:8:numel(xx),'MarkerSize',3);
            line_index=line_index+1; all_curves(line_index)=struct('x',xx.','y',yy.');
        end
    end
    title(a1,sprintf('seed260001 U_{PM}=%d m/s: elevation',cfg.U_PM_mps(u)));
    xlabel(a1,'x (m)'); ylabel(a1,'eta (m), positive upward'); grid(a1,'on'); xlim(a1,[-80 80]);
    title(a2,'Raw/C2 one-sided sampled PSD; values <=1e-28 m^2 omitted');
    set(a2,'YScale','log'); xlabel(a2,'K (rad/m)'); ylabel(a2,'PSD (m^3)'); grid(a2,'on');
    if u==1 && ~isempty(findall(a1,'Type','line')), legend(a1,'Location','northoutside','Orientation','horizontal','FontSize',8); end
end
annotation(fig2,'textbox',[.05 .94 .9 .035],'String', ...
    sprintf('%s | Raw/C2 illustration; C2 is not bandlimited | fixed dx=160/8192 m',cfg.kind), ...
    'EdgeColor','none','HorizontalAlignment','center','FontSize',13,'Interpreter','none');
drawnow; savefig(fig2,fullfile(summary_dir,'surfaces_spectra.fig'));
exportgraphics(fig2,fullfile(summary_dir,'surfaces_spectra.png'),'Resolution',130);
reopened=openfig(fullfile(summary_dir,'surfaces_spectra.fig'),'invisible'); rg=onCleanup(@()close(reopened));
lines=findall(reopened,'Type','line'); assert(numel(lines)==numel(all_curves));
for j=1:numel(lines)
    matched=false;
    for kk=1:numel(all_curves)
        if isequal(lines(j).XData,all_curves(kk).x) && isequal(lines(j).YData,all_curves(kk).y), matched=true; break; end
    end
    assert(matched,'Reopened surface/spectrum curve differs from saved data');
end
receipt=struct('status','SUMMARY_PASS','failed_checks',0,'rows',height(rows),'reopened_curves',curves, ...
    'surface_spectrum_reopened_curves',numel(lines), ...
    'normalization','full-period RMS; no recentering; no height renormalization', ...
    'replication','3 paired seed blocks; 45 configurations are not 45 independent realizations');
pm_sensitivity_io('write',fullfile(summary_dir,'render_receipt.json'),receipt);
end

function label=local_label(embed)
if embed, label='C2'; else, label='raw'; end
end
