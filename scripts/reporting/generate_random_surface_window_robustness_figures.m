function generate_random_surface_window_robustness_figures(validation)
%GENERATE_RANDOM_SURFACE_WINDOW_ROBUSTNESS_FIGURES Compact robustness plots.
v=validation; out=fullfile(v.config.output_dir,'figures');
if ~exist(out,'dir'), mkdir(out); end
t=v.comparison_table;
colors=lines(numel(v.config.windows_m));

f=figure('Visible','off','Color','w','Position',[100 100 1200 760]);
tl=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
metrics={'reflect_tl_error_db','reflect_phase_error_rad','receiver_reflected_center_l2','max_edge5_fraction'};
labels={'|reflect TL error| (dB)','|reflect phase error| (rad)','receiver center complex L2','maximum outer-5% fraction'};
limits=[v.config.tl_limit_db,v.config.phase_limit_rad,v.config.center_l2_limit,v.config.edge5_fraction_limit];
for mm=1:numel(metrics)
    nexttile; hold on; grid on
    for ww=1:numel(v.config.windows_m)
        w=v.config.windows_m(ww); q=t(t.window_nominal_m==w,:);
        x=q.target_hs_m+0.006*(ww-2);
        y=abs(q.(metrics{mm}));
        scatter(x,y,36,colors(ww,:),'filled','DisplayName',sprintf('W=%g m',w));
    end
    yline(limits(mm),'k--','gate','HandleVisibility','off');
    xlabel('target H_s (m)'); ylabel(labels{mm}); set(gca,'YScale','log');
end
legend('Location','best'); title(tl,'4 kHz random-surface window robustness');
exportgraphics(f,fullfile(out,'01_random_surface_gate_scatter.png'),'Resolution',180); close(f);

f=figure('Visible','off','Color','w','Position',[100 100 1050 650]); hold on; grid on
q=t(t.window_nominal_m~=v.config.reference_window_m,:);
gscatter(q.max_edge5_fraction,abs(q.reflect_tl_error_db), ...
    strcat(string(q.sea_state)," / W",string(q.window_nominal_m)));
xline(v.config.edge5_fraction_limit,'k--','edge gate','HandleVisibility','off');
yline(v.config.tl_limit_db,'k:','TL gate','HandleVisibility','off');
set(gca,'XScale','log','YScale','log'); xlabel('maximum outer-5% energy fraction');
ylabel('|reflect TL error| (dB)'); title('Edge validity versus receiver-response error');
exportgraphics(f,fullfile(out,'02_edge_vs_receiver_error.png'),'Resolution',180); close(f);

f=figure('Visible','off','Color','w','Position',[100 100 1100 680]);
rates=v.pass_rate_table; sea=string({v.config.sea_states.name}); data=nan(numel(sea),numel(v.config.windows_m));
for ss=1:numel(sea)
    for ww=1:numel(v.config.windows_m)
        m=strcmp(rates.sea_state,sea(ss)) & rates.window_nominal_m==v.config.windows_m(ww);
        data(ss,ww)=rates.pass_rate(m);
    end
end
b=bar(categorical(sea),data); grid on; ylim([0 1.05]); ylabel('strict pass rate');
legend(compose('W=%g m',v.config.windows_m),'Location','southoutside','Orientation','horizontal');
title('Strict pass rate by sea state');
for ii=1:numel(b), b(ii).FaceColor=colors(ii,:); end
exportgraphics(f,fullfile(out,'03_pass_rates.png'),'Resolution',180); close(f);

f=figure('Visible','off','Color','w','Position',[100 100 1150 720]); hold on; grid on
q=v.energy_table(v.energy_table.window_nominal_m==192,:);
for ss=1:numel(sea)
    s=q(strcmp(q.sea_state,sea(ss)) & q.seed==v.config.seeds(1),:);
    up=s(s.stage=="upward",:); down=s(s.stage=="downward",:);
    plot(100-up.z_m,up.edge5_fraction,'-','LineWidth',1.3,'DisplayName',sea(ss)+" upward");
    plot(100+down.z_m,down.edge5_fraction,'--','LineWidth',1.3,'DisplayName',sea(ss)+" downward");
end
yline(v.config.edge5_fraction_limit,'k:','gate','HandleVisibility','off');
set(gca,'YScale','log'); xlabel('cumulative chain coordinate (m)'); ylabel('outer-5% energy fraction');
title(sprintf('Stage attribution at W=192 m, seed %d',v.config.seeds(1))); legend('Location','best');
exportgraphics(f,fullfile(out,'04_stage_edge_energy.png'),'Resolution',180); close(f);
end
