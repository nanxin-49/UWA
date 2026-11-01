function receipt=render_deepsea3d_65grid_four_v3(project_root)
%RENDER_DEEPSEA3D_65GRID_FOUR Render four saved native points, never solve fields.
arguments
    project_root (1,1) string
end
out=fullfile(project_root,'results','validation','deepsea3d_65grid_four_20261027');
source=fullfile(out,'four_point_blind_metrics.csv');
assert(isfile(source),'deepsea3d65:MissingSavedMetrics');
t=readtable(source);
assert(isequal(t.frequency_hz,[4750;7250;4125;6562.5]));
gates=string(t.prediction_gates_pass);
assert(all(gates=="True" | gates=="true" | gates=="1"));
labels=compose('%g Hz',t.frequency_hz);
colors=[0.10 0.37 0.65;0.10 0.37 0.65;0.75 0.29 0.13;0.75 0.29 0.13];

f=figure('Visible','off','Color','w','Position',[100 100 980 430]);
tl=tiledlayout(f,1,2,'Padding','compact','TileSpacing','compact');
ax=nexttile(tl);hold(ax,'on');grid(ax,'on');
for j=1:height(t)
    plot(ax,t.R_FE_real(j)*1e5,t.R_FE_imag(j)*1e5,'o','Color',colors(j,:), ...
        'MarkerFaceColor',colors(j,:),'MarkerSize',8,'HandleVisibility','off');
    plot(ax,t.Rpred_real(j)*1e5,t.Rpred_imag(j)*1e5,'x','Color',[0 0 0], ...
        'MarkerSize',9,'LineWidth',1.4,'HandleVisibility','off');
    text(ax,t.R_FE_real(j)*1e5+0.12,t.R_FE_imag(j)*1e5,labels(j),'FontSize',9);
end
plot(ax,nan,nan,'ko','MarkerFaceColor','k','DisplayName','FE reference');
plot(ax,nan,nan,'kx','DisplayName','Locked prediction');
xlabel(ax,'Re(R)  [10^{-5} m^{-1}]');ylabel(ax,'Im(R)  [10^{-5} m^{-1}]');
title(ax,'Native reflection complex plane');legend(ax,'Location','best');
xlim(ax,[-6.2 3.8]);ylim(ax,[-5.3 5.0]);axis(ax,'equal');
ax2=nexttile(tl);hold(ax2,'on');grid(ax2,'on');
for j=1:height(t)
    plot(ax2,(t.Rpred_real(j)-t.R_FE_real(j))*1e9, ...
        (t.Rpred_imag(j)-t.R_FE_imag(j))*1e9,'o','Color',colors(j,:), ...
        'MarkerFaceColor',colors(j,:),'MarkerSize',8);
    text(ax2,(t.Rpred_real(j)-t.R_FE_real(j))*1e9+0.08, ...
        (t.Rpred_imag(j)-t.R_FE_imag(j))*1e9,labels(j),'FontSize',9);
end
xline(ax2,0,':');yline(ax2,0,':');
xlabel(ax2,'Re(R_{pred}-R_{FE})  [10^{-9} m^{-1}]');
ylabel(ax2,'Im(R_{pred}-R_{FE})  [10^{-9} m^{-1}]');
title(ax2,'Unscaled complex prediction residual');
xlim(ax2,[-2.1 3.4]);ylim(ax2,[-1.6 1.35]);
savefig(f,fullfile(out,'four_point_reflection_complex_v3.fig'));
exportgraphics(f,fullfile(out,'four_point_reflection_complex_v3.png'),'Resolution',180);
close(f);

f=figure('Visible','off','Color','w','Position',[100 100 850 460]);
ax=axes(f);hold(ax,'on');grid(ax,'on');
x=1:height(t);
semilogy(ax,x,t.absolute_complex_gap_m_inv,'o','LineStyle','none','Color',[0.07 0.34 0.64], ...
    'MarkerFaceColor',[0.07 0.34 0.64],'MarkerSize',8,'DisplayName','|Rpred-RFE|');
semilogy(ax,x,t.FE_empirical_sensitivity_m_inv,'s','LineStyle','none','Color',[0.75 0.29 0.13], ...
    'MarkerFaceColor',[0.75 0.29 0.13],'MarkerSize',8,'DisplayName','FE empirical sensitivity');
semilogy(ax,x,5*t.FE_empirical_sensitivity_m_inv,'^','LineStyle','none','Color',[0.35 0.35 0.35], ...
    'MarkerFaceColor',[0.35 0.35 0.35],'MarkerSize',7,'DisplayName','5 x sensitivity');
set(ax,'YScale','log');ylim(ax,[1e-11 1e-8]);
xticks(ax,x);xticklabels(ax,labels);xlim(ax,[.7 4.3]);
ylabel(ax,'Absolute complex difference  [m^{-1}]');
title(ax,'Four independent points: prediction gap and FE resolution');
legend(ax,'Location','best');
text(ax,1.0,max(t.absolute_complex_gap_m_inv)*1.65,'Routine','Color',[0.07 0.34 0.64]);
text(ax,3.0,max(t.absolute_complex_gap_m_inv)*1.65,'Cancellation stress','Color',[0.75 0.29 0.13]);
savefig(f,fullfile(out,'prediction_vs_FE_sensitivity_v3.fig'));
exportgraphics(f,fullfile(out,'prediction_vs_FE_sensitivity_v3.png'),'Resolution',180);
close(f);
receipt=struct('source',source,'figures',2,'points',height(t));
disp(receipt);
end
