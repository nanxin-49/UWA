function output=render_two_dimensional_s5a_weak_bie_8k(run_dir)
%RENDER_TWO_DIMENSIONAL_S5A_WEAK_BIE_8K Plot the audited 8 kHz A-stage pilot.
arguments
    run_dir (1,:) char
end
audit_path=fullfile(run_dir,'audit_v1','independent_audit.json');
layers_path=fullfile(run_dir,'audit_v1','bie_layers.csv');
out_dir=fullfile(run_dir,'figures_v1');
assert(isfile(audit_path) && isfile(layers_path));
assert(~isfolder(out_dir));
a=jsondecode(fileread(audit_path)); layers=readtable(layers_path);
assert(a.frequency_hz==8000 && a.amplitude_m==.01 && ...
    a.receiver_x_m==0 && a.receiver_z_m==3);
assert(a.reference_gate_1e_4_pass && ...
    ~a.pe_difference_over_5x_resolution && ...
    ~a.fk_difference_over_5x_resolution);
assert(height(layers)==5 && all(layers.node_count==[7680;7680;8960;6400;7472]));
idx=~strcmp(layers.job_id,'main96');
model_ratio=[a.pe_model_relative/(5*a.pe_empirical_resolution), ...
    a.fk_model_relative/(5*a.fk_empirical_resolution)];
assert(all(isfinite(model_ratio)) && all(model_ratio>=0));
mkdir(out_dir);
fig=figure('Color','w','Units','inches','Position',[1 1 10 7], ...
    'Visible','off','Name','8 kHz BIE sensitivity and model resolution');
layout=tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
ax1=nexttile(layout,1);
y=layers.difference_from_main_relative(idx);
semilogy(ax1,1:numel(y),y,'o','MarkerSize',7,'LineWidth',1.4, ...
    'Color',[0.00 0.35 0.65],'MarkerFaceColor',[0.55 0.78 0.92]);
hold(ax1,'on'); yline(ax1,1e-4,'--','Reference gate 10^{-4}', ...
    'Color',[0.72 0.18 0.12],'LineWidth',1.2, ...
    'LabelHorizontalAlignment','left');
hold(ax1,'off'); grid(ax1,'on'); box(ax1,'on');
ax1.XLim=[.5 4.5]; ax1.XTick=1:4;
ax1.XTickLabel=cellstr(layers.job_id(idx));
ax1.YLim=[1e-10 1e-4]; ax1.FontSize=10;
ylabel(ax1,'Complex relative difference from main96');
title(ax1,'Five-layer BIE sensitivity (main96 is the reference)');

ax2=nexttile(layout,2);
b=bar(ax2,1:2,model_ratio,.62,'LineWidth',1.0);
b.FaceColor='flat'; b.CData=[0.00 0.40 0.67;0.50 0.30 0.62];
hold(ax2,'on'); yline(ax2,1,'--','5× resolution boundary', ...
    'Color',[0.72 0.18 0.12],'LineWidth',1.2, ...
    'LabelHorizontalAlignment','left');
hold(ax2,'off'); grid(ax2,'on'); box(ax2,'on');
ax2.XLim=[.4 2.6]; ax2.XTick=1:2;
ax2.XTickLabel={'PE vs BIE','FullFK vs BIE'};
ax2.YLim=[0 1.25]; ax2.FontSize=10;
ylabel(ax2,'Model difference / (5 × empirical resolution)');
title(ax2,'Both ratios below 1: differences are unresolved');
for ii=1:2
    text(ax2,ii,model_ratio(ii)+.035,sprintf('%.3f',model_ratio(ii)), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'FontSize',10,'FontWeight','bold');
end
title(layout,'Fixed 2-D weak sine, P22 (0, 3 m), 8 kHz; unit-source field');

drawnow;
fig_path=fullfile(out_dir,'bie_sensitivity_and_model_resolution.fig');
png_path=fullfile(out_dir,'bie_sensitivity_and_model_resolution.png');
assert(~isfile(fig_path) && ~isfile(png_path));
savefig(fig,fig_path);
exportgraphics(fig,png_path,'Resolution',200);
close(fig);
output=struct('status','FIGURE_SAVED','figure_path',fig_path, ...
    'png_path',png_path,'data_source',audit_path, ...
    'transformation','direct plot of saved BIE layer differences and model differences divided by five times empirical resolution; no interpolation or smoothing', ...
    'resolution_dpi',200,'scope','single 8 kHz weak-sine P22 pilot only');
fid=fopen(fullfile(out_dir,'figure_metadata.json'),'w');
assert(fid>0); cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(output,PrettyPrint=true));
disp(output);
end
