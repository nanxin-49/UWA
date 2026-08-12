function atlas = generate_pe_bellhop_comparison_atlas_vertical(validation_file)
%GENERATE_PE_BELLHOP_COMPARISON_ATLAS_VERTICAL Ten-figure formal atlas.

if nargin<1 || isempty(validation_file), error('validation_file is required.'); end
loaded=load(validation_file,'validation');
if ~isfield(loaded,'validation'), error('MAT must contain validation.'); end
v=loaded.validation; out=v.validation_run_meta.visualization_output_dir;
if ~exist(out,'dir'), mkdir(out); end
names={ ...
    '01_geometry_and_reference_paths.png';'02_phase_reference_and_cir.png'; ...
    '03_path_pdp_comparison.png';'04_delay_residuals_and_small_offset.png'; ...
    '05_path_tl_and_scale_residuals.png';'06_source_directivity_audit.png'; ...
    '07_window_sponge_decoupling.png';'08_bellhop_ray_and_tl_fields.png'; ...
    '09_beam_convergence_and_receiver_slice.png';'10_layered_acceptance_summary.png'};
files=cellfun(@(n)fullfile(out,n),names,'UniformOutput',false);
local_geometry(v,files{1}); local_phase(v,files{2}); local_pdp(v,files{3});
local_delay(v,files{4}); local_tl(v,files{5}); local_source(v,files{6});
local_window(v,files{7}); local_fields(v,files{8}); local_beams(v,files{9});
local_summary(v,files{10});

atlas=struct('schema_version','1.0.0','run_id',v.validation_run_meta.run_id, ...
    'validation_file',validation_file,'figure_files',{files}, ...
    'normalization_policy','Shared reference per comparison; no per-curve normalization.', ...
    'difference_policy','Residual panels show zero-centered signed differences and thresholds.');
atlas.mat_file=fullfile(out,'pe_bellhop_comparison_atlas_data.mat');
atlas.manifest_file=fullfile(out,'pe_bellhop_comparison_atlas_manifest.md');
atlas.summary_file=fullfile(out,'pe_bellhop_comparison_atlas_summary.txt');
atlas_data=struct('schema_version','1.0.0','run_id',atlas.run_id, ...
    'path_table',v.path_table,'small_offset_table',v.small_offset.table, ...
    'source_audit_table',v.source_audit.table, ...
    'window_sponge_table',v.window_audit.table,'acceptance_checks',v.checks, ...
    'amplitude_metrics',v.amplitude_metrics,'phase_closure',v.closure, ...
    'field_convergence_table',v.field_output.convergence_table, ...
    'field_receiver_table',v.field_output.receiver_table);
save(atlas.mat_file,'atlas','atlas_data');
local_manifest(atlas.manifest_file,v,atlas,names);
local_atlas_summary(atlas.summary_file,v,atlas);
for ii=1:numel(files)
    if exist(files{ii},'file')~=2 || dir(files{ii}).bytes==0
        error('Atlas figure is missing or empty: %s',files{ii});
    end
end
end

function local_geometry(v,file)
cfg=v.config; fig=local_fig([100 100 1200 650]); cleanup=onCleanup(@()close(fig));
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; hold on; x=cfg.baseline_offsets_m; colors=lines(numel(x));
for ii=1:numel(x)
    plot([0 x(ii)],[cfg.z_tx_m cfg.z_rx_m],'-','Color',colors(ii,:),'LineWidth',1.8);
    xb=x(ii)*cfg.z_tx_m/(cfg.z_tx_m+cfg.z_rx_m);
    plot([0 xb x(ii)],[cfg.z_tx_m 0 cfg.z_rx_m],'--','Color',colors(ii,:),'LineWidth',1.8);
end
plot(0,cfg.z_tx_m,'kp','MarkerFaceColor',[0.9 0.25 0.1],'MarkerSize',12);
plot(x,cfg.z_rx_m*ones(size(x)),'ks','MarkerFaceColor',[1 0.8 0.1]);
yline(0,'k-','Surface'); yline(cfg.water_depth_m,'k:','Seabed excluded');
set(gca,'YDir','reverse'); xlim([0 10]); ylim([0 cfg.water_depth_m]); grid on;
xlabel('horizontal offset (m)'); ylabel('depth (m)'); title('Direct and image-source paths');
nexttile; hold on; rectangle('Position',[-cfg.width_m/2 0 cfg.width_m cfg.z_tx_m], ...
    'EdgeColor',[0.2 0.4 0.8],'LineWidth',1.8);
sp=cfg.sponge_ratio*cfg.width_m;
patch([-cfg.width_m/2 -cfg.width_m/2+sp -cfg.width_m/2+sp -cfg.width_m/2], ...
    [0 0 cfg.z_tx_m cfg.z_tx_m],[0.8 0.3 0.2],'FaceAlpha',0.2,'EdgeColor','none');
patch([cfg.width_m/2-sp cfg.width_m/2 cfg.width_m/2 cfg.width_m/2-sp], ...
    [0 0 cfg.z_tx_m cfg.z_tx_m],[0.8 0.3 0.2],'FaceAlpha',0.2,'EdgeColor','none');
plot(0,cfg.z_tx_m,'kp','MarkerFaceColor',[0.9 0.25 0.1],'MarkerSize',12);
plot(cfg.baseline_offsets_m,cfg.z_rx_m*ones(size(cfg.baseline_offsets_m)),'ks','MarkerFaceColor',[1 0.8 0.1]);
set(gca,'YDir','reverse'); axis tight; grid on; xlabel('PE x window (m)'); ylabel('depth (m)');
title(sprintf('PE %.0f m window; sponge %.2f m per edge',cfg.width_m,sp));
sgtitle(sprintf('Flat uniform validation geometry; Bellhop fan %.2f to %.1f deg', ...
    cfg.angle_limits_deg(1),cfg.angle_limits_deg(2)));
local_export(fig,file); clear cleanup
end

function local_phase(v,file)
c=v.baseline_cases(2); f=c.pe.f_axis(:)/1000; fig=local_fig([100 100 1250 820]);
cleanup=onCleanup(@()close(fig)); tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile; plot(f,unwrap(angle(c.pe.H_direct_reduced_f)),'-','LineWidth',1.4); hold on;
plot(f,unwrap(angle(c.pe.H_direct_physical_f)),'--','LineWidth',1.4);
plot(f,unwrap(angle(c.pe.H_direct_f)),':','LineWidth',1.8); grid on;
xlabel('frequency (kHz)'); ylabel('unwrapped phase (rad)'); title('Direct representations');
legend('reduced','physical','direct-DSP','Location','best');
nexttile; plot(f,unwrap(angle(c.pe.H_reflect_reduced_f)),'-','LineWidth',1.4); hold on;
plot(f,unwrap(angle(c.pe.H_reflect_physical_f)),'--','LineWidth',1.4);
plot(f,unwrap(angle(c.pe.H_reflect_f)),':','LineWidth',1.8); grid on;
xlabel('frequency (kHz)'); ylabel('unwrapped phase (rad)'); title('Reflected representations');
legend('reduced','physical','direct-DSP','Location','best');
nexttile([1 2]); ref=max([c.direct_cir.power;c.reflection_cir.power;c.bellhop_total_cir.power]);
plot(1000*c.direct_cir.delay_s,10*log10(max(c.direct_cir.power/ref,1e-10)),'-','LineWidth',1.5); hold on;
plot(1000*c.reflection_cir.delay_s,10*log10(max(c.reflection_cir.power/ref,1e-10)),'-','LineWidth',1.5);
plot(1000*c.bellhop_total_cir.delay_s,10*log10(max(c.bellhop_total_cir.power/ref,1e-10)),'k--','LineWidth',1.4);
xline(1000*c.geometry.direct_time_s,':'); xline(1000*c.geometry.surface_time_s,':');
xlim([45 62]); ylim([-80 2]); grid on; xlabel('physical delay (ms)'); ylabel('power / common reference (dB)');
title('x=6 m physical CIR; dotted lines are analytic delays'); legend('PE direct','PE reflected','Bellhop synthesis');
sgtitle('Reduced, physical and direct-DSP phase references'); local_export(fig,file); clear cleanup
end

function local_pdp(v,file)
fig=local_fig([100 100 1300 760]); cleanup=onCleanup(@()close(fig));
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
for ii=1:3
    c=v.baseline_cases(ii); common=max([c.pe_total_cir.power;c.bellhop_total_cir.power]);
    nexttile(ii); plot(1000*c.pe_total_cir.delay_s,10*log10(max(c.pe_total_cir.power/common,1e-8)),'-','LineWidth',1.4); hold on;
    plot(1000*c.bellhop_total_cir.delay_s,10*log10(max(c.bellhop_total_cir.power/common,1e-8)),'--','LineWidth',1.4);
    xlim([45 62]); ylim([-70 2]); grid on; title(sprintf('Total, x=%.0f m',c.geometry.offset_m));
    if ii==1,ylabel('common-reference power (dB)');end
    nexttile(ii+3); path_common=max([c.direct_cir.power;c.reflection_cir.power]);
    plot(1000*c.direct_cir.delay_s,10*log10(max(c.direct_cir.power/path_common,1e-8)),'-','LineWidth',1.3); hold on;
    plot(1000*c.reflection_cir.delay_s,10*log10(max(c.reflection_cir.power/path_common,1e-8)),'-','LineWidth',1.3);
    xline(1000*c.geometry.direct_time_s,':'); xline(1000*c.geometry.surface_time_s,':');
    xlim([45 62]); ylim([-70 2]); grid on; xlabel('delay (ms)'); title('PE separated paths');
    if ii==1,ylabel('common-reference power (dB)');end
end
axes_handles=findobj(fig,'Type','axes');
top_axes=axes_handles(arrayfun(@(a)contains(string(a.Title.String),'Total'),axes_handles));
bottom_axes=axes_handles(arrayfun(@(a)contains(string(a.Title.String),'PE separated'),axes_handles));
if ~isempty(top_axes), legend(top_axes(1),'PE total','Bellhop total','Location','best'); end
if ~isempty(bottom_axes), legend(bottom_axes(1),'PE direct','PE reflected','Location','best'); end
sgtitle('PE/Bellhop PDP with one reference per geometry');
local_export(fig,file); clear cleanup
end

function local_delay(v,file)
t=v.path_table; s=v.small_offset.table; fig=local_fig([100 100 1200 760]); cleanup=onCleanup(@()close(fig));
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile; hold on; m=t.path_name=="direct"; plot(t.offset_m(m),t.pe_minus_analytic_ms(m),'o-','LineWidth',1.4);
plot(t.offset_m(m),t.bellhop_minus_analytic_ms(m),'s--','LineWidth',1.4); m=t.path_name=="surface_reflection";
plot(t.offset_m(m),t.pe_minus_analytic_ms(m),'o-','LineWidth',1.4); plot(t.offset_m(m),t.bellhop_minus_analytic_ms(m),'s--','LineWidth',1.4);
yline(0.25,'r:');yline(-0.25,'r:'); grid on; xlabel('offset (m)');ylabel('delay residual (ms)');title('Baseline residuals');
legend('PE direct','BH direct','PE reflected','BH reflected','Location','best');
paths=["direct","surface_reflection"];
for pp=1:2
    nexttile(pp+1); m=s.path_name==paths(pp); semilogx(max(s.offset_m(m),0.05),s.pe_residual_ms(m),'o-','LineWidth',1.4); hold on;
    valid=m & isfinite(s.bellhop_residual_ms); semilogx(s.offset_m(valid),s.bellhop_residual_ms(valid),'s--','LineWidth',1.4);
    yline(0.25,'r:');yline(-0.25,'r:'); grid on; xlabel('offset (m), x=0 plotted at 0.05');ylabel('residual (ms)');title(paths(pp));
end
nexttile(4); m=isfinite(s.bellhop_time_ms); scatter(s.analytic_time_ms(m),s.pe_time_ms(m),45,'filled'); hold on;
scatter(s.analytic_time_ms(m),s.bellhop_time_ms(m),45,'d'); q=xlim; plot(q,q,'k:'); axis equal; grid on;
xlabel('analytic delay (ms)');ylabel('estimated delay (ms)');title('Small-offset absolute delays');legend('PE','Bellhop','1:1');
sgtitle('Delay residuals and x approaching zero');local_export(fig,file);clear cleanup
end

function local_tl(v,file)
t=v.path_table; fig=local_fig([100 100 1200 760]);cleanup=onCleanup(@()close(fig));
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile; cats=categorical(strcat(string(t.offset_m),'m ',t.path_name)); bar(cats,[t.pe_tl_db t.bellhop_tl_db]);
ylabel('raw TL (dB)');grid on;title('Raw path TL');legend('PE','Bellhop');
nexttile; scaled=-20*log10(max(v.amplitude_metrics.global_scale*t.pe_amplitude,realmin));bar(cats,[scaled t.bellhop_tl_db]);
ylabel('TL (dB)');grid on;title('One global direct-path scale');legend('scaled PE','Bellhop');
nexttile([1 2]); b=bar(cats,t.scaled_tl_residual_db); b.FaceColor=[0.25 0.55 0.8];hold on;
yline(1,'r:');yline(-1,'r:');yline(2,'m--');yline(-2,'m--');grid on;
ylabel('scaled PE - Bellhop TL (dB)');title(sprintf('Residuals: spread %.3f dB, reflected RMS %.3f dB', ...
    v.amplitude_metrics.scale_spread_db,v.amplitude_metrics.reflection_rms_db));
sgtitle(sprintf('Amplitude assessment: %s',v.amplitude_status));local_export(fig,file);clear cleanup
end

function local_source(v,file)
t=v.source_audit.table; fig=local_fig([100 100 1150 720]);cleanup=onCleanup(@()close(fig));
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile; semilogx(t.offset_m,t.pe_amplitude/max(t.pe_amplitude),'-o','LineWidth',1.4);hold on;
semilogx(t.offset_m,t.bellhop_amplitude/max(t.bellhop_amplitude),'--s','LineWidth',1.4);
semilogx(t.offset_m,t.bellhop_source_aware_amplitude/max(t.bellhop_source_aware_amplitude),':d','LineWidth',1.4);grid on;
xlabel('offset (m)');ylabel('normalized direct amplitude');title('Direct-only angular trend');legend('PE','Bellhop point source','Bellhop x Gaussian weight');
nexttile; plot(t.offset_m,t.gaussian_weight,'o-','LineWidth',1.4);grid on;xlabel('offset (m)');ylabel('W(kappa)');title('Gaussian angular-spectrum weight');
nexttile([1 2]); plot(t.offset_m,t.raw_residual_db,'o-','LineWidth',1.4);hold on;plot(t.offset_m,t.source_aware_residual_db,'s--','LineWidth',1.4);
yline(1,'r:');yline(-1,'r:');grid on;xlabel('offset (m)');ylabel('global-scale residual (dB)');
title(sprintf('Raw RMS %.3f dB; source-aware RMS %.3f dB (diagnostic only)',v.source_audit.raw_rms_db,v.source_audit.source_aware_rms_db));
legend('raw point source','source-aware diagnostic');local_export(fig,file);clear cleanup
end

function local_window(v,file)
t=v.window_audit.table; fig=local_fig([100 100 1300 760]);cleanup=onCleanup(@()close(fig));
groups=["sampling","aperture","sponge_location","sponge_strength"];
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
for gg=1:4
    nexttile; m=t.group_name==groups(gg); yyaxis left; bar(categorical(t.case_name(m)),t.phase_rms_rad(m));
ylabel('phase RMS (rad)'); yline(0.15,'r:'); yyaxis right; plot(categorical(t.case_name(m)),t.relative_tl_db(m),'ko-','LineWidth',1.2);
ylabel('relative TL (dB)'); yline(0.5,'m:'); yline(-0.5,'m:');grid on;title(strrep(groups(gg),'_',' '));
    if groups(gg)=="aperture"
        labels=compose('edge %.1f dB',t.edge_level_db(m));
        yyaxis left; yl=ylim;
        for jj=1:sum(m), text(jj,0.92*yl(2),labels(jj),'HorizontalAlignment','center', ...
                'FontSize',8,'Rotation',15); end
    end
end
sgtitle(sprintf('Window/sponge decoupling; aperture edge limit %.0f dB',v.config.edge_level_limit_db));
local_export(fig,file);clear cleanup
end

function local_fields(v,file)
o=v.field_output; cfg=o.config; shade=o.field_runs(1,end).shade;
fig=local_fig([100 100 1450 760]);cleanup=onCleanup(@()close(fig));tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile;hold on;rays=o.ray_data.rays;
for ii=1:numel(rays)
    if rays(ii).bottom_bounce_count>0,col=[0.75 0.2 0.2];elseif rays(ii).top_bounce_count>0,col=[0.2 0.55 0.85];else,col=[0.65 0.65 0.65];end
    plot(rays(ii).range_m,rays(ii).depth_m,'Color',col,'HandleVisibility','off');
end
plot(0,cfg.z_tx_m,'kp','MarkerFaceColor',[0.9 0.2 0.1]);plot(cfg.receiver_offsets_m,cfg.z_rx_m*ones(size(cfg.receiver_offsets_m)),'ks','MarkerFaceColor',[1 0.8 0.1]);
set(gca,'YDir','reverse');xlim(cfg.range_limits_m);ylim([0 cfg.water_depth_m]);grid on;xlabel('range (m)');ylabel('depth (m)');title('Bellhop rays');
fields={o.coherent_tl_db,o.incoherent_tl_db};titles={'Coherent TL','Incoherent TL'};
for jj=1:2
    nexttile;h=imagesc(shade.receiver_range_m,shade.receiver_depth_m,fields{jj});set(h,'AlphaData',isfinite(fields{jj}));
    set(gca,'YDir','reverse');clim([20 80]);colormap(flipud(turbo(256)));colorbar;hold on;
    plot(cfg.receiver_offsets_m,cfg.z_rx_m*ones(size(cfg.receiver_offsets_m)),'ws','MarkerFaceColor',[1 0.8 0.1]);
    xlabel('range (m)');ylabel('depth (m)');title(titles{jj});
end
sgtitle('Bellhop ray topology and shared-scale 4 kHz TL fields');local_export(fig,file);clear cleanup
end

function local_beams(v,file)
o=v.field_output; lowc=o.field_runs(1,1).shade; highc=o.field_runs(1,end).shade;
lowi=o.field_runs(2,1).shade; highi=o.field_runs(2,end).shade;
clow=-20*log10(max(abs(lowc.pressure),realmin));chigh=-20*log10(max(abs(highc.pressure),realmin));
ilow=-20*log10(max(abs(lowi.pressure),realmin));ihigh=-20*log10(max(abs(highi.pressure),realmin));
dc=clow-chigh;di=ilow-ihigh;lim=max([1,local_percentile(abs(dc(isfinite(dc))),99), ...
    local_percentile(abs(di(isfinite(di))),99)]);
fig=local_fig([100 100 1350 780]);cleanup=onCleanup(@()close(fig));tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile;imagesc(highc.receiver_range_m,highc.receiver_depth_m,dc);set(gca,'YDir','reverse');clim([-lim lim]);colorbar;title('Coherent: 5001 - 10001 beams (dB)');
nexttile;imagesc(highi.receiver_range_m,highi.receiver_depth_m,di);set(gca,'YDir','reverse');clim([-lim lim]);colorbar;title('Incoherent: 5001 - 10001 beams (dB)');
colormap(redblue(256));
[~,iz]=min(abs(highc.receiver_depth_m-o.config.z_rx_m));nexttile([1 2]);
plot(highc.receiver_range_m,chigh(iz,:),'-','LineWidth',1.3);hold on;plot(highi.receiver_range_m,ihigh(iz,:),'--','LineWidth',1.3);
plot(o.receiver_table.offset_m,o.receiver_table.arrival_synthesis_tl_db,'kd','MarkerFaceColor','k');
plot(o.receiver_table.offset_m,o.receiver_table.pe_scaled_total_tl_db,'ro');grid on;ylim([20 80]);
xlabel('range (m)');ylabel('TL (dB)');title(sprintf('Receiver-depth slice at z=%.2f m',highc.receiver_depth_m(iz)));
legend('Bellhop coherent','Bellhop incoherent','arrival synthesis','globally scaled PE');
sgtitle('Beam-count convergence and receiver comparison');local_export(fig,file);clear cleanup
end

function local_summary(v,file)
t=v.checks; fig=local_fig([100 100 1200 720]);cleanup=onCleanup(@()close(fig));
ratio=t.value./max(t.limit,eps); ratio(t.limit==0)=double(t.value(t.limit==0)>0); ratio=min(ratio,2);
positions=(1:height(t)).';
barh(positions,ratio,'FaceColor',[0.3 0.55 0.8]);hold on;xline(1,'r--','limit','LineWidth',1.5);
set(gca,'YTick',positions,'YTickLabel',t.check_name,'YDir','reverse', ...
    'TickLabelInterpreter','none');
xlim([0 2]);ylim([0.5 height(t)+0.5]);grid on;xlabel('value / limit (clipped at 2)');
title(sprintf('Core %s; amplitude %s',v.decision,v.amplitude_status),'Interpreter','none');
for ii=1:height(t)
    text(min(ratio(ii)+0.03,1.85),positions(ii),string(t.passed(ii)), ...
        'FontWeight','bold');
end
local_export(fig,file);clear cleanup
end

function fig=local_fig(position)
fig=figure('Visible','off','Color','w','Position',position);
end

function local_export(fig,file)
exportgraphics(fig,file,'Resolution',180);
end

function map=redblue(n)
half=floor(n/2); map=[linspace(0,1,half)' linspace(0,1,half)' ones(half,1); ...
    ones(n-half,1) linspace(1,0,n-half)' linspace(1,0,n-half)'];
end

function value=local_percentile(values,p)
values=sort(values(:));
if isempty(values), value=0; return; end
position=1+(numel(values)-1)*p/100;
lo=floor(position); hi=ceil(position);
if lo==hi
    value=values(lo);
else
    value=values(lo)+(position-lo)*(values(hi)-values(lo));
end
end

function local_manifest(file,v,atlas,names)
fid=fopen(file,'w','n','UTF-8');if fid<0,error('Cannot create %s',file);end;cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'# PE/Bellhop comparison atlas manifest\n\n- run_id: `%s`\n',atlas.run_id);
fprintf(fid,'- physical fields: public `H_*_physical_f`; no manual carrier restoration\n');
fprintf(fid,'- TL field limits: shared 20--80 dB; residual maps are zero-centered\n');
fprintf(fid,'- Bellhop 4 kHz arrivals extended over 3--5 kHz only for delay/PDP diagnostics\n\n');
for ii=1:numel(names),fprintf(fid,'%d. `%s`\n',ii,names{ii});end
fprintf(fid,'\nCore decision `%s`; amplitude `%s`.\n',v.decision,v.amplitude_status);clear cleanup
end

function local_atlas_summary(file,v,atlas)
fid=fopen(file,'w','n','UTF-8');if fid<0,error('Cannot create %s',file);end;cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'PE/Bellhop comparison atlas\nrun_id %s\ncore %s\namplitude %s\nfigures %d\n', ...
    atlas.run_id,v.decision,v.amplitude_status,numel(atlas.figure_files));clear cleanup
end
