function output = generate_fk_optimization_ppt_visuals(output_dir)
%GENERATE_FK_OPTIMIZATION_PPT_VISUALS Render saved FK evidence for slides.
% No PE, FK, BIE, RNG, rank search, or timing benchmark is executed.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if nargin < 1
    output_dir = fullfile(root,'results','visualization','fk_optimization_ppt');
end
figdir = fullfile(output_dir,'figures');
if ~isfolder(figdir), mkdir(figdir); end
v = fullfile(root,'results','validation');
pdir = fullfile(v,'spectral_fast_kirchhoff_prototype');
mdir = fullfile(v,'spectral_fast_kirchhoff_pm_pilot');
bdir = fullfile(v,'full_kirchhoff_pm_validity_boundary');
fdir = fullfile(v,'spectral_fast_kirchhoff_multifrequency_consistency');
P = readtable(fullfile(pdir,'all_metrics.csv'));
S = readtable(fullfile(pdir,'scaling.csv'));
M = readtable(fullfile(mdir,'field_metrics.csv'));
B = readtable(fullfile(bdir,'validity_metrics.csv'));
F = readtable(fullfile(fdir,'dense_metrics.csv'));
A = readtable(fullfile(fdir,'bie_anchors.csv'));
G = readtable(fullfile(fdir,'phase_group_delay.csv'));
C = readtable(fullfile(fdir,'reflected_cir.csv'));
q = load(fullfile(mdir,'seed_260001_fields.mat'),'d'); pm = q.d;
q = load(fullfile(pdir,'weak_fields.mat'),'d'); weak = q.d;
q = load(fullfile(pdir,'strong_height_fields.mat'),'d'); strong = q.d;
assert(height(F)==65 && all(F.rank_in==3 & F.rank_out==12));
assert(all(isfinite(F.weighted)) && max(F.weighted)<1e-3);
blue = [0 .36 .62]; orange = [.82 .32 .05]; purple = [.43 .28 .64];
names = {'ppt_overview_roadmap','ppt_geometry_setup','ppt_incident_field', ...
    'ppt_reflected_field_comparison','ppt_incident_vs_reflected', ...
    'ppt_prototype_accuracy_speed','ppt_pm_pilot_algorithm_vs_physics', ...
    'ppt_pm_validity_boundary','ppt_multifrequency_consistency', ...
    'ppt_rank_evolution','ppt_summary_dashboard'};

% 1: concise opening story, with separate evidence boundaries.
f = canvas('FK optimization: progress and evidence', ...
    'Saved 2-D validation evidence | Production model unchanged');
ax = axes(f,'Position',[.04 .43 .92 .39]); axis(ax,[0 5 0 1]); axis(ax,'off');
stages = {'PE up','FullFK','FastFK','PM validity','Broadband'};
details = {'Frozen incident', 'BIE reference', 'FFT compression', 'Conditional domain', '4-8 kHz / 65 points'};
for j=1:5
    rectangle(ax,'Position',[j-1+.05 .23 .84 .62],'Curvature',.12,'FaceColor',[.94 .97 .99],'EdgeColor',blue);
    text(ax,j-.53,.64,stages{j},'FontSize',24,'FontWeight','bold','HorizontalAlignment','center');
    text(ax,j-.53,.40,details{j},'FontSize',17,'HorizontalAlignment','center');
    if j<5, text(ax,j-.025,.54,'>','FontSize',28,'HorizontalAlignment','center'); end
end
card(f,[.06 .16 .26 .21],sprintf('141x warm\nPrototype, N=4097'),'3/6 ranks; cold 14.7x');
card(f,[.37 .16 .26 .21],'89-112x warm','PM pilot: 3/8 algorithm PASS');
card(f,[.68 .16 .26 .21],sprintf('2.31e-4 max error\nBroadband rank 3/12'),'Fixed canonical surface only');
foot(f,'BIE agreement in tested M99 core; high-cutoff PM fails physics. Next eligible scope: small realization data, not executed.');
emit(f,figdir,names{1});

% 2: physical geometry and honest zoom; no invented ray-trace result.
f = canvas('What reflects at the sea surface?', ...
    'Research fixture: Tx z=100 m, Rx z=3 m | Positive depth points downward');
t = tiledlayout(f,1,2,'Padding','loose','TileSpacing','loose');
ax = nexttile(t); hold(ax,'on');
plot(ax,pm.x_m,pm.eta_m,'Color',blue,'LineWidth',2.3);
plot(ax,[-60 60],[3 3],'--','Color',purple,'LineWidth',2);
plot(ax,0,100,'o','MarkerFaceColor',orange,'MarkerEdgeColor',orange,'MarkerSize',12);
quiver(ax,0,96,-20,-92,0,'Color',orange,'LineWidth',2,'MaxHeadSize',.18);
quiver(ax,-20,1,33,2,0,'Color',purple,'LineWidth',2,'MaxHeadSize',.3);
text(ax,4,93,'Tx','FontSize',22); text(ax,-48,53,'Incident','FontSize',21,'Color',orange);
text(ax,17,12,'Reflected to Rx','FontSize',19,'Color',purple);
set(ax,'YDir','reverse'); xlim(ax,[-60 60]); ylim(ax,[-6 105]);
xlabel(ax,'Lateral position x (m)'); ylabel(ax,'Depth z (m)'); title(ax,'Propagation sketch'); style(ax);
ax = nexttile(t); hold(ax,'on');
plot(ax,pm.x_m,pm.eta_m,'Color',blue,'LineWidth',2.5);
plot(ax,[-60 60],[3 3],'--','Color',purple,'LineWidth',2);
xline(ax,-39.84375,':','M99','LabelVerticalAlignment','middle'); xline(ax,39.84375,':');
set(ax,'YDir','reverse'); ylim(ax,[-.8 3.6]); xlim(ax,[-60 60]);
xlabel(ax,'Lateral position x (m)'); ylabel(ax,'Depth z (m)');
title(ax,'Canonical PM: surface / Rx zoom'); legend(ax,{'PM surface','Rx plane'},'Location','southwest'); style(ax);
foot(f,'Arrows are schematic, not computed rays. C2 surface: plateau |x|<=42 m, support 50 m. Tx100 is not the Tx99 Bellhop installation.');
emit(f,figdir,names{2});

% 3: exact source from saved Fourier coefficients, and saved compressed input.
u = exact_incident(weak); up = pm.refcfg.surface_incident;
f = canvas('Incident: the field arriving at the surface', ...
    '4 kHz | Exact saved-source reconstruction and saved FastFK incident | No phase fitting');
t = tiledlayout(f,2,2,'Padding','loose','TileSpacing','loose');
ds = {weak,pm}; us = {u,up}; labels = {'Weak deterministic surface','Canonical PM surface'};
incident_errors = zeros(2,1);
for j=1:2
    d=ds{j}; ui=us{j}; ui_fast=d.selected.surface_incident;
    incident_errors(j)=norm(ui_fast-ui)/norm(ui);
    ax=nexttile(t,j); hold(ax,'on');
    scale=max(abs(ui)); plot(ax,d.x_m,abs(ui)/scale,'Color',blue,'LineWidth',2.8);
    plot(ax,d.x_m,abs(ui_fast)/scale,'--','Color',orange,'LineWidth',1.8);
    title(ax,labels{j}); ylabel(ax,'Magnitude / exact peak'); xlim(ax,[-40 40]); ylim(ax,[0 1.08]); style(ax);
    legend(ax,{'Incident exact','FastFK input (rank 3)'},'Location','south');
    ax=nexttile(t,j+2); hold(ax,'on');
    plot(ax,d.x_m,unwrap(angle(ui)),'Color',blue,'LineWidth',2.5);
    plot(ax,d.x_m,unwrap(angle(ui_fast)),'--','Color',orange,'LineWidth',1.8);
    xlabel(ax,'Surface lateral position x (m)'); ylabel(ax,'Unwrapped phase (rad)'); xlim(ax,[-40 40]); style(ax);
end
foot(f,'Exact angular-spectrum input is shared by FullFK/BIE. This panel compares its FastFK compression; it is not a new independent PE test.');
emit(f,figdir,names{3});

% 4: native receiving fields, and raw residual phase at a readable scale.
f = canvas('Reflected field: FullFK and FastFK overlap', ...
    '4 kHz | FullFK solid blue, FastFK dashed orange, BIE dotted purple');
t=tiledlayout(f,2,3,'Padding','loose','TileSpacing','loose');
cases={weak,strong,pm}; labels={'Weak / rank 3/6','Strong height / rank 3/6','Canonical PM / rank 3/8'};
for j=1:3
    d=cases{j}; scale=max(abs(d.reference_field(d.mask)));
    ax=nexttile(t,j); hold(ax,'on');
    plot(ax,d.x_m,abs(d.reference_field)/scale,'Color',blue,'LineWidth',3);
    plot(ax,d.x_m,abs(d.selected.receiver_field)/scale,'--','Color',orange,'LineWidth',1.9);
    if j==3
        R=readtable(fullfile(bdir,'case_data','cutoff_n12_receiver_fields.csv'));
        rb=complex(R.BIE_real,R.BIE_imag);
        plot(ax,R.receiver_x_m,abs(rb)/scale,':','Color',purple,'LineWidth',2.2);
    end
    xlim(ax,[-39.84375 39.84375]); ylim(ax,[0 1.12]); title(ax,labels{j}); ylabel(ax,'Normalized magnitude'); style(ax);
    ax=nexttile(t,j+3); hold(ax,'on');
    delta=angle(d.selected.receiver_field.*conj(d.reference_field)); delta(~d.mask)=NaN;
    plot(ax,d.x_m,delta*1e3,'--','Color',orange,'LineWidth',2);
    if j==3
        delta_b=angle(complex(R.Full_real,R.Full_imag).*conj(rb)); delta_b(~logical(R.fixed_mask))=NaN;
        plot(ax,R.receiver_x_m,delta_b*1e3,':','Color',purple,'LineWidth',2);
        legend(ax,{'FastFK-FullFK','FullFK-BIE'},'Location','south');
    end
    xlim(ax,[-39.84375 39.84375]); xlabel(ax,'Receiver lateral position (m)'); ylabel(ax,'Phase difference (mrad)'); style(ax);
end
foot(f,'Each column uses one FullFK normalization; residual phase is unadjusted. M99 / Full-only -40 dB mask; edges are not physically certified.');
emit(f,figdir,names{4});

% 5: common incident reference preserves the amplitude change on reflection.
f=canvas('Canonical PM: before and after reflection', ...
    '4 kHz | One shared amplitude reference: peak exact Incident | Rx depth 3 m');
t=tiledlayout(f,2,2,'Padding','loose','TileSpacing','loose');
scale=max(abs(up));
ax=nexttile(t,[1 2]); plot(ax,pm.x_m,pm.eta_m,'Color',blue,'LineWidth',2.6);
xlim(ax,[-40 40]); ylabel(ax,'Surface depth (m)'); title(ax,'The same fixed random surface'); style(ax);
ax=nexttile(t); plot(ax,pm.x_m,abs(up)/scale,'Color',blue,'LineWidth',2.8);
xlim(ax,[-40 40]); ylim(ax,[0 1.2]); title(ax,'Incident on the true surface'); xlabel(ax,'Surface position x (m)'); ylabel(ax,'Magnitude / Incident peak'); style(ax);
ax=nexttile(t); hold(ax,'on'); plot(ax,pm.x_m,abs(pm.reference_field)/scale,'Color',blue,'LineWidth',2.8);
plot(ax,pm.x_m,abs(pm.selected.receiver_field)/scale,'--','Color',orange,'LineWidth',1.8);
xlim(ax,[-40 40]); ylim(ax,[0 1.2]); title(ax,'Reflected on the Rx plane'); xlabel(ax,'Receiver position x (m)'); ylabel(ax,'Magnitude / Incident peak'); style(ax);
legend(ax,{'FullFK','FastFK'},'Location','south');
foot(f,'Saved 1-D field profiles show redistribution after scattering; no unsaved 2-D field map is inferred.'); emit(f,figdir,names{5});

% 6: logarithmic points avoid misleading truncated bars.
f=canvas('Prototype: accurate and faster', ...
    'FastFK-FullFK algorithm fidelity | Prototype-B (full geometry), padding 2');
t=tiledlayout(f,1,3,'Padding','loose','TileSpacing','loose');
Q=P(strcmp(P.study,'rank') & strcmp(P.geometry,'B') & P.rank_in==3,:);
ax=nexttile(t); semilogy(ax,1:height(Q),Q.complex_l2,'o','Color',blue,'MarkerFaceColor',blue,'MarkerSize',9);
xticks(ax,1:height(Q)); xticklabels(ax,strrep(Q.case_name,'_',' ')); xtickangle(ax,35);
yline(ax,1e-4,'--','1e-4 target'); ylabel(ax,'Weighted relative complex L2'); title(ax,'Six cases / rank 3/6'); style(ax);
ax=nexttile(t); hold(ax,'on');
for c={'weak','strong_height','high_K'}
    Q=P(strcmp(P.study,'outgoing_rank') & strcmp(P.case_name,c{1}) & strcmp(P.geometry,'B'),:);
    plot(ax,Q.rank_out,Q.complex_l2,'-o','LineWidth',2,'DisplayName',strrep(c{1},'_',' '));
end
set(ax,'YScale','log'); xlabel(ax,'Outgoing rank (incoming fixed at 4)'); ylabel(ax,'Weighted relative complex L2'); title(ax,'Rank convergence'); legend(ax,'Location','southwest'); style(ax);
ax=nexttile(t); loglog(ax,S.n_surface,S.full_apply_wall_s,'-o','Color',blue,'LineWidth',2.5); hold(ax,'on');
loglog(ax,S.n_surface,S.fast_apply_wall_s,'--s','Color',orange,'LineWidth',2.5);
xlabel(ax,'N surface = N receiver'); ylabel(ax,'Warm apply wall time (s)'); title(ax,'N=4097: 141x warm'); legend(ax,{'FullFK','FastFK'},'Location','northwest'); style(ax);
xticks(ax,S.n_surface); xticklabels(ax,string(S.n_surface)); ax.XAxis.FontSize=16;
foot(f,'Measured kernel times: FullFK median 2 / FastFK median 5. N4097 cold speedup 14.7x; speed is session-dependent.'); emit(f,figdir,names{6});

% 7: separate model error from compression error, on two axes.
f=canvas('PM pilot: algorithm passes, high-cutoff physics fails', ...
    'Compression error and physical approximation error answer different questions');
t=tiledlayout(f,1,2,'Padding','loose','TileSpacing','loose');
ax=nexttile(t); hold(ax,'on');
for st={'baseline','selected'}
    Q=M(strcmp(M.study,st{1}),:);
    plot(ax,1:height(Q),Q.weighted_complex_l2,'-o','LineWidth',2.5,'MarkerSize',8);
end
set(ax,'YScale','log'); yline(ax,1e-3,'k--','Gate 1e-3'); xticks(ax,1:5); xticklabels(ax,{'Seed 1','Seed 2','Seed 3','Low K','High K'});
ylabel(ax,'Weighted relative complex L2'); title(ax,'FastFK-FullFK / fidelity'); legend(ax,{'3/6 FAIL','3/8 PASS'},'Location','south'); style(ax);
ax=nexttile(t); e=[B.complex_l2(strcmp(B.case_name,'cutoff_n12')); B.complex_l2(strcmp(B.case_name,'cutoff_n19'))];
bar(ax,1:2,e,'FaceColor',purple); yline(ax,1e-3,'k--','Gate 1e-3'); xticks(ax,1:2); xticklabels(ax,{'Canonical','High cutoff'});
ylim(ax,[0 1.5e-3]); ylabel(ax,'Weighted relative complex L2'); title(ax,'FullFK-BIE / physics'); style(ax);
foot(f,'Selected 3/8 warm speedup: 89-112x. Raising FastFK rank cannot repair the FullFK physical-model discrepancy.'); emit(f,figdir,names{7});

% 8: guarded endpoints and all descriptor points, without a fitted law.
f=canvas('PM validity is conditional on illuminated structure', ...
    'One realization / 4 kHz / Tx100-Rx3 / frozen M99 core | No universal cutoff law');
t=tiledlayout(f,1,2,'Padding','loose','TileSpacing','loose');
ax=nexttile(t); Q=B(strcmp(B.family,'cutoff'),:); Q=sortrows(Q,'Kmax_radpm');
plot(ax,Q.Kmax_radpm,Q.complex_l2,'-o','Color',blue,'LineWidth',2.3); hold(ax,'on');
Z=Q(isfinite(Q.numerical_uncertainty),:);
errorbar(ax,Z.Kmax_radpm,Z.complex_l2,Z.numerical_uncertainty,'s','Color',orange,'MarkerSize',10,'LineWidth',2);
yline(ax,1e-3,'k--','Gate 1e-3'); xlabel(ax,'Physical Kmax (rad/m)'); ylabel(ax,'FullFK-BIE weighted complex L2'); title(ax,'Guarded pass / fail bracket'); style(ax);
text(ax,.705,.00072,'0.706858 PASS','HorizontalAlignment','right','FontSize',18);
text(ax,.746,.00131,'0.746128 FAIL','HorizontalAlignment','right','FontSize',18); ylim(ax,[0 .0015]);
ax=nexttile(t); hold(ax,'on');
plot(ax,B.illuminated_curvature_over_k,B.complex_l2,'o','Color',[.55 .55 .55],'MarkerSize',8,'DisplayName','All 15 controls');
Q=B(strcmp(B.family,'amplitude') | strcmp(B.case_name,'cutoff_n19'),:); Q=sortrows(Q,'added_amplitude');
plot(ax,Q.illuminated_curvature_over_k,Q.complex_l2,'s','Color',orange,'MarkerFaceColor',orange,'MarkerSize',9,'DisplayName','Same Kmax / varied amplitude');
yline(ax,1e-3,'k--','HandleVisibility','off'); xlabel(ax,'Illuminated RMS curvature / k'); ylabel(ax,'FullFK-BIE weighted complex L2'); title(ax,'Illuminated structure'); legend(ax,'Location','northwest'); style(ax);
foot(f,'Endpoint bars = summed numerical changes (~6.8e-6), not confidence intervals. Most interior controls are baseline-only.'); emit(f,figdir,names{8});

% 9: primary DSP bridge tables; phase and delay use separate axes.
f=canvas('4-8 kHz: one rank, consistent complex response', ...
    'Canonical seed260001 only | 65 frequencies / uniform rank 3/12 / padding 2');
t=tiledlayout(f,2,2,'Padding','loose','TileSpacing','loose');
ax=nexttile(t); semilogy(ax,F.frequency_hz/1000,F.weighted,'-','Color',blue,'LineWidth',2.6); yline(ax,1e-3,'k--','Gate');
xlabel(ax,'Frequency (kHz)'); ylabel(ax,'FastFK-FullFK weighted L2'); title(ax,'Max error 2.31e-4'); style(ax);
ax=nexttile(t); plot(ax,A.frequency_hz/1000,A.weighted,'o','Color',purple,'MarkerFaceColor',purple,'MarkerSize',10); yline(ax,1e-3,'k--','Gate');
ylim(ax,[0 .0012]); xlim(ax,[3.8 8.2]); xlabel(ax,'Frequency (kHz)'); ylabel(ax,'FullFK-BIE weighted L2'); title(ax,'Three physical anchors'); style(ax);
ax=nexttile(t); plot(ax,G{:,1}/1000,G{:,10}*1e3,'Color',orange,'LineWidth',2.5);
xlabel(ax,'Frequency (kHz)'); ylabel(ax,'FastFK-FullFK phase (mrad)'); title(ax,'Max phase 0.174 mrad / delay 70.2 ns'); style(ax);
ax=nexttile(t); hc=complex(C{:,2},C{:,3}); hf=complex(C{:,4},C{:,5}); scale=max(abs(hc));
plot(ax,C{:,1}*1000,abs(hc)/scale,'-o','Color',blue,'LineWidth',2,'MarkerSize',4); hold(ax,'on');
plot(ax,C{:,1}*1000,abs(hf)/scale,'--','Color',orange,'LineWidth',2);
xlim(ax,[0 16]); xlabel(ax,'Relative delay (ms)'); ylabel(ax,'Magnitude / FullFK peak'); title(ax,'Reflected-only CIR: diagnostic'); legend(ax,{'FullFK','FastFK'},'Location','northeast'); style(ax);
foot(f,'Native-to-DSP bridge applied once. CIR uses rectangular 65-point IFFT; 16 ms periodic window does not certify unaliased long tails.'); emit(f,figdir,names{9});

% 10: rank configuration progression, not continuous optimization evidence.
f=canvas('Uniform ranks grow with the tested scope', ...
    'Rank notation: incoming / outgoing | One shared configuration within each stage');
ax=axes(f,'Position',[.08 .32 .84 .46]); hold(ax,'on');
plot(ax,1:3,[6 8 12],'-s','Color',orange,'LineWidth',3,'MarkerSize',14,'MarkerFaceColor',orange);
plot(ax,1:3,[3 3 3],'--o','Color',blue,'LineWidth',2.5,'MarkerSize',12,'MarkerFaceColor',blue);
xticks(ax,1:3); xticklabels(ax,{'Prototype 3/6','PM pilot 3/8','4-8 kHz 3/12'}); xlim(ax,[.7 3.3]); ylim(ax,[0 14]); yticks(ax,[0 3 6 8 12]); ylabel(ax,'Rank');
legend(ax,{'Outgoing rank','Incoming rank'},'Location','northwest'); style(ax);
card(f,[.12 .13 .34 .12],'3/6 fails on PM','3/8 chosen across five PM cases');
card(f,[.55 .13 .34 .12],'3/8 fails from 6 kHz','3/10 narrowly fails; 3/12 passes band');
foot(f,'Finite rank ladders retain failed trials. Configuration is shared across cases / frequencies; no per-frequency rank tuning.'); emit(f,figdir,names{10});

% 11: large, readable closing dashboard.
f=canvas('FK optimization: established results and limits', ...
    'The successful compression approximates FullFK; its physical domain remains independently checked');
rows={ ...
    'Prototype','Compress deterministic FullFK','Max E 3.37e-5','141x warm / 14.7x cold','3/6 fidelity PASS'; ...
    'PM pilot','Test five fixed PM inputs','Max E 2.48e-4','89-112x warm','3/8 fidelity PASS; physics mixed'; ...
    'PM boundary','Separate physical validity','8.53e-4 pass / 1.20e-3 fail','No Fast benchmark','Conditional single-realization bracket'; ...
    '4-8 kHz','Check one shared rank','Max E 2.31e-4','Min 41.9x warm / 9.83x cold','3/12 PASS; CIR diagnostic'};
heads={'Stage','Goal','Representative error','Measured benefit','Conclusion'};
ax=axes(f,'Position',[.035 .24 .93 .55]); axis(ax,[0 1 0 1]); axis(ax,'off');
widths=[.13 .24 .22 .20 .21]; edges=[0 cumsum(widths)];
for r=0:4
    y=.84-r*.19;
    rectangle(ax,'Position',[0 y .999 .17],'FaceColor', [.95 .97 .99] - .025*mod(r,2),'EdgeColor','none');
    for c=1:5
        if r==0, str=heads{c}; weight='bold'; else, str=rows{r,c}; weight='normal'; end
        text(ax,edges(c)+.009,y+.085,wrap_text(str,round(widths(c)*105)), ...
            'FontSize',18,'FontWeight',weight,'VerticalAlignment','middle','Interpreter','none');
    end
end
card(f,[.08 .10 .84 .10],'Next eligible scope: a small realization dataset','Under tested physics and explicit phase metadata; no new generation in this task');
foot(f,'All E values are relative weighted complex L2. Timings are hardware/session-specific; no Monte Carlo, 3-D or production claim.'); emit(f,figdir,names{11});

% Retain exact plotted incident source data and a compact machine-readable audit.
writetable(table(weak.x_m,real(u),imag(u),real(weak.selected.surface_incident),imag(weak.selected.surface_incident), ...
    'VariableNames',{'x_m','exact_real','exact_imag','Fast_input_real','Fast_input_imag'}),fullfile(output_dir,'weak_incident_source.csv'));
writetable(table(pm.x_m,real(up),imag(up),real(pm.selected.surface_incident),imag(pm.selected.surface_incident), ...
    'VariableNames',{'x_m','exact_real','exact_imag','Fast_input_real','Fast_input_imag'}),fullfile(output_dir,'pm_incident_source.csv'));
output=struct('figure_names',{names},'figure_dir',figdir,'matlab_version',version, ...
    'incident_compression_whole_l2',incident_errors,'multifrequency_max_weighted',max(F.weighted), ...
    'max_phase_rad',max(abs(G{:,10})),'max_delay_difference_s',max(abs(G{:,8}-G{:,9})), ...
    'cir_complex_l2',norm(hf-hc)/norm(hc),'new_numerical_solves',0);
fid=fopen(fullfile(output_dir,'visualization_audit.json'),'w'); clean=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(output,PrettyPrint=true));
disp(output);
end

function f=canvas(title_text,subtitle)
f=figure('Visible','off','Color','white','Position',[30 30 1600 900]);
annotation(f,'textbox',[.04 .925 .92 .06],'String',title_text,'FontSize',28, ...
    'FontName','Arial','FontWeight','bold','EdgeColor','none','Interpreter','none');
annotation(f,'textbox',[.04 .875 .92 .045],'String',subtitle,'FontSize',17, ...
    'FontName','Arial','EdgeColor','none','Interpreter','none','Color',[.25 .25 .25]);
% Leave explicit title/footer margins for all tiled layouts.
set(f,'DefaultAxesFontName','Arial','DefaultTextFontName','Arial');
end

function style(ax)
set(ax,'FontSize',18,'LineWidth',1,'Box','off','TickDir','out','XGrid','on','YGrid','on', ...
    'XMinorGrid','off','YMinorGrid','off','GridAlpha',.12);
if isa(ax.Parent,'matlab.graphics.layout.TiledChartLayout')
    ax.Parent.OuterPosition=[.045 .115 .91 .735];
end
end

function foot(f,str)
annotation(f,'textbox',[.04 .015 .92 .065],'String',str,'FontSize',15,'FontName','Arial', ...
    'EdgeColor','none','Color',[.25 .25 .25],'Interpreter','none');
end

function card(f,pos,headline,detail)
annotation(f,'textbox',pos,'String',{headline,detail},'FontSize',21,'FontName','Arial', ...
    'EdgeColor',[.78 .84 .89],'BackgroundColor',[.96 .98 1],'Margin',12,'Interpreter','none');
end

function emit(f,folder,name)
exportgraphics(f,fullfile(folder,[name '.png']),'Resolution',150,'BackgroundColor','white', ...
    'Units','inches','Width',16,'Height',9,'Padding','figure');
exportgraphics(f,fullfile(folder,[name '.pdf']),'ContentType','vector','BackgroundColor','white', ...
    'Units','inches','Width',16,'Height',9,'Padding','figure');
close(f); fprintf('Rendered %s\n',name);
end

function u=exact_incident(d)
% Reevaluate only the saved discrete source at saved surface points.
% This presentation reconstruction duplicates no scattering/propagation solve.
k=2*pi*4000/1500; kz=sqrt(complex(k^2-d.source_kx.^2,0));
x0=-96.09375; n=numel(d.x_m); u=complex(zeros(n,1));
for first=1:128:n
    rows=first:min(first+127,n);
    phase=exp(1i*(d.x_m(rows)-x0)*d.source_kx.').*exp(-1i*(d.eta_m(rows)-100).*kz.');
    u(rows)=phase*d.source_coefficients/numel(d.source_coefficients);
end
end

function out=wrap_text(str,maxlen)
words=strsplit(str); out=''; line='';
for j=1:numel(words)
    if numel(line)+numel(words{j})+1>maxlen && ~isempty(line)
        out=[out line newline]; %#ok<AGROW>
        line=words{j};
    elseif isempty(line), line=words{j};
    else, line=[line ' ' words{j}]; %#ok<AGROW>
    end
end
out=[out line];
end
