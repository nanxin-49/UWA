function files = render_two_dimensional_closeout_s0_s1(plot_data_path, figure_dir)
%RENDER_TWO_DIMENSIONAL_CLOSEOUT_S0_S1 Draw only from saved plot_data.mat.
assert(isfile(plot_data_path));
assert(~isfolder(figure_dir) || isempty(dir(fullfile(figure_dir,'F*.png'))), ...
    'Refusing to overwrite previously rendered closeout figures.');
if ~isfolder(figure_dir), mkdir(figure_dir); end
t=load(plot_data_path,'plot_data'); d=t.plot_data;
assert(numel(d.line.cases)==3 && isequal(size(d.raw.PE_carrier_restored),[411 45]));
ids=string(d.grid.receiver_id); xx=d.grid.requested_x_m; zz=d.grid.requested_z_m;
colors=[.80 .31 .09;.38 .25 .63;0 .37 .62];

% F01: three different spatial scales stay explicit.
f=canvas('F01 | Fixed geometry and actual saved coordinates');
lay=tiledlayout(f,1,3,'Padding','compact','TileSpacing','compact');
ax=nexttile(lay); hold(ax,'on');
plot(ax,0,100,'^','MarkerFaceColor',colors(1,:),'MarkerSize',10);
plot(ax,xx,zz,'ko','MarkerFaceColor',colors(2,:),'MarkerSize',4);
plot(ax,[-50 50],[0 0],'k-','LineWidth',1.5);
xlim(ax,[-50 50]); ylim(ax,[-2 105]); ylabel(ax,'Depth z (m), downward'); xlabel(ax,'x (m)');
title(ax,'Tx=(0,100), receivers near z=3'); set(ax,'YDir','reverse'); grid(ax,'on');
ax=nexttile(lay); hold(ax,'on'); c=d.line.cases(3);
plot(ax,d.line.x_m,-c.eta_m,'Color',colors(1,:),'LineWidth',1.3);
plot(ax,d.line.x_m,zeros(size(d.line.x_m)),'k:');
xline(ax,-42,'--'); xline(ax,42,'--'); xline(ax,-50,':'); xline(ax,50,':');
xlim(ax,[-52 52]); xlabel(ax,'x (m)'); ylabel(ax,'Physical surface z=-eta (m)');
title(ax,'Strong sine; C2 plateau 42 m, support 50 m'); set(ax,'YDir','reverse'); grid(ax,'on');
ax=nexttile(lay); hold(ax,'on');
for k=1:9
    plot(ax,xx(k),zz(k),'o','MarkerFaceColor',colors(2,:),'Color','k','MarkerSize',7);
    text(ax,xx(k)+.006,zz(k),ids(k),'FontSize',10);
end
xlim(ax,[-.30 .40]); ylim(ax,[2.78 3.22]); xlabel(ax,'x (m)'); ylabel(ax,'z (m), downward');
title(ax,'3 x 3 diagnostic points; E02 x-z window 0.5:0.125:6 m');
set(ax,'YDir','reverse'); grid(ax,'on');
subtitle(lay,'F01: source, physical surface, numerical window and nine exact saved receiver coordinates.');
emit(f,figure_dir,'F01');

% F03: all displayed pixels have saved complex values; Bellhop native scale.
f=canvas('F03 | Strong-sine raw reflected fields, 411 x 45 saved samples');
lay=tiledlayout(f,2,3,'Padding','compact','TileSpacing','compact');
names={'PE_carrier_restored','BIE_native','BH_native'};
titles={'PE Model-0, carrier restored','BIE native','Bellhop X native'};
ceiling=max([abs(d.raw.PE_carrier_restored(:));abs(d.raw.BIE_native(:))]);
bh_ceiling=max(abs(d.raw.BH_native(:)));
for mi=1:3
    u=d.raw.(names{mi}); ax=nexttile(lay,mi);
    imagesc(ax,d.raw.x_m,d.raw.z_m,abs(u).'); set(ax,'YDir','reverse');
    colormap(ax,parula(256));
    if mi==3, clim(ax,[0 bh_ceiling]); else, clim(ax,[0 ceiling]); end
    colorbar(ax); hold(ax,'on'); plot(ax,xx,zz,'ws','MarkerFaceColor','none','MarkerSize',5);
    title(ax,[titles{mi} ' | magnitude']); xlabel(ax,'x (m)'); ylabel(ax,'z (m)');
    ax=nexttile(lay,mi+3); imagesc(ax,d.raw.x_m,d.raw.z_m,angle(u).');
    set(ax,'YDir','reverse'); colormap(ax,hsv(256)); clim(ax,[-pi pi]); colorbar(ax);
    hold(ax,'on'); plot(ax,xx,zz,'ws','MarkerFaceColor','none','MarkerSize',5);
    title(ax,[titles{mi} ' | native phase (rad)']); xlabel(ax,'x (m)'); ylabel(ax,'z (m)');
end
subtitle(lay,sprintf('F03: PE/BIE magnitude common ceiling %.4g; Bellhop native %.4g. Native phase is not cross-method residual.',ceiling,bh_ceiling));
emit(f,figure_dir,'F03');

% F04: saved G on original line, with same-depth point markers only.
f=canvas('F04 | Original M99 receiver line: G magnitude and phase');
lay=tiledlayout(f,3,4,'Padding','compact','TileSpacing','compact');
for ci=1:3
    c=d.line.cases(ci); mask=c.plot_mask; x=d.line.x_m(mask);
    fields={c.G_PE(mask),c.G_BIE(mask),c.G_BH(mask)};
    for col=1:4
        ax=nexttile(lay,(ci-1)*4+col); hold(ax,'on');
        for mi=1:3
            if col==1 || col==3, y=abs(fields{mi});
            else, y=angle(fields{mi}); gaps=[false;abs(diff(y))>pi]; y(gaps)=NaN; end
            plot(ax,x,y,'Color',colors(mi,:),'LineWidth',1.1);
        end
        if col>=3, xlim(ax,[-.5 .5]); else, xlim(ax,[min(x) max(x)]); end
        if col==2 || col==4, ylim(ax,[-pi pi]); end
        for px=[-.1953125 0 .1953125], xline(ax,px,':','Color',[.45 .45 .45]); end
        ylabel(ax,ternary(col==1||col==3,'|G|','arg G (rad)')); xlabel(ax,'x (m), z=3 m');
        labels={'full magnitude','full phase','local magnitude','local phase'};
        title(ax,sprintf('%s | %s',c.name,labels{col}),'Interpreter','none');
        grid(ax,'on');
    end
end
subtitle(lay,'F04: PE orange, BIE purple, Bellhop blue. P21/P22/P23 only; P12/P32 are at other depths. Phase wraps are gaps.');
emit(f,figure_dir,'F04');

% F05: line residuals and explicit missing 3x3 G support.
f=canvas('F05 | PE-BIE G residuals and the incomplete 3 x 3 support');
lay=tiledlayout(f,2,3,'Padding','compact','TileSpacing','compact');
for ci=1:3
    c=d.line.cases(ci); ix=c.plot_mask; A=c.G_PE(ix); B=c.G_BIE(ix); x=d.line.x_m(ix);
    ax=nexttile(lay,ci); hold(ax,'on');
    plot(ax,x,abs(A-B),'Color',colors(1,:),'LineWidth',1.3);
    title(ax,[c.name ' | complex |delta G|'],'Interpreter','none');
    xlabel(ax,'x (m), z=3 m'); grid(ax,'on');
    ax=nexttile(lay,ci+3); hold(ax,'on');
    yyaxis(ax,'left'); plot(ax,x,20*log10(abs(A)./abs(B)),'Color',colors(2,:));
    ylabel(ax,'Magnitude ratio (dB)'); yyaxis(ax,'right');
    plot(ax,x,angle(A.*conj(B)),'Color',colors(3,:)); ylabel(ax,'Wrapped phase (rad)');
    xlabel(ax,'x (m), z=3 m'); grid(ax,'on');
end
subtitle(lay,'F05: top complex residual; bottom signed magnitude dB and wrapped phase on separately labeled axes.');
emit(f,figure_dir,'F05');
f=canvas('F05 continued | Strong-sine G: missing points stay blank');
lay=tiledlayout(f,1,3,'Padding','compact','TileSpacing','compact');
q=d.points(d.points.case_id=="strong_height_low_K" & d.points.pair=="PE_BIE" & d.points.quantity_kind=="G",:);
heat={reshape(q.relative_abs,3,3).',reshape(q.magnitude_db,3,3).',reshape(q.phase_rad,3,3).'};
labels={'Point relative complex error','Signed magnitude difference (dB)','Wrapped phase difference (rad)'};
for k=1:3
    ax=nexttile(lay); h=imagesc(ax,1:3,1:3,heat{k}); h.AlphaData=isfinite(heat{k});
    set(ax,'Color',[.88 .88 .88],'YDir','reverse','XTick',1:3,'XTickLabel',{'left','center','right'}, ...
        'YTick',1:3,'YTickLabel',{'top','middle','bottom'});
    colorbar(ax); title(ax,labels{k});
    for iz=1:3
        for ix=1:3
            id=ids((iz-1)*3+ix);
            if ~isfinite(heat{k}(iz,ix)), text(ax,ix,iz,id+" missing",'HorizontalAlignment','center','FontSize',9);
            else, text(ax,ix,iz,id,'HorizontalAlignment','center','FontSize',10); end
        end
    end
end
subtitle(lay,'F05 grid: saved G exists at z=3 m only. Gray cells are missing; no interpolation.');
emit(f,figure_dir,'F05_grid');

% F06: five fixed points. Mapped complex comparison at z=3; raw amplitudes all five.
f=canvas('F06 | Five fixed points: comparable G versus raw native amplitude');
lay=tiledlayout(f,2,3,'Padding','compact','TileSpacing','compact');
show=[5 2 4 6 8]; case_name="strong_height_low_K";
q=d.points(d.points.case_id==case_name & d.points.pair=="PE_BIE" & d.points.quantity_kind=="G",:);
ax=nexttile(lay); hold(ax,'on');
plot(ax,real(d.line.cases(3).G_BIE(d.line.mask)),imag(d.line.cases(3).G_BIE(d.line.mask)),'.','Color',[.80 .80 .80]);
for k=show
    if q.valid(k)
        plot(ax,q.B_real(k),q.B_imag(k),'o','Color',colors(2,:),'MarkerFaceColor',colors(2,:));
        plot(ax,q.A_real(k),q.A_imag(k),'x','Color',colors(1,:),'MarkerSize',9);
        plot(ax,[q.A_real(k) q.B_real(k)],[q.A_imag(k) q.B_imag(k)],'k-');
        text(ax,q.A_real(k),q.A_imag(k)," "+ids(k));
    end
end
axis(ax,'equal'); grid(ax,'on'); xlabel(ax,'real G'); ylabel(ax,'imag G');
title(ax,'Comparable complex G; three saved points');
ax=nexttile(lay); hold(ax,'on');
plot(ax,1:5,q.A_abs(show),'x-','Color',colors(1,:));
plot(ax,1:5,q.B_abs(show),'o-','Color',colors(2,:));
set(ax,'XTick',1:5,'XTickLabel',cellstr(ids(show))); grid(ax,'on');
ylabel(ax,'|G|'); title(ax,'G magnitude; missing stays blank');
ax=nexttile(lay); hold(ax,'on');
plot(ax,1:5,angle(complex(q.A_real(show),q.A_imag(show))),'x-','Color',colors(1,:));
plot(ax,1:5,angle(complex(q.B_real(show),q.B_imag(show))),'o-','Color',colors(2,:));
set(ax,'XTick',1:5,'XTickLabel',cellstr(ids(show))); grid(ax,'on');
ylabel(ax,'Principal phase of G (rad)'); title(ax,'G phase; missing stays blank');
ax=nexttile(lay); hold(ax,'on');
nf=d.native_fields(d.native_fields.case_id==case_name,:);
for mi=1:2
    values=nan(5,1);
    for j=1:5
        row=nf(nf.receiver_id==ids(show(j)) & nf.method==ternary(mi==1,"PE","BIE"),:);
        values(j)=row.magnitude;
    end
    plot(ax,1:5,values,'-o','Color',colors(mi,:),'LineWidth',1.5);
end
set(ax,'XTick',1:5,'XTickLabel',cellstr(ids(show))); grid(ax,'on');
ylabel(ax,'Native |u|, shared Gaussian source scale'); title(ax,'Raw amplitude; all five saved');
ax=nexttile(lay); hold(ax,'on');
plot(ax,1:5,q.magnitude_db(show),'o-','Color',colors(1,:),'LineWidth',1.5);
set(ax,'XTick',1:5,'XTickLabel',cellstr(ids(show))); grid(ax,'on');
ylabel(ax,'20 log10(|G_PE|/|G_BIE|) dB'); title(ax,'G signed magnitude residual');
ax=nexttile(lay); hold(ax,'on');
plot(ax,1:5,q.phase_rad(show),'o-','Color',colors(3,:),'LineWidth',1.5);
set(ax,'XTick',1:5,'XTickLabel',cellstr(ids(show))); grid(ax,'on');
ylabel(ax,'Wrapped phase residual (rad)'); title(ax,'G phase residual; missing stays blank');
subtitle(lay,'F06: native rough PE/BIE phase is not used as a cross-method error. P12 and P32 lack G denominators.');
emit(f,figure_dir,'F06');

% F07: line, subset, raw-amplitude summaries have separate labels.
f=canvas('F07 | Saved-data summaries and exact magnitude-phase identity');
lay=tiledlayout(f,2,3,'Padding','compact','TileSpacing','compact');
line_rows=d.line_metrics(d.line_metrics.pair=="PE_BIE",:);
point_rows=d.summaries(d.summaries.pair=="PE_BIE" & d.summaries.quantity_kind=="G",:);
ax=nexttile(lay); hold(ax,'on');
bar(ax,1:3,[line_rows.complex_l2 point_rows.normalized_error]);
set(ax,'XTick',1:3,'XTickLabel',{'weak','medium','strong'}); grid(ax,'on');
legend(ax,{'M99 weighted E','3-point E_subset'},'Location','northwest');
ylabel(ax,'Normalized complex error'); title(ax,'Different supports, separate estimands');
ax=nexttile(lay); hold(ax,'on');
bar(ax,1:3,[line_rows.magnitude_fraction line_rows.phase_fraction],'stacked');
set(ax,'XTick',1:3,'XTickLabel',{'weak','medium','strong'}); ylim(ax,[0 1]);
legend(ax,{'Magnitude squared term','Phase squared term'},'Location','southwest');
ylabel(ax,'Fraction of weighted complex-difference square'); title(ax,'Exact algebraic decomposition');
ax=nexttile(lay); hold(ax,'on');
bar(ax,1:3,[point_rows.magnitude_fraction point_rows.phase_fraction],'stacked');
set(ax,'XTick',1:3,'XTickLabel',{'weak','medium','strong'}); ylim(ax,[0 1]);
ylabel(ax,'Fraction of local complex-difference square'); title(ax,'Local 3-point decomposition');
ax=nexttile(lay); hold(ax,'on');
bar(ax,1:3,[point_rows.mean_relative point_rows.median_relative point_rows.max_relative]);
set(ax,'XTick',1:3,'XTickLabel',{'weak','medium','strong'}); grid(ax,'on');
legend(ax,{'mean','median','max'},'Location','northwest');
ylabel(ax,'Pointwise relative G error, three z=3 samples'); title(ax,'Local subset, n=3 of 9');
ax=nexttile(lay); hold(ax,'on');
bar(ax,1:3,[point_rows.mean_absolute point_rows.median_absolute point_rows.max_absolute]);
set(ax,'XTick',1:3,'XTickLabel',{'weak','medium','strong'}); grid(ax,'on');
ylabel(ax,'|delta G|'); title(ax,'Absolute complex difference, n=3');
ax=nexttile(lay); hold(ax,'on');
raw_row=d.summaries(d.summaries.quantity_kind=="raw_native_amplitude_only",:);
bar(ax,[raw_row.mean_relative raw_row.median_relative raw_row.max_relative]);
set(ax,'XTick',1:3,'XTickLabel',{'mean','median','max'}); grid(ax,'on');
ylabel(ax,'Raw PE/BIE amplitude-only relative difference'); title(ax,'All 9 raw amplitudes; no phase ranking');
subtitle(lay,'F07: full 9-point G E_array is unavailable. Native raw amplitude and mapped z=3 raw rows are separate in point_summary.csv.');
emit(f,figure_dir,'F07');
files={'F01','F03','F04','F05','F05_grid','F06','F07'};
end

function f=canvas(title_text)
f=figure('Visible','off','Color','white','Position',[40 40 1500 850],'Name',title_text);
set(f,'DefaultAxesFontSize',11,'DefaultTextFontSize',11);
end
function emit(f,dir,name)
set(f,'PaperPositionMode','auto');
print(f,fullfile(dir,[name '.png']),'-dpng','-r140');
exportgraphics(f,fullfile(dir,[name '.pdf']),'ContentType','vector');
close(f);
end
function v=ternary(test,a,b), if test, v=a; else, v=b; end, end
