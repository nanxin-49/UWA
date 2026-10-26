function summary=render_deepsea3d_5k7k_frequency()
%RENDER_DEEPSEA3D_5K7K_FREQUENCY Five discrete native response points.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
out=fullfile(root,'results','validation','deepsea3d_5k7k_20261026');
a=load(fullfile(out,'five_seven.mat'),'state');
b=load(fullfile(root,'results','validation','deepsea3d_g1_frequency','frequency.mat'),'state');
lock=load(fullfile(root,'results','validation','deepsea_dual_broadband_decision_20261015', ...
    'locked_candidate_and_gate.mat'),'decision','Q','predictions');
evidence=load(fullfile(root,'results','validation','deepsea_dual_broadband_decision_20261015', ...
    'saved_evidence.mat'),'output');
prior={b.state.cases(1),b.state.baseline6,b.state.cases(2)};
cases={prior{1},a.state.cases(1),prior{2},a.state.cases(2),prior{3}};
freq=[4000 5000 6000 7000 8000]';
R=complex(zeros(5,1));R0=R;D=R;H=R;Delta=R;FK=R;N=zeros(5,1);
for j=1:5
    c=cases{j};assert(c.cfg.frequency_hz==freq(j));
    assert(c.reference_ready && c.complete && c.G1_status=="CONDITIONAL_G1_PASS");
    ref=c.runs(c.reference_index).result;
    bg=deepsea3d_flat(c.cfg,0);
    R(j)=ref.reflect_native;R0(j)=bg.reflect_native;D(j)=bg.direct_native;
    Delta(j)=R(j)-R0(j);H(j)=D(j)+R(j);FK(j)=c.fk(2).reflect_native;
    N(j)=ref.order;
    assert(abs(ref.rough_correction-Delta(j))<1e-16);
    assert(abs(ref.total_native-H(j))<1e-16);
end
eta=evidence.output.geometry.surface_height_at_flat_specular_m;
pred=complex(nan(5,1));errR=nan(5,1);errDelta=errR;errH=errR;
lockError=zeros(2,1);
for n=1:2
    j=2*n; f=freq(j); k=n+1;
    cfg=cases{j}.cfg;B=deepsea3d_flat(cfg,eta);
    qmid=(lock.Q(k-1)+lock.Q(k))/2;
    reconstructed=qmid*B.reflect_native;
    saved=complex(lock.predictions(n).reflect_native_real,lock.predictions(n).reflect_native_imag);
    lockError(n)=abs(saved-reconstructed);
    assert(lockError(n)<1e-15);
    pred(j)=saved;
    errR(j)=abs(saved-R(j))/abs(R(j));
    errDelta(j)=abs(saved-R(j))/abs(Delta(j));
    errH(j)=abs(saved-R(j))/abs(H(j));
end
summary=table(freq,N,real(R),imag(R),real(R0),imag(R0),real(D),imag(D), ...
    real(Delta),imag(Delta),real(H),imag(H),real(FK),imag(FK), ...
    real(pred),imag(pred),errR,errDelta,errH, ...
    'VariableNames',{'frequency_hz','FE_order','R_real','R_imag','R0_real','R0_imag', ...
    'D_real','D_imag','Delta_real','Delta_imag','H_real','H_imag', ...
    'FK_R_real','FK_R_imag','locked_pred_R_real','locked_pred_R_imag', ...
    'reconstruction_R_relative','reconstruction_Delta_relative','reconstruction_H_relative'});
dataFile=fullfile(out,'five_frequency_native_summary_v2.csv');
figFile=fullfile(out,'five_frequency_reflection_and_error_v2.fig');
pngFile=fullfile(out,'five_frequency_reflection_and_error_v2.png');
assert(~isfile(dataFile) && ~isfile(figFile) && ~isfile(pngFile));
writetable(summary,dataFile);
fig=figure('Visible','off','Color','w','Position',[100 100 1180 490]);
t=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(t,1);hold(ax,'on');grid(ax,'on');axis(ax,'equal');
plot(ax,real(R([1 3 5]))*1e5,imag(R([1 3 5]))*1e5,'ko', ...
    'MarkerFaceColor',[.2 .2 .2],'MarkerSize',8,'LineStyle','none','DisplayName','4/6/8 kHz reference');
plot(ax,real(R([2 4]))*1e5,imag(R([2 4]))*1e5,'bs', ...
    'MarkerFaceColor',[.1 .45 .8],'MarkerSize',9,'LineStyle','none','DisplayName','5/7 kHz reference');
plot(ax,real(pred([2 4]))*1e5,imag(pred([2 4]))*1e5,'rx', ...
    'MarkerSize',10,'LineWidth',1.7,'LineStyle','none','DisplayName','Locked prediction');
for j=1:5,text(ax,real(R(j))*1e5,imag(R(j))*1e5,sprintf('  %g kHz',freq(j)/1000), ...
        'FontSize',9,'VerticalAlignment','bottom');end
xlim(ax,[-5.8 5.8]);ylim(ax,[-5.8 5.8]);
xlabel(ax,'Re(R) (10^{-5} m^{-1})');ylabel(ax,'Im(R) (10^{-5} m^{-1})');
title(ax,'Native complex reflection: five discrete points');legend(ax,'Location','southoutside');
ax=nexttile(t,2);hold(ax,'on');grid(ax,'on');
plot(ax,freq([2 4])/1000,errR([2 4]),'bs-','LineWidth',1.4,'MarkerFaceColor',[.1 .45 .8], ...
    'DisplayName','|R_{pred}-R|/|R|');
plot(ax,freq([2 4])/1000,errDelta([2 4]),'ko--','LineWidth',1.4, ...
    'DisplayName','|R_{pred}-R|/|\Delta|');
plot(ax,freq([2 4])/1000,errH([2 4]),'rd:','LineWidth',1.4, ...
    'DisplayName','|R_{pred}-R|/|H|');
yline(ax,1e-3,'Color',[.4 .4 .4],'LineStyle','--','DisplayName','Locked 10^{-3} screen');
set(ax,'YScale','log','XTick',[5 7],'XLim',[4.6 7.4]);
xlabel(ax,'Frequency (kHz)');ylabel(ax,'Relative complex reconstruction error');
title(ax,'Independent 5/7 kHz reconstruction checks');
legend(ax,'Location','southoutside');
title(t,'Discrete markers are physical cases; display lines do not certify an interval');
drawnow;savefig(fig,figFile);exportgraphics(fig,pngFile,'Resolution',180);close(fig);
disp(summary(:,{'frequency_hz','FE_order','reconstruction_R_relative', ...
    'reconstruction_Delta_relative','reconstruction_H_relative'}));
disp(lockError);
end
