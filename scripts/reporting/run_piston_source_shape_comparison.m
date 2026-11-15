function receipt=run_piston_source_shape_comparison(run_id,input_mat)
% Reuse saved inner fields; compute only the missing free-field outer strips.
arguments
    run_id (1,:) char
    input_mat (1,:) char
end
assert(~isempty(regexp(run_id,'^[A-Za-z0-9_]+$','once')));
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
out=fullfile(root,'results','validation','piston_source_shape_comparison',run_id);
assert(~isfolder(out)); mkdir(out); diary(fullfile(out,'diary.log'));
cleanup=onCleanup(@() diary('off')); %#ok<NASGU>
t0=tic; old=load(input_mat);
assert(isequal(double(old.x_m),-300:300) && isequal(old.z_m,0:.5:99.5));
assert(isequal(old.cfg.total_widths_m,[.2 .1]) && old.cfg.v0_peak_mps==.001);
x_m=-1000:1:1000; z_m=old.z_m; theta_deg=old.theta_deg;
f_hz=[4000 6000 8000]; radii_m=[97 200 500];
cfg=struct('c_mps',1500,'rho_kgpm3',1025,'v0_mps',.001,'z_source_m',100, ...
    'total_width_m',.1,'nquad',24,'input_mat',input_mat, ...
    'normalization','p(x,z)/p(0,100-R), R=hypot(x,100-z); equal-distance axial reference');
P=complex(zeros(numel(z_m),numel(x_m),3)); N=zeros(size(P));
Dstrip=zeros(numel(theta_deg),3); Dcircle=Dstrip;
pcircle=complex(zeros(numel(theta_deg),3,3));
half3=nan(2,3); nulls=nan(2,3); refine=zeros(1,3); stable=zeros(3,3);
inside=abs(x_m)<=300; outside=~inside; [XO,ZO]=meshgrid(x_m(outside),z_m);
[X,Z]=meshgrid(x_m,z_m); R=hypot(X,100-Z);
q3s=fzero(@(q) sin(q)/q-10^(-3/20),[.1 3]);
q3c=fzero(@(q) 2*besselj(1,q)/q-10^(-3/20),[.1 3]);
qnullc=fzero(@(q)besselj(1,q),[3 4]);
for fi=1:3
    src=struct('frequency_hz',f_hz(fi),'c_mps',1500,'rho_kgpm3',1025, ...
        'v0_mps',.001,'total_width_m',.1,'z_source_m',100);
    P(:,inside,fi)=old.p_field_f(:,:,2,fi);
    P(:,outside,fi)=baffled_strip_piston_field_vertical(XO,ZO,src,24);
    pref=baffled_strip_piston_field_vertical(zeros(size(R)),100-R,src,24);
    N(:,:,fi)=20*log10(abs(P(:,:,fi))./abs(pref));
    ix=[1 250 700 1500 2001]; iz=[1 70 200];
    test=baffled_strip_piston_field_vertical(X(iz,ix),Z(iz,ix),src,96);
    v=P(iz,ix,fi); refine(fi)=max(abs(test(:)-v(:)))/max(abs(test(:)));
    assert(refine(fi)<1e-7);
    ka=pi*.1*f_hz(fi)/1500; u=ka*sind(theta_deg); nz=u~=0;
    ds=ones(size(u)); dc=ds; ds(nz)=sin(u(nz))./u(nz); dc(nz)=2*besselj(1,u(nz))./u(nz);
    Dstrip(:,fi)=ds; Dcircle(:,fi)=dc;
    roots3=[q3s q3c]; roots0=[pi qnullc];
    for j=1:2
        if ka>=roots3(j), half3(j,fi)=asind(roots3(j)/ka); end
        if ka>=roots0(j), nulls(j,fi)=asind(roots0(j)/ka); end
    end
    for ri=1:3
        r=radii_m(ri);
        if ri==1, pcircle(:,ri,fi)=old.p_circle_f(:,2,fi);
        else, pcircle(:,ri,fi)=baffled_strip_piston_field_vertical(r*sind(theta_deg),100-r*cosd(theta_deg),src,96); end
        d=abs(pcircle(:,ri,fi)); d=d/d(theta_deg==0); stable(ri,fi)=max(abs(d-abs(ds)));
    end
end
assert(max(stable(:))<1e-3 && all(isfinite(P(:))));
% Only the 20 cm, 8 kHz reference can exhibit a forward first null here.
src.total_width_m=.2; src.frequency_hz=8000;
P20=complex(zeros(size(X))); P20(:,inside)=old.p_field_f(:,:,1,3);
P20(:,outside)=baffled_strip_piston_field_vertical(XO,ZO,src,24);
pref=baffled_strip_piston_field_vertical(zeros(size(R)),100-R,src,24);
N20=20*log10(abs(P20)./abs(pref));
save(fullfile(out,'data.mat'),'cfg','x_m','z_m','f_hz','theta_deg','P','N','P20','N20', ...
    'Dstrip','Dcircle','radii_m','pcircle','half3','nulls','refine','stable','-v7');
T=table(theta_deg,Dstrip(:,1),Dcircle(:,1),Dstrip(:,2),Dcircle(:,2),Dstrip(:,3),Dcircle(:,3), ...
    'VariableNames',{'theta_deg','strip4','circle4','strip6','circle6','strip8','circle8'});
writetable(T,fullfile(out,'theory_amplitude.csv'));
fig=figure('Visible','off','Color','w','Position',[30 30 1500 570]);
tl=tiledlayout(fig,1,3,'Padding','compact'); title(tl,'10 cm nominal aperture: 2-D strip vs 3-D circular piston');
subtitle(tl,'Different physical dimensions; axial-normalized FAR-FIELD amplitude. No forward first null for either source at 4/6/8 kHz.');
for fi=1:3
    ax=nexttile(tl); hold(ax,'on');
    plot(ax,theta_deg,20*log10(abs(Dstrip(:,fi))),'-','LineWidth',2,'DisplayName','2-D strip, total width 10 cm');
    plot(ax,theta_deg,20*log10(abs(Dcircle(:,fi))),'--','LineWidth',2,'DisplayName','3-D circle, diameter 10 cm');
    yline(ax,-3,'k:','-3 dB','HandleVisibility','off');
    for j=1:2
        if isfinite(half3(j,fi)), plot(ax,[-1 1]*half3(j,fi),[-3 -3],'o','HandleVisibility','off'); end
    end
    title(ax,sprintf('%d kHz | ka = %.3f',f_hz(fi)/1000,pi*.1*f_hz(fi)/1500));
    if fi<3, note='Neither curve reaches -3 dB in forward space';
    else, note=sprintf('-3 dB half-angle: strip %.1f deg, circle %.1f deg',half3(1,fi),half3(2,fi)); end
    text(ax,0,-12,note,'HorizontalAlignment','center','FontSize',9);
    xlabel(ax,'Angle from upward axis (deg)'); ylabel(ax,'Axis-normalized amplitude (dB)');
    xlim(ax,[-90 90]); ylim(ax,[-40 1]); grid(ax,'on');
    if fi==1, legend(ax,'Location','south','FontSize',9); end
end
saveplot(fig,'theory_10cm');
for span=[500 1000]
    fig=figure('Visible','off','Color','w','Position',[30 30 1500 850]);
    tl=tiledlayout(fig,2,3,'Padding','compact','TileSpacing','compact');
    title(tl,sprintf('10 cm total-width 2-D strip | x = +/- %d m | peak v_0 = 1 mm/s',span));
    subtitle(tl,'Top: absolute peak Pa. Bottom: equal-distance axial reference. Display aspect is stretched; slopes do not measure angles.');
    for row=1:2
        for fi=1:3
            ax=nexttile(tl); hold(ax,'on');
            if row==1, V=20*log10(abs(P(:,:,fi))); lim=[-10 60]; lab='dB re 1 Pa (peak)';
            else, V=N(:,:,fi); lim=[-40 0]; lab='20 log10(|p| / |p_axis(R)|), dB'; end
            imagesc(ax,x_m,z_m,V); clim(ax,lim); colormap(ax,parula(256)); cb=colorbar(ax); cb.Label.String=lab;
            set(ax,'YDir','reverse'); xlim(ax,[-span span]); ylim(ax,[0 100]);
            yline(ax,3,'w--','z=3 m','HandleVisibility','off'); plot(ax,0,100,'kv','MarkerFaceColor','w');
            xlabel(ax,'x (m)'); ylabel(ax,'z (m), down'); title(ax,sprintf('%d kHz',f_hz(fi)/1000));
        end
    end
    saveplot(fig,sprintf('strip10_far_%dm',span));
end
fig=figure('Visible','off','Color','w','Position',[30 30 1400 500]); tl=tiledlayout(fig,1,2,'Padding','compact');
title(tl,'8 kHz: equal-distance normalized exact 2-D fields | x = +/-1000 m');
subtitle(tl,'20 cm reference shows a null and weak outer lobe; 10 cm has no forward null. Color floor -40 dB; unequal display scales.');
for wi=1:2
    ax=nexttile(tl); if wi==1,V=N(:,:,3);w=10;else,V=N20;w=20;end
    imagesc(ax,x_m,z_m,V); set(ax,'YDir','reverse'); xlim(ax,[-1000 1000]); ylim(ax,[0 100]); clim(ax,[-40 0]);
    colorbar(ax); title(ax,sprintf('%d cm total width',w)); xlabel(ax,'x (m)');ylabel(ax,'z (m), down');
end
saveplot(fig,'strip_8kHz_width_reference');
fig=figure('Visible','off','Color','w','Position',[30 30 1500 550]);tl=tiledlayout(fig,1,3,'Padding','compact');
title(tl,'10 cm strip: exact equal-radius arcs and far-field sinc');
for fi=1:3
    ax=nexttile(tl); hold(ax,'on');
    styles={'-', '--', ':'};
    for ri=1:3
        d=abs(pcircle(:,ri,fi));d=d/d(theta_deg==0);
        plot(ax,theta_deg,20*log10(d),styles{ri},'LineWidth',1.7,'DisplayName',sprintf('r=%d m',radii_m(ri)));
    end
    ii=1:80:numel(theta_deg);plot(ax,theta_deg(ii),20*log10(abs(Dstrip(ii,fi))),'ko','MarkerSize',3,'DisplayName','far-field sinc');
    xlim(ax,[-90 90]);ylim(ax,[-6 0.2]);grid(ax,'on');title(ax,sprintf('%d kHz',f_hz(fi)/1000));xlabel(ax,'Angle from upward (deg)');ylabel(ax,'Axis-normalized dB');
    if fi==1,legend(ax,'Location','south');end
end
saveplot(fig,'strip10_radius_stability');
receipt=struct('status','COMPLETE_VISUAL_AUDIT_PENDING','elapsed_s',toc(t0),'release',version, ...
    'half_3db_deg',half3,'first_null_deg',nulls,'quadrature_relative',refine,'sinc_max_abs',stable, ...
    'reused_values',numel(z_m)*sum(inside)*3+numel(z_m)*sum(inside),'failed_checks',0);
fid=fopen(fullfile(out,'receipt.json'),'w');fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true));fclose(fid);disp(receipt);
    function saveplot(fig,stem)
        drawnow;savefig(fig,fullfile(out,[stem '.fig']));exportgraphics(fig,fullfile(out,[stem '.png']),'Resolution',160);close(fig);
    end
end
