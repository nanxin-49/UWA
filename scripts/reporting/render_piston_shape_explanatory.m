function receipt=render_piston_shape_explanatory(run_id,input_mat)
% Four explanatory figures; only missing angle/radius samples are computed.
arguments
    run_id (1,:) char
    input_mat (1,:) char
end
assert(~isempty(regexp(run_id,'^[A-Za-z0-9_]+$','once')));
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
out=fullfile(root,'results','visualization','piston_shape_explanatory',run_id);
assert(~isfolder(out));mkdir(out);diary(fullfile(out,'diary.log'));
cleanup=onCleanup(@()diary('off')); %#ok<NASGU>
t0=tic;s=load(input_mat); theta_deg=s.theta_deg; freqs=s.f_hz;
colors=[0 .35 .68;.85 .33 .05];
fig=figure('Visible','off','Color','w','Position',[30 30 1450 580]);
tl=tiledlayout(fig,1,3,'Padding','compact','TileSpacing','compact');
title(tl,'10 cm nominal aperture | strip vs circular piston | normalized far-field amplitude');
subtitle(tl,'Solid: 2-D strip (total width 10 cm). Dashed: 3-D circle (diameter 10 cm). No forward first null at these frequencies.');
for fi=1:3
    ax=nexttile(tl);hold(ax,'on');
    plot(ax,theta_deg,20*log10(abs(s.Dstrip(:,fi))),'-','Color',colors(1,:),'LineWidth',2,'DisplayName','2-D strip');
    plot(ax,theta_deg,20*log10(abs(s.Dcircle(:,fi))),'--','Color',colors(2,:),'LineWidth',2,'DisplayName','3-D circle');
    yline(ax,-3,'k:','-3 dB','HandleVisibility','off');
    title(ax,sprintf('%d kHz',freqs(fi)/1000));xlim(ax,[-90 90]);ylim(ax,[-6 0]);grid(ax,'on');
    xlabel(ax,'Angle from upward axis (deg)');ylabel(ax,'20 log_{10}|D| (dB)');
    if fi==3
        text(ax,0,-5.4,'-3 dB half-angles: 56.0 / 74.4 deg','HorizontalAlignment','center','FontSize',10);
    else
        text(ax,0,-5.4,'No forward -3 dB crossing','HorizontalAlignment','center','FontSize',10);
    end
    if fi==1,legend(ax,'Location','southwest');end
end
saveplot(fig,'01_theory_cartesian');
fig=figure('Visible','off','Color','w','Position',[30 30 1450 650]);
tl=tiledlayout(fig,1,3,'Padding','compact','TileSpacing','compact');
title(tl,'10 cm nominal aperture | forward polar directivity (linear amplitude radius)');
subtitle(tl,'Upward = 0 deg. Radius = |D|, not dB. Dotted semicircle = -3 dB threshold; sound continues outside this threshold.');
for fi=1:3
    ax=polaraxes(tl); ax.Layout.Tile=fi;hold(ax,'on');
    polarplot(ax,deg2rad(theta_deg),abs(s.Dstrip(:,fi)),'-','Color',colors(1,:),'LineWidth',2,'DisplayName','2-D strip');
    polarplot(ax,deg2rad(theta_deg),abs(s.Dcircle(:,fi)),'--','Color',colors(2,:),'LineWidth',2,'DisplayName','3-D circle');
    polarplot(ax,deg2rad(theta_deg),10^(-3/20)*ones(size(theta_deg)),'k:','HandleVisibility','off');
    ax.ThetaZeroLocation='top';ax.ThetaDir='clockwise';thetalim(ax,[-90 90]);rlim(ax,[0 1]);rticks(ax,[0 .25 .5 .75 1]);
    title(ax,sprintf('%d kHz',freqs(fi)/1000));if fi==1,legend(ax,'Location','southoutside');end
end
saveplot(fig,'02_theory_polar');
radius_m=unique([logspace(log10(.5),log10(500),65) 97 200 500]);
[TH,RR]=meshgrid(theta_deg,radius_m);P=complex(zeros([size(TH) 2]));
normalized_dB=zeros(size(P)); reuse_count=0;
for wi=1:2
    src=struct('frequency_hz',8000,'c_mps',1500,'rho_kgpm3',1025, ...
        'v0_mps',.001,'total_width_m',[.1 .2]*[wi==1;wi==2],'z_source_m',100);
    missing=true(size(radius_m));
    if wi==1
        for ri=1:3
            j=find(radius_m==s.radii_m(ri));P(j,:,wi)=s.pcircle(:,ri,3).';
            missing(j)=false;reuse_count=reuse_count+numel(theta_deg);
        end
    end
    P(missing,:,wi)=baffled_strip_piston_field_vertical(RR(missing,:).*sind(TH(missing,:)), ...
        100-RR(missing,:).*cosd(TH(missing,:)),src,96);
    axis_pressure=baffled_strip_piston_field_vertical(zeros(size(radius_m)),100-radius_m,src,96);
    normalized_dB(:,:,wi)=20*log10(max(abs(P(:,:,wi))./abs(axis_pressure(:)),realmin));
end
fig=figure('Visible','off','Color','w','Position',[30 30 1100 740]);
tl=tiledlayout(fig,1,1,'Padding','compact');
title(tl,'10 cm 2-D strip | 8 kHz | exact aperture pressure');
subtitle(tl,'Equal-distance axial normalization: 20 log10(|p(r,theta)|/|p(r,0)|). Distance uses log scale.');
ax=nexttile(tl);map(ax,normalized_dB(:,:,1),[-6 0]);hold(ax,'on');
[C,h]=contour(ax,theta_deg,radius_m,normalized_dB(:,:,1),[-4 -3 -2 -1],'k','LineWidth',1.1);
clabel(C,h,'FontSize',10,'Color','k');title(ax,'Labeled contours are amplitude levels, not beam boundaries');
saveplot(fig,'03_strip10_angle_distance');
fig=figure('Visible','off','Color','w','Position',[30 30 1450 740]);
tl=tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');
title(tl,'8 kHz | exact 2-D strip fields | same normalization and -40 to 0 dB color scale');
subtitle(tl,'10 cm: broad main lobe, no forward null. 20 cm: deep minimum near +/-69.64 deg, weak outer-lobe segment. Distance uses log scale.');
for wi=1:2
    ax=nexttile(tl);map(ax,normalized_dB(:,:,wi),[-40 0]);hold(ax,'on');
    title(ax,sprintf('%d cm total strip width',[10 20]*[wi==1;wi==2]));
    if wi==2
        null_angle=asind(1500/(8000*.2));
        for sg=[-1 1]
            plot(ax,sg*null_angle*[1 1],[10 500],'w:','LineWidth',1.1);
        end
        text(ax,0,.7,'White guides: FAR-FIELD first-null angle only', ...
            'HorizontalAlignment','center','FontSize',9,'BackgroundColor','w');
    end
end
saveplot(fig,'04_strip_width_angle_distance');
plot_cfg=struct('radius_m',radius_m,'angles_deg',theta_deg,'normalization','exact equal-radius axis pressure', ...
    'quadrature_order',96,'near_radius_min_m',.5,'new_physics_validation',false);
save(fullfile(out,'plot_data.mat'),'plot_cfg','P','normalized_dB','radius_m','theta_deg','reuse_count','-v7');
receipt=struct('status','RENDERED_VISUAL_REVIEW_PENDING','elapsed_s',toc(t0),'release',version, ...
    'reused_circle_samples',reuse_count,'figure_count',4,'new_model_validation',false);
fid=fopen(fullfile(out,'receipt.json'),'w');fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true));fclose(fid);
    function map(ax,V,limits)
        contourf(ax,theta_deg,radius_m,min(max(V,limits(1)),limits(2)),linspace(limits(1),limits(2),81),'LineStyle','none');
        set(ax,'YScale','log','YDir','normal','FontSize',11);xlim(ax,[-90 90]);ylim(ax,[.5 500]);
        clim(ax,limits);colormap(ax,parula(256));cb=colorbar(ax);cb.Label.String='Equal-radius normalized amplitude (dB)';
        xlabel(ax,'Angle from upward axis (deg)');ylabel(ax,'Distance from source r (m), log scale');
    end
    function saveplot(fig,stem)
        drawnow;savefig(fig,fullfile(out,[stem '.fig']));exportgraphics(fig,fullfile(out,[stem '.png']),'Resolution',160);close(fig);
    end
end
