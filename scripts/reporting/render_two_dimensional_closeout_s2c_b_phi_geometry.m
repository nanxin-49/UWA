function summary = render_two_dimensional_closeout_s2c_b_phi_geometry(root, out, source_mat)
%RENDER_TWO_DIMENSIONAL_CLOSEOUT_S2C_B_PHI_GEOMETRY Analytic M05a only.
% Reads the frozen S2c-A geometry; performs no acoustic field evaluation.
arguments
    root (1,:) char
    out (1,:) char
    source_mat (1,:) char = ''
end
assert(~isfolder(out),'Refusing to overwrite an existing figure folder');
if isempty(source_mat)
    prior=fullfile(root,'results','comparisons','two_dimensional_channel_closeout','20261023_s2c_a_1642');
    g=readtable(fullfile(prior,'protocol','candidate_geometry.csv'),TextType='string');
    a=g(g.case_id=="BASE_M",:); b=g(g.case_id=="PHI_PLUS90",:);
else
    frozen=load(source_mat,'a','b'); a=frozen.a; b=frozen.b;
    assert(all(a.case_id=="BASE_M") && all(b.case_id=="PHI_PLUS90"));
end
assert(height(a)==12001 && height(b)==12001);
assert(isequal(a.x_m,b.x_m));
x=a.x_m;
assert(max(abs(x-(-60:.01:60)'))<1e-10);
assert(all(a.eta_m(abs(x)>=50)==0) && all(b.eta_m(abs(x)>=50)==0));
assert(all(a.slope(abs(x)>=50)==0) && all(b.slope(abs(x)>=50)==0));
assert(all(a.curvature_per_m(abs(x)>=50)==0) && all(b.curvature_per_m(abs(x)>=50)==0));
i=find(x==0,1);
center=[a.eta_m(i),a.slope(i),a.curvature_per_m(i); ...
    b.eta_m(i),b.slope(i),b.curvature_per_m(i)];
assert(max(abs(center-[0,.005,0;.05,0,-.0005]),[],'all')<1e-12);
summary=table(["BASE_M";"PHI_PLUS90"],center(:,1),center(:,2),center(:,3), ...
    [max(abs(a.eta_m));max(abs(b.eta_m))], ...
    [max(abs(a.slope));max(abs(b.slope))], ...
    [max(abs(a.curvature_per_m));max(abs(b.curvature_per_m))], ...
    'VariableNames',{'case_id','center_eta_m','center_slope','center_curvature_per_m', ...
    'max_abs_eta_m','max_abs_slope','max_abs_curvature_per_m'});
mkdir(out);
writetable([a;b],fullfile(out,'geometry.csv'));
writetable(summary,fullfile(out,'geometry_summary.csv'));
save(fullfile(out,'geometry.mat'),'a','b','summary','-v7');
f=figure('Visible','off','Color','w','Position',[80 80 1100 900]);
t=tiledlayout(f,3,1,'TileSpacing','compact','Padding','compact');
names={'Surface height \eta_z (m)','Slope d\eta/dx','Curvature d^2\eta/dx^2 (m^{-1})'};
vars={'eta_m','slope','curvature_per_m'};
for k=1:3
    ax=nexttile(t,k); hold(ax,'on');
    h1=plot(ax,x,a.(vars{k}),'Color',[.12 .35 .75],'LineWidth',1.5);
    h2=plot(ax,x,b.(vars{k}),'Color',[.8 .24 .14],'LineWidth',1.5);
    v=xline(ax,0,':k','Receiver x=0','LabelVerticalAlignment','bottom');
    v.HandleVisibility='off';
    for q=[-50 -42 42 50]
        v=xline(ax,q,':','Color',[.6 .6 .6]); v.HandleVisibility='off';
    end
    ylabel(ax,names{k}); xlim(ax,[-60 60]); grid(ax,'on');
    if k==1
        h=plot(ax,0,center(1,1),'o','Color',[.12 .35 .75],'MarkerFaceColor',[.12 .35 .75]); h.HandleVisibility='off';
        h=plot(ax,0,center(2,1),'o','Color',[.8 .24 .14],'MarkerFaceColor',[.8 .24 .14]); h.HandleVisibility='off';
        legend(ax,[h1 h2],{'BASE_M: sine','PHI_PLUS90: cosine'},'Location','northeast','Interpreter','none');
    end
end
xlabel(nexttile(t,3),'Surface coordinate x (m)');
title(t,'M05a | frozen C2 geometry only; no PHI acoustic field');
drawnow;
savefig(f,fullfile(out,'M05a_geometry_only.fig'));
exportgraphics(f,fullfile(out,'M05a_geometry_only.png'),'Resolution',180);
close(f);
end
