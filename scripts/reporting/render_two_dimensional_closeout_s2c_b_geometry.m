function files = render_two_dimensional_closeout_s2c_b_geometry(p0_mat,figure_dir)
%RENDER_TWO_DIMENSIONAL_CLOSEOUT_S2C_B_GEOMETRY Geometry-only M04a.
arguments
    p0_mat (1,:) char
    figure_dir (1,:) char
end
assert(isfile(p0_mat) && ~isfolder(figure_dir),'Use saved P0 MAT and a new figure directory');
t=load(p0_mat,'gx','eta0','d10','d20','eta1','d11','d21');
assert(all(t.eta1(abs(t.gx)>=50)==0));
mkdir(figure_dir);
f=figure('Visible','off','Color','w','Position',[80 80 1500 1050]);
layout=tiledlayout(f,3,1,'TileSpacing','compact','Padding','compact');
names={'Physical surface height \eta_z (m)','Surface slope d\eta_z/dx','Surface curvature d^2\eta_z/dx^2 (1/m)'};
base={t.eta0,t.d10,t.d20}; higher={t.eta1,t.d11,t.d21};
for j=1:3
    ax=nexttile(layout); hold(ax,'on');
    plot(ax,t.gx,base{j},'Color',[0 .35 .7],'LineWidth',1.5,'DisplayName','BASE_M K=0.10');
    plot(ax,t.gx,higher{j},'Color',[.8 .2 .05],'LineWidth',1.35,'LineStyle','--','DisplayName','K_HIGH K=0.20');
    for edge=[-50 -42 42 50]
        xline(ax,edge,':','Color',[.4 .4 .4],'LineWidth',.8,'HandleVisibility','off');
    end
    xlim(ax,[-60 60]); grid(ax,'on'); ylabel(ax,names{j},'Interpreter','tex');
    if j==1
        title(ax,'M04a PRE-EXPERIMENT GEOMETRY ONLY | C2 plateau |x|<=42 m, transition 42-50 m');
        legend(ax,'Location','northoutside','Orientation','horizontal','Interpreter','none');
    end
    if j==3, xlabel(ax,'surface x (m)'); end
end
drawnow;
png=fullfile(figure_dir,'M04a_geometry_only.png'); fig=fullfile(figure_dir,'M04a_geometry_only.fig');
assert(~isfile(png) && ~isfile(fig));
savefig(f,fig); exportgraphics(f,png,'Resolution',220); close(f);
files=struct('png',png,'fig',fig,'input_mat',p0_mat);
end
