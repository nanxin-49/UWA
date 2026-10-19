function output = generate_pe_bellhop_bie_xz_visuals(run_dir, figure_dir)
%GENERATE_PE_BELLHOP_BIE_XZ_VISUALS Raw complex fields on a true x-z grid.
% Pure saved-data renderer. Never invokes propagation, BIE, or Bellhop.
if nargin<2, figure_dir=fullfile(run_dir,'figures'); end
if ~isfolder(figure_dir), mkdir(figure_dir); end
t=load(fullfile(run_dir,'plot_data','xz_data.mat'),'data'); d=t.data;
assert(isequal(size(d.PE_carrier_restored),[411 45]));
methods={'PE_carrier_restored','BH_native','BIE_native'};
titles={{'PE Model-0','Analytic carrier restored'}, ...
    {'Bellhop','Native line-source units'},{'BIE','Gaussian source peak = 1'}};
names={'XZ01_raw_magnitude','XZ02_raw_phase','XZ03_complex_field','XZ04_raw_profiles'};
shared_max=max([abs(d.PE_carrier_restored(:));abs(d.BIE_native(:))]);
bh_max=max(abs(d.BH_native(:))); ceilings=[shared_max bh_max shared_max];
palette=[.43 .28 .64;0 .36 .62;.82 .32 .05];

f=canvas('Raw reflected field: amplitude in the x-z plane', ...
    'A=0.20 m | 4 kHz | 411 x 45 computed points | No rough/flat ratio');
lay=tiledlayout(f,1,3,'Padding','loose','TileSpacing','loose');
for j=1:3
    ax=nexttile(lay); imagesc(ax,d.x_m,d.z_m,abs(d.(methods{j})).');
    colormap(ax,parula(256)); clim(ax,[0 ceilings(j)]); cb=colorbar(ax,'southoutside');
    cb.Label.String='|u| (native source units)'; title(ax,titles{j},'FontSize',16); field_axes(ax);
end
foot(f,'PE/BIE share one amplitude scale; Bellhop retains its own source scale. Axes are expanded vertically. No fitted scaling.');
emit(f,figure_dir,names{1});

f=canvas('Raw reflected field: principal phase in the x-z plane', ...
    'Color = arg(u), rad | All panels use the same cyclic phase scale');
lay=tiledlayout(f,1,3,'Padding','loose','TileSpacing','loose');
for j=1:3
    ax=nexttile(lay); imagesc(ax,d.x_m,d.z_m,angle(d.(methods{j})).');
    colormap(ax,hsv(256)); clim(ax,[-pi pi]); cb=colorbar(ax,'southoutside');
    cb.Ticks=[-pi -pi/2 0 pi/2 pi]; cb.TickLabels={'-pi','-pi/2','0','pi/2','pi'};
    cb.Label.String='Principal phase (rad)'; title(ax,titles{j},'FontSize',16); field_axes(ax);
end
foot(f,'Native source phase references are retained; PE restores exp(+i*k0*(100+z)). These maps do not define a new phase-error ranking.');
emit(f,figure_dir,names{2});

f=canvas('Complex reflected field: phase and amplitude together', ...
    'Hue = principal phase | Brightness = linear |u| / displayed amplitude ceiling');
rgb=struct;
for j=1:3
    ax=axes(f,'Position',[.11+(j-1)*.28 .35 .22 .42]); u=d.(methods{j});
    hue=(angle(u)+pi)/(2*pi); value=min(abs(u)/ceilings(j),1);
    rgb.(methods{j})=hsv2rgb(cat(3,hue,ones(size(hue)),value));
    image(ax,d.x_m,d.z_m,permute(rgb.(methods{j}),[2 1 3]));
    title(ax,titles{j},'FontSize',16); field_axes(ax);
end
% A real color wheel, with amplitude varying radially, discloses both encodings.
ax=axes(f,'Position',[.43 .09 .14 .14]);
[wx,wy]=meshgrid(linspace(-1,1,151)); rr=hypot(wx,wy); hh=(atan2(wy,wx)+pi)/(2*pi);
wheel=hsv2rgb(cat(3,hh,ones(size(hh)),min(rr,1)));
for c=1:3, plane=wheel(:,:,c); plane(rr>1)=1; wheel(:,:,c)=plane; end
image(ax,[-1 1],[-1 1],wheel); set(ax,'YDir','normal'); axis(ax,'image'); axis(ax,'off');
text(ax,1.12,0,'0','FontSize',11); text(ax,-1.55,0,'+/-pi','FontSize',11);
text(ax,-.2,1.15,'pi/2','FontSize',11); text(ax,-.3,-1.18,'-pi/2','FontSize',11);
annotation(f,'textbox',[.10 .12 .32 .07],'String',sprintf('PE/BIE ceiling %.5g\nBellhop ceiling %.5g (native units)',shared_max,bh_max), ...
    'FontSize',14,'EdgeColor','none');
annotation(f,'textbox',[.61 .10 .32 .10],'String',{'Phase follows the hue angle', ...
    'Radius shows amplitude / ceiling','Center = zero; rim = ceiling'}, ...
    'FontSize',13,'EdgeColor','none');
foot(f,'Linear brightness, no gamma enhancement. The color wheel is a legend; every field pixel comes from a computed grid sample.');
emit(f,figure_dir,names{3});

f=canvas('Raw amplitude profiles at three physical depths', ...
    'Top: BIE purple solid / PE orange dash-dot | Bottom: Bellhop blue dashed, native units');
lay=tiledlayout(f,2,3,'Padding','loose','TileSpacing','loose'); depths=[1 3 6];
for j=1:3
    zi=find(abs(d.z_m-depths(j))<1e-12); assert(isscalar(zi));
    ax=nexttile(lay,j); hold(ax,'on');
    plot(ax,d.x_m,abs(d.BIE_native(:,zi)),'-','Color',palette(1,:),'LineWidth',2,'DisplayName','BIE');
    plot(ax,d.x_m,abs(d.PE_carrier_restored(:,zi)),'-.','Color',palette(3,:),'LineWidth',2,'DisplayName','PE Model-0');
    title(ax,sprintf('Depth z = %g m',depths(j))); ylabel(ax,'|u| (Gaussian source units)');
    ylim(ax,[0 shared_max*1.04]); style(ax);
    ax=nexttile(lay,j+3);
    plot(ax,d.x_m,abs(d.BH_native(:,zi)),'--','Color',palette(2,:),'LineWidth',2);
    ylabel(ax,'|u| (Bellhop native units)'); xlabel(ax,'Lateral position x (m)');
    ylim(ax,[0 bh_max*1.04]); style(ax);
end
foot(f,'Raw magnitudes, zero-based axes, shared limits within each row. Bellhop source amplitude is not calibrated to the Gaussian peak.');
emit(f,figure_dir,names{4});

transforms=struct('phase_rad',d.phase_rad,'magnitude',d.magnitude,'complex_rgb',rgb, ...
    'ceilings',ceilings,'phase_limits_rad',[-pi pi],'brightness_rule','linear abs(u)/ceiling; no gamma', ...
    'fields_use_ratio',false,'interpolation','none: imagesc displays computed samples');
save(fullfile(figure_dir,'plot_transforms.mat'),'transforms','-v7.3');
metadata=table(string(names(:)),repmat(2400,4,1),repmat(1350,4,1), ...
    'VariableNames',{'name','width','height'});
for j=1:4
    info=imfinfo(fullfile(figure_dir,[names{j} '.png']));
    assert(info.Width==2400 && info.Height==1350);
end
writetable(metadata,fullfile(figure_dir,'image_metadata.csv'));
output=struct('PNG_PDF_pairs',4,'new_solver_runs',0,'shared_PE_BIE_ceiling',shared_max, ...
    'BH_native_ceiling',bh_max,'figure_names',{names});
fid=fopen(fullfile(figure_dir,'render_audit.json'),'w'); cl=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(output,PrettyPrint=true)); clear cl
write_report(d,run_dir,names,ceilings); disp(output);
end

function write_report(d,out,names,ceilings)
fid=fopen(fullfile(out,'xz_field_report.md'),'w','n','UTF-8'); cl=onCleanup(@()fclose(fid));
fprintf(fid,'# PE–Bellhop–BIE：原始二维反射声场\n\n');
fprintf(fid,'这组新增图直接画复声场u(x,z)，不画粗糙/平面比G，不用误差指标替代声场。固定4 kHz、c=1500 m/s、Tx100 m、Gaussian sigma=.3 m、Bellhop X线源；强海面A=.20 m、K=.10 rad/m、C2平台42/支撑50 m。原有三案例比值图册保留，仍用于其既定定量比较。\n\n');
fprintf(fid,'|网格或计算|设置|\n|---|---|\n|横向|411点，x=%.8g至%.8g m，dx=.1953125 m|\n|深度|z=.5:.125:6 m，共45层|\n|真实二维采样|每种方法18495个位置，无图像拉伸伪造深度，无空间插值|\n|PE|复用海面反射初始场，FFT传播到新深度|\n|BIE|复用已保存4480节点边界密度，仅评价Green积分，边界方程新求解0次|\n|Bellhop|唯一新增10001波束二维网格作业，步长.05 m、4097点内部壁|\n\n',min(d.x_m),max(d.x_m));
fprintf(fid,'## 原始值、表示及色标\n\n');
fprintf(fid,'PE原始包络完整保存在PE_envelope；为显示波前，以解析因子exp(+i k0(100+z))恢复载波，保存PE_carrier_restored。该转换不改变幅值。BIE为原生exp(-i omega t)场；Bellhop为既定conj(raw SHD)原生场，原始SHD同样保留。没有按案例拟合幅度、扣相位、峰值归一化或构造G来冒充原始场。\n\n');
fprintf(fid,'原图册的固定表示转换作用于G：G_BIE=conj(u_BIE_rough/u_BIE_flat)，G_BH同理；它不能直接当作原始u的完整表示转换。这里保留原Model-0海面相位屏与各实现的场，载波恢复不重定义原相位屏。因此本组原始相位图与原G相位图是不同对象，不能把原G误差数值套到这些原始场上。\n\n');
fprintf(fid,'平面控制中，PE载波恢复后与BIE原生平面场相对复差%.8g；Bellhop/BIE原始平面中心幅值比%.8g、相位差%.8g rad，仅是源标度/参考相位诊断，没有把它作为修正系数应用。新彩图不能把这些原始源差别当作模型误差。\n\n', ...
    d.source_diagnostic.PE_restored_flat_vs_BIE_l2,d.source_diagnostic.BH_over_BIE_flat_center_abs,d.source_diagnostic.BH_over_BIE_flat_center_phase_rad);
fprintf(fid,'PE/BIE幅值共用[0,%.8g]；Bellhop原始源标度使用独立[0,%.8g]，均标为源单位，不声称Pa或绝对TL。全部相位图共用[-pi,pi]循环色标，保留各自原始源相位。相位图用来读空间波前，不能直接据色差重新评判原比值实验的误差大小。\n\n',ceilings(1),ceilings(2));
captions={'直接查看三种方法预测的反射场幅值在横向与深度上的分布。', ...
    '相位颜色展示波前条纹；循环色标的两端代表同一相位。', ...
    '同一像素的颜色给出相位，明暗给出幅值，暗色表示弱场而不是缺失。', ...
    '在1、3、6 m深度直接查看原始幅值曲线，避免仅凭彩图判断差异。'};
for j=1:4
    fprintf(fid,'## 图%d\n\n![二维声场图%d](figures/%s.png)\n\n%s\n\n[PDF](figures/%s.pdf)\n\n',j,j,names{j},captions{j},names{j});
    if j==1
        fprintf(fid,'每列对应一种方法，横轴x、纵轴深度z，深度向下。每个色块对应实际计算网格；为读清近海面结构，图轴在竖直方向作了显示上的展开，不表示实际几何比例。PE/BIE可按同一色标读取幅度；Bellhop应按自身色条读值，不能跨方法用颜色深浅估计绝对声压比。幅值的主体是同一声束包络，模型间的细小差别可能仍不显著，这不是缺少差异或应增加对比度的理由。\n\n');
    elseif j==2
        fprintf(fid,'相位是完整复声场的主值arg(u)，单位rad，不是相对于平面反射的相位增量。条纹展示相位随x和z变化；从+pi绕回-pi的颜色边界不代表场不连续。Bellhop保留其原始源参考相位，PE恢复了解析载波，因此应先看各自波前结构，不把整幅图的整体色移当作反射模型误差。\n\n');
    elseif j==3
        fprintf(fid,'色相h=(arg(u)+pi)/(2pi)，饱和度1，明度v=abs(u)/色条上限，线性映射、不作gamma增强。圆形图例中半径表示幅值，周向颜色表示相位。PE/BIE共用亮度上限；Bellhop上限单列注明。所有RGB数组、幅相数组与转换公式都保存在plot_transforms.mat，原始复数在xz_data.mat和CSV。\n\n');
    else
        fprintf(fid,'每列一个深度，上排直接叠加PE与BIE的原始幅值，下排显示Bellhop原始幅值。纵轴均从零开始，同一排使用共同范围；没有除以平面场。不同深度的曲线来自真实场点评价，3 m截线还与原图册原始作业核对。这里展示反射分量，不包含直达波与反射波叠加后的总场。\n\n');
    end
end
fprintf(fid,'## 数值核对与边界\n\n- PE新网格3 m截线与保存值相对差：%.8g。\n- BIE新评价3 m截线与保存值相对差：%.8g。\n- Bellhop新网格3 m截线与原10001波束作业相对差：%.8g。\n- 保存BIE的12/14 ppw两级密度在44个新域探针上的场变化：%.8g；这是本轮场评价敏感性，不是整个新域严格误差界。\n- 新Bellhop作业实际wall时间：%.6g s；不是方法间速度基准。\n\n', ...
    d.checks.PE_z3_closure,d.checks.BIE_z3_closure,d.checks.BH_z3_closure,d.checks.BIE_12_vs_14_ppw_probe_l2,d.BH_solver_wall_seconds);
fprintf(fid,'这是0.5–6 m近海面区域的二维反射场，未展示100 m全传播域、直达/反射总场、PM或三维。波前相似与原始幅值相似不替代误差检查，原生相位表示的差异也不能直接认定为物理模型误差。源模型、生产API、冻结求解器未改。\n\n');
fprintf(fid,'## 保存与独立重绘\n\n`plot_data/xz_data.mat`与`xz_fields.csv`保存原始复数、坐标、全部幅相与表示约定；`source_data`保存源文件哈希、数值边界密度、Bellhop输入输出与专用程序。纯绘图入口只读取xz_data.mat：\n\n```matlab\ngenerate_pe_bellhop_bie_xz_visuals(run_dir);\n```\n\n');
fprintf(fid,'工作流采用matlab与scientific-visualization技能。软件流程归因：Kassis T., Agarwal V., He Y., Patel D., Brueckner A.M. (2026), [Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents](https://doi.org/10.48550/arXiv.2609.00065)。已核对当前arXiv记录；此引用不是声学证据。\n');
end

function f=canvas(head,sub)
f=figure('Visible','off','Color','white','Position',[20 20 1600 900]);
set(f,'DefaultAxesFontName','Arial','DefaultTextFontName','Arial');
annotation(f,'textbox',[.04 .93 .92 .055],'String',head,'FontSize',25,'FontWeight','bold','EdgeColor','none','Interpreter','none');
annotation(f,'textbox',[.04 .875 .92 .04],'String',sub,'FontSize',16,'EdgeColor','none','Interpreter','none');
end

function field_axes(ax)
set(ax,'YDir','reverse'); xlabel(ax,'Lateral position x (m)'); ylabel(ax,'Depth z (m)'); style(ax);
end

function style(ax)
set(ax,'FontSize',15,'LineWidth',1,'Box','on','TickDir','out');
if isa(ax.Parent,'matlab.graphics.layout.TiledChartLayout'), ax.Parent.OuterPosition=[.045 .16 .91 .66]; end
end

function foot(f,str)
annotation(f,'textbox',[.04 .015 .92 .04],'String',str,'FontSize',13,'EdgeColor','none','Interpreter','none');
end

function emit(f,folder,name)
drawnow;
exportgraphics(f,fullfile(folder,[name '.png']),'Resolution',150,'BackgroundColor','white', ...
    'Units','inches','Width',16,'Height',9,'Padding','figure');
close(f);
% Encode the reviewed PNG as one PDF page. A second graphics export of the
% original tiled layout can reflow/clip labels in R2025b, even in image mode.
pixels=imread(fullfile(folder,[name '.png']));
page=figure('Visible','off','Color','white','Position',[20 20 1600 900]);
ax=axes(page,'Position',[0 0 1 1]); image(ax,pixels);
set(ax,'Visible','off','Position',[0 0 1 1], ...
    'XLim',[.5 size(pixels,2)+.5],'YLim',[.5 size(pixels,1)+.5]);
drawnow;
exportgraphics(page,fullfile(folder,[name '.pdf']),'ContentType','image','Resolution',150,'BackgroundColor','white', ...
    'Units','inches','Width',16,'Height',9,'Padding','figure');
close(page); fprintf('Rendered %s\n',name);
end
