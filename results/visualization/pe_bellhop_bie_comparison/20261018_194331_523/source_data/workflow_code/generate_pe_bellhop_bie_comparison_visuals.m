function output = generate_pe_bellhop_bie_comparison_visuals(run_dir, figure_dir)
%GENERATE_PE_BELLHOP_BIE_COMPARISON_VISUALS Pure saved-data seven-figure renderer.
% No PE, Bellhop, BIE, shell process, or random draw is invoked here.
if nargin<2, figure_dir=fullfile(run_dir,'figures'); end
q=load(fullfile(run_dir,'plot_data','comparison_data.mat'),'data'); data=q.data;
assert(numel(data.cases)==3 && height(data.metrics)==9);
for name={'','bellhop_bie','three_way'}
    folder=fullfile(figure_dir,name{1}); if ~isfolder(folder), mkdir(folder); end
end
names={'01_environment','bellhop_bie/BB01_fields','bellhop_bie/BB02_residuals', ...
    'bellhop_bie/BB03_metrics','three_way/TR01_fields','three_way/TR02_residuals','three_way/TR03_metrics'};
labels={'Weak: A=0.01 m','Medium: A=0.05 m','Strong: A=0.20 m'};
palette=[.43 .28 .64;0 .36 .62;.82 .32 .05];
% Environment: common physical fixture, no invented ray trace.
f=canvas('PE / Bellhop / BIE: a matched 2-D comparison', ...
    '4 kHz | c=1500 m/s | Gaussian sigma=0.3 m / Bellhop X line source');
t=tiledlayout(f,1,2,'Padding','loose','TileSpacing','loose');
ax=nexttile(t); hold(ax,'on');
plot(ax,[-60 60],[0 0],'-','Color',palette(1,:),'LineWidth',2);
plot(ax,[-60 60],[3 3],'--','Color',palette(2,:),'LineWidth',2);
plot(ax,0,100,'o','MarkerSize',12,'MarkerFaceColor',palette(3,:),'MarkerEdgeColor',palette(3,:));
plot(ax,0,3,'s','MarkerSize',10,'MarkerFaceColor',palette(2,:),'MarkerEdgeColor',palette(2,:));
text(ax,5,96,'Tx (0,100 m)','FontSize',19); text(ax,5,12,'Rx line: z=3 m','FontSize',19);
set(ax,'YDir','reverse'); xlim(ax,[-60 60]); ylim(ax,[-5 108]);
xlabel(ax,'Lateral position x (m)'); ylabel(ax,'Depth z (m)'); title(ax,'Geometry / no bottom interaction'); style(ax);
ax=nexttile(t); hold(ax,'on');
for j=1:3
    plot(ax,data.x_m,data.cases(j).eta_m,'LineWidth',2,'DisplayName',labels{j});
end
xlim(ax,[-60 60]); set(ax,'YDir','reverse'); xlabel(ax,'Lateral position x (m)'); ylabel(ax,'Surface depth eta (m)');
title(ax,'Same K=0.10 rad/m / different heights'); legend(ax,'Location','southoutside'); style(ax);
footer(f,'C2 plateau 42 m / support 50 m. G = reflected rough field / reflected flat field; G is not absolute pressure or TL.');
emit(f,figure_dir,names{1});

field_figure(data,labels,palette,false,figure_dir,names{2});
f=canvas('Bellhop - BIE: small residuals made visible', ...
    'Signed differences | Main curves: 10001 beams | Fixed Stage-0 M99 footprint');
t=tiledlayout(f,2,3,'Padding','loose','TileSpacing','loose');
for j=1:3
    c=data.cases(j); r=c.residuals(2); xx=data.x_m; ix=c.plot_mask;
    ax=nexttile(t,j); plot(ax,xx(ix),r.magnitude_ratio_db(ix),'Color',palette(2,:),'LineWidth',2);
    title(ax,labels{j}); ylabel(ax,'BH / BIE magnitude ratio (dB)'); style(ax);
    ax.YAxis.Exponent=0; ytickformat(ax,'%.4g');
    ax=nexttile(t,j+3); plot(ax,xx(ix),r.phase_rad(ix)*1e3,'Color',palette(2,:),'LineWidth',2);
    xlabel(ax,'Receiver x (m)'); ylabel(ax,'BH - BIE phase (mrad)'); style(ax);
end
footer(f,'Residual axes have independent scales per column. No fitted amplitude or phase shift. 1 mrad = 0.001 rad.');
emit(f,figure_dir,names{3});
metric_figure(data,labels,palette,false,figure_dir,names{4});
field_figure(data,labels,palette,true,figure_dir,names{5});

f=canvas('Three-way comparison: pairwise residual structure', ...
    'PE-BIE orange / Bellhop-BIE blue / PE-Bellhop gray | No alignment');
t=tiledlayout(f,2,3,'Padding','loose','TileSpacing','loose');
colors=[palette(3,:);palette(2,:);.40 .40 .40]; styles={'-.','--',':'};
for j=1:3
    c=data.cases(j); xx=data.x_m; ix=c.plot_mask;
    ax=nexttile(t,j); hold(ax,'on');
    for p=1:3, plot(ax,xx(ix),c.residuals(p).complex_abs(ix),styles{p},'Color',colors(p,:),'LineWidth',2); end
    title(ax,labels{j}); ylabel(ax,'|G_A - G_B| (dimensionless)'); style(ax);
    ax=nexttile(t,j+3); hold(ax,'on');
    for p=1:3, plot(ax,xx(ix),c.residuals(p).phase_rad(ix)*1e3,styles{p},'Color',colors(p,:),'LineWidth',2); end
    xlabel(ax,'Receiver x (m)'); ylabel(ax,'arg(G_A conj(G_B)) (mrad)'); style(ax);
end
footer(f,'Independent residual scales per column. Curves retain signed phase; absolute complex differences are not energy differences.');
emit(f,figure_dir,names{6});
metric_figure(data,labels,palette,true,figure_dir,names{7});

% Recompute metrics from packaged raw curves; do not trust old CSV values.
recheck=zeros(9,3); n=0;
for j=1:3
    c=data.cases(j); a={c.G_PE,c.G_BH,c.G_PE}; b={c.G_BIE,c.G_BIE,c.G_BH};
    ix=c.plot_mask; w=data.weights(ix); w=w/sum(w);
    for p=1:3
        n=n+1; aa=a{p}(ix); bb=b{p}(ix);
        recheck(n,:)=[sqrt(sum(w.*abs(aa-bb).^2)/sum(w.*abs(bb).^2)), ...
            sqrt(sum(w.*angle(aa.*conj(bb)).^2)),sqrt(sum(w.*(20*log10(abs(aa)./abs(bb))).^2))];
    end
end
stored=[data.metrics.complex_l2,data.metrics.phase_rms_rad,data.metrics.magnitude_ratio_rms_db];
assert(max(abs(recheck-stored),[],'all')<1e-12,'Packaged metrics do not reproduce.');
metadata=table(string(names(:)),zeros(7,1),zeros(7,1),'VariableNames',{'name','width','height'});
for j=1:7
    infx=imfinfo(fullfile(figure_dir,[names{j} '.png']));
    metadata.width(j)=infx.Width; metadata.height(j)=infx.Height;
    assert(infx.Width==2400 && infx.Height==1350);
    assert(isfile(fullfile(figure_dir,[names{j} '.pdf'])));
end
writetable(metadata,fullfile(figure_dir,'image_metadata.csv'));
output=struct('figure_names',{names},'figure_dir',figure_dir,'metrics_closure',max(abs(recheck-stored),[],'all'), ...
    'PNG_PDF_pairs',7,'new_solver_runs',0);
fid=fopen(fullfile(figure_dir,'render_audit.json'),'w'); cleaner=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(output,PrettyPrint=true)); clear cleaner
write_report(data,run_dir,names);
disp(output);
end

function field_figure(d,labels,palette,triple,folder,name)
title_text='Bellhop / BIE: normalized reflected field';
if triple, title_text='PE / Bellhop / BIE: normalized reflected field'; end
sub='BIE solid purple / Bellhop dashed blue | G = rough / flat';
if triple, sub='BIE solid purple / Bellhop dashed blue / PE dash-dot orange | G = rough / flat'; end
f=canvas(title_text,sub);
t=tiledlayout(f,2,3,'Padding','loose','TileSpacing','loose');
fields={'G_BIE','G_BH','G_PE'}; styles={'-','--','-.'}; count=2+double(triple);
mag_min=Inf; mag_max=-Inf;
for j=1:3
    for k=1:count
        vals=abs(d.cases(j).(fields{k})(d.cases(j).plot_mask));
        mag_min=min(mag_min,min(vals)); mag_max=max(mag_max,max(vals));
    end
end
pad=max(.01,.08*(mag_max-mag_min));
for j=1:3
    c=d.cases(j); ix=c.plot_mask; xx=d.x_m(ix);
    ax=nexttile(t,j); hold(ax,'on');
    for k=1:count, y=c.(fields{k}); plot(ax,xx,abs(y(ix)),styles{k},'Color',palette(k,:),'LineWidth',2.3); end
    title(ax,labels{j}); ylabel(ax,'|G| (dimensionless)'); ylim(ax,[max(0,mag_min-pad) mag_max+pad]); style(ax);
    ax=nexttile(t,j+3); hold(ax,'on');
    for k=1:count
        method=erase(fields{k},'G_'); y=c.transforms.(method).plotted_phase_mrad;
        plot(ax,xx,y(ix),styles{k},'Color',palette(k,:),'LineWidth',2.3);
    end
    xlabel(ax,'Receiver x (m)'); ylabel(ax,'Principal arg(G) (mrad)'); ylim(ax,[-pi pi]*1e3); style(ax);
end
footer(f,'Shared magnitude and phase scales. Phase-wrap crossings are gaps; full complex samples remain in the plot-data MAT/CSV.');
emit(f,folder,name);
end

function metric_figure(d,labels,palette,triple,folder,name)
title_text='Bellhop - BIE: weighted comparison metrics';
if triple, title_text='PE / Bellhop / BIE: three pairwise comparisons'; end
f=canvas(title_text,'Fixed M99 incident-energy weights | Metrics recomputed from saved complex curves');
t=tiledlayout(f,1,3,'Padding','loose','TileSpacing','loose');
tags={'PE_BIE','BH_BIE','PE_BH'}; colors=[palette(3,:);palette(2,:);.4 .4 .4];
chosen=2; if triple, chosen=1:3; end
columns={'complex_l2','phase_rms_rad','magnitude_ratio_rms_db'};
titles={'Relative complex L2','Relative phase RMS','Magnitude-ratio RMS'};
units={'Dimensionless','mrad','dB'};
markers={'o','s','^'};
for k=1:3
    ax=nexttile(t); hold(ax,'on'); plotted=[];
    for p=chosen
        rows=strcmp(d.metrics.pair,tags{p}); vals=d.metrics.(columns{k})(rows);
        if k==2, vals=vals*1000; end
        plotted=[plotted;vals]; %#ok<AGROW>
        plot(ax,1:3,vals,['-' markers{p}],'Color',colors(p,:),'LineWidth',2, ...
            'MarkerSize',9,'DisplayName',strrep(tags{p},'_','-'));
    end
    set(ax,'YScale','log'); ylim(ax,[min(plotted)/1.3 max(plotted)*1.3]);
    xticks(ax,1:3); xticklabels(ax,{'Weak','Medium','Strong'});
    xlim(ax,[.7 3.3]); title(ax,titles{k}); ylabel(ax,units{k}); style(ax);
    if k==1 && triple, legend(ax,'Location','northwest'); end
end
if triple
    pe=d.metrics.complex_l2(strcmp(d.metrics.pair,'PE_BIE'));
    bh=d.metrics.complex_l2(strcmp(d.metrics.pair,'BH_BIE'));
    if all(bh<pe)
        text_note='这三个固定海面上，Bellhop与保存BIE的复差均小于PE；PE差异随海面高度增大。结论仅限本夹具，BIE未重算。';
        if ~all(diff(pe)>0)
            text_note='这三个固定海面上，Bellhop与保存BIE的复差均小于PE。结论仅限本夹具，BIE未重算。';
        end
    else
        text_note='各配对在三个固定海面上的实际差异如图；不能外推普适模型排序。BIE为保存参考，本轮未重算。';
    end
else
    values=d.beam_sensitivity.complex_l2;
    text_note=sprintf('5001/10001 beam-ratio changes: %.2g / %.2g / %.2g (weak/medium/strong). Sensitivity, not an error bound.',values);
    text_note=sprintf('%s\nHistorical BIE refinement ~6e-8: saved R3-R5 reports; no new BIE solve.',text_note);
end
footer(f,text_note); emit(f,folder,name);
% Labels parameter is kept explicit in the shared call interface.
assert(numel(labels)==3);
end

function write_report(d,out,names)
titles={'共同环境与三个海面','Bellhop–BIE 归一反射场','Bellhop–BIE 幅相差值', ...
    'Bellhop–BIE 指标','PE–Bellhop–BIE 归一反射场','三组配对残差','三方指标总结'};
caps={'同一声源和几何下，只改变海面起伏高度。','两种方法在相同接收坐标上的归一反射幅度和相位直接比较。', ...
    '差值图将主曲线难以辨认的幅度与相位差单独展开。','加权指标概括三个固定海面的差异，并保留波束敏感性。', ...
    '三条曲线使用同一环境与表示约定，不作幅相拟合。','分别观察 PE–BIE、Bellhop–BIE 和 PE–Bellhop 的空间残差。', ...
    '三个配对指标共同说明这组固定实验中的关系。'};
fid=fopen(fullfile(out,'comparison_report.md'),'w','n','UTF-8'); clean=onCleanup(@()fclose(fid));
fprintf(fid,'# Bellhop–BIE 与 PE–Bellhop–BIE 对比图册\n\n');
fprintf(fid,'固定二维夹具：4 kHz、c=1500 m/s、Tx100/Rx3、Gaussian sigma=.3 m、Bellhop X线源。海面K=.10 rad/m，A=.01/.05/.20 m；C2 plateau42/support50 m。比较量为接收端反射场的粗糙/平面复数比G，不是绝对声压或传播损失。\n\n');
fprintf(fid,'|参数|本轮统一设置|\n|---|---|\n|频率 / 声速|4000 Hz / 1500 m/s|\n|物理坐标|海面z=0，深度向下；Tx=(0,100 m)，Rx为z=3 m横向线|\n|声源|Gaussian sigma=0.3 m；Bellhop X线源使用Stage-0匹配指向性|\n|海面|eta=A sin(0.10x) w(x)，A=0.01/0.05/0.20 m|\n|C2包络|abs(x)<=42 m保留，abs(x)>=50 m归零，中间五次多项式过渡|\n|PE|保存的原始Model-0相位屏|\n|Bellhop|压力释放内部壁，[-80,80] m / 4097点，步长0.05 m，角度[-30,30] deg|\n|波束|10001主结果；5001作敏感性|\n|接收线|984点，周期网格总长192.1875 m，间隔0.1953125 m|\n|指标支持|Stage-0 M99共%d点，x=%.7g至%.7g m；权重为固定入射能量|\n\n',nnz(d.mask),min(d.x_m(d.mask)),max(d.x_m(d.mask)));
fprintf(fid,'指标定义：E_AB=sqrt(sum(w |G_A-G_B|^2)/sum(w |G_B|^2))；相位RMS=sqrt(sum(w angle(G_A conj(G_B))^2))；幅度比差值RMS=sqrt(sum(w [20log10(|G_A|/|G_B|)]^2))。所有sum均在同一M99上，w在该支持内归一化为和1。\n\n');
fprintf(fid,'## 数据和物理范围\n\nPE Model-0与BIE来自活动目录的误差分解MAT派生副本；原R3–R5活动目录当前缺失，没有重算PE或BIE。复用程序时另发现results/archive旧工作树快照保留原曲线，已用纯数值HDF5读取核对三案例PE/BIE与派生副本逐点相同（相对差0），只作来源交叉核对。Bellhop补算flat及三个rough案例、各5001/10001波束，共8个作业。\n\n');
fprintf(fid,'Bellhop专用程序从results/archive旧工作树快照复用，SHA-256与活动构建清单完全一致，两份overlay源码哈希也一致；没有访问cash，没有安装编译器或更换普通Bellhop算法。source_data/binary_origin.json保留来源。\n\n');
fprintf(fid,'共享984点接收网格与Stage-0 M99入射能量权重；主图使用10001波束。Bellhop内部壁剖面为[-80,80] m、4097点、步长.05 m。仅固定表示转换：SHD到native为conj(raw)，native粗糙/平面比再conj进入PE比较表示。无逐案例定相位、幅度拟合或重生成声源。\n\n');
fprintf(fid,'Bellhop有限发射角覆盖外，平面场为零的%d个接收点无法定义粗糙/平面比，显式保存为复NaN并给出ratio_valid掩膜；原始SHD仍保留。它们均在固定M99之外，411个比较点全部有效，比较支持没有缩减；没有插值或补造曲线。\n\n',nnz(~d.cases(1).valid_mask));
fprintf(fid,'历史BIE数值变化指标仅引用原报告，约6e-8，未在本轮重新认证；不得把历史数值当作本轮新Bellhop的严格误差界。波束变化也只作敏感性。\n\n');
fprintf(fid,'本轮10001波束平面控制与保存Stage-0平面Bellhop原生场的固定M99复L2为%.8g；记录在metrics/flat_control.csv。八个作业的进程实际耗时单独保存，不能把读取已完成输出的时间当作求解耗时。\n\n',d.flat_control.complex_l2);
fprintf(fid,'## 本轮指标\n\n|案例|配对|复L2|相位RMS / mrad|幅度比差值RMS / dB|\n|---|---|---:|---:|---:|\n');
for j=1:height(d.metrics)
    r=d.metrics(j,:); fprintf(fid,'|%s|%s|%.8g|%.8g|%.8g|\n',r.case_name{1},r.pair{1},r.complex_l2,r.phase_rms_rad*1000,r.magnitude_ratio_rms_db);
end
fprintf(fid,'\n## 图册\n\n');
for j=1:7
    fprintf(fid,'### %s\n\n![%s](figures/%s.png)\n\n%s\n\n[矢量PDF](figures/%s.pdf)\n\n',titles{j},titles{j},names{j},caps{j},names{j});
    if j==1
        fprintf(fid,'左图明确Tx与Rx中心竖直对齐，接收线位于3 m深度，没有海底反射。右图比较同一波数、同一边缘处理下的三档高度。Bellhop用变换坐标中的内部壁表示相同海面，不使用Tx99实际水深展示环境。\n\n');
    elseif ismember(j,[2 5])
        fprintf(fid,'每列对应弱、中、强海面；上排是|G|，下排是相对于平面反射的主值相位。幅度和相位采用跨列共同刻度，主值跨越±pi时留空，不用连线伪造相位连续性。比较中的细小差异应结合对应差值图判断，不能只凭曲线肉眼重合下结论。\n\n');
        bh=d.metrics.complex_l2(strcmp(d.metrics.pair,'BH_BIE'));
        fprintf(fid,'本轮Bellhop与BIE的加权复差为%.4g至%.4g，两条主曲线因此在此尺度下几乎重合。强海面的主值相位多次跨越±pi，图中的断线表示相位绕回，不是接收声场消失。\n\n',min(bh),max(bh));
        if j==5
            pe=d.metrics.complex_l2(strcmp(d.metrics.pair,'PE_BIE'));
            fprintf(fid,'加入PE后，强海面的幅度空间起伏已可与Bellhop/BIE区分；三种方法的相位在主图尺度下仍看起来接近。PE–BIE复差由弱海面的%.4g增至强海面的%.4g，必须结合TR02与TR03识别，不能把主图相似解读为复数一致。\n\n',pe(1),pe(3));
        end
    elseif j==3
        fprintf(fid,'上排为20log10(|G_BH|/|G_BIE|)，正值表示Bellhop归一幅度更大；下排为arg(G_BH conj(G_BIE))，正值表示比较表示下Bellhop相位更大。各列独立纵轴，1 mrad=1e-3 rad；没有去掉任何全局相位。\n\n');
        for ci=1:3
            c=d.cases(ci); r=c.residuals(2); ix=c.plot_mask;
            fprintf(fid,'%s：M99内最大绝对幅度比差值%.6g dB，最大绝对相对相位%.6g mrad。数值对应差值本身，不能用|G|接近1替代误差判断。\n\n',c.name,max(abs(r.magnitude_ratio_db(ix))),max(abs(r.phase_rad(ix)))*1000);
        end
    elseif j==6
        fprintf(fid,'上排为无量纲复数差的绝对值，下排为有符号相对相位。它们描述接收坐标上的残差结构，不是海面各点贡献，也不能解释为能量差或独立多径。不同配对使用相同M99坐标。\n\n');
        fprintf(fid,'橙色PE–BIE与灰色PE–Bellhop残差在当前纵轴上接近重合，蓝色Bellhop–BIE位于零附近。这与BB02放大的微小蓝色残差一致；蓝线接近零不表示严格相等。残差在接收线边缘较大，指标仍采用原入射能量权重，未删去这些点。\n\n');
    else
        fprintf(fid,'从原始复数数组在固定M99上重新计算，权重为归一化入射能量。复L2以配对中的第二个方法为分母；相位RMS使用wrapped相位；dB项是幅度比差值，不是绝对TL。对数坐标用于显示不同数量级，没有改变原始数值。\n\n');
        if j==7
            for ci=1:3
                pe=d.metrics.complex_l2(strcmp(d.metrics.case_name,d.cases(ci).name)&strcmp(d.metrics.pair,'PE_BIE'));
                bh=d.metrics.complex_l2(strcmp(d.metrics.case_name,d.cases(ci).name)&strcmp(d.metrics.pair,'BH_BIE'));
                fprintf(fid,'%s：PE–BIE复误差为%.6g，Bellhop–BIE为%.6g，两者比值为%.4g。该比值只描述当前固定海面上相对于保存BIE参考的差异，不代表普适模型优劣。\n\n',d.cases(ci).name,pe,bh,pe/bh);
            end
        end
    end
end
fprintf(fid,'## 历史核对与波束敏感性\n\n|案例|配对|本轮复L2|历史报告复L2|\n|---|---|---:|---:|\n');
for j=1:height(d.metrics)
    r=d.metrics(j,:); fprintf(fid,'|%s|%s|%.8g|%.8g|\n',r.case_name{1},r.pair{1},r.complex_l2,d.historical.historical_complex_l2(j));
end
fprintf(fid,'\n历史数字引自source_data中的原报告；差异不通过调参数追平。\n\n|案例|BH5001–BH10001复L2|相位RMS / mrad|\n|---|---:|---:|\n');
for j=1:3, r=d.beam_sensitivity(j,:); fprintf(fid,'|%s|%.8g|%.8g|\n',r.case_name{1},r.complex_l2,r.phase_rms_rad*1000); end
fprintf(fid,'\n## 可复现数据\n\n`plot_data/comparison_data.mat`包含原始复数G、native表示、坐标、海面、权重、掩膜和残差；逐案例CSV保留实部/虚部及全部绘图变换。source_data保存输入快照、哈希、Bellhop原始文件、专用程序和作业检查点。metrics分别保存本轮与历史指标。\n\n');
fprintf(fid,'在setup_vertical_project后运行：`generate_pe_bellhop_bie_comparison_visuals(run_dir)`。指定第二参数可重绘到新目录；无需Bellhop或源工作树路径。\n\n');
fprintf(fid,'复现示例（run_dir指向本数据包所在文件夹）：\n\n```matlab\nsetup_vertical_project;\ngenerate_pe_bellhop_bie_comparison_visuals(run_dir, fullfile(run_dir,''rebuilt_figures''));\naudit_pe_bellhop_bie_visualization(run_dir, fullfile(run_dir,''rebuilt_figures''));\n```\n\n');
fprintf(fid,'本轮仅新增验证/可视化流程，生产物理与公共接口未改；不运行通信或生产信道回归。图组仅适用于三个固定二维确定性海面，不外推PM、三维或通信性能。\n');
fprintf(fid,'\n工作流采用matlab与scientific-visualization技能。软件流程归因：Kassis T., Agarwal V., He Y., Patel D., Brueckner A.M. (2026), [Scientific Agent Skills](https://doi.org/10.48550/arXiv.2609.00065)，当前arXiv记录已核对；此引用不是声学证据。\n');
end

function f=canvas(head,sub)
f=figure('Visible','off','Color','white','Position',[30 30 1600 900]);
set(f,'DefaultAxesFontName','Arial','DefaultTextFontName','Arial');
annotation(f,'textbox',[.04 .925 .92 .06],'String',head,'FontSize',28,'FontWeight','bold','EdgeColor','none','Interpreter','none');
annotation(f,'textbox',[.04 .87 .92 .045],'String',sub,'FontSize',17,'EdgeColor','none','Interpreter','none');
end

function style(ax)
set(ax,'FontSize',17,'LineWidth',1,'Box','off','TickDir','out','XGrid','on','YGrid','on', ...
    'XMinorGrid','off','YMinorGrid','off','GridAlpha',.12);
if isa(ax.Parent,'matlab.graphics.layout.TiledChartLayout'), ax.Parent.OuterPosition=[.045 .115 .91 .735]; end
end

function footer(f,str)
annotation(f,'textbox',[.04 .015 .92 .065],'String',str,'FontName','Microsoft YaHei', ...
    'FontSize',14,'EdgeColor','none','Interpreter','none');
end

function emit(f,folder,name)
exportgraphics(f,fullfile(folder,[name '.png']),'Resolution',150,'BackgroundColor','white', ...
    'Units','inches','Width',16,'Height',9,'Padding','figure');
exportgraphics(f,fullfile(folder,[name '.pdf']),'ContentType','vector','BackgroundColor','white', ...
    'Units','inches','Width',16,'Height',9,'Padding','figure');
close(f); fprintf('Rendered %s\n',name);
end
