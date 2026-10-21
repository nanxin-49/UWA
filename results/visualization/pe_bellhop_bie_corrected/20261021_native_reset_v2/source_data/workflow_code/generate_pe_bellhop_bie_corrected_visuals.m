function output=generate_pe_bellhop_bie_corrected_visuals(run_dir,figure_dir)
%GENERATE_PE_BELLHOP_BIE_CORRECTED_VISUALS Replot corrected saved-data atlas.
% The packaged workflow_code directory permits redraw without source runs.
if nargin<1 || isempty(run_dir), error('An explicit saved run directory is required.'); end
if nargin<2, figure_dir=fullfile(run_dir,'figures'); end
workflow=fullfile(run_dir,'source_data','workflow_code');
assert(isfile(fullfile(workflow,'generate_pe_bellhop_bie_comparison_visuals.m')));
assert(isfile(fullfile(workflow,'generate_pe_bellhop_bie_xz_visuals.m')));
old_path=path; cl_path=onCleanup(@()path(old_path)); addpath(workflow,'-begin');
s=load(fullfile(run_dir,'plot_data','comparison_data.mat'),'data'); d=s.data;
t=load(fullfile(run_dir,'plot_data','xz_data.mat'),'data'); xz=t.data;
assert(strcmp(d.schema_version,'corrected-native-v1'));
assert(strcmp(xz.schema_version,'corrected-native-v1'));
assert(isequal(d.x_m(d.mask),xz.x_m));
assert(all(isfinite(real(xz.PE_carrier_restored(:))) & ...
    isfinite(imag(xz.PE_carrier_restored(:)))));
figs={ '01_environment','bellhop_bie/BB01_fields', ...
    'bellhop_bie/BB02_residuals','bellhop_bie/BB03_metrics', ...
    'three_way/TR01_fields','three_way/TR02_residuals', ...
    'three_way/TR03_metrics','XZ01_raw_magnitude','XZ02_raw_phase', ...
    'XZ03_complex_field','XZ04_raw_profiles','R01_height_sign_reset'};
ready=all(cellfun(@(name) isfile(fullfile(figure_dir,[name '.png'])) && ...
    isfile(fullfile(figure_dir,[name '.pdf'])),figs));
if ready
    line_output=struct('metrics_closure',0);
    xz_output=struct('BH_native_ceiling',max(abs(xz.BH_native(:))));
else
    line_output=generate_pe_bellhop_bie_comparison_visuals(run_dir,figure_dir);
    xz_output=generate_pe_bellhop_bie_xz_visuals(run_dir,figure_dir);
end
assert(line_output.metrics_closure<1e-12);
closure=recheck_metrics(d,run_dir);
make_reset_figure(d,figure_dir);
write_report(d,xz,run_dir,figure_dir);
metadata=table(string(figs(:)),zeros(numel(figs),1),zeros(numel(figs),1), ...
    'VariableNames',{'name','width','height'});
for j=1:numel(figs)
    png=fullfile(figure_dir,[figs{j} '.png']);
    pdf=fullfile(figure_dir,[figs{j} '.pdf']);
    assert(isfile(png) && isfile(pdf)); info=imfinfo(png);
    assert(info.Width==2400 && info.Height==1350);
    metadata.width(j)=info.Width; metadata.height(j)=info.Height;
end
writetable(metadata,fullfile(figure_dir,'image_metadata.csv'));
output=struct('run_dir',run_dir,'figure_dir',figure_dir, ...
    'PNG_PDF_pairs',12,'new_solver_runs',0, ...
    'line_metrics_closure',max(line_output.metrics_closure,closure), ...
    'xz_PE_z3_closure',xz.checks.PE_z3_closure, ...
    'BH_native_ceiling',xz_output.BH_native_ceiling);
fid=fopen(fullfile(run_dir,'render_audit.json'),'w','n','UTF-8');
assert(fid>=0); cl=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(output,PrettyPrint=true)); clear cl
clear cl_path
disp(output);
end

function closure=recheck_metrics(d,run_dir)
stored=readtable(fullfile(run_dir,'metrics','corrected_vs_legacy.csv'));
assert(height(stored)==9); n=0; closure=0;
for j=1:3
    c=d.cases(j); mask=c.plot_mask; w=d.weights(mask); w=w/sum(w);
    A={c.G_PE,c.G_BH,c.G_PE}; B={c.G_BIE,c.G_BIE,c.G_BH};
    for p=1:3
        n=n+1; a=A{p}(mask); b=B{p}(mask);
        E=sqrt(sum(w.*abs(a-b).^2)/sum(w.*abs(b).^2));
        phase=sqrt(sum(w.*angle(a.*conj(b)).^2));
        db=sqrt(sum(w.*(20*log10(abs(a)./abs(b))).^2));
        closure=max(closure,max(abs([E-stored.complex_l2(n), ...
            E-d.metrics.complex_l2(n),phase-stored.phase_rms_rad(n), ...
            db-stored.magnitude_ratio_rms_db(n)])));
    end
end
assert(closure<1e-12);
assert(norm(d.cases(3).G_PE(d.mask)-d.cases(3).G_PE_legacy(d.mask))>1);
end

function make_reset_figure(d,folder)
f=figure('Visible','off','Color','white','Position',[20 20 1600 900]);
annotation(f,'textbox',[.04 .93 .92 .05],'String', ...
    'Height-sign reset: what changed on the same receiver line', ...
    'FontSize',23,'FontWeight','bold','EdgeColor','none');
annotation(f,'textbox',[.04 .88 .92 .04],'String', ...
    'BIE native: purple | PE corrected -2k eta: orange | PE legacy +2k eta: gray', ...
    'FontSize',15,'EdgeColor','none');
lay=tiledlayout(f,2,3,'Padding','loose','TileSpacing','loose');
lay.OuterPosition=[.06 .10 .90 .72];
labels={'Weak A=0.01 m','Medium A=0.05 m','Strong A=0.20 m'};
colors=[.43 .28 .64;.82 .32 .05;.35 .35 .35];
for j=1:3
    c=d.cases(j); mask=c.plot_mask; x=d.x_m(mask);
    ax=nexttile(lay,j); hold(ax,'on');
    plot(ax,x,abs(c.G_BIE(mask))-1,'-','Color',colors(1,:), ...
        'LineWidth',2,'DisplayName','BIE native');
    plot(ax,x,abs(c.G_PE(mask))-1,'-.','Color',colors(2,:), ...
        'LineWidth',2,'DisplayName','PE corrected');
    plot(ax,x,abs(c.G_PE_legacy(mask))-1,':','Color',colors(3,:), ...
        'LineWidth',2,'DisplayName','PE legacy');
    yline(ax,0,'Color',[.7 .7 .7],'HandleVisibility','off');
    title(ax,labels{j}); ylabel(ax,'|G|-1 (dimensionless)');
    set(ax,'FontSize',14,'Box','on');
    ax=nexttile(lay,j+3); hold(ax,'on');
    new_phase=angle(c.G_PE(mask).*conj(c.G_BIE(mask)))*1e3;
    old_phase=angle(c.G_PE_legacy(mask).*conj(c.G_BIE_legacy(mask)))*1e3;
    plot(ax,x,new_phase,'-.','Color',colors(2,:),'LineWidth',2, ...
        'DisplayName','corrected PE - native BIE');
    plot(ax,x,old_phase,':','Color',colors(3,:),'LineWidth',2, ...
        'DisplayName','legacy PE - conjugated BIE');
    yline(ax,0,'Color',[.7 .7 .7],'HandleVisibility','off');
    xlabel(ax,'Receiver x (m)'); ylabel(ax,'Signed phase difference (mrad)');
    set(ax,'FontSize',14,'Box','on');
end
annotation(f,'textbox',[.04 .02 .92 .035],'String', ...
    'G = rough/flat. Top shows amplitude deviation from 1; bottom compares distinct representation pairs. No fitting. 1 mrad = 0.001 rad.', ...
    'FontSize',12,'EdgeColor','none','Interpreter','none');
drawnow; png=fullfile(folder,'R01_height_sign_reset.png');
exportgraphics(f,png,'Resolution',150,'BackgroundColor','white', ...
    'Units','inches','Width',16,'Height',9,'Padding','figure'); close(f);
pixels=imread(png); page=figure('Visible','off','Color','white', ...
    'Position',[20 20 1600 900]);
ax=axes(page,'Position',[0 0 1 1]); image(ax,pixels);
set(ax,'Visible','off','Position',[0 0 1 1]);
exportgraphics(page,fullfile(folder,'R01_height_sign_reset.pdf'), ...
    'ContentType','image','BackgroundColor','white'); close(page);
end

function write_report(d,xz,out,figdir)
tab=readtable(fullfile(out,'metrics','corrected_vs_legacy.csv'));
decomp=readtable(fullfile(out,'metrics','complex_error_decomposition.csv'));
fid=fopen(fullfile(out,'corrected_report.md'),'w','n','UTF-8');
assert(fid>=0); cl=onCleanup(@()fclose(fid));
fprintf(fid,'# PE–Bellhop–BIE 高度符号修正后的可视化\n\n');
fprintf(fid,'本包与原图册分开保存。固定4 kHz、c=1500 m/s、Tx深100 m、Rx深3 m、Gaussian σ=0.3 m、三个正弦海面A=0.01/0.05/0.20 m，K=0.10 rad/m。海面在|x|≤42 m保持原形，50 m以外归零。\n\n');
fprintf(fid,'海面几何采用z=η（深度向下），时间约定exp(−iωt)。PE相位屏由旧的exp(+2ikη)改为exp(−2ikη)，并从同一保存入射场FFT传播到接收线与真实x–z网格。BIE、Bellhop沿用保存原生复场；三方G都用原生表示。没有重新传播声源、求解BIE或运行Bellhop，也没有改变生产API。\n\n');
fprintf(fid,'## 先看直观图\n\n');
figs={'R01_height_sign_reset','three_way/TR01_fields','three_way/TR02_residuals', ...
    'three_way/TR03_metrics','XZ01_raw_magnitude','XZ02_raw_phase', ...
    'XZ03_complex_field','XZ04_raw_profiles'};
captions={'同一接收线的旧PE、修正PE与原生BIE；上排是|G|−1，下排是两种各自配对的相位差。', ...
    '三方法直接比较G的幅值和主值相位；修正PE幅值调制方向与BIE一致。', ...
    '三组配对的点位复差和有符号相位差；展示差异出现的位置，纵轴各案例独立。', ...
    '固定M99权重下的复L2、相位RMS及dB幅值比RMS；曲线趋势限于这三个海面。', ...
    '强海面真实x–z原始反射场幅值；PE/BIE同色标，Bellhop独立源标度。', ...
    '同一x–z复场的主值相位；色环首尾相接，色差不能单独当作误差。', ...
    '色相表示相位、明暗表示幅值；暗区相位虽存储，但不宜凭颜色解读。', ...
    'z=1、3、6 m原始幅值剖面；上排PE/BIE共标度，下排Bellhop原生单位。'};
for j=1:numel(figs)
    fprintf(fid,'### %s\n\n![%s](figures/%s.png)\n\n%s\n\n', ...
        figs{j},captions{j},figs{j},captions{j});
end
fprintf(fid,'Bellhop–BIE子图BB01–BB03和共同环境图01仍在figures/目录；该配对同时转为原生表示，其复差数值不变。彩图画的是原始u(x,z)，其Bellhop源标度与Gaussian峰值尚未标定，不能把三列像素亮度直接当作跨方法误差。接收线G=粗糙反射场/平面反射场，是无量纲比值，不是绝对声压或传播损失。\n\n');
fprintf(fid,'## 固定接收线指标（本轮与旧值分列）\n\n');
fprintf(fid,'|海面|配对|本轮复L2|本轮幅值相对L2|本轮相位RMS (mrad)|本轮dB幅值比RMS|旧复L2|\n|---|---|---:|---:|---:|---:|---:|\n');
for j=1:height(tab)
    fprintf(fid,'|%s|%s|%.8g|%.8g|%.8g|%.8g|%.8g|\n', ...
        char(string(tab.case_name(j))),char(string(tab.pair(j))),tab.complex_l2(j), ...
        tab.magnitude_relative_l2(j),tab.phase_rms_rad(j)*1e3, ...
        tab.magnitude_ratio_rms_db(j),tab.legacy_complex_l2(j));
end
fprintf(fid,'\n复差精确拆分：E² = 加权(|A|−|B|)²/加权|B|² + 加权[2|A||B|(1−cos Δφ)]/加权|B|²。强海面PE–BIE的两项为%.8g与%.8g，相位项占%.5f%%；旧8.0908%%包含先前高度/表示失配，不能继续当作修正后的模型误差。\n\n', ...
    decomp.magnitude_term(7),decomp.phase_term(7),decomp.phase_fraction(7)*100);
fprintf(fid,'修正后强海面PE–BIE：复L2 %.5f%%，幅值相对L2 %.5f%%，相位RMS %.5f mrad。仍有残余相位差；本次不把它归因给单一机制，也不宣称模型普遍精度。Bellhop–BIE复L2 %.6f%%。\n\n', ...
    100*tab.complex_l2(7),100*tab.magnitude_relative_l2(7), ...
    1000*tab.phase_rms_rad(7),100*tab.complex_l2(8));
fprintf(fid,'## 检查与复现\n\n');
fprintf(fid,'平面PE分母复现相对L2 %.3g；强海面二维z=3 m截线闭合 %.3g。每种方法的x–z实数/虚数保存在plot_data/xz_data.mat和xz_fields.csv；三案例G、掩膜、权重和旧PE保存在comparison_data.mat与逐案例CSV。source_data/source_hashes.csv记录三个原始MAT的路径及SHA256；旧图册、原BIE/Bellhop输出保持不变。\n\n', ...
    d.correction.flat_relative_l2,xz.checks.PE_z3_closure);
fprintf(fid,'在MATLAB执行 `setup_vertical_project; generate_pe_bellhop_bie_corrected_visuals(''本包绝对目录'',''新图目录'');`，只读取本包MAT和随包保存的绘图源码；不运行声学求解。生成12对2400×1350 PNG/PDF。\n\n');
fprintf(fid,'适用边界：固定单频、共同M99接收区、三档确定性海面；BIE为历史派生副本，没有本轮收敛复核；5001/10001波束差是敏感性而非误差上界；二维彩图仅强海面、近海面z=0.5–6 m。\n');
clear cl
% Replace reports emitted by the legacy figure functions with one pointer.
for name={'comparison_report.md','xz_field_report.md'}
    f=fopen(fullfile(out,name{1}),'w','n','UTF-8'); assert(f>=0);
    fprintf(f,'# 本轮已重置高度符号\n\n请阅读 [修正后的统一报告](corrected_report.md)。\n'); fclose(f);
end
assert(isfolder(figdir));
end
