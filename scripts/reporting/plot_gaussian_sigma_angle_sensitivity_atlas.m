function analysis = plot_gaussian_sigma_angle_sensitivity_atlas
%PLOT_GAUSSIAN_SIGMA_ANGLE_SENSITIVITY_ATLAS
% Replot saved Gaussian-sigma PE/FK/BIE data; never runs a field solver.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
base_dir = fullfile(root,'results','validation','pe_bie_full_kirchhoff_gaussian_sigma_sensitivity');
kz_dir = fullfile(root,'results','validation','pe_bie_full_kirchhoff_gaussian_sigma_kzaware');
out_dir = fullfile(root,'results','visualization','pe_bie_gaussian_sigma_angle_sensitivity');
fig_dir = fullfile(out_dir,'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
base_path = fullfile(base_dir,'gaussian_sigma_sensitivity.mat');
kz_path = fullfile(kz_dir,'gaussian_sigma_kzaware_comparison.mat');
assert(exist(base_path,'file')==2 && exist(kz_path,'file')==2, ...
    'Required saved sensitivity MAT files were not found. No computations were run.');
base_loaded=load(base_path,'result'); base=base_loaded.result;
kz_loaded=load(kz_path,'result'); kz=kz_loaded.result;
analysis=local_build_analysis(base,kz,base_path,kz_path);
local_write_outputs(analysis,root,out_dir,fig_dir);
fprintf('Saved-data plotting completed: %s\n',out_dir);
disp(analysis.summary(:,{'scene','sigma_m','theta99_deg','n_common', ...
    'Model0_Ec_common','Model1_Ec_common','FullKirchhoff_Ec_common','improvement_common'}));
end

function a=local_build_analysis(base,kz,base_path,kz_path)
assert(isequal(base.config.sigmas,kz.config.sigmas),'Sigma nodes differ between source MAT files.');
assert(numel(base.cases)==15 && numel(kz.cases)==15,'Expected complete 3-scene x 5-sigma MAT data.');
scene_names={base.config.cases.name}; sigma_nodes=base.config.sigmas(:).';
summary_rows=repmat(local_summary_row(),15,1); cases=repmat(local_case_row(),15,1); index=0;
theta_by_sigma=nan(3,numel(sigma_nodes),3);
for ic=1:numel(scene_names)
    for is=1:numel(sigma_nodes)
        index=index+1; scene=scene_names{ic}; sig=sigma_nodes(is);
        ib=find(strcmp({base.cases.case_name},scene) & abs([base.cases.sigma_m]-sig)<1e-12);
        ik=find(strcmp({kz.cases.case_name},scene) & abs([kz.cases.sigma_m]-sig)<1e-12);
        assert(isscalar(ib)&&isscalar(ik),'Missing/ambiguous saved fields for %s sigma %.6g.',scene,sig);
        d0=base.cases(ib); d1=kz.cases(ik);
        assert(isequal(d0.x_m(:),d1.x_m(:)),'Receiver grids differ for %s sigma %.6g.',scene,sig);
        assert(isequal(d0.footprint.m99.mask(:),d1.footprint.m99.mask(:)), ...
            'Saved M99 footprint differs between MAT sources for %s sigma %.6g.',scene,sig);
        assert(max(abs(d0.receiver_PE(:)-d1.receiver_PE(:)))<1e-10*max(1,max(abs(d0.receiver_PE(:)))), ...
            'Model-0 saved fields disagree for %s sigma %.6g.',scene,sig);
        fields=struct('Model0',d1.receiver_PE(:),'Model1',d1.receiver_Model1(:), ...
            'FullKirchhoff',d1.receiver_FK(:),'BIE',d1.receiver_BIE(:));
        base_mask=d1.footprint.m99.mask(:);
        finite_mask=isfinite(fields.Model0)&isfinite(fields.Model1)& ...
            isfinite(fields.FullKirchhoff)&isfinite(fields.BIE);
        % Apply the -40 dB gate within the saved M99 footprint; receiver-tail
        % outliers outside M99 must not determine the formal threshold.
        common_peak=max([abs(fields.Model0(base_mask));abs(fields.Model1(base_mask)); ...
            abs(fields.FullKirchhoff(base_mask));abs(fields.BIE(base_mask))]);
        threshold=10^(-40/20)*common_peak;
        common_mask=base_mask&finite_mask&abs(fields.Model0)>=threshold& ...
            abs(fields.Model1)>=threshold&abs(fields.FullKirchhoff)>=threshold&abs(fields.BIE)>=threshold;
        assert(any(common_mask),'Common four-field mask is empty for %s sigma %.6g.',scene,sig);
        w=d1.footprint.energy_weights(:); w=w(common_mask); w=w/max(sum(w),realmin);
        metrics=struct(); names={'Model0','Model1','FullKirchhoff'};
        for jj=1:numel(names), metrics.(names{jj})=local_metrics(fields.(names{jj}),fields.BIE,common_mask,w); end
        % Preserve the original pairwise-mask values separately for auditing.
        old=struct('Model0_Ec',d1.metrics_Model0_BIE.complex_l2, ...
            'Model1_Ec',d1.metrics_Model1_BIE.complex_l2, ...
            'FullKirchhoff_Ec',d1.metrics_FK_BIE.complex_l2, ...
            'Model0_improvement',d1.metrics_Model0_BIE.complex_l2/ ...
            max(d1.metrics_Model1_BIE.complex_l2,realmin));
        theta=d1.spectrum; theta_by_sigma(ic,is,:)=[theta.theta90_deg theta.theta95_deg theta.theta99_deg];
        s=local_summary_row(); s.scene=scene; s.A_m=d1.A_m; s.K_radpm=d1.K_radpm; s.sigma_m=sig;
        s.theta90_deg=theta.theta90_deg; s.theta95_deg=theta.theta95_deg; s.theta99_deg=theta.theta99_deg;
        s.n_common=sum(common_mask); s.mask_x_min_m=min(d1.x_m(common_mask)); s.mask_x_max_m=max(d1.x_m(common_mask));
        s.Model0_Ec_common=metrics.Model0.complex_l2; s.Model1_Ec_common=metrics.Model1.complex_l2;
        s.FullKirchhoff_Ec_common=metrics.FullKirchhoff.complex_l2;
        s.Model0_mag_common=metrics.Model0.magnitude_l2; s.Model1_mag_common=metrics.Model1.magnitude_l2;
        s.FullKirchhoff_mag_common=metrics.FullKirchhoff.magnitude_l2;
        s.Model0_phase_common=metrics.Model0.phase_rms; s.Model1_phase_common=metrics.Model1.phase_rms;
        s.FullKirchhoff_phase_common=metrics.FullKirchhoff.phase_rms;
        s.improvement_common=metrics.Model0.complex_l2/max(metrics.Model1.complex_l2,realmin);
        s.old_Model0_Ec_pairwise=old.Model0_Ec; s.old_Model1_Ec_pairwise=old.Model1_Ec;
        s.old_FullKirchhoff_Ec_pairwise=old.FullKirchhoff_Ec; s.old_improvement_pairwise=old.Model0_improvement;
        summary_rows(index)=s;
        cr=local_case_row(); cr.scene=scene; cr.sigma_m=sig; cr.x_m=d1.x_m(:); cr.fields=fields;
        cr.base_M99=base_mask; cr.common_mask=common_mask; cr.energy_weights=d1.footprint.energy_weights(:);
        cr.metrics=metrics; cr.spectrum=d1.spectrum; cr.old_pairwise=old; cases(index)=cr;
    end
end
for jj=1:3, assert(max(max(abs(theta_by_sigma(:,:,jj)-theta_by_sigma(1,:,jj)),[],2))<1e-10, ...
        'Angular-spectrum quantiles unexpectedly differ by surface case.'); end
a=struct('schema_version','1.0.0','source_mat_files',{{base_path,kz_path}}, ...
    'sigma_nodes_m',sigma_nodes,'scene_names',{scene_names},'summary',struct2table(summary_rows), ...
    'cases',cases,'mask_definition',['New common mask: saved sigma-specific M99 footprint intersected ', ...
    'with finite values and all four fields (Model-0, Model-1, Full Kirchhoff, BIE) above -40 dB ', ...
    'relative to the maximum of those four fields; saved M99 energy weights renormalized on the intersection.'], ...
    'old_mask_definition',base.mask_definition,'comparison_convention',base.comparison_convention, ...
    'model1_comparison_convention',kz.comparison_convention,'theta_quantiles_deg',squeeze(theta_by_sigma(1,:,:)));
end

function m=local_metrics(u,ref,mask,w)
u=u(mask); ref=ref(mask); phase=angle(u.*conj(ref)); au=abs(u); ar=abs(ref);
m=struct('complex_l2',sqrt(sum(w.*abs(u-ref).^2)/max(sum(w.*ar.^2),realmin)), ...
    'magnitude_l2',sqrt(sum(w.*(au-ar).^2)/max(sum(w.*ar.^2),realmin)), ...
    'phase_rms',sqrt(sum(w.*phase.^2)));
end

function r=local_summary_row()
r=struct('scene','','A_m',NaN,'K_radpm',NaN,'sigma_m',NaN,'theta90_deg',NaN, ...
    'theta95_deg',NaN,'theta99_deg',NaN,'n_common',NaN,'mask_x_min_m',NaN,'mask_x_max_m',NaN, ...
    'Model0_Ec_common',NaN,'Model1_Ec_common',NaN,'FullKirchhoff_Ec_common',NaN, ...
    'Model0_mag_common',NaN,'Model1_mag_common',NaN,'FullKirchhoff_mag_common',NaN, ...
    'Model0_phase_common',NaN,'Model1_phase_common',NaN,'FullKirchhoff_phase_common',NaN, ...
    'improvement_common',NaN,'old_Model0_Ec_pairwise',NaN,'old_Model1_Ec_pairwise',NaN, ...
    'old_FullKirchhoff_Ec_pairwise',NaN,'old_improvement_pairwise',NaN);
end
function r=local_case_row()
r=struct('scene','','sigma_m',NaN,'x_m',[],'fields',[],'base_M99',[],'common_mask',[], ...
    'energy_weights',[],'metrics',[],'spectrum',[],'old_pairwise',[]);
end

function local_write_outputs(a,root,out_dir,fig_dir)
T=a.summary; csv=fullfile(out_dir,'gaussian_sigma_common_mask_summary.csv'); mat=fullfile(out_dir,'gaussian_sigma_angle_sensitivity_visualization.mat');
writetable(T,csv); save(mat,'a','-v7.3');
local_plot_angle(a,fig_dir); local_plot_ec(a,fig_dir); local_plot_improvement(a,fig_dir);
local_plot_decomposition(a,fig_dir,'phase'); local_plot_decomposition(a,fig_dir,'magnitude');
local_plot_heatmap(a,fig_dir); local_plot_representatives(a,fig_dir);
 local_write_report(a,root,out_dir,csv,mat);
end

function local_plot_angle(a,fd)
sig=a.sigma_nodes_m; q=a.theta_quantiles_deg; f=local_figure([8.2 5.3]); hold on;
plot(sig,q(:,1),'o-','Color',[0.00 0.447 0.741],'LineWidth',1.8,'MarkerFaceColor',[0.00 0.447 0.741]);
plot(sig,q(:,2),'s--','Color',[0.85 0.325 0.098],'LineWidth',1.8,'MarkerFaceColor',[0.85 0.325 0.098]);
plot(sig,q(:,3),'d-.','Color',[0.466 0.674 0.188],'LineWidth',1.8,'MarkerFaceColor',[0.466 0.674 0.188]);
grid on; xlabel('$\sigma$ (m)','Interpreter','latex'); ylabel('Incidence angle (deg)','Interpreter','latex');
title('Gaussian source angular-spectrum quantiles','Interpreter','none'); legend('$\theta_{90}$','$\theta_{95}$','$\theta_{99}$','Interpreter','latex','Location','northeast');
local_save_figure(f,fd,'01_incident_angle_quantiles');
end

function local_plot_ec(a,fd)
T=a.summary; scenes=a.scene_names; colors=local_model_colors(); f=local_figure([12.2 4.5]); tiledlayout(1,3,'TileSpacing','compact');
vals=[T.Model0_Ec_common;T.Model1_Ec_common;T.FullKirchhoff_Ec_common]; ymin=10^floor(log10(min(vals)*0.75)); ymax=10^ceil(log10(max(vals)*1.4));
for ic=1:3
    nexttile; q=T(strcmp(T.scene,scenes{ic}),:); hold on;
    loglog(q.sigma_m,q.Model0_Ec_common,'o-','Color',colors(1,:),'LineWidth',1.8,'MarkerFaceColor',colors(1,:));
    loglog(q.sigma_m,q.Model1_Ec_common,'s--','Color',colors(2,:),'LineWidth',1.8,'MarkerFaceColor',colors(2,:));
    loglog(q.sigma_m,q.FullKirchhoff_Ec_common,'d-.','Color',colors(3,:),'LineWidth',1.8,'MarkerFaceColor',colors(3,:));
    set(gca,'XScale','log','YScale','log','YLim',[ymin ymax]); grid on; title(local_scene_title(scenes{ic}),'Interpreter','none');
    xlabel('$\sigma$ (m)','Interpreter','latex'); if ic==1,ylabel('$E_c$ vs. BIE','Interpreter','latex');end
    if ic==1,legend('Model-0 PE','Model-1 $k_z$-aware','Full Kirchhoff','Interpreter','latex','Location','southwest');end
end
sgtitle('Common-mask normalized complex error','Interpreter','none'); local_save_figure(f,fd,'02_complex_error_by_sigma');
end

function local_plot_improvement(a,fd)
T=a.summary; scenes=a.scene_names; colors=local_model_colors(); f=local_figure([8.2 5.3]); hold on;
for ic=1:3
    q=T(strcmp(T.scene,scenes{ic}),:); q=sortrows(q,'theta99_deg');
    plot(q.theta99_deg,q.improvement_common,'o-','Color',colors(ic,:),'LineWidth',1.8,'MarkerFaceColor',colors(ic,:), ...
        'DisplayName',local_scene_title(scenes{ic}));
end
yline(1,'k--','Improvement=1','LabelHorizontalAlignment','left'); grid on;
xlabel('$\theta_{99}$ (deg)','Interpreter','latex'); ylabel('$E_c(\mathrm{Model\!-!0},\mathrm{BIE})/E_c(\mathrm{Model\!-!1},\mathrm{BIE})$','Interpreter','latex');
title('$k_z$-aware improvement vs. incident angular width','Interpreter','latex'); legend('Location','best','Interpreter','none');
local_save_figure(f,fd,'03_kz_aware_improvement_vs_theta99');
end

function local_plot_decomposition(a,fd,kind)
T=a.summary; scenes=a.scene_names; colors=local_model_colors(); f=local_figure([12.2 4.5]); tiledlayout(1,3,'TileSpacing','compact');
if strcmp(kind,'phase'), fields={'Model0_phase_common','Model1_phase_common','FullKirchhoff_phase_common'}; yl='Phase RMS (rad)'; stem='04_phase_error_by_sigma'; ttl='Phase error relative to BIE';
else, fields={'Model0_mag_common','Model1_mag_common','FullKirchhoff_mag_common'}; yl='Magnitude relative L2'; stem='05_magnitude_error_by_sigma'; ttl='Magnitude error relative to BIE'; end
vv=[]; for ic=1:3,q=T(strcmp(T.scene,scenes{ic}),:);vv=[vv;q{:,fields}(:)];end %#ok<AGROW>
ymin=10^floor(log10(min(vv(vv>0))*0.75)); ymax=10^ceil(log10(max(vv)*1.4));
for ic=1:3
    nexttile; q=T(strcmp(T.scene,scenes{ic}),:); hold on;
    loglog(q.sigma_m,q.(fields{1}),'o-','Color',colors(1,:),'LineWidth',1.8,'MarkerFaceColor',colors(1,:));
    loglog(q.sigma_m,q.(fields{2}),'s--','Color',colors(2,:),'LineWidth',1.8,'MarkerFaceColor',colors(2,:));
    loglog(q.sigma_m,q.(fields{3}),'d-.','Color',colors(3,:),'LineWidth',1.8,'MarkerFaceColor',colors(3,:));
    set(gca,'XScale','log','YScale','log','YLim',[ymin ymax]); grid on; title(local_scene_title(scenes{ic}),'Interpreter','none');
    xlabel('$\sigma$ (m)','Interpreter','latex'); if ic==1,ylabel(yl,'Interpreter','latex');end
    if ic==1,legend('Model-0 PE','Model-1 $k_z$-aware','Full Kirchhoff','Interpreter','latex','Location','southwest');end
end
sgtitle(ttl,'Interpreter','none'); local_save_figure(f,fd,stem);
end

function local_plot_heatmap(a,fd)
T=a.summary; names=cell(height(T),1);
for i=1:height(T)
    switch T.scene{i}
        case 'weak_low_K', short='Weak/low-K';
        case 'strong_height_low_K', short='Strong-height/low-K';
        case 'weak_high_K', short='Weak/high-K';
        otherwise, short=T.scene{i};
    end
    names{i}=sprintf('%s | $\\sigma=%.1f$ m',short,T.sigma_m(i));
end
E=[T.Model0_Ec_common T.Model1_Ec_common T.FullKirchhoff_Ec_common]; I=T.improvement_common;
f=local_figure([14.0 8.0]); tiledlayout(1,2,'TileSpacing','compact');
nexttile; imagesc(log10(E)); set(gca,'YDir','reverse','XTick',1:3,'XTickLabel',{'Model-0','Model-1','Full K'}, ...
    'YTick',1:height(T),'YTickLabel',names,'TickLabelInterpreter','latex'); colorbar; clim([floor(log10(min(E(:)))) ceil(log10(max(E(:))))]);
title('$\log_{10}(E_c)$, common mask','Interpreter','latex');
for i=1:size(E,1),for j=1:size(E,2),text(j,i,sprintf('%.2g',E(i,j)),'HorizontalAlignment','center','FontSize',8,'Color','k');end,end
nexttile; imagesc(I); set(gca,'YDir','reverse','XTick',1,'XTickLabel',{'M0/M1'}, ...
    'YTick',1:height(T),'YTickLabel',names,'TickLabelInterpreter','latex'); colorbar; clim([min(I) max(I)]);
title('Improvement ratio, common mask','Interpreter','none');
for i=1:numel(I),text(1,i,sprintf('%.2f',I(i)),'HorizontalAlignment','center','FontSize',8,'Color','k');end
sgtitle('All 15 scene–sigma comparisons','Interpreter','none'); local_save_figure(f,fd,'06_all_case_metric_heatmap');
end

function local_plot_representatives(a,fd)
spec={{'weak_low_K',0.3},{'strong_height_low_K',1},{'weak_high_K',5}}; colors=local_model_colors(); biecolor=[0.18 0.18 0.18];
for k=1:numel(spec)
    scene=spec{k}{1}; sig=spec{k}{2}; idx=find(strcmp({a.cases.scene},scene)&abs([a.cases.sigma_m]-sig)<1e-12); assert(isscalar(idx)); d=a.cases(idx);
    flds=d.fields; x=d.x_m; mask=d.common_mask; xm=x(mask); f=local_figure([11.2 8.0]); tiledlayout(2,2,'TileSpacing','compact');
    nexttile; hold on; plot(xm,abs(flds.Model0(mask)),'-','Color',colors(1,:),'LineWidth',1.4);
    plot(xm,abs(flds.Model1(mask)),'--','Color',colors(2,:),'LineWidth',1.4);
    plot(xm,abs(flds.FullKirchhoff(mask)),'-.','Color',colors(3,:),'LineWidth',1.4); plot(xm,abs(flds.BIE(mask)),':','Color',biecolor,'LineWidth',1.6);
    grid on; ylabel('Normalized magnitude','Interpreter','latex'); title('Receiver-line amplitude','Interpreter','none');
    legend('Model-0 PE','Model-1 $k_z$-aware','Full Kirchhoff','Helmholtz BIE','Interpreter','latex','Location','best');
    nexttile; hold on; plot(xm,angle(flds.Model0(mask)),'-','Color',colors(1,:),'LineWidth',1.4);
    plot(xm,angle(flds.Model1(mask)),'--','Color',colors(2,:),'LineWidth',1.4);
    plot(xm,angle(flds.FullKirchhoff(mask)),'-.','Color',colors(3,:),'LineWidth',1.4); plot(xm,angle(flds.BIE(mask)),':','Color',biecolor,'LineWidth',1.6);
    grid on; ylabel('Wrapped phase (rad)','Interpreter','latex'); title('Receiver-line wrapped phase','Interpreter','none');
    nexttile; hold on; plot(xm,20*log10(max(abs(flds.Model0(mask)),realmin)./max(abs(flds.BIE(mask)),realmin)),'-','Color',colors(1,:),'LineWidth',1.4);
    plot(xm,20*log10(max(abs(flds.Model1(mask)),realmin)./max(abs(flds.BIE(mask)),realmin)),'--','Color',colors(2,:),'LineWidth',1.4);
    grid on; ylabel('Magnitude difference vs BIE (dB)','Interpreter','latex'); title('Pointwise magnitude difference','Interpreter','none');
    legend('Model-0/BIE','Model-1/BIE','Location','best');
    nexttile; hold on; plot(xm,angle(flds.Model0(mask).*conj(flds.BIE(mask))),'-','Color',colors(1,:),'LineWidth',1.4);
    plot(xm,angle(flds.Model1(mask).*conj(flds.BIE(mask))),'--','Color',colors(2,:),'LineWidth',1.4);
    grid on; ylabel('Wrapped phase difference (rad)','Interpreter','latex'); title('Pointwise phase difference','Interpreter','none'); legend('Model-0/BIE','Model-1/BIE','Location','best');
    for p=1:4,nexttile(p);xlabel('Receiver transverse coordinate $x$ (m)','Interpreter','latex');end
    sgtitle(sprintf('%s, sigma=%.1f m; common mask n=%d, x=[%.2f, %.2f] m', ...
        scene,sig,sum(mask),min(xm),max(xm)),'Interpreter','none');
    stem=sprintf('07_receiver_%s_sigma_%g',scene,sig); local_save_figure(f,fd,stem);
end
end

function local_write_report(a,root,out_dir,csv,mat)
p=fullfile(root,'reports','pe_bie_gaussian_sigma_angle_sensitivity_visualization_report.md'); fid=fopen(p,'w','n','UTF-8'); assert(fid>=0); cl=onCleanup(@()fclose(fid)); %#ok<NASGU>
T=a.summary; fprintf(fid,'# Gaussian 入射角谱与 PE/FK/BIE 误差可视化报告\n\n');
fprintf(fid,'状态：**只做已有数据核查、指标重算与绘图；没有执行 PE、Full Kirchhoff 或 BIE 求解，也没有修改任何物理模型。**\n\n');
fprintf(fid,'## 1. 实验设置和数据来源\n\n');
fprintf(fid,'读取 `results/validation/pe_bie_full_kirchhoff_gaussian_sigma_sensitivity/gaussian_sigma_sensitivity.mat` 与 `results/validation/pe_bie_full_kirchhoff_gaussian_sigma_kzaware/gaussian_sigma_kzaware_comparison.mat`。预期 3 个海面场景 × 5 个 sigma 共 15 组；MAT 中已找到全部 15 组四种场（Model-0、Model-1、Full Kirchhoff、BIE）的 receiver 复场。\n\n');
fprintf(fid,'环境为 2D/one-transverse，4 kHz，均匀 `c=1500 m/s`，Tx=100 m、Rx=3 m、压力释放海面，sigma=`0.3, 0.5, 1, 2, 5 m`。场景：weak_low_K (`A=0.05 m, K=0.10 rad/m`)、strong_height_low_K (`A=0.20 m, K=0.10 rad/m`)、weak_high_K (`A=0.02 m, K=0.47 rad/m`)。\n\n');
fprintf(fid,'这些 MAT 中的接收复场是无量纲 rough/flat reflected-response ratio，并按既有比较约定处理 FK/BIE 的共轭；不是未归一化的绝对声压。源数据与原仿真结果均保持未改。\n\n');
fprintf(fid,'### 评价 mask 与指标定义\n\n已有 sensitivity / kz-aware 报告使用 sigma-specific M99 与逐模型对 BIE 的 pairwise -40 dB 门槛；不同模型 Ec 可能来自不同的接收点集合，不应横向视为完全同采样区域。本报告另外从已存复场重算统一指标：每个 sigma 使用其已存 M99 mask，取所有四个复场均有限且各自幅度不低于四场共同峰值的 -40 dB 的交集，并在交集上重新归一化原 M99 energy weights。mask 随 sigma 仍可能变化，故跨 sigma 是各自有效 footprint 上的一致口径，不是固定空间区域的严格比较。\n\n');
fprintf(fid,'统一复场误差 `E_c=sqrt(sum(w|u-u_BIE|^2)/sum(w|u_BIE|^2))`；幅度误差 `E_a=sqrt(sum(w(|u|-|u_BIE|)^2)/sum(w|u_BIE|^2))`；phase RMS 为 `sqrt(sum(w angle(u conj(u_BIE))^2))`。完整的新旧 15 组指标见 CSV。\n\n');
fprintf(fid,'图件均导出为矢量 PDF 和 300 dpi PNG，存于 `%s`。\n\n',local_rel(root,fullfile(out_dir,'figures')));
fprintf(fid,'## 2. 入射角谱随 sigma 的变化\n\n'); local_fig_caption(fid,'01_incident_angle_quantiles','图1：同一 Gaussian 源的传播角谱 theta90、theta95、theta99 累积能量角度分位数随 sigma 增大而下降；三个海面场景使用相同源谱，因此角度分位数相同。');
fprintf(fid,'数值从宽角端至窄角端，theta99 为 `%.4g` 到 `%.4g deg`；theta90 为 `%.4g` 到 `%.4g deg`。角度统计只计传播谱分量，曲线没有插值或平滑。\n\n',T.theta99_deg(1),T.theta99_deg(5),T.theta90_deg(1),T.theta90_deg(5));
fprintf(fid,'## 3. Model-0、Model-1、Full Kirchhoff 总体精度\n\n'); local_fig_caption(fid,'02_complex_error_by_sigma','图2（核心）：三联图比较三个场景的 Model-0、kz-aware Model-1 和 Full Kirchhoff 相对 BIE 的统一 mask 复场误差；共用 log y 轴和一致的模型色/线型。');
fprintf(fid,'在 weak_low_K、sigma=0.3 m，统一 mask 下 Model-1 相对 Model-0 的 Ec 改善为 `%.3g` 倍；三场景、sigma=5 m 时 Ec(M0)/Ec(M1) 分别为 `%.3g / %.3g / %.3g`。Full Kirchhoff 的同口径误差范围为 `%.3g` 至 `%.3g`。\n\n', ...
    local_ratio(T,'weak_low_K',0.3),local_ratio(T,'weak_low_K',5),local_ratio(T,'strong_height_low_K',5),local_ratio(T,'weak_high_K',5),min(T.FullKirchhoff_Ec_common),max(T.FullKirchhoff_Ec_common));
fprintf(fid,'## 4. kz-aware 改善倍数\n\n'); local_fig_caption(fid,'03_kz_aware_improvement_vs_theta99','图3：改善倍数定义为同一 common mask 上 Ec(Model-0,BIE)/Ec(Model-1,BIE)；水平线 1 表示没有改善，低于 1 表示 Model-1 略差。');
fprintf(fid,'在宽角 sigma=0.3 m，三个 case 的改善倍数依次为 `%.3g, %.3g, %.3g`；在 sigma=5 m 分别为 `%.3g, %.3g, %.3g`。宽角下收益具有场景依赖性；缩窄角谱后多数比值趋近 1，不应据此推断误差全都消失。\n\n', ...
    local_ratio(T,'weak_low_K',.3),local_ratio(T,'strong_height_low_K',.3),local_ratio(T,'weak_high_K',.3), ...
    local_ratio(T,'weak_low_K',5),local_ratio(T,'strong_height_low_K',5),local_ratio(T,'weak_high_K',5));
fprintf(fid,'## 5. 相位与幅度误差分解\n\n'); local_fig_caption(fid,'04_phase_error_by_sigma','图4：三联图展示两种 PE 与 Full Kirchhoff 的加权相位 RMS；对数纵轴突出相位误差数量级。');
local_fig_caption(fid,'05_magnitude_error_by_sigma','图5：三联图展示幅度相对 L2；重点可检查 weak_high_K 随 sigma 增大后的幅度误差，而不把相位改善误认为幅度也改善。');
fprintf(fid,'weak_low_K 中随角谱变窄，Model-0/Model-1 phase RMS 整体显著下降，支持有限角度项的存在。strong_height_low_K 仍保留较大残差。weak_high_K 在窄角端仍有可见幅度误差；因此其剩余差异不能仅用入射角宽度解释。该组数据本身不能唯一判定残差是局部边界近似、有限窗/离散误差或其他模型因素。\n\n');
fprintf(fid,'## 6. 代表性接收场分析\n\n');
reps={{'weak_low_K',.3},{'strong_height_low_K',1},{'weak_high_K',5}};
for i=1:numel(reps)
    scene=reps{i}{1}; sig=reps{i}{2}; stem=sprintf('07_receiver_%s_sigma_%g',scene,sig);
    local_fig_caption(fid,stem,sprintf('代表性接收线图：幅度和 wrapped phase 对比四种方法，下排给出 Model-0/BIE、Model-1/BIE 的逐点幅度(dB)与 wrapped phase 差；只画 common M99/-40 dB 有效点，坐标 gaps 表示被 mask 的样本。(%s, sigma=%.1f m)',scene,sig));
end
fprintf(fid,'相位差采用 `angle(u_model conj(u_BIE))`，限制在 `[-pi,pi]`；没有对曲线做 unwrap，以免把真实的 wrapped residual 改写成累计相位。每个图标题列出 mask 样本数和接收横向范围。\n\n');
fprintf(fid,'## 7. 全部 15 组总结热力表\n\n'); local_fig_caption(fid,'06_all_case_metric_heatmap','图6：左侧颜色表示 log10(Ec)，单元格标注原始 Ec；右侧以独立色阶显示 improvement，单元格标注 M0/M1 倍数。每行是一个场景和 sigma。');
fprintf(fid,'完整数据可读 `results/visualization/pe_bie_gaussian_sigma_angle_sensitivity/gaussian_sigma_common_mask_summary.csv`；其中同时保留 common-mask 新指标与旧 pairwise-mask 指标，并有 `n_common` 及 mask x 边界供审计。\n\n');
fprintf(fid,'## 8. 结论与适用范围\n\n');
fprintf(fid,'1. 宽角谱端 weak_low_K 的 kz-aware 改善约 `%.2f` 倍；strong-height 约 `%.2f` 倍，weak-high-K 约 `%.2f` 倍，体现明显场景依赖。\n',local_ratio(T,'weak_low_K',.3),local_ratio(T,'strong_height_low_K',.3),local_ratio(T,'weak_high_K',.3));
fprintf(fid,'2. sigma 增大时传播角谱变窄，Model-0 中有限入射角相关的相位误差下降；两种 PE 的总误差逐渐接近，表示 kz 修正的增益收缩到很小，并不代表它们都与 BIE 完全一致。\n');
fprintf(fid,'3. sigma=5 m 时，common-mask Ec(M0/M1/FK) 的逐场景数值如下：\n\n| case | M0 | M1 | Full Kirchhoff |\n|---|---:|---:|---:|\n');
for i=1:3,q=T(strcmp(T.scene,a.scene_names{i})&T.sigma_m==5,:);fprintf(fid,'| %s | %.4g | %.4g | %.4g |\n',q.scene{1},q.Model0_Ec_common,q.Model1_Ec_common,q.FullKirchhoff_Ec_common);end
fprintf(fid,'\n4. weak_low_K 的宽角误差主要有可被 kz-aware 削弱的相位部分；strong-height 仍有相位残差；weak_high_K 窄角端的幅度残差仍不可忽略。\n');
fprintf(fid,'5. 现有数据足以支持开展 receiver-depth sensitivity test：已有 15 个源宽度/场景的完整复场和明确 mask，可将接收深度作为下一项受控维度；但需要锁定统一指标定义并逐深度记录有效 mask，避免跨深度误读为同一采样区域。该建议是可行性判断，不表示本轮已执行 depth sweep。\n\n');
fprintf(fid,'### 原有指标与本轮统一指标的区分\n\n');
fprintf(fid,'原始报告的 pairwise mask 指标来自原始 MAT 中 `metrics_Model0_BIE`、`metrics_Model1_BIE`、`metrics_FK_BIE`；本报告的曲线、表和改善倍数全部由保存的四种原始复场在同一个四场交集 mask 上重算。CSV 中的 `old_*_pairwise` 列仅供数值核对，不能与 common-mask 列混为同一评价口径。原始报告的 mask 定义为：%s\n\n',a.old_mask_definition);
fprintf(fid,'## 文件\n\n- 绘图脚本：`scripts/reporting/plot_gaussian_sigma_angle_sensitivity_atlas.m`\n- 完整 15 组 CSV：`%s`\n- 可审计 MAT：`%s`\n- 数据源：`%s`、`%s`\n',local_rel(root,csv),local_rel(root,mat),local_rel(root,a.source_mat_files{1}),local_rel(root,a.source_mat_files{2}));
end

function local_fig_caption(fid,stem,caption)
if endsWith(strtrim(caption),'。'), suffix=''; else, suffix='。'; end
fprintf(fid,'![%s](../results/visualization/pe_bie_gaussian_sigma_angle_sensitivity/figures/%s.png)\n\n**%s%s**\n\n[PDF 矢量版](../results/visualization/pe_bie_gaussian_sigma_angle_sensitivity/figures/%s.pdf)\n\n',stem,stem,caption,suffix,stem);
end
function val=local_ratio(T,scene,sig)
q=T(strcmp(T.scene,scene)&abs(T.sigma_m-sig)<1e-12,:); assert(isscalar(q.sigma_m)); val=q.improvement_common;
end
function f=local_figure(inches)
f=figure('Visible','off','Color','w','Units','inches','Position',[1 1 inches]);
set(f,'PaperUnits','inches','PaperPosition',[0 0 inches],'PaperSize',inches);
end
function local_save_figure(f,fd,stem)
set(findall(f,'Type','axes'),'FontName','Arial','FontSize',9,'LineWidth',0.8,'TickLabelInterpreter','latex');
exportgraphics(f,fullfile(fd,[stem '.png']),'Resolution',300);
exportgraphics(f,fullfile(fd,[stem '.pdf']),'ContentType','vector'); close(f);
end
function c=local_model_colors()
c=[0.000 0.447 0.741;0.850 0.325 0.098;0.466 0.674 0.188];
end
function s=local_scene_title(name)
switch name
    case 'weak_low_K', s='Weak height / low K';
    case 'strong_height_low_K', s='Strong height / low K';
    case 'weak_high_K', s='Weak height / high K';
    otherwise, s=name;
end
end
function p=local_rel(root,p)
p=strrep(p,[root filesep],''); p=strrep(p,filesep,'/');
end
