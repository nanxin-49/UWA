function receipt=render_munk_ffp_5000m(run_id,audit_id,version_suffix)
%RENDER_MUNK_FFP_5000M Five saved-data figure groups plus offline HTML.
arguments
    run_id (1,:) char = '20261112_full_01'
    audit_id (1,:) char = 'audit_01'
    version_suffix (1,:) char = ''
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'scripts','validation','support'),'-begin');
d=fullfile(root,'results','validation','munk_ffp_5000m',run_id); a=fullfile(d,audit_id);
assert(isempty(regexp(version_suffix,'[^A-Za-z0-9_-]','once')));
v=fullfile(root,'results','visualization','munk_ffp_5000m',[run_id version_suffix]);
assert(~isfolder(v),'Preserve existing visualization attempts'); mkdir(v); t=tic; diary(fullfile(v,'diary.log'));
receipt=struct('stage','render','status','RUNNING','elapsed_s',0,'exit_status',NaN);
try
    audit=jsondecode(fileread(fullfile(a,'receipt.json')));
    if isfield(audit,'final_job'), last=audit.final_job; else, last=6; end
    s=load(fullfile(d,sprintf('J%02d.mat',last)),'output'); o=s.output;
    qual=audit.qualification;
    source='C:/Users/ASUS/.codex/worktrees/a28b/Explain/results/validation/bellhop_fig12b/user_supplied_fig12b.png';
    assert(isfile(source),'Paper raster missing; do not replace with a fabricated image');
    copyfile(source,fullfile(v,'paper_supplied_original.png')); paper=imread(source); crop=paper(31:610,131:700,:);
    f=local_fig('F01'); h=(.35*1600)*(580/570)/1000;
    ax1=axes(f,'Position',[.075 .2 .35 h]); local_field(ax1,o,[0 5000]);
    title(ax1,{'Independent layered-wave / exact J0','10 kHz; source 4950 m; no axis filling'});
    cb=colorbar(ax1); cb.Direction='reverse'; cb.Position=[.435 .2 .014 h]; cb.Label.String='TL=-20 log10 |4 pi G| (dB)';
    ax2=axes(f,'Position',[.56 .2 .35 h]); image(ax2,[-5 5],[0 5000],crop); set(ax2,'YDir','reverse');
    xlim(ax2,[-5 5]); ylim(ax2,[0 5000]); xlabel(ax2,'Range x (m)'); ylabel(ax2,'Depth z (m)');
    title(ax2,{'Supplied paper Fig.12(b) raster','Original color/plot-box aspect preserved'});
    colormap(ax2,flipud(jet(256))); clim(ax2,[20 70]); cb=colorbar(ax2); cb.Direction='reverse'; cb.Position=[.92 .2 .014 h]; cb.Label.String='Source raster scale (dB)';
    ax1.Position=[.075 .2 .35 h]; ax2.Position=[.56 .2 .35 h];
    annotation(f,'textbox',[.05 .03 .9 .075],'String',{['Numeric status: ' qual], ...
        'Raster comparison is qualitative. Gray = exact zero or true source singularity; black x marks source. No fitted gain.'}, ...
        'EdgeColor','none','HorizontalAlignment','center','Interpreter','none','FontSize',13);
    local_export(f,v,'F01_full_field_vs_paper');

    axis_idx=o.range_m==0; p=o.pressure_1m(:,axis_idx); z=o.depth_m;
    cum=[0;cumsum(diff(o.cfg.layer_edges_m)./o.cfg.layer_c_mps)];
    delay=abs(interp1(o.cfg.layer_edges_m,cum,z)-interp1(o.cfg.layer_edges_m,cum,o.cfg.source_depth_m));
    referenced=angle(p.*exp(-1i*2*pi*o.cfg.frequency_hz*delay));
    f=local_fig('F02'); tl=tiledlayout(f,1,3,'TileSpacing','compact','Padding','compact');
    ax=nexttile(tl); semilogx(ax,abs(p),z,'-'); set(ax,'YDir','reverse'); grid(ax,'on'); xlabel(ax,'|p|, unit-1m source reference'); ylabel(ax,'Depth z (m)');
    title(ax,'Direct r=0 calculation; true source is NaN'); yline(ax,4950,':k','source');
    ax=nexttile(tl); plot(ax,angle(p),z,'.','MarkerSize',3); set(ax,'YDir','reverse'); grid(ax,'on');
    xlabel(ax,'Wrapped phase (rad)'); ylabel(ax,'Depth z (m)'); title(ax,'Actual raw phase; no global 5m unwrap'); xlim(ax,[-pi pi]);
    ax=nexttile(tl); plot(ax,referenced,z,'.','MarkerSize',3); set(ax,'YDir','reverse'); grid(ax,'on');
    xlabel(ax,'Phase relative to fixed one-pass delay (rad)'); ylabel(ax,'Depth z (m)'); title(ax,'Display reference only; not estimated group delay'); xlim(ax,[-pi pi]);
    title(tl,['F02 axial complex-pressure profiles: ' qual],'Interpreter','none');
    local_export(f,v,'F02_axis_magnitude_phase');

    f=local_fig('F03'); tl=tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
    ax=nexttile(tl); local_field(ax,o,[4900 5000]); title(ax,'Actual 0.01 m near-bottom receiver samples'); colorbar(ax);
    ax=nexttile(tl); local_field(ax,o,[4970 4972]); title(ax,'4970--4972 m zoom; no display interpolation'); colorbar(ax);
    ax=nexttile(tl); [~,ix]=min(abs(o.range_m-1)); k=z>=4970 & z<=4980;
    plot(ax,z(k),-20*log10(abs(o.pressure_1m(k,ix)))); grid(ax,'on'); xlabel(ax,'Depth z (m)'); ylabel(ax,'TL (dB)'); title(ax,'x=1 m depth interference, raw samples');
    ax=nexttile(tl); k=z>=4948 & z<=4952;
    semilogy(ax,z(k),abs(p(k)),'-'); grid(ax,'on'); xlabel(ax,'Depth z (m)'); ylabel(ax,'|p(r=0)|');
    title(ax,'Source neighborhood: excluded source is a gap'); xline(ax,4950,':k','true point source');
    title(tl,'F03 local resolution: maximum dz < lambda/10; finite source size is not assumed');
    local_export(f,v,'F03_near_bottom_resolution');

    pre=fullfile(root,'results','validation','munk_ffp_5000m','20261112_preflight_01');
    sf=load(fullfile(pre,'uniform_free.mat'),'uniform_free','truth'); sr=load(fullfile(pre,'uniform_release.mat'),'uniform_release','mirror');
    f=local_fig('F04'); tl=tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
    long_reference=load(fullfile(d,'long_uniform_controls.mat'),'free','release','exact_free','exact_release');
    sf.uniform_free=long_reference.free; sf.truth=long_reference.exact_free;
    sr.uniform_release=long_reference.release; sr.mirror=long_reference.exact_release;
    long_receipt=jsondecode(fileread(fullfile(d,'long_uniform_receipt.json')));
    [~,iz]=min(abs(sf.uniform_free.depth_m-4900)); x=sf.uniform_free.range_m;
    ax=nexttile(tl); plot(ax,x,real(sf.truth(iz,:)),'ko',x,real(sf.uniform_free.green(iz,:)),'b+', ...
        x,imag(sf.truth(iz,:)),'ks',x,imag(sf.uniform_free.green(iz,:)),'rx');
    grid(ax,'on'); xlabel(ax,'Range x (m)'); ylabel(ax,'Uniform Green, source z=4950m / receiver z=4900m');
    legend(ax,{'Exact Re','Numeric Re','Exact Im','Numeric Im'}); title(ax,'Existing analytic point source, no subtraction');
    ax=nexttile(tl); plot(ax,x,real(sr.mirror(iz,:)),'ko',x,real(sr.uniform_release.green(iz,:)),'b+', ...
        x,imag(sr.mirror(iz,:)),'ks',x,imag(sr.uniform_release.green(iz,:)),'rx');
    grid(ax,'on'); xlabel(ax,'Range x (m)'); ylabel(ax,'Pressure-release mirror Green'); title(ax,'Existing image solution, no fitted phase/gain');
    ax=nexttile(tl); semilogy(ax,x,abs(sf.uniform_free.green(iz,:)-sf.truth(iz,:)),'o-', ...
        x,abs(sr.uniform_release.green(iz,:)-sr.mirror(iz,:)),'s-'); grid(ax,'on'); xlabel(ax,'Range x (m)'); ylabel(ax,'Absolute complex error'); legend(ax,{'Free','Mirror'});
    ax=nexttile(tl); axis(ax,'off'); metrics=audit.metrics;
    text(ax,.02,.95,{sprintf('5000m uniform max per-depth complex error: %.5g',long_receipt.metrics.max_per_depth_complex_rel), ...
        sprintf('Short-scale uniform/mirror errors: %.5g / %.5g',metrics.uniform_free_rel,metrics.uniform_mirror_rel), ...
        sprintf('5000m uniform maximum phase error: %.5g rad',long_receipt.metrics.max_phase_error_rad), ...
        sprintf('Surface residual: %.5g',metrics.surface_rel), ...
        sprintf('Final vertical one-pass time error: %.5g s',metrics.max_onepass_delay_error_s)}, ...
        'Units','normalized','VerticalAlignment','top','FontSize',15,'Interpreter','none');
    title(tl,'F04 existing Green/image references at actual 5000m scale: full kernel integrated without analytic subtraction');
    local_export(f,v,'F04_uniform_reference_errors');

    cv=readtable(fullfile(a,'convergence.csv')); f=local_fig('F05'); tl=tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
    comparisons=cellstr(unique(string(cv.comparison),'stable')); nc=numel(comparisons); regions={'upper_column','near_bottom','axis','off_axis'};
    labels={'2.5 to 1.25m','K80 to 100','early .25 to .125','offset .2 to .1','1.25 to .625m'};
    if nc==6, labels{6}='final .125 to .0625'; end
    vars={'complex_rel_l2','phase_p95_rad','tl_p95_db'}; limits=[1e-3 .01 .05];
    for j=1:3
        ax=nexttile(tl); hold(ax,'on');
        for r=1:4
            values=zeros(1,nc);
            for c=1:nc, row=strcmp(cv.comparison,comparisons{c}) & strcmp(cv.region,regions{r}); values(c)=cv.(vars{j})(row); end
            display_values=max(values,1e-12); % display floor only; CSV retains exact zeros
            plot(ax,1:nc,display_values,'-o','DisplayName',regions{r});
        end
        set(ax,'YScale','log'); yline(ax,limits(j),'--k','required gate','HandleVisibility','off'); grid(ax,'on');
        xticks(ax,1:nc); xticklabels(ax,labels); xtickangle(ax,15); title(ax,[vars{j} ' (zero errors at 1e-12 display floor)'],'Interpreter','none');
        if j==1, legend(ax,'Location','best','Interpreter','none'); end
    end
    ax=nexttile(tl); axis(ax,'off');
    text(ax,.02,.95,{['Qualification: ' qual],sprintf('Failed required checks: %d',audit.failed_checks), ...
        sprintf('Maximum local dz: %.9g m',metrics.max_patch_dz_m),sprintf('Local lambda/10: %.9g m',metrics.lambda_min_m/10), ...
        'Early layer/quadrature failures retained; finer final pairs required.', ...
        'Each error uses the same region and explicit reference denominator.', ...
        'Symmetry is also a consequence of the J0 operator; it alone is not accuracy.'}, ...
        'Units','normalized','VerticalAlignment','top','FontSize',14,'Interpreter','none');
    title(tl,'F05 independent convergence tests: truncation, quadrature, contour and vertical layers');
    local_export(f,v,'F05_convergence_qualification');

    % Numeric structure comparison: actual axis and r=0.1/1/5m at named depths.
    rows=struct([]);
    for zz=[1000 3000 4000 4900 4970 5000]
        [~,zi]=min(abs(o.depth_m-zz)); [~,i0]=min(abs(o.range_m));
        for xx=[0 .1 1 5]
            [~,xi]=min(abs(o.range_m-xx)); pval=o.pressure_1m(zi,xi);
            rows=[rows struct('depth_m',o.depth_m(zi),'range_m',o.range_m(xi), ...
                'real_pressure_1m',real(pval),'imag_pressure_1m',imag(pval),'tl_db',-20*log10(abs(pval)), ...
                'axis_relative_brightness_db',20*log10(abs(o.pressure_1m(zi,i0))/abs(pval)))]; %#ok<AGROW>
        end
    end
    writetable(struct2table(rows),fullfile(v,'structure_metrics.csv'));
    stems={'F01_full_field_vs_paper','F02_axis_magnitude_phase','F03_near_bottom_resolution','F04_uniform_reference_errors','F05_convergence_qualification'};
    captions={'全场：独立三维点源轴对称切面与用户论文栅格并排。原图只裁绘图区、保留原始副本与比例；对比只定性，不拟合源级。原生计算直接包含r=0。全场概览只显示真实0:5:5000m均匀子集，近底使用0.01m网格；概览不能解析全深度细驻波条纹。', ...
        '中央轴：真实复声压幅度、原始wrapped相位及固定one-pass参考相位。单频相位不等于群时延；全局5m样本不做错误解缠。源点本身是NaN标记，其他轴点直接积分。', ...
        '近底：4900–5000m按0.01m采样，保留4970–4972m放大与深度剖面；所有图来自原始复场，不复制邻轴、不做显示插值。', ...
        '均匀基准：复用已有Green、镜像和Weyl；本图显示5000m/源4950m同尺度对照，短尺度对照也保留；数值控制未用解析奇点减除，因此不是把解析值原样返回的恒等检查。', ...
        '收敛：独立改变K截断、积分面板宽度、数值轮廓和垂向层厚，在上部/近底/轴/离轴分别报告。小误差不自动证明论文未知参数一致。'};
    fid=fopen(fullfile(v,'index.html'),'w','n','UTF-8'); clean=onCleanup(@()fclose(fid));
    fprintf(fid,'<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>5000m Munk 独立波动参考</title><style>body{font:17px system-ui;max-width:1500px;margin:24px auto;padding:0 20px;background:#f4f5f7;color:#15202b}article{background:white;border:1px solid #ccd2d8;padding:20px;margin:24px 0}img{width:100%%;height:auto}.status{padding:16px;background:#fff0ce}a{color:#064ca4}</style>');
    fprintf(fid,'<h1>5000m Munk：分层Green / 精确J0波数积分</h1><p class="status">数值资格：<b>%s</b>；必需失败检查%d项。源与底质仍为明确假设，不能据此声称论文原始物理参数被复原。</p>',qual,audit.failed_checks);
    fprintf(fid,'<p>run_id=%s；现有MATLAB MCP串行，10kHz，水深5000m/源4950m。真实轴点直接计算；点源本身排除。没有Bellhop/PE调用或CIR/通信扩展。</p>',run_id);
    fprintf(fid,'<p><a href="../../../../reports/munk_ffp_5000m_design_2026-11-12.md">模型与事前设计</a> | <a href="../../../../reports/munk_ffp_5000m_resource_budget_2026-11-12.md">获批预算范围</a> | <a href="../../../validation/munk_ffp_5000m/%s/%s/convergence.csv">完整误差CSV</a> | <a href="../../../validation/munk_ffp_5000m/%s/%s/axis_complex.csv">原始轴复场CSV</a> | <a href="paper_supplied_original.png">完整用户原图</a> | <a href="structure_metrics.csv">结构指标</a></p>',run_id,audit_id,run_id,audit_id);
    for j=1:5
        fprintf(fid,'<article><h2>%s</h2><p>%s</p><a href="%s.png"><img src="%s.png" alt="%s" loading="lazy"></a><p><a href="%s.fig">可编辑FIG</a> | <a href="%s.png">PNG</a></p></article>',stems{j},captions{j},stems{j},stems{j},captions{j},stems{j},stems{j});
    end
    fprintf(fid,'<p>原始复Green、pressure_1m=4πGreen、完整配置与六组结果在validation同run目录。TL单位为指定单位点源参考，不是校准Pa；示例流体底cp1600/rho1.8/衰减.8dB每波长不是论文原底质。局部网格高分辨率不等于实物源区可信。</p></html>');
    clear clean;
    receipt.status='FIVE_FIGURE_GROUPS_AND_HTML_SAVED'; receipt.exit_status=0; receipt.qualification=qual;
    receipt.figure_stems=stems;
catch err
    if exist('f','var') && isgraphics(f) && strcmp(f.Tag,'munk_ffp_5000m_owned'), close(f); end
    receipt.status='RENDER_FAILED'; receipt.exit_status=1; receipt.error=getReport(err,'extended','hyperlinks','off');
end
receipt.elapsed_s=toc(t); receipt.visualization_dir=v; local_json(fullfile(v,'receipt.json'),receipt); diary off;
fprintf('%s\n',jsonencode(receipt));
end
function f=local_fig(name)
f=figure('Visible','off','Tag','munk_ffp_5000m_owned','Name',name,'Color','w','Units','pixels','Position',[20 20 1600 1000]);
set(f,'PaperUnits','inches','PaperPosition',[0 0 16 10],'PaperSize',[16 10],'PaperPositionMode','manual');
end
function local_field(ax,o,limits)
rows=find(o.depth_m>=limits(1) & o.depth_m<=limits(2)); zplot=o.depth_m(rows);
if numel(zplot)>2 && max(diff(zplot))-min(diff(zplot))>1e-8
    desired=(ceil(limits(1)/5)*5:5:floor(limits(2)/5)*5).';
    [present,rows]=ismember(desired,o.depth_m); assert(all(present),'Missing actual global receiver samples');
    zplot=o.depth_m(rows);
end
p=o.pressure_1m(rows,:); tl=-20*log10(abs(p)); valid=isfinite(tl); tl(~valid)=NaN;
im=imagesc(ax,o.range_m,zplot,tl); im.AlphaData=valid; set(ax,'YDir','reverse','Color',[.78 .78 .78]);
colormap(ax,flipud(jet(256))); clim(ax,[20 70]); xlim(ax,[-5 5]); ylim(ax,limits);
xlabel(ax,'Range x (m)'); ylabel(ax,'Depth z (m)'); hold(ax,'on'); plot(ax,0,4950,'kx','MarkerSize',8,'LineWidth',1.5);
end
function local_export(f,v,stem)
assert(~isfile(fullfile(v,[stem '.fig'])) && ~isfile(fullfile(v,[stem '.png'])));
drawnow; savefig(f,fullfile(v,[stem '.fig'])); print(f,fullfile(v,[stem '.png']),'-dpng','-r100'); close(f);
end
function local_json(path,value)
fid=fopen(path,'w'); assert(fid>=0); cleaner=onCleanup(@()fclose(fid)); fprintf(fid,'%s',jsonencode(value,'PrettyPrint',true));
end
