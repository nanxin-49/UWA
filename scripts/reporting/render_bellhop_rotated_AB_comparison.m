function receipt=render_bellhop_rotated_AB_comparison(A_file,B_file,output_dir)
%RENDER_BELLHOP_ROTATED_AB_COMPARISON Saved native A/B grids, no interpolation.
assert(~isfolder(output_dir)); [~,ram]=memory; assert(ram.PhysicalMemory.Available>=1.5*2^30,'Comparison figure RAM gate'); mkdir(output_dir); start=tic;
A=load(A_file,'total','envelope','plan'); B=load(B_file,'total','envelope','plan');
assert(isequal(A.plan.depths_m,B.plan.depths_m)&&isequal(A.plan.x_m,B.plan.x_m));
x=B.plan.x_m(:).'; z=B.plan.depths_m(:); delta=abs(A.total-B.total)./B.envelope;
checks=table(strings(0,1),false(0,1),zeros(0,1),'VariableNames',{'check','pass','max_difference'});
fig=figure('Visible','off','Color','w','Position',[30 30 1500 1000],'Tag','AB_comparison_owned'); tiles=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for j=1:4
    ax=nexttile(tiles); if mod(j,2)==1, field='A'; P=A.total; label='A: Q=1.25m'; else, field='B'; P=B.total; label='B: Q=0.625m'; end
    plot_strips(ax,-20*log10(abs(P)),[0 5000],field); if j>2, xlim(ax,[-1 1]); label=[label ' | central x=-1 to1m']; end
    clim(ax,[20 70]); colormap(ax,parula(256)); cb=colorbar(ax); cb.Label.String='TL (dB, fixed1m pressure reference)'; title(ax,label);
end
title(tiles,{'A/B coherent N16 field | only Q spacing changes','5000m Munk, Tx4950m,10kHz | C/D and FFP comparison not run'},'FontSize',16);
save_pair(fig,'F01_AB_full_and_central');
fig=figure('Visible','off','Color','w','Position',[30 30 1650 750],'Tag','AB_comparison_owned'); tiles=tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
for j=1:3
    ax=nexttile(tiles);
    if j==1, C=-20*log10(abs(A.total)); field='A'; label='A near-bottom'; limits=[20 70];
    elseif j==2, C=-20*log10(abs(B.total)); field='B'; label='B near-bottom'; limits=[20 70];
    else, C=delta; field='delta'; label='|pA-pB| / B envelope'; limits=[0 .02]; end
    plot_strips(ax,C,[4900 5000],field); clim(ax,limits); colormap(ax,parula(256)); cb=colorbar(ax); if j<3, cb.Label.String='TL (dB)'; else, cb.Label.String='Envelope-normalized difference'; end; title(ax,label);
end
title(tiles,{'A/B native0.01m near-bottom samples | same401 beams, step1m and receiver grid','Difference panel limit0.02 is the original comparison gate; not FFP error'},'FontSize',16);
save_pair(fig,'F02_AB_near_bottom');
fig=figure('Visible','off','Color','w','Position',[30 30 1450 650],'Tag','AB_comparison_owned'); tiles=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact'); ix=find(x==0);
ax=nexttile(tiles); semilogy(ax,z,abs(A.total(:,ix)),'--','DisplayName','A Q1.25'); hold(ax,'on'); semilogy(ax,z,abs(B.total(:,ix)),'-','DisplayName','B Q0.625'); grid(ax,'on'); xlim(ax,[0 5000]); xlabel(ax,'Physical depth z (m)'); ylabel(ax,'|p| / fixed1m reference'); legend(ax,'Location','best'); title(ax,'Actual axis samples; true source remains masked');
ax=nexttile(tiles); phase=angle(A.total(:,ix).*conj(B.total(:,ix))); phase(abs(B.total(:,ix))<.05*B.envelope(:,ix))=NaN; plot(ax,z,phase); grid(ax,'on'); xlim(ax,[0 5000]); xlabel(ax,'Physical depth z (m)'); ylabel(ax,'arg(pA conj(pB)) (rad)'); title(ax,'Effective phase comparison; original gate0.05rad');
title(tiles,'A/B Q-sampling comparison; no fitted phase, gain or coordinate resampling','FontSize',16); save_pair(fig,'F03_AB_axis');
receipt=struct('status','AB_SAVED_FIELD_FIGURES','failed_checks',nnz(~checks.pass),'checks_count',height(checks),'elapsed_s',toc(start),'A_file',A_file,'A_sha256',black_file_hash(A_file), ...
    'B_file',B_file,'B_sha256',black_file_hash(B_file),'output_dir',output_dir,'native_grid',[10981 201],'new_physics_solves',0, ...
    'encoding','Original uniform depth strips with complete coordinate vectors; no data resampling, source fitting or smoothing');
writetable(checks,fullfile(output_dir,'figure_checks.csv')); write(fullfile(output_dir,'checkpoint.json'),jsonencode(receipt,PrettyPrint=true));
html=['<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>B组与A组对照</title><style>body{font:16px/1.6 sans-serif;max-width:1650px;margin:24px auto;padding:0 20px}img{max-width:100%;height:auto}</style>' ...
'<h1>B组Q加密：A/B完整复场对照</h1><p>仅Q间隔1.25→0.625m；5000m Munk、源4950m、10kHz、N16、401反射束、原接收网格和步长不变。C/D及FFP新计算未执行，模型假设不等于论文原参数。</p>' ...
'<h2>全深度及中央区域</h2><img src="F01_AB_full_and_central.png" alt="A和B的N16相干传播损失并排显示，全深度和中央横向区域。">' ...
'<h2>近底干涉与复数差</h2><img src="F02_AB_near_bottom.png" alt="A/B在4900至5000m的原生0.01m接收网格，以及归一到B分量包络的复场差。">' ...
'<h2>轴线复场</h2><img src="F03_AB_axis.png" alt="A/B轴幅度和有效相位差，保留真实源点掩膜和弱场缺口。">' ...
'<p>TL=-20log10|p|，固定1m参考；20–70dB仅显示范围，原复数和超限值保存在MAT。没有插值、平滑、补轴或任意幅相拟合。全深度5m只作概览；各图同名FIG可编辑。阶段通过与否以中文报告和原始CSV门槛为准。</p></html>']; write(fullfile(output_dir,'index.html'),html); disp(jsonencode(receipt));
    function plot_strips(ax,C,limits,field)
        hold(ax,'on'); rows={find(z<=4900),find(z>=4900)};
        for m=1:2
            r=rows{m}; assert(max(abs(diff(z(r))-mean(diff(z(r)))))<1e-9);
            im=imagesc(ax,x,z(r),C(r,:)); im.Tag=sprintf('%s_strip%d',field,m); im.UserData=struct('field',field,'rows',r); im.AlphaData=isfinite(C(r,:)); im.Interpolation='nearest';
        end
        set(ax,'YDir','reverse','Color',[.9 .9 .9],'FontSize',12); xlim(ax,[-5 5]); ylim(ax,limits); xlabel(ax,'Physical x (m)'); ylabel(ax,'Physical depth z (m)');
    end
    function save_pair(fig,stem)
        fpath=fullfile(output_dir,[stem '.fig']); ppath=fullfile(output_dir,[stem '.png']); assert(~isfile(fpath)&&~isfile(ppath)); drawnow; savefig(fig,fpath); exportgraphics(fig,ppath,'Resolution',140); close(fig);
        h=openfig(fpath,'new','invisible'); cleanup=onCleanup(@()close(h)); images=findall(h,'Type','image');
        for m=1:numel(images)
            u=images(m).UserData; r=u.rows;
            if strcmp(u.field,'A'), expected=-20*log10(abs(A.total(r,:))); elseif strcmp(u.field,'B'), expected=-20*log10(abs(B.total(r,:))); else, expected=delta(r,:); end
            pass=isequaln(images(m).CData,expected)&&isequal(images(m).XData(:),x(:))&&isequal(images(m).YData(:),z(r));
            use=isfinite(expected); diff=max(abs(images(m).CData(use)-expected(use))); checks(end+1,:)={string([stem '_' images(m).Tag]),pass,diff};
        end
        info=imfinfo(ppath); checks(end+1,:)={string([stem '_PNG']),info.Width>1000&&info.Height>500,0};
    end
end
function write(path,text)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); c=onCleanup(@()fclose(fid)); fprintf(fid,'%s',text);
end