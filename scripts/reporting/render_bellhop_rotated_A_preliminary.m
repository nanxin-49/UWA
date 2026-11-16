function receipt=render_bellhop_rotated_A_preliminary(run_dir,output_dir)
%RENDER_BELLHOP_ROTATED_A_PRELIMINARY Saved A data; native uniform strips only.
assert(~isfolder(output_dir)); [~,ram]=memory; assert(ram.PhysicalMemory.Available>=2^30,'Plot RAM gate'); mkdir(output_dir); start=tic;
s=load(fullfile(run_dir,'assembly','result.mat'),'total','direct','surface_family','bottom_family','envelope','plan');
coarse=load(fullfile(run_dir,'assembly','sum_through_order12.mat'),'total'); x=s.plan.x_m(:).'; z=s.plan.depths_m(:);
assert(numel(z)==10981 && numel(x)==201 && nnz(s.plan.source_singularity_mask)==1);
fields={'direct','surface_family','bottom_family','total'}; labels={'Direct charts','All orders starting at surface','All orders starting at bottom','Coherent total through N=16'};
checks=table(strings(0,1),false(0,1),zeros(0,1),'VariableNames',{'check','pass','max_difference'});
fig=figure('Visible','off','Color','w','Position',[30 30 1500 1000],'Tag','A_preliminary_owned');
tiles=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for j=1:4
    ax=nexttile(tiles); plot_strips(ax,-20*log10(abs(s.(fields{j}))),[0 5000],fields{j});
    clim(ax,[20 70]); colormap(ax,parula(256)); cb=colorbar(ax); cb.Label.String='TL (dB, fixed 1m pressure reference)'; title(ax,labels{j});
end
title(tiles,{'A preliminary finite-order field | 5000m Munk, Tx4950m,10kHz','N12-to16 gates passed; Q/step/beam convergence and FFP comparison pending'},'FontSize',16);
save_pair(fig,'F01_A_components_and_total');
fig=figure('Visible','off','Color','w','Position',[30 30 1500 750],'Tag','A_preliminary_owned'); tiles=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(tiles); plot_strips(ax,-20*log10(abs(s.total)),[4900 5000],'total_near'); clim(ax,[20 70]); colormap(ax,parula(256)); cb=colorbar(ax); cb.Label.String='TL (dB)'; title(ax,'Native 0.01m near-bottom samples');
ax=nexttile(tiles); delta=abs(s.total-coarse.total)./s.envelope; plot_strips(ax,delta,[4900 5000],'tail_near'); clim(ax,[0 2e-5]); colormap(ax,parula(256)); cb=colorbar(ax); cb.Label.String='|p16-p12| / component envelope'; title(ax,{'N12-to16 change (not FFP error)','Original limit=0.001; displayed range0 to2e-5'});
title(tiles,{'A preliminary near-bottom field | same physical x,z and model','True point source remains masked; no smoothing, interpolation or source fitting'},'FontSize',16);
save_pair(fig,'F02_A_near_bottom_and_order_change');
fig=figure('Visible','off','Color','w','Position',[30 30 1450 650],'Tag','A_preliminary_owned'); tiles=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ix=find(x==0); ax=nexttile(tiles); semilogy(ax,z,abs(coarse.total(:,ix)),'--','DisplayName','N12'); hold(ax,'on'); semilogy(ax,z,abs(s.total(:,ix)),'-','DisplayName','N16'); grid(ax,'on'); xlabel(ax,'Physical depth z (m)'); ylabel(ax,'|p| / fixed1m reference'); xlim(ax,[0 5000]); legend(ax,'Location','best'); title(ax,'Actual axis samples; source point omitted');
ax=nexttile(tiles); phase=angle(coarse.total(:,ix).*conj(s.total(:,ix))); phase(abs(s.total(:,ix))<.05*s.envelope(:,ix))=NaN; plot(ax,z,phase,'Color',[.1 .35 .65]); grid(ax,'on'); xlim(ax,[0 5000]); xlabel(ax,'Physical depth z (m)'); ylabel(ax,'arg(p12 conj(p16)) (rad)'); title(ax,'Effective-phase samples; original limit0.01rad');
title(tiles,'A axis order comparison; global5m samples give an overview, not all fringes','FontSize',16);
save_pair(fig,'F03_A_axis_order_comparison');
receipt=struct('status','A_PRELIMINARY_FIGURES_COMPLETE','elapsed_s',toc(start),'output_dir',output_dir,'failed_checks',nnz(~checks.pass), ...
    'field_file',fullfile(run_dir,'assembly','result.mat'),'field_sha256',black_file_hash(fullfile(run_dir,'assembly','result.mat')), ...
    'native_grid',[10981 201],'image_encoding','Two overlapping native uniform strips z<=4900 and z>=4900; shared4900row has identical value; no interpolation/resampling', ...
    'source_mask_count',1,'TL_limits_dB',[20 70],'qualification','A-only finite-order preliminary; B/C/D and FFP comparison unrun');
writetable(checks,fullfile(output_dir,'figure_checks.csv')); fid=fopen(fullfile(output_dir,'checkpoint.json'),'w','n','UTF-8'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid);
html=['<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>A组初步完整场</title><style>body{font:16px/1.6 sans-serif;max-width:1500px;margin:24px auto;padding:0 20px}img{max-width:100%;height:auto}</style>' ...
'<h1>A组：N16初步完整接收网格场</h1><p>5000m Munk，Tx4950m，10kHz，全向三维点源和示例流体底。A的截断与原始数据审计通过；B/C/D和FFP共点比较尚未执行，非论文最终复现。</p>' ...
'<h2>直达与反射分量族的复数相干总场</h2><p>两反射族分别按最先命中的海面/底面分类，均包含后续混合反射；不是只含单次反射。TL=-20log10|p|，20至70dB仅显示范围，原始复声压未裁剪。</p><img src="F01_A_components_and_total.png" alt="直达、先海面反射族、先海底反射族及N16相干总场，真实源点掩膜。">' ...
'<h2>近底原始0.01m网格</h2><img src="F02_A_near_bottom_and_order_change.png" alt="近底N16传播损失和N12至N16的包络归一变化，后者不是FFP误差。">' ...
'<h2>轴线截断对照</h2><img src="F03_A_axis_order_comparison.png" alt="轴线上原始N12和N16幅度及有效点相位变化，保留源点和弱场缺口。">' ...
'<p>各图同名FIG可编辑；原始复数MAT位于相邻validation结果的assembly/result.mat。没有平滑、插值、任意幅相拟合；全深度5m网格只作概览，未解析所有波长尺度驻波。</p></html>'];
fid=fopen(fullfile(output_dir,'index.html'),'w','n','UTF-8'); fprintf(fid,'%s',html); fclose(fid); disp(jsonencode(receipt));
    function plot_strips(ax,C,limits,field)
        hold(ax,'on'); rows={find(z<=4900),find(z>=4900)};
        for k=1:2
            iz=rows{k}; dz=diff(z(iz)); assert(max(abs(dz-dz(1)))<1e-9);
            im=imagesc(ax,x,z(iz),C(iz,:)); im.Tag=sprintf('%s_strip%d',field,k); im.UserData=struct('field',field,'rows',iz);
            im.AlphaData=isfinite(C(iz,:)); im.Interpolation='nearest';
        end
        set(ax,'YDir','reverse','Color',[.9 .9 .9],'FontSize',12); xlim(ax,[-5 5]); ylim(ax,limits); xlabel(ax,'Physical x (m)'); ylabel(ax,'Physical depth z (m)');
    end
    function save_pair(fig,stem)
        figpath=fullfile(output_dir,[stem '.fig']); pngpath=fullfile(output_dir,[stem '.png']); assert(~isfile(figpath)&&~isfile(pngpath));
        drawnow; savefig(fig,figpath); exportgraphics(fig,pngpath,'Resolution',140); close(fig);
        % Reopen trusted FIG and independently bind image colors to saved data.
        h=openfig(figpath,'new','invisible'); cleanup=onCleanup(@()close(h)); images=findall(h,'Type','image');
        for k=1:numel(images)
            u=images(k).UserData; if isempty(u), continue; end
            field=u.field; rows=u.rows;
            if strcmp(field,'tail_near'), expected=abs(s.total(rows,:)-coarse.total(rows,:))./s.envelope(rows,:);
            else, if strcmp(field,'total_near'), field='total'; end; expected=-20*log10(abs(s.(field)(rows,:))); end
            actual=images(k).CData; finite=isfinite(expected); difference=max(abs(actual(finite)-expected(finite))); if isempty(difference), difference=0; end
            pass=isequaln(actual,expected) && isequal(images(k).YData(:),z(rows)) && isequal(images(k).XData(:),x(:));
            checks(end+1,:)={string([stem '_' images(k).Tag]),pass,difference}; %#ok<AGROW>
        end
        info=imfinfo(pngpath); checks(end+1,:)={string([stem '_PNG_nonempty']),info.Width>1000&&info.Height>500,0}; %#ok<AGROW>
    end
end