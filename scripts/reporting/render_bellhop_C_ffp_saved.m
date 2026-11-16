function receipt=render_bellhop_C_ffp_saved(project_root,comparison_dir,output_dir)
% Saved-coordinate visualization only. No solver, interpolation or fitting.
assert(~isfolder(output_dir)); mkdir(output_dir); start=tic; [~,ram]=memory; assert(ram.PhysicalMemory.Available>=1.5*2^30,'Saved-field graphics RAM gate');
s=load(fullfile(comparison_dir,'comparison.mat')); x=s.x(:).'; z=s.z(:); C=s.p_bellhop; F=s.p_ffp; E=s.envelope; clear s
checks=table(strings(0,1),false(0,1),'VariableNames',{'check','pass'}); figures={}; D=abs(C-F)./E; phase=angle(C.*conj(F)); phase(abs(C)<.05*E|abs(F)<.05*E)=NaN;
fig=owned([1700 1000]); t=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for j=1:4
 ax=nexttile(t); g=1+mod(j-1,2); if g==1, P=C; label='Bellhop C N16'; else, P=F; label='FFP J07'; end
 strips(ax,-20*log10(abs(P)),[0 5000],g,'TL'); if j>2, xlim(ax,[-1 1]); label=[label ' | central x=-1..1m']; end; title(ax,label); clim(ax,[20 70]); colorbar(ax);
end
title(t,{'Coherent TL (dB), fixed1m reference | identical native receiving coordinates','Top: full x=-5..5m; bottom: central view | no interpolation or fitted gain/phase'},'FontSize',15); pair(fig,'F01_C_FFP_full_central');
fig=owned([1700 1000]); t=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for j=1:4
 ax=nexttile(t); g=1+mod(j-1,2); if g==1, P=C; label='Bellhop C'; else, P=F; label='FFP J07'; end
 if j<=2, limits=[4900 5000]; else, limits=[4970 4972]; end
 strips(ax,-20*log10(abs(P)),limits,g,'TL'); title(ax,label); clim(ax,[20 70]); colorbar(ax);
end
title(t,{'Near-bottom coherent TL (dB) | original0.01m samples','Top4900..5000m; bottom4970..4972m | same20..70dB color limits'},'FontSize',15); pair(fig,'F02_C_FFP_near_bottom');
fig=owned([1500 650]); t=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact'); ix=find(x==0); ax=nexttile(t);
a=semilogy(ax,z,abs(C(:,ix)),'-','Color',[.1 .1 .1],'LineWidth',1,'DisplayName','Bellhop C'); a.UserData=struct('field',1,'kind','amplitude'); hold(ax,'on'); b=semilogy(ax,z,abs(F(:,ix)),'--','Color',[0 .35 .8],'DisplayName','FFP J07'); b.UserData=struct('field',2,'kind','amplitude'); xlabel(ax,'Depth z (m)'); ylabel(ax,'|p| / fixed1m reference'); xlim(ax,[0 5000]); grid(ax,'on'); legend(ax,[a b],{'Bellhop C','FFP J07'},'Visible','on','Location','northwest'); title(ax,'Raw axis amplitudes; true source remains masked');
ax=nexttile(t); h=plot(ax,z,phase(:,ix)); h.UserData=struct('field',0,'kind','phase'); xlabel(ax,'Depth z (m)'); ylabel(ax,'arg(pC conj(pFFP)) (rad)'); xlim(ax,[0 5000]); grid(ax,'on'); title(ax,'Effective phase: both |p|>=0.05 C envelope'); title(t,'Axis field | no phase or gain alignment','FontSize',15); pair(fig,'F03_C_FFP_axis');
fig=owned([1500 800]); t=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact'); ax=nexttile(t); strips(ax,D,[0 5000],0,'difference'); clim(ax,[0 .001]); colorbar(ax); title(ax,'|pC-pFFP| / C component envelope'); ax=nexttile(t); strips(ax,phase,[0 5000],0,'phase'); clim(ax,[-.015 .015]); colorbar(ax); title(ax,'Effective phase difference (rad)'); title(t,'Native-point differences | display limits are not acceptance gates','FontSize',15); pair(fig,'F04_C_FFP_differences');
% Reuse real pre-existing direct rays; no new trace was started.
family=fullfile(project_root,'results','validation','bellhop_rotated_fig12b'); saved=load(fullfile(family,'20261114_s1_rays_01','result.mat'),'runs');
fig=owned([1500 900]); t=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact'); ray_geometry_checks=0;
for g=1:2
 a=saved.runs{g}; sg=1; if strcmp(a.cfg.direction,'up'), sg=-1; end
 ax=nexttile(t,2*g-1); hold(ax,'on'); bx=nexttile(t,2*g); hold(bx,'on');
 for k=1:numel(a.rays)
  rq=a.rays(k).rq; keep=rq(:,1)>=0 & rq(:,1)<=a.cfg.range_limit_m; r=rq(keep,1); q=rq(keep,2)-a.cfg.transverse_origin_m; zz=a.cfg.source_z_m+sg*r; xx=-sg*q;
  plot(ax,r/1000,q); plot(bx,xx,zz); assert(max(abs(sg*(zz-a.cfg.source_z_m)-r))<1e-8); ray_geometry_checks=ray_geometry_checks+1;
 end
 set(ax,'YDir','reverse'); set(bx,'YDir','reverse'); xlabel(ax,'Chart range r (km)'); ylabel(ax,'q-q0 (m)'); xlabel(bx,'Physical x (m)'); ylabel(bx,'Physical depth z (m)'); grid(ax,'on'); grid(bx,'on'); title(ax,[a.cfg.direction ' direct | rotated']); title(bx,[a.cfg.direction ' direct | physical']);
end
title(t,{'Saved native direct rays: historical S1, Q2.5m, step0.5m','Actual .ray nodes; not C-field production ray output | reflected full trajectories unavailable'},'FontSize',15); pair(fig,'F05_saved_direct_rays');
% High-order reflection logs contain true hit states, not complete trajectories.
file=fullfile(family,'20261115_performance_review_01','optimized_clock_case01','case.perfwall'); fid=fopen(file,'rb','ieee-le'); assert(fid>=0); bytes=fread(fid,inf,'*uint8'); fclose(fid); assert(mod(numel(bytes),420)==0); raw=reshape(bytes,420,[]); ints=reshape(typecast(reshape(raw(1:28,:),[],1),'int32'),7,[]); geom=reshape(typecast(reshape(raw(29:end,:),[],1),'double'),49,[]); assert(size(geom,2)==1313);
fig=owned([1500 700]); t=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact'); ax=nexttile(t); hold(ax,'on'); bx=nexttile(t); hold(bx,'on');
for beam=[1 51 101]
 columns=(beam-1)*13+(1:13); bounce=double(ints(4,columns)+ints(5,columns)); assert(isequal(bounce,1:13)); r=geom(8,columns); q=geom(9,columns); xx=(-1).^(bounce-1).*q; zz=5000*(mod(bounce,2)==0); assert(max(abs(r-(4950+(bounce-1)*5000)))<1e-8);
 scatter(ax,r/1000,q,36,bounce,'filled'); scatter(bx,xx,zz,36,bounce,'filled');
end
xlabel(ax,'Unfolded range r (km)'); ylabel(ax,'Incident q (m)'); title(ax,'Actual wall-hit states, selected3 of101 rays'); xlabel(bx,'Physical x at hit (m)'); ylabel(bx,'Physical depth z (m)'); set(bx,'YDir','reverse'); ylim(bx,[-200 5200]); title(bx,'13 reflections: surface odd, bottom even'); colorbar(ax); colorbar(bx); clim(ax,[1 13]); clim(bx,[1 13]); grid(ax,'on'); grid(bx,'on'); title(t,{'Saved order13 reflection-hit records, Windows qualified performance case','Markers only: no missing path segment invented; full first/high-order ray curves remain a gap'},'FontSize',15); pair(fig,'F06_saved_reflection_hits');
[v,im]=max(D(:),[],'omitnan'); [iz,jx]=ind2sub(size(D),im); [p,ip]=max(abs(phase(:)),[],'omitnan'); [zp,xp]=ind2sub(size(phase),ip); [~,neighbor]=min(abs(x-.1));
receipt=struct('status','SAVED_C_FFP_COMPARISON_FIGURES_RAY_GAP_RETAINED','failed_checks',nnz(~checks.pass),'figure_checks',height(checks),'saved_direct_ray_map_checks',ray_geometry_checks,'new_Bellhop',0,'new_FFP',0,'figures',{figures},'max_normalized_difference',v,'max_difference_location_m',[x(jx) z(iz)],'max_effective_phase_rad',p,'max_phase_location_m',[x(xp) z(zp)],'axis_neighbor_actual_x_m',x(neighbor),'elapsed_s',toc(start),'output_dir',output_dir);
writetable(checks,fullfile(output_dir,'figure_checks.csv')); fid=fopen(fullfile(output_dir,'checkpoint.json'),'w','n','UTF-8'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid);
html=['<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>Bellhop C与FFP共点比较</title><style>body{max-width:1700px;margin:24px auto;padding:0 20px;font:16px/1.65 sans-serif}img{max-width:100%;height:auto}.warning{padding:15px;background:#fff0d5}</style><h1>Bellhop C / 已有FFP J07：原网格共点复场</h1><p>5000m水深、源4950m、10kHz、全向3D点源与同一流体底假设。10981×201点完全相同；C=-conj(native)，FFP=4πGreen；没有增益/相位拟合、插值或平滑。TL固定20–70dB，超界值仍保存在原MAT。</p><p class="warning">本轮不重算A/B/C或FFP。新R模式在RAM启动门停止、PID0，未重试；下面完整声线仅为既有S1直达记录，高阶图只显示真实反射命中点，不能当成完整反射轨迹。</p>'];
for k=1:numel(figures), stem=figures{k}; html=[html '<h2>' stem '</h2><img src="' stem '.png" alt="' stem '"><p><a href="' stem '.fig">可编辑FIG</a></p>']; end %#ok<AGROW>
html=[html '<p>原始共点差CSV及报告保存在本实验comparison目录。C/FFP方法差异、弱场支持和中央亮带结论以中文报告为准；D与新FFP未执行。本页无需联网。</p></html>']; fid=fopen(fullfile(output_dir,'index.html'),'w','n','UTF-8'); fprintf(fid,'%s',html); fclose(fid); disp(jsonencode(receipt));
 function fig=owned(sz), fig=figure('Visible','off','Color','w','Position',[20 20 sz]); end
 function strips(ax,V,limits,g,kind)
  hold(ax,'on'); if limits(1)>=4900, rr={find(z>=limits(1)&z<=limits(2))}; else, rr={find(z<=4900),find(z>=4900)}; end
  for m=1:numel(rr), rows=rr{m}; h=imagesc(ax,x,z(rows),V(rows,:)); h.UserData=struct('field',g,'kind',kind,'rows',rows); h.AlphaData=~isnan(V(rows,:)); h.Interpolation='nearest'; end
  set(ax,'YDir','reverse','Color',[.85 .85 .85],'FontSize',11); xlim(ax,[-5 5]); ylim(ax,limits); xlabel(ax,'Physical x (m)'); ylabel(ax,'Physical depth z (m)'); colormap(ax,parula(256));
 end
 function pair(fig,stem)
  drawnow; savefig(fig,fullfile(output_dir,[stem '.fig'])); exportgraphics(fig,fullfile(output_dir,[stem '.png']),'Resolution',130); close(fig); figures{end+1}=stem;
  if numel(figures)<=4
   h=openfig(fullfile(output_dir,[stem '.fig']),'new','invisible'); close_fig=onCleanup(@()close(h)); ims=findall(h,'Type','image');
   for u=1:numel(ims)
    m=ims(u).UserData; rows=m.rows; if strcmp(m.kind,'TL'), if m.field==1, expected=-20*log10(abs(C(rows,:))); else, expected=-20*log10(abs(F(rows,:))); end; elseif strcmp(m.kind,'difference'), expected=D(rows,:); else, expected=phase(rows,:); end
    pass=isequaln(ims(u).CData,expected)&&isequal(ims(u).XData(:),x(:))&&isequal(ims(u).YData(:),z(rows))&&isequal(logical(ims(u).AlphaData),~isnan(expected)); checks(end+1,:)={string([stem '_image' num2str(u)]),pass}; %#ok<AGROW>
   end
   if isempty(ims)
    ls=findall(h,'Type','line'); for u=1:numel(ls), m=ls(u).UserData; if ~isstruct(m)||~isfield(m,'kind'), continue; end; if strcmp(m.kind,'amplitude'), if m.field==1, expected=abs(C(:,ix)); else, expected=abs(F(:,ix)); end; else, expected=phase(:,ix); end; checks(end+1,:)={string([stem '_line' num2str(u)]),isequal(ls(u).XData(:),z)&&isequaln(ls(u).YData(:),expected)}; end %#ok<AGROW>
   end
  end
 end
end
