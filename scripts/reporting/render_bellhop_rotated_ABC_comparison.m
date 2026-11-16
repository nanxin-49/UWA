function receipt=render_bellhop_rotated_ABC_comparison(A_file,B_file,C_file,output_dir)
% Saved native A/B/C fields; no solve, interpolation or amplitude/phase fitting.
assert(~isfolder(output_dir)); [~,ram]=memory; assert(ram.PhysicalMemory.Available>=1.5*2^30,'Figure RAM gate'); mkdir(output_dir); start=tic;
F={load(A_file,'total','envelope','plan'),load(B_file,'total','envelope','plan'),load(C_file,'total','envelope','plan')};
x=F{3}.plan.x_m(:).'; z=F{3}.plan.depths_m(:); labels={'A Q1.25 h1','B Q0.625 h1','C Q0.625 h0.5'};
for j=1:2, assert(isequal(F{j}.plan.depths_m,F{3}.plan.depths_m)&&isequal(F{j}.plan.x_m,F{3}.plan.x_m)); end
fig=figure('Visible','off','Color','w','Position',[20 20 1800 1000]); t=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
for j=1:6
 ax=nexttile(t); g=mod(j-1,3)+1; strips(ax,-20*log10(abs(F{g}.total)),[0 5000],g); if j>3, xlim(ax,[-1 1]); end
 clim(ax,[20 70]); colorbar(ax); title(ax,labels{g});
end
title(t,{'A/B/C coherent N16 TL (dB) | same receiving coordinates and source reference','Top: full field; bottom: central x=-1 to1m | 5000m Munk, Tx4950m,10kHz'},'FontSize',16); pair(fig,'F01_ABC_full_central');
fig=figure('Visible','off','Color','w','Position',[20 20 1800 750]); t=tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
for g=1:3, ax=nexttile(t); strips(ax,-20*log10(abs(F{g}.total)),[4900 5000],g); clim(ax,[20 70]); colorbar(ax); title(ax,labels{g}); end
title(t,{'A/B/C near-bottom coherent TL (dB) | native0.01m depth spacing','No receiver resampling, smoothing, fitted gain or phase'},'FontSize',16); pair(fig,'F02_ABC_near_bottom');
fig=figure('Visible','off','Color','w','Position',[20 20 1500 650]); t=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact'); ix=find(x==0);
ax=nexttile(t); hold(ax,'on'); for g=1:3, h=semilogy(ax,z,abs(F{g}.total(:,ix)),'DisplayName',labels{g}); h.UserData=struct('group',g,'kind','amplitude'); end
set(ax,'YScale','log'); grid(ax,'on'); xlim(ax,[0 5000]); xlabel(ax,'Physical depth (m)'); ylabel(ax,'|p| / fixed1m reference'); legend(ax,'Location','best'); title(ax,'Actual axis samples; true source remains masked');
ax=nexttile(t); phase=angle(F{2}.total(:,ix).*conj(F{3}.total(:,ix))); phase(abs(F{3}.total(:,ix))<.05*F{3}.envelope(:,ix))=NaN;
h=plot(ax,z,phase); h.UserData=struct('group',3,'kind','BC_phase'); grid(ax,'on'); xlim(ax,[0 5000]); xlabel(ax,'Physical depth (m)'); ylabel(ax,'arg(pB conj(pC)) (rad)'); title(ax,'Effective B/C phase; original gate0.05rad');
title(t,'A/B/C axis amplitude and B/C phase | no phase or gain fit','FontSize',16); pair(fig,'F03_ABC_axis');
receipt=struct('status','ABC_SAVED_FIELD_FIGURES','failed_checks',0,'elapsed_s',toc(start),'A_file',A_file,'B_file',B_file,'C_file',C_file,'output_dir',output_dir,'new_physics_solves',0,'native_grid',[10981 201]);
write(fullfile(output_dir,'checkpoint.json'),jsonencode(receipt,PrettyPrint=true));
html=['<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>A/B/C复场对照</title><style>body{font:16px/1.6 sans-serif;max-width:1800px;margin:24px auto;padding:0 20px}img{width:100%;height:auto}</style>' ...
'<h1>C组步长收敛：A/B/C完整复场</h1><p>A：Q1.25m、步长1m；B：Q0.625m、步长1m；C：Q0.625m、步长0.5m。近源图表原有步长同步减半。其余模型、N16、声束和原接收网格不变。源是理想3D点源，条带宽度不适用。D与新FFP未执行，Munk及底参数是假设，尚非论文最终物理验收。</p>' ...
'<h2>全场与中央区域</h2><img src="F01_ABC_full_central.png" alt="ABC全场和中央区域传播损失">' ...
'<h2>近底干涉</h2><img src="F02_ABC_near_bottom.png" alt="ABC原生0.01m近底接收网格传播损失">' ...
'<h2>轴线幅度与有效相位差</h2><img src="F03_ABC_axis.png" alt="ABC轴幅度及BC有效相位差">' ...
'<p>TL=-20log10|p|，固定1m参考，色标20至70dB仅用于显示。全深度以原5m与0.01m均匀分区分别绘制；无插值、平滑或坐标修改。真源点保持缺失，有效相位保留弱场缺口。同名FIG可编辑；数值结论以报告和CSV为准。</p></html>']; write(fullfile(output_dir,'index.html'),html); disp(jsonencode(receipt));
 function strips(ax,v,limits,g)
  hold(ax,'on'); rr={find(z<=4900),find(z>=4900)};
  for m=1:2, r=rr{m}; h=imagesc(ax,x,z(r),v(r,:)); h.UserData=struct('group',g,'rows',r); h.AlphaData=isfinite(v(r,:)); h.Interpolation='nearest'; end
  set(ax,'YDir','reverse','Color',[.9 .9 .9],'FontSize',12); xlim(ax,[-5 5]); ylim(ax,limits); xlabel(ax,'Physical x (m)'); ylabel(ax,'Physical depth (m)'); colormap(ax,parula(256));
 end
 function pair(fig,stem)
  drawnow; savefig(fig,fullfile(output_dir,[stem '.fig'])); exportgraphics(fig,fullfile(output_dir,[stem '.png']),'Resolution',140); close(fig);
 end
end
function write(path,value)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); c=onCleanup(@()fclose(fid)); fprintf(fid,'%s',value);
end
