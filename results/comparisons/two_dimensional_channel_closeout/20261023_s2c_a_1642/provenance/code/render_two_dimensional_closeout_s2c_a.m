function render_two_dimensional_closeout_s2c_a(mat_path,out)
% Read-only rendering of saved geometry and real complex receiver fields.
assert(~isfolder(out),'Use a fresh rendering directory'); mkdir(out);
S=load(mat_path,'plot_data'); d=S.plot_data; x=d.x(1:411); U=d.U;
col=[0 .38 .70;.82 .32 .08;.15 .55 .28]; labs={'PE','BIE','FullFK complete'};
f=newfig(1120,1080); t=tiledlayout(f,4,1,TileSpacing='compact',Padding='compact');
c=4; g=d.geometry(d.geometry.case_name==d.names(c),:);
a=nexttile(t); plot(a,g.x_m,g.eta_m,'k',LineWidth=1.4); ylabel(a,'\eta [m]'); title(a,'M00 强档总览 | A=0.20 m, K=0.1 rad/m, \phi=0, z_{Rx}=3 m');
a=nexttile(t); hold(a,'on');
for m=1:3, plot(a,x,abs(U(1:411,c,m)),Color=col(m,:),LineWidth=1.3); end
ylabel(a,'|u_{rough}| [unit source]'); legend(a,labs,Location='eastoutside');
a=nexttile(t); hold(a,'on');
for m=[1 3], phaseplot(a,x,angle(U(1:411,c,m).*conj(U(1:411,c,2))),col(m,:)); end
ylabel(a,'wrapped \Delta\phi [rad]'); legend(a,{'PE-BIE','FullFK-BIE'},Location='eastoutside');
a=nexttile(t); hold(a,'on'); a.YScale='log';
for m=[1 3], semilogy(a,x,max(abs(U(1:411,c,m)-U(1:411,c,2)),realmin),Color=col(m,:),LineWidth=1.3); end
ylabel(a,'|u_m-u_{BIE}| [unit source]'); xlabel(a,'receiver x [m]');
title(a,'原生复差；参考近零=0；全线两差异已分辨；对齐x不表示局部因果');
finish(f,t,[-40.1 40.1]); store(f,out,'M00_overview');
f=newfig(1150,1020); t=tiledlayout(f,3,1,TileSpacing='compact',Padding='compact');
vars={'eta_m','slope','curvature_per_m'}; yl={'\eta [m]','d\eta/dx [m/m]','d^2\eta/dx^2 [1/m]'};
for j=1:3
 a=nexttile(t); hold(a,'on');
 for c=2:4
  g=d.geometry(d.geometry.case_name==d.names(c),:); plot(a,g.x_m,g.(vars{j}),Color=col(c-1,:),LineWidth=1.3);
 end
 ylabel(a,yl{j}); xline(a,-42,':',HandleVisibility='off'); xline(a,42,':',HandleVisibility='off'); xline(a,-50,'--',HandleVisibility='off'); xline(a,50,'--',HandleVisibility='off');
 if j==1, title(a,'M01 C2正弦几何 | K=0.1 rad/m, \phi=0 | 平台\leq42 m，过渡42–50 m'); legend(a,{'A=.01 m','A=.05 m','A=.20 m'},Location='eastoutside'); end
 if j==2, title(a,'x=0高度/曲率为0，坡度分别为0.001 / 0.005 / 0.020；并非平海面'); end
 if j==3
  h=d.illumination(d.illumination.case_name=="flat",:); yyaxis(a,'right'); plot(a,h.x_m,abs(complex(h.incident_real,h.incident_imag)).^2/max(abs(complex(h.incident_real,h.incident_imag)).^2),'k--'); ylabel(a,'保存平面入射强度/max [1]'); ylim(a,[0 1.1]); yyaxis(a,'left'); xlabel(a,'surface x [m]');
  title(a,'入射强度仅用保存GL节点；不映射为接收x的局部贡献');
 end
end
finish(f,t,[-55 55]); store(f,out,'M01_geometry_all');
for c=2:4
 stems=["weak","medium","strong"]; amplitudes=[.01 .05 .2]; stem=char(stems(c-1));
 r=U(1:411,c,2); flat=U(1:411,1,2); du=r-flat;
 f=newfig(1240,1000); t=tiledlayout(f,2,2,TileSpacing='compact',Padding='compact');
 a=nexttile(t); plot(a,x,abs(flat),'k--',x,abs(r),Color=col(2,:),LineWidth=1.3); ylabel(a,'|u_{BIE}| [unit source]'); legend(a,{'flat','rough'}); title(a,'数值BIE参考的海面响应；不是模型误差');
 a=nexttile(t); hold(a,'on'); phaseplot(a,x,angle(flat),[.3 .3 .3]); phaseplot(a,x,angle(r),col(2,:)); ylabel(a,'arg(u_{BIE}) wrapped [rad]'); legend(a,{'flat','rough'}); title(a,'绝对相位在\pm\pi处断线');
 a=nexttile(t); hold(a,'on');
 for m=1:3, plot(a,x,abs(U(1:411,c,m)),Color=col(m,:),LineWidth=1.3); end
 ylabel(a,'|u_{rough}| [unit source]'); legend(a,labs); title(a,'三方法原生幅度；没有rough/flat归一化');
 a=nexttile(t); plot(a,x,abs(du),'k',LineWidth=1.3); ylabel(a,'|u_{BIE,rough}-u_{BIE,flat}| [unit source]'); title(a,'物理增量：相消与相长共同决定，非PE/FK误差');
 title(t,sprintf('M02a %s | A=%.2f m | 411真实z=3 m样本 | native exp(-i omega t)',stem,amplitudes(c-1)),Interpreter='tex');
 finish(f,t,[-40.1 40.1]); store(f,out,['M02a_' stem '_response']);
 f=newfig(1230,1250); t=tiledlayout(f,3,2,TileSpacing='compact',Padding='compact');
 for p=1:6
  a=nexttile(t); hold(a,'on');
  for m=[1 3]
   u=U(1:411,c,m); dp=angle(u.*conj(r)); aa=(abs(u)-abs(r)).^2; pp=4*abs(u).*abs(r).*sin(dp/2).^2;
   switch p
    case 1, plot(a,x,20*log10(abs(u)./abs(r)),Color=col(m,:));
    case 2, phaseplot(a,x,dp,col(m,:));
    case 3, semilogy(a,x,max(abs(u-r),realmin),Color=col(m,:));
    case 4, semilogy(a,x,max(aa,realmin),Color=col(m,:));
    case 5, semilogy(a,x,max(pp,realmin),Color=col(m,:));
    case 6, semilogy(a,x,max(aa+pp,realmin),Color=col(m,:));
   end
  end
  ys={'20log_{10}(|u_m|/|u_{BIE}|) [dB]','wrapped \Delta\phi [rad]','|e_m| [unit source]','(|u_m|-|u_{BIE}|)^2 [unit source^2]','4|u_m||u_{BIE}|sin^2(\Delta\phi/2) [unit source^2]','amplitude^2 + phase^2 = |e_m|^2'};
  if p>=3, a.YScale='log'; end
  ylabel(a,ys{p}); legend(a,{'PE-BIE','FullFK-BIE'},Location='best');
  if p==1, title(a,'相对BIE的幅度差'); elseif p==2, title(a,'相位差不作常数拟合'); elseif p==3, title(a,'绝对复残差；微小值使用独立纵轴'); elseif p==4, title(a,'精确幅度平方项'); elseif p==5, title(a,'精确相位平方项'); else, title(a,'平方恒等式；显示下限realmin仅供对数轴'); end
 end
 title(t,sprintf('M02b %s | 模型误差与精确平方分解 | 参考近零=0',stem));
 finish(f,t,[-40.1 40.1]); store(f,out,['M02b_' stem '_errors']);
 f=newfig(1260,1470); t=tiledlayout(f,4,2,TileSpacing='compact',Padding='compact');
 for k=1:2
  a=nexttile(t); hold(a,'on'); ii=(1:411)'; if k==2, ii=d.local; end
  for m=[1 3]
   semilogy(a,x(ii),max(abs(U(ii,c,m)-U(ii,c,2)),realmin),Color=col(m,:),LineWidth=1.3);
   jj=d.ix9(d.grid.z_m==3); scatter(a,d.x(jj),abs(U(jj,c,m)-U(jj,c,2)),35,col(m,:),'o');
  end
  a.YScale='log';
  ylabel(a,'|e| [unit source]'); xlabel(a,'receiver x [m]'); xline(a,0,':',HandleVisibility='off'); grid(a,'on');
  if k==1, xlim(a,[-40.1 40.1]); title(a,'全线411点；q=1/411'); else, xlim(a,[-1 1]); title(a,sprintf('中心|x|<=1 m：%d真样本；圆点仅P21/P22/P23',numel(ii))); end
  legend(a,{'PE-BIE','PE z3 grid','FullFK-BIE','FK z3 grid'},Location='best');
  if k==1
   for m=[1 3]
    rs=d.resolutions(d.resolutions.case_name==d.names(c)&d.resolutions.candidate==d.methods(m)&d.resolutions.reference=="BIE"&d.resolutions.support_id=="M99_z3_411",:);
    yline(a,5*rs.numerical_resolution_estimate*sqrt(mean(abs(U(1:411,c,2)).^2)),'--',Color=col(m,:),DisplayName=char(d.methods(m)+" 5x support RMS floor"));
   end
  end
 end
 for p=1:3
  for m=[1 3]
   a=nexttile(t); ii=d.ix9; u=U(ii,c,m); r9=U(ii,c,2);
   if p==1
    v=20*log10(abs(u)./abs(r9)); cap='amplitude difference [dB]';
   elseif p==2
    v=angle(u.*conj(r9)); cap='wrapped phase difference [rad]';
   else
    v=abs(u-r9)./abs(r9); cap='|e| / |u_{BIE}| [1]';
   end
   imagesc(a,1:3,1:3,reshape(v,3,3)'); a.YDir='reverse';
   allu=d.U(d.ix9,2:4,m); allr=d.U(d.ix9,2:4,2);
   if p==1
    scale=max(abs(20*log10(abs(allu(:))./abs(allr(:))))); clim(a,[-scale scale]);
   elseif p==2
    scale=max(abs(angle(allu(:).*conj(allr(:))))); clim(a,[-scale scale]);
   else
    scale=max(abs(allu(:)-allr(:))./abs(allr(:))); clim(a,[0 scale]);
   end
   a.XTick=1:3; a.XTickLabel=compose('%.6g',unique(d.grid.x_m));
   a.YTick=1:3; a.YTickLabel=compose('%.3f',unique(d.grid.z_m)); xlabel(a,'real receiver x [m]'); ylabel(a,'depth z [m]');
   cb=colorbar(a); cb.Label.String=cap; colormap(a,parula(128));
   for j=1:9, text(a,mod(j-1,3)+1,ceil(j/3),sprintf('%s\n%.2g',d.grid.receiver_id(j),v(j)),HorizontalAlignment='center',FontSize=10,Color='k',BackgroundColor='w',Margin=.5); end
   rs=d.resolutions(d.resolutions.case_name==d.names(c)&d.resolutions.candidate==d.methods(m)&d.resolutions.reference=="BIE"&d.resolutions.support_id=="P11-P33",:);
   if p==1, title(a,sprintf('%s 9点 E=%.3g | %s',d.methods(m),rs.relative_difference,resolve_label(rs.resolved_above_5x_estimate))); else, title(a,sprintf('%s | %s',d.methods(m),cap)); end
  end
 end
 title(t,sprintf('M03 %s | 全线/中心/九点独立支持；九点q=1/9；6个离深度点从未投到线上',stem));
 style(f); store(f,out,['M03_' stem '_supports']);
end
end

function s=resolve_label(flag)
if flag, s='resolved >5x empirical'; else, s='UNRESOLVED (9-point support)'; end
end
function f=newfig(w,h)
f=figure(Visible='off',Color='w',Units='pixels',Position=[40 40 w h]);
end
function phaseplot(a,x,y,col)
% Insert a NaN separator WITHOUT deleting either real sample at a wrap.
jump=abs(diff(y))>pi; pos=(1:numel(x))'+cumsum([0;jump(:)]);
xx=nan(numel(x)+sum(jump),1); yy=xx; xx(pos)=x; yy(pos)=y;
plot(a,xx,yy,Color=col,LineWidth=1.1);
end
function finish(f,t,xlimits)
style(f); aa=findall(f,'Type','axes');
for k=1:numel(aa), xlim(aa(k),xlimits); xline(aa(k),0,':',Color=[.5 .5 .5],HandleVisibility='off'); grid(aa(k),'on'); if isempty(aa(k).XLabel.String), xlabel(aa(k),'x [m]'); end, end
t.TileSpacing='compact';
end
function style(f)
aa=findall(f,'Type','axes'); for k=1:numel(aa), aa(k).FontName='Microsoft YaHei'; aa(k).FontSize=11; aa(k).Box='on'; end
tt=findall(f,'Type','text'); for k=1:numel(tt), tt(k).FontName='Microsoft YaHei'; end
end
function store(f,out,stem)
drawnow; savefig(f,fullfile(out,[stem '.fig'])); exportgraphics(f,fullfile(out,[stem '.png']),Resolution=150); close(f);
end
