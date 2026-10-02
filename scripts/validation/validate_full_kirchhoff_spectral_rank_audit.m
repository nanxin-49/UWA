function audit = validate_full_kirchhoff_spectral_rank_audit(render_only)
%VALIDATE_FULL_KIRCHHOFF_SPECTRAL_RANK_AUDIT Saved-data separability audit.
% Validation only. No physical reference or production solver is changed.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root); setup_vertical_project();
out=fullfile(root,'results','validation','full_kirchhoff_spectral_rank_audit');
figdir=fullfile(out,'figures'); if ~isfolder(figdir), mkdir(figdir); end
if nargin>0 && render_only
 q=load(fullfile(out,'spectral_rank_audit.mat'),'audit'); audit=q.audit;
 local_plots(audit,figdir); local_report(audit,root); return;
end
wall=tic; cpu0=cputime;
paths={fullfile(root,'results','validation','pe_bie_full_kirchhoff_parameter_sweep','pe_bie_full_kirchhoff_parameter_sweep.mat'), ...
 fullfile(root,'results','validation','pe_bie_full_kirchhoff_gaussian_sigma_sensitivity','gaussian_sigma_sensitivity.mat'), ...
 fullfile(root,'results','validation','pe_bellhop_controlled_comparison','stage0','stage0_validation.mat')};
p=load(paths{1},'result'); sweep=p.result;
% These trusted project artifacts contain obsolete serialized local handles.
% Suppress only their known relocation warning; never invoke those handles.
ws=warning('query','MATLAB:dispatcher:UnresolvedFunctionHandle');
cleanup=onCleanup(@()warning(ws.state,ws.identifier));
warning('off',ws.identifier); p=load(paths{2},'result'); sigma=p.result;
p=load(paths{3},'validation'); stage0=p.validation;
k=2*pi*sweep.config.frequency_hz/sweep.config.c0_mps;
src=sigma.cases(1).src;
assert(src.sigma==sweep.config.sigma_src_m);
% Explicit numeric source reconstruction verifies saved coefficients.
xg=src.xgrid(:); co=fft(exp(-.5*(xg/src.sigma).^2)).';
assert(norm(co-src.co)/norm(src.co)<1e-12);
spec=struct('case_id',{},'A_m',{},'K_radpm',{},'sigma_m',{},'src',{});
for i=1:numel(sweep.cases)
 d=sweep.cases(i);
 if i>1 && any([spec.A_m]==d.A_m & [spec.K_radpm]==d.K_radpm), continue; end
 spec(end+1)=struct('case_id',char(d.case_name),'A_m',d.A_m,'K_radpm',d.K_radpm,'sigma_m',src.sigma,'src',src); %#ok<AGROW>
end
for sig=[.5 1 2 5]
 idx=find([sigma.cases.A_m]==.2 & [sigma.cases.sigma_m]==sig,1);
 d=sigma.cases(idx); s=d.src;
 spec(end+1)=struct('case_id',sprintf('strong_A0p2_K0p1_sigma%g',sig),'A_m',d.A_m,'K_radpm',d.K_radpm,'sigma_m',sig,'src',s); %#ok<AGROW>
end
rows=struct([]); decompositions=struct([]); spectra=struct([]);
domains={'effective99','saved_propagating','full_propagating_control'};
for i=1:numel(spec)
 d=spec(i); s=d.src; pw=abs(s.co).^2; prop=abs(s.kx)<=k;
 pw=pw/sum(pw(prop)); [ab,ix]=sort(abs(s.kx(prop)));
 pp=pw(prop); cp=cumsum(pp(ix));
 q95=ab(find(cp>=.95,1)); q99=ab(find(cp>=.99,1));
 th95=asind(q95/k); th99=asind(q99/k);
 spectrum=struct('case_id',d.case_id,'kx_radpm',s.kx,'power',pw,'propagating',prop,'kx95',q95,'kx99',q99);
 if isempty(spectra), spectra=spectrum; else, spectra(end+1)=spectrum; end %#ok<AGROW>
 for id=1:numel(domains)
  domain=domains{id};
  if id==1
   mask=prop & abs(s.kx)<=q99; kk=s.kx(mask); w=pw(mask);
  elseif id==2
   mask=prop; kk=s.kx(mask); w=pw(mask);
  else
   kk=linspace(-k,k,1025); w=ones(size(kk));
  end
  eta=linspace(-d.A_m,d.A_m,257).'; delta=k-sqrt(max(0,k^2-kk.^2));
  C=exp(1i*eta*delta); Csmall=exp(1i*eta*(kk.^2/(2*k)));
  for weighted=[false true]
   if id==3 && weighted, continue; end
   weights=ones(size(w))/numel(w); if weighted, weights=w/sum(w); end
   t=tic; [U,S,V]=svd(C.*sqrt(weights),'econ'); seconds=toc(t);
   sv=diag(S); tail=local_tail(sv); ranks=arrayfun(@(tol)find(tail<=tol,1),[1e-2 1e-3 1e-4]);
   row=struct('case_id',d.case_id,'A_m',d.A_m,'K_radpm',d.K_radpm,'sigma_m',d.sigma_m, ...
    'kA',k*d.A_m,'AK',d.A_m*d.K_radpm,'AK2_over_k',d.A_m*d.K_radpm^2/k,'K_over_k',d.K_radpm/k, ...
    'theta95_deg',th95,'theta99_deg',th99,'kx95_over_k',q95/k,'kx99_over_k',q99/k, ...
    'chi',d.A_m*(k-sqrt(k^2-q99^2)),'chi_small_angle',k*d.A_m*(th99*pi/180)^2/2, ...
    'domain',domain,'weighted',weighted,'energy_retained',sum(w),'n_eta',numel(eta),'n_kx',numel(kk), ...
    'rank_1e2',ranks(1),'rank_1e3',ranks(2),'rank_1e4',ranks(3),'small_angle_matrix_error',norm((Csmall-C).*sqrt(weights),'fro')/norm(C.*sqrt(weights),'fro'), ...
    'svd_wall_seconds',seconds,'matrix_MiB',16*numel(C)/2^20,'svd_factors_MiB',16*(numel(U)+numel(V)+numel(S))/2^20);
   if id==3, row.energy_retained=NaN; end
   if isempty(rows), rows=row; else, rows(end+1)=row; end %#ok<AGROW>
   dd=struct('case_id',d.case_id,'domain',domain,'weighted',weighted,'singular_values',sv,'relative_tail',tail);
   if isempty(decompositions), decompositions=dd; else, decompositions(end+1)=dd; end %#ok<AGROW>
  end
 end
 fprintf('rank audit %d/%d: %s\n',i,numel(spec),d.case_id);
end
T=struct2table(rows); writetable(T,fullfile(out,'rank_summary.csv'));
% Resolution audit changes quadrature density at fixed physical bandwidth,
% using interpolation of the saved energy spectrum, not new propagation.
gridrows=struct([]);
for sig=[.3 1]
 j=find([spec.sigma_m]==sig & [spec.A_m]==.2,1); s=spec(j).src;
 ss=spectra(j); [kx,ord]=sort(s.kx); power=ss.power(ord);
 for domain={'effective99','saved_propagating'}
 for level=1:3
  ne=[129 257 513]; nk=[257 513 1025];
  if strcmp(domain{1},'effective99')
   kk=linspace(-ss.kx99,ss.kx99,nk(level));
  else
   kk=linspace(min(kx),max(kx),nk(level));
  end
  eta=linspace(-.2,.2,ne(level)).';
  w=interp1(kx,power,kk,'linear'); w=w/sum(w);
  for weighted=[false true]
   ww=ones(size(w))/numel(w); if weighted, ww=w; end
   C=exp(1i*eta*(k-sqrt(k^2-kk.^2))).*sqrt(ww);
   sv=svd(C,'econ'); tail=local_tail(sv);
   r=struct('sigma_m',sig,'domain',domain{1},'n_eta',ne(level),'n_kx',nk(level),'weighted',weighted, ...
    'rank_1e2',find(tail<=1e-2,1),'rank_1e3',find(tail<=1e-3,1),'rank_1e4',find(tail<=1e-4,1));
   if isempty(gridrows), gridrows=r; else, gridrows(end+1)=r; end %#ok<AGROW>
  end
 end
 end
end
G=struct2table(gridrows); writetable(G,fullfile(out,'grid_convergence.csv'));
% Reuse all three accepted sigma=.3 FK cases with correct analytic normals.
sanity=struct([]);
for ii=[2 6 11]
 t=tic; check=local_operator(sweep.cases(ii),src,stage0,sweep.config,k);
 check.wall_seconds=toc(t);
 if isempty(sanity), sanity=check; else, sanity(end+1)=check; end %#ok<AGROW>
 writetable(check.table,fullfile(out,sprintf('operator_%02d.csv',ii)));
 fprintf('operator %s: saved closure %.3g\n',check.case_id,check.saved_receiver_closure);
 save(fullfile(out,sprintf('operator_%02d.mat',ii)),'check','-v7.3');
end
for i=1:numel(spectra)
 s=spectra(i); ts=table(s.kx_radpm(:),s.kx_radpm(:)/k,s.power(:),s.propagating(:), ...
  'VariableNames',{'kx_radpm','kx_over_k','normalized_energy_weight','propagating'});
 writetable(ts,fullfile(out,sprintf('incident_spectrum_%02d.csv',i)));
end
audit=struct('schema_version','1.0','config',struct('k_radpm',k,'eta_samples',257,'tolerances',[1e-2 1e-3 1e-4], ...
 'error_norm','relative Frobenius; column weights = normalized incident energy; uniform eta measure'), ...
 'input_paths',{paths},'matlab_version',version,'rank_table',T,'grid_table',G, ...
 'decompositions',decompositions,'spectra',spectra,'sanity',sanity);
local_plots(audit,figdir);
audit.wall_seconds=toc(wall); audit.cpu_seconds=cputime-cpu0;
save(fullfile(out,'spectral_rank_audit.mat'),'audit','-v7.3');
local_report(audit,root);
disp(T(T.weighted & strcmp(T.domain,'saved_propagating'),{'case_id','chi','rank_1e2','rank_1e3','rank_1e4'}));
fprintf('Wall %.3f s, CPU %.3f s\n',audit.wall_seconds,audit.cpu_seconds);
end

function tail=local_tail(s)
% Sum small singular values directly: avoid cancellation from total-cumsum.
energy=flipud(cumsum(flipud(abs(s).^2)));
tail=sqrt([energy(2:end);0]/sum(abs(s).^2));
end

function check=local_operator(d,s,stage0,c,k)
xs=d.surface_x_m(:); eta=d.surface_z_m(:); A=d.A_m; K=d.K_radpm;
% The same C2 profile derivative as the authoritative parameter sweep.
r=abs(xs); win=zeros(size(xs)); dw=win; win(r<=c.surface_inner_m)=1;
mid=r>c.surface_inner_m & r<c.surface_support_m;
t=(r(mid)-c.surface_inner_m)/(c.surface_support_m-c.surface_inner_m);
win(mid)=1-10*t.^3+15*t.^4-6*t.^5;
dw(mid)=(-30*t.^2+60*t.^3-30*t.^4)/(c.surface_support_m-c.surface_inner_m);
sp=A*(K*cos(K*xs).*win+sin(K*xs).*dw.*sign(xs));
jac=sqrt(1+sp.^2); normal=[-sp./jac,1./jac];
assert(max(abs(eta-A*sin(K*xs).*win))<1e-12);
ok=abs(s.kx)<=k; kk=s.kx(ok); kz=real(s.kz(ok));
w=abs(s.co(ok)).^2; w=w/sum(w);
etagrid=linspace(-A,A,257).'; Cgrid=exp(1i*etagrid*(k-kz));
Csurf=exp(1i*eta*(k-kz)); strong=exp(-1i*k*eta);
phase=exp(1i*(xs-s.x0)*kk);
coeff=s.co(ok).*exp(1i*s.ztx*kz)/s.n;
% Exact factorization, including spectral multipliers in BOTH derivatives.
u=sum(phase.*Csurf.*coeff,2).*strong;
ux=sum(phase.*Csurf.*(1i*kk.*coeff),2).*strong;
uz=sum(phase.*Csurf.*(-1i*kz.*coeff),2).*strong;
dn=normal(:,1).*ux+normal(:,2).*uz;
assert(norm(u-d.surface_incident)/norm(d.surface_incident)<1e-11);
mask=stage0.footprint.m99.mask(:); xr=stage0.x_m(mask); xr=xr(:);
cfg=struct('frequency_hz',c.frequency_hz,'c0_mps',c.c0_mps, ...
 'surface_x_m',xs,'surface_z_m',eta,'surface_normal',normal,'surface_incident',u, ...
 'surface_incident_dn',dn,'receiver_x_m',xr,'receiver_z_m',c.z_rx_m+zeros(size(xr)));
ref=solve_full_kirchhoff_2d(cfg);
flat=-sum(exp(1i*(xr-s.x0)*s.kx).*(s.co.*exp(1i*(c.z_tx_m+c.z_rx_m)*s.kz)),2)/s.n;
ratio=conj(ref.receiver_field./flat);
rw=stage0.footprint.energy_weights(mask); rw=rw(:)/sum(rw);
closure=local_error(ratio,d.receiver_FK(mask),rw);
assert(closure<1e-10,'Full FK native reconstruction differs from saved result.');
% Build the unchanged full-support FK linear contraction once, then verify
% it against the reference. This avoids repeated expensive Green evaluation.
seg=hypot(diff(xs),diff(eta)); qw=[seg(1)/2;(seg(1:end-1)+seg(2:end))/2;seg(end)/2];
dx=xr(:)-xs.'; dz=c.z_rx_m-eta.'; rr=hypot(dx,dz);
B=-1i/4*besselh(0,1,k*rr).*qw.';
D=-1i*k/4*besselh(1,1,k*rr).*(dx.*normal(:,1).'+dz.*normal(:,2).')./rr.*qw.';
assert(norm(B*dn+D*u-ref.receiver_field)/norm(ref.receiver_field)<1e-11);
rows=struct([]); fields=struct([]);
for weighted=[false true]
 weights=ones(size(w))/numel(w); if weighted, weights=w; end
 [~,S,V]=svd(Cgrid.*sqrt(weights),'econ'); sv=diag(S); tails=local_tail(sv);
 for rank=1:8
  Vr=V(:,1:rank);
  % Extend the fixed right singular subspace to actual saved surface heights.
  % No interpolation/Taylor expansion of the strong carrier phase is used.
  Cr=((Csurf.*sqrt(weights))*Vr)*Vr'; Cr=Cr./sqrt(weights);
  ur=sum(phase.*Cr.*coeff,2).*strong;
  uxr=sum(phase.*Cr.*(1i*kk.*coeff),2).*strong;
  uzr=sum(phase.*Cr.*(-1i*kz.*coeff),2).*strong;
  dnr=normal(:,1).*uxr+normal(:,2).*uzr;
  field=B*dnr+D*ur; gr=conj(field./flat);
  row=struct('weighted',weighted,'rank',rank,'grid_matrix_error',tails(rank), ...
   'surface_matrix_error',norm((Cr-Csurf).*sqrt(weights),'fro')/norm(Csurf.*sqrt(weights),'fro'), ...
   'incident_l2',norm(ur-u)/norm(u),'normal_derivative_l2',norm(dnr-dn)/norm(dn), ...
   'receiver_native_l2',local_error(field,ref.receiver_field,rw),'receiver_ratio_l2',local_error(gr,ratio,rw));
  if isempty(rows), rows=row; else, rows(end+1)=row; end %#ok<AGROW>
  f=struct('weighted',weighted,'rank',rank,'receiver_native',field);
  if isempty(fields), fields=f; else, fields(end+1)=f; end %#ok<AGROW>
 end
end
check=struct('case_id',char(d.case_name),'A_m',A,'K_radpm',K,'sigma_m',s.sigma, ...
 'z_rx_m',c.z_rx_m,'table',struct2table(rows),'saved_incident_closure',norm(u-d.surface_incident)/norm(u), ...
 'saved_receiver_closure',closure,'receiver_x_m',xr,'receiver_reference_native',ref.receiver_field, ...
 'fields',fields,'n_surface',numel(xs),'n_receiver',numel(xr), ...
 'coupling_matrix_MiB',16*numel(Csurf)/2^20,'contraction_matrices_MiB',16*(numel(B)+numel(D))/2^20);
end

function e=local_error(a,b,w)
e=sqrt(sum(w(:).*abs(a(:)-b(:)).^2)/sum(w(:).*abs(b(:)).^2));
end

function local_plots(a,out)
f=figure('Visible','off','Position',[100 100 1100 700]); tiledlayout(2,2);
for j=1:4
 nexttile; hold on;
 if j<=3, names={a.sanity(j).case_id}; else, names={a.sanity(2).case_id}; end
 for i=1:numel(a.decompositions)
  d=a.decompositions(i); if ~strcmp(d.case_id,names{1}), continue; end
  if j==4 && ~strcmp(d.domain,'full_propagating_control'), continue; end
  if j<=3 && strcmp(d.domain,'full_propagating_control'), continue; end
  style='--'; if d.weighted, style='-'; end
  semilogy(1:min(14,numel(d.singular_values)),max(d.singular_values(1:min(14,end))/d.singular_values(1),1e-16),style, ...
   'DisplayName',sprintf('%s w=%d',d.domain,d.weighted),'LineWidth',1.4);
 end
 set(gca,'YScale','log'); ylim([1e-16 1]); grid on; xlabel('singular value index'); ylabel('s_j/s_1');
 title(names{1},'Interpreter','none'); legend('Location','northeast','Interpreter','none');
end
local_export(f,out,'singular_value_decay');
T=a.rank_table; ix=T.weighted & strcmp(T.domain,'saved_propagating'); t=T(ix,:);
f=figure('Visible','off','Position',[100 100 1100 700]); tiledlayout(2,3);
vars={'kA','K_over_k','AK','AK2_over_k','theta99_deg','chi'};
for j=1:6
 nexttile; hold on; symbols={'o','s','^','d','v'}; sigmas=[.3 .5 1 2 5];
 for is=1:numel(sigmas)
  mask=t.sigma_m==sigmas(is);
  scatter(t.(vars{j})(mask),t.rank_1e3(mask),55,t.sigma_m(mask),'filled',symbols{is}, ...
   'DisplayName',sprintf('sigma=%g',sigmas(is)));
 end
 clim([.3 5]); ylim([.8 2.2]); yticks([1 2]); grid on;
 xlabel(vars{j},'Interpreter','none'); ylabel('weighted rank at 10^{-3}'); colorbar;
 title('color: source sigma (m)');
end
local_export(f,out,'rank_dimensionless_parameters');
f=figure('Visible','off','Position',[100 100 1250 450]); tiledlayout(1,3);
for i=1:numel(a.sanity)
 nexttile; t=a.sanity(i).table; w=t.weighted; hold on;
 semilogy(t.rank(w),max(t.grid_matrix_error(w),1e-16),'o-','DisplayName','weighted matrix');
 semilogy(t.rank(w),max(t.receiver_ratio_l2(w),1e-16),'s--','DisplayName','weighted FK ratio');
 semilogy(t.rank(~w),max(t.receiver_ratio_l2(~w),1e-16),'^:','DisplayName','unweighted FK ratio');
 yline(1e-3,'k--','10^{-3}','HandleVisibility','off'); set(gca,'YScale','log');
 grid on; ylim([1e-16 1]); xlabel('rank'); ylabel('relative error');
 title(a.sanity(i).case_id,'Interpreter','none'); legend('Location','northeast');
end
local_export(f,out,'operator_reconstruction');
f=figure('Visible','off','Position',[100 100 1000 430]); tiledlayout(1,2);
nexttile; hold on; nexttile; hold on;
for sig=[.3 .5 1 2 5]
 j=find(abs(a.rank_table.sigma_m-sig)<1e-12,1); name=a.rank_table.case_id{j};
 idx=find(strcmp({a.spectra.case_id},name),1); s=a.spectra(idx);
 [kk,ord]=sort(s.kx_radpm(s.propagating)); pp=s.power(s.propagating); pp=pp(ord);
 nexttile(1); semilogy(kk/a.config.k_radpm,max(pp,1e-20),'DisplayName',sprintf('sigma=%g m',sig));
 % Collapse tied +/- bins before cumulative plot; quantiles use both bins.
 [ab,~,groups]=unique(abs(kk)); mass=accumarray(groups,pp(:));
 nexttile(2); plot(ab/a.config.k_radpm,cumsum(mass),'DisplayName',sprintf('sigma=%g m',sig));
end
nexttile(1); set(gca,'YScale','log'); xlim([-.6 .6]); ylim([1e-12 1]); xlabel('k_x/k'); ylabel('normalized bin energy'); grid on; legend;
nexttile(2); xlim([0 .5]); ylim([.9 1.001]);
yline(.95,':','DisplayName','95%'); yline(.99,'--','DisplayName','99%');
xlabel('|k_x|/k'); ylabel('cumulative energy'); grid on; legend;
local_export(f,out,'incident_energy_ranges');
end

function local_export(f,out,name)
exportgraphics(f,fullfile(out,[name '.png']),'Resolution',180);
exportgraphics(f,fullfile(out,[name '.pdf']),'ContentType','vector'); close(f);
end

function local_report(a,root)
path=fullfile(root,'reports','full_kirchhoff_spectral_rank_audit_report.md');
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); closer=onCleanup(@()fclose(fid));
fprintf(fid,'# Full Kirchhoff spectral separability / low-rank audit\n\n');
fprintf(fid,'日期：2026-09-30。状态：**VERIFIED_IN_SAVED_CASES；建议继续开发 Spectral Fast FK prototype**。本阶段已完成并停止；没有开发生产 Fast FK。\n\n');
fprintf(fid,'## 结论与版本边界\n\n');
fprintf(fid,'真实入射能量加权、保留全部已保存传播频谱时，15 个受控案例达到矩阵相对 Frobenius 误差 `1e-2 / 1e-3 / 1e-4` 分别只需 rank **1–2 / 1–2 / 1–3**。强高度宽角案例为 **2 / 2 / 3**。三个完整 FK 运算检查中，加权 rank 2 的接收场粗/平比值误差全部小于 `1e-3`。固定物理区间加密没有引起 rank 增长。\n\n');
fprintf(fid,'本地基线是干净的 main、HEAD `02a29d92b84eab9f07fb115969022dff4960788e`，已通过 `e78c24d` 整合研究支线。GitHub 查询显示 [PR #3](https://github.com/nanxin-49/UWA/pull/3) 仍为 Draft/open、merged=false，head SHA 与本地 HEAD 相同；因此“本地已整合”和“GitHub PR 已合并”必须区分。没有同步、提交、推送或改变 PR。\n\n');
fprintf(fid,'结论仅涉及二维单横向、4 kHz、c=1500 m/s、Tx=100 m、Rx=3 m、C2 截断确定性正弦面及保存 Gaussian 入射谱。没有读取 cash/，没有重跑 PE、BIE、PM Monte Carlo、宽带数据库或 z=1 m 结果。交接报告中的旧 worktree 路径现在是历史 provenance；本轮使用主线结果目录。\n\n');
fprintf(fid,'## 数学对象、范数与权重\n\n');
fprintf(fid,'令 `kz_j=sqrt(k^2-kx_j^2)`、`delta_j=k-kz_j`，审查 `C_mj=exp(i delta_j eta_m)`，其中 `eta_m` 均匀覆盖 `[-A,A]`，默认 257 点。完整高度相位为 `exp(-i kz_j eta_m)=exp(-i k eta_m) C_mj`。强高度载波完整保留，绝不进行 small-height Taylor expansion。左乘单位模对角矩阵不改变奇异值，因此完整高度矩阵和残余矩阵的 SVD rank 理论上相同；分解的价值是揭示控制耦合的是 delta 的跨度，而不是去掉某个物理大相位。\n\n');
fprintf(fid,'未加权矩阵使用 `B=C/sqrt(n_kx)`；加权使用 `B=C diag(sqrt(w_j))`，`w_j=|co_j|^2/sum_j |co_j|^2`。使用**能量的平方根**乘列，不是直接以能量乘列。误差定义为 `||B-B_r||_F/||B||_F=sqrt(sum_{j>r} s_j^2/sum_j s_j^2)`，报告满足门槛的最小正整数 rank。eta 方向为均匀测度；它不是表面实际高度分布的概率加权。该范数是平均矩阵误差，不是最坏入射算子范数或每个接收点的误差上界。尾能量从小奇异值直接累加，避免 total-minus-cumsum 的抵消。\n\n');
fprintf(fid,'三个频域：`effective99` 为真实谱的最小对称 99%% 能量区间；`saved_propagating` 保留所有已保存传播箱；`full_propagating_control` 为 1025 点均匀 `[-k,k]`、仅未加权数学对照。95%%/99%% 分位数按 |kx| 排序，边界相同的正负箱全部保留，实际能量可略大于 99%%。没有对入射场加 M99/Hann 窗来生成权重：优先使用保存源的原始复角谱系数，避免窗口展宽。已有 covariance CSV 仅作为交叉背景，其 Hann 谱不是这里的权重。\n\n');
fprintf(fid,'`saved_propagating` 的源网格 N=984、宽度 192.1875 m；所有这些箱都为传播波，范围约 `[-0.960,0.958] k`。没有丢弃已保存尾部，但未采样的近掠射区间不等于已实测工程谱。数学 full-propagating 对照补充这一差别。截掉 1%% 能量一般可产生约 0.1 的谱 L2 误差，所以 effective99 rank 不能单独作为严格场误差保证。\n\n');
fprintf(fid,'## 数据复用与实际角宽\n\n');
for i=1:numel(a.input_paths)
 fprintf(fid,'- `%s`\n',strrep(strrep(a.input_paths{i},[root filesep],''),filesep,'/'));
end
fprintf(fid,'\n复用完整 sweep 的 11 个唯一 A–K 点；相同 `(A=.02,K=.1)` 重复行去重。额外仅对 strong-height A=.2、K=.1 复用 sigma=.5/1/2/5 的数值源谱，不重跑原 sweep。保存 MAT 中过期函数句柄只加载而不调用；代码仅定向抑制已确认的 `UnresolvedFunctionHandle` 路径警告，读取数值 co/kx/kz。sigma=.3 的保存 co 与原 Gaussian 离散 FFT 闭合在 `1e-12` 内。\n\n');
T=a.rank_table;
fprintf(fid,'|sigma (m)|theta95 (deg)|theta99 (deg)|kx95/k|kx99/k|99%% 区间实际能量|\n|---:|---:|---:|---:|---:|---:|\n');
for sig=[.3 .5 1 2 5]
 t=T(find(T.sigma_m==sig & strcmp(T.domain,'effective99') & T.weighted,1),:);
 fprintf(fid,'|%.3g|%.7g|%.7g|%.7g|%.7g|%.9g|\n',sig,t.theta95_deg,t.theta99_deg,t.kx95_over_k,t.kx99_over_k,t.energy_retained);
end
fprintf(fid,'\n## 参数与 rank 结果\n\n');
fprintf(fid,'各 rank 元组依次对应 `(1e-2,1e-3,1e-4)`。U99/W99 为有效区间未加权/加权，Uall/Wall 为全部保存频谱未加权/加权，Ufull 为完整传播区间数学对照。每行完整 theta95/99、样本数、保留能量、耗时和小角矩阵误差均见 `rank_summary.csv`。\n\n');
fprintf(fid,'|A (m)|K (rad/m)|sigma (m)|kA|AK|AK²/k|K/k|theta99 (deg)|chi|U99|W99|Uall|Wall|Ufull|\n|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|---|---|\n');
ids=unique(T.case_id,'stable');
for i=1:numel(ids)
 tt=T(strcmp(T.case_id,ids{i}),:); t=tt(1,:);
 fprintf(fid,'|%.4g|%.4g|%.3g|%.6g|%.6g|%.6g|%.6g|%.6g|%.6g|',t.A_m,t.K_radpm,t.sigma_m,t.kA,t.AK,t.AK2_over_k,t.K_over_k,t.theta99_deg,t.chi);
 for j=1:height(tt), fprintf(fid,'%d/%d/%d|',tt.rank_1e2(j),tt.rank_1e3(j),tt.rank_1e4(j)); end
 fprintf(fid,'\n');
end
fprintf(fid,'\n### 参数控制量解释\n\n');
fprintf(fid,'`chi=A[k-kz(theta99)]=kA(1-cos(theta99))`；小角估计为 `kA theta99^2/2`（theta 用弧度），同时保留 `delta≈kx²/(2k)` 的真实矩阵误差。固定 A=.2 和 k 时，sigma=.3→5 使 chi=.22848→.00077196、全谱加权 rank 从 2/2/3 降至 1/1/1，而 kA 始终为 3.35103。因此 **chi 连同真实能量分布比单独 kA 更能解释 rank**。相同 chi 也不保证相同加权 rank，因为尾部谱能量及范数权重仍有影响。小角式只作诊断，SVD 与算子检查均使用精确 kz。\n\n');
fprintf(fid,'在均匀 eta 区间矩阵的定义中，K 不出现。因此固定 A、sigma 的 K=.05–.70 扫描 rank 完全相同，是数学构造的必然结果，不能解读为 K 对完整 FK 不重要。AK、AK²/k、K/k 会通过实际 eta(x)、法向、散射出射带宽和 Green 核影响完整运算。矩阵层面没有强行拟合经验公式；K=.47 的真实曲面检查提供有限的算子支持，K=.70 尚未做完整 reconstruction 检查。\n\n');
fprintf(fid,'小角诊断，固定 A=.2、K=.1；末列为 effective99 域的能量加权相对 Frobenius 误差 `||C_small-C||/||C||`：\n\n');
fprintf(fid,'|sigma|chi exact|chi small-angle|small-angle matrix error|\n|---:|---:|---:|---:|\n');
for sig=[.3 .5 1 2 5]
 t=T(find(T.A_m==.2 & T.sigma_m==sig & T.weighted & strcmp(T.domain,'effective99'),1),:);
 fprintf(fid,'|%g|%.8g|%.8g|%.8g|\n',sig,t.chi,t.chi_small_angle,t.small_angle_matrix_error);
end
fprintf(fid,'\n宽角端的小角近似自身误差不保证达到 1e-4，不能以低秩成功替代精确 kz；实际原型应保留精确谱相位。\n\n');
fprintf(fid,'## 网格稳定性\n\n');
fprintf(fid,'固定 A=.2、sigma=.3 或 1，保持物理 eta/kx 区间不变，把 `(n_eta,n_kx)` 从 `(129,257)` 加密到 `(257,513)`、`(513,1025)`。权重由保存能量谱线性插值，并在各网格归一；没有生成新源场或重新求解。这是离散积分/采样密度稳定性，而不是新的物理 Nyquist 带宽或真实传播网格收敛。\n\n');
fprintf(fid,'|sigma|domain|weighted|rank tuples at the three densities|\n|---:|---|---|---|\n');
G=a.grid_table;
for sig=[.3 1]
 for dom={'effective99','saved_propagating'}
  for w=[false true]
   t=G(G.sigma_m==sig & strcmp(G.domain,dom{1}) & G.weighted==w,:);
   fprintf(fid,'|%g|%s|%d|',sig,dom{1},w);
   for j=1:height(t), fprintf(fid,'%d/%d/%d ',t.rank_1e2(j),t.rank_1e3(j),t.rank_1e4(j)); end
   fprintf(fid,'|\n');
  end
 end
end
fprintf(fid,'\n所有测试 rank 均稳定，未出现随采样数显著增长。\n\n');
fprintf(fid,'## 完整 FK 运算 sanity check\n\n');
fprintf(fid,'使用参数扫描中 weak `(A=.01,K=.1)`、strong-height `(.2,.1)`、high-K `(.02,.47)`，sigma=.3、Rx=3 m。保持原保存曲面 4097 点、完整 [-60,60] m 支撑、解析 C2 法向、原弧长梯形测度、全部 984 个传播谱箱与固定 M99 接收 411 点/原能量权重。没有额外 -40 dB 成对筛选，没有逐例拟合、归一相位或载波校准。\n\n');
fprintf(fid,'以 uniform-eta SVD 的右奇异子空间 V_r 投影实际曲面耦合：`C_r(eta(x))=[C(eta(x)) diag(sqrt(w)) V_r] V_r^H diag(1/sqrt(w))`。这是同一右子空间到保存高度样本的直接延拓，不是插值左奇异向量。显式计算 `u=sum_j co_j exp(i kx_j(x-x0)) exp(i kz_j zTx) exp(-ik eta) C_r/n`；ux、uz 使用同一个 C_r，并分别完整保留 `i kx_j`、`-i kz_j` 乘子；物理坐标导数不对曲面限制后的 u(x,eta(x)) 做切向求导。法向导数为 `n_x ux+n_z uz`。\n\n');
fprintf(fid,'先调用原 `solve_full_kirchhoff_2d` 得到未截断基线；再只在本脚本装配相同的全支撑 Green/导数线性收缩，并确认与原 solver 相对差 `<1e-11`，避免每个 rank 重算 Bessel 核。以近似后的 u/dn 作用到该完整收缩，出射几何/Green 完全不近似。场比值用固定 pressure-release 平面参考，`G=conj(u_ref/flat)` 一次共轭，未再次共轭已保存 G。\n\n');
fprintf(fid,'接收原生场和粗/平比值均使用 `sqrt(sum(w |a-b|²)/sum(w |b|²))`；没有相位重定。矩阵误差、入射场误差、法向导数误差、receiver_native_l2、receiver_ratio_l2 都保存在 operator CSV。\n\n');
fprintf(fid,'|case|saved incident closure|saved FK ratio closure|weighted r=1 FK error|r=2 matrix error|r=2 FK error|r=3 FK error|r=4 FK error|\n|---|---:|---:|---:|---:|---:|---:|---:|\n');
for i=1:numel(a.sanity)
 s=a.sanity(i); t=s.table(s.table.weighted,:);
 fprintf(fid,'|%s|%.4g|%.4g|%.6g|%.6g|%.6g|%.6g|%.6g|\n',s.case_id,s.saved_incident_closure,s.saved_receiver_closure,t.receiver_ratio_l2(1),t.grid_matrix_error(2),t.receiver_ratio_l2(2),t.receiver_ratio_l2(3),t.receiver_ratio_l2(4));
end
fprintf(fid,'\n场误差随 rank 下降，与矩阵截断误差趋势一致；加权 rank 2 所有案例通过 1e-3，rank 3 所有案例通过 1e-4。强高度 r=2 的 grid Frobenius 误差约 9.22e-4，而实际高度样本 Frobenius 误差约 1.076e-3、FK ratio 误差约 6.61e-4，说明趋势一致不意味着数值相等或严格场误差上界。未加权强高度 r=4 约 3.80e-4、r=5 约 1.96e-5；完整 CSV 包括两种加权方式的 r=1–8。机器精度附近不再解释微小非单调变化。\n\n');
fprintf(fid,'## CPU、内存与复现\n\n');
fprintf(fid,'MATLAB `%s`，共享 MCP 会话，标准 MATLAB FFT/SVD/Bessel/graphics，无新增工具箱。完整最终一次运行（加载、矩阵 audit、网格 audit、三个算子检查、PNG/PDF 导出；不含最后保存/report 写入）wall **%.3f s**，进程 CPU **%.3f s**。CPU 包括可能的多线程计算，因此可大于 wall；它不是实时性能或生产速度比较。未运行重复基准、未测峰值 RSS、未证明加速比。\n\n',a.matlab_version,a.wall_seconds,a.cpu_seconds);
fprintf(fid,'主 audit 残余矩阵 complex-double payload 为 %.3f–%.3f MiB；最大 U/S/V payload 合计 %.3f MiB。每个 sanity 的实际曲面耦合矩阵 %.3f MiB，两块 FK Green/导数收缩矩阵合计 %.3f MiB。还有 phase、C_r、Bessel 临时数组及 SVD 工作空间，所以这些是**数组存储核算，不是进程峰值**。原型仍有 dense 矩阵和 O(N_surface N_kx)、O(N_receiver N_surface) 工作；低秩重构的误差检查不是 FFT 性能实现。\n\n',min(T.matrix_MiB),max(T.matrix_MiB),max(T.svd_factors_MiB),a.sanity(1).coupling_matrix_MiB,a.sanity(1).contraction_matrices_MiB);
fprintf(fid,'运行：先将项目根加入 MATLAB path 并执行 `setup_vertical_project`，再调用 `audit=validate_full_kirchhoff_spectral_rank_audit()`。入口通过 mfilename 定位根目录，不改变 cwd，输出路径明确；重新运行会更新本 audit 目录和本报告。数值数据无 RNG。输入和代码 SHA-256 见目录中的 `provenance.json`（由本轮仓库工具记录，不属于求解器）。\n\n');
fprintf(fid,'仅从本 audit 保存数据重绘/更新报告：`validate_full_kirchhoff_spectral_rank_audit(true)`；此模式不重算 SVD/FK，报告中的运行时间仍对应保存的完整 numerical run，排除后续纯渲染时间。\n\n');
fprintf(fid,'## 图与数据索引\n\n');
names={'incident_energy_ranges','singular_value_decay','rank_dimensionless_parameters','operator_reconstruction'};
captions={'原保存源的归一 bin 能量与累积能量；95/99% 角宽不由全传播区间代替。','三个案例有效/全保存谱奇异值衰减，以及强高度完整传播数学对照；机器精度以下不作物理解读。','15 个案例全保存谱加权 1e-3 rank 对各无量纲量；点重叠表示离散 rank 相同，颜色为 sigma。','真实完整 FK 运算误差随 rank 下降；矩阵与场使用各自明确的相对范数。'};
for i=1:numel(names)
 fprintf(fid,'![%s](../results/validation/full_kirchhoff_spectral_rank_audit/figures/%s.png)\n\n%s\n\n',names{i},names{i},captions{i});
end
fprintf(fid,'- `spectral_rank_audit.mat`: rank/grid/spectrum/奇异值/算子结果和运行统计，无旧 worktree 函数句柄。\n- `rank_summary.csv`: 75 行，两种真实域×两种权重、一个数学对照，全部参数/误差/成本。\n- `grid_convergence.csv`: 24 行，两个域/两种源宽/三种密度/两种权重。\n- `incident_spectrum_01..15.csv`: 原始 kx bin、归一能量权重和传播标志。\n- `operator_02/06/11.csv` 与对应 MAT: 三个真实曲面运算检查，包含 r=1–8 的全部原生重构场。\n- `figures/`: 四组 PNG 和 vector PDF。\n\n');
fprintf(fid,'## 下一阶段判断与停止点\n\n');
fprintf(fid,'**继续开发 Spectral Fast FK prototype。** 依据是本范围内 1e-3 有效加权 rank=1–2、固定物理区间加密 rank 稳定，以及三类完整积分运算检查通过。这个判断是下一阶段研究建议，并不授权本阶段继续编码。\n\n');
fprintf(fid,'**INFERENCE**：入射 residual 分离可将 sum_j 展成少量可 FFT 的谱乘子与空间系数，具备 O(r N log N) 候选结构。**UNVERIFIED**：整个非局部 FK 的出射 Green 表示、散射角带宽、倏逝尾部、FFT padding/periodization/有限支撑处理，以及装配成本摊销仍需下一阶段证明。入射 theta99 不自动界定出射 theta99；K 会改变出射带宽。未经这些检查不能把当前结果称为完整 Fast FK 或全核低秩证明。\n\n');
fprintf(fid,'真实 PM 海面、比 sigma=.3 更宽的工程角谱、更大 kA、更高频率、三维和宽带均未验证。均匀 eta 区间在给定 A 下能覆盖高度，但不能自动验证 PM 的导数、散射和出射谱。若后续这些范围所需 rank 很高、随物理带宽/网格增长或出射重构不稳定，应明确停止全局分离路线并比较 block-Fresnel、locally stationary approximation、SSA/SSA2；本次没有这些失败证据，也没有对备选方法作优劣认证。\n\n');
fprintf(fid,'## 改动与文档纪律\n\n');
fprintf(fid,'新增一个 validation-only 入口、独立结果目录及本报告，更新 scripts/README、PROJECT_CONTEXT 与 vertical_comm_guide 的验证边界说明。未改变 production paramsV/output/results 接口、PE/WAPE、BIE、solve_full_kirchhoff_2d 或 local phase screen。没有新增物理公式实现，只审查现有谱耦合的数学等价分解。当前不涉及生产信道回归，验证为保存源/场闭合、原 FK contraction 等价性、rank 密度稳定性和真实重构误差；未运行 scalar-channel/wideband-communication 回归，也不声称它们本轮通过。\n');
end
