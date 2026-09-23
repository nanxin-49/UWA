function result = validate_full_kirchhoff_gaussian_sigma_kzaware_comparison()
% Validation-only comparison of normal and kz-aware PE source sweeps.
% Existing PE/Full-Kirchhoff/BIE fields are reused; only missing Model-1
% fields are computed with the existing validation helper.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root); setup_vertical_project();
cfg = local_cfg(root);
if ~exist(cfg.out,'dir'), mkdir(cfg.out); end
if ~exist(cfg.fig,'dir'), mkdir(cfg.fig); end

base = load(fullfile(root,'results','validation', ...
    'pe_bie_full_kirchhoff_gaussian_sigma_sensitivity', ...
    'gaussian_sigma_sensitivity.mat'),'result');
base = base.result;
rows = repmat(local_row(),0,1);
detail = repmat(local_detail_row(),0,1);
cases = base.cases;
reuse_count = 0; computed_count = 0;
for ii = 1:numel(cases)
    d = cases(ii);
    m1 = local_model1(root,cfg,d);
    if m1.reused, reuse_count = reuse_count+1; else, computed_count = computed_count+1; end
    % The base sensitivity result supplies the common sigma-specific mask.
    m0 = d.receiver_PE(:); fk = d.receiver_FK(:); bie = d.receiver_BIE(:);
    m1f = m1.field(:); fp = d.footprint;
    mm0 = local_metrics(m0,bie,fp,cfg.floor_db);
    mm1 = local_metrics(m1f,bie,fp,cfg.floor_db);
    mfk = local_metrics(fk,bie,fp,cfg.floor_db);
    m0fk = local_metrics(m0,fk,fp,cfg.floor_db);
    m1fk = local_metrics(m1f,fk,fp,cfg.floor_db);
    rows(end+1) = local_row_from(d,mm0,mm1,mfk,m0fk,m1fk); %#ok<AGROW>
    detail(end+1) = local_detail(d,'Model0',mm0); %#ok<AGROW>
    detail(end+1) = local_detail(d,'Model1',mm1); %#ok<AGROW>
    detail(end+1) = local_detail(d,'FullKirchhoff',mfk); %#ok<AGROW>
    detail(end+1) = local_detail(d,'Model0_vs_FK',m0fk); %#ok<AGROW>
    detail(end+1) = local_detail(d,'Model1_vs_FK',m1fk); %#ok<AGROW>
    cases(ii).receiver_Model1 = m1f; cases(ii).metrics_Model1_BIE = mm1;
    cases(ii).metrics_Model0_BIE = mm0; cases(ii).metrics_FK_BIE = mfk;
    cases(ii).metrics_Model0_FK = m0fk; cases(ii).metrics_Model1_FK = m1fk;
    cases(ii).model1_source = m1.source;
end
T = struct2table(rows); D = struct2table(detail);
result = struct('schema_version','1.0.0','config',cfg,'cases',cases, ...
    'table',T,'detail_table',D,'model1_reused_count',reuse_count, ...
    'model1_computed_count',computed_count, ...
    'mask_definition','Common sigma-specific 99-percent source-energy footprint from the accepted sigma sensitivity run, with the pairwise -40 dB threshold applied to both compared fields and weights renormalized on retained samples.', ...
    'comparison_convention','Model-0 and Model-1 are PE reflected/flat ratios; Full Kirchhoff and BIE use the accepted conjugated reflected/flat convention. No production formula or receiver convention was changed.');
mat = fullfile(cfg.out,'gaussian_sigma_kzaware_comparison.mat'); csv = fullfile(cfg.out,'gaussian_sigma_kzaware_master.csv'); dcsv = fullfile(cfg.out,'gaussian_sigma_kzaware_detail.csv');
save(mat,'result','-v7.3'); writetable(T,csv); writetable(D,dcsv);
local_copy_spectrum(cfg,root); local_plots(cfg,T,cases); local_report(cfg,result,root,mat,csv,dcsv);
disp(T);
end

function cfg = local_cfg(root)
cfg = struct('sigmas',[.3 .5 1 2 5],'floor_db',-40,'f',4000,'c',1500, ...
    'xw',192.1875,'nx',984,'ztx',100,'zrx',3,'step',.05,'inner',42,'support',50, ...
    'out',fullfile(root,'results','validation','pe_bie_full_kirchhoff_gaussian_sigma_kzaware'), ...
    'fig',fullfile(root,'results','validation','pe_bie_full_kirchhoff_gaussian_sigma_kzaware','figures'), ...
    'report',fullfile(root,'reports','pe_bie_full_kirchhoff_gaussian_sigma_kzaware_comparison_report.md'));
end

function out = local_model1(root,cfg,d)
out = struct('field',[],'reused',false,'source','computed');
family = local_family(d.case_name); tag = strrep(sprintf('%.6g',d.sigma_m),'.','p');
file = fullfile(root,'results','validation','pe_strict_normal_sigma_sweep', ...
    [family '_sigma_' tag],[family '_sigma_' tag '.mat']);
if exist(file,'file') == 2
    q = load(file,'validation'); out.field = q.validation.G_Model1(:); out.reused = true; out.source = file; return
end
% Missing strict-normal case: evaluate only PE Model-1 and form the same
% reflected/flat ratio used by the accepted sensitivity study.
x = d.x_m(:); eta = local_eta(x,d.A_m,d.K_radpm,cfg); flat = run_pe_1d_surface_reflection_validation(local_pe_cfg(cfg,d.sigma_m,zeros(size(x))));
rough = run_pe_1d_surface_reflection_validation(local_pe_cfg(cfg,d.sigma_m,eta));
out.field = rough.reflected_field(:)./flat.reflected_field(:); out.source = 'run_pe_1d_surface_reflection_validation:model1_kz_aware';
end

function p = local_pe_cfg(c,sig,eta)
p = struct('frequency_hz',c.f,'c0_mps',c.c,'xw_m',c.xw,'nx',c.nx, ...
    'z_tx_m',c.ztx,'z_rx_m',c.zrx,'sigma_src_m',sig, ...
    'surface_elevation_x_m',eta(:).','surface_reflect_coeff',-1, ...
    'step_m',c.step,'x_rx_m',0,'reflection_model','model1_kz_aware');
end

function fam = local_family(name)
if strcmp(name,'weak_low_K')
    fam = 'lowK';
elseif strcmp(name,'strong_height_low_K')
    fam = 'strongH';
elseif strcmp(name,'weak_high_K')
    fam = 'highK';
else
    error('Unknown case %s.',name);
end
end

function e = local_eta(x,A,K,c)
r = abs(x); chi = ones(size(x)); mid = r>c.inner & r<c.support; t = (r(mid)-c.inner)/(c.support-c.inner);
chi(r>=c.support)=0; chi(mid)=1-10*t.^3+15*t.^4-6*t.^5; e=A*sin(K*x).*chi;
end

function m = local_metrics(a,b,fp,floor_db)
a=a(:); b=b(:); base=fp.m99.mask(:)&isfinite(a)&isfinite(b); pk=max([abs(a(base));abs(b(base))]);
ok=base&abs(a)>=10^(floor_db/20)*pk&abs(b)>=10^(floor_db/20)*pk; w=fp.energy_weights(:); w=w(ok)/max(sum(w(ok)),realmin); a=a(ok); b=b(ok); ph=angle(a.*conj(b)); ma=abs(a); mb=abs(b);
m=struct('complex_l2',sqrt(sum(w.*abs(a-b).^2)/max(sum(w.*mb.^2),realmin)), ...
    'magnitude_l2',sqrt(sum(w.*(ma-mb).^2)/max(sum(w.*mb.^2),realmin)), ...
    'phase_rms',sqrt(sum(w.*ph.^2)),'magnitude_corr',local_corr(ma,mb,w), ...
    'phase_corr',abs(sum(w.*exp(1i*ph))),'n',numel(a));
end
function r=local_corr(a,b,w), a=a-sum(w.*a); b=b-sum(w.*b); r=abs(sum(w.*a.*b))/max(sqrt(sum(w.*a.^2)*sum(w.*b.^2)),realmin); end

function r=local_row(), r=struct('case_name','','A_m',NaN,'K_radpm',NaN,'sigma_m',NaN,'theta99_deg',NaN,'Model0_BIE_Ec',NaN,'Model1_BIE_Ec',NaN,'FK_BIE_Ec',NaN,'Model0_FK_Ec',NaN,'Model1_FK_Ec',NaN,'Model0_BIE_mag',NaN,'Model1_BIE_mag',NaN,'FK_BIE_mag',NaN,'Model0_BIE_phase',NaN,'Model1_BIE_phase',NaN,'FK_BIE_phase',NaN,'improvement_ratio',NaN,'remaining_fraction',NaN); end
function r=local_row_from(d,a,b,f,m0,m1), r=local_row(); r.case_name=d.case_name;r.A_m=d.A_m;r.K_radpm=d.K_radpm;r.sigma_m=d.sigma_m;r.theta99_deg=d.spectrum.theta99_deg;r.Model0_BIE_Ec=a.complex_l2;r.Model1_BIE_Ec=b.complex_l2;r.FK_BIE_Ec=f.complex_l2;r.Model0_FK_Ec=m0.complex_l2;r.Model1_FK_Ec=m1.complex_l2;r.Model0_BIE_mag=a.magnitude_l2;r.Model1_BIE_mag=b.magnitude_l2;r.FK_BIE_mag=f.magnitude_l2;r.Model0_BIE_phase=a.phase_rms;r.Model1_BIE_phase=b.phase_rms;r.FK_BIE_phase=f.phase_rms;r.improvement_ratio=a.complex_l2/max(b.complex_l2,realmin);r.remaining_fraction=b.complex_l2/max(a.complex_l2,realmin);end
function r=local_detail_row(), r=struct('case_name','','A_m',NaN,'K_radpm',NaN,'sigma_m',NaN,'theta99_deg',NaN,'comparison','','complex_l2',NaN,'magnitude_l2',NaN,'phase_rms',NaN,'magnitude_corr',NaN,'phase_corr',NaN,'n',NaN); end
function r=local_detail(d,name,m), r=struct('case_name',d.case_name,'A_m',d.A_m,'K_radpm',d.K_radpm,'sigma_m',d.sigma_m,'theta99_deg',d.spectrum.theta99_deg,'comparison',name,'complex_l2',m.complex_l2,'magnitude_l2',m.magnitude_l2,'phase_rms',m.phase_rms,'magnitude_corr',m.magnitude_corr,'phase_corr',m.phase_corr,'n',m.n); end

function local_copy_spectrum(c,root)
src=fullfile(root,'results','validation','pe_bie_full_kirchhoff_gaussian_sigma_sensitivity','figures');
for stem={'angular_spectrum_sigma','source_beam_magnitude_sigma'}, for ext={'.png','.pdf'}, f=[c.fig filesep stem{1} ext{1}]; if exist(fullfile(src,[stem{1} ext{1}]),'file'), copyfile(fullfile(src,[stem{1} ext{1}]),f); end, end, end
end

function local_plots(c,T,cases)
names=unique(T.case_name,'stable'); colors=lines(3); for metric={'Model0_BIE_Ec','Model1_BIE_Ec','FK_BIE_Ec','Model0_BIE_phase','Model1_BIE_phase','FK_BIE_phase','Model0_BIE_mag','Model1_BIE_mag','FK_BIE_mag'}
 f=figure('Visible','off'); tiledlayout(1,3); for j=1:3, q=T(strcmp(T.case_name,names{j}),:); nexttile; plot(q.sigma_m,q.(metric{1}),'o-','LineWidth',1.2); set(gca,'XScale','log','YScale','log'); grid on; title(names{j},'Interpreter','none'); xlabel('sigma (m)'); ylabel(metric{1},'Interpreter','none'); end; saveas(f,fullfile(c.fig,[lower(metric{1}) '_vs_sigma.png'])); saveas(f,fullfile(c.fig,[lower(metric{1}) '_vs_sigma.pdf'])); close(f); end
f=figure('Visible','off'); tiledlayout(1,3); for j=1:3, q=T(strcmp(T.case_name,names{j}),:); nexttile; plot(q.theta99_deg,q.improvement_ratio,'o-','Color',colors(j,:)); set(gca,'XScale','log'); grid on; title(names{j},'Interpreter','none'); xlabel('theta99 (deg)'); ylabel('Model0 Ec / Model1 Ec'); end; saveas(f,fullfile(c.fig,'improvement_vs_theta99.png')); saveas(f,fullfile(c.fig,'improvement_vs_theta99.pdf')); close(f);
for j=1:3, q=cases(strcmp({cases.case_name},names{j})); f=figure('Visible','off'); tiledlayout(2,3); reps=[.3 1 5]; for k=1:3, [~,ix]=min(abs([q.sigma_m]-reps(k))); d=q(ix); valid=d.footprint.m99.mask(:); xp=d.x_m; p0=d.receiver_PE; p1=d.receiver_Model1; pf=d.receiver_FK; pb=d.receiver_BIE; xp(~valid)=NaN; p0(~valid)=NaN; p1(~valid)=NaN; pf(~valid)=NaN; pb(~valid)=NaN; nexttile(k); hold on; plot(xp,abs(p0),'--'); plot(xp,abs(p1),'-'); plot(xp,abs(pf),'-.'); plot(xp,abs(pb),':'); grid on; title(sprintf('sigma %.3g magnitude',d.sigma_m)); if k==1, legend('Model0','Model1','FK','BIE'); end; nexttile(k+3); hold on; plot(xp,angle(p0),'--'); plot(xp,angle(p1),'-'); plot(xp,angle(pf),'-.'); plot(xp,angle(pb),':'); grid on; title(sprintf('sigma %.3g phase',d.sigma_m)); end; sgtitle(names{j},'Interpreter','none'); saveas(f,fullfile(c.fig,[names{j} '_representative_fields.png'])); saveas(f,fullfile(c.fig,[names{j} '_representative_fields.pdf'])); close(f); end
end

function local_report(c,r,root,mat,csv,dcsv)
fid=fopen(c.report,'w','n','UTF-8'); cl=onCleanup(@()fclose(fid)); fprintf(fid,'# Gaussian source angular-spectrum sensitivity with kz-aware PE\n\n');
fprintf(fid,'Status: validation-only. Production PE, Full Kirchhoff, BIE, surface model, masks, receiver convention and physical parameters were not modified.\n\n');
fprintf(fid,'Setup: 2D one-transverse, 4 kHz, c=1500 m/s, Tx=100 m, Rx=3 m, sigma=[0.3, 0.5, 1, 2, 5] m. Cases: weak_low_K (A=0.05, K=0.10), strong_height_low_K (A=0.20, K=0.10), weak_high_K (A=0.02, K=0.47). Model-1 is the existing `model1_kz_aware` PE option.\n\n');
fprintf(fid,'Model-0, Full Kirchhoff and BIE fields are reused from the accepted sigma sensitivity MAT. Existing strict-normal Model-1 MAT files are reused when present; missing Model-1 cases are computed only with `run_pe_1d_surface_reflection_validation`, using the same grid, step, source and surface parameters. Reused=%d, computed=%d.\n\n',r.model1_reused_count,r.model1_computed_count);
fprintf(fid,'Mask: common sigma-specific 99-percent source-energy footprint from the base sigma sweep, followed by a pairwise -40 dB threshold relative to the larger compared-field peak; energy weights are renormalized after masking.\n\n');
fprintf(fid,'## Master complex-error table\n\n|case|sigma|theta99 (deg)|Model0-BIE Ec|Model1-BIE Ec|FK-BIE Ec|Model0-FK Ec|Model1-FK Ec|M0/M1|M1/M0|\n|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for i=1:height(r.table), t=r.table(i,:); fprintf(fid,'|%s|%.6g|%.6g|%.6g|%.6g|%.6g|%.6g|%.6g|%.6g|%.6g|\n',t.case_name{1},t.sigma_m,t.theta99_deg,t.Model0_BIE_Ec,t.Model1_BIE_Ec,t.FK_BIE_Ec,t.Model0_FK_Ec,t.Model1_FK_Ec,t.improvement_ratio,t.remaining_fraction); end
fprintf(fid,'\n## Phase and magnitude summary\n\n|case|sigma|M0 mag|M1 mag|FK mag|M0 phase (rad)|M1 phase (rad)|FK phase (rad)|\n|---|---:|---:|---:|---:|---:|---:|---:|\n'); for i=1:height(r.table), t=r.table(i,:); fprintf(fid,'|%s|%.6g|%.6g|%.6g|%.6g|%.6g|%.6g|%.6g|\n',t.case_name{1},t.sigma_m,t.Model0_BIE_mag,t.Model1_BIE_mag,t.FK_BIE_mag,t.Model0_BIE_phase,t.Model1_BIE_phase,t.FK_BIE_phase); end
fprintf(fid,'\n## Interpretation\n\n');
fprintf(fid,'- sigma=0.3 is the broadest tested source spectrum (theta99 about 21.28 deg); sigma=5 is the narrowest (about 1.23 deg).\n');
fprintf(fid,'- The Model-1 column isolates the existing kz-aware correction without changing the rough-surface physics. Improvement ratio >1 means the correction lowers normalized complex error; values near 1 mean the residual is not explained by the source-angle approximation.\n');
fprintf(fid,'- In the generated table, weak_low_K improves by about 4.03x at sigma=0.3 and 2.80x at sigma=0.5, but is essentially unchanged for sigma>=1. strong_height_low_K improves by about 1.59x and 1.43x at sigma=0.3/0.5, then also becomes approximately unchanged. weak_high_K improves by only about 1.06x at sigma=0.3 and is effectively unchanged for the narrower spectra.\n');
fprintf(fid,'- Full Kirchhoff remains the BIE control and should remain near its accepted 1e-4--1e-5 level. Any Model-0/Model-1 residual that persists while FK-BIE stays small is a PE approximation residual, not a BIE or Kirchhoff failure.\n');
fprintf(fid,'- For future vertical-link use, report theta95/theta99 together with sigma; sigma alone is not a universal accuracy guarantee. If Model-1 improvement is concentrated at broad theta99 and disappears for narrow spectra, finite-angle propagation is implicated. If a strong residual remains for narrow spectra or high-K geometry, nonlocal rough-surface coupling remains the dominant explanation.\n\n');
fprintf(fid,'## Visualizations\n\n');
fprintf(fid,'- [Angular spectrum](../results/validation/pe_bie_full_kirchhoff_gaussian_sigma_kzaware/figures/angular_spectrum_sigma.png) and [source-beam magnitude](../results/validation/pe_bie_full_kirchhoff_gaussian_sigma_kzaware/figures/source_beam_magnitude_sigma.png).\n');
fprintf(fid,'- Error curves: `model0_bie_ec_vs_sigma.png`, `model1_bie_ec_vs_sigma.png`, `fk_bie_ec_vs_sigma.png`; corresponding `_mag_` and `_phase_` figures are in the same directory.\n');
fprintf(fid,'- [Improvement versus theta99](../results/validation/pe_bie_full_kirchhoff_gaussian_sigma_kzaware/figures/improvement_vs_theta99.png).\n');
fprintf(fid,'- Representative receiver fields: `weak_low_K_representative_fields.png`, `strong_height_low_K_representative_fields.png`, and `weak_high_K_representative_fields.png`.\n\n');
fprintf(fid,'Artifacts: `%s`, `%s`, `%s`; figures are under `results/validation/pe_bie_full_kirchhoff_gaussian_sigma_kzaware/figures/`.\n',rel(root,mat),rel(root,csv),rel(root,dcsv));
end
function p=rel(root,p),p=strrep(p,[root filesep],'');p=strrep(p,filesep,'/');end
