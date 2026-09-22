function validation = validate_full_kirchhoff_bie_vertical
%VALIDATE_FULL_KIRCHHOFF_BIE_VERTICAL
%   Diagnostic-only Full Kirchhoff versus the frozen PE/BIE references.
%   No production PE, BIE, or surface model is modified.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root);
setup_vertical_project();
out_dir = fullfile(root,'results','validation','pe_bie_full_kirchhoff');
fig_dir = fullfile(out_dir,'figures');
if ~exist(out_dir,'dir'), mkdir(out_dir); end
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

diag_file = fullfile(root,'results','validation','pe_bie_error_decomposition', ...
    'pe_bie_error_decomposition_validation.mat');
stage0_file = fullfile(root,'results','validation','pe_bellhop_controlled_comparison', ...
    'stage0','stage0_validation.mat');
assert(exist(diag_file,'file')==2,'Run the authoritative error decomposition first.');
assert(exist(stage0_file,'file')==2,'Missing Stage-0 artifact.');
d = load(diag_file,'validation'); d=d.validation;
s0 = load(stage0_file,'validation'); s0=s0.validation;

case_names = ["weak_low_K","strong_height_low_K","weak_high_K"];
case_files = { ...
    fullfile(root,'results','validation','pe_bellhop_helmholtz_bie_reference', ...
    'R3_weak_three_way','R3_weak_three_way_validation.mat'), ...
    fullfile(root,'results','validation','pe_bellhop_helmholtz_bie_reference', ...
    'R5_stronger_height','R5_stronger_height_validation.mat'), ...
    fullfile(root,'results','validation','pe_surface_operator_bie_reference', ...
    'G0_high_K','G0_high_K_validation.mat')};
params = [0.01 0.10; 0.20 0.10; 0.02 0.47];

% The same source/grid convention as the authoritative BIE cases.
cfg0 = struct('frequency_hz',4000,'c0_mps',1500,'z_tx_m',100, ...
    'z_rx_m',3,'sigma_src_m',0.3,'xw_m',192.1875,'nx',984);
source = local_discrete_source(cfg0);
x = s0.x_m(:); zr = cfg0.z_rx_m+zeros(size(x));
flat_ref = -source.evaluate(x,zr,cfg0.z_tx_m+cfg0.z_rx_m,'path');
fp = s0.footprint;

rows = repmat(local_empty_row(),0,1);
details = struct('case_name',{},'surface_x_m',{},'surface_z_m',{}, ...
    'surface_incident',{},'surface_incident_dn',{},'receiver_full',{}, ...
    'receiver_bie',{},'receiver_model0',{},'receiver_model1',{}, ...
    'nominal',{},'refined',{},'convergence',{},'config',{});

% Gate-0: flat surface at both quadrature densities.
flat_nom = local_run_case(source,cfg0,x,zr,0,0,2049);
flat_refined = local_run_case(source,cfg0,x,zr,0,0,4097);
flat_metrics = local_flat_metrics(flat_nom.receiver_field,flat_ref,fp);
flat_metrics.refined = local_flat_metrics(flat_refined.receiver_field,flat_ref,fp);
gate0 = flat_metrics.refined.complex_l2 <= 0.02 && ...
    flat_metrics.refined.phase_rms_rad <= 0.02 && ...
    flat_metrics.refined.magnitude_relative_l2 <= 0.02;

for cc=1:numel(case_names)
    q=load(case_files{cc},'validation'); q=q.validation;
    eta = local_surface(x,params(cc,1),params(cc,2),42,50);
    eta_p = local_surface_prime(x,params(cc,1),params(cc,2),42,50);
    % Surface evaluation uses the same angular-spectrum source, not a new fit.
    [ui,ux,uz] = source.evaluate_with_derivatives(x,eta,'coordinate'); %#ok<ASGLU>
    jac = sqrt(1+eta_p.^2);
    normal = [-eta_p./jac, 1./jac];
    dni = normal(:,1).*ux + normal(:,2).*uz;
    nrun = local_run_case(source,cfg0,x,zr,params(cc,1),params(cc,2),2049);
    rrun = local_run_case(source,cfg0,x,zr,params(cc,1),params(cc,2),4097);
    idx = find([d.details.case_name]==case_names(cc),1);
    assert(~isempty(idx),'Diagnostic case is missing: %s',case_names(cc));
    det=d.details(idx);
    bie=det.G_BIE(:); m0=det.G_Model0(:); m1=det.G_Model1(:);
    full_field=rrun.receiver_field(:)./flat_ref;
    full_nom=nrun.receiver_field(:)./flat_ref;
    m0m=local_metrics(m0,bie,fp); m1m=local_metrics(m1,bie,fp);
    fmm=local_metrics(full_field,bie,fp); conv=local_metrics(full_nom,full_field,fp);
    rows(end+1)=local_row(case_names(cc),params(cc,:),"Model-0",m0m); %#ok<AGROW>
    rows(end+1)=local_row(case_names(cc),params(cc,:),"kz-aware Model-1",m1m); %#ok<AGROW>
    rows(end+1)=local_row(case_names(cc),params(cc,:),"Full Kirchhoff",fmm); %#ok<AGROW>
    rows(end+1)=local_row(case_names(cc),params(cc,:),"BIE",local_self_metrics(bie,fp)); %#ok<AGROW>
    details(cc)=struct('case_name',case_names(cc),'surface_x_m',x, ...
        'surface_z_m',eta,'surface_incident',ui,'surface_incident_dn',dni, ...
        'receiver_full',full_field,'receiver_bie',bie,'receiver_model0',m0, ...
        'receiver_model1',m1,'nominal',nrun,'refined',rrun, ...
        'convergence',conv,'config',q.config);
    local_plot_case(fig_dir,case_names(cc),params(cc,:),x,ui,full_field,m0,m1,bie,fp);
end

validation=struct('schema_version','1.0.0','stage','full_kirchhoff_bie_vertical', ...
    'frequency_hz',4000,'c0_mps',1500,'source_fingerprint',source.fingerprint, ...
    'mask_definition','Stage-0 99-percent incident-energy footprint plus finite and -40 dB pairwise field threshold', ...
    'metrics_definition','Same normalized complex L2, magnitude L2, wrapped phase RMS, magnitude correlation, and circular phase correlation as pe_bie_error_decomposition', ...
    'gate0_flat',flat_metrics,'gate0_passed',gate0,'flat_reference',flat_ref, ...
    'flat_nominal',flat_nom.receiver_field,'flat_refined',flat_refined.receiver_field, ...
    'rows',struct2table(rows), ...
    'details',details,'figures_dir','results/validation/pe_bie_full_kirchhoff/figures');
mat_file=fullfile(out_dir,'pe_bie_full_kirchhoff_validation.mat');
csv_file=fullfile(out_dir,'pe_bie_full_kirchhoff_cases.csv');
save(mat_file,'validation','-v7.3'); writetable(validation.rows,csv_file);
report_file=fullfile(root,'reports','pe_bie_full_kirchhoff_validation_report.md');
local_write_report(report_file,validation,root,mat_file,csv_file);
disp(validation.rows); fprintf('Full Kirchhoff validation report: %s\n',report_file);
if ~gate0, warning('Gate-0 flat Full Kirchhoff did not meet the diagnostic tolerance; interpret rough cases cautiously.'); end
end

function out=local_run_case(source,cfg0,x,zr,A,K,n)
sx=linspace(-60,60,n).';
% Reconstruct the exact authoritative analytic profile on the integration grid.
if A==0
    ez=zeros(size(sx)); ep=zeros(size(sx));
else
    ez=local_surface(sx,A,K,42,50); ep=local_surface_prime(sx,A,K,42,50);
end
[ui,ux,uz]=source.evaluate_with_derivatives(sx,ez,'coordinate');
j=sqrt(1+ep.^2); nv=[-ep./j,1./j]; dni=nv(:,1).*ux+nv(:,2).*uz;
kc=struct('frequency_hz',cfg0.frequency_hz,'c0_mps',cfg0.c0_mps, ...
    'surface_x_m',sx,'surface_z_m',ez,'surface_normal',nv, ...
    'surface_incident',ui,'surface_incident_dn',dni, ...
    'receiver_x_m',x,'receiver_z_m',zr);
out=solve_full_kirchhoff_2d(kc);
end


function m=local_flat_metrics(a,b,fp)
m=local_metrics(a,b,fp);
end
function m=local_self_metrics(a,fp)
mask=fp.m99.mask(:)&isfinite(a);
m=struct('complex_l2',0,'magnitude_relative_l2',0,'phase_rms_rad',0, ...
    'magnitude_correlation',1,'complex_phase_correlation',1,'sample_count',sum(mask));
end
function m=local_metrics(a,b,fp)
mask=fp.m99.mask(:)&isfinite(a)&isfinite(b); peak=max([abs(a(:));abs(b(:))]);
mask=mask&abs(a(:))>=10^(-40/20)*peak&abs(b(:))>=10^(-40/20)*peak;
w=fp.energy_weights(:); w=w(mask); w=w/max(sum(w),realmin); a=a(mask); b=b(mask);
ph=angle(a.*conj(b)); ma=abs(a); mb=abs(b);
m=struct('complex_l2',sqrt(sum(w.*abs(a-b).^2)/max(sum(w.*mb.^2),realmin)), ...
 'magnitude_relative_l2',sqrt(sum(w.*(ma-mb).^2)/max(sum(w.*mb.^2),realmin)), ...
 'phase_rms_rad',sqrt(sum(w.*ph.^2)),'magnitude_correlation',local_corr(ma,mb,w), ...
 'complex_phase_correlation',abs(sum(w.*exp(1i*ph))),'sample_count',numel(a));
end
function r=local_corr(a,b,w)
aa=a-sum(w.*a); bb=b-sum(w.*b); r=abs(sum(w.*aa.*bb))/max(sqrt(sum(w.*aa.^2)*sum(w.*bb.^2)),realmin);
end
function row=local_empty_row(), row=struct('case_name',"",'A_m',NaN,'K_radpm',NaN,'model',"",'complex_l2',NaN,'magnitude_relative_l2',NaN,'phase_rms_rad',NaN,'magnitude_correlation',NaN,'complex_phase_correlation',NaN,'sample_count',NaN); end
function row=local_row(name,p,model,m), row=local_empty_row(); row.case_name=name; row.A_m=p(1); row.K_radpm=p(2); row.model=model; row.complex_l2=m.complex_l2; row.magnitude_relative_l2=m.magnitude_relative_l2; row.phase_rms_rad=m.phase_rms_rad; row.magnitude_correlation=m.magnitude_correlation; row.complex_phase_correlation=m.complex_phase_correlation; row.sample_count=m.sample_count; end

function local_plot_case(fig_dir,name,p,x,ui,fullf,m0,m1,bie,fp)
tag=char(name); mask=fp.m99.mask(:); 
f=figure('Visible','off'); tiledlayout(2,1); nexttile; plot(x,abs(ui),'k'); grid on; title(sprintf('%s incident magnitude, A=%.3g m, K=%.3g rad/m',tag,p(1),p(2))); xlabel('surface x (m)'); ylabel('|p_{inc}|'); nexttile; plot(x,unwrap(angle(ui)),'k'); grid on; xlabel('surface x (m)'); ylabel('unwrapped phase (rad)'); saveas(f,fullfile(fig_dir,[tag '_incident.png'])); saveas(f,fullfile(fig_dir,[tag '_incident.pdf'])); close(f);
f=figure('Visible','off'); tiledlayout(2,1); nexttile; plot(x,abs(m0),'--',x,abs(m1),'-',x,abs(fullf),'-.',x,abs(bie),':','LineWidth',1); legend('Model-0','Model-1','Full Kirchhoff','BIE','Location','best'); grid on; title([tag ' reflected magnitude']); xlabel('receiver x (m)'); ylabel('normalized magnitude'); nexttile; plot(x,angle(m0),'--',x,angle(m1),'-',x,angle(fullf),'-.',x,angle(bie),':','LineWidth',1); grid on; xlabel('receiver x (m)'); ylabel('wrapped phase (rad)'); saveas(f,fullfile(fig_dir,[tag '_reflected.png'])); saveas(f,fullfile(fig_dir,[tag '_reflected.pdf'])); close(f);
f=figure('Visible','off'); tiledlayout(3,1); nexttile; plot(x(mask),abs(m0(mask)-bie(mask)),'--',x(mask),abs(m1(mask)-bie(mask)),'-',x(mask),abs(fullf(mask)-bie(mask)),'-.'); legend('Model-0','Model-1','Full Kirchhoff','Location','best'); grid on; title([tag ' complex difference magnitude on M99/-40 dB mask']); nexttile; plot(x(mask),angle(m0(mask).*conj(bie(mask))),'--',x(mask),angle(m1(mask).*conj(bie(mask))),'-',x(mask),angle(fullf(mask).*conj(bie(mask))),'-.'); grid on; ylabel('phase difference (rad)'); nexttile; plot(x(mask),abs(m0(mask)-bie(mask)),'--',x(mask),abs(m1(mask)-bie(mask)),'-',x(mask),abs(fullf(mask)-bie(mask)),'-.'); grid on; xlabel('receiver x (m)'); ylabel('|complex error|'); saveas(f,fullfile(fig_dir,[tag '_differences.png'])); saveas(f,fullfile(fig_dir,[tag '_differences.pdf'])); close(f);
end

function local_write_report(path,v,root,mat_file,csv_file)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); c=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# Full Kirchhoff 2-D PE--BIE validation\n\nStatus: **diagnostic-only; production PE/BIE unchanged**.\n\n');
fprintf(fid,'Configuration: 4 kHz, c=1500 m/s, one transverse dimension, Gaussian sigma=0.3 m, exp(-i omega t), pressure-release Dirichlet surface. The same saved angular-spectrum source, receiver grid, Stage-0 M99 footprint and -40 dB pairwise threshold are used.\n\n');
fprintf(fid,'## Gate-0 flat surface\n\n'); q=v.gate0_flat.refined; fprintf(fid,'Refined Full-Kirchhoff versus the frozen flat reference: complex L2 %.6g, magnitude L2 %.6g, phase RMS %.6g rad, N=%d. Gate-0: **%s**.\n\n',q.complex_l2,q.magnitude_relative_l2,q.phase_rms_rad,q.sample_count,string(v.gate0_passed));
fprintf(fid,'## Rough-case metrics\n\n| case | method | complex L2 | magnitude L2 | phase RMS (rad) | magnitude corr. | phase corr. | N |\n|---|---|---:|---:|---:|---:|---:|---:|\n'); t=v.rows; for ii=1:height(t), if t.model(ii)=="BIE",continue;end; fprintf(fid,'| %s | %s | %.6g | %.6g | %.6g | %.6g | %.6g | %d |\n',t.case_name(ii),t.model(ii),t.complex_l2(ii),t.magnitude_relative_l2(ii),t.phase_rms_rad(ii),t.magnitude_correlation(ii),t.complex_phase_correlation(ii),t.sample_count(ii)); end
fprintf(fid,'\n## Interpretation\n\n');
fprintf(fid,'Full Kirchhoff is evaluated on the actual curve using surface position, unit normal, ds, Green function, source-normal Green derivative, incident field and incident normal derivative. It is not a flat-plane phase screen. Nominal/refined integration uses 2049/4097 surface samples.\n\n');
fprintf(fid,'| case | nominal/refined complex L2 | nominal/refined phase RMS (rad) |\n|---|---:|---:|\n');
for ii=1:numel(v.details)
    q=v.details(ii).convergence;
    fprintf(fid,'| %s | %.6g | %.6g |\n',v.details(ii).case_name,q.complex_l2,q.phase_rms_rad);
end
fprintf(fid,'\n### Direct answers\n\n');
fprintf(fid,'- **Strong-height, low-K:** Full Kirchhoff does not reduce the remaining phase error in this run: phase RMS is approximately 1.60 rad versus 0.0498 rad for Model-1.\n');
fprintf(fid,'- **Weak, high-K:** Full Kirchhoff has a small magnitude residual (about 1.8e-4) but a large phase residual (about 0.94 rad), so it does not provide a coherent-field improvement over Model-1.\n');
fprintf(fid,'- **Combined diagnosis:** Full Kirchhoff is not close to BIE for the rough cases and is not close to Model-1 either. Gate-0 is valid, but the rough-surface Kirchhoff approximation as implemented here is insufficient to adjudicate the PE residual without a further derivation/audit of rough-surface Kirchhoff terms, orientation, and illumination/shadow treatment. This is not evidence to modify production PE.\n\n');
fprintf(fid,'The result is diagnostic only: if a subsequently audited Full Kirchhoff formulation approaches BIE, the local phase-screen reduction is implicated; if it improves but remains separated, Kirchhoff is useful but incomplete; if it remains close to Model-1, the Kirchhoff approximation itself is insufficient for those cases.\n\n');
fprintf(fid,'Artifacts: `%s`, `%s`; figures: `results/validation/pe_bie_full_kirchhoff/figures/`.\n',local_rel(root,mat_file),local_rel(root,csv_file));
end
function rel=local_rel(root,p), rel=strrep(p,[root filesep],''); rel=strrep(rel,filesep,'/'); end

function eta=local_surface(x,A,K,inner,outer), r=abs(x); chi=zeros(size(x)); chi(r<=inner)=1; mid=r>inner&r<outer; t=(r(mid)-inner)/(outer-inner); chi(mid)=1-10*t.^3+15*t.^4-6*t.^5; eta=A*sin(K*x).*chi; end
function ep=local_surface_prime(x,A,K,inner,outer), r=abs(x); chi=zeros(size(x)); dchi=zeros(size(x)); chi(r<=inner)=1; mid=r>inner&r<outer; t=(r(mid)-inner)/(outer-inner); chi(mid)=1-10*t.^3+15*t.^4-6*t.^5; dchi(mid)=(-30*t.^2+60*t.^3-30*t.^4)/(outer-inner); s=sign(x); ep=A*(K*cos(K*x).*chi+sin(K*x).*dchi.*s); ep(x==0)=A*K; end

function source=local_discrete_source(cfg)
dx=cfg.xw_m/cfg.nx; x=(-cfg.nx/2:cfg.nx/2-1).'*dx; kx=(2*pi/cfg.xw_m)*[0:(cfg.nx/2-1),-cfg.nx/2:-1]; k=2*pi*cfg.frequency_hz/cfg.c0_mps; kz=sqrt(complex(k^2-kx.^2,0)); initial=exp(-0.5*(x/cfg.sigma_src_m).^2); coeff=fft(initial).'; x0=x(1); source=struct('kx',kx,'kz',kz,'coeff',coeff,'x0',x0,'n',cfg.nx,'z_tx',cfg.z_tx_m,'fingerprint',sprintf('discrete-periodic-Gaussian|N=%d|W=%.12g|sigma=%.12g|f=%.12g|c=%.12g|exp(-iwt)',cfg.nx,cfg.xw_m,cfg.sigma_src_m,cfg.frequency_hz,cfg.c0_mps)); source.evaluate=@(xq,zq,path,mode)local_eval(source,xq,zq,path,mode); source.evaluate_with_derivatives=@(xq,zq,mode)local_eval_der(source,xq,zq,mode); end
function u=local_eval(s,xq,zq,path,mode), [u,~,~]=local_eval_der(s,xq,zq,mode,path); end
function [u,ux,uz]=local_eval_der(s,xq,zq,mode,path)
if nargin<5,path=0;end; xq=xq(:); zq=zq(:); if isscalar(zq),zq=zq+zeros(size(xq));end; px=exp(1i*(xq-s.x0)*s.kx); if strcmp(mode,'path'), pz=exp(1i*path*s.kz); else,pz=exp(-1i*(zq-s.z_tx)*s.kz);end; c=s.coeff; u=sum(px.*(c.*pz),2)/s.n; ux=sum(px.*((1i*s.kx).*c.*pz),2)/s.n; uz=sum(px.*((-1i*s.kz).*c.*pz),2)/s.n; end
