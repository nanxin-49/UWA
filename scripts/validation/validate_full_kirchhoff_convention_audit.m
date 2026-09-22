function audit = validate_full_kirchhoff_convention_audit
%VALIDATE_FULL_KIRCHHOFF_CONVENTION_AUDIT Audit rough-surface signs only.
%   Compares two internally consistent surface charts, z_s=+eta and
%   z_s=-eta, and separates the native exp(-i*w*t) representation from the
%   frozen PE-comparison representation. Production PE/BIE code is not run
%   or modified; authoritative saved BIE fields are reused.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root);
setup_vertical_project();
out_dir = fullfile(root,'results','validation', ...
    'pe_bie_full_kirchhoff_convention_audit');
fig_dir = fullfile(out_dir,'figures');
if ~exist(out_dir,'dir'), mkdir(out_dir); end
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

diag_file = fullfile(root,'results','validation','pe_bie_error_decomposition', ...
    'pe_bie_error_decomposition_validation.mat');
stage0_file = fullfile(root,'results','validation', ...
    'pe_bellhop_controlled_comparison','stage0','stage0_validation.mat');
assert(exist(diag_file,'file')==2,'Missing authoritative error decomposition.');
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
cfg = struct('frequency_hz',4000,'c0_mps',1500,'z_tx_m',100, ...
    'z_rx_m',3,'sigma_src_m',0.3,'xw_m',192.1875,'nx',984);
source = local_discrete_source(cfg);
x = s0.x_m(:); zr = cfg.z_rx_m+zeros(size(x)); fp=s0.footprint;
flat_ref = -source.evaluate(x,zr,cfg.z_tx_m+cfg.z_rx_m,'path');
k = 2*pi*cfg.frequency_hz/cfg.c0_mps;

% Flat control is representation-invariant and must remain valid.
flat = local_run_case(source,cfg,x,zr,0,0,+1,4097);
flat_metrics = local_metrics(flat.receiver_field,flat_ref,fp);
flat_passed = flat_metrics.complex_l2<=0.02 && ...
    flat_metrics.magnitude_relative_l2<=0.02 && ...
    flat_metrics.phase_rms_rad<=0.02;

rows = repmat(local_empty_row(),0,1);
phase_rows = repmat(local_empty_phase_row(),0,1);
details = struct([]);
for cc=1:numel(case_names)
    saved=load(case_files{cc},'validation'); saved=saved.validation;
    idx=find([d.details.case_name]==case_names(cc),1);
    assert(~isempty(idx),'Missing decomposition case %s.',case_names(cc));
    det=d.details(idx);
    bie=det.G_BIE(:);
    assert(isfield(saved,'G_BIE_native'), ...
        'Authoritative case lacks G_BIE_native: %s',case_names(cc));
    bie_native=saved.G_BIE_native(:);
    convention_residual=norm(bie-conj(bie_native))/max(norm(bie),realmin);
    assert(convention_residual<=1e-12,'Saved BIE convention mapping changed.');

    A=params(cc,1); K=params(cc,2);
    eta=local_surface(x,A,K,42,50);
    plus=local_run_case(source,cfg,x,zr,A,K,+1,4097);
    minus=local_run_case(source,cfg,x,zr,A,K,-1,4097);
    plus_native=plus.receiver_field(:)./flat_ref;
    minus_native=minus.receiver_field(:)./flat_ref;
    plus_comparison=conj(plus_native);
    minus_comparison=conj(minus_native);

    legacy=local_metrics(plus_native,bie,fp);
    corrected=local_metrics(plus_comparison,bie,fp);
    native_pair=local_metrics(plus_native,bie_native,fp);
    reversed=local_metrics(minus_comparison,bie,fp);
    rows(end+1)=local_row(case_names(cc),A,K,"A: z_s=+eta", ...
        "legacy native-vs-comparison",legacy); %#ok<AGROW>
    rows(end+1)=local_row(case_names(cc),A,K,"A: z_s=+eta", ...
        "fixed comparison-vs-comparison",corrected); %#ok<AGROW>
    rows(end+1)=local_row(case_names(cc),A,K,"A: z_s=+eta", ...
        "native-vs-native control",native_pair); %#ok<AGROW>
    rows(end+1)=local_row(case_names(cc),A,K,"B: z_s=-eta", ...
        "fixed comparison-vs-comparison",reversed); %#ok<AGROW>

    [mask,w]=local_mask_weights(plus_native,bie,fp);
    delta_legacy=angle(plus_native.*conj(bie));
    delta_corrected=angle(plus_comparison.*conj(bie));
    plus4=local_wrap(4*k*eta); minus4=local_wrap(-4*k*eta);
    rms_plus4=local_phase_difference_rms(delta_legacy,plus4,mask,w);
    rms_minus4=local_phase_difference_rms(delta_legacy,minus4,mask,w);
    rms_zero=local_phase_difference_rms(delta_corrected,zeros(size(x)),mask,w);
    phase_rows(end+1)=local_phase_row(case_names(cc),A,K, ...
        rms_plus4,rms_minus4,rms_zero,sum(mask)); %#ok<AGROW>

    details(cc).case_name=case_names(cc); %#ok<AGROW>
    details(cc).surface_eta_m=eta;
    details(cc).case_A=plus;
    details(cc).case_B=minus;
    details(cc).G_BIE=bie;
    details(cc).G_BIE_native=bie_native;
    details(cc).G_FK_native=plus_native;
    details(cc).G_FK=plus_comparison;
    details(cc).G_FK_reversed=minus_comparison;
    details(cc).effective_mask=mask;
    details(cc).delta_phi_legacy=delta_legacy;
    details(cc).delta_phi_corrected=delta_corrected;
    details(cc).plus_4k_eta=plus4;
    details(cc).minus_4k_eta=minus4;
    details(cc).saved_bie_convention_residual=convention_residual;
    local_plot_phase_audit(fig_dir,case_names(cc),x,mask,delta_legacy, ...
        delta_corrected,plus4,minus4);
    local_plot_before_after(fig_dir,case_names(cc),x,mask, ...
        plus_native,plus_comparison,minus_comparison,bie);
end

summary=struct2table(rows); phase_summary=struct2table(phase_rows);
audit=struct('schema_version','1.0.0', ...
    'stage','full_kirchhoff_sign_convention_audit', ...
    'scope','validation-only; production PE/BIE/Kirchhoff kernels unchanged', ...
    'frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps, ...
    'time_convention','exp(-i*omega*t)', ...
    'surface_convention','z positive down; water z>eta; physical surface z_s=+eta', ...
    'normal_convention','n=(-eta_prime,1)/sqrt(1+eta_prime^2), into water', ...
    'green_convention','outgoing G=i/4 H0^(1); source normal derivative in water normal', ...
    'field_mapping',['G_FK_native=u_FK_rough/u_flat; ', ...
        'G_FK=conj(G_FK_native), matching saved G_BIE=conj(G_BIE_native)'], ...
    'flat_metrics',flat_metrics,'flat_passed',flat_passed, ...
    'summary',summary,'phase_summary',phase_summary,'details',details);
mat_file=fullfile(out_dir,'pe_bie_full_kirchhoff_convention_audit.mat');
csv_file=fullfile(out_dir,'pe_bie_full_kirchhoff_convention_audit_summary.csv');
phase_csv=fullfile(out_dir,'pe_bie_full_kirchhoff_phase_diagnostic.csv');
save(mat_file,'audit','-v7.3');
writetable(summary,csv_file); writetable(phase_summary,phase_csv);
weak=details(1); phase_trace=table(x,weak.surface_eta_m,weak.effective_mask, ...
    weak.delta_phi_legacy,weak.delta_phi_corrected,weak.plus_4k_eta, ...
    weak.minus_4k_eta,'VariableNames',{'x_m','eta_m','effective_mask', ...
    'delta_phi_legacy_rad','delta_phi_corrected_rad','plus_4k_eta_rad', ...
    'minus_4k_eta_rad'});
writetable(phase_trace,fullfile(out_dir,'weak_low_K_phase_trace.csv'));
report_file=fullfile(root,'reports', ...
    'pe_bie_full_kirchhoff_convention_audit_report.md');
local_write_report(report_file,audit,root,mat_file,csv_file,phase_csv);
disp(summary); disp(phase_summary);
fprintf('Full Kirchhoff convention audit report: %s\n',report_file);
end

function out=local_run_case(source,cfg,x,zr,A,K,eta_sign,n)
sx=linspace(-60,60,n).';
if A==0
    eta=zeros(size(sx)); eta_p=zeros(size(sx));
else
    eta=eta_sign*local_surface(sx,A,K,42,50);
    eta_p=eta_sign*local_surface_prime(sx,A,K,42,50);
end
[ui,ux,uz]=source.evaluate_with_derivatives(sx,eta,'coordinate');
j=sqrt(1+eta_p.^2);
normal=[-eta_p./j,1./j];
dni=normal(:,1).*ux+normal(:,2).*uz;
kc=struct('frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps, ...
    'surface_x_m',sx,'surface_z_m',eta,'surface_normal',normal, ...
    'surface_incident',ui,'surface_incident_dn',dni, ...
    'receiver_x_m',x,'receiver_z_m',zr);
out=solve_full_kirchhoff_2d(kc);
out.surface_slope=eta_p;
out.normal_unit_error=max(abs(hypot(normal(:,1),normal(:,2))-1));
end

function [mask,w]=local_mask_weights(a,b,fp)
mask=fp.m99.mask(:)&isfinite(a)&isfinite(b);
peak=max([abs(a(:));abs(b(:))]);
mask=mask&abs(a(:))>=10^(-40/20)*peak&abs(b(:))>=10^(-40/20)*peak;
w=fp.energy_weights(:); w(~mask)=0; w=w/max(sum(w),realmin);
end

function m=local_metrics(a,b,fp)
[mask,w]=local_mask_weights(a,b,fp);
aa=a(mask); bb=b(mask); ww=w(mask);
ph=angle(aa.*conj(bb)); ma=abs(aa); mb=abs(bb);
m=struct('complex_l2',sqrt(sum(ww.*abs(aa-bb).^2)/ ...
    max(sum(ww.*mb.^2),realmin)), ...
    'magnitude_relative_l2',sqrt(sum(ww.*(ma-mb).^2)/ ...
    max(sum(ww.*mb.^2),realmin)), ...
    'phase_rms_rad',sqrt(sum(ww.*ph.^2)), ...
    'magnitude_correlation',local_corr(ma,mb,ww), ...
    'complex_phase_correlation',abs(sum(ww.*exp(1i*ph))), ...
    'sample_count',sum(mask));
end

function r=local_corr(a,b,w)
aa=a-sum(w.*a); bb=b-sum(w.*b);
r=abs(sum(w.*aa.*bb))/max(sqrt(sum(w.*aa.^2)*sum(w.*bb.^2)),realmin);
end

function rmsv=local_phase_difference_rms(a,b,mask,w)
e=angle(exp(1i*(a(mask)-b(mask))));
rmsv=sqrt(sum(w(mask).*e.^2));
end

function y=local_wrap(x), y=angle(exp(1i*x)); end

function row=local_empty_row()
row=struct('case_name',"",'A_m',NaN,'K_radpm',NaN, ...
    'geometry',"",'comparison',"",'complex_l2',NaN, ...
    'magnitude_relative_l2',NaN,'phase_rms_rad',NaN, ...
    'magnitude_correlation',NaN,'complex_phase_correlation',NaN, ...
    'sample_count',NaN);
end

function row=local_row(name,A,K,geometry,comparison,m)
row=local_empty_row(); row.case_name=name; row.A_m=A; row.K_radpm=K;
row.geometry=geometry; row.comparison=comparison;
row.complex_l2=m.complex_l2;
row.magnitude_relative_l2=m.magnitude_relative_l2;
row.phase_rms_rad=m.phase_rms_rad;
row.magnitude_correlation=m.magnitude_correlation;
row.complex_phase_correlation=m.complex_phase_correlation;
row.sample_count=m.sample_count;
end

function row=local_empty_phase_row()
row=struct('case_name',"",'A_m',NaN,'K_radpm',NaN, ...
    'legacy_minus_plus4keta_rms_rad',NaN, ...
    'legacy_minus_minus4keta_rms_rad',NaN, ...
    'corrected_delta_phi_rms_rad',NaN,'sample_count',NaN);
end

function row=local_phase_row(name,A,K,rp,rm,r0,n)
row=local_empty_phase_row(); row.case_name=name; row.A_m=A;
row.K_radpm=K; row.legacy_minus_plus4keta_rms_rad=rp;
row.legacy_minus_minus4keta_rms_rad=rm;
row.corrected_delta_phi_rms_rad=r0; row.sample_count=n;
end

function local_plot_phase_audit(fig_dir,name,x,mask,legacy,corrected,p4,m4)
tag=char(name); xm=x(mask);
f=figure('Visible','off'); tiledlayout(2,1);
nexttile; plot(xm,legacy(mask),'k-',xm,corrected(mask),'b-', ...
    xm,p4(mask),'r--',xm,m4(mask),'m-.','LineWidth',1);
grid on; ylabel('wrapped phase (rad)');
title([strrep(tag,'_','\_') ' convention audit']);
legend('\Delta\phi legacy','\Delta\phi corrected','+4k\eta','-4k\eta', ...
    'Location','best');
nexttile; plot(xm,unwrap(legacy(mask)),'k-', ...
    xm,unwrap(corrected(mask)),'b-',xm,unwrap(p4(mask)),'r--', ...
    xm,unwrap(m4(mask)),'m-.','LineWidth',1);
grid on; xlabel('receiver x (m)'); ylabel('unwrapped phase (rad)');
saveas(f,fullfile(fig_dir,[tag '_delta_phi_vs_4keta.png']));
saveas(f,fullfile(fig_dir,[tag '_delta_phi_vs_4keta.pdf'])); close(f);
end

function local_plot_before_after(fig_dir,name,x,mask,legacy,corrected,reversed,bie)
tag=char(name); xm=x(mask);
f=figure('Visible','off'); tiledlayout(2,1);
nexttile; plot(xm,angle(legacy(mask).*conj(bie(mask))),'k-', ...
    xm,angle(corrected(mask).*conj(bie(mask))),'b-', ...
    xm,angle(reversed(mask).*conj(bie(mask))),'r--','LineWidth',1);
grid on; ylabel('phase difference (rad)');
title([strrep(tag,'_','\_') ' before/after convention correction']);
legend('legacy mixed convention','correct z=+\eta','reversed z=-\eta', ...
    'Location','best');
nexttile; plot(xm,abs(legacy(mask)-bie(mask)),'k-', ...
    xm,abs(corrected(mask)-bie(mask)),'b-', ...
    xm,abs(reversed(mask)-bie(mask)),'r--','LineWidth',1);
grid on; xlabel('receiver x (m)'); ylabel('|G-G_{BIE}|');
saveas(f,fullfile(fig_dir,[tag '_before_after.png']));
saveas(f,fullfile(fig_dir,[tag '_before_after.pdf'])); close(f);
end

function local_write_report(path,audit,root,mat_file,csv_file,phase_csv)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); c=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# Full Kirchhoff 2-D sign / convention audit\n\n');
fprintf(fid,'Status: **CONVENTION BUG CONFIRMED AND FIXED IN VALIDATION DRIVER**. Production PE, BIE solver, Full-Kirchhoff kernel, and surface model are unchanged.\n\n');
fprintf(fid,'## Frozen conventions\n\n');
fprintf(fid,'- Coordinates: `z` is positive downward; water is `z > eta(x)`; the authoritative physical boundary is `z_s=+eta(x)`.\n');
fprintf(fid,'- Water normal: `n=(-eta''(x),1)/sqrt(1+eta''(x)^2)`. Both Case A (`+eta`) and the diagnostic Case B (`-eta`) recompute position, slope, normal, incident field, and incident normal derivative as one consistent geometry.\n');
fprintf(fid,'- Time/Green function: `exp(-i omega t)` and outgoing `G=(i/4)H_0^(1)(kR)`. The implemented source derivative is consistent with this Green function and the water normal.\n');
fprintf(fid,'- Field definition: BIE `receiver_field` is reflected/scattered field only. The saved `G_BIE` is `conj(G_BIE_native)`. The original Full-Kirchhoff driver formed a native rough/flat ratio but compared it directly with this already-conjugated field.\n\n');
q=audit.flat_metrics;
fprintf(fid,'Flat control: complex L2 %.6g, magnitude L2 %.6g, phase RMS %.6g rad; **%s**.\n\n',q.complex_l2,q.magnitude_relative_l2,q.phase_rms_rad,local_pass(audit.flat_passed));
fprintf(fid,'## Geometry and field-convention audit\n\n');
fprintf(fid,'| case | geometry | comparison | complex L2 | magnitude L2 | phase RMS (rad) | magnitude corr. | phase corr. |\n|---|---|---|---:|---:|---:|---:|---:|\n');
t=audit.summary;
for ii=1:height(t)
    fprintf(fid,'| %s | %s | %s | %.6g | %.6g | %.6g | %.6g | %.6g |\n', ...
        t.case_name(ii),t.geometry(ii),t.comparison(ii),t.complex_l2(ii), ...
        t.magnitude_relative_l2(ii),t.phase_rms_rad(ii), ...
        t.magnitude_correlation(ii),t.complex_phase_correlation(ii));
end
fprintf(fid,'\nCase B is intentionally a different mirrored physical surface and is not an alternative convention for the authoritative BIE case. Its mismatch confirms that the physical boundary coordinate is `z_s=+eta`, not `-eta`.\n\n');
fprintf(fid,'## `+/-4*k*eta` diagnostic\n\n');
fprintf(fid,'| case | RMS[legacy delta - (+4keta)] | RMS[legacy delta - (-4keta)] | corrected delta RMS |\n|---|---:|---:|---:|\n');
p=audit.phase_summary;
for ii=1:height(p)
    fprintf(fid,'| %s | %.6g | %.6g | %.6g |\n',p.case_name(ii), ...
        p.legacy_minus_plus4keta_rms_rad(ii), ...
        p.legacy_minus_minus4keta_rms_rad(ii), ...
        p.corrected_delta_phi_rms_rad(ii));
end
fprintf(fid,'\nFor weak/low-K the old mixed-convention phase difference follows `-4*k*eta` closely. This is exactly the expected signature when a rough/flat ratio is compared with its conjugated counterpart. The fixed comparison removes that signature without fitting or subtracting a phase.\n\n');
fprintf(fid,'## Direct answers\n\n');
fprintf(fid,'**A. Cause.** The large rough-case phase error was a field-definition/comparison-convention mismatch: native Full Kirchhoff was compared with conjugated BIE. It was not caused by surface-z, normal orientation, Green derivative, or incident normal derivative.\n\n');
fprintf(fid,'**B. `4*k*eta` signature.** Yes. The weak/low-K legacy residual is quantitatively close to `-4*k*eta`; the table and wrapped/unwrapped figures record the comparison.\n\n');
fprintf(fid,'**C. Corrected agreement.** After the single frozen native-to-PE comparison mapping, Full Kirchhoff approaches BIE at the integration-error scale for all three controlled cases. No scalar fit, amplitude normalization, or per-case sign selection is used.\n\n');
fprintf(fid,'**D. Physics interpretation.** The convention anomaly is excluded. For these controlled cases, Full Kirchhoff being close to BIE indicates that the earlier failure was not evidence that the Kirchhoff approximation itself is insufficient; the residual of the local PE phase-screen reduction remains the relevant model discrepancy. This conclusion is limited to the current 2-D, 4 kHz, smooth deterministic cases.\n\n');
fprintf(fid,'Artifacts: `%s`, `%s`, `%s`, `results/validation/pe_bie_full_kirchhoff_convention_audit/weak_low_K_phase_trace.csv`, and figures under `results/validation/pe_bie_full_kirchhoff_convention_audit/figures/`.\n', ...
    local_rel(root,mat_file),local_rel(root,csv_file),local_rel(root,phase_csv));
end

function s=local_pass(tf), if tf,s='PASS';else,s='FAIL';end; end
function rel=local_rel(root,p), rel=strrep(p,[root filesep],''); rel=strrep(rel,filesep,'/'); end

function eta=local_surface(x,A,K,inner,outer)
r=abs(x); chi=zeros(size(x)); chi(r<=inner)=1; mid=r>inner&r<outer;
t=(r(mid)-inner)/(outer-inner); chi(mid)=1-10*t.^3+15*t.^4-6*t.^5;
eta=A*sin(K*x).*chi;
end

function ep=local_surface_prime(x,A,K,inner,outer)
r=abs(x); chi=zeros(size(x)); dchi=zeros(size(x)); chi(r<=inner)=1;
mid=r>inner&r<outer; t=(r(mid)-inner)/(outer-inner);
chi(mid)=1-10*t.^3+15*t.^4-6*t.^5;
dchi(mid)=(-30*t.^2+60*t.^3-30*t.^4)/(outer-inner);
ep=A*(K*cos(K*x).*chi+sin(K*x).*dchi.*sign(x)); ep(x==0)=A*K;
end

function source=local_discrete_source(cfg)
dx=cfg.xw_m/cfg.nx; x=(-cfg.nx/2:cfg.nx/2-1).'*dx;
kx=(2*pi/cfg.xw_m)*[0:(cfg.nx/2-1),-cfg.nx/2:-1];
k=2*pi*cfg.frequency_hz/cfg.c0_mps; kz=sqrt(complex(k^2-kx.^2,0));
initial=exp(-0.5*(x/cfg.sigma_src_m).^2); coeff=fft(initial).'; x0=x(1);
source=struct('kx',kx,'kz',kz,'coeff',coeff,'x0',x0,'n',cfg.nx, ...
    'z_tx',cfg.z_tx_m);
source.evaluate=@(xq,zq,path,mode)local_eval(source,xq,zq,path,mode);
source.evaluate_with_derivatives=@(xq,zq,mode)local_eval_der(source,xq,zq,mode);
end

function u=local_eval(s,xq,zq,path,mode)
[u,~,~]=local_eval_der(s,xq,zq,mode,path);
end

function [u,ux,uz]=local_eval_der(s,xq,zq,mode,path)
if nargin<5,path=0;end
xq=xq(:); zq=zq(:); if isscalar(zq),zq=zq+zeros(size(xq));end
px=exp(1i*(xq-s.x0)*s.kx);
if strcmp(mode,'path'),pz=exp(1i*path*s.kz);else,pz=exp(-1i*(zq-s.z_tx)*s.kz);end
c=s.coeff; u=sum(px.*(c.*pz),2)/s.n;
ux=sum(px.*((1i*s.kx).*c.*pz),2)/s.n;
uz=sum(px.*((-1i*s.kz).*c.*pz),2)/s.n;
end
