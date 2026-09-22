function result = validate_full_kirchhoff_reduced_operator_stage1
%VALIDATE_FULL_KIRCHHOFF_REDUCED_OPERATOR_STAGE1
% Validation-only ablation and localization audit of the accepted 2-D
% Full-Kirchhoff model. Production PE, BIE, and surface physics are unchanged.

root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root);setup_vertical_project();
out_dir=fullfile(root,'results','validation','pe_bie_full_kirchhoff_reduced_stage1');
fig_dir=fullfile(out_dir,'figures');
if ~exist(out_dir,'dir'),mkdir(out_dir);end
if ~exist(fig_dir,'dir'),mkdir(fig_dir);end

sweep_file=fullfile(root,'results','validation','pe_bie_full_kirchhoff_parameter_sweep', ...
    'pe_bie_full_kirchhoff_parameter_sweep.mat');
stage0_file=fullfile(root,'results','validation','pe_bellhop_controlled_comparison', ...
    'stage0','stage0_validation.mat');
assert(exist(sweep_file,'file')==2,'Run the accepted Full-Kirchhoff sweep first.');
s=load(sweep_file,'result');s=s.result;
q=load(stage0_file,'validation');s0=q.validation;
cfg=s.config;source=local_source(cfg);x=s0.x_m(:);zr=cfg.z_rx_m+zeros(size(x));
fp=s0.footprint;flat_ref=-source.evaluate(x,zr,cfg.z_tx_m+cfg.z_rx_m,'path');

names=["weak_low_K","strong_height_low_K","weak_high_K"];
indices=[2 6 11];
radii=[0.5 1 2 4 8 16 32 64 Inf];
summary_rows=repmat(local_row(),0,1);
radius_rows=repmat(local_radius_row(),0,1);
details=repmat(local_detail(),numel(indices),1);

% PS--FK sweep metrics reuse the already accepted fields; no solver rerun.
sweep_rows=repmat(local_sweep_row(),numel(s.cases),1);
for ii=1:numel(s.cases)
    c=s.cases(ii);m=local_metrics(c.receiver_PE,c.receiver_FK,fp,cfg.phase_floor_db);
    sweep_rows(ii)=struct('family',c.family,'A_m',c.A_m,'K_radpm',c.K_radpm, ...
        'slope_AK',c.geometry.slope_AK,'curvature_AK2',c.geometry.curvature_AK2, ...
        'complex_l2',m.complex_l2,'magnitude_relative_l2',m.magnitude_relative_l2, ...
        'phase_rms_rad',m.phase_rms_rad,'magnitude_correlation',m.magnitude_correlation, ...
        'complex_phase_correlation',m.complex_phase_correlation);
end

for cc=1:numel(indices)
    base=s.cases(indices(cc));A=base.A_m;K=base.K_radpm;
    sx=base.surface_x_m(:);sz=base.surface_z_m(:);
    sp=local_surface_prime(sx,A,K,cfg.surface_inner_m,cfg.surface_support_m);
    jac=sqrt(1+sp.^2);normal=[-sp./jac,1./jac];flat_normal=[zeros(size(sp)),ones(size(sp))];
    [ui,ux,uz]=source.evaluate_with_derivatives(sx,sz,'coordinate');
    dni=normal(:,1).*ux+normal(:,2).*uz;dni_flat=uz;
    common=struct('frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps, ...
        'surface_x_m',sx,'surface_z_m',sz,'surface_incident',ui, ...
        'receiver_x_m',x,'receiver_z_m',zr);

    full_cfg=common;full_cfg.surface_normal=normal;full_cfg.surface_incident_dn=dni;
    full=solve_full_kirchhoff_2d(full_cfg);g_full=conj(full.receiver_field(:)./flat_ref);
    reconstruction=local_metrics(g_full,base.receiver_FK,fp,cfg.phase_floor_db);
    assert(reconstruction.complex_l2<1e-12,'Full-field reconstruction changed for %s.',names(cc));

    normal_cfg=common;normal_cfg.surface_normal=flat_normal;normal_cfg.surface_incident_dn=dni_flat;
    no_normal=solve_full_kirchhoff_2d(normal_cfg);g_no_normal=conj(no_normal.receiver_field(:)./flat_ref);
    jac_cfg=full_cfg;jac_cfg.quadrature_measure="parameter_x";
    no_jac=solve_full_kirchhoff_2d(jac_cfg);g_no_jac=conj(no_jac.receiver_field(:)./flat_ref);
    % Surface-height/path ablation with consistent incident evaluation.
    [ui0,~,uz0]=source.evaluate_with_derivatives(sx,zeros(size(sx)),'coordinate');
    flat_surface_cfg=common;flat_surface_cfg.surface_z_m=zeros(size(sx));
    flat_surface_cfg.surface_normal=flat_normal;flat_surface_cfg.surface_incident=ui0;
    flat_surface_cfg.surface_incident_dn=uz0;
    no_height=solve_full_kirchhoff_2d(flat_surface_cfg);
    g_no_height=conj(no_height.receiver_field(:)./flat_ref);

    m_ps=local_metrics(base.receiver_PE,g_full,fp,cfg.phase_floor_db);
    m_normal=local_metrics(g_no_normal,g_full,fp,cfg.phase_floor_db);
    m_jac=local_metrics(g_no_jac,g_full,fp,cfg.phase_floor_db);
    summary_rows(end+1)=local_make_row(names(cc),A,K,"PE phase screen",m_ps); %#ok<AGROW>
    summary_rows(end+1)=local_make_row(names(cc),A,K,"flat normal",m_normal); %#ok<AGROW>
    summary_rows(end+1)=local_make_row(names(cc),A,K,"ds=dx",m_jac); %#ok<AGROW>
    m_height=local_metrics(g_no_height,g_full,fp,cfg.phase_floor_db);
    summary_rows(end+1)=local_make_row(names(cc),A,K,"flat surface z=0",m_height); %#ok<AGROW>
    summary_rows(end+1)=local_make_row(names(cc),A,K,"Full Kirchhoff",local_self_metrics(g_full,fp)); %#ok<AGROW>

    localized=complex(zeros(numel(x),numel(radii)));
    correction_error=zeros(size(radii));
    delta_full=g_full-base.receiver_PE;
    for rr=1:numel(radii)
        lc=full_cfg;lc.local_radius_m=radii(rr);
        lr=solve_full_kirchhoff_2d(lc);localized(:,rr)=conj(lr.receiver_field(:)./flat_ref);
        ml=local_metrics(localized(:,rr),g_full,fp,cfg.phase_floor_db);
        correction_error(rr)=local_relative_correction(localized(:,rr)-base.receiver_PE,delta_full,fp);
        radius_rows(end+1)=local_make_radius_row(names(cc),A,K,radii(rr),ml,correction_error(rr)); %#ok<AGROW>
    end
    spatial=local_spatial_diagnostics(x,base.receiver_PE,g_full,sz,sx,sp,fp,cfg.phase_floor_db);
    details(cc)=struct('case_name',names(cc),'A_m',A,'K_radpm',K,'x_m',x, ...
        'surface_x_m',sx,'surface_z_m',sz,'receiver_PS',base.receiver_PE, ...
        'receiver_FK',g_full,'receiver_BIE',base.receiver_BIE, ...
        'receiver_flat_normal',g_no_normal,'receiver_no_jacobian',g_no_jac, ...
        'receiver_flat_surface',g_no_height, ...
        'localized_radius_m',radii,'receiver_localized',localized, ...
        'metrics_PS_FK',m_ps,'metrics_flat_normal_FK',m_normal, ...
        'metrics_no_jacobian_FK',m_jac,'correction_relative_error',correction_error, ...
        'spatial_diagnostics',spatial,'reconstruction',reconstruction);
    local_plot_case(fig_dir,details(cc),fp);
    local_plot_ablation(fig_dir,details(cc));
    local_plot_radius(fig_dir,details(cc),radius_rows);
end

sweep_table=struct2table(sweep_rows);summary_table=struct2table(summary_rows);
radius_table=struct2table(radius_rows);
local_plot_sweep(fig_dir,sweep_table);
result=struct('schema_version','1.0.0','stage','full_kirchhoff_reduced_operator_stage1', ...
    'scope','validation-only diagnostic; production PE/BIE/surface unchanged', ...
    'config',cfg,'representative_case_names',names,'local_radius_m',radii, ...
    'mask_definition','Fixed Stage-0 M99 and Full-Kirchhoff-reference -40 dB mask with saved incident-energy weights', ...
    'summary',summary_table,'radius_scan',radius_table,'sweep_PS_FK',sweep_table, ...
    'details',details);
mat_file=fullfile(out_dir,'pe_bie_full_kirchhoff_reduced_stage1.mat');
summary_csv=fullfile(out_dir,'pe_bie_full_kirchhoff_reduced_stage1_summary.csv');
radius_csv=fullfile(out_dir,'pe_bie_full_kirchhoff_reduced_stage1_radius_scan.csv');
sweep_csv=fullfile(out_dir,'pe_bie_full_kirchhoff_reduced_stage1_sweep.csv');
save(mat_file,'result','-v7.3');writetable(summary_table,summary_csv);writetable(radius_table,radius_csv);writetable(sweep_table,sweep_csv);
report=fullfile(root,'reports','pe_bie_full_kirchhoff_reduced_operator_stage1_report.md');
local_write_report(report,result,root,mat_file,summary_csv,radius_csv,sweep_csv);
disp(summary_table);fprintf('Reduced Full-Kirchhoff Stage-1 report: %s\n',report);
end

function s=local_source(c)
dx=c.xw_m/c.nx;x=(-c.nx/2:c.nx/2-1).'*dx;kx=(2*pi/c.xw_m)*[0:(c.nx/2-1),-c.nx/2:-1];k=2*pi*c.frequency_hz/c.c0_mps;kz=sqrt(complex(k^2-kx.^2,0));coeff=fft(exp(-0.5*(x/c.sigma_src_m).^2)).';x0=x(1);
s=struct('kx',kx,'kz',kz,'coeff',coeff,'x0',x0,'n',c.nx,'z_tx',c.z_tx_m);
s.evaluate=@(xq,zq,path,mode)local_eval(s,xq,zq,path,mode);s.evaluate_with_derivatives=@(xq,zq,mode)local_der(s,xq,zq,mode);
end
function u=local_eval(s,xq,zq,path,mode),xq=xq(:);zq=zq(:);if isscalar(zq),zq=zq+zeros(size(xq));end;px=exp(1i*(xq-s.x0)*s.kx);if strcmp(mode,'path'),pz=exp(1i*path*s.kz);else,pz=exp(-1i*(zq-path)*s.kz);end;u=sum(px.*(s.coeff.*pz),2)/s.n;end
function [u,ux,uz]=local_der(s,xq,zq,mode),assert(strcmp(mode,'coordinate'));xq=xq(:);zq=zq(:);px=exp(1i*(xq-s.x0)*s.kx);pz=exp(-1i*(zq-s.z_tx)*s.kz);u=sum(px.*(s.coeff.*pz),2)/s.n;ux=sum(px.*((1i*s.kx).*s.coeff.*pz),2)/s.n;uz=sum(px.*((-1i*s.kz).*s.coeff.*pz),2)/s.n;end
function ep=local_surface_prime(x,A,K,inner,outer),r=abs(x);chi=zeros(size(x));dchi=zeros(size(x));chi(r<=inner)=1;mid=r>inner&r<outer;t=(r(mid)-inner)/(outer-inner);chi(mid)=1-10*t.^3+15*t.^4-6*t.^5;dchi(mid)=(-30*t.^2+60*t.^3-30*t.^4)/(outer-inner);ep=A*(K*cos(K*x).*chi+sin(K*x).*dchi.*sign(x));ep(x==0)=A*K;end

function [mask,w]=local_mask(ref,fp,floor_db),ref=ref(:);mask=fp.m99.mask(:)&isfinite(ref);peak=max(abs(ref));mask=mask&abs(ref)>=10^(floor_db/20)*peak;w=fp.energy_weights(:);w=w(mask);w=w/max(sum(w),realmin);end
function m=local_metrics(a,b,fp,floor_db),a=a(:);b=b(:);[mask,~]=local_mask(b,fp,floor_db);mask=mask&isfinite(a);wf=fp.energy_weights(:);w=wf(mask);w=w/max(sum(w),realmin);a=a(mask);b=b(mask);ph=angle(a.*conj(b));ma=abs(a);mb=abs(b);m=struct('complex_l2',sqrt(sum(w.*abs(a-b).^2)/max(sum(w.*abs(b).^2),realmin)),'magnitude_relative_l2',sqrt(sum(w.*(ma-mb).^2)/max(sum(w.*mb.^2),realmin)),'phase_rms_rad',sqrt(sum(w.*ph.^2)),'magnitude_correlation',local_corr(ma,mb,w),'complex_phase_correlation',abs(sum(w.*exp(1i*ph))),'sample_count',numel(a));end
function m=local_self_metrics(a,fp),[mask,~]=local_mask(a,fp,-40);m=struct('complex_l2',0,'magnitude_relative_l2',0,'phase_rms_rad',0,'magnitude_correlation',1,'complex_phase_correlation',1,'sample_count',sum(mask));end
function r=local_corr(a,b,w),aa=a-sum(w.*a);bb=b-sum(w.*b);r=abs(sum(w.*aa.*bb))/max(sqrt(sum(w.*aa.^2)*sum(w.*bb.^2)),realmin);end
function e=local_relative_correction(test,ref,fp),mask=fp.m99.mask(:)&isfinite(test)&isfinite(ref);w=fp.energy_weights(:);w=w(mask);w=w/max(sum(w),realmin);e=sqrt(sum(w.*abs(test(mask)-ref(mask)).^2)/max(sum(w.*abs(ref(mask)).^2),realmin));end
function d=local_spatial_diagnostics(x,ps,fk,eta,sx,slope,fp,floor_db),[mask,w]=local_mask(fk,fp,floor_db);err=abs(ps(mask)-fk(mask));e2=w.*err.^2;r=abs(x(mask));[rs,ord]=sort(r);cs=cumsum(e2(ord))/max(sum(e2),realmin);r95=rs(find(cs>=0.95,1));eh=interp1(sx,abs(eta),x(mask),'linear',0);es=interp1(sx,abs(slope),x(mask),'linear',0);d=struct('error_r95_m',r95,'error_height_correlation',local_unweighted_corr(err,eh),'error_slope_correlation',local_unweighted_corr(err,es),'max_abs_phase_difference_rad',max(abs(angle(ps(mask).*conj(fk(mask))))),'max_abs_magnitude_difference_db',max(abs(20*log10(max(abs(ps(mask)),realmin)./max(abs(fk(mask)),realmin)))));end
function r=local_unweighted_corr(a,b),a=a(:);b=b(:);if std(a)==0||std(b)==0,r=NaN;else,q=corrcoef(a,b);r=q(1,2);end;end

function r=local_row(),r=struct('case_name',"",'A_m',NaN,'K_radpm',NaN,'model',"",'complex_l2',NaN,'magnitude_relative_l2',NaN,'phase_rms_rad',NaN,'magnitude_correlation',NaN,'complex_phase_correlation',NaN,'sample_count',NaN);end
function r=local_make_row(n,A,K,model,m),r=local_row();r.case_name=n;r.A_m=A;r.K_radpm=K;r.model=model;r.complex_l2=m.complex_l2;r.magnitude_relative_l2=m.magnitude_relative_l2;r.phase_rms_rad=m.phase_rms_rad;r.magnitude_correlation=m.magnitude_correlation;r.complex_phase_correlation=m.complex_phase_correlation;r.sample_count=m.sample_count;end
function r=local_radius_row(),r=struct('case_name',"",'A_m',NaN,'K_radpm',NaN,'radius_m',NaN,'complex_l2',NaN,'magnitude_relative_l2',NaN,'phase_rms_rad',NaN,'magnitude_correlation',NaN,'complex_phase_correlation',NaN,'correction_relative_error',NaN);end
function r=local_make_radius_row(n,A,K,L,m,c),r=local_radius_row();r.case_name=n;r.A_m=A;r.K_radpm=K;r.radius_m=L;r.complex_l2=m.complex_l2;r.magnitude_relative_l2=m.magnitude_relative_l2;r.phase_rms_rad=m.phase_rms_rad;r.magnitude_correlation=m.magnitude_correlation;r.complex_phase_correlation=m.complex_phase_correlation;r.correction_relative_error=c;end
function r=local_sweep_row(),r=struct('family',"",'A_m',NaN,'K_radpm',NaN,'slope_AK',NaN,'curvature_AK2',NaN,'complex_l2',NaN,'magnitude_relative_l2',NaN,'phase_rms_rad',NaN,'magnitude_correlation',NaN,'complex_phase_correlation',NaN);end
function r=local_detail(),r=struct('case_name',"",'A_m',NaN,'K_radpm',NaN,'x_m',[],'surface_x_m',[],'surface_z_m',[],'receiver_PS',[],'receiver_FK',[],'receiver_BIE',[],'receiver_flat_normal',[],'receiver_no_jacobian',[],'receiver_flat_surface',[],'localized_radius_m',[],'receiver_localized',[],'metrics_PS_FK',[],'metrics_flat_normal_FK',[],'metrics_no_jacobian_FK',[],'correction_relative_error',[],'spatial_diagnostics',[],'reconstruction',[]);end

function local_plot_case(fd,d,fp),[mask,~]=local_mask(d.receiver_FK,fp,-40);x=d.x_m;tag=char(d.case_name);f=figure('Visible','off');tiledlayout(2,1);nexttile;plot(x,abs(d.receiver_PS),'--',x,abs(d.receiver_FK),'-',x,abs(d.receiver_BIE),':','LineWidth',1);grid on;legend('PE phase screen','Full Kirchhoff','BIE');ylabel('magnitude');title(sprintf('%s: A=%.3g m, K=%.3g rad/m',tag,d.A_m,d.K_radpm));nexttile;plot(x,angle(d.receiver_PS),'--',x,angle(d.receiver_FK),'-',x,angle(d.receiver_BIE),':','LineWidth',1);grid on;ylabel('wrapped phase');xlabel('receiver x (m)');saveas(f,fullfile(fd,[tag '_fields.png']));saveas(f,fullfile(fd,[tag '_fields.pdf']));close(f);f=figure('Visible','off');tiledlayout(4,1);nexttile;plot(x(mask),abs(d.receiver_PS(mask)-d.receiver_FK(mask)));grid on;ylabel('|PS-FK|');nexttile;plot(x(mask),angle(d.receiver_PS(mask).*conj(d.receiver_FK(mask))));grid on;ylabel('wrapped phase');nexttile;plot(x(mask),unwrap(angle(d.receiver_PS(mask).*conj(d.receiver_FK(mask)))));grid on;ylabel('unwrapped phase');nexttile;plot(x(mask),20*log10(max(abs(d.receiver_PS(mask)),realmin)./max(abs(d.receiver_FK(mask)),realmin)));grid on;ylabel('magnitude diff (dB)');xlabel('receiver x (m)');saveas(f,fullfile(fd,[tag '_PS_FK_error.png']));saveas(f,fullfile(fd,[tag '_PS_FK_error.pdf']));close(f);end
function local_plot_ablation(fd,d),x=d.x_m;tag=char(d.case_name);f=figure('Visible','off');tiledlayout(2,1);nexttile;plot(x,abs(d.receiver_FK),'-',x,abs(d.receiver_flat_normal),'--',x,abs(d.receiver_no_jacobian),':',x,abs(d.receiver_flat_surface),'-.');grid on;legend('Full','flat normal','ds=dx','flat surface z=0');ylabel('magnitude');title([tag ' physical-term ablations']);nexttile;plot(x,angle(d.receiver_FK),'-',x,angle(d.receiver_flat_normal),'--',x,angle(d.receiver_no_jacobian),':',x,angle(d.receiver_flat_surface),'-.');grid on;ylabel('wrapped phase');xlabel('receiver x (m)');saveas(f,fullfile(fd,[tag '_ablations.png']));saveas(f,fullfile(fd,[tag '_ablations.pdf']));close(f);end
function local_plot_radius(fd,d,rows),t=rows([rows.case_name]==d.case_name);L=[t.radius_m].';Lf=L;Lf(isinf(Lf))=128;f=figure('Visible','off');semilogy(Lf,[t.complex_l2].','o-',Lf,[t.magnitude_relative_l2].','s-',Lf,[t.phase_rms_rad].','^-','LineWidth',1.1);grid on;legend('complex L2','magnitude L2','phase RMS');xlabel('local radius L (m); 128 denotes Inf');ylabel('error relative to Full Kirchhoff');title([char(d.case_name) ' localization scan']);saveas(f,fullfile(fd,[char(d.case_name) '_radius_scan.png']));saveas(f,fullfile(fd,[char(d.case_name) '_radius_scan.pdf']));close(f);end
function local_plot_sweep(fd,t),h=t(t.family=="height",:);k=t(t.family=="wavenumber",:);f=figure('Visible','off');tiledlayout(1,2);nexttile;semilogy(h.A_m,h.complex_l2,'o-',h.A_m,h.magnitude_relative_l2,'s-',h.A_m,h.phase_rms_rad,'^-');grid on;xlabel('A (m)');ylabel('PS-FK error');legend('complex','magnitude','phase');title('K=0.10 rad/m');nexttile;semilogy(k.K_radpm,k.complex_l2,'o-',k.K_radpm,k.magnitude_relative_l2,'s-',k.K_radpm,k.phase_rms_rad,'^-');grid on;xlabel('K (rad/m)');legend('complex','magnitude','phase');title('A=0.02 m');saveas(f,fullfile(fd,'PS_FK_parameter_trends.png'));saveas(f,fullfile(fd,'PS_FK_parameter_trends.pdf'));close(f);end

function local_write_report(path,r,root,mat,summary_csv,radius_csv,sweep_csv)
fid=fopen(path,'w','n','UTF-8');assert(fid>=0);cl=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# Full Kirchhoff to reduced-operator Stage-1 analysis\n\nStatus: **validation-only diagnostic; production PE/BIE/surface models unchanged**.\n\n');
fprintf(fid,'Configuration and convention are inherited from the accepted 4 kHz deterministic Full-Kirchhoff sweep. Metrics use a fixed Stage-0 M99 footprint, a Full-Kirchhoff-reference -40 dB mask, and saved incident-energy weights. No scalar fitting or case-dependent phase adjustment is used.\n\n');
fprintf(fid,'## PE phase screen versus Full Kirchhoff\n\n| case | complex L2 | magnitude L2 | phase RMS rad | magnitude corr. | phase corr. | error r95 m | corr(error,|eta|) | corr(error,|slope|) |\n|---|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for i=1:numel(r.details),d=r.details(i);m=d.metrics_PS_FK;z=d.spatial_diagnostics;fprintf(fid,'| %s | %.6g | %.6g | %.6g | %.6g | %.6g | %.6g | %.6g | %.6g |\n',d.case_name,m.complex_l2,m.magnitude_relative_l2,m.phase_rms_rad,m.magnitude_correlation,m.complex_phase_correlation,z.error_r95_m,z.error_height_correlation,z.error_slope_correlation);end
fprintf(fid,'\n## Physical-term ablations relative to Full Kirchhoff\n\n| case | model | complex L2 | magnitude L2 | phase RMS rad | magnitude corr. | phase corr. |\n|---|---|---:|---:|---:|---:|---:|\n');t=r.summary;for i=1:height(t),if t.model(i)=="Full Kirchhoff",continue;end;fprintf(fid,'| %s | %s | %.6g | %.6g | %.6g | %.6g | %.6g |\n',t.case_name(i),t.model{i},t.complex_l2(i),t.magnitude_relative_l2(i),t.phase_rms_rad(i),t.magnitude_correlation(i),t.complex_phase_correlation(i));end
fprintf(fid,'\n`flat normal` retains the actual surface position and arc-length measure but consistently uses `n=(0,1)` in both Green and incident normal derivatives. `ds=dx` retains the curved normal and actual surface position but removes only the arc-length Jacobian. `flat surface z=0` removes the physical surface height/path perturbation and recomputes the incident field consistently at z=0.\n\n');
fprintf(fid,'## Receiver-centered finite-aperture integral\n\n| case | L m | complex L2 | magnitude L2 | phase RMS rad | correction relative error |\n|---|---:|---:|---:|---:|---:|\n');t=r.radius_scan;for i=1:height(t),fprintf(fid,'| %s | %.6g | %.6g | %.6g | %.6g | %.6g |\n',t.case_name(i),t.radius_m(i),t.complex_l2(i),t.magnitude_relative_l2(i),t.phase_rms_rad(i),t.correction_relative_error(i));end
fprintf(fid,'\nThe finite-aperture model keeps all accepted Full-Kirchhoff terms but sets contributions with `|x_r-x_s|>L` to zero. Its correction diagnostic compares `(u_L-u_PS)` with the complete `(u_FK-u_PS)`; it is not a fitted convolution kernel.\n\n');
fprintf(fid,'## Interpretation and ranking\n\n');
for i=1:numel(r.details),d=r.details(i);ps=d.metrics_PS_FK;complex=[d.metrics_flat_normal_FK.complex_l2,d.metrics_no_jacobian_FK.complex_l2];t=r.radius_scan(r.radius_scan.case_name==d.case_name,:);target=0.5*ps.complex_l2;j=find(t.complex_l2<=target,1);if isempty(j),Ltxt='not reached';else,Ltxt=sprintf('%.6g m',t.radius_m(j));end;fprintf(fid,'- **%s:** PS/FK `Ec=%.6g`; flat-normal `Ec=%.6g`; no-Jacobian `Ec=%.6g`; first aperture below half the PS error: **%s**.\n',d.case_name,ps.complex_l2,complex(1),complex(2),Ltxt);end
fprintf(fid,'\nThe physical importance ordering is determined from the tables rather than assumed. A term whose ablation error is far below PS/FK cannot explain the phase-screen failure by itself. The aperture needed to beat the phase screen measures the minimum nonlocal support for a reduced operator under this receiver-centered truncation.\n\n');
fprintf(fid,'The flat-surface ablation shows that physical surface height/path phase is essential, but it does **not** mean PE omits height phase: PE already contains the local `2*k*eta` proxy. The residual diagnosis is that this local proxy does not reproduce the complete surface-to-receiver path phase and nonlocal coherent integration. Curved-normal and Jacobian corrections are too small to close that residual on their own.\n\n');
fprintf(fid,'A reduced short-range Kirchhoff integral is the preferred next prototype if a finite `L` consistently beats the phase screen across all three cases. A phase-screen-plus-small-kernel form should only be implemented after extracting a case-independent kernel from this aperture study; directly inserting the measured `u_FK-u_PS` would be tautological. SSA remains out of scope because complete Full Kirchhoff already agrees with BIE in the tested envelope.\n\n');
fprintf(fid,'Artifacts: `%s`, `%s`, `%s`, `%s`; figures: `results/validation/pe_bie_full_kirchhoff_reduced_stage1/figures/`.\n',local_rel(root,mat),local_rel(root,summary_csv),local_rel(root,radius_csv),local_rel(root,sweep_csv));
end
function p=local_rel(root,x),p=strrep(x,[root filesep],'');p=strrep(p,filesep,'/');end
