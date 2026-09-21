function validation = validate_pe_surface_kz_aware_operator()
%VALIDATE_PE_SURFACE_KZ_AWARE_OPERATOR G2 validation-only Model-1 test.

root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup_vertical_project();
out_dir=fullfile(root,'results','validation','pe_surface_operator_bie_reference', ...
    'G2_kz_aware');
if ~exist(out_dir,'dir'),mkdir(out_dir);end

s0=local_load_validation(fullfile(root,'results','validation', ...
    'pe_bellhop_controlled_comparison','stage0','stage0_validation.mat'));
r1=local_load_validation(fullfile(root,'results','validation', ...
    'pe_bellhop_helmholtz_bie_reference','R1_flat','R1_flat_validation.mat'));
g1=load(fullfile(root,'results','validation','pe_surface_operator_bie_reference', ...
    'G1_normal_approximation','G1_normal_approximation_audit.mat'),'audit');
assert(s0.passed && r1.passed && g1.audit.passed,'G2 prerequisite gate failed.');

paths={ ...
    fullfile(root,'results','validation','pe_bellhop_helmholtz_bie_reference', ...
        'R3_weak_three_way','R3_weak_three_way_validation.mat'); ...
    fullfile(root,'results','validation','pe_bellhop_helmholtz_bie_reference', ...
        'R4_region_II','R4_region_II_validation.mat'); ...
    fullfile(root,'results','validation','pe_bellhop_helmholtz_bie_reference', ...
        'R5_stronger_height','R5_stronger_height_validation.mat'); ...
    fullfile(root,'results','validation','pe_surface_operator_bie_reference', ...
        'G0_high_K','G0_high_K_validation.mat')};
names=["weak_low_K";"region_II_low_K";"strong_height_low_K";"weak_high_K"];
flat_ref=s0.pe.reflected_field(:); x=s0.pe.x_m(:); fp=local_footprint(x,flat_ref);

flat_cfg=local_cfg(s0.config,zeros(size(x)));
flat_cfg.reflection_model='model0_normal'; flat0=run_pe_1d_surface_reflection_validation(flat_cfg);
flat_cfg.reflection_model='model1_kz_aware'; flat1=run_pe_1d_surface_reflection_validation(flat_cfg);
flat_model_difference=norm(flat1.reflected_field(:)-flat0.reflected_field(:))/ ...
    norm(flat0.reflected_field(:));

n=numel(paths); A=zeros(n,1);K=A;
E0=A;E1=A;EBH=A;phase0=A;phase1=A;phaseBH=A;tl0=A;tl1=A;tlBH=A;
rho0=A;rho1=A;rhoBH=A;aligned0=A;aligned1=A;alignedBH=A;legacy_reproduction=A;
for ii=1:n
    v=local_load_validation(paths{ii});
    assert(v.passed,'Authoritative case did not pass: %s',paths{ii});
    eta=local_surface(x,v.config.amplitude_m,v.config.wavenumber_radpm, ...
        v.config.surface_taper_inner_m,v.config.surface_support_m);
    cfg=local_cfg(s0.config,eta);
    cfg.reflection_model='model0_normal'; p0=run_pe_1d_surface_reflection_validation(cfg);
    cfg.reflection_model='model1_kz_aware'; p1=run_pe_1d_surface_reflection_validation(cfg);
    G0=p0.reflected_field(:)./flat_ref;
    G1=p1.reflected_field(:)./flat_ref;
    legacy_reproduction(ii)=norm(G0(fp.mask)-v.G_PE(fp.mask))/norm(v.G_PE(fp.mask));
    m0=local_metrics(G0,v.G_BIE(:),fp);
    m1=local_metrics(G1,v.G_BIE(:),fp);
    mb=v.metrics.BH_BIE;
    A(ii)=v.config.amplitude_m;K(ii)=v.config.wavenumber_radpm;
    [E0(ii),phase0(ii),tl0(ii),rho0(ii),aligned0(ii)]=local_unpack(m0);
    [E1(ii),phase1(ii),tl1(ii),rho1(ii),aligned1(ii)]=local_unpack(m1);
    [EBH(ii),phaseBH(ii),tlBH(ii),rhoBH(ii),alignedBH(ii)]=local_unpack(mb);
end
cases=table(names,A,K,E0,E1,EBH,phase0,phase1,phaseBH,tl0,tl1,tlBH, ...
    rho0,rho1,rhoBH,aligned0,aligned1,alignedBH,legacy_reproduction, ...
    E0./E1,phase0./phase1, ...
    'VariableNames',{'case_name','A_m','K_radpm','Model0_E_G','Model1_E_G', ...
    'Bellhop_E_G','Model0_phase_rms_rad','Model1_phase_rms_rad', ...
    'Bellhop_phase_rms_rad','Model0_TL_rms_db','Model1_TL_rms_db', ...
    'Bellhop_TL_rms_db','Model0_rho_shape','Model1_rho_shape', ...
    'Bellhop_rho_shape','Model0_E_aligned','Model1_E_aligned', ...
    'Bellhop_E_aligned','legacy_reproduction_error','E_improvement_factor', ...
    'phase_improvement_factor'});

checks=struct( ...
    'flat_not_degraded',flat_model_difference<=1e-10, ...
    'legacy_reproduced',max(legacy_reproduction)<=1e-12, ...
    'weak_preserved',E1(1)<=E0(1), ...
    'region_II_significant',E1(2)<=0.5*E0(2) && phase1(2)<=0.5*phase0(2), ...
    'strong_height_significant',E1(3)<=0.75*E0(3) && phase1(3)<=0.75*phase0(3), ...
    'consistent_all_cases',all(E1<E0) && all(phase1<phase0), ...
    'finite',all(isfinite(cases{:,2:end}),'all'));
checks.all=all(cell2mat(struct2cell(checks)));
high_K_limited=E1(4)>0.75*E0(4) || phase1(4)>0.75*phase0(4);
if checks.all
    conclusion='NORMAL_APPROXIMATION_CONFIRMED';
else
    conclusion='IMPROVEMENT_NOT_CONFIRMED';
end
validation=struct('schema_version','1.0.0','stage','G2_kz_aware', ...
    'model0','R0*exp(+i*2*k*eta)*psi_inc', ...
    'model1','componentwise R0*exp(+i*2*kz(kx)*eta)*Psi_inc(kx)', ...
    'no_fitted_parameters',true,'flat_model_difference',flat_model_difference, ...
    'flat_BIE_metrics',r1.main.metrics,'flat_Bellhop_metrics',s0.metrics.flat_internal, ...
    'cases',cases,'checks',checks,'high_K_limited',high_K_limited, ...
    'conclusion',conclusion,'passed',checks.all);
mat_file=fullfile(out_dir,'G2_kz_aware_validation.mat');
csv_file=fullfile(out_dir,'G2_kz_aware_cases.csv');
save(mat_file,'validation','-v7.3');writetable(cases,csv_file);
report_file=fullfile(root,'reports','pe_surface_kz_aware_G2_validation.md');
local_write_report(report_file,validation,mat_file,csv_file);
fprintf('G2 %s high_K_limited=%d\n',conclusion,high_K_limited);disp(cases);
if ~validation.passed,error('G2 gate failed; G3 remains locked.');end
end

function v=local_load_validation(path)
assert(exist(path,'file')==2,'Missing authoritative artifact: %s',path);
q=load(path,'validation');v=q.validation;
end

function cfg=local_cfg(base,eta)
cfg=struct('frequency_hz',4000,'c0_mps',1500,'xw_m',192.1875, ...
    'nx',984,'z_tx_m',100,'z_rx_m',3,'sigma_src_m',0.3, ...
    'surface_elevation_x_m',eta(:).','surface_reflect_coeff',-1, ...
    'step_m',base.pe_step_m,'x_rx_m',0);
end

function eta=local_surface(x,A,K,inner,outer)
r=abs(x);chi=zeros(size(x));chi(r<=inner)=1;mid=r>inner&r<outer;
t=(r(mid)-inner)/(outer-inner);chi(mid)=1-10*t.^3+15*t.^4-6*t.^5;
eta=A*sin(K*x).*chi;
end

function fp=local_footprint(x,field)
e=abs(field(:)).^2;[~,i0]=min(abs(x));r=abs(x-x(i0));[rs,ord]=sort(r);
c=cumsum(e(ord))/sum(e);r99=rs(find(c>=.99,1));mask=r<=r99;
fp=struct('mask',mask,'weights',e(mask)/sum(e(mask)));
end

function m=local_metrics(a,b,fp)
a=a(:);b=b(:);mask=fp.mask;w=fp.weights;
phase=angle(a.*conj(b));tl=20*log10(max(abs(a),realmin)./max(abs(b),realmin));
S=sum(w.*a(mask).*conj(b(mask)));den=sqrt(max(sum(w.*abs(a(mask)).^2)*sum(w.*abs(b(mask)).^2),realmin));
c=S/den;phi0=angle(S);
m=struct('E_G',sqrt(sum(w.*abs(a(mask)-b(mask)).^2)/max(sum(w.*abs(b(mask)).^2),realmin)), ...
    'phase_rms_rad',sqrt(sum(w.*phase(mask).^2)), ...
    'tl_rms_db',sqrt(sum(w.*tl(mask).^2)), ...
    'rho_shape',abs(c),'E_aligned',sqrt(sum(w.*abs(a(mask)-exp(1i*phi0)*b(mask)).^2)/ ...
    max(sum(w.*abs(a(mask)).^2),realmin)));
end

function [E,p,t,r,a]=local_unpack(m)
E=m.E_G;p=m.phase_rms_rad;t=m.tl_rms_db;r=m.rho_shape;a=m.E_aligned;
end

function local_write_report(path,v,mat_file,csv_file)
fid=fopen(path,'w','n','UTF-8');if fid<0,error('Cannot write %s.',path);end
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# PE surface kz-aware Model-1 G2 validation\n\n');
fprintf(fid,'Status: **%s**; high-K limited: **%d**.\n\n',v.conclusion,v.high_K_limited);
fprintf(fid,'- Model-0: `%s`.\n- Model-1: `%s`.\n',v.model0,v.model1);
fprintf(fid,'- No BIE-fitted parameter; production PE is unchanged.\n');
fprintf(fid,'- Flat Model-0/Model-1 relative field difference: `%.9g`.\n',v.flat_model_difference);
fprintf(fid,'- Independent flat BIE error: `%.9g`; flat Bellhop validation error: `%.9g`.\n\n', ...
    v.flat_BIE_metrics.complex_l2_m99,v.flat_Bellhop_metrics.l2_m99);
fprintf(fid,'| case | M0 E | M1 E | BH E | M0 phase | M1 phase | BH phase | M0 TL | M1 TL | BH TL | M0 rho | M1 rho | M0 aligned | M1 aligned | E gain | phase gain |\n');
fprintf(fid,'|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n');
t=v.cases;
for ii=1:height(t)
    fprintf(fid,'| %s | %.6g | %.6g | %.6g | %.6g | %.6g | %.6g | %.6g | %.6g | %.6g | %.8g | %.8g | %.6g | %.6g | %.4g | %.4g |\n', ...
        t.case_name(ii),t.Model0_E_G(ii),t.Model1_E_G(ii),t.Bellhop_E_G(ii), ...
        t.Model0_phase_rms_rad(ii),t.Model1_phase_rms_rad(ii),t.Bellhop_phase_rms_rad(ii), ...
        t.Model0_TL_rms_db(ii),t.Model1_TL_rms_db(ii),t.Bellhop_TL_rms_db(ii), ...
        t.Model0_rho_shape(ii),t.Model1_rho_shape(ii),t.Model0_E_aligned(ii), ...
        t.Model1_E_aligned(ii),t.E_improvement_factor(ii),t.phase_improvement_factor(ii));
end
fprintf(fid,'\n## Gates\n\n');names=fieldnames(v.checks);
for ii=1:numel(names),fprintf(fid,'- `%s`: %s\n',names{ii},string(v.checks.(names{ii})));end
fprintf(fid,'\nIf high-K remains limited while G2 passes, the Goal proceeds to G3 local-slope coupling without fitting Model-1.\n\n');
fprintf(fid,'Artifacts: `%s`, `%s`.\n',mat_file,csv_file);
end
