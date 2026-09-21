function validation=validate_pe_surface_angle_slope_operator()
%VALIDATE_PE_SURFACE_ANGLE_SLOPE_OPERATOR G3 local-specular phase audit.

root=fileparts(fileparts(fileparts(mfilename('fullpath'))));setup_vertical_project();
out_dir=fullfile(root,'results','validation','pe_surface_operator_bie_reference', ...
    'G3_angle_slope');if ~exist(out_dir,'dir'),mkdir(out_dir);end
s0=local_load(fullfile(root,'results','validation','pe_bellhop_controlled_comparison', ...
    'stage0','stage0_validation.mat'));
g2=local_load(fullfile(root,'results','validation','pe_surface_operator_bie_reference', ...
    'G2_kz_aware','G2_kz_aware_validation.mat'));
assert(s0.passed&&g2.passed&&g2.high_K_limited,'G3 prerequisite gate failed.');
paths={ ...
    fullfile(root,'results','validation','pe_bellhop_helmholtz_bie_reference','R3_weak_three_way','R3_weak_three_way_validation.mat'); ...
    fullfile(root,'results','validation','pe_bellhop_helmholtz_bie_reference','R4_region_II','R4_region_II_validation.mat'); ...
    fullfile(root,'results','validation','pe_bellhop_helmholtz_bie_reference','R5_stronger_height','R5_stronger_height_validation.mat'); ...
    fullfile(root,'results','validation','pe_surface_operator_bie_reference','G0_high_K','G0_high_K_validation.mat')};
names=["weak_low_K";"region_II_low_K";"strong_height_low_K";"weak_high_K"];
x=s0.pe.x_m(:);flat=s0.pe.reflected_field(:);fp=local_footprint(x,flat);

cfg=local_cfg(s0.config,zeros(size(x)),zeros(size(x)));
cfg.reflection_model='model0_normal';f0=run_pe_1d_surface_reflection_validation(cfg);
cfg.reflection_model='model2_angle_slope';f2=run_pe_1d_surface_reflection_validation(cfg);
flat_difference=norm(f2.reflected_field(:)-f0.reflected_field(:))/norm(f0.reflected_field(:));

n=numel(paths);A=zeros(n,1);K=A;E0=A;E1=A;E2=A;EB=A;
p0=A;p1=A;p2=A;pB=A;t0=A;t1=A;t2=A;tB=A;r0=A;r1=A;r2=A;rB=A;
a0=A;a1=A;a2=A;aB=A;nonreturn=A;minimum_kzr=A;
for ii=1:n
    v=local_load(paths{ii});assert(v.passed,'Prerequisite case failed.');
    [eta,slope]=local_surface(x,v.config.amplitude_m,v.config.wavenumber_radpm, ...
        v.config.surface_taper_inner_m,v.config.surface_support_m);
    cfg=local_cfg(s0.config,eta,slope);
    cfg.reflection_model='model0_normal';q0=run_pe_1d_surface_reflection_validation(cfg);
    cfg.reflection_model='model1_kz_aware';q1=run_pe_1d_surface_reflection_validation(cfg);
    cfg.reflection_model='model2_angle_slope';q2=run_pe_1d_surface_reflection_validation(cfg);
    G0=q0.reflected_field(:)./flat;G1=q1.reflected_field(:)./flat;G2=q2.reflected_field(:)./flat;
    m0=local_metrics(G0,v.G_BIE(:),fp);m1=local_metrics(G1,v.G_BIE(:),fp);
    m2=local_metrics(G2,v.G_BIE(:),fp);mb=v.metrics.BH_BIE;
    A(ii)=v.config.amplitude_m;K(ii)=v.config.wavenumber_radpm;
    [E0(ii),p0(ii),t0(ii),r0(ii),a0(ii)]=local_unpack(m0);
    [E1(ii),p1(ii),t1(ii),r1(ii),a1(ii)]=local_unpack(m1);
    [E2(ii),p2(ii),t2(ii),r2(ii),a2(ii)]=local_unpack(m2);
    [EB(ii),pB(ii),tB(ii),rB(ii),aB(ii)]=local_unpack(mb);
    nonreturn(ii)=q2.reflection_meta.nonreturning_component_count;
    minimum_kzr(ii)=q2.reflection_meta.minimum_reflected_kz_radpm;
end
cases=table(names,A,K,E0,E1,E2,EB,p0,p1,p2,pB,t0,t1,t2,tB,r0,r1,r2,rB, ...
    a0,a1,a2,aB,E1./E2,p1./p2,minimum_kzr,nonreturn, ...
    'VariableNames',{'case_name','A_m','K_radpm','Model0_E_G','Model1_E_G', ...
    'Model2_E_G','Bellhop_E_G','Model0_phase_rms_rad','Model1_phase_rms_rad', ...
    'Model2_phase_rms_rad','Bellhop_phase_rms_rad','Model0_TL_rms_db', ...
    'Model1_TL_rms_db','Model2_TL_rms_db','Bellhop_TL_rms_db', ...
    'Model0_rho_shape','Model1_rho_shape','Model2_rho_shape','Bellhop_rho_shape', ...
    'Model0_E_aligned','Model1_E_aligned','Model2_E_aligned','Bellhop_E_aligned', ...
    'M1_to_M2_E_gain','M1_to_M2_phase_gain','minimum_reflected_kz_radpm', ...
    'nonreturning_component_count'});
checks=struct('flat_not_degraded',flat_difference<=1e-10, ...
    'low_K_not_degraded',all(E2(1:3)<=1.1*E1(1:3)), ...
    'high_K_significant',E2(4)<=0.75*E1(4)&&p2(4)<=0.75*p1(4), ...
    'high_K_better_than_Model0',E2(4)<E0(4)&&p2(4)<p0(4), ...
    'returning_branch',all(nonreturn==0)&&all(minimum_kzr>0), ...
    'finite',all(isfinite(cases{:,2:end}),'all'));
checks.all=all(cell2mat(struct2cell(checks)));
if checks.all,conclusion='ANGLE_SLOPE_MECHANISM_CONFIRMED';else,conclusion='NONLOCAL_EFFECT_REQUIRED';end
validation=struct('schema_version','1.0.0','stage','G3_angle_slope', ...
    'formula',['kzr=((1-s^2)kzi+2*s*kxi)/(1+s^2); ', ...
    'R0*exp(+i*(kzi+kzr)*eta) per incident component'], ...
    'no_fitted_parameters',true,'flat_model_difference',flat_difference, ...
    'cases',cases,'checks',checks,'conclusion',conclusion,'passed',checks.all);
mat_file=fullfile(out_dir,'G3_angle_slope_validation.mat');
csv_file=fullfile(out_dir,'G3_angle_slope_cases.csv');save(mat_file,'validation','-v7.3');writetable(cases,csv_file);
report_file=fullfile(root,'reports','pe_surface_angle_slope_G3_validation.md');
local_report(report_file,validation,mat_file,csv_file);
fprintf('G3 %s\n',conclusion);disp(cases);
if ~validation.passed,error('G3 gate failed: nonlocal-effect review is required before G5.');end
end

function v=local_load(path)
assert(exist(path,'file')==2,'Missing artifact: %s',path);q=load(path,'validation');v=q.validation;
end
function cfg=local_cfg(base,eta,slope)
cfg=struct('frequency_hz',4000,'c0_mps',1500,'xw_m',192.1875,'nx',984, ...
    'z_tx_m',100,'z_rx_m',3,'sigma_src_m',0.3,'surface_elevation_x_m',eta(:).', ...
    'surface_slope_x',slope(:).','surface_reflect_coeff',-1, ...
    'step_m',base.pe_step_m,'x_rx_m',0);
end
function [eta,deta]=local_surface(x,A,K,inner,outer)
r=abs(x);chi=zeros(size(x));dchi=zeros(size(x));inside=r<=inner;chi(inside)=1;
mid=r>inner&r<outer;t=(r(mid)-inner)/(outer-inner);
chi(mid)=1-10*t.^3+15*t.^4-6*t.^5;
dchi(mid)=(-30*t.^2+60*t.^3-30*t.^4)/(outer-inner).*sign(x(mid));
eta=A*sin(K*x).*chi;deta=A*(K*cos(K*x).*chi+sin(K*x).*dchi);
end
function fp=local_footprint(x,f)
e=abs(f(:)).^2;[~,i0]=min(abs(x));r=abs(x-x(i0));[rs,o]=sort(r);c=cumsum(e(o))/sum(e);
mask=r<=rs(find(c>=.99,1));fp=struct('mask',mask,'weights',e(mask)/sum(e(mask)));
end
function m=local_metrics(a,b,fp)
a=a(:);b=b(:);mask=fp.mask;w=fp.weights;ph=angle(a.*conj(b));tl=20*log10(max(abs(a),realmin)./max(abs(b),realmin));
S=sum(w.*a(mask).*conj(b(mask)));den=sqrt(max(sum(w.*abs(a(mask)).^2)*sum(w.*abs(b(mask)).^2),realmin));c=S/den;phi=angle(S);
m=struct('E_G',sqrt(sum(w.*abs(a(mask)-b(mask)).^2)/max(sum(w.*abs(b(mask)).^2),realmin)), ...
    'phase_rms_rad',sqrt(sum(w.*ph(mask).^2)),'tl_rms_db',sqrt(sum(w.*tl(mask).^2)), ...
    'rho_shape',abs(c),'E_aligned',sqrt(sum(w.*abs(a(mask)-exp(1i*phi)*b(mask)).^2)/max(sum(w.*abs(a(mask)).^2),realmin)));
end
function [E,p,t,r,a]=local_unpack(m)
E=m.E_G;p=m.phase_rms_rad;t=m.tl_rms_db;r=m.rho_shape;a=m.E_aligned;
end
function local_report(path,v,mat_file,csv_file)
fid=fopen(path,'w','n','UTF-8');
if fid<0,error('Cannot write report.');end
c=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# PE surface angle--slope Model-2 G3 validation\n\nStatus: **%s**.\n\n',v.conclusion);
fprintf(fid,'- Formula: `%s`.\n- No fitted parameter; production PE/marching unchanged.\n',v.formula);
fprintf(fid,'- Flat Model-0/Model-2 relative difference: `%.9g`.\n\n',v.flat_model_difference);
fprintf(fid,'| case | M0 E | M1 E | M2 E | BH E | M0 phase | M1 phase | M2 phase | BH phase | M1/M2 E | M1/M2 phase | min kzr | nonreturn |\n');
fprintf(fid,'|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n');t=v.cases;
for ii=1:height(t)
    fprintf(fid,'| %s | %.6g | %.6g | %.6g | %.6g | %.6g | %.6g | %.6g | %.6g | %.4g | %.4g | %.6g | %d |\n', ...
        t.case_name(ii),t.Model0_E_G(ii),t.Model1_E_G(ii),t.Model2_E_G(ii),t.Bellhop_E_G(ii), ...
        t.Model0_phase_rms_rad(ii),t.Model1_phase_rms_rad(ii),t.Model2_phase_rms_rad(ii), ...
        t.Bellhop_phase_rms_rad(ii),t.M1_to_M2_E_gain(ii),t.M1_to_M2_phase_gain(ii), ...
        t.minimum_reflected_kz_radpm(ii),t.nonreturning_component_count(ii));
end
fprintf(fid,'\n## Gates\n\n');n=fieldnames(v.checks);
for ii=1:numel(n)
    fprintf(fid,'- `%s`: %s\n',n{ii},string(v.checks.(n{ii})));
end
fprintf(fid,'\nArtifacts: `%s`, `%s`.\n',mat_file,csv_file);
end
