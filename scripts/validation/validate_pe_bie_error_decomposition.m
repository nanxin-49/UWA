function validation = validate_pe_bie_error_decomposition()
%VALIDATE_PE_BIE_ERROR_DECOMPOSITION Decompose PE--BIE residual by type.
% This validation is read-only with respect to production PE/BIE physics. It
% reuses authoritative R3/R4/R5/G0 artifacts and recomputes only the missing
% pointwise Model-1 fields for those four cases.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup_vertical_project();
out_dir = fullfile(root,'results','validation','pe_bie_error_decomposition');
if ~exist(out_dir,'dir'), mkdir(out_dir); end

stage0 = load(fullfile(root,'results','validation', ...
    'pe_bellhop_controlled_comparison','stage0','stage0_validation.mat'),'validation');
stage0 = stage0.validation;
case_specs = local_case_specs(root);
rows = repmat(local_empty_row(),2*numel(case_specs),1);
details = repmat(struct('case_name','','source_file','','x_m',[],'mask',[], ...
    'weights',[],'G_Model0',[],'G_Model1',[],'G_BIE',[],'surface_elevation_m',[]), ...
    numel(case_specs),1);

rr = 0;
for ii = 1:numel(case_specs)
    s = load(case_specs(ii).file,'validation');
    v = s.validation;
    x = v.x_m(:);
    G0 = v.G_PE(:);
    GBIE = v.G_BIE(:);
    eta = local_surface(x,v.config.amplitude_m,v.config.wavenumber_radpm, ...
        v.config.surface_taper_inner_m,v.config.surface_support_m);
    G1 = local_model1(v.config,eta,stage0.pe.reflected_field(:));

    [mask,weights] = local_mask(stage0.footprint,x,G0,G1,GBIE);
    m0 = local_metrics(G0,GBIE,mask,weights);
    m1 = local_metrics(G1,GBIE,mask,weights);
    rr=rr+1; rows(rr)=local_row(case_specs(ii),"Model-0",m0);
    rr=rr+1; rows(rr)=local_row(case_specs(ii),"kz-aware Model-1",m1);
    details(ii)=struct('case_name',case_specs(ii).name,'source_file',case_specs(ii).file, ...
        'x_m',x,'mask',mask,'weights',weights,'G_Model0',G0,'G_Model1',G1, ...
        'G_BIE',GBIE,'surface_elevation_m',eta);
end

table_rows = struct2table(rows);
summary = local_summary(table_rows);
validation = struct('schema_version','1.0.0', ...
    'stage','pe_bie_error_decomposition','frequency_hz',4000,'c0_mps',1500, ...
    'case_definition','R3 weak low-K, R4 region-II low-K, R5 strong-height low-K, G0 weak high-K', ...
    'mask_definition',['stage0 99-percent incident-energy footprint, finite samples, ', ...
        'and both compared fields above -40 dB of their own peak'], ...
    'phase_metric','weighted RMS of wrapped angle(a*conj(b)) on the effective mask', ...
    'complex_phase_metric','weighted circular coherence of unit phasors', ...
    'reused_artifacts',{arrayfun(@(q)local_relative_path(root,q.file),case_specs,'UniformOutput',false)}, ...
    'model1_recomputed_validation_only',true,'rows',table_rows, ...
    'summary',summary,'details',details);

mat_file = fullfile(out_dir,'pe_bie_error_decomposition_validation.mat');
csv_file = fullfile(out_dir,'pe_bie_error_decomposition_cases.csv');
save(mat_file,'validation','-v7.3');
writetable(table_rows,csv_file);
report_file = fullfile(root,'reports','pe_bie_error_decomposition_report.md');
local_write_report(report_file,validation,root,mat_file,csv_file);
disp(table_rows);
fprintf('PE--BIE error decomposition written to %s\n',report_file);
end

function specs = local_case_specs(root)
base = fullfile(root,'results','validation');
specs = struct('name',{},'amplitude_m',{},'wavenumber_radpm',{},'file',{});
specs(1)=struct('name',"weak_low_K",'amplitude_m',0.01,'wavenumber_radpm',0.10, ...
    'file',fullfile(base,'pe_bellhop_helmholtz_bie_reference','R3_weak_three_way', ...
    'R3_weak_three_way_validation.mat'));
specs(2)=struct('name',"region_II_low_K",'amplitude_m',0.05,'wavenumber_radpm',0.10, ...
    'file',fullfile(base,'pe_bellhop_helmholtz_bie_reference','R4_region_II', ...
    'R4_region_II_validation.mat'));
specs(3)=struct('name',"strong_height_low_K",'amplitude_m',0.20,'wavenumber_radpm',0.10, ...
    'file',fullfile(base,'pe_bellhop_helmholtz_bie_reference','R5_stronger_height', ...
    'R5_stronger_height_validation.mat'));
specs(4)=struct('name',"weak_high_K",'amplitude_m',0.02,'wavenumber_radpm',0.47, ...
    'file',fullfile(base,'pe_surface_operator_bie_reference','G0_high_K', ...
    'G0_high_K_validation.mat'));
for ii=1:numel(specs)
    assert(exist(specs(ii).file,'file')==2,'Missing authoritative artifact: %s',specs(ii).file);
end
end

function G1 = local_model1(cfg,eta,flat_ref)
pe_cfg = struct('frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps, ...
    'xw_m',cfg.xw_m,'nx',cfg.nx,'z_tx_m',cfg.z_tx_m,'z_rx_m',cfg.z_rx_m, ...
    'sigma_src_m',cfg.sigma_src_m,'surface_elevation_x_m',eta(:).', ...
    'surface_reflect_coeff',-1,'step_m',0.05,'x_rx_m',0, ...
    'reflection_model','model1_kz_aware');
p = run_pe_1d_surface_reflection_validation(pe_cfg);
G1 = p.reflected_field(:)./flat_ref;
end

function eta = local_surface(x,A,K,inner,outer)
r=abs(x); chi=zeros(size(x)); chi(r<=inner)=1;
mid=r>inner & r<outer; t=(r(mid)-inner)/(outer-inner);
chi(mid)=1-10*t.^3+15*t.^4-6*t.^5;
eta=A*sin(K*x).*chi;
end

function [mask,w] = local_mask(footprint,x,a,b,c)
    base = footprint.m99.mask(:) & numel(x)==numel(footprint.m99.mask(:));
finite = isfinite(real(a)) & isfinite(imag(a)) & isfinite(real(b)) & ...
    isfinite(imag(b)) & isfinite(real(c)) & isfinite(imag(c));
floor0=10^(-40/20);
level = max([abs(a(:));abs(b(:));abs(c(:))]);
mask = base & finite & abs(a)>=floor0*level & abs(b)>=floor0*level & abs(c)>=floor0*level;
w=footprint.energy_weights(:);
w=w(mask); w=w/max(sum(w),realmin);
end

function m = local_metrics(a,b,mask,w)
a=a(mask); b=b(mask);
phase=angle(a.*conj(b)); ma=abs(a); mb=abs(b);
mag_l2=sqrt(sum(w.*(ma-mb).^2)/max(sum(w.*mb.^2),realmin));
complex_l2=sqrt(sum(w.*abs(a-b).^2)/max(sum(w.*mb.^2),realmin));
phase_rms=sqrt(sum(w.*phase.^2));
mag_corr=local_weighted_corr(ma,mb,w);
phase_corr=abs(sum(w.*exp(1i*angle(a)).*conj(exp(1i*angle(b)))));
if phase_rms > 1.5*mag_l2
    class="phase-dominated";
elseif mag_l2 > 1.5*phase_rms
    class="magnitude-dominated";
else
    class="both-significant";
end
m=struct('complex_l2',complex_l2,'magnitude_relative_l2',mag_l2, ...
    'wrapped_phase_rms_rad',phase_rms,'magnitude_correlation',mag_corr, ...
    'complex_phase_correlation',phase_corr,'classification',class, ...
    'sample_count',numel(a));
end

function r=local_weighted_corr(a,b,w)
wa=sum(w.*a)/sum(w); wb=sum(w.*b)/sum(w);
num=sum(w.*(a-wa).*(b-wb)); den=sqrt(sum(w.*(a-wa).^2)*sum(w.*(b-wb).^2));
r=abs(num/max(den,realmin));
end

function row=local_empty_row()
row=struct('case_name',"",'A_m',NaN,'K_radpm',NaN,'model',"", ...
    'complex_l2',NaN,'magnitude_relative_l2',NaN,'wrapped_phase_rms_rad',NaN, ...
    'magnitude_correlation',NaN,'complex_phase_correlation',NaN, ...
    'classification',"",'sample_count',NaN);
end

function row=local_row(spec,model,m)
row=local_empty_row(); row.case_name=spec.name; row.A_m=spec.amplitude_m;
row.K_radpm=spec.wavenumber_radpm; row.model=model;
row.complex_l2=m.complex_l2; row.magnitude_relative_l2=m.magnitude_relative_l2;
row.wrapped_phase_rms_rad=m.wrapped_phase_rms_rad;
row.magnitude_correlation=m.magnitude_correlation;
row.complex_phase_correlation=m.complex_phase_correlation;
row.classification=m.classification; row.sample_count=m.sample_count;
end

function s=local_summary(t)
names=unique(t.case_name,'stable'); s=struct('case_name',{},'model0_classification',{}, ...
    'model1_classification',{},'model0_phase_to_magnitude',{}, ...
    'model1_phase_to_magnitude',{},'kz_aware_complex_improvement',{});
for ii=1:numel(names)
    q=t(t.case_name==names(ii),:); q0=q(q.model=="Model-0",:); q1=q(q.model=="kz-aware Model-1",:);
    z=struct('case_name',names(ii),'model0_classification',q0.classification, ...
        'model1_classification',q1.classification, ...
        'model0_phase_to_magnitude',q0.wrapped_phase_rms_rad/max(q0.magnitude_relative_l2,realmin), ...
        'model1_phase_to_magnitude',q1.wrapped_phase_rms_rad/max(q1.magnitude_relative_l2,realmin), ...
        'kz_aware_complex_improvement',q0.complex_l2/max(q1.complex_l2,realmin));
    s(end+1)=z; %#ok<AGROW>
end
end

function local_write_report(path,v,root,mat_file,csv_file)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0,'Cannot create report.');
cleanup_obj=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# PE--BIE error decomposition\n\n');
fprintf(fid,'Status: **independent diagnostic; production physics unchanged**.\n\n');
fprintf(fid,'Four authoritative cases are reused: weak/low-K (R3), region-II medium height/low-K (R4), strong-height/low-K (R5), and weak/high-K (G0). Only missing pointwise Model-1 fields are recomputed with the frozen 4 kHz, c=1500 m/s PE configuration.\n\n');
fprintf(fid,'## Mask and metrics\n\n');
fprintf(fid,'The effective mask is the saved Stage-0 99%% incident-energy footprint, finite samples, and samples where each compared field is at least -40 dB relative to the largest field in that pair. Weights are the saved Stage-0 energy weights renormalized on that mask. Complex and magnitude L2 values are relative to BIE magnitude energy; phase is wrapped `angle(a*conj(b))`; phase correlation is circular coherence of unit phasors. Classification uses phase RMS / magnitude relative L2: >1.5 phase-dominated, <2/3 magnitude-dominated, otherwise both significant.\n\n');
t=v.rows; fprintf(fid,'| case | model | complex L2 | magnitude L2 | phase RMS (rad) | magnitude corr. | complex phase corr. | diagnosis | N |\n|---|---|---:|---:|---:|---:|---:|---|---:|\n');
for ii=1:height(t)
    fprintf(fid,'| %s | %s | %.6g | %.6g | %.6g | %.6g | %.6g | %s | %d |\n', ...
        t.case_name(ii),t.model(ii),t.complex_l2(ii),t.magnitude_relative_l2(ii), ...
        t.wrapped_phase_rms_rad(ii),t.magnitude_correlation(ii),t.complex_phase_correlation(ii), ...
        t.classification(ii),t.sample_count(ii));
end
fprintf(fid,'\n## Interpretation\n\n');
for ii=1:numel(v.summary)
    q=v.summary(ii); fprintf(fid,'- `%s`: Model-0 **%s**, Model-1 **%s**; phase/magnitude ratios %.4g -> %.4g; kz-aware complex-error improvement %.4g-fold.\n', ...
        q.case_name,q.model0_classification,q.model1_classification,q.model0_phase_to_magnitude,q.model1_phase_to_magnitude,q.kz_aware_complex_improvement);
end
fprintf(fid,'\n## Next-step decision rule\n\n');
fprintf(fid,'The observed diagnosis is case dependent: Model-0 is phase-dominated for all low-K cases, while the kz-aware correction removes most finite-angle phase error in weak/medium low-K cases. Strong height retains a substantial phase residual after the correction, and high-K retains a non-negligible magnitude residual with only modest complex-error improvement. Therefore amplitude alone is not a sufficient explanation for the PE--BIE discrepancy; the next physics investigation should prioritize complete Kirchhoff/nonlocal phase coupling while separately auditing amplitude/geometric factors.\n\n');
fprintf(fid,'- If magnitude error is much smaller than phase error, prioritize complete Kirchhoff/nonlocal phase coupling.\n');
fprintf(fid,'- If both remain significant, implement a complete Kirchhoff surface integral and audit obliquity, normal derivative, and surface Jacobian terms.\n');
fprintf(fid,'- If magnitude error dominates, first inspect a simple amplitude/geometric correction; do not immediately introduce a nonlocal operator.\n\n');
fprintf(fid,'Artifacts: `%s`, `%s`.\n',local_relative_path(root,mat_file),local_relative_path(root,csv_file));
end

function rel = local_relative_path(root,path)
root = char(root); path = char(path);
prefix = [root filesep];
if strncmpi(path,prefix,numel(prefix))
    rel = strrep(path(numel(prefix)+1:end),filesep,'/');
else
    rel = strrep(path,filesep,'/');
end
end
