function audit = validate_pe_surface_normal_approximation_audit()
%VALIDATE_PE_SURFACE_NORMAL_APPROXIMATION_AUDIT G1 angle-phase mechanism audit.
%   Read-only with respect to the accepted PE/Bellhop/BIE artifacts. The
%   componentwise finite-angle field is a diagnostic prediction only; it is
%   not exposed as a PE reflection model by this stage.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup_vertical_project();
out_dir = fullfile(root,'results','validation', ...
    'pe_surface_operator_bie_reference','G1_normal_approximation');
if ~exist(out_dir,'dir'), mkdir(out_dir); end

case_files = { ...
    fullfile(root,'results','validation','pe_bellhop_helmholtz_bie_reference', ...
        'R3_weak_three_way','R3_weak_three_way_validation.mat'); ...
    fullfile(root,'results','validation','pe_bellhop_helmholtz_bie_reference', ...
        'R4_region_II','R4_region_II_validation.mat'); ...
    fullfile(root,'results','validation','pe_bellhop_helmholtz_bie_reference', ...
        'R5_stronger_height','R5_stronger_height_validation.mat'); ...
    fullfile(root,'results','validation','pe_surface_operator_bie_reference', ...
        'G0_high_K','G0_high_K_validation.mat')};
case_names = ["weak_low_K";"region_II_low_K";"strong_height_low_K";"weak_high_K"];

stage0_file = fullfile(root,'results','validation', ...
    'pe_bellhop_controlled_comparison','stage0','stage0_validation.mat');
assert(exist(stage0_file,'file')==2,'Authoritative Stage-0 artifact is missing.');
s0 = load(stage0_file,'validation'); s0 = s0.validation;
assert(s0.passed,'Authoritative Stage-0 prerequisite did not pass.');
flat_ref = s0.pe.reflected_field(:);
x = s0.pe.x_m(:);
fp = local_footprint(x,flat_ref);

first = load(case_files{1},'validation'); first = first.validation;
pe0 = first.PE;
k = pe0.k0_rad_per_m;
kx = pe0.kx_rad_per_m(:);
kz = sqrt(max(k^2-kx.^2,0));
incident_k = fft(pe0.surface_incident_field(:));
spectral_weight = abs(incident_k).^2;
spectral_weight = spectral_weight/sum(spectral_weight);
beta = 2*(k-kz);
theta = asin(min(abs(kx)/k,1));

spectrum = struct( ...
    'kx_rms_radpm',sqrt(sum(spectral_weight.*kx.^2)), ...
    'theta_rms_deg',rad2deg(sqrt(sum(spectral_weight.*theta.^2))), ...
    'kz_mean_radpm',sum(spectral_weight.*kz), ...
    'kz_over_k_mean',sum(spectral_weight.*kz)/k, ...
    'beta_mean_radpm',sum(spectral_weight.*beta), ...
    'beta_rms_radpm',sqrt(sum(spectral_weight.*beta.^2)), ...
    'abs_kx_quantiles_radpm',local_weighted_quantile(abs(kx),spectral_weight,[0.5 0.9 0.95 0.99]), ...
    'abs_theta_quantiles_deg',rad2deg(local_weighted_quantile(theta,spectral_weight,[0.5 0.9 0.95 0.99])));

ncase = numel(case_files);
A = zeros(ncase,1); K = A; eta_rms = A; observed_phase = A;
mean_kz_prediction = A; component_prediction = A;
angle_correction_rms = A; phase_correlation = A; phase_projection = A;
angle_model_phase_residual = A; angle_model_E_G = A; improvement_factor = A;
for ii = 1:ncase
    q = load(case_files{ii},'validation'); v = q.validation;
    assert(v.passed,'Authoritative prerequisite failed: %s',case_files{ii});
    assert(max(abs(v.x_m(:)-x))<=1e-12,'Receiver grid mismatch in %s.',case_files{ii});
    A(ii) = v.config.amplitude_m;
    K(ii) = v.config.wavenumber_radpm;
    eta = local_surface(x,A(ii),K(ii),v.config.surface_taper_inner_m, ...
        v.config.surface_support_m);
    w = fp.weights;
    eta_rms(ii) = sqrt(sum(w.*eta(fp.mask).^2));
    mean_kz_prediction(ii) = spectrum.beta_mean_radpm*eta_rms(ii);
    component_prediction(ii) = spectrum.beta_rms_radpm*eta_rms(ii);

    angle_surface = local_angle_aware_surface_reflection( ...
        v.PE.surface_incident_field(:),eta,kz);
    angle_receiver = local_uniform_march(angle_surface,kz,k, ...
        v.config.z_rx_m,v.config.bellhop_step_m);
    G_angle = complex(zeros(size(x)));
    valid = abs(flat_ref)>0;
    G_angle(valid) = angle_receiver(valid)./flat_ref(valid);

    actual = angle(v.G_PE(:).*conj(v.G_BIE(:)));
    predicted = angle(v.G_PE(:).*conj(G_angle));
    observed_phase(ii) = v.metrics.PE_BIE.phase_rms_rad;
    angle_correction_rms(ii) = sqrt(sum(w.*predicted(fp.mask).^2));
    phase_correlation(ii) = sum(w.*actual(fp.mask).*predicted(fp.mask))/ ...
        sqrt(max(sum(w.*actual(fp.mask).^2)*sum(w.*predicted(fp.mask).^2),realmin));
    phase_projection(ii) = sum(w.*actual(fp.mask).*predicted(fp.mask))/ ...
        max(sum(w.*predicted(fp.mask).^2),realmin);
    m = local_metrics(G_angle,v.G_BIE(:),fp);
    angle_model_phase_residual(ii) = m.phase_rms_rad;
    angle_model_E_G(ii) = m.E_G;
    improvement_factor(ii) = observed_phase(ii)/max(m.phase_rms_rad,realmin);
end

cases = table(case_names,A,K,A.*K,A.*K.^2,eta_rms,observed_phase, ...
    mean_kz_prediction,component_prediction,angle_correction_rms, ...
    phase_correlation,phase_projection,angle_model_phase_residual, ...
    angle_model_E_G,improvement_factor, ...
    'VariableNames',{'case_name','A_m','K_radpm','max_slope_nominal', ...
    'max_curvature_nominal_per_m','eta_rms_weighted_m','PE_BIE_phase_rms_rad', ...
    'mean_kz_phase_prediction_rad','component_phase_prediction_rad', ...
    'angle_correction_rms_rad','phase_error_correlation','phase_projection', ...
    'diagnostic_angle_model_phase_residual_rad','diagnostic_angle_model_E_G', ...
    'phase_improvement_factor'});

low_k = 1:3;
checks = struct( ...
    'authoritative_inputs',true, ...
    'finite',all(isfinite(cases{:,2:end}),'all'), ...
    'low_K_direction',median(cases.phase_error_correlation(low_k))>=0.9, ...
    'low_K_major_fraction',median(cases.angle_correction_rms_rad(low_k)./ ...
        cases.PE_BIE_phase_rms_rad(low_k))>=0.5, ...
    'multi_case_improvement',all(cases.phase_improvement_factor>=1.1));
checks.all = all(cell2mat(struct2cell(checks)));
if checks.all
    conclusion = 'NORMAL_APPROXIMATION_MECHANISM_SUPPORTED_FOR_G2_TEST';
else
    conclusion = 'NORMAL_APPROXIMATION_NOT_PRIMARY__STOP_BEFORE_G2';
end

audit = struct('schema_version','1.0.0','stage','G1_normal_approximation', ...
    'convention','exp(+i*2*k*eta) Model-0; diagnostic finite-angle component phase exp(+i*2*kz*eta)', ...
    'spectrum',spectrum,'cases',cases,'checks',checks,'conclusion',conclusion, ...
    'diagnostic_only',true,'passed',checks.all);
mat_file = fullfile(out_dir,'G1_normal_approximation_audit.mat');
csv_file = fullfile(out_dir,'G1_normal_approximation_cases.csv');
save(mat_file,'audit','-v7.3'); writetable(cases,csv_file);
report_file = fullfile(root,'reports','pe_surface_normal_approximation_G1_audit.md');
local_write_report(report_file,audit,mat_file,csv_file);
fprintf('G1 %s\n',conclusion);
disp(cases);
if ~audit.passed
    error('G1 gate failed; G2 remains locked.');
end
end

function eta = local_surface(x,A,K,inner,outer)
r = abs(x); chi = zeros(size(x));
chi(r<=inner) = 1;
mid = r>inner & r<outer;
t = (r(mid)-inner)/(outer-inner);
chi(mid) = 1-10*t.^3+15*t.^4-6*t.^5;
eta = A*sin(K*x).*chi;
end

function reflected = local_angle_aware_surface_reflection(incident,eta,kz)
n = numel(incident);
indices = (0:n-1).';
fourier = exp(1i*2*pi*(indices*indices.')/n)/n;
phase = exp(1i*2*(eta(:)*kz(:).'));
reflected = -sum(fourier.*phase.*fft(incident(:)).',2);
end

function out = local_uniform_march(in,kz,k,distance,step)
nstep = max(1,ceil(abs(distance)/step));
ds = abs(distance)/nstep;
factor = exp(-1i*nstep*ds*(k-kz(:)));
out = ifft(factor.*fft(in(:)));
end

function fp = local_footprint(x,field)
energy = abs(field(:)).^2;
[~,axis_index] = min(abs(x));
radius = abs(x-x(axis_index));
[rs,ord] = sort(radius);
cumulative = cumsum(energy(ord))/sum(energy);
r99 = rs(find(cumulative>=0.99,1));
mask = radius<=r99;
weights = energy(mask)/sum(energy(mask));
fp = struct('mask',mask,'weights',weights,'radius99_m',r99);
end

function m = local_metrics(test,ref,fp)
mask = fp.mask; w = fp.weights;
test = test(:); ref = ref(:);
phase = angle(test.*conj(ref));
m = struct('E_G',sqrt(sum(w.*abs(test(mask)-ref(mask)).^2)/ ...
    max(sum(w.*abs(ref(mask)).^2),realmin)), ...
    'phase_rms_rad',sqrt(sum(w.*phase(mask).^2)));
end

function q = local_weighted_quantile(x,w,p)
[x,ord] = sort(x(:)); w = w(ord); c = cumsum(w)/sum(w);
q = zeros(size(p));
for ii=1:numel(p), q(ii)=x(find(c>=p(ii),1)); end
end

function local_write_report(path,audit,mat_file,csv_file)
fid=fopen(path,'w','n','UTF-8');if fid<0,error('Cannot write %s.',path);end
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# PE surface normal-approximation G1 audit\n\n');
fprintf(fid,'Status: **%s**\n\n',audit.conclusion);
s=audit.spectrum;
fprintf(fid,'## Frozen incident angular spectrum\n\n');
fprintf(fid,'- `kx_rms = %.9g rad/m`; `theta_rms = %.9g deg`.\n',s.kx_rms_radpm,s.theta_rms_deg);
fprintf(fid,'- `<kz>/k = %.12g`; `<2(k-kz)> = %.9g rad/m`; RMS `= %.9g rad/m`.\n', ...
    s.kz_over_k_mean,s.beta_mean_radpm,s.beta_rms_radpm);
fprintf(fid,'- `|kx|` q50/q90/q95/q99: `%.6g / %.6g / %.6g / %.6g rad/m`.\n',s.abs_kx_quantiles_radpm);
fprintf(fid,'- `|theta|` q50/q90/q95/q99: `%.6g / %.6g / %.6g / %.6g deg`.\n\n',s.abs_theta_quantiles_deg);
fprintf(fid,'## Case comparison\n\n');
fprintf(fid,'The componentwise diagnostic uses `exp(+i 2 kz eta)` before the unchanged receiver march. It is not fitted and is not yet a selectable PE model.\n\n');
fprintf(fid,'| case | A | K | PE-BIE phase RMS | mean-kz prediction | component prediction | exact angle correction | correlation | projection | angle residual | E_G | improvement |\n');
fprintf(fid,'|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n');
t=audit.cases;
for ii=1:height(t)
    fprintf(fid,'| %s | %.4g | %.4g | %.6g | %.6g | %.6g | %.6g | %.5g | %.5g | %.6g | %.6g | %.4g |\n', ...
        t.case_name(ii),t.A_m(ii),t.K_radpm(ii),t.PE_BIE_phase_rms_rad(ii), ...
        t.mean_kz_phase_prediction_rad(ii),t.component_phase_prediction_rad(ii), ...
        t.angle_correction_rms_rad(ii),t.phase_error_correlation(ii), ...
        t.phase_projection(ii),t.diagnostic_angle_model_phase_residual_rad(ii), ...
        t.diagnostic_angle_model_E_G(ii),t.phase_improvement_factor(ii));
end
fprintf(fid,'\n## Gate\n\n');
fprintf(fid,['The frozen mechanism gate requires median low-K error-vector correlation ', ...
    '>=0.9, median low-K correction/error RMS >=0.5, and phase-RMS improvement ', ...
    '>=1.1 in every case.\n\n']);
names=fieldnames(audit.checks);
for ii=1:numel(names),fprintf(fid,'- `%s`: %s\n',names{ii},string(audit.checks.(names{ii})));end
fprintf(fid,'\nInterpretation: G1 tests whether the omitted finite-angle phase has the correct scale and error-vector direction across frozen cases. High-K residual after this diagnostic is evidence for a separate slope-coupling test, not permission to fit a coefficient.\n\n');
fprintf(fid,'Artifacts: `%s`, `%s`.\n',mat_file,csv_file);
end
