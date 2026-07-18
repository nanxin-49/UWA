function [deltaG_xy_fm, meta] = sample_kirchhoff_kstat_joint_frequency_vertical( ...
    spec, f_axis_hz, c0_mps, reflect_coeff, n_realizations, rng_seed, mode)
%SAMPLE_KIRCHHOFF_KSTAT_JOINT_FREQUENCY_VERTICAL Joint-frequency K-stat screen.
%   This validation/prototype interface generates zero-mean random phase
%   screens deltaG(x,y,f) on the periodic PM grid. In joint mode it uses the
%   analytic cross-frequency covariance and pseudo-covariance of a Gaussian
%   Kirchhoff phase screen. Independent mode retains only the per-frequency
%   covariance diagonal as an engineering comparison.

arguments
    spec (1,1) struct
    f_axis_hz (:,1) double {mustBeFinite, mustBePositive}
    c0_mps (1,1) double {mustBeFinite, mustBePositive}
    reflect_coeff (1,1) double {mustBeFinite}
    n_realizations (1,1) double {mustBeInteger, mustBePositive}
    rng_seed (1,1) double {mustBeFinite}
    mode (1,:) char {mustBeMember(mode, {'joint', 'independent'})}
end
required = {'W_eta_kstat','dkx_rad_per_m','dky_rad_per_m', ...
    'dx_m','dy_m','nx','ny'};
for ii = 1:numel(required)
    if ~isfield(spec, required{ii})
        error('sample_kirchhoff_kstat_joint_frequency_vertical:MissingField', ...
            'spec.%s is required.', required{ii});
    end
end

f_axis_hz = f_axis_hz(:);
F = numel(f_axis_hz);
ny = spec.ny;
nx = spec.nx;
Nxy = nx * ny;
dxdy = spec.dx_m * spec.dy_m;
fourier_area_scale = spec.dkx_rad_per_m * spec.dky_rad_per_m / (2*pi)^2;
C_eta = real(ifft2(spec.W_eta_kstat)) * Nxy * fourier_area_scale;
sigma2 = max(real(C_eta(1,1)), 0);
alpha_f = 4*pi*f_axis_hz/c0_mps;
R0_f = repmat(complex(reflect_coeff), F, 1);
R_coh_f = R0_f .* exp(-0.5 .* alpha_f.^2 .* sigma2);

C0 = complex(zeros(F, F));
P0 = complex(zeros(F, F));
for ii = 1:F
    for jj = 1:F
        exponent0 = -0.5 * (alpha_f(ii)^2 + alpha_f(jj)^2) * sigma2;
        product_alpha = alpha_f(ii) * alpha_f(jj);
        C0(ii,jj) = R0_f(ii) * conj(R0_f(jj)) * ...
            (exp(exponent0 + product_alpha * sigma2) - exp(exponent0));
        P0(ii,jj) = R0_f(ii) * R0_f(jj) * ...
            (exp(exponent0 - product_alpha * sigma2) - exp(exponent0));
    end
end

rng(mod(round(rng_seed), 2^32), 'twister');
tic_build = tic;
negative_eigenvalue_count = 0;
negative_eigenvalue_abs_sum = 0;
positive_eigenvalue_sum = 0;
max_augmented_size = 0;

if strcmp(mode, 'joint')
    % Store compact single-precision spectra during factorization. Factors
    % are not retained: all requested realizations are drawn in one pass.
    S_cov = complex(zeros(ny, nx, F, F, 'single'));
    S_pseudo = complex(zeros(ny, nx, F, F, 'single'));
    for ii = 1:F
        for jj = ii:F
            exponent0 = -0.5 * (alpha_f(ii)^2 + alpha_f(jj)^2) * sigma2;
            product_alpha = alpha_f(ii) * alpha_f(jj);
            C_xy = R0_f(ii) * conj(R0_f(jj)) .* ...
                (exp(exponent0 + product_alpha .* C_eta) - exp(exponent0));
            P_xy = R0_f(ii) * R0_f(jj) .* ...
                (exp(exponent0 - product_alpha .* C_eta) - exp(exponent0));
            Sc = fft2(C_xy) * dxdy;
            Sp = fft2(P_xy) * dxdy;
            S_cov(:,:,ii,jj) = single(Sc);
            S_cov(:,:,jj,ii) = single(conj(Sc));
            S_pseudo(:,:,ii,jj) = single(Sp);
            S_pseudo(:,:,jj,ii) = single(Sp);
        end
    end
    spectrum_fft_count = F * (F + 1); % covariance + pseudo-covariance.
    spectrum_storage_bytes = numel(S_cov) * 8 + numel(S_pseudo) * 8;
else
    S_diag = zeros(ny, nx, F, 'single');
    for ii = 1:F
        exponent0 = -alpha_f(ii)^2 * sigma2;
        C_xy = abs(R0_f(ii))^2 .* ...
            (exp(exponent0 + alpha_f(ii)^2 .* C_eta) - exp(exponent0));
        S_diag(:,:,ii) = single(real(fft2(C_xy) * dxdy));
    end
    spectrum_fft_count = F;
    spectrum_storage_bytes = numel(S_diag) * 4;
end
build_spectra_time_s = toc(tic_build);

tic_sample = tic;
X_k_fm = complex(zeros(ny, nx, F, n_realizations));
q_scale = Nxy / dxdy;
if strcmp(mode, 'independent')
    for ii = 1:F
        q = q_scale .* max(double(S_diag(:,:,ii)), 0);
        z = (randn(ny,nx,n_realizations) + 1i*randn(ny,nx,n_realizations)) / sqrt(2);
        X_k_fm(:,:,ii,:) = reshape(sqrt(q), ny,nx,1,1) .* reshape(z, ny,nx,1,n_realizations);
    end
else
    visited = false(ny, nx);
    for iy = 1:ny
        jy = mod(-(iy-1), ny) + 1;
        for ix = 1:nx
            if visited(iy,ix)
                continue
            end
            jx = mod(-(ix-1), nx) + 1;
            visited(iy,ix) = true;
            visited(jy,jx) = true;

            Ck = q_scale .* squeeze(double(S_cov(iy,ix,:,:)));
            Cminus = q_scale .* squeeze(double(S_cov(jy,jx,:,:)));
            Pk = q_scale .* squeeze(double(S_pseudo(iy,ix,:,:)));
            Ck = 0.5 * (Ck + Ck');
            Cminus = 0.5 * (Cminus + Cminus');

            if iy == jy && ix == jx
                % Self-conjugate FFT bins need a real 2F-dimensional draw.
                RR = 0.5 * real(Ck + Pk);
                II = 0.5 * real(Ck - Pk);
                RI = 0.5 * (imag(Pk) - imag(Ck));
                augmented = [RR, RI; RI.', II];
                augmented = real(0.5 * (augmented + augmented.'));
                [factor, eig_meta] = local_psd_factor(augmented);
                draw = factor * randn(size(factor,2), n_realizations);
                X_k_fm(iy,ix,:,:) = reshape( ...
                    draw(1:F,:) + 1i*draw(F+1:end,:), 1,1,F,n_realizations);
            else
                % y=[X(k);conj(X(-k))] has covariance [C P;P^H conj(C-)].
                augmented = [Ck, Pk; Pk', conj(Cminus)];
                augmented = 0.5 * (augmented + augmented');
                [factor, eig_meta] = local_psd_factor(augmented);
                z = (randn(size(factor,2), n_realizations) + ...
                    1i*randn(size(factor,2), n_realizations)) / sqrt(2);
                draw = factor * z;
                X_k_fm(iy,ix,:,:) = reshape(draw(1:F,:), 1,1,F,n_realizations);
                X_k_fm(jy,jx,:,:) = reshape(conj(draw(F+1:end,:)), 1,1,F,n_realizations);
            end
            negative_eigenvalue_count = negative_eigenvalue_count + eig_meta.negative_count;
            negative_eigenvalue_abs_sum = negative_eigenvalue_abs_sum + eig_meta.negative_abs_sum;
            positive_eigenvalue_sum = positive_eigenvalue_sum + eig_meta.positive_sum;
            max_augmented_size = max(max_augmented_size, size(augmented,1));
        end
    end
end

deltaG_xy_fm = complex(zeros(ny, nx, F, n_realizations));
for mm = 1:n_realizations
    for ii = 1:F
        deltaG_xy_fm(:,:,ii,mm) = ifft2(X_k_fm(:,:,ii,mm));
    end
end
sample_time_s = toc(tic_sample);

meta = struct();
meta.model = ['kirchhoff_kstat_', mode, '_frequency_phase_screen_v1'];
meta.mode = mode;
meta.f_axis_hz = f_axis_hz;
meta.alpha_f_rad_per_m = alpha_f;
meta.reflect_coeff = reflect_coeff;
meta.sigma_eta2_m2 = sigma2;
meta.R_coh_f = R_coh_f;
meta.C_zero_lag_f = C0;
meta.P_zero_lag_f = P0;
meta.pseudo_to_cov_fro_ratio_zero_lag = norm(P0, 'fro') / max(norm(C0, 'fro'), eps);
meta.n_realizations = n_realizations;
meta.rng_seed = mod(round(rng_seed), 2^32);
meta.spectrum_fft_count = spectrum_fft_count;
meta.realization_ifft_count = F * n_realizations;
meta.spectrum_storage_bytes = spectrum_storage_bytes;
meta.output_storage_bytes = numel(deltaG_xy_fm) * 16;
meta.build_spectra_time_s = build_spectra_time_s;
meta.sample_time_s = sample_time_s;
meta.negative_eigenvalue_count = negative_eigenvalue_count;
meta.negative_eigenvalue_abs_sum = negative_eigenvalue_abs_sum;
meta.positive_eigenvalue_sum = positive_eigenvalue_sum;
meta.negative_to_positive_eigenvalue_ratio = negative_eigenvalue_abs_sum / ...
    max(positive_eigenvalue_sum, eps);
meta.max_augmented_factor_size = max_augmented_size;
meta.covariance_formula = ['C_Gij(rho)=R0_i*conj(R0_j)*exp(-0.5*(alpha_i^2+', ...
    'alpha_j^2)*sigma^2)*(exp(alpha_i*alpha_j*C_eta(rho))-1)'];
meta.pseudocovariance_formula = ['P_Gij(rho)=R0_i*R0_j*exp(-0.5*(alpha_i^2+', ...
    'alpha_j^2)*sigma^2)*(exp(-alpha_i*alpha_j*C_eta(rho))-1)'];
meta.fft_normalization = ['S=fft2(C)*dx*dy; E[X(k)X(k)^H]=', ...
    'numel(grid)/(dx*dy)*S(k); deltaG=ifft2(X).'];
end

function [factor, meta] = local_psd_factor(A)
A = 0.5 * (A + A');
[V, D] = eig(A, 'vector');
D = real(D);
scale = max(max(abs(D)), 1);
tol = 1e-11 * scale;
negative = D < -tol;
positive = D > tol;
meta = struct( ...
    'negative_count', nnz(negative), ...
    'negative_abs_sum', sum(abs(D(negative))), ...
    'positive_sum', sum(D(positive)));
D = max(D, 0);
keep = D > tol;
if any(keep)
    factor = V(:,keep) .* reshape(sqrt(D(keep)).', 1, []);
else
    factor = zeros(size(A,1), 0, 'like', A);
end
end
