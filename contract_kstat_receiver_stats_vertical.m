function [stats, meta] = contract_kstat_receiver_stats_vertical( ...
    pm_spec, projection, options)
%CONTRACT_KSTAT_RECEIVER_STATS_VERTICAL Propagate joint K-stat C/P to one Rx.
%   Validation-only dense or PM-periodic FFT contraction. The current
%   deltaG convention already includes R0; projection weights must not.

arguments
    pm_spec (1,1) struct
    projection (1,1) struct
    options (1,1) struct = struct()
end
options = local_options(options);
local_validate_inputs(pm_spec, projection, options);

F = numel(projection.f_axis_hz);
ny = pm_spec.ny;
nx = pm_spec.nx;
Nxy = nx*ny;
mapping = projection.pm_mapping;
a_pm_xy_f = complex(zeros(ny,nx,F));
a_pm_xy_f(mapping.iy,mapping.ix,:) = projection.a_pe_xy_f;

inside_error = max(abs(reshape( ...
    a_pm_xy_f(mapping.iy,mapping.ix,:)-projection.a_pe_xy_f,[],1)));
outside_mask = true(ny,nx);
outside_mask(mapping.iy,mapping.ix) = false;
outside_values = a_pm_xy_f(repmat(outside_mask,1,1,F));
if isempty(outside_values)
    outside_max_abs = 0;
else
    outside_max_abs = max(abs(outside_values));
end
if inside_error ~= 0 || outside_max_abs ~= 0
    error('contract_kstat_receiver_stats_vertical:EmbeddingFailure', ...
        'PM zero embedding does not exactly match the cached central crop.');
end

fourier_area_scale = pm_spec.dkx_rad_per_m*pm_spec.dky_rad_per_m/(2*pi)^2;
C_eta_xy = real(ifft2(pm_spec.W_eta_kstat))*Nxy*fourier_area_scale;
sigma_eta2 = max(real(C_eta_xy(1,1)),0);
alpha_f = 4*pi*projection.f_axis_hz(:)/projection.c0_mps;
R0_f = repmat(complex(projection.reflect_coeff),F,1);
R_coh_f = R0_f.*exp(-0.5*alpha_f.^2*sigma_eta2);

C_H = complex(zeros(F,F));
P_H = complex(zeros(F,F));
pair_count = F*(F+1)/2;
pair_time_s = zeros(pair_count,1);
pair_index = 0;
total_timer = tic;

if strcmp(options.method,'fft')
    A_pm_k_f = fft2(a_pm_xy_f);
    B_pm_k_f = fft2(conj(a_pm_xy_f));
    lag_index = [];
else
    A_pm_k_f = complex(zeros(0,0,0));
    B_pm_k_f = complex(zeros(0,0,0));
    lag_index = local_periodic_lag_index(ny,nx);
end

for ii = 1:F
    for jj = ii:F
        pair_timer = tic;
        [C_lag_xy,P_lag_xy] = local_joint_lag_pair( ...
            C_eta_xy,sigma_eta2,alpha_f(ii),alpha_f(jj),R0_f(ii),R0_f(jj));
        if strcmp(options.method,'fft')
            C_hat = fft2(C_lag_xy);
            P_hat = fft2(P_lag_xy);
            Ai = A_pm_k_f(:,:,ii);
            Aj = A_pm_k_f(:,:,jj);
            Bj = B_pm_k_f(:,:,jj);
            term_c = conj(Ai).*C_hat.*Aj;
            term_p = conj(Ai).*P_hat.*Bj;
            C_value = sum(term_c(:))/Nxy;
            P_value = sum(term_p(:))/Nxy;
        else
            C_matrix = reshape(C_lag_xy(lag_index),Nxy,Nxy);
            P_matrix = reshape(P_lag_xy(lag_index),Nxy,Nxy);
            ai = reshape(a_pm_xy_f(:,:,ii),Nxy,1);
            aj = reshape(a_pm_xy_f(:,:,jj),Nxy,1);
            C_value = ai'*(C_matrix*aj);
            P_value = ai'*(P_matrix*conj(aj));
        end
        C_H(ii,jj) = C_value;
        C_H(jj,ii) = conj(C_value);
        P_H(ii,jj) = P_value;
        P_H(jj,ii) = P_value;
        pair_index = pair_index+1;
        pair_time_s(pair_index) = toc(pair_timer);
    end
end
total_s = toc(total_timer);

hermitian_error_before = norm(C_H-C_H','fro')/max(norm(C_H,'fro'),eps);
pseudo_symmetry_error_before = norm(P_H-P_H.','fro')/max(norm(P_H,'fro'),eps);
C_H = 0.5*(C_H+C_H');
P_H = 0.5*(P_H+P_H.');
augmented = [C_H,P_H;conj(P_H),conj(C_H)];
augmented = 0.5*(augmented+augmented');
augmented_eigenvalues = real(eig(augmented));

mu_scatter_f = complex(zeros(F,1));
mu_total_f = projection.H_direct_f(:)+projection.H_ref_coh_f(:)+mu_scatter_f;
stats = struct();
stats.method = options.method;
stats.f_axis_hz = projection.f_axis_hz(:);
stats.mu_scatter_f = mu_scatter_f;
stats.mu_total_f = mu_total_f;
stats.C_H = C_H;
stats.P_H = P_H;
stats.E_abs_H2_f = abs(mu_scatter_f).^2+real(diag(C_H));
stats.R_coh_f = R_coh_f;
stats.sigma_eta2_m2 = sigma_eta2;
stats.augmented_covariance_eigenvalues = augmented_eigenvalues;
stats.augmented_covariance_min_eigenvalue = min(augmented_eigenvalues);
stats.deltaG_definition = projection.deltaG_definition;
stats.weight_definition = projection.weight_definition;

dense_block_bytes = 2*Nxy^2*16;
full_dense_cp_bytes = 2*F^2*Nxy^2*16;
meta = struct();
meta.method = options.method;
meta.total_s = total_s;
meta.pair_count = pair_count;
meta.mean_pair_s = mean(pair_time_s);
meta.max_pair_s = max(pair_time_s);
meta.pair_time_s = pair_time_s;
meta.embedding_inside_max_abs_error = inside_error;
meta.embedding_outside_max_abs = outside_max_abs;
meta.hermitian_relative_error_before_symmetrize = hermitian_error_before;
meta.pseudo_symmetry_relative_error_before_symmetrize = pseudo_symmetry_error_before;
meta.a_pm_spatial_bytes = numel(a_pm_xy_f)*16;
meta.a_pm_fft_bytes = numel(A_pm_k_f)*16+numel(B_pm_k_f)*16;
meta.streamed_dense_cp_block_bytes = dense_block_bytes;
meta.forbidden_full_dense_cp_bytes = full_dense_cp_bytes;
meta.streamed_lag_pair_bytes = 2*Nxy*16;
meta.fft_scaling = 'MATLAB unnormalized fft2; Parseval contraction uses 1/(nx*ny)';
meta.lag_origin = '[1,1], unshifted FFT order';
meta.pseudo_k_coupling = 'B_j=fft2(conj(a_j)); equivalent K/-K coupling is retained explicitly';
meta.memory_snapshot_bytes = local_memory_used_bytes();
meta.options = options;
end

function options = local_options(options)
defaults = struct('method','fft','max_dense_points',32^2);
names = fieldnames(defaults);
for ii = 1:numel(names)
    if ~isfield(options,names{ii})
        options.(names{ii}) = defaults.(names{ii});
    end
end
if isstring(options.method)
    options.method = char(options.method);
end
options.method = lower(strtrim(options.method));
if ~ismember(options.method,{'dense','fft'})
    error('contract_kstat_receiver_stats_vertical:Method', ...
        'options.method must be dense or fft.');
end
validateattributes(options.max_dense_points,{'numeric'}, ...
    {'scalar','integer','positive','finite'});
end

function local_validate_inputs(pm_spec, projection, options)
pm_required = {'nx','ny','dx_m','dy_m','dkx_rad_per_m','dky_rad_per_m','W_eta_kstat'};
for ii = 1:numel(pm_required)
    if ~isfield(pm_spec,pm_required{ii})
        error('contract_kstat_receiver_stats_vertical:InvalidPMSpec', ...
            'pm_spec.%s is required.',pm_required{ii});
    end
end
projection_required = {'f_axis_hz','a_pe_xy_f','H_direct_f','H_ref_coh_f', ...
    'pm_mapping','c0_mps','reflect_coeff','deltaG_definition','weight_definition'};
for ii = 1:numel(projection_required)
    if ~isfield(projection,projection_required{ii})
        error('contract_kstat_receiver_stats_vertical:InvalidProjection', ...
            'projection.%s is required.',projection_required{ii});
    end
end
if ~isa(projection.a_pe_xy_f,'double') || isa(projection.a_pe_xy_f,'gpuArray')
    error('contract_kstat_receiver_stats_vertical:CPUdoubleRequired', ...
        'projection weights must be CPU double.');
end
if ~isequal(projection.pm_mapping.pm_grid,[pm_spec.ny,pm_spec.nx])
    error('contract_kstat_receiver_stats_vertical:MappingMismatch', ...
        'Projection mapping and PM grid dimensions differ.');
end
if size(projection.a_pe_xy_f,3) ~= numel(projection.f_axis_hz)
    error('contract_kstat_receiver_stats_vertical:FrequencyMismatch', ...
        'Projection weight frequency dimension is inconsistent.');
end
if strcmp(options.method,'dense') && pm_spec.nx*pm_spec.ny > options.max_dense_points
    error('contract_kstat_receiver_stats_vertical:DenseGridTooLarge', ...
        'Dense contraction is limited to options.max_dense_points.');
end
end

function [C_lag_xy,P_lag_xy] = local_joint_lag_pair( ...
    C_eta_xy,sigma_eta2,alpha_i,alpha_j,R0_i,R0_j)
exponent0 = -0.5*(alpha_i^2+alpha_j^2)*sigma_eta2;
product_alpha = alpha_i*alpha_j;
C_lag_xy = R0_i*conj(R0_j).*( ...
    exp(exponent0+product_alpha.*C_eta_xy)-exp(exponent0));
P_lag_xy = R0_i*R0_j.*( ...
    exp(exponent0-product_alpha.*C_eta_xy)-exp(exponent0));
end

function lag_index = local_periodic_lag_index(ny,nx)
Nxy = nx*ny;
[iy,ix] = ind2sub([ny,nx],1:Nxy);
lag_y = mod(iy(:)-iy(:).',ny)+1;
lag_x = mod(ix(:)-ix(:).',nx)+1;
lag_index = uint32(sub2ind([ny,nx],lag_y,lag_x));
end

function value = local_memory_used_bytes()
value = NaN;
try
    info = memory;
    value = info.MemUsedMATLAB;
catch
end
end
