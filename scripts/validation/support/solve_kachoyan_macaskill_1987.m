function out = solve_kachoyan_macaskill_1987(cfg)
%SOLVE_KACHOYAN_MACASKILL_1987 Validation-only implementation of the
% dimensionless collocation and reflection-coefficient method of
% Kachoyan and Macaskill, JASA 82, 1720-1726 (1987).
%
% The surface is z = hhat*xi(x), x is dimensionless, and khat=k*L.
% cfg.x and cfg.xi contain N uniformly spaced samples. The optional
% periodic_replicas field evaluates the kernel over repeated periods for
% the sinusoidal Table-II cases. No production channel code calls this file.

arguments
    cfg (1,1) struct
end
required = {'x','xi','xi_prime','xi_second','dx','khat','hhat', ...
    'theta_deg','periodic_replicas','mode','z_samples'};
for ii = 1:numel(required)
    if ~isfield(cfg,required{ii})
        error('solve_kachoyan_macaskill_1987:MissingField', ...
            'cfg.%s is required.',required{ii});
    end
end
x = cfg.x(:); xi = cfg.xi(:); xi_prime = cfg.xi_prime(:); xi_second = cfg.xi_second(:);
N = numel(x);
if any([numel(xi),numel(xi_prime),numel(xi_second)] ~= N)
    error('Surface arrays must all have length N.');
end
if N < 8 || mod(N,2) ~= 0 || any(~isfinite([x;xi;xi_prime;xi_second]))
    error('x and surface arrays must be finite with even N >= 8.');
end
if cfg.dx <= 0 || cfg.khat <= 0 || cfg.hhat < 0
    error('dx and khat must be positive and hhat must be nonnegative.');
end
if abs(median(diff(x))-cfg.dx) > 1e-9*max(1,cfg.dx) || ...
        max(abs(diff(x)-cfg.dx)) > 1e-7*max(1,cfg.dx)
    error('x must be uniformly spaced with cfg.dx.');
end
alpha0 = cosd(cfg.theta_deg);
gamma0 = sind(cfg.theta_deg);
if gamma0 <= 0
    error('theta_deg must lie in (0,180) degrees from the horizontal.');
end
T = N*cfg.dx;
if cfg.periodic_replicas < 0 || cfg.periodic_replicas ~= floor(cfg.periodic_replicas)
    error('periodic_replicas must be a nonnegative integer.');
end
if isfield(cfg,'window') && ~isempty(cfg.window)
    window = cfg.window(:);
    if numel(window) ~= N, error('cfg.window must have N elements.'); end
else
    window = ones(N,1);
end
if any(~isfinite(window)) || any(window < 0)
    error('cfg.window must be finite and nonnegative.');
end

% Eq. (5) kernel and its local diagonal limit. The coefficient below is
% the dimensionless Hankel H_1^(1) form; the limit is the curvature term
% stated below Eq. (5) in the paper and avoids evaluating H_1 at zero.
Ksum = complex(zeros(N,N));
rep_range = -cfg.periodic_replicas:cfg.periodic_replicas;
for qq = rep_range
    source_x = x.' + qq*T;
    dxmat = source_x - x;
    dximat = xi.' - xi;
    rho2 = dxmat.^2 + cfg.hhat^2*dximat.^2;
    rho = sqrt(rho2);
    num = xi.' - xi - dxmat.*xi_prime;
    Kq = complex(zeros(N,N));
    nz = rho > 0;
    Kq(nz) = -0.5i*pi*cfg.khat*(cfg.hhat^2) .* ...
        (num(nz)./rho(nz)) .* besselh(1,1,cfg.khat*rho(nz));
    if qq == 0
        diag_lim = 0.5*(cfg.hhat^2)*xi_second ./ ...
            (1 + cfg.hhat^2*xi_prime.^2);
        Kq(1:N+1:end) = diag_lim;
    end
    Ksum = Ksum + Kq .* window.';
end

psi_inc = -(gamma0 + alpha0*cfg.hhat*xi_prime) .* ...
    exp(1i*cfg.khat*(alpha0*x - gamma0*cfg.hhat*xi));
A = eye(N) + cfg.dx*Ksum;
rhs = 2*psi_inc;
psi_full = A\rhs;
psi_ka = rhs;
residual_full = norm(A*psi_full-rhs,2)/max(norm(rhs,2),realmin);
residual_ka = norm(A*psi_ka-rhs,2)/max(norm(rhs,2),realmin);

% Eq. (8): the normalized transform uses the + sign, hence MATLAB ifft.
psi_bar_full = exp(-1i*cfg.khat*alpha0*x).*psi_full;
psi_bar_ka = exp(-1i*cfg.khat*alpha0*x).*psi_ka;
psi_hat_full = ifft(psi_bar_full);
psi_hat_ka = ifft(psi_bar_ka);
% Independent transformed-system solve in the DFT variables.
F = ifft(eye(N)); Finv = fft(eye(N)); D = diag(exp(-1i*cfg.khat*alpha0*x));
B = F*D*A*(D\Finv); rhs_hat_direct = F*D*rhs;
psi_hat_direct = B\rhs_hat_direct;
dft_system_error = norm(psi_hat_direct-psi_hat_full,2)/max(norm(psi_hat_full,2),realmin);
mode_index = [0:(N/2-1), -N/2:-1].';
q_mode = alpha0 + 2*pi*mode_index/(cfg.khat*T);
gamma_mode = complex(nan(N,1));
propagating = abs(q_mode) < 1;
gamma_mode(propagating) = sqrt(max(0,1-q_mode(propagating).^2));
prop_indices = find(propagating);
R_full = local_reflection_coefficients(psi_hat_full,cfg.hhat*xi,cfg.khat,gamma_mode,propagating);
R_ka = local_reflection_coefficients(psi_hat_ka,cfg.hhat*xi,cfg.khat,gamma_mode,propagating);
E_full = real(sum(abs(R_full(propagating)).^2 .* gamma_mode(propagating))/gamma0);
E_ka = real(sum(abs(R_ka(propagating)).^2 .* gamma_mode(propagating))/gamma0);

z = cfg.z_samples(:).';
if any(~isfinite(z)) || any(z <= max(cfg.hhat*xi))
    error('All z_samples must satisfy z > hhat*max(xi), as required by Appendix A.');
end
x_eval = x.';
ps_full = local_reconstruct(R_full,propagating,q_mode,gamma_mode, ...
    cfg.khat,x_eval,z);
ps_ka = local_reconstruct(R_ka,propagating,q_mode,gamma_mode, ...
    cfg.khat,x_eval,z);
out = struct();
out.schema_version = '1.0.0';
out.paper = 'Kachoyan and Macaskill JASA 82(5), 1720-1726 (1987)';
out.formulation = 'dimensionless Fredholm second-kind collocation plus Eq. (15)-(17) reflection modes';
out.coordinate_convention = 'x horizontal; z positive into water; theta from horizontal';
out.time_convention = 'exp(-i omega t)';
out.x = x; out.xi = xi; out.xi_prime = xi_prime; out.xi_second = xi_second;
out.dx = cfg.dx; out.period_length = T; out.khat = cfg.khat; out.hhat = cfg.hhat;
out.theta_deg = cfg.theta_deg; out.alpha0 = alpha0; out.gamma0 = gamma0;
out.periodic_replicas = cfg.periodic_replicas; out.window = window;
out.kernel_matrix = Ksum; out.system_matrix = A; out.rhs = rhs;
out.psi_inc = psi_inc; out.psi_full = psi_full; out.psi_ka = psi_ka;
out.psi_hat_full = psi_hat_full; out.psi_hat_ka = psi_hat_ka;
out.linear_residual_full = residual_full; out.linear_residual_ka = residual_ka; out.dft_system_error = dft_system_error;
out.mode_index = mode_index; out.q_mode = q_mode; out.gamma_mode = gamma_mode;
out.propagating_mask = propagating; out.propagating_indices = prop_indices;
out.propagating_count = numel(prop_indices);
out.R_full = R_full; out.R_ka = R_ka; out.energy_full = E_full; out.energy_ka = E_ka;
out.z_samples = z; out.ps_full = ps_full; out.ps_ka = ps_ka;
out.nearfield_intensity_full = abs(ps_full).^2;
out.nearfield_intensity_ka = abs(ps_ka).^2;
out.finite = all(isfinite([real(psi_full);imag(psi_full);real(R_full(propagating)); ...
    imag(R_full(propagating));real(ps_full(:));imag(ps_full(:))]));
out.mode_formula = 'q_r=alpha0+2*pi*r/(khat*N*dx), gamma_r=sqrt(1-q_r^2)';
out.energy_formula = 'E=Re(sum(|R_r|^2*gamma_r/gamma0)) over propagating modes';
out.diagonal_limit_formula = '0.5*hhat^2*xi''''/(1+hhat^2*xi''^2)';
end

function R = local_reflection_coefficients(psi_hat,xi,khat,gamma_mode,mask)
N = numel(psi_hat);
R = complex(nan(N,1));
for rr = find(mask(:)).'
    gamma = gamma_mode(rr);
    C = fft(exp(-1i*khat*gamma*xi))/N;
    % C_{r-m} with indices interpreted modulo N.
    idx = mod((0:N-1) - (rr-1),N) + 1;
    R(rr) = sum(psi_hat(:).*C(idx))/(2*gamma);
end
end

function ps = local_reconstruct(R,mask,q_mode,gamma_mode,khat,x,z)
ps = complex(zeros(numel(z),numel(x)));
for rr = find(mask(:)).'
    phase_x = exp(1i*khat*q_mode(rr)*x);
    phase_z = exp(1i*khat*gamma_mode(rr)*z(:));
    ps = ps + (phase_z*phase_x)*R(rr);
end
end
