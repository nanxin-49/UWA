function result = run_pe_1d_surface_reflection_validation(cfg)
%RUN_PE_1D_SURFACE_REFLECTION_VALIDATION Validation-only x-z PE bridge.
%   This helper mirrors the production split-step operator in one transverse
%   dimension. It is not used by the public channel API.

if nargin < 1 || ~isstruct(cfg), error('cfg must be a struct.'); end
required = {'frequency_hz','c0_mps','xw_m','nx','z_tx_m','z_rx_m', ...
    'sigma_src_m','surface_elevation_x_m','surface_reflect_coeff','step_m'};
for ii = 1:numel(required)
    if ~isfield(cfg, required{ii}), error('Missing cfg.%s.', required{ii}); end
end
if mod(cfg.nx,2) ~= 0 || cfg.nx < 4, error('nx must be an even integer >= 4.'); end
x = (-0.5*cfg.xw_m) + (0:cfg.nx-1) * (cfg.xw_m/cfg.nx);
eta = double(cfg.surface_elevation_x_m(:).');
if numel(eta) ~= cfg.nx || any(~isfinite(eta))
    error('surface_elevation_x_m must contain nx finite values.');
end
if any(~isfinite([cfg.frequency_hz,cfg.c0_mps,cfg.xw_m,cfg.z_tx_m,cfg.z_rx_m, ...
        cfg.sigma_src_m,cfg.surface_reflect_coeff,cfg.step_m])) || ...
        cfg.frequency_hz <= 0 || cfg.c0_mps <= 0 || cfg.xw_m <= 0 || ...
        cfg.z_tx_m <= cfg.z_rx_m || cfg.sigma_src_m <= 0 || cfg.step_m <= 0
    error('Invalid 1-D PE configuration.');
end

kx = (2*pi/cfg.xw_m) * [0:(cfg.nx/2-1), -cfg.nx/2:-1];
k0 = 2*pi*cfg.frequency_hz/cfg.c0_mps;
kz = sqrt(complex(k0^2-kx.^2,0));
source = exp(-0.5*(x/cfg.sigma_src_m).^2);
psi_surface_inc = local_march(source, cfg.z_tx_m, 0, kz, k0, cfg.step_m);
if ~isfield(cfg,'reflection_model') || isempty(cfg.reflection_model)
    reflection_model = 'model0_normal';
else
    reflection_model = lower(char(cfg.reflection_model));
end
switch reflection_model
    case 'model0_normal'
        phase_screen = cfg.surface_reflect_coeff .* exp(1i*2*k0*eta);
        psi_surface_ref = phase_screen .* psi_surface_inc;
        reflection_formula = 'R0*exp(+i*2*k0*eta(x))*psi_inc(x)';
        reflection_meta = struct('minimum_reflected_kz_radpm',NaN, ...
            'nonreturning_component_count',0);
    case 'model1_kz_aware'
        phase_screen = complex(zeros(0,0));
        psi_surface_ref = local_kz_aware_reflection( ...
            psi_surface_inc,eta,kz,cfg.surface_reflect_coeff);
        reflection_formula = ['R0*sum_k Psi_inc(k)*exp(+i*kx*x)', ...
            '*exp(+i*2*kz(k)*eta(x))/N'];
        reflection_meta = struct('minimum_reflected_kz_radpm',min(real(kz)), ...
            'nonreturning_component_count',sum(real(kz)<=0));
    case 'model2_angle_slope'
        if ~isfield(cfg,'surface_slope_x')
            error('model2_angle_slope requires cfg.surface_slope_x.');
        end
        slope = double(cfg.surface_slope_x(:).');
        if numel(slope)~=cfg.nx || any(~isfinite(slope))
            error('surface_slope_x must contain nx finite values.');
        end
        phase_screen = complex(zeros(0,0));
        [psi_surface_ref,reflection_meta] = local_angle_slope_reflection( ...
            psi_surface_inc,eta,slope,kx,kz,cfg.surface_reflect_coeff);
        reflection_formula = ['R0*sum_k Psi_inc(k)*exp(+i*kx*x)', ...
            '*exp(+i*(kz_i+kz_r)*eta(x))/N; ', ...
            'kz_r=((1-s^2)*kz_i+2*s*kx_i)/(1+s^2)'];
    otherwise
        error('Unknown validation-only reflection_model %s.',reflection_model);
end
psi_direct = local_march(source, cfg.z_tx_m, cfg.z_rx_m, kz, k0, cfg.step_m);
psi_reflected = local_march(psi_surface_ref, 0, cfg.z_rx_m, kz, k0, cfg.step_m);
dx = cfg.xw_m/cfg.nx;
[~,ix_rx] = min(abs(x-cfg.x_rx_m));
result = struct('schema_version','1.0.0','operator','1-D split-step square-root PE', ...
    'frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps,'k0_rad_per_m',k0, ...
    'x_m',x(:).','kx_rad_per_m',kx(:).','dx_m',dx, ...
    'source_field',source(:).','surface_elevation_x_m',eta(:).', ...
    'phase_screen',phase_screen(:).','surface_incident_field',psi_surface_inc(:).', ...
    'surface_reflected_field',psi_surface_ref(:).', ...
    'direct_field',psi_direct(:).','reflected_field',psi_reflected(:).', ...
    'x_rx_m',cfg.x_rx_m,'ix_rx',ix_rx, ...
    'direct_receiver',psi_direct(ix_rx),'reflected_receiver',psi_reflected(ix_rx), ...
    'step_m',cfg.step_m,'direct_range_m',cfg.z_tx_m-cfg.z_rx_m, ...
    'reflected_range_m',cfg.z_tx_m+cfg.z_rx_m, ...
    'reflection_model',reflection_model,'reflection_formula',reflection_formula, ...
    'reflection_meta',reflection_meta, ...
    'phase_reference','raw reduced PE envelope');
end

function [reflected,meta] = local_angle_slope_reflection( ...
        incident,eta,slope,kx,kz,reflect_coeff)
n=numel(incident);incident_k=fft(incident(:)).';indices=0:n-1;
reflected=complex(zeros(n,1));block=128;
minimum_kzr=inf;nonreturning=0;
for first=1:block:n
    rows=first:min(first+block-1,n);
    s=slope(rows).';
    denominator=1+s.^2;
    kzr=((1-s.^2).*kz(:).'+2*s.*kx(:).')./denominator;
    minimum_kzr=min(minimum_kzr,min(real(kzr),[],'all'));
    nonreturning=nonreturning+sum(real(kzr)<=0,'all');
    fourier=exp(1i*2*pi*((rows-1).'*indices)/n)/n;
    height_phase=exp(1i*eta(rows).' .*(kz(:).'+kzr));
    reflected(rows)=reflect_coeff*sum(fourier.*height_phase.*incident_k,2);
end
reflected=reflected.';
meta=struct('minimum_reflected_kz_radpm',minimum_kzr, ...
    'nonreturning_component_count',nonreturning);
end

function reflected = local_kz_aware_reflection(incident,eta,kz,reflect_coeff)
n = numel(incident);
incident_k = fft(incident(:)).';
indices = 0:n-1;
reflected = complex(zeros(n,1));
block = 128;
for first = 1:block:n
    rows = first:min(first+block-1,n);
    fourier = exp(1i*2*pi*((rows-1).'*indices)/n)/n;
    height_phase = exp(1i*2*(eta(rows).'*kz(:).'));
    reflected(rows) = reflect_coeff*sum(fourier.*height_phase.*incident_k,2);
end
reflected = reflected.';
end

function psi_end = local_march(psi_start, z_start, z_end, kz, k0, step_m)
if abs(z_end-z_start) <= eps
    psi_end = psi_start;
    return
end
n_step = max(1,ceil(abs(z_end-z_start)/step_m));
ds = abs((z_end-z_start)/n_step);
fr = exp(-1i*0.5*ds*(k0-kz));
psi_k = fft(psi_start);
for jj = 1:n_step
    psi_k = fr .* fft(ifft(fr .* psi_k));
end
psi_end = ifft(psi_k);
end
