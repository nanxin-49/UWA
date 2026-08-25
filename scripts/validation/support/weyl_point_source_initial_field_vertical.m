function [psi_init_xy,meta] = weyl_point_source_initial_field_vertical(X,Y,f_hz,cfg)
%WEYL_POINT_SOURCE_INITIAL_FIELD_VERTICAL Discrete Weyl point-source plane.
% Validation-only source. The continuous transverse transform convention is
% F(kx,ky)=int int f(x,y) exp(-i(kx*x+ky*y)) dx dy. For
% G=exp(i*k*R)/(4*pi*R), F_G=i*exp(i*kz*s0)/(2*kz). The axial carrier
% exp(i*k*s0) is removed to match the PE reduced-field convention.

arguments
    X double
    Y double
    f_hz (1,1) double {mustBeFinite,mustBePositive}
    cfg (1,1) struct
end
if ~isequal(size(X),size(Y)), error('X and Y must have identical sizes.'); end
required={'virtual_source_distance_m','c0','xw','yw','nx','ny','x_tx','y_tx'};
for ii=1:numel(required)
    if ~isfield(cfg,required{ii}), error('Missing cfg.%s.',required{ii}); end
end
s0=cfg.virtual_source_distance_m;
validateattributes(s0,{'numeric'},{'scalar','finite','positive'});
if isfield(cfg,'weyl_regularization_np_per_m')
    eta=cfg.weyl_regularization_np_per_m;
else
    eta=0;
end
validateattributes(eta,{'numeric'},{'scalar','finite','nonnegative'});
if isfield(cfg,'weyl_subcell_order')
    subcell_order=cfg.weyl_subcell_order;
else
    subcell_order=4;
end
validateattributes(subcell_order,{'numeric'},{'scalar','integer','positive','<=',16});

k=2*pi*f_hz/cfg.c0+1i*eta;
k_ref=real(k);
kx=(2*pi/cfg.xw)*[0:(cfg.nx/2-1),-cfg.nx/2:-1];
ky=(2*pi/cfg.yw)*[0:(cfg.ny/2-1),-cfg.ny/2:-1];
[KX0,KY0]=meshgrid(kx,ky);
dkx=2*pi/cfg.xw; dky=2*pi/cfg.yw;
F_continuous=complex(zeros(size(KX0)));
min_abs_kz=Inf; propagating_fraction=0;
sub_offsets=((1:subcell_order)-0.5)/subcell_order-0.5;
for ax=1:subcell_order
    for ay=1:subcell_order
        KX=KX0+sub_offsets(ax)*dkx;
        KY=KY0+sub_offsets(ay)*dky;
        kz=sqrt(complex(k.^2-KX.^2-KY.^2));
        flip_mask=imag(kz)<0 | (abs(imag(kz))<eps(k_ref) & real(kz)<0);
        kz(flip_mask)=-kz(flip_mask);
        source_shift=exp(-1i*(KX*cfg.x_tx+KY*cfg.y_tx));
        F_continuous=F_continuous+1i*exp(1i*(kz-k_ref)*s0).*source_shift./(2*kz);
        min_abs_kz=min(min_abs_kz,min(abs(kz(:))));
        propagating_fraction=propagating_fraction+nnz(real(kz)>0)/numel(kz);
    end
end
F_continuous=F_continuous/(subcell_order^2);
grid_origin_shift=exp(-1i*(KX0*cfg.xw/2+KY0*cfg.yw/2));
dx=cfg.xw/cfg.nx; dy=cfg.yw/cfg.ny;
psi_init_xy=ifft2(F_continuous.*grid_origin_shift/(dx*dy));
meta=struct('reference','weyl_virtual_source_axial_phase_removed', ...
    'virtual_source_distance_m',s0,'regularization_np_per_m',eta, ...
    'continuous_spectrum','i*exp(i*(kz-kref)*s0)/(2*kz)', ...
    'spectral_cell_quadrature','uniform midpoint tensor rule', ...
    'subcell_order',subcell_order,'dft_scale','1/(dx*dy)', ...
    'carrier_to_source_reference',exp(1i*k_ref*s0), ...
    'min_subcell_abs_kz',min_abs_kz, ...
    'propagating_subcell_fraction',propagating_fraction/(subcell_order^2));
end
