function out=deepsea3d_shape_spectral_check(cfg,domain_halfwidth,spacing)
%DEEPSEA3D_SHAPE_SPECTRAL_CHECK Independent Fourier half-space propagation.
% This cross-check solves the same FIRST-ORDER boundary data, on a periodic
% computational domain. Domain growth must diagnose image contamination.
clock_run=tic; cpu0=cputime;
n=2*ceil(domain_halfwidth/spacing); h=2*domain_halfwidth/n;
x=-domain_halfwidth+(0:n-1)*h; k=2*pi*cfg.frequency_hz/cfg.c0_mps;
b=complex(zeros(n,n));
for first=1:32:n
    rows=first:min(first+31,n); [xx,yy]=meshgrid(x,x(rows));
    eta=deepsea3d_surface(xx,yy,cfg.profile)-cfg.profile.height_m;
    pts=[xx(:)+cfg.center_m(1),yy(:)+cfg.center_m(2),ones(numel(xx),1)*cfg.profile.height_m];
    [~,d]=deepsea3d_point_field(pts,cfg.tx_m,k);
    b(rows,:)=reshape(-2*eta(:).*d(:,3),numel(rows),n);
end
B=fft2(b); clear b;
kv=2*pi/(2*domain_halfwidth)*[0:n/2-1 -n/2:-1];
offset=cfg.rx_m(1:2)-cfg.center_m;
phase_x=exp(1i*kv*(offset(1)+domain_halfwidth));
value=0;
for first=1:32:n
    rows=first:min(first+31,n); kz=sqrt(complex(k^2-kv.^2-kv(rows).'.^2));
    phase_y=exp(1i*kv(rows).'*(offset(2)+domain_halfwidth));
    kernel=exp(1i*kz*(cfg.rx_m(3)-cfg.profile.height_m)).*phase_y.*phase_x;
    value=value+sum(B(rows,:).*kernel,'all')/n^2;
end
out=struct('rough_correction',value,'nodes',n^2,'spacing_m',h, ...
    'wall_s',toc(clock_run),'cpu_s',cputime-cpu0, ...
    'named_fft_array_MiB',16*n^2/2^20);
end
