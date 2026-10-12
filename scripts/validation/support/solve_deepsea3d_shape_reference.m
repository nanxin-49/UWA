function out = solve_deepsea3d_shape_reference(cfg)
%SOLVE_DEEPSEA3D_SHAPE_REFERENCE Independent first Dirichlet shape derivative.
% Different formulation/quadrature from FK, flat background h0.
% Boundary data b=-eta*partial_z(inc+ref_flat)=-2*eta*partial_z inc.
% Half-space Rayleigh kernel K=-2*partial_R G*(z_rx-h0)/R.
% delta_u=int K*b dxdy. Evanescent components retained by real-space kernel.
[x,w]=deepsea3d_quadrature(cfg.half_width_m,cfg.spacing_m,'gauss4');
k=2*pi*cfg.frequency_hz/cfg.c0_mps; h0=cfg.profile.height_m;
clock_run=tic; cpu0=cputime; delta=0;
for first=1:32:numel(x)
    rows=first:min(first+31,numel(x)); [xx,yy]=meshgrid(x,x(rows));
    weights=w(rows).'*w;
    eta=deepsea3d_surface(xx,yy,cfg.profile)-h0;
    pts=[xx(:)+cfg.center_m(1),yy(:)+cfg.center_m(2),ones(numel(xx),1)*h0];
    [~,di]=deepsea3d_point_field(pts,cfg.tx_m,k);
    R=sqrt(sum((pts-cfg.rx_m).^2,2));
    rayleigh=-2*exp(1i*k*R)./(4*pi*R).*(1i*k-1./R)*(cfg.rx_m(3)-h0)./R;
    delta=delta+sum(weights(:).*rayleigh.*(-2*eta(:).*di(:,3)));
end
bg=deepsea3d_flat(cfg,h0);
out=struct('rough_correction',delta,'reflect_native',bg.reflect_native+delta, ...
    'node_count',numel(x)^2,'wall_s',toc(clock_run),'cpu_s',cputime-cpu0, ...
    'scope','first-order shape derivative only, not finite-roughness BIE');
end
