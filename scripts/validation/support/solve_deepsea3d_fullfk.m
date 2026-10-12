function out = solve_deepsea3d_fullfk(cfg)
%SOLVE_DEEPSEA3D_FULLFK Validation-only single-Rx 3D Kirchhoff integral.
% -int (G*d_n uinc + d_n G*uinc)dS, n=(-eta_x,-eta_y,1)/J.
% For compact perturbations of height h0, add exact infinite flat background
% and integrate the DIFFERENCE of full rough and full flat densities.
% This is exact background subtraction inside the Kirchhoff approximation.
arguments
    cfg (1,1) struct
end
assert(cfg.rx_m(3)<cfg.tx_m(3) && cfg.rx_m(3)>cfg.profile.height_m+sum(abs(cfg.profile.amplitudes_m)), ...
    'deepsea3d:Clearance','Invalid conservative receiver clearance.');
[x,w]=deepsea3d_quadrature(cfg.half_width_m,cfg.spacing_m,cfg.quadrature);
k=2*pi*cfg.frequency_hz/cfg.c0_mps; bg=deepsea3d_flat(cfg,cfg.profile.height_m);
total=0; raw=0; rawflat=0; term1=0; term2=0; absolute_sum=0;
clock_run=tic; cpu0=cputime;
for first=1:32:numel(x)
    rows=first:min(first+31,numel(x));
    [xx,yy]=meshgrid(x,x(rows)); weights=w(rows).'*w;
    [z,ex,ey]=deepsea3d_surface(xx,yy,cfg.profile);
    p=[xx(:)+cfg.center_m(1),yy(:)+cfg.center_m(2),z(:)];
    % n*dS=(-eta_x,-eta_y,1)dxdy: exact analytic normal/Jacobian cancellation.
    [ui,du]=deepsea3d_point_field(p,cfg.tx_m,k);
    [G,dG]=deepsea3d_point_field(p,cfg.rx_m,k); % gradient w.r.t. surface point
    normal_area=[-ex(:),-ey(:),ones(numel(z),1)];
    t1=-G.*sum(du.*normal_area,2); t2=-ui.*sum(dG.*normal_area,2);
    p(:,3)=cfg.profile.height_m;
    [uf,df]=deepsea3d_point_field(p,cfg.tx_m,k);
    [gf,dgf]=deepsea3d_point_field(p,cfg.rx_m,k);
    f1=-gf.*df(:,3); f2=-uf.*dgf(:,3);
    ww=weights(:); delta=(t1-f1)+(t2-f2);
    total=total+sum(ww.*delta); raw=raw+sum(ww.*(t1+t2));
    rawflat=rawflat+sum(ww.*(f1+f2));
    term1=term1+sum(ww.*(t1-f1)); term2=term2+sum(ww.*(t2-f2));
    absolute_sum=absolute_sum+sum(ww.*abs(delta));
end
pr=bg.reflect_native+total;
out=struct('reflect_native',pr,'direct_native',bg.direct_native, ...
    'total_native',pr+bg.direct_native,'rough_correction',total, ...
    'term1_correction',term1,'term2_correction',term2, ...
    'raw_truncated_reflect',raw,'raw_truncated_flat',rawflat, ...
    'flat_exact',bg.reflect_native,'quadrature_identity_residual',abs(total-(raw-rawflat)), ...
    'cancellation_condition',absolute_sum/max(abs(total),realmin), ...
    'node_count',numel(x)^2,'wall_s',toc(clock_run),'cpu_s',cputime-cpu0, ...
    'estimated_block_payload_MiB',numel(x)*32*30*16/2^20, ...
    'time_convention','exp(-i omega t)','normal','into water', ...
    'integration','infinite flat + compact full-density difference');
end
