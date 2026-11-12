function out = deepsea3d_flat(cfg, height)
%DEEPSEA3D_FLAT Independent infinite pressure-release image solution.
k=2*pi*cfg.frequency_hz/cfg.c0_mps;
tx=cfg.tx_m; rx=cfg.rx_m;
assert(tx(3)>height && rx(3)>height,'deepsea3d:Clearance','Endpoints in water.');
image=tx; image(3)=2*height-tx(3);
Ld=norm(rx-tx); Lr=norm(rx-image);
pd=exp(1i*k*Ld)/(4*pi*Ld); pr=-exp(1i*k*Lr)/(4*pi*Lr);
origin=exp(1i*k*Ld);
out=struct('direct_native',pd,'reflect_native',pr,'total_native',pd+pr, ...
    'direct_dsp',conj(pd)*origin,'reflect_dsp',conj(pr)*origin, ...
    'total_dsp',conj(pd+pr)*origin,'Ld_m',Ld,'Lr_m',Lr, ...
    'excess_delay_s',(Lr-Ld)/cfg.c0_mps);
end
