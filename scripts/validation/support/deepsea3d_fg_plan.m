function p=deepsea3d_fg_plan(cfg,opts)
%DEEPSEA3D_FG_PLAN Same collocation coefficient interpolants for paired runs.
L=opts.domain_halfwidth_m;n=2*ceil(L/opts.spacing_m);h=2*L/n;
x=-L+(0:n-1)*h;[xx,yy]=meshgrid(x,x);
g=deepsea3d_tfe_geometry(xx,yy,cfg.profile);original=g;
for key=["eta","ex","ey","lap"]
    g.(key)=real(deepsea3d_fg_resize(g.(key),n));
end
g.slope2=g.ex.^2+g.ey.^2;
pad=struct();m=2*n;
for key=["eta","ex","ey","lap"]
    pad.(key)=real(deepsea3d_fg_resize(g.(key),m));
end
pad.slope2=pad.ex.^2+pad.ey.^2;
[s,D,D2]=deepsea3d_tfe_chebyshev(opts.vertical_degree,opts.top_depth_m);
kv=pi/L*[0:n/2-1 -n/2:-1];k=2*pi*cfg.frequency_hz/cfg.c0_mps;
gamma=sqrt(complex(k^2-kv.^2-kv.'.^2));
p=struct('n',n,'h',h,'L',L,'x',x,'xx',xx,'yy',yy,'geo',g,'original_geo',original, ...
    'pad_geo',pad,'s',s,'D',D,'D2',D2,'kv',kv,'k',k,'b',opts.top_depth_m, ...
    'gamma',gamma,'modal',deepsea3d_tfe_modal_plan(D,D2,gamma), ...
    'geometry_Nyquist_relative',norm(g.eta-original.eta,'fro')/norm(original.eta,'fro'));
end
