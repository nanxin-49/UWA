function [eta, ex, ey] = deepsea3d_surface(x, y, profile)
%DEEPSEA3D_SURFACE Deterministic two-horizontal-coordinate compact surface.
% C2 radial geometric embedding; not an absorbing window or a PM sea.
r=hypot(x,y); a=profile.plateau_m; b=profile.support_m;
assert(b>a && a>0,'deepsea3d:Profile','Invalid embedding radii.');
w=ones(size(r)); wr=zeros(size(r)); w(r>=b)=0;
j=r>a & r<b; t=(r(j)-a)/(b-a);
w(j)=1-10*t.^3+15*t.^4-6*t.^5;
wr(j)=(-30*t.^2+60*t.^3-30*t.^4)/(b-a);
wx=zeros(size(r)); wy=wx; j=r>0;
wx(j)=wr(j).*x(j)./r(j); wy(j)=wr(j).*y(j)./r(j);
v=zeros(size(r)); vx=v; vy=v;
for ii=1:numel(profile.amplitudes_m)
    arg=profile.wavevectors_radpm(ii,1)*x+profile.wavevectors_radpm(ii,2)*y+profile.phases_rad(ii);
    A=profile.amplitudes_m(ii); kv=profile.wavevectors_radpm(ii,:);
    v=v+A*cos(arg); vx=vx-A*kv(1)*sin(arg); vy=vy-A*kv(2)*sin(arg);
end
eta=profile.height_m+w.*v; ex=wx.*v+w.*vx; ey=wy.*v+w.*vy;
end
