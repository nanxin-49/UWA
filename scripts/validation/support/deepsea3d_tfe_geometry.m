function geo=deepsea3d_tfe_geometry(x,y,p)
%DEEPSEA3D_TFE_GEOMETRY Analytic C2 surface including horizontal Laplacian.
[eta,ex,ey]=deepsea3d_surface(x,y,p);
r=hypot(x,y);a=p.plateau_m;b=p.support_m;
w=ones(size(r));wr=zeros(size(r));wrr=wr;w(r>=b)=0;
j=r>a & r<b;t=(r(j)-a)/(b-a);
w(j)=1-10*t.^3+15*t.^4-6*t.^5;
wr(j)=(-30*t.^2+60*t.^3-30*t.^4)/(b-a);
wrr(j)=(-60*t+180*t.^2-120*t.^3)/(b-a)^2;
v=zeros(size(r));vx=v;vy=v;lapv=v;
for j=1:numel(p.amplitudes_m)
    q=p.wavevectors_radpm(j,:);arg=q(1)*x+q(2)*y+p.phases_rad(j);
    A=p.amplitudes_m(j);v=v+A*cos(arg);vx=vx-A*q(1)*sin(arg);vy=vy-A*q(2)*sin(arg);
    lapv=lapv-A*sum(q.^2)*cos(arg);
end
wx=zeros(size(r));wy=wx;lapw=wrr;j=r>0;
wx(j)=wr(j).*x(j)./r(j);wy(j)=wr(j).*y(j)./r(j);lapw(j)=lapw(j)+wr(j)./r(j);
geo=struct('eta',eta,'ex',ex,'ey',ey,'lap',w.*lapv+2*wx.*vx+2*wy.*vy+lapw.*v, ...
    'slope2',ex.^2+ey.^2);
end
