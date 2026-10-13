function tests=test_deepsea3d_fg_products
tests=functiontests(localfunctions);
end
function testScalingComplexModes(t)
n=32;x=2*pi*(0:n-1)/n;[xx,yy]=meshgrid(x,x);U=(2+3i)*exp(1i*(7*xx-11*yy));
V=deepsea3d_fg_resize(U,64);back=deepsea3d_fg_resize(V,n);
verifyLessThan(t,norm(U-back,'fro')/norm(U,'fro'),1e-13);
verifyLessThan(t,abs(V(1,1)-(2+3i)),1e-13);
end
function testCubicNeedsDoublePadding(t)
n=32;x=2*pi*(0:n-1)/n;[xx,yy]=meshgrid(x,x);U=cos(12*(xx+yy));expected=.75*U;
small=deepsea3d_fg_resize(deepsea3d_fg_resize(U,48).^3,n);
correct=deepsea3d_fg_resize(deepsea3d_fg_resize(U,64).^3,n);
verifyGreaterThan(t,norm(small-expected,'fro')/norm(expected,'fro'),.3);
verifyLessThan(t,norm(correct-expected,'fro')/norm(expected,'fro'),1e-12);
end
function testNyquistIsExplicit(t)
n=32;[i,j]=meshgrid(0:n-1,0:n-1);U=(-1).^i+(-1).^j;
verifyLessThan(t,norm(deepsea3d_fg_resize(U,64),'fro'),1e-12);
end
function testOperatorsAndTopAgainstAnalyticOversampling(t)
n=32;b=.8;k=50;x=-pi+(0:n-1)*2*pi/n;[xx,yy]=meshgrid(x,x);
[s,D,D2]=deepsea3d_tfe_chebyshev(4,b);kv=[0:n/2-1 -n/2:-1];
profile=struct('height_m',0,'plateau_m',100,'support_m',101, ...
    'amplitudes_m',.01,'wavevectors_radpm',[12 12],'phases_rad',0);
g=deepsea3d_tfe_geometry(xx,yy,profile);gp=struct();
for key=["eta","ex","ey","lap"],gp.(key)=real(deepsea3d_fg_resize(g.(key),2*n));end
gp.slope2=gp.ex.^2+gp.ey.^2;
p=struct('n',n,'s',s,'D',D,'D2',D2,'k',k,'kv',kv,'b',b,'geo',g,'pad_geo',gp, ...
    'gamma',sqrt(complex(k^2-kv.^2-kv.'.^2)));
v=cos(12*(xx+yy));W=v(:)*s.'.^2;[A,B]=deepsea3d_fg_operators(W,p,"galerkin");
x4=-pi+(0:4*n-1)*2*pi/(4*n);[x4,y4]=meshgrid(x4,x4);q=12*(x4+y4);
e=.01*cos(q)/b;ex=-.12*sin(q);lap=-2.88*cos(q);slope=2*ex.^2;
for j=1:numel(s)
    r=1-s(j)/b;H=(k^2-288)*cos(q)*s(j)^2;S=2*s(j)*cos(q);SS=2*cos(q);
    mix=2*ex.*(-24*s(j)*sin(q));
    a=-2*e.*H-2*r*mix-lap*r.*S;
    c=e.^2.*H+slope*r^2.*SS+2*e*r.*mix+(e.*lap-2*slope/b)*r.*S;
    a=deepsea3d_fg_resize(a,n);c=deepsea3d_fg_resize(c,n);
    verifyLessThan(t,norm(A(:,j)-a(:))+norm(B(:,j)-c(:)),1e-9);
end
top=deepsea3d_fg_top_product(v*b^2,p,"galerkin");
expected=-.005i*sqrt(k^2-288)*b*ones(n);
verifyLessThan(t,norm(top-expected,'fro'),1e-10);
end
function testCollocationProjectionMatchesDeclaredSpace(t)
n=32;x=-pi+(0:n-1)*2*pi/n;[xx,yy]=meshgrid(x,x);b=.8;
[s,D,D2]=deepsea3d_tfe_chebyshev(4,b);kv=[0:n/2-1 -n/2:-1];
profile=struct('height_m',0,'plateau_m',100,'support_m',101, ...
    'amplitudes_m',.01,'wavevectors_radpm',[12 0],'phases_rad',0);
g=deepsea3d_tfe_geometry(xx,yy,profile);
p=struct('n',n,'s',s,'D',D,'D2',D2,'k',50,'kv',kv,'b',b,'geo',g);
v=cos(8*xx);[A,B]=deepsea3d_fg_operators(v(:)*s.'.^2,p,"collocation");
for j=1:numel(s)
    F=fft2(reshape(A(:,j)+B(:,j),n,n));
    verifyLessThan(t,norm(F(n/2+1,:))+norm(F(:,n/2+1)),1e-10);
end
end
