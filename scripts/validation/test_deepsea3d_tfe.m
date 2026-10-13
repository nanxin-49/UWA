function tests=test_deepsea3d_tfe
tests=functiontests(localfunctions);
end
function testVerticalDerivativeAndRadiation(t)
[s,D,D2]=deepsea3d_tfe_chebyshev(40,.8);
verifyLessThan(t,norm(D*s-1),1e-10);
gamma=[sqrt(25^2-1),1i*12,0,pi/.8];p=deepsea3d_tfe_modal_plan(D,D2,gamma);
verifyGreaterThanOrEqual(t,numel(p.fallback),1);
exact=exp(1i*s*gamma);F=zeros(size(exact));
W=deepsea3d_tfe_modal_solve(p,F,ones(1,4),zeros(1,4));
verifyLessThan(t,norm(W-exact,'fro')/norm(exact,'fro'),1e-9);
verifyLessThan(t,norm(D(end,:)*W-1i*gamma.*W(end,:)),1e-8);
end
function testModalForcedSolveAgainstCoupledSystem(t)
[s,D,D2]=deepsea3d_tfe_chebyshev(24,.8);
gamma=[sqrt(25^2-17),1i*9];p=deepsea3d_tfe_modal_plan(D,D2,gamma);
F=sin(3*s)*[1 2i];g=[.3i .1];h=[.2 .7i];
W=deepsea3d_tfe_modal_solve(p,F,g,h);
for j=1:2
    A=D2+gamma(j)^2*eye(size(D));A(1,:)=0;A(1,1)=1;
    A(end,:)=D(end,:);A(end,end)=A(end,end)-1i*gamma(j);
    rhs=F(:,j);rhs(1)=g(j);rhs(end)=h(j);exact=A\rhs;
    verifyLessThan(t,norm(W(:,j)-exact)/norm(exact),1e-9);
end
end
function testMappedHelmholtzOnExactPlaneWave(t)
n=32;L=8;b=.8;k=15;[s,D,D2]=deepsea3d_tfe_chebyshev(32,b);
kv=pi/L*[0:n/2-1 -n/2:-1];x=-L+(0:n-1)*2*L/n;[xx,yy]=meshgrid(x,x);
p=struct('height_m',0,'plateau_m',100,'support_m',101,'amplitudes_m',.01, ...
    'wavevectors_radpm',[2*pi/L pi/L],'phases_rad',.2);
geo=deepsea3d_tfe_geometry(xx,yy,p);
qx=pi/L;qy=2*pi/L;gamma=sqrt(k^2-qx^2-qy^2);
W0=exp(1i*(qx*xx(:)+qy*yy(:)+gamma*s.'));
A=1i*gamma*geo.eta(:).*(1-s.'/b);W1=A.*W0;W2=A.^2/2.*W0;
[Z0,Z1,Z2]=deepsea3d_tfe_operators(W0,geo,s,D,D2,kv,k,b);
[Y0,Y1,~]=deepsea3d_tfe_operators(W1,geo,s,D,D2,kv,k,b);
[X0,~,~]=deepsea3d_tfe_operators(W2,geo,s,D,D2,kv,k,b);
verifyLessThan(t,norm(Z0,'fro')/norm(W0,'fro')/k^2,1e-8);
verifyLessThan(t,norm(Y0+Z1,'fro')/norm(W1,'fro')/k^2,1e-8);
verifyLessThan(t,norm(X0+Y1+Z2,'fro')/norm(W2,'fro')/k^2,1e-8);
end
