function tests=test_deepsea3d_finite_reference
tests=functiontests(localfunctions);
end
function testFiniteShiftPlaneWave(t)
% Finite-amplitude exact reflection, not a first-order test or FK comparison.
n=16; k=2*pi*6000/1500; a=.03; z=5; N=14;
eta=ones(n,n)*a; kv=2*pi/16*[0:n/2-1 -n/2:-1];
kz=sqrt(complex(k^2-kv.^2-kv.'.^2));
forcing=cell(N+1,1); forcing{1}=zeros(n,n);
for m=1:N
    forcing{m+1}=ones(n,n)*a^m*((-1i*k)^m-(1i*k)^m)/factorial(m);
end
[~,terms]=deepsea3d_fe_recursion(eta,1i*kz,forcing,exp(1i*kz*z));
exact=(1-exp(-2i*k*a))*exp(1i*k*z);
verifyLessThan(t,abs(sum(terms)-exact)/abs(exact),1e-8);
end
function testFiniteObliqueShift(t)
n=16; k=2*pi*6000/1500; a=.025; z=5; N=14;
kv=2*pi/16*[0:n/2-1 -n/2:-1]; [x,y]=meshgrid(0:n-1,0:n-1);
mode=exp(1i*(kv(2)*x+kv(3)*y)); alpha=sqrt(k^2-kv(2)^2-kv(3)^2);
kz=sqrt(complex(k^2-kv.^2-kv.'.^2));
forcing=cell(N+1,1); forcing{1}=zeros(n,n);
for m=1:N
    forcing{m+1}=mode*a^m*((-1i*alpha)^m-(1i*alpha)^m)/factorial(m);
end
[~,terms]=deepsea3d_fe_recursion(ones(n,n)*a,1i*kz,forcing,exp(1i*kz*z));
exact=(1-exp(-2i*alpha*a))*exp(1i*alpha*z);
verifyLessThan(t,abs(sum(terms)-exact)/abs(exact),1e-8);
end
function testExactPointForcing(t)
cfg=struct('frequency_hz',6000,'c0_mps',1500,'tx_m',[0 0 1500], ...
    'rx_m',[0 0 5],'center_m',[0 0],'profile',struct('height_m',0,'plateau_m',8, ...
    'support_m',12,'amplitudes_m',[.02 .01],'wavevectors_radpm',[.4 0;.2 .6],'phases_rad',[.2 .7]));
opts=struct('domain_halfwidth_m',24,'spacing_m',.25,'order',14,'residual_order',32);
a=deepsea3d_dirichlet_fe(cfg,opts);
verifyLessThan(t,a.forcing_series_relative,1e-8);
verifyLessThan(t,a.boundary_residual_relative,1e-8);
end
