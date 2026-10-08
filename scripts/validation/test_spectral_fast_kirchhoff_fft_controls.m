function tests=test_spectral_fast_kirchhoff_fft_controls
% Independent physics/control tests, not production-channel regressions.
tests=functiontests(localfunctions);
end

function testFlatGreenAndSourceDerivatives(testCase)
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root); setup_vertical_project();
n=513; x=linspace(-60,60,n)'; eta=zeros(n,1);
k=2*pi*4000/1500; nx=984; W=192.1875;
xg=(-nx/2:nx/2-1)'*W/nx; kx=2*pi/W*[0:nx/2-1,-nx/2:-1]';
kz=sqrt(complex(k^2-kx.^2,0)); coeff=fft(exp(-.5*(xg/.3).^2));
src=struct('kx',kx,'kz',kz,'coeff',coeff,'period_m',W,'n_original',nx,'x0_m',xg(1),'z_tx_m',100);
% Independent direct sum, including the physical source normal derivative.
phase=exp(1i*(x-xg(1))*kx.').*exp(1i*100*kz.');
u=phase*coeff/nx; dn=phase*(-1i*kz.*coeff)/nx;
normal=[zeros(n,1),ones(n,1)];
c=struct('frequency_hz',4000,'c0_mps',1500,'surface_x_m',x,'surface_z_m',eta, ...
 'surface_normal',normal,'height_bound_m',.01,'receiver_depth_m',3, ...
 'spectral_source',src,'max_rank',6,'height_nodes',33,'rank',6,'geometry','B','padding_factor',2);
plan=prepare_spectral_fast_kirchhoff_2d(c); f=solve_spectral_fast_kirchhoff_2d(c,plan);
q=struct('frequency_hz',4000,'c0_mps',1500,'surface_x_m',x,'surface_z_m',eta, ...
 'surface_normal',normal,'surface_incident',u,'surface_incident_dn',dn, ...
 'receiver_x_m',x,'receiver_z_m',3+zeros(n,1));
ref=solve_full_kirchhoff_2d(q);
verifyLessThan(testCase,norm(f.surface_incident-u)/norm(u),1e-11);
verifyLessThan(testCase,norm(f.surface_incident_dn-dn)/norm(dn),1e-11);
verifyLessThan(testCase,norm(f.receiver_field-ref.receiver_field)/norm(ref.receiver_field),1e-10);
c.padding_factor=4; padded=solve_spectral_fast_kirchhoff_2d(c,plan);
verifyLessThan(testCase,norm(padded.receiver_field-f.receiver_field)/norm(f.receiver_field),1e-12);
c.padding_factor=1; wrapped=solve_spectral_fast_kirchhoff_2d(c,plan);
verifyGreaterThan(testCase,norm(wrapped.receiver_field-f.receiver_field)/norm(f.receiver_field),1e-4);
end

function testHeightInterpolationComplexPolynomial(testCase)
nodes=.2*cos(pi*(0:32)'/32); values=nodes.^7+1i*nodes.^3;
query=linspace(-.2,.2,93)'; exact=query.^7+1i*query.^3;
actual=fast_fk_height_interpolate(nodes,values,query);
verifyLessThan(testCase,max(abs(actual-exact)),1e-14);
at_nodes=fast_fk_height_interpolate(nodes,values,nodes);
verifyEqual(testCase,at_nodes,values);
end
