function tests=test_pm_validity_metrics
%TEST_PM_VALIDITY_METRICS Exact error decomposition and source closure.
tests=functiontests(localfunctions);
end

function testPurePhaseIsNotAligned(testCase)
b=[1;2i;-.5]; phi=.003; a=b*exp(1i*phi); w=[.2;.3;.5];
m=pm_validity_field_metrics(a,b,w,true(3,1));
verifyEqual(testCase,m.complex_l2,2*sin(phi/2),'AbsTol',1e-14);
verifyEqual(testCase,m.phase_mean_rad,phi,'AbsTol',1e-14);
verifyEqual(testCase,m.phase_rms_rad,phi,'AbsTol',1e-14);
verifyLessThan(testCase,m.magnitude_l2,1e-14);
verifyLessThan(testCase,m.decomposition_closure,1e-14);
end

function testMixedMagnitudePhaseIdentity(testCase)
b=[1;2i;-.5;-.2i]; a=b.*[1.01;.98;1.03;.99].*exp(1i*[.01;-.02;.004;.008]);
m=pm_validity_field_metrics(a,b,[.2;.3;.1;.4],true(4,1));
verifyEqual(testCase,m.complex_l2^2,m.magnitude_l2^2+m.phase_contribution_l2^2,'AbsTol',1e-14);
end

function testFrozenSourceAndPhysicalDerivatives(testCase)
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
q=load(fullfile(root,'results','validation','spectral_fast_kirchhoff_pm_pilot','spectral_fast_fk_pm_pilot.mat'),'result');
s=q.result.source; d=q.result.cases(1); x=d.x_m; eta=d.eta_m;
[u,~,~]=evaluate_pm_validity_incident(s,x,eta);
verifyLessThan(testCase,norm(u-d.refcfg.surface_incident)/norm(u),1e-13);
x=[-8;0;7]; z=[.1;-.2;.3]; h=1e-5;
[~,ux,uz]=evaluate_pm_validity_incident(s,x,z);
up=evaluate_pm_validity_incident(s,x+h,z); um=evaluate_pm_validity_incident(s,x-h,z);
verifyLessThan(testCase,norm((up-um)/(2*h)-ux)/norm(ux),1e-7);
up=evaluate_pm_validity_incident(s,x,z+h); um=evaluate_pm_validity_incident(s,x,z-h);
verifyLessThan(testCase,norm((up-um)/(2*h)-uz)/norm(uz),1e-7);
end
