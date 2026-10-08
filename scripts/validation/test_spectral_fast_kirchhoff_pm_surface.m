function tests=test_spectral_fast_kirchhoff_pm_surface
%TEST_SPECTRAL_FAST_KIRCHHOFF_PM_SURFACE Analytic PM embedding invariants.
tests=functiontests(localfunctions);
end

function testAnalyticDerivatives(testCase)
p=struct('k_rad_per_m',[.17;.43;.71],'cos_coeff_m',[.2;-.04;.03], ...
 'sin_coeff_m',[-.1;.08;.02]);
x=[-49.2;-46.1;-42.5;-30;0;35;42.8;47.3;49.6]; h=1e-4;
[~,s,c]=evaluate_pm_pilot_surface(p,x);
[up,sp]=evaluate_pm_pilot_surface(p,x+h);
[um,sm]=evaluate_pm_pilot_surface(p,x-h);
verifyLessThan(testCase,max(abs((up-um)/(2*h)-s)),1e-8);
verifyLessThan(testCase,max(abs((sp-sm)/(2*h)-c)),1e-8);
end

function testCentralPMAndFlatExterior(testCase)
p=struct('k_rad_per_m',[.17;.43],'cos_coeff_m',[.2;-.04],'sin_coeff_m',[-.1;.08]);
x=[-42;-20;0;20;42]; [eta,s,c]=evaluate_pm_pilot_surface(p,x);
phase=x*p.k_rad_per_m.';
verifyEqual(testCase,eta,cos(phase)*p.cos_coeff_m+sin(phase)*p.sin_coeff_m,'AbsTol',1e-14);
verifyEqual(testCase,s,sin(phase)*(-p.k_rad_per_m.*p.cos_coeff_m)+cos(phase)*(p.k_rad_per_m.*p.sin_coeff_m),'AbsTol',1e-14);
verifyEqual(testCase,c,cos(phase)*(-p.k_rad_per_m.^2.*p.cos_coeff_m)+sin(phase)*(-p.k_rad_per_m.^2.*p.sin_coeff_m),'AbsTol',1e-14);
[eta,s,c]=evaluate_pm_pilot_surface(p,[-60;-50;50;60]);
verifyEqual(testCase,[eta;s;c],zeros(12,1),'AbsTol',1e-14);
end

function testUnmodifiedFourierSurface(testCase)
p=struct('k_rad_per_m',[.17;.43],'cos_coeff_m',[.2;-.04],'sin_coeff_m',[-.1;.08]);
x=[-60;-45;0;45;60]; [eta,s,c]=evaluate_pm_pilot_surface(p,x,false);
phase=x*p.k_rad_per_m.';
verifyEqual(testCase,eta,cos(phase)*p.cos_coeff_m+sin(phase)*p.sin_coeff_m,'AbsTol',1e-14);
verifyEqual(testCase,s,sin(phase)*(-p.k_rad_per_m.*p.cos_coeff_m)+cos(phase)*(p.k_rad_per_m.*p.sin_coeff_m),'AbsTol',1e-14);
verifyEqual(testCase,c,cos(phase)*(-p.k_rad_per_m.^2.*p.cos_coeff_m)+sin(phase)*(-p.k_rad_per_m.^2.*p.sin_coeff_m),'AbsTol',1e-14);
end
