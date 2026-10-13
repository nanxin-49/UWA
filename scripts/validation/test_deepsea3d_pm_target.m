function tests=test_deepsea3d_pm_target
tests=functiontests(localfunctions);
end
function testSpectrumAndParseval(t)
s=build_deepsea3d_pm_target();d=s.parent_diagnostics;
verifyLessThan(t,abs(d.ensemble_Hm0_m-.5),1e-12);
verifyLessThan(t,d.parseval_relative,1e-12);
alpha=s.alpha;kp=s.kp_radpm;
m0=integral(@(k)alpha./(2*k.^3).*exp(-1.25*(kp./k).^2),1e-5,Inf);
verifyLessThan(t,abs(m0-alpha/(5*kp^2))/m0,1e-8);
verifyLessThan(t,d.height_global_bound_m,5);
end
function testComplexToRealSurfaceAndGradient(t)
s=build_deepsea3d_pm_target();m=s.baseline_mask;
c=complex(s.coefficient_real(m),s.coefficient_imag(m));
xy=[0 0;1.1 -2.3;-3.2 .7];
phase=exp(1i*(s.kx_radpm(m)*xy(:,1).'+s.ky_radpm(m)*xy(:,2).'));
exact=2*real(sum(c.*phase,1)).';
dx=2*real(sum(1i*s.kx_radpm(m).*c.*phase,1)).';
dy=2*real(sum(1i*s.ky_radpm(m).*c.*phase,1)).';
[eta,ex,ey]=deepsea3d_surface(xy(:,1),xy(:,2),s.profile);
verifyLessThan(t,norm(eta-exact),1e-12);
verifyLessThan(t,norm(ex-dx)+norm(ey-dy),1e-12);
end
function testDeterministicAndSharedCoefficients(t)
s=build_deepsea3d_pm_target();b=build_deepsea3d_pm_target();
verifyEqual(t,s.coefficient_real,b.coefficient_real);
verifyEqual(t,s.coefficient_imag,b.coefficient_imag);
verifyTrue(t,all(~s.baseline_mask | s.short_extended_mask));
verifyTrue(t,all(~s.short_removed_mask | s.baseline_mask));
verifyGreaterThan(t,s.short_extended_diagnostics.ensemble_Hm0_m,.5);
verifyLessThan(t,s.short_removed_diagnostics.ensemble_Hm0_m,.5);
verifyGreaterThan(t,abs(s.parent_diagnostics.realized_Hm0_m-.5),1e-4);
end
