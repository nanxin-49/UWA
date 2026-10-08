function tests=test_spectral_fast_kirchhoff_multifrequency_phase
%TEST_SPECTRAL_FAST_KIRCHHOFF_MULTIFREQUENCY_PHASE Independent phase controls.
tests=functiontests(localfunctions);
end

function testKnownDelayAndPhysicalSign(t)
f=(4000:62.5:8000)'; g=struct('z_tx',100,'z_rx',3,'c0',1500,'z_surface',0);
r=struct('direct_f',ones(size(f)),'reflect_fm',-ones(size(f)));
[a,m]=apply_pe_channel_phase_reference_vertical(f,r,g,'direct_dsp');
[b,~]=apply_pe_channel_phase_reference_vertical(f,r,g,'absolute_physical');
verifyEqual(t,m.relative_delay_s,.004,'AbsTol',1e-16);
phi=unwrap(angle(a.reflect_fm)); tau=-gradient(phi,f)/(2*pi);
verifyEqual(t,tau,.004+zeros(size(f)),'AbsTol',1e-14);
% Equivalent multiplication orders round the ~3500-radian argument by ulps.
phase_roundoff=8*eps(max(2*pi*f*103/1500));
verifyEqual(t,b.reflect_fm,-exp(1i*2*pi*f*103/1500),'AbsTol',phase_roundoff);
verifyEqual(t,a.total_fm,a.direct_f+a.reflect_fm,'AbsTol',1e-14);
h=ifft(a.reflect_fm); [~,ii]=max(abs(h));
verifyLessThan(t,abs((ii-1)/(65*62.5)-.004),1/(65*62.5));
end

function testSparseAliasCannotCertifyDelay(t)
f=(4000:1000:8000)'; h=-exp(-1i*2*pi*f*.004);
tau=-gradient(unwrap(angle(h)),f)/(2*pi);
verifyLessThan(t,max(abs(tau)),1e-14);
verifyEqual(t,1/mean(diff(f)),.001,'AbsTol',1e-16);
end

function testTransversePhaseRequiresWholeEnvelopeConversion(t)
f=(5990:.1:6010)'; k=2*pi*f/1500; kz=sqrt(k.^2-4^2);
p=-exp(1i*kz*103); reduced=p.*exp(-1i*k*103);
g=struct('z_tx',100,'z_rx',3,'c0',1500);
a=apply_pe_channel_phase_reference_vertical(f, ...
 struct('direct_f',ones(size(f)),'reflect_fm',reduced),g,'direct_dsp');
native=conj(p).*exp(1i*k*97);
adapted=apply_pe_channel_phase_reference_vertical(f, ...
 struct('direct_f',ones(size(f)),'reflect_fm',conj(reduced)),g,'direct_dsp');
verifyEqual(t,adapted.reflect_fm,native,'AbsTol',1e-12);
tau_native=-gradient(unwrap(angle(native)),f)/(2*pi);
tau_project=-gradient(unwrap(angle(a.reflect_fm)),f)/(2*pi);
mid=101; expected=103/1500*k(mid)/kz(mid)-97/1500;
verifyEqual(t,tau_native(mid),expected,'AbsTol',1e-10);
verifyEqual(t,tau_project(mid),.008-expected,'AbsTol',1e-10);
verifyGreaterThan(t,abs(tau_native(mid)-tau_project(mid)),1e-3);
end

function testFixedSurfaceNestedGrid(t)
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
p=load_fixed_pm_profile_for_pe_bellhop_validation(fullfile(root,'results','validation', ...
 'bellhop_internal_pm_fixed_realization','fixed_pm_fourier_coefficients.csv'));
x=linspace(-60,60,4097)'; xf=linspace(-60,60,8193)';
[eta,sl,sec]=evaluate_pm_pilot_surface(p,x);
[ef,sf,cf]=evaluate_pm_pilot_surface(p,xf);
verifyEqual(t,eta,ef(1:2:end)); verifyEqual(t,sl,sf(1:2:end)); verifyEqual(t,sec,cf(1:2:end));
verifyEqual(t,max(p.k_rad_per_m),12*2*pi/160,'AbsTol',1e-14);
end

function testActualPhaseAssemblerWithKnownDelay(t)
f=(4000:62.5:8000)'; extra=.0003;
physical=-exp(1i*2*pi*f*(103/1500+extra));
d=struct([]);
for i=1:65
 d(i).frequency_hz=f(i);
 d(i).H_full_dsp=-exp(-1i*2*pi*f(i)*(.004+extra));
 d(i).H_fast_dsp=d(i).H_full_dsp;
 d(i).H_full_reduced=-exp(1i*2*pi*f(i)*extra); d(i).H_fast_reduced=d(i).H_full_reduced;
 d(i).H_full_physical=physical(i); d(i).H_fast_physical=physical(i);
end
p=audit_fk_multifrequency_phase(d);
verifyEqual(t,p.full_dsp_group_delay_s,.004+extra+zeros(65,1),'AbsTol',1e-14);
verifyEqual(t,p.full_physical_group_delay_s,-103/1500-extra+zeros(65,1),'AbsTol',1e-14);
verifyEqual(t,p.cir_full,ifft(-exp(-1i*2*pi*f*(.004+extra))),'AbsTol',1e-14);
[~,i]=max(abs(p.cir_full)); verifyLessThan(t,abs(p.delay_s(i)-.004-extra),1/(65*62.5));
end
