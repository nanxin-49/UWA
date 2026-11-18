function audit = audit_d0_incident_case(case_dir,cfg)
%AUDIT_D0_INCIDENT_CASE Independent saved field checks; writes failures first.
v=load(fullfile(case_dir,'field.mat'),'data'); d=v.data; s=d.source;
checks=struct('name',{},'value',{},'limit',{},'passed',{});
g=cfg.gates; x=cfg.audit_x_m(:); z=zeros(size(x));
[~,ix]=ismember(x,d.x_m); assert(all(ix>0));
actual=[d.p_f(ix) d.px_f(ix) d.pz_f(ix)];
[p,px,pz]=baffled_strip_piston_field_vertical(x,z,s,128); ref=[p px pz];
add('aperture64_128_pressure',rel(actual(:,1),p),g.source_relative);
add('aperture64_128_gradient',rel(actual(:,2:3),ref(:,2:3)),g.source_relative);
AS16=d0_incident_spectrum_reference(x,z,s,16,8);
AS32=d0_incident_spectrum_reference(x,z,s,32,4);
add('AS16_32_pressure',rel(AS16.total(:,1),AS32.total(:,1)),g.source_relative);
add('AS16_32_gradient',rel(AS16.total(:,2:3),AS32.total(:,2:3)),g.source_relative);
add('aperture_AS_pressure',rel(p,AS32.total(:,1)),g.source_relative);
add('aperture_AS_gradient',rel(ref(:,2:3),AS32.total(:,2:3)),g.source_relative);
% Regression against the unmodified qualified adaptive P0 integrator at short distance.
xx=[-.2;0;.2]; zz=s.z_source_m-.1;
[p0,x0,z0]=piston_independent_angular_spectrum_vertical(xx,zz,s);
b0=d0_incident_spectrum_reference(xx,zz,s,32,4);
add('segmented_legacy_P0',rel(b0.total,[p0 x0 z0]),g.source_relative);
for h=[1e-4 5e-5]
    f=@(xx,zz)baffled_strip_piston_field_vertical(xx,zz,s,128);
    fx=(-f(x+2*h,z)+8*f(x+h,z)-8*f(x-h,z)+f(x-2*h,z))/(12*h);
    fz=(-f(x,z+2*h)+8*f(x,z+h)-8*f(x,z-h)+f(x,z-2*h))/(12*h);
    add(sprintf('five_point_gradient_h_%g',h),rel([fx fz],ref(:,2:3)),g.gradient_relative);
end
T=readtable(fullfile(case_dir,'points.csv'));
add('MAT_CSV',rel(T{:,:},d.points{:,:}),g.csv_relative);
pw=d.plane_wave;
add('plane_wave_error_decomposition',max(abs(pw.complex_relative_L2.^2- ...
    pw.amplitude_component_relative.^2-pw.phase_component_relative.^2)),g.fft_relative);
for j=1:numel(d.window_spectra)
    w=d.window_spectra{j};
    add(sprintf('FFT_Parseval_%d',j),w.parseval_relative,g.fft_relative);
    add(sprintf('FFT_inverse_%d',j),w.fft_inverse_relative,g.fft_relative);
    sel=unique(round(linspace(1,numel(w.kx_rad_per_m),11)));
    xi=d.x_m>=-w.halfwidth_m & d.x_m<w.halfwidth_m;
    xx=d.x_m(xi); f=d.p_f(xi);
    if w.hann, f=f.*(.5*(1+cos(pi*xx/w.halfwidth_m))); end
    direct=cfg.dx_m*(exp(-1i*w.kx_rad_per_m(sel)*xx.')*f);
    add(sprintf('FFT_direct_DFT_%d',j),norm(direct-w.spectrum_f(sel))/norm(w.spectrum_f),g.dft_relative);
end
add('adjacent_phase_sampling',max(abs(angle(d.p_f(2:end).*conj(d.p_f(1:end-1))))),g.max_phase_step_rad);
if isfile(fullfile(case_dir,'fine.mat'))
    v=load(fullfile(case_dir,'fine.mat'),'fine'); fine=v.fine;
    add('fine_common_pressure',rel(fine.p_f(1:2:end),d.p_f),g.common_point_relative);
    add('fine_common_gradients',rel([fine.px_f(1:2:end) fine.pz_f(1:2:end)],[d.px_f d.pz_f]),g.common_point_relative);
    predicted=interp1(d.x_m,abs(d.p_f),fine.x_m,'spline').* ...
        exp(1i*interp1(d.x_m,unwrap(angle(d.p_f)),fine.x_m,'spline'));
    add('fine_amplitude_phase_spline',rel(predicted,fine.p_f),g.sampling_relative);
    fc=d0_incident_diagnostics(fine.x_m,fine.p_f,fine.px_f,fine.pz_f,s,cfg);
    for j=[2 4] % Hann only; rectangular edge dependence retained as diagnostic.
        a=d.window_spectra{j}; b=fc.window_spectra{j};
        [~,ids]=ismember(round(a.kx_rad_per_m,10),round(b.kx_rad_per_m,10));
        sel=abs(a.kx_rad_per_m)<.95*(2*pi*s.frequency_hz/s.c_mps);
        add(sprintf('fine_Hann_FFT_%d',j),rel(a.spectrum_f(sel),b.spectrum_f(ids(sel))),g.sampling_relative);
    end
end
r=load(fullfile(case_dir,'rough.mat'),'rough'); r=r.rough;
add('rough_MAT_CSV',rel(readmatrix(fullfile(case_dir,'rough.csv')),r.points{:,:}),g.csv_relative);
% Independent scalar Fourier/C2 geometry reconstruction from frozen coefficients.
c=readtable(fullfile(case_dir,'..','..','inputs','k0','coefficients.csv'));
xr=r.x_m; raw=zeros(size(xr)); slope=raw;
for j=1:height(c)
    q=c.k_rad_per_m(j); aa=c.cos_depth_coeff_m(j); bb=c.sin_depth_coeff_m(j);
    raw=raw+aa*cos(q*xr)+bb*sin(q*xr);
    slope=slope-aa*q*sin(q*xr)+bb*q*cos(q*xr);
end
chi=ones(size(xr)); cp=zeros(size(xr)); mid=abs(xr)>42; t=(abs(xr(mid))-42)/8;
chi(mid)=1-10*t.^3+15*t.^4-6*t.^5;
cp(mid)=(-30*t.^2+60*t.^3-30*t.^4).*sign(xr(mid))/8;
add('K0_geometry_depth_m',max(abs(r.z_m-raw.*chi)),g.geometry_absolute_m);
add('K0_geometry_slope',max(abs(r.depth_slope-(slope.*chi+raw.*cp))),g.geometry_absolute_m);
sel=unique(round(linspace(1,numel(xr),13))); xx=xr(sel); zz=r.z_m(sel);
h=1e-4; f=@(xx,zz)baffled_strip_piston_field_vertical(xx,zz,s,128);
fx=(-f(xx+2*h,zz)+8*f(xx+h,zz)-8*f(xx-h,zz)+f(xx-2*h,zz))/(12*h);
fz=(-f(xx,zz+2*h)+8*f(xx,zz+h)-8*f(xx,zz-h)+f(xx,zz-2*h))/(12*h);
normal=(-r.depth_slope(sel).*fx+fz)./hypot(1,r.depth_slope(sel));
add('rough_normal_finite_difference',rel(normal,r.dn_water_f(sel)),g.gradient_relative);
ar=d0_incident_spectrum_reference(xx,zz,s,32,4);
add('rough_AS_pressure',rel(r.p_f(sel),ar.total(:,1)),g.source_relative);
add('rough_AS_gradient',rel([r.px_f(sel) r.pz_f(sel)],ar.total(:,2:3)),g.source_relative);
% Physical representation round-trip using v2; each point has explicit vertical carrier.
phase_error=0;
for j=1:numel(p)
    geom=struct('z_tx',s.z_source_m,'z_rx',0,'c0',s.c_mps);
    reduced=struct('direct_f',p(j)*exp(-1i*2*pi*s.frequency_hz*s.z_source_m/s.c_mps));
    [back,meta]=apply_pe_channel_phase_reference_vertical(s.frequency_hz,reduced,geom,'absolute_physical');
    assert(strcmp(meta.phase_contract_version,'2.0.0') && ~meta.conjugate_input);
    phase_error=max(phase_error,abs(back.direct_f-p(j)));
end
add('phase_v2_roundtrip_Pa',phase_error,g.phase_absolute_Pa);
save(fullfile(case_dir,'independent_spectrum.mat'),'AS16','AS32','x','s','-v7');
audit=struct('status','AUDITED_PASS','checks',checks,'failed_checks',sum(~[checks.passed]), ...
    'branch_pressure_relative',norm(AS32.evanescent(:,1))/norm(AS32.total(:,1)), ...
    'evanescent_semantics','coherent branch norm on audit points, not energy fraction', ...
    'next_allowed_step','next D0 depth only; no reflection solve');
if audit.failed_checks>0, audit.status='FAILED_STOPPED'; end
pm_sensitivity_io('write',fullfile(case_dir,'audit.json'),audit);
writetable(struct2table(checks),fullfile(case_dir,'audit_checks.csv'));
assert(audit.failed_checks==0,'D0 numerical gates failed; inspect saved audit_checks.csv');
    function add(name,value,limit)
        checks(end+1)=struct('name',name,'value',value,'limit',limit,'passed',isfinite(value)&&value<limit);
    end
end
function value=rel(a,b)
assert(isequal(size(a),size(b)) && isequal(isnan(a),isnan(b)),'Shape/NaN mask mismatch');
ii=~isnan(b); value=norm(a(ii)-b(ii))/max(norm(b(ii)),realmin);
end
