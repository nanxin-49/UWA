function checks=audit_deepsea3d_g1_q10_artifacts()
%AUDIT_DEEPSEA3D_G1_Q10_ARTIFACTS Audit geometry, native CSV and numerical gates.
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
out=fullfile(root,'results','validation','deepsea3d_g1_q10');
a=load(fullfile(out,'q10.mat'),'state');s=a.state;
assert(s.complete && s.global_G1=="OPEN" && ~s.WP2_allowed);
assert(isequal(s.cfg.tx_m,[150.5 0 1500]) && isequal(s.cfg.rx_m,[0 0 5]));
restored=s.cfg;restored.tx_m=s.baseline.cfg.tx_m;assert(isequal(restored,s.baseline.cfg));
assert(s.flat_control.passed && s.flat_control.q==.1 && s.flat_control.specular_x_m==.5);
assert(s.cfg.center_m(1)==0 && isequal(s.cfg.profile,s.baseline.cfg.profile));
x=-12:.0625:12;[xx,yy]=meshgrid(x,x);
[eta0,dx0,dy0]=deepsea3d_surface(xx,yy,s.baseline.cfg.profile);
[eta1,dx1,dy1]=deepsea3d_surface(xx,yy,s.cfg.profile);
surface_difference=max(abs([eta0(:)-eta1(:);dx0(:)-dx1(:);dy0(:)-dy1(:)]));assert(surface_difference==0);
checks=struct('complete',true,'geometry_only_change',true,'fixed_surface_difference',surface_difference, ...
    'flat_control_pass',s.flat_control.passed,'reference_ready',s.reference_ready, ...
    'q10_G1',s.G1_status,'global_G1',"OPEN",'WP2_allowed',false);
if ~s.reference_ready
    assert(isempty(s.fk) && s.G1_status=="INCONCLUSIVE");
else
    g=readtable(fullfile(out,'reference_guards.csv'));assert(all(g{1,:}<1e-4));
    r=s.runs(s.reference_index).result;o=s.observer;
    assert(all(isfinite([r.terms;r.partials;r.boundary_coefficient_norms])));
    assert(all(diff(r.boundary_coefficient_norms(end-4:end))<0) && all(diff(abs(r.terms(end-4:end)))<0));
    assert(o.term_match_error<1e-12 && o.norm_match_error<1e-12 && r.forcing_series_relative<1e-8);
    assert(norm(o.terms-r.terms)/norm(r.terms)<1e-12);
    fk=readtable(fullfile(out,'fk_runs.csv'));refs=readtable(fullfile(out,'reference_runs.csv'));
    assert(height(fk)==5 && height(refs)==numel(s.runs));
    F=fk.reflect_real(2)+1i*fk.reflect_imag(2);D=fk.delta_real(2)+1i*fk.delta_imag(2);
    R=refs.reflect_real(s.reference_index)+1i*refs.reflect_imag(s.reference_index);
    I=refs.delta_real(s.reference_index)+1i*refs.delta_imag(s.reference_index);
    assert(abs(F-s.fk(2).reflect_native)/abs(F)<1e-13 && abs(I-r.rough_correction)/abs(I)<1e-13);
    ratio=[F/R,D/I];native_errors=abs(ratio-1);phases=angle(ratio);rho=abs(ratio);
    identity=max(abs(native_errors.^2-((rho-1).^2+4*rho.*sin(phases/2).^2)));assert(identity<1e-16);
    c=summarize_deepsea3d_q10_comparison(s);assert(isequaln(c,s.comparison));
    primary=readtable(fullfile(out,'geometry_comparison.csv'));
    assert(max(abs(native_errors-[primary.reflected_complex(2),primary.increment_complex(2)]))<1e-13);
    assert(max(abs(phases-[primary.reflected_phase_rad(2),primary.increment_phase_rad(2)]))<1e-13);
    old=readtable(fullfile(root,'results','validation','deepsea3d_alpha25_domain','physical_metrics.csv'));
    assert(abs(primary.reflected_complex(1)-old.reflected_complex)<1e-14);
    assert(abs(primary.increment_complex(1)-old.increment_complex)<1e-14);
    envelopes=readtable(fullfile(out,'error_envelopes.csv'));
    assert(isequal(envelopes.ranges_disjoint,[c.envelopes.ranges_disjoint].'));
    num=s.fk_guards;
    classified=classify_deepsea3d_amplitude(true,max([num.grid,num.aperture_12_16,num.aperture_16_24]),num.gauss,c.metrics_q10);
    assert(classified==s.G1_status);
    % Counterfactual guard: an unresolved reference cannot become physical FAIL.
    unsafe=struct('reflected_complex',.2,'increment_complex',.3);
    assert(classify_deepsea3d_amplitude(false,0,0,unsafe)=="INCONCLUSIVE");
    assert(classify_deepsea3d_amplitude(true,.002,0,unsafe)=="INCONCLUSIVE");
    checks.CSV_native_pressure_and_phase_matches=true;checks.amplitude_phase_identity_abs_error=identity;
    checks.q0_saved_metrics_unchanged=true;checks.empirical_envelopes_recomputed=true;
    checks.reference_failures_never_physical_fail=true;
end
fid=fopen(fullfile(out,'saved_data_checks.json'),'w');assert(fid>=0);
clean=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(checks,PrettyPrint=true));disp(checks);
end
