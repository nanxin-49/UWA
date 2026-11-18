function receipt = review_d0_saved(source_dir,out)
%REVIEW_D0_SAVED Recheck saved arrays only. No source/propagation/SSH calls.
arguments
    source_dir (1,:) char
    out (1,:) char
end
assert(~isfolder(out)); mkdir(out); t0=tic;
diary(fullfile(out,'review_diary.log')); dg=onCleanup(@()diary('off'));
repo=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repo,'scripts','experiments'),fullfile(repo,'scripts','reporting'));
manifest=jsondecode(fileread(fullfile(source_dir,'retrieval','manifest.json')));
pm_sensitivity_io('verify',source_dir,manifest);
cfg=jsondecode(fileread(fullfile(source_dir,'config.json')));
checks=struct('depth_m',{},'name',{},'value',{},'limit',{},'passed',{});
summaries=cell(5,1); rough_stats=cell(5,1); spectral=cell(0,1); pw_tables=cell(5,1);
saved=cell(5,1); original_fig_lines=0;
for j=1:5
    depth=cfg.depths_m(j); folder=fullfile(source_dir,'output',sprintf('depth_%04dm',depth));
    v=load(fullfile(folder,'field.mat'),'data'); d=v.data; saved{j}=d;
    x=d.x_m; p=d.p_f; px=d.px_f; pz=d.pz_f; k=2*pi*4000/1500; center=4401;
    assert(numel(x)==8801 && x(center)==0 && d.source.z_source_m==depth);
    check('MAT_CSV_points',relative(readmatrix(fullfile(folder,'points.csv')),d.points{:,:}),1e-12);
    check('complex_center_normalization',relative(d.center_normalized_f,p/p(center)),1e-12);
    check('real_RMS_normalization',relative(d.rms_normalized_f,p/sqrt(mean(abs(p).^2))),1e-12);
    phase_unwrapped=unwrap(angle(p));
    check('relative_phase',max(abs(d.relative_phase_rad-(phase_unwrapped-phase_unwrapped(center)))),1e-10);
    % Independent derivative: 16th-order centered stencil on saved p(x), interior only.
    offsets=-8:8; weights=zeros(17,1);
    for mm=1:8
        weights(9+mm)=(-1)^(mm+1)*factorial(8)^2/(mm*factorial(8-mm)*factorial(8+mm));
        weights(9-mm)=-weights(9+mm);
    end
    ids=(9:numel(x)-8)'; derivative=complex(zeros(size(ids)));
    for q=1:17, derivative=derivative+weights(q)*p(ids+offsets(q))/.05; end
    check('saved_pressure_dx_16th',relative(derivative,px(ids)),1e-6);
    a=load(fullfile(folder,'independent_spectrum.mat'));
    [~,ix]=ismember(a.x,x); assert(all(ix>0));
    check('stored_AS_pressure',relative(p(ix),a.AS32.total(:,1)),1e-7);
    check('stored_AS_gradient',relative([px(ix) pz(ix)],a.AS32.total(:,2:3)),1e-7);
    check('stored_AS_branch_sum',relative(a.AS32.total,a.AS32.propagating+a.AS32.evanescent),1e-12);
    check('phase_kx_from_gradient',relative(d.local_kx_rad_per_m,imag(px./p)),1e-12);
    check('local_direction',max(abs(d.local_angle_deg-atan2d(imag(px./p),-imag(pz./p)))),1e-10);
    if isfile(fullfile(folder,'fine.mat'))
        f=load(fullfile(folder,'fine.mat'),'fine'); f=f.fine;
        check('fine_shared_pressure',relative(f.p_f(1:2:end),p),1e-12);
        check('fine_shared_gradient',relative([f.px_f(1:2:end) f.pz_f(1:2:end)],[px pz]),1e-12);
    end
    r=load(fullfile(folder,'rough.mat'),'rough'); r=r.rough;
    check('rough_MAT_CSV',relative(readmatrix(fullfile(folder,'rough.csv')),r.points{:,:}),1e-12);
    check('water_normal_contraction',relative(r.dn_water_f,(-r.depth_slope.*r.px_f+r.pz_f)./hypot(1,r.depth_slope)),1e-12);
    if j==1, frozen_geometry=[r.x_m r.z_m r.depth_slope];
    else, check('identical_rough_geometry',max(abs([r.x_m r.z_m r.depth_slope]-frozen_geometry),[],'all'),1e-14); end
    theta=atan2d(imag(r.px_f./r.p_f),-imag(r.pz_f./r.p_f));
    normal_tilt=atan(r.depth_slope)*180/pi;
    local_incidence=abs(theta-normal_tilt);
    rough_stats{j}=struct('depth_m',depth,'point_count',numel(r.x_m),'min_x_m',min(r.x_m),'max_x_m',max(r.x_m), ...
        'pressure_rms_Pa',sqrt(mean(abs(r.p_f).^2)),'dn_water_rms_Pa_per_m',sqrt(mean(abs(r.dn_water_f).^2)), ...
        'center_ray_abs_angle_max_deg',atand(50/depth),'normal_incidence_min_deg',min(local_incidence), ...
        'normal_incidence_max_deg',max(local_incidence),'normal_tilt_rms_deg',sqrt(mean(normal_tilt.^2)));
    Q=d.plane_wave; check('plane_wave_CSV',relative(readmatrix(fullfile(folder,'plane_wave.csv')),Q{:,:}),1e-12);
    regenerated=zeros(height(Q),5);
    for q=1:height(Q)
        x0=Q.anchor_x_m(q); [~,anchor]=min(abs(x-x0)); sel=abs(x-x0)<=Q.patch_width_m(q)/2;
        plane=p(anchor)*exp(1i*imag(px(anchor)/p(anchor))*(x(sel)-x0));
        diffphase=angle(p(sel).*conj(plane));
        E=norm(p(sel)-plane)/norm(p(sel)); EA=norm(abs(p(sel))-abs(plane))/norm(p(sel));
        EP=sqrt(sum(4*abs(p(sel)).*abs(plane).*sin(diffphase/2).^2))/norm(p(sel));
        regenerated(q,:)=[E max(abs(diffphase)) EA EP abs(E^2-EA^2-EP^2)];
    end
    check('plane_wave_complex_error',relative(Q.complex_relative_L2,regenerated(:,1)),1e-12);
    check('plane_wave_amplitude_phase_decomposition',max(regenerated(:,5)),1e-12);
    Q.depth_m=repmat(depth,height(Q),1); pw_tables{j}=Q;
    s=d.continuous_spectrum;
    if j==1, continuous_modulus=abs(s.prop_density_f);
    else, check('propagating_spectrum_modulus_invariant',relative(abs(s.prop_density_f),continuous_modulus),1e-12); end
    edge=x==200; support=x==50;
    span=abs(x)<=200; rms_target=sqrt(mean(abs(p(span)).^2));
    phase=d.relative_phase_rad; excess=hypot(depth,x)-depth;
    geomphase=k*excess;
    sph= sqrt(depth./hypot(depth,x)); u=k*.05*x./hypot(depth,x); sincweight=ones(size(u));
    ii=u~=0; sincweight(ii)=sin(u(ii))./u(ii);
    asymptotic_ratio=sph.*sincweight;
    summaries{j}=struct('depth_m',depth,'center_peak_Pa',abs(p(center)),'rms_full_line_Pa',sqrt(mean(abs(p).^2)), ...
        'rms_target_200_Pa',rms_target,'amplitude_ratio_200',abs(p(edge)/p(center)), ...
        'amplitude_drop_200_percent',100*(1-abs(p(edge)/p(center))), ...
        'amplitude_ratio_220',abs(p(end)/p(center)),'phase_200_rad',phase(edge),'phase_220_rad',phase(end), ...
        'phase_50_rad',phase(support),'phase_curvature_axis_rad_per_m2',d.phase_curvature_rad_per_m2(center), ...
        'geometric_angle_200_deg',atand(200/depth),'local_angle_200_deg',d.local_angle_deg(edge), ...
        'wavefront_radius_axis_m',d.metrics.axis_phase_front_radius_m, ...
        'full_line_axial_plane_error',norm(p-p(center))/norm(p), ...
        'target_200_axial_plane_error',norm(p(span)-p(center))/norm(p(span)), ...
        'phase_geometry_max_abs_rad',max(abs(phase-geomphase)), ...
        'strip_farfield_amplitude_shape_max_abs',max(abs(abs(p/p(center))-asymptotic_ratio)), ...
        'point_source_spreading_edge',sph(edge),'strip_directivity_edge',sincweight(edge), ...
        'pressure_symmetry_relative',relative(p,flipud(p)), ...
        'evanescent_branch_pressure_norm_ratio',norm(a.AS32.evanescent(:,1))/norm(a.AS32.total(:,1)));
    for q=1:4
        w=d.window_spectra{q}; ii=x>=-w.halfwidth_m & x<w.halfwidth_m; xx=x(ii); pp=p(ii); n=numel(xx);
        taper=ones(n,1); if w.hann, taper=.5*(1+cos(pi*xx/w.halfwidth_m)); end
        probes=unique(round(linspace(1,n,31)));
        direct=.05*(exp(-1i*w.kx_rad_per_m(probes)*xx.')*(pp.*taper));
        check(sprintf('FFT_direct_DFT_%d',q),norm(direct-w.spectrum_f(probes))/norm(w.spectrum_f),1e-10);
        csv=readmatrix(fullfile(folder,sprintf('window_spectrum_%d.csv',q)));
        check(sprintf('FFT_MAT_CSV_%d',q),relative(csv(:,1:5),[w.kx_rad_per_m real(w.spectrum_f) imag(w.spectrum_f) abs(w.spectrum_f) w.phase_rad]),1e-12);
        power=abs(w.spectrum_f).^2; [sortedq,order]=sort(abs(w.kx_rad_per_m)); cumulative=cumsum(power(order))/sum(power);
        q95=sortedq(find(cumulative>=.95,1));
        spectral{end+1}=struct('depth_m',depth,'halfwidth_m',w.halfwidth_m,'hann',w.hann, ...
            'kx95_abs_rad_per_m',q95,'mapped_angle95_deg',asind(min(1,q95/k)), ...
            'window_outside_band_power_fraction',sum(power(abs(w.kx_rad_per_m)>k))/sum(power)); %#ok<AGROW>
    end
end
summaries=vertcat(summaries{:}); rough_stats=vertcat(rough_stats{:}); spectral=vertcat(spectral{:});
% Verify all 60 original FIG lines against saved arrays, preserving multiplicity.
stems={'F01_absolute_pressure','F02_normalized_amplitude','F03_relative_phase','F04_complex_angular_spectrum','F05_incidence_angles','F06_gradients_curvature'};
for fi=1:6
    expected=struct('x',{},'y',{});
    for j=1:5
        d=saved{j};
        switch fi
            case 1, xx={d.x_m}; yy={abs(d.p_f)};
            case 2, xx={d.x_m,d.x_m}; yy={abs(d.center_normalized_f),abs(d.rms_normalized_f)};
            case 3, xx={d.x_m}; yy={d.relative_phase_rad};
            case 4
                w=d.window_spectra{4}; q=d.continuous_spectrum;
                xx={w.kx_rad_per_m,w.kx_rad_per_m,q.prop_kx_rad_per_m,q.prop_kx_rad_per_m};
                yy={abs(w.spectrum_f),w.phase_rad,abs(q.prop_density_f),q.prop_phase_unwrapped_rad};
            case 5, xx={d.x_m,d.x_m}; yy={d.geometric_angle_deg,d.local_angle_deg};
            case 6, xx={d.x_m,d.x_m}; yy={hypot(abs(d.px_f),abs(d.pz_f)),d.phase_curvature_rad_per_m2};
        end
        for q=1:numel(xx), expected(end+1)=struct('x',xx{q}(:),'y',yy{q}(:)); end %#ok<AGROW>
    end
    f=openfig(fullfile(source_dir,'output','figures',[stems{fi} '.fig']),'new','invisible'); fg=onCleanup(@()close(f));
    lines=findall(f,'Type','line'); used=false(size(expected));
    for q=1:numel(lines)
        match=0;
        for z=1:numel(expected)
            if ~used(z) && isequaln(lines(q).XData(:),expected(z).x) && isequaln(lines(q).YData(:),expected(z).y), match=z; break; end
        end
        assert(match>0,'Unexpected original FIG line'); used(match)=true;
    end
    assert(all(used) && numel(lines)==numel(expected)); original_fig_lines=original_fig_lines+numel(lines); clear fg
end
writetable(struct2table(summaries),fullfile(out,'depth_summary_review.csv'));
writetable(struct2table(rough_stats),fullfile(out,'rough_summary_review.csv'));
writetable(struct2table(spectral),fullfile(out,'window_spectrum_diagnostics.csv'));
writetable(vertcat(pw_tables{:}),fullfile(out,'plane_wave_review.csv'));
writetable(struct2table(checks),fullfile(out,'local_saved_checks.csv'));
pm_sensitivity_io('write',fullfile(out,'analysis.json'),struct('depths',summaries,'rough',rough_stats,'window_spectra',spectral));
assert(all([checks.passed]),'Saved-data audit failed');
figure_receipt=render_d0_saved_review(saved,out);
pm_sensitivity_io('verify',source_dir,manifest);
receipt=struct('status','SAVED_DATA_REVIEW_PASS','elapsed_s',toc(t0),'matlab_release',version('-release'), ...
    'checks',numel(checks),'failed_checks',sum(~[checks.passed]),'original_FIG_lines',original_fig_lines, ...
    'new_figures',figure_receipt,'source_manifest_files',numel(manifest), ...
    'source_SHA_unchanged',true,'new_source_or_solver_calls',0,'SSH_calls',0, ...
    'scope','Stored-array audits, finite differences, FFT/DFT, geometry and posthoc analytic interpretation only');
pm_sensitivity_io('write',fullfile(out,'receipt.json'),receipt); disp(jsonencode(receipt));
    function check(name,value,limit)
        checks(end+1)=struct('depth_m',depth,'name',name,'value',value,'limit',limit,'passed',isfinite(value)&&value<limit);
    end
end
function v=relative(a,b)
assert(isequal(size(a),size(b)) && isequal(isnan(a),isnan(b)),'Array shape/NaN mask differs');
ii=~isnan(b); v=norm(a(ii)-b(ii))/max(norm(b(ii)),realmin);
end
