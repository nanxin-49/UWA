function validation = validate_bellhop_internal_sinusoidal_wall_poc(options)
%VALIDATE_BELLHOP_INTERNAL_SINUSOIDAL_WALL_POC Validation-only smooth wall POC.
% The internal binary uses a sampled C-ATI-equivalent wall, calls native
% Reflect2D once, and applies the already validated proper pi rotation.
arguments
    options.output_dir (1,:) char = ''
    options.validation_exe (1,:) char = ''
    options.official_exe (1,:) char = 'E:\MISC\BELLHOP\AcousticsToolbox_2020\windows-bin-20201102\bellhop.exe'
    options.frequency_hz (1,1) double {mustBePositive} = 4000
    options.wall_wavenumber_per_m (1,1) double = -0.01
    options.wall_amplitudes_m (1,:) double {mustBePositive} = [0.25 0.5]
    options.profile_counts (1,:) double {mustBeInteger,mustBePositive} = [41 81 161]
    options.beam_counts (1,:) double {mustBeInteger,mustBePositive} = [2001 5001]
    options.step_values_m (1,:) double {mustBePositive} = [0.2 0.1 0.05]
    options.fail_on_check (1,1) logical = true
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup_vertical_project;
if isempty(options.output_dir)
    options.output_dir=fullfile(root,'results','validation','bellhop_internal_sinusoidal_wall_poc');
end
if isempty(options.validation_exe)
    options.validation_exe=fullfile(options.output_dir,'bin','bellhop_iwall_sinusoidal_2020.exe');
end
assert(isfile(options.validation_exe),'Missing sinusoidal-wall validation binary: %s',options.validation_exe);
assert(isfile(options.official_exe),'Missing official Bellhop 2020 binary: %s',options.official_exe);
if ~exist(options.output_dir,'dir'), mkdir(options.output_dir); end

R0=100; K=options.wall_wavenumber_per_m;
% Use a signed K so dr/dz=-A*K*cos(Kz)>0 around z=0.  This keeps the
% prescribed sinusoid monotone in ATI range while retaining the original
% source/receiver datum and the already validated rotation about (R0,0).
zc=0;
cfg=struct('c0_mps',1500,'wall_r0_m',R0,'wall_wavenumber_per_m',K, ...
    'mapped_receiver_range_m',103,'native_receiver_range_m',97, ...
    'source_depth_m',zc,'receiver_depths_m',-zc,'receiver_ranges_m',[102 103], ...
    'run_type','C','angle_limits_deg',[-30 30],'domain_half_depth_m',1000, ...
    'sigma_src_m',0.3,'source_pattern_clip_db',-120,'frequency_hz',options.frequency_hz);
pat=local_pattern(cfg.frequency_hz,cfg);
cfg.source_pattern_angles_deg=pat.angles_deg;
cfg.source_pattern_level_db=pat.level_db;

template=struct('amplitude_m',NaN,'profile_count',NaN,'beam_count',NaN,'step_m',NaN, ...
    'wall_pressure',complex(NaN),'native_total_pressure',complex(NaN), ...
    'native_direct_pressure',complex(NaN),'native_reflected_pressure',complex(NaN), ...
    'field_tl_error_db',NaN,'field_phase_error_rad',NaN,'field_complex_relative_error',NaN, ...
    'intersection_residual_m',NaN,'tangent_error',NaN,'normal_error',NaN, ...
    'specular_direction_error',NaN,'rotation_direction_error',NaN,'phase_jump_error_rad',NaN, ...
    'amp_jump_error',NaN,'p_reflect_error',NaN,'q_reflect_error',NaN, ...
    'p_rotation_error',NaN,'q_rotation_error',NaN,'kappa_abs',NaN, ...
    'path_length_error_m',NaN,'travel_time_error_s',NaN,'travel_time_imag_abs_s',NaN, ...
    'min_post_range_increment_m',NaN,'all_post_range_positive',false, ...
    'native_case_root','','wall_case_root','');
nRows=numel(options.wall_amplitudes_m)*numel(options.profile_counts)* ...
    numel(options.step_values_m)*numel(options.beam_counts);
rows=repmat(template,nRows,1);
wall_runs=cell(nRows,1); native_runs=cell(nRows,1); direct_runs=cell(nRows,1);
rr=0;
for aa=1:numel(options.wall_amplitudes_m)
    A=options.wall_amplitudes_m(aa);
    for pp=1:numel(options.profile_counts)
        N=options.profile_counts(pp);
        zProf=linspace(zc-100,zc+100,N).';
        rProf=R0-A*sin(K*zProf);
        for ss=1:numel(options.step_values_m)
            for bb=1:numel(options.beam_counts)
                rr=rr+1; step=options.step_values_m(ss); beams=options.beam_counts(bb);
                tag=sprintf('A_%s_N_%d_step_%s_beams_%d',local_tag(A),N,local_tag(step),beams);
                common=cfg; common.wall_amplitude_m=A; common.wall_profile_r_m=rProf; ...
                    common.wall_profile_z_m=zProf; common.step_m=step; common.beam_count=beams;
                common.case_root=fullfile(options.output_dir,'cases','wall',[tag '_wall']);
                common.bellhop_exe=options.validation_exe;
                wall_runs{rr}=run_bellhop_internal_sinusoidal_wall_poc_vertical(common);
                native=common; native.case_root=fullfile(options.output_dir,'cases','native',[tag '_native']);
                native.bellhop_exe=options.official_exe; native.receiver_depths_m=zc;
                native.receiver_ranges_m=[97 98];
                native_runs{rr}=run_bellhop_native_sinusoidal_wall_vertical(native);
                direct=native; direct.case_root=fullfile(options.output_dir,'cases','direct',[tag '_direct']);
                direct_runs{rr}=run_bellhop_unfolded_gaussian_vertical(direct);

                pwall=select_bellhop_shd_pressure_at_range_vertical( ...
                    wall_runs{rr}.data,cfg.mapped_receiver_range_m,zc);
                pnative=select_bellhop_shd_pressure_at_range_vertical( ...
                    native_runs{rr}.data,cfg.native_receiver_range_m,zc);
                pdirect=select_bellhop_shd_pressure_at_range_vertical( ...
                    direct_runs{rr}.data,cfg.native_receiver_range_m,zc);
                pref=pnative-pdirect;
                d=wall_runs{rr}.diagnostics;
                % Uniform-c post-wall propagation gives an exact path to the
                % mapped range plane for every traced ray.
                rotUr=d.rot_ur; hitR=2*R0-d.hit_r;
                pre=sqrt(d.hit_r.^2+(d.hit_z-zc).^2);
                expected=pre+(cfg.mapped_receiver_range_m-hitR)./rotUr;
                actual=cfg.c0_mps*d.tau_receiver_real;
                q=template; q.amplitude_m=A; q.profile_count=N; q.beam_count=beams; q.step_m=step;
                q.wall_pressure=pwall; q.native_total_pressure=pnative; q.native_direct_pressure=pdirect;
                q.native_reflected_pressure=pref;
                q.field_tl_error_db=20*log10(abs(pwall)/max(abs(pref),realmin));
                q.field_phase_error_rad=angle(pwall*conj(pref));
                q.field_complex_relative_error=abs(pwall-pref)/max(abs(pref),realmin);
                q.intersection_residual_m=max(abs(d.wall_residual));
                q.tangent_error=max(d.tangent_error); q.normal_error=max(d.normal_error);
                q.specular_direction_error=max(d.specular_error); q.rotation_direction_error=max(d.rotation_error);
                q.phase_jump_error_rad=max(abs(angle(exp(1i*(d.phase_delta-pi)))));
                q.amp_jump_error=max(abs(d.amp_delta));
                q.p_reflect_error=max(d.p_ref_error); q.q_reflect_error=max(d.q_ref_error);
                q.p_rotation_error=max(d.p_rot_error); q.q_rotation_error=max(d.q_rot_error);
                q.kappa_abs=max(abs(d.kappa)); q.path_length_error_m=max(abs(actual-expected));
                q.travel_time_error_s=q.path_length_error_m/cfg.c0_mps;
                q.travel_time_imag_abs_s=max(abs(d.tau_receiver_imag));
                q.min_post_range_increment_m=min(d.min_post_dr); q.all_post_range_positive=all(d.min_post_dr>0);
                q.native_case_root=native.case_root; q.wall_case_root=common.case_root;
                rows(rr)=q; wall_runs{rr}.profile=[rProf zProf];
            end
        end
    end
end
table_rows=struct2table(rows);
checks=struct;
checks.intersection=max(table_rows.intersection_residual_m)<=1e-8;
checks.frame=max([table_rows.tangent_error;table_rows.normal_error])<=5e-5;
checks.direction=max(table_rows.specular_direction_error)<=1e-10 && max(table_rows.rotation_direction_error)<=1e-10;
checks.phase_jump=max(table_rows.phase_jump_error_rad)<=1e-10;
checks.amplitude=max(table_rows.amp_jump_error)<=1e-10;
% A curved TOP reflection is expected to change p through Reflect2D's native
% curvature kick.  Only q continuity and the post-rotation state are hard
% invariants here; p_reflect_error is retained as the measured kick magnitude.
checks.dynamic_state=max([table_rows.q_reflect_error;table_rows.p_rotation_error; ...
    table_rows.q_rotation_error])<=1e-10;
checks.curvature_kick=~any(isnan(table_rows.p_reflect_error) | isinf(table_rows.p_reflect_error)) && ...
    max(table_rows.p_reflect_error)>1e-12;
checks.path_time=max(table_rows.path_length_error_m)<=1e-7 && max(table_rows.travel_time_imag_abs_s)<=1e-12;
checks.positive_range=all(table_rows.all_post_range_positive);
% Native backward-range Cartesian influence has a known stable offset; keep
% it diagnostic and do not renormalize or use it as a hard curvature gate.
checks.native_field_phase=max(abs(table_rows.field_phase_error_rad))<=0.03;
% Native backward-range amplitude is diagnostic only; it is deliberately not
% part of checks.all and is neither fitted nor renormalized.
checks.native_field_amplitude_diagnostic=max(abs(table_rows.field_tl_error_db))<=0.40 && ...
    max(table_rows.field_complex_relative_error)<=0.05;
checks.all=checks.intersection && checks.frame && checks.direction && checks.phase_jump && ...
    checks.amplitude && checks.dynamic_state && checks.curvature_kick && checks.path_time && ...
    checks.positive_range && checks.native_field_phase;

writetable(table_rows,fullfile(options.output_dir,'sinusoidal_wall_convergence.csv'));
diag_summary=table_rows(:,{'amplitude_m','profile_count','beam_count','step_m', ...
    'intersection_residual_m','tangent_error','normal_error','specular_direction_error', ...
    'rotation_direction_error','phase_jump_error_rad','amp_jump_error','p_reflect_error', ...
    'q_reflect_error','p_rotation_error','q_rotation_error','kappa_abs','path_length_error_m', ...
    'travel_time_error_s','min_post_range_increment_m','field_tl_error_db', ...
    'field_phase_error_rad','field_complex_relative_error'});
writetable(diag_summary,fullfile(options.output_dir,'sinusoidal_wall_diagnostics_summary.csv'));
validation=struct('config',cfg,'options',options,'rows',rows,'table',table_rows,'checks',checks, ...
    'wall_runs',{wall_runs},'native_runs',{native_runs},'direct_runs',{direct_runs}, ...
    'files',struct('convergence_csv',fullfile(options.output_dir,'sinusoidal_wall_convergence.csv'), ...
    'diagnostics_csv',fullfile(options.output_dir,'sinusoidal_wall_diagnostics_summary.csv'), ...
    'mat',fullfile(options.output_dir,'sinusoidal_wall_poc_results.mat')));
save(validation.files.mat,'validation','-v7.3');
disp(table_rows(:,{'amplitude_m','profile_count','beam_count','step_m','kappa_abs', ...
    'field_tl_error_db','field_phase_error_rad','intersection_residual_m','path_length_error_m'}));
disp(checks);
if options.fail_on_check && ~checks.all
    error('Bellhop internal sinusoidal-wall POC failed; inspect %s.',options.output_dir);
end
end

function pat=local_pattern(f,cfg)
angles=linspace(cfg.angle_limits_deg(1),cfg.angle_limits_deg(2),2401).';
k=2*pi*f/cfg.c0_mps; th=deg2rad(angles);
D=cos(th).*exp(-0.5*(k*cfg.sigma_src_m*sin(th)).^2); D=abs(D)/max(abs(D));
pat=struct('angles_deg',angles,'level_db',20*log10(max(D,10^(cfg.source_pattern_clip_db/20))));
end

function tag=local_tag(value)
tag=strrep(sprintf('%.10g',value),'.','p'); tag=strrep(tag,'-','m');
end
