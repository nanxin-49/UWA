function validation = validate_bellhop_internal_tilted_wall_poc(options)
%VALIDATE_BELLHOP_INTERNAL_TILTED_WALL_POC Validate a straight internal wall.
% The validation binary traces to r=R0+a*z, invokes Bellhop's native
% Reflect2D, rotates the post-wall branch by pi, and leaves the influence
% routine untouched.  The companion native case uses the same line as ATI.
arguments
    options.output_dir (1,:) char = ''
    options.validation_exe (1,:) char = ''
    options.official_exe (1,:) char = ''
    options.frequency_hz (1,1) double {mustBePositive} = 4000
    options.wall_slope (1,1) double = 0.005
    options.beam_counts (1,:) double {mustBeInteger,mustBePositive} = [2001 5001 10001]
    options.step_values_m (1,:) double {mustBePositive} = [0.2 0.1 0.05]
    options.fail_on_check (1,1) logical = true
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup_vertical_project;
if isempty(options.official_exe), options.official_exe=getenv('BELLHOP_EXE'); end
if isempty(options.official_exe)
    error('Set options.official_exe or BELLHOP_EXE to the AcousticsToolbox 2020 bellhop.exe.');
end
if isempty(options.output_dir)
    options.output_dir=fullfile(root,'results','validation','bellhop_internal_tilted_wall_poc');
end
if isempty(options.validation_exe)
    options.validation_exe=fullfile(options.output_dir,'bin','bellhop_iwall_tilted_2020.exe');
end
assert(isfile(options.validation_exe),'Missing tilted-wall validation binary: %s',options.validation_exe);
assert(isfile(options.official_exe),'Missing official Bellhop 2020 binary: %s',options.official_exe);
assert(abs(options.wall_slope-0.005)<=1e-12,'This POC fixes slope a=0.005.');
if ~exist(options.output_dir,'dir'), mkdir(options.output_dir); end

cfg=struct('c0_mps',1500,'wall_r0_m',100,'wall_slope',options.wall_slope, ...
    'mapped_receiver_range_m',103,'native_receiver_range_m',97, ...
    'source_depth_m',0,'receiver_depths_m',0,'receiver_ranges_m',[102 103], ...
    'run_type','C','angle_limits_deg',[-30 30],'domain_half_depth_m',1000, ...
    'sigma_src_m',0.3,'source_pattern_clip_db',-120,'frequency_hz',options.frequency_hz);
pat=local_pattern(cfg.frequency_hz,cfg);
cfg.source_pattern_angles_deg=pat.angles_deg;
cfg.source_pattern_level_db=pat.level_db;
[target_alpha_deg,~,target_path_m]=local_target_ray(cfg);

template=struct('beam_count',NaN,'step_m',NaN,'wall_pressure',complex(NaN), ...
    'native_total_pressure',complex(NaN),'native_direct_pressure',complex(NaN), ...
    'native_reflected_pressure',complex(NaN),'field_tl_error_db',NaN, ...
    'field_phase_error_rad',NaN,'field_complex_relative_error',NaN, ...
    'intersection_residual_m',NaN,'intersection_point_error_m',NaN, ...
    'tangent_error',NaN,'normal_error',NaN,'specular_direction_error',NaN, ...
    'rotation_direction_error',NaN,'phase_jump_error_rad',NaN,'amp_jump_error',NaN, ...
    'p_reflect_error',NaN,'q_reflect_error',NaN,'p_rotation_error',NaN,'q_rotation_error',NaN, ...
    'kappa_abs_max',NaN,'path_length_error_m',NaN,'travel_time_error_s',NaN, ...
    'target_travel_time_error_s',NaN, ...
    'travel_time_imag_abs_s',NaN,'min_post_range_increment_m',NaN,'all_post_range_positive',false, ...
    'target_ray_angle_deg',NaN,'target_ray_angle_error_deg',NaN,'target_path_length_m',NaN, ...
    'native_case_root','','wall_case_root','');
rows=repmat(template,numel(options.step_values_m)*numel(options.beam_counts),1);
wall_runs=cell(size(rows)); native_runs=cell(size(rows)); direct_runs=cell(size(rows));
rr=0; native_arrival_case=[];
for ss=1:numel(options.step_values_m)
    for bb=1:numel(options.beam_counts)
        rr=rr+1; step=options.step_values_m(ss); beams=options.beam_counts(bb);
        tag=sprintf('step_%s_beams_%d',local_tag(step),beams);
        common=cfg; common.step_m=step; common.beam_count=beams;
        common.case_root=fullfile(options.output_dir,'cases','wall',[tag '_wall']);
        common.bellhop_exe=options.validation_exe;
        wall_runs{rr}=run_bellhop_internal_tilted_wall_poc_vertical(common);

        native=common; native.case_root=fullfile(options.output_dir,'cases','native',[tag '_native']);
        native.bellhop_exe=options.official_exe;
        native.receiver_ranges_m=[cfg.native_receiver_range_m 98];
        native_runs{rr}=run_bellhop_native_tilted_wall_vertical(native);

        direct=native; direct.case_root=fullfile(options.output_dir,'cases','direct',[tag '_direct']);
        direct_runs{rr}=run_bellhop_unfolded_gaussian_vertical(direct);

        pwall=select_bellhop_shd_pressure_at_range_vertical( ...
            wall_runs{rr}.data,cfg.mapped_receiver_range_m,cfg.receiver_depths_m(1));
        pnative=select_bellhop_shd_pressure_at_range_vertical( ...
            native_runs{rr}.data,cfg.native_receiver_range_m,cfg.receiver_depths_m(1));
        pdirect=select_bellhop_shd_pressure_at_range_vertical( ...
            direct_runs{rr}.data,cfg.native_receiver_range_m,cfg.receiver_depths_m(1));
        pref=pnative-pdirect;
        d=wall_runs{rr}.diagnostics;
        alpha=deg2rad(d.alpha_deg);
        ux=cos(alpha); uz=sin(alpha);
        h=(cfg.wall_r0_m-cfg.wall_slope*0)./(ux-cfg.wall_slope*uz);
        hit_r=h.*ux; hit_z=h.*uz;
        % tau is sampled where the transformed ray crosses the mapped range
        % plane. Only the target fan ray reaches the physical receiver point.
        ref_ur=ux-2*(ux-cfg.wall_slope*uz)/(1+cfg.wall_slope^2);
        rot_ur=-ref_ur;
        hit_r_rot=2*cfg.wall_r0_m-hit_r;
        expected_len=h+(cfg.mapped_receiver_range_m-hit_r_rot)./rot_ur;
        path_len=c0_from_cfg(cfg).*d.tau_receiver_real;

        q=template; q.beam_count=beams; q.step_m=step; q.wall_pressure=pwall;
        q.native_total_pressure=pnative; q.native_direct_pressure=pdirect;
        q.native_reflected_pressure=pref;
        q.field_tl_error_db=20*log10(abs(pwall)/max(abs(pref),realmin));
        q.field_phase_error_rad=angle(pwall*conj(pref));
        q.field_complex_relative_error=abs(pwall-pref)/max(abs(pref),realmin);
        q.intersection_residual_m=max(abs(d.line_residual));
        q.intersection_point_error_m=max(hypot(d.hit_r-hit_r,d.hit_z-hit_z));
        q.tangent_error=max(d.tangent_error); q.normal_error=max(d.normal_error);
        q.specular_direction_error=max(d.specular_error);
        q.rotation_direction_error=max(d.rotation_error);
        q.phase_jump_error_rad=max(abs(angle(exp(1i*(d.phase_delta-pi)))));
        q.amp_jump_error=max(abs(d.amp_delta));
        q.p_reflect_error=max(d.p_ref_error); q.q_reflect_error=max(d.q_ref_error);
        q.p_rotation_error=max(d.p_rot_error); q.q_rotation_error=max(d.q_rot_error);
        q.kappa_abs_max=max(abs(d.kappa));
        q.path_length_error_m=max(abs(path_len-expected_len));
        q.travel_time_error_s=q.path_length_error_m/cfg.c0_mps;
        [~,i0]=min(abs(d.alpha_deg-target_alpha_deg));
        q.target_travel_time_error_s=d.tau_receiver_real(i0)-target_path_m/cfg.c0_mps;
        q.travel_time_imag_abs_s=max(abs(d.tau_receiver_imag));
        q.min_post_range_increment_m=min(d.min_post_dr);
        q.all_post_range_positive=all(d.min_post_dr>0);
        q.target_ray_angle_deg=target_alpha_deg;
        q.target_ray_angle_error_deg=d.alpha_deg(i0)-target_alpha_deg;
        q.target_path_length_m=target_path_m;
        q.native_case_root=native.case_root; q.wall_case_root=common.case_root;
        rows(rr)=q;

        % One arrivals run at the finest setting supplies an independent
        % native path/time/amplitude/phase audit without changing the C field.
        if ss==numel(options.step_values_m) && bb==numel(options.beam_counts)
            ar=native; ar.run_type='A'; ar.case_root=fullfile(options.output_dir,'cases','native', ...
                [tag '_native_arrivals']);
            native_arrival_case=run_bellhop_native_tilted_wall_vertical(ar);
        end
    end
end
table_rows=struct2table(rows);

checks=struct;
checks.intersection=max(table_rows.intersection_residual_m)<=1e-9 && max(table_rows.intersection_point_error_m)<=1e-8;
checks.frame=max([table_rows.tangent_error;table_rows.normal_error])<=1e-12;
checks.direction=max(table_rows.specular_direction_error)<=1e-12 && max(table_rows.rotation_direction_error)<=1e-12;
checks.phase_jump=max(table_rows.phase_jump_error_rad)<=1e-12;
checks.amplitude=max(table_rows.amp_jump_error)<=1e-12;
checks.flat_curvature=max(table_rows.kappa_abs_max)==0;
checks.dynamic_state=max([table_rows.p_reflect_error;table_rows.q_reflect_error; ...
    table_rows.p_rotation_error;table_rows.q_rotation_error])<=1e-12;
% The target-ray crossing is sampled from a finite beam fan.  Keep its small
% interpolation error as a diagnostic, but do not let fan discretization mask
% the exact segment/path-time regression gate.
checks.path_time=max(table_rows.path_length_error_m)<=1e-8 && ...
    max(abs(table_rows.travel_time_error_s))<=1e-10 && ...
    max(abs(table_rows.target_travel_time_error_s))<=2e-7 && ...
    max(table_rows.travel_time_imag_abs_s)<=1e-12;
checks.positive_range=all(table_rows.all_post_range_positive);
% Native Bellhop evaluates a backward physical branch in its original range
% chart. Its coherent Cartesian amplitude therefore has a repeatable small
% range-spreading offset relative to the proper-rotation chart; geometry and
% phase remain strict gates.
checks.native_field_phase=max(abs(table_rows.field_phase_error_rad))<=0.02;
checks.native_field_amplitude=max(abs(table_rows.field_tl_error_db))<=0.35 && ...
    max(table_rows.field_complex_relative_error)<=0.04;
checks.native_field=checks.native_field_phase && checks.native_field_amplitude;
checks.all=all(structfun(@(x)logical(x),checks));

writetable(table_rows,fullfile(options.output_dir,'tilted_wall_convergence.csv'));
diag_summary=table_rows(:,{'beam_count','step_m','intersection_residual_m','intersection_point_error_m', ...
    'tangent_error','normal_error','specular_direction_error','rotation_direction_error', ...
    'phase_jump_error_rad','amp_jump_error','p_reflect_error','q_reflect_error','p_rotation_error', ...
    'q_rotation_error','path_length_error_m','travel_time_error_s','target_travel_time_error_s', ...
    'target_ray_angle_error_deg', ...
    'min_post_range_increment_m', ...
    'field_tl_error_db','field_phase_error_rad','field_complex_relative_error'});
writetable(diag_summary,fullfile(options.output_dir,'tilted_wall_diagnostics_summary.csv'));
validation=struct('config',cfg,'options',options,'rows',rows,'table',table_rows, ...
    'checks',checks,'wall_runs',{wall_runs},'native_runs',{native_runs}, ...
    'direct_runs',{direct_runs},'native_arrival_case',native_arrival_case, ...
    'files',struct('convergence_csv',fullfile(options.output_dir,'tilted_wall_convergence.csv'), ...
    'diagnostics_csv',fullfile(options.output_dir,'tilted_wall_diagnostics_summary.csv'), ...
    'mat',fullfile(options.output_dir,'tilted_wall_poc_results.mat')));
save(validation.files.mat,'validation','-v7.3');
disp(table_rows(:,{'beam_count','step_m','field_tl_error_db','field_phase_error_rad', ...
    'field_complex_relative_error','intersection_residual_m','intersection_point_error_m', ...
    'path_length_error_m','travel_time_error_s','target_travel_time_error_s', ...
    'min_post_range_increment_m'}));
disp(checks);
if options.fail_on_check && ~checks.all
    error('Bellhop internal tilted-wall POC failed; inspect %s.',options.output_dir);
end
end

function [alpha_deg,hit,path_len]=local_target_ray(cfg)
% Mirror-source construction gives the exact specular ray from (0,0) to the
% physical receiver (97,0) for r=R0+a*z.
n=[1;-cfg.wall_slope]/sqrt(1+cfg.wall_slope^2);
wall0=[cfg.wall_r0_m;0]; src=[0;0]; rx=[cfg.native_receiver_range_m;0];
image=src-2*dot(src-wall0,n)*n;
v=rx-image; lam=-dot(n,image-wall0)/dot(n,v);
hit=image+lam*v;
alpha_deg=rad2deg(atan2(hit(2),hit(1)));
path_len=norm(hit-src)+norm(rx-hit);
end

function pat=local_pattern(f,cfg)
angles=linspace(cfg.angle_limits_deg(1),cfg.angle_limits_deg(2),2401).';
k=2*pi*f/cfg.c0_mps; th=deg2rad(angles);
D=cos(th).*exp(-0.5*(k*cfg.sigma_src_m*sin(th)).^2);
D=abs(D)/max(abs(D));
pat=struct('angles_deg',angles,'level_db',20*log10(max(D,10^(cfg.source_pattern_clip_db/20))));
end

function tag=local_tag(value)
tag=strrep(sprintf('%.10g',value),'.','p'); tag=strrep(tag,'-','m');
end

function c=c0_from_cfg(cfg)
c=cfg.c0_mps;
end
