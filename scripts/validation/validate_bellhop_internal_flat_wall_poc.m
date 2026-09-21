function validation = validate_bellhop_internal_flat_wall_poc(options)
%VALIDATE_BELLHOP_INTERNAL_FLAT_WALL_POC Validate r=100 m native reflection.
arguments
    options.output_dir (1,:) char = ''
    options.validation_exe (1,:) char = ''
    options.official_exe (1,:) char = ''
    options.frequency_hz (1,1) double {mustBePositive} = 4000
    options.beam_counts (1,:) double {mustBeInteger,mustBePositive} = [2001 5001 10001]
    options.step_values_m (1,:) double {mustBePositive} = [0.2 0.1 0.05]
    options.fail_on_check (1,1) logical = true
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup_vertical_project;
if isempty(options.official_exe), options.official_exe=resolve_bellhop_exe_vertical(); end
if isempty(options.output_dir)
    options.output_dir=fullfile(root,'results','validation','bellhop_internal_flat_wall_poc');
end
if isempty(options.validation_exe)
    options.validation_exe=fullfile(options.output_dir,'bin','bellhop_iwall_flat_2020.exe');
end
assert(isfile(options.validation_exe),'Missing validation binary: %s',options.validation_exe);
assert(isfile(options.official_exe),'Missing official Bellhop 2020 binary: %s',options.official_exe);

cfg=struct('c0_mps',1500,'wall_range_m',100,'mapped_receiver_range_m',103, ...
    'source_depth_m',0,'receiver_depths_m',0,'receiver_ranges_m',[102 103], ...
    'run_type','C','angle_limits_deg',[-30 30],'domain_half_depth_m',1000, ...
    'sigma_src_m',0.3,'source_pattern_clip_db',-120,'frequency_hz',options.frequency_hz);
pat=local_pattern(cfg.frequency_hz,cfg);
cfg.source_pattern_angles_deg=pat.angles_deg;
cfg.source_pattern_level_db=pat.level_db;

template=struct('beam_count',NaN,'step_m',NaN,'wall_pressure',complex(NaN), ...
    'reference_pressure',complex(NaN),'expected_pressure',complex(NaN), ...
    'tl_error_db',NaN,'phase_error_rad',NaN,'complex_relative_error',NaN, ...
    'intersection_residual_m',NaN,'specular_direction_error',NaN, ...
    'rotation_direction_error',NaN,'phase_jump_error_rad',NaN,'amp_jump_error',NaN, ...
    'p_reflect_error',NaN,'q_reflect_error',NaN,'p_rotation_error',NaN,'q_rotation_error',NaN, ...
    'kappa_abs_max',NaN,'travel_time_s',NaN,'travel_time_error_s',NaN, ...
    'travel_time_imag_abs_s',NaN,'min_post_range_increment_m',NaN,'all_post_range_positive',false, ...
    'wall_case_root','','reference_case_root','');
rows=repmat(template,numel(options.step_values_m)*numel(options.beam_counts),1);
wall_runs=cell(size(rows)); reference_runs=cell(size(rows)); rr=0;
for ss=1:numel(options.step_values_m)
    for bb=1:numel(options.beam_counts)
        rr=rr+1; step=options.step_values_m(ss); beams=options.beam_counts(bb);
        tag=sprintf('step_%s_beams_%d',local_tag(step),beams);
        common=cfg; common.step_m=step; common.beam_count=beams;
        common.case_root=fullfile(options.output_dir,'cases','reference',[tag '_reference']);
        common.bellhop_exe=options.official_exe;
        reference_runs{rr}=run_bellhop_unfolded_gaussian_vertical(common);
        wall_cfg=common; wall_cfg.case_root=fullfile(options.output_dir,'cases','wall',[tag '_wall']);
        wall_cfg.bellhop_exe=options.validation_exe;
        wall_runs{rr}=run_bellhop_internal_flat_wall_poc_vertical(wall_cfg);

        pref=select_bellhop_shd_pressure_at_range_vertical( ...
            reference_runs{rr}.data,cfg.mapped_receiver_range_m,cfg.receiver_depths_m(1));
        pwall=select_bellhop_shd_pressure_at_range_vertical( ...
            wall_runs{rr}.data,cfg.mapped_receiver_range_m,cfg.receiver_depths_m(1));
        expected=-pref;
        d=wall_runs{rr}.diagnostics;
        [~,i0]=min(abs(d.alpha_deg));
        spec_err=max(hypot(d.ref_ur+d.inc_ur,d.ref_uz-d.inc_uz));
        rot_err=max(hypot(d.rot_ur+d.ref_ur,d.rot_uz+d.ref_uz));
        q=template; q.beam_count=beams; q.step_m=step; q.wall_pressure=pwall;
        q.reference_pressure=pref; q.expected_pressure=expected;
        q.tl_error_db=20*log10(abs(pwall)/abs(expected));
        q.phase_error_rad=angle(pwall*conj(expected));
        q.complex_relative_error=abs(pwall-expected)/abs(expected);
        q.intersection_residual_m=max(abs(d.residual));
        q.specular_direction_error=spec_err; q.rotation_direction_error=rot_err;
        q.phase_jump_error_rad=max(abs(d.phase_delta-pi));
        q.amp_jump_error=max(abs(d.amp_delta));
        q.p_reflect_error=max(d.p_ref_error); q.q_reflect_error=max(d.q_ref_error);
        q.p_rotation_error=max(d.p_rot_error); q.q_rotation_error=max(d.q_rot_error);
        q.kappa_abs_max=max(abs(d.kappa)); q.travel_time_s=d.tau_receiver_real(i0);
        q.travel_time_error_s=q.travel_time_s-103/1500;
        q.travel_time_imag_abs_s=abs(d.tau_receiver_imag(i0));
        q.min_post_range_increment_m=min(d.min_post_dr);
        q.all_post_range_positive=all(d.min_post_dr>0);
        q.wall_case_root=wall_cfg.case_root; q.reference_case_root=common.case_root;
        rows(rr)=q;
    end
end
table_rows=struct2table(rows);

checks=struct;
checks.intersection=max(abs(table_rows.intersection_residual_m))<=1e-9;
checks.direction=max(table_rows.specular_direction_error)<=1e-12 && max(table_rows.rotation_direction_error)<=1e-12;
checks.phase_jump=max(table_rows.phase_jump_error_rad)<=1e-12;
checks.amplitude=max(table_rows.amp_jump_error)<=1e-12;
checks.flat_curvature=max(table_rows.kappa_abs_max)==0;
checks.dynamic_state=max([table_rows.p_reflect_error;table_rows.q_reflect_error; ...
    table_rows.p_rotation_error;table_rows.q_rotation_error])<=1e-12;
checks.travel_time=max(abs(table_rows.travel_time_error_s))<=1e-10 && ...
    max(table_rows.travel_time_imag_abs_s)<=1e-12;
checks.positive_range=all(table_rows.all_post_range_positive);
checks.field=max(abs(table_rows.tl_error_db))<=0.05 && ...
    max(abs(table_rows.phase_error_rad))<=0.01 && max(table_rows.complex_relative_error)<=0.01;
checks.all=all(structfun(@(x)logical(x),checks));

if ~exist(options.output_dir,'dir'), mkdir(options.output_dir); end
writetable(table_rows,fullfile(options.output_dir,'flat_wall_convergence.csv'));
diag_summary=table_rows(:,{'beam_count','step_m','intersection_residual_m','specular_direction_error', ...
    'rotation_direction_error','phase_jump_error_rad','amp_jump_error','p_reflect_error','q_reflect_error', ...
    'p_rotation_error','q_rotation_error','travel_time_error_s','min_post_range_increment_m', ...
    'tl_error_db','phase_error_rad','complex_relative_error'});
writetable(diag_summary,fullfile(options.output_dir,'flat_wall_diagnostics_summary.csv'));
validation=struct('config',cfg,'options',options,'rows',rows,'table',table_rows, ...
    'checks',checks,'wall_runs',{wall_runs},'reference_runs',{reference_runs}, ...
    'files',struct('convergence_csv',fullfile(options.output_dir,'flat_wall_convergence.csv'), ...
    'diagnostics_csv',fullfile(options.output_dir,'flat_wall_diagnostics_summary.csv'), ...
    'mat',fullfile(options.output_dir,'flat_wall_poc_results.mat')));
save(validation.files.mat,'validation','-v7.3');
disp(table_rows(:,{'beam_count','step_m','tl_error_db','phase_error_rad','complex_relative_error', ...
    'intersection_residual_m','travel_time_error_s','min_post_range_increment_m'}));
disp(checks);
if options.fail_on_check && ~checks.all
    error('Bellhop internal flat-wall POC failed; inspect %s.',options.output_dir);
end
end

function pat=local_pattern(f,cfg)
angles=linspace(cfg.angle_limits_deg(1),cfg.angle_limits_deg(2),2401).';
k=2*pi*f/cfg.c0_mps; th=deg2rad(angles);
D=cos(th).*exp(-0.5*(k*cfg.sigma_src_m*sin(th)).^2);
D=abs(D)/max(abs(D));
pat=struct('angles_deg',angles, ...
    'level_db',20*log10(max(D,10^(cfg.source_pattern_clip_db/20))));
end

function tag=local_tag(value)
tag=strrep(sprintf('%.10g',value),'.','p');
tag=strrep(tag,'-','m');
end
