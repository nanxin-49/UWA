function validation = validate_bellhop_internal_pm_wall(options)
%VALIDATE_BELLHOP_INTERNAL_PM_WALL Fixed-seed 1-D PM wall validation.
% The same sampled eta(s) is written to native C-ATI and to the rotated
% validation-only internal wall.  No PE or rough-wall scattering formula is
% used; Bellhop's native Reflect2D and InfluenceGeoHatCart remain authoritative.
arguments
    options.output_dir (1,:) char = ''
    options.validation_exe (1,:) char = ''
    options.official_exe (1,:) char = 'E:\MISC\BELLHOP\AcousticsToolbox_2020\windows-bin-20201102\bellhop.exe'
    options.frequency_hz (1,1) double {mustBePositive} = 4000
    options.wind_speed_mps (1,1) double {mustBePositive} = 6
    options.profile_span_m (1,1) double {mustBePositive} = 160
    options.profile_origin_m (1,1) double = NaN
    options.profile_counts (1,:) double {mustBeInteger,mustBePositive} = [129 257]
    options.beam_counts (1,:) double {mustBeInteger,mustBePositive} = [2001 5001]
    options.step_values_m (1,:) double {mustBePositive} = [0.2 0.1]
    options.seed (1,1) double {mustBeFinite} = 260001
    options.fail_on_check (1,1) logical = true
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup_vertical_project;
if isempty(options.output_dir)
    options.output_dir=fullfile(root,'results','validation','bellhop_internal_pm_wall_poc');
end
if isempty(options.validation_exe)
    options.validation_exe=fullfile(options.output_dir,'bin','bellhop_iwall_pm_2020.exe');
end
assert(isfile(options.validation_exe),'Missing PM validation binary: %s',options.validation_exe);
assert(isfile(options.official_exe),'Missing official Bellhop 2020 binary: %s',options.official_exe);
if ~exist(options.output_dir,'dir'), mkdir(options.output_dir); end

R0=100; zc=0;
cfg=struct('c0_mps',1500,'wall_r0_m',R0,'mapped_receiver_range_m',103, ...
    'native_receiver_range_m',97,'source_depth_m',zc,'receiver_depths_m',zc, ...
    'receiver_ranges_m',[102 103],'run_type','C','angle_limits_deg',[-30 30], ...
    'domain_half_depth_m',1000,'sigma_src_m',0.3,'source_pattern_clip_db',-120, ...
    'frequency_hz',options.frequency_hz);
pat=local_pattern(cfg.frequency_hz,cfg);
cfg.source_pattern_angles_deg=pat.angles_deg;
cfg.source_pattern_level_db=pat.level_db;

template=struct('seed',NaN,'wind_speed_mps',NaN,'profile_count',NaN,'beam_count',NaN, ...
    'step_m',NaN,'rms_height_m',NaN,'rms_slope',NaN,'max_slope',NaN, ...
    'rms_curvature_per_m',NaN,'max_curvature_per_m',NaN,'grazing_hit_count',NaN, ...
    'near_vertex_hit_count',NaN,'rejected_nonmonotone_count',NaN, ...
    'successful_wall_reflection_fraction',NaN,'wall_pressure',complex(NaN), ...
    'native_source_clearance_m',NaN,'native_field_nonzero',false, ...
    'profile_nonmonotone_segment_count',NaN,'profile_min_abs_dr_m',NaN, ...
    'native_total_pressure',complex(NaN),'native_direct_pressure',complex(NaN), ...
    'native_reflected_pressure',complex(NaN),'field_tl_error_db',NaN, ...
    'field_phase_error_rad',NaN,'field_complex_relative_error',NaN, ...
    'intersection_residual_m',NaN,'tangent_error',NaN,'normal_error',NaN, ...
    'curvature_error_per_m',NaN,'specular_direction_error',NaN, ...
    'rotation_direction_error',NaN,'phase_jump_error_rad',NaN,'amp_jump_error',NaN, ...
    'p_reflect_error',NaN,'q_reflect_error',NaN,'p_rotation_error',NaN, ...
    'q_rotation_error',NaN,'path_length_error_m',NaN,'travel_time_error_s',NaN, ...
    'travel_time_imag_abs_s',NaN,'min_post_range_increment_m',NaN, ...
    'all_post_range_positive',false,'wall_case_root','','native_case_root','');
nRows=numel(options.profile_counts)*numel(options.step_values_m)*numel(options.beam_counts);
rows=repmat(template,nRows,1); wall_runs=cell(nRows,1); native_runs=cell(nRows,1); direct_runs=cell(nRows,1);
profile_meta=cell(numel(options.profile_counts),1); rr=0;
for pp=1:numel(options.profile_counts)
    N=options.profile_counts(pp);
    [s,eta,pm]=sample_fixed_pm_profile_vertical(wind_speed_mps=options.wind_speed_mps, ...
        span_m=options.profile_span_m,sample_count=N,seed=options.seed,origin_m=options.profile_origin_m);
    profile_meta{pp}=pm;
    rRot=R0-eta; zRot=s; rNat=R0+s; zNat=eta;
    for ss=1:numel(options.step_values_m)
        for bb=1:numel(options.beam_counts)
            rr=rr+1; step=options.step_values_m(ss); beams=options.beam_counts(bb);
            tag=sprintf('seed_%d_N_%d_step_%s_beams_%d',pm.seed,N,local_tag(step),beams);
            common=cfg; common.wall_seed=pm.seed; common.wall_profile_r_m=rRot; ...
                common.wall_profile_z_m=zRot; common.step_m=step; common.beam_count=beams;
            common.case_root=fullfile(options.output_dir,'cases','wall',[tag '_wall']);
            common.bellhop_exe=options.validation_exe;
            wall_runs{rr}=run_bellhop_internal_pm_wall_poc_vertical(common);
            native=common; native.case_root=fullfile(options.output_dir,'cases','native',[tag '_native']);
            native.bellhop_exe=options.official_exe; native.receiver_depths_m=zc;
            native.receiver_ranges_m=[97 98]; native.wall_profile_r_m=rNat; native.wall_profile_z_m=zNat;
            native_runs{rr}=run_bellhop_native_pm_wall_vertical(native);
            direct=native; direct.case_root=fullfile(options.output_dir,'cases','direct',[tag '_direct']);
            direct_runs{rr}=run_bellhop_unfolded_gaussian_vertical(direct);

            pwall=select_bellhop_shd_pressure_at_range_vertical(wall_runs{rr}.data,103,zc);
            pnative=select_bellhop_shd_pressure_at_range_vertical(native_runs{rr}.data,97,zc);
            pdirect=select_bellhop_shd_pressure_at_range_vertical(direct_runs{rr}.data,97,zc);
            pref=pnative-pdirect; d=wall_runs{rr}.diagnostics;
            geom=local_geometry_audit(d,s,eta,R0);
            hitR=2*R0-d.hit_r; pre=sqrt(d.hit_r.^2+d.hit_z.^2);
            expected=pre+(103-hitR)./d.rot_ur; actual=cfg.c0_mps*d.tau_receiver_real;
            q=template; q.seed=pm.seed; q.wind_speed_mps=pm.wind_speed_mps; q.profile_count=N;
            q.beam_count=beams; q.step_m=step; q.rms_height_m=pm.rms_height_m;
            q.rms_slope=pm.rms_slope; q.max_slope=pm.max_slope;
            q.rms_curvature_per_m=pm.rms_curvature_per_m; q.max_curvature_per_m=pm.max_curvature_per_m;
            q.native_source_clearance_m=cfg.source_depth_m-eta(1);
            drot=diff(rRot);
            q.profile_nonmonotone_segment_count=sum(drot(1:end-1).*drot(2:end)<0);
            q.profile_min_abs_dr_m=min(abs(drot));
            q.grazing_hit_count=sum(abs(geom.incident_normal_component)<1e-3);
            q.near_vertex_hit_count=sum(geom.near_vertex);
            q.rejected_nonmonotone_count=sum(d.min_post_dr<=0);
            q.successful_wall_reflection_fraction=height(d)/beams;
            q.wall_pressure=pwall; q.native_total_pressure=pnative; q.native_direct_pressure=pdirect;
            q.native_reflected_pressure=pref; q.field_tl_error_db=20*log10(abs(pwall)/max(abs(pref),realmin));
            q.field_phase_error_rad=angle(pwall*conj(pref));
            q.field_complex_relative_error=abs(pwall-pref)/max(abs(pref),realmin);
            q.native_field_nonzero=abs(pref)>1e-6*max(abs(pnative),realmin);
            q.intersection_residual_m=max(abs(d.wall_residual)); q.tangent_error=max(geom.tangent_error);
            q.normal_error=max(geom.normal_error); q.curvature_error_per_m=max(abs(geom.curvature_error));
            q.specular_direction_error=max(d.specular_error); q.rotation_direction_error=max(d.rotation_error);
            q.phase_jump_error_rad=max(abs(angle(exp(1i*(d.phase_delta-pi))))); q.amp_jump_error=max(abs(d.amp_delta));
            q.p_reflect_error=max(d.p_ref_error); q.q_reflect_error=max(d.q_ref_error);
            q.p_rotation_error=max(d.p_rot_error); q.q_rotation_error=max(d.q_rot_error);
            q.path_length_error_m=max(abs(actual-expected)); q.travel_time_error_s=q.path_length_error_m/cfg.c0_mps;
            q.travel_time_imag_abs_s=max(abs(d.tau_receiver_imag)); q.min_post_range_increment_m=min(d.min_post_dr);
            q.all_post_range_positive=all(d.min_post_dr>0); q.wall_case_root=common.case_root; q.native_case_root=native.case_root;
            rows(rr)=q;
        end
    end
end

table_rows=struct2table(rows);
checks=struct;
checks.intersection=max(table_rows.intersection_residual_m)<=1e-8;
checks.geometry=max([table_rows.tangent_error;table_rows.normal_error])<=5e-5;
checks.curvature=max(table_rows.curvature_error_per_m)<=5e-5;
checks.direction=max(table_rows.specular_direction_error)<=1e-10 && max(table_rows.rotation_direction_error)<=1e-10;
checks.phase_jump=max(table_rows.phase_jump_error_rad)<=1e-10;
checks.amplitude=max(table_rows.amp_jump_error)<=1e-10;
checks.dynamic_state=max([table_rows.q_reflect_error;table_rows.p_rotation_error;table_rows.q_rotation_error])<=1e-10;
% Curvature is intentionally not required to vanish for PM walls.  A finite
% p kick is the expected native Reflect2D signature and is retained as a
% diagnostic hard gate against silently dropping kappa.
checks.curvature_kick=~any(isnan(table_rows.p_reflect_error) | isinf(table_rows.p_reflect_error)) && ...
    max(table_rows.p_reflect_error)>1e-12;
checks.wall_pathologies=all(table_rows.grazing_hit_count==0) && ...
    all(table_rows.near_vertex_hit_count==0) && all(table_rows.rejected_nonmonotone_count==0);
profile_groups=unique(table_rows.profile_count);
kick_by_profile=arrayfun(@(p) max(table_rows.p_reflect_error(table_rows.profile_count==p)),profile_groups);
checks.curvature_kick_sampling=numel(kick_by_profile)<2 || ...
    max(kick_by_profile)/max(min(kick_by_profile),realmin)<=10;
checks.path_time=max(table_rows.path_length_error_m)<=1e-7 && max(table_rows.travel_time_imag_abs_s)<=1e-12;
checks.positive_range=all(table_rows.all_post_range_positive) && all(table_rows.rejected_nonmonotone_count==0);
checks.reflection_fraction=all(table_rows.successful_wall_reflection_fraction>=0.999);
checks.native_field_phase=max(abs(table_rows.field_phase_error_rad))<=0.10;
checks.native_field_nonzero=all(table_rows.native_field_nonzero);
checks.all=checks.intersection && checks.geometry && checks.curvature && checks.direction && checks.phase_jump && ...
    checks.amplitude && checks.dynamic_state && checks.curvature_kick && checks.wall_pathologies && ...
    checks.curvature_kick_sampling && checks.path_time && checks.positive_range && ...
    checks.reflection_fraction && checks.native_field_phase && checks.native_field_nonzero;

csv_file=fullfile(options.output_dir,'pm_wall_convergence.csv'); writetable(table_rows,csv_file);
diag_summary=table_rows(:,{'seed','wind_speed_mps','profile_count','beam_count','step_m','rms_height_m', ...
    'rms_slope','max_slope','rms_curvature_per_m','max_curvature_per_m','grazing_hit_count', ...
    'near_vertex_hit_count','rejected_nonmonotone_count','successful_wall_reflection_fraction', ...
    'native_source_clearance_m','native_field_nonzero', ...
    'profile_nonmonotone_segment_count','profile_min_abs_dr_m', ...
    'intersection_residual_m','tangent_error','normal_error','curvature_error_per_m', ...
    'phase_jump_error_rad','amp_jump_error','p_reflect_error','q_reflect_error','path_length_error_m', ...
    'travel_time_error_s','min_post_range_increment_m','field_tl_error_db','field_phase_error_rad', ...
    'field_complex_relative_error'});
diag_file=fullfile(options.output_dir,'pm_wall_diagnostics_summary.csv'); writetable(diag_summary,diag_file);
profile_json=fullfile(options.output_dir,'pm_profile_meta.json');
fid=fopen(profile_json,'w'); if fid<0, error('Cannot create %s.',profile_json); end
cleanup=onCleanup(@()fclose(fid)); fprintf(fid,'%s',jsonencode(profile_meta)); clear cleanup
validation=struct('config',cfg,'options',options,'table',table_rows,'checks',checks, ...
    'profile_meta',{profile_meta},'files',struct('convergence_csv',csv_file,'diagnostics_csv',diag_file,'profile_json',profile_json));
save(fullfile(options.output_dir,'pm_wall_validation_results.mat'),'validation','-v7.3');
disp(table_rows(:,{'profile_count','beam_count','step_m','rms_height_m','max_slope','max_curvature_per_m', ...
    'profile_nonmonotone_segment_count','profile_min_abs_dr_m','native_source_clearance_m', ...
    'field_tl_error_db','field_phase_error_rad','intersection_residual_m','curvature_error_per_m'}));
disp(checks);
if options.fail_on_check && ~checks.all
    error('Bellhop fixed-seed PM wall validation failed; inspect %s.',options.output_dir);
end
end

function audit=local_geometry_audit(d,s,eta,R0)
Q=[0 -1;1 0]; n=height(d); audit=struct('tangent_error',zeros(n,1), ...
    'normal_error',zeros(n,1),'curvature_error',zeros(n,1), ...
    'incident_normal_component',zeros(n,1),'near_vertex',false(n,1));
for ii=1:n
    sm=d.hit_z(ii); ev=interp1(s,eta,sm,'linear','extrap');
    nat=local_polyline_geometry([R0+s(:) eta(:)].',[R0+sm ev]);
    rot=local_polyline_geometry([R0-eta(:) s(:)].',[d.hit_r(ii) d.hit_z(ii)]);
    audit.tangent_error(ii)=norm([d.wall_t_r(ii);d.wall_t_z(ii)]-Q*nat.t);
    audit.normal_error(ii)=norm([d.wall_n_r(ii);d.wall_n_z(ii)]-Q*nat.n);
    audit.curvature_error(ii)=d.kappa(ii)-nat.kappa;
    audit.incident_normal_component(ii)=dot([d.inc_ur(ii);d.inc_uz(ii)],[d.wall_n_r(ii);d.wall_n_z(ii)]);
    audit.near_vertex(ii)=min(abs(sm-s))<=1e-8;
    if norm(rot.t-Q*nat.t)>1e-8 || norm(rot.n-Q*nat.n)>1e-8
        error('Shared PM profile geometry lost proper-rotation covariance.');
    end
end
end

function geom=local_polyline_geometry(points,hit)
x=points(:,1); z=points(:,2); n=numel(x); seg=diff(points,1,1); len=vecnorm(seg,2,2); seg=seg./len;
node=zeros(n,2); node(1,:)=[1 0]; node(n,:)=[1 0]; node(2:n-1,:)=0.5*(seg(1:n-2,:)+seg(2:n-1,:));
% Bellhop's Dss override uses one Dx value per node: the left extension
% contributes zero at the first sample, then each finite-segment slope.
dx=[0;diff(z)./diff(x)];
kappa=(dx(2:end)-dx(1:end-1))./diff(x).*seg(:,1).^3;
segvec=points(2:end,:)-points(1:end-1,:);
q=hit-points(1:end-1,:);
lam=sum(q.*segvec,2)./sum(segvec.^2,2);
valid=lam>=-1e-8 & lam<=1+1e-8; dist=inf(n-1,1);
dist(valid)=abs(q(valid,1).*seg(valid,2)-q(valid,2).*seg(valid,1)); [~,i]=min(dist);
lambda=min(1,max(0,lam(i))); t=(1-lambda)*node(i,:)+lambda*node(i+1,:);
% Keep Bellhop's raw averaged node frame (it is not renormalized before
% Reflect2D); the normal follows the same TOP convention.
nn=[t(2) -t(1)];
geom=struct('t',t(:),'n',nn(:),'kappa',kappa(i));
end

function pat=local_pattern(f,cfg)
angles=linspace(cfg.angle_limits_deg(1),cfg.angle_limits_deg(2),2401).'; k=2*pi*f/cfg.c0_mps; th=deg2rad(angles);
D=cos(th).*exp(-0.5*(k*cfg.sigma_src_m*sin(th)).^2); D=abs(D)/max(abs(D));
pat=struct('angles_deg',angles,'level_db',20*log10(max(D,10^(cfg.source_pattern_clip_db/20))));
end
function tag=local_tag(value)
tag=strrep(sprintf('%.10g',value),'.','p'); tag=strrep(tag,'-','m');
end
