function audit = audit_bellhop_fig12b_artifacts(output_dir)
%AUDIT_BELLHOP_FIG12B_ARTIFACTS Check saved native records, not physics PASS.
arguments
    output_dir (1,:) char
end
s=load(fullfile(output_dir,'plan.mat'),'plan');
audit=struct('case_count',numel(s.plan),'artifact_checks_passed',false, ...
    'physical_reproduction_passed',false,'nonfinite_count',0, ...
    'storage_warning_count',0,'zero_axis_forced_by_native_code',true);
for ii=1:numel(s.plan)
    t=load(fullfile(output_dir,[s.plan(ii).tag '.mat']),'output'); o=t.output;
    assert(isfile(o.env_file) && isfile(o.prt_file));
    assert(o.exit_status==0);
    audit.storage_warning_count=audit.storage_warning_count+o.storage_warnings;
    if isempty(o.pressure), continue; end
    assert(isequal(size(o.pressure),[numel(o.depth_m) numel(o.range_m)]));
    assert(max(abs(o.range_m-o.cfg.range_m))<1e-5);
    assert(max(abs(o.depth_m-o.cfg.depth_m(:)))<1e-3);
    audit.nonfinite_count=audit.nonfinite_count+nnz(~isfinite(o.pressure));
    if any(o.range_m==0), assert(o.zero_axis_all_zero); end
end
assert(audit.nonfinite_count==0);
a=load(fullfile(output_dir,'initial_n50_h10.mat'),'output');
b=load(fullfile(output_dir,'left_only.mat'),'output');
audit.independent_left_right_relative_l2= ...
    norm(a.output.pressure-fliplr(b.output.pressure),'fro')/norm(a.output.pressure,'fro');
assert(audit.independent_left_right_relative_l2<1e-5);
s=load(fullfile(output_dir,'strict_vertical_field.mat'),'output');
audit.strict_vertical_field_zero_fraction=s.output.zero_fraction;
assert(audit.strict_vertical_field_zero_fraction==1);
s=load(fullfile(output_dir,'strict_vertical_ray.mat'),'output'); ray=s.output.rays;
audit.strict_vertical_ray_max_abs_range_m=max(abs(ray.rz_m(:,1)));
audit.strict_vertical_ray_top_bounces=ray.top_bounces;
audit.strict_vertical_ray_bottom_bounces=ray.bottom_bounces;
assert(audit.strict_vertical_ray_max_abs_range_m<1e-6);
assert(ray.top_bounces>0 && ray.bottom_bounces>0);
t=readtable(fullfile(output_dir,'analytic_control.csv'),'TextType','string');
mask=t.tag=="control_narrow_2000" & t.samples>0;
audit.narrow_control_max_p95_tl_error_db=max(t.p95_tl_error_db(mask));
assert(audit.narrow_control_max_p95_tl_error_db<1e-3);
s=load(fullfile(output_dir,'surface_patch.mat'),'output');
audit.surface_pressure_max=max(abs(s.output.pressure(1,:)));
audit.surface_pressure_median=median(abs(s.output.pressure(1,:)));
audit.artifact_checks_passed=true;
fid=fopen(fullfile(output_dir,'artifact_audit.json'),'w'); assert(fid>=0);
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(audit,'PrettyPrint',true));
clear cleanup
disp(audit);
end
