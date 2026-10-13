function summary=audit_deepsea3d_fg_artifacts()
%AUDIT_DEEPSEA3D_FG_ARTIFACTS Reconstruct paired conclusions from saved fields.
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
dir=fullfile(root,'results','validation','deepsea3d_tfe_dealiasing');
b=load(fullfile(root,'results','validation','deepsea3d_tfe_prototype','prototype.mat'),'state');
names=["dual_fg","quarter_fg","pm_fg","pm_collocation_fixed"];
initial_wall=0;
if isfile(fullfile(dir,'pm_collocation.mat')),names=[names,"pm_collocation"];end
outputs=struct();wall=0;cpu=0;nyqmax=0;rows=struct([]);mem=0;
for name=names
    a=load(fullfile(dir,name+".mat"),'output','snapshots','cfg');o=a.output;
    restored=a.cfg;restored.profile.amplitudes_m=restored.profile.amplitudes_m/o.alpha;
    if name=="dual_fg",assert(isequaln(restored,b.state.dual_cfg));else,assert(isequaln(restored,b.state.pm_cfg));end
    assert(o.options.order<=32 && a.cfg.profile.support_m==12 && o.options.spacing_m==.25);
    T=readtable(fullfile(dir,name+".csv"));assert(height(T)==numel(o.checkpoints));
    for j=1:numel(o.checkpoints)
        q=o.checkpoints(j);r=q.independent4;
        assert(abs(complex(T.delta_real(j),T.delta_imag(j))-q.rough_correction) ...
            <1e-13*max(abs(q.rough_correction),realmin));
        assert(max(q.linear_residual)<1e-8);
        assert(abs(q.independent3.retained_polynomial_mapped_PDE-r.retained_polynomial_mapped_PDE)<1e-10);
        n=round(sqrt(size(a.snapshots(j).W,1)));F=fft2(reshape(a.snapshots(j).W,n,n,[]));
        e=norm(reshape(F(n/2+1,:,:),[],1))+norm(reshape(F(:,n/2+1,:),[],1));
        nyq=e/norm(F(:));
        if name~="pm_collocation",assert(nyq<1e-12);nyqmax=max(nyqmax,nyq);end
        row=struct('job',name,'N',q.N,'delta_real',real(q.rough_correction), ...
            'delta_imag',imag(q.rough_correction),'physical_PDE',r.physical_PDE, ...
            'retained_physical_PDE',r.projected_physical_PDE,'boundary',r.true_boundary, ...
            'DtN',r.true_DtN,'lasttwo',abs(sum(q.terms(end-1:end)))/abs(q.rough_correction), ...
            'last_norm',q.volume_coefficient_norms(end),'NYquist_relative',nyq, ...
            'valid_paired_space',name~="pm_collocation");
        if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
    end
    outputs.(name)=o;wall=wall+o.wall_s;cpu=cpu+o.cpu_s;
    if name=="pm_collocation",initial_wall=o.wall_s;end
    mem=max(mem,max(o.process_before_MiB,o.process_after_MiB));
end
d=outputs.dual_fg.checkpoints(end);quarter=outputs.quarter_fg.checkpoints(end);
assert(local_good(d) && local_good(quarter));
dual_error=abs(d.rough_correction-b.state.dual_reference.rough_correction)/abs(b.state.dual_reference.rough_correction);
assert(dual_error<1e-4);
fg=outputs.pm_fg;colloc=outputs.pm_collocation_fixed;
for j=1:3,assert(fg.checkpoints(j).N==colloc.checkpoints(j).N && ~local_good(fg.checkpoints(j)));end
q=fg.checkpoints(end);c=colloc.checkpoints(end);
% Amplitude homogeneity is an exact algebraic property, tested by the
% independently executed quarter-amplitude case, not a fabricated field.
n=(1:16)';prediction=fg.checkpoints(1).terms.*.25.^n;
homogeneity=norm(prediction-quarter.terms)/norm(quarter.terms);
norm_homogeneity=norm(fg.checkpoints(1).volume_coefficient_norms.*.25.^n-quarter.volume_coefficient_norms) ...
    /norm(quarter.volume_coefficient_norms);
assert(homogeneity<1e-10 && norm_homogeneity<1e-10);
a=load(fullfile(dir,'pm_fg_geometry_probe.mat'),'probe');probe=a.probe;
assert(probe.four.geometry_operator_difference<1e-4 && probe.four.physical_PDE>1e-4);
summary=struct('controls_local_passed',true,'dual_frozen_reference_difference',dual_error, ...
    'quarter_homogeneity_pressure_error',homogeneity,'quarter_homogeneity_volume_error',norm_homogeneity, ...
    'valid_space_Nyquist_max',nyqmax,'initial_comparator_valid',false, ...
    'coefficient_norm_reduction_N32',c.volume_coefficient_norms(end)/q.volume_coefficient_norms(end), ...
    'PDE_ratio_dealias_to_collocation_N32',q.independent4.physical_PDE/c.independent4.physical_PDE, ...
    'DtN_ratio_dealias_to_collocation_N32',q.independent4.true_DtN/c.independent4.true_DtN, ...
    'retained_residual_energy_fraction', ...
    (q.independent4.projected_physical_PDE/q.independent4.physical_PDE)^2, ...
    'PM_reference','INCONCLUSIVE','global_G1','OPEN','WP2_allowed',false, ...
    'new_solves',numel(names),'support24_solves',0,'FullFK_solves',0,'higher_order_solves',0, ...
    'next_reference','bounded direct t=1 mapped Helmholtz, matrix-free; BIE fallback');
cost=struct('new_solve_wall_s',wall,'new_solve_CPU_s',cpu,'saved_field_probe_wall_s',probe.wall_s, ...
    'all_saved_kernel_wall_s',wall+probe.wall_s,'maximum_named_workspace_MiB',384.75, ...
    'maximum_process_snapshot_MiB',mem,'initial_invalid_run_wall_s',initial_wall, ...
    'scope','includes initial invalid comparator and all residual evaluations; excludes IO/tests/report; no peak instrumentation');
assert(cost.all_saved_kernel_wall_s<600);
writetable(struct2table(rows),fullfile(dir,'all_checkpoints.csv'));
local_json(fullfile(dir,'artifact_assertions.json'),summary);local_json(fullfile(dir,'cost_summary.json'),cost);
disp(summary);disp(cost);
end
function good=local_good(q)
a=q.independent3;b=q.independent4;
v=[a.physical_PDE,a.true_boundary,a.true_DtN,b.physical_PDE,b.true_boundary,b.true_DtN];
good=all(v<1e-4) && abs(sum(q.terms(end-1:end)))/abs(q.rough_correction)<1e-4;
end
function local_json(file,data)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(data,PrettyPrint=true));
end
