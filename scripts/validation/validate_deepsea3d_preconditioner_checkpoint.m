function output=validate_deepsea3d_preconditioner_checkpoint(job)
%VALIDATE_DEEPSEA3D_PRECONDITIONER_CHECKPOINT Single-route capped new cases.
arguments
    job (1,1) string {mustBeMember(job,["flat","dual","quarter","pm"])}
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
outdir=fullfile(root,'results','validation','deepsea_two_track_20261014');
file=fullfile(outdir,"right_"+job+".mat");assert(~isfile(file),'Checkpoint:Exists','No overwrite.');
a=load(fullfile(root,'results','validation','deepsea3d_tfe_prototype','prototype.mat'),'state');
cfg=a.state.pm_cfg;
switch job
    case "flat",cfg=a.state.dual_cfg;cfg.profile.amplitudes_m(:)=0;
    case "dual",cfg=a.state.dual_cfg;
    case "quarter",cfg.profile.amplitudes_m=cfg.profile.amplitudes_m*.25;
    case "pm"
        for name=["flat","dual","quarter"]
            saved=load(fullfile(outdir,"right_"+name+".mat"),'output');
            assert(saved.output.linear_pass && saved.output.initial_physical_pass,'Checkpoint:Control','Controls must pass.');
        end
end
spent=0;
for name=["flat","dual","quarter","pm"]
    f=fullfile(outdir,"right_"+name+".mat");
    if isfile(f)
        a=load(f,'output');spent=spent+a.output.wall_s;
    end
end
assert(spent+600<3600,'Checkpoint:TotalBudget','Insufficient remaining budget.');
opts=struct('domain_halfwidth_m',24,'spacing_m',.25,'vertical_degree',32,'top_depth_m',.8, ...
    'restart',12,'max_cycles',8,'krylov_tolerance',1e-10,'krylov_wall_limit_s',560, ...
    'max_operator_calls',640,'max_workspace_MiB',737.5);
[output,W]=solve_deepsea3d_right_neumann(cfg,opts);output.job=job;
if job=="pm"
    output.physics_tier="PM_unvalidated_diagnostic";
elseif job=="quarter"
    output.physics_tier="scaled_PM_diagnostic_control";
else
    output.physics_tier="accepted_control_regression";
end
save(file,'output','W','cfg','-v7');
summary=rmfield(output,{'rough_correction','direct_native','reflect_native','total_native'});
summary.rough_correction_real=real(output.rough_correction);summary.rough_correction_imag=imag(output.rough_correction);
fid=fopen(fullfile(outdir,"right_"+job+".json"),'w');assert(fid>=0);cleaner=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(summary,PrettyPrint=true));
fprintf('%s fresh %.4g truePDE %.4g boundary %.4g DtN %.4g wall %.3fs\n', ...
    job,output.fresh_preconditioned_residual,output.independent4.physical_PDE, ...
    output.independent4.true_boundary,output.independent4.true_DtN,output.wall_s);
end
