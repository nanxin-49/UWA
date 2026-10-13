function output=validate_deepsea3d_mapped_direct(job)
%VALIDATE_DEEPSEA3D_MAPPED_DIRECT Bounded staged matrix-free reference pilot.
arguments
    job (1,1) string
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
dir=fullfile(root,'results','validation','deepsea3d_mapped_direct');file=fullfile(dir,job+".mat");
assert(~isfile(file),'Mapped:Exists','Do not overwrite prior results.');
a=load(fullfile(root,'results','validation','deepsea3d_tfe_prototype','prototype.mat'),'state');
cfg=a.state.pm_cfg;alpha=1;cycles=8;time=180;calls=100;restart=8;
switch job
    case "flat",cfg=a.state.dual_cfg;cfg.profile.amplitudes_m(:)=0;
    case "dual",cfg=a.state.dual_cfg;
    case "quarter",alpha=.25;
    case {"pm","pm_restart12"}
        for name=["flat","dual","quarter"]
            b=load(fullfile(dir,name+".mat"),'output');
            assert(b.output.linear_pass && b.output.initial_physical_pass,'Mapped:Controls','Controls must pass before original PM.');
        end
        cycles=12;time=300;calls=160;
        if job=="pm_restart12",restart=12;cycles=8;end
    otherwise,error('Mapped:Job','Undefined bounded job.');
end
cfg.profile.amplitudes_m=cfg.profile.amplitudes_m*alpha;
spent=0;
for previous=["flat","dual","quarter","pm","pm_restart12"]
    oldfile=fullfile(dir,previous+".mat");
    if isfile(oldfile),saved=load(oldfile,'output');spent=spent+saved.output.wall_s;end
end
remaining=900-spent;assert(remaining>60,'Mapped:StudyBudget','No budget for another case and residual audit.');
time=min(time,remaining-60);
opts=struct('domain_halfwidth_m',24,'spacing_m',.25,'vertical_degree',32,'top_depth_m',.8, ...
    'restart',restart,'max_cycles',cycles,'krylov_tolerance',1e-10,'krylov_wall_limit_s',time, ...
    'max_operator_calls',calls,'max_workspace_MiB',900);
[output,W]=solve_deepsea3d_mapped_direct(cfg,opts);output.job=job;output.alpha=alpha;
save(file,'output','W','cfg','-v7');
if ~isempty(output.cycles),writetable(struct2table(output.cycles),fullfile(dir,job+"_cycles.csv"));end
q=output.independent4;
row=struct('job',job,'alpha',alpha,'delta_real',real(output.rough_correction),'delta_imag',imag(output.rough_correction), ...
    'fresh_preconditioned',output.fresh_preconditioned_residual,'discrete_PDE',output.discrete.mapped_PDE, ...
    'discrete_Dirichlet',output.discrete.Dirichlet,'discrete_DtN',output.discrete.DtN, ...
    'physical_PDE',q.physical_PDE,'true_boundary',q.true_boundary,'true_DtN',q.true_DtN, ...
    'linear_pass',output.linear_pass,'initial_physical_pass',output.initial_physical_pass, ...
    'wall_s',output.wall_s,'cpu_s',output.cpu_s);
writetable(struct2table(row),fullfile(dir,job+".csv"));
end
