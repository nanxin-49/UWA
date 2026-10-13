function output=validate_deepsea3d_tfe_dealiasing(job)
%VALIDATE_DEEPSEA3D_TFE_DEALIASING Bounded independent operator experiment.
arguments
    job (1,1) string
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
out=fullfile(root,'results','validation','deepsea3d_tfe_dealiasing');
file=fullfile(out,job+".mat");assert(~isfile(file),'FG:Exists','Do not overwrite a saved experiment.');
if isfile(fullfile(out,'artifact_assertions.json')) && any(job==["half_fg","pm_fine_fg"])
    error('FG:StudyStopped','Completed diagnosis stops expansion; use a separately scoped future study.');
end
a=load(fullfile(root,'results','validation','deepsea3d_tfe_prototype','prototype.mat'),'state');
s=a.state;kind="galerkin";alpha=1;h=.25;N=32;snaps=[16 24 32];cfg=s.pm_cfg;
switch job
    case "dual_fg",cfg=s.dual_cfg;N=16;snaps=16;
    case "quarter_fg",alpha=.25;N=16;snaps=16;
    case "pm_collocation",error('FG:RetiredComparator','Use pm_collocation_fixed; original data are retained read-only.');
    case "pm_collocation_fixed",kind="collocation";
    case "pm_fg"
    case "half_fg",alpha=.5;
    case "pm_fine_fg",h=.125;
    otherwise,error('FG:Job','Undefined bounded experiment.');
end
cfg.profile.amplitudes_m=cfg.profile.amplitudes_m*alpha;
opts=struct('domain_halfwidth_m',24,'spacing_m',h,'vertical_degree',32,'top_depth_m',.8, ...
    'order',N,'checkpoint_orders',snaps,'product_kind',kind,'max_workspace_MiB',900);
% Base-grid runs stay far below this ceiling; fine-grid execution requires a
% separate capacity check and is not automatically scheduled.
[before,sys]=memory;
estimate=12*16*(2*ceil(24/h))^2*33/2^20+18*16*(4*2*ceil(24/h))^2/2^20;
assert(estimate<=opts.max_workspace_MiB,'FG:MemoryBudget','Named-array budget exceeded before allocation.');
assert(estimate<.85*sys.PhysicalMemory.Available/2^20,'FG:PhysicalMemory','Insufficient available RAM for this bounded run.');
[output,snapshots]=solve_deepsea3d_tfe_galerkin(cfg,opts);[after,~]=memory;
output.job=job;output.alpha=alpha;output.process_before_MiB=before.MemUsedMATLAB/2^20;
output.process_after_MiB=after.MemUsedMATLAB/2^20;
save(file,'output','snapshots','cfg','-v7');
rows=struct([]);
for j=1:numel(output.checkpoints)
    q=output.checkpoints(j);r=q.independent4;
    row=struct('N',q.N,'delta_real',real(q.rough_correction),'delta_imag',imag(q.rough_correction), ...
        'physical_PDE',r.physical_PDE,'retained_PDE',r.retained_true_mapped_PDE, ...
        'poly_PDE',r.coefficient_interpolant_PDE,'poly_retained_PDE',r.retained_polynomial_mapped_PDE, ...
        'alias_operator',r.collocation_alias_operator,'boundary',r.true_boundary,'DtN',r.true_DtN, ...
        'tailtwo',abs(sum(q.terms(end-1:end)))/abs(q.rough_correction), ...
        'last_coefficient_norm',q.volume_coefficient_norms(end));
    if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
end
writetable(struct2table(rows),fullfile(out,job+".csv"));
end
