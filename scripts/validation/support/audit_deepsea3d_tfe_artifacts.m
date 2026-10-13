function summary=audit_deepsea3d_tfe_artifacts()
%AUDIT_DEEPSEA3D_TFE_ARTIFACTS Independent saved-data and gate reconstruction.
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
out=fullfile(root,'results','validation','deepsea3d_tfe_prototype');
a=load(fullfile(out,'prototype.mat'),'state');s=a.state;
pm=load(fullfile(root,'results','validation','deepsea3d_pm_feasibility','feasibility.mat'),'state');
assert(isequaln(s.pm_cfg,pm.state.cfg) && isequaln(s.frozen_sea,pm.state.sea));
assert(s.frozen_sea.seed==261013 && s.pm_cfg.frequency_hz==6000);
assert(s.WP2_allowed==false && s.FK_solves==0 && s.global_G1=="OPEN");
R=readtable(fullfile(out,'runs.csv'),TextType='string');assert(height(R)==numel(s.runs));
T=readtable(fullfile(out,'terms.csv'),TextType='string');
for j=1:numel(s.runs)
    q=s.runs(j);r=q.result;
    assert(R.name(j)==q.name);
    assert(abs(complex(R.delta_real(j),R.delta_imag(j))-r.rough_correction)<1e-16);
    tt=T(T.name==q.name,:);assert(height(tt)==r.order);
    assert(norm(complex(tt.term_real,tt.term_imag)-r.terms)<1e-14);
    assert(abs(sum(r.terms)-r.rough_correction)<1e-17);
    assert(all(isfinite([r.terms;r.volume_coefficient_norms;r.linear_residual])));
    assert(r.min_map_jacobian>0 && max(r.linear_residual)<1e-8);
end
names=string({s.runs.name});get=@(name)s.runs(find(names==name,1)).result;
d=get("dual_order");d20=get("dual_order20");fine=get("dual_horizontal");
guards=struct('order16_to20',local_change(d20,d), ...
    'horizontal',local_change(fine,d),'vertical',local_change(get("dual_vertical"),d), ...
    'domain',local_change(get("dual_domain"),d),'transparent_plane',local_change(get("dual_top"),d), ...
    'frozen_reference',local_change(fine,s.dual_reference));
dual_names=["dual_order","dual_order20","dual_horizontal","dual_vertical","dual_top","dual_domain"];
dual_local=true;
for name=dual_names,dual_local=dual_local && local_good(get(name));end
dual_ready=dual_local && all(structfun(@(v)v<1e-4,guards));
assert(dual_ready,'TFE:DualGate','Dual reference verification unresolved.');
p16=get("pm16");p24=get("pm24");p32=get("pm32");p40=get("pm40");
pm_guards=struct('order16_to24',local_change(p24,p16),'order24_to32',local_change(p32,p24), ...
    'order32_to40',local_change(p40,p32),'horizontal',local_change(get("pm_horizontal"),p32), ...
    'vertical',local_change(get("pm_vertical"),p32),'domain',local_change(get("pm_domain"),p32), ...
    'transparent_plane',local_change(get("pm_top"),p32));
assert(~local_good(p32) && ~local_good(p40) && ~local_good(get("pm_horizontal")));
scaled=get("pm_scaled");
summary=struct('dual_reference_ready',dual_ready,'dual_guards',guards, ...
    'PM_reference_ready',false,'PM_guards',pm_guards,'scaled_local_guards',local_good(scaled), ...
    'scaled_formal_reference_ready',false,'support24_runs',0,'G1_status','INCONCLUSIVE', ...
    'global_G1','OPEN','WP2_allowed',false,'FK_solves',0, ...
    'surface_loaded_not_regenerated',true,'saved_native_CSV_and_terms',true, ...
    'new_rough_solves',numel(s.runs),'flat_controls',1);
assert(summary.scaled_local_guards);
cost=struct('rough_TFE_wall_s',sum(arrayfun(@(q)q.result.wall_s,s.runs)), ...
    'rough_TFE_cpu_s',sum(arrayfun(@(q)q.result.cpu_s,s.runs)), ...
    'flat_wall_s',s.flat.result.wall_s,'flat_cpu_s',s.flat.result.cpu_s, ...
    'maximum_volume_MiB',max(arrayfun(@(q)q.result.volume_MiB,s.runs)), ...
    'maximum_workspace_estimate_MiB',max(arrayfun(@(q)q.result.estimated_workspace_MiB,s.runs)), ...
    'maximum_process_snapshot_MiB',max([s.runs.process_before_MiB,s.runs.process_after_MiB]), ...
    'scope','saved kernel/plan/forcing/residual evaluations, excludes IO/tests/report; named-array estimates and snapshots are not peaks');
local_json(fullfile(out,'artifact_assertions.json'),summary);
local_json(fullfile(out,'cost_summary.json'),cost);
writetable(struct2table(guards),fullfile(out,'dual_guards.csv'));
writetable(struct2table(pm_guards),fullfile(out,'PM_guards.csv'));
disp(summary);disp(cost);
end
function x=local_change(a,b)
x=abs(a.rough_correction-b.rough_correction)/abs(b.rough_correction);
end
function good=local_good(r)
v=[r.boundary_residual_relative,r.offgrid_boundary_residual_relative, ...
    r.PDE_residual_relative,r.offgrid_PDE_residual_relative,r.transparent_residual_relative];
tail=abs(r.partials(end)-r.partials(end-2))/max(abs(r.rough_correction),realmin);
good=all(isfinite(v)) && all(v<1e-4) && tail<1e-4 && max(r.linear_residual)<1e-8 ...
    && all(diff(r.volume_coefficient_norms(end-4:end))<0) ...
    && all(diff(abs(r.terms(end-4:end)))<0);
end
function local_json(file,x)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(x,PrettyPrint=true));
end
