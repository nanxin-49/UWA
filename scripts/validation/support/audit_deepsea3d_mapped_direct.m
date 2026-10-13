function summary=audit_deepsea3d_mapped_direct()
%AUDIT_DEEPSEA3D_MAPPED_DIRECT Separate Krylov/discrete/physical gate evidence.
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
dir=fullfile(root,'results','validation','deepsea3d_mapped_direct');
a=load(fullfile(root,'results','validation','deepsea3d_tfe_prototype','prototype.mat'),'state');base=a.state;
names=["flat","dual","quarter","pm","pm_restart12"];
outputs=struct();rows=struct([]);wall=0;cpu=0;snapshot=0;
for name=names
    a=load(fullfile(dir,name+".mat"),'output','W','cfg');o=a.output;
    expected=base.pm_cfg;
    if name=="dual" || name=="flat",expected=base.dual_cfg;end
    expected.profile.amplitudes_m=expected.profile.amplitudes_m*o.alpha;
    if name=="flat",expected.profile.amplitudes_m(:)=0;end
    assert(isequaln(a.cfg,expected));
    assert(o.options.domain_halfwidth_m==24 && o.options.spacing_m==.25 ...
        && o.options.vertical_degree==32 && o.options.top_depth_m==.8);
    assert(o.options.restart*o.options.max_cycles<=96 && o.options.max_workspace_MiB==900);
    p=deepsea3d_fg_plan(a.cfg,o.options);n=p.n;
    points=[p.xx(:)+a.cfg.center_m(1),p.yy(:)+a.cfg.center_m(2),p.original_geo.eta(:)];
    image=a.cfg.tx_m;image(3)=-image(3);
    flat=deepsea3d_point_field(points,a.cfg.tx_m,p.k)-deepsea3d_point_field(points,image,p.k);
    g=deepsea3d_fg_resize(-reshape(flat,n,n),n);
    % Independent audit recomputes the unpreconditioned discrete equations.
    r=deepsea3d_mapped_discrete_residual(a.W,g,p);
    assert(abs(r.mapped_PDE-o.discrete.mapped_PDE)<1e-13);
    assert(abs(r.Dirichlet-o.discrete.Dirichlet)<1e-13 && abs(r.DtN-o.discrete.DtN)<1e-13);
    linear=all([o.fresh_preconditioned_residual,r.mapped_PDE,r.Dirichlet,r.DtN]<1e-8);
    physical=all([o.independent3.physical_PDE,o.independent3.true_boundary,o.independent3.true_DtN, ...
        o.independent4.physical_PDE,o.independent4.true_boundary,o.independent4.true_DtN]<1e-4);
    assert(linear==o.linear_pass && physical==o.initial_physical_pass);
    if any(name==["flat","dual","quarter"]),assert(linear && physical);else,assert(~linear && ~physical);end
    assert(strcmp(o.PM_status,'INCONCLUSIVE') && strcmp(o.global_G1,'OPEN') && ~o.WP2_allowed);
    T=readtable(fullfile(dir,name+".csv"));
    assert(abs(complex(T.delta_real(1),T.delta_imag(1))-o.rough_correction)<1e-13*max(abs(o.rough_correction),realmin) ...
        || o.rough_correction==0);
    bg=deepsea3d_flat(a.cfg,0);
    assert(abs(o.reflect_native-bg.reflect_native-o.rough_correction)<1e-17);
    assert(abs(o.total_native-o.reflect_native-o.direct_native)<1e-17);
    if ~isempty(o.cycles)
        assert(o.cycles(end).operator_calls<=o.options.max_operator_calls);
        assert(o.cycles(end).wall_s<=o.options.krylov_wall_limit_s+5);
    end
    row=struct('job',name,'fresh_B',o.fresh_preconditioned_residual,'discrete_PDE',r.mapped_PDE, ...
        'true_PDE3',o.independent3.physical_PDE,'true_PDE4',o.independent4.physical_PDE, ...
        'true_boundary',o.independent4.true_boundary,'true_DtN',o.independent4.true_DtN, ...
        'linear_pass',linear,'initial_physical_pass',physical,'delta_real',real(o.rough_correction), ...
        'delta_imag',imag(o.rough_correction),'wall_s',o.wall_s);
    if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
    outputs.(name)=o;wall=wall+o.wall_s;cpu=cpu+o.cpu_s;
    snapshot=max(snapshot,max(o.process_before_MiB,o.process_after_MiB));
end
d=outputs.dual;dual_error=abs(d.rough_correction-base.dual_reference.rough_correction)/abs(base.dual_reference.rough_correction);
assert(dual_error<1e-4);
cert=certify_deepsea3d_mapped_geometry(base.pm_cfg,.8);
local_json(fullfile(dir,'continuous_map_certificate.json'),cert);
old=load(fullfile(root,'results','validation','deepsea3d_tfe_dealiasing','pm_fg.mat'),'output');
taylor=old.output.checkpoints(end).independent4;pm=outputs.pm_restart12;
summary=struct('controls_passed',true,'dual_frozen_reference_difference',dual_error, ...
    'PM_krylov_converged',false,'PM_discrete_passed',false,'PM_physical_passed',false, ...
    'restart8_to12_pressure_change',abs(pm.rough_correction-outputs.pm.rough_correction)/abs(pm.rough_correction), ...
    'PDE_improvement_over_Taylor32',taylor.physical_PDE/pm.independent4.physical_PDE, ...
    'DtN_improvement_over_Taylor32',taylor.true_DtN/pm.independent4.true_DtN, ...
    'PM_reference','INCONCLUSIVE','global_G1','OPEN','WP2_allowed',false, ...
    'no_Taylor_source',true,'no_fitting',true,'new_channel_solves',5, ...
    'FullFK_solves',0,'support24_solves',0,'grid_refinement_solves',0);
cost=struct('saved_kernel_wall_s',wall,'saved_kernel_CPU_s',cpu, ...
    'maximum_named_workspace_MiB',max(arrayfun(@(v)outputs.(v).estimated_workspace_MiB,names)), ...
    'maximum_process_snapshot_MiB',snapshot,'scope', ...
    'plans/source/GMRES/discrete and independent physical residuals; excludes IO/tests/audit/report; not peak monitoring');
assert(wall<900);
writetable(struct2table(rows),fullfile(dir,'summary.csv'));
local_json(fullfile(dir,'artifact_assertions.json'),summary);local_json(fullfile(dir,'cost_summary.json'),cost);
disp(summary);disp(cost);
end
function local_json(file,value)
fid=fopen(file,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(value,PrettyPrint=true));
end
