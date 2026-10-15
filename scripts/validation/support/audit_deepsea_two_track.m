function audit=audit_deepsea_two_track()
%AUDIT_DEEPSEA_TWO_TRACK Recheck new fields against frozen operators/fixtures.
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
outdir=fullfile(root,'results','validation','deepsea_two_track_20261014');timer=tic;
prototype=load(fullfile(root,'results','validation','deepsea3d_tfe_prototype','prototype.mat'),'state');
baseline=load(fullfile(root,'results','validation','deepsea3d_mapped_direct','pm_restart12.mat'),'output');
rows=struct([]);spent=0;
for name=["flat","dual","quarter","pm"]
    a=load(fullfile(outdir,"right_"+name+".mat"),'output','W','cfg');o=a.output;
    expected=prototype.state.pm_cfg;
    if name=="dual" || name=="flat",expected=prototype.state.dual_cfg;end
    if name=="flat",expected.profile.amplitudes_m(:)=0;end
    if name=="quarter",expected.profile.amplitudes_m=expected.profile.amplitudes_m*.25;end
    assert(isequaln(a.cfg,expected),'Audit:Surface','Exactly frozen coefficient/config match.');
    p=deepsea3d_fg_plan(a.cfg,o.options);n=p.n;
    points=[p.xx(:)+a.cfg.center_m(1),p.yy(:)+a.cfg.center_m(2),p.original_geo.eta(:)];
    image=a.cfg.tx_m;image(3)=-image(3);
    flat=deepsea3d_point_field(points,a.cfg.tx_m,p.k)-deepsea3d_point_field(points,image,p.k);
    g=deepsea3d_fg_resize(-reshape(flat,n,n),n);
    r=deepsea3d_mapped_discrete_residual(a.W,g,p);
    assert(max(abs([r.mapped_PDE,r.Dirichlet,r.DtN]- ...
        [o.discrete.mapped_PDE,o.discrete.Dirichlet,o.discrete.DtN]))<1e-13);
    assert(o.operator_calls<=640 && o.wall_s<=600 && ~o.WP2_formal && strcmp(o.G1_target_PM,'OPEN'));
    assert(o.options.restart*o.options.max_cycles==96);
    if name~="pm",assert(o.linear_pass && o.initial_physical_pass); ...
    else,assert(~o.linear_pass && ~o.initial_physical_pass);end
    row=struct('job',name,'fresh_B',o.fresh_preconditioned_residual,'discrete_PDE',r.mapped_PDE, ...
        'true_PDE',o.independent4.physical_PDE,'true_Dirichlet',o.independent4.true_boundary, ...
        'true_DtN',o.independent4.true_DtN,'B_calls',o.operator_calls,'wall_s',o.wall_s, ...
        'linear_pass',o.linear_pass,'physical_pass',o.initial_physical_pass);
    if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
    spent=spent+o.wall_s;
    if name=="dual"
        fixture=readtable(fullfile(root,'results','validation','deepsea3d_g1_frequency','frequency_summary.csv'));
        ref=complex(fixture.reference_delta_real(2),fixture.reference_delta_imag(2));
        dual_error=abs(o.rough_correction-ref)/abs(ref);assert(dual_error<1e-4);
    end
end
assert(spent<3600);
b=load(fullfile(outdir,'singlelink_readiness.mat'),'output','responses','dualmeta');B=b.output;
assert(B.provenance_rejections==5 && all(B.flat_two_path_bit_errors==0) && all(B.dual_scalar_bit_errors==0));
assert(~B.WP2_formal && strcmp(B.G1_target_PM,'OPEN'));
for j=1:numel(b.responses)
    s=b.responses(j);
    assert(max(abs(s.total_native-s.direct_native-s.reflect_native))<1e-15);
    assert(norm(s.absolute_DSP-conj(s.total_native))<1e-15);
    assert(norm(s.relative_DSP-s.absolute_DSP.*exp(1i*2*pi*s.f_axis_hz*s.tau_direct_s))<1e-15);
end
audit=struct('status','PASS','A_rows',rows,'A_new_physical_kernel_wall_s',spent, ...
    'baseline_fresh_B',baseline.output.fresh_preconditioned_residual, ...
    'baseline_true_PDE',baseline.output.independent4.physical_PDE, ...
    'baseline_true_DtN',baseline.output.independent4.true_DtN, ...
    'dual_increment_reference_error',dual_error, ...
    'A_decision','STOP_ONE_INTERVENTION_FAILED; matched 3D BIE feasibility only', ...
    'B_status',B.status,'linear_interpolation_513_1e3_pass',false, ...
    'generic_broadband_grid_accepted',false,'PM_linear_solved',false, ...
    'reference_candidate',false,'reference_ready',false,'G1_target_PM','OPEN','WP2_formal',false, ...
    'audit_wall_s',toc(timer));
writetable(struct2table(rows),fullfile(outdir,'A_comparison.csv'));
fid=fopen(fullfile(outdir,'data_audit.json'),'w');assert(fid>=0);cleaner=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(audit,PrettyPrint=true));
fprintf('Data audit PASS; A total kernel %.3fs; dual relative error%.4g; audit %.3fs\n',spent,dual_error,audit.audit_wall_s);
end
