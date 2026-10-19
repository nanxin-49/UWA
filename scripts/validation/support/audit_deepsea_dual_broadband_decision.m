function audit=audit_deepsea_dual_broadband_decision()
%AUDIT_DEEPSEA_DUAL_BROADBAND_DECISION Saved-data arithmetic, not physics tests.
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
outdir=fullfile(root,'results','validation','deepsea_dual_broadband_decision_20261015');
file=fullfile(outdir,'numerical_audit.json');assert(~isfile(file),'dualDecision:Exists','No overwrite.');
a=load(fullfile(outdir,'saved_evidence.mat'),'output','raw');
b=load(fullfile(outdir,'locked_candidate_and_gate.mat'),'decision');g=b.decision;
assert(~g.memory_gate_pass && strcmp(g.expansion_decision,'STOP_RESOURCE_GATE_NO_NEW_PHYSICS'));
assert(g.new_physical_frequencies==0 && g.new_physical_kernel_wall_s==0);
assert(g.estimated_extra_working_set_MiB>.7*g.physical_available_MiB);
assert(~a.output.continuous_band_certified && ~a.output.WP2_formal && strcmp(a.output.G1_target_PM,'OPEN'));
rows=struct([]);identity=0;csv_consistency=0;
for j=1:3
    x=a.raw(j);m=a.output.metrics(j);r=x.reference;
    names={'flat_direct','flat_reflect','flat_total','rough_reflect','rough_increment','rough_total', ...
        'FK_reflect','FK_increment','FK_total'};
    values=[x.flat.direct_native,x.flat.reflect_native,x.flat.total_native,r.reflect_native, ...
        r.rough_correction,r.total_native,x.fullfk.reflect_native,x.fullfk.rough_correction,x.fullfk.total_native];
    for z=1:numel(names)
        row=struct('frequency_hz',x.cfg.frequency_hz,'component',string(names{z}), ...
            'native_real',real(values(z)),'native_imag',imag(values(z)), ...
            'native_abs',abs(values(z)),'native_phase_wrapped_rad',angle(values(z)));
        if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
    end
    identity=max(identity,abs(r.total_native-x.flat.direct_native-r.reflect_native));
    identity=max(identity,abs(r.reflect_native-x.flat.reflect_native-r.rough_correction));
    assert(m.fk_reflected_model_error<.01 && m.fk_increment_model_error<.01);
    assert(m.effect_over_ref_sensitivity>1e4 && m.effect_over_model_difference>4e3);
    % Frozen CSV independently stores all three complex components per solve.
    if x.cfg.frequency_hz==6000
        folder=fullfile(root,'results','validation','deepsea3d_g1_azimuth');
    else
        folder=fullfile(root,'results','validation','deepsea3d_g1_frequency',sprintf('f%d',x.cfg.frequency_hz));
    end
    t=readtable(fullfile(folder,'reference_runs.csv'));f=readtable(fullfile(folder,'fk_runs.csv'));
    refvalue=complex(t.total_real(end),t.total_imag(end));fkvalue=complex(f.total_real(2),f.total_imag(2));
    csv_consistency=max([csv_consistency,abs(refvalue-r.total_native),abs(fkvalue-x.fullfk.total_native)]);
end
assert(identity<1e-16 && csv_consistency<1e-16);
writetable(struct2table(rows),fullfile(outdir,'native_pressures.csv'));
T=struct2table(a.output.cost_candidates);
assert(all(T.one_midpoint_total_points==2*T.source_points-1));
assert(all(T.new_points_after_reusing_three==T.one_midpoint_total_points-3));
geometry=a.output.geometry;
raw_required=ceil(4000/geometry.raw_reflected_linear_df_for_1e3_hz)+1;
% Archive file sizes are an observed storage baseline, not a saved-field budget.
storage=struct([]);
for f=[4000 8000]
    folder=fullfile(root,'results','validation','deepsea3d_g1_frequency',sprintf('f%d',f));
    files=dir(fullfile(folder,'*'));files=files(~[files.isdir]);
    row=struct('frequency_hz',f,'per_frequency_CSV_JSON_files',numel(files),'bytes',sum([files.bytes]));
    if isempty(storage),storage=row;else,storage(end+1)=row;end %#ok<AGROW>
end
audit=struct('saved_native_identities_max_abs',identity,'MAT_CSV_max_abs_difference',csv_consistency, ...
    'saved_frequency_count',3,'new_physical_frequency_count',0,'new_physical_kernel_wall_s',0, ...
    'old_conditional_physics_threshold',.01,'old_reference_threshold',1e-4,'old_FK_numerical_threshold',1e-3, ...
    'raw_relative_DSP_reflected_linear_1e3_sufficient_uniform_source_count',raw_required, ...
    'raw_linear_count_scope','analytic single reflected exponential only; no rough-band certification', ...
    'observed_export_storage',storage,'decision','B_CURRENT_DISCRETE_ONLY', ...
    'reason','Memory gate rejected optional two-point validation; candidate bandwidth unsupported by independent holdouts', ...
    'no_new_reference_failure_claim',true,'G1_target_PM','OPEN','WP2_formal',false,'passed',true);
fid=fopen(file,'w');assert(fid>=0);clean=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(audit,PrettyPrint=true));disp(audit);
end
