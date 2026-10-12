function checks = audit_fk_single_receiver_methods_artifacts()
%AUDIT_FK_SINGLE_RECEIVER_METHODS_ARTIFACTS Saved-evidence audit, no solver run.
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
outdir=fullfile(root,'results','validation','fk_single_receiver_methods');
loaded=load(fullfile(outdir,'evaluation.mat'),'result'); r=loaded.result;
assert(r.pass && r.references_unchanged);
assert(isequal(r.frequencies,[4000 6000 8000]) && isequal(r.seeds,260001:260012));
assert(height(r.accuracy)==102 && height(r.timing)==45);
assert(height(r.guards)==75 && height(r.plan_controls)==141 && height(r.cancellation)==6);
assert(all(r.guards.pass) && all(r.accuracy.relative_error<1e-3) && all(r.accuracy.finite_metrics));
assert(r.hardening.freeze_eligible && r.hardening.pressure_unchanged);
local_metrics_audit(r.accuracy); local_metrics_audit(r.cancellation);
assert(all(r.plan_controls.pass) && all(r.plan_controls.value<=r.plan_controls.tolerance));
negative=~cellfun(@isempty,r.plan_controls.expected_id);
assert(all(strcmp(r.plan_controls.expected_id(negative),r.plan_controls.actual_id(negative))));
assert(all(isfinite([r.cancellation.rank5_relative_error; r.cancellation.nodes65_relative_error; ...
    r.cancellation.quadrature4097_8193_relative;r.cancellation.cancellation_ratio])));
assert(isequal(r.cancellation.original_relative_gate_pass,r.cancellation.relative_error<1e-3));
recomputed=abs(complex(r.accuracy.approx_real,r.accuracy.approx_imag)- ...
    complex(r.accuracy.full_real,r.accuracy.full_imag));
assert(max(abs(recomputed-r.accuracy.absolute_error))<1e-15);
for j=1:numel(r.reference_paths)
    assert(strcmp(local_hash(r.reference_paths{j}),r.reference_hashes{j}),'Reference drift.');
end
for j=1:numel(r.code_paths)
    assert(strcmp(local_hash(r.code_paths{j}),r.code_hashes{j}),'Execution code drift; rerun evaluation.');
end
for j=1:numel(r.artifact_names)
    assert(strcmp(local_hash(fullfile(outdir,r.artifact_names{j})),r.artifact_hashes{j}),'CSV drift.');
end
manifest=jsondecode(fileread(fullfile(outdir,'freeze_manifest.json')));
assert(strcmp(manifest.evaluation_sha256,local_hash(fullfile(outdir,'evaluation.mat'))));
assert(strcmp(manifest.previous_evaluation_sha256,local_hash(fullfile(outdir,'pre_hardening','evaluation.mat'))));
assert(manifest.relative_complex_gate==1e-3 && manifest.product_identity_gate==1e-12 && ...
    manifest.node_refinement_gate==1e-5 && manifest.rank_refinement_gate==1e-3 && ...
    manifest.quadrature_refinement_gate==1e-3 && manifest.reuse_identity_gate==1e-12);
assert(manifest.freeze_eligible==r.hardening.freeze_eligible && ...
    manifest.cancellation_relative_gate_pass==r.hardening.cancellation_relative_gate_pass);
for j=1:numel(r.profiles)
    t=readtable(fullfile(outdir,sprintf('coefficients_%d.csv',r.seeds(j))));
    p=r.profiles{j};
    assert(max(abs(t.cos_coeff_m-p.cos_coeff_m))<1e-14);
    assert(max(abs(t.sin_coeff_m-p.sin_coeff_m))<1e-14);
end
for f=r.frequencies
    a=r.accuracy(r.accuracy.frequency_hz==f,:); assert(height(a)==34);
    assert(sum(strcmp(a.method,'fastin'))==17 && sum(strcmp(a.method,'composite'))==17);
    ctl=r.plan_controls(r.plan_controls.frequency_hz==f,:);
    assert(height(ctl)==47 && sum(strcmp(ctl.method,'verifier'))==7);
    for method={'fastin','composite'}
        ct=ctl(strcmp(ctl.method,method{1}),:);
        assert(height(ct)==20 && numel(unique(ct.case_name))==20);
        expected={'height','source','source_kx','normalization','period','origin','frequency','sound_speed', ...
            'dispersion','grid','Tx','bound','nodes','rank','Rx_x','Rx_z','metadata','surface','kind','missing'};
        assert(isequal(sort(ct.case_name),sort(expected(:))));
    end
    for m={'full','fastin','composite'}
        t=r.timing(r.timing.frequency_hz==f & strcmp(r.timing.method,m{1}),:);
        assert(height(t)==5 && all(isfinite(t.per_realization_s)) && all(isfinite(t.cold_s)) && ...
            all(t.per_realization_s>0) && all(t.cold_s>0));
        assert(max(abs(t.cold_s-t.setup_s-t.first_apply_s))<1e-14);
    end
end
checks=struct('pass',true,'accuracy_rows',height(r.accuracy),'timing_rows',height(r.timing), ...
    'guard_rows',height(r.guards),'plan_control_rows',height(r.plan_controls), ...
    'cancellation_rows',height(r.cancellation),'primary_freeze_eligible',r.hardening.freeze_eligible, ...
    'cancellation_relative_gate_pass',r.hardening.cancellation_relative_gate_pass, ...
    'reference_hashes_unchanged',true,'execution_hashes_match',true,'artifact_hashes_match',true);
disp(checks);
end
function local_metrics_audit(t)
for j=1:height(t)
    ref=complex(t.full_real(j),t.full_imag(j)); a=complex(t.approx_real(j),t.approx_imag(j));
    m=fk_single_receiver_error_metrics(ref,a,t.diagnostic_scale(j));
    names=fieldnames(m);
    for k=1:numel(names)
        name=names{k}; actual=t.(name)(j); expected=m.(name);
        assert((isnan(actual) && isnan(expected)) || ...
            (isfinite(actual) && isfinite(expected) && abs(double(actual)-double(expected))<1e-14), ...
            'fk_validation:MetricAuditMismatch','Metric does not match complex fields: %s.',name);
    end
end
end
function h=local_hash(path)
fid=fopen(path,'rb'); assert(fid>=0); closer=onCleanup(@()fclose(fid));
bytes=fread(fid,Inf,'*uint8'); md=java.security.MessageDigest.getInstance('SHA-256'); md.update(bytes);
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
