function receipt=summarize_r2_qualified_closeout(case_root)
%SUMMARIZE_R2_QUALIFIED_CLOSEOUT Qualified three-seed saved-data pairing.
pre=jsondecode(fileread(fullfile(case_root,'source_preflight.json')));
root=fileparts(fileparts(fileparts(mfilename('fullpath')))); oldpath=path; pg=onCleanup(@()path(oldpath));
addpath(fullfile(root,'scripts','validation','support'),'-begin');
work_dirs=cell(1,3); registry=table; baseline=table;
for j=1:3
    seed=260000+j; rec=pre.sources([pre.sources.seed]==seed);
    if j==1
        prior=jsondecode(fileread(rec.prior_formal_review.path));
        assert(strcmp(prior.status,'SAVED_RESULT_REVIEW_PASS') && prior.failed_checks==0);
        assert(strcmp(r1_bie_manual_io('sha',rec.prior_formal_review.path),rec.prior_formal_review.sha256));
        work_dirs{j}=rec.work; qualification="FORMAL_PRIOR_REVIEW_REVALIDATED_SHA";
        audit_sha=rec.prior_formal_review.sha256; evidence=string(rec.prior_formal_review.path);
    else
        task=fullfile(case_root,sprintf('seed%d',seed));
        seal=local_verify(task,[rec.manifest_sha256 ':LOCAL-SAVED-CLOSEOUT']);
        assert(strcmp(seal.status,'AUDITED_PASS'));
        n=r1_bie_manual_io('read',fullfile(task,'closeout_receipt.json'));
        assert(strcmp(n.status,'FORMAL_LOCAL_SAVED_CLOSEOUT_PASS') && n.failed_checks==0);
        work_dirs{j}=fullfile(task,'work'); qualification="FORMAL_THIS_LOCAL_SAVED_AUDIT";
        audit_sha=r1_bie_manual_io('sha',fullfile(task,'committed.json')); evidence=string(fullfile(task,'committed.json'));
    end
    registry=[registry;table(seed,qualification,string(evidence),string(audit_sha),string(rec.original_status), ...
        'VariableNames',{'seed','qualification','evidence_path','evidence_sha256','historical_remote_status'})]; %#ok<AGROW>
    old=readtable(fullfile(rec.work,'model_metrics.csv')); old=old(ismember(string(old.method),["PS","FullFK"]),:);
    old.seed=repmat(seed,height(old),1); old.snapshot_qualification=repmat("ORIGINAL_SAVED_SNAPSHOT_NOT_A_NEW_AUDIT",height(old),1);
    old.historical_remote_status=repmat(string(rec.original_status),height(old),1); baseline=[baseline;old]; %#ok<AGROW>
end
out=fullfile(case_root,'three_seed_summary'); assert(~isfolder(out)); timer=tic;
raw=summarize_r2_pm_fields(work_dirs,out); assert(raw.failed_checks==0);
M=readtable(fullfile(out,'seed_metrics.csv')); P=readtable(fullfile(out,'point_errors.csv'));
assert(height(M)==6 && height(P)==2466);
M.qualification=strings(6,1); M.historical_remote_status=strings(6,1); M.formal_pass=true(6,1);
for j=1:3
    q=load(fullfile(work_dirs{j},'common.mat'),'x'); a=load(fullfile(work_dirs{j},'comparison.mat'));
    for method=["PS","FullFK"]
        m=1; if method=="FullFK", m=3; end
        row=find(M.seed==260000+j & string(M.method)==method); assert(isscalar(row));
        e=a.U(:,m)-a.U(:,2); E=norm(e)/norm(a.U(:,2));
        uncertainty=a.punc; if m==3, uncertainty=a.func; end
        limit=5*norm(a.bunc+uncertainty)/norm(a.U(:,2));
        assert(abs(M.complex_relative_L2(row)-E)<1e-12 && abs(M.five_times_empirical_resolution(row)-limit)<1e-12 && M.resolved(row)==(E>limit));
        ii=P.seed==260000+j & string(P.method)==method;
        assert(isequal(P.x_m(ii),q.x) && max(abs(P.residual_abs_Pa(ii)-abs(e)))<1e-11);
        assert(max(abs(P.amplitude_difference_Pa(ii)-(abs(a.U(:,m))-abs(a.U(:,2)))))<1e-11);
        assert(max(abs(P.amplitude_error_energy_Pa2(ii)+P.phase_error_energy_Pa2(ii)-abs(e).^2))<1e-10);
        M.qualification(row)=registry.qualification(j); M.historical_remote_status(row)=registry.historical_remote_status(j);
        old=baseline(baseline.seed==260000+j & string(baseline.method)==method,:);
        assert(abs(old.complex_relative_L2-E)<1e-12,'r2summary:ChangedMetric','Original saved metric changed');
    end
end
M.complex_relative_L2_percent=100*M.complex_relative_L2; M.phase_energy_percent=100*M.phase_energy_fraction;
writetable(M,fullfile(out,'qualified_metrics.csv')); writetable(registry,fullfile(out,'qualification_registry.csv'));
writetable(baseline,fullfile(out,'original_stored_metrics.csv'));
R=readtable(fullfile(out,'disjoint_energy_regions.csv'));
for seed=260001:260003
    for m=["PS","FullFK"]
        ii=R.seed==seed & string(R.method)==m;
        assert(sum(R.points(ii))==411 && abs(sum(R.residual_energy_fraction(ii))-1)<1e-12);
    end
end
render_r2_stable_saved(work_dirs,fullfile(out,'figures'));
figures=audit_r2_final_figures(work_dirs,fullfile(out,'figures'),fullfile(out,'figure_audit'));
assert(figures.failed_checks==0);
receipt=struct('status','QUALIFIED_THREE_SEED_SAVED_SUMMARY_PASS','failed_checks',0,'elapsed_s',toc(timer), ...
    'formal_rows',height(M),'snapshot_only_rows',height(baseline),'point_rows',height(P),'energy_regions',height(R), ...
    'figure_audit',figures,'new_acoustic_solves',0,'SSH_connections',0, ...
    'qualification_note','R1 prior formal review is SHA revalidated; both R2 seeds independently audited and locally sealed in this task. Original snapshot rows carry no new audit claim.', ...
    'scope','Three fixed 4kHz PM realizations and 411-point 2D receiver line; no population inference, R3 or broadband');
r1_bie_manual_io('write',fullfile(out,'qualified_receipt.json'),receipt);
fprintf('THREE_SEED formal_rows=%d failed=%d elapsed_s=%.6f\n',height(M),receipt.failed_checks,receipt.elapsed_s);
end
function s=local_verify(root,signature)
% The approved v2 SHA primitive handles empty diary files; v1 does not.
s=r1_bie_manual_io('read',fullfile(root,'committed.json'));
assert(s.schema_version==1 && strcmp(s.signature,signature) && ~isempty(s.files));
assert(numel(unique({s.files.path}))==numel(s.files));
for file_index=1:numel(s.files)
    p=r1_bie_manual_io('path',root,s.files(file_index).path);
    assert(strcmp(r1_bie_continuation_io('sha',p),s.files(file_index).sha256),'r2summary:Seal','Local sealed artifact drift');
end
end
