function receipt=validate_r1_bie_c4_preparation(root,out)
%VALIDATE_R1_BIE_C4_PREPARATION Local acceptance. Never invokes C4.
% Audits saved C3 density during migration; BIE/PS/FullFK solves are sentinels.
arguments
    root (1,:) char
    out (1,:) char
end
assert(~isfolder(out),'Use a new validation output directory'); mkdir(out);
previous_threads=maxNumCompThreads(1); thread_guard=onCleanup(@()maxNumCompThreads(previous_threads));
diary(fullfile(out,'diary.log')); diary_guard=onCleanup(@()diary('off'));
timer=tic; passed=0; failures=struct('name',{},'identifier',{},'message',{});
root=char(java.io.File(root).getCanonicalPath()); out=char(java.io.File(out).getCanonicalPath());
code=fullfile(root,'code'); previous_path=path; path_guard=onCleanup(@()path(previous_path)); addpath(code,'-begin');
sentinels=fullfile(out,'solve_sentinels'); mkdir(sentinels);
for name={'solve_helmholtz_bie_halfplane_vertical','solve_helmholtz_bie_halfplane_cached_vertical', ...
        'run_pe_1d_surface_reflection_validation','closeout_s3a_pm_fullfk'}
    fid=fopen(fullfile(sentinels,[name{1} '.m']),'w'); assert(fid>=0);
    logpath=strrep(fullfile(out,'forbidden_solve_calls.txt'),'''','''''');
    fprintf(fid,'function varargout=%s(varargin)\nfid=fopen(''%s'',''a''); fprintf(fid,''CALLED\\n''); fclose(fid); error(''r1v2:ForbiddenSolve'',''No density/PS/FK solve allowed in preparation'');\nend\n',name{1},logpath);
    fclose(fid);
end
addpath(sentinels,'-begin'); rehash;
try
    assert(strcmp(which('solve_helmholtz_bie_halfplane_vertical'),fullfile(sentinels,'solve_helmholtz_bie_halfplane_vertical.m')));
    assert(~isfile(fullfile(root,'jobs','bie_window','data.mat')));
    r1_bie_continuation_io('verify_package',root); r1_bie_continuation_io('verify_origin',root); passed=passed+1;
    policy=r1_bie_continuation_io('read',fullfile(root,'policy.json'));
    x=linspace(-40,40,411)'; z=3*ones(411,1); one=ones(411,1);
    for fine=[.5 1 2]
        fields={one+1e-7,one,one+1e-7,one+fine*1e-7,one};
        g=r1_bie_continuation_gates(fields,x,z,policy);
        assert(strcmp(g.status,'PASS') && g.grid_diagnostic.grid_nonmonotonic==(g.pairs(2).relative>g.pairs(1).relative));
        if fine==2, assert(g.grid_diagnostic.grid_nonmonotonic && strcmp(g.legacy_policy_status,'FAIL')); end
        passed=passed+1;
    end
    equal=r1_bie_continuation_gates({one,one,one,one},x,z,policy);
    assert(~equal.grid_diagnostic.grid_nonmonotonic && equal.grid_diagnostic.first==0 && ...
        contains(jsonencode(equal.grid_diagnostic),'"ratio":null')); passed=passed+1;
    ref=zeros(411,1); ref(1)=1; exact=ref; exact(2)=1e-4;
    g=r1_bie_continuation_gates({exact,ref,ref,ref,ref},x,z,policy);
    assert(g.pairs(1).relative==1e-4 && strcmp(g.status,'FAIL')); passed=passed+1;
    g=r1_bie_continuation_gates({ref+1e-3,ref,ref,ref,ref},x,z,policy);
    assert(strcmp(g.status,'FAIL')); passed=passed+1;
    bad=one; bad(7)=NaN;
    local_reject(@()r1_bie_continuation_gates({bad,one,one,one},x,z,policy)); passed=passed+1;
    local_reject(@()r1_bie_continuation_gates({one,one,one},x,z,policy)); passed=passed+1;
    local_reject(@()r1_bie_continuation_gates({one,one,one,one},x,4*ones(411,1),policy)); passed=passed+1;
    bad_policy=policy; bad_policy.gates.sensitivity=1.1e-4;
    local_reject(@()r1_bie_continuation_io('policy',bad_policy)); passed=passed+1;
    bad_policy=policy; bad_policy.grid_monotonicity_mode='ignore';
    local_reject(@()r1_bie_continuation_io('policy',bad_policy)); passed=passed+1;
    % v2 state fixtures use nonphysical files; they cannot be loaded as fields.
    fixture=fullfile(out,'state_fixture'); mkdir(fixture); sig='nonphysical-v2-fixture';
    assert(strcmp(r1_bie_continuation_io('c4_state',fixture,sig),'SOLVE_C4_ONCE')); passed=passed+1;
    payload={'data.mat','points.csv','config.json','timings.json','solve_invocation.json'};
    for j=1:numel(payload), r1_bie_continuation_io('write',fullfile(fixture,payload{j}),struct('nonphysical',true)); end
    local_reject(@()r1_bie_continuation_io('c4_state',fixture,sig)); passed=passed+1;
    r1_bie_continuation_io('seal',fixture,'solved.json',sig,payload,'SOLVED_UNAUDITED');
    assert(strcmp(r1_bie_continuation_io('c4_state',fixture,sig),'AUDIT_ONLY')); passed=passed+1;
    local_reject(@()r1_bie_continuation_io('c4_state',fixture,[sig 'changed'])); passed=passed+1;
    r1_bie_continuation_io('write',fullfile(fixture,'audit.json'),struct('status','PASS','failed_checks',0));
    r1_bie_continuation_io('write',fullfile(fixture,'gates_v2.json'),struct('status','PASS'));
    r1_bie_continuation_io('write',fullfile(fixture,'audit_reconstruction.mat'),struct('nonphysical',true));
    r1_bie_continuation_io('seal',fixture,'committed.json',sig, ...
        [payload {'solved.json','audit.json','audit_reconstruction.mat','gates_v2.json'}],'AUDITED_PASS');
    assert(strcmp(r1_bie_continuation_io('c4_state',fixture,sig),'SKIP_AUDITED')); passed=passed+1;
    local_reject(@()r1_bie_continuation_io('seal',fixture,'committed.json',sig,payload,'AUDITED_PASS')); passed=passed+1;
    local_reject(@()r1_bie_continuation_io('verify_seal',fixture,'committed.json',sig,{'missing.json'})); passed=passed+1;
    r1_bie_continuation_io('write',fullfile(fixture,'data.mat'),struct('tampered',true));
    local_reject(@()r1_bie_continuation_io('c4_state',fixture,sig)); passed=passed+1;
    % Source SHA check before and after migration; no original writes permitted.
    origin_before=r1_bie_continuation_io('verify_origin',root);
    migration=migrate_r1_bie_checkpoints(root); passed=passed+1;
    repeated=migrate_r1_bie_checkpoints(root); assert(strcmp(repeated.status,'SKIP_MIGRATION_AUDITED')); passed=passed+1;
    origin_after=r1_bie_continuation_io('verify_origin',root); assert(isequal(origin_before,origin_after)); passed=passed+1;
    chain=r1_bie_continuation_io('verify_migration',root);
    dry=run_r1_bie_c4_continuation(root); assert(strcmp(dry.status,'DRY_RUN_PASS') && dry.new_acoustic_solves==0); passed=passed+1;
    passed=passed+local_final_fixture(root,out,policy);
    source=fileread(fullfile(code,'run_r1_bie_c4_continuation.m'));
    assert(isscalar(regexp(source,'boundary=solve_helmholtz_bie_halfplane_vertical\(cfg\)','match')) && ...
        ~contains(source,'for stage=0:4')); passed=passed+1;
    assert(~isfile(fullfile(out,'forbidden_solve_calls.txt')) && ~isfile(fullfile(root,'jobs','bie_window','data.mat'))); passed=passed+1;
    files={'run_r1_bie_c4_continuation.m','migrate_r1_bie_checkpoints.m','r1_bie_continuation_io.m', ...
        'r1_bie_continuation_gates.m','audit_r1_bie_continuation_final.m','validate_r1_bie_c4_preparation.m'};
    analyzer=struct('file',{},'issues',{}); messages=0;
    for j=1:numel(files)
        issues=checkcode(fullfile(code,files{j}),'-id'); analyzer(j)=struct('file',files{j},'issues',issues);
        messages=messages+numel(issues);
    end
    r1_bie_continuation_io('write',fullfile(out,'code_analyzer.json'),analyzer);
    assert(messages==0,'r1v2:Analyzer','Review and fix all new MATLAB Code Analyzer messages before readiness'); passed=passed+1;
    receipt=struct('status','PASS','failed_checks',0,'passed_checks',passed,'elapsed_s',toc(timer), ...
        'new_acoustic_solves',0,'C3_saved_density_audit',migration,'MATLAB_release',version,'max_threads',maxNumCompThreads, ...
        'manifest_sha256',chain.manifest_sha256,'policy_sha256',chain.policy_sha256, ...
        'C4_state','SOLVE_C4_ONCE','C4_not_run',true,'code_analyzer_messages',messages, ...
        'unrun_checks',{{'Real C4 solve and final five-stage numeric/PNG audit','Remote Linux launcher lifecycle'}});
    r1_bie_continuation_io('write',fullfile(out,'receipt.json'),receipt);
    dest=fullfile(root,'preparation_validation'); assert(~isfolder(dest)); mkdir(dest);
    copyfile(fullfile(out,'receipt.json'),fullfile(dest,'receipt.json')); copyfile(fullfile(out,'code_analyzer.json'),fullfile(dest,'code_analyzer.json'));
    r1_bie_continuation_io('write',fullfile(root,'local_validation.json'),receipt);
    r1_bie_continuation_io('seal',root,'local_validation_commit.json', ...
        [chain.predecessor ':VALIDATED:' chain.manifest_sha256 ':' chain.policy_sha256], ...
        {'local_validation.json','preparation_validation/receipt.json','preparation_validation/code_analyzer.json'},'LOCAL_ACCEPTANCE_PASS');
    diary('off'); clear diary_guard path_guard; disp(jsonencode(receipt));
catch ME
    failures(end+1)=struct('name','preparation','identifier',ME.identifier,'message',ME.message);
    receipt=struct('status','FAIL','failed_checks',numel(failures),'passed_checks',passed, ...
        'failures',failures,'elapsed_s',toc(timer),'new_acoustic_solves',0);
    r1_bie_continuation_io('write',fullfile(out,'receipt.json'),receipt);
    fprintf(2,'%s\n',getReport(ME,'extended','hyperlinks','off')); diary('off'); rethrow(ME);
end
end
function count=local_final_fixture(root,out,policy)
% Nonphysical window fixture tests final policy/metrics, never C4 convergence.
pkg=fullfile(out,'FINAL_TEST_ONLY_PACKAGE'); copyfile(root,pkg);
work=fullfile(out,'FINAL_TEST_ONLY_WORK'); mkdir(work);
r1_bie_continuation_io('write',fullfile(work,'SYNTHETIC_NOT_A_C4_RESULT.json'), ...
    struct('synthetic_window_field','C1 field reused only to test postprocessing','scientific_result',false));
copyfile(fullfile(root,'common.mat'),fullfile(work,'common.mat'));
copyfile(fullfile(root,'input'),fullfile(work,'input'));
ids={'bie_coarse','bie_main','bie_self','bie_space','bie_window'}; fields=cell(5,1);
for j=1:5
    from=ids{min(j,4)}; if j==5, from='bie_main'; end
    d=fullfile(work,'jobs',ids{j}); mkdir(d);
    copyfile(fullfile(root,'jobs',from,'data.mat'),fullfile(d,'data.mat'));
    b=load(fullfile(d,'data.mat'),'field'); fields{j}=b.field;
end
q=load(fullfile(root,'common.mat'));
gate=r1_bie_continuation_gates(fields,q.x,q.z,policy);
gate.policy_sha256=r1_bie_continuation_io('sha',fullfile(root,'policy.json'));
d=fullfile(pkg,'jobs','bie_window'); mkdir(d);
r1_bie_continuation_io('write',fullfile(d,'gates_v2.json'),gate);
run_piston_pm_r1_job(work,'analyze');
r=audit_r1_bie_continuation_final(work,pkg);
assert(strcmp(r.status,'PASS') && r.grid_nonmonotonic && r.reference_stage==3); count=1;
% A mismatched stage/final policy cannot be accepted.
gate.policy_sha256=repmat('0',1,64); r1_bie_continuation_io('write',fullfile(d,'gates_v2.json'),gate);
local_reject(@()audit_r1_bie_continuation_final(work,pkg)); count=count+1;
gate.policy_sha256=r1_bie_continuation_io('sha',fullfile(root,'policy.json'));
r1_bie_continuation_io('write',fullfile(d,'gates_v2.json'),gate);
M=readtable(fullfile(work,'model_metrics.csv')); M.complex_relative_L2(1)=M.complex_relative_L2(1)+1e-6;
writetable(M,fullfile(work,'model_metrics.csv'));
local_reject(@()audit_r1_bie_continuation_final(work,pkg)); count=count+1;
% Window sensitivity above 1e-4 must still hard-stop the final audit.
field=fields{2}*1.001; save(fullfile(work,'jobs','bie_window','data.mat'),'field','-v7');
local_reject(@()audit_r1_bie_continuation_final(work,pkg)); count=count+1;
end
function local_reject(callback)
failed=false;
try
    callback();
catch
    failed=true;
end
assert(failed,'r1v2:ExpectedRejection','Invalid evidence unexpectedly accepted');
end
