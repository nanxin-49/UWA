function receipt=validate_r2_pm_preparation(package_dir,output_dir)
%VALIDATE_R2_PM_PREPARATION Solver-free runtime acceptance; synthetic state only.
arguments
    package_dir (1,:) char
    output_dir (1,:) char
end
assert(~isfolder(output_dir)); mkdir(output_dir); sentinel=fullfile(output_dir,'sentinels'); mkdir(sentinel);
diary(fullfile(output_dir,'stdout.log')); diary_guard=onCleanup(@()diary('off'));
original_path=path; guard=onCleanup(@()path(original_path));
blocked={'solve_helmholtz_bie_halfplane_vertical','solve_helmholtz_bie_halfplane_cached_vertical', ...
    'baffled_strip_piston_field_vertical','piston_independent_angular_spectrum_vertical', ...
    'closeout_s3a_pm_fullfk','run_pe_1d_surface_reflection_validation','besselh'};
for j=1:numel(blocked)
    fid=fopen(fullfile(sentinel,[blocked{j} '.m']),'w'); assert(fid>=0);
    fprintf(fid,'function varargout=%s(varargin)\nerror(''r2:AcousticForbidden'',''Forbidden during preparation'');\nend\n',blocked{j}); fclose(fid);
end
addpath(fullfile(package_dir,'code'),'-begin'); addpath(sentinel,'-begin');
try
    for j=1:numel(blocked), assert(startsWith(which(blocked{j}),sentinel)); end
    assert(~isfolder(fullfile(package_dir,'jobs')) || isfile(fullfile(package_dir,'manifest.json')));
    r=run_r2_pm_seed(package_dir,'dry-run'); assert(r.new_acoustic_solves==0);
    fixture=fullfile(output_dir,'synthetic_checkpoint'); mkdir(fixture); sig='synthetic:seed:manifest:predecessor';
    assert(strcmp(r2_pm_io('resume_state',fixture,sig),'SOLVE'));
    files={'data.mat','points.csv','config.json','timings.json'};
    for j=1:numel(files), local_text(fullfile(fixture,files{j}),'synthetic; never an acoustic result'); end
    local_reject(@()r2_pm_io('resume_state',fixture,sig));
    r2_pm_io('seal',fixture,'solved.json',sig,files,'SOLVED_UNAUDITED');
    assert(strcmp(r2_pm_io('resume_state',fixture,sig),'AUDIT_ONLY'));
    r2_pm_io('write',fullfile(fixture,'audit.json'),struct('status','PASS','failed_checks',0));
    r2_pm_io('write',fullfile(fixture,'gates.json'),struct('status','PASS'));
    r2_pm_io('seal',fixture,'committed.json',sig,[files {'solved.json','audit.json','gates.json'}],'AUDITED_PASS');
    assert(strcmp(r2_pm_io('resume_state',fixture,sig),'SKIP_AUDITED'));
    local_reject(@()r2_pm_io('resume_state',fixture,[sig ':changed_predecessor']));
    original=fileread(fullfile(fixture,'data.mat')); local_text(fullfile(fixture,'data.mat'),'tampered');
    local_reject(@()r2_pm_io('resume_state',fixture,sig)); local_text(fullfile(fixture,'data.mat'),original);
    policy=r2_pm_io('read',fullfile(package_dir,'policy.json'));
    x=(-205:205)'*.1953125; z=3*ones(411,1); one=ones(411,1);
    F={one+1e-7,one,one+1e-7,one+1.1e-7,one};
    g=r1_bie_continuation_gates(F,x,z,policy);
    assert(strcmp(g.status,'PASS') && g.grid_diagnostic.grid_nonmonotonic && strcmp(g.legacy_policy_status,'FAIL'));
    F{5}=one+.001; g=r1_bie_continuation_gates(F,x,z,policy); assert(strcmp(g.status,'FAIL'));
    policy.gates.sensitivity=1e-3;
    local_reject(@()r1_bie_continuation_gates(F,x,z,policy));
    files=dir(fullfile(package_dir,'code','*.m')); analyzer=struct; warnings=0;
    for j=1:numel(files)
        issues=checkcode(fullfile(files(j).folder,files(j).name),'-id');
        analyzer.(erase(files(j).name,'.m'))=issues; warnings=warnings+numel(issues);
    end
    save(fullfile(output_dir,'checks.mat'),'r','analyzer','blocked','-v7');
    receipt=struct('status','SOLVER_FREE_DRY_RUN_PASS','failed_checks',0,'new_acoustic_solves',0, ...
        'SSH_connections',0,'MATLAB_version',version,'manifest_sha256',r.manifest_sha256, ...
        'blocked_acoustic_entrypoints',numel(blocked),'Code_Analyzer_messages',warnings, ...
        'synthetic_checks',{{'fresh','orphan refusal','audit-only','skip audited','predecessor drift', ...
        'artifact drift','v2 nonmonotonic diagnostic','sensitivity failure','policy tamper refusal'}});
    r2_pm_io('write',fullfile(output_dir,'receipt.json'),receipt); disp(jsonencode(receipt));
catch ME
    r2_pm_io('write',fullfile(output_dir,'failure.json'),struct('status','FAIL','error',ME.message,'identifier',ME.identifier));
    fprintf(2,'%s\n',getReport(ME,'extended','hyperlinks','off')); rethrow(ME);
end
end
function local_reject(fn)
rejected=false; try, fn(); catch, rejected=true; end
assert(rejected,'r2:RejectionTest','Invalid synthetic state accepted');
end
function local_text(path,value)
fid=fopen(path,'w'); assert(fid>=0); c=onCleanup(@()fclose(fid)); fprintf(fid,'%s',value);
end
