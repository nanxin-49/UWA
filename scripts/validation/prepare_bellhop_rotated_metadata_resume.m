function receipt = prepare_bellhop_rotated_metadata_resume(old_suite_file,latest_run_dir,validation_dir,disabled_dir,output_dir)
%PREPARE_BELLHOP_ROTATED_METADATA_RESUME Explicit cache requalification.
% Rebuild every cached input with new I/O code. Byte identity plus unchanged
% executable/pressure kernels proves reuse. Preserve original source versions
% and checkpoints; new receipts are compatibility qualifications, not solves.
arguments
    old_suite_file (1,:) char
    latest_run_dir (1,:) char
    validation_dir (1,:) char
    disabled_dir (1,:) char
    output_dir (1,:) char
end
assert(~isfolder(output_dir)); mkdir(output_dir); t=tic;
old=load(old_suite_file,'suite'); suite=old.suite;
v=jsondecode(fileread(fullfile(validation_dir,'checkpoint.json')));
d=jsondecode(fileread(fullfile(disabled_dir,'checkpoint.json')));
assert(v.failed_checks==0 && d.failed_checks==0);
ledger=jsondecode(fileread(fullfile(latest_run_dir,'checkpoint.json')));
assert(strcmp(ledger.status,'FAILED'));
source_control=fullfile(latest_run_dir,suite.labels{1});
parent=jsondecode(fileread(fullfile(source_control,'checkpoint.json')));
assert(parent.completed_jobs==31,'Expected7 direct plus up/down1..12; inspect any changed boundary.');
allowed={'write_bellhop_rotated_field_options','write_bellhop_rotated_multifield_q','run_bellhop_rotated_field_component'};
changed={};
for j=1:numel(suite.code_versions)
    path=suite.code_versions(j).path; current=black_file_hash(path);
    if ~strcmpi(current,suite.code_versions(j).sha256)
        [~,name]=fileparts(path); assert(ismember(name,allowed),'Unexpected source change.'); changed{end+1}=name; %#ok<AGROW>
    end
    suite.code_versions(j).sha256=current;
end
assert(numel(changed)==3);
helper=which('repair_bellhop_rotated_range_metadata');
suite.code_versions(end+1)=struct('path',helper,'sha256',black_file_hash(helper));
suite.source_approved_suite_sha256=black_file_hash(old_suite_file);
suite.status='METADATA_IO_FIX_SAME_PHYSICS_PLANS_BUDGET';
suite.software_fix='ENV float-label collision repair only; exact physical receiver coordinates remain byte-identical';
suite.metadata_validation=validation_dir; suite.disabled_validation=disabled_dir;
for j=1:numel(suite.plan_files), assert(strcmpi(black_file_hash(suite.plan_files{j}),suite.plan_sha256{j})); end
new_suite_file=fullfile(output_dir,'suite.mat'); save(new_suite_file,'suite'); new_suite_sha=black_file_hash(new_suite_file);
authorization=jsondecode(fileread(ledger.approval_file));
assert(authorization.user_approved && strcmpi(authorization.suite_sha256,suite.source_approved_suite_sha256));
authorization.basis_authorization_file=ledger.approval_file;
authorization.basis_authorization_sha256=black_file_hash(ledger.approval_file);
authorization.suite_sha256=new_suite_sha;
authorization.max_wall_s=authorization.max_wall_s-30;
authorization.conservative_unmeasured_preparation_reserve_s=30;
authorization.software_fix_scope='Goal-authorized local metadata I/O repair; same four physical plans, same binary, same RAM/disk gates.30s reserved from original upper bound; no scope/budget increase.';
authorization.allowed_plan_sha256=suite.plan_sha256;
new_auth_file=fullfile(output_dir,'authorization_software_fix.json');
write(new_auth_file,jsonencode(authorization,PrettyPrint=true));
reuse=fullfile(output_dir,'reuse'); mkdir(reuse); new_control=fullfile(reuse,suite.labels{1}); mkdir(new_control);
plan_data=load(suite.plan_files{1},'plan'); jobs=plan_data.plan.jobs;
checks=table(strings(0,1),false(0,1),zeros(0,1),zeros(0,1),'VariableNames',{'check','pass','value','limit'});
source_names={'run_bellhop_rotated_field_component','run_bellhop_rotated_local_case', ...
    'write_bellhop_rotated_multifield_q','write_bellhop_rotated_munk_q','write_bellhop_rotated_horizontal_q', ...
    'write_bellhop_rotated_field_options','read_bellhop_rotated_field','munk_profile_sound_speed','repair_bellhop_rotated_range_metadata'};
current_versions=struct('path',{},'sha256',{});
for j=1:numel(source_names)
    path=which(source_names{j}); current_versions(j)=struct('path',path,'sha256',black_file_hash(path));
end
for j=1:parent.completed_jobs
    job=jobs(j); cfg=job.cfg; prior_path=fullfile(source_control,job.tag,'checkpoint.json');
    prior=jsondecode(fileread(prior_path)); assert(strcmp(prior.status,'COMPLETE'));
    assert(strcmpi(prior.executable_sha256,suite.executable_sha256));
    assert(strcmp(jsonencode(prior.cfg),jsonencode(cfg)));
    for b=1:numel(prior.code_versions)
        path=prior.code_versions(b).path; [~,name]=fileparts(path);
        if ~ismember(name,allowed), assert(strcmpi(black_file_hash(path),prior.code_versions(b).sha256),'Non-I/O pressure source drift.'); end
    end
    for b=1:numel(prior.bindings), assert(strcmpi(black_file_hash(prior.bindings(b).path),prior.bindings(b).sha256),'Original cache SHA mismatch.'); end
    proof_dir=fullfile(output_dir,'input_proofs',job.tag); mkdir(proof_dir); stem=fullfile(proof_dir,'case');
    if isfield(cfg,'reflection_order')
        write_bellhop_rotated_multifield_q(stem,cfg);
    elseif strcmp(cfg.chart_kind,'horizontal')
        write_bellhop_rotated_horizontal_q(stem,cfg,cfg.horizontal_side);
    else
        meta=write_bellhop_rotated_munk_q(stem,cfg);
        fid=fopen([stem '.iwcov'],'w'); fprintf(fid,'2\n0\n0\n0\n0\n'); fclose(fid);
        write_bellhop_rotated_field_options(stem,cfg.receiver_range_m,meta.q_receiver_m,point3d=true,fast_q=true);
    end
    for ext={'.env','.ssp','.iwcov','.rotq','.r25d'}
        if ~isfile([stem ext{1}]), continue; end
        same=strcmpi(black_file_hash([stem ext{1}]),black_file_hash([prior.raw_root ext{1}]));
        checks(end+1,:)={string([job.tag ext{1}]),same,double(~same),0}; %#ok<AGROW>
        assert(same,'Cache input changed; do not migrate this contribution.');
    end
    migrated=prior;
    migrated.original_execution_code_versions=prior.code_versions;
    migrated.original_checkpoint_path=prior_path; migrated.original_checkpoint_sha256=black_file_hash(prior_path);
    migrated.code_versions=current_versions;
    migrated.checkpoint_role='BYTE_IDENTICAL_INPUT_AND_BINARY_COMPATIBILITY_QUALIFICATION_NOT_A_NEW_SOLVE';
    migrated.qualification_input_root=stem;
    migrated.qualification_preparer_sha256=black_file_hash([mfilename('fullpath') '.m']);
    migrated.metadata_helper_sha256=black_file_hash(helper);
    target=fullfile(new_control,job.tag); mkdir(target);
    write(fullfile(target,'checkpoint.json'),jsonencode(migrated,PrettyPrint=true));
end
parent.compatibility_qualification='Input bytes identical; original pressure binary and non-I/O sources unchanged; original source versions preserved per component';
write(fullfile(new_control,'checkpoint.json'),jsonencode(parent,PrettyPrint=true));
qualification_s=toc(t);
ledger.source_ledger_path=fullfile(latest_run_dir,'checkpoint.json');
ledger.source_ledger_sha256=black_file_hash(ledger.source_ledger_path);
ledger.source_approved_suite_sha256=ledger.suite_sha256; ledger.suite_sha256=new_suite_sha;
ledger.status='METADATA_COMPATIBILITY_RESUME_ANCHOR_NOT_A_NUMERICAL_RUN'; ledger.result_dir=reuse;
ledger.cumulative_elapsed_s=ledger.cumulative_elapsed_s+v.elapsed_s+d.elapsed_s+qualification_s;
ledger.extra_measured_metadata_validation_s=v.elapsed_s+d.elapsed_s;
ledger.extra_measured_cache_qualification_s=qualification_s;
ledger.disk_roots=unique([cellstr(string(ledger.disk_roots(:)));{output_dir}]);
write(fullfile(reuse,'checkpoint.json'),jsonencode(ledger,PrettyPrint=true));
writetable(checks,fullfile(output_dir,'checks.csv'));
receipt=struct('status','CACHE_INPUT_COMPATIBILITY_QUALIFIED_NO_NEW_FULL_SOLVE','failed_checks',nnz(~checks.pass), ...
    'checks_count',height(checks),'cached_components',parent.completed_jobs,'qualification_elapsed_s',qualification_s, ...
    'new_suite_file',new_suite_file,'new_suite_sha256',new_suite_sha,'new_authorization_file',new_auth_file, ...
    'resume_dir',reuse,'cumulative_measured_before_resume_s',ledger.cumulative_elapsed_s, ...
    'remaining_solver_audit_budget_s',authorization.max_wall_s-ledger.cumulative_elapsed_s, ...
    'unmeasured_preparation_reserved_s',30,'unchanged_physical_plan_sha256',{suite.plan_sha256}, ...
    'unchanged_executable_sha256',suite.executable_sha256,'changed_source_roles',{changed},'output_dir',output_dir);
write(fullfile(output_dir,'checkpoint.json'),jsonencode(receipt,PrettyPrint=true)); disp(jsonencode(receipt));
end
function write(path,content)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); cleaner=onCleanup(@()fclose(fid)); fprintf(fid,'%s',content);
end