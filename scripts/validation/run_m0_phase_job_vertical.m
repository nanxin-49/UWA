function receipt=run_m0_phase_job_vertical(stage,run_dir)
%RUN_M0_PHASE_JOB_VERTICAL Serial, bounded M0 work in the existing session.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root); setup_vertical_project();
assert(java.io.File(run_dir).isAbsolute(),'Absolute run directory required.');
assert(startsWith(char(java.io.File(run_dir).getCanonicalPath()), ...
    fullfile(root,'results','validation','m0_phase_contract')),'M0 result path required.');
if ~isfolder(run_dir), mkdir(run_dir); end
budget_file=fullfile(run_dir,'budget_state.json');
if isfile(budget_file), budget=jsondecode(fileread(budget_file));
else, budget=struct('limit_s',600,'used_s',15,'preflight_conservative_s',15); end
assert(budget.used_s<600,'M0 MATLAB budget exhausted.');
job_dir=fullfile(run_dir,stage); assert(~isfolder(job_dir),'Never overwrite a job.'); mkdir(job_dir);
timer=tic; old_threads=maxNumCompThreads(1);
restore_threads=onCleanup(@()maxNumCompThreads(old_threads)); %#ok<NASGU>
diary(fullfile(job_dir,'stdout.log')); stop_diary=onCleanup(@()diary('off')); %#ok<NASGU>
receipt=struct('run_id',string(java.io.File(run_dir).getName()),'job_id',stage, ...
    'status','FAILED','elapsed_s',0,'failed_checks',1,'result_dir',job_dir, ...
    'execution_mode','existing_local_matlab_mcp','matlab_release',version, ...
    'exit_code',1,'exit_code_kind','MATLAB exception status; shared process is not exited', ...
    'stdout_log',fullfile(job_dir,'stdout.log'),'stderr_log',fullfile(job_dir,'stderr.log'), ...
    'checkpoint_path',fullfile(job_dir,'checkpoint.json'),'threads',maxNumCompThreads);
payload=struct();
try
    [~,ram_before]=memory;
    receipt.available_physical_before_bytes=ram_before.PhysicalMemory.Available;
    receipt.estimated_peak_additional_bytes=200*2^20;
    assert(receipt.estimated_peak_additional_bytes<=0.7*ram_before.PhysicalMemory.Available, ...
        'M0 RAM gate: preserve 30 percent of currently available physical RAM.');
    payload=validate_m0_static_phase_vertical(stage,run_dir);
    assert(payload.pass,'M0 stage checks failed.');
    assert(toc(timer)<=180,'M0 individual job exceeded 180 seconds.');
    receipt.status='COMPLETE'; receipt.failed_checks=0; receipt.exit_code=0;
    local_json(receipt.stderr_log,struct('error','none'));
catch exception
    receipt.error_identifier=exception.identifier; receipt.error_message=exception.message;
    fid=fopen(receipt.stderr_log,'w'); fprintf(fid,'%s\n',getReport(exception,'extended','hyperlinks','off')); fclose(fid);
    fprintf('M0 FAILED %s: %s\n',exception.identifier,exception.message);
end
save(fullfile(job_dir,'payload.mat'),'payload','-v7');
receipt.completed_steps={stage};
receipt.config=struct('phase_contract_version','2.0.0','serial',true, ...
    'single_job_limit_s',180,'total_budget_s',600,'reserved_available_RAM_fraction',.3);
manifest_path=fullfile(run_dir,'code_manifest_review.json');
if ~isfile(manifest_path), manifest_path=fullfile(run_dir,'code_manifest.json'); end
manifest=jsondecode(fileread(manifest_path));
for j=1:numel(manifest), manifest(j).sha256=local_sha(manifest(j).path); end
receipt.code_versions_sha256=manifest;
receipt.input_versions_sha256=struct('path',fullfile(run_dir,'checkpoint_initial.json'), ...
    'sha256',local_sha(fullfile(run_dir,'checkpoint_initial.json')));
receipt.elapsed_s=toc(timer); budget.used_s=budget.used_s+receipt.elapsed_s;
receipt.cumulative_elapsed_s=budget.used_s;
[~,ram_after]=memory; receipt.available_physical_after_bytes=ram_after.PhysicalMemory.Available;
receipt.resource_gate_pass=ram_after.PhysicalMemory.Available>=0.3*receipt.available_physical_before_bytes;
if ~receipt.resource_gate_pass || budget.used_s>600
    receipt.status='FAILED'; receipt.failed_checks=receipt.failed_checks+1; receipt.exit_code=1;
end
receipt.next_allowed_step='Next authorized M0 job only; no M1 solve.';
receipt.audit_status_and_path='Independent audit required before M0 closure.';
if isfield(payload,'checks')
    receipt.failed_or_unrun_checks=payload.checks(~[payload.checks.pass]);
    receipt.failed_checks=max(receipt.failed_checks,numel(receipt.failed_or_unrun_checks));
else
    receipt.failed_or_unrun_checks=struct('name','job did not return validation checks');
end
receipt.resource_usage=struct('elapsed_s',receipt.elapsed_s,'cumulative_s',budget.used_s, ...
    'budget_s',600,'job_limit_s',180,'threads',receipt.threads, ...
    'available_before_bytes',receipt.available_physical_before_bytes, ...
    'available_after_bytes',receipt.available_physical_after_bytes);
clear stop_diary
files=dir(job_dir); artifacts=struct('path',{},'sha256',{});
for j=1:numel(files)
    if files(j).isdir, continue; end
    ap=fullfile(job_dir,files(j).name);
    artifacts(end+1)=struct('path',ap,'sha256',local_sha(ap)); %#ok<AGROW>
end
receipt.artifacts_sha256=artifacts;
local_json(budget_file,budget); local_json(receipt.checkpoint_path,receipt);
short=rmfield(receipt,{'code_versions_sha256','artifacts_sha256'});
fprintf('M0_RECEIPT %s\n',jsonencode(short));
end
function value=local_sha(path)
md=java.security.MessageDigest.getInstance('SHA-256');
fid=fopen(path,'rb'); c=onCleanup(@()fclose(fid)); %#ok<NASGU>
bytes=fread(fid,Inf,'*uint8');
if ~isempty(bytes), md.update(typecast(bytes,'int8')); end
value=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
function local_json(path,value)
fid=fopen(path,'w','n','UTF-8'); assert(fid>0); c=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(value,PrettyPrint=true));
end
