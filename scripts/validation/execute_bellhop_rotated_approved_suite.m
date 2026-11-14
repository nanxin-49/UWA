function receipt = execute_bellhop_rotated_approved_suite(suite_file,executable,authorization_file,output_dir,controller_dir,options)
%EXECUTE_BELLHOP_ROTATED_APPROVED_SUITE Bounded MCP stdout/error persistence.
% The frozen suite enforces its shared solver/audit budget and resource gates.
% This wrapper starts no process when the approved RAM gate is not satisfied.
arguments
    suite_file (1,:) char
    executable (1,:) char
    authorization_file (1,:) char
    output_dir (1,:) char
    controller_dir (1,:) char
    options.reuse_dir (1,:) char=''
end
assert(~isfolder(controller_dir)); mkdir(controller_dir);
a=jsondecode(fileread(authorization_file));
assert(a.user_approved && strcmp(a.scope,'rotated_fig12b_full_local'));
assert(strcmpi(black_file_hash(suite_file),a.suite_sha256));
s=load(suite_file,'suite'); drift=0;
for j=1:numel(s.suite.code_versions)
    drift=drift+~strcmpi(black_file_hash(s.suite.code_versions(j).path),s.suite.code_versions(j).sha256);
end
assert(drift==0,'Frozen suite execution source drift.');
[~,mem]=memory;
receipt=struct('status','NOT_STARTED_RESOURCE_GATE','controller_dir',controller_dir, ...
    'authorization_file',authorization_file,'authorization_sha256',black_file_hash(authorization_file), ...
    'suite_file',suite_file,'suite_sha256',a.suite_sha256,'source_drift_count',drift, ...
    'available_bytes',mem.PhysicalMemory.Available,'required_bytes',a.minimum_available_bytes, ...
    'output_dir',output_dir,'solver_dispatch_attempted',false,'checked_utc',char(datetime('now','TimeZone','UTC')));
if mem.PhysicalMemory.Available<a.minimum_available_bytes
    write(fullfile(controller_dir,'checkpoint.json'),jsonencode(receipt,PrettyPrint=true)); disp(jsonencode(receipt)); return
end
receipt.status='DISPATCHING_APPROVED_SUITE'; receipt.solver_dispatch_attempted=true;
write(fullfile(controller_dir,'checkpoint.json'),jsonencode(receipt,PrettyPrint=true));
t=tic; failure=''; suite_result=[];
invoke=@()run_bellhop_rotated_field_suite(suite_file,executable,output_dir,authorization_file,reuse_dir=options.reuse_dir); %#ok<NASGU> Called by fixed evalc expression.
stdout=evalc('try; suite_result=invoke(); catch exception; failure=getReport(exception,"extended","hyperlinks","off"); fprintf("%s\n",failure); end');
write(fullfile(controller_dir,'stdout.log'),stdout); write(fullfile(controller_dir,'stderr.log'),failure);
receipt.wrapper_elapsed_s=toc(t);
receipt.suite_checkpoint=fullfile(output_dir,'checkpoint.json');
if isempty(failure)
    receipt.status=suite_result.status; receipt.cumulative_solver_audit_s=suite_result.cumulative_elapsed_s;
else
    receipt.status='FAILED_INSPECT_DURABLE_CHECKPOINT_BEFORE_ANY_RETRY'; receipt.error=failure;
    if isfile(receipt.suite_checkpoint)
        terminal=jsondecode(fileread(receipt.suite_checkpoint));
        receipt.suite_status=terminal.status;
        if isfield(terminal,'cumulative_elapsed_s'), receipt.cumulative_solver_audit_s=terminal.cumulative_elapsed_s; end
    end
end
write(fullfile(controller_dir,'checkpoint.json'),jsonencode(receipt,PrettyPrint=true));
small=rmfield(receipt,intersect(fieldnames(receipt),{'error'})); disp(jsonencode(small));
end
function write(path,content)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); cleaner=onCleanup(@()fclose(fid)); fprintf(fid,'%s',content);
end