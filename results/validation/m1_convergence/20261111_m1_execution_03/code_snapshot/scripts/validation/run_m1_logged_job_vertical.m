function receipt=run_m1_logged_job_vertical(job_id,run_dir)
%RUN_M1_LOGGED_JOB_VERTICAL Persist evalc output despite outer MCP timeouts.
% Executes one job only, synchronously in the existing MATLAB process.
if ~isfolder(run_dir), mkdir(run_dir); end
dispatch_path=fullfile(run_dir,[job_id '_dispatch.json']);
assert(~isfile(dispatch_path),'Existing dispatch: inspect original job; never duplicate.');
t=tic;
dispatch=struct('job_id',job_id,'status','DISPATCHED_CAPTURE_PENDING', ...
    'matlab_pid',feature('getpid'),'entrypoint',mfilename, ...
    'wrapper_sha256',sha([mfilename('fullpath') '.m']));
jsonfile(dispatch_path,dispatch);
% A literal, reviewed nested capture is needed: MCP itself captures command
% output, and diary alone produced empty logs in the first completed pilot.
captured=evalc('receipt=run_m1_convergence_job_vertical(job_id,run_dir);');
job_dir=fullfile(run_dir,job_id);
file=fullfile(job_dir,'stdout_captured.log');
fid=fopen(file,'w','n','UTF-8'); assert(fid>0); fprintf(fid,'%s',captured); fclose(fid);
assert(~isempty(captured),'Explicit stdout capture unexpectedly empty.');
copyfile(receipt.checkpoint_path,fullfile(job_dir,'checkpoint_solver.json'));
receipt.diary_log=receipt.stdout_log; receipt.stdout_log=file;
receipt.stdout_capture='COMPLETE_EXPLICIT_EVALC';
receipt.wrapper_sha256=dispatch.wrapper_sha256;
receipt.wrapper_elapsed_s=toc(t);
receipt.artifacts_sha256(end+1)=struct('path','stdout_captured.log','sha256',sha(file));
receipt.artifacts_sha256(end+1)=struct('path','checkpoint_solver.json','sha256',sha(fullfile(job_dir,'checkpoint_solver.json')));
jsonfile(receipt.checkpoint_path,receipt);
dispatch.status=receipt.status; dispatch.elapsed_s=receipt.wrapper_elapsed_s;
dispatch.checkpoint=receipt.checkpoint_path; dispatch.stdout_log=file; dispatch.stdout_sha256=sha(file);
jsonfile(dispatch_path,dispatch);
fprintf('M1_LOGGED_RECEIPT job=%s status=%s elapsed=%.3fs capture_bytes=%d checkpoint=%s\n', ...
    job_id,receipt.status,receipt.elapsed_s,numel(unicode2native(captured,'UTF-8')),receipt.checkpoint_path);
end
function value=sha(file)
md=java.security.MessageDigest.getInstance('SHA-256'); fid=fopen(file,'rb'); c=onCleanup(@()fclose(fid)); %#ok<NASGU>
while ~feof(fid), bytes=fread(fid,2^20,'*uint8'); if ~isempty(bytes), md.update(typecast(bytes,'int8')); end, end
value=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
function jsonfile(file,s)
fid=fopen(file,'w','n','UTF-8'); fprintf(fid,'%s',jsonencode(s,PrettyPrint=true)); fclose(fid);
end
