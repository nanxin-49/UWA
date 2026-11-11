function receipt=run_m1_convergence_stage_vertical(run_dir,stage)
%RUN_M1_CONVERGENCE_STAGE_VERTICAL Serial dispatch with durable boundaries.
% P4 is dispatched separately to verify the pilot and capture mechanism.
assert(ismember(stage,{'references','remaining_pe','sponge_remaining'}));
checkpoint=fullfile(run_dir,[stage '_dispatch_checkpoint.json']);
assert(~isfile(checkpoint),'Stage already dispatched: inspect checkpoint, never relaunch blindly.');
p4=jsondecode(fileread(fullfile(run_dir,'P4','checkpoint.json')));
assert(strcmp(p4.status,'COMPLETE') && p4.elapsed_s<=900 && isfield(p4,'stdout_capture'));
switch stage
    case 'references'
        jobs={'A-C4','A-C6','A-C8','A-W4','A-W6','A-W8','A-G4','A-G6','A-G8','A-R4'};
    case 'remaining_pe'
        reference_path=fullfile(run_dir,'reference_audit_v2','checkpoint.json');
        if ~isfile(reference_path), reference_path=fullfile(run_dir,'reference_audit','checkpoint.json'); end
        refs=jsondecode(fileread(reference_path));
        assert(strcmp(refs.status,'AUDIT_PASS'),'Reference convergence/edge qualification failed.');
        jobs={'P6','P8','Z4','Z8'};
    case 'sponge_remaining'
        sp4=jsondecode(fileread(fullfile(run_dir,'SP4','checkpoint.json')));
        assert(strcmp(sp4.status,'COMPLETE') && sp4.elapsed_s<=900 && isfield(sp4,'stdout_capture'));
        jobs={'SP6','SP8'};
end
receipt=struct('status','RUNNING','stage',stage,'matlab_pid',feature('getpid'), ...
    'completed_jobs',{{}},'next_job','','stopped_reason','', ...
    'code_sha256',local_sha([mfilename('fullpath') '.m']));
write_checkpoint();
for j=1:numel(jobs)
    budget=jsondecode(fileread(fullfile(run_dir,'budget.json')));
    base=p4.elapsed_s;
    previous={'P6','P8','Z4','Z8','SP4','SP6','SP8'}; weights=[1.5 2 2 4 1 1.5 2];
    for k=1:numel(previous)
        file=fullfile(run_dir,previous{k},'checkpoint.json');
        if isfile(file), old=jsondecode(fileread(file)); base=max(base,old.elapsed_s/weights(k)); end
    end
    % Conservative remaining scope: all outstanding approved PE plus
    % reference/audit/packaging allowance. Never start if it cannot fit.
    weight_left=0;
    for k=1:numel(previous)
        if ~isfile(fullfile(run_dir,previous{k},'checkpoint.json')), weight_left=weight_left+weights(k); end
    end
    predicted=1.35*base*weight_left+1200;
    receipt.remaining_prediction_s=predicted; receipt.budget_used_s=budget.used_s;
    if budget.used_s+predicted>14400
        receipt.status='STOPPED_BUDGET_PREDICTION'; receipt.stopped_reason='Updated remaining scope cannot fit authorized cumulative budget.'; write_checkpoint(); return
    end
    receipt.next_job=jobs{j}; write_checkpoint();
    r=run_m1_logged_job_vertical(jobs{j},run_dir);
    receipt.completed_jobs{end+1}=jobs{j}; write_checkpoint();
    if ~strcmp(r.status,'COMPLETE') || r.pilot_stop
        receipt.status='STOPPED_JOB_GATE'; receipt.stopped_reason=[jobs{j} ': ' r.status]; write_checkpoint(); return
    end
end
receipt.status='COMPLETE'; receipt.next_job=''; write_checkpoint();
fprintf('M1_STAGE_RECEIPT %s\n',jsonencode(receipt));
    function write_checkpoint()
        fid=fopen(checkpoint,'w','n','UTF-8'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid);
    end
end
function value=local_sha(file)
md=java.security.MessageDigest.getInstance('SHA-256'); fid=fopen(file,'rb'); c=onCleanup(@()fclose(fid)); %#ok<NASGU>
while ~feof(fid), bytes=fread(fid,2^20,'*uint8'); if ~isempty(bytes), md.update(typecast(bytes,'int8')); end, end
value=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
