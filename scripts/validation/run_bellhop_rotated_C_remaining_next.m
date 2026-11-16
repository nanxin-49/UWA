function receipt=run_bellhop_rotated_C_remaining_next(run_dir)
%RESUME_BELLHOP_ROTATED_A_NEXT Dispatch at most one missing formal component.
s=load(fullfile(run_dir,'preflight.mat')); p=s.receipt; plan=s.plan; start=tic;
cpfile=fullfile(run_dir,'checkpoint.json');
if isfile(cpfile)
    receipt=jsondecode(fileread(cpfile)); assert(strcmp(receipt.status,'READY_NEXT_C_COMPONENT'),'Stop state forbids automatic retry.');
else
    receipt=struct('status','READY_NEXT_C_COMPONENT','result_dir',run_dir,'components',p.components,'completed_count',0, ...
        'new_native_count',0,'new_native_elapsed_s',0,'resume_elapsed_s',0,'prior_used_and_reserved_s',p.prior_accounted_s+p.preflight_elapsed_s, ...
        'plan_sha256',p.plan_sha256,'executable_sha256',p.executable_sha256,'scope','C only; no D/FFP or commits');
end
j=find(cellfun(@isempty,{receipt.components.path}),1); assert(~isempty(j),'No missing A component.');
try
    assert(strcmpi(black_file_hash(p.plan_file),p.plan_sha256) && strcmpi(black_file_hash(p.executable),p.executable_sha256));
    remain=3600-receipt.prior_used_and_reserved_s-receipt.resume_elapsed_s; assert(remain>60,'Shared time gate with60s audit reserve');
    [~,ram]=memory; receipt.latest_available_bytes=ram.PhysicalMemory.Available; assert(ram.PhysicalMemory.Available>=p.budget.minimum_available_bytes,'RAM launch gate');
    bytes=0; for k=1:numel(p.disk_roots), d=dir(fullfile(p.disk_roots{k},'**','*')); bytes=bytes+sum([d(~[d.isdir]).bytes]); end
    receipt.output_bytes=bytes; assert(bytes<p.budget.max_output_bytes,'Shared output-size gate');
    receipt.status='RUNNING_ONE_C_COMPONENT'; receipt.current_tag=plan.jobs(j).tag; write(cpfile,receipt);
    out=fullfile(run_dir,'C_Q625_h05_B401',plan.jobs(j).tag);
    job=plan.jobs(j); Bcp=jsondecode(fileread(p.B_components(j).checkpoint_path)); job.B_raw_root=Bcp.raw_root;
    if isfield(receipt,'input_verifier_recovery_tag') && strcmp(job.tag,receipt.input_verifier_recovery_tag)
        job.previous_unstarted_input_root=fullfile(out,'case'); out=[out '_input02'];
    end
    if isfield(receipt,'authorized_RAM_recovery_tag') && strcmp(job.tag,receipt.authorized_RAM_recovery_tag)
        assert(~isfield(job,'previous_unstarted_input_root')); job.previous_unstarted_input_root=fullfile(out,'case'); out=[out receipt.authorized_RAM_recovery_suffix];
    end
    result=run_bellhop_rotated_REAL8_C_component(p.executable,job,out,min(120,remain-60));
    receipt.components(j).path=fullfile(out,'component.mat'); receipt.components(j).checkpoint_path=fullfile(out,'checkpoint.json');
    receipt.components(j).component_sha256=black_file_hash(receipt.components(j).path); receipt.components(j).checkpoint_sha256=black_file_hash(receipt.components(j).checkpoint_path);
    receipt.new_native_count=receipt.new_native_count+1; receipt.new_native_elapsed_s=receipt.new_native_elapsed_s+result.receipt.process.elapsed_s;
    receipt.completed_count=receipt.completed_count+1; receipt.resume_elapsed_s=receipt.resume_elapsed_s+toc(start);
    receipt.cumulative_accounted_s=receipt.prior_used_and_reserved_s+receipt.resume_elapsed_s;
    receipt.status='READY_NEXT_C_COMPONENT'; if receipt.completed_count==39, receipt.status='ALL39_COMPONENTS_COMPLETE_ASSEMBLY_PENDING'; end
    write(cpfile,receipt); fprintf('C_PROGRESS %d/39 new%d native%.3fs cumulative%.3fs remaining%.3fs status%s\n',receipt.completed_count,receipt.new_native_count,receipt.new_native_elapsed_s,receipt.cumulative_accounted_s,3600-receipt.cumulative_accounted_s,receipt.status);
catch exception
    receipt.resume_elapsed_s=receipt.resume_elapsed_s+toc(start); receipt.cumulative_accounted_s=receipt.prior_used_and_reserved_s+receipt.resume_elapsed_s;
    receipt.status='STOPPED_NO_RETRY'; receipt.error=getReport(exception,'extended','hyperlinks','off'); write(cpfile,receipt); fprintf('%s\n',receipt.error);
end
end
function write(path,receipt)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); c=onCleanup(@()fclose(fid)); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true));
end
