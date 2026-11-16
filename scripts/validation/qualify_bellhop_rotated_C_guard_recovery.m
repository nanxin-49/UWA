function receipt=qualify_bellhop_rotated_C_guard_recovery(run_dir)
%QUALIFY_BELLHOP_ROTATED_C_GUARD_RECOVERY No native process/recalculation.
start=tic; statefile=fullfile(run_dir,'checkpoint.json'); state=jsondecode(fileread(statefile));
assert(strcmp(state.status,'STOPPED_NO_RETRY') && state.completed_count==2 && state.new_native_count==2 && strcmp(state.current_tag,'direct_up_part3'));
pre=load(fullfile(run_dir,'preflight.mat'),'plan','receipt'); j=3; job=pre.plan.jobs(j); Bcp=jsondecode(fileread(pre.receipt.B_components(j).checkpoint_path));
folder=fullfile(run_dir,'C_Q625_h05_B401',job.tag); assert(~isfile(fullfile(folder,'pid.txt')) && ~isfile(fullfile(folder,'case.shd')),'Blocked input was already executed');
[checks,proof]=check_bellhop_rotated_C_input_contract(fullfile(folder,'case'),Bcp.raw_root,job.cfg); assert(all(checks.pass));
for k=1:2
    c=state.components(k); cp=jsondecode(fileread(c.checkpoint_path)); assert(strcmp(cp.status,'COMPLETE')&&cp.failed_checks==0);
    assert(strcmpi(black_file_hash(c.path),c.component_sha256)); assert(strcmpi(black_file_hash(c.checkpoint_path),c.checkpoint_sha256));
end
out=fullfile(run_dir,'input_guard_qualification01'); assert(~isfolder(out)); mkdir(out); writetable(checks,fullfile(out,'checks.csv'));
receipt=struct('status','PRE_SOLVER_GUARD_CHECKER_QUALIFIED_NO_SOLVE','checks_count',height(checks),'failed_checks',0,'elapsed_s',toc(start), ...
    'input_contract',proof,'native_processes_launched',0,'completed_C_recomputed',0,'blocked_input_has_no_PID',true, ...
    'scope','Only derived SSP guard follows unchanged original writer formula beyond unchanged Beam.Box; no physical nodes/c/receivers changed');
write(fullfile(out,'checkpoint.json'),receipt); copyfile(statefile,fullfile(run_dir,'checkpoint_input_guard_stop01.json'));
state.status='READY_NEXT_C_COMPONENT'; state.input_verifier_recovery_tag=job.tag; state.input_verifier_qualification=fullfile(out,'checkpoint.json');
state.resume_elapsed_s=state.resume_elapsed_s+receipt.elapsed_s; state.cumulative_accounted_s=state.prior_used_and_reserved_s+state.resume_elapsed_s;
state.historical_input_verifier_error=state.error; state=rmfield(state,'error'); write(statefile,state); disp(jsonencode(receipt));
end
function write(path,value)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); c=onCleanup(@()fclose(fid)); fprintf(fid,'%s',jsonencode(value,PrettyPrint=true));
end