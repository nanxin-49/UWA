function receipt=recover_bellhop_rotated_B_metadata_only(run_dir)
%RECOVER_BELLHOP_ROTATED_B_METADATA_ONLY No solve/retry; saved successful job.
start=tic; state_path=fullfile(run_dir,'checkpoint.json'); state=jsondecode(fileread(state_path));
assert(strcmp(state.status,'STOPPED_NO_RETRY') && state.completed_count==0 && strcmp(state.current_tag,'direct_up_part1'));
folder=fullfile(run_dir,'B_Q625_h1_B401','direct_up_part1'); oldcp=fullfile(folder,'checkpoint.json');
old=jsondecode(fileread(oldcp)); assert(contains(old.error,'case.rotq') && old.process.exit_code==0 && strcmp(strtrim(fileread(fullfile(folder,'state.txt'))),'TERMINAL'));
checks=readtable(fullfile(folder,'checks.csv')); assert(all(checks.pass));
d=load(fullfile(folder,'component.mat')); assert(d.receipt.failed_checks==0 && strcmp(d.receipt.status,'COMPLETE'));
s=load(fullfile(run_dir,'preflight.mat'),'plan'); job=s.plan.jobs(1); assert(strcmp(jsonencode(job.cfg),jsonencode(d.cfg)));
raw=read_bellhop_rotated_field(old.raw_root); P=complex(zeros(size(d.pressure_ffp))); P(d.cfg.receiver_depth_order,d.cfg.receiver_x_order)=(-conj(raw.pressure)).';
assert(isequaln(P,d.pressure_ffp),'Saved raw physical pressure mismatch');
assert(~isfile([old.raw_root '.rotq']) && job.order==0,'Optional wall file scope differs');
receipt=d.receipt; receipt.bindings=struct('path',{},'sha256',{});
for ext={'.env','.ssp','.iwcov','.r25d','.shd','.prt','.r25diag','.mwdiag','.covdiag','.rr64input','.rr64eval','.rr64scale'}
    path=[old.raw_root ext{1}]; assert(isfile(path)); receipt.bindings(end+1)=struct('path',path,'sha256',black_file_hash(path)); %#ok<AGROW>
end
path=fullfile(folder,'component.mat'); receipt.bindings(end+1)=struct('path',path,'sha256',black_file_hash(path));
receipt.metadata_recovery=struct('original_failed_checkpoint',oldcp,'original_failed_checkpoint_sha256',black_file_hash(oldcp), ...
    'reason','Optional rotq applies only reflection, not direct job; all saved native outputs and numerical checks pass', ...
    'new_native_processes',0,'raw_pressure_max_difference',0,'checks',height(checks));
out=fullfile(folder,'completion_audit01'); assert(~isfolder(out)); mkdir(out); newcp=fullfile(out,'checkpoint.json'); write(newcp,receipt);
copyfile(state_path,fullfile(run_dir,'checkpoint_metadata_stop01.json'));
state.components(1).path=path; state.components(1).component_sha256=black_file_hash(path);
state.components(1).checkpoint_path=newcp; state.components(1).checkpoint_sha256=black_file_hash(newcp);
state.completed_count=1; state.new_native_count=1; state.new_native_elapsed_s=d.receipt.process.elapsed_s;
state.resume_elapsed_s=state.resume_elapsed_s+toc(start); state.cumulative_accounted_s=state.prior_used_and_reserved_s+state.resume_elapsed_s;
state.status='READY_NEXT_B_COMPONENT'; state.metadata_recovery=receipt.metadata_recovery; write(state_path,state);
fprintf('B_METADATA_ONLY_RECOVERY native_retries0 completed1/39 native%.4fs pressure_diff0 checks%d\n',state.new_native_elapsed_s,height(checks));
end
function write(path,value)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); c=onCleanup(@()fclose(fid)); fprintf(fid,'%s',jsonencode(value,PrettyPrint=true));
end