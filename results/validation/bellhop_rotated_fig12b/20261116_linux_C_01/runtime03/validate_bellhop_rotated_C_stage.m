function receipt=validate_bellhop_rotated_C_stage(run_dir)
%VALIDATE_BELLHOP_ROTATED_B_STAGE Saved-field gates only; no C/D/FFP dispatch.
start=tic; data=load(fullfile(run_dir,'preflight.mat'),'receipt'); p=data.receipt;
cpfile=fullfile(run_dir,'checkpoint.json'); state=jsondecode(fileread(cpfile)); assert(strcmp(state.status,'ASSEMBLED_C_RAW_AUDIT_PENDING'));
receipt=struct('status','RUNNING_C_STAGE_AUDITS','run_dir',run_dir,'failed_checks',0,'completed_stage','none');
diary(fullfile(run_dir,'postprocess_diary.log')); cleanup=onCleanup(@()diary('off'));
try
    assert(strcmpi(black_file_hash(p.B_field_file),p.B_field_sha256),'B source changed');
    [~,ram]=bellhop_linux_memory(); assert(ram.PhysicalMemory.Available>=1.5*2^30,'C postprocess RAM gate');
    receipt.raw=audit_bellhop_rotated_C_raw(fullfile(run_dir,'assembly'),fullfile(run_dir,'raw_audit01')); receipt.completed_stage='raw'; write(fullfile(run_dir,'checkpoint_numeric.json'),receipt);
    assert(receipt.raw.failed_checks==0,'C raw reconstruction gate');
    receipt.tail=assess_bellhop_rotated_field_change(fullfile(run_dir,'assembly','sum_through_order12.mat'),fullfile(run_dir,'assembly','result.mat'),fullfile(run_dir,'N12_to_N16_audit01'),label='C_N12_to_N16',relative_tolerance=.001,envelope_tolerance=.001,phase_tolerance_rad=.01);
    receipt.completed_stage='tail'; write(fullfile(run_dir,'checkpoint_numeric.json'),receipt); assert(receipt.tail.failed_checks==0,'C order convergence gate');
    receipt.change=assess_bellhop_rotated_field_change(p.B_field_file,fullfile(run_dir,'assembly','result.mat'),fullfile(run_dir,'B_to_C_audit01'),label='B_to_C_step_only',relative_tolerance=.02,envelope_tolerance=.02,phase_tolerance_rad=.05);
    receipt.completed_stage='B_to_C'; write(fullfile(run_dir,'checkpoint_numeric.json'),receipt); assert(receipt.change.failed_checks==0,'B-toC original numerical gate');
    receipt.boundary=check_bellhop_rotated_A_boundaries(run_dir); receipt.completed_stage='boundary'; assert(receipt.boundary.failed_checks==0,'C top pressure release gate');
    receipt.status='C_STAGE_GATES_PASS_STOPPED_FOR_REVIEW';
catch exception
    receipt.status='STOPPED_C_STAGE_FAILED_NO_RETRY'; receipt.failed_checks=1; receipt.error=getReport(exception,'extended','hyperlinks','off'); fprintf('%s\n',receipt.error);
end
receipt.elapsed_s=toc(start); write(fullfile(run_dir,'checkpoint_numeric.json'),receipt);
state.resume_elapsed_s=state.resume_elapsed_s+receipt.elapsed_s; state.cumulative_accounted_s=state.prior_used_and_reserved_s+state.resume_elapsed_s; state.status=receipt.status;
state.remaining_budget_s=3600-state.cumulative_accounted_s; state.numeric_checkpoint=fullfile(run_dir,'checkpoint_numeric.json'); state.next_allowed_step='Saved-field figures/report and manual review only. No D/FFP/new solves or commit/push/merge.';
if isfield(state,'error'), state.historical_metadata_error=state.error; state=rmfield(state,'error'); end
write(cpfile,state); fprintf('C_NUMERIC %s cumulative%.3fs remaining%.3fs\n',receipt.status,state.cumulative_accounted_s,state.remaining_budget_s);
end
function write(path,value)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); c=onCleanup(@()fclose(fid)); fprintf(fid,'%s',jsonencode(value,PrettyPrint=true));
end