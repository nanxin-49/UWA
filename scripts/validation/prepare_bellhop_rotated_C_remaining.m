function receipt=prepare_bellhop_rotated_C_remaining(project_root,run_dir,resume_dir)
% User-authorized RAM-stop recovery; no native process in this function.
assert(~isfolder(resume_dir)); mkdir(resume_dir); start=tic; diary(fullfile(resume_dir,'preflight_diary.log')); guard=onCleanup(@()diary('off'));
state_file=fullfile(run_dir,'checkpoint.json'); state=jsondecode(fileread(state_file)); assert(strcmp(state.status,'STOPPED_NO_RETRY')&&state.completed_count==30&&state.new_native_count==30);
assert(abs(state.cumulative_accounted_s-2580.8344009)<1e-6); assert(contains(state.error,'Immediate native launch RAM gate')); s=load(fullfile(run_dir,'preflight.mat'),'receipt','plan'); p=s.receipt;
[mu,ms]=memory; receipt=struct('status','PREFLIGHT_RUNNING_NO_SOLVE','matlab_version',version,'matlab_pid',feature('getpid'),'matlab_memory_bytes',mu.MemUsedMATLAB,'available_bytes',ms.PhysicalMemory.Available,'total_RAM_bytes',ms.PhysicalMemory.Total,'temporary_variables_released',{{}},'memory_diagnosis','Existing shared session has no large task-owned base arrays. Application/native/cache and other applications dominate; no unowned variables, figures or processes were cleared.','previous_accounted_s',state.cumulative_accounted_s,'previous_remaining_s',3600-state.cumulative_accounted_s);
assert(strcmpi(black_file_hash(p.plan_file),p.plan_sha256)&&strcmpi(black_file_hash(p.executable),p.executable_sha256));
manifest_paths={fullfile(p.A_dir,'sha256_manifest.csv'),fullfile(p.B_dir,'sha256_manifest.csv'),fullfile(run_dir,'sha256_manifest.csv')}; counts=zeros(1,3);
for g=1:3
 m=readtable(manifest_paths{g},'TextType','string','Delimiter',','); counts(g)=height(m);
 for k=1:height(m), assert(strcmpi(black_file_hash(char(m.path(k))),m.sha256(k)),'Pre-resume manifest drift'); end
end
for group={p.A_dir,p.B_dir,run_dir}
 old=jsondecode(fileread(fullfile(group{1},'checkpoint.json')));
 for k=1:old.completed_count
  c=old.components(k); assert(strcmpi(black_file_hash(c.path),c.component_sha256)&&strcmpi(black_file_hash(c.checkpoint_path),c.checkpoint_sha256));
 end
end
assert(strcmpi(black_file_hash(fullfile(p.A_dir,'checkpoint_stage.json')),p.A_stage_checkpoint_sha256)&&strcmpi(black_file_hash(fullfile(p.B_dir,'checkpoint_stage.json')),p.B_stage_checkpoint_sha256));
failed=fullfile(run_dir,'C_Q625_h05_B401','down_order12'); assert(isfolder(failed)&&isempty(dir(fullfile(failed,'*pid*')))&&~isfile(fullfile(failed,'case.shd'))); failedcp=jsondecode(fileread(fullfile(failed,'checkpoint.json'))); assert(contains(failedcp.error,'Immediate native launch RAM gate')&&~isfield(failedcp,'process'));
bytes=0; for k=1:numel(p.disk_roots), d=dir(fullfile(p.disk_roots{k},'**','*')); bytes=bytes+sum([d(~[d.isdir]).bytes]); end
[~,ms]=memory; disk=java.io.File(run_dir); receipt.current_available_bytes=ms.PhysicalMemory.Available; receipt.free_disk_bytes=disk.getUsableSpace(); receipt.output_bytes=bytes; receipt.manifest_counts=counts; receipt.completed_C_SHA_revalidated=30; receipt.A_B_checkpoint_components_revalidated=78;
assert(ms.PhysicalMemory.Available>=p.budget.minimum_available_bytes,'Resume RAM launch gate'); assert(bytes<p.budget.max_output_bytes&&disk.getUsableSpace()>=p.budget.minimum_free_disk_bytes,'Resume disk gate');
mkdir(fullfile(resume_dir,'history')); for name={'checkpoint.json','checkpoint_stage.json','checkpoint_delivery.json'}, path=fullfile(run_dir,name{1}); copyfile(path,fullfile(resume_dir,'history',name{1})); end
receipt.status='READY_AUTHORIZED_C_REMAINING9'; receipt.elapsed_s=toc(start); receipt.diagnostic_and_preparation_reserve_s=20; receipt.failed_checks=0; receipt.preflight_only_no_native=true; receipt.resume_dir=resume_dir;
state.resume_elapsed_s=state.resume_elapsed_s+receipt.elapsed_s+20; state.cumulative_accounted_s=state.prior_used_and_reserved_s+state.resume_elapsed_s; state.remaining_budget_s=3600-state.cumulative_accounted_s; assert(state.remaining_budget_s>500,'Insufficient shared budget for9 jobs and audits');
state.historical_RAM_error=state.error; state=rmfield(state,'error'); state.authorized_RAM_recovery_tag='down_order12'; state.authorized_RAM_recovery_suffix='_resume_20261116_01'; state.recovery_evidence_dir=resume_dir; state.reused_completed_C_count=30; state.new_native_at_recovery=30; state.status='READY_NEXT_C_COMPONENT'; state.next_allowed_step='Only missing C9 serial under original resource gates. Stop on any new failure; no D/FFP or commits.'; write(state_file,state); receipt.remaining_original_budget_s=state.remaining_budget_s; write(fullfile(resume_dir,'checkpoint_preflight.json'),receipt); fprintf('C_RECOVERY_PREFLIGHT SHA manifests[%s] C30 AB78 PASS elapsed%.3fs RAM%.3fGiB remaining%.3fs\n',num2str(counts),receipt.elapsed_s,ms.PhysicalMemory.Available/2^30,state.remaining_budget_s);
end
function write(path,value)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); guard=onCleanup(@()fclose(fid)); fprintf(fid,'%s',jsonencode(value,PrettyPrint=true));
end
