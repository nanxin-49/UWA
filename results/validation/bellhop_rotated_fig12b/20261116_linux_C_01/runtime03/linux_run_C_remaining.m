function receipt=linux_run_C_remaining(project_root,run_dir)
% Sole C9 dispatcher after qualification. Original artifacts remain immutable.
addpath(fullfile(project_root,'scripts','validation'),fullfile(project_root,'scripts','validation','support'),fullfile(project_root,'scripts','reporting')); addpath(fullfile(run_dir,'runtime03'),'-begin');
q=jsondecode(fileread(fullfile(run_dir,'checkpoint_qualification03.json'))); assert(strcmp(q.status,'QUALIFICATION_PASS_C9_ALLOWED')&&q.failed_checks==0); budget_initial=q.remaining_original_budget_s; epoch=str2double(getenv('BL_JOB_START_EPOCH')); reserve=60;
out=fullfile(run_dir,'C_resume'); assert(~isfolder(out)); mkdir(out); diary(fullfile(run_dir,'C_diary.log')); guard=onCleanup(@()diary('off'));
old_dir=fullfile(project_root,'results','validation','bellhop_rotated_fig12b','20261115_C_run_01'); old=jsondecode(fileread(fullfile(old_dir,'checkpoint.json'))); assert(old.completed_count==30&&old.new_native_count==30);
s=load(fullfile(old_dir,'preflight.mat'),'receipt','plan'); plan=s.plan; p=bellhop_linux_map(s.receipt); p.executable=fullfile(run_dir,'bin','bellhop_rotated_receiver64_2020'); p.executable_sha256=black_file_hash(p.executable); assert(strcmp(p.executable_sha256,q.executable_sha256)); receipt=p; save(fullfile(out,'preflight.mat'),'receipt','plan');
state=bellhop_linux_map(old); state.status='RUNNING_REMOTE_PREFLIGHT'; state.components=bellhop_linux_map(old.components); state.result_dir=out; state.reused_completed_C_count=30; state.new_native_at_recovery=30; state.new_native_elapsed_this_round_s=0; state.prior_used_and_reserved_s=3600-budget_initial; state.scope='C remaining9 only; no D/FFP/commit/push';
try
 assert(strcmpi(black_file_hash(p.plan_file),p.plan_sha256));
 for k=1:30, c=state.components(k); assert(strcmpi(black_file_hash(c.path),c.component_sha256)&&strcmpi(black_file_hash(c.checkpoint_path),c.checkpoint_sha256)); end
 for group={p.A_dir,p.B_dir}
  cp=jsondecode(fileread(fullfile(group{1},'checkpoint.json')));
  for k=1:39, c=cp.components(k); assert(strcmpi(black_file_hash(c.path),c.component_sha256)&&strcmpi(black_file_hash(c.checkpoint_path),c.checkpoint_sha256)); end
 end
 assert(strcmpi(black_file_hash(fullfile(p.A_dir,'checkpoint_stage.json')),p.A_stage_checkpoint_sha256)&&strcmpi(black_file_hash(fullfile(p.B_dir,'checkpoint_stage.json')),p.B_stage_checkpoint_sha256));
 state.status='READY_NEXT_C_COMPONENT'; sync();
 for j=31:39
  assert(remaining()>reserve,'Shared budget60s audit reserve'); resource_gate(); job=plan.jobs(j); bc=jsondecode(fileread(p.B_components(j).checkpoint_path)); job.B_raw_root=bellhop_linux_path(bc.raw_root);
  if j==31, assert(strcmp(job.tag,'down_order12')); job.previous_unstarted_input_root=fullfile(old_dir,'C_Q625_h05_B401','down_order12','case'); assert(~isfile([job.previous_unstarted_input_root '.shd'])); end
  dest=fullfile(out,'C_Q625_h05_B401',job.tag); state.current_tag=job.tag; state.status='RUNNING_ONE_C_COMPONENT'; sync();
  result=run_bellhop_rotated_REAL8_C_component(p.executable,job,dest,min(120,remaining()-reserve));
  state.components(j).path=fullfile(dest,'component.mat'); state.components(j).checkpoint_path=fullfile(dest,'checkpoint.json'); state.components(j).component_sha256=black_file_hash(state.components(j).path); state.components(j).checkpoint_sha256=black_file_hash(state.components(j).checkpoint_path);
  state.completed_count=j; state.new_native_count=j; state.new_native_elapsed_s=state.new_native_elapsed_s+result.receipt.process.elapsed_s; state.new_native_elapsed_this_round_s=state.new_native_elapsed_this_round_s+result.receipt.process.elapsed_s; state.status='READY_NEXT_C_COMPONENT'; sync(); fprintf('REMOTE_C_PROGRESS %d/39 native_new%.3fs remaining%.3fs\n',j,state.new_native_elapsed_this_round_s,remaining());
 end
 state.status='ALL39_COMPONENTS_COMPLETE_ASSEMBLY_PENDING'; sync(); resource_gate(); assert(remaining()>reserve); assembly=assemble_bellhop_rotated_C_run(out); %#ok<NASGU>
 state=jsondecode(fileread(fullfile(out,'checkpoint.json'))); sync(); number=validate_bellhop_rotated_C_stage(out); state=jsondecode(fileread(fullfile(out,'checkpoint.json'))); sync(); assert(strcmp(number.status,'C_STAGE_GATES_PASS_STOPPED_FOR_REVIEW')&&number.failed_checks==0,'C numerical gate');
 resource_gate(); assert(remaining()>30,'Shared budget before figures'); figs=fullfile(run_dir,'ABC_figures'); fig=render_bellhop_rotated_ABC_comparison(p.A_field_file,p.B_field_file,fullfile(out,'assembly','result.mat'),figs); %#ok<NASGU>
 fig_audit=audit_bellhop_rotated_ABC_figures(p.A_field_file,p.B_field_file,fullfile(out,'assembly','result.mat'),figs,fullfile(out,'figure_audit01')); assert(fig_audit.failed_checks==0); state.status='C_COMPLETE_ALL_STAGE_GATES_PASS_STOPPED_FOR_REVIEW'; state.figure_dir=figs; sync();
 receipt=struct('status',state.status,'C_completed',39,'C_reused',30,'new_formal_C',9,'new_formal_C_native_s',state.new_native_elapsed_this_round_s,'C_all_native_s',state.new_native_elapsed_s,'numeric',number,'figure_audit',fig_audit,'field_sha256',black_file_hash(fullfile(out,'assembly','result.mat')),'binary_sha256',p.executable_sha256,'original_plan_sha256',p.plan_sha256,'remaining_original_budget_s',remaining(),'accounted_cumulative_s',3600-remaining(),'scope',state.scope,'unrun',{{'D','new FFP','final independent paper qualification'}});
catch e
 state.status='STOPPED_REMOTE_C_NO_RETRY'; state.error=getReport(e,'extended','hyperlinks','off'); sync(); receipt=struct('status',state.status,'C_completed',state.completed_count,'new_formal_C',state.completed_count-30,'remaining_original_budget_s',remaining(),'accounted_cumulative_s',3600-remaining(),'error',state.error); fprintf('%s\n',state.error);
end
receipt.MATLAB_version=version; receipt.A_B_recomputed=0; receipt.C_completed_recomputed=0; receipt.windows_checkpoint_overwritten=false; write(fullfile(run_dir,'checkpoint_C_stage.json'),receipt); sync(); fprintf('REMOTE_C_FINAL %s count%d remaining%.3fs\n',receipt.status,state.completed_count,remaining());
 function r=remaining(), r=budget_initial-(posixtime(datetime('now','TimeZone','UTC'))-epoch); end
 function sync()
  state.resume_elapsed_s=budget_initial-remaining(); state.cumulative_accounted_s=3600-remaining(); state.remaining_budget_s=remaining(); state.next_allowed_step='Only currently authorized missing C tasks; after complete or failure manual review. No D/FFP or commits.'; write(fullfile(out,'checkpoint.json'),state);
 end
 function resource_gate()
  [~,ram]=bellhop_linux_memory(); assert(ram.PhysicalMemory.Available>=1.5*2^30,'Remote RAM launch gate'); drive=java.io.File(out); assert(drive.getUsableSpace()>=20*2^30,'Remote free disk gate'); d=dir(fullfile(run_dir,'**','*')); bytes=sum([d(~[d.isdir]).bytes]); assert(10.630438*2^30+bytes<15*2^30,'Logical original output budget; transferred replicas excluded from new output'); state.original_logical_output_GiB=10.630438+bytes/2^30;
 end
end
function write(path,value)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); guard=onCleanup(@()fclose(fid)); fprintf(fid,'%s',jsonencode(value,PrettyPrint=true));
end
