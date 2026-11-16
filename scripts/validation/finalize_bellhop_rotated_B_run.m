function receipt=finalize_bellhop_rotated_B_run(run_dir,figure_dir)
%FINALIZE_BELLHOP_ROTATED_B_RUN Stop after B; preserve A and original budget.
start=tic; s=load(fullfile(run_dir,'preflight.mat'),'receipt','plan'); p=s.receipt; state_path=fullfile(run_dir,'checkpoint.json'); state=jsondecode(fileread(state_path));
assert(strcmp(state.status,'B_STAGE_GATES_PASS_STOPPED_FOR_REVIEW') && state.completed_count==39 && state.new_native_count==39);
number=jsondecode(fileread(fullfile(run_dir,'checkpoint_numeric.json'))); assert(number.failed_checks==0);
figaudit=jsondecode(fileread(fullfile(run_dir,'figure_audit01','checkpoint.json'))); assert(figaudit.failed_checks==0);
fig=jsondecode(fileread(fullfile(figure_dir,'checkpoint.json'))); assert(fig.failed_checks==0);
Amanifest=readtable(fullfile(p.A_dir,'sha256_manifest.csv'),'TextType','string');
for k=1:height(Amanifest), assert(strcmpi(black_file_hash(char(Amanifest.path(k))),Amanifest.sha256(k)),'A artifact changed'); end
Astate=jsondecode(fileread(fullfile(p.A_dir,'checkpoint.json')));
for k=1:39
    c=Astate.components(k); assert(strcmpi(black_file_hash(c.path),c.component_sha256)); assert(strcmpi(black_file_hash(c.checkpoint_path),c.checkpoint_sha256));
end
assert(strcmpi(black_file_hash(p.A_field_file),p.A_field_sha256)); assert(strcmpi(black_file_hash(fullfile(p.A_dir,'checkpoint_stage.json')),p.A_stage_checkpoint_sha256));
original25=dir(fullfile(fileparts(run_dir),'20261114_full_01','A_Q125_h1_B401','*','component.mat')); assert(numel(original25)==25);
for k=1:25
    path=fullfile(original25(k).folder,original25(k).name); cp=jsondecode(fileread(fullfile(original25(k).folder,'checkpoint.json'))); ix=find(strcmp({cp.bindings.path},path)); assert(isscalar(ix)); assert(strcmpi(black_file_hash(path),cp.bindings(ix).sha256));
end
timing=table(strings(39,1),zeros(39,1),zeros(39,1),zeros(39,1),zeros(39,1),'VariableNames',{'tag','order','native_s','native_peak_MiB','runtime_min_available_GiB'}); component_checks=0;
for k=1:39
    c=state.components(k); cp=jsondecode(fileread(c.checkpoint_path)); assert(strcmp(cp.status,'COMPLETE')&&cp.failed_checks==0);
    assert(strcmpi(black_file_hash(c.path),c.component_sha256)); assert(strcmpi(black_file_hash(c.checkpoint_path),c.checkpoint_sha256));
    assert(strcmp(strtrim(fileread(fullfile(fileparts(c.path),'state.txt'))),'TERMINAL'));
    checks=readtable(fullfile(fileparts(c.path),'checks.csv')); assert(all(checks.pass)); component_checks=component_checks+height(checks);
    timing(k,:)={string(c.tag),s.plan.jobs(k).order,cp.process.elapsed_s,cp.process.peak_working_set_bytes/2^20,cp.process.minimum_observed_available_bytes/2^30};
end
writetable(timing,fullfile(run_dir,'formal_component_timing.csv'));
state.resume_elapsed_s=state.resume_elapsed_s+fig.elapsed_s+figaudit.elapsed_s+toc(start); state.cumulative_accounted_s=state.prior_used_and_reserved_s+state.resume_elapsed_s;
assert(state.cumulative_accounted_s<=3600); state.remaining_budget_s=3600-state.cumulative_accounted_s;
roots=p.disk_roots; roots{end+1}=figure_dir; roots{end+1}=fullfile(fileparts(run_dir),'20261115_B_run_01'); bytes=0;
for k=1:numel(roots), d=dir(fullfile(roots{k},'**','*')); bytes=bytes+sum([d(~[d.isdir]).bytes]); end
assert(bytes<=15*2^30); state.output_bytes=bytes; state.status='B_COMPLETE_ALL_STAGE_GATES_PASS_STOPPED_FOR_REVIEW'; state.figure_dir=figure_dir;
state.next_allowed_step='Manual review only. No C,D,FFP solve,commit,push or merge.'; write(state_path,state);
receipt=struct('status',state.status,'formal_B_components',39,'A_components_recomputed',0,'original25_SHA_unchanged',true,'A_all39_and_manifest277_unchanged',true, ...
    'B_native_elapsed_s',sum(timing.native_s),'maximum_native_job_s',max(timing.native_s),'maximum_sampled_native_working_set_MiB',max(timing.native_peak_MiB), ...
    'minimum_observed_runtime_available_GiB',min(timing.runtime_min_available_GiB),'component_checks',component_checks, ...
    'raw_reconstruction_checks',number.raw.checks_count,'B_tail_failed_regions',number.tail.failed_checks,'A_to_B_failed_regions',number.change.failed_checks, ...
    'figure_checks',figaudit.checks_count,'figure_failed',figaudit.failed_checks,'top_pressure_over_envelope',number.boundary.top_pressure_release_max_over_envelope, ...
    'accounted_cumulative_s',state.cumulative_accounted_s,'remaining_original_budget_s',state.remaining_budget_s, ...
    'output_GiB',bytes/2^30,'remaining_original_output_GiB',15-bytes/2^30,'figure_dir',figure_dir,'field_file',fullfile(run_dir,'assembly','result.mat'), ...
    'field_sha256',black_file_hash(fullfile(run_dir,'assembly','result.mat')),'plan_sha256',p.plan_sha256,'binary_sha256',p.executable_sha256, ...
    'metadata_recovery','First direct job: optional rotq hash finalization recovered from saved pressure;0 native retry; old failure retained', ...
    'plot_path_correction','9 task-owned artifacts moved from default compare workspace to target worktree; SHA identical; original receipt and relocation log retained', ...
    'unrun',{{'C','D','FFP new calculation','full independentFFP comparison','final paper reproduction qualification'}}, ...
    'scope','B Q=.625m only, frozen step/beams/physics/native receiving grid/REAL8; serial, original gates retained');
write(fullfile(run_dir,'checkpoint_stage.json'),receipt); disp(jsonencode(receipt));
end
function write(path,value)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); c=onCleanup(@()fclose(fid)); fprintf(fid,'%s',jsonencode(value,PrettyPrint=true));
end