function receipt=seal_bellhop_rotated_C_stop(run_dir)
% Read-only evidence audit and stop receipt. Never dispatches Bellhop.
start=tic; state_path=fullfile(run_dir,'checkpoint.json'); state=jsondecode(fileread(state_path)); assert(strcmp(state.status,'STOPPED_NO_RETRY')&&state.completed_count==30);
s=load(fullfile(run_dir,'preflight.mat'),'receipt'); p=s.receipt; assert(contains(state.error,'Immediate native launch RAM gate'));
failed_dir=fullfile(run_dir,'C_Q625_h05_B401',state.current_tag); assert(~isfile(fullfile(failed_dir,'pid.txt'))&&~isfile(fullfile(failed_dir,'case.shd')));
assert(isempty(dir(fullfile(failed_dir,'*pid*'))));
manifest_counts=zeros(1,2); groups={p.A_dir,p.B_dir};
for g=1:2
 m=readtable(fullfile(groups{g},'sha256_manifest.csv'),'TextType','string'); manifest_counts(g)=height(m);
 for k=1:height(m), assert(strcmpi(black_file_hash(char(m.path(k))),m.sha256(k)),'A/B manifest drift'); end
 old=jsondecode(fileread(fullfile(groups{g},'checkpoint.json')));
 for k=1:39, c=old.components(k); assert(strcmpi(black_file_hash(c.path),c.component_sha256)); assert(strcmpi(black_file_hash(c.checkpoint_path),c.checkpoint_sha256)); end
end
original25=dir(fullfile(fileparts(run_dir),'20261114_full_01','A_Q125_h1_B401','*','component.mat')); assert(numel(original25)==25);
for k=1:25
 path=fullfile(original25(k).folder,original25(k).name); cp=jsondecode(fileread(fullfile(original25(k).folder,'checkpoint.json'))); ix=find(strcmp({cp.bindings.path},path)); assert(isscalar(ix)&&strcmpi(black_file_hash(path),cp.bindings(ix).sha256));
end
timing=table(strings(30,1),zeros(30,1),zeros(30,1),zeros(30,1),'VariableNames',{'tag','native_s','peak_MiB','runtime_min_available_GiB'}); checks_count=0; input_checks_count=0; binding_count=0;
for k=1:30
 c=state.components(k); cp=jsondecode(fileread(c.checkpoint_path)); assert(strcmp(cp.status,'COMPLETE')&&cp.failed_checks==0); assert(strcmpi(black_file_hash(c.path),c.component_sha256)&&strcmpi(black_file_hash(c.checkpoint_path),c.checkpoint_sha256));
 for j=1:numel(cp.bindings), assert(strcmpi(black_file_hash(cp.bindings(j).path),cp.bindings(j).sha256),'C binding drift'); binding_count=binding_count+1; end
 for j=1:numel(cp.code_versions)
  v=cp.code_versions(j); actual=v.path;
  if ~strcmpi(black_file_hash(actual),v.sha256), [~,n,e]=fileparts(actual); actual=fullfile(run_dir,'executed_sources_initial',[n e]); end
  assert(isfile(actual)&&strcmpi(black_file_hash(actual),v.sha256),'Executed C source unavailable');
 end
 d=fileparts(c.path); t=readtable(fullfile(d,'checks.csv')); assert(all(t.pass)); checks_count=checks_count+height(t);
 t=readtable(fullfile(d,'input_contract_checks.csv'),'Delimiter',',','ReadVariableNames',true,'VariableNamingRule','preserve'); assert(all(t.pass)); input_checks_count=input_checks_count+height(t);
 assert(strcmp(strtrim(fileread(fullfile(d,'state.txt'))),'TERMINAL'));
 timing(k,:)={string(c.tag),cp.process.elapsed_s,cp.process.peak_working_set_bytes/2^20,cp.process.minimum_observed_available_bytes/2^30};
end
writetable(timing,fullfile(run_dir,'completed_component_timing.csv'));
files=dir(fullfile(run_dir,'**','*')); files=files(~[files.isdir]); manifest=table(strings(0,1),strings(0,1),'VariableNames',{'path','sha256'});
for k=1:numel(files)
 path=fullfile(files(k).folder,files(k).name);
 if strcmp(files(k).folder,run_dir)&&ismember(files(k).name,{'checkpoint.json','checkpoint_stage.json','sha256_manifest.csv'}), continue; end
 manifest(end+1,:)={string(path),string(black_file_hash(path))};
end
writetable(manifest,fullfile(run_dir,'sha256_manifest.csv'));
bytes=0; for k=1:numel(p.disk_roots), d=dir(fullfile(p.disk_roots{k},'**','*')); bytes=bytes+sum([d(~[d.isdir]).bytes]); end
[~,ram]=memory; drive=java.io.File(run_dir); elapsed=toc(start); state.resume_elapsed_s=state.resume_elapsed_s+elapsed; state.cumulative_accounted_s=state.prior_used_and_reserved_s+state.resume_elapsed_s; state.remaining_budget_s=3600-state.cumulative_accounted_s; state.output_bytes=bytes;
state.next_allowed_step='Manual review only; no automatic retry. Reuse30 completed C components; resolve RAM and review resource-stopped down_order12 before any future authorization. No D/FFP/commit/push/merge.';
write(state_path,state);
receipt=struct('status','C_PARTIAL_30_OF39_RAM_STOPPED_FOR_REVIEW','completed_C',30,'remaining_C',9,'native_elapsed_s',sum(timing.native_s),'maximum_native_job_s',max(timing.native_s),'maximum_native_peak_MiB',max(timing.peak_MiB),'minimum_observed_runtime_available_GiB',min(timing.runtime_min_available_GiB), ...
 'component_checks',checks_count,'input_contract_checks',input_checks_count,'C_binding_SHA_checks',binding_count,'C_manifest_files',height(manifest),'A_B_manifest_counts',manifest_counts,'A_B_all78_SHA_unchanged',true,'original25_SHA_unchanged',true,'A_B_recomputed',0,'native_retries',0, ...
 'stop_tag',state.current_tag,'stop_before_native',true,'failed_launch_threshold_GiB',1.5,'failed_launch_actual_available_bytes','unknown: original assertion did not persist the immediate RAM value','controller_previous_available_GiB',state.latest_available_bytes/2^30, ...
 'seal_current_available_GiB',ram.PhysicalMemory.Available/2^30,'output_GiB',bytes/2^30,'free_disk_GiB',drive.getUsableSpace()/2^30,'elapsed_s',elapsed,'accounted_cumulative_s',state.cumulative_accounted_s,'remaining_original_budget_s',state.remaining_budget_s, ...
 'plan_sha256',p.plan_sha256,'binary_sha256',p.executable_sha256,'unrun',{{'down_order12','up/down order13-16','C complete assembly','C N12-to16','B-toC complete field','ABC complete field figures','D','new FFP'}},'execution_mode','Existing shared MATLAB MCP; explicit target project_path');
write(fullfile(run_dir,'checkpoint_stage.json'),receipt); disp(jsonencode(receipt));
end
function write(path,value)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); c=onCleanup(@()fclose(fid)); fprintf(fid,'%s',jsonencode(value,PrettyPrint=true));
end
