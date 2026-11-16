function receipt=audit_bellhop_rotated_linux_return(project_root,returned_root,output_dir)
%AUDIT_BELLHOP_ROTATED_LINUX_RETURN Transfer/immutable/component evidence only.
% No propagation, resampling, coordinate migration or large-field reconstruction.
assert(~isfolder(output_dir)); mkdir(output_dir); start=tic; checks=table(strings(0,1),false(0,1),'VariableNames',{'check','pass'});
manifest=splitlines(string(fileread(fullfile(returned_root,'return_SHA256SUMS')))); verified=0;
for k=1:numel(manifest)
 line=char(manifest(k)); if isempty(strtrim(line)), continue; end
 sha=line(1:64); relative=strtrim(line(67:end)); path=fullfile(returned_root,strrep(relative,'/',filesep)); assert(strcmpi(black_file_hash(path),sha),'Returned SHA mismatch: %s',relative); verified=verified+1;
end
add('all_returned_SHA',true); state=jsondecode(fileread(fullfile(returned_root,'C_resume','checkpoint.json'))); stage=jsondecode(fileread(fullfile(returned_root,'checkpoint_C_stage.json'))); add('all39_and_reuse30',state.completed_count==39&&state.reused_completed_C_count==30); add('remote_numerical_pass',strcmp(stage.status,'C_COMPLETE_ALL_STAGE_GATES_PASS_STOPPED_FOR_REVIEW'));
component_checks=0; input_checks=0; bindings=0; timing=table(strings(39,1),zeros(39,1),zeros(39,1),zeros(39,1),'VariableNames',{'tag','native_s','peak_MiB','min_available_GiB'});
for k=1:39
 c=state.components(k); cpfile=local_path(c.checkpoint_path); matfile=local_path(c.path); assert(strcmpi(black_file_hash(cpfile),c.checkpoint_sha256)&&strcmpi(black_file_hash(matfile),c.component_sha256)); cp=jsondecode(fileread(cpfile)); assert(strcmp(cp.status,'COMPLETE')&&cp.failed_checks==0);
 d=fileparts(matfile); t=readtable(fullfile(d,'checks.csv'),'Delimiter',',','VariableNamingRule','preserve'); assert(all(t.pass)); component_checks=component_checks+height(t); t=readtable(fullfile(d,'input_contract_checks.csv'),'Delimiter',',','VariableNamingRule','preserve'); assert(all(t.pass)); input_checks=input_checks+height(t);
 if k>=31
  for j=1:numel(cp.bindings), assert(strcmpi(black_file_hash(local_path(cp.bindings(j).path)),cp.bindings(j).sha256)); bindings=bindings+1; end
 end
 timing(k,:)={string(c.tag),cp.process.elapsed_s,cp.process.peak_working_set_bytes/2^20,cp.process.minimum_observed_available_bytes/2^30};
end
add('C39_SHA_and_all_case_checks',true);
family=fullfile(project_root,'results','validation','bellhop_rotated_fig12b');
for group={'20261115_A_resume_01','20261115_B_run_02'}
 old=jsondecode(fileread(fullfile(family,group{1},'checkpoint.json')));
 for k=1:39, c=old.components(k); assert(strcmpi(black_file_hash(c.path),c.component_sha256)&&strcmpi(black_file_hash(c.checkpoint_path),c.checkpoint_sha256)); end
end
add('A_B_all78_original_SHA_unchanged',true);
original=dir(fullfile(family,'20261114_full_01','A_Q125_h1_B401','*','component.mat')); assert(numel(original)==25);
for k=1:25
 path=fullfile(original(k).folder,original(k).name); cp=jsondecode(fileread(fullfile(original(k).folder,'checkpoint.json'))); ix=find(strcmp({cp.bindings.path},path)); assert(isscalar(ix)&&strcmpi(black_file_hash(path),cp.bindings(ix).sha256));
end
add('original25_SHA_unchanged',true); add('Linux_qualification18_checks_pass',jsondecode(fileread(fullfile(returned_root,'checkpoint_qualification03.json'))).failed_checks==0);
field=fullfile(returned_root,'C_resume','assembly','result.mat'); add('final_complex_field_SHA',strcmpi(black_file_hash(field),stage.field_sha256));
receipt=struct('status','LINUX_C_RETURNED_EVIDENCE_AUDIT','failed_checks',nnz(~checks.pass),'checks_count',height(checks),'returned_SHA_files',verified,'component_checks',component_checks,'input_contract_checks',input_checks,'new_component_binding_SHA',bindings,'A_B_all78_SHA_unchanged',true,'original25_SHA_unchanged',true,'completed_C_recomputed',0,'new_formal_C',9,'elapsed_s',toc(start),'remote_checkpoint',stage,'field_sha256',black_file_hash(field),'new_physics_solves',0);
writetable(checks,fullfile(output_dir,'checks.csv')); writetable(timing,fullfile(output_dir,'component_timing.csv')); fid=fopen(fullfile(output_dir,'checkpoint.json'),'w','n','UTF-8'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid); fprintf('LINUX_RETURN SHA%d C39 AB78 original25 PASS case_checks%d input%d elapsed%.3fs\n',verified,component_checks,input_checks,receipt.elapsed_s);
 function add(name,pass), checks(end+1,:)={string(name),logical(pass)}; end
 function out=local_path(in)
  remote='/home/misc/cwx/bellhop_rotated_C_transfer_20261116_01/Explain'; run=[remote '/results/validation/bellhop_rotated_fig12b/20261116_linux_C_01'];
  if startsWith(in,run), out=[returned_root strrep(in(numel(run)+1:end),'/',filesep)]; elseif startsWith(in,remote), out=[project_root strrep(in(numel(remote)+1:end),'/',filesep)]; else, out=in; end
 end
end
