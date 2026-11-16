function receipt=prepare_bellhop_rotated_C_run(project_root,output_dir)
%PREPARE_BELLHOP_ROTATED_C_RUN Only frozen C; preserve accepted A/B evidence.
assert(~isfolder(output_dir)); mkdir(output_dir); start=tic;
f=fullfile(project_root,'results','validation','bellhop_rotated_fig12b'); A_dir=fullfile(f,'20261115_A_resume_01'); B_dir=fullfile(f,'20261115_B_run_02');
A=jsondecode(fileread(fullfile(A_dir,'checkpoint_stage.json'))); B=jsondecode(fileread(fullfile(B_dir,'checkpoint_stage.json')));
assert(strcmp(B.status,'B_COMPLETE_ALL_STAGE_GATES_PASS_STOPPED_FOR_REVIEW'));
manifest_counts=zeros(1,2);
for group={A_dir,B_dir}
    manifest=readtable(fullfile(group{1},'sha256_manifest.csv'),'TextType','string');
    for k=1:height(manifest), assert(strcmpi(black_file_hash(char(manifest.path(k))),manifest.sha256(k)),'Prior A/B manifest drift'); end
    old=jsondecode(fileread(fullfile(group{1},'checkpoint.json'))); assert(old.completed_count==39);
    for k=1:39
        c=old.components(k); assert(strcmpi(black_file_hash(c.path),c.component_sha256)); assert(strcmpi(black_file_hash(c.checkpoint_path),c.checkpoint_sha256));
    end
    if strcmp(group{1},A_dir), manifest_counts(1)=height(manifest); else, manifest_counts(2)=height(manifest); B_components=old.components; end
end
assert(strcmpi(black_file_hash(A.field_file),A.field_sha256)&&strcmpi(black_file_hash(B.field_file),B.field_sha256));
C_file=fullfile(f,'20261114_full_suite_plan_03','C_Q625_h05_B401.mat'); assert(strcmpi(black_file_hash(C_file),'cce53c33f09104675eed7ef6131b6b9d6154259b25baaeffa83570d331947e97'));
c=load(C_file,'plan'); plan=c.plan; b=load(fullfile(f,'20261114_full_suite_plan_03','B_Q625_h1_B401.mat'),'plan');
assert(numel(plan.jobs)==39 && isequal(plan.depths_m,b.plan.depths_m) && isequal(plan.x_m,b.plan.x_m));
for k=1:39
    bc=b.plan.jobs(k); cc=plan.jobs(k); assert(strcmp(bc.tag,cc.tag));
    assert(bc.cfg.ssp_spacing_m==.625 && cc.cfg.ssp_spacing_m==.625 && cc.cfg.step_m==bc.cfg.step_m/2);
    assert(isequaln(rmfield(bc.cfg,'step_m'),rmfield(cc.cfg,'step_m')),'C changes something besides step');
    assert(isequal(bc.z_indices,cc.z_indices)&&isequal(bc.x_indices,cc.x_indices)&&isequal(bc.coverage,cc.coverage));
end
old=load(fullfile(B_dir,'preflight.mat'),'receipt'); auth=old.receipt.budget; executable=old.receipt.executable;
assert(strcmpi(black_file_hash(executable),'4d3b99f4365f778eff40b52714794eeee627b1f5b54af51baeba8b0b811ecdf2'));
[~,ram]=memory; disk_roots=old.receipt.disk_roots; disk_roots{end+1}=B.figure_dir; disk_roots{end+1}=fullfile(f,'20261115_B_run_01'); disk_roots{end+1}=output_dir;
bytes=0; for k=1:numel(disk_roots), d=dir(fullfile(disk_roots{k},'**','*')); bytes=bytes+sum([d(~[d.isdir]).bytes]); end
assert(bytes<auth.max_output_bytes); drive=java.io.File(output_dir); assert(drive.getUsableSpace()>=auth.minimum_free_disk_bytes);
components=repmat(struct('tag','','path','','checkpoint_path','','checkpoint_sha256','','component_sha256','','order',0),1,39);
for k=1:39, components(k).tag=plan.jobs(k).tag; components(k).order=plan.jobs(k).order; end
receipt=struct('status','PREPARED_C_ONLY_NO_SOLVE','plan_file',C_file,'plan_sha256',black_file_hash(C_file), ...
    'A_dir',A_dir,'A_field_file',A.field_file,'A_field_sha256',A.field_sha256,'B_dir',B_dir,'B_field_file',B.field_file,'B_field_sha256',B.field_sha256, ...
    'AB_all78_and_manifests_preserved',true,'manifest_counts',manifest_counts,'B_components',B_components, ...
    'A_stage_checkpoint_sha256',black_file_hash(fullfile(A_dir,'checkpoint_stage.json')),'B_stage_checkpoint_sha256',black_file_hash(fullfile(B_dir,'checkpoint_stage.json')), ...
    'executable',executable,'executable_sha256',black_file_hash(executable),'prior_accounted_s',B.accounted_cumulative_s+20, ...
    'new_preparation_reserve_s',20,'preflight_elapsed_s',toc(start),'budget',auth,'disk_roots',{disk_roots},'output_bytes',bytes, ...
    'available_RAM_bytes',ram.PhysicalMemory.Available,'RAM_start_gate_met',ram.PhysicalMemory.Available>=auth.minimum_available_bytes,'components',components, ...
    'estimated_native_s',[500 800],'estimated_postprocess_s',[180 300],'estimated_new_output_GiB',[4 4.5], ...
    'user_scope','C39 only; Q.625 unchanged and all original per-chart steps halved; no A/B recompute,D,FFP,production change,commit,push,merge');
save(fullfile(output_dir,'preflight.mat'),'receipt','plan'); write(fullfile(output_dir,'checkpoint_preflight.json'),receipt);
disp(jsonencode(rmfield(receipt,{'B_components','components','budget','disk_roots'})));
end
function write(path,value)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); c=onCleanup(@()fclose(fid)); fprintf(fid,'%s',jsonencode(value,PrettyPrint=true));
end