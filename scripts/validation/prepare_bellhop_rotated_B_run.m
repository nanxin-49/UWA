function receipt=prepare_bellhop_rotated_B_run(project_root,output_dir)
%PREPARE_BELLHOP_ROTATED_B_RUN Only frozen B; preserve all accepted A evidence.
assert(~isfolder(output_dir)); mkdir(output_dir); start=tic;
f=fullfile(project_root,'results','validation','bellhop_rotated_fig12b'); A_dir=fullfile(f,'20261115_A_resume_01');
A=jsondecode(fileread(fullfile(A_dir,'checkpoint_stage.json'))); assert(strcmp(A.status,'A_COMPLETE_STAGE_GATES_PASS_STOPPED_FOR_REVIEW'));
assert(strcmpi(black_file_hash(A.field_file),A.field_sha256));
manifest=readtable(fullfile(A_dir,'sha256_manifest.csv'),'TextType','string');
for k=1:height(manifest), assert(strcmpi(black_file_hash(char(manifest.path(k))),manifest.sha256(k)),'A manifest drift'); end
Astate=jsondecode(fileread(fullfile(A_dir,'checkpoint.json')));
for k=1:numel(Astate.components)
    c=Astate.components(k); assert(strcmpi(black_file_hash(c.path),c.component_sha256)); assert(strcmpi(black_file_hash(c.checkpoint_path),c.checkpoint_sha256));
end
B_file=fullfile(f,'20261114_full_suite_plan_03','B_Q625_h1_B401.mat'); assert(strcmpi(black_file_hash(B_file),'899a9b805c95646f43264584de14512b4d2017cec51517578625f590720625e6'));
b=load(B_file,'plan'); plan=b.plan; a=load(fullfile(f,'20261114_full_suite_plan_03','A_Q125_h1_B401.mat'),'plan');
assert(numel(plan.jobs)==39 && isequal(plan.depths_m,a.plan.depths_m) && isequal(plan.x_m,a.plan.x_m));
for k=1:39
    ac=a.plan.jobs(k); bc=plan.jobs(k); assert(strcmp(ac.tag,bc.tag));
    assert(ac.cfg.ssp_spacing_m==1.25 && bc.cfg.ssp_spacing_m==.625);
    q_fields={'ssp_spacing_m'};
    if isfield(ac.cfg,'horizontal_ssp_spacing_m')
        assert(ac.cfg.horizontal_ssp_spacing_m==1.25 && bc.cfg.horizontal_ssp_spacing_m==.625);
        q_fields{end+1}='horizontal_ssp_spacing_m';
    end
    assert(isequaln(rmfield(ac.cfg,q_fields),rmfield(bc.cfg,q_fields)),'B changes something besides Q');
    assert(isequal(ac.z_indices,bc.z_indices) && isequal(ac.x_indices,bc.x_indices) && isequal(ac.coverage,bc.coverage));
end
old=load(fullfile(A_dir,'preflight.mat'),'receipt'); auth=old.receipt.budget;
executable=old.receipt.executable; assert(strcmpi(black_file_hash(executable),'4d3b99f4365f778eff40b52714794eeee627b1f5b54af51baeba8b0b811ecdf2'));
[~,ram]=memory; disk_roots=old.receipt.disk_roots; disk_roots{end+1}=A.figure_dir; disk_roots{end+1}=output_dir;
bytes=0; for k=1:numel(disk_roots), d=dir(fullfile(disk_roots{k},'**','*')); bytes=bytes+sum([d(~[d.isdir]).bytes]); end
assert(bytes<auth.max_output_bytes); drive=java.io.File(output_dir); assert(drive.getUsableSpace()>=auth.minimum_free_disk_bytes);
components=repmat(struct('tag','','path','','checkpoint_path','','checkpoint_sha256','','component_sha256','','order',0),1,39);
for k=1:39, components(k).tag=plan.jobs(k).tag; components(k).order=plan.jobs(k).order; end
receipt=struct('status','PREPARED_B_ONLY_NO_SOLVE','plan_file',B_file,'plan_sha256',black_file_hash(B_file), ...
    'A_dir',A_dir,'A_field_file',A.field_file,'A_field_sha256',A.field_sha256,'A_all39_and_manifest277_preserved',true, ...
    'A_stage_checkpoint_sha256',black_file_hash(fullfile(A_dir,'checkpoint_stage.json')),'executable',executable,'executable_sha256',black_file_hash(executable), ...
    'prior_accounted_s',A.accounted_elapsed_s+40,'new_preparation_reserve_s',40,'preflight_elapsed_s',toc(start),'budget',auth, ...
    'disk_roots',{disk_roots},'output_bytes',bytes,'available_RAM_bytes',ram.PhysicalMemory.Available,'components',components, ...
    'user_scope','Only B39 tasks. Q1.25 to.625m only. No C/D/FFP/new production physics/commit/push/merge.', ...
    'scope_equivalence','All39cfg differ only ssp_spacing_m; receiving grid/order/coverage exactly unchanged');
save(fullfile(output_dir,'preflight.mat'),'receipt','plan'); write(fullfile(output_dir,'checkpoint_preflight.json'),receipt);
receipt_view=rmfield(receipt,{'components','budget','disk_roots'}); disp(jsonencode(receipt_view));
end
function write(path,value)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); c=onCleanup(@()fclose(fid)); fprintf(fid,'%s',jsonencode(value,PrettyPrint=true));
end