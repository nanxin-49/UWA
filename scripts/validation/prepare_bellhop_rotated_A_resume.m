function receipt=prepare_bellhop_rotated_A_resume(project_root,output_dir)
%PREPARE_BELLHOP_ROTATED_A_RESUME Verify immutable plan/cache; no native solve.
assert(~isfolder(output_dir)); mkdir(output_dir); start=tic;
family=fullfile(project_root,'results','validation','bellhop_rotated_fig12b');
plan_file=fullfile(family,'20261114_full_suite_plan_03','A_Q125_h1_B401.mat');
assert(strcmpi(black_file_hash(plan_file),'d3566fc6d4e2bc163826ea74eebe61d7a806a2ebe6598ed1db31b97c6742558e'));
x=load(plan_file,'plan'); plan=x.plan; assert(numel(plan.jobs)==39);
old_group=fullfile(family,'20261114_full_03','A_Q125_h1_B401');
executable=fullfile(family,'20261115_performance_review_01','optimized_plain_build01','bin','bellhop_rotated_receiver64_2020.exe');
assert(strcmpi(black_file_hash(executable),'4d3b99f4365f778eff40b52714794eeee627b1f5b54af51baeba8b0b811ecdf2'));
proof=jsondecode(fileread(fullfile(family,'20261115_performance_review_01','saved_audit01','checkpoint.json'))); assert(proof.failed_checks==0);
components=repmat(struct('tag','','path','','checkpoint_path','','checkpoint_sha256','','component_sha256','','order',0),1,39);
completed=[]; remaining=[]; bindings=struct('path',{},'sha256',{});
for j=1:numel(plan.jobs)
    job=plan.jobs(j); cpfile=fullfile(old_group,job.tag,'checkpoint.json');
    components(j).tag=job.tag; components(j).order=job.order;
    if isfile(cpfile)
        cp=jsondecode(fileread(cpfile));
        if strcmp(cp.status,'COMPLETE')
            assert(strcmp(jsonencode(cp.cfg),jsonencode(job.cfg)),'Frozen config mismatch.');
            assert(strcmpi(cp.executable_sha256,'93ab5754e45525ab0bff0e5490035c74ffedac83472e9c298b3690eb9928c107'));
            for k=1:numel(cp.bindings)
                assert(strcmpi(black_file_hash(cp.bindings(k).path),cp.bindings(k).sha256),'Cached output SHA mismatch.');
                bindings(end+1)=cp.bindings(k); %#ok<AGROW>
            end
            path=fullfile(fileparts(cpfile),'component.mat');
            components(j).path=path; components(j).checkpoint_path=cpfile;
            components(j).checkpoint_sha256=black_file_hash(cpfile); components(j).component_sha256=black_file_hash(path);
            completed(end+1)=j; %#ok<AGROW>
            continue
        end
        if isfile(fullfile(fileparts(cpfile),'pid.txt'))
            assert(strcmp(strtrim(fileread(fullfile(fileparts(cpfile),'state.txt'))),'TERMINAL'),'Prior process not terminal.');
        end
    end
    remaining(end+1)=j; %#ok<AGROW>
end
assert(numel(completed)==31 && numel(remaining)==8 && all([plan.jobs(remaining).order]>=13));
original=dir(fullfile(family,'20261114_full_01','A_Q125_h1_B401','*','component.mat')); assert(numel(original)==25);
for k=1:numel(original)
    path=fullfile(original(k).folder,original(k).name); cp=jsondecode(fileread(fullfile(original(k).folder,'checkpoint.json')));
    index=find(strcmp({cp.bindings.path},path)); assert(isscalar(index)); assert(strcmpi(black_file_hash(path),cp.bindings(index).sha256));
end
original_auth=fullfile(family,'20261114_full_authorization_01','authorization.json'); auth=jsondecode(fileread(original_auth));
assert(auth.max_wall_s==3600 && auth.max_job_s==120 && auth.minimum_available_bytes==1.5*2^30 && auth.minimum_runtime_available_bytes==256*2^20);
ledger=jsondecode(fileread(fullfile(family,'20261114_metadata_fix_01','checkpoint.json')));
used=ledger.cumulative_measured_before_resume_s; accounting=struct('path',{},'elapsed_s',{});
% Count old metadata qualification once; add subsequent measured diagnostic work.
extra={'20261114_receiver64_minimal_01/case01/checkpoint_preflight.json','20261114_receiver64_minimal_01/case02/checkpoint_preflight.json', ...
'20261114_receiver64_minimal_01/audit01/checkpoint.json','20261114_receiver64_minimal_01/saved_audit01/checkpoint.json', ...
'20261115_performance_review_01/baseline_profile_case02/checkpoint.json','20261115_performance_review_01/baseline_clock_case01/checkpoint.json', ...
'20261115_performance_review_01/optimized_clock_case01/checkpoint.json','20261115_performance_review_01/timing_pair1_baseline/checkpoint.json', ...
'20261115_performance_review_01/timing_pair1_optimized/checkpoint.json','20261115_performance_review_01/timing_pair2_optimized/checkpoint.json', ...
'20261115_performance_review_01/saved_audit01/checkpoint.json'};
for k=1:numel(extra)
    path=fullfile(family,extra{k}); cp=jsondecode(fileread(path)); used=used+cp.elapsed_s;
    accounting(end+1)=struct('path',path,'elapsed_s',cp.elapsed_s); %#ok<AGROW>
end
used=used+30+60; % retain old30s reserve and add60s for unmeasured preparation/checks
[~,ram]=memory; drive=java.io.File(output_dir); assert(drive.getUsableSpace()>=auth.minimum_free_disk_bytes,'Disk free-space gate.');
roots={'20261114_full_01','20261114_full_02','20261114_full_03','20261114_resume_anchor_01','20261114_metadata_fix_01', ...
'20261114_receiver64_minimal_01','20261115_performance_review_01'};
disk_roots=cellfun(@(name)fullfile(family,name),roots,'UniformOutput',false); disk_roots{end+1}=output_dir;
bytes=0; for k=1:numel(disk_roots), files=dir(fullfile(disk_roots{k},'**','*')); bytes=bytes+sum([files(~[files.isdir]).bytes]); end
assert(bytes<auth.max_output_bytes,'Original cumulative disk-output gate.');
receipt=struct('status','PREPARED_A_ONLY_NO_SOLVE','result_dir',output_dir,'plan_file',plan_file,'plan_sha256',black_file_hash(plan_file), ...
    'executable',executable,'executable_sha256',black_file_hash(executable),'completed_count',31,'remaining_count',8, ...
    'remaining_tags',{{plan.jobs(remaining).tag}},'original25_SHA_unchanged',true,'components',components, ...
    'authorization_basis',original_auth,'authorization_basis_sha256',black_file_hash(original_auth), ...
    'user_scope_20261115','Resume only original formal A missing components with validated REAL8/plane overlay; no B/C/D, commit/push', ...
    'budget',auth,'prior_accounting',accounting,'prior_measured_and_reserved_s',used,'preflight_elapsed_s',toc(start), ...
    'remaining_budget_before_preflight_s',3600-used,'disk_roots',{disk_roots},'output_bytes',bytes, ...
    'available_RAM_bytes',ram.PhysicalMemory.Available,'RAM_start_gate_met',ram.PhysicalMemory.Available>=auth.minimum_available_bytes, ...
    'estimated_new_native_s',[90 150],'estimated_validation_and_assembly_s',[60 180], ...
    'estimated_native_peak_MiB',[60 120],'estimated_MATLAB_extra_peak_MiB',[300 600],'estimated_extra_output_GiB',[1 2], ...
    'native_processes_started',0,'stop_rule','Any resource/physics/SHA failure stops dispatch; no retries or threshold changes; B/C/D require review');
save(fullfile(output_dir,'preflight.mat'),'receipt','plan','components','completed','remaining','bindings');
fid=fopen(fullfile(output_dir,'checkpoint_preflight.json'),'w','n','UTF-8'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid);
disp(jsonencode(rmfield(receipt,{'components','budget','prior_accounting','disk_roots'})));
end