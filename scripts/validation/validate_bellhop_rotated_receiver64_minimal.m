function receipt=validate_bellhop_rotated_receiver64_minimal(preflight_file,baseline_root,output_dir)
%VALIDATE_BELLHOP_ROTATED_RECEIVER64_MINIMAL Exactly one authorized native job.
arguments
    preflight_file (1,:) char
    baseline_root (1,:) char
    output_dir (1,:) char
end
assert(~isfolder(output_dir)); mkdir(output_dir); started=tic;
s=load(preflight_file); assert(s.receipt.failed_checks==0);
root=s.receipt.case_root; executable=s.receipt.executable;
assert(strcmpi(black_file_hash(executable),s.receipt.executable_sha256));
assert(strcmpi(black_file_hash(s.receipt.plan_file),s.receipt.plan_sha256));
assert(~isfile([root '.shd']) && ~isfile([root '.rr64input']),'No repeat of an existing job.');
for ext={'.ssp','.iwcov','.rotq','.r25d'}
    assert(strcmpi(black_file_hash([root ext{1}]),black_file_hash([baseline_root ext{1}])),'Baseline physical input differs: %s',ext{1});
end
newlines=splitlines(string(fileread([root '.env']))); oldlines=splitlines(string(fileread([baseline_root '.env'])));
nonempty=find(strlength(strtrim(newlines))>0); changes=find(newlines~=oldlines);
assert(isequal(changes,nonempty(end-4)),'Only baseline ENV metadata range row may differ.');
assert(isequal(newlines,splitlines(string(fileread([root '.env'])))));
[~,mem]=memory; available_before=mem.PhysicalMemory.Available;
assert(available_before>=1.5*2^30,'Original1.5GiB launch gate.');
family=fileparts(fileparts(fileparts(root)));
protected_roots={fullfile(family,'20261114_full_03'),fullfile(family,'20261114_resume_anchor_01')};
protected_paths={s.receipt.plan_file};
for j=1:numel(protected_roots)
    files=dir(fullfile(protected_roots{j},'**','component.mat'));
    for k=1:numel(files), protected_paths{end+1}=fullfile(files(k).folder,files(k).name); end %#ok<AGROW>
end
for ext={'.env','.r25d','.prt'}, protected_paths{end+1}=[s.receipt.original_failed_root ext{1}]; end
before_hash=cellfun(@black_file_hash,protected_paths,'UniformOutput',false);
save(fullfile(output_dir,'before.mat'),'protected_paths','before_hash','available_before');
checks=table(strings(0,1),false(0,1),zeros(0,1),zeros(0,1),'VariableNames',{'check','pass','value','limit'});
process=[]; field=[]; status='FAILED'; failure=''; native_count=0;
try
    native_count=1;
    process=run_bellhop_rotated_local_case(executable,root,20,allow_task_termination=true,minimum_runtime_available_bytes=256*2^20);
    add('native_exit_zero',abs(process.exit_code),0);
    add('stderr_empty',strlength(string(process.stderr)),0);
    fid=fopen([root '.rr64input'],'rb','ieee-le'); assert(fid>=0);
    header=fread(fid,4,'int32'); input64=fread(fid,10981,'double'); evaluated64=fread(fid,10981,'double');
    q64=fread(fid,3,'double'); metadata32=fread(fid,10981,'float32'); tail=fread(fid,inf,'uint8'); fclose(fid);
    add('runtime_count_and_REAL8_kinds',double(~isequal(header,[10981;3;8;8])),0);
    add('input_audit_complete',numel(tail),0);
    add('runtime_ENV64_bits_match_preflight',bitdiff(input64,s.env_m64),0);
    add('runtime_effective_ranges_bitwise_frozen_MAT',bitdiff(evaluated64,s.r),0);
    add('runtime_effective_columns_bitwise_frozen_MAT',bitdiff(q64,s.q),0);
    add('runtime_input64_strict',nnz(diff(input64)<=0),0);
    add('runtime_evaluation64_strict',nnz(diff(evaluated64)<=0),0);
    add('runtime_native_metadata_strict',nnz(diff(metadata32)<=0),0);
    add('runtime_native_metadata_matches_planned',max(abs(metadata32-double(s.new_metadata))),0);
    fid=fopen([root '.rr64eval'],'rb','ieee-le'); h=fread(fid,2,'int32'); bytes=fread(fid,inf,'*uint8'); fclose(fid);
    assert(mod(numel(bytes),24)==0); rec=reshape(bytes,24,[]);
    ir=double(typecast(reshape(rec(1:4,:),[],1),'int32')); iz=double(typecast(reshape(rec(5:8,:),[],1),'int32'));
    ar=typecast(reshape(rec(9:16,:),[],1),'double'); aq=typecast(reshape(rec(17:24,:),[],1),'double');
    assert(all(ir>=1 & ir<=10981 & iz>=1 & iz<=3)); multiplicity=accumarray((ir-1)*3+iz,1,[10981*3 1]);
    add('actual_evaluation_header',double(~isequal(h,[10981;3])),0);
    add('all32943_evaluation_pairs_once',nnz(multiplicity~=1),0);
    add('actual_evaluation_range_bits',bitdiff(ar,s.r(ir)),0);
    add('actual_evaluation_column_bits',bitdiff(aq,s.q(iz)),0);
    fid=fopen([root '.rr64scale'],'rb','ieee-le'); hs=fread(fid,1,'int32'); bytes=fread(fid,inf,'*uint8'); fclose(fid);
    assert(mod(numel(bytes),12)==0); rec=reshape(bytes,12,[]);
    sr=double(typecast(reshape(rec(1:4,:),[],1),'int32')); scale_r=typecast(reshape(rec(5:12,:),[],1),'double');
    add('scale_header',abs(hs-10981),0);
    add('all10981_scale_ranges_once',nnz(accumarray(sr,1,[10981 1])~=1),0);
    add('scale_uses_exact_double_bits',bitdiff(scale_r,s.r(sr)),0);
    mw=readmatrix([root '.mwdiag'],'FileType','text','NumHeaderLines',1);
    rd=readmatrix([root '.r25diag'],'FileType','text','NumHeaderLines',1);
    add('all101_rays_logged',abs(size(mw,1)-101)+abs(size(rd,1)-101),0);
    add('all13_walls_and_7top6bottom',nnz(mw(:,2)~=13 | mw(:,3)~=13 | mw(:,4)~=7 | mw(:,5)~=6),0);
    add('ray_capacity_original_margin',max(mw(:,6)),.95*600000);
    add('all_rays_reach_final_receiver',max(s.r)-min(mw(:,7)),0);
    add('ray_unit_tangent_original_gate',max(rd(:,8)),1e-8);
    field=read_bellhop_rotated_field(root); baseline=read_bellhop_rotated_field(baseline_root);
    add('output_exact_range_binding',bitdiff(field.receiver_range_m,s.r),0);
    add('output_exact_column_binding',bitdiff(field.receiver_depth_m,s.q),0);
    add('pressure_count',abs(numel(field.pressure)-32943),0);
    add('complex_pressure_all_finite',nnz(~isfinite(field.pressure)),0);
    add('complex_pressure_bitwise_baseline',double(~isequal(field.pressure,baseline.pressure)),0);
    add('complex_pressure_max_absolute_difference',max(abs(field.pressure-baseline.pressure),[],'all'),0);
    relative_l2=norm(field.pressure(:)-baseline.pressure(:))/norm(baseline.pressure(:));
    physical_x_m=-s.q; physical_z_m=s.cfg_minimal.receiver_depths_m(:);
    add('physical_depth_forward_mapping_bitwise',bitdiff(64950+physical_z_m,s.r),0);
    after_hash=cellfun(@black_file_hash,protected_paths,'UniformOutput',false);
    add('all_existing_components_and_failure_inputs_unchanged',nnz(~strcmpi(before_hash,after_hash)),0);
    status='COMPLETE';
catch exception
    failure=getReport(exception,'extended','hyperlinks','off'); relative_l2=NaN;
end
[~,mem]=memory;
receipt=struct('status',status,'failed_checks',nnz(~checks.pass)+~strcmp(status,'COMPLETE'), ...
    'checks_count',height(checks),'elapsed_s',toc(started),'bellhop_processes_launched',native_count, ...
    'case_root',root,'output_dir',output_dir,'preflight_file',preflight_file,'baseline_root',baseline_root, ...
    'executable_sha256',black_file_hash(executable),'complex_relative_l2',relative_l2, ...
    'available_before_bytes',available_before,'available_after_bytes',mem.PhysicalMemory.Available, ...
    'protected_files_count',numel(protected_paths),'process',process,'matlab_release',version, ...
    'scope','Only up_order13,101beams,3original columns,10981original ranges; full A/order14-16/B/C/D remain paused', ...
    'next_allowed_step','Manual review; no new native solve authorized');
writetable(checks,fullfile(output_dir,'checks.csv'));
save(fullfile(output_dir,'result.mat'),'receipt','checks','field','process','protected_paths','before_hash','s');
write(fullfile(output_dir,'checkpoint.json'),jsonencode(receipt,PrettyPrint=true)); write(fullfile(output_dir,'error.log'),failure);
disp(jsonencode(receipt)); disp(checks(~checks.pass,:)); if ~isempty(failure), fprintf('%s\n',failure); end
    function add(name,value,limit)
        checks(end+1,:)={string(name),isfinite(value)&&value<=limit,value,limit};
    end
end
function n=bitdiff(a,b)
assert(numel(a)==numel(b)); n=nnz(typecast(a(:),'uint64')~=typecast(b(:),'uint64'));
end
function write(path,content)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); cleanup=onCleanup(@()fclose(fid)); fprintf(fid,'%s',content);
end