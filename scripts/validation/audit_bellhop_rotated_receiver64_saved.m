function receipt=audit_bellhop_rotated_receiver64_saved(run_dir,output_dir)
%AUDIT_BELLHOP_ROTATED_RECEIVER64_SAVED Independent saved-byte audit, no solver.
arguments
    run_dir (1,:) char
    output_dir (1,:) char
end
assert(~isfolder(output_dir)); mkdir(output_dir); t=tic;
a=load(fullfile(run_dir,'audit01','result.mat')); assert(a.receipt.failed_checks==0);
root=a.receipt.case_root; p=load(a.s.receipt.plan_file,'plan');
j=find(strcmp({p.plan.jobs.tag},'up_order13')); cfg=p.plan.jobs(j).cfg;
frozen_r=cfg.receiver_range_m(:); q=a.s.q; n=10981; nz=3;
checks=table(strings(0,1),false(0,1),zeros(0,1),zeros(0,1),'VariableNames',{'check','pass','value','limit'});
fid=fopen([root '.rr64input'],'rb','ieee-le'); header=fread(fid,4,'int32'); env64=fread(fid,n,'double'); eval64=fread(fid,n,'double'); q64=fread(fid,nz,'double'); fclose(fid);
add('saved_REAL8_header',double(~isequal(header,[n;nz;8;8])),0);
add('saved_effective_range_bits_vs_frozen_plan',bits(eval64,frozen_r),0);
add('saved_effective_column_bits',bits(q64,q),0);
fid=fopen([root '.rr64eval'],'rb','ieee-le'); h=fread(fid,2,'int32');
fseek(fid,8,'bof'); ir=fread(fid,inf,'int32',20); fseek(fid,12,'bof'); iz=fread(fid,inf,'int32',20);
fseek(fid,16,'bof'); er=fread(fid,inf,'double',16); fseek(fid,24,'bof'); eq=fread(fid,inf,'double',16); fclose(fid);
add('saved_evaluation_header',double(~isequal(h,[n;nz])),0);
add('saved_evaluation_coverage',nnz(accumarray((ir-1)*nz+iz,1,[n*nz 1])~=1),0);
add('saved_evaluation_range_bit_mismatches',bits(er,frozen_r(ir)),0);
add('saved_evaluation_transverse_bit_mismatches',bits(eq,q(iz)),0);
fid=fopen([root '.rr64scale'],'rb','ieee-le'); hs=fread(fid,1,'int32');
fseek(fid,4,'bof'); si=fread(fid,inf,'int32',8); fseek(fid,8,'bof'); se=fread(fid,inf,'double',4); fclose(fid);
add('saved_scale_header',abs(hs-n),0);
add('saved_scale_coverage',nnz(accumarray(si,1,[n 1])~=1),0);
add('saved_scale_bit_mismatches',bits(se,frozen_r(si)),0);
% Independently decode fixed Bellhop2020 SHD records; no heuristic/field helper.
[native_pressure,metadata_range,metadata_q]=raw_shd([root '.shd'],n,nz);
[old_pressure,~,~]=raw_shd([a.receipt.baseline_root '.shd'],n,nz);
add('raw_SHD_complex32_bitwise_vs_existing_baseline',double(~isequal(typecast(native_pressure(:),'uint32'),typecast(old_pressure(:),'uint32'))),0);
pressure_complex=complex(double(native_pressure(1:2:end,:)),double(native_pressure(2:2:end,:))).';
add('raw_SHD_vs_saved_MAT_pressure',max(abs(pressure_complex-a.field.pressure),[],'all'),0);
add('all_complex_pressure_finite',nnz(~isfinite(pressure_complex)),0);
add('SHD_metadata_strict_order',nnz(diff(metadata_range)<=0),0);
add('unchanged_two_ULP_range_metadata_gate',max(abs(metadata_range-frozen_r)./double(eps(single(metadata_range)))),2);
add('SHD_transverse_binding',max(abs(metadata_q-q)),0);
% Restore only the frozen index map; no sort, interpolation or r-offset inverse.
z_native=cfg.nominal_receiver_depths_m(cfg.receiver_depth_order).'; z_native=z_native(:);
x_native=-q; original_x_indices=zeros(1,nz);
for k=1:nz, hit=find(cfg.nominal_receiver_x_m==x_native(k)); assert(isscalar(hit)); original_x_indices(k)=hit; end
add('physical_z_forward_transform_exact_bits',bits(64950+z_native,frozen_r),0);
add('physical_x_original_index_mapping',max(abs(cfg.nominal_receiver_x_m(original_x_indices).'+q)),0);
physical_binding=struct('native_r_m',frozen_r,'native_q_m',q,'physical_z_m',z_native,'physical_x_m',x_native, ...
    'frozen_depth_indices',cfg.receiver_depth_order,'frozen_x_indices',original_x_indices, ...
    'native_pressure',pressure_complex,'mapping','p(native_q,native_r) retains original frozen indices; no interpolation');
family=fileparts(run_dir);
% Check original25 components against their original checkpoint hashes.
original_group=fullfile(family,'20261114_full_01','A_Q125_h1_B401'); files=dir(fullfile(original_group,'*','component.mat'));
add('original25_component_count',abs(numel(files)-25),0); historical_mismatches=0; component_files=cell(numel(files),1); component_hashes=component_files;
for k=1:numel(files)
    path=fullfile(files(k).folder,files(k).name); cp=jsondecode(fileread(fullfile(files(k).folder,'checkpoint.json')));
    match=find(strcmp({cp.bindings.path},path)); assert(isscalar(match)); actual=black_file_hash(path);
    historical_mismatches=historical_mismatches+~strcmpi(actual,cp.bindings(match).sha256); component_files{k}=path; component_hashes{k}=actual;
end
add('original25_checkpoint_SHA256_mismatches',historical_mismatches,0);
add('failed_original_ENV_SHA256_unchanged',double(~strcmpi(black_file_hash([a.s.receipt.original_failed_root '.env']),'e4746cb9b98c2d5b72df02fb1faac70bf0fa4a8d9408dad49e0f954d23d12527')),0);
add('failed_original_r25d_SHA256_unchanged',double(~strcmpi(black_file_hash([a.s.receipt.original_failed_root '.r25d']),'1845078c5423b885cc0f3f19fe642bda86deac459050ba71fcb04685851bc6f0')),0);
oldprt=fileread([a.s.receipt.original_failed_root '.prt']);
add('old_fatal_input_record_retained',double(~contains(oldprt,'Receiver ranges are not monotonically increasing')),0);
% Bind the actual compiled copies to retained reviewed overlay sources.
source_files=dir(fullfile(run_dir,'overlay02','*.f90')); compiled_mismatches=0;
for k=1:numel(source_files)
    name=source_files(k).name; module='Bellhop'; if strcmp(name,'SourceReceiverPositions.f90'), module='misc'; end
    source=fullfile(source_files(k).folder,name); compiled=fullfile(run_dir,'build01','build','AcousticsToolbox_2020',module,name);
    compiled_mismatches=compiled_mismatches+~strcmpi(black_file_hash(source),black_file_hash(compiled));
end
add('compiled_source_copies_match_reviewed_overlay',compiled_mismatches,0);
add('native_binary_SHA256_still_matches_preflight',double(~strcmpi(black_file_hash(a.s.receipt.executable),a.s.receipt.executable_sha256)),0);
for k=1:numel(a.protected_paths), add(sprintf('protected_before_after_%02d',k),double(~strcmpi(black_file_hash(a.protected_paths{k}),a.before_hash{k})),0); end
phase_difference=angle(pressure_complex.*conj(complex(double(old_pressure(1:2:end,:)),double(old_pressure(2:2:end,:))).'));
receipt=struct('status','SAVED_DATA_AUDIT_COMPLETE','failed_checks',nnz(~checks.pass),'checks_count',height(checks), ...
    'elapsed_s',toc(t),'bellhop_processes_launched',0,'original25_preserved',historical_mismatches==0&&numel(files)==25, ...
    'max_complex_difference',max(abs(pressure_complex-a.field.pressure),[],'all'), ...
    'max_phase_difference_rad',max(abs(phase_difference),[],'all'), ...
    'min_ENV64_gap_m',min(diff(env64)),'min_actual64_gap_m',min(diff(eval64)), ...
    'max_ENV64_unit_roundtrip_m',max(abs(env64-frozen_r)), ...
    'max_SHD_metadata_quantization_m',max(abs(metadata_range-frozen_r)), ...
    'pressure_min_magnitude',min(abs(pressure_complex),[],'all'),'pressure_max_magnitude',max(abs(pressure_complex),[],'all'), ...
    'run_dir',run_dir,'output_dir',output_dir,'next_allowed_step','Manual review only; no full A,order14-16,B/C/D');
writetable(checks,fullfile(output_dir,'checks.csv'));
save(fullfile(output_dir,'audit.mat'),'receipt','checks','physical_binding','component_files','component_hashes');
fid=fopen(fullfile(output_dir,'checkpoint.json'),'w','n','UTF-8'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid);
disp(jsonencode(receipt)); disp(checks(~checks.pass,:));
    function add(name,value,limit)
        checks(end+1,:)={string(name),isfinite(value)&&value<=limit,value,limit};
    end
end
function n=bits(a,b)
assert(numel(a)==numel(b)); n=nnz(typecast(a(:),'uint64')~=typecast(b(:),'uint64'));
end
function [raw,mr,mq]=raw_shd(path,n,nz)
fid=fopen(path,'rb','ieee-le'); assert(fid>=0); cleanup=onCleanup(@()fclose(fid)); rb=4*fread(fid,1,'int32');
fseek(fid,8*rb,'bof'); mq=fread(fid,nz,'float32');
fseek(fid,9*rb,'bof'); mr=fread(fid,n,'float32');
raw=zeros(2*n,nz,'single');
for j=1:nz, fseek(fid,(9+j)*rb,'bof'); raw(:,j)=fread(fid,2*n,'*float32'); end
end