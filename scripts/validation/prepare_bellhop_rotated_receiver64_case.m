function receipt = prepare_bellhop_rotated_receiver64_case(plan_file,failed_root,overlay_dir,executable,output_dir)
%PREPARE_BELLHOP_ROTATED_RECEIVER64_CASE No-solver coordinate and chain gate.
% Only requested beam/column reductions. ENV range text remains untouched.
arguments
    plan_file (1,:) char
    failed_root (1,:) char
    overlay_dir (1,:) char
    executable (1,:) char
    output_dir (1,:) char
end
assert(~isfolder(output_dir)); assert(isfile(executable)); mkdir(output_dir); t=tic;
p=load(plan_file,'plan'); i=find(strcmp({p.plan.jobs.tag},'up_order13')); assert(isscalar(i)); cfg=p.plan.jobs(i).cfg;
assert(cfg.reflection_order==13 && strcmp(cfg.direction,'up'));
r=cfg.receiver_range_m(:); assert(isa(r,'double') && numel(r)==10981 && all(isfinite(r)) && all(diff(r)>0));
source_env=fileread([failed_root '.env']); lines=splitlines(string(source_env)); nonempty=find(strlength(strtrim(lines))>0);
original_range_line=lines(nonempty(end-4));
fid=fopen([failed_root '.r25d'],'r'); assert(fid>=0); close_side=onCleanup(@()fclose(fid));
flags=fscanf(fid,'%d',2); counts=fscanf(fid,'%d',2); original_r=fscanf(fid,'%f',counts(1)); full_q=fscanf(fid,'%f',counts(2)); clear close_side
assert(counts(1)==10981 && flags(1)==1 && bits_equal(original_r,r));
q=[-1;0;1]; selected=zeros(3,1);
for j=1:3, found=find(full_q==q(j)); assert(isscalar(found)); selected(j)=found; end
q=full_q(selected);
lines(nonempty(end-7))="3"; lines(nonempty(end-6))=string([sprintf('%.17g ',q) '/']);
lines(nonempty(end-2))="101";
allowed_lines=nonempty([end-7 end-6 end-2]); changed=find(lines~=splitlines(string(source_env)));
assert(all(ismember(changed,allowed_lines)) && lines(nonempty(end-4))==original_range_line);
stem=fullfile(output_dir,'case'); write([stem '.env'],char(join(lines,newline)));
for ext={'.ssp','.iwcov','.rotq'}, copyfile([failed_root ext{1}],[stem ext{1}]); end
fid=fopen([stem '.r25d'],'w'); assert(fid>=0); closer=onCleanup(@()fclose(fid));
fprintf(fid,'%d %d\n%d %d\n',flags,counts(1),3); fprintf(fid,'%.17g ',r); fprintf(fid,'\n'); fprintf(fid,'%.17g ',q); fprintf(fid,'\n'); clear closer
fid=fopen([stem '.r25d'],'r'); fscanf(fid,'%d',2); fscanf(fid,'%d',2); planned_r=fscanf(fid,'%f',numel(r)); planned_q=fscanf(fid,'%f',3); fclose(fid);
text_km=sscanf(char(original_range_line),'%f'); env_m64=1000*text_km;
old_metadata=single(single(text_km).*single(1000)); new_metadata=single(env_m64);
reader=fileread(fullfile(overlay_dir,'SourceReceiverPositions.f90'));
state=fileread(fullfile(overlay_dir,'RotatedFieldState.f90'));
environment=fileread(fullfile(overlay_dir,'ReadEnvironmentBell.f90'));
ip=fileread(fullfile(overlay_dir,'influence.f90'));
a=strfind(ip,'  SUBROUTINE InfluenceGeoHatCart('); b=strfind(ip,'  END SUBROUTINE InfluenceGeoHatCart'); hat=ip(a:b);
a=strfind(reader,'    IF(exactInput) THEN'); b=strfind(reader,'    ELSE'); b=b(b>a); exact_reader=reader(a:b(1));
checks=table(strings(0,1),false(0,1),zeros(0,1),zeros(0,1),'VariableNames',{'check','pass','value','limit'});
add('all10981_original_points_retained',abs(numel(planned_r)-10981),0);
add('original_double_all_finite',double(any(~isfinite(r))),0);
add('original_double_strictly_increasing',nnz(diff(r)<=0),0);
add('frozen_MAT_to_original_r25d_bitwise',double(~bits_equal(original_r,r)),0);
add('frozen_MAT_to_planned_r25d_bitwise',double(~bits_equal(planned_r,r)),0);
add('selected_original_columns_bitwise',double(~bits_equal(planned_q,q)),0);
add('no_physical_grid_or_order_change',double(~isequal(p.plan.depths_m,cfg.nominal_receiver_depths_m)),0);
add('ENV_range_text_unchanged',double(lines(nonempty(end-4))~=original_range_line),0);
add('only_beam_and_column_lines_changed',nnz(~ismember(changed,allowed_lines)),0);
add('old_reader_duplicate_pairs_reproduced',abs(nnz(diff(old_metadata)==0)-166),0);
add('new_reader_metadata_strict',nnz(diff(new_metadata)<=0),0);
add('existing_two_ULP_metadata_gate',max(abs(double(new_metadata)-r)./double(eps(new_metadata))),2);
add('ENV64_to_sidecar_conversion_rounding',max(abs(env_m64-r)./eps(env_m64)),2);
add('REAL8_read_and_1000d0_present',double(~(contains(reader,'REAL(KIND=8), ALLOCATABLE :: rangesKm64(:)')&&contains(reader,'ReceiverRangeInput64=1000d0*rangesKm64'))),0);
add('exact_reader_no_sort_or_subtab',double(contains(exact_reader,'CALL Sort')||contains(exact_reader,'CALL SubTab')),0);
add('raw_strict_and_legacy_metadata_strict_retained',double(~(contains(exact_reader,'rangesKm64(2:)<=')&&contains(reader,'.NOT. monotonic( Pos%rr'))),0);
add('early_exact_read_route',double(~contains(environment,'CALL ReadRcvrRanges(ValidationExactRanges)')),0);
add('exact_route_cannot_use_R4_fallback',double(~contains(state,'Exact path cannot use float32 fallback')),0);
add('field_receiver_R8_declaration',double(~contains(hat,'REAL    (KIND=8) :: x_ray')),0);
add('field_no_PosRr_or_PosRz_coordinates',double(contains(hat,'Pos%Rr')||contains(hat,'Pos%Rz')),0);
add('field_and_scale_runtime_bit_audits_enabled',double(~(contains(hat,'CALL RecordExactReceiverEvaluation')&&contains(ip,'CALL RecordExactReceiverScale'))),0);
add('no_single_conversion_of_ReceiverRange64',double(contains(ip,'SNGL(ReceiverRange64')||contains(state,'REAL(ReceiverRange64,KIND=4)')),0);
for ext={'.ssp','.iwcov','.rotq'}
    add(['unchanged_input' ext{1}],double(~strcmpi(black_file_hash([stem ext{1}]),black_file_hash([failed_root ext{1}]))),0);
end
cfg_minimal=cfg; cfg_minimal.beam_count=101; cfg_minimal.nominal_receiver_x_m=[-1 0 1];
receipt=struct('status','NO_SOLVER_PREFLIGHT_COMPLETE','failed_checks',nnz(~checks.pass),'checks_count',height(checks), ...
    'elapsed_s',toc(t),'case_root',stem,'executable',executable,'executable_sha256',black_file_hash(executable), ...
    'plan_file',plan_file,'plan_sha256',black_file_hash(plan_file),'original_failed_root',failed_root, ...
    'original_range_count',numel(r),'beam_count',101,'column_count',3,'q_columns_m',q.', ...
    'physical_x_columns_m',-q.','ENV_range_text_unchanged',true,'original_coordinate_bits_preserved',bits_equal(planned_r,r), ...
    'native_old_duplicates',nnz(diff(old_metadata)==0),'native_new_duplicates',nnz(diff(new_metadata)==0), ...
    'maximum_ENV64_unit_roundtrip_m',max(abs(env_m64-r)),'bellhop_processes_launched',0, ...
    'scope','Gate before the single authorized13th-order101-beam3-column solve; no full A/B/C/D or other order');
writetable(checks,fullfile(output_dir,'preflight_checks.csv'));
save(fullfile(output_dir,'preflight.mat'),'receipt','checks','r','q','selected','env_m64','old_metadata','new_metadata','cfg_minimal');
write(fullfile(output_dir,'checkpoint_preflight.json'),jsonencode(receipt,PrettyPrint=true)); disp(jsonencode(receipt)); disp(checks(~checks.pass,:));
    function add(name,value,limit)
        checks(end+1,:)={string(name),isfinite(value)&&value<=limit,value,limit};
    end
end
function yes=bits_equal(a,b)
yes=isequal(typecast(a(:),'uint64'),typecast(b(:),'uint64'));
end
function write(path,content)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); cleaner=onCleanup(@()fclose(fid)); fprintf(fid,'%s',content);
end