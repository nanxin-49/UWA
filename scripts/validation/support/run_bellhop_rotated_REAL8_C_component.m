function output=run_bellhop_rotated_REAL8_C_component(executable,job,output_dir,max_seconds)
%RUN_BELLHOP_ROTATED_REAL8_B_COMPONENT Formal frozen C job, no label repair.
cfg=job.cfg; assert(~isfolder(output_dir)); [~,ram]=memory; assert(ram.PhysicalMemory.Available>=1.5*2^30,'RAM launch gate');
assert(cfg.ssp_spacing_m==.625);
mkdir(output_dir); root=fullfile(output_dir,'case'); start=tic;
receipt=struct('status','FAILED','cfg',cfg,'result_dir',output_dir,'executable_sha256',black_file_hash(executable),'matlab_release',version,'raw_root',root);
checks=table(strings(0,1),false(0,1),zeros(0,1),zeros(0,1),'VariableNames',{'check','pass','value','limit'});
diary(fullfile(output_dir,'diary.log')); closeDiary=onCleanup(@()diary('off'));
try
    [meta,kind]=write_bellhop_rotated_REAL8_C_inputs(root,cfg);
    [input_proof,input_contract]=check_bellhop_rotated_C_input_contract(root,job.B_raw_root,cfg);
    writetable(input_proof,fullfile(output_dir,'input_contract_checks.csv'));
    write(fullfile(output_dir,'input_contract.json'),jsonencode(input_contract,PrettyPrint=true));
    if isfield(job,'previous_unstarted_input_root')
        for ext={'.env','.ssp','.iwcov','.rotq','.r25d'}
            old=[job.previous_unstarted_input_root ext{1}]; current=[root ext{1}]; assert(isfile(old)==isfile(current));
            if isfile(old), assert(strcmpi(black_file_hash(old),black_file_hash(current)),'Recovery changed previously frozen C input'); end
        end
    end
    [~,launchRAM]=memory; assert(launchRAM.PhysicalMemory.Available>=1.5*2^30,'Immediate native launch RAM gate');
    process=run_bellhop_rotated_local_case(executable,root,min(120,max_seconds),allow_task_termination=true,minimum_runtime_available_bytes=256*2^20);
    d=readmatrix([root '.r25diag'],'FileType','text','CommentStyle','#'); d=d(~all(isnan(d),2),:);
    if strcmp(kind,'reflection')
        wall=readmatrix([root '.mwdiag'],'FileType','text','CommentStyle','#'); wall=wall(~all(isnan(wall),2),:);
        add('all_requested_rays',abs(size(d,1)-cfg.beam_count)+abs(size(wall,1)-cfg.beam_count),0);
        add('exact_physical_bounce_counts',nnz(d(:,9)~=cfg.expected_surface_bounces | d(:,10)~=cfg.expected_bottom_bounces),0);
        add('requested_wall_order',nnz(wall(:,3)~=cfg.reflection_order),0);
        add('capacity_original_strict_safety_margin',double(max(wall(:,6))>=.95*600000),0);
    else
        wall=[]; add('all_requested_rays',abs(size(d,1)-cfg.beam_count),0);
        add('no_artificial_boundary_reflections',nnz(d(:,9:10)),0);
    end
    add('all_receivers_reached',max(cfg.receiver_range_m)-min(d(:,2)),0);
    add('end_unit_tangent_original_gate',max(d(:,8)),1e-8);
    raw=read_bellhop_rotated_field(root); r=cfg.receiver_range_m(:); q=meta.q_receiver_m(:); nr=numel(r); nq=numel(q);
    add('all_native_pressure_finite',nnz(~isfinite(raw.pressure)),0);
    required=true(size(raw.pressure));
    if job.order==0
        required=job.coverage(cfg.receiver_depth_order,cfg.receiver_x_order);
        if ~strcmp(kind,'horizontal'), required=required.'; end
    end
    add('finite_nonzero_required_chart_coverage',nnz(abs(raw.pressure(required))==0),0);
    fid=fopen([root '.rr64input'],'rb','ieee-le'); assert(fid>=0); h=fread(fid,4,'int32'); env64=fread(fid,nr,'double'); exact64=fread(fid,nr,'double'); q64=fread(fid,nq,'double'); fclose(fid);
    add('REAL8_header',double(~isequal(h,[nr;nq;8;8])),0);
    add('frozen_effective_range_bits',bitdiff(exact64,r),0); add('frozen_effective_transverse_bits',bitdiff(q64,q),0);
    add('strict_real8_range_order',nnz(diff(exact64)<=0),0); add('ENV64_unit_roundoff_existing_gate',max(abs(env64-r)./eps(env64)),2);
    fid=fopen([root '.rr64eval'],'rb','ieee-le'); h=fread(fid,2,'int32'); assert(isequal(h,[nr;nq]));
    seen=false(nr*nq,1); count=0; mismatches=0;
    while true
        bytes=fread(fid,24*32768,'*uint8'); if isempty(bytes), break; end; assert(mod(numel(bytes),24)==0);
        rec=reshape(bytes,24,[]); ir=double(typecast(reshape(rec(1:4,:),[],1),'int32')); iq=double(typecast(reshape(rec(5:8,:),[],1),'int32'));
        ar=typecast(reshape(rec(9:16,:),[],1),'double'); aq=typecast(reshape(rec(17:24,:),[],1),'double');
        assert(all(ir>=1 & ir<=nr & iq>=1 & iq<=nq)); mismatches=mismatches+bitdiff(ar,r(ir))+bitdiff(aq,q(iq));
        seen((ir-1)*nq+iq)=true; count=count+numel(ir);
    end
    fclose(fid); add('all_actual_evaluation_bits',mismatches,0); add('evaluation_records_unique',count-nnz(seen),0); add('required_evaluation_pair_coverage',nnz(required(:) & ~seen),0); clear seen bytes rec
    fid=fopen([root '.rr64scale'],'rb','ieee-le'); hs=fread(fid,1,'int32'); bytes=fread(fid,inf,'*uint8'); fclose(fid);
    assert(hs==nr && numel(bytes)==12*nr); rec=reshape(bytes,12,[]); ir=double(typecast(reshape(rec(1:4,:),[],1),'int32')); ar=typecast(reshape(rec(5:12,:),[],1),'double');
    add('all_scale_coordinate_bits',bitdiff(ar,r(ir)),0); add('all_scale_indices_once',nnz(accumarray(ir,1,[nr 1])~=1),0);
    add('stderr_empty',strlength(string(process.stderr)),0); assert(all(checks.pass),'Formal component numerical gate failed.');
    pressure_ffp=complex(zeros(numel(cfg.nominal_receiver_depths_m),numel(cfg.nominal_receiver_x_m)));
    if strcmp(kind,'horizontal')
        pressure_ffp(cfg.receiver_depth_order,cfg.receiver_x_order)=-conj(raw.pressure);
    else
        pressure_ffp(cfg.receiver_depth_order,cfg.receiver_x_order)=(-conj(raw.pressure)).';
    end
    receipt.status='COMPLETE'; receipt.failed_checks=0; receipt.process=process; receipt.elapsed_s=toc(start); receipt.kind=kind;
    receipt.max_ray_end_unit_error=max(d(:,8)); receipt.max_ray_range_m=max(d(:,2)); receipt.available_before_bytes=ram.PhysicalMemory.Available;
    code={'run_bellhop_rotated_REAL8_C_component','check_bellhop_rotated_C_input_contract','write_bellhop_rotated_REAL8_C_inputs','write_bellhop_rotated_munk_q','write_bellhop_rotated_field_options','read_bellhop_rotated_field','run_bellhop_rotated_local_case'};
    receipt.code_versions=struct('path',{},'sha256',{}); for k=1:numel(code), path=which(code{k}); receipt.code_versions(k)=struct('path',path,'sha256',black_file_hash(path)); end
    save(fullfile(output_dir,'component.mat'),'pressure_ffp','cfg','meta','receipt','d','wall','-v7.3');
    bindings=struct('path',{},'sha256',{});
    for ext={'.env','.ssp','.iwcov','.rotq','.r25d','.shd','.prt','.r25diag','.mwdiag','.covdiag','.rr64input','.rr64eval','.rr64scale'}
        path=[root ext{1}]; if ~isfile(path), continue; end; bindings(end+1)=struct('path',path,'sha256',black_file_hash(path)); %#ok<AGROW>
    end
    path=fullfile(output_dir,'component.mat'); bindings(end+1)=struct('path',path,'sha256',black_file_hash(path)); receipt.bindings=bindings;
    output=struct('pressure_ffp',pressure_ffp,'receipt',receipt,'meta',meta);
catch exception
    receipt.status='FAILED';
    receipt.elapsed_s=toc(start); receipt.failed_checks=nnz(~checks.pass)+1; receipt.error=getReport(exception,'extended','hyperlinks','off');
    writetable(checks,fullfile(output_dir,'checks.csv')); write(fullfile(output_dir,'checkpoint.json'),jsonencode(receipt,PrettyPrint=true)); rethrow(exception)
end
writetable(checks,fullfile(output_dir,'checks.csv')); write(fullfile(output_dir,'checkpoint.json'),jsonencode(receipt,PrettyPrint=true));
fprintf('FORMAL_C %s order%d points%d native%.4fs audit_total%.4fs failed%d\n',job.tag,job.order,numel(pressure_ffp),process.elapsed_s,receipt.elapsed_s,receipt.failed_checks);
    function add(name,value,limit)
        checks(end+1,:)={string(name),isfinite(value)&&value<=limit,value,limit};
    end
end
function n=bitdiff(a,b)
assert(numel(a)==numel(b)); n=nnz(typecast(a(:),'uint64')~=typecast(b(:),'uint64'));
end
function write(path,text)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); c=onCleanup(@()fclose(fid)); fprintf(fid,'%s',text);
end