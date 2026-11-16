function receipt = audit_bellhop_rotated_A_resume_raw(input_dir,output_dir,options)
%AUDIT_BELLHOP_ROTATED_FIELD_ASSEMBLY Reconstruct every chart from raw SHD.
arguments
    input_dir (1,:) char
    output_dir (1,:) char
    options.ffp_file (1,:) char=''
    options.ffp_expected_sha256 (1,:) char=''
end
assert(~isfolder(output_dir)); mkdir(output_dir); t=tic;
saved=load(fullfile(input_dir,'result.mat'),'plan','direct','total','order_fields','envelope','component_sources'); plan=saved.plan;
checks=table(strings(0,1),false(0,1),zeros(0,1),zeros(0,1), ...
    'VariableNames',{'check','pass','value','limit'});
P0=complex(nan(size(saved.direct))); P=P0; envelope=zeros(size(P0)); counts=zeros(size(P0));
bindings=struct('path',{},'sha256',{});
for j=1:numel(plan.jobs)
    job=plan.jobs(j); root=fullfile(input_dir,job.tag,'case');
    component_path=saved.component_sources(j).path;
    component=load(component_path,'pressure_ffp','receipt');
    if isfield(component.receipt,'raw_root'), root=component.receipt.raw_root; end
    cpfile=saved.component_sources(j).checkpoint_path;
    if isfile(cpfile)
        jobcp=jsondecode(fileread(cpfile));
        if isfield(jobcp,'bindings')
            for h=1:numel(jobcp.bindings)
                assert(strcmpi(black_file_hash(jobcp.bindings(h).path),jobcp.bindings(h).sha256),'Component artifact SHA drift.');
            end
        end
    end
    raw=read_bellhop_rotated_field(root);
    cfg=job.cfg; A=complex(zeros(numel(job.z_indices),numel(job.x_indices)));
    if isfield(cfg,'chart_kind') && strcmp(cfg.chart_kind,'horizontal')
        A(cfg.receiver_depth_order,cfg.receiver_x_order)=-conj(raw.pressure);
    else
        A(cfg.receiver_depth_order,cfg.receiver_x_order)=(-conj(raw.pressure)).';
    end
    add([job.tag '_raw_physical_binding'],max(abs(component.pressure_ffp-A),[],'all'),0);
    for ext={'.r25d','.shd','.env','.ssp','.iwcov'}
        path=[root ext{1}]; bindings(end+1)=struct('path',path,'sha256',black_file_hash(path)); %#ok<AGROW>
    end
    if job.order==0
        sub=P0(job.z_indices,job.x_indices); sub(job.coverage)=A(job.coverage);
        P0(job.z_indices,job.x_indices)=sub;
        counts(job.z_indices,job.x_indices)=counts(job.z_indices,job.x_indices)+job.coverage;
        P=P0; envelope=abs(P0);
    else
        P=P+A; envelope=envelope+abs(A);
        if j==numel(plan.jobs) || plan.jobs(j+1).order~=job.order
            partial=load(fullfile(input_dir,sprintf('sum_through_order%d.mat',job.order)),'total');
            valid=~plan.source_singularity_mask;
            add(sprintf('sum_order%d_raw_complex',job.order),max(abs(P(valid)-partial.total(valid))),0);
        end
    end
end
valid=~plan.source_singularity_mask;
add('direct_chart_coverage_exactly_once',max(abs(counts(valid)-1)),0);
add('source_coverage_zero',sum(counts(plan.source_singularity_mask)),0);
add('only_true_source_mask',double(~isequal(isnan(P),plan.source_singularity_mask)),0);
add('direct_reconstruction',max(abs(P0(valid)-saved.direct(valid))),0);
add('total_reconstruction',max(abs(P(valid)-saved.total(valid))),0);
if isfield(saved,'envelope'), add('envelope_reconstruction',max(abs(envelope(valid)-saved.envelope(valid))),0); end
comparison=struct('status','NOT_REQUESTED');
if ~isempty(options.ffp_file)
    [~,mem]=memory; assert(mem.PhysicalMemory.Available>=2^30,'Reference read RAM gate.');
    sha=black_file_hash(options.ffp_file);
    assert(~isempty(options.ffp_expected_sha256) && strcmpi(sha,options.ffp_expected_sha256),'FFP SHA binding failed.');
    ffp=load(options.ffp_file,'output');
    [dz,iz]=min(abs(ffp.output.depth_m(:)-plan.depths_m),[],1);
    [dx,ix]=min(abs(ffp.output.range_m(:)-plan.x_m),[],1);
    assert(all(dz<1e-9) && all(dx<1e-9),'No interpolation of unmatched FFP points.');
    reference=ffp.output.pressure_1m(iz,ix); mask=valid & isfinite(reference);
    comparison=struct('status','DIAGNOSTIC_ONLY_PENDING_CONVERGENCE', ...
        'ffp_sha256',sha,'coordinate_max_difference_m',max([dz dx]), ...
        'complex_relative_L2',norm(P(mask)-reference(mask))/norm(reference(mask)), ...
        'maximum_envelope_normalized_difference',max(abs(P(mask)-reference(mask))./envelope(mask)), ...
        'reflection_order',plan.settings.max_order);
    save(fullfile(output_dir,'reference_subset.mat'),'reference','comparison','envelope','plan');
end
writetable(checks,fullfile(output_dir,'checks.csv'));
receipt=struct('status','COMPLETE_RAW_RECONSTRUCTION_ONLY','failed_checks',nnz(~checks.pass), ...
    'checks_count',height(checks),'elapsed_s',toc(t),'result_dir',output_dir,'comparison',comparison);
save(fullfile(output_dir,'result.mat'),'receipt','checks','bindings','envelope');
total=P; direct=P0; save(fullfile(output_dir,'reconstructed_field.mat'),'total','direct','envelope','plan','-v7.3');
write(fullfile(output_dir,'checkpoint.json'),jsonencode(receipt,PrettyPrint=true));
write(fullfile(output_dir,'input_sha256.json'),jsonencode(bindings,PrettyPrint=true));
disp(jsonencode(receipt));
    function add(name,value,limit)
        checks(end+1,:)={string(name),isfinite(value)&&value<=limit,value,limit};
    end
end
function write(path,content)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); cleaner=onCleanup(@()fclose(fid)); fprintf(fid,'%s',content);
end