function receipt = run_bellhop_rotated_field_assembly(plan_file,executable,output_dir,options)
%RUN_BELLHOP_ROTATED_FIELD_ASSEMBLY Streaming coherent exact-grid assembly.
% Full mode requires a durable user-approval receipt bound to the plan SHA.
arguments
    plan_file (1,:) char
    executable (1,:) char
    output_dir (1,:) char
    options.scope (1,:) char {mustBeMember(options.scope,{'lightweight_assembly_pilot','approved_full_local'})}='lightweight_assembly_pilot'
    options.reuse_dir (1,:) char=''
    options.approval_file (1,:) char=''
    options.max_wall_s (1,1) double {mustBePositive}=60
    options.max_job_s (1,1) double {mustBePositive}=20
    options.minimum_available_bytes (1,1) double {mustBePositive}=512*2^20
end
assert(~isfolder(output_dir)); data=load(plan_file,'plan'); plan=data.plan;
plan_sha=black_file_hash(plan_file); allow_termination=false; runtime_ram=0;
if strcmp(options.scope,'lightweight_assembly_pilot')
    assert(numel(plan.depths_m)*numel(plan.x_m)<=10000 && plan.settings.max_order<=3 && ...
        options.max_wall_s<=60 && options.max_job_s<=20,'Pilot scope cannot authorize a complete-grid run.');
else
    assert(isfile(options.approval_file),'Full local run needs explicit durable user approval.');
    approval=jsondecode(fileread(options.approval_file));
    assert(strcmp(approval.scope,'rotated_fig12b_full_local') && approval.user_approved && ...
        any(strcmpi(string(approval.allowed_plan_sha256),plan_sha)), 'Plan is not covered by user approval.');
    assert(options.max_wall_s<=approval.max_wall_s && options.max_job_s<=approval.max_job_s && ...
        options.minimum_available_bytes>=approval.minimum_available_bytes,'Approval resource bounds exceeded.');
    allow_termination=approval.allow_task_termination; runtime_ram=approval.minimum_runtime_available_bytes;
end
[~,mem]=memory; assert(mem.PhysicalMemory.Available>=options.minimum_available_bytes,'RAM gate.');
if ~isempty(options.reuse_dir)
    prior=jsondecode(fileread(fullfile(options.reuse_dir,'checkpoint.json')));
    assert(strcmpi(prior.plan_sha256,plan_sha) && strcmpi(prior.executable_sha256,black_file_hash(executable)),'Assembly reuse contract drift.');
end
mkdir(output_dir); t=tic; nz=numel(plan.depths_m); nx=numel(plan.x_m);
direct=complex(nan(nz,nx)); total=direct; envelope=abs(direct); coverage_count=zeros(nz,nx,'uint8');
runs=struct('tag',{},'path',{},'elapsed_s',{},'order',{}); order_fields=struct('order',{},'path',{});
receipt=struct('status','FAILED','scope',options.scope,'result_dir',output_dir, ...
    'plan_file',plan_file,'plan_sha256',plan_sha,'executable_sha256',black_file_hash(executable), ...
    'matlab_release',version,'options',options);
try
    for j=1:numel(plan.jobs)
        job=plan.jobs(j); remain=options.max_wall_s-toc(t);
        assert(remain>0,'Assembly cumulative time gate.');
        d=fullfile(output_dir,job.tag);
        reuse='';
        if ~isempty(options.reuse_dir)
            candidate=fullfile(options.reuse_dir,job.tag,'checkpoint.json');
            if isfile(candidate)
                prior=jsondecode(fileread(candidate));
                if strcmp(prior.status,'COMPLETE')
                    reuse=fileparts(candidate);
                elseif isfile(fullfile(fileparts(candidate),'pid.txt'))
                    state=fullfile(fileparts(candidate),'state.txt');
                    assert(isfile(state) && startsWith(strtrim(fileread(state)),'TERMINAL'), ...
                        'Prior task PID may still be live; inspect before retry.');
                end
            end
        end
        result=run_bellhop_rotated_field_component(executable,job.cfg,d, ...
            max_seconds=min(options.max_job_s,remain),minimum_available_bytes=options.minimum_available_bytes, ...
            allow_task_termination=allow_termination,minimum_runtime_available_bytes=runtime_ram,reuse_dir=reuse);
        P=result.pressure_ffp;
        if job.order==0
            ix=job.x_indices; iz=job.z_indices; mask=job.coverage;
            c=coverage_count(iz,ix); assert(all(c(mask)==0),'Direct chart double assignment.');
            assert(all(isfinite(P(mask))) && all(abs(P(mask))>0),'Direct chart coverage missing.');
            sub=direct(iz,ix); sub(mask)=P(mask); direct(iz,ix)=sub;
            c(mask)=c(mask)+1; coverage_count(iz,ix)=c;
        else
            if ~any([runs.order]>0)
                assert(all(coverage_count(~plan.source_singularity_mask)==1),'Direct coverage gap.');
                assert(all(coverage_count(plan.source_singularity_mask)==0),'True source singularity incorrectly filled.');
                total=direct; envelope=abs(direct);
                save(fullfile(output_dir,'direct.mat'),'direct','coverage_count','plan','-v7.3');
            end
            assert(all(isfinite(P),'all') && all(abs(P)>0,'all'),'Reflected component coverage missing.');
            total=total+P; envelope=envelope+abs(P);
            if j==numel(plan.jobs) || plan.jobs(j+1).order~=job.order
                n=job.order; path=fullfile(output_dir,sprintf('sum_through_order%d.mat',n));
                save(path,'total','envelope','n','plan','-v7.3');
                order_fields(end+1)=struct('order',n,'path',path); %#ok<AGROW>
            end
        end
        runs(end+1)=struct('tag',job.tag,'path',d,'elapsed_s',result.receipt.elapsed_s,'order',job.order); %#ok<AGROW>
        clear result P
        receipt.completed_jobs=numel(runs); receipt.elapsed_s=toc(t); receipt.runs=runs;
        write(fullfile(output_dir,'checkpoint.json'),jsonencode(receipt,PrettyPrint=true));
    end
    assert(isequal(isnan(total),plan.source_singularity_mask),'Only the true source point may be masked.');
    receipt.status='COMPLETE_FINITE_ORDER_FIELD_NOT_QUALIFIED'; receipt.failed_checks=0;
    receipt.elapsed_s=toc(t); receipt.runs=runs; receipt.order_fields=order_fields;
    receipt.unrun_checks={'independent field convergence','FFP comparison','raw-data/figure reconstruction audits'};
    save(fullfile(output_dir,'result.mat'),'receipt','total','direct','envelope','coverage_count','order_fields','plan','-v7.3');
catch exception
    receipt.elapsed_s=toc(t); receipt.error=getReport(exception,'extended','hyperlinks','off');
    write(fullfile(output_dir,'checkpoint.json'),jsonencode(receipt,PrettyPrint=true)); rethrow(exception)
end
write(fullfile(output_dir,'checkpoint.json'),jsonencode(receipt,PrettyPrint=true)); disp(jsonencode(rmfield(receipt,{'runs','order_fields'})));
end
function write(path,content)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); cleaner=onCleanup(@()fclose(fid)); fprintf(fid,'%s',content);
end