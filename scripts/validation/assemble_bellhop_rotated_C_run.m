function receipt=assemble_bellhop_rotated_C_run(run_dir)
%ASSEMBLE_BELLHOP_ROTATED_A_RESUME Reuse39 saved physical components, no solve.
start=tic; s=load(fullfile(run_dir,'preflight.mat'),'plan','receipt'); plan=s.plan; p=s.receipt;
cpfile=fullfile(run_dir,'checkpoint.json'); state=jsondecode(fileread(cpfile)); assert(state.completed_count==39 && strcmp(state.status,'ALL39_COMPONENTS_COMPLETE_ASSEMBLY_PENDING'));
assert(state.cumulative_accounted_s<3500); [~,ram]=memory; assert(ram.PhysicalMemory.Available>=1.5*2^30,'Assembly RAM gate.');
out=fullfile(run_dir,'assembly'); assert(~isfolder(out)); mkdir(out);
nz=numel(plan.depths_m); nx=numel(plan.x_m); direct=complex(nan(nz,nx)); coverage_count=zeros(nz,nx,'uint8');
surface_family=complex(zeros(nz,nx)); bottom_family=surface_family; total=[]; envelope=[];
component_sources=state.components; order_fields=struct('order',{},'path',{}); mask=~plan.source_singularity_mask;
for j=1:numel(plan.jobs)
    job=plan.jobs(j); source=component_sources(j);
    assert(strcmp(source.tag,job.tag) && strcmpi(black_file_hash(source.path),source.component_sha256));
    assert(strcmpi(black_file_hash(source.checkpoint_path),source.checkpoint_sha256));
    d=load(source.path,'pressure_ffp','cfg'); assert(strcmp(jsonencode(d.cfg),jsonencode(job.cfg))); P=d.pressure_ffp;
    if job.order==0
        sub=direct(job.z_indices,job.x_indices); c=coverage_count(job.z_indices,job.x_indices);
        assert(all(c(job.coverage)==0) && all(isfinite(P(job.coverage))) && all(abs(P(job.coverage))>0));
        sub(job.coverage)=P(job.coverage); c(job.coverage)=c(job.coverage)+1;
        direct(job.z_indices,job.x_indices)=sub; coverage_count(job.z_indices,job.x_indices)=c;
    else
        if isempty(total)
            assert(all(coverage_count(mask)==1) && coverage_count(plan.source_singularity_mask)==0);
            total=direct; envelope=abs(direct);
            save(fullfile(out,'direct.mat'),'direct','coverage_count','plan','-v7.3');
        end
        assert(all(isfinite(P),'all') && all(abs(P)>0,'all'));
        total=total+P; envelope=envelope+abs(P);
        if strcmp(job.cfg.direction,'up'), surface_family=surface_family+P; else, bottom_family=bottom_family+P; end
        if j==numel(plan.jobs) || plan.jobs(j+1).order~=job.order
            n=job.order; path=fullfile(out,sprintf('sum_through_order%d.mat',n)); assert(~isfile(path));
            save(path,'total','envelope','n','plan','-v7.3'); order_fields(end+1)=struct('order',n,'path',path); %#ok<AGROW>
        end
    end
end
assert(isequal(isnan(total),plan.source_singularity_mask));
clear d P
family_error=max(abs((direct(mask)+surface_family(mask)+bottom_family(mask))-total(mask)));
% Different summation association is recorded, not used to replace original order.
receipt=struct('status','FINITE_ORDER_C_COMPLETE_NOT_YET_QUALIFIED','failed_checks',0,'result_dir',out,'run_dir',run_dir, ...
    'elapsed_s',toc(start),'completed_components',39,'reused_components',0,'new_components',39, ...
    'original_A_not_recomputed',true,'family_sum_association_max_difference',family_error, ...
    'family_definition','surface_family=first-surface-launch all orders; bottom_family=first-bottom-launch all orders, mixed reflections included', ...
    'component_sources',component_sources,'order_fields',order_fields,'plan_sha256',p.plan_sha256, ...
    'next_allowed_step','Independent raw SHD reconstruction, C N12-to16 and B-toC gates; D/FFP not authorized');
save(fullfile(out,'result.mat'),'receipt','total','direct','surface_family','bottom_family','envelope','coverage_count','order_fields','plan','component_sources','-v7.3');
receipt.elapsed_s=toc(start); write(fullfile(out,'checkpoint.json'),receipt);
state.assembly_path=out; state.status='ASSEMBLED_C_RAW_AUDIT_PENDING'; state.resume_elapsed_s=state.resume_elapsed_s+receipt.elapsed_s; state.cumulative_accounted_s=state.prior_used_and_reserved_s+state.resume_elapsed_s; write(cpfile,state);
fprintf('C_ASSEMBLED39 elapsed%.3fs cumulative%.3fs\n',receipt.elapsed_s,state.cumulative_accounted_s);
end
function write(path,receipt)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); c=onCleanup(@()fclose(fid)); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true));
end