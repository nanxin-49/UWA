function receipt = run_bellhop_rotated_field_suite(suite_file,executable,output_dir,approval_file,options)
%RUN_BELLHOP_ROTATED_FIELD_SUITE Approved local four-control streaming suite.
% Shares one wall-time/disk budget; each failed numerical gate stops dispatch.
arguments
    suite_file (1,:) char
    executable (1,:) char
    output_dir (1,:) char
    approval_file (1,:) char
    options.reuse_dir (1,:) char=''
end
assert(~isfolder(output_dir));
authorization=jsondecode(fileread(approval_file));
assert(authorization.user_approved && strcmp(authorization.scope,'rotated_fig12b_full_local'));
assert(strcmpi(authorization.suite_sha256,black_file_hash(suite_file)),'Suite approval SHA mismatch.');
s=load(suite_file,'suite'); suite=s.suite;
assert(strcmpi(black_file_hash(executable),suite.executable_sha256),'Binary drift.');
for k=1:numel(suite.code_versions)
    assert(strcmpi(black_file_hash(suite.code_versions(k).path),suite.code_versions(k).sha256),'Execution code drift.');
end
for k=1:numel(suite.plan_files)
    assert(any(strcmpi(string(authorization.allowed_plan_sha256),black_file_hash(suite.plan_files{k}))));
end
assert(numel(suite.plan_files)==4,'Frozen single-variable controls require four plans.');
[~,ram]=memory; assert(ram.PhysicalMemory.Available>=authorization.minimum_available_bytes,'Suite RAM gate.');
parent=fileparts(output_dir); assert(isfolder(parent)); drive=java.io.File(parent);
assert(drive.getUsableSpace()>=authorization.minimum_free_disk_bytes,'Suite disk startup gate.');
used_before_s=0; disk_roots={output_dir};
if ~isempty(options.reuse_dir)
    old=jsondecode(fileread(fullfile(options.reuse_dir,'checkpoint.json')));
    assert(strcmpi(old.suite_sha256,black_file_hash(suite_file)),'Suite reuse revision differs.');
    assert(~strcmp(old.status,'RUNNING'),'Prior suite may still be active; inspect before resume.');
    used_before_s=old.cumulative_elapsed_s;
    disk_roots=unique([cellstr(string(old.disk_roots(:)));{output_dir}]);
    assert(used_before_s<authorization.max_wall_s,'Original budget exhausted; no automatic budget increase.');
end
mkdir(output_dir); t=tic; results=cell(1,4); audits=cell(1,4); changes=cell(1,3); tails=cell(1,4);
receipt=struct('status','RUNNING','result_dir',output_dir,'suite_file',suite_file, ...
    'suite_sha256',black_file_hash(suite_file),'used_before_s',used_before_s, ...
    'approval_file',approval_file,'approval_sha256',black_file_hash(approval_file),'matlab_release',version);
receipt.cumulative_elapsed_s=used_before_s; receipt.disk_roots=disk_roots;
write(fullfile(output_dir,'checkpoint.json'),jsonencode(receipt,PrettyPrint=true));
try
    for k=1:4
        remaining=authorization.max_wall_s-used_before_s-toc(t)-60;
        assert(remaining>0,'Shared suite wall-time gate; reserve 60s for numerical audits.');
        out=fullfile(output_dir,suite.labels{k}); prior_dir=''; complete_prior=false;
        if ~isempty(options.reuse_dir)
            candidate=fullfile(options.reuse_dir,suite.labels{k},'checkpoint.json');
            if isfile(candidate)
                prior_dir=fileparts(candidate); previous=jsondecode(fileread(candidate));
                assert(strcmpi(previous.plan_sha256,black_file_hash(suite.plan_files{k})),'Prior plan differs.');
                complete_prior=strcmp(previous.status,'COMPLETE_FINITE_ORDER_FIELD_NOT_QUALIFIED');
            end
        end
        if complete_prior
            out=prior_dir; results{k}=previous;
        else
            results{k}=run_bellhop_rotated_field_assembly(suite.plan_files{k},executable,out, ...
                scope='approved_full_local',approval_file=approval_file,max_wall_s=remaining, ...
                max_job_s=authorization.max_job_s,minimum_available_bytes=authorization.minimum_available_bytes,reuse_dir=prior_dir);
        end
        tails{k}=assess_bellhop_rotated_field_change(fullfile(out,'sum_through_order12.mat'),fullfile(out,'result.mat'), ...
            fullfile(output_dir,[suite.labels{k} '_order_convergence']),label='N12_to_N16', ...
            relative_tolerance=.001,envelope_tolerance=.001,phase_tolerance_rad=.01);
        assert(tails{k}.failed_checks==0,'Reflection-order convergence failed.');
        audits{k}=audit_bellhop_rotated_field_assembly(out,fullfile(output_dir,[suite.labels{k} '_raw_audit']));
        assert(audits{k}.failed_checks==0,'Raw reconstruction failed.');
        if k>1
            previous=fullfile(results{k-1}.result_dir,'result.mat');
            changes{k-1}=assess_bellhop_rotated_field_change(previous,fullfile(out,'result.mat'), ...
                fullfile(output_dir,[suite.labels{k} '_change']),label=suite.control_labels{k-1});
            assert(changes{k-1}.failed_checks==0,'Single-variable numerical convergence failed.');
        end
        disk_bytes=0;
        for d=1:numel(disk_roots)
            used=dir(fullfile(disk_roots{d},'**','*')); used=used(~[used.isdir]); disk_bytes=disk_bytes+sum([used.bytes]);
        end
        assert(disk_bytes<=authorization.max_output_bytes,'Suite output-size gate.');
        assert(used_before_s+toc(t)<=authorization.max_wall_s,'Suite shared wall-time gate.');
        receipt.completed_controls=k; receipt.elapsed_s=toc(t); receipt.cumulative_elapsed_s=used_before_s+toc(t); receipt.output_bytes=disk_bytes;
        receipt.results=results; receipt.raw_audits=audits; receipt.order_audits=tails; receipt.changes=changes;
        write(fullfile(output_dir,'checkpoint.json'),jsonencode(receipt,PrettyPrint=true));
    end
    receipt.status='NUMERICAL_SUITE_COMPLETE_FFP_AND_FIGURES_PENDING';
    receipt.failed_checks=0; receipt.elapsed_s=toc(t); receipt.cumulative_elapsed_s=used_before_s+toc(t);
catch exception
    receipt.status='FAILED';
    receipt.elapsed_s=toc(t); receipt.cumulative_elapsed_s=used_before_s+toc(t); receipt.error=getReport(exception,'extended','hyperlinks','off');
    write(fullfile(output_dir,'checkpoint.json'),jsonencode(receipt,PrettyPrint=true)); rethrow(exception)
end
write(fullfile(output_dir,'checkpoint.json'),jsonencode(receipt,PrettyPrint=true));
fprintf('SUITE %s elapsed=%.3fs\n',receipt.status,receipt.elapsed_s);
end
function write(path,content)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); cleaner=onCleanup(@()fclose(fid)); fprintf(fid,'%s',content);
end