function receipt=audit_d1_saved_tx100(origin,out)
%AUDIT_D1_SAVED_TX100 Full saved five-layer audit and fixed FIG pipeline, no solve.
arguments
    origin (1,:) char
    out (1,:) char
end
assert(~isfolder(out)); mkdir(out); diary(fullfile(out,'diary.log')); cleanup=onCleanup(@()diary('off')); t0=tic;
previous_threads=maxNumCompThreads; maxNumCompThreads(1); restore_threads=onCleanup(@()maxNumCompThreads(previous_threads));
rows=jsondecode(fileread(fullfile(origin,'retrieval','manifest.json'))); pm_sensitivity_io('verify',origin,rows);
fp=pm_sensitivity_io('sha',fullfile(origin,'manifest.json'));
assert(strcmp(fp,'6221c13050770c38392be06d9675328673585a96000b5ba00e399235441ac727'));
cp=jsondecode(fileread(fullfile(origin,'output','checkpoint.json')));
assert(strcmp(cp.status,'FAILED_STOPPED') && cp.new_BIE_solves==4 && numel(cp.completed_steps)==5);
assert(contains(strjoin(string(cp.failed_or_unrun_checks),' '),'Reopened PNG pixel mismatch'));
src=fullfile(origin,'output','depth_0100m'); work=fullfile(out,'audited_work'); mkdir(work);
for name={'common.mat','comparison.mat','model_metrics.csv','receiver_fields_and_errors.csv','policy.json','source_geometry_audit.json','source_geometry_audit.mat'}
    copyfile(fullfile(src,name{1}),fullfile(work,name{1}));
end
copyfile(fullfile(src,'input'),fullfile(work,'input')); copyfile(fullfile(src,'jobs'),fullfile(work,'jobs'));
ids={'bie_coarse','bie_main','bie_self','bie_space','bie_window'}; recipes=[60 10 96;60 12 96;60 12 64;60 14 96;70 12 96];
for j=1:5
    stage=fullfile(work,'jobs',ids{j}); seal=jsondecode(fileread(fullfile(stage,'checkpoint.json')));
    assert(strcmp(seal.status,'AUDITED_PASS')); pm_sensitivity_io('verify',stage,seal.artifacts_sha256);
    original_seal_sha=pm_sensitivity_io('sha',fullfile(stage,'checkpoint.json'));
    audit_r1_bie_manual_stage(stage,fullfile(work,'common.mat'),recipes(j,:));
    seal.origin_checkpoint_sha256=original_seal_sha; seal.audit_execution='local saved-density re-audit; original preserved';
    seal.artifacts_sha256=pm_sensitivity_io('manifest',stage,{seal.artifacts_sha256.path});
    pm_sensitivity_io('write',fullfile(stage,'checkpoint.json'),seal);
end
audit_d1_models(work); diary(fullfile(out,'diary.log')); audit_d1_final(work);
render_d1_depth_results(work); diary(fullfile(out,'diary.log'));
fig=jsondecode(fileread(fullfile(work,'figures','figure_audit.json'))); assert(fig.failed_checks==0 && fig.reopened_curves==24);
pm_sensitivity_io('verify',origin,rows);
receipt=struct('status','SAVED_TX100_DEPTH_AUDIT_PASS','failed_checks',0,'elapsed_s',toc(t0),'new_BIE_solves',0, ...
    'origin_manifest_sha256',fp,'origin_retrieval_sha_files',numel(rows),'completed_BIE_layers',5, ...
    'remaining_BIE_solves',5,'figure_audit',fig,'audited_work_dir',work,'matlab_version',version, ...
    'next_allowed_step','User manually starts new package importing all five Tx100 layers; no Tx100 re-solve');
pm_sensitivity_io('write',fullfile(out,'receipt.json'),receipt);
listing=dir(fullfile(work,'**','*')); paths={'receipt.json'};
for j=1:numel(listing)
    if ~listing(j).isdir
        p=fullfile(listing(j).folder,listing(j).name); paths{end+1}=strrep(p(numel(out)+2:end),'\','/'); %#ok<AGROW>
    end
end
repo=fileparts(fileparts(fileparts(mfilename('fullpath'))));
pm_sensitivity_io('write',fullfile(out,'checkpoint.json'),struct('status','AUDITED_PASS', ...
    'code_versions_sha256',pm_sensitivity_io('manifest',repo,{'scripts/reporting/render_d1_depth_results.m','scripts/validation/audit_d1_saved_tx100.m'}), ...
    'artifacts_sha256',pm_sensitivity_io('manifest',out,paths),'failed_or_unrun_checks',{{'Tx1490 and native Linux continuation not run'}}, ...
    'resource_usage',struct('new_BIE_solves',0,'elapsed_s',toc(t0)),'next_allowed_step',receipt.next_allowed_step));
fprintf('%s\n',jsonencode(receipt));
end
