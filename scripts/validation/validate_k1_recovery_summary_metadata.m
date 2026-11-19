function receipt=validate_k1_recovery_summary_metadata(master,out)
%VALIDATE_K1_RECOVERY_SUMMARY_METADATA Test-only seals/arrays, renderer mocked.
% Does not consume a graphics attempt: no figures or PNGs are constructed.
assert(~isfolder(out)); mkdir(out); timer=tic;
oldpath=path; pg=onCleanup(@()path(oldpath)); addpath(fullfile(master,'code'),'-begin');
fixture=fullfile(out,'TEST_ONLY_package'); mkdir(fixture);
copyfile(fullfile(master,'anchor'),fullfile(fixture,'anchor'));
cfg=k1_pm_io('read',fullfile(master,'config.json')); cfg.kind='TEST_ONLY_K1_SUMMARY';
k1_pm_io('write',fullfile(fixture,'config.json'),cfg);
ids={'bie_coarse','bie_main','bie_self','bie_space','bie_window'};
for j=1:6
    source=fullfile(master,'cases',cfg.cases(j).case_id);
    case_root=fullfile(fixture,'cases',cfg.cases(j).case_id); mkdir(case_root);
    for name={'code','input'}, copyfile(fullfile(source,name{1}),fullfile(case_root,name{1})); end
    for name={'policy.json','manifest.json'}, copyfile(fullfile(source,name{1}),fullfile(case_root,name{1})); end
    fp=k1_pm_io('sha',fullfile(case_root,'manifest.json'));
    model=fullfile(case_root,'model_stage'); mkdir(fullfile(model,'jobs','models'));
    copyfile(fullfile(fixture,'anchor','common.mat'),fullfile(model,'common.mat'));
    k1_pm_io('write',fullfile(model,'jobs','models','receipt.json'),struct('status','COMPLETE','failed_checks',0,'test_only',true));
    k1_pm_io('write',fullfile(model,'solve_invocation.json'),struct('invocations',1,'stage','MODELS','test_only',true));
    sig=[fp ':MODELS:' fp];
    k1_pm_io('seal',model,'solved.json',sig,{'common.mat','jobs/models/receipt.json','solve_invocation.json'},'SOLVED_UNAUDITED');
    k1_pm_io('seal',model,'committed.json',sig,{'common.mat','jobs/models/receipt.json','solved.json','solve_invocation.json'},'AUDITED_PASS');
    pred=k1_pm_io('sha',fullfile(model,'committed.json'));
    for k=1:5
        dest=fullfile(case_root,'jobs',ids{k}); mkdir(dest); synthetic_payload='TEST_ONLY_NO_DENSITY'; save(fullfile(dest,'data.mat'),'synthetic_payload','-v7');
        k1_pm_io('write',fullfile(dest,'audit.json'),struct('status','PASS','failed_checks',0,'test_only',true));
        k1_pm_io('write',fullfile(dest,'gates.json'),struct('status','PASS','test_only',true));
        k1_pm_io('write',fullfile(dest,'solve_invocation.json'),struct('invocations',1,'test_only',true));
        sig=[pred ':C' num2str(k-1) ':' fp];
        k1_pm_io('seal',dest,'committed.json',sig,{'data.mat','audit.json','gates.json','solve_invocation.json'},'AUDITED_PASS');
        pred=k1_pm_io('sha',fullfile(dest,'committed.json'));
    end
    final=fullfile(case_root,'finalization'); work=fullfile(final,'attempts','attempt_0001','work'); mkdir(work);
    for name={'common.mat','comparison.mat','model_metrics.csv','receiver_fields_and_errors.csv','independent_final_metrics_audit.json'}
        copyfile(fullfile(fixture,'anchor',name{1}),fullfile(work,name{1}));
    end
    mkdir(fullfile(work,'audit')); copyfile(fullfile(fixture,'anchor','audit','receipt.json'),fullfile(work,'audit','receipt.json'));
    k1_pm_io('write',fullfile(final,'pointer.json'),struct('relative_work_dir','attempts/attempt_0001/work'));
    listing=dir(fullfile(work,'**','*')); files={'pointer.json'};
    for k=1:numel(listing)
        if ~listing(k).isdir, files{end+1}=strrep(erase(fullfile(listing(k).folder,listing(k).name),[final filesep]),'\','/'); end %#ok<AGROW>
    end
    k1_pm_io('seal',final,'committed.json',[pred ':FINAL:' fp],files,'AUDITED_PASS');
    k1_pm_io('write',fullfile(case_root,'checkpoint.json'),struct('status','K1_CASE_COMPLETE','failed_checks',0,'test_only',true));
end
files={'config.json','anchor/anchor_proof.json'}; rows=struct('path',{},'sha256',{});
for j=1:numel(files), rows(j)=struct('path',files{j},'sha256',k1_pm_io('sha',fullfile(fixture,files{j}))); end
k1_pm_io('write',fullfile(fixture,'manifest.json'),struct('files',rows));
guards=fullfile(out,'NO_GRAPHICS_STUB'); mkdir(guards);
fid=fopen(fullfile(guards,'render_k1_pm_summary_fields.m'),'w'); assert(fid>=0);
fprintf(fid,['function r=render_k1_pm_summary_fields(records,out)\n' ...
    'assert(numel(records)==7 && all([records.test_only]));\n' ...
    'r=struct(''failed_checks'',0,''test_only'',true,''curves_checked'',56);\nend\n']); fclose(fid);
addpath(guards,'-begin'); rehash; assert(startsWith(which('render_k1_pm_summary_fields'),guards));
rejected=false;
try, summarize_k1_recovery_cases(fixture,fullfile(out,'real_mode_reject')); catch, rejected=true; end
assert(rejected);
r=summarize_k1_recovery_cases(fixture,fullfile(out,'TEST_ONLY_summary'),'test-only');
assert(r.model_rows==14 && r.point_rows==5754 && r.new_BIE_invocations==0 && r.imported_BIE_invocations==0);
T=readtable(fullfile(out,'TEST_ONLY_summary','seven_condition_metrics.csv'));
assert(nnz(string(T.BIE_origin)=="imported_audited_N7")==2 && nnz(string(T.BIE_origin)=="qualified_R1_anchor")==2);
assert(nnz(string(T.BIE_origin)=="this_recovery_independent")==10);
receipt=struct('status','RECOVERY_SUMMARY_TEST_ONLY_METADATA_PASS','failed_checks',0, ...
    'model_rows',14,'point_rows',5754,'renderer_mocked_no_visual_acceptance',true, ...
    'new_acoustic_solves',0,'real_mode_fixture_rejected',true,'elapsed_s',toc(timer));
k1_pm_io('write',fullfile(out,'receipt.json'),receipt); disp(receipt);
end
