function receipt=validate_k1_pm_summary_pipeline(master,out)
%VALIDATE_K1_PM_SUMMARY_PIPELINE Synthetic seals + historical arrays, no solves.
assert(~isfolder(out)); mkdir(out); timer=tic;
oldpath=path; pg=onCleanup(@()path(oldpath)); addpath(fullfile(master,'code'),'-begin');
fixture=fullfile(out,'TEST_ONLY_package'); copyfile(master,fixture);
cfg=k1_pm_io('read',fullfile(fixture,'config.json')); cfg.kind='TEST_ONLY_K1_SUMMARY';
k1_pm_io('write',fullfile(fixture,'config.json'),cfg);
ids={'bie_coarse','bie_main','bie_self','bie_space','bie_window'};
for j=1:6
    case_root=fullfile(fixture,'cases',cfg.cases(j).case_id); fp=k1_pm_io('sha',fullfile(case_root,'manifest.json'));
    model=fullfile(case_root,'model_stage'); mkdir(fullfile(model,'jobs','models'));
    copyfile(fullfile(fixture,'anchor','common.mat'),fullfile(model,'common.mat'));
    k1_pm_io('write',fullfile(model,'jobs','models','receipt.json'),struct('status','COMPLETE','failed_checks',0,'test_only',true));
    k1_pm_io('write',fullfile(model,'solve_invocation.json'),struct('invocations',1,'stage','MODELS','test_only',true));
    sig=[fp ':MODELS:' fp]; k1_pm_io('seal',model,'solved.json',sig,{'common.mat','jobs/models/receipt.json','solve_invocation.json'},'SOLVED_UNAUDITED');
    k1_pm_io('seal',model,'committed.json',sig,{'common.mat','jobs/models/receipt.json','solved.json','solve_invocation.json'},'AUDITED_PASS');
    pred=k1_pm_io('sha',fullfile(model,'committed.json'));
    for k=1:5
        dest=fullfile(case_root,'jobs',ids{k}); mkdir(dest); synthetic_payload='TEST_ONLY_NOT_A_DENSITY'; save(fullfile(dest,'data.mat'),'synthetic_payload','-v7');
        k1_pm_io('write',fullfile(dest,'audit.json'),struct('status','PASS','failed_checks',0,'test_only',true));
        k1_pm_io('write',fullfile(dest,'gates.json'),struct('status','PASS','test_only',true));
        k1_pm_io('write',fullfile(dest,'solve_invocation.json'),struct('invocations',1,'stage',sprintf('C%d',k-1),'test_only',true));
        sig=[pred ':C' num2str(k-1) ':' fp];
        k1_pm_io('seal',dest,'committed.json',sig,{'data.mat','audit.json','gates.json','solve_invocation.json'},'AUDITED_PASS');
        pred=k1_pm_io('sha',fullfile(dest,'committed.json'));
    end
    final=fullfile(case_root,'finalization'); work=fullfile(final,'attempts','attempt_0001','work'); mkdir(work);
    for file={'common.mat','comparison.mat','model_metrics.csv','receiver_fields_and_errors.csv','independent_final_metrics_audit.json'}
        copyfile(fullfile(fixture,'anchor',file{1}),fullfile(work,file{1}));
    end
    mkdir(fullfile(work,'audit')); copyfile(fullfile(fixture,'anchor','audit','receipt.json'),fullfile(work,'audit','receipt.json'));
    k1_pm_io('write',fullfile(final,'pointer.json'),struct('relative_work_dir','attempts/attempt_0001/work','status','PASS','test_only',true));
    listing=dir(fullfile(work,'**','*')); paths={'pointer.json'};
    for k=1:numel(listing)
        if ~listing(k).isdir, paths{end+1}=strrep(erase(fullfile(listing(k).folder,listing(k).name),[final filesep]),'\','/'); end %#ok<AGROW>
    end
    k1_pm_io('seal',final,'committed.json',[pred ':FINAL:' fp],paths,'AUDITED_PASS');
    k1_pm_io('write',fullfile(case_root,'checkpoint.json'),struct('status','K1_CASE_COMPLETE','failed_checks',0,'test_only',true));
end
% Bind only immutable snapshot inputs, not fake outputs, to the fixture master.
m=k1_pm_io('read',fullfile(fixture,'manifest.json'));
ix=find(strcmp({m.files.path},'config.json')); assert(isscalar(ix)); m.files(ix).sha256=k1_pm_io('sha',fullfile(fixture,'config.json'));
k1_pm_io('write',fullfile(fixture,'manifest.json'),m);
rejected=false;
try
    summarize_k1_pm_cases(fixture,fullfile(out,'real_mode_must_reject'));
catch
    rejected=true;
end
assert(rejected);
r=summarize_k1_pm_cases(fixture,fullfile(out,'TEST_ONLY_summary'),'test-only');
assert(strcmp(r.status,'K1_SUMMARY_TEST_ONLY_PASS') && r.model_rows==14 && r.point_rows==5754 && r.new_physical_cases==0 && r.new_BIE_invocations==0);
receipt=struct('status','K1_SUMMARY_PIPELINE_TEST_ONLY_PASS','failed_checks',0,'model_rows',r.model_rows, ...
    'point_rows',r.point_rows,'curves',r.figure_curves,'real_mode_fixture_rejected',true, ...
    'new_acoustic_solves',0,'new_K1_physical_results',0,'elapsed_s',toc(timer));
k1_pm_io('write',fullfile(out,'receipt.json'),receipt); disp(jsonencode(receipt));
end
