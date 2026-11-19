function receipt=audit_d1_recovery_origin(origin,out)
%AUDIT_D1_RECOVERY_ORIGIN Audit saved C0 in a new workspace; no BIE/PS/FK solve.
arguments
    origin (1,:) char
    out (1,:) char
end
pm_sensitivity_io('absolute',origin); pm_sensitivity_io('absolute',out);
assert(~isfolder(out)); mkdir(out); diary(fullfile(out,'diary.log')); cleanup=onCleanup(@()diary('off')); t0=tic;
previous_threads=maxNumCompThreads; maxNumCompThreads(1); thread_cleanup=onCleanup(@()maxNumCompThreads(previous_threads));
m=jsondecode(fileread(fullfile(origin,'retrieval','manifest.json')));
pm_sensitivity_io('verify',origin,m);
fp=pm_sensitivity_io('sha',fullfile(origin,'manifest.json'));
assert(strcmp(fp,'74b1a47e6bd59ece8589510b7287a3b3de00dc6bf3774142f99f2081be19d81a'));
cp=jsondecode(fileread(fullfile(origin,'output','checkpoint.json')));
assert(strcmp(cp.status,'FAILED_STOPPED') && cp.new_BIE_solves==1 && isempty(cp.completed_steps));
assert(contains(strjoin(string(cp.failed_or_unrun_checks),' '),'MATLAB:assertion:failed'));
work=fullfile(out,'audited_work'); mkdir(work); mkdir(fullfile(work,'input'));
src=fullfile(origin,'output','depth_0100m'); dest=fullfile(work,'jobs','bie_coarse'); mkdir(dest);
copyfile(fullfile(src,'common.mat'),fullfile(work,'common.mat'));
for name={'coefficients.csv','receivers.csv','frozen_config.json'}
    copyfile(fullfile(src,'input',name{1}),fullfile(work,'input',name{1}));
end
copyfile(fullfile(origin,'inputs','policy.json'),fullfile(work,'policy.json'));
for name={'data.mat','solve_invocation.json','timings.json'}
    copyfile(fullfile(src,'jobs','bie_coarse',name{1}),fullfile(dest,name{1}));
end
copyfile(fullfile(src,'jobs','bie_coarse','points.csv'),fullfile(out,'original_points.csv'));
q=load(fullfile(work,'common.mat')); b=load(fullfile(dest,'data.mat'));
rx=readtable(fullfile(work,'input','receivers.csv'));
assert(isequal(q.x,rx.x_m) && all(q.z==6));
assert(q.cfg.source.z_source_m==100 && q.cfg.source.total_width_m==.1 && q.cfg.seed==260001);
before=readtable(fullfile(out,'original_points.csv'));
write_d1_csv(table(q.x,q.z,real(b.field),imag(b.field),'VariableNames',{'x_m','z_m','real_Pa','imag_Pa'}),fullfile(dest,'points.csv'));
after=readtable(fullfile(dest,'points.csv'));
assert(isequal(after.x_m,q.x) && isequal(after.z_m,q.z),'D1 exact coordinate round-trip');
assert(max(abs(after.real_Pa+1i*after.imag_Pa-b.field))<1e-11);
% The existing independent stage auditor keeps every original assertion.
audit_r1_bie_manual_stage(dest,fullfile(work,'common.mat'),[60 10 96]);
audit=jsondecode(fileread(fullfile(dest,'audit.json'))); assert(audit.failed_checks==0);
% Include extreme finite doubles, signed small numbers, masked phase NaNs,
% string method labels and logical masks in a separate serialization check.
v=[q.x;real(b.field);imag(b.field);pi;-pi;realmin;realmax;eps;NaN];
T=table(v,repmat("FullFK",numel(v),1),isfinite(v),'VariableNames',{'value','method','phase_valid'});
write_d1_csv(T,fullfile(out,'roundtrip.csv')); back=readtable(fullfile(out,'roundtrip.csv'));
assert(isequaln(back.value,v) && all(string(back.method)=="FullFK") && isequal(logical(back.phase_valid),T.phase_valid));
pm_sensitivity_io('verify',origin,m);
receipt=struct('status','SAVED_C0_INDEPENDENT_AUDIT_PASS','failed_checks',0,'elapsed_s',toc(t0), ...
    'origin_manifest_sha256',fp,'origin_retrieval_sha_files',numel(m),'new_BIE_solves',0, ...
    'original_coordinate_changed_count',nnz(before.x_m~=q.x),'original_coordinate_max_difference_m',max(abs(before.x_m-q.x)), ...
    'fixed_coordinate_max_difference_m',max(abs(after.x_m-q.x)),'audit',audit, ...
    'audited_work_dir',work,'matlab_version',version,'matlab_pid',feature('getpid'), ...
    'next_allowed_step','Freeze recovery package importing this audited C0; user manually starts remaining nine BIE solves');
pm_sensitivity_io('write',fullfile(out,'receipt.json'),receipt);
paths={'audited_work/common.mat','audited_work/jobs/bie_coarse/data.mat','audited_work/jobs/bie_coarse/points.csv', ...
    'audited_work/jobs/bie_coarse/audit.json','audited_work/jobs/bie_coarse/audit_reconstruction.mat','roundtrip.csv','receipt.json'};
pm_sensitivity_io('write',fullfile(out,'checkpoint.json'),struct('status','AUDITED_PASS','completed_steps',{{'origin SHA','CSV fix','C0 independent audit','roundtrip regression'}}, ...
    'origin_input_versions_sha256',m,'code_versions_sha256',pm_sensitivity_io('manifest',fileparts(fileparts(fileparts(mfilename('fullpath')))), ...
    {'scripts/validation/support/write_d1_csv.m','scripts/validation/audit_d1_recovery_origin.m','scripts/validation/support/audit_r1_bie_manual_stage.m'}), ...
    'artifacts_sha256',pm_sensitivity_io('manifest',out,paths),'resource_usage',struct('new_BIE_solves',0,'elapsed_s',toc(t0)), ...
    'failed_or_unrun_checks',{{'C1-C4 and Tx1490 not run'}},'next_allowed_step',receipt.next_allowed_step));
fprintf('%s\n',jsonencode(receipt));
end
