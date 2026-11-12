function receipt=run_munk_ffp_full(stage,run_id)
%RUN_MUNK_FFP_FULL Frozen six-job full-field/reference matrix, serial only.
% init is preparation only. next requires a user-budget authorization receipt.
arguments
    stage (1,:) char = 'status'
    run_id (1,:) char = '20261112_full_01'
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
assert(isempty(regexp(run_id,'[^A-Za-z0-9_-]','once')));
addpath(fullfile(root,'scripts','validation','support'),'-begin');
d=fullfile(root,'results','validation','munk_ffp_5000m',run_id);
if strcmp(stage,'init')
    assert(~isfolder(d),'Never overwrite run directory'); mkdir(d);
    plan=struct([]); h=[2.5 1.25 1.25 1.25 1.25 .625];
    K=[80 80 100 100 100 100]; panel=[.25 .25 .25 .125 .125 .125]; a=[.2 .2 .2 .2 .1 .1];
    for j=1:6
        c=make_munk_5000m_ffp_config(h(j)); c.kappa_max_radpm=K(j); c.panel_width_radpm=panel(j); c.contour_offset_radpm=a(j);
        c.depth_m=unique([(0:5:5000)';(4900:.01:5000)']); c.max_wall_s=90;
        plan=[plan struct('job_id',sprintf('J%02d',j),'config',c)]; %#ok<AGROW>
    end
    save(fullfile(d,'frozen_plan.mat'),'plan','-v7');
    receipt=struct('stage','init','status','PREPARED_AWAITING_USER_RESOURCE_BUDGET', ...
        'run_id',run_id,'proposed_total_budget_s',900,'per_job_budget_s',90, ...
        'required_start_available_bytes',2^30,'minimum_runtime_available_bytes',256*2^20, ...
        'job_count',6,'full_grid_size',[numel(plan(1).config.depth_m) numel(plan(1).config.range_m)], ...
        'next_allowed_step','Obtain explicit user budget authorization; verify RAM before next.');
    local_json(fullfile(d,'plan_receipt.json'),receipt); fprintf('%s\n',jsonencode(receipt)); return;
end
assert(isfolder(d),'Initialize a new isolated run'); s=load(fullfile(d,'frozen_plan.mat'),'plan'); plan=s.plan;
if strcmp(stage,'status')
    complete=false(1,6); for j=1:6, complete(j)=isfile(fullfile(d,[plan(j).job_id '.mat'])); end
    receipt=struct('run_id',run_id,'completed_jobs',complete,'budget_authorized',isfile(fullfile(d,'authorization.json')));
    fprintf('%s\n',jsonencode(receipt)); return;
end
assert(strcmp(stage,'next'),'Supported stages init/status/next');
authfile=fullfile(d,'authorization.json');
assert(isfile(authfile),'munk_ffp:BudgetPending','Full-field matrix awaits explicit user resource budget.');
auth=jsondecode(fileread(authfile)); assert(auth.approved_total_budget_s>=900 && strcmp(auth.approved_by,'user'));
used=0; recs=dir(fullfile(d,'J*_receipt.json'));
for j=1:numel(recs), old=jsondecode(fileread(fullfile(d,recs(j).name))); used=used+old.elapsed_s; end
assert(used+90<=auth.approved_total_budget_s,'munk_ffp:TotalBudget','No remaining per-job budget.');
[~,mem]=memory; assert(mem.PhysicalMemory.Available>=2^30,'munk_ffp:StartRAM','At least 1GiB available RAM required to begin full-field job.');
job=0;
for j=1:6
    path=fullfile(d,[plan(j).job_id '.mat']); rec=fullfile(d,[plan(j).job_id '_receipt.json']);
    if isfile(path), continue; end
    assert(~isfile(rec),'Previous attempt failed; inspect its checkpoint before any new retry'); job=j; break;
end
assert(job>0,'All six planned jobs already complete; audit saved results next.');
tag=plan(job).job_id; receiptfile=fullfile(d,[tag '_receipt.json']);
diary(fullfile(d,[tag '_diary.log'])); t=tic;
receipt=struct('job_id',tag,'stage','full_field','status','RUNNING','exit_status',NaN,'elapsed_s',0, ...
    'matlab_release',version,'pid',feature('getpid'),'execution_mode','existing local MATLAB MCP', ...
    'starting_available_bytes',mem.PhysicalMemory.Available,'failed_or_unrun_checks',{{}}, ...
    'qualification','UNVERIFIED until complete saved-data audit');
try
    output=solve_layered_fluid_ffp(plan(job).config);
    save(fullfile(d,[tag '.mat']),'output','-v7');
    receipt.status='NUMERICAL_OUTPUT_SAVED_NOT_YET_QUALIFIED'; receipt.exit_status=0;
    receipt.metrics=struct('spectral_nodes',output.meta.spectral_node_count,'solve_wall_s',output.elapsed_s, ...
        'finite_receiver_count',nnz(isfinite(output.green)),'source_singular_count',nnz(output.source_singularity_mask), ...
        'grid_size',size(output.green));
catch err
    receipt.status='FAILED'; receipt.exit_status=1; receipt.failed_or_unrun_checks={err.identifier,err.message};
    fid=fopen(fullfile(d,[tag '_error.log']),'w'); fprintf(fid,'%s',getReport(err,'extended','hyperlinks','off')); fclose(fid);
end
receipt.elapsed_s=toc(t); receipt.prior_full_job_s=used;
local_json(receiptfile,receipt); diary off; fprintf('%s\n',jsonencode(receipt));
if receipt.exit_status~=0, error('munk_ffp:JobFailed','Read saved checkpoint; no automatic retry.'); end
end
function local_json(path,value)
fid=fopen(path,'w'); assert(fid>=0); cleaner=onCleanup(@()fclose(fid)); fprintf(fid,'%s',jsonencode(value,'PrettyPrint',true));
end
