function receipt=run_r1_bie_c4_continuation(root,mode,session_id)
%RUN_R1_BIE_C4_CONTINUATION v2 batch: exactly one possible C4 solve call.
% No C0--C3 solve branch. No implicit migration. Default never solves.
arguments
    root (1,:) char
    mode (1,:) char = 'dry-run'
    session_id (1,:) char = ''
end
assert(ismember(mode,{'dry-run','run'}));
chain=r1_bie_continuation_io('verify_migration',root);
policy=r1_bie_continuation_io('read',fullfile(root,'policy.json'));
receipt=struct('status','DRY_RUN_PASS','mode',mode,'new_acoustic_solves',0, ...
    'C4_solve_invocations_this_session',0,'legacy_completed_solves',3,'reference_stage',3, ...
    'manifest_sha256',chain.manifest_sha256,'policy_sha256',chain.policy_sha256);
dest=fullfile(root,'jobs','bie_window');
sig=[chain.predecessor ':C4:' chain.manifest_sha256 ':' chain.policy_sha256];
if isfolder(dest), state=r1_bie_continuation_io('c4_state',dest,sig); else, state='SOLVE_C4_ONCE'; end
receipt.C4_state=state;
if strcmp(mode,'dry-run'), disp(jsonencode(receipt)); return; end
assert(isunix && ~ismac && strcmp(version('-release'),'2025b'));
assert(maxNumCompThreads==1,'r1v2:Threads','Use the singleCompThread launcher');
manifest=r1_bie_continuation_io('read',fullfile(root,'manifest.json')); [~,task]=fileparts(root);
assert(strcmp(task,manifest.remote_task_name) && ~isempty(regexp(session_id,'^[A-Za-z0-9_-]+$','once')));
assert(strcmp(getenv('R1_C4_ACTIVE_SESSION'),session_id) && ...
    strcmp(strtrim(fileread(fullfile(root,'control','active_session.txt'))),session_id));
validation=r1_bie_continuation_io('read',fullfile(root,'local_validation.json'));
assert(strcmp(validation.status,'PASS') && validation.failed_checks==0 && ...
    strcmp(validation.manifest_sha256,chain.manifest_sha256) && strcmp(validation.policy_sha256,chain.policy_sha256));
runtime=r1_bie_continuation_io('read',fullfile(root,'sessions',session_id,'runtime.json'));
baseline=runtime.RAM_baseline_bytes;
assert(baseline==max(runtime.initial_available_bytes,runtime.starting_available_bytes) && runtime.reserve_fraction==.3);
local_guard(baseline); receipt.status='RUNNING'; receipt.session_id=session_id;
receipt.matlab_release=version; receipt.resource_before=local_resource(baseline); receipt.failed_checks=0;
q=load(fullfile(root,'common.mat')); ids={'bie_coarse','bie_main','bie_self','bie_space','bie_window'};
try
    if ~isfolder(dest), mkdir(dest); end
    if ~strcmp(state,'SKIP_AUDITED')
        attempt=local_attempt(dest); diary(fullfile(attempt,'diary.log')); cleanup=onCleanup(@()diary('off'));
        local_stage(root,'C4'); r1_bie_continuation_io('write',fullfile(attempt,'resource_before.json'),local_resource(baseline));
        if strcmp(state,'SOLVE_C4_ONCE')
            cfg=local_cfg(q,baseline); local_guard(baseline);
            % Persist BEFORE invoking the solver: an interrupted call cannot auto-retry.
            invocation=struct('stage',4,'invocations',1,'session_id',session_id, ...
                'manifest_sha256',chain.manifest_sha256,'policy_sha256',chain.policy_sha256);
            assert(~isfile(fullfile(dest,'solve_invocation.json')));
            r1_bie_continuation_io('write',fullfile(dest,'solve_invocation.json'),invocation);
            receipt.C4_solve_invocations_this_session=1; local_checkpoint(root,receipt,baseline,'C4_SOLVING');
            timer=tic;
            boundary=solve_helmholtz_bie_halfplane_vertical(cfg);
            timings=boundary.cost_meta.timings; timings.solver_call_s=toc(timer);
            assert(boundary.finite && all(boundary.linear_residual<1e-10) && all(boundary.offgrid_boundary_residual<1e-4));
            assert(boundary.node_count==4480 && numel(boundary.boundary_check_x_m)==476);
            field=q.flat+boundary.receiver_field; incident=cfg.incident_fn(boundary.node_x_m,boundary.node_z_m);
            numeric_config=rmfield(cfg,{'eta_fn','eta_prime_fn','incident_fn','execution_guard_fn'});
            save(fullfile(dest,'data.mat'),'boundary','field','incident','numeric_config','-v7');
            writetable(table(q.x,q.z,real(field),imag(field),'VariableNames', ...
                {'x_m','z_m','real_Pa','imag_Pa'}),fullfile(dest,'points.csv'));
            r1_bie_continuation_io('write',fullfile(dest,'config.json'),struct('stage','C4','recipe',[70 12 96], ...
                'source',q.cfg.source,'seed',260001,'receivers',411,'cache_enabled',true, ...
                'manifest_sha256',chain.manifest_sha256,'policy_sha256',chain.policy_sha256));
            r1_bie_continuation_io('write',fullfile(dest,'timings.json'),timings);
            r1_bie_continuation_io('seal',dest,'solved.json',sig, ...
                {'data.mat','points.csv','config.json','timings.json','solve_invocation.json'},'SOLVED_UNAUDITED');
            receipt.new_acoustic_solves=1;
        end
        local_guard(baseline); r1_bie_continuation_io('verify_migration',root);
        for name={'audit.json','audit_reconstruction.mat','gates_v2.json'}
            if isfile(fullfile(dest,name{1})), copyfile(fullfile(dest,name{1}),fullfile(attempt,['previous_' name{1}])); end
        end
        audit_r1_bie_manual_stage(dest,fullfile(root,'common.mat'),[70 12 96]);
        fields=cell(5,1);
        for j=1:5, b=load(fullfile(root,'jobs',ids{j},'data.mat'),'field'); fields{j}=b.field; end
        gates=r1_bie_continuation_gates(fields,q.x,q.z,policy); gates.policy_sha256=chain.policy_sha256;
        r1_bie_continuation_io('write',fullfile(dest,'gates_v2.json'),gates);
        assert(strcmp(gates.status,'PASS'),'r1v2:Convergence','One of the original four sensitivity gates failed');
        r1_bie_continuation_io('write',fullfile(attempt,'resource_after.json'),local_resource(baseline));
        diary('off'); clear cleanup;
        r1_bie_continuation_io('verify_migration',root);
        r1_bie_continuation_io('seal',dest,'committed.json',sig, ...
            {'data.mat','points.csv','config.json','timings.json','solve_invocation.json','solved.json', ...
            'audit.json','audit_reconstruction.mat','gates_v2.json'},'AUDITED_PASS');
    end
    r1_bie_continuation_io('c4_state',dest,sig);
    local_checkpoint(root,receipt,baseline,'C4_AUDITED_PASS'); local_stage(root,'FINAL');
    final=fullfile(root,'finalization'); if ~isfolder(final), mkdir(final); end
    fsig=[r1_bie_continuation_io('sha',fullfile(dest,'committed.json')) ':FINAL:' chain.manifest_sha256 ':' chain.policy_sha256];
    if isfile(fullfile(final,'committed.json'))
        seal=r1_bie_continuation_io('verify_seal',final,'committed.json',fsig,{'pointer.json'});
        assert(strcmp(seal.status,'AUDITED_PASS'));
    else
        attempt=local_attempt(final); work=fullfile(attempt,'work'); mkdir(work);
        copyfile(fullfile(root,'common.mat'),fullfile(work,'common.mat')); copyfile(fullfile(root,'input'),fullfile(work,'input'));
        for j=1:5
            d=fullfile(work,'jobs',ids{j}); mkdir(d); copyfile(fullfile(root,'jobs',ids{j},'data.mat'),fullfile(d,'data.mat'));
        end
        diary(fullfile(attempt,'diary.log')); cleanup=onCleanup(@()diary('off'));
        run_piston_pm_r1_job(work,'analyze'); diary(fullfile(attempt,'diary.log'));
        local_guard(baseline); audit_piston_pm_r1(work); diary(fullfile(attempt,'diary.log'));
        audit_r1_bie_continuation_final(work,root); local_guard(baseline);
        render_piston_pm_r1(work,'figures_final'); local_report(work,chain);
        diary('off'); clear cleanup;
        r1_bie_continuation_io('verify_migration',root); r1_bie_continuation_io('c4_state',dest,sig);
        pointer=struct('status','PASS','relative_work_dir',local_relative(final,work));
        r1_bie_continuation_io('write',fullfile(final,'pointer.json'),pointer);
        listing=dir(fullfile(work,'**','*')); paths={'pointer.json'};
        for j=1:numel(listing)
            if ~listing(j).isdir, paths{end+1}=local_relative(final,fullfile(listing(j).folder,listing(j).name)); end %#ok<AGROW>
        end
        r1_bie_continuation_io('seal',final,'committed.json',fsig,paths,'AUDITED_PASS');
    end
    receipt.status='R1_REFERENCE_COMPLETE_WAITING_REVIEW'; receipt.next_allowed_step='STOP; no R2';
    receipt.C4_completed_solves_total=1; receipt.historical_completed_solves_including_C4=4;
    local_stage(root,'DONE'); local_checkpoint(root,receipt,baseline,'COMPLETE'); disp(jsonencode(receipt));
catch ME
    receipt.status='FAILED_STOPPED'; receipt.failed_checks=1; receipt.error_identifier=ME.identifier; receipt.error=ME.message;
    receipt.next_allowed_step='Review failure; never repeat C0-C3 or an invoked C4 solve';
    local_stage(root,'STOPPED'); local_checkpoint(root,receipt,baseline,'FAIL');
    failure=fullfile(root,'sessions',session_id,'failure.json');
    r1_bie_continuation_io('write',failure,receipt); fprintf(2,'%s\n',getReport(ME,'extended','hyperlinks','off')); rethrow(ME);
end
end
function cfg=local_cfg(q,baseline)
cfg=struct('frequency_hz',4000,'c0_mps',1500,'half_width_m',70, ...
    'points_per_wavelength',12,'panel_order',8,'self_quadrature_order',96, ...
    'halfplane_h_m',-2,'window_plateau_ratio',.85,'eta_fn',@(x)closeout_s3a_pm_surface(q.profile,x), ...
    'eta_prime_fn',@(x)local_slope(q.profile,x),'incident_fn',@(x,z) ...
    baffled_strip_piston_field_vertical(x,z,q.cfg.source)-baffled_strip_piston_field_vertical(x,-z,q.cfg.source), ...
    'receiver_x_m',q.x,'receiver_z_m',q.z,'enable_static_cache',true,'execution_guard_fn',@()local_guard(baseline));
end
function slope=local_slope(profile,x), [~,slope]=closeout_s3a_pm_surface(profile,x); end
function a=local_attempt(root)
base=fullfile(root,'attempts'); if ~isfolder(base), mkdir(base); end
j=1; a=fullfile(base,sprintf('attempt_%04d',j));
while isfolder(a), j=j+1; a=fullfile(base,sprintf('attempt_%04d',j)); end
mkdir(a);
end
function r=local_resource(baseline)
t=fileread('/proc/meminfo'); m=regexp(t,'MemAvailable:\s+(\d+)\s+kB','tokens','once'); assert(~isempty(m));
r=struct('available_bytes',str2double(m{1})*1024,'RAM_baseline_bytes',baseline,'reserve_fraction',.3);
end
function local_guard(baseline)
r=local_resource(baseline); assert(r.available_bytes>=.3*baseline,'r1v2:RAM','RAM reserve fell below 30 percent');
end
function local_checkpoint(root,receipt,baseline,state)
receipt.stage_state=state; receipt.resource=local_resource(baseline);
receipt.execution_mode='User-started R2025b singleCompThread';
receipt.artifacts_sha256='jobs/*/migration_commit.json; jobs/bie_window/committed.json; finalization/committed.json';
r1_bie_continuation_io('write',fullfile(root,'checkpoint.json'),receipt);
end
function local_stage(root,label)
fid=fopen(fullfile(root,'control','stage_id.txt'),'w'); assert(fid>=0); c=onCleanup(@()fclose(fid)); fprintf(fid,'%s\n',label);
end
function rel=local_relative(root,path)
prefix=[root filesep]; assert(strncmp(path,prefix,numel(prefix))); rel=strrep(path(numel(prefix)+1:end),'\','/');
end
function local_report(work,chain)
M=readtable(fullfile(work,'model_metrics.csv'));
fid=fopen(fullfile(work,'report.md'),'w','n','UTF-8'); assert(fid>=0); c=onCleanup(@()fclose(fid));
fprintf(fid,'# R1-BIE v2 C4 completion\n\n10 cm / 4 kHz / seed260001 / 411 actual receivers. Reference=C3.\n\n');
fprintf(fid,'C0-C3 imported without new density solves. C3 was re-audited during migration. Only C4 solved here.\n');
fprintf(fid,'Original strict-policy FAIL and 8.18%% nonmonotonicity preserved in origin_evidence; v2 records it as diagnostic.\n');
fprintf(fid,'Policy SHA: %s. All four original <1e-4 gates and independent audits passed.\n\n',chain.policy_sha256);
fprintf(fid,'|Method|Complex L2|Five-times empirical resolution|Resolved|Peak x (m)|\n|---|---:|---:|---|---:|\n');
for j=[1 3]
    fprintf(fid,'|%s|%.12g|%.12g|%d|%.9g|\n',string(M.method(j)),M.complex_relative_L2(j), ...
        M.five_times_empirical_resolution(j),M.resolved(j),M.peak_x_m(j));
end
fprintf(fid,'\nEmpirical sensitivity budget, not a rigorous error bound. Unresolved does not prove equality. STOP for review; no R2.\n');
end
