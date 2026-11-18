function receipt=run_r2_pm_seed(root,mode,session_id)
%RUN_R2_PM_SEED One frozen seed, models then five independent BIE solves.
% Default is geometry/SHA only. Real execution requires the Linux lock owner.
arguments
    root (1,:) char
    mode (1,:) char = 'dry-run'
    session_id (1,:) char = ''
end
assert(ismember(mode,{'dry-run','run'}));
m=r2_pm_io('verify_manifest',root);
fp=r2_pm_io('sha',fullfile(root,'manifest.json'));
policy=r2_pm_io('read',fullfile(root,'policy.json'));
recipes=[60 10 96;60 12 96;60 12 64;60 14 96;70 12 96];
ids={'bie_coarse','bie_main','bie_self','bie_space','bie_window'};
receipt=struct('status','DRY_RUN_INPUT_CHECKS_PASS','seed',m.seed,'frequency_hz',4000, ...
    'strip_width_m',.1,'receivers',411,'manifest_sha256',fp,'recipes',recipes, ...
    'new_acoustic_solves',0,'failed_checks',0,'mode',mode);
if strcmp(mode,'dry-run'), disp(jsonencode(receipt)); return; end
assert(strcmp(version('-release'),'2025b'));
[~,name]=fileparts(root); assert(strcmp(name,m.remote_task_name));
assert(~isempty(regexp(session_id,'^[A-Za-z0-9_-]+$','once')));
assert(strcmp(getenv('R2_ACTIVE_SESSION'),session_id),'r2:Launcher','Use the Linux launcher');
assert(strcmp(strtrim(fileread(fullfile(root,'control','active_session.txt'))),session_id));
runtime=r2_pm_io('read',fullfile(root,'sessions',session_id,'shared_runtime.json'));
baseline=runtime.RAM_baseline_bytes; local_guard(baseline);
receipt.status='RUNNING'; receipt.session_id=session_id; receipt.completed_steps={};
receipt.matlab_version=version; receipt.matlab_pid=feature('getpid'); receipt.RAM_baseline_bytes=baseline;
receipt.solve_counter_scope='this MATLAB session; durable per-stage invocation markers are authoritative';
receipt.counter_definition='new_acoustic_solves counts BIE solver invocations; MODELS is separately marked (six PS and four FullFK evaluations maximum)';
receipt.input_code_versions='manifest.json'; receipt.policy_sha256=m.policy_sha256;
receipt.stdout_log=['sessions/' session_id '/stdout.log']; receipt.stderr_log=['sessions/' session_id '/stderr.log'];
pred=fp;
try
    % Models have their own sealed result; a failed/unsealed invocation cannot
    % be repeated by Continue. No R1 common.mat or scattering density imported.
    dest=fullfile(root,'model_stage'); if ~isfolder(dest), mkdir(dest); end
    sig=[pred ':MODELS:' fp]; local_status(root,receipt,'MODELS','CHECKING');
    if isfile(fullfile(dest,'committed.json'))
        seal=r2_pm_io('verify_seal',dest,'committed.json',sig); assert(strcmp(seal.status,'AUDITED_PASS'));
    elseif isfile(fullfile(dest,'solved.json'))
        r2_pm_io('verify_seal',dest,'solved.json',sig);
        local_models_check(dest,m.seed);
        r2_pm_io('seal',dest,'committed.json',sig,{'common.mat','jobs/models/receipt.json','solved.json','solve_invocation.json'},'AUDITED_PASS');
    else
        local_mark_invocation(dest,m.seed,'MODELS',fp);
        mkdir(fullfile(dest,'input')); copyfile(fullfile(root,'input','*'),fullfile(dest,'input'));
        local_guard(baseline); run_r2_piston_pm_job(dest,'models');
        r2_pm_io('seal',dest,'solved.json',sig,{'common.mat','jobs/models/receipt.json','solve_invocation.json'},'SOLVED_UNAUDITED');
        local_models_check(dest,m.seed);
        r2_pm_io('seal',dest,'committed.json',sig,{'common.mat','jobs/models/receipt.json','solved.json','solve_invocation.json'},'AUDITED_PASS');
    end
    local_verify(root,fp); local_models_check(dest,m.seed);
    q=load(fullfile(dest,'common.mat')); pred=r2_pm_io('sha',fullfile(dest,'committed.json'));
    receipt.completed_steps{end+1}='MODELS';
    for stage=0:4
        local_verify(root,fp); local_guard(baseline);
        label=sprintf('C%d',stage); dest=fullfile(root,'jobs',ids{stage+1});
        if ~isfolder(dest), mkdir(dest); end
        sig=[pred ':' label ':' fp]; state=r2_pm_io('resume_state',dest,sig);
        local_status(root,receipt,label,state);
        if ~strcmp(state,'SKIP_AUDITED')
            attempt=local_attempt(dest); diary(fullfile(attempt,'diary.log'));
            cleanup=onCleanup(@()diary('off'));
            r2_pm_io('write',fullfile(attempt,'resource_before.json'),local_resource(baseline));
            if strcmp(state,'SOLVE')
                local_mark_invocation(dest,m.seed,label,fp);
                cfg=local_cfg(q,recipes(stage+1,:)); cfg.enable_static_cache=true;
                cfg.execution_guard_fn=@()local_guard(baseline); timer=tic;
                receipt.new_acoustic_solves=receipt.new_acoustic_solves+1;
                boundary=solve_helmholtz_bie_halfplane_vertical(cfg);
                timings=boundary.cost_meta.timings; timings.solver_call_s=toc(timer);
                field=q.flat+boundary.receiver_field;
                incident=cfg.incident_fn(boundary.node_x_m,boundary.node_z_m);
                numeric_config=rmfield(cfg,{'eta_fn','eta_prime_fn','incident_fn','execution_guard_fn'});
                save(fullfile(dest,'data.mat'),'boundary','field','incident','numeric_config','-v7');
                writetable(table(q.x,q.z,real(field),imag(field),'VariableNames', ...
                    {'x_m','z_m','real_Pa','imag_Pa'}),fullfile(dest,'points.csv'));
                r2_pm_io('write',fullfile(dest,'config.json'),struct('seed',m.seed,'source',q.cfg.source, ...
                    'recipe',recipes(stage+1,:),'phase_contract_version','2.0.0','manifest_sha256',fp));
                r2_pm_io('write',fullfile(dest,'timings.json'),timings);
                r2_pm_io('seal',dest,'solved.json',sig, ...
                    {'data.mat','points.csv','config.json','timings.json','solve_invocation.json'},'SOLVED_UNAUDITED');
            end
            % Re-audits go into a fresh attempt; failed audit products survive.
            for item={'data.mat','points.csv'}, copyfile(fullfile(dest,item{1}),fullfile(attempt,item{1})); end
            audit_r1_bie_manual_stage(attempt,fullfile(root,'model_stage','common.mat'),recipes(stage+1,:));
            for item={'audit.json','audit_reconstruction.mat'}, copyfile(fullfile(attempt,item{1}),fullfile(dest,item{1})); end
            gates=local_gates(root,ids,stage,q,policy);
            r2_pm_io('write',fullfile(attempt,'gates.json'),gates);
            r2_pm_io('write',fullfile(dest,'gates.json'),gates);
            assert(strcmp(gates.status,'PASS'),'r2:Convergence','Original sensitivity gate failed');
            r2_pm_io('write',fullfile(attempt,'resource_after.json'),local_resource(baseline));
            diary('off'); clear cleanup;
            local_verify(root,fp);
            r2_pm_io('seal',dest,'committed.json',sig, ...
                {'data.mat','points.csv','config.json','timings.json','solved.json','solve_invocation.json', ...
                'audit.json','audit_reconstruction.mat','gates.json'},'AUDITED_PASS');
        end
        pred=r2_pm_io('sha',fullfile(dest,'committed.json'));
        receipt.completed_steps{end+1}=label;
        local_status(root,receipt,label,'AUDITED_PASS');
    end
    final=fullfile(root,'finalization'); if ~isfolder(final), mkdir(final); end
    sig=[pred ':FINAL:' fp]; local_status(root,receipt,'FINAL','CHECKING');
    if isfile(fullfile(final,'committed.json'))
        seal=r2_pm_io('verify_seal',final,'committed.json',sig); assert(strcmp(seal.status,'AUDITED_PASS'));
        pointer=r2_pm_io('read',fullfile(final,'pointer.json'));
    else
        attempt=local_attempt(final); work=fullfile(attempt,'work'); mkdir(work);
        copyfile(fullfile(root,'model_stage','common.mat'),fullfile(work,'common.mat'));
        copyfile(fullfile(root,'policy.json'),fullfile(work,'policy.json'));
        for j=1:5
            target=fullfile(work,'jobs',ids{j}); mkdir(target);
            copyfile(fullfile(root,'jobs',ids{j},'data.mat'),fullfile(target,'data.mat'));
        end
        diary(fullfile(attempt,'diary.log')); cleanup=onCleanup(@()diary('off'));
        run_r2_piston_pm_job(work,'analyze'); diary(fullfile(attempt,'diary.log'));
        audit_piston_pm_r1(work); diary(fullfile(attempt,'diary.log'));
        audit_r2_pm_final(work); render_r2_pm_seed(work,'figures_final');
        summarize_r2_pm_fields({work},fullfile(work,'spatial_summary'));
        diary('off'); clear cleanup; local_guard(baseline); local_verify(root,fp);
        pointer=struct('status','PASS','relative_work_dir',local_relative(final,work),'seed',m.seed);
        r2_pm_io('write',fullfile(final,'pointer.json'),pointer);
        listing=dir(fullfile(work,'**','*')); paths={'pointer.json'};
        for j=1:numel(listing)
            if ~listing(j).isdir, paths{end+1}=local_relative(final,fullfile(listing(j).folder,listing(j).name)); end %#ok<AGROW>
        end
        r2_pm_io('seal',final,'committed.json',sig,paths,'AUDITED_PASS');
    end
    receipt.status='R2_SEED_COMPLETE_WAITING_REVIEW'; receipt.final_work_dir=['finalization/' pointer.relative_work_dir];
    receipt.next_allowed_step='STOP; user reviews this seed before manually starting the other seed or saved-data summary';
    local_status(root,receipt,'DONE','COMPLETE'); disp(jsonencode(receipt));
catch ME
    diary('off'); receipt.status='FAILED_STOPPED'; receipt.failed_checks=1;
    receipt.error=ME.message; receipt.error_identifier=ME.identifier;
    if exist('attempt','var') && isfolder(attempt)
        r2_pm_io('write',fullfile(attempt,'failure.json'),struct('status','FAIL','error',ME.message, ...
            'identifier',ME.identifier,'time',local_time()));
    end
    receipt.next_allowed_step='Inspect logs and SHA; manual audit-only recovery for sealed solves; never automatically retry a solve';
    local_status(root,receipt,'STOPPED','FAIL');
    fprintf(2,'%s\n',getReport(ME,'extended','hyperlinks','off')); rethrow(ME);
end
end
function local_models_check(dest,seed)
q=load(fullfile(dest,'common.mat')); r=r2_pm_io('read',fullfile(dest,'jobs','models','receipt.json'));
assert(strcmp(r.status,'COMPLETE') && r.failed_checks==0 && q.cfg.seed==seed);
parent=fileparts(dest); frozen=r2_pm_io('read',fullfile(parent,'input','frozen_config.json'));
assert(isequal(orderfields(q.cfg.source),orderfields(frozen.source)) && ...
    isequal(q.profile,readtable(fullfile(parent,'input','coefficients.csv'))) && ...
    isequal(q.x,(-205:205)'*.1953125) && isequal(q.z,3*ones(411,1)));
assert(all(q.ps_values(1:2)<1e-4) && q.ps_values(3)<1e-10 && all(q.fk_values<1e-5));
assert(isequal(q.cfg.ps_recipe,[128 .05 .1875;256 .05 .1875;512 .05 .1875;1024 .05 .1875;1024 .025 .1875;1024 .025 .09375]));
assert(isequal(q.cfg.fk_recipe,[400 50;800 50;960 60;1440 90]));
end
function local_verify(root,fp)
assert(strcmp(r2_pm_io('sha',fullfile(root,'manifest.json')),fp),'r2:ManifestDrift','Manifest changed during run');
r2_pm_io('verify_manifest',root);
end
function cfg=local_cfg(q,rc)
cfg=struct('frequency_hz',4000,'c0_mps',1500,'half_width_m',rc(1), ...
    'points_per_wavelength',rc(2),'panel_order',8,'self_quadrature_order',rc(3), ...
    'halfplane_h_m',-2,'window_plateau_ratio',.85,'eta_fn',@(x)closeout_s3a_pm_surface(q.profile,x), ...
    'eta_prime_fn',@(x)local_slope(q.profile,x), ...
    'incident_fn',@(x,z)baffled_strip_piston_field_vertical(x,z,q.cfg.source)- ...
    baffled_strip_piston_field_vertical(x,-z,q.cfg.source),'receiver_x_m',q.x,'receiver_z_m',q.z);
end
function y=local_slope(p,x)
[~,y]=closeout_s3a_pm_surface(p,x);
end
function g=local_gates(root,ids,stage,q,policy)
if stage>=3
    fields=cell(stage+1,1);
    for j=1:stage+1, b=load(fullfile(root,'jobs',ids{j},'data.mat'),'field'); fields{j}=b.field; end
    g=r1_bie_continuation_gates(fields,q.x,q.z,policy);
else
    pairs=[1 2;3 2]; v=[];
    for j=1:stage
        A=load(fullfile(root,'jobs',ids{pairs(j,1)},'data.mat'),'field');
        B=load(fullfile(root,'jobs',ids{pairs(j,2)},'data.mat'),'field');
        value=norm(A.field-B.field)/norm(B.field); edges=round(linspace(1,411,5));
        for k=1:4, ii=edges(k):edges(k+1); value=max(value,norm(A.field(ii)-B.field(ii))/norm(B.field(ii))); end
        v(end+1)=value; %#ok<AGROW>
    end
    status='PASS'; if ~all(isfinite(v) & v<1e-4), status='FAIL'; end
    g=struct('status',status,'sensitivity_gate',1e-4,'pair_values',v);
end
end
function local_mark_invocation(dest,seed,label,fp)
p=fullfile(dest,'solve_invocation.json');
assert(~isfile(p),'r2:NoRetry','Unsealed prior invocation exists; stop for review: %s',p);
r2_pm_io('write',p,struct('seed',seed,'stage',label,'invocations',1,'manifest_sha256',fp,'time',local_time()));
end
function p=local_attempt(root)
j=1; p=fullfile(root,'attempts',sprintf('attempt_%04d',j));
while isfolder(p), j=j+1; p=fullfile(root,'attempts',sprintf('attempt_%04d',j)); end
mkdir(p);
end
function r=local_resource(baseline)
t=fileread('/proc/meminfo'); m=regexp(t,'MemAvailable:\s+(\d+)\s+kB','tokens','once'); assert(~isempty(m));
group=getenv('R2_RAM_GROUP'); assert(strcmp(group,'/home/misc/cwx/r2_pm_20261117_03_control'));
shared=str2double(strtrim(fileread(fullfile(group,'RAM_baseline_bytes.txt')))); assert(isfinite(shared) && shared>0);
baseline=max(baseline,shared);
r=struct('time',local_time(),'available_bytes',str2double(m{1})*1024,'baseline_bytes',baseline,'reserve_fraction',.3);
end
function local_guard(baseline)
r=local_resource(baseline); assert(r.available_bytes>=.3*baseline,'r2:RAM','30 percent RAM reserve breached');
end
function local_status(root,r,stage,state)
r.stage=stage; r.stage_state=state; r.updated_at=local_time(); r.resource_usage=local_resource(r.RAM_baseline_bytes);
r.session_elapsed_s=posixtime(datetime('now','TimeZone','UTC'))- ...
    str2double(strtrim(fileread(fullfile(root,'sessions',r.session_id,'start_epoch.txt'))));
r.audit_and_sha_paths='jobs/*/{audit.json,gates.json,solved.json,committed.json}; finalization/committed.json';
fid=fopen(fullfile(root,'control','stage_id.txt'),'w'); assert(fid>=0); fprintf(fid,'%s\n',stage); fclose(fid);
r2_pm_io('write',fullfile(root,'checkpoint.json'),r);
end
function s=local_time()
s=char(datetime('now','TimeZone','Asia/Shanghai','Format','yyyy-MM-dd HH:mm:ss Z'));
end
function s=local_relative(root,path)
assert(startsWith(path,[root filesep])); s=strrep(path(numel(root)+2:end),'\','/');
end
