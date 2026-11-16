function receipt=run_r1_bie_manual(root,mode,session_id)
%RUN_R1_BIE_MANUAL Manual, resumable C1--C4 batch. Default is solver-free.
% Use the packaged Linux launcher for 'run'; it supplies the exclusive lock,
% one singleCompThread process, stdout/stderr and task-owned RAM watchdog.
arguments
    root (1,:) char
    mode (1,:) char = 'dry-run'
    session_id (1,:) char = ''
end
assert(ismember(mode,{'dry-run','run'}));
manifest=r1_bie_manual_io('verify_manifest',root);
fingerprint=r1_bie_manual_io('sha',fullfile(root,'manifest.json'));
q=load(fullfile(root,'common.mat')); s=q.cfg.source;
assert(s.total_width_m==.1 && s.frequency_hz==4000 && s.z_source_m==100 && ...
    s.c_mps==1500 && s.rho_kgpm3==1025 && s.v0_mps==.001);
assert(numel(q.x)==411 && all(q.z==3) && numel(unique(q.x))==411);
assert(height(q.profile)==12 && abs(max(q.profile.k_rad_per_m)-.471238898038469)<1e-14);
assert(strcmp(q.cfg.phase_contract_version,'2.0.0'));
p=r1_bie_manual_io('read',fullfile(root,'parent','checkpoint.json'));
assert(p.failed_checks==0 && p.C0_equivalence.density_relative==0 && ...
    p.C0_equivalence.reflected_relative==0 && p.C0_equivalence.node_geometry_exact);
old_audit=r1_bie_manual_io('read',fullfile(root,'parent','C0_audit_receipt.json'));
assert(strcmp(old_audit.status,'PASS') && old_audit.failed_checks==0);
C0=load(fullfile(root,'parent','C0_data.mat'));
assert(C0.boundary.node_count==3200 && C0.boundary.half_width_m==60 && ...
    C0.boundary.points_per_wavelength_requested==10 && C0.boundary.self_quadrature_order==96);
assert(isequal(C0.boundary.receiver_x_m,q.x) && isequal(C0.boundary.receiver_z_m,q.z));
assert(isequal(C0.field,q.flat+C0.boundary.receiver_field) && C0.boundary.finite);
recipes=[60 10 96;60 12 96;60 12 64;60 14 96;70 12 96];
ids={'bie_coarse','bie_main','bie_self','bie_space','bie_window'};
receipt=struct('status','DRY_RUN_INPUT_CHECKS_PASS','mode',mode,'schema_version',1, ...
    'manifest_sha256',fingerprint,'stage_ids',{{'C0_import','C1','C2','C3','C4'}}, ...
    'recipes',recipes,'C0_solve_forbidden',true,'cumulative_time_cap',false,'cost_prediction_gate',false, ...
    'new_acoustic_solves',0,'solve_counter_scope','this MATLAB session only','execution_policy',manifest.execution_policy);
if strcmp(mode,'dry-run')
    receipt.resume_selftest=local_resume_selftest(root);
    fprintf('%s\n',jsonencode(receipt)); return;
end
assert(strcmp(version('-release'),'2025b'));
[~,task_name]=fileparts(root); assert(strcmp(task_name,manifest.remote_task_name),'Remote task directory does not match package');
assert(~isempty(regexp(session_id,'^[A-Za-z0-9_-]+$','once')));
assert(strcmp(strtrim(fileread(fullfile(root,'control','active_session.txt'))),session_id));
assert(strcmp(getenv('R1_MANUAL_ACTIVE_SESSION'),session_id),'Use the exclusive-lock Linux launcher');
session=fullfile(root,'sessions',session_id); assert(isfolder(session));
runtime=r1_bie_manual_io('read',fullfile(session,'runtime.json'));
baseline=runtime.RAM_baseline_bytes; local_guard(baseline);
receipt.status='RUNNING'; receipt.session_id=session_id; receipt.matlab_release=version;
receipt.execution_mode='User-started SSH MATLAB R2025b -singleCompThread'; receipt.failed_checks=0;
receipt.config=struct('source',s,'seed',260001,'receivers',411,'receiver_z_m',3, ...
    'recipes',recipes,'phase_contract_version','2.0.0','RAM_reserve_fraction',.3);
receipt.next_allowed_step='Next fixed stage only after predecessor SHA and independent audit PASS';
receipt.matlab_pid=feature('getpid'); receipt.RAM_baseline_bytes=baseline;
receipt.stdout_log=['sessions/' session_id '/stdout.log']; receipt.stderr_log=['sessions/' session_id '/stderr.log'];
receipt.completed_steps={}; receipt.started_at=local_time(); receipt.new_acoustic_solves=0;
predecessor=fingerprint;
try
    for stage=0:4
        local_verify_package(root,fingerprint);
        label=sprintf('C%d',stage); dest=fullfile(root,'jobs',ids{stage+1});
        if ~isfolder(dest), mkdir(dest); end
        signature=[predecessor ':' label ':' fingerprint];
        state=r1_bie_manual_io('resume_state',dest,signature);
        local_status(root,receipt,label,state);
        if strcmp(state,'SKIP_AUDITED')
            fprintf('%s SKIP_AUDITED (no solve, no audit repeated)\n',label);
        else
            attempt=local_attempt(dest); diary(fullfile(attempt,'diary.log'));
            attempt_guard=onCleanup(@()diary('off'));
            try
                r1_bie_manual_io('write',fullfile(attempt,'resource_before.json'),local_resource(baseline));
                if strcmp(state,'SOLVE')
                    if stage==0
                        % The only C0 path is an import. No C0 solver call.
                        copyfile(fullfile(root,'parent','C0_data.mat'),fullfile(dest,'data.mat'));
                        copyfile(fullfile(root,'parent','C0_points.csv'),fullfile(dest,'points.csv'));
                        numeric_config=C0.numeric_config; boundary=C0.boundary; %#ok<NASGU>
                        timings=struct('imported',true,'source_C0_function_s',p.C0_function_s);
                    else
                        cfg=local_cfg(q,recipes(stage+1,:)); cfg.enable_static_cache=true;
                        cfg.execution_guard_fn=@()local_guard(baseline);
                        local_guard(baseline); timer=tic;
                        boundary=solve_helmholtz_bie_halfplane_vertical(cfg);
                        timings=boundary.cost_meta.timings; timings.solver_call_s=toc(timer);
                        field=q.flat+boundary.receiver_field;
                        incident=cfg.incident_fn(boundary.node_x_m,boundary.node_z_m);
                        numeric_config=rmfield(cfg,{'eta_fn','eta_prime_fn','incident_fn','execution_guard_fn'});
                        temp=fullfile(dest,'data_pending.mat');
                        if isfile(temp), movefile(temp,fullfile(attempt,'previous_unsealed_pending.mat')); end
                        save(temp,'boundary','field','incident','numeric_config','-v7');
                        movefile(temp,fullfile(dest,'data.mat'));
                        writetable(table(q.x,q.z,real(field),imag(field),'VariableNames', ...
                            {'x_m','z_m','real_Pa','imag_Pa'}),fullfile(dest,'points.csv'));
                        receipt.new_acoustic_solves=receipt.new_acoustic_solves+1;
                    end
                    r1_bie_manual_io('write',fullfile(dest,'config.json'),local_config_record(q,recipes(stage+1,:),label,fingerprint));
                    r1_bie_manual_io('write',fullfile(dest,'timings.json'),timings);
                    r1_bie_manual_io('seal',dest,'solved.json',signature, ...
                        {'data.mat','points.csv','config.json','timings.json'},'SOLVED_UNAUDITED');
                end
                local_guard(baseline);
                local_verify_package(root,fingerprint);
                % Preserve failed audit evidence before an audit-only restart.
                for name={'audit.json','gates.json','audit_reconstruction.mat'}
                    if isfile(fullfile(dest,name{1})), copyfile(fullfile(dest,name{1}),fullfile(attempt,['previous_' name{1}])); end
                end
                audit_r1_bie_manual_stage(dest,fullfile(root,'common.mat'),recipes(stage+1,:));
                gates=local_gates(root,ids,stage);
                r1_bie_manual_io('write',fullfile(dest,'gates.json'),gates);
                r1_bie_manual_io('write',fullfile(attempt,'resource_after.json'),local_resource(baseline));
                diary('off'); clear attempt_guard;
                local_verify_package(root,fingerprint);
                r1_bie_manual_io('seal',dest,'committed.json',signature, ...
                    {'data.mat','points.csv','config.json','timings.json','audit.json','audit_reconstruction.mat','gates.json'},'AUDITED_PASS');
            catch ME
                r1_bie_manual_io('write',fullfile(attempt,'failure.json'), ...
                    struct('status','FAIL','error',ME.message,'identifier',ME.identifier,'time',local_time()));
                fprintf(2,'%s\n',getReport(ME,'extended','hyperlinks','off')); diary('off'); clear attempt_guard;
                rethrow(ME);
            end
        end
        predecessor=r1_bie_manual_io('sha',fullfile(dest,'committed.json'));
        receipt.completed_steps{end+1}=[label ' AUDITED_PASS'];
        local_status(root,receipt,label,'AUDITED_PASS');
    end
    local_guard(baseline); local_verify_package(root,fingerprint); local_stage_file(root,'FINAL');
    final=fullfile(root,'finalization'); if ~isfolder(final), mkdir(final); end
    signature=[predecessor ':FINAL:' fingerprint];
    if isfile(fullfile(final,'committed.json'))
        r1_bie_manual_io('verify_seal',final,'committed.json',signature);
        pointer=r1_bie_manual_io('read',fullfile(final,'pointer.json'));
    else
        attempt=local_attempt(final); work=fullfile(attempt,'work'); mkdir(work);
        local_prepare_final(root,work,ids); diary(fullfile(attempt,'diary.log'));
        try
            run_piston_pm_r1_job(work,'analyze');
            diary(fullfile(attempt,'diary.log'));
            local_guard(baseline); audit_piston_pm_r1(work);
            diary(fullfile(attempt,'diary.log'));
            audit_r1_bie_manual_final(work);
            local_guard(baseline); render_piston_pm_r1(work,'figures_final');
            local_final_report(work);
            diary('off');
        catch ME
            fprintf(2,'%s\n',getReport(ME,'extended','hyperlinks','off')); diary('off');
            r1_bie_manual_io('write',fullfile(attempt,'failure.json'),struct('error',ME.message,'identifier',ME.identifier)); rethrow(ME);
        end
        pointer=struct('status','PASS','relative_work_dir',strrep(local_relative_path(final,work),'\','/'));
        r1_bie_manual_io('write',fullfile(final,'pointer.json'),pointer);
        local_verify_package(root,fingerprint);
        listing=dir(fullfile(work,'**','*')); paths={'pointer.json'};
        for j=1:numel(listing)
            if ~listing(j).isdir
                path=fullfile(listing(j).folder,listing(j).name);
                paths{end+1}=local_relative_path(final,path); %#ok<AGROW>
            end
        end
        r1_bie_manual_io('seal',final,'committed.json',signature,paths,'AUDITED_PASS');
    end
    receipt.status='R1_REFERENCE_COMPLETE_WAITING_REVIEW';
    receipt.final_report=['finalization/' pointer.relative_work_dir '/report.md'];
    receipt.completed_steps{end+1}='Independent final audits, resolution and PNG/FIG reconstruction PASS';
    receipt.next_allowed_step='STOP and review R1; no R2 authorization';
    local_status(root,receipt,'DONE','COMPLETE'); disp(jsonencode(receipt));
catch ME
    receipt.status='FAILED_STOPPED'; receipt.failed_checks=1; receipt.error=ME.message; receipt.error_identifier=ME.identifier;
    receipt.next_allowed_step='Manual restart after reviewing failure; audited stages must not be recalculated';
    local_status(root,receipt,'STOPPED','FAIL');
    fprintf(2,'%s\n',getReport(ME,'extended','hyperlinks','off')); rethrow(ME);
end
end

function cfg=local_cfg(q,recipe)
cfg=struct('frequency_hz',4000,'c0_mps',1500,'half_width_m',recipe(1), ...
    'points_per_wavelength',recipe(2),'panel_order',8,'self_quadrature_order',recipe(3), ...
    'halfplane_h_m',-2,'window_plateau_ratio',.85,'eta_fn',@(x)closeout_s3a_pm_surface(q.profile,x), ...
    'eta_prime_fn',@(x)local_slope(q.profile,x), ...
    'incident_fn',@(x,z)baffled_strip_piston_field_vertical(x,z,q.cfg.source)- ...
    baffled_strip_piston_field_vertical(x,-z,q.cfg.source),'receiver_x_m',q.x,'receiver_z_m',q.z);
end
function local_verify_package(root,fingerprint)
assert(strcmp(r1_bie_manual_io('sha',fullfile(root,'manifest.json')),fingerprint), ...
    'r1manual:InputDrift','Package manifest changed during execution');
r1_bie_manual_io('verify_manifest',root);
end
function slope=local_slope(profile,x)
[~,slope]=closeout_s3a_pm_surface(profile,x);
end
function c=local_config_record(q,recipe,label,fingerprint)
c=struct('stage',label,'L_m',recipe(1),'ppw',recipe(2),'qa',recipe(3),'panel_order',8, ...
    'h_m',-2,'window_plateau_ratio',.85,'source',q.cfg.source,'seed',260001,'receivers',411, ...
    'phase_contract_version','2.0.0','cache_enabled',true,'manifest_sha256',fingerprint);
end
function gates=local_gates(root,ids,stage)
gates=struct('status','PASS','stage',sprintf('C%d',stage),'sensitivity_gate',1e-4, ...
    'pairs',struct('first',{},'second',{},'relative',{}));
if stage==0, return; end
pairs=[1 2;3 2;2 4;2 5];
for j=1:stage
    a=fullfile(root,'jobs',ids{pairs(j,1)},'data.mat'); b=fullfile(root,'jobs',ids{pairs(j,2)},'data.mat');
    if isfile(a) && isfile(b)
        A=load(a,'field'); B=load(b,'field'); v=local_worst(A.field,B.field);
        gates.pairs(end+1)=struct('first',ids{pairs(j,1)},'second',ids{pairs(j,2)},'relative',v);
        if v>=1e-4
            gates.status='FAIL'; r1_bie_manual_io('write',fullfile(root,'jobs',ids{stage+1},'gates.json'),gates);
        end
        assert(v<1e-4,'r1manual:Convergence','Reference sensitivity failed, pair %d: %.9g',j,v);
    end
end
if stage>=3
    A=load(fullfile(root,'jobs',ids{1},'data.mat'),'field'); B=load(fullfile(root,'jobs',ids{2},'data.mat'),'field');
    C=load(fullfile(root,'jobs',ids{4},'data.mat'),'field');
    gates.grid_nonmonotonic=local_worst(B.field,C.field)>local_worst(A.field,B.field);
    if gates.grid_nonmonotonic
        gates.status='FAIL'; r1_bie_manual_io('write',fullfile(root,'jobs',ids{stage+1},'gates.json'),gates);
    end
    assert(~gates.grid_nonmonotonic, ...
        'r1manual:Nonmonotonic','Grid sensitivity nonmonotonic; stop for review');
end
end
function v=local_worst(a,b)
v=norm(a-b)/max(norm(b),realmin); edges=round(linspace(1,numel(a),5));
for j=1:4, ii=edges(j):edges(j+1); v=max(v,norm(a(ii)-b(ii))/max(norm(b(ii)),realmin)); end
end
function a=local_attempt(root)
base=fullfile(root,'attempts'); if ~isfolder(base), mkdir(base); end
j=1; a=fullfile(base,sprintf('attempt_%04d',j));
while isfolder(a), j=j+1; a=fullfile(base,sprintf('attempt_%04d',j)); end
mkdir(a);
end
function local_guard(baseline)
r=local_resource(baseline);
assert(r.available_bytes>=.3*baseline,'r1manual:RAM','Available RAM fell below 30 percent reserve');
end
function r=local_resource(baseline)
t=fileread('/proc/meminfo'); m=regexp(t,'MemAvailable:\s+(\d+)\s+kB','tokens','once'); assert(~isempty(m));
r=struct('time',local_time(),'available_bytes',str2double(m{1})*1024,'baseline_bytes',baseline,'reserve_fraction',.3);
end
function t=local_time()
t=char(datetime('now','TimeZone','Asia/Shanghai','Format','yyyy-MM-dd HH:mm:ss Z'));
end
function local_stage_file(root,stage)
fid=fopen(fullfile(root,'control','stage_id.txt'),'w'); assert(fid>=0); fprintf(fid,'%s\n',stage); fclose(fid);
end
function local_status(root,receipt,stage,state)
local_stage_file(root,stage); receipt.stage=stage; receipt.stage_state=state; receipt.updated_at=local_time();
receipt.session_elapsed_s=posixtime(datetime('now','TimeZone','UTC'))- ...
    str2double(strtrim(fileread(fullfile(root,'sessions',receipt.session_id,'start_epoch.txt'))));
receipt.resource_usage=local_resource(receipt.RAM_baseline_bytes);
receipt.input_versions_sha256='manifest.json'; receipt.code_versions_sha256='manifest.json';
receipt.artifacts_sha256='jobs/*/committed.json and finalization/committed.json';
receipt.audit_status_and_path='Per-stage audit.json/gates.json/committed.json; finalization/pointer.json after final audits';
receipt.failed_or_unrun_checks='See completed_steps, stage_state and preserved attempt/session failures';
r1_bie_manual_io('write',fullfile(root,'checkpoint.json'),receipt);
end
function local_prepare_final(root,work,ids)
copyfile(fullfile(root,'common.mat'),fullfile(work,'common.mat')); mkdir(fullfile(work,'input'));
copyfile(fullfile(root,'input','*'),fullfile(work,'input'));
for j=1:numel(ids)
    dest=fullfile(work,'jobs',ids{j}); mkdir(dest);
    copyfile(fullfile(root,'jobs',ids{j},'data.mat'),fullfile(dest,'data.mat'));
end
end
function rel=local_relative_path(root,path)
prefix=[root filesep]; assert(strncmp(path,prefix,numel(prefix))); rel=path(numel(prefix)+1:end);
end
function local_final_report(work)
M=readtable(fullfile(work,'model_metrics.csv')); a=load(fullfile(work,'comparison.mat'),'bie_values');
fid=fopen(fullfile(work,'report.md'),'w','n','UTF-8'); assert(fid>=0); c=onCleanup(@()fclose(fid));
fprintf(fid,'# R1-BIE C1–C4 completion\n\nFixed 10 cm strip, 4 kHz, PM seed260001, 411 actual receivers at z=3 m.\n\n');
fprintf(fid,'C0 was imported, never re-solved. C1–C4 and independent stage/final audits passed.\n\n');
fprintf(fid,'| Method | Complex relative L2 | Five-times empirical resolution | Resolved | Peak x (m) |\n|---|---:|---:|---|---:|\n');
for j=[1 3]
    fprintf(fid,'| %s | %.12g | %.12g | %d | %.9g |\n',string(M.method(j)), ...
        M.complex_relative_L2(j),M.five_times_empirical_resolution(j),M.resolved(j),M.peak_x_m(j));
end
fprintf(fid,'\nBIE sensitivities: %.9g, %.9g, %.9g, %.9g (each <1e-4).\n',a.bie_values);
fprintf(fid,'\nSee receiver_fields_and_errors.csv, independent_reference_resolution.csv, audit/ and figures_final/.\n');
fprintf(fid,'\nEmpirical numerical reference for this case, not a mathematical truth bound. STOP for result review; R2 remains unauthorized.\n');
end
function out=local_resume_selftest(root)
% Filesystem-only fixtures: no solver, Green function or acoustic evaluation.
base=fullfile(root,'dry_run_tests'); if ~isfolder(base), mkdir(base); end
test=local_attempt(base); stage=fullfile(test,'fixture_only'); mkdir(stage);
fid=fopen(fullfile(stage,'data.txt'),'w'); fprintf(fid,'NONPHYSICAL RESUME FIXTURE\n'); fclose(fid);
r1_bie_manual_io('write',fullfile(stage,'audit.json'),struct('status','PASS','failed_checks',0));
r1_bie_manual_io('write',fullfile(stage,'gates.json'),struct('status','PASS'));
sig='FIXTURE_NOT_A_RESEARCH_STAGE';
assert(strcmp(r1_bie_manual_io('resume_state',stage,sig),'SOLVE'));
r1_bie_manual_io('seal',stage,'solved.json',sig,{'data.txt'},'SOLVED_UNAUDITED');
assert(strcmp(r1_bie_manual_io('resume_state',stage,sig),'AUDIT_ONLY'));
r1_bie_manual_io('seal',stage,'committed.json',sig,{'data.txt','audit.json','gates.json'},'AUDITED_PASS');
assert(strcmp(r1_bie_manual_io('resume_state',stage,sig),'SKIP_AUDITED'));
stale=false;
try
    r1_bie_manual_io('resume_state',stage,'changed');
catch ME
    stale=strcmp(ME.identifier,'r1manual:StaleCheckpoint');
end
assert(stale);
fid=fopen(fullfile(stage,'data.txt'),'a'); fprintf(fid,'tamper\n'); fclose(fid);
tamper=false;
try
    r1_bie_manual_io('resume_state',stage,sig);
catch ME
    tamper=strcmp(ME.identifier,'r1manual:ArtifactDrift');
end
assert(tamper);
orphan=fullfile(test,'orphan_fixture'); mkdir(orphan);
fid=fopen(fullfile(orphan,'data.mat'),'w'); fprintf(fid,'NOT A MAT FILE - ORPHAN STATE FIXTURE'); fclose(fid);
orphan_rejected=false;
try
    r1_bie_manual_io('resume_state',orphan,sig);
catch ME
    orphan_rejected=strcmp(ME.identifier,'r1manual:OrphanData');
end
assert(orphan_rejected);
out=struct('status','PASS','fixture_only',true,'acoustic_solves',0, ...
    'fresh_solve_state',true,'audit_only_recovery',true,'audited_skip',true, ...
    'changed_predecessor_rejected',stale,'tampered_artifact_rejected',tamper,'unsealed_orphan_rejected',orphan_rejected);
r1_bie_manual_io('write',fullfile(test,'result.json'),out);
end
