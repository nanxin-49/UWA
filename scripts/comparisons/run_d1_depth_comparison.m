function receipt=run_d1_depth_comparison(root,mode)
%RUN_D1_DEPTH_COMPARISON Isolated two-depth, Rx6m, fixed PM reflection screen.
arguments
    root (1,:) char
    mode (1,:) char = 'dry-run'
end
pm_sensitivity_io('absolute',root); assert(ismember(mode,{'dry-run','audit-import','execute'}));
m=jsondecode(fileread(fullfile(root,'manifest.json')));
pm_sensitivity_io('verify',root,m.files);
c=jsondecode(fileread(fullfile(root,'config.json')));
assert(strcmp(c.kind,'D1A_FIXED_PM_DEPTH_COMPARISON') && isequal(c.depths_m(:),[100;1490]));
assert(c.receiver_z_m==6 && c.receiver_halfwidth_m==36 && c.receivers==411);
binding=pm_sensitivity_io('sha',fullfile(root,'manifest.json'));
receipt=struct('status','DRY_RUN_STATIC_INPUTS_ONLY','new_BIE_solves',0,'failed_or_unrun_checks',{{'MATLAB numerical runtime'}}, ...
    'manifest_sha256',binding,'next_allowed_step','User reviews and manually starts isolated remote package');
if strcmp(mode,'dry-run'), disp(jsonencode(receipt)); return; end
if strcmp(mode,'audit-import')
    assert(isfield(c,'recovery')); verify_d1_saved_import(root,fullfile(root,'imported','depth_0100m'),c.recovery);
    receipt.status='IMPORTED_DATA_AND_MODEL_GATES_VERIFIED_NO_SOLVE'; disp(jsonencode(receipt)); return;
end
assert(strcmp(version('-release'),'2025b') && maxNumCompThreads==1);
assert(strcmp(getenv('D1_SUPERVISED_MANIFEST_SHA256'),binding),'Use supervised launcher');
out=fullfile(root,'output'); assert(~isfolder(out),'One shot only; inspect existing task before recovery'); mkdir(out);
diary(fullfile(out,'diary.log')); cleanup=onCleanup(@()diary('off')); t0=tic;
receipt.status='RUNNING'; receipt.completed_steps={}; receipt.failed_or_unrun_checks={'D1 stages incomplete'};
receipt.imported_BIE_solves=0;
if isfield(c,'recovery'), receipt.imported_BIE_solves=10-c.recovery.remaining_BIE_solves; end
receipt.config=c; receipt.code_versions_sha256=m.files; receipt.input_versions_sha256=m.files;
receipt.execution_mode='remote R2025b single process single computational thread'; receipt.matlab_release=version;
receipt.resource_usage=struct('reserve_fraction',.3,'runtime_receipt','../sessions/*/runtime.json','resource_log','../sessions/*/resource.tsv');
ids={'bie_coarse','bie_main','bie_self','bie_space','bie_window'};
recipes=[60 10 96;60 12 96;60 12 64;60 14 96;70 12 96];
try
    for depth=c.depths_m(:).'
        work=fullfile(out,sprintf('depth_%04dm',depth));
        imported=isfield(c,'recovery') && depth==100;
        imported_depth=imported && strcmp(c.recovery.kind,'IMPORT_AUDITED_TX100_NO_RESOLVE');
        if imported
            copyfile(fullfile(root,'imported','depth_0100m'),work);
            verify_d1_saved_import(root,work,c.recovery);
        else
            mkdir(work); mkdir(fullfile(work,'input'));
        end
        pm_sensitivity_io('write',fullfile(work,'started.json'),struct('depth_m',depth,'utc',char(datetime('now','TimeZone','UTC'))));
        copyfile(fullfile(root,'inputs','coefficients.csv'),fullfile(work,'input','coefficients.csv'));
        copyfile(fullfile(root,'inputs','receivers.csv'),fullfile(work,'input','receivers.csv'));
        copyfile(fullfile(root,'inputs','policy.json'),fullfile(work,'policy.json'));
        frozen=struct('source',c.source,'seed',260001,'modes',12,'U_PM_mps',6,'case_id',sprintf('Tx%d_Rx6',depth),'family','original12_no_resampling');
        frozen.source.z_source_m=depth;
        pm_sensitivity_io('write',fullfile(work,'input','frozen_config.json'),frozen);
        if ~imported
            local_source_preflight(work,frozen.source);
            local_guard(); run_d1_model_job(work,'models'); diary(fullfile(out,'diary.log'));
        end
        local_seal(work,'MODELS',binding); q=load(fullfile(work,'common.mat')); fields=cell(5,1);
        for j=1:5
            local_guard(); pm_sensitivity_io('verify',root,m.files);
            dest=fullfile(work,'jobs',ids{j});
            reuse=imported && (j==1 || imported_depth);
            if reuse
                saved_sha=pm_sensitivity_io('sha',fullfile(dest,'data.mat'));
                expected=c.recovery.C0_data_sha256;
                if imported_depth, expected=c.recovery.layer_data_sha256{j}; end
                assert(strcmp(saved_sha,expected));
                old=load(fullfile(dest,'data.mat')); field=old.field; fields{j}=field;
                % Re-audit this same saved density. Never call the solver for C0.
            else
            assert(~isfolder(dest),'D1:ExistingStage','Refuse to repeat an existing stage'); mkdir(dest);
            pm_sensitivity_io('write',fullfile(dest,'solve_invocation.json'),struct('depth_m',depth,'stage',j-1,'calls',1,'manifest_sha256',binding));
            cfg=struct('frequency_hz',4000,'c0_mps',1500,'half_width_m',recipes(j,1), ...
                'points_per_wavelength',recipes(j,2),'panel_order',8,'self_quadrature_order',recipes(j,3), ...
                'halfplane_h_m',-2,'window_plateau_ratio',.85,'eta_fn',@(x)closeout_s3a_pm_surface(q.profile,x), ...
                'eta_prime_fn',@(x)local_slope(q.profile,x),'incident_fn',@(x,z) ...
                baffled_strip_piston_field_vertical(x,z,q.cfg.source)-baffled_strip_piston_field_vertical(x,-z,q.cfg.source), ...
                'receiver_x_m',q.x,'receiver_z_m',q.z,'enable_static_cache',true,'execution_guard_fn',@local_guard);
            receipt.new_BIE_solves=receipt.new_BIE_solves+1; boundary=solve_helmholtz_bie_halfplane_vertical(cfg);
            field=q.flat+boundary.receiver_field; fields{j}=field;
            incident=cfg.incident_fn(boundary.node_x_m,boundary.node_z_m);
            numeric_config=rmfield(cfg,{'eta_fn','eta_prime_fn','incident_fn','execution_guard_fn'});
            save(fullfile(dest,'data.mat'),'boundary','field','incident','numeric_config','-v7');
            pm_sensitivity_io('write',fullfile(dest,'timings.json'),boundary.cost_meta.timings);
            write_d1_csv(table(q.x,q.z,real(field),imag(field),'VariableNames',{'x_m','z_m','real_Pa','imag_Pa'}),fullfile(dest,'points.csv'));
            end
            audit_r1_bie_manual_stage(dest,fullfile(work,'common.mat'),recipes(j,:));
            if reuse, assert(strcmp(pm_sensitivity_io('sha',fullfile(dest,'data.mat')),saved_sha)); end
            if j>=4
                policy=jsondecode(fileread(fullfile(work,'policy.json'))); g=d1_bie_gates(fields(1:j),q.x,q.z,policy);
                pm_sensitivity_io('write',fullfile(dest,'gates.json'),g); assert(strcmp(g.status,'PASS'),'D1 BIE sensitivity gate');
            elseif j>1
                if j==2, v=local_worst(fields{1},fields{2}); else, v=local_worst(fields{3},fields{2}); end
                pm_sensitivity_io('write',fullfile(dest,'gates.json'),struct('worst_sensitivity',v,'gate',1e-4)); assert(v<1e-4);
            end
            local_seal(dest,sprintf('C%d',j-1),binding);
            receipt.completed_steps{end+1}=sprintf('Tx%d_C%d',depth,j-1); local_checkpoint();
        end
        if ~imported_depth, run_d1_model_job(work,'analyze'); diary(fullfile(out,'diary.log')); end
        audit_d1_models(work); diary(fullfile(out,'diary.log')); audit_d1_final(work);
        if isfield(c,'final_figures_mode')
            pm_sensitivity_io('write',fullfile(work,'figures_pending.json'),struct('status','LOCAL_FIGURE_AUDIT_REQUIRED','required_differing_pixel_channels',0));
            local_seal(work,'DEPTH_NUMERICAL_AUDITED',binding);
        else
            render_d1_depth_results(work); diary(fullfile(out,'diary.log')); local_guard();
            local_seal(work,'DEPTH_AUDITED',binding);
        end
        receipt.completed_steps{end+1}=sprintf('Tx%d_FINAL',depth); local_checkpoint();
    end
    a=readtable(fullfile(out,'depth_0100m','model_metrics.csv')); a.tx_depth_m=100*ones(height(a),1);
    b=readtable(fullfile(out,'depth_1490m','model_metrics.csv')); b.tx_depth_m=1490*ones(height(b),1);
    write_d1_csv([a;b],fullfile(out,'depth_comparison.csv'));
    methods=string(a.method([1 3])); change=b.complex_relative_L2([1 3])-a.complex_relative_L2([1 3]);
    resolution_bound=a.five_times_empirical_resolution([1 3])+b.five_times_empirical_resolution([1 3]);
    write_d1_csv(table(methods,change,resolution_bound,abs(change)>resolution_bound,'VariableNames', ...
        {'method','E1490_minus_E100','conservative_change_resolution','change_resolved'}),fullfile(out,'depth_change.csv'));
    fid=fopen(fullfile(out,'report.md'),'w','n','UTF-8'); assert(fid>=0);
    fprintf(fid,'# D1a fixed-surface depth comparison\n\n0.10 m 2D strip, 4 kHz; Tx100/1490 m, Rx6 m; original U6/N12/260001/L160/C2.\n\n');
    fprintf(fid,'411 direct points on [-36,36] m. BIE C3 empirical reference, C4 window check. No R1/R2 density imported.\n\n');
    fprintf(fid,'|Tx m|Method|Complex L2|Amplitude L2|Phase RMS rad|5x resolution|Resolved|Phase valid|\n|--:|--|--:|--:|--:|--:|--|--:|\n');
    T=[a;b]; for j=1:height(T), fprintf(fid,'|%g|%s|%.12g|%.12g|%.12g|%.12g|%d|%d|\n',T.tx_depth_m(j),string(T.method(j)),T.complex_relative_L2(j),T.amplitude_relative_L2(j),T.phase_rms_rad(j),T.five_times_empirical_resolution(j),T.resolved(j),T.phase_valid_points(j)); end
    fprintf(fid,'\nPer-depth common/comparison MAT, model_metrics/receiver_fields_and_errors CSV, audits and figures contain full fields, masks and spatial errors. Inspect them before physical interpretation. One paired realization; receiver points are not independent sea realizations.\n\nUniform lossless water, upward ideal baffled source, sound-soft surface on flat exterior. No bottom/downward radiation/repeated surface-bottom paths/structure. This is not a complete real deep-sea multipath model. D1b is not authorized.\n'); fclose(fid);
    pm_sensitivity_io('verify',root,m.files);
    receipt.status='D1A_COMPLETE_WAITING_REVIEW'; receipt.failed_or_unrun_checks={}; receipt.audit_status_and_path='Both depths stage audits + independent final metrics + figure audits passed';
    if isfield(c,'final_figures_mode')
        receipt.status='D1A_NUMERICS_COMPLETE_LOCAL_FIGURES_REQUIRED';
        receipt.failed_or_unrun_checks={'Local PNG/FIG reconstruction and zero-pixel audit required'};
        receipt.audit_status_and_path='Both depth numerical audits passed; local figure gate pending';
    end
    receipt.next_allowed_step='STOP; user reviews D1a before any D1b design authorization'; local_checkpoint();
catch ME
    receipt.status='FAILED_STOPPED'; receipt.failed_or_unrun_checks={ME.identifier,ME.message}; local_checkpoint();
    receipt.next_allowed_step='STOP: inspect checkpoint and preserve every solved stage; never repeat a completed solve'; local_checkpoint();
    fprintf(2,'%s\n',getReport(ME,'extended','hyperlinks','off')); rethrow(ME);
end
    function local_checkpoint()
        receipt.elapsed_s=toc(t0); receipt.artifacts_sha256=local_files(out);
        if isfield(c,'recovery'), receipt.cumulative_elapsed_s=receipt.elapsed_s+c.recovery.prior_supervisor_elapsed_s; end
        pm_sensitivity_io('write',fullfile(out,'checkpoint.json'),receipt);
    end
end
function local_guard()
base=str2double(getenv('D1_RAM_BASELINE_BYTES')); assert(isfinite(base) && base>0);
t=regexp(fileread('/proc/meminfo'),'MemAvailable:\s+(\d+)\s+kB','tokens','once');
assert(str2double(t{1})*1024>=.3*base,'D1:RAM','30 percent initial available RAM reserve breached');
end
function y=local_slope(p,x)
[~,y]=closeout_s3a_pm_surface(p,x);
end
function local_source_preflight(work,s)
p=readtable(fullfile(work,'input','coefficients.csv'));
x=linspace(-70,70,14001)'; [z,slope]=closeout_s3a_pm_surface(p,x);
assert(all(isfinite([z;slope])) && min(z)>-2 && max(z)<6,'Surface and auxiliary-plane clearance');
assert(all(z(abs(x)>=50)==0) && all(slope(abs(x)>=50)==0),'Original compact support');
xx=[-50;-36;-20;0;20;36;50]; zz=closeout_s3a_pm_surface(p,xx);
[u,ux,uz]=baffled_strip_piston_field_vertical(xx,zz,s,128);
A=d0_incident_spectrum_reference(xx,zz,s,16,8); B=d0_incident_spectrum_reference(xx,zz,s,32,4);
pr=norm(u-B.total(:,1))/norm(B.total(:,1)); gr=norm([ux uz]-B.total(:,2:3),'fro')/norm(B.total(:,2:3),'fro');
qr=norm(A.total-B.total,'fro')/norm(B.total,'fro');
pm_sensitivity_io('write',fullfile(work,'source_geometry_audit.json'),struct('source',s, ...
    'source_pressure_relative',pr,'source_gradient_relative',gr,'independent_quadrature_relative',qr, ...
    'minimum_surface_z_m',min(z),'maximum_surface_z_m',max(z),'gate',1e-7,'geometry_scope','sampled original compact C2; solver checks all numerical nodes'));
save(fullfile(work,'source_geometry_audit.mat'),'xx','zz','u','ux','uz','A','B','-v7');
assert(all([pr gr qr]<1e-7),'Independent source/gradient quadrature gate');
end
function v=local_worst(a,b)
v=norm(a-b)/max(norm(b),realmin); edges=round(linspace(1,411,5));
for j=1:4, ii=edges(j):edges(j+1); v=max(v,norm(a(ii)-b(ii))/max(norm(b(ii)),realmin)); end
end
function f=local_files(root)
d=dir(fullfile(root,'**','*')); paths={};
for j=1:numel(d)
    if ~d(j).isdir && ~strcmp(d(j).name,'diary.log') && ...
            ~(strcmp(d(j).name,'checkpoint.json') && strcmp(d(j).folder,root)) && ~endsWith(d(j).name,'.pending')
        p=fullfile(d(j).folder,d(j).name); paths{end+1}=strrep(p(numel(root)+2:end),'\','/'); %#ok<AGROW>
    end
end
f=pm_sensitivity_io('manifest',root,paths);
end
function local_seal(root,label,binding)
f=local_files(root); pm_sensitivity_io('verify',root,f);
input_root=root;
while ~isfile(fullfile(input_root,'manifest.json'))
    next=fileparts(input_root); assert(~strcmp(next,input_root),'Missing isolated manifest root'); input_root=next;
end
m=jsondecode(fileread(fullfile(input_root,'manifest.json'))); c=jsondecode(fileread(fullfile(input_root,'config.json')));
resource=struct('baseline_available_bytes',str2double(getenv('D1_RAM_BASELINE_BYTES')), ...
    'reserve_fraction',.3,'measured_meminfo',fileread('/proc/meminfo'),'time_utc',char(datetime('now','TimeZone','UTC')), ...
    'full_time_RSS_receipt','sessions/*/resource.tsv; sessions/*/exit.json');
status='AUDITED_PASS'; unrun={};
if strcmp(label,'MODELS'), status='NUMERICAL_GATES_PASS_PENDING_INDEPENDENT_AUDIT'; unrun={'Independent source, PS Green, FullFK and MAT/CSV audits at depth finalization'}; end
if strcmp(label,'DEPTH_NUMERICAL_AUDITED'), unrun={'Local PNG/FIG reconstruction and zero-pixel audit required'}; end
pm_sensitivity_io('write',fullfile(root,'checkpoint.json'),struct('status',status,'stage',label, ...
    'manifest_sha256',binding,'config',c,'input_code_versions_sha256',m.files,'resource_usage',resource, ...
    'execution_mode','remote R2025b single computational thread','matlab_release',version, ...
    'artifacts_sha256',f,'failed_or_unrun_checks',{unrun}, ...
    'next_allowed_step','Next frozen stage only; one-shot failure requires manual review'));
end
