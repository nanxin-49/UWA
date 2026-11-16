function receipt=run_r1_bie_cache_continuation(root)
%RUN_R1_BIE_CACHE_CONTINUATION Authorized, budget-gated C0--C4 continuation.
% Prior 439 s is immutable. All new process overhead shares the 661 s cap.
arguments
    root (1,:) char
end
assert(isfolder(root)); clock0=tic;
diary(fullfile(root,'execution','matlab_diary.log')); cleanup=onCleanup(@()diary('off')); %#ok<NASGU>
epoch=str2double(strtrim(fileread(fullfile(root,'execution','process_start_epoch.txt'))));
baseline=str2double(strtrim(fileread(fullfile(root,'execution','starting_available_bytes.txt'))));
receipt=struct('run_id','20261116_r1_bie_cache_01','status','RUNNING','failed_checks',0, ...
    'completed_steps',{{}},'matlab_release',version,'pid',feature('getpid'), ...
    'prior_process_s',439,'cumulative_cap_s',1100,'new_process_cap_s',661, ...
    'config',struct('source_width_m',.1,'frequency_hz',4000,'seed',260001, ...
    'receivers',411,'receiver_z_m',3,'cache_default',false,'singleCompThread',true), ...
    'code_versions_sha256','snapshot.sha256','input_versions_sha256','snapshot.sha256', ...
    'next_allowed_step','S0 only, followed by C0 on equivalence PASS');
local_checkpoint(root,receipt,epoch);
try
assert(strcmp(version('-release'),'2025b'));
q=load(fullfile(root,'common.mat')); source=q.cfg.source;
assert(source.total_width_m==.1 && source.frequency_hz==4000 && numel(q.x)==411 && all(q.z==3));
stage_deadline=epoch+40; local_deadline(root,stage_deadline,'S0');
code_analysis=struct;
for name={'solve_helmholtz_bie_halfplane_vertical','solve_helmholtz_bie_halfplane_cached_vertical','run_r1_bie_cache_continuation'}
    code_analysis.(name{1})=checkcode(which(name{1}),'-id');
end
save(fullfile(root,'equivalence','code_analysis.mat'),'code_analysis','-v7');
test_rows=table; case_id=0;
for ppw=[.2 .4]
    for qa=[64 96]
        local_guard(epoch,stage_deadline,baseline);
        case_id=case_id+1; cfg=local_cfg(q,60,ppw,qa);
        cfg.return_system_for_validation=true;
        original=solve_helmholtz_bie_halfplane_original_r1(cfg);
        implicit=solve_helmholtz_bie_halfplane_vertical(cfg);
        cfg.enable_static_cache=false; disabled=solve_helmholtz_bie_halfplane_vertical(cfg);
        assert(isequaln(original,implicit) && isequaln(original,disabled),'Disabled/default path changed');
        cfg.enable_static_cache=true;
        cfg.execution_guard_fn=@()local_guard(epoch,stage_deadline,baseline);
        cached=solve_helmholtz_bie_halfplane_vertical(cfg);
        differences=[local_relative(cached.validation_system.A,original.validation_system.A), ...
            local_relative(cached.validation_system.rhs,original.validation_system.rhs), ...
            local_relative(cached.density,original.density), ...
            local_relative(cached.receiver_field,original.receiver_field), ...
            norm(cached.boundary_total_residual-original.boundary_total_residual)/ ...
            max(norm(local_input(original.boundary_check_x_m,q.profile,source)),realmin)];
        assert(all(differences<=[1e-12 1e-12 1e-11 1e-11 1e-10]),'Small cache equivalence gate');
        test_rows=[test_rows;table(case_id,original.node_count,ppw,qa,true, ...
            differences(1),differences(2),differences(3),differences(4),differences(5), ...
            'VariableNames',{'case_id','nodes','ppw','qa','implicit_disabled_exact', ...
            'A_relative','rhs_relative','density_relative','receiver_relative','boundary_difference_over_drive'})]; %#ok<AGROW>
        save(fullfile(root,'equivalence',sprintf('case_%d.mat',case_id)), ...
            'original','implicit','disabled','cached','differences','-v7');
        local_guard(epoch,stage_deadline,baseline);
    end
end
writetable(test_rows,fullfile(root,'equivalence','metrics.csv'));
receipt.completed_steps{end+1}='S0 four small original/implicit/disabled/cached cases PASS';
receipt.S0_process_elapsed_s=local_elapsed(epoch);
assert(receipt.S0_process_elapsed_s<=40,'S0 40-second cap');
local_checkpoint(root,receipt,epoch);
% Full C0 calibration is the only high-cost solve allowed before cost PASS.
stage_deadline=min(epoch+661,local_now()+120); local_deadline(root,stage_deadline,'C0');
[r,tC0]=local_solve(q,root,'bie_coarse',60,10,96,epoch,stage_deadline,baseline);
old=load(fullfile(root,'parent','C0_data.mat'));
equivalence=struct('density_relative',local_relative(r.boundary.density,old.boundary.density), ...
    'correction_relative',local_relative(r.boundary.receiver_field,old.boundary.receiver_field), ...
    'reflected_relative',local_relative(r.field,old.field), ...
    'node_geometry_exact',isequaln(r.boundary.node_x_m,old.boundary.node_x_m) && ...
    isequaln(r.boundary.node_z_m,old.boundary.node_z_m) && isequaln(r.boundary.node_normal,old.boundary.node_normal));
local_json(fullfile(root,'equivalence','C0_full.json'),equivalence);
assert(equivalence.density_relative<=1e-11 && equivalence.correction_relative<=1e-11 && ...
    equivalence.reflected_relative<=1e-11 && equivalence.node_geometry_exact,'C0 historical equivalence gate');
receipt.completed_steps{end+1}='C0 full configuration historical replay PASS';
receipt.C0_equivalence=equivalence; receipt.C0_function_s=tC0;
receipt.C0_phase_timings=r.boundary.cost_meta.timings;
receipt.C0_speedup_vs_old_function=254.689257/tC0;
remaining_factors=[1.2^3 1.2^3 1.4^3 1.4^3];
remaining_names={'bie_main','bie_self','bie_space','bie_window'};
recipes=[60 12 96;60 12 64;60 14 96;70 12 96];
forecast=1.25*tC0*sum(remaining_factors)+110;
receipt.cost_gate=struct('elapsed_new_process_s',local_elapsed(epoch),'remaining_cap_s',661-local_elapsed(epoch), ...
    'C1_C4_and_closeout_forecast_s',forecast,'margin_factor',1.25,'closeout_reserve_s',110, ...
    'pass',local_elapsed(epoch)+forecast<=661,'model','calibrated N^3; no qa64 speed credit');
local_checkpoint(root,receipt,epoch);
if ~receipt.cost_gate.pass
    receipt.status='STOPPED_COST_GATE_AFTER_C0';
    receipt.unrun_checks={'C1 grid12','C2 quadrature64','C3 grid14','C4 window70','qualified BIE model comparison'};
    receipt.next_allowed_step='STOP and review actual C0 timings; no C1-C4 dispatched';
else
    % Revalidate all remaining cost and reference gates at each boundary.
    calibration=tC0;
    for j=1:4
        forecast=1.25*calibration*sum(remaining_factors(j:end))+110;
        admission=struct('job',remaining_names{j},'elapsed_new_s',local_elapsed(epoch), ...
            'all_remaining_and_closeout_forecast_s',forecast,'pass',local_elapsed(epoch)+forecast<=661);
        local_json(fullfile(root,'execution',sprintf('admission_C%d.json',j)),admission);
        if ~admission.pass
            receipt.status='STOPPED_UPDATED_COST_GATE'; receipt.next_allowed_step='Review saved boundary; no retry or dropped levels'; break;
        end
        stage_deadline=epoch+661-110; local_deadline(root,stage_deadline,sprintf('C%d',j));
        rc=recipes(j,:);
        [~,dt]=local_solve(q,root,remaining_names{j},rc(1),rc(2),rc(3),epoch,stage_deadline,baseline);
        calibration=max(calibration,dt/remaining_factors(j));
        local_reference_gates(root,j);
        receipt.completed_steps{end+1}=sprintf('C%d residual and available convergence checks PASS',j);
        local_checkpoint(root,receipt,epoch);
    end
    if isfile(fullfile(root,'jobs','bie_window','data.mat')) && strcmp(receipt.status,'RUNNING')
        receipt.status='REFERENCE_LAYERS_COMPLETE';
    end
end
% Independent C0 audit always executes. It does not re-solve an operator.
local_deadline(root,epoch+661,'saved_data_closeout'); local_guard(epoch,epoch+661,baseline);
local_audit_C0(root,q);
receipt.completed_steps{end+1}='Independent C0 saved-density receiver and MAT/CSV audit PASS';
if strcmp(receipt.status,'REFERENCE_LAYERS_COMPLETE')
    run_piston_pm_r1_job(root,'analyze');
    audit_piston_pm_r1(root); render_piston_pm_r1(root,'figures_final');
    receipt.status='R1_REFERENCE_AND_COMPARISON_COMPLETE';
    receipt.completed_steps{end+1}='Independent reference/comparison and final figure audits PASS';
    receipt.next_allowed_step='STOP for R1 result review; R2 remains unauthorized';
end
receipt.elapsed_function_s=toc(clock0); local_checkpoint(root,receipt,epoch); disp(jsonencode(receipt));
catch ME
    receipt.status='FAILED_OR_STOPPED'; receipt.failed_checks=1; receipt.error=ME.message;
    receipt.error_identifier=ME.identifier; receipt.elapsed_function_s=toc(clock0);
    receipt.next_allowed_step='STOP; inspect original job/checkpoint, no automatic retry';
    local_checkpoint(root,receipt,epoch); fprintf(2,'%s\n',getReport(ME,'extended','hyperlinks','off')); rethrow(ME);
end
end

function [r,elapsed]=local_solve(q,root,name,L,ppw,qa,epoch,deadline,baseline)
dest=fullfile(root,'jobs',name); assert(~isfolder(dest)); mkdir(dest);
local_guard(epoch,deadline,baseline); cfg=local_cfg(q,L,ppw,qa); cfg.enable_static_cache=true;
cfg.execution_guard_fn=@()local_guard(epoch,deadline,baseline);
clock=tic; boundary=solve_helmholtz_bie_halfplane_vertical(cfg); elapsed=toc(clock);
field=q.flat+boundary.receiver_field;
incident=cfg.incident_fn(boundary.node_x_m,boundary.node_z_m);
numeric_config=rmfield(cfg,{'eta_fn','eta_prime_fn','incident_fn','execution_guard_fn'});
save(fullfile(dest,'data.mat'),'boundary','field','incident','numeric_config','-v7');
writetable(table(q.x,q.z,real(field),imag(field),'VariableNames',{'x_m','z_m','real_Pa','imag_Pa'}),fullfile(dest,'points.csv'));
metrics=struct('node_count',boundary.node_count,'linear_residual',max(boundary.linear_residual), ...
    'offgrid_boundary_residual',max(boundary.offgrid_boundary_residual),'boundary_midpoints',numel(boundary.boundary_check_x_m));
rr=struct('status','DATA_SAVED','function_s',elapsed,'timings',boundary.cost_meta.timings,'metrics',metrics);
local_json(fullfile(dest,'receipt.json'),rr);
assert(boundary.finite && all(boundary.linear_residual<1e-10) && all(boundary.offgrid_boundary_residual<1e-4),'BIE residual gate');
rr.status='RESIDUAL_GATES_PASS'; local_json(fullfile(dest,'receipt.json'),rr);
r=struct('field',field,'boundary',boundary);
end

function cfg=local_cfg(q,L,ppw,qa)
cfg=struct('frequency_hz',4000,'c0_mps',1500,'half_width_m',L, ...
    'points_per_wavelength',ppw,'panel_order',8,'self_quadrature_order',qa, ...
    'halfplane_h_m',-2,'window_plateau_ratio',.85, ...
    'eta_fn',@(x)closeout_s3a_pm_surface(q.profile,x), ...
    'eta_prime_fn',@(x)local_slope(q.profile,x), ...
    'incident_fn',@(x,z)baffled_strip_piston_field_vertical(x,z,q.cfg.source)- ...
    baffled_strip_piston_field_vertical(x,-z,q.cfg.source), ...
    'receiver_x_m',q.x,'receiver_z_m',q.z);
end
function p=local_input(x,profile,source)
z=closeout_s3a_pm_surface(profile,x);
p=baffled_strip_piston_field_vertical(x,z,source)-baffled_strip_piston_field_vertical(x,-z,source);
end
function s=local_slope(profile,x)
[~,s]=closeout_s3a_pm_surface(profile,x);
end
function local_reference_gates(root,stage)
pairs={'bie_coarse','bie_main';'bie_self','bie_main';'bie_main','bie_space';'bie_main','bie_window'};
for j=1:stage
    a=fullfile(root,'jobs',pairs{j,1},'data.mat'); b=fullfile(root,'jobs',pairs{j,2},'data.mat');
    if isfile(a) && isfile(b)
        x=load(a,'field'); y=load(b,'field'); value=local_worst(x.field,y.field);
        local_json(fullfile(root,'execution',sprintf('convergence_pair_%d.json',j)),struct('relative',value,'gate',1e-4));
        assert(value<1e-4,'BIE convergence gate');
    end
end
if stage>=3
    a=load(fullfile(root,'jobs','bie_coarse','data.mat'),'field');
    b=load(fullfile(root,'jobs','bie_main','data.mat'),'field');
    c=load(fullfile(root,'jobs','bie_space','data.mat'),'field');
    first=local_worst(a.field,b.field); second=local_worst(b.field,c.field);
    assert(second<=first,'Nonmonotonic grid sensitivity; stop for review');
end
end
function local_audit_C0(root,q)
dest=fullfile(root,'C0_independent_audit'); mkdir(dest); b=load(fullfile(root,'jobs','bie_coarse','data.mat'));
bn=b.boundary; n=8; j=(1:n-1)'; beta=j./sqrt(4*j.^2-1); [V,D]=eig(diag(beta,1)+diag(beta,-1));
[g,ix]=sort(diag(D)); w=2*(V(1,ix).^2)';
edges=linspace(-60,60,bn.panel_count+1)'; centers=(edges(1:end-1)+edges(2:end))/2;
half=60/bn.panel_count; xs=reshape((centers+half*g.').',[],1); qw=repmat(half*w,bn.panel_count,1);
[zs,sp]=closeout_s3a_pm_surface(q.profile,xs); jac=hypot(1,sp); nx=-sp./jac; nz=1./jac;
r=abs(xs)/60; window=ones(size(r)); window(r>=1)=0; ii=r>.85 & r<1;
t=(r(ii)-.85)/.15; left=exp(-1./t); right=exp(-1./(1-t)); window(ii)=right./(left+right);
field=q.flat; k=2*pi*4000/1500;
for first=1:32:411
    ii=first:min(first+31,411); dx=q.x(ii)-xs.'; dz=q.z(ii)-zs.'; dzI=q.z(ii)+4+zs.';
    R=hypot(dx,dz); RI=hypot(dx,dzI);
    S=1i/4*(besselh(0,1,k*R)-besselh(0,1,k*RI));
    K=1i*k/4*(besselh(1,1,k*R).*(dx.*nx.'+dz.*nz.')./R- ...
        besselh(1,1,k*RI).*(dx.*nx.'-dzI.*nz.')./RI);
    field(ii)=field(ii)+(K-1i*k*S)*(qw.*jac.*window.*bn.density);
end
T=readtable(fullfile(root,'jobs','bie_coarse','points.csv')); csv=T.real_Pa+1i*T.imag_Pa;
relative=local_relative(field,b.field); csv_abs=max(abs(csv-b.field));
assert(relative<1e-11 && csv_abs<1e-11,'Independent C0 reconstruction failed');
save(fullfile(dest,'reconstruction.mat'),'field','-v7');
local_json(fullfile(dest,'receipt.json'),struct('status','PASS','failed_checks',0, ...
    'density_receiver_relative',relative,'csv_max_abs_Pa',csv_abs,'reference_convergence_qualified',false));
end
function local_guard(epoch,deadline,baseline)
assert(local_now()<min(epoch+661,deadline),'R1:TimeGate','Runtime deadline reached');
text=fileread('/proc/meminfo'); m=regexp(text,'MemAvailable:\s+(\d+)\s+kB','tokens','once');
assert(~isempty(m) && str2double(m{1})*1024>=.3*baseline,'R1:RAMGate','30 percent RAM reserve breached');
end
function local_deadline(root,deadline,stage)
fid=fopen(fullfile(root,'execution','stage_deadline_epoch.txt'),'w'); fprintf(fid,'%.6f\n',deadline); fclose(fid);
local_json(fullfile(root,'execution','current_stage.json'),struct('stage',stage,'deadline_epoch',deadline));
end
function local_checkpoint(root,r,epoch)
r.new_process_elapsed_s=local_elapsed(epoch); r.cumulative_process_s=439+r.new_process_elapsed_s;
r.resource_usage=struct('prior_process_s',439,'new_process_elapsed_s',r.new_process_elapsed_s,'cap_s',1100);
r.artifacts_sha256='final_outputs.sha256'; r.audit_status_and_path='See equivalence/ and C0_independent_audit/';
r.failed_or_unrun_checks='Status and completed_steps list the executed subset; full reference unqualified until C1-C4 all pass';
local_json(fullfile(root,'checkpoint.json'),r);
end
function t=local_now()
t=posixtime(datetime('now','TimeZone','UTC'));
end
function t=local_elapsed(epoch)
t=local_now()-epoch;
end
function v=local_relative(a,b)
v=norm(a(:)-b(:))/max(norm(b(:)),realmin);
end
function v=local_worst(a,b)
v=local_relative(a,b); edges=round(linspace(1,numel(a),5));
for j=1:4, ii=edges(j):edges(j+1); v=max(v,local_relative(a(ii),b(ii))); end
end
function local_json(p,s)
fid=fopen(p,'w'); assert(fid>0); c=onCleanup(@()fclose(fid)); fprintf(fid,'%s',jsonencode(s,PrettyPrint=true));
end
