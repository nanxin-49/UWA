function receipt = run_two_dimensional_closeout_s2c_b_bie_job(run_dir,job_id,available_physical_bytes)
%RUN_TWO_DIMENSIONAL_CLOSEOUT_S2C_B_BIE_JOB One guarded K_HIGH BIE solve.
arguments
    run_dir (1,:) char
    job_id (1,:) char
    available_physical_bytes (1,1) double {mustBePositive}
end
ids={'main96','self64','space14','window60','window70'};
levels=find(strcmp(ids,job_id)); assert(isscalar(levels),'Only five frozen K_HIGH levels');
assert(isfolder(run_dir) && contains(run_dir,'20261023_s2c_b_khigh_'));
files=cellfun(@(q)fullfile(run_dir,'jobs',['bie_' q '.mat']),ids,'UniformOutput',false);
assert(~isfile(files{levels}),'No duplicate K_HIGH solve');
for j=1:levels-1, assert(isfile(files{j}),'Required earlier reference level absent'); end
for j=levels+1:5, assert(~isfile(files{j}),'Unexpected later reference level'); end
gate=closeout_s2c_b_effective_gate(run_dir);
authorization=jsondecode(fileread(fullfile(run_dir,'protocol','authorization.json')));
assert(strcmp(authorization.authorized_case,'K_HIGH') && authorization.explicit_user_authorization);
root=fileparts(fileparts(fileparts(fileparts(run_dir))));
prior=fullfile(root,'results/comparisons/two_dimensional_channel_closeout/20261023_s2c_a_1642');
formal=fullfile(root,'results/comparisons/two_dimensional_channel_closeout/20261023_s2b_1305/formal_reference_v2');
cfg_path=fullfile(prior,'protocol','K_HIGH_configuration.json');
assert(strcmp(closeout_s2b_hash(cfg_path),authorization.frozen_configuration_sha256));
src_path=fullfile(prior,'protocol','frozen_source.mat');
srcdata=load(src_path,'src'); src=srcdata.src;
f=load(fullfile(formal,'reference_data.mat'),'x','z'); x=f.x(:); z=f.z(:);
assert(numel(x)==417 && numel(z)==417 && all(z(1:411)==3));
nodes=[3840 3840 4480 3200 3736];
L=[60 60 60 60 70]; ppw=[12 12 14 10 10]; qa=[96 64 96 96 96];
historical=[242.9793715 276.6899304 289.741466 ...
    242.9793715*(3200/3840)^2 242.9793715*(3736/3840)^2];
upper=1.2*historical;
actual=zeros(1,5);
for j=1:levels-1
    q=load(files{j},'receipt'); assert(strcmp(q.receipt.status,'COMPLETE'));
    actual(j)=q.receipt.wall_seconds;
end
if levels>1
    for j=2:5
        if j>=levels
            scale=(nodes(j)/nodes(1))^2;
            if j==2, scale=historical(2)/historical(1); end
            upper(j)=max(upper(j),1.2*actual(1)*scale);
        end
    end
end
upper(1:levels-1)=actual(1:levels-1);
projected=sum(upper)+gate.other_physics_reserve_s+gate.failure_contingency_s;
remaining=gate.core_remaining_s-projected;
assert(projected<=gate.subbudget_s && remaining>=gate.S3A_reserve_s, ...
    'Whole five-level cost projection failed before next BIE job');
assert(all(upper(levels:end)<gate.single_job_cap_s));
peak=4*16*nodes.^2+128*1024^2;
assert(all(peak(levels:end)<=.7*available_physical_bytes), ...
    'Current physical RAM cannot support all remaining frozen BIE levels');
assert(gate.charged_prior_seconds+projected<=7200);
solver=fullfile(root,'scripts/validation/support/solve_helmholtz_bie_halfplane_vertical.m');
surface=fullfile(root,'scripts/validation/support/closeout_s2c_b_surface.m');
incident=fullfile(root,'scripts/validation/support/fk_single_receiver_incident.m');
solver_hash=closeout_s2b_hash(solver);
preflight=struct('job_id',job_id,'status','PREFLIGHT_PASS','new_case','K_HIGH', ...
    'available_physical_bytes',available_physical_bytes,'measurement','Windows Memory Available Bytes immediately before MCP job', ...
    'estimated_increment_bytes',peak(levels),'maximum_remaining_increment_bytes',max(peak(levels:end)), ...
    'reserve_fraction',.3,'projected_all_K_HIGH_core_s',projected, ...
    'authorized_subbudget_s',gate.subbudget_s,'budget_amendment_sha256',gate.budget_amendment_sha256, ...
    'five_level_upper_s',upper,'remaining_after_conservative_plan_s',remaining, ...
    'prior_charged_core_s',gate.charged_prior_seconds, ...
    'solver_sha256',solver_hash,'surface_sha256',closeout_s2b_hash(surface), ...
    'incident_sha256',closeout_s2b_hash(incident),'source_sha256',closeout_s2b_hash(src_path), ...
    'configuration_sha256',authorization.frozen_configuration_sha256, ...
    'timestamp',char(datetime('now','TimeZone','Asia/Shanghai')));
closeout_s2b_json(fullfile(run_dir,'jobs',['bie_' job_id '_preflight.json']),preflight);
cfg=struct('frequency_hz',4000,'c0_mps',1500,'half_width_m',L(levels), ...
    'points_per_wavelength',ppw(levels),'panel_order',8,'self_quadrature_order',qa(levels), ...
    'halfplane_h_m',-2,'window_plateau_ratio',.85, ...
    'eta_fn',@(xx) closeout_s2c_b_surface(xx,.05,.2,0), ...
    'eta_prime_fn',@(xx) local_slope(xx), ...
    'incident_fn',@(xx,zz) fk_single_receiver_incident(src,xx,zz), ...
    'receiver_x_m',x,'receiver_z_m',z);
cfg_numeric=rmfield(cfg,{'eta_fn','eta_prime_fn','incident_fn'});
assert(nodes(levels)==8*ceil((2*L(levels))/(1500/4000/ppw(levels)*8)));
clock=tic;
try
    boundary=solve_helmholtz_bie_halfplane_vertical(cfg);
catch ME
    seconds=toc(clock);
    failed=struct('job_id',job_id,'status','FAILED','wall_seconds',seconds, ...
        'conservative_charge_seconds',max(150,seconds),'message',ME.message, ...
        'available_physical_bytes',available_physical_bytes,'timestamp',char(datetime('now','TimeZone','Asia/Shanghai')));
    closeout_s2b_json(fullfile(run_dir,'jobs',['bie_' job_id '_failure.json']),failed);
    rethrow(ME);
end
seconds=toc(clock);
assert(seconds<900,'Job exceeded single-job cap; stop all subsequent jobs');
assert(boundary.finite && boundary.node_count==nodes(levels));
receipt=struct('job_id',job_id,'status','COMPLETE','case_id','K_HIGH', ...
    'node_count',boundary.node_count,'self_quadrature_order',qa(levels), ...
    'half_width_m',L(levels),'points_per_wavelength',ppw(levels), ...
    'wall_seconds',seconds,'receiver_values',numel(boundary.receiver_field), ...
    'offgrid_boundary_residual',boundary.offgrid_boundary_residual, ...
    'linear_residual',boundary.linear_residual, ...
    'available_physical_bytes',available_physical_bytes, ...
    'estimated_increment_bytes',peak(levels),'reserve_fraction',.3, ...
    'authorized_subbudget_s',gate.subbudget_s,'budget_amendment_sha256',gate.budget_amendment_sha256, ...
    'solver_sha256',solver_hash,'source_sha256',closeout_s2b_hash(src_path), ...
    'configuration_sha256',authorization.frozen_configuration_sha256, ...
    'field_evaluation_included_in_wall_seconds',true, ...
    'cumulative_core_charged_seconds',gate.charged_prior_seconds+sum(actual)+seconds, ...
    'timestamp',char(datetime('now','TimeZone','Asia/Shanghai')));
save(files{levels},'boundary','src','cfg_numeric','receipt','-v7.3');
closeout_s2b_json(fullfile(run_dir,'jobs',['bie_' job_id '.json']),receipt);
disp(receipt);
end

function d1=local_slope(x)
[~,d1]=closeout_s2c_b_surface(x,.05,.2,0);
end
