function receipt = run_two_dimensional_closeout_s2c_b_fullfk_job(run_dir,available_physical_bytes)
%RUN_TWO_DIMENSIONAL_CLOSEOUT_S2C_B_FULLFK_JOB Complete K_HIGH GL8 reference.
arguments
    run_dir (1,:) char
    available_physical_bytes (1,1) double {mustBePositive}
end
assert(isfolder(run_dir) && contains(run_dir,'20261023_s2c_b_khigh_'));
file=fullfile(run_dir,'jobs','fullfk_khigh.mat'); assert(~isfile(file));
root=fileparts(fileparts(fileparts(fileparts(run_dir))));
prior=fullfile(root,'results/comparisons/two_dimensional_channel_closeout/20261023_s2c_a_1642');
formal=fullfile(root,'results/comparisons/two_dimensional_channel_closeout/20261023_s2b_1305/formal_reference_v2');
gate=closeout_s2c_b_effective_gate(run_dir);
ids={'main96','self64','space14','window60','window70'};
bie_seconds=0;
for j=1:5
    p=fullfile(run_dir,'jobs',['bie_' ids{j} '.mat']); assert(isfile(p));
    q=load(p,'receipt'); assert(strcmp(q.receipt.status,'COMPLETE')); bie_seconds=bie_seconds+q.receipt.wall_seconds;
end
pe_file=fullfile(run_dir,'jobs','pe_khigh.mat'); assert(isfile(pe_file));
q=load(pe_file,'receipt'); assert(strcmp(q.receipt.status,'COMPLETE')); pe_seconds=q.receipt.wall_seconds;
assert(bie_seconds+pe_seconds+200<=gate.subbudget_s && ...
    gate.core_remaining_s-bie_seconds-pe_seconds-200>=gate.S3A_reserve_s);
estimated_increment=256*1024^2;
assert(estimated_increment<=.7*available_physical_bytes,'Current FullFK physical RAM gate failed');
src_file=fullfile(prior,'protocol','frozen_source.mat'); s=load(src_file,'src');
t=load(fullfile(formal,'reference_data.mat'),'x','z','flat_exact');
x=t.x(:); z=t.z(:); flat_exact=t.flat_exact(:);
assert(numel(x)==417 && numel(z)==417);
helper=fullfile(root,'scripts/validation/support/closeout_s2c_b_fullfk_analytic_exterior.m');
preflight=struct('job_id','fullfk_khigh','status','PREFLIGHT_PASS', ...
    'available_physical_bytes',available_physical_bytes, ...
    'estimated_increment_bytes',estimated_increment,'reserve_fraction',.3, ...
    'authorized_subbudget_s',gate.subbudget_s,'budget_amendment_sha256',gate.budget_amendment_sha256, ...
    'prior_BIE_seconds',bie_seconds,'prior_PE_seconds',pe_seconds, ...
    'source_sha256',closeout_s2b_hash(src_file),'helper_sha256',closeout_s2b_hash(helper), ...
    'timestamp',char(datetime('now','TimeZone','Asia/Shanghai')));
closeout_s2b_json(fullfile(run_dir,'jobs','fullfk_khigh_preflight.json'),preflight);
clock=tic; wall=zeros(4,1);
try
    tt=tic; low=closeout_s2c_b_fullfk_analytic_exterior(s.src,x,z,400,50,flat_exact); wall(1)=toc(tt);
    tt=tic; high=closeout_s2c_b_fullfk_analytic_exterior(s.src,x,z,800,50,flat_exact); wall(2)=toc(tt);
    tt=tic; extent60=closeout_s2c_b_fullfk_analytic_exterior(s.src,x,z,960,60,flat_exact); wall(3)=toc(tt);
    tt=tic; extent90=closeout_s2c_b_fullfk_analytic_exterior(s.src,x,z,1440,90,flat_exact); wall(4)=toc(tt);
catch ME
    seconds=toc(clock);
    failed=struct('job_id','fullfk_khigh','status','FAILED','wall_seconds',seconds, ...
        'conservative_charge_seconds',max(150,seconds),'message',ME.message);
    closeout_s2b_json(fullfile(run_dir,'jobs','fullfk_khigh_failure.json'),failed);
    rethrow(ME);
end
seconds=toc(clock); assert(seconds<900);
assert(all(isfinite(high.receiver_field)) && all(high.local_radius_m==Inf));
component_max=max(abs(high.receiver_field-high.term_incident_dn-high.term_green_dn));
quadrature_relative=norm(high.receiver_field-low.receiver_field)/norm(high.receiver_field);
L60_L90_relative=norm(extent60.receiver_field-extent90.receiver_field)/norm(high.receiver_field);
compact_L60_relative=norm(high.receiver_field-extent60.receiver_field)/norm(high.receiver_field);
payload=struct('u_native',high.receiver_field,'x_m',x,'z_m',z, ...
    'GL8_400',low,'GL8_800',high,'extent60',extent60,'extent90',extent90, ...
    'quadrature_relative',quadrature_relative,'L60_L90_relative',L60_L90_relative, ...
    'compact_L60_relative',compact_L60_relative,'component_max_abs',component_max, ...
    'version','FullFK_complete_analytic_flat_exterior_GL8_v1');
receipt=struct('job_id','fullfk_khigh','status','COMPLETE','wall_seconds',seconds, ...
    'job_wall_seconds',wall,'available_physical_bytes',available_physical_bytes, ...
    'estimated_increment_bytes',estimated_increment,'reserve_fraction',.3, ...
    'authorized_subbudget_s',gate.subbudget_s,'budget_amendment_sha256',gate.budget_amendment_sha256, ...
    'quadrature_relative',quadrature_relative,'L60_L90_relative',L60_L90_relative, ...
    'component_max_abs',component_max,'local_radius_m',Inf, ...
    'cumulative_core_charged_seconds',gate.charged_prior_seconds+bie_seconds+pe_seconds+seconds, ...
    'timestamp',char(datetime('now','TimeZone','Asia/Shanghai')));
save(file,'payload','receipt','-v7.3');
closeout_s2b_json(fullfile(run_dir,'jobs','fullfk_khigh.json'),receipt); disp(receipt);
end
