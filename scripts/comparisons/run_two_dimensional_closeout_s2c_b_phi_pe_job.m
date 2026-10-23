function receipt = run_two_dimensional_closeout_s2c_b_phi_pe_job(run_dir,available_physical_bytes)
%RUN_TWO_DIMENSIONAL_CLOSEOUT_S2C_B_PHI_PE_JOB Frozen-source PHI_PLUS90 PE control.
arguments
    run_dir (1,:) char
    available_physical_bytes (1,1) double {mustBePositive}
end
assert(isfolder(run_dir) && contains(run_dir,'20261023_s2c_b_phi_plus90_01'));
file=fullfile(run_dir,'jobs','pe_phi.mat'); assert(~isfile(file));
ids={'main96','self64','space14','window60','window70'};
root=fileparts(fileparts(fileparts(fileparts(run_dir))));
prior=fullfile(root,'results/comparisons/two_dimensional_channel_closeout/20261023_s2c_a_1642');
formal=fullfile(root,'results/comparisons/two_dimensional_channel_closeout/20261023_s2b_1305/formal_reference_v2');
gate=closeout_s2c_b_phi_gate(run_dir);
bie_seconds=0;
for j=1:5
    p=fullfile(run_dir,'jobs',['bie_' ids{j} '.mat']); assert(isfile(p),'Complete five-level BIE first');
    q=load(p,'receipt'); assert(strcmp(q.receipt.status,'COMPLETE')); bie_seconds=bie_seconds+q.receipt.wall_seconds;
end
assert(bie_seconds+300<=gate.subbudget_s && gate.core_remaining_s-bie_seconds-300>=gate.S3A_reserve_s);
estimated_increment=128*1024^2;
assert(estimated_increment<=.7*available_physical_bytes,'Current PE physical RAM gate failed');
src_file=fullfile(prior,'protocol','frozen_source.mat'); s=load(src_file,'src'); src=s.src;
t=load(fullfile(formal,'reference_data.mat'),'x','z','U'); x=t.x(:); z=t.z(:);
assert(numel(x)==417 && numel(src.coeff)==984);
W=192.1875; x0=(-W/2)+(0:983)'*(W/984); source=ifft(src.coeff(:));
assert(max(abs(source-exp(-.5*(x0/.3).^2)))<1e-12);
eta=closeout_s2c_b_surface(x0,.05,.1,pi/2);
cfg=struct('frequency_hz',4000,'c0_mps',1500,'xw_m',W,'nx',984, ...
    'z_tx_m',100,'z_rx_m',3,'sigma_src_m',.3,'surface_reflect_coeff',-1, ...
    'step_m',.05,'x_rx_m',0,'reflection_model','model0_normal', ...
    'source_field_override',source,'eta_z_m',eta);
preflight=struct('job_id','pe_phi','status','PREFLIGHT_PASS', ...
    'available_physical_bytes',available_physical_bytes, ...
    'estimated_increment_bytes',estimated_increment,'reserve_fraction',.3, ...
    'authorized_subbudget_s',gate.subbudget_s,'budget_amendment_sha256',gate.budget_amendment_sha256, ...
    'source_sha256',closeout_s2b_hash(src_file), ...
    'adapter_sha256',closeout_s2b_hash(fullfile(root,'scripts/validation/support/run_pe_physical_surface_model0_validation.m')), ...
    'surface_sha256',closeout_s2b_hash(fullfile(root,'scripts/validation/support/closeout_s2c_b_surface.m')), ...
    'timestamp',char(datetime('now','TimeZone','Asia/Shanghai')));
closeout_s2b_json(fullfile(run_dir,'jobs','pe_phi_preflight.json'),preflight);
clock=tic;
try
    flat_cfg=cfg; flat_cfg.eta_z_m=zeros(size(eta));
    flat_run=run_pe_physical_surface_model0_validation(flat_cfg);
    flat_native=flat_run.reflected_native(:);
    [ok,ixflat]=ismember(x(1:411),x0); assert(all(ok));
    flat_relative=norm(flat_native(ixflat)-t.U(1:411,1,2))/norm(t.U(1:411,1,2));
    assert(flat_relative<1e-10,'Flat same-source PE/BIE control failed');
    u=complex(zeros(417,1)); half=u;
    depths=unique(z).';
    for zz=depths
        cfg.z_rx_m=zz; a=run_pe_physical_surface_model0_validation(cfg);
        rows=z==zz; [ok,ii]=ismember(x(rows),x0); assert(all(ok));
        v=a.reflected_native(:); u(rows)=v(ii);
        cfg_half=cfg; cfg_half.step_m=.025;
        h=run_pe_physical_surface_model0_validation(cfg_half);
        vh=h.reflected_native(:); half(rows)=vh(ii);
    end
    probes=readtable(fullfile(prior,'protocol','PE_frozen_probes.csv'),TextType='string');
    assert(height(probes)==13);
    coeff=complex(zeros(1968,1)); coeff(1:492)=src.coeff(1:492)*2;
    coeff(end-491:end)=src.coeff(493:end)*2;
    fine_source=ifft(coeff); assert(max(abs(fine_source(1:2:end)-source))<1e-12);
    xf=(-W/2)+(0:1967)'*(W/1968);
    fine_cfg=cfg; fine_cfg.nx=1968; fine_cfg.source_field_override=fine_source;
    fine_cfg.eta_z_m=closeout_s2c_b_surface(xf,.05,.1,pi/2);
    fine_cfg.step_m=.05;
    doubled=complex(zeros(13,1));
    for zz=unique(probes.z_m).'
        fine_cfg.z_rx_m=zz; a=run_pe_physical_surface_model0_validation(fine_cfg);
        rows=probes.z_m==zz; [ok,ii]=ismember(probes.x_m(rows),xf);
        assert(all(ok)); v=a.reflected_native(:); doubled(rows)=v(ii);
    end
    [ok,ip]=ismember([probes.x_m,probes.z_m],[x,z],'rows'); assert(all(ok));
    primary_probes=u(ip); half_probes=half(ip);
catch ME
    seconds=toc(clock);
    failed=struct('job_id','pe_phi','status','FAILED','wall_seconds',seconds, ...
        'conservative_charge_seconds',max(150,seconds),'message',ME.message);
    closeout_s2b_json(fullfile(run_dir,'jobs','pe_phi_failure.json'),failed);
    rethrow(ME);
end
seconds=toc(clock); assert(seconds<900);
assert(all(isfinite([u;half;doubled])));
payload=struct('x_m',x,'z_m',z,'u_native',u,'half_step_native',half, ...
    'double_grid_probe_native',doubled,'probe_ids',probes.id, ...
    'probe_x_m',probes.x_m,'probe_z_m',probes.z_m, ...
    'primary_probe_native',primary_probes,'half_step_probe_native',half_probes, ...
    'flat_z3_relative_to_BIE',flat_relative, ...
    'double_grid_source','zero-padded SAME frozen N984 complex spectrum; no new Gaussian source', ...
    'surface_version','C2 PHI_PLUS90 A0.05 K0.1 phi+pi/2 physical depth', ...
    'representation','native exp(-i omega t); reflected unit-source response');
receipt=struct('job_id','pe_phi','status','COMPLETE','wall_seconds',seconds, ...
    'available_physical_bytes',available_physical_bytes, ...
    'estimated_increment_bytes',estimated_increment,'reserve_fraction',.3, ...
    'authorized_subbudget_s',gate.subbudget_s,'budget_amendment_sha256',gate.budget_amendment_sha256, ...
    'flat_control_relative',flat_relative,'primary_points',417, ...
    'half_step_points',417,'double_grid_probe_points',13, ...
    'bie_seconds_prior',bie_seconds,'cumulative_core_charged_seconds', ...
    gate.charged_prior_seconds+bie_seconds+seconds, ...
    'timestamp',char(datetime('now','TimeZone','Asia/Shanghai')));
save(file,'payload','receipt','-v7.3');
closeout_s2b_json(fullfile(run_dir,'jobs','pe_phi.json'),receipt); disp(receipt);
end
