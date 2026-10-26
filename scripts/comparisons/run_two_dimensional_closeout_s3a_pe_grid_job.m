function receipt = run_two_dimensional_closeout_s3a_pe_grid_job(root,run_dir,available_physical_bytes)
%RUN_TWO_DIMENSIONAL_CLOSEOUT_S3A_PE_GRID_JOB Same frozen spectrum N984->N1968.
arguments
    root (1,:) char
    run_dir (1,:) char
    available_physical_bytes (1,1) double {mustBePositive}
end
assert(isfolder(root) && isfolder(run_dir) && startsWith(run_dir,root));
file=fullfile(run_dir,'jobs','pe_grid.mat'); assert(~isfile(file));
budget=jsondecode(fileread(fullfile(run_dir,'protocol','budget.json')));
assert(ismember(budget.seed,[260001 260002]) && ...
    contains(run_dir,sprintf('s3a_seed%d',budget.seed)));
ids={'main96','self64','space14','window60','window70','pe','fullfk'};
spent=0;
for j=1:7
    a=load(fullfile(run_dir,'jobs',[ids{j} '.mat']),'receipt');
    assert(strcmp(a.receipt.status,'COMPLETE')); spent=spent+a.receipt.wall_seconds;
end
if isfield(budget,'upper_pe_grid_s'), grid_upper=budget.upper_pe_grid_s;
else, grid_upper=30; end
assert(spent+grid_upper+budget.failure_reserve_s<budget.subbudget_s);
estimated_increment=256*1024^2;
assert(estimated_increment<.7*available_physical_bytes);
q=load(fullfile(run_dir,'jobs','pe.mat'),'payload');
t=load(fullfile(run_dir,'input','reference_data.mat'),'x','z');
s=load(fullfile(run_dir,'input','frozen_source.mat'),'src'); src=s.src;
p=readtable(fullfile(run_dir,'input',sprintf('seed_%d_coefficients.csv',budget.seed)));
W=192.1875; n=1968; xf=(-W/2)+(0:n-1)'*(W/n);
coeff=complex(zeros(n,1)); coeff(1:492)=2*src.coeff(1:492);
coeff(end-491:end)=2*src.coeff(493:end);
source=ifft(coeff); original=ifft(src.coeff(:));
assert(max(abs(source(1:2:end)-original))<1e-12);
eta=closeout_s3a_pm_surface(p,xf);
cfg=struct('frequency_hz',4000,'c0_mps',1500,'xw_m',W,'nx',n, ...
    'z_tx_m',100,'z_rx_m',3,'sigma_src_m',.3,'surface_reflect_coeff',-1, ...
    'step_m',.05,'x_rx_m',0,'reflection_model','model0_normal', ...
    'source_field_override',source,'eta_z_m',eta);
pre=struct('job_id','pe_grid','available_physical_bytes',available_physical_bytes, ...
    'estimated_increment_bytes',estimated_increment,'reserve_fraction',.3, ...
    'same_frozen_source_spectrum',true,'same_surface_coefficients',true, ...
    'same_receiver_points',417);
local_json(fullfile(run_dir,'jobs','pe_grid_preflight.json'),pre);
clock=tic;
try
    u=complex(zeros(417,1));
    for zz=unique(t.z(:)).'
        cfg.z_rx_m=zz; a=run_pe_physical_surface_model0_validation(cfg);
        rows=t.z(:)==zz; [ok,ii]=ismember(t.x(rows),xf); assert(all(ok));
        v=a.reflected_native(:); u(rows)=v(ii);
    end
    seconds=toc(clock); assert(seconds<900);
    assert(all(isfinite(u)) && numel(u)==417);
    payload=struct('u_native',u,'x_m',t.x(:),'z_m',t.z(:), ...
        'source_spectrum_N984_to_N1968','zero-padded exact same N984 coefficients', ...
        'version','PE-mode0 S2a-depth-positive-native-v1 grid sensitivity');
    receipt=struct('job_id','pe_grid','status','COMPLETE','wall_seconds',seconds, ...
        'available_physical_bytes',available_physical_bytes, ...
        'estimated_increment_bytes',estimated_increment, ...
        'relative_to_primary',norm(u-q.payload.u_native(:))/norm(q.payload.u_native(:)), ...
        'receiver_values',417);
    save(file,'payload','receipt','-v7.3');
    local_json(fullfile(run_dir,'jobs','pe_grid.json'),receipt);
catch ME
    failure=struct('job_id','pe_grid','status','FAILED','wall_seconds',toc(clock), ...
        'message',ME.message,'identifier',ME.identifier);
    local_json(fullfile(run_dir,'jobs','pe_grid_failure.json'),failure);
    rethrow(ME);
end
disp(receipt);
end

function local_json(path,value)
fid=fopen(path,'w'); assert(fid>0); c=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(value));
end
