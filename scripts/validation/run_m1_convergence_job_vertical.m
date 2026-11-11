function receipt=run_m1_convergence_job_vertical(job_id,run_dir)
%RUN_M1_CONVERGENCE_JOB_VERTICAL One serial M1 job in the attached MATLAB.
% No dispatch loop: caller inspects each durable checkpoint before next job.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
assert(java.io.File(run_dir).isAbsolute(),'Absolute result path required.');
canonical=char(java.io.File(run_dir).getCanonicalPath());
assert(startsWith(canonical,[fullfile(root,'results','validation','m1_convergence') filesep]));
setup_vertical_project();
if ~isfolder(run_dir), mkdir(run_dir); end
budget_path=fullfile(run_dir,'budget.json');
if isfile(budget_path), budget=jsondecode(fileread(budget_path));
else, budget=struct('used_s',15,'limit_s',14400,'initial_control_allowance_s',15); end
assert(budget.used_s<budget.limit_s,'M1 cumulative budget exhausted.');
job_dir=fullfile(run_dir,job_id);
assert(~isfolder(job_dir),'Never overwrite or restart an existing job.');
mkdir(job_dir); t=tic;
old_threads=maxNumCompThreads(1);
cleanup_threads=onCleanup(@()maxNumCompThreads(old_threads)); %#ok<NASGU>
diary(fullfile(job_dir,'stdout.log'));
cleanup_diary=onCleanup(@()diary('off')); %#ok<NASGU>
cfg=local_job(job_id);
receipt=struct('run_id',string(java.io.File(run_dir).getName()), ...
    'job_id',job_id,'status','RUNNING','matlab_pid',feature('getpid'), ...
    'matlab_version',version,'threads',maxNumCompThreads,'config',cfg, ...
    'elapsed_s',0,'exit_code',NaN,'exit_code_kind','job exception status; shared MATLAB process is not exited', ...
    'stdout_log',fullfile(job_dir,'stdout.log'),'stderr_log',fullfile(job_dir,'stderr.log'), ...
    'checkpoint_path',fullfile(job_dir,'checkpoint.json'),'failed_checks',NaN);
receipt.source_sha256=local_sources(root);
receipt.available_before_bytes=0;
receipt.min_available_bytes=Inf;
receipt.max_matlab_memory_bytes=0;
receipt.last_observation=struct();
receipt.estimated_extra_peak_bytes=cfg.estimated_extra_peak_bytes;
trace=struct('direct',[],'incident',[],'reflect',[]);
roi=struct();
payload=struct();
local_json(fullfile(job_dir,'state.json'),receipt);
try
    assert(contains(version,'R2025b') && contains(version,'Update 6'),'Expected R2025b Update 6.');
    [u,ram]=memory;
    receipt.available_before_bytes=ram.PhysicalMemory.Available;
    receipt.min_available_bytes=receipt.available_before_bytes;
    receipt.max_matlab_memory_bytes=u.MemUsedMATLAB;
    disk=java.io.File(run_dir).getUsableSpace(); receipt.disk_before_bytes=double(disk);
    assert(disk>=2*2^30,'Disk gate: less than 2 GiB available.');
    assert(cfg.estimated_extra_peak_bytes<=.7*receipt.available_before_bytes,'RAM reserve preflight failed.');
    if cfg.n>=1966, assert(receipt.available_before_bytes>=1.5*2^30,'Large AS requires at least 1.5 GiB available.'); end
    if startsWith(job_id,'SP')
        gate=fullfile(run_dir,'stage1_audit','checkpoint.json');
        assert(isfile(gate),'Stage 1 independent acceptance checkpoint missing.');
        g=jsondecode(fileread(gate)); assert(strcmp(g.status,'STAGE1_QUALIFIED'),'Stage 1 not qualified.');
    end
    local_json(fullfile(job_dir,'state.json'),receipt);
    x=(-cfg.n/2:cfg.n/2-1)*cfg.dx_m; [X,Y]=meshgrid(x,x);
    center=hypot(X,Y)<=2;
    if cfg.grid_stride==2
        common=false(cfg.n); common(1:2:end,1:2:end)=true; center=center & common;
    end
    roi_idx=find(center); coords=[X(roi_idx),Y(roi_idx)];
    edge=abs(X)>=.45*cfg.width_m | abs(Y)>=.45*cfg.width_m;
    boundary=false(cfg.n); boundary([1 end],:)=true; boundary(:,[1 end])=true;
    clear X Y center
    if strcmp(cfg.kind,'AS')
        source=exp(-((x(:).^2)+(x.^2))/(2*.3^2));
        fields=complex(zeros(numel(roi_idx),3));
        edge_rows=zeros(3,4);
        for j=1:3
            L=[97 100 103];
            [field,meta]=exact_angular_spectrum_one_step_vertical(source,cfg.width_m,cfg.width_m,cfg.f_hz,1500,L(j));
            fields(:,j)=field(roi_idx);
            [ef,bd]=edge_stats(field);
            edge_rows(j,:)=[L(j),ef,bd,max(abs(field(:)))];
            observe_resources();
            local_json(fullfile(job_dir,'state.json'),receipt);
            clear field
        end
        payload=struct('config',cfg,'coords_m',coords,'native_direct',fields(:,1), ...
            'native_incident',fields(:,2),'native_reflect',-fields(:,3), ...
            'edge_table',array2table(edge_rows,'VariableNames',{'distance_m','edge5_fraction','boundary_db','peak'}), ...
            'as_meta',meta,'source_expression','exp(-(x.^2+y.^2)/(2*.3^2))', ...
            'source_peak',max(source(:)),'grid_x_m',x);
        clear source fields
        writetable(payload.edge_table,fullfile(job_dir,'edge_metrics.csv'));
        boundary_pass=all(edge_rows(:,2)<1e-4 & edge_rows(:,3)<=-40);
    elseif strcmp(cfg.kind,'CONTROL')
        p=local_params(cfg); p.save_mode='rx_only';
        a=vertical_channel_model(p);
        p.validation_step_observer=@observe;
        b=vertical_channel_model(p);
        required={'h_direct','h_reflect','h_total','H_f','f_axis','idx_f_ref'};
        dif=0;
        for j=1:numel(required), dif=max(dif,max(abs(a.(required{j})(:)-b.(required{j})(:)))); end
        assert(dif==0,'Observer changed channel outputs.');
        p.enable_surface_reflection=false; p.validation_step_observer=[];
        c=vertical_channel_model(p);
        assert(abs(c.H_f-a.H_direct_f)<1e-12 && c.h_reflect==0);
        payload=struct('disabled_enabled_max_difference',dif,'direct_only_closure',abs(c.H_f-a.H_direct_f), ...
            'note','Reduced-grid implementation control; not a qualification case. No broadband scan or new communication experiment.');
        boundary_pass=true;
    else
        p=local_params(cfg); p.validation_step_observer=@observe;
        out=vertical_channel_model(p);
        p.validation_step_observer=[]; out.config.validation_step_observer=[];
        assert(strcmp(out.phase_reference_meta.phase_contract_version,'2.0.0'));
        payload=struct('config',cfg,'params',p,'output',out,'coords_m',coords, ...
            'native_direct',roi.direct,'native_incident',roi.incident,'native_reflect',roi.reflect, ...
            'trace',trace);
        names={'direct','incident','reflect'}; all_rows=[];
        for j=1:3
            rows=trace.(names{j}); expected=ceil([97 100 3]/(cfg.stepz_lamb*1500/cfg.f_hz));
            assert(size(rows,1)==expected(j)+1 && isequal(rows(:,1),(0:expected(j)).'),'Missing propagation steps.');
            tab=array2table(rows,'VariableNames',{'step','z_m','edge5_fraction','boundary_db','peak','total_power'});
            writetable(tab,fullfile(job_dir,[names{j} '_every_step.csv']));
            all_rows=[all_rows;rows]; %#ok<AGROW>
        end
        boundary_pass=all(all_rows(:,3)<1e-4 & all_rows(:,4)<=-40);
        c=[out.H_f-out.H_direct_f-out.H_reflect_f;out.h_total-out.H_f(out.idx_f_ref); ...
            out.h_direct-out.H_direct_f(out.idx_f_ref);out.h_reflect-out.H_reflect_f(out.idx_f_ref)];
        payload.interface_closure=max(abs(c)); assert(payload.interface_closure<=1e-12);
    end
    if isfield(payload,'native_direct')
        k=2*pi*cfg.f_hz/1500;
        payload.native_total=payload.native_direct+payload.native_reflect;
        payload.physical_direct=payload.native_direct*exp(1i*k*97);
        payload.physical_reflect=payload.native_reflect*exp(1i*k*103);
        payload.physical_total=payload.physical_direct+payload.physical_reflect;
        payload.dsp_direct=conj(payload.native_direct);
        payload.dsp_reflect=conj(payload.native_reflect)*exp(-1i*k*6);
        payload.dsp_total=payload.dsp_direct+payload.dsp_reflect;
        payload.phase_contract_version='2.0.0';
        payload.phase_contract_relative_error=max(abs(payload.dsp_total-conj(payload.physical_total)*exp(1i*k*97)))/max(abs(payload.native_direct));
        assert(payload.phase_contract_relative_error<=1e-12,'Phase v2 relation failed.');
    end
    payload.boundary_pass=boundary_pass;
    save(fullfile(job_dir,'payload.mat'),'payload','-v7.3');
    if isfield(payload,'native_direct'), local_field_csv(fullfile(job_dir,'fields_ieee754.csv'),payload); end
    observe_resources();
    receipt.status='COMPLETE'; receipt.failed_checks=double(~boundary_pass);
    if ~boundary_pass, receipt.status='FAILED_BOUNDARY_GATE'; end
    receipt.exit_code=double(~boundary_pass);
    receipt.boundary_pass=boundary_pass;
    fid=fopen(receipt.stderr_log,'w'); fclose(fid);
catch exception
    receipt.status='FAILED'; receipt.exit_code=1; receipt.failed_checks=1;
    receipt.error_identifier=exception.identifier; receipt.error_message=exception.message;
    fid=fopen(receipt.stderr_log,'w'); fprintf(fid,'%s\n',getReport(exception,'extended','hyperlinks','off')); fclose(fid);
    save(fullfile(job_dir,'partial.mat'),'trace','roi','cfg','-v7.3');
    fprintf('M1_FAILURE %s: %s\n',exception.identifier,exception.message);
end
receipt.elapsed_s=toc(t); budget.used_s=budget.used_s+receipt.elapsed_s;
receipt.cumulative_s=budget.used_s;
receipt.pilot_stop=ismember(job_id,{'P4','SP4'}) && receipt.elapsed_s>900;
receipt.next_allowed_step='Inspect job/checkpoint and independent audit before dispatch; no M2.';
if receipt.elapsed_s>3600 || budget.used_s>14400, receipt.status='FAILED_TIME_GATE'; receipt.failed_checks=receipt.failed_checks+1; end
local_json(budget_path,budget);
local_json(fullfile(job_dir,'state.json'),receipt);
clear cleanup_diary
files=dir(job_dir); artifacts=struct('path',{},'sha256',{});
for j=1:numel(files)
    if files(j).isdir || ismember(files(j).name,{'checkpoint.json'}), continue; end
    artifacts(end+1)=struct('path',files(j).name,'sha256',local_sha(fullfile(job_dir,files(j).name))); %#ok<AGROW>
end
receipt.artifacts_sha256=artifacts;
local_json(receipt.checkpoint_path,receipt);
fprintf('M1_RECEIPT %s\n',jsonencode(rmfield(receipt,{'source_sha256','artifacts_sha256'})));

    function [ef,bd,peak,power]=edge_stats(field)
        pwr=abs(field).^2; power=sum(pwr(:)); peak=max(abs(field(:)));
        ef=sum(pwr(edge))/power;
        bd=20*log10(max(abs(field(boundary)))/peak);
    end
    function observe(field,state)
        [ef,bd,peak,power]=edge_stats(field);
        name=state.path;
        trace.(name)(state.step+1,:)=[state.step,state.z_m,ef,bd,peak,power];
        if state.step==state.n_steps, roi.(name)=field(roi_idx); end
        if mod(state.step,16)==0 || state.step==state.n_steps
            receipt.last_observation=state; observe_resources();
            receipt.elapsed_s=toc(t); local_json(fullfile(job_dir,'state.json'),receipt);
        end
    end
    function observe_resources()
        [um,sm]=memory;
        receipt.min_available_bytes=min(receipt.min_available_bytes,sm.PhysicalMemory.Available);
        receipt.max_matlab_memory_bytes=max(receipt.max_matlab_memory_bytes,um.MemUsedMATLAB);
        assert(sm.PhysicalMemory.Available>=.3*receipt.available_before_bytes,'RAM reserve fell below 30 percent of initial available RAM.');
        assert(toc(t)<=3600,'Single-job 3600 s time gate.');
        assert(budget.used_s+toc(t)<=14400,'Cumulative 14400 s time gate.');
    end
end

function cfg=local_job(id)
cfg=struct('id',id,'kind','PE','f_hz',4000,'n',984,'width_m',192.1875, ...
    'dx_m',.1953125,'grid_stride',1,'stepz_lamb',.5,'sponge_ratio',0, ...
    'alpha_max_np_per_m',0,'estimated_extra_peak_bytes',512*2^20);
if strcmp(id,'CONTROL'), cfg.kind='CONTROL'; cfg.n=64; cfg.width_m=12.5; cfg.estimated_extra_peak_bytes=64*2^20; return; end
if startsWith(id,'A-')
    cfg.kind='AS'; group=id(3); cfg.f_hz=str2double(id(4:end))*1000;
    switch group
        case 'C'
        case 'W', cfg.n=1310; cfg.width_m=255.859375;
        case 'G', cfg.n=1968; cfg.dx_m=.09765625; cfg.grid_stride=2; cfg.estimated_extra_peak_bytes=2^30;
        case 'R', assert(cfg.f_hz==4000); cfg.n=1966; cfg.width_m=383.984375; cfg.estimated_extra_peak_bytes=2^30;
        otherwise, error('Unknown AS job.');
    end
else
    cfg.f_hz=str2double(id(end))*1000;
    assert(ismember(id,{'P4','P6','P8','Z4','Z8','SP4','SP6','SP8'}));
    if startsWith(id,'Z'), cfg.stepz_lamb=.25; end
    if startsWith(id,'SP'), cfg.sponge_ratio=.12; cfg.alpha_max_np_per_m=.15; end
end
assert(ismember(cfg.f_hz,[4000,6000,8000]) && cfg.n*cfg.dx_m==cfg.width_m);
end
function p=local_params(c)
p=struct('f0',c.f_hz,'enable_wideband',false,'c0',1500,'z_max',100, ...
    'z_tx',100,'z_rx',3,'x_tx',0,'y_tx',0,'x_rx',0,'y_rx',0, ...
    'nx',c.n,'ny',c.n,'xw',c.width_m,'yw',c.width_m,'sigma_src_m',.3, ...
    'stepz_lamb',c.stepz_lamb,'show_figures',false,'enforce_1_over_R',false, ...
    'save_mode','slice','nout',1,'use_gpu',false,'enable_bubbles',false, ...
    'env_mode','uniform','enable_surface_reflection',true, ...
    'surface_boundary_model','kirchhoff_spatial','surface_reflect_coeff',-1, ...
    'surface_phase_mode','normal','sea_hs_target',0,'sea_seed',910100, ...
    'surface_elevation_override_xy',zeros(c.n),'alpha_max_np_per_m',c.alpha_max_np_per_m, ...
    'sponge_ratio',c.sponge_ratio,'validation_allow_extended_sponge_ratio',true, ...
    'validation_allow_extended_window',true,'channel_phase_reference','direct_dsp');
end
function local_field_csv(file,s)
fid=fopen(file,'w'); cleaner=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'x_hex,y_hex,d_real_hex,d_imag_hex,r_real_hex,r_imag_hex\n');
v=[s.coords_m,real(s.native_direct),imag(s.native_direct),real(s.native_reflect),imag(s.native_reflect)];
for j=1:size(v,1)
    h=num2hex(v(j,:).'); fprintf(fid,'%s,%s,%s,%s,%s,%s\n',h(1,:),h(2,:),h(3,:),h(4,:),h(5,:),h(6,:));
end
end
function source=local_sources(root)
old=jsondecode(fileread(fullfile(root,'results','validation','m0_phase_contract','20261110_m0_native_dsp_v2_01','code_versions_final.json')));
paths={old(endsWith({old.path},'.m')).path};
paths=[paths,{'scripts/validation/run_m1_convergence_job_vertical.m','scripts/validation/support/exact_angular_spectrum_one_step_vertical.m'}];
source=struct('path',{},'sha256',{});
for j=1:numel(paths), source(end+1)=struct('path',paths{j},'sha256',local_sha(fullfile(root,paths{j}))); end %#ok<AGROW>
end
function value=local_sha(file)
md=java.security.MessageDigest.getInstance('SHA-256'); fid=fopen(file,'rb'); cleaner=onCleanup(@()fclose(fid)); %#ok<NASGU>
while ~feof(fid), bytes=fread(fid,2^20,'*uint8'); if ~isempty(bytes), md.update(typecast(bytes,'int8')); end, end
value=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
function local_json(file,s)
temp=[file '.writing']; fid=fopen(temp,'w','n','UTF-8'); assert(fid>0);
fprintf(fid,'%s',jsonencode(s,PrettyPrint=true)); fclose(fid); movefile(temp,file,'f');
end
