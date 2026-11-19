function receipt = run_d1b_scale_pilot(package_root,out,manifest_sha,mode)
%RUN_D1B_SCALE_PILOT Manual-only 1536 control -> 4096 audit -> 8192 pilot.
arguments
    package_root (1,1) string
    out (1,1) string
    manifest_sha (1,1) string = ""
    mode (1,1) string {mustBeMember(mode,["dry-run","execute"])} = "dry-run"
end
cfg=jsondecode(fileread(fullfile(package_root,'config.json')));
assert(strcmp(cfg.scope,'D1B_SCALE_PILOT_ONLY') && ~cfg.full_d1b_enabled);
assert(isequal(cfg.nodes,[4096;8192]) && cfg.control_nodes==1536 && cfg.workers==2 && cfg.threads==1);
assert(cfg.max_wall_s==1800 && cfg.reserve_fraction==.3 && cfg.rss_cap_GiB==12);
data=load(fullfile(package_root,'inputs','profile.mat'),'profile','line_x','line_z');
assert(numel(data.profile.k)==157 && data.profile.plateau==220 && data.profile.support==240);
assert(isequal(data.line_x,(-1139:1139)'*(200/1139)) && isequal(data.line_z,6*ones(2279,1)));
receipt=struct('status','DRY_RUN_NO_POOL_NO_SOLVE','scope',cfg.scope,'nodes',cfg.nodes, ...
    'SSH_connections',0,'new_solves',0,'full_22400_enabled',false);
if mode=="dry-run", return, end
assert(isunix && strcmp(version('-release'),'2025b') && maxNumCompThreads==1);
assert(strcmp(getenv('D1B_SCALE_MANIFEST_SHA'),manifest_sha) && strlength(manifest_sha)==64);
assert(isfolder(out) && ~isfile(fullfile(out,'checkpoint.json')) && isempty(gcp('nocreate')));
t0=tic; receipt=struct('status','RUNNING','config',cfg,'manifest_sha256',manifest_sha, ...
    'matlab_version',version,'main_pid',feature('getpid'),'completed',{{}},'failed_checks',0);
try
    product=ver('parallel'); [profiles,default_profile]=parallel.listProfiles;
    [ok,reason]=license('checkout','Distrib_Computing_Toolbox');
    d1b_scale_json(fullfile(out,'pct.json'),struct('product',product,'profiles',{profiles}, ...
        'default_profile',default_profile,'checkout',ok,'message',reason));
    assert(ok==1,'d1b:License','PCT checkout failed: %s',reason);
    assert(~isempty(product) && any(strcmp(profiles,'Processes')),'d1b:Profile','Processes profile unavailable');
    src=struct('frequency_hz',4000,'c_mps',1500,'rho_kgpm3',1025,'v0_mps',.001, ...
        'total_width_m',.1,'z_source_m',1490);
    local_checkpoint(out,receipt,t0,'SERIAL_CONTROL');
    local_case(1536,0,'control_serial',"",0);
    local_gate('control_serial');
    local_checkpoint(out,receipt,t0,'POOL_START');
    monitor=struct('monitor_dir',out,'case_id','pool'); d1b_scale_phase(monitor,'startup');
    c=parcluster('Processes'); c.NumThreads=1;
    storage=fullfile(out,'pct_jobs_2'); mkdir(storage); c.JobStorageLocation=char(storage);
    pool=parpool(c,2,AttachedFiles={which('d1b_scale_bie'),which('d1b_tiles'), ...
        which('d1b_pm_surface'),which('d1b_boundary_audit'),which('baffled_strip_piston_field_vertical')});
    spmd
        registration=[feature('getpid') maxNumCompThreads];
    end
    ids=[registration{1};registration{2}];
    assert(numel(unique(ids(:,1)))==2 && all(ids(:,2)==1) && contains(class(pool),'ProcessPool'));
    pool_info=struct('worker_pid',ids(:,1),'compute_threads',ids(:,2),'profile',c.Profile,'pool_class',class(pool));
    d1b_scale_json(fullfile(out,'pool.json'),pool_info);
    d1b_scale_json(fullfile(out,'active_pool.json'),pool_info);
    d1b_scale_phase(monitor,'IDLE');
    local_case(1536,2,'control_process2',fullfile(out,'control_serial','data.mat'),ids(:,1));
    local_gate('control_process2');
    receipt.completed{end+1}='CONTROL_EQUIVALENCE_PASS';
    local_checkpoint(out,receipt,t0,'CONTROL_PASS');
    for nodes=cfg.nodes.'
        if nodes==8192, local_gate('N4096'); end
        name=sprintf('N%d',nodes); local_checkpoint(out,receipt,t0,name);
        local_case(nodes,2,name,"",ids(:,1)); local_gate(name);
        receipt.completed{end+1}=[name '_INDEPENDENT_PASS'];
        local_checkpoint(out,receipt,t0,[name '_PASS']);
    end
    monitor=struct('monitor_dir',out,'case_id','closeout'); d1b_scale_phase(monitor,'pool_close');
    delete(pool);
    d1b_scale_json(fullfile(out,'active_pool.json'),struct('worker_pid',[]));
    d1b_scale_phase(monitor,'IDLE');
    receipt.status='PILOT_NUMERICS_PASS_REVIEW_REQUIRED';
    local_checkpoint(out,receipt,t0,'NUMERICS_DONE');
catch err
    receipt.status='FAILED_STOPPED'; receipt.failed_checks=1;
    receipt.error=struct('identifier',err.identifier,'message',err.message);
    local_checkpoint(out,receipt,t0,'FAILED');
    % Keep phase receipt and original error; supervisor owns process cleanup.
    rethrow(err)
end

    function local_case(n,w,id,reference,worker_ids)
        folder=fullfile(out,id); assert(~isfolder(folder)); mkdir(folder);
        L=n*.375/(2*12);
        q=struct('frequency_hz',4000,'c0_mps',1500,'half_width_m',L, ...
            'points_per_wavelength',12,'panel_order',8,'self_quadrature_order',96, ...
            'halfplane_h_m',-2,'window_plateau_ratio',.85, ...
            'eta_fn',@(a) d1b_pm_surface(data.profile,a),'eta_prime_fn',@(a) local_slope(data.profile,a), ...
            'incident_fn',@(a,b) baffled_strip_piston_field_vertical(a,b,src)-baffled_strip_piston_field_vertical(a,-b,src), ...
            'receiver_x_m',data.line_x,'receiver_z_m',data.line_z, ...
            'return_system_for_validation',true,'monitor_dir',out,'case_id',id,'worker_ids',worker_ids);
        start=tic; b=d1b_scale_bie(q,w); assert(b.node_count==n);
        d1b_scale_phase(q,'save_system');
        A=b.validation_system.A; rhs=b.validation_system.rhs; b=rmfield(b,'validation_system');
        save(fullfile(folder,'data.mat'),'b','A','rhs','-v7.3'); clear A rhs b
        audit=audit_d1b_scale_saved(fullfile(folder,'data.mat'),data.profile,src,w,reference,q);
        d1b_scale_phase(q,'seal');
        d1b_scale_json(fullfile(folder,'seal.json'),struct('status','AUDITED_PASS', ...
            'data_sha256',audit.data_sha256,'audit_sha256',d1b_scale_hash(fullfile(folder,'audit.json')), ...
            'independent_sha256',d1b_scale_hash(fullfile(folder,'independent.mat')), ...
            'manifest_sha256',manifest_sha,'nodes',n,'total_case_s',toc(start)));
        d1b_scale_phase(q,'IDLE');
    end
    function local_gate(id)
        folder=fullfile(out,id); seal=jsondecode(fileread(fullfile(folder,'seal.json')));
        audit=jsondecode(fileread(fullfile(folder,'audit.json')));
        assert(strcmp(seal.status,'AUDITED_PASS') && strcmp(audit.status,'PASS_INDEPENDENT_SAVED_AUDIT') && audit.failed_checks==0);
        assert(strcmp(seal.manifest_sha256,manifest_sha) && strcmp(seal.data_sha256,d1b_scale_hash(fullfile(folder,'data.mat'))));
        assert(strcmp(seal.audit_sha256,d1b_scale_hash(fullfile(folder,'audit.json'))));
        assert(strcmp(seal.independent_sha256,d1b_scale_hash(fullfile(folder,'independent.mat'))));
    end
end

function s=local_slope(p,x)
[~,s]=d1b_pm_surface(p,x);
end
function local_checkpoint(out,receipt,t0,stage)
receipt.active_stage=stage; receipt.elapsed_s=toc(t0);
receipt.next_allowed_step='STOP after 8192; no 22400 or extra workers/retries';
d1b_scale_json(fullfile(out,'checkpoint.json'),receipt);
end
