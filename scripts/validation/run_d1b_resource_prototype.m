function run_d1b_resource_prototype(package_root,out,manifest_sha)
%RUN_D1B_RESOURCE_PROTOTYPE Entry called once by the user-started supervisor.
arguments
    package_root (1,1) string
    out (1,1) string
    manifest_sha (1,1) string
end
assert(isunix && strcmp(version('-release'),'2025b') && maxNumCompThreads==1);
assert(strcmp(getenv('D1B_MANIFEST_SHA'),manifest_sha) && strlength(manifest_sha)==64);
assert(isfolder(out) && ~isfile(fullfile(out,'checkpoint.json')));
cfg = jsondecode(fileread(fullfile(package_root,'config.json')));
assert(strcmp(cfg.scope,'D1B_RESOURCE_PROTOTYPE_ONLY') && ~cfg.full_d1b_enabled);
assert(isequal(cfg.prototype_half_widths_m,[6;12;24]) && cfg.max_nodes==1536 && cfg.repeats==2);
assert(isequal(cfg.pool_workers,[0;2;4]));
t0 = tic;
receipt = struct('status','RUNNING','manifest_sha256',manifest_sha,'config',cfg,'failed_checks',0, ...
    'matlab_version',version,'main_pid',feature('getpid'),'completed',{{}}, ...
    'next_allowed_step','STOP after resource prototype; user review required');
local_checkpoint(out,receipt,t0,'SURFACES');
try
    surface_dir = fullfile(out,'surfaces');
    surface = build_d1b_pm_design(package_root,surface_dir);
    assert(strcmp(surface.status,'PASS_SURFACE_ONLY'));
    receipt.completed{end+1} = 'SURFACES';
    local_checkpoint(out,receipt,t0,'PCT_CHECKOUT');
    product = ver('parallel');
    assert(~isempty(product),'D1b:PCTInstall','PCT missing');
    [profiles,default_profile] = parallel.listProfiles;
    assert(any(strcmp(profiles,'Processes')),'D1b:Profile','Local Processes profile missing');
    [ok,reason] = license('checkout','Distrib_Computing_Toolbox');
    local_json(fullfile(out,'pct.json'),struct('product',product,'profiles',{profiles}, ...
        'default_profile',default_profile,'checkout',ok,'message',reason));
    assert(ok==1,'D1b:PCTLicense','PCT checkout failed: %s',reason);
    receipt.completed{end+1} = 'PCT_CHECKOUT';
    data = load(fullfile(surface_dir,cfg.prototype_case,'profile.mat'),'profile','line_x','line_z');
    p = data.profile; x = data.line_x; z = data.line_z;
    assert(numel(x)==2279 && isequal(x([1 end]),[-200;200]));
    src = struct('frequency_hz',4000,'c_mps',1500,'rho_kgpm3',1025,'v0_mps',.001, ...
        'total_width_m',.1,'z_source_m',1490);
    fk_src = struct('frequency_hz',4000,'c0_mps',1500,'incident_field_fn', ...
        @(a,b) baffled_strip_piston_field_vertical(a,b,src));
    flat = -baffled_strip_piston_field_vertical(x,-z,src);
    refs = fullfile(out,'serial_references'); mkdir(refs);
    for workers = cfg.pool_workers.'
        label = sprintf('process_%d',workers);
        if workers==0, label = 'serial'; end
        local_checkpoint(out,receipt,t0,[label '_startup']);
        folder = fullfile(out,label); mkdir(folder);
        [pool,ids,info] = local_pool(out,workers);
        cleanup = onCleanup(@() local_close_pool(pool,out));
        local_json(fullfile(folder,'pool.json'),info);
        for repeat = 1:cfg.repeats
            for half_width = cfg.prototype_half_widths_m.'
                id = sprintf('%s_r%d_L%02d',label,repeat,half_width);
                local_checkpoint(out,receipt,t0,id);
                ref = fullfile(refs,sprintf('L%02d',half_width));
                result = local_bie(src,p,x,z,half_width,workers,pool,ids,ref,workers==0 && repeat==1);
                local_json(fullfile(folder,[id '.json']),result);
                receipt.completed{end+1} = id;
                local_checkpoint(out,receipt,t0,id);
            end
            id = sprintf('%s_r%d_receivers',label,repeat);
            local_checkpoint(out,receipt,t0,id);
            ref = fullfile(refs,'receivers.mat');
            result = local_receivers(src,fk_src,p,x,z,flat,workers,pool,ids,ref,workers==0 && repeat==1);
            local_json(fullfile(folder,[id '.json']),result);
            receipt.completed{end+1} = id;
            local_checkpoint(out,receipt,t0,id);
        end
        clear cleanup
    end
    receipt.status = 'PROTOTYPE_NUMERICS_PASS';
    local_checkpoint(out,receipt,t0,'DONE');
catch err
    receipt.status = 'FAILED_STOPPED';
    receipt.failed_checks = 1;
    receipt.error = struct('identifier',err.identifier,'message',err.message);
    local_checkpoint(out,receipt,t0,'FAILED');
    rethrow(err)
end
end

function [pool,ids,info] = local_pool(out,workers)
pool = []; ids = 0;
info = struct('workers',workers,'startup_s',0,'worker_pid',[],'compute_threads',[],'profile','serial');
if workers==0, return, end
assert(isempty(gcp('nocreate')));
cluster = parcluster('Processes');
storage = fullfile(out,sprintf('pct_jobs_%d',workers)); mkdir(storage);
cluster.JobStorageLocation = char(storage);
cluster.NumThreads = 1;
t0 = tic;
attached = {which('d1b_tiles'),which('d1b_bie_prototype'),which('d1b_boundary_audit'), ...
    which('d1b_pm_surface'),which('d1b_ps_prototype'),which('d1b_fullfk_prototype'), ...
    which('run_d1b_resource_prototype'),which('baffled_strip_piston_field_vertical')};
pool = parpool(cluster,workers,AttachedFiles=attached);
try
    traffic_token = ticBytes(pool);
    spmd
        register = [feature('getpid') maxNumCompThreads];
    end
    rows = zeros(workers,2);
    for w = 1:workers, rows(w,:) = register{w}; end
    assert(size(rows,1)==workers && numel(unique(rows(:,1)))==workers && all(rows(:,2)==1));
    ids = rows(:,1);
    startup_traffic = local_bytes_end(pool,traffic_token);
    info = struct('workers',workers,'startup_s',toc(t0),'worker_pid',ids, ...
        'compute_threads',rows(:,2),'profile',cluster.Profile,'pool_class',class(pool), ...
        'startup_bytes_by_worker',startup_traffic);
    local_json(fullfile(out,'active_pool.json'),info);
catch err
    delete(pool);
    rethrow(err)
end
end

function local_close_pool(pool,out)
if isempty(pool), return, end
local_json(fullfile(out,'active_pool.json'),struct('worker_pid',[]));
delete(pool);
end

function result = local_bie(src,p,x,z,L,workers,pool,ids,ref,make_ref)
cfg = struct('frequency_hz',4000,'c0_mps',1500,'half_width_m',L, ...
    'points_per_wavelength',12,'panel_order',8,'self_quadrature_order',96, ...
    'halfplane_h_m',-2,'window_plateau_ratio',.85, ...
    'eta_fn',@(a) d1b_pm_surface(p,a),'eta_prime_fn',@(a) local_slope(p,a), ...
    'incident_fn',@(a,b) baffled_strip_piston_field_vertical(a,b,src)-baffled_strip_piston_field_vertical(a,-b,src), ...
    'receiver_x_m',x,'receiver_z_m',z,'return_system_for_validation',true);
assert(2*L*12/.375<=1536);
token = local_bytes_begin(pool);
t0 = tic; b = d1b_bie_prototype(cfg,workers); bie_wall_s = toc(t0);
traffic_bie = local_bytes_end(pool,token);
A = b.validation_system.A; rhs = b.validation_system.rhs;
b = rmfield(b,'validation_system');
assert(b.finite && all(isfinite(A),'all') && b.linear_residual<1e-10);
token = local_bytes_begin(pool);
t0 = tic; [ind,ind_owner,ind_relative] = d1b_boundary_audit(b,p,src,workers); audit_s = toc(t0);
traffic_audit = local_bytes_end(pool,token);
assert(all(isfinite(ind)) && isfinite(ind_relative) && isfinite(b.offgrid_boundary_residual));
owner = b.worker_rows;
density = b.density; rx = b.receiver_field; boundary = b.boundary_total_residual;
if make_ref
    save(ref+'.mat','A','rhs','density','rx','boundary','ind','-v7.3');
end
f = matfile(ref+'.mat');
% The serial matrix stays on disk. Compare one 32-row slab at a time.
row_error = zeros(size(A,1),1); scale = row_error; peak = row_error; ref_peak = row_error;
for first = 1:32:size(A,1)
    ii = first:min(first+31,size(A,1));
    baseline = f.A(ii,:); diff = A(ii,:)-baseline;
    row_error(ii) = vecnorm(diff,2,2); scale(ii) = vecnorm(baseline,2,2);
    peak(ii) = max(abs(diff),[],2); ref_peak(ii) = max(abs(baseline),[],2);
end
matrix = local_worker_error(ids,owner.assembly,row_error,scale,peak,ref_peak,1e-11);
result = struct('status','PASS_PROTOTYPE_EQUIVALENCE','failed_checks',0,'workers',workers, ...
    'node_count',b.node_count,'half_width_m',L,'timings_s',b.cost_meta.timings, ...
    'bie_wall_s',bie_wall_s,'independent_audit_s',audit_s,'linear_relative',b.linear_residual, ...
    'midpoint_relative',b.offgrid_boundary_residual,'independent_boundary_relative',ind_relative, ...
    'midpoint_scientific_gate_pass',b.offgrid_boundary_residual<1e-4, ...
    'independent_scientific_gate_pass',ind_relative<1e-4, ...
    'BIE_bytes_by_worker',traffic_bie,'audit_bytes_by_worker',traffic_audit,'matrix_checks',matrix, ...
    'rhs_checks',local_compare(0,zeros(size(rhs,1),1),rhs,f.rhs,1e-11), ...
    'density_checks',local_compare(0,zeros(size(rhs,1),1),density,f.density,1e-10), ...
    'receiver_checks',local_compare(ids,owner.receiver,rx,f.rx,1e-11), ...
    'midpoint_checks',local_compare(ids,owner.boundary,boundary,f.boundary,1e-10), ...
    'independent_checks',local_compare(ids,ind_owner,ind,f.ind,1e-10));
end

function result = local_receivers(src,fk_src,p,x,z,flat,workers,pool,ids,ref,make_ref)
token = local_bytes_begin(pool);
[ps,po,ps_setup,ps_rx] = d1b_ps_prototype(src,p,x,z,flat,workers);
traffic_ps = local_bytes_end(pool,token);
token = local_bytes_begin(pool);
[fk,fo,fk_setup,fk_rx,terms] = d1b_fullfk_prototype(fk_src,p,x,z,480,240,flat,workers);
traffic_fk = local_bytes_end(pool,token);
if make_ref, save(ref,'ps','fk','terms','-v7'); end
ref = load(ref,'ps','fk','terms');
result = struct('status','PASS_PROTOTYPE_EQUIVALENCE','failed_checks',0,'workers',workers, ...
    'receiver_count',numel(x),'PS_setup_s',ps_setup,'PS_receive_s',ps_rx, ...
    'FullFK_setup_s',fk_setup,'FullFK_receive_s',fk_rx,'PS_bytes_by_worker',traffic_ps, ...
    'FullFK_bytes_by_worker',traffic_fk,'PS_checks',local_compare(ids,po,ps,ref.ps,1e-11), ...
    'FullFK_checks',local_compare(ids,fo,fk,ref.fk,1e-11), ...
    'FullFK_term1_checks',local_compare(ids,fo,terms(:,1),ref.terms(:,1),1e-11), ...
    'FullFK_term2_checks',local_compare(ids,fo,terms(:,2),ref.terms(:,2),1e-11));
end

function checks = local_compare(ids,owner,actual,baseline,tolerance)
assert(isequal(size(actual),size(baseline)) && all(isfinite(actual),'all') && all(isfinite(baseline),'all'));
e = abs(actual-baseline);
checks = local_worker_error(ids,owner,vecnorm(e,2,2),vecnorm(baseline,2,2), ...
    max(e,[],2),max(abs(baseline),[],2),tolerance);
end

function checks = local_worker_error(ids,owner,e,scale,peak,ref_peak,tolerance)
assert(isequal(sort(unique(owner)),sort(ids(:))),'Missing or unexpected worker row ownership');
checks = repmat(struct('pid',0,'rows',0,'relative_L2',0,'max_scaled_abs',0),numel(ids),1);
for j = 1:numel(ids)
    rows = owner==ids(j);
    rel = norm(e(rows))/max(norm(scale(rows)),realmin);
    absolute = max(peak(rows))/max(max(ref_peak),1);
    assert(rel<tolerance && absolute<tolerance,'D1b:Equivalence','PID %.0f failed: %.4g, %.4g',ids(j),rel,absolute);
    checks(j) = struct('pid',ids(j),'rows',nnz(rows),'relative_L2',rel,'max_scaled_abs',absolute);
end
end

function token = local_bytes_begin(pool)
token = [];
if ~isempty(pool), token = ticBytes(pool); end
end

function traffic = local_bytes_end(pool,token)
traffic = [];
if ~isempty(pool)
    bytes = tocBytes(pool,token);
    assert(isnumeric(bytes) && isequal(size(bytes),[pool.NumWorkers 2]) && ...
        all(isfinite(bytes),'all') && all(bytes>=0,'all'), ...
        'D1b:TrafficShape','Expected one row per worker and two columns from tocBytes');
    traffic = struct('columns',{{'BytesSentToWorkers','BytesReceivedFromWorkers'}}, ...
        'rows',(1:pool.NumWorkers).','bytes',bytes);
end
end

function s = local_slope(p,x)
[~,s] = d1b_pm_surface(p,x);
end

function local_checkpoint(out,receipt,t0,stage)
receipt.active_stage = stage; receipt.elapsed_s = toc(t0);
local_json(fullfile(out,'checkpoint.json'),receipt);
end

function local_json(path,value)
path = string(path); tmp = path + ".pending";
fid = fopen(tmp,'w'); assert(fid>0);
clean = onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
clear clean
assert(movefile(tmp,path,'f'));
end
