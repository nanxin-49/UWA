function d1b_scale_phase(cfg,stage)
%D1B_SCALE_PHASE Client-only phase/traffic receipt; no numerical operations.
persistent prior
if ~isfield(cfg,'monitor_dir'), return, end
folder = char(cfg.monitor_dir);
assert(isfolder(folder));
if ~isempty(prior)
    assert(strcmp(prior.folder,folder));
    rec = struct('case_id',prior.case_id,'phase',prior.stage,'wall_s',toc(prior.clock), ...
        'main_pid',feature('getpid'),'bytes_by_worker',[]);
    if ~isempty(prior.pool)
        bytes = tocBytes(prior.pool,prior.token);
        assert(isnumeric(bytes) && isequal(size(bytes),[2 2]) && all(isfinite(bytes),'all') && all(bytes>=0,'all'));
        rec.bytes_by_worker = bytes;
    end
    fid = fopen(fullfile(folder,'phase_events.jsonl'),'a'); assert(fid>0);
    fprintf(fid,'%s\n',jsonencode(rec)); fclose(fid);
end
prior = [];
body = struct('case_id',cfg.case_id,'phase',char(stage),'main_pid',feature('getpid'));
tmp = fullfile(folder,'phase.pending');
fid = fopen(tmp,'w'); assert(fid>0); fprintf(fid,'%s\n',jsonencode(body)); fclose(fid);
assert(movefile(tmp,fullfile(folder,'phase.json'),'f'));
if strcmp(stage,'IDLE'), return, end
pool = []; token = [];
% pool_close measures deletion itself, so no live pool remains at phase end.
if ~strcmp(stage,'pool_close'), pool = gcp('nocreate'); end
if ~isempty(pool), assert(pool.NumWorkers==2); token = ticBytes(pool); end
prior = struct('folder',folder,'case_id',cfg.case_id,'stage',char(stage), ...
    'clock',tic,'pool',pool,'token',token);
end
