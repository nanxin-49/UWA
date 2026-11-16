function output = run_bellhop_rotated_field_component(executable,cfg,output_dir,options)
%RUN_BELLHOP_ROTATED_FIELD_COMPONENT One bounded exact-coordinate component.
% Returns a physical z-by-x matrix. Orchestration controls scope/budget; this
% helper does not authorize a full run or a retry of a still-live process.
arguments
    executable (1,:) char
    cfg (1,1) struct
    output_dir (1,:) char
    options.max_seconds (1,1) double {mustBePositive}=20
    options.minimum_available_bytes (1,1) double {mustBePositive}=512*2^20
    options.reuse_dir (1,:) char=''
    options.allow_task_termination (1,1) logical=false
    options.minimum_runtime_available_bytes (1,1) double {mustBeNonnegative}=0
end
assert(~isfolder(output_dir));
[~,before]=memory; assert(before.PhysicalMemory.Available>=options.minimum_available_bytes,'RAM gate.');
mkdir(output_dir); stem=fullfile(output_dir,'case'); started=tic;
diary(fullfile(output_dir,'diary.log')); cleaner=onCleanup(@()diary('off'));
receipt=struct('status','FAILED','result_dir',output_dir,'cfg',cfg, ...
    'executable_sha256',black_file_hash(executable),'matlab_release',version);
try
    if isfield(cfg,'reflection_order')
        kind='reflection'; meta=write_bellhop_rotated_multifield_q(stem,cfg);
    elseif strcmp(cfg.chart_kind,'horizontal')
        kind='horizontal'; meta=write_bellhop_rotated_horizontal_q(stem,cfg,cfg.horizontal_side);
    else
        kind='direct'; meta=write_bellhop_rotated_munk_q(stem,cfg);
        fid=fopen([stem '.iwcov'],'w'); assert(fid>=0); fprintf(fid,'2\n0\n0\n0\n0\n'); fclose(fid);
        write_bellhop_rotated_field_options(stem,cfg.receiver_range_m,meta.q_receiver_m,point3d=true,fast_q=true);
    end
    source_names={'run_bellhop_rotated_field_component','run_bellhop_rotated_local_case', ...
        'write_bellhop_rotated_multifield_q','write_bellhop_rotated_munk_q','write_bellhop_rotated_horizontal_q', ...
        'write_bellhop_rotated_field_options','read_bellhop_rotated_field','munk_profile_sound_speed', ...
        'repair_bellhop_rotated_range_metadata'};
    code_versions=struct('path',{},'sha256',{});
    for b=1:numel(source_names)
        path=which(source_names{b}); code_versions(b)=struct('path',path,'sha256',black_file_hash(path));
    end
    read_stem=stem;
    if ~isempty(options.reuse_dir)
        old=jsondecode(fileread(fullfile(options.reuse_dir,'checkpoint.json')));
        assert(strcmp(old.status,'COMPLETE') && strcmpi(old.executable_sha256,receipt.executable_sha256));
        assert(strcmp(jsonencode(old.cfg),jsonencode(cfg)) && strcmp(jsonencode(old.code_versions),jsonencode(code_versions)),'Reuse config/code drift.');
        for b=1:numel(old.bindings)
            assert(strcmpi(black_file_hash(old.bindings(b).path),old.bindings(b).sha256),'Reuse artifact SHA mismatch.');
        end
        for ext={'.env','.ssp','.iwcov','.rotq','.r25d'}
            if ~isfile([stem ext{1}]), continue; end
            assert(strcmpi(black_file_hash([stem ext{1}]),black_file_hash([old.raw_root ext{1}])),'Reuse input drift.');
        end
        read_stem=old.raw_root;
        process=struct('status','REUSED_VERIFIED_TERMINAL','source_checkpoint',fullfile(options.reuse_dir,'checkpoint.json'), ...
            'elapsed_s',0,'peak_working_set_bytes',NaN);
    else
        process=run_bellhop_rotated_local_case(executable,stem,options.max_seconds, ...
            allow_task_termination=options.allow_task_termination, ...
            minimum_runtime_available_bytes=options.minimum_runtime_available_bytes);
    end
    raw=read_bellhop_rotated_field(read_stem);
    d=readmatrix([read_stem '.r25diag'],'FileType','text','CommentStyle','#'); d=d(~all(isnan(d),2),:);
    assert(height(d)==cfg.beam_count && all(isfinite(raw.pressure),'all'));
    if strcmp(kind,'reflection')
        assert(all(d(:,9)==cfg.expected_surface_bounces & d(:,10)==cfg.expected_bottom_bounces), ...
            'Physical/artificial bounce counts disagree.');
        wall=readmatrix([read_stem '.mwdiag'],'FileType','text','CommentStyle','#'); wall=wall(~all(isnan(wall),2),:);
        assert(all(wall(:,3)==cfg.reflection_order) && max(wall(:,6))<.95*600000,'Reflection count/capacity gate.');
    else
        assert(all(d(:,9:10)==0,'all'),'Artificial boundary reflection.'); wall=[];
    end
    assert(all(d(:,2)>=max(cfg.receiver_range_m)),'Receiver range not reached.');
    pressure_ffp=complex(zeros(numel(cfg.nominal_receiver_depths_m),numel(cfg.nominal_receiver_x_m)));
    if strcmp(kind,'horizontal')
        pressure_ffp(cfg.receiver_depth_order,cfg.receiver_x_order)=-conj(raw.pressure);
    else
        pressure_ffp(cfg.receiver_depth_order,cfg.receiver_x_order)=(-conj(raw.pressure)).';
    end
    [~,after]=memory;
    receipt.raw_root=read_stem; receipt.code_versions=code_versions;
    receipt.status='COMPLETE'; receipt.elapsed_s=toc(started); receipt.failed_checks=0;
    receipt.kind=kind; receipt.process=process;
    receipt.available_before_bytes=before.PhysicalMemory.Available;
    receipt.available_after_bytes=after.PhysicalMemory.Available;
    receipt.max_ray_end_unit_error=max(d(:,8)); receipt.max_ray_range_m=max(d(:,2));
    save(fullfile(output_dir,'component.mat'),'pressure_ffp','cfg','meta','receipt','d','wall','-v7.3');
    bindings=struct('path',{},'sha256',{});
    for ext={'.env','.ssp','.iwcov','.rotq','.r25d','.shd','.prt','.r25diag','.mwdiag','.covdiag'}
        path=[read_stem ext{1}]; if ~isfile(path), continue; end
        bindings(end+1)=struct('path',path,'sha256',black_file_hash(path)); %#ok<AGROW>
    end
    path=fullfile(output_dir,'component.mat');
    bindings(end+1)=struct('path',path,'sha256',black_file_hash(path));
    receipt.bindings=bindings;
    output=struct('pressure_ffp',pressure_ffp,'receipt',receipt,'cfg',cfg,'meta',meta,'d',d,'wall',wall);
catch exception
    receipt.elapsed_s=toc(started); receipt.error=getReport(exception,'extended','hyperlinks','off');
    write(fullfile(output_dir,'checkpoint.json'),jsonencode(receipt,PrettyPrint=true));
    rethrow(exception)
end
write(fullfile(output_dir,'checkpoint.json'),jsonencode(receipt,PrettyPrint=true));
fprintf('COMPONENT %s points=%d elapsed=%.3fs RAM_before=%.3fGiB\n',kind,numel(pressure_ffp),receipt.elapsed_s,receipt.available_before_bytes/2^30);
end
function write(path,content)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); cleaner=onCleanup(@()fclose(fid)); fprintf(fid,'%s',content);
end