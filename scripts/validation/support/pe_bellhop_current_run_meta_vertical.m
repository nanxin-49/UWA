function meta = pe_bellhop_current_run_meta_vertical(root_dir,options)
%PE_BELLHOP_CURRENT_RUN_META_VERTICAL Create a reproducible Bellhop run ID.

arguments
    root_dir (1,:) char
    options (1,1) struct = struct()
end
if ~isfield(options,'run_id') || isempty(options.run_id)
    options.run_id=['bellhop_current_' char(datetime('now','Format','yyyyMMdd_HHmmss'))];
end
if ~isfield(options,'bellhop_exe') || isempty(options.bellhop_exe)
    options.bellhop_exe=getenv('BELLHOP_EXE');
end
if isempty(options.bellhop_exe) || exist(options.bellhop_exe,'file')~=2
    error('pe_bellhop_current_run_meta_vertical:MissingBellhop', ...
        'Set BELLHOP_EXE to a stable absolute bellhop.exe path.');
end
bellhop_exe=char(java.io.File(options.bellhop_exe).getCanonicalPath());
if ~local_is_absolute(bellhop_exe)
    error('Bellhop executable path must be absolute: %s',bellhop_exe);
end
temp_roots={tempdir,getenv('TEMP'),getenv('TMP')};
for ii=1:numel(temp_roots)
    if ~isempty(temp_roots{ii}) && local_is_under(bellhop_exe,temp_roots{ii})
        error('pe_bellhop_current_run_meta_vertical:TemporaryBellhop', ...
            'Formal Bellhop validation rejects temporary executables: %s',bellhop_exe);
    end
end

files={ ...
    'setup_vertical_project.m'; ...
    fullfile('src','channel','apply_pe_channel_phase_reference_vertical.m'); ...
    fullfile('src','channel','build_channel_cir_vertical.m'); ...
    'vertical_channel_model.m'; ...
    fullfile('src','channel','vertical_channel_model_impl.m'); ...
    fullfile('src','propagation','vertical_wape_propagator.m'); ...
    fullfile('scripts','validation','support','run_bellhop_flat_surface_arrivals_vertical.m'); ...
    fullfile('scripts','validation','support','pe_bellhop_current_run_meta_vertical.m'); ...
    fullfile('scripts','validation','validate_pe_phase_convention_uniform_vertical.m'); ...
    fullfile('scripts','validation','validate_pe_bellhop_flat_surface_current_vertical.m'); ...
    fullfile('scripts','reporting','generate_pe_bellhop_comparison_atlas_vertical.m')};
code_hash=local_hash_sources(root_dir,files);
[git_status,revision]=system(sprintf('git -C "%s" rev-parse HEAD',root_dir));
if git_status~=0, revision='unavailable'; else, revision=strtrim(revision); end
[~,dirty]=system(sprintf('git -C "%s" status --porcelain --untracked-files=no',root_dir));

validation_root=fullfile(root_dir,'results','validation', ...
    'pe_bellhop_flat_surface_current',options.run_id);
visual_root=fullfile(root_dir,'results','visualization', ...
    'pe_bellhop_flat_surface_current',options.run_id);
if ~exist(validation_root,'dir'), mkdir(validation_root); end
if ~exist(visual_root,'dir'), mkdir(visual_root); end

bellhop_info=dir(bellhop_exe);
meta=struct( ...
    'schema_version','1.0.0', ...
    'run_id',options.run_id, ...
    'created_at',char(datetime('now','TimeZone','local', ...
        'Format','yyyy-MM-dd''T''HH:mm:ssXXX')), ...
    'code_revision',revision, ...
    'code_dirty',~isempty(strtrim(dirty)), ...
    'code_fingerprint',code_hash, ...
    'code_fingerprint_files',{files}, ...
    'matlab_version',version, ...
    'geometry',struct('z_tx',80,'z_rx',10,'z_surface',0,'c0',1500), ...
    'frequency_axis_hz',linspace(3000,5000,33).', ...
    'seed_definition',struct('deterministic',true,'random_surface',false), ...
    'phase_reference','direct_dsp', ...
    'bellhop_executable',bellhop_exe, ...
    'bellhop_sha256',local_hash_file(bellhop_exe), ...
    'bellhop_bytes',bellhop_info.bytes, ...
    'validation_output_dir',validation_root, ...
    'visualization_output_dir',visual_root, ...
    'status','ACTIVE');

state_dir=fileparts(validation_root);
save(fullfile(state_dir,'current_run.mat'),'meta');
save(fullfile(validation_root,'validation_run_meta.mat'),'meta');
end

function tf=local_is_absolute(path)
if ispc
    tf=~isempty(regexp(path,'^[A-Za-z]:[\\/]|^\\\\','once'));
else
    tf=startsWith(path,'/');
end
end

function tf=local_is_under(path,parent)
path=lower(char(java.io.File(path).getCanonicalPath()));
parent=lower(char(java.io.File(parent).getCanonicalPath()));
if ~endsWith(parent,filesep), parent=[parent filesep]; end
tf=startsWith(path,parent);
end

function hash=local_hash_sources(root_dir,files)
md=java.security.MessageDigest.getInstance('SHA-256');
for ii=1:numel(files)
    path=fullfile(root_dir,files{ii});
    if exist(path,'file')~=2
        error('Fingerprint source is missing: %s',path);
    end
    md.update(typecast(uint8(files{ii}),'int8'));
    fid=fopen(path,'rb'); cleanup=onCleanup(@()fclose(fid));
    bytes=fread(fid,Inf,'*uint8'); clear cleanup
    md.update(typecast(bytes,'int8'));
end
hash=local_digest(md);
end

function hash=local_hash_file(path)
md=java.security.MessageDigest.getInstance('SHA-256');
fid=fopen(path,'rb'); cleanup=onCleanup(@()fclose(fid));
while true
    bytes=fread(fid,1024*1024,'*uint8');
    if isempty(bytes), break; end
    md.update(typecast(bytes,'int8'));
end
clear cleanup
hash=local_digest(md);
end

function hash=local_digest(md)
digest=typecast(md.digest(),'uint8');
hash=lower(reshape(dec2hex(digest,2).',1,[]));
end
