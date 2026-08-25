function preparation = prepare_pe_phase_reference_release_candidate_vertical()
%PREPARE_PE_PHASE_REFERENCE_RELEASE_CANDIDATE_VERTICAL Audit and archive old outputs.
root_dir=fileparts(fileparts(fileparts(mfilename('fullpath')))); addpath(root_dir); setup_vertical_project();
f64=linspace(4000,8000,64).';
seed_definition=struct('phase',810001,'adjoint','script_fixed_seeds', ...
    'u5','1100001..1900001','u8','2100001..2900001','two_node','script_fixed_seeds');
run_meta=pe_phase_release_run_meta_vertical(root_dir,struct('create',true, ...
    'frequency_axis_hz',f64,'seed_definition',seed_definition));
run_dir=fullfile(root_dir,'results','validation','pe_phase_release_candidate',run_meta.run_id);
if ~exist(run_dir,'dir'), mkdir(run_dir); end
pre_audit=audit_phase_reference_artifacts_vertical(root_dir,run_meta,run_dir);

stamp=char(datetime('now','Format','yyyyMMdd_HHmmss'));
archive_root=fullfile(root_dir,'results','archive','pe_phase_reference_pre_rc',stamp);
items=local_archive_items(root_dir);
moves=struct('source',{},'destination',{},'bytes',{});
for ii=1:numel(items)
    src=items{ii}; if ~exist(src,'file'), continue; end
    rel=erase(src,[root_dir filesep]); dst=fullfile(archive_root,rel);
    local_assert_under(src,root_dir); local_assert_under(dst,archive_root);
    if exist(dst,'file'), error('Archive destination already exists: %s',dst); end
    parent=fileparts(dst); if ~exist(parent,'dir'), mkdir(parent); end
    info=dir(src); [ok,msg]=movefile(src,dst,'f');
    if ~ok, error('Failed to archive %s: %s',src,msg); end
    moves(end+1)=struct('source',rel,'destination',erase(dst,[root_dir filesep]), ... %#ok<AGROW>
        'bytes',info.bytes);
end
migration_audit=local_migrate_legacy_models(root_dir,archive_root,run_dir,run_meta);
schema_version='1.0.0'; phase_reference_meta=local_phase_meta(run_meta.geometry);
validation_run_meta=run_meta;
preparation=struct('schema_version',schema_version,'phase_reference_meta',phase_reference_meta, ...
    'validation_run_meta',validation_run_meta,'archive_root',archive_root, ...
    'audit_file',fullfile(run_dir,'phase_artifact_audit.mat'),'moves',moves, ...
    'migration_audit',migration_audit);
save(fullfile(run_dir,'release_candidate_preparation.mat'),'preparation', ...
    'schema_version','phase_reference_meta','validation_run_meta','-v7.3');
if ~isempty(moves), writetable(struct2table(moves),fullfile(run_dir,'archive_manifest.csv')); end
local_write_manifest(preparation,fullfile(run_dir,'archive_manifest.md'));
fprintf('Prepared %s; archived %d phase-sensitive files to %s\n', ...
    run_meta.run_id,numel(moves),archive_root);
end

function audit=local_migrate_legacy_models(root_dir,archive_root,run_dir,run_meta)
rows=struct('source',{},'kind',{},'migrated',{},'max_closure_error',{},'message',{});
if exist(archive_root,'dir'), d=dir(fullfile(archive_root,'**','*.mat')); else, d=[]; end
for ii=1:numel(d)
    if d(ii).isdir || ~contains(lower(d(ii).name),'conditional'), continue; end
    src=fullfile(d(ii).folder,d(ii).name); vars={whos('-file',src).name};
    kind=''; if any(strcmp(vars,'model')), kind='model'; elseif any(strcmp(vars,'library')), kind='library'; end
    if isempty(kind), continue; end
    try
        q=load(src,kind); old=q.(kind);
        if ~isfield(old,'schema_version') || ~startsWith(old.schema_version,'1.'), continue; end
        [new,meta]=upgrade_conditional_channel_phase_vertical(old,run_meta.geometry);
        closure=local_migration_closure(old,new);
        rel=erase(src,[archive_root filesep]); dst=fullfile(run_dir,'migrated_legacy_models',rel);
        parent=fileparts(dst); if ~exist(parent,'dir'), mkdir(parent); end
        schema_version=new.schema_version; phase_reference_meta=new.phase_reference_meta;
        validation_run_meta=run_meta;
        if strcmp(kind,'model'), model=new; save(dst,'model','schema_version','phase_reference_meta','validation_run_meta','meta','closure','-v7.3');
        else, library=new; save(dst,'library','schema_version','phase_reference_meta','validation_run_meta','meta','closure','-v7.3'); end
        rows(end+1)=struct('source',erase(src,[root_dir filesep]),'kind',kind, ... %#ok<AGROW>
            'migrated',true,'max_closure_error',closure,'message','audit-only copy');
    catch exception
        rows(end+1)=struct('source',erase(src,[root_dir filesep]),'kind',kind, ... %#ok<AGROW>
            'migrated',false,'max_closure_error',Inf,'message',exception.message);
    end
end
audit=struct('rows',rows,'migrated_count',sum([rows.migrated]), ...
    'failed_count',sum(~[rows.migrated]));
if ~isempty(rows), writetable(struct2table(rows),fullfile(run_dir,'legacy_migration_audit.csv')); end
end

function e=local_migration_closure(old,new)
if isfield(old,'conditions')
    e=0; for ii=1:numel(old.conditions), e=max(e,local_migration_closure(old.conditions(ii).model,new.conditions(ii).model)); end
    return
end
d=new.phase_reference_meta.reflect_dsp_factor_f(:); s0=old.stats; s1=new.stats;
T=[diag(real(d)),-diag(imag(d));diag(imag(d)),diag(real(d))];
errors=[local_rel(s1.mu_scatter_f,d.*s0.mu_scatter_f), ...
    local_rel(s1.C_scatter_f,d.*s0.C_scatter_f.*conj(d.')), ...
    local_rel(s1.P_scatter_f,d.*s0.P_scatter_f.*d.'), ...
    local_rel(s1.eigenvectors,d.*s0.eigenvectors), ...
    local_rel(s1.augmented_mean,T*s0.augmented_mean), ...
    local_rel(s1.augmented_covariance,T*s0.augmented_covariance*T.'), ...
    local_rel(s1.augmented_eigenvectors,T*s0.augmented_eigenvectors)];
e=max(errors);
end

function e=local_rel(a,b), e=norm(a(:)-b(:))/max(norm(b(:)),eps); end

function items=local_archive_items(root_dir)
v=fullfile(root_dir,'results','validation'); vis=fullfile(root_dir,'results','visualization');
whole_dirs={'adjoint_pe_receiver_projection','conditional_channel_library', ...
    'pe_channel_phase_reference','pe_phase_convention_uniform','two_node_communication'};
items={};
for ii=1:numel(whole_dirs), items=[items;local_files(fullfile(v,whole_dirs{ii}))]; end %#ok<AGROW>
items=[items;local_files(fullfile(vis,'pe_propagation_atlas'))];
for name={'u5_conditional_channel_f64','u8_conditional_channel_f64'}
    files=local_files(fullfile(v,name{1}));
    for jj=1:numel(files)
        q=lower(files{jj});
        keep=contains(q,'joint_pe_cache.mat') || contains(q,'streaming_joint_pe_cache.mat') || ...
            contains(q,'joint_streaming_model_u5.mat') || contains(q,'joint_streaming_model_u5_meta.mat') || ...
            contains(q,'active_k_audit.txt') || contains(q,'u8_raw_pm_aperture_audit.');
        if ~keep, items{end+1,1}=files{jj}; end %#ok<AGROW>
    end
end
patterns={'lfm_*','validate_kstat_vs_kdomain_lfm_*','validate_comm_link_minimal_vertical_*', ...
    'validate_public_channel_modes_vertical_result.mat'};
for pp=1:numel(patterns)
    d=dir(fullfile(v,patterns{pp}));
    for jj=1:numel(d), if ~d(jj).isdir, items{end+1,1}=fullfile(d(jj).folder,d(jj).name); end, end %#ok<AGROW>
end
cached_public=fullfile(v,'cached_joint_kstat_pe_receiver','cached_public_consistency_double.mat');
if isfile(cached_public), items{end+1,1}=cached_public; end
items=unique(items,'stable');
end

function files=local_files(folder)
files={}; if ~exist(folder,'dir'), return; end
d=dir(fullfile(folder,'**','*'));
for ii=1:numel(d), if ~d(ii).isdir, files{end+1,1}=fullfile(d(ii).folder,d(ii).name); end, end %#ok<AGROW>
end

function local_assert_under(path,parent)
p=char(java.io.File(path).getCanonicalPath()); root=char(java.io.File(parent).getCanonicalPath());
if ~(strcmpi(p,root) || startsWith(lower(p),[lower(root) filesep]))
    error('Path escapes intended root: %s',path);
end
end

function m=local_phase_meta(g)
f=[4000;8000]; [~,m]=apply_pe_channel_phase_reference_vertical(f, ...
    struct('direct_f',ones(2,1),'reflect_fm',ones(2,1)),g,'direct_dsp');
end

function local_write_manifest(p,path)
fid=fopen(path,'w'); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'# PE phase-reference pre-RC archive\n\n');
fprintf(fid,'- run_id: `%s`\n- archive: `%s`\n- moved files: %d\n\n', ...
    p.validation_run_meta.run_id,strrep(p.archive_root,'\','/'),numel(p.moves));
fprintf(fid,'Old files were moved, not deleted or overwritten. Relative paths are preserved.\n');
end
