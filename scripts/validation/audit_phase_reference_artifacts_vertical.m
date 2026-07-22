function audit = audit_phase_reference_artifacts_vertical(root_dir,run_meta,out_dir)
%AUDIT_PHASE_REFERENCE_ARTIFACTS_VERTICAL Read-only registered MAT audit.
if nargin<1 || isempty(root_dir), root_dir=fileparts(fileparts(fileparts(mfilename('fullpath')))); end
if nargin<2 || isempty(run_meta), run_meta=pe_phase_release_run_meta_vertical(root_dir); end
if nargin<3 || isempty(out_dir)
    out_dir=fullfile(root_dir,'results','validation','pe_phase_release_candidate',run_meta.run_id);
end
if ~exist(out_dir,'dir'), mkdir(out_dir); end
roots={fullfile(root_dir,'results','validation'),fullfile(root_dir,'results','visualization')};
rows=struct('relative_path',{},'bytes',{},'classification',{},'reason',{}, ...
    'schema_version',{},'phase_reference',{},'has_geometry',{},'has_run_meta',{});
for rr=1:numel(roots)
    if ~exist(roots{rr},'dir'), continue; end
    listing=dir(fullfile(roots{rr},'**','*.mat'));
    for ii=1:numel(listing)
        path=fullfile(listing(ii).folder,listing(ii).name);
        if contains(path,[filesep 'archive' filesep]) || contains(path,[filesep run_meta.run_id filesep]), continue; end
        rel=erase(path,[root_dir filesep]); vars=whos('-file',path); names={vars.name};
        [kind,reason]=local_classify(rel,names);
        [schema,phase,has_geometry,has_run]=local_metadata(path,names);
        rows(end+1)=struct('relative_path',rel,'bytes',listing(ii).bytes, ... %#ok<AGROW>
            'classification',kind,'reason',reason,'schema_version',schema, ...
            'phase_reference',phase,'has_geometry',has_geometry,'has_run_meta',has_run);
    end
end
audit=struct('schema_version','1.0.0','created_at',char(datetime('now')), ...
    'validation_run_meta',run_meta,'rows',rows);
table_rows=struct2table(rows);
writetable(table_rows,fullfile(out_dir,'phase_artifact_audit.csv'));
fid=fopen(fullfile(out_dir,'phase_artifact_audit.md'),'w'); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'# PE phase-reference artifact audit\n\n');
fprintf(fid,'- run_id: `%s`\n- read-only scan: yes\n- files: %d\n\n',run_meta.run_id,numel(rows));
fprintf(fid,'| File | Classification | Bytes | Reason |\n|---|---:|---:|---|\n');
for ii=1:numel(rows)
    fprintf(fid,'| `%s` | `%s` | %d | %s |\n',strrep(rows(ii).relative_path,'\','/'), ...
        rows(ii).classification,rows(ii).bytes,strrep(rows(ii).reason,'|','/'));
end
save(fullfile(out_dir,'phase_artifact_audit.mat'),'audit','-v7.3');
end

function [kind,reason]=local_classify(rel,names)
q=lower(strrep(rel,'\','/'));
neutral=contains(q,'joint_pe_cache') || contains(q,'streaming_joint_pe_cache') || ...
    contains(q,'joint_streaming_model') || contains(q,'raw_pm_aperture_audit');
if neutral
    kind='phase_neutral_reusable'; reason='PE/PM/joint spatial operator data; phase transform is receiver-side only';
elseif contains(q,'external') && ~any(strcmp(names,'phase_reference_meta'))
    kind='external_dsp_assumed'; reason='explicit external input without project phase metadata';
elseif contains(q,'conditional_channel_model') || contains(q,'conditional_channel_library')
    kind='migratable_statistics'; reason='conditional mean/C/P/EVD fields can be rotated when geometry is reliable';
elseif any(ismember(names,{'H_direct_f','H_reflect_f','H_ref_sca_f','ensemble'}))
    kind='migratable_components'; reason='separate channel components may be migrated only with reliable geometry';
else
    kind='regenerate_required'; reason='derived result, figure payload, or total-only semantics';
end
end

function [schema,phase,has_geometry,has_run]=local_metadata(path,names)
schema=''; phase=''; has_geometry=false; has_run=any(strcmp(names,'validation_run_meta'));
load_names=intersect(names,{'schema_version','phase_reference_meta','validation_run_meta'});
if ~isempty(load_names)
    try
        s=load(path,load_names{:});
        if isfield(s,'schema_version'), schema=char(string(s.schema_version)); end
        if isfield(s,'phase_reference_meta') && isfield(s.phase_reference_meta,'target_reference')
            phase=s.phase_reference_meta.target_reference;
        end
        if isfield(s,'validation_run_meta')
            has_geometry=isfield(s.validation_run_meta,'geometry');
        end
    catch
        schema='unreadable_metadata';
    end
end
end
