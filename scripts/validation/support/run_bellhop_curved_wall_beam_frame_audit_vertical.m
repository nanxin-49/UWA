function result = run_bellhop_curved_wall_beam_frame_audit_vertical(cfg)
%RUN_BELLHOP_CURVED_WALL_BEAM_FRAME_AUDIT_VERTICAL Run one audit case.
% The audit binary is the Bellhop 2020 validation overlay with optional
% internal-wall and InfluenceGeoHatCart logging.  No production source is used.
arguments
    cfg (1,1) struct
end
required={'bellhop_exe','case_root','audit_mode','audit_receiver_range_m', ...
    'audit_receiver_depth_m'};
for ii=1:numel(required)
    if ~isfield(cfg,required{ii}), error('Missing cfg.%s.',required{ii}); end
end
if exist(cfg.bellhop_exe,'file')~=2
    error('Curved-wall audit executable is missing: %s',cfg.bellhop_exe);
end
out_dir=fileparts(cfg.case_root);
if ~exist(out_dir,'dir'), mkdir(out_dir); end
audit_file=[cfg.case_root '.audit'];
fid=fopen(audit_file,'w');
if fid<0, error('Cannot create %s.',audit_file); end
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n%.17g\n%.17g\n',char(cfg.audit_mode),cfg.audit_receiver_range_m,cfg.audit_receiver_depth_m);
clear cleanup

mode=lower(char(cfg.audit_mode));
switch mode
    case 'rotated'
        base=run_bellhop_internal_sinusoidal_wall_poc_vertical(cfg);
    case 'native'
        base=run_bellhop_native_sinusoidal_wall_vertical(cfg);
    otherwise
        error('audit_mode must be rotated or native.');
end

reflect_file=[cfg.case_root '.iwaudit_reflect'];
rotation_file=[cfg.case_root '.iwaudit_rotation'];
contribution_file=[cfg.case_root '.iwaudit_contribution'];
for path={reflect_file,rotation_file,contribution_file}
    if exist(path{1},'file')~=2, error('Audit binary did not create %s.',path{1}); end
end
result=base;
result.audit=struct('mode',mode, ...
    'reflection',local_read_audit_table(reflect_file), ...
    'rotation',local_read_audit_table(rotation_file), ...
    'contribution',local_read_audit_table(contribution_file));
result.files.audit=audit_file;
result.files.reflection_audit=reflect_file;
result.files.rotation_audit=rotation_file;
result.files.contribution_audit=contribution_file;
end

function table_out=local_read_audit_table(path)
fid=fopen(path,'r');
if fid<0, error('Cannot read %s.',path); end
cleanup=onCleanup(@()fclose(fid));
header={}; modes={}; values=[];
while true
    line=fgetl(fid);
    if ~ischar(line), break; end
    line=strtrim(line);
    if isempty(line), continue; end
    tokens=strsplit(line);
    if tokens{1}(1)=='#'
        header=tokens(2:end);
        continue;
    end
    if isempty(header), error('Missing header in %s.',path); end
    modes{end+1,1}=tokens{1}; %#ok<AGROW>
    row=str2double(tokens(2:end));
    if any(isnan(row)), error('Non-numeric audit row in %s.',path); end
    values(end+1,1:numel(row))=row; %#ok<AGROW>
end
clear cleanup
if isempty(header)
    error('Empty audit table: %s.',path);
end
if isempty(values)
    values=zeros(0,numel(header)-1);
end
if size(values,2)~=numel(header)-1
    error('Audit column mismatch in %s: header=%d data=%d.',path,numel(header),size(values,2));
end
table_out=array2table(values,'VariableNames',header(2:end));
table_out.mode=modes;
table_out=table_out(:,[end 1:end-1]);
end
