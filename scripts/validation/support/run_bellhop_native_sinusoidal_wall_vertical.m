function result = run_bellhop_native_sinusoidal_wall_vertical(cfg)
%RUN_BELLHOP_NATIVE_SINUSOIDAL_WALL_VERTICAL Run an equivalent native ATI case.
% The physical wall is represented by Bellhop's original C-ATI top boundary,
% using the exact same sampled profile supplied to the internal-wall binary.
arguments
    cfg (1,1) struct
end
required={'bellhop_exe','case_root','run_type','frequency_hz','c0_mps', ...
    'source_depth_m','receiver_depths_m','receiver_ranges_m','beam_count', ...
    'angle_limits_deg','step_m','domain_half_depth_m','source_pattern_angles_deg', ...
    'source_pattern_level_db','wall_profile_r_m','wall_profile_z_m'};
for ii=1:numel(required)
    if ~isfield(cfg,required{ii}), error('Missing cfg.%s.',required{ii}); end
end
if exist(cfg.bellhop_exe,'file')~=2
    error('Native Bellhop executable is missing: %s',cfg.bellhop_exe);
end
out_dir=fileparts(cfg.case_root); if ~exist(out_dir,'dir'), mkdir(out_dir); end
env_file=[cfg.case_root '.env']; ati_file=[cfg.case_root '.ati'];
bty_file=[cfg.case_root '.bty']; sbp_file=[cfg.case_root '.sbp'];
run_type=upper(char(cfg.run_type));
if ~ismember(run_type(1),['C','A']), error('Native comparison supports C or A runs.'); end
source_geometry='R';
if isfield(cfg,'source_geometry'), source_geometry=upper(char(cfg.source_geometry)); end
if ~isscalar(source_geometry) || ~ismember(source_geometry,['R','X'])
    error('source_geometry must be R (point) or X (line).');
end

fid=fopen(sbp_file,'w'); if fid<0, error('Cannot create %s.',sbp_file); end
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%d\n',numel(cfg.source_pattern_angles_deg));
for ii=1:numel(cfg.source_pattern_angles_deg)
    fprintf(fid,'%.12g %.12g\n',cfg.source_pattern_angles_deg(ii),cfg.source_pattern_level_db(ii));
end
clear cleanup

% The ATI points are the same profile samples as the internal wall. Bellhop
% extends the first/last segment to infinity, just as the validation overlay.
r_ati=cfg.wall_profile_r_m(:).'; z_ati=cfg.wall_profile_z_m(:).';
if any(diff(r_ati)<=0)
    error('Native ATI profile ranges must be strictly increasing.');
end
fid=fopen(ati_file,'w'); if fid<0, error('Cannot create %s.',ati_file); end
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'''C''\n%d\n',numel(r_ati));
fprintf(fid,'%.17g %.17g\n',[r_ati(:)/1000 z_ati(:)].');
clear cleanup

d=cfg.domain_half_depth_m;
fid=fopen(env_file,'w'); if fid<0, error('Cannot create %s.',env_file); end
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'''Native sinusoidal ATI wall comparison''\n%.17g\n1\n''CVW *''\n',cfg.frequency_hz);
fprintf(fid,'2 0 %.17g\n-%.17g %.17g /\n%.17g %.17g /\n',d,d,cfg.c0_mps,d,cfg.c0_mps);
fprintf(fid,'''A'' 0\n%.17g %.17g 0 1 0 /\n',d,cfg.c0_mps);
fprintf(fid,'1\n%.17g /\n',cfg.source_depth_m);
fprintf(fid,'%d\n',numel(cfg.receiver_depths_m));
fprintf(fid,'%.17g ',cfg.receiver_depths_m); fprintf(fid,'/\n');
fprintf(fid,'%d\n',numel(cfg.receiver_ranges_m));
fprintf(fid,'%.17g ',cfg.receiver_ranges_m/1000); fprintf(fid,'/\n');
fprintf(fid,'''%s *%s''\n%d\n%.17g %.17g /\n%.17g %.17g %.17g\n', ...
    run_type(1),source_geometry,cfg.beam_count,cfg.angle_limits_deg,cfg.step_m,d+1, ...
    max(cfg.receiver_ranges_m)*1.05/1000);
clear cleanup

% Make the matched flat lower half-space explicit so no bottom reflection is
% part of the native comparison.
fid=fopen(bty_file,'w'); if fid<0, error('Cannot create %s.',bty_file); end
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'''L''\n2\n-1 %.17g\n1 %.17g\n',d,d);
clear cleanup

for ext={'.arr','.shd','.ray','.prt'}
    target=[cfg.case_root ext{1}]; if exist(target,'file')==2, delete(target); end
end
old=pwd; cleanup=onCleanup(@()cd(old)); cd(out_dir);
[~,name]=fileparts(cfg.case_root);
prt_file=[cfg.case_root '.prt'];
if run_type(1)=='C'
    data_file=[cfg.case_root '.shd'];
else
    data_file=[cfg.case_root '.arr'];
end
[status,command_output]=system(sprintf('"%s" "%s"',cfg.bellhop_exe,name));
clear cleanup
if status~=0, error('Native Bellhop failed for %s: %s',name,command_output); end
if exist(prt_file,'file')~=2, error('Native Bellhop did not create %s.',prt_file); end
if run_type(1)=='C'
    if exist(data_file,'file')~=2, error('Native Bellhop did not create %s.',data_file); end
    data=read_bellhop_shd_unfolded_vertical(data_file);
else
    if exist(data_file,'file')~=2, error('Native Bellhop did not create %s.',data_file); end
    data=read_bellhop_arrivals_ascii_vertical(data_file);
end
result=struct('config',cfg,'data',data,'command_output',command_output, ...
    'files',struct('env',env_file,'ati',ati_file,'bty',bty_file,'sbp',sbp_file, ...
    'data',data_file,'prt',prt_file));
end
