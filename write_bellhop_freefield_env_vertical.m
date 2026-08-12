function env_file = write_bellhop_freefield_env_vertical(case_root,cfg)
%WRITE_BELLHOP_FREEFIELD_ENV_VERTICAL Write an omni-style free-space case.
% The upper and lower acousto-elastic halfspaces match the homogeneous
% fluid, following the Acoustic Toolbox tests/BeamPattern/omni.env case.

arguments
    case_root (1,:) char
    cfg (1,1) struct
end
required={'frequency_hz','c0_mps','source_depth_m','receiver_depths_m', ...
    'receiver_ranges_m','run_type','beam_count','angle_limits_deg', ...
    'step_m','domain_half_depth_m'};
for ii=1:numel(required)
    if ~isfield(cfg,required{ii}), error('Missing cfg.%s.',required{ii}); end
end
if ~ismember(upper(char(cfg.run_type)),{'R','A','C','I'})
    error('run_type must be R, A, C, or I.');
end
if any(abs(cfg.receiver_depths_m)>cfg.domain_half_depth_m)
    error('Receiver depths must lie inside the symmetric free-space domain.');
end
env_file=[case_root '.env'];
fid=fopen(env_file,'w');
if fid<0, error('Cannot create Bellhop environment: %s',env_file); end
cleanup=onCleanup(@()fclose(fid));
d=cfg.domain_half_depth_m;
fprintf(fid,'''Point source in free space - PE validation''\n');
fprintf(fid,'%.12g\n1\n''CAF''\n',cfg.frequency_hz);
fprintf(fid,'%.12g %.12g 0.0 1.0 /\n',-d,cfg.c0_mps);
fprintf(fid,'%.12g 0.0 %.12g\n',cfg.c0_mps,d);
fprintf(fid,'%.12g %.12g 0.0 1.0 /\n',-d,cfg.c0_mps);
fprintf(fid,'%.12g /\n''A'' 0.0\n/\n',d);
fprintf(fid,'1\n%.12g /\n',cfg.source_depth_m);
fprintf(fid,'%d\n',numel(cfg.receiver_depths_m));
fprintf(fid,'%.12g ',cfg.receiver_depths_m); fprintf(fid,'/\n');
fprintf(fid,'%d\n',numel(cfg.receiver_ranges_m));
fprintf(fid,'%.12g ',cfg.receiver_ranges_m/1000); fprintf(fid,'/\n');
fprintf(fid,'''%s''\n%d\n%.12g %.12g /\n',upper(char(cfg.run_type)), ...
    cfg.beam_count,cfg.angle_limits_deg);
rbox_km=max(cfg.receiver_ranges_m)*1.05/1000;
fprintf(fid,'%.12g %.12g %.12g\n',cfg.step_m,d+1,rbox_km);
clear cleanup
end
