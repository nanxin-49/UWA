function [env_file,sbp_file] = write_bellhop_unfolded_gaussian_env_vertical(case_root,cfg)
%WRITE_BELLHOP_UNFOLDED_GAUSSIAN_ENV_VERTICAL Write a matched-halfspace
% Bellhop environment for the unfolded (vertical-to-range) comparison.
% The source pattern is written in the Bellhop .sbp format and is centered
% on the positive range axis.  The artificial depth coordinate is one PE
% transverse coordinate; the matched halfspaces prevent image reflections.

arguments
    case_root (1,:) char
    cfg (1,1) struct
end
required={'frequency_hz','c0_mps','source_depth_m','receiver_depths_m', ...
    'receiver_ranges_m','beam_count','angle_limits_deg','step_m', ...
    'domain_half_depth_m','source_pattern_angles_deg','source_pattern_level_db'};
for ii=1:numel(required)
    if ~isfield(cfg,required{ii}), error('Missing cfg.%s.',required{ii}); end
end
run_type=upper(char(cfg.run_type));
if isempty(run_type), run_type='C'; end
if ~ismember(run_type(1),['C','A','R','I'])
    error('run_type must start with C, A, R, or I.');
end
if numel(cfg.source_pattern_angles_deg)~=numel(cfg.source_pattern_level_db) || ...
        any(~isfinite(cfg.source_pattern_angles_deg(:))) || any(~isfinite(cfg.source_pattern_level_db(:)))
    error('The .sbp angle and level vectors must be finite and have equal length.');
end
if any(abs(cfg.receiver_depths_m)>cfg.domain_half_depth_m)
    error('Receiver depths must lie inside the matched free-space domain.');
end
out_dir=fileparts(case_root); if ~exist(out_dir,'dir'), mkdir(out_dir); end
env_file=[case_root '.env']; sbp_file=[case_root '.sbp'];
fid=fopen(sbp_file,'w'); if fid<0, error('Cannot create %s.',sbp_file); end
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%d\n',numel(cfg.source_pattern_angles_deg));
for ii=1:numel(cfg.source_pattern_angles_deg)
    fprintf(fid,'%.12g %.12g\n',cfg.source_pattern_angles_deg(ii),cfg.source_pattern_level_db(ii));
end
clear cleanup

fid=fopen(env_file,'w'); if fid<0, error('Cannot create %s.',env_file); end
cleanup=onCleanup(@()fclose(fid)); d=cfg.domain_half_depth_m;
fprintf(fid,'''Unfolded Gaussian PE--Bellhop validation''\n');
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
% Bellhop reads a source beam pattern when the third RunType character is *.
fprintf(fid,'''%s *''\n%d\n%.12g %.12g /\n',run_type(1), ...
    cfg.beam_count,cfg.angle_limits_deg);
rbox_km=max(cfg.receiver_ranges_m)*1.05/1000;
fprintf(fid,'%.12g %.12g %.12g\n',cfg.step_m,d+1,rbox_km);
clear cleanup
end
