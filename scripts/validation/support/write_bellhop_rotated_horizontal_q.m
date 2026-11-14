function meta = write_bellhop_rotated_horizontal_q(case_root,cfg,side)
%WRITE_BELLHOP_ROTATED_HORIZONTAL_Q Direct source-plane chart, natural c(q).
% right: r=x, q=z-zs; left: r=-x, q=-(z-zs), both det=+1.
arguments
    case_root (1,:) char
    cfg (1,1) struct
    side (1,1) double {mustBeMember(side,[-1 1])}
end
assert(~isfile([case_root '.env'])); folder=fileparts(case_root);
if ~isfolder(folder), mkdir(folder); end
q=(-20:.25:20).'; z=cfg.source_z_m+side*q;
c=munk_profile_sound_speed(z,cfg.munk_parameters);
if cfg.uniform_control, c(:)=cfg.uniform_c_mps; end
r=[0;cfg.range_limit_m+1];
fid=fopen([case_root '.ssp'],'w'); assert(fid>=0); close_file=onCleanup(@()fclose(fid));
fprintf(fid,'2\n%.17g %.17g\n',r/1000);
fprintf(fid,'%.17g %.17g\n',[c c].'); clear close_file
fid=fopen([case_root '.env'],'w'); assert(fid>=0); close_file=onCleanup(@()fclose(fid));
fprintf(fid,'''Rotated Munk source-plane horizontal chart; direct only''\n');
fprintf(fid,'%.17g\n1\n''QVW''\n%d 0 %.17g\n',cfg.frequency_hz,numel(q),q(end));
fprintf(fid,'%.17g %.17g 0 1 0 0 /\n',[q c].');
fprintf(fid,'''R'' 0\n1\n0 /\n%d\n',numel(cfg.receiver_q_m));
fprintf(fid,'%.17g ',cfg.receiver_q_m); fprintf(fid,'/\n%d\n',numel(cfg.receiver_range_m));
fprintf(fid,'%.17g ',cfg.receiver_range_m/1000); fprintf(fid,'/\n');
fprintf(fid,'''C  R''\n%d\n%.17g %.17g /\n',cfg.beam_count,cfg.angle_limits_deg);
fprintf(fid,'%.17g 21 %.17g\n',cfg.step_m,cfg.range_limit_m/1000); clear close_file
fid=fopen([case_root '.iwcov'],'w'); assert(fid>=0); fprintf(fid,'4\n0\n0\n0\n0\n'); fclose(fid);
write_bellhop_rotated_field_options(case_root,cfg.receiver_range_m,cfg.receiver_q_m, ...
    point3d=cfg.point3d,fast_q=true);
meta=struct('side',side,'physical_depth_m',z,'q_nodes_m',q,'c_mps',c, ...
    'physical_x_m',side*cfg.receiver_range_m,'physical_z_m',cfg.source_z_m+side*cfg.receiver_q_m, ...
    'source_z_m',cfg.source_z_m,'component','direct source-plane chart only', ...
    'source_singularity','not evaluated; all chart ranges positive');
end
