function meta = write_bellhop_rotated_munk_q(case_root,cfg)
%WRITE_BELLHOP_ROTATED_MUNK_Q Write a direct-only proper-rotation Q SSP case.
% Physical (x,z) maps to (r,q)=(zs-z,q0+x) up and (z-zs,q0-x) down.
% Both maps have determinant +1. No physical surface/bottom is present here.
arguments
    case_root (1,:) char
    cfg (1,1) struct
end
assert(~isfile([case_root '.env']),'Refusing to overwrite an existing case.');
assert(ismember(cfg.run_type,{'A','C','R'}),'Unsupported run type.');
assert(ismember(cfg.source_geometry,{'R','X'}),'Invalid source geometry.');
assert(mod(cfg.beam_count,2)==1 && cfg.beam_count>=3,'Use an odd beam fan including the axis.');
d=cfg.artificial_halfwidth_m; q0=cfg.transverse_origin_m;
assert(d>max(abs(cfg.receiver_x_m)) && q0>=d,'Invalid artificial domain.');
r=unique([0:cfg.ssp_spacing_m:cfg.range_limit_m,cfg.range_limit_m]).';
sgn=1; if strcmp(cfg.direction,'up'), sgn=-1; end
z=cfg.source_z_m+sgn*r;
assert(all(z>=0 & z<=cfg.water_depth_m),'Q SSP extends outside physical water.');
c=munk_profile_sound_speed(z,cfg.munk_parameters);
if cfg.uniform_control, c(:)=cfg.uniform_c_mps; end
physical_sample_count=numel(r);
% Numeric guard lies beyond Beam.Box.r; no receiver/physical result uses it.
r=[r;cfg.range_limit_m+max(4*cfg.step_m,1)];
z=[z;z(end)]; c=[c;c(end)];
assert(all(cfg.receiver_range_m>0 & cfg.receiver_range_m<cfg.range_limit_m));
out_dir=fileparts(case_root); if ~isfolder(out_dir), mkdir(out_dir); end
fid=fopen([case_root '.ssp'],'w'); assert(fid>=0); closer=onCleanup(@()fclose(fid));
fprintf(fid,'%d\n',numel(r)); fprintf(fid,'%.17g ',r/1000); fprintf(fid,'\n');
% Quad reads one row per artificial depth, with Nr range columns.
fprintf(fid,'%.17g ',c); fprintf(fid,'\n'); fprintf(fid,'%.17g ',c); fprintf(fid,'\n');
clear closer
fid=fopen([case_root '.env'],'w'); assert(fid>=0); closer=onCleanup(@()fclose(fid));
fprintf(fid,'''Rotated Munk direct Q diagnostic; no physical reflections''\n');
fprintf(fid,'%.17g\n1\n''QVW''\n2 0 %.17g\n',cfg.frequency_hz,q0+d);
fprintf(fid,'%.17g %.17g 0 1 0 0 /\n',q0-d,c(1));
fprintf(fid,'%.17g %.17g 0 1 0 0 /\n',q0+d,c(1));
fprintf(fid,'''R'' 0\n1\n%.17g /\n%d\n',q0,numel(cfg.receiver_x_m));
% q=q0-sgn*x is the transverse coordinate in the proper rotation.
q=q0-sgn*cfg.receiver_x_m;
assert(all(diff(q)>=0),'Receiver transverse coordinates must increase in Bellhop chart.');
fprintf(fid,'%.17g ',q); fprintf(fid,'/\n%d\n',numel(cfg.receiver_range_m));
fprintf(fid,'%.17g ',cfg.receiver_range_m/1000); fprintf(fid,'/\n');
fprintf(fid,'''%s  %s''\n%d\n%.17g %.17g /\n',cfg.run_type,cfg.source_geometry, ...
    cfg.beam_count,cfg.angle_limits_deg);
fprintf(fid,'%.17g %.17g %.17g\n',cfg.step_m,q0+d+1,cfg.range_limit_m/1000);
clear closer
meta=struct('range_m',r,'physical_depth_m',z,'c_mps',c,'q_receiver_m',q, ...
    'physical_sample_count',physical_sample_count,'numeric_guard_only',true,'physical_to_chart_matrix',[0 sgn;-sgn 0], ...
    'physical_source',[0 cfg.source_z_m],'chart_source',[0 q0], ...
    'components','direct-only; no total field','ssp_file',[case_root '.ssp'], ...
    'env_file',[case_root '.env'],'source_geometry',cfg.source_geometry);
end