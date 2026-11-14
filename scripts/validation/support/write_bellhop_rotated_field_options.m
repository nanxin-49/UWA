function path = write_bellhop_rotated_field_options(case_root,ranges_m,depths_m,options)
%WRITE_BELLHOP_ROTATED_FIELD_OPTIONS Exact evaluator coordinates, no resampling.
arguments
    case_root (1,:) char
    ranges_m double {mustBeFinite,mustBePositive}
    depths_m double {mustBeFinite}
    options.point3d (1,1) logical = true
    options.fast_q (1,1) logical = true
end
assert(isvector(ranges_m) && isvector(depths_m));
assert(all(diff(ranges_m(:))>0) && all(diff(depths_m(:))>0));
path=[case_root '.r25d']; assert(~isfile(path),'Refusing to overwrite coordinate options.');
fid=fopen(path,'w'); assert(fid>=0); closer=onCleanup(@()fclose(fid));
fprintf(fid,'%d %d\n%d %d\n',options.point3d,options.fast_q,numel(ranges_m),numel(depths_m));
fprintf(fid,'%.17g ',ranges_m); fprintf(fid,'\n'); fprintf(fid,'%.17g ',depths_m); fprintf(fid,'\n');
end