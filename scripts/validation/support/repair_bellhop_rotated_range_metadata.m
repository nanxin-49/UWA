function meta = repair_bellhop_rotated_range_metadata(case_root,ranges_m)
%REPAIR_BELLHOP_ROTATED_RANGE_METADATA Canonical float32 ENV labels only.
% Bellhop reads km into real32 then multiplies by1000 in real32. At >65km
% this can collapse .01m labels although exact physical r25d samples differ.
% Keep native monotonicity and the existing2ULP metadata binding; pressure
% coordinates remain byte-identical in r25d. Only C/R point exact-mode input.
arguments
    case_root (1,:) char
    ranges_m double {mustBeFinite,mustBePositive}
end
r=ranges_m(:); assert(all(diff(r)>0));
km=single(r/1000); before=single(km.*single(1000));
meta=struct('status','NO_REPAIR_REQUIRED','duplicate_count_before',nnz(diff(before)<=0), ...
    'physical_receiver_coordinates_changed',false);
if all(diff(before)>0), return; end
side=[case_root '.r25d']; assert(isfile(side),'Exact coordinate sidecar required.');
fid=fopen(side,'r'); assert(fid>=0); closer=onCleanup(@()fclose(fid));
flags=fscanf(fid,'%d',2); counts=fscanf(fid,'%d',2); exact=fscanf(fid,'%f',counts(1)); clear closer
assert(flags(1)==1 && counts(1)==numel(r) && isequal(exact,r),'Point/exact coordinate binding mismatch.');
side_sha=black_file_hash(side);
file=[case_root '.env']; original=fileread(file); lines=splitlines(string(original));
nonempty=find(strlength(strtrim(lines))>0);
assert(numel(nonempty)>=6);
mode_line=nonempty(end-3); range_line=nonempty(end-4); count_line=nonempty(end-5);
assert(~isempty(regexp(strtrim(lines(mode_line)),"^'C\s+R'\s*$",'once')),'Only exact geometric point C/R fields supported.');
assert(str2double(strtrim(lines(count_line)))==numel(r));
written=sscanf(char(lines(range_line)),'%f');
assert(numel(written)==numel(r) && max(abs(written-r/1000))<1e-10,'Unexpected ENV range labels.');
fixed=before; steps=zeros(size(r));
for j=2:numel(r)
    while fixed(j)<=fixed(j-1)
        bits=typecast(km(j),'uint32'); km(j)=typecast(bits+uint32(1),'single');
        fixed(j)=single(km(j)*single(1000)); steps(j)=steps(j)+1;
        assert(steps(j)<=4,'Native label allocation would exceed a local quantization correction.');
    end
end
error_m=abs(double(fixed)-r); tolerance=max(1e-9,2*double(eps(fixed)));
assert(all(diff(fixed)>0) && all(error_m<=tolerance),'Cannot preserve strict metadata order within existing2ULP bound.');
backup=[case_root '.env.before_range_metadata_repair']; assert(~isfile(backup)); copyfile(file,backup);
lines(range_line)=string([sprintf('%.17g ',double(km)) '/']);
updated=char(join(lines,newline));
fid=fopen(file,'w','n','UTF-8'); assert(fid>=0); closer=onCleanup(@()fclose(fid)); fprintf(fid,'%s',updated); clear closer
assert(strcmp(black_file_hash(side),side_sha),'Physical receiver sidecar changed.');
meta=struct('status','NATIVE_RANGE_LABEL_QUANTIZATION_REPAIRED','duplicate_count_before',nnz(diff(before)<=0), ...
    'duplicate_count_after',nnz(diff(fixed)<=0),'adjusted_labels',nnz(steps), ...
    'maximum_error_m',max(error_m),'maximum_error_over_ulp',max(error_m./double(eps(fixed))), ...
    'physical_receiver_coordinates_changed',false,'r25d_sha256',side_sha, ...
    'env_before_sha256',black_file_hash(backup),'env_after_sha256',black_file_hash(file), ...
    'scope','Float32 metadata labels only; native strict monotonicity and2ULP binding kept; no pressure interpolation, gain or phase fit');
save([case_root '.range_metadata.mat'],'r','km','before','fixed','steps','meta');
fid=fopen([case_root '.range_metadata.json'],'w'); assert(fid>=0); closer=onCleanup(@()fclose(fid)); fprintf(fid,'%s',jsonencode(meta,PrettyPrint=true));
end