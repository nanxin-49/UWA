function shade = read_bellhop_rotated_field(case_root)
%READ_BELLHOP_ROTATED_FIELD Read pressure with both exact coordinate bindings.
% A monotonicity-only SHD range heuristic is unsafe when NRz > NRr.
% Native float32 coordinate metadata is retained; pressure is never resampled.
arguments
    case_root (1,:) char
end
case_root=bellhop_linux_path(case_root); side=[case_root '.r25d']; assert(isfile(side));
fid=fopen(side,'r'); assert(fid>=0); close_side=onCleanup(@()fclose(fid));
flags=fscanf(fid,'%d',2); counts=fscanf(fid,'%d',2);
r=fscanf(fid,'%f',counts(1)); q=fscanf(fid,'%f',counts(2)); clear close_side
assert(numel(r)==counts(1) && numel(q)==counts(2));
file=[case_root '.shd']; fid=fopen(file,'rb'); assert(fid>=0); closer=onCleanup(@()fclose(fid));
record_words=fread(fid,1,'int32'); assert(isscalar(record_words) && record_words>0);
rb=4*record_words; title=deblank(fread(fid,80,'*char').');
fseek(fid,2*rb,'bof'); frequency=fread(fid,1,'float32'); sizes=fread(fid,6,'int32').';
attenuation=fread(fid,1,'float32');
assert(numel(sizes)==6 && sizes(1)==1 && sizes(4)==1);
assert(sizes(5)==numel(q) && sizes(6)==numel(r));
matches=[]; range_metadata=[]; depth_metadata=[];
for index=[8 9 10]
    fseek(fid,index*rb,'bof'); mr=fread(fid,numel(r),'float32');
    fseek(fid,(index-1)*rb,'bof'); mq=fread(fid,numel(q),'float32');
    if numel(mr)~=numel(r) || numel(mq)~=numel(q), continue; end
    range_ok=all(abs(r-mr)<=max(1e-9,2*double(eps(single(mr)))));
    depth_ok=all(abs(q-mq)<=max(1e-9,2*double(eps(single(mq)))));
    if range_ok && depth_ok
        matches(end+1)=index;
        range_metadata=mr; depth_metadata=mq;
    end
end
assert(isscalar(matches),'SHD records do not uniquely match both bound coordinate arrays.');
range_record=matches;
pressure=complex(zeros(numel(q),numel(r)));
for j=1:numel(q)
    fseek(fid,(range_record+j)*rb,'bof'); raw=fread(fid,2*numel(r),'float32');
    assert(numel(raw)==2*numel(r),'Incomplete SHD pressure row.');
    pressure(j,:)=raw(1:2:end)+1i*raw(2:2:end);
end
shade=struct('title',title,'frequency_hz',frequency,'attenuation',attenuation, ...
    'counts',sizes,'pressure',pressure,'receiver_range_metadata_m',range_metadata, ...
    'receiver_depth_metadata_m',depth_metadata,'receiver_range_m',r,'receiver_depth_m',q, ...
    'range_record',range_record,'record_bytes',rb,'file',file, ...
    'evaluation_coordinate_precision','float64 from bound r25d; SHD metadata float32 retained', ...
    'point3d_enabled',logical(flags(1)),'fast_q_enabled',logical(flags(2)), ...
    'exact_coordinate_sha256',black_file_hash(side));
end
