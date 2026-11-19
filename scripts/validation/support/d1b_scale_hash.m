function value = d1b_scale_hash(path)
%D1B_SCALE_HASH Bounded-memory SHA-256, including multi-GiB saved matrices.
fid = fopen(path,'rb'); assert(fid>0); cleaner = onCleanup(@() fclose(fid));
digest = java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid)
    block = fread(fid,1024*1024,'*uint8');
    if ~isempty(block), digest.update(typecast(block,'int8')); end
end
value = lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
end
