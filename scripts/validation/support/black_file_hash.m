function value = black_file_hash(filename, algorithm)
%BLACK_FILE_HASH Stream a named file without loading its contents into RAM.
if nargin < 2, algorithm = 'SHA-256'; end
digest = java.security.MessageDigest.getInstance(algorithm);
fid = fopen(filename, 'rb');
assert(fid >= 0, 'black:HashOpen', 'Cannot open %s', filename);
cleanup = onCleanup(@() fclose(fid));
while true
    bytes = fread(fid, 1024*1024, '*uint8');
    if isempty(bytes), break; end
    digest.update(typecast(bytes, 'int8'));
end
bytes = typecast(digest.digest(), 'uint8');
value = lower(reshape(dec2hex(bytes, 2).', 1, []));
end
