function h=closeout_s2b_hash(path)
md=java.security.MessageDigest.getInstance('SHA-256');
fid=fopen(path,'rb'); assert(fid>=0); cl=onCleanup(@()fclose(fid));
while ~feof(fid), md.update(fread(fid,1048576,'*uint8')); end
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
