function d1b_scale_json(path,value)
%D1B_SCALE_JSON Atomically replace a task-owned JSON receipt.
path = string(path); tmp = path+".pending";
fid = fopen(tmp,'w','n','UTF-8'); assert(fid>0);
cleaner = onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true)); clear cleaner
assert(movefile(tmp,path,'f'));
end
