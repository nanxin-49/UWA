function closeout_s2b_json(path,value)
assert(~isfile(path),'Refusing overwrite');
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); cl=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(value,PrettyPrint=true));
end
