function out=bellhop_linux_map(in)
out=in;
if isstruct(in)
 names=fieldnames(in); for j=1:numel(in), for k=1:numel(names), out(j).(names{k})=bellhop_linux_map(in(j).(names{k})); end; end
elseif iscell(in)
 for k=1:numel(in), out{k}=bellhop_linux_map(in{k}); end
elseif ischar(in)||isstring(in)
 out=bellhop_linux_path(in);
end
end
