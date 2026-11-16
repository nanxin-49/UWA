function out=bellhop_linux_path(in)
out=in; prefix='C:\Users\ASUS\.codex\worktrees\bellhop_rotated_fig12b\Explain'; replacement='/home/misc/cwx/bellhop_rotated_C_transfer_20261116_01/Explain';
if ischar(in)&&startsWith(in,prefix), out=[replacement strrep(in(numel(prefix)+1:end),'\','/')]; elseif isstring(in), out=string(arrayfun(@(v)bellhop_linux_path(char(v)),in,'UniformOutput',false)); end
end
