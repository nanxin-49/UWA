function [user,system]=bellhop_linux_memory()
t=fileread('/proc/meminfo'); available=str2double(regexp(t,'MemAvailable:\s*(\d+)','tokens','once')); total=str2double(regexp(t,'MemTotal:\s*(\d+)','tokens','once')); system=struct('PhysicalMemory',struct('Available',1024*available,'Total',1024*total)); user=struct();
end
