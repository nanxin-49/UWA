function output=run_bellhop_rotated_local_case(executable,case_root,max_seconds,options)
% Linux validation adapter; only this Process handle can be terminated.
arguments
 executable (1,:) char
 case_root (1,:) char
 max_seconds (1,1) double {mustBePositive}=20
 options.allow_task_termination (1,1) logical=false
 options.minimum_runtime_available_bytes (1,1) double {mustBeNonnegative}=0
end
assert(isfile(executable)&&isfile([case_root '.env'])); [~,ram]=bellhop_linux_memory(); assert(ram.PhysicalMemory.Available>=1.5*2^30,'Linux launch RAM gate');
reserve=max(options.minimum_runtime_available_bytes,.3*ram.PhysicalMemory.Available); [case_dir,stem]=fileparts(case_root); args=java.util.ArrayList(); args.add('/bin/bash'); args.add('-c'); args.add('printf ''%s\n'' "$$" > pid.txt; exec "$1" "$2"'); args.add('bellhop-task-owned'); args.add(executable); args.add(stem); builder=java.lang.ProcessBuilder(args); builder.directory(java.io.File(case_dir)); builder.redirectOutput(java.io.File(fullfile(case_dir,'stdout.log'))); builder.redirectError(java.io.File(fullfile(case_dir,'stderr.log')));
started=tic; proc=builder.start(); pid_path=fullfile(case_dir,'pid.txt'); while ~isfile(pid_path) && proc.isAlive() && toc(started)<2, pause(.01); end; if ~isfile(pid_path), proc.destroyForcibly(); proc.waitFor(); error('Task PID receipt not created'); end; pid=str2double(strtrim(fileread(pid_path))); assert(isfinite(pid)&&pid>0); if isfile(sprintf('/proc/%d/stat',pid)), write('process_identity.txt',fileread(sprintf('/proc/%d/stat',pid))); end
peak=0; minimum=ram.PhysicalMemory.Available; abort=''; last=-inf;
while proc.isAlive()
 if toc(started)-last>=.5
  [~,r]=bellhop_linux_memory(); minimum=min(minimum,r.PhysicalMemory.Available); last=toc(started); if minimum<reserve, abort='RUNTIME_RAM_RESERVE_GATE'; end
  status=sprintf('/proc/%d/status',pid); if isfile(status), v=regexp(fileread(status),'VmHWM:\s*(\d+)','tokens','once'); if ~isempty(v), peak=max(peak,1024*str2double(v{1})); end; end
 end
 if toc(started)>=min(120,max_seconds), abort='OBSERVATION_TIME_GATE'; end
 if ~isempty(abort)
  if options.allow_task_termination, proc.destroyForcibly(); proc.waitFor(); break; else, write('state.txt',['LIVE_' abort]); error('Task PID%d alive; no retry',pid); end
 end
 pause(.02);
end
exit_code=double(proc.waitFor()); elapsed=toc(started); stdout=fileread(fullfile(case_dir,'stdout.log')); stderr=fileread(fullfile(case_dir,'stderr.log')); write('exit_code.txt',num2str(exit_code)); write('state.txt','TERMINAL');
output=struct('pid',pid,'exit_code',exit_code,'elapsed_s',elapsed,'stdout',stdout,'stderr',stderr,'case_root',case_root,'peak_working_set_bytes',peak,'cpu_time_s',NaN,'abort_reason',abort,'minimum_observed_available_bytes',minimum,'runtime_reserve_bytes',reserve,'memory_counter_scope','Linux /proc VmHWM .5s samples; final tail excluded');
fid=fopen(fullfile(case_dir,'process_receipt.json'),'w'); fprintf(fid,'%s',jsonencode(output,PrettyPrint=true)); fclose(fid);
assert(isempty(abort),'Linux task-owned child stopped: %s',abort); assert(exit_code==0,'Nonzero Bellhop exit'); prt=lower(fileread([case_root '.prt'])); assert(~contains(prt,'fatal')&&~contains(prt,'insufficient storage'),'Bellhop PRT fatal/storage');
 function write(name,value), fid=fopen(fullfile(case_dir,name),'w','n','UTF-8'); assert(fid>=0); c=onCleanup(@()fclose(fid)); fprintf(fid,'%s',value); end
end
