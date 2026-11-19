function output=run_bellhop_rotated_local_case(executable,case_root,max_seconds,options)
% Task-owned Linux adapter: identity is saved by the shell BEFORE exec.
arguments
 executable (1,:) char
 case_root (1,:) char
 max_seconds (1,1) double {mustBePositive}=120
 options.allow_task_termination (1,1) logical=false
 options.minimum_runtime_available_bytes (1,1) double=256*2^20
end
assert(isfile(executable)&&isfile([case_root '.iwcov']));
[~,ram]=bellhop_linux_memory(); available=ram.PhysicalMemory.Available;
assert(available>=1.5*2^30,'Launch RAM gate');
reserve=max(options.minimum_runtime_available_bytes,.3*available);
[case_dir,stem]=fileparts(case_root); args=java.util.ArrayList();
args.add('/bin/bash'); args.add(fullfile(fileparts(mfilename('fullpath')),'native_identity_exec.sh')); args.add(executable); args.add(stem);
builder=java.lang.ProcessBuilder(args); builder.directory(java.io.File(case_dir));
builder.redirectOutput(java.io.File(fullfile(case_dir,'stdout.log'))); builder.redirectError(java.io.File(fullfile(case_dir,'stderr.log')));
started=tic; proc=builder.start(); minimum=available; samples=0; abort='';
while proc.isAlive()
 [~,ram]=bellhop_linux_memory(); minimum=min(minimum,ram.PhysicalMemory.Available); samples=samples+1;
 if ram.PhysicalMemory.Available<reserve, abort='RUNTIME_RAM_RESERVE_GATE'; end
 if toc(started)>=min(120,max_seconds), abort='NATIVE_TIME_GATE'; end
 if ~isempty(abort)
  assert(options.allow_task_termination,'No permission to terminate owned Process'); proc.destroyForcibly(); break;
 end
 pause(.02);
end
% Always reap and save the real exit code before optional diagnostics.
exit_code=double(proc.waitFor()); elapsed=toc(started);
fid=fopen(fullfile(case_dir,'exit_code.txt'),'w'); fprintf(fid,'%d\n',exit_code); fclose(fid);
pid=NaN; if isfile(fullfile(case_dir,'pid.txt')), pid=str2double(strtrim(fileread(fullfile(case_dir,'pid.txt')))); end
[~,ram]=bellhop_linux_memory(); minimum=min(minimum,ram.PhysicalMemory.Available);
output=struct('pid',pid,'exit_code',exit_code,'elapsed_s',elapsed,'stdout',fileread(fullfile(case_dir,'stdout.log')),'stderr',fileread(fullfile(case_dir,'stderr.log')),'case_root',case_root,'abort_reason',abort,'minimum_observed_available_bytes',minimum,'runtime_reserve_bytes',reserve,'launch_available_bytes',available,'live_RAM_samples',samples,'identity_before_exec',isfile(fullfile(case_dir,'process_identity.txt')),'memory_scope','MemAvailable before/20ms live/after; sampled, not continuous');
fid=fopen(fullfile(case_dir,'process_receipt.json'),'w'); fprintf(fid,'%s',jsonencode(output,PrettyPrint=true)); fclose(fid);
assert(isempty(abort),'Native resource/time failure'); assert(exit_code==0,'Native Bellhop nonzero exit');
assert(output.identity_before_exec&&isfinite(pid),'Native identity receipt missing');
prt=lower(fileread([case_root '.prt'])); assert(~contains(prt,'fatal')&&~contains(prt,'insufficient storage'),'PRT fatal');
assert(isfile([case_root '.shd']),'C mode SHD missing');
end