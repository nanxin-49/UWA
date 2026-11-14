function output = run_bellhop_rotated_local_case(executable,case_root,max_seconds,options)
%RUN_BELLHOP_ROTATED_LOCAL_CASE Run a single task-owned external Bellhop case.
% The caller writes inputs. This helper never changes MATLAB's current folder.
arguments
    executable (1,:) char
    case_root (1,:) char
    max_seconds (1,1) double {mustBePositive} = 20
    options.allow_task_termination (1,1) logical=false
    options.minimum_runtime_available_bytes (1,1) double {mustBeNonnegative}=0
end
assert(isfile(executable) && isfile([case_root '.env']));
case_dir=fileparts(case_root); [~,stem]=fileparts(case_root);
proc=System.Diagnostics.Process(); proc.StartInfo.FileName=executable;
proc.StartInfo.Arguments=stem; proc.StartInfo.WorkingDirectory=case_dir;
proc.StartInfo.UseShellExecute=false; proc.StartInfo.CreateNoWindow=true;
proc.StartInfo.RedirectStandardOutput=true; proc.StartInfo.RedirectStandardError=true;
started=tic; assert(proc.Start(),'Bellhop did not start.');
stdout_task=proc.StandardOutput.ReadToEndAsync(); stderr_task=proc.StandardError.ReadToEndAsync();
pid=proc.Id; text_file('pid.txt',num2str(pid));
peak_working_set_bytes=NaN; cpu_time_s=NaN; abort_reason=''; last_ram_check=-inf; minimum_observed_available_bytes=inf;
while true
    % Windows may refuse memory queries after exit; sample while task is live.
    try
        proc.Refresh();
        peak_working_set_bytes=max(peak_working_set_bytes,double(proc.PeakWorkingSet64),'omitnan');
        cpu_time_s=double(proc.TotalProcessorTime.TotalSeconds);
    catch
        % A process may exit between Refresh and the counter query.
    end
    if proc.WaitForExit(int32(20)), break; end
    if options.minimum_runtime_available_bytes>0 && toc(started)-last_ram_check>=.5
        [~,ram]=memory; last_ram_check=toc(started);
        minimum_observed_available_bytes=min(minimum_observed_available_bytes,ram.PhysicalMemory.Available);
        if ram.PhysicalMemory.Available<options.minimum_runtime_available_bytes, abort_reason='RUNTIME_RAM_GATE'; end
    end
    if toc(started)>=max_seconds, abort_reason='OBSERVATION_TIME_GATE'; end
    if ~isempty(abort_reason)
        if options.allow_task_termination
            if ~proc.HasExited, proc.Kill(); end
            proc.WaitForExit(); break
        else
            text_file('state.txt',['LIVE_ON_' abort_reason ': inspect recorded PID before retry; no automatic kill']);
            error('rotated:LiveProcess','Bellhop task PID %d is still live in %s.',pid,case_dir);
        end
    end
end
proc.WaitForExit();
stdout=char(stdout_task.Result); stderr=char(stderr_task.Result);
text_file('stdout.log',stdout); text_file('stderr.log',stderr);
text_file('exit_code.txt',num2str(proc.ExitCode));
if isempty(abort_reason), text_file('state.txt','TERMINAL'); else, text_file('state.txt',['TERMINAL_TASK_ABORT_' abort_reason]); end
output=struct('pid',pid,'exit_code',proc.ExitCode,'elapsed_s',toc(started), ...
    'stdout',stdout,'stderr',stderr,'case_root',case_root, ...
    'peak_working_set_bytes',peak_working_set_bytes, ...
    'cpu_time_s',cpu_time_s, ...
    'memory_counter_scope','20ms live-process peak counter samples; final unobserved tail excluded', ...
    'abort_reason',abort_reason,'minimum_observed_available_bytes',minimum_observed_available_bytes);
assert(isempty(abort_reason),'Task-owned Bellhop stopped at %s in %s.',abort_reason,case_dir);
assert(proc.ExitCode==0,'Nonzero Bellhop exit code; inspect %s.',case_dir);
prt=fileread([case_root '.prt']);
assert(~contains(lower(prt),'fatal') && ~contains(lower(prt),'insufficient storage'), ...
    'Bellhop fatal/storage failure; inspect %s.',case_dir);
    function text_file(name,value)
        fid=fopen(fullfile(case_dir,name),'w','n','UTF-8'); assert(fid>=0);
        close_file=onCleanup(@()fclose(fid));
        fprintf(fid,'%s',value);
    end
end