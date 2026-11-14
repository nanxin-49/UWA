function output = run_bellhop_rotated_local_case(executable,case_root,max_seconds)
%RUN_BELLHOP_ROTATED_LOCAL_CASE Run a single task-owned external Bellhop case.
% The caller writes inputs. This helper never changes MATLAB's current folder.
arguments
    executable (1,:) char
    case_root (1,:) char
    max_seconds (1,1) double {mustBePositive} = 20
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
if ~proc.WaitForExit(int32(1000*max_seconds))
    text_file('state.txt','LIVE_ON_OBSERVATION_TIMEOUT: inspect recorded PID before retry; no automatic kill');
    error('rotated:LiveProcess','Bellhop task PID %d is still live in %s.',pid,case_dir);
end
proc.WaitForExit();
stdout=char(stdout_task.Result); stderr=char(stderr_task.Result);
text_file('stdout.log',stdout); text_file('stderr.log',stderr);
text_file('exit_code.txt',num2str(proc.ExitCode));
text_file('state.txt','TERMINAL');
output=struct('pid',pid,'exit_code',proc.ExitCode,'elapsed_s',toc(started), ...
    'stdout',stdout,'stderr',stderr,'case_root',case_root, ...
    'peak_working_set_bytes',double(proc.PeakWorkingSet64), ...
    'cpu_time_s',double(proc.TotalProcessorTime.TotalSeconds));
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