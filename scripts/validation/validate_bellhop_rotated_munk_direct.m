function receipt = validate_bellhop_rotated_munk_direct(options)
%VALIDATE_BELLHOP_ROTATED_MUNK_DIRECT Small Q-SSP axis time/convergence gate.
% Does not qualify rotated point spreading, reflections, total field, or Fig12b.
arguments
    options.bellhop_exe (1,:) char
    options.output_dir (1,:) char
    options.refinements (:,2) double {mustBePositive} = [10 2;5 1;2.5 .5]
    options.max_elapsed_s (1,1) double {mustBePositive} = 120
    options.minimum_available_bytes (1,1) double {mustBePositive} = 512*2^20
end
assert(isfile(options.bellhop_exe),'Missing user-selected Bellhop executable.');
assert(~isfolder(options.output_dir),'Use a new output directory.');
mkdir(options.output_dir); diary(fullfile(options.output_dir,'diary.log'));
diary_cleanup=onCleanup(@()diary('off'));
clock_start=tic; status='FAILED'; exit_code=1;
checks=table(strings(0,1),false(0,1),zeros(0,1),zeros(0,1), ...
    'VariableNames',{'check','pass','value','limit'});
configs={}; runs={}; err_text='';
[~,mem_before]=memory;
try
    assert(mem_before.PhysicalMemory.Available>=options.minimum_available_bytes,'RAM preflight failed.');
    directions={'up','down'};
    refinements=options.refinements;
    for id=1:2
        for ir=1:size(refinements,1)
            assert(toc(clock_start)<options.max_elapsed_s,'Elapsed time gate.');
            [~,mem_now]=memory;
            assert(mem_now.PhysicalMemory.Available>=options.minimum_available_bytes,'RAM runtime gate.');
            cfg=make_bellhop_rotated_munk_config(directions{id});
            cfg.ssp_spacing_m=refinements(ir,1); cfg.step_m=refinements(ir,2);
            cfg.beam_count=101; % centered axis ray; not a field beam-density claim
            tag=sprintf('%s_refine%d',directions{id},ir);
            case_dir=fullfile(options.output_dir,tag); mkdir(case_dir);
            case_root=fullfile(case_dir,'case');
            meta=write_bellhop_rotated_munk_q(case_root,cfg);
            proc=System.Diagnostics.Process(); proc.StartInfo.FileName=options.bellhop_exe;
            proc.StartInfo.Arguments='case'; proc.StartInfo.WorkingDirectory=case_dir;
            proc.StartInfo.UseShellExecute=false; proc.StartInfo.CreateNoWindow=true;
            proc.StartInfo.RedirectStandardOutput=true; proc.StartInfo.RedirectStandardError=true;
            remaining_ms=max(1,floor(1000*(options.max_elapsed_s-toc(clock_start))));
            assert(proc.Start(),'Bellhop did not start.');
            % Read async to avoid pipe-buffer deadlocks. Only this task-owned process is used.
            stdout_task=proc.StandardOutput.ReadToEndAsync();
            stderr_task=proc.StandardError.ReadToEndAsync();
            finished=proc.WaitForExit(int32(remaining_ms));
            if ~finished
                % Do not kill or retry a potentially live external process.
                local_text(fullfile(case_dir,'live_pid.txt'),num2str(proc.Id));
                error('rotated:LiveProcess','Time gate: task PID %d remains live; inspect before retry.',proc.Id);
            end
            proc.WaitForExit(); stdout=char(stdout_task.Result); stderr=char(stderr_task.Result);
            local_text(fullfile(case_dir,'stdout.log'),stdout);
            local_text(fullfile(case_dir,'stderr.log'),stderr);
            local_text(fullfile(case_dir,'exit_code.txt'),num2str(proc.ExitCode));
            assert(proc.ExitCode==0,'Bellhop exit code %d; inspect stderr/prt.',proc.ExitCode);
            prt=fileread([case_root '.prt']);
            assert(~contains(lower(prt),'fatal') && ~contains(lower(prt),'insufficient storage'), ...
                'Bellhop PRT reports a fatal/storage failure.');
            arr=read_bellhop_arrivals_ascii_vertical([case_root '.arr']);
            expected_z=cfg.source_z_m+(1-2*strcmp(cfg.direction,'up'))*cfg.receiver_range_m;
            tau=munk_vertical_travel_time_reference(expected_z,cfg.source_z_m,cfg.munk_parameters);
            delay=nan(size(tau)); counts=zeros(size(tau));
            for k=1:numel(tau)
                hit=abs(arr.receiver_range_m-cfg.receiver_range_m(k))<1e-3 & ...
                    abs(arr.source_angle_deg)<1e-5;
                counts(k)=nnz(hit);
                if counts(k)==1, delay(k)=arr.delay_s(hit); end
            end
            no_outer=all(arr.top_bounce_count==0 & arr.bottom_bounce_count==0);
            map_matrix=meta.physical_to_chart_matrix;
            add([tag '_proper_rotation'],abs(det(map_matrix)-1),1e-14);
            add([tag '_axis_arrival_count'],max(abs(counts-1)),0);
            add([tag '_outer_bounce_count'],double(~no_outer),0);
            dt=max(abs(delay-tau));
            add([tag '_axis_time_error_s'],dt,0.05/(2*pi*cfg.frequency_hz));
            np=meta.physical_sample_count;
            mid=.5*(meta.range_m(1:np-1)+meta.range_m(2:np));
            zmid=cfg.source_z_m+(1-2*strcmp(cfg.direction,'up'))*mid;
            c_linear=.5*(meta.c_mps(1:np-1)+meta.c_mps(2:np));
            c_exact=munk_profile_sound_speed(zmid,cfg.munk_parameters);
            interp_error=max(abs(c_linear-c_exact));
            configs{end+1}=cfg; %#ok<AGROW>
            runs{end+1}=struct('tag',tag,'meta',meta,'arrivals',arr,'axis_delay_s',delay, ...
                'reference_delay_s',tau,'time_error_s',dt,'ssp_midpoint_error_mps',interp_error); %#ok<AGROW>
            save(fullfile(case_dir,'result.mat'),'cfg','meta','arr','tau','delay','interp_error');
            writetable(arr,fullfile(case_dir,'arrivals.csv'));
            fprintf('%s time_error_s=%.9g SSP_midpoint_error_mps=%.9g\n',tag,dt,interp_error);
        end
    end
    status='COMPLETE'; exit_code=0;
catch exception
    err_text=getReport(exception,'extended','hyperlinks','off');
    fprintf(2,'%s\n',err_text);
end
writetable(checks,fullfile(options.output_dir,'checks.csv'));
local_text(fullfile(options.output_dir,'error.log'),err_text);
local_text(fullfile(options.output_dir,'exit_code.txt'),num2str(exit_code));
[~,mem_after]=memory;
unrun={'off-axis ray/path audit','beam-density and complex-field convergence','rotated 3D spreading', ...
    'surface/bottom/high-order reflections','independent FFP comparison','figures and figure audits'};
if ~strcmp(status,'COMPLETE'), unrun=[{'remaining direct-stage checks'},unrun]; end
receipt=struct('run_id',string(options.output_dir),'status',status,'elapsed_s',toc(clock_start), ...
    'failed_checks',nnz(~checks.pass)+double(~strcmp(status,'COMPLETE')),'unrun_checks',{unrun},'result_dir',options.output_dir, ...
    'matlab_release',version,'execution_mode','existing MATLAB MCP + local official Bellhop', ...
    'available_bytes_before',mem_before.PhysicalMemory.Available, ...
    'available_bytes_after',mem_after.PhysicalMemory.Available, ...
    'maximum_elapsed_s',options.max_elapsed_s,'minimum_available_bytes',options.minimum_available_bytes, ...
    'completed_steps',{cellfun(@(r)r.tag,runs,'UniformOutput',false)}, ...
    'next_allowed_step','Independent direct audit; then Q overlay reflection development. No large/remote run authorized.', ...
    'audit_status','pending independent replay; not physical qualification');
save(fullfile(options.output_dir,'result.mat'),'receipt','checks','configs','runs');
local_text(fullfile(options.output_dir,'checkpoint.json'),jsonencode(receipt,PrettyPrint=true));
disp(jsonencode(receipt));
    function add(name,value,limit)
        checks(end+1,:)={string(name),isfinite(value)&&value<=limit,value,limit};
    end
end
function local_text(path,content)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); c=onCleanup(@()fclose(fid));
fprintf(fid,'%s',content);
end