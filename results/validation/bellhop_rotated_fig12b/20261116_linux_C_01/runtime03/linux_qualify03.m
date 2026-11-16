function receipt=linux_qualify03(project_root,run_dir)
% Fixed original101-beam,10981x3; retain original zero-difference gates.
addpath(fullfile(project_root,'scripts','validation','support')); addpath(fullfile(run_dir,'runtime03'),'-begin'); start=tic; diary(fullfile(run_dir,'qualification03_diary.log')); clean=onCleanup(@()diary('off'));
checks=table(strings(0,1),false(0,1),zeros(0,1),zeros(0,1),'VariableNames',{'check','pass','value','limit'}); process=[]; exception_text='';
baseline=fullfile(project_root,'results','validation','bellhop_rotated_fig12b','20261115_performance_review_01','timing_pair1_optimized','case'); out=fullfile(run_dir,'minimal03'); assert(~isfolder(out)); mkdir(out); root=fullfile(out,'case');
exe=fullfile(run_dir,'bin','bellhop_rotated_receiver64_2020');
try
 assert(license('test','MATLAB')==1); env=splitlines(string(fileread([baseline '.env']))); nonempty=find(strlength(strtrim(env))>0); assert(str2double(env(nonempty(end-2)))==101 && str2double(env(nonempty(end-7)))==3 && str2double(env(nonempty(end-5)))==10981);
 for ext={'.env','.ssp','.iwcov','.rotq','.r25d'}, copyfile([baseline ext{1}],[root ext{1}]); add(['input_SHA' ext{1}],double(~strcmpi(black_file_hash([root ext{1}]),black_file_hash([baseline ext{1}]))),0); end
 process=run_bellhop_rotated_local_case(exe,root,120,allow_task_termination=true,minimum_runtime_available_bytes=256*2^20);
 for ext={'.rr64input','.rr64eval','.rr64scale'}, add(['REAL8_actual_chain_exact_bytes' ext{1}],double(~strcmpi(black_file_hash([root ext{1}]),black_file_hash([baseline ext{1}]))),0); end
 a=read_bellhop_rotated_field(root); b=read_bellhop_rotated_field(baseline); add('exact_receiver_ranges',double(~isequal(a.receiver_range_m,b.receiver_range_m)),0); add('exact_receiver_transverse',double(~isequal(a.receiver_depth_m,b.receiver_depth_m)),0); add('all_pressure_finite',nnz(~isfinite(a.pressure)),0);
 add('original_complex_pressure_max_difference',max(abs(a.pressure-b.pressure),[],'all'),0); relative_L2=norm(a.pressure(:)-b.pressure(:))/norm(b.pressure(:)); add('original_complex_pressure_exact_values',nnz(a.pressure~=b.pressure),0);
 mw=readmatrix([root '.mwdiag'],'FileType','text','NumHeaderLines',1); rd=readmatrix([root '.r25diag'],'FileType','text','NumHeaderLines',1); add('all101_rays',abs(size(mw,1)-101)+abs(size(rd,1)-101),0); add('13_walls_7surface_6bottom',nnz(mw(:,2)~=13|mw(:,3)~=13|mw(:,4)~=7|mw(:,5)~=6),0); add('unit_tangent_original',max(rd(:,8)),1e-8); add('original_capacity_margin',max(mw(:,6)),.95*600000); add('stderr_empty',strlength(string(process.stderr)),0);
catch e
 exception_text=getReport(e,'extended','hyperlinks','off'); relative_L2=NaN;
end
elapsed=toc(start); build=0; job_elapsed=posixtime(datetime('now','TimeZone','UTC'))-str2double(getenv('BL_JOB_START_EPOCH')); used=build+job_elapsed;
status='QUALIFICATION_FAILED_C_NOT_STARTED'; if isempty(exception_text)&&all(checks.pass), status='QUALIFICATION_PASS_C9_ALLOWED'; end
receipt=struct('status',status,'failed_checks',nnz(~checks.pass)+~isempty(exception_text),'checks_count',height(checks),'elapsed_s',elapsed,'process',process,'relative_L2',relative_L2,'original_pressure_tolerance',0,'original_coordinate_tolerance',0,'matlab_version',version,'executable_sha256',black_file_hash(exe),'build_elapsed_s',build,'MATLAB_job_wall_s',job_elapsed,'prior_budget_remaining_s',726.081105536377,'charged_this_round_s',used,'remaining_original_budget_s',726.081105536377-used,'C_completed',30,'new_formal_C',0,'baseline',baseline,'root',root,'error',exception_text);
writetable(checks,fullfile(run_dir,'qualification03_checks.csv')); save(fullfile(run_dir,'qualification03.mat'),'receipt','checks'); fid=fopen(fullfile(run_dir,'checkpoint_qualification03.json'),'w','n','UTF-8'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid); fprintf('LINUX_QUAL %s failed%d relL2%.9g remaining%.3fs\n',status,receipt.failed_checks,relative_L2,receipt.remaining_original_budget_s); disp(checks(~checks.pass,:)); if ~isempty(exception_text), fprintf('%s\n',exception_text); end
 function add(name,value,limit)
  checks(end+1,:)={string(name),isfinite(value)&&value<=limit,value,limit};
 end
end
