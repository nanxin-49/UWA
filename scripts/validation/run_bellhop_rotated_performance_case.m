function receipt=run_bellhop_rotated_performance_case(executable,frozen_root,output_dir)
%RUN_BELLHOP_ROTATED_PERFORMANCE_CASE One fixed101beam/3column13th-order case.
arguments
    executable (1,:) char
    frozen_root (1,:) char
    output_dir (1,:) char
end
assert(~isfolder(output_dir)); mkdir(output_dir); root=fullfile(output_dir,'case'); t=tic;
env=splitlines(string(fileread([frozen_root '.env']))); nonempty=find(strlength(strtrim(env))>0);
assert(str2double(env(nonempty(end-2)))==101 && str2double(env(nonempty(end-7)))==3);
assert(str2double(env(nonempty(end-5)))==10981);
fid=fopen([frozen_root '.iwcov'],'r'); assert(fid>=0); mode=fscanf(fid,'%d',1); fclose(fid); assert(mode==5);
[~,mem]=memory; assert(mem.PhysicalMemory.Available>=1.5*2^30,'RAM launch gate');
versions=struct('path',{},'sha256',{});
for ext={'.env','.ssp','.iwcov','.rotq','.r25d'}
    source=[frozen_root ext{1}]; destination=[root ext{1}]; copyfile(source,destination);
    versions(end+1)=struct('path',destination,'sha256',black_file_hash(source)); %#ok<AGROW>
    assert(strcmpi(versions(end).sha256,black_file_hash(destination)));
end
status='FAILED'; failure=''; process=[]; checks=table(strings(0,1),false(0,1),zeros(0,1),zeros(0,1),'VariableNames',{'check','pass','value','limit'});
try
    process=run_bellhop_rotated_local_case(executable,root,20,allow_task_termination=true,minimum_runtime_available_bytes=256*2^20);
    for ext={'.shd','.rr64input','.rr64eval','.rr64scale','.mwdiag','.r25diag'}
        add(['exact_saved_bytes' ext{1}],double(~strcmpi(black_file_hash([root ext{1}]),black_file_hash([frozen_root ext{1}]))),0);
    end
    field=read_bellhop_rotated_field(root); old=read_bellhop_rotated_field(frozen_root);
    add('finite_complex_pressure',nnz(~isfinite(field.pressure)),0);
    add('complex_pressure_max_difference',max(abs(field.pressure-old.pressure),[],'all'),0);
    add('stderr_empty',strlength(string(process.stderr)),0);
    status='COMPLETE';
catch exception
    failure=getReport(exception,'extended','hyperlinks','off');
end
receipt=struct('status',status,'failed_checks',nnz(~checks.pass)+~strcmp(status,'COMPLETE'),'elapsed_s',toc(t), ...
    'process',process,'frozen_root',frozen_root,'root',root,'input_sha256',versions,'executable',executable, ...
    'executable_sha256',black_file_hash(executable),'scope','Fixed order13,101beams,10981x3; no complete A or14-16/B/C/D', ...
    'matlab_release',version,'next_step','Read profiling/audit evidence before another bounded comparison');
writetable(checks,fullfile(output_dir,'checks.csv'));
save(fullfile(output_dir,'result.mat'),'receipt','checks');
write(fullfile(output_dir,'checkpoint.json'),jsonencode(receipt,PrettyPrint=true)); write(fullfile(output_dir,'error.log'),failure);
disp(jsonencode(struct('status',status,'failed_checks',receipt.failed_checks,'process',process,'output_dir',output_dir))); disp(checks(~checks.pass,:));
if ~isempty(failure), fprintf('%s\n',failure); end
    function add(name,value,limit)
        checks(end+1,:)={string(name),isfinite(value)&&value<=limit,value,limit};
    end
end
function write(path,text)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); c=onCleanup(@()fclose(fid)); fprintf(fid,'%s',text);
end