function receipt = validate_bellhop_rotated_range_metadata(executable,source_checkpoint,output_dir)
%VALIDATE_BELLHOP_ROTATED_RANGE_METADATA Actual13-wall exact-coordinate pilot.
% Same physical r25d samples, different strictly ordered ENV float labels.
% Two101-beam/3-column cases, <=60s cumulative, <=20s native process each.
arguments
    executable (1,:) char
    source_checkpoint (1,:) char
    output_dir (1,:) char
end
assert(~isfolder(output_dir)); mkdir(output_dir); t=tic;
s=jsondecode(fileread(source_checkpoint)); old=s.cfg;
cfg=make_bellhop_rotated_multifield_config('up',13,old.nominal_receiver_depths_m(:).',[-1 0 1]);
cfg.ssp_spacing_m=old.ssp_spacing_m; cfg.step_m=old.step_m; cfg.beam_count=101;
cfg.angle_limits_deg=old.angle_limits_deg(:).';
checks=table(strings(0,1),false(0,1),zeros(0,1),zeros(0,1),'VariableNames',{'check','pass','value','limit'});
fields=cell(1,2); processes=cell(1,2); status='FAILED'; failure='';
try
for j=1:2
    [~,mem]=memory; assert(mem.PhysicalMemory.Available>=512*2^20 && toc(t)<60,'Pilot resource gate.');
    folder=fullfile(output_dir,sprintf('labels%d',j)); mkdir(folder); stem=fullfile(folder,'case');
    write_bellhop_rotated_multifield_q(stem,cfg);
    repaired=jsondecode(fileread([stem '.range_metadata.json']));
    add(sprintf('case%d_strict_native_metadata',j),repaired.duplicate_count_after,0);
    add(sprintf('case%d_existing_metadata_ULP_bound',j),repaired.maximum_error_over_ulp,2);
    if j==2
        r=cfg.receiver_range_m(:); km=single((r+.5*double(eps(single(r))))/1000);
        label=single(km.*single(1000));
        for n=2:numel(r)
            while label(n)<=label(n-1)
                km(n)=typecast(typecast(km(n),'uint32')+uint32(1),'single'); label(n)=single(km(n)*single(1000));
            end
        end
        assert(all(diff(label)>0) && all(abs(double(label)-r)<=2*double(eps(label))));
        text=splitlines(string(fileread([stem '.env']))); nonempty=find(strlength(strtrim(text))>0);
        row=nonempty(end-4); text(row)=string([sprintf('%.17g ',double(km)) '/']);
        fid=fopen([stem '.env'],'w'); fprintf(fid,'%s',join(text,newline)); fclose(fid);
        save([stem '.alternate_metadata.mat'],'r','km','label');
    end
    processes{j}=run_bellhop_rotated_local_case(executable,stem,min(20,60-toc(t)));
    fields{j}=read_bellhop_rotated_field(stem);
    add(sprintf('case%d_actual_range_binding',j),max(abs(fields{j}.receiver_range_m-cfg.receiver_range_m(:))),0);
    add(sprintf('case%d_raw_metadata_strict',j),double(any(diff(fields{j}.receiver_range_metadata_m)<=0)),0);
    add(sprintf('case%d_raw_metadata_existing_bound',j), ...
        max(abs(fields{j}.receiver_range_metadata_m-cfg.receiver_range_m(:))./double(eps(single(fields{j}.receiver_range_metadata_m)))),2);
end
add('different_metadata_actually_used',double(isequal(fields{1}.receiver_range_metadata_m,fields{2}.receiver_range_metadata_m)),0);
add('physical_coordinates_identical',max(abs(fields{1}.receiver_range_m-fields{2}.receiver_range_m)),0);
add('complex_pressure_identical',max(abs(fields{1}.pressure-fields{2}.pressure),[],'all'),0);
add('source_binary_unchanged',double(~strcmpi(black_file_hash(executable),s.executable_sha256)),0);
status='COMPLETE';
catch exception
failure=getReport(exception,'extended','hyperlinks','off');
end
receipt=struct('status',status,'failed_checks',nnz(~checks.pass)+~strcmp(status,'COMPLETE'), ...
    'elapsed_s',toc(t),'output_dir',output_dir,'scope','Native13-wall input-label monotonicity and exact pressure invariance only; no complete field or physical FFP qualification', ...
    'executable_sha256',black_file_hash(executable),'source_checkpoint',source_checkpoint);
writetable(checks,fullfile(output_dir,'checks.csv'));
save(fullfile(output_dir,'result.mat'),'receipt','cfg','checks','fields','processes');
write(fullfile(output_dir,'checkpoint.json'),jsonencode(receipt,PrettyPrint=true));
write(fullfile(output_dir,'error.log'),failure); disp(jsonencode(receipt)); disp(checks(~checks.pass,:)); if ~isempty(failure), fprintf('%s\n',failure); end
    function add(name,value,limit)
        checks(end+1,:)={string(name),isfinite(value)&&value<=limit,value,limit};
    end
end
function write(path,content)
fid=fopen(path,'w'); assert(fid>=0); cleaner=onCleanup(@()fclose(fid)); fprintf(fid,'%s',content);
end