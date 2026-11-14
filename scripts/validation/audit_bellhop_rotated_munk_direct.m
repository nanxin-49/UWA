function receipt = audit_bellhop_rotated_munk_direct(input_dirs,output_dir)
%AUDIT_BELLHOP_ROTATED_MUNK_DIRECT Independently replay saved axis evidence.
% Integrates the actual written piecewise-linear Q table analytically, and a
% separately expressed continuous Munk formula; also checks raw ARR/CSV binding.
arguments
    input_dirs (1,:) cell
    output_dir (1,:) char
end
assert(~isfolder(output_dir)); mkdir(output_dir);
diary(fullfile(output_dir,'diary.log')); cleanup=onCleanup(@()diary('off'));
t=tic; checks=table(strings(0,1),false(0,1),zeros(0,1),zeros(0,1), ...
    'VariableNames',{'check','pass','value','limit'});
bindings=struct('path',{},'sha256',{}); details={}; error_text=''; status='FAILED';
try
[~,mem]=memory; assert(mem.PhysicalMemory.Available>=512*2^20,'RAM preflight gate.');
for j=1:numel(input_dirs)
    receipt_path=fullfile(input_dirs{j},'checkpoint.json');
    cp=jsondecode(fileread(receipt_path));
    assert(strcmp(cp.status,'COMPLETE'),'Input run is not complete.');
    bind(receipt_path);
    files=dir(fullfile(input_dirs{j},'**','arrivals.csv'));
    assert(~isempty(files),'No raw axis ARR cases in %s.',input_dirs{j});
    for k=1:numel(files)
        folder=files(k).folder; stem=fullfile(folder,'case');
        names={files(k).name,'result.mat','case.env','case.ssp','case.arr','case.prt','exit_code.txt','stdout.log','stderr.log'};
        for n=1:numel(names), bind(fullfile(folder,names{n})); end
        data=load(fullfile(folder,'result.mat'),'cfg','arr','tau','delay');
        cfg=data.cfg;
        a=read_bellhop_arrivals_ascii_vertical([stem '.arr']);
        csv=readtable(fullfile(folder,'arrivals.csv'));
        assert(isequal(a.Properties.VariableNames,csv.Properties.VariableNames));
        % CSV is a second serialization; compare every retained real/complex
        % field separately (complex table columns are represented as text).
        numeric_names={'receiver_range_m','receiver_depth_m','delay_s','source_angle_deg', ...
            'receiver_angle_deg','top_bounce_count','bottom_bounce_count','amplitude','phase_deg'};
        difference=0;
        for n=1:numel(numeric_names)
            v=numeric_names{n}; difference=max(difference,max(abs(a.(v)-csv.(v))));
        end
        tag=sprintf('run%d_case%d',j,k);
        add([tag '_raw_csv_difference'],difference,1e-12);
        fid=fopen([stem '.ssp'],'r'); assert(fid>=0); closer=onCleanup(@()fclose(fid));
        nr=fscanf(fid,'%d',1); ranges=1000*fscanf(fid,'%f',nr);
        c_rows=fscanf(fid,'%f',[nr 2]); clear closer
        assert(size(c_rows,2)==2 && size(c_rows,1)==nr);
        add([tag '_transverse_ssp_difference_mps'],max(abs(diff(c_rows,1,2))),1e-12);
        tau_q=nan(size(cfg.receiver_range_m)); delay=tau_q; tau_munk=tau_q;
        for n=1:numel(delay)
            r=cfg.receiver_range_m(n);
            use=ranges<r; rr=[ranges(use);r];
            cc=interp1(ranges,c_rows(:,1),rr,'linear');
            u=diff(cc)./cc(1:end-1); ratio=ones(size(u));
            nz=u~=0; ratio(nz)=log1p(u(nz))./u(nz);
            tau_q(n)=sum(diff(rr)./cc(1:end-1).*ratio);
            pp=cfg.munk_parameters;
            sign_z=1-2*strcmp(cfg.direction,'up');
            z=cfg.source_z_m+sign_z*r;
            inverse_c=@(zz)1./(pp.c0_mps*(1+pp.epsilon*(2*(zz-pp.axis_depth_m)/pp.scale_m-1+ ...
                exp(-2*(zz-pp.axis_depth_m)/pp.scale_m))));
            tau_munk(n)=abs(integral(inverse_c,cfg.source_z_m,z,'RelTol',1e-12,'AbsTol',1e-13));
            hit=abs(a.receiver_range_m-r)<1e-3 & abs(a.source_angle_deg)<1e-5;
            assert(nnz(hit)==1,'Expected one centered ray arrival.');
            delay(n)=a.delay_s(hit);
        end
        add([tag '_continuous_reference_replay_s'],max(abs(tau_munk-data.tau)),1e-12);
        add([tag '_axis_delay_binding_s'],max(abs(delay-data.delay)),1e-15);
        add([tag '_native_to_Q_time_s'],max(abs(delay-tau_q)),.05/(2*pi*cfg.frequency_hz));
        add([tag '_native_to_Munk_time_s'],max(abs(delay-tau_munk)),.05/(2*pi*cfg.frequency_hz));
        add([tag '_artificial_bounces'],max(a.top_bounce_count+a.bottom_bounce_count),0);
        details{end+1}=struct('case_dir',folder,'delay_s',delay,'tau_Q_exact_s',tau_q, ...
            'tau_Munk_s',tau_munk,'cfg',cfg);
    end
end
status='COMPLETE';
catch exception
    error_text=getReport(exception,'extended','hyperlinks','off');
end
writetable(checks,fullfile(output_dir,'checks.csv'));
receipt=struct('status',status,'elapsed_s',toc(t),'failed_checks',nnz(~checks.pass)+~strcmp(status,'COMPLETE'), ...
    'execution_mode','existing MATLAB MCP; saved data only','matlab_release',version, ...
    'qualification','Direct axis Q/Munk time and serialization only; no total-field qualification', ...
    'result_dir',output_dir);
save(fullfile(output_dir,'result.mat'),'checks','details','receipt','bindings');
text(fullfile(output_dir,'input_sha256.json'),jsonencode(bindings,PrettyPrint=true));
text(fullfile(output_dir,'error.log'),error_text);
text(fullfile(output_dir,'checkpoint.json'),jsonencode(receipt,PrettyPrint=true));
disp(jsonencode(receipt));
    function bind(path)
        bindings(end+1)=struct('path',path,'sha256',black_file_hash(path));
    end
    function add(name,value,limit)
        checks(end+1,:)={string(name),isfinite(value)&&value<=limit,value,limit};
    end
end
function text(path,value)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',value);
end