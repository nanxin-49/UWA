function receipt=audit_m1_convergence_vertical(run_dir,mode)
%AUDIT_M1_CONVERGENCE_VERTICAL Independently reopen M1 MAT/CSV and identity.
% Does not call production propagation or the one-step reference helper.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
assert(java.io.File(run_dir).isAbsolute());
assert(startsWith(char(java.io.File(run_dir).getCanonicalPath()), ...
    [fullfile(root,'results','validation','m1_convergence') filesep]));
out_dir=fullfile(run_dir,mode); assert(~isfolder(out_dir),'Audit directory already exists.'); mkdir(out_dir);
audit_kind=regexprep(mode,'_v\d+$','');
t=tic; threads=maxNumCompThreads(1); restore=onCleanup(@()maxNumCompThreads(threads)); %#ok<NASGU>
diary(fullfile(out_dir,'stdout.log')); dc=onCleanup(@()diary('off')); %#ok<NASGU>
checks=struct('name',{},'value',{},'limit',{},'pass',{}); rows=struct([]);
budget=jsondecode(fileread(fullfile(run_dir,'budget.json')));
receipt=struct('status','FAILED','elapsed_s',0,'failed_checks',1,'matlab_pid',feature('getpid'), ...
    'matlab_version',version,'exit_code',1,'exit_code_kind','job exception status, not process exit', ...
    'audit_code_sha256',sha([mfilename('fullpath') '.m']),'result_dir',out_dir);
try
    if strcmp(audit_kind,'stage1_audit')
        ids={'P4','P6','P8','Z4','Z8','A-C4','A-C6','A-C8','A-W4','A-W6','A-W8','A-G4','A-G6','A-G8','A-R4'};
    elseif strcmp(audit_kind,'reference_audit')
        ids={'P4','A-C4','A-C6','A-C8','A-W4','A-W6','A-W8','A-G4','A-G6','A-G8','A-R4'};
    elseif strcmp(audit_kind,'stage2_audit')
        ids={'P4','P6','P8','SP4','SP6','SP8','A-C4','A-C6','A-C8'};
    else
        assert(startsWith(mode,'audit_')); ids={extractAfter(mode,6)};
    end
    for j=1:numel(ids)
        id=char(ids{j}); job_dir=fullfile(run_dir,id);
        state=jsondecode(fileread(fullfile(job_dir,'checkpoint.json')));
        for k=1:numel(state.source_sha256)
            src=state.source_sha256(k);
            add([id '_source_' num2str(k)],~strcmpi(src.sha256,sha(fullfile(root,src.path))),0);
        end
        for k=1:numel(state.artifacts_sha256)
            a=state.artifacts_sha256(k);
            add([id '_artifact_' num2str(k)],~strcmpi(a.sha256,sha(fullfile(job_dir,a.path))),0);
        end
        add([id '_job_exit'],state.exit_code,0);
        add([id '_stdout_captured'],~isfield(state,'stdout_capture'),0);
        add([id '_single_thread'],abs(state.threads-1),0);
        add([id '_runtime_s'],state.elapsed_s,3600);
        q=load(fullfile(job_dir,'payload.mat'),'payload'); q=q.payload;
        cfg=q.config;
        [u,ram]=memory; %#ok<ASGLU>
        receipt.min_available_bytes=ram.PhysicalMemory.Available;
        add([id '_coordinate_grid'],max(abs(q.coords_m/.1953125-round(q.coords_m/.1953125)),[],'all'),0);
        add([id '_physical_component_closure'],max(abs(q.physical_total-q.physical_direct-q.physical_reflect)),1e-12);
        k0=2*pi*cfg.f_hz/1500;
        expected=conj(q.physical_total)*exp(1i*k0*97);
        add([id '_v2_direct_scaled'],max(abs(q.dsp_total-expected))/max(abs(q.physical_direct)),1e-12);
        options=detectImportOptions(fullfile(job_dir,'fields_ieee754.csv'),'TextType','string');
        options=setvartype(options,options.VariableNames,'string');
        csv=readtable(fullfile(job_dir,'fields_ieee754.csv'),options);
        a=csv{:,:}; numbers=zeros(size(a));
        for k=1:6, numbers(:,k)=hex2num(char(a(:,k))); end
        expected=[q.coords_m,real(q.native_direct),imag(q.native_direct),real(q.native_reflect),imag(q.native_reflect)];
        add([id '_MAT_CSV_bitwise'],max(abs(numbers-expected),[],'all'),0);
        if strcmp(cfg.kind,'AS')
            estimate=cfg.estimated_extra_peak_bytes;
            assert(estimate<=.7*ram.PhysicalMemory.Available,'Audit RAM resource gate.');
            if cfg.n>=1966, assert(ram.PhysicalMemory.Available>=1.5*2^30); end
            x=(-cfg.n/2:cfg.n/2-1)*cfg.dx_m; [X,Y]=meshgrid(x,x);
            k=(2*pi/cfg.width_m)*[0:cfg.n/2-1,-cfg.n/2:-1];
            % Preserve the stated left-to-right k0^2-kx^2-ky^2 evaluation.
            % Reassociation changes FFT roundoff at noise-floor boundaries.
            kz=sqrt(complex((k0^2-k.^2)-k(:).^2,0));
            spectrum=fft2(exp(-(X.^2+Y.^2)/.18));
            [~,ix]=ismember(q.coords_m(:,1),x); [~,iy]=ismember(q.coords_m(:,2),x);
            index=sub2ind([cfg.n,cfg.n],iy,ix);
            targets={q.native_direct,q.native_incident,-q.native_reflect};
            L=[97 100 103]; edges=zeros(3,2);
            edge=abs(X)>=.45*cfg.width_m | abs(Y)>=.45*cfg.width_m;
            bd=false(cfg.n); bd([1 end],:)=true; bd(:,[1 end])=true;
            for k=1:3
                field=ifft2(spectrum.*exp(1i*L(k)*(kz-k0)));
                add([id '_independent_AS_L' num2str(L(k))],norm(field(index)-targets{k})/norm(targets{k}),1e-10);
                power=abs(field).^2;
                edges(k,:)=[sum(power(edge))/sum(power(:)),20*log10(max(abs(field(bd)))/max(abs(field(:))))];
                add([id '_edge_L' num2str(L(k))],edges(k,1),1e-4,true);
                add([id '_boundary_db_L' num2str(L(k))],edges(k,2),-40);
                clear field power
            end
            add([id '_edge_table_reconstruction'],max(abs(edges-[q.edge_table.edge5_fraction,q.edge_table.boundary_db]),[],'all'),1e-9);
            clear X Y kz spectrum edge bd
        else
            paths={'direct','incident','reflect'};
            for k=1:3
                name=paths{k}; tr=q.trace.(name);
                lengths=[97 100 3]; count=ceil(lengths(k)/(cfg.stepz_lamb*1500/cfg.f_hz));
                add([id '_' name '_step_coverage'],~isequal(tr(:,1),(0:count).'),0);
                add([id '_' name '_edge_max'],max(tr(:,3)),1e-4,true);
                add([id '_' name '_boundary_max_db'],max(tr(:,4)),-40);
                table=readtable(fullfile(job_dir,[name '_every_step.csv']));
                finite=isfinite(tr); from_csv=table{:,:};
                add([id '_' name '_trace_csv'],max(abs(from_csv(finite)-tr(finite))),1e-9);
            end
            out=q.output;
            c=[out.h_total-out.H_f(out.idx_f_ref);out.h_direct-out.H_direct_f(out.idx_f_ref); ...
                out.h_reflect-out.H_reflect_f(out.idx_f_ref);out.H_f-out.H_direct_f-out.H_reflect_f];
            add([id '_public_scalar_closure'],max(abs(c)),1e-12);
            ix=find(all(q.coords_m==0,2));
            add([id '_scalar_roi_direct'],abs(q.native_direct(ix)-out.H_direct_reduced_f),1e-12);
            add([id '_scalar_roi_reflect'],abs(q.native_reflect(ix)-out.H_reflect_reduced_f),1e-12);
            % Recompute end-plane boundary metrics from retained raw 2D data.
            field=out.psifinal_xy; [X,Y]=meshgrid(out.x,out.y);
            edge=abs(X)>=.45*cfg.width_m | abs(Y)>=.45*cfg.width_m;
            bd=false(cfg.n); bd([1 end],:)=true; bd(:,[1 end])=true;
            power=abs(field).^2;
            add([id '_direct_final_edge_reconstruction'],abs(sum(power(edge))/sum(power(:))-q.trace.direct(end,3)),1e-12);
            add([id '_direct_final_boundary_reconstruction'],abs(20*log10(max(abs(field(bd)))/max(abs(field(:))))-q.trace.direct(end,4)),1e-10);
            clear field out X Y power edge bd
        end
        clear q
        assert(toc(t)<=3600 && budget.used_s+toc(t)<=14400,'Audit time gate.');
    end
    if strcmp(audit_kind,'reference_audit')
        compare('P4','A-C4','PE_AS',1e-10);
        for f=[4 6 8]
            compare(['A-C' num2str(f)],['A-W' num2str(f)],'window',1e-3);
            compare(['A-C' num2str(f)],['A-G' num2str(f)],'grid',1e-3);
        end
        compare('A-W4','A-R4','reference_window',1e-3);
    elseif strcmp(audit_kind,'stage1_audit')
        for f=[4 6 8]
            compare(['P' num2str(f)],['A-C' num2str(f)],'PE_AS',1e-10);
            compare(['A-C' num2str(f)],['A-W' num2str(f)],'window',1e-3);
            compare(['A-C' num2str(f)],['A-G' num2str(f)],'grid',1e-3);
        end
        compare('A-W4','A-R4','reference_window',1e-3);
        for f=[4 8], compare(['P' num2str(f)],['Z' num2str(f)],'step',1e-3); end
    elseif strcmp(audit_kind,'stage2_audit')
        for f=[4 6 8]
            compare(['SP' num2str(f)],['P' num2str(f)],'sponge_no_sponge',1e-3);
            compare(['SP' num2str(f)],['A-C' num2str(f)],'sponge_AS',1e-3);
        end
    end
    receipt.status='AUDIT_PASS';
    receipt.failed_checks=sum(~[checks.pass]);
    if receipt.failed_checks>0, receipt.status='AUDIT_FAILED';
    elseif strcmp(audit_kind,'stage1_audit'), receipt.status='STAGE1_QUALIFIED';
    elseif strcmp(audit_kind,'stage2_audit'), receipt.status='SPONGE_QUALIFIED'; end
    receipt.exit_code=double(receipt.failed_checks>0);
    fid=fopen(fullfile(out_dir,'stderr.log'),'w'); fclose(fid);
catch ex
    receipt.error_identifier=ex.identifier; receipt.error_message=ex.message;
    receipt.failed_checks=max(1,sum(~[checks.pass]));
    fid=fopen(fullfile(out_dir,'stderr.log'),'w'); fprintf(fid,'%s\n',getReport(ex,'extended','hyperlinks','off')); fclose(fid);
    fprintf('M1_AUDIT_FAILURE %s\n',ex.message);
end
receipt.elapsed_s=toc(t); budget.used_s=budget.used_s+receipt.elapsed_s; receipt.cumulative_s=budget.used_s;
receipt.failed_or_unrun_checks=checks(~[checks.pass]);
save(fullfile(out_dir,'audit.mat'),'checks','rows','receipt');
writetable(struct2table(checks),fullfile(out_dir,'checks.csv'));
if ~isempty(rows), writetable(struct2table(rows),fullfile(out_dir,'comparisons.csv')); end
jsonfile(fullfile(run_dir,'budget.json'),budget);
clear dc
jsonfile(fullfile(out_dir,'checkpoint.json'),receipt);
fprintf('M1_AUDIT_RECEIPT %s\n',jsonencode(receipt));

    function add(name,value,limit,strict)
        if nargin<4, strict=false; end
        pass=isfinite(value) && value<=limit;
        if value==-Inf, pass=true; end
        if strict, pass=isfinite(value) && value<limit; end
        checks(end+1)=struct('name',name,'value',value,'limit',limit,'pass',pass);
    end
    function compare(a_id,b_id,kind,limit)
        a=load(fullfile(run_dir,a_id,'payload.mat'),'payload'); a=a.payload;
        b=load(fullfile(run_dir,b_id,'payload.mat'),'payload'); b=b.payload;
        [ok,index]=ismember(a.coords_m,b.coords_m,'rows'); assert(all(ok),'Non-nested coordinates.');
        component={'direct','reflect','total'};
        for radius=[.5 1 2]
            mask=hypot(a.coords_m(:,1),a.coords_m(:,2))<=radius;
            coordinates=a.coords_m(mask,:);
            scale=b.physical_direct(index(mask));
            for cc=1:3
                name=component{cc}; aa=a.(['physical_' name])(mask); bb=b.(['physical_' name])(index(mask));
                denominator=norm(bb); if strcmp(name,'total'), denominator=norm(scale); end
                l2=norm(aa-bb)/denominator;
                valid=abs(bb)>max(abs(bb))*1e-12;
                amp=abs(20*log10(abs(aa(valid))./abs(bb(valid))));
                phase=angle(aa(valid).*conj(bb(valid)));
                [max_amp,where]=max(amp); vcoords=coordinates(valid,:);
                [max_phase,phase_where]=max(abs(phase));
                [max_complex,wi]=max(abs(aa-bb)/max(abs(scale)));
                tag=sprintf('%s_%s_%s_%s_rho%g',kind,a_id,b_id,name,radius);
                add([tag '_complex_l2'],l2,limit);
                add([tag '_max_amplitude_db'],max_amp,.05);
                add([tag '_phase_rms_rad'],sqrt(mean(phase.^2)),.01);
                row=struct('kind',kind,'candidate',a_id,'reference',b_id, ...
                    'frequency_hz',a.config.f_hz,'component',name,'radius_m',radius, ...
                    'complex_l2',l2,'max_amplitude_db',max_amp,'amplitude_x_m',vcoords(where,1), ...
                    'amplitude_y_m',vcoords(where,2),'phase_rms_rad',sqrt(mean(phase.^2)), ...
                    'max_phase_rad',max_phase,'phase_x_m',vcoords(phase_where,1),'phase_y_m',vcoords(phase_where,2), ...
                    'max_direct_scaled_abs',max_complex,'complex_x_m',coordinates(wi,1), ...
                    'complex_y_m',coordinates(wi,2),'excluded_points',sum(~valid), ...
                    'denominator',string(conditional(strcmp(name,'total'),'direct_L2','component_L2')));
                if isempty(rows), rows=row; else, rows(end+1)=row; end %#ok<AGROW>
            end
        end
    end
end
function y=conditional(c,a,b)
if c, y=a; else, y=b; end
end
function value=sha(file)
md=java.security.MessageDigest.getInstance('SHA-256'); fid=fopen(file,'rb'); c=onCleanup(@()fclose(fid)); %#ok<NASGU>
while ~feof(fid), bytes=fread(fid,2^20,'*uint8'); if ~isempty(bytes), md.update(typecast(bytes,'int8')); end, end
value=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
function jsonfile(file,s)
fid=fopen(file,'w','n','UTF-8'); fprintf(fid,'%s',jsonencode(s,PrettyPrint=true)); fclose(fid);
end
