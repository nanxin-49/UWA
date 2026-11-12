function receipt = run_piston_source_quicklook(run_id)
%RUN_PISTON_SOURCE_QUICKLOOK Exact 2-D free-field visualization data only.
% No production-source integration or P0 qualification is performed.
arguments
    run_id (1,:) char
end
assert(~isempty(regexp(run_id,'^[A-Za-z0-9_]+$','once')));
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
out=fullfile(root,'results','validation','piston_source_quicklook',run_id);
assert(~isfolder(out),'Use a unique run ID.'); mkdir(out);
t0=tic;
cfg=struct('c_mps',1500,'rho_kgpm3',1025,'v0_mps',.001, ...
    'z_tx_m',100,'z_rx_m',3,'frequencies_hz',[4000 6000 8000], ...
    'total_widths_m',[.20 .10 .05],'field_nquad',24,'near_nquad',96, ...
    'phasor','exp(-i omega t); outgoing H0 first kind', ...
    'velocity_definition','peak, uniform, in phase, upward', ...
    'historical_5cm','circular radius .05 m; strip half-width .05 m => total width .10 m', ...
    'qualification','QUICKLOOK_ONLY; P0 NOT EXECUTED', ...
    'random_seed','none; deterministic quadrature');
cfg.tolerances=struct('quadrature_complex_relative',1e-7, ...
    'angular_spectrum_complex_relative',1e-7, ...
    'farfield_normalized_absolute',1e-3,'csv_relative',1e-12);
[u0,m0]=memory;
resource=struct('available_start_B',m0.PhysicalMemory.Available, ...
    'matlab_start_B',u0.MemUsedMATLAB,'reserve_B',.3*m0.PhysicalMemory.Available, ...
    'estimated_additional_peak_B',180*1024^2,'available_min_sampled_B',m0.PhysicalMemory.Available, ...
    'cpu','existing local MATLAB; serial; no pool or thread-count change', ...
    'pid',feature('getpid'),'release',version,'peak_measurement','boundary samples, not a peak profiler');
assert(resource.available_start_B-resource.estimated_additional_peak_B>resource.reserve_B);
x_m=-120:.5:120; z_m=0:.5:99.5;
[X,Z]=meshgrid(x_m,z_m);
x_near_m=-.4:.005:.4; s_near_m=.01:.005:1;
[XN,SN]=meshgrid(x_near_m,s_near_m);
x_rx_m=(-120:.05:120).';
theta_deg=(-89.5:.25:89.5).'; circle_radius_m=97;
save(fullfile(out,'config.mat'),'cfg','resource');
write_state('RUNNING',0);
try
    logmsg('RUN_ID=%s PID=%d VERSION=%s\n',run_id,resource.pid,version);
    for fi=1:3
        f_hz=cfg.frequencies_hz(fi); tf=tic;
        p_field_f=complex(zeros([size(X),3]));
        p_near_f=complex(zeros([size(XN),3]));
        p_rx_f=complex(zeros(numel(x_rx_m),3));
        p_circle_f=complex(zeros(numel(theta_deg),3));
        D_sinc=zeros(numel(theta_deg),3);
        for wi=1:3
            [~,ms]=memory; resource.available_min_sampled_B=min(resource.available_min_sampled_B,ms.PhysicalMemory.Available);
            assert(ms.PhysicalMemory.Available>resource.reserve_B,'RAM reserve gate.');
            W=cfg.total_widths_m(wi);
            p_field_f(:,:,wi)=piston_strip_pressure_quicklook(X,100-Z,W,f_hz,cfg,cfg.field_nquad);
            p_near_f(:,:,wi)=piston_strip_pressure_quicklook(XN,SN,W,f_hz,cfg,cfg.near_nquad);
            p_rx_f(:,wi)=piston_strip_pressure_quicklook(x_rx_m,97+zeros(size(x_rx_m)),W,f_hz,cfg,cfg.field_nquad);
            p_circle_f(:,wi)=piston_strip_pressure_quicklook(circle_radius_m*sind(theta_deg), ...
                circle_radius_m*cosd(theta_deg),W,f_hz,cfg,cfg.field_nquad);
            q=pi*f_hz/cfg.c_mps*W*sind(theta_deg);
            d=ones(size(q)); nz=q~=0; d(nz)=sin(q(nz))./q(nz); D_sinc(:,wi)=d;
        end
        assert(all(isfinite(p_field_f(:))) && all(isfinite(p_near_f(:))));
        T=table(x_rx_m,3+zeros(size(x_rx_m)),'VariableNames',{'x_m','z_m'});
        for wi=1:3
            tag=sprintf('W%02dcm',round(100*cfg.total_widths_m(wi)));
            p=p_rx_f(:,wi); pa=abs(p); p_axis=pa(x_rx_m==0);
            T.([tag '_real_Pa'])=real(p); T.([tag '_imag_Pa'])=imag(p);
            T.([tag '_abs_peak_Pa'])=pa;
            T.([tag '_axis_normalized_dB'])=20*log10(pa/p_axis);
        end
        T.ratio20_over_5_dB=20*log10(abs(p_rx_f(:,1))./abs(p_rx_f(:,3)));
        T.ratio20_over_10_dB=20*log10(abs(p_rx_f(:,1))./abs(p_rx_f(:,2)));
        writetable(T,fullfile(out,sprintf('receiver_%dkHz.csv',f_hz/1000)));
        elapsed_frequency_s=toc(tf);
        save(fullfile(out,sprintf('field_%dkHz.mat',f_hz/1000)), ...
            'cfg','f_hz','x_m','z_m','p_field_f','x_near_m','s_near_m','p_near_f', ...
            'x_rx_m','p_rx_f','theta_deg','circle_radius_m','p_circle_f','D_sinc','elapsed_frequency_s','-v7');
        logmsg('F_HZ=%d ELAPSED_S=%.3f AXIS_PEAK_PA_W20_W10_W5=%s\n', ...
            f_hz,elapsed_frequency_s,mat2str(abs(p_rx_f(x_rx_m==0,:)),9));
        write_state('RUNNING',fi);
    end
    resource.elapsed_s=toc(t0); [u1,m1]=memory;
    resource.matlab_end_B=u1.MemUsedMATLAB; resource.available_end_B=m1.PhysicalMemory.Available;
    save(fullfile(out,'config.mat'),'cfg','resource');
    write_state('DATA_COMPLETE_AUDIT_PENDING',3);
    fid=fopen(fullfile(out,'exit_code.txt'),'w'); fprintf(fid,'0\n'); fclose(fid);
    receipt=struct('run_id',run_id,'status','DATA_COMPLETE_AUDIT_PENDING', ...
        'elapsed_s',toc(t0),'result_dir',out,'failed_checks',0,'unrun_checks',{{'independent audit','figure reconstruction','visual review','P0 qualification (out of scope)'}});
    logmsg('%s\n',jsonencode(receipt));
catch ME
    fid=fopen(fullfile(out,'stderr.log'),'w'); fprintf(fid,'%s\n',getReport(ME,'extended','hyperlinks','off')); fclose(fid);
    fid=fopen(fullfile(out,'exit_code.txt'),'w'); fprintf(fid,'1\n'); fclose(fid);
    write_state('FAILED',0); rethrow(ME);
end
    function logmsg(varargin)
        msg=sprintf(varargin{:}); fprintf('%s',msg);
        fid=fopen(fullfile(out,'stdout.log'),'a'); fprintf(fid,'%s',msg); fclose(fid);
    end
    function write_state(status,n)
        checkpoint=struct('run_id',run_id,'status',status,'completed_frequency_count',n, ...
            'completed_steps',{{'historical_definition_review','configuration_frozen'}}, ...
            'config',cfg,'resource_usage',resource,'execution_mode','existing local MATLAB MCP', ...
            'matlab_release',version,'elapsed_s',toc(t0),'code_versions_sha256','pending external manifest', ...
            'input_versions_sha256','pending external manifest','artifacts_sha256','pending external manifest', ...
            'failed_or_unrun_checks',{{'independent audit pending','P0 not authorized in this quicklook'}}, ...
            'next_allowed_step','audit and render this quicklook only; stop before P0/P1');
        fid=fopen(fullfile(out,'checkpoint_data.json'),'w'); fprintf(fid,'%s',jsonencode(checkpoint,PrettyPrint=true)); fclose(fid);
    end
end
