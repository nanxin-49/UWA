function receipt = run_piston_p0p1_remote_job(job_id,run_root)
%RUN_PISTON_P0P1_REMOTE_JOB Isolated remote, single-thread job dispatcher.
arguments
    job_id (1,:) char
    run_root (1,:) char
end
out=fullfile(run_root,'jobs',job_id); assert(~isfolder(out)); mkdir(out);
t0=tic; receipt=struct('job_id',job_id,'status','RUNNING','matlab_version',version, ...
    'pid',feature('getpid'),'execution_mode','remote -batch -singleCompThread', ...
    'failed_checks',0,'elapsed_s',0,'result_dir',out);
write_json(fullfile(out,'receipt.json'),receipt);
try
    switch job_id
        case 'prepare'
            cfg=struct('c_mps',1500,'rho_kgpm3',1025,'v0_mps',.001,'z_source_m',100, ...
                'frequencies_hz',[4000 6000 8000],'widths_m',[.2 .1], ...
                'launch_offsets_m',[1 2],'quadrature_order',64,'phase_contract_version','2.0.0', ...
                'rx_x_m',(-40:.025:40).','rx_z_m',3,'surface_A_m',.01);
            cfg.gates=struct('source_relative',1e-7,'gradient_relative',1e-6, ...
                'phase_absolute',1e-10,'pe_input_gradient_relative',1e-4, ...
                'pe_window_relative',1e-4,'pe_grid_relative',1e-4,'pe_step_relative',1e-10, ...
                'physical_free_flat_relative',1e-3,'bie_convergence_relative',1e-4, ...
                'bie_boundary_relative',1e-4,'fk_convergence_relative',1e-5, ...
                'receiver_interleaving_relative',1e-12,'resolved_factor',5);
            p_direct=complex(zeros(numel(cfg.rx_x_m),2)); p_mirror=p_direct;
            for wi=1:2
                source=local_source(cfg,4000,cfg.widths_m(wi));
                p_direct(:,wi)=baffled_strip_piston_field_vertical(cfg.rx_x_m,3,source);
                p_mirror(:,wi)=-baffled_strip_piston_field_vertical(cfg.rx_x_m,-3,source);
            end
            save(fullfile(run_root,'common.mat'),'cfg','p_direct','p_mirror','-v7');
            receipt.status='COMMON_SOURCE_PREPARED';
        case 'p0_source'
            q=load(fullfile(run_root,'common.mat')); cfg=q.cfg; metrics=[];
            for f=cfg.frequencies_hz
                for W=cfg.widths_m
                    s=local_source(cfg,f,W);
                    x=[0;.13;-.2;.2]; z=100-[.02;.04;.1;.4];
                    [p,px,pz]=baffled_strip_piston_field_vertical(x,z,s,128);
                    [a,ax,az]=piston_independent_angular_spectrum_vertical(x,z,s);
                    as_error=norm([p-a;px-ax;pz-az])/norm([a;ax;az]);
                    h=1e-4; fn=@(xx,zz)baffled_strip_piston_field_vertical(xx,zz,s,128);
                    dx=(-fn(x+2*h,z)+8*fn(x+h,z)-8*fn(x-h,z)+fn(x-2*h,z))/(12*h);
                    dz=(-fn(x,z+2*h)+8*fn(x,z+h)-8*fn(x,z-h)+fn(x,z-2*h))/(12*h);
                    grad_error=norm([dx-px;dz-pz])/norm([px;pz]);
                    % The integrated normal velocity is independently checked
                    % from the angular spectrum at a small positive distance.
                    d=1e-6; k=2*pi*f/s.c_mps; a=atan(W/(2*d));
                    scale=s.rho_kgpm3*2*pi*f*s.v0_mps/2;
                    vzgrad=scale*k*d*integral(@(t)besselh(1,1,k*d./cos(t))./cos(t), ...
                        -a,a,'RelTol',1e-10,'AbsTol',1e-10);
                    vz=vzgrad/(1i*s.rho_kgpm3*2*pi*f);
                    velocity_error=abs(vz+s.v0_mps)/s.v0_mps;
                    phase_error=local_phase_control(s);
                    metrics=[metrics;f W as_error grad_error velocity_error phase_error]; %#ok<AGROW>
                end
            end
            T=array2table(metrics,'VariableNames',{'frequency_hz','total_width_m','AS_pressure_gradient_relative','finite_difference_gradient_relative','near_face_velocity_relative','phase_v2_absolute'});
            writetable(T,fullfile(out,'metrics.csv')); save(fullfile(out,'data.mat'),'T','cfg');
            assert(all(T.AS_pressure_gradient_relative<1e-7) && all(T.finite_difference_gradient_relative<1e-6) && all(T.near_face_velocity_relative<1e-3) && all(T.phase_v2_absolute<1e-10));
            receipt.status='P0_SOURCE_CHECKS_PASS';
        case {'p0_pe','p0_bie','p0_fullfk'}
            q=load(fullfile(run_root,'common.mat')); cfg=q.cfg; metrics=[]; inputs=cell(3,2);
            for fi=1:3
                for wi=1:2
                    s=local_source(cfg,cfg.frequencies_hz(fi),cfg.widths_m(wi));
                    x=(-40:.05:40).'; [z,~]=closeout_s2b_surface(x,.01);
                    if strcmp(job_id,'p0_pe'), z=99+zeros(size(x)); end
                    [p,px,pz]=baffled_strip_piston_field_vertical(x,z,s);
                    if strcmp(job_id,'p0_fullfk')
                        src=local_fk_source(s); [u,ux,uz]=fk_single_receiver_incident(src,x,z);
                    elseif strcmp(job_id,'p0_bie')
                        incident=@(xx,zz)baffled_strip_piston_field_vertical(xx,zz,s);
                        u=incident(x,z); ux=px; uz=pz;
                    else
                        % Actual PE source array, explicit carrier removal.
                        [src,pcfg]=local_pe_input(cfg,s,512,.025,1,0,.1875);
                        ix=round((x-src.x(1))/.025)+1;
                        u=pcfg.source_field_override(ix).'*exp(1i*2*pi*s.frequency_hz/s.c_mps);
                        k=2*pi*s.frequency_hz/s.c_mps;
                        kx=2*pi/512*[0:numel(src.x)/2-1 -numel(src.x)/2:-1];
                        spectrum=fft(src.p.*src.taper);
                        dpx=ifft(1i*kx.*spectrum);
                        dpz=ifft(-1i*sqrt(complex(k^2-kx.^2,0)).*spectrum);
                        ux=reshape(dpx(ix),[],1); uz=reshape(dpz(ix),[],1);
                    end
                    err=norm(u-p)/norm(p);
                    gradient_error=norm([ux-px;uz-pz])/norm([px;pz]);
                    inputs{fi,wi}=struct('source',s,'x_m',x,'z_m',z,'p',u,'px',ux,'pz',uz);
                    metrics=[metrics;s.frequency_hz s.total_width_m err gradient_error]; %#ok<AGROW>
                end
            end
            T=array2table(metrics,'VariableNames',{'frequency_hz','total_width_m','actual_input_relative','actual_gradient_relative'});
            writetable(T,fullfile(out,'metrics.csv')); save(fullfile(out,'data.mat'),'T','cfg','inputs');
            assert(all(T.actual_input_relative<1e-7) && all(T.actual_gradient_relative<1e-4)); receipt.status='P0_MODEL_INPUT_PASS';
        case 'gaussian_regression'
            run_new_ubuntu_minimal_channel_baseline(fullfile(out,'direct'),'direct_4k');
            run_new_ubuntu_minimal_channel_baseline(fullfile(out,'reflect'),'reflect_4k');
            % Compare implicit vs explicit default source under the same API.
            d=load(fullfile(out,'reflect','channel_result.mat'));
            par=d.paramsV; par.source_mode='gaussian'; b=vertical_channel_model(par);
            err=max(abs([b.H_f(:)-d.output.H_f(:);b.H_direct_f(:)-d.output.H_direct_f(:);b.H_reflect_f(:)-d.output.H_reflect_f(:)]));
            assert(err==0); receipt.gaussian_default_max_abs=err;
            receipt.status='GAUSSIAN_DEFAULT_REGRESSION_PASS';
        otherwise
            receipt=local_p1(job_id,run_root,out,receipt);
    end
    receipt.elapsed_s=toc(t0); write_json(fullfile(out,'receipt.json'),receipt);
    fprintf('%s\n',jsonencode(receipt));
catch ME
    receipt.status='FAILED_OR_INCONCLUSIVE'; receipt.failed_checks=1; receipt.elapsed_s=toc(t0);
    receipt.error=ME.message; write_json(fullfile(out,'receipt.json'),receipt);
    fid=fopen(fullfile(out,'matlab_error.txt'),'w'); fprintf(fid,'%s',getReport(ME,'extended','hyperlinks','off')); fclose(fid);
    rethrow(ME);
end
end
function s=local_source(cfg,f,W)
s=struct('frequency_hz',f,'c_mps',cfg.c_mps,'rho_kgpm3',cfg.rho_kgpm3, ...
    'v0_mps',cfg.v0_mps,'total_width_m',W,'z_source_m',cfg.z_source_m);
end
function src=local_fk_source(s)
src=struct('frequency_hz',s.frequency_hz,'c0_mps',s.c_mps, ...
    'incident_field_fn',@(x,z)baffled_strip_piston_field_vertical(x,z,s));
end
function e=local_phase_control(s)
k=2*pi*s.frequency_hz/s.c_mps;
p=baffled_strip_piston_field_vertical(0,3,s); r=-baffled_strip_piston_field_vertical(0,-3,s);
red=struct('direct_f',p*exp(-1i*k*97),'reflect_coh_f',r*exp(-1i*k*103));
g=struct('z_tx',100,'z_rx',3,'c0',s.c_mps);
[a,m]=apply_pe_channel_phase_reference_vertical(s.frequency_hz,red,g,'absolute_physical');
b=apply_pe_channel_phase_reference_vertical(s.frequency_hz,red,g,'direct_dsp');
e=max(abs([a.direct_f-p;a.reflect_fm-r;b.total_fm-conj(p+r)*exp(1i*k*97)]));
assert(strcmp(m.phase_contract_version,'2.0.0'));
end
function [src,pcfg]=local_pe_input(cfg,s,W,dx,offset,A,step)
N=round(W/dx); assert(mod(N,2)==0); x=-W/2+(0:N-1)*dx;
[p,px,pz]=baffled_strip_piston_field_vertical(x,100-offset,s);
r=abs(x)/(W/2); taper=ones(size(x)); mid=r>.7;
t=min(1,(r(mid)-.7)/.3); taper(mid)=1-10*t.^3+15*t.^4-6*t.^5;
k=2*pi*s.frequency_hz/s.c_mps;
pcfg=struct('frequency_hz',s.frequency_hz,'c0_mps',s.c_mps,'xw_m',W,'nx',N, ...
    'z_tx_m',100-offset,'z_rx_m',3,'sigma_src_m',.3,'surface_reflect_coeff',-1, ...
    'step_m',step,'x_rx_m',0,'reflection_model','model0_normal', ...
    'source_field_override',p.*taper*exp(-1i*k*offset),'eta_z_m',closeout_s2b_surface(x,A));
src=struct('source',s,'x',x,'p',p,'px',px,'pz',pz,'taper',taper, ...
    'numerical_plane_z_m',100-offset,'physical_plane_z_m',100,'offset_m',offset, ...
    'launch_reduced_formula','p_physical(z_launch)*exp(-i*k*offset)', ...
    'taper_semantics','C2 numerical launch window; unity in central 70 percent, independently converge W/dx/offset');
end
function receipt=local_p1(job_id,run_root,out,receipt)
q=load(fullfile(run_root,'common.mat')); cfg=q.cfg;
% P1 dispatch is added after the common-source pilot; never default a solve.
error('piston:P1NotConfigured','Unknown or not-yet-configured job %s; no solve dispatched.',job_id);
end
function write_json(p,v)
fid=fopen(p,'w'); assert(fid>0); fprintf(fid,'%s',jsonencode(v,PrettyPrint=true)); fclose(fid);
end
