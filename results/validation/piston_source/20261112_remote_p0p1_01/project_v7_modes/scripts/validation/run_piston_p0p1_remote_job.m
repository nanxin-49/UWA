function receipt = run_piston_p0p1_remote_job(job_id,run_root)
%RUN_PISTON_P0P1_REMOTE_JOB Isolated remote, single-thread job dispatcher.
arguments
    job_id (1,:) char
    run_root (1,:) char
end
if contains(job_id,'_pe_cont_') || startsWith(job_id,'gate_cont_')
    receipt=run_piston_matched_p1_job(job_id,run_root);
    return;
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
        case {'p0_pe','p0_pe_refined','p0_pe_refined_recovery01','p0_pe_refined_window32k','p0_bie','p0_fullfk'}
            q=load(fullfile(run_root,'common.mat')); cfg=q.cfg; metrics=[]; inputs=cell(3,2);
            for fi=1:3
                for wi=1:2
                    s=local_source(cfg,cfg.frequencies_hz(fi),cfg.widths_m(wi));
                    x=(-40:.05:40).'; [z,~]=closeout_s2b_surface(x,.01);
                    if startsWith(job_id,'p0_pe'), z=99+zeros(size(x)); end
                    [p,px,pz]=baffled_strip_piston_field_vertical(x,z,s);
                    if strcmp(job_id,'p0_fullfk')
                        src=local_fk_source(s); [u,ux,uz]=fk_single_receiver_incident(src,x,z);
                    elseif strcmp(job_id,'p0_bie')
                        incident=@(xx,zz)baffled_strip_piston_field_vertical(xx,zz,s);
                        u=incident(x,z); ux=px; uz=pz;
                    else
                        % Actual PE source array, explicit carrier removal.
                        window=512; if startsWith(job_id,'p0_pe_refined'), window=8192; end
                        if strcmp(job_id,'p0_pe_refined_window32k'), window=32768; end
                        [src,pcfg]=local_pe_input(cfg,s,window,.025,1,0,.1875);
                        ix=round((x-src.x(1))/.025)+1;
                        u=pcfg.source_field_override(ix).'*exp(1i*2*pi*s.frequency_hz/s.c_mps);
                        k=2*pi*s.frequency_hz/s.c_mps;
                        kx=2*pi/window*[0:numel(src.x)/2-1 -numel(src.x)/2:-1];
                        spectrum=fft(src.p.*src.taper);
                        dpx=ifft(1i*kx.*spectrum);
                        dpz=ifft(-1i*sqrt(complex(k^2-kx.^2,0)).*spectrum);
                        ux=reshape(dpx(ix),[],1); uz=reshape(dpz(ix),[],1);
                    end
                    err=norm(u-p)/norm(p);
                    gradient_error=norm([ux-px;uz-pz])/norm([px;pz]);
                    inputs{fi,wi}=struct('source',s,'x_m',x,'z_m',z,'p',u,'px',ux,'pz',uz);
                    if startsWith(job_id,'p0_pe')
                        inputs{fi,wi}.numerical_window_m=window;
                        inputs{fi,wi}.dx_m=.025;
                        inputs{fi,wi}.source_offset_m=1;
                        inputs{fi,wi}.pressure_quadrature_order=24;
                    end
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
        case 'p0_finalize'
            selected_refined='p0_pe_refined';
            selection_file=fullfile(run_root,'p0_recovery_selection.json');
            if isfile(fullfile(run_root,'p0_active_selection_v5.json'))
                selection_file=fullfile(run_root,'p0_active_selection_v5.json');
            end
            if isfile(selection_file)
                selection=jsondecode(fileread(selection_file)); selected_refined=selection.active_refined_job_id;
                assert(~isempty(regexp(selected_refined,'^p0_pe_refined[A-Za-z0-9_]*$','once')));
            end
            required={'prepare','p0_source',selected_refined,'p0_bie','p0_fullfk','gaussian_regression'};
            for j=1:numel(required)
                rr=jsondecode(fileread(fullfile(run_root,'jobs',required{j},'receipt.json')));
                assert(rr.failed_checks==0 && (contains(rr.status,'PASS') || strcmp(rr.status,'COMMON_SOURCE_PREPARED')));
            end
            q=load(fullfile(run_root,'common.mat')); cfg=q.cfg;
            offset_error=0; launch_quadrature_error=0;
            for f=cfg.frequencies_hz
                k=2*pi*f/cfg.c_mps;
                for W=cfg.widths_m
                    s=local_source(cfg,f,W);
                    probes=[-4096;-100;-40;-.1;0;.1;40;100;4096];
                    for offset=cfg.launch_offsets_m
                        pa=baffled_strip_piston_field_vertical(probes,100-offset,s,24);
                        pb=baffled_strip_piston_field_vertical(probes,100-offset,s,64);
                        launch_quadrature_error=max(launch_quadrature_error,max(abs(pa-pb)./abs(pb)));
                        [src,pcfg]=local_pe_input(cfg,s,512,.025,offset,0,.1875);
                        mask=abs(src.x)<=40;
                        reconstructed=pcfg.source_field_override(mask)*exp(1i*k*offset);
                        offset_error=max(offset_error,norm(reconstructed-src.p(mask))/norm(src.p(mask)));
                    end
                end
            end
            % Regression of the unchanged frozen-Fourier FullFK fallback.
            fsrc=struct('x0_m',-.5,'z_tx_m',100,'kx',[-2;0;3], ...
                'kz',sqrt((2*pi*4000/1500)^2-[-2;0;3].^2),'coeff',[1+2i;3;.4i],'n_original',3);
            x=[-1;0;2]; z=[3;4;5]; [a,ax,az]=fk_single_receiver_incident(fsrc,x,z);
            E=exp(1i*(x-fsrc.x0_m)*fsrc.kx.').*exp(-1i*(z-100)*fsrc.kz.');
            fallback_error=max(abs([a-E*fsrc.coeff/3;ax-E*(1i*fsrc.kx.*fsrc.coeff)/3;az-E*(-1i*fsrc.kz.*fsrc.coeff)/3]));
            assert(offset_error<1e-10 && fallback_error<1e-12 && launch_quadrature_error<1e-7);
            receipt.launch_offset_phase_relative=offset_error; receipt.frozen_fourier_fallback_absolute=fallback_error;
            receipt.launch_24_vs_64_max_relative=launch_quadrature_error;
            ql_path=fullfile(fileparts(fileparts(fileparts(fileparts(run_root)))),'inputs','quicklook_field_4kHz.mat');
            if isfile(ql_path)
                old=load(ql_path,'x_rx_m','p_rx_f'); ix=old.x_rx_m>=-40 & old.x_rx_m<=40;
                reuse_error=norm(q.p_direct(1:2:end,:)-old.p_rx_f(ix,1:2),'fro')/norm(old.p_rx_f(ix,1:2),'fro');
                assert(reuse_error<1e-7); receipt.quicklook_reuse_relative=reuse_error;
            else
                error('piston:MissingPriorEvidence','Required verified quicklook pressure evidence missing.');
            end
            project=fileparts(fileparts(fileparts(mfilename('fullpath'))));
            code_analysis=checkcode(fullfile(project,'scripts','validation','run_piston_p0p1_remote_job.m'),'-id');
            save(fullfile(out,'code_analysis.mat'),'code_analysis'); receipt.code_analysis_messages=numel(code_analysis);
            receipt.status='P0_PASS';
        case 'review_code_recovery01'
            project=fileparts(fileparts(fileparts(mfilename('fullpath'))));
            names={'src/propagation/baffled_strip_piston_field_vertical.m', ...
                'scripts/validation/run_piston_p0p1_remote_job.m', ...
                'scripts/validation/audit_piston_p0p1_remote.m', ...
                'scripts/reporting/render_piston_p0p1_remote.m', ...
                'scripts/validation/support/piston_independent_angular_spectrum_vertical.m'};
            analysis=cell(size(names));
            for j=1:numel(names), analysis{j}=checkcode(fullfile(project,names{j}),'-id'); end
            save(fullfile(out,'code_analysis.mat'),'names','analysis');
            receipt.code_files_checked=numel(names); receipt.message_counts=cellfun(@numel,analysis);
            write_json(fullfile(out,'code_analysis.json'),struct('files',{names},'messages',{analysis}));
            receipt.status='CODE_ANALYZER_EXECUTED_REVIEW_MESSAGES';
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
function [src,pcfg]=local_pe_input(~,s,W,dx,offset,A,step)
N=round(W/dx); assert(mod(N,2)==0); x=-W/2+(0:N-1)*dx;
% At >=1 m clearance, GL24 is cross-checked against GL64 in P0. Avoid
% recomputing unused full-window derivative arrays during every PE sweep.
p=baffled_strip_piston_field_vertical(x,100-offset,s,24); px=[]; pz=[];
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
    'pressure_quadrature_order',24,'gradient_arrays','not duplicated; P0 analytic and FFT gradients saved separately', ...
    'taper_semantics','C2 numerical launch window; unity in central 70 percent, independently converge W/dx/offset');
end
function receipt=local_p1(job_id,run_root,out,receipt)
q=load(fullfile(run_root,'common.mat')); cfg=q.cfg;
if startsWith(job_id,'gate_')
    receipt=local_gate(job_id(6:end),run_root,out,q,receipt); return;
end
parts=regexp(job_id,'^(free|flat|weak)_(pe|bie|fullfk)_l(\d+)$','tokens','once');
assert(~isempty(parts),'Unknown P1 job ID.');
stage=parts{1}; method=parts{2}; level=str2double(parts{3});
A=0; if strcmp(stage,'weak'), A=.01; end
result=struct('stage',stage,'method',method,'level',level,'frequency_hz',4000, ...
    'widths_m',cfg.widths_m,'x_m',cfg.rx_x_m,'z_m',3,'phase_contract_version','2.0.0', ...
    'representation','absolute_physical peak Pa; exp(-i omega t)', ...
    'direct',q.p_direct,'reflect',complex(zeros(size(q.p_direct))));
if strcmp(method,'pe')
    recipe=[512 .05 1 .1875;1024 .05 1 .1875;2048 .05 1 .1875;4096 .05 1 .1875; ...
        4096 .025 1 .1875;4096 .025 1 .09375;4096 .025 2 .09375];
    assert(level<=size(recipe,1)); rc=recipe(level,:); W=rc(1); dx=rc(2); offset=rc(3); step=rc(4);
    result.numeric_config=struct('window_m',W,'dx_m',dx,'source_offset_m',offset,'step_m',step);
    fields=cell(1,2); launches=fields;
    for wi=1:2
        s=local_source(cfg,4000,cfg.widths_m(wi));
        [launch,pcfg]=local_pe_input(cfg,s,W,dx,offset,A,step);
        if strcmp(stage,'free'), pcfg.surface_reflect_coeff=0; end
        pe=run_pe_physical_surface_model0_validation(pcfg);
        k=2*pi*4000/1500;
        pd=pe.direct_field*exp(1i*k*97);
        pr=pe.reflected_field*exp(1i*k*103);
        factor=round(dx/.025); assert(abs(dx/factor-.025)<1e-12);
        pd_fine=interpft(pd,numel(pd)*factor); pr_fine=interpft(pr,numel(pr)*factor);
        ix=round((cfg.rx_x_m+W/2)/.025)+1;
        result.direct(:,wi)=reshape(pd_fine(ix),[],1);
        result.reflect(:,wi)=reshape(pr_fine(ix),[],1);
        fields{wi}=pe; launches{wi}=launch;
    end
    result.pe_fields=fields; result.launch_planes=launches;
elseif strcmp(stage,'free')
    assert(level==1); result.numeric_config=struct('quadrature_order',128);
    for wi=1:2
        s=local_source(cfg,4000,cfg.widths_m(wi));
        if strcmp(method,'bie')
            result.direct(:,wi)=baffled_strip_piston_field_vertical(cfg.rx_x_m,3,s,128);
        else
            result.direct(:,wi)=fk_single_receiver_incident(local_fk_source(s),cfg.rx_x_m,3);
        end
    end
elseif strcmp(method,'bie')
    % One boundary matrix and factorization serve both aperture right sides.
    recipe=[60 3 32;60 6 32;60 6 64;80 6 64;100 6 64;100 8 64;120 8 64;100 8 96];
    assert(level<=size(recipe,1)); rc=recipe(level,:);
    bc=struct('frequency_hz',4000,'c0_mps',1500,'half_width_m',rc(1), ...
        'points_per_wavelength',rc(2),'panel_order',8,'self_quadrature_order',rc(3), ...
        'halfplane_h_m',-2,'window_plateau_ratio',.85, ...
        'eta_fn',@(x)closeout_s2b_surface(x,A),'eta_prime_fn',@(x)local_slope(x,A), ...
        'incident_fn',@(x,z)local_bie_input(x,z,cfg,A), ...
        'receiver_x_m',cfg.rx_x_m,'receiver_z_m',3);
    boundary=solve_helmholtz_bie_halfplane_vertical(bc);
    if A==0
        result.reflect=boundary.receiver_field;
    else
        % Exact flat exterior + compactly driven correction, same PDE.
        result.reflect=q.p_mirror+boundary.receiver_field;
    end
    result.numeric_config=rmfield(bc,{'eta_fn','eta_prime_fn','incident_fn'});
    result.boundary=boundary;
    result.actual_incident_nodes=local_bie_input(boundary.node_x_m,boundary.node_z_m,cfg,A);
    result.incident_semantics='flat: actual incident pressure; weak: incident+exact-flat-reflection boundary residual, correction solved by the same BIE';
elseif strcmp(method,'fullfk')
    panels=[200 400 800 1600]; assert(level<=4);
    result.numeric_config=struct('panels',panels(level),'gauss_order',8,'flat_exterior','analytic infinite plane');
    details=cell(1,2);
    for wi=1:2
        s=local_source(cfg,4000,cfg.widths_m(wi)); src=local_fk_source(s);
        v=closeout_s2b_fullfk_analytic_exterior(src,A,cfg.rx_x_m,3+zeros(size(cfg.rx_x_m)),panels(level),q.p_mirror(:,wi));
        v.source=s; result.reflect(:,wi)=v.receiver_field; details{wi}=v;
    end
    result.fullfk=details;
    if strcmp(stage,'flat')
        probes=(-40:10:40).'; L=200*2^(level-1);
        numerical=local_flat_fk_probes(cfg,probes,L);
        exact=complex(zeros(size(numerical)));
        for wi=1:2
            exact(:,wi)=-baffled_strip_piston_field_vertical(probes,-3,local_source(cfg,4000,cfg.widths_m(wi)));
        end
        result.flat_independent_probe=struct('x_m',probes,'L_m',L,'field',numerical,'exact',exact,'relative',norm(numerical-exact,'fro')/norm(exact,'fro'));
    end
end
result.total=result.direct+result.reflect;
result.H_direct_dsp=conj(result.direct)*exp(1i*2*pi*4000*97/1500);
result.H_reflect_dsp=conj(result.reflect)*exp(1i*2*pi*4000*97/1500);
result.H_total_dsp=result.H_direct_dsp+result.H_reflect_dsp;
assert(all(isfinite(result.total(:))));
save(fullfile(out,'data.mat'),'result','cfg','-v7');
T=local_point_table(result); writetable(T,fullfile(out,'points.csv'));
receipt.status='P1_DATA_COMPLETE_AUDIT_PENDING'; receipt.stage=stage; receipt.method=method; receipt.level=level;
if isfield(result,'boundary')
    receipt.node_count=result.boundary.node_count;
    receipt.linear_residual_max=max(result.boundary.linear_residual);
    receipt.boundary_residual_max=max(result.boundary.offgrid_boundary_residual);
end
end
function slope=local_slope(x,A)
[~,slope]=closeout_s2b_surface(x,A);
end
function p=local_bie_input(x,z,cfg,A)
x=x(:); z=z(:)+zeros(size(x)); p=complex(zeros(numel(x),2));
for wi=1:2
    s=local_source(cfg,4000,cfg.widths_m(wi));
    p(:,wi)=baffled_strip_piston_field_vertical(x,z,s);
    if A~=0, p(:,wi)=p(:,wi)-baffled_strip_piston_field_vertical(x,-z,s); end
end
end
function p=local_flat_fk_probes(cfg,x,L)
% Independent finite-extent two-term planar quadrature, no mirror supplied.
n=8; j=(1:n-1)'; beta=j./sqrt(4*j.^2-1); [V,D]=eig(diag(beta,1)+diag(beta,-1));
[g,ix]=sort(diag(D)); weights=2*(V(1,ix).^2)';
panels=ceil(4*L); half=L/panels; centers=-L+half+(0:panels-1)'*2*half;
xs=reshape((centers+half*g.').',[],1); qw=repmat(half*weights,panels,1);
k=2*pi*4000/1500; p=complex(zeros(numel(x),2));
for wi=1:2
    [u,~,uz]=baffled_strip_piston_field_vertical(xs,0,local_source(cfg,4000,cfg.widths_m(wi)));
    R=hypot(x-xs.',3); G=1i/4*besselh(0,1,k*R);
    Gn=1i*k/4*besselh(1,1,k*R).*3./R;
    p(:,wi)=-sum((G.*uz.'+Gn.*u.').*qw.',2);
end
end
function T=local_point_table(r)
n=numel(r.x_m); T=table(repmat(r.x_m,2,1),repelem(r.widths_m(:),n), ...
    real(r.direct(:)),imag(r.direct(:)),real(r.reflect(:)),imag(r.reflect(:)), ...
    real(r.total(:)),imag(r.total(:)),abs(r.reflect(:)), ...
    'VariableNames',{'x_m','total_width_m','direct_real_Pa','direct_imag_Pa', ...
    'reflect_real_Pa','reflect_imag_Pa','total_real_Pa','total_imag_Pa','reflect_abs_Pa'});
end
function receipt=local_gate(stage,root,out,q,receipt)
% Gates are frozen in common.mat and never retuned after seeing results.
models={'pe','bie','fullfk'}; selected=cell(1,3); rows=[];
for mi=1:3
    model=models{mi}; files=dir(fullfile(root,'jobs',[stage '_' model '_l*'],'data.mat'));
    assert(~isempty(files),'Missing model data.'); vals=cell(size(files)); levels=zeros(size(files));
    for ii=1:numel(files)
        v=load(fullfile(files(ii).folder,files(ii).name),'result'); vals{ii}=v.result; levels(ii)=v.result.level;
    end
    [levels,ord]=sort(levels); vals=vals(ord); best=vals{end}; selected{mi}=best;
    field='reflect'; if strcmp(stage,'free'), field='direct'; end
    for ii=2:numel(vals)
        a=vals{ii-1}.(field); b=vals{ii}.(field);
        rows=[rows;mi levels(ii-1) levels(ii) norm(a-b,'fro')/max(norm(b,'fro'),realmin)]; %#ok<AGROW>
    end
end
M=array2table(rows,'VariableNames',{'method_id','level_before','level_after','relative_L2'});
writetable(M,fullfile(out,'successive_convergence.csv'));
pe=selected{1}; bie=selected{2}; fk=selected{3}; checks=struct;
if strcmp(stage,'free')
    checks.pe_physical=norm(pe.direct-q.p_direct,'fro')/norm(q.p_direct,'fro')<q.cfg.gates.physical_free_flat_relative;
    checks.bie_physical=norm(bie.direct-q.p_direct,'fro')/norm(q.p_direct,'fro')<1e-7;
    checks.fk_physical=norm(fk.direct-q.p_direct,'fro')/norm(q.p_direct,'fro')<1e-7;
else
    if strcmp(stage,'flat')
        checks.pe_mirror=norm(pe.reflect-q.p_mirror,'fro')/norm(q.p_mirror,'fro')<q.cfg.gates.physical_free_flat_relative;
        checks.bie_mirror=norm(bie.reflect-q.p_mirror,'fro')/norm(q.p_mirror,'fro')<q.cfg.gates.physical_free_flat_relative;
        checks.fullfk_independent_plane=fk.flat_independent_probe.relative<1e-3;
    end
    checks.bie_residual=max(bie.boundary.offgrid_boundary_residual)<q.cfg.gates.bie_boundary_relative;
    checks.fk_convergence=local_pair(root,stage,'fullfk',2,3,field)<q.cfg.gates.fk_convergence_relative;
    checks.bie_quadrature=local_pair(root,stage,'bie',6,8,field)<q.cfg.gates.bie_convergence_relative;
    checks.bie_grid=local_pair(root,stage,'bie',5,6,field)<q.cfg.gates.bie_convergence_relative;
    checks.bie_window=local_pair(root,stage,'bie',6,7,field)<q.cfg.gates.bie_convergence_relative;
end
checks.pe_window=local_pair(root,stage,'pe',3,4,field)<q.cfg.gates.pe_window_relative;
checks.pe_grid=local_pair(root,stage,'pe',4,5,field)<q.cfg.gates.pe_grid_relative;
checks.pe_step=local_pair(root,stage,'pe',5,6,field)<q.cfg.gates.pe_step_relative;
checks.pe_launch_offset=local_pair(root,stage,'pe',6,7,field)<q.cfg.gates.pe_grid_relative;
receipt.checks=checks; receipt.failed_checks=nnz(~cell2mat(struct2cell(checks)));
receipt.status='P1_STAGE_PASS'; if receipt.failed_checks>0, receipt.status='P1_STAGE_INCONCLUSIVE'; end
save(fullfile(out,'stage_summary.mat'),'selected','checks','M','-v7');
write_json(fullfile(out,'checks.json'),checks);
end
function v=local_pair(root,stage,method,a,b,field)
% Explicit file paths keep the audit independent of caller working directory.
x=load(fullfile(root,'jobs',sprintf('%s_%s_l%d',stage,method,a),'data.mat'),'result');
y=load(fullfile(root,'jobs',sprintf('%s_%s_l%d',stage,method,b),'data.mat'),'result');
A=x.result.(field); B=y.result.(field);
v=max(vecnorm(A-B)./max(vecnorm(B),realmin));
% Keep edge segments visible instead of letting a whole-line mean hide them.
N=size(B,1); edges=round(linspace(1,N,5));
for j=1:4
    ix=edges(j):edges(j+1);
    v=max(v,max(vecnorm(A(ix,:)-B(ix,:))./max(vecnorm(B(ix,:)),realmin)));
end
end
function write_json(p,v)
fid=fopen(p,'w'); assert(fid>0); fprintf(fid,'%s',jsonencode(v,PrettyPrint=true)); fclose(fid);
end
