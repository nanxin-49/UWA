function receipt = run_k1_recovery_model_job(run_dir,job_id)
%RUN_K1_PISTON_PM_JOB Frozen R2 formula, K0 coefficient-defined U_PM/N adapter.
% Research only. No production defaults, sparse receivers, or FastFK.
arguments
    run_dir (1,:) char
    job_id (1,:) char
end
out=fullfile(run_dir,'jobs',job_id); assert(~isfolder(out)); mkdir(out);
diary(fullfile(out,'diary.log')); cleanup=onCleanup(@()diary('off'));
t0=tic; receipt=struct('job_id',job_id,'status','RUNNING','failed_checks',0, ...
    'pid',feature('getpid'),'matlab_version',version);
try
assert(strcmp(version('-release'),'2025b'));
if strcmp(job_id,'models')
    k1_recovery_event(run_dir,'MODELS','BEGIN');
    k1_recovery_event(run_dir,'INPUT','BEGIN');
    frozen=jsondecode(fileread(fullfile(run_dir,'input','frozen_config.json')));
    profile=readtable(fullfile(run_dir,'input','coefficients.csv'));
    receivers=readtable(fullfile(run_dir,'input','receivers.csv'));
    x=receivers.x_m; z=receivers.z_m;
    assert(numel(x)==411 && all(z==3) && numel(unique(x))==411);
    assert(height(profile)==frozen.modes && abs(max(profile.k_rad_per_m)-frozen.modes*2*pi/160)<1e-14);
    source=struct('frequency_hz',4000,'c_mps',1500,'rho_kgpm3',1025, ...
        'v0_mps',.001,'total_width_m',.1,'z_source_m',100);
    cfg=struct('source',source,'seed',frozen.seed,'period_m',160, ...
        'surface_plateau_halfwidth_m',42,'surface_support_halfwidth_m',50, ...
        'phase_contract_version','2.0.0','representation','absolute_physical peak Pa; exp(-i omega t)', ...
        'phase_mask_fraction_of_BIE_RMS',.01,'bie_gate',1e-4,'ps_gate',1e-4, ...
        'ps_step_gate',1e-10,'fk_gate',1e-5,'linear_gate',1e-10, ...
        'receiver_evaluation','direct at 411 physical coordinates; no receiver interpolation');
    k1_recovery_event(run_dir,'INPUT','END');
    k1_recovery_event(run_dir,'DIRECT','BEGIN');
    direct=baffled_strip_piston_field_vertical(x,z,source);
    k1_recovery_event(run_dir,'DIRECT','END');
    k1_recovery_event(run_dir,'FLAT','BEGIN');
    flat=-baffled_strip_piston_field_vertical(x,-z,source);
    k1_recovery_event(run_dir,'FLAT','END');
    recipe=[128 .05 .1875;256 .05 .1875;512 .05 .1875;1024 .05 .1875;1024 .025 .1875;1024 .025 .09375];
    PS=cell(size(recipe,1),1);
    for j=1:numel(PS)
        k1_recovery_event(run_dir,sprintf('PS%d',j),'BEGIN');
        PS{j}=local_ps(x,flat,source,profile,recipe(j,:));
        k1_recovery_event(run_dir,sprintf('PS%d',j),'END');
    end
    ps_values=[local_worst(PS{3}.field,PS{4}.field),local_worst(PS{4}.field,PS{5}.field),local_worst(PS{5}.field,PS{6}.field)];
    src=struct('frequency_hz',4000,'c0_mps',1500,'incident_field_fn', ...
        @(xx,zz)baffled_strip_piston_field_vertical(xx,zz,source));
    FK=cell(4,1); recipes=[400 50;800 50;960 60;1440 90];
    for j=1:4
        k1_recovery_event(run_dir,sprintf('FullFK%d',j),'BEGIN');
        FK{j}=closeout_s3a_pm_fullfk(src,profile,x,z,recipes(j,1),recipes(j,2),flat);
        k1_recovery_event(run_dir,sprintf('FullFK%d',j),'END');
    end
    fk_values=[local_worst(FK{1}.receiver_field,FK{2}.receiver_field), ...
        local_worst(FK{2}.receiver_field,FK{3}.receiver_field),local_worst(FK{3}.receiver_field,FK{4}.receiver_field)];
    sx=linspace(-50,50,2001)'; [eta,slope,curvature]=closeout_s3a_pm_surface(profile,sx);
    assert(all(z>closeout_s3a_pm_surface(profile,x)) && min(eta)>-2);
    cfg.ps_recipe=recipe; cfg.fk_recipe=recipes; cfg.U_PM_mps=frozen.U_PM_mps; cfg.modes=frozen.modes; cfg.case_id=frozen.case_id; cfg.physical_Kmax_rad_per_m=frozen.modes*2*pi/160; cfg.family=frozen.family;
    k1_recovery_event(run_dir,'COMMON_SAVE','BEGIN');
    save(fullfile(run_dir,'common.mat'),'cfg','profile','x','z','direct','flat','PS','FK','ps_values','fk_values','sx','eta','slope','curvature','-v7');
    k1_recovery_event(run_dir,'COMMON_SAVE','END');
    k1_recovery_event(run_dir,'MODEL_GATES','BEGIN');
    receipt.metrics=struct('ps_window',ps_values(1),'ps_grid',ps_values(2),'ps_step',ps_values(3), ...
        'fk_quadrature',fk_values(1),'fk_extent',max(fk_values(2:3)),'surface_min_z_m',min(eta),'surface_max_z_m',max(eta));
    assert(all(ps_values(1:2)<cfg.ps_gate) && ps_values(3)<cfg.ps_step_gate,'PS convergence gate');
    assert(all(fk_values<cfg.fk_gate),'FullFK convergence gate');
    k1_recovery_event(run_dir,'MODEL_GATES','END');
    k1_recovery_event(run_dir,'MODELS','END');
elseif strcmp(job_id,'analyze')
    q=load(fullfile(run_dir,'common.mat')); ids={'bie_coarse','bie_main','bie_self','bie_space','bie_window'};
    reference_qualified=strcmp(job_id,'analyze');
    if reference_qualified
        B=cell(5,1); for j=1:5, B{j}=load(fullfile(run_dir,'jobs',ids{j},'data.mat')); end
        bie=B{4}.field;
        bie_values=[local_worst(B{1}.field,B{2}.field),local_worst(B{2}.field,B{4}.field), ...
            local_worst(B{3}.field,B{2}.field),local_worst(B{2}.field,B{5}.field)];
        assert(all(bie_values<q.cfg.bie_gate),'BIE convergence gate');
        bunc=abs(B{1}.field-B{2}.field)+abs(B{2}.field-B{4}.field)+abs(B{3}.field-B{2}.field)+abs(B{2}.field-B{5}.field);
        reference_job='bie_space';
    else
        B=load(fullfile(run_dir,'jobs','bie_coarse','data.mat')); bie=B.field;
        bie_values=NaN(1,4); bunc=NaN(411,1); reference_job='bie_coarse';
    end
    U=[q.PS{6}.field bie q.FK{2}.receiver_field];
    punc=abs(q.PS{3}.field-q.PS{4}.field)+abs(q.PS{4}.field-q.PS{5}.field)+abs(q.PS{5}.field-q.PS{6}.field);
    func=abs(q.FK{1}.receiver_field-q.FK{2}.receiver_field)+abs(q.FK{2}.receiver_field-q.FK{3}.receiver_field)+abs(q.FK{3}.receiver_field-q.FK{4}.receiver_field);
    scale=norm(bie)/sqrt(411); mask=abs(bie)>=q.cfg.phase_mask_fraction_of_BIE_RMS*scale;
    points=table; metrics=table; names=["PS" "BIE" "FullFK"];
    for m=1:3
        e=U(:,m)-bie; amp=abs(U(:,m))-abs(bie); phase=angle(U(:,m).*conj(bie)); phase(~mask)=NaN;
        points=[points;table(repmat(names(m),411,1),q.x,q.z,real(U(:,m)),imag(U(:,m)), ...
            real(e),imag(e),abs(e),abs(e)/scale,amp,phase,mask, ...
            'VariableNames',{'method','x_m','z_m','real_Pa','imag_Pa','residual_real_Pa','residual_imag_Pa', ...
            'residual_abs_Pa','residual_over_BIE_RMS','amplitude_difference_Pa','phase_difference_rad','phase_valid'})]; %#ok<AGROW>
        num=bunc; if m==1, num=num+punc; elseif m==3, num=num+func; end
        E=norm(e)/norm(bie); res=5*norm(num)/norm(bie); [peak,jj]=max(abs(e));
        metrics=[metrics;table(names(m),E,res,E>res,peak,q.x(jj),norm(amp)/norm(bie),nnz(mask), ...
            'VariableNames',{'method','complex_relative_L2','five_times_empirical_resolution','resolved', ...
            'peak_residual_Pa','peak_x_m','amplitude_relative_L2','phase_valid_points'})]; %#ok<AGROW>
    end
    rough_flat=norm(bie-q.flat)/norm(q.flat);
    % Historical Gaussian source pairing is outside the R2 scope.
    writetable(points,fullfile(run_dir,'receiver_fields_and_errors.csv'));
    writetable(metrics,fullfile(run_dir,'model_metrics.csv'));
    save(fullfile(run_dir,'comparison.mat'),'U','bie_values','bunc','punc','func','scale','mask','metrics','rough_flat', ...
        'reference_qualified','reference_job','-v7');
    receipt.metrics=struct('PS_BIE_L2',metrics.complex_relative_L2(1),'FullFK_BIE_L2',metrics.complex_relative_L2(3), ...
        'BIE_rough_flat_L2',rough_flat,'BIE_worst_sensitivity',max(bie_values),'phase_valid_points',nnz(mask));
    receipt.reference_qualified=reference_qualified;
    if ~reference_qualified, receipt.unrun_checks={'BIE grid sensitivity','BIE self-quadrature sensitivity','BIE window sensitivity','five-times reference resolution'}; end
else
    error('Unknown K1 job');
end
receipt.status='COMPLETE'; receipt.elapsed_s=toc(t0); local_json(fullfile(out,'receipt.json'),receipt);
fprintf('%s\n',jsonencode(receipt));
catch ME
    receipt.status='FAILED_OR_STOPPED'; receipt.failed_checks=1; receipt.elapsed_s=toc(t0); receipt.error=ME.message;
    local_json(fullfile(out,'receipt.json'),receipt); fprintf(2,'%s\n',getReport(ME,'extended','hyperlinks','off')); rethrow(ME);
end
end
function r=local_ps(rx,flat,s,profile,rc)
W=rc(1); dx=rc(2); N=round(W/dx); xx=-W/2+(0:N-1)'*dx;
k=2*pi*s.frequency_hz/s.c_mps; eta=closeout_s3a_pm_surface(profile,xx);
inc=baffled_strip_piston_field_vertical(xx,0,s); delta=(1-exp(-2i*k*eta)).*inc;
kx=2*pi/W*[0:N/2-1 -N/2:-1]'; kz=sqrt(complex(k^2-kx.^2,0));
spec=fft(delta)/N; ns=ceil(3/rc(3)); ds=3/ns; half=exp(1i*(kz-k)*ds/2);
for j=1:ns, spec=half.*(half.*spec); end
spec=spec*exp(1i*k*3); field=flat;
% Direct Fourier evaluation, including evanescent spectrum, at actual rx.
for a=1:32:numel(rx), ii=a:min(a+31,numel(rx)); field(ii)=field(ii)+exp(1i*(rx(ii)+W/2)*kx.')*spec; end
keep=abs(xx)<=50;
r=struct('field',field,'window_m',W,'dx_m',dx,'step_m',ds,'kx',kx,'propagated_spectrum',spec, ...
    'surface_x_m',xx(keep),'surface_z_m',eta(keep),'incident_Pa',inc(keep),'delta_Pa',delta(keep), ...
    'operator','homogeneous square-root PE exp(i(kz-k)ds) with physical carrier', ...
    'reflection_model','PS: exact flat background plus (1-exp(-2ik eta_z))*p_inc(x,0)');
end
function v=local_worst(a,b)
v=norm(a-b)/max(norm(b),realmin); edges=round(linspace(1,numel(a),5));
for j=1:4, ii=edges(j):edges(j+1); v=max(v,norm(a(ii)-b(ii))/max(norm(b(ii)),realmin)); end
end
function local_json(p,v)
fid=fopen(p,'w'); assert(fid>0); c=onCleanup(@()fclose(fid)); fprintf(fid,'%s',jsonencode(v,PrettyPrint=true));
end
