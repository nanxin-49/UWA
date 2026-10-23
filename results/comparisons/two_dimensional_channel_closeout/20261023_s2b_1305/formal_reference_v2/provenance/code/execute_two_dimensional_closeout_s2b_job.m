function receipt=execute_two_dimensional_closeout_s2b_job(out,job_id,available_physical_bytes)
% One checkpointed P2 batch, existing MATLAB session only.
assert(~isempty(regexp(job_id,'^(bie[123]|flat|pe[1234]|fk[1234])$','once')));
dest=fullfile(out,'checkpoints'); if ~isfolder(dest), mkdir(dest); end
file=fullfile(dest,[job_id '.mat']); assert(~isfile(file),'Job already checkpointed');
t=load(fullfile(out,'data/p1.mat')); plan=t.plan;
p=load(fullfile(out,'data/p0_data.mat'),'fields'); f=p.fields;
assert(strcmp(plan.status,'P1_INPUT_RESOURCE_GATE_PASS_DATA_EVALUATION_ONLY'));
previous=dir(fullfile(dest,'*.mat')); spent=plan.prior_recorded_seconds+plan.prior_unmetered_reserve_seconds;
for j=1:numel(previous), old=load(fullfile(dest,previous(j).name),'receipt'); spent=spent+old.receipt.wall_seconds; end
estimated_peak=96*1024^2; assert(estimated_peak<=.7*available_physical_bytes,'Physical RAM gate');
assert(spent+900<7200,'Cumulative budget gate');
if ~strcmp(job_id,'bie3'), assert(isfile(fullfile(dest,'bie3.mat')),'Representative density job must precede batches'); end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
infile=fullfile(root,'results/visualization/pe_bellhop_bie_xz/20261019_115214_348/source_data/numeric_inputs.mat');
s=load(infile,'inputs'); in=s.inputs; src=f.source;
[xx,zz]=meshgrid([-.1953125 0 .1953125],[2.875 3 3.125]);
xx=reshape(xx.',[],1); zz=reshape(zz.',[],1); ids=string(compose('P%d%d',repelem((1:3)',3),repmat((1:3)',3,1)));
% Frozen probes: nine points plus each method-pair's worst line location
% and minimum reference amplitude. Geometry is selected from saved P0 only.
if strcmp(job_id,'flat'), c=1; else, c=str2double(job_id(end)); if startsWith(job_id,'bie'), c=c+1; end, end
extra=[]; pairs={'PE','BIE';'FullFK','BIE';'PE','FullFK'};
for j=1:3
    a=f.(pairs{j,1}){c}; b=f.(pairs{j,2}){c};
    [~,im]=max(abs(a-b)./abs(b)); [~,imin]=min(abs(b)); extra=[extra;im;imin]; %#ok<AGROW>
end
extra=unique(extra); xp=[xx;f.x_m(extra)]; zp=[zz;3+zeros(size(extra))];
A=[0 .01 .05 .2]; payload=struct('case_index',c,'x_m',xp,'z_m',zp,'receiver_ids',ids, ...
    'line_probe_index',extra,'source',src,'quantity_kind','u_reflected_native', ...
    'representation','native exp(-i omega t); outgoing exp(+ikR)', ...
    'source_id','GAUSS_X_4K_UNIT_PEAK_N984','unit','unit-source response; not Pa');
ticjob=tic;
if startsWith(job_id,'bie')
    ss=t.snapshots{c-1}; names=fieldnames(ss); vals=struct;
    for j=1:numel(names)
        levels=ss.(names{j}); U=complex(zeros(numel(xp),numel(levels)));
        for k=1:numel(levels), U(:,k)=evaluate_saved_halfplane_bie_field(levels{k},xp,zp); end
        vals.(names{j})=U;
    end
    payload.level_fields=vals; payload.u=vals.spatial_runs(:,end);
    [ok,ix]=ismember(xp(zp==3),f.x_m); assert(all(ok));
    closure=norm(payload.u(zp==3)-f.BIE{c}(ix))/norm(f.BIE{c}(ix));
    assert(closure<1e-11); payload.saved_line_closure=closure;
elseif strcmp(job_id,'flat')
    payload.u=-fk_single_receiver_incident(src,xp,-zp);
    [ok,ix]=ismember(xp(zp==3),f.x_m); assert(all(ok));
    assert(norm(payload.u(zp==3)-f.BIE{1}(ix))/norm(f.BIE{1}(ix))<1e-11);
elseif startsWith(job_id,'pe')
    xall=in.xfull_m(:); k=2*pi*4000/1500;
    eta=closeout_s2b_surface(xall,A(c));
    corrected_start=-exp(-2i*k*eta).*in.pe.surface_incident_field(:);
    payload.corrected_start=corrected_start; payload.u=complex(zeros(numel(xp),1));
    for iz=unique(zp).'
        q=ifft(fft(corrected_start).*exp(1i*iz*(src.kz-k)))*exp(1i*k*(100+iz));
        rows=zp==iz; [ok,ix]=ismember(xp(rows),xall); assert(all(ok)); payload.u(rows)=q(ix);
    end
    [ok,ix]=ismember(xp(zp==3),f.x_m); assert(all(ok));
    assert(norm(payload.u(zp==3)-f.PE{c}(ix))/norm(f.PE{c}(ix))<1e-10);
    cfg=struct('frequency_hz',4000,'c0_mps',1500,'xw_m',192.1875,'nx',984, ...
        'z_tx_m',100,'z_rx_m',3,'sigma_src_m',.3,'surface_reflect_coeff',-1, ...
        'step_m',.025,'x_rx_m',0,'reflection_model','model0_normal', ...
        'source_field_override',in.pe.source_field(:),'eta_z_m',eta);
    q=run_pe_physical_surface_model0_validation(cfg);
    halfstep=q.reflected_native(:); payload.step_half_line=halfstep(in.original_M99);
    payload.step_half_difference=norm(payload.step_half_line-f.PE{c})/norm(f.PE{c});
    % Double numerical grid, retaining the SAME frozen Fourier source.
    cfg.nx=1968; cfg.step_m=.05; xf=-192.1875/2+(0:1967)'*192.1875/1968;
    coeff=complex(zeros(1968,1)); coeff(1:492)=src.coeff(1:492)*2;
    coeff(end-491:end)=src.coeff(493:end)*2;
    cfg.source_field_override=ifft(coeff); cfg.eta_z_m=closeout_s2b_surface(xf,A(c));
    assert(max(abs(cfg.source_field_override(1:2:end)-in.pe.source_field(:)))<1e-12);
    ug=complex(zeros(size(xp)));
    for iz=unique(zp).'
        cfg.z_rx_m=iz; q=run_pe_physical_surface_model0_validation(cfg);
        rows=zp==iz; [ok,ix]=ismember(xp(rows),xf); assert(all(ok)); ur=q.reflected_native(:); ug(rows)=ur(ix);
    end
    payload.grid_double_u=ug; payload.grid_double_source_note='zero padded frozen N984 spectrum, same W; no source replacement';
elseif startsWith(job_id,'fk')
    sizes=[2049 4097 8193 12289]; L=[60 60 60 90]; runs=cell(4,1);
    for j=1:4
        xs=linspace(-L(j),L(j),sizes(j))'; [eta,ep]=closeout_s2b_surface(xs,A(c));
        nv=[-ep./hypot(1,ep),1./hypot(1,ep)];
        [ui,ux,uz]=fk_single_receiver_incident(src,xs,eta);
        cfg=struct('frequency_hz',4000,'c0_mps',1500,'surface_x_m',xs,'surface_z_m',eta, ...
            'surface_normal',nv,'surface_incident',ui,'surface_incident_dn',nv(:,1).*ux+nv(:,2).*uz, ...
            'receiver_x_m',xp,'receiver_z_m',zp,'quadrature_measure','arc_length','local_radius_m',Inf);
        runs{j}=solve_full_kirchhoff_2d(cfg);
    end
    payload.runs=runs; payload.u=runs{2}.receiver_field;
    [ok,ix]=ismember(xp(zp==3),f.x_m); assert(all(ok));
    assert(norm(payload.u(zp==3)-f.FullFK{c}(ix))/norm(f.FullFK{c}(ix))<1e-10);
end
seconds=toc(ticjob); assert(seconds<900,'Job exceeded limit; no further jobs');
receipt=struct('job_id',job_id,'wall_seconds',seconds,'available_physical_bytes',available_physical_bytes, ...
    'memory_measurement','Windows Memory Available Bytes immediately before MCP call', ...
    'estimated_peak_increment_bytes',estimated_peak,'reserve_fraction',.3, ...
    'cumulative_charged_seconds',spent+seconds,'n_probe_values',numel(xp), ...
    'timestamp',char(datetime('now','TimeZone','Asia/Shanghai')), ...
    'field_evaluation_is_new_physics',true,'new_boundary_solves',0);
save(file,'payload','receipt','-v7.3');
closeout_s2b_json(fullfile(dest,[job_id '.json']),receipt); disp(receipt);
end
