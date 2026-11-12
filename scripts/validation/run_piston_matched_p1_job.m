function receipt = run_piston_matched_p1_job(job_id,run_root)
%RUN_PISTON_MATCHED_P1_JOB Uniform-medium PE with complete incident exterior.
% Existing truncated-periodic trials are retained as failed diagnostics.
arguments
    job_id (1,:) char
    run_root (1,:) char
end
out=fullfile(run_root,'jobs',job_id); assert(~isfolder(out)); mkdir(out); t0=tic;
receipt=struct('job_id',job_id,'status','RUNNING','pid',feature('getpid'),'matlab_version',version,'failed_checks',0);
local_json(fullfile(out,'receipt.json'),receipt);
try
    q=load(fullfile(run_root,'common.mat')); cfg=q.cfg;
    parts=regexp(job_id,'^(free|flat|weak)_pe_cont_l(\d+)$','tokens','once');
    if ~isempty(parts)
        stage=parts{1}; level=str2double(parts{2});
        result=struct('stage',stage,'method','pe_cont','level',level,'frequency_hz',4000, ...
            'widths_m',cfg.widths_m,'x_m',cfg.rx_x_m,'z_m',3, ...
            'representation','absolute_physical peak Pa, exp(-i omega t)', ...
            'phase_contract_version','2.0.0','direct',q.p_direct,'reflect',complex(zeros(size(q.p_direct))));
        if strcmp(stage,'free')
            recipe=[128 .1875 1;256 .1875 1;512 .1875 1;512 .09375 1;512 .09375 2; ...
                1024 .09375 2;2048 .09375 2;1024 .1875 2;1024 .09375 1];
            rc=recipe(level,:); info=cell(1,2);
            for wi=1:2
                source=local_source(cfg,cfg.widths_m(wi));
                [result.direct(:,wi),info{wi}]=piston_continuous_pe_march(cfg.rx_x_m,3,source,rc(1),rc(2),rc(3));
            end
            result.numeric_config=struct('modes_per_branch',rc(1),'step_m',rc(2),'offset_m',rc(3));
            result.continuum_meta=info;
        else
            base=load(fullfile(run_root,'jobs','free_pe_cont_l7','data.mat'),'result');
            result.direct=base.result.direct;
            recipe=[128 .05 .1875;256 .05 .1875;512 .05 .1875;1024 .05 .1875; ...
                1024 .025 .1875;1024 .025 .09375];
            rc=recipe(level,:); W=rc(1); dx=rc(2); N=round(W/dx);
            x=-W/2+(0:N-1)'*dx; A=0; if strcmp(stage,'weak'), A=.01; end
            [eta,slope]=closeout_s2b_surface(x,A);
            result.numeric_config=struct('window_m',W,'dx_m',dx,'step_m',rc(3),'incident_exterior','continuous full-aperture background');
            delta_fields=cell(1,2); surface_inputs=delta_fields;
            for wi=1:2
                source=local_source(cfg,cfg.widths_m(wi));
                incident=baffled_strip_piston_field_vertical(x,0,source);
                k=2*pi*4000/1500;
                % Exact flat background plus compact Model-0 screen change.
                delta=(1-exp(-2i*k*eta)).*incident;
                kz=sqrt(complex(k^2-(2*pi/W*[0:N/2-1 -N/2:-1]').^2,0));
                spectrum=fft(delta); ns=max(1,ceil(3/rc(3))); ds=3/ns;
                half=exp(1i*(kz-k)*ds/2);
                for j=1:ns, spectrum=half.*(half.*spectrum); end
                field=ifft(spectrum)*exp(1i*k*3);
                factor=round(dx/.025); field_fine=interpft(field,N*factor);
                ix=round((cfg.rx_x_m+W/2)/.025)+1;
                result.reflect(:,wi)=q.p_mirror(:,wi)+field_fine(ix);
                delta_fields{wi}=field;
                surface_inputs{wi}=struct('x_m',x,'physical_surface_z_m',eta,'surface_slope',slope, ...
                    'incident_at_reference_plane_Pa',incident,'delta_screen_pressure_Pa',delta);
            end
            result.compact_delta_fields=delta_fields; result.surface_inputs=surface_inputs;
        end
        result.total=result.direct+result.reflect;
        result.H_direct_dsp=conj(result.direct)*exp(1i*2*pi*4000*97/1500);
        result.H_reflect_dsp=conj(result.reflect)*exp(1i*2*pi*4000*97/1500);
        result.H_total_dsp=result.H_direct_dsp+result.H_reflect_dsp;
        save(fullfile(out,'data.mat'),'result','cfg','-v7'); local_points(result,out);
        receipt.status='MATCHED_PE_DATA_COMPLETE';
    elseif startsWith(job_id,'gate_cont_')
        parts=regexp(job_id,'^gate_cont_(free|flat|weak)','tokens','once'); assert(~isempty(parts));
        stage=parts{1}; local_gate(stage,run_root,out,q,receipt);
        r=jsondecode(fileread(fullfile(out,'receipt.json'))); receipt=r;
    else
        error('Unknown matched P1 job.');
    end
    receipt.elapsed_s=toc(t0); local_json(fullfile(out,'receipt.json'),receipt); fprintf('%s\n',jsonencode(receipt));
catch ME
    receipt.status='FAILED_OR_INCONCLUSIVE'; receipt.failed_checks=1; receipt.error=ME.message; receipt.elapsed_s=toc(t0);
    local_json(fullfile(out,'receipt.json'),receipt); fid=fopen(fullfile(out,'matlab_error.txt'),'w'); fprintf(fid,'%s',getReport(ME,'extended','hyperlinks','off')); fclose(fid); rethrow(ME);
end
end
function s=local_source(cfg,W)
s=struct('frequency_hz',4000,'c_mps',1500,'rho_kgpm3',1025,'v0_mps',.001,'total_width_m',W,'z_source_m',cfg.z_source_m);
end
function local_points(r,out)
n=numel(r.x_m);
T=table(repmat(r.x_m,2,1),repelem(r.widths_m(:),n),real(r.direct(:)),imag(r.direct(:)), ...
    real(r.reflect(:)),imag(r.reflect(:)),real(r.total(:)),imag(r.total(:)),abs(r.reflect(:)), ...
    'VariableNames',{'x_m','total_width_m','direct_real_Pa','direct_imag_Pa','reflect_real_Pa','reflect_imag_Pa','total_real_Pa','total_imag_Pa','reflect_abs_Pa'});
writetable(T,fullfile(out,'points.csv'));
end
function local_gate(stage,root,out,q,receipt)
if strcmp(stage,'free'), bestlevel=7; field='direct'; else, bestlevel=6; field='reflect'; end
a=load(fullfile(root,'jobs',sprintf('%s_pe_cont_l%d',stage,bestlevel),'data.mat'),'result'); pe=a.result;
checks=struct; values=struct;
if strcmp(stage,'free')
    values.spectral_grid=local_difference(root,stage,6,7,field);
    values.step=local_difference(root,stage,6,8,field);
    values.launch_offset=local_difference(root,stage,6,9,field);
    values.physical=local_worst(pe.direct,q.p_direct);
    checks.spectral_grid=values.spectral_grid<q.cfg.gates.pe_grid_relative;
    checks.step=values.step<q.cfg.gates.pe_step_relative;
    checks.launch_offset=values.launch_offset<q.cfg.gates.pe_grid_relative;
    checks.physical=values.physical<q.cfg.gates.physical_free_flat_relative;
else
    values.window=local_difference(root,stage,3,4,field);
    values.grid=local_difference(root,stage,4,5,field);
    values.step=local_difference(root,stage,5,6,field);
    checks.window=values.window<q.cfg.gates.pe_window_relative;
    checks.grid=values.grid<q.cfg.gates.pe_grid_relative;
    checks.step=values.step<q.cfg.gates.pe_step_relative;
    if strcmp(stage,'flat')
        values.physical=local_worst(pe.reflect,q.p_mirror); checks.mirror=values.physical<q.cfg.gates.physical_free_flat_relative;
    end
end
b=dir(fullfile(root,'jobs',[stage '_bie_l*'],'data.mat')); f=dir(fullfile(root,'jobs',[stage '_fullfk_l*'],'data.mat'));
bie=local_latest(b); fk=local_latest(f); selected={pe,bie,fk};
if strcmp(stage,'free')
    checks.bie=local_worst(bie.direct,q.p_direct)<1e-7; checks.fk=local_worst(fk.direct,q.p_direct)<1e-7;
else
    checks.bie_residual=max(bie.boundary.offgrid_boundary_residual)<q.cfg.gates.bie_boundary_relative;
    values.bie_grid=local_existing_pair(root,stage,'bie',5,6,field);
    values.bie_window=local_existing_pair(root,stage,'bie',6,7,field);
    values.bie_quadrature=local_existing_pair(root,stage,'bie',6,8,field);
    values.fk_quadrature=local_existing_pair(root,stage,'fullfk',2,3,field);
    checks.bie_grid=values.bie_grid<q.cfg.gates.bie_convergence_relative;
    checks.bie_window=values.bie_window<q.cfg.gates.bie_convergence_relative;
    checks.bie_quadrature=values.bie_quadrature<q.cfg.gates.bie_convergence_relative;
    checks.fk_quadrature=values.fk_quadrature<q.cfg.gates.fk_convergence_relative;
    if strcmp(stage,'flat')
        values.bie_mirror=local_worst(bie.reflect,q.p_mirror); checks.bie_mirror=values.bie_mirror<1e-3;
        checks.fk_independent_plane=fk.flat_independent_probe.relative<1e-3;
    end
end
receipt.values=values; receipt.checks=checks; receipt.failed_checks=nnz(~cell2mat(struct2cell(checks)));
receipt.status='P1_STAGE_PASS'; if receipt.failed_checks>0, receipt.status='P1_STAGE_INCONCLUSIVE'; end
M=table; save(fullfile(out,'stage_summary.mat'),'selected','checks','values','M','-v7');
local_json(fullfile(out,'receipt.json'),receipt); local_json(fullfile(out,'checks.json'),checks);
end
function r=local_latest(files)
assert(~isempty(files)); levels=zeros(size(files));
for j=1:numel(files), a=load(fullfile(files(j).folder,files(j).name),'result'); levels(j)=a.result.level; end
[~,ix]=max(levels); a=load(fullfile(files(ix).folder,files(ix).name),'result'); r=a.result;
end
function v=local_difference(root,stage,a,b,field)
v=local_existing_pair(root,stage,'pe_cont',a,b,field);
end
function v=local_existing_pair(root,stage,method,a,b,field)
x=load(fullfile(root,'jobs',sprintf('%s_%s_l%d',stage,method,a),'data.mat'),'result');
y=load(fullfile(root,'jobs',sprintf('%s_%s_l%d',stage,method,b),'data.mat'),'result');
v=local_worst(x.result.(field),y.result.(field));
end
function v=local_worst(A,B)
v=max(vecnorm(A-B)./max(vecnorm(B),realmin)); edges=round(linspace(1,size(A,1),5));
for j=1:4, ix=edges(j):edges(j+1); v=max(v,max(vecnorm(A(ix,:)-B(ix,:))./max(vecnorm(B(ix,:)),realmin))); end
end
function local_json(p,v)
fid=fopen(p,'w'); assert(fid>0); fprintf(fid,'%s',jsonencode(v,PrettyPrint=true)); fclose(fid);
end
