function [checks,details]=audit_m0_static_phase_vertical(run_dir)
%AUDIT_M0_STATIC_PHASE_VERTICAL Reopen saved MAT/CSV, reconstruct independently.
% No propagation, external solver, renderer, or production phase helper used.
checks=struct('name',{},'value',{},'limit',{},'pass',{});
jobs={'A_algebra','B_zero','B_positive','B_negative','C_receipt_recovery', ...
    'C_direct_only','D_statistics','E_cir_comm_boundary_fixed','F_external_empirical', ...
    'G_contract_edges','H_metadata_guards','static_analysis_final','dashboard'};
details=struct('reopened_jobs',{jobs},'reconstructed_physical_cases',5, ...
    'historical_failures_retained',{{'C_gaussian receipt SHA failure','E_cir_comm legacy head BER failure'}});
export_dir=fullfile(run_dir,'lossless_metrics');
if ~isfolder(export_dir), mkdir(export_dir); end
for j=1:numel(jobs)
    q=load(fullfile(run_dir,jobs{j},'payload.mat'),'payload'); saved=q.payload;
    % MATLAB writetable's display decimals are not a lossless interchange.
    % Preserve those files; encode IEEE-754 hex in a new CSV, then decode it
    % and require exact equality. No physics/statistics tolerance is changed.
    export=fullfile(export_dir,[jobs{j} '.csv']);
    if ~isfile(export)
        fid=fopen(export,'w'); cleanup=onCleanup(@()fclose(fid));
        fprintf(fid,'name,value_decimal,value_ieee754_hex,limit_ieee754_hex,pass\n');
        for k=1:numel(saved.checks)
            c=saved.checks(k); fprintf(fid,'%s,%.17g,%s,%s,%d\n', ...
                c.name,c.value,num2hex(c.value),num2hex(c.limit),c.pass);
        end
        clear cleanup
    end
    import=detectImportOptions(export,'TextType','string','VariableNamingRule','preserve');
    import=setvartype(import,{'value_ieee754_hex','limit_ieee754_hex'},'string');
    csv=readtable(export,import);
    exact=hex2num(char(csv.value_ieee754_hex));
    checks(end+1)=gate([jobs{j} '_MAT_CSV'],max(abs(exact-[saved.checks.value]')),0); %#ok<AGROW>
    checks(end+1)=gate([jobs{j} '_recorded_gates'],sum(~[saved.checks.pass]),0); %#ok<AGROW>
end
q=load(fullfile(run_dir,'A_algebra','payload.mat'),'payload'); a=q.payload.data;
omega=2*pi*a.f; E=diag(exp(-1i*omega*(6/1500)));
native_total=a.native.direct_f+a.native.reflect_coh_f+a.native.reflect_sca_fm;
expected=conj(a.native.direct_f)+E*conj(a.native.reflect_coh_f+a.native.reflect_sca_fm);
checks(end+1)=gate('independent_complex_mapping',max(abs(expected-a.dsp.total_fm),[],'all'),1e-12);
checks(end+1)=gate('native_components_not_mutated',max(abs(native_total-(a.physical.direct_f.*exp(-1i*omega*(97/1500))+ ...
    a.physical.reflect_fm.*exp(-1i*omega*(103/1500)))),[],'all'),1e-12);
physical_jobs={'B_zero','B_positive','B_negative','C_receipt_recovery','C_direct_only'};
for j=1:numel(physical_jobs)
    q=load(fullfile(run_dir,physical_jobs{j},'payload.mat'),'payload'); a=q.payload.data; out=a.out;
    f=4062.5; k=2*pi*f/1500;
    if j<=3
        heights=[0,.01,-.01]; eta_z=-heights(j);
        pd=exp(1i*k*97); pr=-exp(1i*k*(103-2*eta_z));
    else
        n=256; w=50; sigma=.3; coordinates=(-n/2:n/2-1)*(w/n);
        [X,Y]=meshgrid(coordinates,coordinates); source=exp(-(X.^2+Y.^2)/(2*sigma^2));
        A=fft2(source); K=2*pi/w*[0:n/2-1,-n/2:-1]; [KX,KY]=meshgrid(K,K);
        kz=sqrt(complex(k^2-KX.^2-KY.^2,0));
        % Explicit inverse-DFT receiver evaluation; do not call the march or
        % reuse stored reference fields / their metadata.
        [I,J]=ndgrid(0:n-1,0:n-1); receiver_phase=(-1).^(I+J);
        pd=sum(A.*exp(1i*97*kz).*receiver_phase,'all')/n^2;
        pr=-sum(A.*exp(1i*103*kz).*receiver_phase,'all')/n^2;
        if j==5, pr=0; end
    end
    dsp=conj(pd+pr)*exp(1i*k*97);
    checks(end+1)=gate([physical_jobs{j} '_reconstructed_DSP'],abs(dsp-out.H_f)/max(abs(pd),eps),1e-10); %#ok<AGROW>
    checks(end+1)=gate([physical_jobs{j} '_reconstructed_physical'],abs(pd+pr-out.H_physical_f)/max(abs(pd),eps),1e-10); %#ok<AGROW>
end
q=load(fullfile(run_dir,'D_statistics','payload.mat'),'payload'); s=q.payload.data;
f=s.cfg.f_axis_hz; D=diag(exp(-1i*2*pi*f*.004)); X=D*conj(s.native_samples);
mu=sum(X,2)/size(X,2); Z=X-mu; C=Z*Z'/(size(X,2)-1); P=Z*Z.'/(size(X,2)-1);
checks(end+1)=gate('independent_sample_C',rel(C,s.migrated.stats.C_scatter_f),1e-12);
checks(end+1)=gate('independent_sample_P',rel(P,s.migrated.stats.P_scatter_f),1e-12);
checks(end+1)=gate('independent_old_DSP_migration',rel(C,s.old_migrated.stats.C_scatter_f),1e-12);
checks(end+1)=gate('independent_native_preservation',rel(s.forward.H_direct_reduced_f,s.flat.H_direct_reduced_f),1e-12);
q=load(fullfile(run_dir,'E_cir_comm_boundary_fixed','payload.mat'),'payload'); e=q.payload.data;
N=numel(e.f); time=(0:N-1)'/(N*62.5);
% Direct DFT synthesis, independent of build_channel_cir_vertical.
h=exp(1i*2*pi*(0:N-1)'*(0:N-1)/N)*e.H/N;
checks(end+1)=gate('independent_CIR_DFT',max(abs(h-e.cir.h_tau)),1e-12);
[~,peak]=max(abs(e.reflect_cir.h_tau));
checks(end+1)=gate('independent_arrival_s',abs(time(peak)-.004),1e-12);
checks(end+1)=gate('independent_final_production_BER',e.production_comm.noiseless_BER,0);
% Retain and independently explain the earlier communication failure using
% exactly the same stored physical response, symbols and observations.
old=load(fullfile(run_dir,'E_cir_comm','payload.mat'),'payload'); old=old.payload.data;
checks(end+1)=gate('same_physical_channel_after_boundary_fix',rel(e.production.H_f,old.production.H_f),0);
checks(end+1)=gate('old_failure_preserved',double(old.production_comm.noiseless_BER<=0),0);
% Resource and genuine status receipts; no process exit inferred for MCP.
budget=jsondecode(fileread(fullfile(run_dir,'budget_state.json')));
checks(end+1)=gate('MATLAB_budget_s',budget.used_s,600);
for j=1:numel(jobs)
    receipt=jsondecode(fileread(fullfile(run_dir,jobs{j},'checkpoint.json')));
    checks(end+1)=gate([jobs{j} '_job_time_s'],receipt.elapsed_s,180); %#ok<AGROW>
    checks(end+1)=gate([jobs{j} '_resource_gate'],~receipt.resource_gate_pass,0); %#ok<AGROW>
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
names={'vertical_channel_model','vertical_channel_model_impl','vertical_wape_propagator', ...
    'gaussian_source_initial_field_vertical','pm_surface_boundary_model', ...
    'apply_pe_channel_phase_reference_vertical'};
relative={'vertical_channel_model.m','src/channel/vertical_channel_model_impl.m', ...
    'src/propagation/vertical_wape_propagator.m','src/propagation/gaussian_source_initial_field_vertical.m', ...
    'src/surface/pm_surface_boundary_model.m','src/channel/apply_pe_channel_phase_reference_vertical.m'};
for j=1:numel(names)
    checks(end+1)=gate([names{j} '_path'],~strcmpi(which(names{j}),fullfile(root,relative{j})),0); %#ok<AGROW>
end
render=jsondecode(fileread(fullfile(root,'results','visualization','project_progress','render_receipt.json')));
checks(end+1)=gate('dashboard_references',render.failed_checks,0);
required_docs={'AGENTS.md','PROJECT_CONTEXT.md','vertical_comm_guide.md','src/README.md', ...
    'scripts/README.md','reports/project_progress.json','reports/project_progress.md', ...
    'reports/m0_phase_contract_closeout_2026-11-10.md'};
for j=1:numel(required_docs)
    checks(end+1)=gate(['document_' num2str(j)],~isfile(fullfile(root,required_docs{j})),0); %#ok<AGROW>
end
details.budget_at_audit=budget; details.total_checks=numel(checks);
end
function c=gate(name,value,limit)
c=struct('name',name,'value',double(value),'limit',limit,'pass',isfinite(value)&&value<=limit);
end
function value=rel(a,b)
value=norm(a(:)-b(:))/max(norm(b(:)),eps);
end
