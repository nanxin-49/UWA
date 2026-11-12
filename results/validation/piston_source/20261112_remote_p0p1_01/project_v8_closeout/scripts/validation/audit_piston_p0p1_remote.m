function receipt = audit_piston_p0p1_remote(task_root)
%AUDIT_PISTON_P0P1_REMOTE Independent saved-data reconstruction; remote only.
arguments
    task_root (1,:) char
end
run='20261112_remote_p0p1_01';
vr=fullfile(task_root,'results','validation','piston_source',run);
cr=fullfile(task_root,'results','comparisons','piston_2d_single_frequency',run);
out=fullfile(vr,'independent_audit'); assert(~isfolder(out)); mkdir(out);
t0=tic; rows=[]; resolution_rows=[]; sampling_rows=[];
refined_job='p0_pe_refined';
selection_file=fullfile(vr,'p0_recovery_selection.json');
if isfile(fullfile(vr,'p0_active_selection_v5.json'))
    selection_file=fullfile(vr,'p0_active_selection_v5.json');
end
if isfile(selection_file)
    selection=jsondecode(fileread(selection_file)); refined_job=selection.active_refined_job_id;
end
p0jobs={'p0_source',refined_job,'p0_bie','p0_fullfk'};
for j=1:numel(p0jobs)
    d=fullfile(vr,'jobs',p0jobs{j}); s=load(fullfile(d,'data.mat'),'T'); t=readtable(fullfile(d,'metrics.csv'));
    err=max(abs(t{:,:}-s.T{:,:}),[],'all'); assert(err<1e-9);
    r=jsondecode(fileread(fullfile(d,'receipt.json'))); assert(r.failed_checks==0 && contains(r.status,'PASS'));
    if j>=2
        stored=load(fullfile(d,'data.mat'),'inputs');
        for fi=1:3
            for wi=1:2
                v=stored.inputs{fi,wi};
                [p,px,pz]=baffled_strip_piston_field_vertical(v.x_m,v.z_m,v.source,128);
                assert(norm(v.p-p)/norm(p)<1e-7);
                assert(norm([v.px-px;v.pz-pz])/norm([px;pz])<1e-4);
            end
        end
    end
end
p0=jsondecode(fileread(fullfile(vr,'jobs','p0_finalize','receipt.json')));
assert(strcmp(p0.status,'P0_PASS'));
common=load(fullfile(vr,'common.mat')); gradient_rows=[]; launch_rows=[];
for frequency_hz=common.cfg.frequencies_hz
    for width_m=common.cfg.widths_m
        source=struct('frequency_hz',frequency_hz,'c_mps',1500,'rho_kgpm3',1025, ...
            'v0_mps',.001,'total_width_m',width_m,'z_source_m',100);
        launch_x=[-40;-20;-.2;0;.2;20;40];
        for offset=[1 2]
            continuous=piston_continuous_pe_march(launch_x,100-offset,source,1024,.09375,offset);
            aperture=baffled_strip_piston_field_vertical(launch_x,100-offset,source,128);
            pressure_error=norm(continuous-aperture)/norm(aperture);
            launch_rows=[launch_rows;frequency_hz width_m offset pressure_error]; %#ok<AGROW>
            assert(pressure_error<1e-7);
        end
        x=[-50;-42;-20;0;20;42;50];
        [surface_z,slope]=closeout_s2b_surface(x,.01);
        for plane=1:3
            if plane==1
                z=surface_z;
            elseif plane==2
                z=3+zeros(size(x));
            else
                z=98+zeros(size(x));
            end
            [~,px,pz]=baffled_strip_piston_field_vertical(x,z,source,128);
            h=1e-4;
            f=@(xx,zz)baffled_strip_piston_field_vertical(xx,zz,source,128);
            fx=(-f(x+2*h,z)+8*f(x+h,z)-8*f(x-h,z)+f(x-2*h,z))/(12*h);
            fz=(-f(x,z+2*h)+8*f(x,z+h)-8*f(x,z-h)+f(x,z-2*h))/(12*h);
            relative=norm([fx-px;fz-pz])/norm([px;pz]);
            normal=norm((-slope.*fx+fz)./hypot(1,slope)-(-slope.*px+pz)./hypot(1,slope))/norm([px;pz]);
            gradient_rows=[gradient_rows;frequency_hz width_m plane relative normal]; %#ok<AGROW>
            assert(relative<1e-6 && normal<1e-6);
        end
    end
end
gradient_table=array2table(gradient_rows,'VariableNames',{'frequency_hz','total_width_m','plane_id','gradient_relative','water_normal_gradient_relative'});
writetable(gradient_table,fullfile(out,'independent_boundary_receiver_launch_gradients.csv'));
launch_table=array2table(launch_rows,'VariableNames',{'frequency_hz','total_width_m','launch_offset_m','continuous_initial_pressure_relative'});
writetable(launch_table,fullfile(out,'continuous_source_launch_pressure.csv'));
stages={'free','flat','weak'}; stages_done={}; comparison_rows=[]; point_tables=cell(1,3);
for si=1:3
    stage=stages{si}; gate_job=['gate_cont_' stage '_r2'];
    gate=fullfile(cr,'jobs',gate_job,'receipt.json');
    if ~isfile(gate), continue; end
    gr=jsondecode(fileread(gate));
    s=load(fullfile(cr,'jobs',gate_job,'stage_summary.mat'),'selected','checks');
    selected=s.selected; check_values=cell2mat(struct2cell(s.checks));
    stages_done{end+1}=stage; %#ok<AGROW>
    for mi=1:3
        r=selected{mi}; d=fullfile(cr,'jobs',sprintf('%s_%s_l%d',stage,r.method,r.level));
        T=readtable(fullfile(d,'points.csv'));
        v=[T.direct_real_Pa+1i*T.direct_imag_Pa T.reflect_real_Pa+1i*T.reflect_imag_Pa T.total_real_Pa+1i*T.total_imag_Pa];
        csv_err=max(abs(v-[r.direct(:) r.reflect(:) r.total(:)]),[],'all');
        closure=max(abs(r.total-r.direct-r.reflect),[],'all');
        dsp_err=max(abs(r.H_total_dsp-conj(r.total)*exp(1i*2*pi*4000*97/1500)),[],'all');
        coarse=1:2:numel(r.x_m);
        field='reflect'; if strcmp(stage,'free'), field='direct'; end
        fine=r.(field); interpolated=interp1(r.x_m(coarse),fine(coarse,:),r.x_m,'spline');
        sampling=max(vecnorm(interpolated-fine)./max(vecnorm(fine),realmin));
        coarse_rows=ismember(T.x_m,r.x_m(coarse));
        writetable(T(coarse_rows,:),fullfile(out,[stage '_' r.method '_receiver_dx0p05.csv']));
        writetable(T,fullfile(out,[stage '_' r.method '_receiver_dx0p025.csv']));
        for wi=1:2
            [coarse_max,cix]=max(abs(fine(coarse,wi)));
            [fine_max,fix]=max(abs(fine(:,wi)));
            coarse_min=min(abs(fine(coarse,wi))); fine_min=min(abs(fine(:,wi)));
            sampling_rows=[sampling_rows;si mi wi .05 .025 coarse_max fine_max coarse_min fine_min ...
                r.x_m(coarse(cix)) r.x_m(fix) sampling]; %#ok<AGROW>
        end
        density_error=0; term_error=0;
        if isfield(r,'boundary')
            [density_error,weights]=local_bie_reconstruct(r.boundary);
            save(fullfile(out,[stage '_bie_arc_weights.mat']),'weights');
        end
        if isfield(r,'fullfk')
            for wi=1:2
                v=r.fullfk{wi}; term_error=max(term_error,max(abs(v.receiver_field-v.term_incident_dn-v.term_green_dn)));
            end
        end
        fail=csv_err>1e-10 || closure>1e-10 || dsp_err>1e-10 || sampling>1e-4 || density_error>1e-10 || term_error>1e-10;
        rows=[rows;si mi csv_err closure dsp_err sampling density_error term_error fail]; %#ok<AGROW>
    end
    % Physical references for free/flat; converged BIE for the weak surface.
    q=load(fullfile(vr,'common.mat'));
    if ~strcmp(stage,'free')
        empirical=zeros(3,2);
        pairs={{[3 4],[4 5],[5 6]},{[5 6],[6 7],[6 8]},{[2 3]}};
        for mi=1:3
            for pi=1:numel(pairs{mi})
                ids=pairs{mi}{pi};
                a=load(fullfile(cr,'jobs',sprintf('%s_%s_l%d',stage,selected{mi}.method,ids(1)),'data.mat'),'result');
                b=load(fullfile(cr,'jobs',sprintf('%s_%s_l%d',stage,selected{mi}.method,ids(2)),'data.mat'),'result');
                empirical(mi,:)=empirical(mi,:)+vecnorm(a.result.reflect-b.result.reflect)./max(vecnorm(selected{2}.reflect),realmin);
            end
        end
        for mi=1:3
            for wi=1:2
                model_delta=norm(selected{mi}.reflect(:,wi)-selected{2}.reflect(:,wi))/norm(selected{2}.reflect(:,wi));
                threshold=5*(empirical(mi,wi)+empirical(2,wi));
                resolution_rows=[resolution_rows;si mi wi empirical(mi,wi) empirical(2,wi) model_delta threshold model_delta>threshold all(check_values)]; %#ok<AGROW>
            end
        end
    end
    for wi=1:2
        for component={'direct','reflect','total'}
            name=component{1};
            if strcmp(stage,'free'), reference=q.p_direct(:,wi); if ~strcmp(name,'direct'), continue; end
            elseif strcmp(stage,'flat')
                if strcmp(name,'direct'), reference=q.p_direct(:,wi);
                elseif strcmp(name,'reflect'), reference=q.p_mirror(:,wi);
                else, reference=q.p_direct(:,wi)+q.p_mirror(:,wi); end
            else, reference=selected{2}.(name)(:,wi); end
            for mi=1:3
                value=selected{mi}.(name)(:,wi); delta=value-reference;
                relative=norm(delta)/max(norm(reference),realmin);
                max_abs=max(abs(delta)); valid=abs(reference)>=1e-4*max(abs(reference));
                adb=NaN(size(value)); phase=adb;
                adb(valid)=20*log10(abs(value(valid))./abs(reference(valid)));
                phase(valid)=angle(value(valid).*conj(reference(valid)))*180/pi;
                P=table(si+zeros(size(value)),mi+zeros(size(value)),q.cfg.widths_m(wi)+zeros(size(value)),q.cfg.rx_x_m, ...
                    real(reference),imag(reference),real(value),imag(value),abs(delta),adb,phase,valid, ...
                    'VariableNames',{'stage_id','model_id','total_width_m','x_m','reference_real_Pa','reference_imag_Pa', ...
                    'model_real_Pa','model_imag_Pa','absolute_error_Pa','amplitude_delta_dB','phase_delta_deg','valid_for_ratio'});
                writetable(P,fullfile(out,sprintf('%s_%s_%s_W%02d.csv',stage,selected{mi}.method,name,100*q.cfg.widths_m(wi))));
                for segment=0:4
                    mask=true(size(value));
                    if segment>0
                        edges=[-40 -20 0 20 40]; mask=q.cfg.rx_x_m>=edges(segment)&q.cfg.rx_x_m<=edges(segment+1);
                    end
                    seg_rel=norm(delta(mask))/max(norm(reference(mask)),realmin);
                    comparison_rows=[comparison_rows;si mi wi find(strcmp({'direct','reflect','total'},name)) segment relative max_abs seg_rel nnz(~valid)]; %#ok<AGROW>
                end
            end
        end
    end
    point_tables{si}=struct('stage',stage,'gate_status',gr.status,'all_gate_checks',all(check_values));
end
A=array2table(rows,'VariableNames',{'stage_id','model_id','MAT_CSV_max_abs_Pa','path_closure_max_abs_Pa', ...
    'DSP_v2_max_abs','receiver_spline_relative','BIE_density_reconstruction_relative','FullFK_term_closure_abs','failed'});
writetable(A,fullfile(out,'data_audit.csv'));
M=array2table(comparison_rows,'VariableNames',{'stage_id','model_id','width_id','component_id','segment_id', ...
    'whole_line_relative_L2','max_absolute_Pa','segment_relative_L2','weak_reference_point_count'});
writetable(M,fullfile(out,'model_comparison_metrics.csv'));
R=array2table(resolution_rows,'VariableNames',{'stage_id','model_id','width_id','method_empirical_uncertainty', ...
    'BIE_reference_empirical_uncertainty','method_vs_BIE_relative_L2','five_times_combined_uncertainty', ...
    'difference_above_empirical_resolution','reference_and_all_stage_gates_pass'});
writetable(R,fullfile(out,'empirical_resolution.csv'));
S=array2table(sampling_rows,'VariableNames',{'stage_id','model_id','width_id','coarse_dx_m','fine_dx_m', ...
    'coarse_max_peak_Pa','fine_max_peak_Pa','coarse_min_peak_Pa','fine_min_peak_Pa', ...
    'coarse_max_x_m','fine_max_x_m','coarse_spline_vs_actual_fine_relative'});
writetable(S,fullfile(out,'receiver_sampling_metrics.csv'));
receipt=struct('status','INDEPENDENT_SAVED_DATA_AUDIT_COMPLETE','elapsed_s',toc(t0), ...
    'P0_status',p0.status,'stages_audited',{stages_done},'failed_checks',sum(A.failed), ...
    'stage_gates',{point_tables},'claim_limit','Empirical numerical convergence; no P2 or 6/8 kHz propagation, PM, broadband, CIR or BER.');
save(fullfile(out,'audit.mat'),'receipt','A','M','R','S');
fid=fopen(fullfile(out,'receipt.json'),'w'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid);
fprintf('%s\n',jsonencode(receipt)); assert(receipt.failed_checks==0);
end
function [relative,weighted]=local_bie_reconstruct(b)
% Fresh GL nodes/weights and direct Green kernels, using saved density only.
n=b.panel_order; j=(1:n-1)'; beta=j./sqrt(4*j.^2-1); [V,D]=eig(diag(beta,1)+diag(beta,-1));
[~,ix]=sort(diag(D)); gw=2*(V(1,ix).^2)';
arc=(b.half_width_m/b.panel_count)*repmat(gw,b.panel_count,1)./b.node_normal(:,2);
weighted=arc.*b.window;
indices=unique(round(linspace(1,numel(b.receiver_x_m),17)));
x=b.receiver_x_m(indices); z=b.receiver_z_m(indices); xs=b.node_x_m.'; zs=b.node_z_m.';
nx=b.node_normal(:,1).'; nz=b.node_normal(:,2).'; k=b.k_radpm;
dx=x-xs; dz=z-zs; R=hypot(dx,dz);
dzm=z-(2*b.halfplane_h_m-zs); Rm=hypot(dx,dzm);
S=1i/4*(besselh(0,1,k*R)-besselh(0,1,k*Rm));
K=1i*k/4*(besselh(1,1,k*R).*(dx.*nx+dz.*nz)./R ...
    -besselh(1,1,k*Rm).*(dx.*nx-dzm.*nz)./Rm);
value=(K-1i*k*S)*(weighted.*b.density);
reference=b.receiver_field(indices,:);
relative=norm(value-reference,'fro')/max(norm(reference,'fro'),realmin);
end
