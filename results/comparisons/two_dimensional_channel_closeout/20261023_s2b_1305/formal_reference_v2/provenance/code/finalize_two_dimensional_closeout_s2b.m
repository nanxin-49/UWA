function summary=finalize_two_dimensional_closeout_s2b(out)
% Assemble saved fields only; no physical evaluation in this entrypoint.
assert(~isfile(fullfile(out,'data/s2b_data.mat')));
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
p=load(fullfile(out,'data/p0_data.mat')); f=p.fields;
t=load(fullfile(out,'data/p1.mat')); plan=t.plan;
grid=t.coverage(1:9,{'receiver_id','x_m','z_m'}); x=grid.x_m; z=grid.z_m;
methods={'PE','BIE','FullFK'}; cases=f.case_names; U=complex(nan(9,4,3)); provenance=strings(9,4,3);
strongfile=fullfile(root,'results/visualization/pe_bellhop_bie_corrected/20261021_native_reset_v2/plot_data/xz_data.mat');
st=load(strongfile,'data'); strong=st.data;
oldpath='20261023_s2a_final/data/s2a_data.mat';
for c=1:4
    for m=1:3
        for i=1:9
            if z(i)==3
                ix=find(abs(f.x_m-x(i))<1e-12); assert(isscalar(ix)); U(i,c,m)=f.(methods{m}){c}(ix); provenance(i,c,m)=oldpath;
            elseif c==4 && m<=2
                ix=find(strong.x_m==x(i)); iz=find(strong.z_m==z(i)); assert(isscalar(ix)&&isscalar(iz));
                if m==1, U(i,c,m)=strong.PE_carrier_restored(ix,iz); else, U(i,c,m)=strong.BIE_native(ix,iz); end
                provenance(i,c,m)=string(strongfile);
            end
        end
    end
end
assert(nnz(isfinite(U))==48); oldU=U;
receipts={}; jobs=struct; max_recheck=0;
for id={'bie3','flat','bie1','bie2','pe1','pe2','pe3','pe4','fk1','fk2','fk3','fk4'}
    a=load(fullfile(out,'checkpoints',[id{1} '.mat'])); jobs.(id{1})=a.payload;
    a.receipt.job_id=string(a.receipt.job_id); a.receipt.memory_measurement=string(a.receipt.memory_measurement);
    a.receipt.timestamp=string(a.receipt.timestamp); receipts{end+1,1}=struct2table(a.receipt); %#ok<AGROW>
    if strcmp(id{1},'flat'), c=1;m=2;
    elseif startsWith(id{1},'bie'), c=str2double(id{1}(end))+1;m=2;
    elseif startsWith(id{1},'pe'), c=str2double(id{1}(end));m=1;
    else, c=str2double(id{1}(end));m=3;
    end
    val=a.payload.u(1:9); present=isfinite(U(:,c,m));
    if any(present), max_recheck=max(max_recheck,max(abs(val(present)-U(present,c,m)))); end
    missing=~present; U(missing,c,m)=val(missing); provenance(missing,c,m)=string(['checkpoints/' id{1} '.mat']);
end
assert(all(isfinite(U),'all') && nnz(isfinite(U))==108 && max_recheck<1e-10);
lineage={}; points={}; rows={}; diagrows={}; crossrows={}; sens={}; biesens={};
pairs=[1 2;3 2;1 3];
for c=1:4
    for m=1:3
        for i=1:9
            u=U(i,c,m); lineage{end+1,1}=table(string(cases{c}),string(methods{m}),4000,grid.receiver_id(i),x(i),z(i), ...
                real(u),imag(u),real(u),imag(u),isfinite(oldU(i,c,m)),provenance(i,c,m), ...
                "u_reflected_native","GAUSS_X_4K_UNIT_PEAK_N984","physical_map", ...
                "S2a-depth-positive-native-v1","native exp(-i omega t); outgoing exp(+ikR)", ...
                "unit-source response; not Pa","no calibration", ...
                'VariableNames',{'case_name','method','frequency_hz','receiver_id','x_m','z_m', ...
                'u_native_real','u_native_imag','u_mapped_real','u_mapped_imag','reused_existing', ...
                'provenance','quantity_kind','source_id','source_tier','source_map_version','phase_convention','unit','calibration'}); %#ok<AGROW>
            if c>1
                flat=U(i,1,m); g=u/flat; delta=u-flat;
                diagrows{end+1,1}=table(string(cases{c}),string(methods{m}),grid.receiver_id(i),real(g),imag(g),real(delta),imag(delta),abs(flat)^2, ...
                    'VariableNames',{'case_name','method','receiver_id','G_real','G_imag','delta_real','delta_imag','flat_denominator_squared'}); %#ok<AGROW>
            end
        end
    end
    for j=1:3
        ca=pairs(j,1); rb=pairs(j,2); a=U(:,c,ca); b=U(:,c,rb);
        m=closeout_s2b_metrics(a,b,ones(9,1)/9);
        row=addvars(struct2table(m),string(cases{c}),string(methods{ca}),string(methods{rb}), ...
            "spatial_uniform_9","P11-P33","E_array",grid.receiver_id(m.worst_index), ...
            'Before',1,'NewVariableNames',{'case_name','candidate','reference','weight_kind','support_id','metric_kind','worst_receiver'});
        rows{end+1,1}=row; %#ok<AGROW>
        for i=1:9
            d=abs(a(i)-b(i)); ph=angle(a(i)*conj(b(i)));
            points{end+1,1}=table(string(cases{c}),string(methods{ca}),string(methods{rb}),grid.receiver_id(i),x(i),z(i), ...
                real(a(i)),imag(a(i)),real(b(i)),imag(b(i)),d,d/abs(b(i)),20*log10(abs(a(i))/abs(b(i))),ph, ...
                (abs(a(i))-abs(b(i)))^2,4*abs(a(i))*abs(b(i))*sin(ph/2)^2,abs(b(i))^2, ...
                'VariableNames',{'case_name','candidate','reference','receiver_id','x_m','z_m','a_real','a_imag','b_real','b_imag', ...
                'absolute_error','relative_error','db_signed','phase_signed_rad','magnitude_squared','phase_squared','denominator_squared'}); %#ok<AGROW>
        end
        if c>1
            ef=U(:,1,ca)-U(:,1,rb); er=a-b; ed=er-ef; den=mean(abs(U(:,1,rb)).^2);
            terms=[mean(abs(ef).^2),mean(abs(ed).^2),2*real(mean(ef.*conj(ed))),mean(abs(er).^2)]/den;
            crossrows{end+1,1}=table(string(cases{c}),string(methods{ca}),string(methods{rb}),den,terms(1),terms(2),terms(3),terms(4), ...
                terms(4)-sum(terms(1:3)),'VariableNames',{'case_name','candidate','reference','flat_denominator','flat_term','delta_term','cross_term','rough_term','closure'}); %#ok<AGROW>
        end
    end
    pe=jobs.(sprintf('pe%d',c)); fk=jobs.(sprintf('fk%d',c));
    for j=1:4
        if j==1, a=pe.step_half_line; b=f.PE{c}; label='PE_step_0.05_to_0.025_line411';
        elseif j==2, a=pe.grid_double_u; b=pe.u; label='PE_N984_to1968_same_spectrum_probes';
        elseif j==3, a=fk.runs{3}.receiver_field; b=fk.runs{2}.receiver_field; label='FullFK_L60_4097_to8193';
        else, a=fk.runs{4}.receiver_field; b=fk.runs{3}.receiver_field; label='FullFK_L60_to90_same_fine_spacing';
        end
        sens{end+1,1}=table(string(cases{c}),string(label),numel(a),norm(a-b)/norm(a),max(abs(a-b)), ...
            'VariableNames',{'case_name','check_kind','n_probes','relative_l2','max_abs_change'}); %#ok<AGROW>
    end
    % Preserve both consecutive quadrature changes, never infer convergence from one difference.
    a=fk.runs{2}.receiver_field; b=fk.runs{1}.receiver_field;
    sens{end+1,1}=table(string(cases{c}),"FullFK_L60_2049_to4097",numel(a),norm(a-b)/norm(a),max(abs(a-b)), ...
        'VariableNames',{'case_name','check_kind','n_probes','relative_l2','max_abs_change'}); %#ok<AGROW>
    if c>1
        bie=jobs.(sprintf('bie%d',c-1)); lv=bie.level_fields;
        sp=abs(lv.spatial_runs(:,end)-lv.spatial_runs(:,end-1));
        wi=max(abs(diff(lv.window_runs,1,2)),[],2);
        qu=nan(size(sp)); if isfield(lv,'quadrature_runs'), qu=abs(diff(lv.quadrature_runs,1,2)); end
        for i=1:numel(sp)
            envelope=sp(i)+wi(i); if isfinite(qu(i)), envelope=envelope+qu(i); end
            biesens{end+1,1}=table(string(cases{c}),i,bie.x_m(i),bie.z_m(i),sp(i),wi(i),qu(i),envelope,abs(bie.u(i)), ...
                "empirical_known_components_only","REFERENCE_LIMITED", ...
                'VariableNames',{'case_name','probe_index','x_m','z_m','spatial_abs_change','window_abs_change','quadrature_abs_change', ...
                'known_component_envelope','reference_amplitude','uncertainty_kind','reference_status'}); %#ok<AGROW>
        end
    end
end
fields=vertcat(lineage{:}); array_metrics=vertcat(rows{:}); point_metrics=vertcat(points{:});
diagnostics=vertcat(diagrows{:}); increment=vertcat(crossrows{:}); numerical=vertcat(sens{:}); reference=vertcat(biesens{:});
cost=vertcat(receipts{:});
summary=struct('status','S2B_DATA_COMPLETE_REFERENCE_LIMITED','P0_pass',p.audit.pass,'coverage',108,'target',108, ...
    'new_values',60,'reused_values',48,'max_existing_recheck_abs',max_recheck, ...
    'new_measured_core_seconds',sum(cost.wall_seconds),'S2a_measured_core_seconds',plan.prior_recorded_seconds, ...
    'historical_unmetered_reserve_seconds',120,'failed_attempt_conservative_charge_seconds',30, ...
    'charged_total_seconds',sum(cost.wall_seconds)+plan.prior_recorded_seconds+150, ...
    'reference_limitation','Missing original execution-time BIE code hashes and medium/strong independent self-quadrature layers; saved spatial/window sensitivities are not rigorous bounds.', ...
    'G1_target_PM','OPEN','WP2_formal',false,'new_BIE_solves',0,'new_Bellhop_jobs',0,'S3_S6_executed',false);
save(fullfile(out,'data/s2b_data.mat'),'U','oldU','provenance','grid','fields','array_metrics','point_metrics', ...
    'diagnostics','increment','numerical','reference','cost','summary','-v7.3');
for name={'fields','array_metrics','point_metrics','diagnostics','increment','numerical','reference','cost'}
    switch name{1}
        case 'fields', tab=fields;
        case 'array_metrics', tab=array_metrics;
        case 'point_metrics', tab=point_metrics;
        case 'diagnostics', tab=diagnostics;
        case 'increment', tab=increment;
        case 'numerical', tab=numerical;
        case 'reference', tab=reference;
        case 'cost', tab=cost;
    end
    writetable(tab,fullfile(out,'data',[name{1} '.csv']));
end
closeout_s2b_json(fullfile(out,'data/summary.json'),summary);
plot_data=struct('U',U,'grid',grid,'line',f,'line_metrics',p.metrics,'array_metrics',array_metrics, ...
    'strong',struct('x_m',strong.x_m,'z_m',strong.z_m,'PE',strong.PE_carrier_restored,'BIE',strong.BIE_native), ...
    'summary',summary);
save(fullfile(out,'data/plot_data.mat'),'plot_data','-v7.3'); disp(summary); disp(array_metrics(:,[1:3 8 25]));
end
