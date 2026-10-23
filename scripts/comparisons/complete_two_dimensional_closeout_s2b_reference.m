function summary=complete_two_dimensional_closeout_s2b_reference(out,available_physical_bytes)
% S2b reference closure: known-version BIE layers and independent full FK.
dest=fullfile(out,'formal_reference_v2'); assert(~isfolder(dest));
assert(192*1024^2<=.7*available_physical_bytes,'Physical RAM gate');
initial=jsondecode(fileread(fullfile(out,'data/cost.json'))); spent=initial.charged_total_seconds;
files=dir(fullfile(out,'reference_followup_01/bie_*.mat')); bcost=cell(numel(files),1);
for j=1:numel(files)
    s=load(fullfile(files(j).folder,files(j).name),'receipt'); spent=spent+s.receipt.wall_seconds;
    bcost{j}=s.receipt;
end
assert(spent+120<7200,'Core budget'); mkdir(dest);
p=load(fullfile(out,'data/p0_data.mat')); old=load(fullfile(out,'data/p1.mat'));
d=load(fullfile(out,'data/s2b_data.mat')); f=p.fields; grid=d.grid; src=f.source;
off=find(grid.z_m~=3); x=[f.x_m;grid.x_m(off)]; z=[f.z_m;grid.z_m(off)]; n=numel(x);
ix9=zeros(9,1);
for i=1:9, ix9(i)=find(x==grid.x_m(i)&z==grid.z_m(i)); end
U=complex(zeros(n,4,3));
for m=1:3
    for c=1:4, U(:,c,m)=[f.(method(m)){c};d.U(off,c,m)]; end
end
U_legacy=U; variants={}; reference_rows=cell(3,1); checks=cell(3,1); fk_runs=cell(4,2); costs={};
clock=tic;
tt=tic; flat_exact=-fk_single_receiver_incident(src,x,-z); sec=toc(tt);
costs{end+1}=costrow('flat angular spectrum',1,sec);
assert(norm(flat_exact-U(:,1,2))/norm(U(:,1,2))<1e-11);
for c=2:4
    ids={'weak64','medium64','strong64'}; q64=load(fullfile(out,'reference_followup_01',['bie_' ids{c-1} '.mat']),'boundary');
    if c==2
        q96=old.snapshots{1}.spatial_runs{end}; source_old=old.snapshots{1}.quadrature_runs{1};
        density_match=norm(q64.boundary.density-source_old.density)/norm(source_old.density);
    else
        ids96={'medium96','strong96'}; q=load(fullfile(out,'reference_followup_01',['bie_' ids96{c-2} '.mat']),'boundary'); q96=q.boundary;
        source_old=old.snapshots{c-1}.spatial_runs{end}; density_match=norm(q96.density-source_old.density)/norm(source_old.density);
    end
    assert(density_match<1e-10,'Current-version density reproduction failed');
    tt=tic; u64=evaluate_saved_halfplane_bie_field(q64.boundary,x,z); u96=evaluate_saved_halfplane_bie_field(q96,x,z); sec=toc(tt);
    costs{end+1}=costrow('BIE saved-density whole support evaluation',c,sec); %#ok<AGROW>
    assert(norm(u96-U(:,c,2))/norm(U(:,c,2))<1e-11);
    cp=load(fullfile(out,'checkpoints',sprintf('bie%d.mat',c-1)),'payload'); lv=cp.payload.level_fields; ss=old.snapshots{c-1};
    sp=saved_support(ss.spatial_runs{end-1},lv.spatial_runs(:,end-1),f.x_m,off);
    w1=saved_support(ss.window_runs{1},lv.window_runs(:,1),f.x_m,off);
    w2=saved_support(ss.window_runs{2},lv.window_runs(:,2),f.x_m,off);
    w3=saved_support(ss.window_runs{3},lv.window_runs(:,3),f.x_m,off);
    variation=abs(u96-sp)+max(abs([w2-w1,w3-w2]),[],2)+abs(u64-u96);
    variants{end+1}=field_variant(c,'BIE_primary',U(:,c,2),x,z); %#ok<AGROW>
    names={'BIE_q64_current','BIE_q96','BIE_spatial_previous','BIE_window60','BIE_window70','BIE_window80'};
    vals={u64,u96,sp,w1,w2,w3};
    for j=1:6, variants{end+1}=field_variant(c,names{j},vals{j},x,z); end %#ok<AGROW>
    reference_rows{c-1}=table(repmat(string(f.case_names{c}),n,1),(1:n)',x,z,abs(u96-sp),max(abs([w2-w1,w3-w2]),[],2),abs(u64-u96), ...
        variation,abs(u96),'VariableNames',{'case_name','sample_index','x_m','z_m','spatial_abs_change','window_abs_change','self_quad_abs_change','resolution_estimate','reference_amplitude'});
    checks{c-1}=table(string(f.case_names{c}),density_match,q64.boundary.offgrid_boundary_residual,q96.offgrid_boundary_residual, ...
        norm(u64-u96)/norm(u96),'VariableNames',{'case_name','current_saved_density_relative','q64_boundary_residual','q96_boundary_residual','self_quad_support_relative'});
end
% Flat controls use the analytic source mirror, not a roughness normalization.
for c=1:4
    A=[0 .01 .05 .2];
    for j=1:2
        tt=tic; fk_runs{c,j}=closeout_s2b_fullfk_analytic_exterior(src,A(c),x,z,400*j,flat_exact); sec=toc(tt);
        costs{end+1}=costrow(sprintf('Complete FullFK Gauss8 %d panels',400*j),c,sec); %#ok<AGROW>
        variants{end+1}=field_variant(c,sprintf('FullFK_complete_GL%d',400*j),fk_runs{c,j}.receiver_field,x,z); %#ok<AGROW>
    end
    % The independent finite-integral subtraction was computed previously.
    tail=flat_exact+U_legacy(:,c,3)-U_legacy(:,1,3);
    variants{end+1}=field_variant(c,'FullFK_analytic_exterior_trapezoid4097',tail,x,z); %#ok<AGROW>
    variants{end+1}=field_variant(c,'FullFK_finite_L60_4097',U_legacy(:,c,3),x,z); %#ok<AGROW>
    U(:,c,3)=fk_runs{c,2}.receiver_field;
end
seconds=toc(clock); assert(seconds<900 && spent+seconds<7200);
reference=vertcat(reference_rows{:}); density_checks=vertcat(checks{:}); sensitivity_fields=vertcat(variants{:});
line_metrics={}; array_metrics={}; fields_line={}; fields_core={}; point_metrics={}; diagnostics={}; increments={}; numerical=cell(4,1); resolutions={};
pairs=[1 2;3 2;1 3]; weights={ones(411,1)/411,f.weights}; kinds={'spatial_uniform_411','legacy_incident_energy_M99'};
for c=1:4
    f.FullFK{c}=U(1:411,c,3);
    for m=1:3
        fields_line{end+1}=raw_table(f.case_names{c},method(m),compose('L%03d',(1:411)'),x(1:411),z(1:411),U(1:411,c,m)); %#ok<AGROW>
        fields_core{end+1}=raw_table(f.case_names{c},method(m),grid.receiver_id,grid.x_m,grid.z_m,U(ix9,c,m)); %#ok<AGROW>
    end
    if c==1, err=zeros(n,1);
    else, err=reference.resolution_estimate(reference.case_name==string(f.case_names{c})); end
    for j=1:3
        ca=pairs(j,1); rb=pairs(j,2);
        for k=1:2
            metric=closeout_s2b_metrics(U(1:411,c,ca),U(1:411,c,rb),weights{k});
            row=metric_row(metric,c,ca,rb,kinds{k},'M99_z3_411',f.case_names); row.metric_kind="E_line";
            line_metrics{end+1}=row; %#ok<AGROW>
        end
        metric=closeout_s2b_metrics(U(ix9,c,ca),U(ix9,c,rb),ones(9,1)/9);
        row=metric_row(metric,c,ca,rb,'spatial_uniform_9','P11-P33',f.case_names);
        row.metric_kind="E_array"; row.worst_receiver=grid.receiver_id(metric.worst_index);
        array_metrics{end+1}=row; %#ok<AGROW>
        a=U(ix9,c,ca); b=U(ix9,c,rb); ph=angle(a.*conj(b));
        point_metrics{end+1}=table(repmat(string(f.case_names{c}),9,1),repmat(string(method(ca)),9,1),repmat(string(method(rb)),9,1), ...
            grid.receiver_id,abs(a-b),abs(a-b)./abs(b),20*log10(abs(a)./abs(b)),ph,abs(b),abs(b).^2, ...
            'VariableNames',{'case_name','candidate','reference','receiver_id','absolute_error','relative_error','db_signed','phase_signed_rad','reference_abs','denominator_squared'}); %#ok<AGROW>
        for support=1:2
            if support==1, ix=(1:411)'; label='M99_z3_411'; else, ix=ix9; label='P11-P33'; end
            model_relative=norm(U(ix,c,ca)-U(ix,c,rb))/norm(U(ix,c,rb));
            fk_error=abs(fk_runs{c,2}.receiver_field-fk_runs{c,1}.receiver_field);
            den=norm(U(ix,c,rb)); bie_resolution=0; pe_resolution=0; fk_resolution=0;
            if ca==2 || rb==2, bie_resolution=norm(err(ix))/den; end
            if ca==1 || rb==1
                pe=load(fullfile(out,'checkpoints',sprintf('pe%d.mat',c)),'payload');
                % Grid sensitivity is measured on the frozen probes only;
                % its max is conservatively carried over this support.
                pe_bound=max(abs(pe.payload.grid_double_u-pe.payload.u))+max(abs(pe.payload.step_half_line-f.PE{c}));
                pe_resolution=sqrt(numel(ix))*pe_bound/den;
            end
            if ca==3 || rb==3, fk_resolution=norm(fk_error(ix))/den; end
            roundoff=256*eps*norm(flat_exact(ix))/den;
            resolution=bie_resolution+pe_resolution+fk_resolution+roundoff;
            resolutions{end+1}=table(string(f.case_names{c}),string(method(ca)),string(method(rb)),string(label),model_relative, ...
                bie_resolution,pe_resolution,fk_resolution,roundoff,resolution, ...
                model_relative>5*resolution,'VariableNames',{'case_name','candidate','reference','support_id','relative_difference', ...
                'BIE_resolution_estimate','PE_sampled_resolution_estimate','FullFK_resolution_estimate','roundoff_allowance','numerical_resolution_estimate','resolved_above_5x_estimate'}); %#ok<AGROW>
        end
        if c>1
            ef=U(ix9,1,ca)-U(ix9,1,rb); er=a-b; ed=er-ef; den=mean(abs(U(ix9,1,rb)).^2);
            terms=[mean(abs(ef).^2),mean(abs(ed).^2),2*real(mean(ef.*conj(ed))),mean(abs(er).^2)]/den;
            increments{end+1}=table(string(f.case_names{c}),string(method(ca)),string(method(rb)),den,terms(1),terms(2),terms(3),terms(4), ...
                terms(4)-sum(terms(1:3)),'VariableNames',{'case_name','candidate','reference','denominator','flat_term','delta_term','cross_term','rough_term','closure'}); %#ok<AGROW>
        end
    end
    for m=1:3
        if c>1
            u=U(ix9,c,m); flat=U(ix9,1,m); gg=u./flat; du=u-flat;
            diagnostics{end+1}=table(repmat(string(f.case_names{c}),9,1),repmat(string(method(m)),9,1),grid.receiver_id,real(gg),imag(gg),real(du),imag(du), ...
                'VariableNames',{'case_name','method','receiver_id','G_real','G_imag','delta_real','delta_imag'}); %#ok<AGROW>
        end
    end
    tail=flat_exact+U_legacy(:,c,3)-U_legacy(:,1,3);
    numerical{c}=table(string(f.case_names{c}),norm(fk_runs{c,2}.receiver_field-fk_runs{c,1}.receiver_field)/norm(fk_runs{c,2}.receiver_field), ...
        max(abs(fk_runs{c,2}.receiver_field-tail)),max(abs(fk_runs{c,2}.receiver_field-fk_runs{c,2}.term_incident_dn-fk_runs{c,2}.term_green_dn)), ...
        'VariableNames',{'case_name','FullFK_GL400_to800_relative','FullFK_independent_tail_max_abs','FullFK_two_term_closure_max_abs'});
end
line_metrics=vertcat(line_metrics{:}); array_metrics=vertcat(array_metrics{:}); fields_line=vertcat(fields_line{:}); fields_core=vertcat(fields_core{:});
point_metrics=vertcat(point_metrics{:}); diagnostics=vertcat(diagnostics{:}); increments=vertcat(increments{:}); numerical=vertcat(numerical{:}); resolutions=vertcat(resolutions{:}); costs=vertcat(costs{:});
summary=struct('status','S2B_REFERENCE_CLOSURE_PENDING_INDEPENDENT_AUDIT','coverage',108,'P0_preserved',true, ...
    'reference_closure_seconds',seconds,'prior_charged_seconds',spent,'charged_total_seconds',spent+seconds, ...
    'new_BIE_boundary_solves',numel(files),'formal_FullFK_version','complete two-term integral with analytic exact flat exterior, Gauss8/800panels', ...
    'resolution_kind','empirical matched-support spatial+window+self-quadrature variations, not rigorous error bound', ...
    'G1_target_PM','OPEN','WP2_formal',false,'production_solver_changes',false);
save(fullfile(dest,'reference_data.mat'),'U','U_legacy','ix9','x','z','grid','src','flat_exact','reference','density_checks','sensitivity_fields', ...
    'line_metrics','array_metrics','fields_line','fields_core','point_metrics','diagnostics','increments','numerical','resolutions','costs','bcost','fk_runs','summary','-v7.3');
tabs={'reference','density_checks','sensitivity_fields','line_metrics','array_metrics','fields_line','fields_core','point_metrics','diagnostics','increments','numerical','resolutions','costs'};
values={reference,density_checks,sensitivity_fields,line_metrics,array_metrics,fields_line,fields_core,point_metrics,diagnostics,increments,numerical,resolutions,costs};
for j=1:numel(tabs), writetable(values{j},fullfile(dest,[tabs{j} '.csv'])); end
closeout_s2b_json(fullfile(dest,'summary_pending_audit.json'),summary);
pd=load(fullfile(out,'data/plot_data.mat'),'plot_data'); plot_data=pd.plot_data;
plot_data.U=U(ix9,:,:); plot_data.line=f; plot_data.line_metrics=line_metrics; plot_data.array_metrics=array_metrics;
plot_data.status_caption='Same-source complete reflected fields; numerical resolution assessed separately; no physical PASS';
save(fullfile(dest,'plot_data.mat'),'plot_data','-v7.3'); disp(summary); disp(density_checks); disp(numerical);
end

function u=saved_support(b,probe,line_x,off)
[ok,ix]=ismember(line_x,b.receiver_x_m); assert(all(ok)); u=[b.receiver_field(ix);probe(off)];
end
function row=field_variant(c,name,u,x,z)
row=table(repmat(c,numel(u),1),repmat(string(name),numel(u),1),(1:numel(u))',x,z,real(u),imag(u), ...
    'VariableNames',{'case_index','variant','sample_index','x_m','z_m','u_real','u_imag'});
end
function row=costrow(name,c,seconds)
row=table(string(name),c,seconds,'VariableNames',{'job_kind','case_index','wall_seconds'});
end
function row=metric_row(m,c,ca,rb,weight,support,cases)
row=addvars(struct2table(m),string(cases{c}),string(method(ca)),string(method(rb)),string(weight),string(support), ...
    'Before',1,'NewVariableNames',{'case_name','candidate','reference','weight_kind','support_id'});
end
function row=raw_table(c,m,ids,x,z,u)
n=numel(u); row=table(repmat(string(c),n,1),repmat(string(m),n,1),string(ids),x,z,real(u),imag(u), ...
    repmat("native exp(-i omega t); outgoing exp(+ikR)",n,1),repmat("GAUSS_X_4K_UNIT_PEAK_N984",n,1), ...
    repmat("unit-source response; not Pa",n,1),repmat("physical_map; no calibration; no rough/flat normalization",n,1), ...
    'VariableNames',{'case_name','method','receiver_id','x_m','z_m','u_real','u_imag','representation','source_id','unit','mapping'});
row.u_native_real=real(u); row.u_native_imag=imag(u); row.u_mapped_real=real(u); row.u_mapped_imag=imag(u);
row.frequency_hz=4000+zeros(n,1); row.source_tier=repmat("physical_map",n,1);
row.source_map_version=repmat("S2a-depth-positive-native-v1",n,1);
row.reference_id=repmat("BIE_native_same_source_requalified",n,1); row.validity=isfinite(u);
row.quantity_kind=repmat("u_reflected_native",n,1);
if strcmp(m,'FullFK'), version="FullFK_complete_analytic_flat_exterior_GL8_v1";
else, version="S2a_native_retained_and_reference_requalified";
end
row.field_version=repmat(version,n,1);
end
function name=method(m)
names={'PE','BIE','FullFK'}; name=names{m};
end
