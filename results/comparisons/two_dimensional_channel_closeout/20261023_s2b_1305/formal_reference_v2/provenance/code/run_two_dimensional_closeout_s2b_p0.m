function result=run_two_dimensional_closeout_s2b_p0(run_id)
% Pure saved-data P0. No propagation, incident, layer or FK evaluations.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
assert(~isempty(regexp(run_id,'^[A-Za-z0-9_-]+$','once')));
out=fullfile(root,'results/comparisons/two_dimensional_channel_closeout',run_id);
assert(~isfolder(out),'Refusing overwrite'); mkdir(out); mkdir(fullfile(out,'data'));
old=fullfile(root,'results/comparisons/two_dimensional_channel_closeout/20261023_s2a_final');
t=load(fullfile(old,'data/s2a_data.mat'),'fields','metrics','cost','checks');
fields=t.fields; assert(numel(fields.x_m)==411 && all(fields.z_m==3));
T=readtable(fullfile(old,'data/line_fields.csv'),'TextType','string');
assert(height(T)==4932); methods={'PE','BIE','FullFK'};
max_field_error=0; max_weight_error=0;
for c=1:4
    for j=1:3
        s=T(T.case_name==string(fields.case_names{c}) & T.method==string(methods{j}),:);
        s=sortrows(s,'line_index'); assert(height(s)==411);
        assert(isequal(s.x_m,fields.x_m) && all(s.z_m==3));
        max_field_error=max(max_field_error,max(abs(complex(s.u_real,s.u_imag)-fields.(methods{j}){c})));
        max_weight_error=max(max_weight_error,max(abs(s.weight-fields.weights)));
    end
end
assert(max_field_error<1e-14 && max_weight_error<1e-14);
q=ones(411,1)/411; w=fields.weights(:); assert(abs(sum(w)-1)<1e-12);
pairs={'PE','BIE';'FullFK','BIE';'PE','FullFK'}; rows=cell(24,1); n=0; max_csv_error=0; max_legacy_error=0;
for c=1:4
    for p=1:3
        a=fields.(pairs{p,1}){c}; b=fields.(pairs{p,2}){c};
        for k=1:2
            weights={q,w}; kinds={'spatial_uniform_411','legacy_incident_energy_M99'};
            m=closeout_s2b_metrics(a,b,weights{k}); n=n+1;
            row=addvars(struct2table(m),string(fields.case_names{c}),string(pairs{p,1}),string(pairs{p,2}), ...
                string(kinds{k}),"M99_z3_411",'Before',1,'NewVariableNames',{'case_name','candidate','reference','weight_kind','support_id'});
            rows{n}=row;
            % Independent CSV path: direct cartesian and cosine identities.
            sa=sortrows(T(T.case_name==string(fields.case_names{c})&T.method==string(pairs{p,1}),:),'line_index');
            sb=sortrows(T(T.case_name==string(fields.case_names{c})&T.method==string(pairs{p,2}),:),'line_index');
            ac=complex(sa.u_real,sa.u_imag); bc=complex(sb.u_real,sb.u_imag);
            if k==1, wc=ones(height(sa),1)/height(sa); else, wc=sa.weight/sum(sa.weight); end
            denom=sum(wc.*(sb.u_real.^2+sb.u_imag.^2));
            ec=sum(wc.*((sa.u_real-sb.u_real).^2+(sa.u_imag-sb.u_imag).^2))/denom;
            magc=sum(wc.*(abs(ac)-abs(bc)).^2)/denom;
            phc=sum(wc.*2.*(abs(ac).*abs(bc)-real(ac.*conj(bc))))/denom;
            max_csv_error=max([max_csv_error,abs(sqrt(ec)-m.complex_l2),abs(ec-magc-phc)]);
            if k==2
                ix=(c-1)*3+p;
                for name={'complex_l2','magnitude_l2','phase_rms_rad','db_rms','magnitude_squared_term','phase_squared_term'}
                    max_legacy_error=max(max_legacy_error,abs(m.(name{1})-t.metrics.(name{1})(ix)));
                end
            end
        end
    end
end
metrics=vertcat(rows{:}); assert(max_csv_error<1e-12 && max_legacy_error<1e-12);
weights=table((1:411)',fields.x_m,fields.z_m,q,w,'VariableNames',{'line_index','x_m','z_m','spatial_uniform_411','legacy_incident_energy_M99'});
% Formula controls include cancellation, zero, near-zero and missing samples.
aa={[1;2i],[2;4i],exp(.2i)*[1;2i],[1;-1],[0;0],[1;2],[0;0],[1;1],[1;NaN]};
bb={[1;2i],[1;2i],[1;2i],[1;1],[1;2],[0;0],[0;0],[1e-30;1],[1;1]};
tests=cell(numel(aa),1);
for j=1:numel(aa)
    m=closeout_s2b_metrics(aa{j},bb{j},[.5;.5]);
    tests{j}=m; if isfinite(m.decomposition_residual), assert(abs(m.decomposition_residual)<1e-12); end
end
assert(tests{1}.complex_l2==0 && abs(tests{2}.complex_l2-1)<eps);
assert(abs(tests{3}.complex_l2-2*sin(.1))<1e-14 && tests{5}.n_phase_valid==0);
assert(isnan(tests{6}.complex_l2) && isnan(tests{7}.complex_l2));
assert(tests{8}.n_near_null==1 && tests{9}.n_valid==1);
audit=struct('pass',true,'physical_field_computations',0,'csv_rows',height(T), ...
    'metric_rows',height(metrics),'formula_controls',numel(aa),'max_field_csv_mat_error',max_field_error, ...
    'max_weight_csv_mat_error',max_weight_error,'max_independent_csv_error',max_csv_error, ...
    'max_legacy_reproduction_error',max_legacy_error,'near_null_reference_uncertainty','unknown');
save(fullfile(out,'data/p0_data.mat'),'fields','metrics','weights','audit','tests','-v7.3');
writetable(metrics,fullfile(out,'data/line_metrics_both_weights.csv'));
writetable(weights,fullfile(out,'data/line_weights.csv'));
copyfile(fullfile(old,'data/line_fields.csv'),fullfile(out,'data/line_fields_source.csv'));
fid=fopen(fullfile(out,'data/p0_audit.json'),'w'); fprintf(fid,'%s',jsonencode(audit,PrettyPrint=true)); fclose(fid);
result=struct('run_dir',out,'audit',audit); disp(audit); disp(metrics(:,[1:5 19]));
end
