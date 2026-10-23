function audit=audit_two_dimensional_closeout_s2b_tail(out)
% Saved-data-only analytic exterior test; never replaces frozen fields.
dest=fullfile(out,'reference_followup_01'); assert(~isfolder(dest)); mkdir(dest);
p=load(fullfile(out,'data/p0_data.mat'),'fields'); f=p.fields;
d=load(fullfile(out,'data/s2b_data.mat'),'U','grid');
fl=load(fullfile(out,'checkpoints/fk1.mat'),'payload');
mirror=load(fullfile(out,'checkpoints/flat.mat'),'payload');
U_complete=complex(zeros(9,4,4)); rows={}; line_complete=cell(4,1); line_rows=cell(4,1);
for c=1:4
    r=load(fullfile(out,'checkpoints',sprintf('fk%d.mat',c)),'payload');
    for j=1:4
        rough=r.payload.runs{j}; flat=fl.payload.runs{j};
        assert(isequal(rough.surface_x_m,flat.surface_x_m));
        assert(all(rough.receiver_x_m(1:9)==d.grid.x_m) && all(rough.receiver_z_m(1:9)==d.grid.z_m));
        assert(isinf(rough.local_radius_m)&&isinf(flat.local_radius_m));
        ext=abs(rough.surface_x_m)>=50;
        assert(all(rough.surface_z_m(ext)==0) && all(rough.surface_normal(ext,1)==0));
        assert(max(abs(rough.surface_incident(ext)-flat.surface_incident(ext)))<1e-14);
        assert(max(abs(rough.surface_incident_dn(ext)-flat.surface_incident_dn(ext)))<1e-14);
        U_complete(:,c,j)=mirror.payload.u(1:9)+rough.receiver_field(1:9)-flat.receiver_field(1:9);
        for i=1:9
            u=U_complete(i,c,j);
            rows{end+1,1}=table(string(f.case_names{c}),j,d.grid.receiver_id(i),d.grid.x_m(i),d.grid.z_m(i), ...
                real(u),imag(u),real(rough.receiver_field(i)),imag(rough.receiver_field(i)), ...
                real(mirror.payload.u(i)-flat.receiver_field(i)),imag(mirror.payload.u(i)-flat.receiver_field(i)), ...
                'VariableNames',{'case_name','variant','receiver_id','x_m','z_m','u_complete_real','u_complete_imag', ...
                'u_finite_real','u_finite_imag','analytic_exterior_real','analytic_exterior_imag'}); %#ok<AGROW>
        end
    end
    line_complete{c}=f.BIE{1}+f.FullFK{c}-f.FullFK{1};
    m=closeout_s2b_metrics(line_complete{c},f.BIE{c},ones(411,1)/411);
    line_rows{c}=addvars(struct2table(m),string(f.case_names{c}),"FullFK_analytic_exterior_candidate", ...
        "spatial_uniform_411","M99_z3_411",'Before',1,'NewVariableNames',{'case_name','method','weight_kind','support_id'});
end
grid_delta=max(abs(U_complete(:,:,3)-U_complete(:,:,2)),[],'all');
window_delta=max(abs(U_complete(:,:,4)-U_complete(:,:,3)),[],'all');
assert(window_delta<1e-13);
audit=struct('saved_data_only',true,'new_field_evaluations',0,'same_source',true, ...
    'rough_flat_exterior_identical_for_abs_x_ge_50',true,'max_complete_fine_window_abs_difference',window_delta, ...
    'max_complete_4097_to8193_abs_difference',grid_delta, ...
    'formula','u_complete = analytic_flat + FK_rough_finite - FK_flat_finite', ...
    'status','ANALYTIC_EXTERIOR_IDENTITY_VERIFIED_ON_SAVED_ARRAYS_NOT_YET_PROMOTED', ...
    'interpretation','Exact flat exterior added analytically, not fitted; original finite-window fields retained. Does not resolve BIE qualification by itself.');
fields=vertcat(rows{:}); line_metrics=vertcat(line_rows{:});
save(fullfile(dest,'tail_identity.mat'),'U_complete','line_complete','line_metrics','audit','-v7.3');
writetable(fields,fullfile(dest,'nine_point_variants.csv'));
writetable(line_metrics,fullfile(dest,'line_candidate_metrics.csv'));
closeout_s2b_json(fullfile(dest,'audit.json'),audit); disp(audit);
end
