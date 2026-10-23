function audit = audit_two_dimensional_closeout_s2c_b_p1_csv(run_dir)
%AUDIT_TWO_DIMENSIONAL_CLOSEOUT_S2C_B_P1_CSV Independent saved CSV arithmetic.
% This does not qualify the one-level K_HIGH BIE numerical reference.
arguments
    run_dir (1,:) char
end
T=readtable(fullfile(run_dir,'P1','representative_BIE_physical_response_UNQUALIFIED.csv'));
S=readtable(fullfile(run_dir,'P1','representative_BIE_supports_UNQUALIFIED.csv'),TextType='string');
assert(height(T)==417 && height(S)==5 && all(T.sample_index(:)==(1:417)'));
u0=complex(T.BASE_M_BIE_real,T.BASE_M_BIE_imag);
u1=complex(T.K_HIGH_BIE_main96_real,T.K_HIGH_BIE_main96_imag);
du=u1-u0;
point_error=max(abs(du-complex(T.delta_K_real,T.delta_K_imag)));
magnitude_error=max(abs([abs(u0)-T.BASE_M_BIE_abs;abs(u1)-T.K_HIGH_BIE_main96_abs; ...
    abs(du)-T.delta_K_abs]));
db_error=max(abs(20*log10(abs(u1)./abs(u0))-T.amplitude_change_dB));
phase_error=max(abs(angle(exp(1i*(angle(u1.*conj(u0))-T.wrapped_phase_change_rad)))));
a2=(abs(u1)-abs(u0)).^2;
p2=4*abs(u1).*abs(u0).*sin(angle(u1.*conj(u0))/2).^2;
identity_error=max(abs(abs(du).^2-a2-p2));
assert(point_error<1e-13 && magnitude_error<1e-13 && db_error<1e-12);
assert(phase_error<1e-12 && identity_error<1e-12);
root=fileparts(fileparts(fileparts(fileparts(run_dir))));
f=load(fullfile(root,'results/comparisons/two_dimensional_channel_closeout', ...
    '20261023_s2b_1305/formal_reference_v2/reference_data.mat'),'ix9');
supports={(1:411)',find(abs(T.x_m(1:411))<=1),f.ix9(:), ...
    f.ix9([2 4 5 6 8]),find(abs(T.x_m(1:411))>=30)};
names=["line_411","center_11","nine_9","cross_5","edge_104"];
metric_error=0;
for j=1:5
    ii=supports{j}; q=ones(numel(ii),1)/numel(ii);
    den=sum(q.*abs(u0(ii)).^2);
    E=sqrt(sum(q.*abs(du(ii)).^2)/den);
    [v,k]=max(abs(du(ii))); k=ii(k);
    row=S(j,:);
    assert(row.support_id==names(j) && row.n==numel(ii));
    assert(~row.K_HIGH_reference_resolution_qualified);
    checks=[den-row.BASE_M_denominator,E-row.BIE_response_L2_provisional, ...
        mean(abs(du(ii)))-row.mean_abs_delta,median(abs(du(ii)))-row.median_abs_delta, ...
        v-row.max_abs_delta,T.x_m(k)-row.worst_x_m,T.z_m(k)-row.worst_z_m];
    metric_error=max(metric_error,max(abs(checks)));
end
assert(metric_error<1e-13);
audit=struct('pass',true,'source','saved CSV only plus frozen support indices', ...
    'field_rows',height(T),'support_rows',height(S),'point_complex_max_abs',point_error, ...
    'magnitude_max_abs',magnitude_error,'amplitude_db_max_abs',db_error, ...
    'wrapped_phase_max_abs',phase_error,'amplitude_phase_square_identity_max_abs',identity_error, ...
    'support_metrics_max_abs',metric_error,'reference_qualified',false, ...
    'no_new_acoustic_calculation',true);
closeout_s2b_json(fullfile(run_dir,'audit','P1_csv_independent_audit.json'),audit);
disp(audit);
end
