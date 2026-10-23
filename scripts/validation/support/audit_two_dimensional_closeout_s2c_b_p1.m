function audit = audit_two_dimensional_closeout_s2c_b_p1(run_dir)
%AUDIT_TWO_DIMENSIONAL_CLOSEOUT_S2C_B_P1 Audit one BIE receipt and stop gate.
% Saved-data work only: no new acoustic evaluation or integration.
arguments
    run_dir (1,:) char
end
root=fileparts(fileparts(fileparts(fileparts(run_dir))));
formal=fullfile(root,'results/comparisons/two_dimensional_channel_closeout', ...
    '20261023_s2b_1305/formal_reference_v2/reference_data.mat');
saved=load(fullfile(run_dir,'jobs','bie_main96.mat'),'boundary','receipt');
reported=jsondecode(fileread(fullfile(run_dir,'jobs','bie_main96.json')));
preflight=jsondecode(fileread(fullfile(run_dir,'jobs','bie_main96_preflight.json')));
gate=jsondecode(fileread(fullfile(run_dir,'protocol','P1_initial_gate.json')));
assert(strcmp(saved.receipt.status,'COMPLETE') && strcmp(reported.status,'COMPLETE'));
assert(strcmp(saved.receipt.job_id,'main96') && strcmp(reported.job_id,'main96'));
assert(saved.receipt.wall_seconds==reported.wall_seconds);
assert(saved.receipt.node_count==3840 && saved.boundary.node_count==3840);
assert(saved.boundary.self_quadrature_order==96 && saved.boundary.half_width_m==60);
assert(saved.boundary.points_per_wavelength_requested==12);
assert(saved.boundary.finite && saved.boundary.linear_residual<1e-10);
assert(saved.boundary.offgrid_boundary_residual<1e-6);
assert(numel(saved.boundary.receiver_field)==417);
assert(all(isfinite(real(saved.boundary.receiver_field))) && ...
    all(isfinite(imag(saved.boundary.receiver_field))));
assert(saved.receipt.wall_seconds<gate.single_job_cap_s);
assert(preflight.available_physical_bytes==saved.receipt.available_physical_bytes);
assert(preflight.estimated_increment_bytes<=.7*preflight.available_physical_bytes);
assert(preflight.maximum_remaining_increment_bytes<=.7*preflight.available_physical_bytes);
assert(strcmp(saved.receipt.solver_sha256, ...
    closeout_s2b_hash(fullfile(root,'scripts/validation/support/solve_helmholtz_bie_halfplane_vertical.m'))));
assert(strcmp(saved.receipt.configuration_sha256, ...
    jsondecode(fileread(fullfile(run_dir,'protocol','authorization.json'))).frozen_configuration_sha256));
f=load(formal,'U','x','z','ix9');
assert(max(abs(saved.boundary.receiver_x_m(:)-f.x(:)))<1e-12);
assert(max(abs(saved.boundary.receiver_z_m(:)-f.z(:)))<1e-12);
u0=f.U(:,3,2); u1=saved.boundary.receiver_field(:); du=u1-u0;
assert(isequal(size(u1),[417 1]));
point=table((1:417)',f.x(:),f.z(:),real(u0),imag(u0),real(u1),imag(u1), ...
    real(du),imag(du),abs(u0),abs(u1),abs(du), ...
    20*log10(abs(u1)./abs(u0)),angle(u1.*conj(u0)), ...
    'VariableNames',{'sample_index','x_m','z_m','BASE_M_BIE_real','BASE_M_BIE_imag', ...
    'K_HIGH_BIE_main96_real','K_HIGH_BIE_main96_imag','delta_K_real','delta_K_imag', ...
    'BASE_M_BIE_abs','K_HIGH_BIE_main96_abs','delta_K_abs','amplitude_change_dB', ...
    'wrapped_phase_change_rad'});
dest=fullfile(run_dir,'P1'); if ~isfolder(dest), mkdir(dest); end
writetable(point,fullfile(dest,'representative_BIE_physical_response_UNQUALIFIED.csv'));
center=find(abs(f.x(1:411))<=1); edge=find(abs(f.x(1:411))>=30);
five=f.ix9([2 4 5 6 8]); supports={(1:411)',center,f.ix9(:),five(:),edge};
names=["line_411","center_11","nine_9","cross_5","edge_104"];
rows=table;
for j=1:5
    ii=supports{j}; n=numel(ii); q=ones(n,1)/n;
    den=sum(q.*abs(u0(ii)).^2);
    E=sqrt(sum(q.*abs(du(ii)).^2)/den);
    [worst,k]=max(abs(du(ii))); k=ii(k);
    rows=[rows;table(names(j),n,den,E,mean(abs(du(ii))),median(abs(du(ii))), ...
        worst,f.x(k),f.z(k),false, ...
        'VariableNames',{'support_id','n','BASE_M_denominator','BIE_response_L2_provisional', ...
        'mean_abs_delta','median_abs_delta','max_abs_delta','worst_x_m','worst_z_m', ...
        'K_HIGH_reference_resolution_qualified'})]; %#ok<AGROW>
end
assert(isequal(rows.n,[411;11;9;5;104]));
writetable(rows,fullfile(dest,'representative_BIE_supports_UNQUALIFIED.csv'));
actual=saved.receipt.wall_seconds;
nodes=[3840 3840 4480 3200 3736];
historical=gate.historical_stage_basis_s(:)';
upper=gate.upper_stage_s(:)'; upper(1)=actual;
upper(2)=max(upper(2),1.2*actual*historical(2)/historical(1));
for j=3:5
    upper(j)=max(upper(j),1.2*actual*(nodes(j)/nodes(1))^2);
end
bie_total=sum(upper);
projected=bie_total+gate.other_physics_reserve_s+gate.failure_contingency_s;
assert(projected>gate.subbudget_s);
assert(gate.core_remaining_s-projected>=gate.S3A_reserve_s);
remaining=gate.subbudget_s-actual;
assert(sum(upper(2:5))+gate.other_physics_reserve_s+gate.failure_contingency_s>remaining);
reassessment=struct('decision','P1_STOP_COST_GATE_AFTER_REPRESENTATIVE', ...
    'actual_new_acoustic_seconds',actual,'remaining_S2c_B_subbudget_s',remaining, ...
    'BIE_five_reforecast_s',bie_total,'other_physics_reserve_s',gate.other_physics_reserve_s, ...
    'failure_contingency_s',gate.failure_contingency_s,'total_reforecast_s',projected, ...
    'subbudget_s',gate.subbudget_s,'over_subbudget_s',projected-gate.subbudget_s, ...
    'S3A_balance_after_reforecast_s',gate.core_remaining_s-projected, ...
    'S3A_reserve_s',gate.S3A_reserve_s,'remaining_reference_jobs_unrun',4, ...
    'BIE_level_upper_s',upper,'reference_qualified',false, ...
    'new_PE_jobs',0,'new_FullFK_jobs',0,'source','Independent read of saved MAT/JSON and frozen formal BASE_M');
closeout_s2b_json(fullfile(run_dir,'protocol','P1_after_representative.json'),reassessment);
audit=struct('pass',true,'status','ONE_BIE_SOLVE_SAVED_REFERENCE_INCOMPLETE', ...
    'MAT_JSON_receipt_equal',true,'field_417_finite',true,'same_receivers_as_BASE_M',true, ...
    'same_frozen_configuration_hash',true,'solver_hash_pass',true, ...
    'actual_seconds',actual,'max_abs_delta_K_provisional',max(abs(du)), ...
    'numerical_reference_qualified',false,'stop_cost_gate',true, ...
    'P1_reassessment_sha256',closeout_s2b_hash(fullfile(run_dir,'protocol','P1_after_representative.json')));
closeout_s2b_json(fullfile(run_dir,'audit','P1_independent_audit.json'),audit);
disp(reassessment); disp(audit);
end
