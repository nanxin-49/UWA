function checks=audit_deepsea3d_pm_feasibility_artifacts()
%AUDIT_DEEPSEA3D_PM_FEASIBILITY_ARTIFACTS Recompute saved-data contracts.
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
out=fullfile(root,'results','validation','deepsea3d_pm_feasibility');
a=load(fullfile(out,'feasibility.mat'),'state');s=a.state;
assert(s.complete && numel(s.runs)==9 && numel(s.trace)==2);
assert(~s.reference_ready && ~s.support_controlled && s.FK_solves==0);
assert(s.global_G1=="OPEN" && s.G1_status=="INCONCLUSIVE" && ~s.WP2_allowed);
sea=build_deepsea3d_pm_target();assert(isequaln(sea,s.sea));
assert(isequal(s.cfg.tx_m,[0 150.5 1500]) && isequal(s.cfg.rx_m,[0 0 5]));
assert(abs(norm(s.cfg.tx_m(1:2)-s.cfg.rx_m(1:2))/(1500+5)-.1)<1e-15);
T=readtable(fullfile(out,'parent_coefficients.csv'));
assert(height(T)==numel(sea.k_radpm) && nnz(T.baseline)==652);
assert(norm(T.c_real_m-sea.coefficient_real)+norm(T.c_imag_m-sea.coefficient_imag)<1e-13);
assert(norm(T.kx_radpm-sea.kx_radpm)+norm(T.ky_radpm-sea.ky_radpm)<1e-12);
R=readtable(fullfile(out,'runs.csv'),TextType='string');
assert(height(R)==9);
for j=1:9
    r=s.runs(j).result;
    assert(R.name(j)==s.runs(j).name);
    assert(abs(complex(R.delta_real(j),R.delta_imag(j))-r.rough_correction)<1e-17);
    assert(abs(complex(R.reflect_real(j),R.reflect_imag(j))-r.reflect_native)<1e-17);
    assert(abs(complex(R.total_real(j),R.total_imag(j))-r.total_native)<1e-17);
    assert(abs(r.total_native-r.reflect_native-s.flat_control.direct_native)<1e-17);
    assert(all(isfinite([r.terms(:);r.boundary_coefficient_norms(:)])));
end
assert(s.guards.N16_to_N20>1e-4 && s.guards.boundary>1e-4 ...
    && s.runs(3).result.forcing_series_relative>1e-8);
for j=1:2
    d=s.trace(j);den=norm(d.flat_native);
    assert(abs(norm(d.flat_native+d.exact_trace_native)/den-d.exact_trace_residual_relative)<1e-12);
    assert(d.exact_trace_residual_relative>1 && d.term_match_error<1e-12 && d.norm_match_error<1e-12);
end
% Independently recompute mirror solution and native-to-DSP representation.
k=2*pi*6000/1500;Ld=norm(s.cfg.tx_m-s.cfg.rx_m);
image=s.cfg.tx_m;image(3)=-image(3);Lr=norm(image-s.cfg.rx_m);
pd=exp(1i*k*Ld)/(4*pi*Ld);pr=-exp(1i*k*Lr)/(4*pi*Lr);
assert(abs(pd-s.flat_control.direct_native)+abs(pr-s.flat_control.reflect_native)<1e-16);
assert(abs(conj(pd+pr)*exp(1i*k*Ld)-s.flat_control.total_dsp)<1e-16);
checks=struct('surface_rebuild_exact',true,'coefficients_CSV',true,'native_pressure_CSV',true, ...
    'trace_replay',true,'flat_phase',true,'numerical_failure_classification',true, ...
    'G1_status','INCONCLUSIVE','global_G1','OPEN','WP2_allowed',false,'FK_solves',0);
fid=fopen(fullfile(out,'independent_artifact_audit.json'),'w');assert(fid>=0);
cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(checks,PrettyPrint=true));
disp(checks);
end
