function evidence=assert_spectral_fast_kirchhoff_pm_pilot_artifacts()
%ASSERT_SPECTRAL_FAST_KIRCHHOFF_PM_PILOT_ARTIFACTS Verify saved pilot evidence.
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
out=fullfile(root,'results','validation','spectral_fast_kirchhoff_pm_pilot');
q=load(fullfile(out,'spectral_fast_fk_pm_pilot.mat'),'result'); r=q.result;
assert(numel(r.profiles)==5 && numel(r.cases)==5 && numel(r.bie)==2);
assert(height(r.metrics)==10 && numel(r.rank_diagnostics)==45);
assert(numel(r.untapered_cases)==5 && height(r.untapered_metrics)==5);
base=readtable(r.profiles(1).coeff_file);
high=readtable(r.profiles(5).coeff_file); low=readtable(r.profiles(4).coeff_file);
assert(isequal(base{:,:},high{1:height(base),:}));
assert(isequal(base{1:height(low),:},low{:,:}));
assert(r.padding_2_4_l2<1e-12 && r.grid.surface_sample_closure<1e-12);
assert(r.grid.full_weighted_l2<1e-4 && r.grid.full_whole_l2<1e-4);
assert(r.height_nodes_33_65_l2<1e-4);
assert(isfield(r,'highK_full_quadrature'));
q=r.highK_full_quadrature;
assert(q.n_surface==8193 && q.convergence_weighted_l2<1e-4);
actual=norm(q.fine_full_field-q.bie_field)/norm(q.bie_field);
assert(abs(actual-q.fine_full_bie_whole_l2)<1e-12);
for i=1:5
 d=r.cases(i); raw=r.untapered_cases(i); idx=find(strcmp(r.metrics.case_name,d.name) & strcmp(r.metrics.study,'selected'));
 expected=norm(d.selected.receiver_field-d.reference_field)/norm(d.reference_field);
 assert(abs(expected-r.metrics.whole_grid_complex_l2(idx))<1e-12);
 assert(d.selected.rank==r.selected_rank_in && d.selected.outgoing_rank==r.selected_rank_out);
 assert(d.selected.dense_receiver_surface_entries==0 && d.selected.linear_convolution && ~d.selected.sponge);
 [eta,slope]=evaluate_pm_pilot_surface(r.profiles(i),d.x_m);
 assert(max(abs(d.eta_m-eta))<1e-12 && max(abs(d.slope-slope))<1e-12);
 [eta,slope]=evaluate_pm_pilot_surface(r.profiles(i),raw.x_m,false);
 assert(max(abs(raw.eta_m-eta))<1e-12 && max(abs(raw.slope-slope))<1e-12);
 assert(max(abs(raw.eta_m))<r.common_height_bound_m);
 assert(max(r.profiles(i).k_rad_per_m)<pi/mean(diff(d.x_m))/50);
 s=d.spectrum;
 assert(abs(norm(s.FastFK-s.FK)/norm(s.FK)-s.complex_l2)<1e-12);
 assert(abs(s.propagating_energy_fraction+s.near_evanescent_energy_fraction+s.far_evanescent_energy_fraction-1)<1e-12);
 assert(exist(fullfile(out,[d.name '_reflected_spectrum.csv']),'file')==2);
end
for i=2:3
 p=r.profiles(i); tbl=readtable(p.coeff_file); stream=RandStream('mt19937ar','Seed',p.seed);
 sd=sqrt(tbl.Sk_m3*(2*pi/160));
 a=sd.*randn(stream,height(tbl),1); b=sd.*randn(stream,height(tbl),1);
 assert(max(abs(a-tbl.cos_coeff_m))<1e-14 && max(abs(b-tbl.sin_coeff_m))<1e-14);
end
for i=1:2
 b=r.bie(i); assert(numel(b.runs)==2+(i==2));
 assert(all(isfinite([b.runs.weighted_complex_l2])) && all([b.runs.linear_residual]<1e-10));
end
evidence=struct('artifact_completeness_pass',true,'algorithm_pass',r.algorithm_pass, ...
 'untapered_pass',r.untapered_pass,'physics_pass',r.physics_pass,'pilot_pass',r.pilot_pass);
fprintf('PM_PILOT_ARTIFACT_ASSERTIONS_PASS\n'); disp(evidence);
end
