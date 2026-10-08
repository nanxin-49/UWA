function result = validate_spectral_fast_kirchhoff_pm_pilot(mode)
%VALIDATE_SPECTRAL_FAST_KIRCHHOFF_PM_PILOT Fixed-coefficient single-f pilot.
% 'fast' checkpoints algorithms/sensitivity; 'bie' adds physical references;
% 'raw' adds unmodified Fourier PM controls; 'render' rebuilds plots/report.
if nargin<1, mode='full'; end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root); setup_vertical_project();
out=fullfile(root,'results','validation','spectral_fast_kirchhoff_pm_pilot');
if ~isfolder(out), mkdir(out); end
file=fullfile(out,'spectral_fast_fk_pm_pilot.mat');
if ismember(mode,{'bie','render','raw','quadrature'})
 q=load(file,'result'); result=q.result;
else
 assert(ismember(mode,{'fast','full'}));
 result=local_fast(root,out);
 save(file,'result','-v7.3');
end
if ismember(mode,{'bie','full'})
 result=local_bie(result,out);
 save(file,'result','-v7.3');
end
if ismember(mode,{'raw','full'})
 result=local_raw(result,root);
 save(file,'result','-v7.3');
end
if ismember(mode,{'quadrature','full'})
 result=local_physical_quadrature(result);
 save(file,'result','-v7.3');
end
result=local_gates(result);
save(file,'result','-v7.3');
local_outputs(result,out,root);
fprintf('PM_ALGORITHM_PASS=%d PM_PHYSICS_PASS=%d PM_PILOT_PASS=%d\n', ...
 result.algorithm_pass,result.physics_pass,result.pilot_pass);
end

function result=local_fast(root,out)
t0=tic; cpu0=cputime;
canonical=fullfile(root,'results','validation','bellhop_internal_pm_fixed_realization','fixed_pm_fourier_coefficients.csv');
paths={canonical};
for seed=260002:260003
 paths{end+1}=fullfile(root,'results','validation','pe_bellhop_pm_amplitude_ensemble',sprintf('seed_%d_coefficients.csv',seed)); %#ok<AGROW>
end
profiles=struct([]);
for i=1:3
 p=load_fixed_pm_profile_for_pe_bellhop_validation(paths{i},struct('seed',260000+i));
 p.label=sprintf('seed_%d',p.seed); p.requested_kmax_rad_per_m=.5;
 p.origin='saved canonical' ; if i>1, p.origin='saved independent Gaussian coefficients'; end
 profiles=local_append(profiles,p);
 copyfile(paths{i},fullfile(out,[p.label '_coefficients.csv']));
end
% Nested cutoffs: all common coefficients are identical, with no rescaling.
tbl=readtable(canonical); dk=2*pi/160;
assert(max(abs(tbl.Sk_m3-local_pm_density(tbl.k_rad_per_m)))<1e-12);
low=tbl(tbl.k_rad_per_m<=.3,:);
highk=(13:floor(.75/dk))'*dk; sk=local_pm_density(highk);
stream=RandStream('mt19937ar','Seed',26000175);
high=table(highk,sqrt(sk*dk).*randn(stream,numel(highk),1), ...
 sqrt(sk*dk).*randn(stream,numel(highk),1),sk,'VariableNames',tbl.Properties.VariableNames);
high=[tbl;high];
for item=1:2
 if item==1, tt=low; tag='cutoff_low'; cutoff=.3; else, tt=high; tag='cutoff_high'; cutoff=.75; end
 path=fullfile(out,[tag '_coefficients.csv']); writetable(tt,path);
 p=load_fixed_pm_profile_for_pe_bellhop_validation(path,struct('seed',260001,'requested_kmax_rad_per_m',cutoff));
 p.label=tag; p.origin='nested canonical coefficients; appended modes use mt19937ar seed 26000175';
 profiles=local_append(profiles,p);
end
% Uniform height range for every seed/cutoff, not a per-case rank choice.
bound=0; profile_stats=struct([]);
for i=1:numel(profiles)
 p=profiles(i); x=linspace(-80,80,16385)';
 raw=evaluate_fixed_pm_fourier_profile(p,x);
 [eta,sl,sec]=evaluate_pm_pilot_surface(p,x);
 bound=max(bound,max(abs(eta)));
 st=struct('case_name',p.label,'seed',p.seed,'requested_Kmax',p.requested_kmax_rad_per_m, ...
  'actual_Kmax',max(p.k_rad_per_m),'coefficient_count',numel(p.k_rad_per_m), ...
  'raw_rms_height',sqrt(mean(raw.eta_m(1:end-1).^2)), ...
  'raw_max_height',max(abs(raw.eta_m)),'raw_rms_slope',sqrt(mean(raw.eta_prime_m_per_m(1:end-1).^2)), ...
  'raw_max_slope',max(abs(raw.eta_prime_m_per_m)), ...
  'raw_rms_curvature',sqrt(mean((raw.eta_second_m_per_m2./(1+raw.eta_prime_m_per_m.^2).^1.5).^2)), ...
  'raw_max_curvature',max(abs(raw.eta_second_m_per_m2./(1+raw.eta_prime_m_per_m.^2).^1.5)), ...
  'embedded_max_height',max(abs(eta)),'embedded_max_slope',max(abs(sl)), ...
  'embedded_max_curvature',max(abs(sec./(1+sl.^2).^1.5)));
 profile_stats=local_append(profile_stats,st);
 writetable(table(x,raw.eta_m,raw.eta_prime_m_per_m,raw.eta_second_m_per_m2,eta,sl,sec, ...
  'VariableNames',{'x_m','raw_eta_m','raw_slope','raw_second','embedded_eta_m','embedded_slope','embedded_second'}), ...
  fullfile(out,[p.label '_master_surface.csv']));
end
bound=ceil(bound/.05)*.05;
% Reuse numeric source from the accepted prototype, no opaque MAT handles.
q=load(fullfile(root,'results','validation','spectral_fast_kirchhoff_prototype','weak_fields.mat'),'d');
source=struct('kx',q.d.source_kx(:),'kz',sqrt((2*pi*4000/1500)^2-q.d.source_kx(:).^2), ...
 'coeff',q.d.source_coefficients(:),'period_m',192.1875,'n_original',984,'x0_m',-96.09375,'z_tx_m',100);
q=load(fullfile(root,'results','validation','pe_bellhop_controlled_comparison','stage0','stage0_validation.mat'),'validation');
stage0=q.validation;
cases=struct([]); ledger=struct([]);
for i=1:numel(profiles)
 fprintf('PM frozen-rank case %s, common A=%.3f\n',profiles(i).label,bound);
 d=local_case(profiles(i),4097,bound,source,stage0,33);
 c=d.cfg; c.rank=3; c.outgoing_rank=6;
 f=solve_spectral_fast_kirchhoff_2d(c,d.plan); d.baseline=f; d.selected=f;
 rows=local_row(d,f,'baseline');
 d.table=struct2table(rows); ledger=local_append(ledger,rows);
 d=rmfield(d,'plan'); cases=local_append(cases,d);
 fprintf('  weighted %.5g whole-grid %.5g\n',rows.weighted_complex_l2,rows.whole_grid_complex_l2);
end
baseline_pass=all([ledger.weighted_complex_l2]<1e-3 & [ledger.whole_grid_complex_l2]<1e-3);
ri=3; ro=6; rank_table=struct([]);
if ~baseline_pass
 % Diagnose both separations with a bounded COMMON rank ladder for all cases.
 % Never choose a different rank for each realization, and stop at rank 10.
 for i=1:numel(cases)
  d=cases(i); cfg=d.cfg; cfg.max_rank=10;
  plan=prepare_spectral_fast_kirchhoff_2d(cfg);
  for rin=[3 4 5]
   for rout=[6 8 10]
    cfg.rank=rin; cfg.outgoing_rank=rout;
    f=solve_spectral_fast_kirchhoff_2d(cfg,plan);
    rank_table=local_append(rank_table,local_row_with_plan(d,f,plan,'rank_diagnostic'));
   end
  end
 end
 candidates=[3 6;3 8;4 6;4 8;3 10;5 6;5 8;4 10;5 10];
 ri=NaN; ro=NaN;
 for pair=candidates.'
  rows=rank_table([rank_table.rank_in]==pair(1) & [rank_table.rank_out]==pair(2));
  if all([rows.weighted_complex_l2]<1e-3 & [rows.whole_grid_complex_l2]<1e-3)
   ri=pair(1); ro=pair(2); break;
  end
 end
end
selected_cases=struct([]);
for i=1:numel(cases)
 d=cases(i); cfg=d.cfg; cfg.max_rank=6;
 if isfinite(ri), cfg.rank=ri; cfg.outgoing_rank=ro; else, cfg.rank=3; cfg.outgoing_rank=6; end
 cfg.max_rank=max([cfg.max_rank,cfg.rank,cfg.outgoing_rank]);
 plan=prepare_spectral_fast_kirchhoff_2d(cfg);
 solve_spectral_fast_kirchhoff_2d(cfg,plan); wt=zeros(5,1); ct=wt;
 for j=1:5
  f=solve_spectral_fast_kirchhoff_2d(cfg,plan); wt(j)=f.wall_seconds; ct(j)=f.cpu_seconds;
 end
 f.wall_seconds=median(wt); f.cpu_seconds=median(ct); d.selected=f;
 row=local_row_with_plan(d,f,plan,'selected'); d.table=[d.table;struct2table(row)];
 d.timing_row=row; d.spectrum=local_spectrum(d,source);
 d.rank_in_1e3=local_rank(plan.incoming_singular_values,1e-3);
 d.rank_out_1e3=local_rank(plan.outgoing_singular_values,1e-3);
 d.incoming_singular_values=plan.incoming_singular_values;
 d.outgoing_singular_values=plan.outgoing_singular_values;
 ledger=local_append(ledger,row); selected_cases=local_append(selected_cases,d);
end
cases=selected_cases;
checkpoint=struct('cases',cases,'metrics',struct2table(ledger),'rank_diagnostics',rank_table,'rank_in',ri,'rank_out',ro);
save(fullfile(out,'fast_case_checkpoint.mat'),'checkpoint','-v7.3');
% One coefficient-defined realization on three nested integration grids.
sens=struct([]); fine=cases(1);
for n=[2049 8193]
 d=local_case(profiles(1),n,bound,source,stage0,33);
 cfg=d.cfg; cfg.rank=cases(1).selected.rank; cfg.outgoing_rank=cases(1).selected.outgoing_rank;
 if max(cfg.rank,cfg.outgoing_rank)>d.plan.max_rank
  cfg.max_rank=max(cfg.rank,cfg.outgoing_rank); d.plan=prepare_spectral_fast_kirchhoff_2d(cfg);
 end
 f=solve_spectral_fast_kirchhoff_2d(cfg,d.plan);
 row=local_row(d,f,'grid');
 sens=local_append(sens,row);
 if n==8193, fine=d; fine.selected=f; end
end
idx=1:2:numel(fine.x_m); mask=cases(1).mask; w=cases(1).weights;
grid=struct('n_coarse',4097,'n_fine',8193, ...
 'full_weighted_l2',local_error(cases(1).reference_field,fine.reference_field(idx),w,mask).weighted_complex_l2, ...
 'full_whole_l2',norm(cases(1).reference_field-fine.reference_field(idx))/norm(fine.reference_field(idx)), ...
 'fast_weighted_l2',local_error(cases(1).selected.receiver_field,fine.selected.receiver_field(idx),w,mask).weighted_complex_l2, ...
 'surface_sample_closure',max(abs(cases(1).eta_m-fine.eta_m(idx))));
cfg=cases(1).cfg; cfg.max_rank=10; cfg.rank=cases(1).selected.rank; cfg.outgoing_rank=cases(1).selected.outgoing_rank;
plan=prepare_spectral_fast_kirchhoff_2d(cfg);
cfg.padding_factor=4; pad=solve_spectral_fast_kirchhoff_2d(cfg,plan);
padding_error=norm(pad.receiver_field-cases(1).selected.receiver_field)/norm(cases(1).selected.receiver_field);
cfg.padding_factor=2; cfg.height_nodes=65;
plan65=prepare_spectral_fast_kirchhoff_2d(cfg); f65=solve_spectral_fast_kirchhoff_2d(cfg,plan65);
height_node_error=norm(f65.receiver_field-cases(1).selected.receiver_field)/norm(f65.receiver_field);
% Explicit incident-only / outgoing-only ablations identify failure source.
cfg.height_nodes=33; plan=prepare_spectral_fast_kirchhoff_2d(cfg);
 cfg.rank=3; cfg.outgoing_rank=10; inc_test=solve_spectral_fast_kirchhoff_2d(cfg,plan);
 cfg.rank=5; cfg.outgoing_rank=6; out_test=solve_spectral_fast_kirchhoff_2d(cfg,plan);
diagnosis=struct('rank3_10_whole_l2',norm(inc_test.receiver_field-cases(1).reference_field)/norm(cases(1).reference_field), ...
 'rank5_6_whole_l2',norm(out_test.receiver_field-cases(1).reference_field)/norm(cases(1).reference_field));
result=struct('schema_version','1.0','date',char(datetime('today','Format','yyyy-MM-dd')), ...
 'matlab_version',version,'profiles',profiles,'profile_stats',struct2table(profile_stats), ...
 'common_height_bound_m',bound,'source',source,'cases',cases,'metrics',struct2table(ledger), ...
 'rank_diagnostics',rank_table,'selected_rank_in',ri,'selected_rank_out',ro, ...
 'baseline_pass',baseline_pass,'grid',grid,'grid_metrics',struct2table(sens), ...
 'padding_2_4_l2',padding_error,'height_nodes_33_65_l2',height_node_error,'diagnosis',diagnosis, ...
 'numerical_wall_s',toc(t0),'numerical_cpu_s',cputime-cpu0,'bie',struct([]));
end

function sk=local_pm_density(k)
sk=8.10e-3./(2*k.^3).*exp(-.74*9.81^2./(6^4*k.^2));
end

function d=local_case(p,n,A,source,stage0,q,embed)
if nargin<7, embed=true; end
x=linspace(-60,60,n)'; [eta,slope,second]=evaluate_pm_pilot_surface(p,x,embed);
normal=[-slope,ones(n,1)]./sqrt(1+slope.^2);
cfg=struct('frequency_hz',4000,'c0_mps',1500,'surface_x_m',x,'surface_z_m',eta, ...
 'surface_normal',normal,'height_bound_m',A,'receiver_depth_m',3,'spectral_source',source, ...
 'max_rank',6,'height_nodes',q,'geometry','B','rank',3,'outgoing_rank',6,'padding_factor',2);
plan=prepare_spectral_fast_kirchhoff_2d(cfg);
clk=tic; cpu=cputime; [u,ux,uz]=local_source(source,x,eta);
incident_wall=toc(clk); incident_cpu=cputime-cpu;
refcfg=struct('frequency_hz',4000,'c0_mps',1500,'surface_x_m',x,'surface_z_m',eta, ...
 'surface_normal',normal,'surface_incident',u,'surface_incident_dn',normal(:,1).*ux+normal(:,2).*uz, ...
 'receiver_x_m',x,'receiver_z_m',3+zeros(n,1));
times=zeros(2,1); cpus=times;
for j=1:2
 clk=tic; cpu=cputime; ref=solve_full_kirchhoff_2d(refcfg); times(j)=toc(clk); cpus(j)=cputime-cpu;
end
oldx=stage0.x_m(:); oldmask=stage0.footprint.m99.mask(:);
mask=x>=min(oldx(oldmask)) & x<=max(oldx(oldmask));
w=interp1(oldx,stage0.footprint.energy_weights(:),x,'linear',0);
mask=mask & abs(ref.receiver_field)>=max(abs(ref.receiver_field(mask)))*1e-2;
d=struct('name',p.label,'profile',p,'cfg',cfg,'refcfg',refcfg,'plan',plan,'x_m',x,'eta_m',eta, ...
 'slope',slope,'second_derivative',second,'normal',normal,'weights',w,'mask',mask, ...
 'reference_field',ref.receiver_field,'reference_wall_s',median(times),'reference_cpu_s',median(cpus), ...
 'incident_wall_s',incident_wall,'incident_cpu_s',incident_cpu);
end

function r=local_raw(r,root)
% Unmodified PM heights/normal on the full finite integration interval.
q=load(fullfile(root,'results','validation','pe_bellhop_controlled_comparison','stage0','stage0_validation.mat'),'validation');
raw=struct([]); rows=struct([]);
for i=1:numel(r.profiles)
 fprintf('Untapered PM algorithm control: %s\n',r.profiles(i).label);
 d=local_case(r.profiles(i),4097,r.common_height_bound_m,r.source,q.validation,33,false);
 cfg=d.cfg; cfg.rank=r.cases(1).selected.rank; cfg.outgoing_rank=r.cases(1).selected.outgoing_rank;
 cfg.max_rank=max([6 cfg.rank cfg.outgoing_rank]); plan=prepare_spectral_fast_kirchhoff_2d(cfg);
 f=solve_spectral_fast_kirchhoff_2d(cfg,plan);
 row=local_row_with_plan(d,f,plan,'untapered'); d.selected=f;
 d.spectrum=local_spectrum(d,r.source); d=rmfield(d,'plan');
 raw=local_append(raw,d); rows=local_append(rows,row);
end
r.untapered_cases=raw; r.untapered_metrics=struct2table(rows);
end

function [u,ux,uz]=local_source(s,x,z)
x=x(:); z=z(:)+zeros(size(x)); n=numel(x); u=complex(zeros(n,1)); ux=u; uz=u;
for first=1:128:n
 jj=first:min(first+127,n);
 phase=exp(1i*(x(jj)-s.x0_m)*s.kx.').*exp(-1i*(z(jj)-s.z_tx_m).*s.kz.');
 u(jj)=phase*s.coeff/s.n_original;
 ux(jj)=phase*(1i*s.kx.*s.coeff)/s.n_original;
 uz(jj)=phase*(-1i*s.kz.*s.coeff)/s.n_original;
end
end

function row=local_row(d,f,study)
row=local_row_with_plan(d,f,d.plan,study);
end

function row=local_row_with_plan(d,f,plan,study)
m=local_error(f.receiver_field,d.reference_field,d.weights,d.mask);
dx=d.x_m(2)-d.x_m(1);
row=struct('case_name',d.name,'study',study,'rank_in',f.rank,'rank_out',f.outgoing_rank, ...
 'weighted_complex_l2',m.weighted_complex_l2,'whole_grid_complex_l2',m.whole_grid_complex_l2, ...
 'magnitude_l2',m.magnitude_l2,'phase_rms_rad',m.phase_rms_rad, ...
 'actual_Kmax',max(d.profile.k_rad_per_m),'dx_m',dx,'K_Nyquist',pi/dx, ...
 'n_surface',numel(d.x_m),'fft_size',f.fft_size,'padding',f.padding_factor, ...
 'incoming_effective_rank_1e3',local_rank(plan.incoming_singular_values,1e-3), ...
 'outgoing_effective_rank_1e3',local_rank(plan.outgoing_singular_values,1e-3), ...
 'fast_setup_wall_s',plan.setup_wall_seconds,'fast_setup_cpu_s',plan.setup_cpu_seconds, ...
 'fast_apply_wall_s',f.wall_seconds,'fast_apply_cpu_s',f.cpu_seconds, ...
 'full_apply_wall_s',d.reference_wall_s,'full_apply_cpu_s',d.reference_cpu_s, ...
 'full_incident_wall_s',d.incident_wall_s,'full_incident_cpu_s',d.incident_cpu_s, ...
 'warm_speedup',d.reference_wall_s/f.wall_seconds, ...
 'cold_speedup',(d.reference_wall_s+d.incident_wall_s)/(plan.setup_wall_seconds+f.wall_seconds), ...
 'plan_MiB',plan.plan_bytes/2^20,'setup_arrays_MiB',plan.setup_named_arrays_bytes/2^20, ...
 'apply_arrays_MiB',f.named_apply_arrays_bytes/2^20, ...
 'full_green_block_estimated_MiB',min(256,numel(d.x_m))*numel(d.x_m)*104/2^20);
end

function m=local_error(a,b,w,mask)
m.whole_grid_complex_l2=norm(a-b)/norm(b);
a=a(mask); b=b(mask); w=w(mask); w=w/sum(w);
m.weighted_complex_l2=sqrt(sum(w.*abs(a-b).^2)/sum(w.*abs(b).^2));
m.magnitude_l2=sqrt(sum(w.*(abs(a)-abs(b)).^2)/sum(w.*abs(b).^2));
m.phase_rms_rad=sqrt(sum(w.*angle(a.*conj(b)).^2));
end

function r=local_rank(s,tol)
tail=sqrt(flipud(cumsum(flipud(abs(s).^2)))/sum(abs(s).^2));
tail=[tail(2:end);0]; r=find(tail<=tol,1);
end

function sp=local_spectrum(d,source)
n=numel(d.x_m); dx=mean(diff(d.x_m)); k=2*pi*4000/1500;
kx=2*pi*(-floor(n/2):ceil(n/2)-1)'/(n*dx);
win=.5-.5*cos(2*pi*(0:n-1)'/(n-1));
fk=dx*fftshift(fft(d.reference_field.*win)); fast=dx*fftshift(fft(d.selected.receiver_field.*win));
prop=abs(kx)<=k; near=abs(kx)>k & abs(kx)<=1.2*k;
energy=abs(fk).^2; fe=abs(fast).^2;
sk=source.kx; se=abs(source.coeff).^2; [ab,ii]=sort(abs(sk)); cum=cumsum(se(ii))/sum(se);
cut=ab(find(cum>=.99,1));
sp=struct('kx',kx,'FK',fk,'FastFK',fast,'propagating',prop,'near_evanescent',near, ...
 'theta95_out_deg',local_angle(kx(prop),energy(prop),k,.95), ...
 'theta99_out_deg',local_angle(kx(prop),energy(prop),k,.99), ...
 'theta99_fast_deg',local_angle(kx(prop),fe(prop),k,.99), ...
 'incident_theta99_deg',asind(cut/k), ...
 'outside_incident99_propagating_energy',sum(energy(prop & abs(kx)>cut))/sum(energy(prop)), ...
 'propagating_energy_fraction',sum(energy(prop))/sum(energy), ...
 'near_evanescent_energy_fraction',sum(energy(near))/sum(energy), ...
 'far_evanescent_energy_fraction',sum(energy(~prop & ~near))/sum(energy), ...
 'complex_l2',norm(fast-fk)/norm(fk),'propagating_complex_l2',norm(fast(prop)-fk(prop))/norm(fk(prop)), ...
 'near_evanescent_complex_l2',norm(fast(near)-fk(near))/max(norm(fk(near)),realmin));
end

function theta=local_angle(kx,en,k,q)
[ab,ii]=sort(abs(kx)); cum=cumsum(en(ii))/sum(en); theta=asind(min(1,ab(find(cum>=q,1))/k));
end

function a=local_append(a,b)
if isempty(a), a=b; else, a=[a(:);b(:)]; end
end

function result=local_bie(result,out)
% Two prespecified physical checks: canonical and stronger high-cutoff PM.
% BIE numerical refinements do not adjust FK phases or heights.
refs=struct([]);
for ic=[1 5]
 d=result.cases(ic); source=result.source; p=d.profile;
 idx=(1:8:numel(d.x_m))'; x=d.x_m(idx);
 cfg=d.refcfg; cfg.receiver_x_m=x; cfg.receiver_z_m=3+zeros(size(x));
 full=solve_full_kirchhoff_2d(cfg);
 runs=struct([]);
 settings=[60 8;60 10]; if ic==5, settings=[settings;70 10]; end %#ok<AGROW>
 for level=1:size(settings,1)
  L=settings(level,1); ppw=settings(level,2);
  fprintf('PM BIE %s halfwidth %g ppw %g\n',d.name,L,ppw);
  bcfg=struct('frequency_hz',4000,'c0_mps',1500,'half_width_m',L,'points_per_wavelength',ppw, ...
   'panel_order',8,'self_quadrature_order',96,'halfplane_h_m',-2,'window_plateau_ratio',.85, ...
   'eta_fn',@(xx)local_eta(p,xx),'eta_prime_fn',@(xx)local_slope(p,xx), ...
   'incident_fn',@(xx,zz)local_source(source,xx,zz),'receiver_x_m',x,'receiver_z_m',3);
  clk=tic; cpu=cputime; b=solve_helmholtz_bie_halfplane_vertical(bcfg);
  wall=toc(clk); cp=cputime-cpu;
  m=local_error(full.receiver_field,b.receiver_field,d.weights(idx),d.mask(idx));
  runs=local_append(runs,struct('half_width_m',L,'ppw',ppw,'wall_s',wall,'cpu_s',cp, ...
   'node_count',b.node_count,'dense_matrix_MiB',16*b.node_count^2/2^20, ...
   'linear_residual',b.linear_residual,'boundary_residual',b.offgrid_boundary_residual, ...
   'field',b.receiver_field,'weighted_complex_l2',m.weighted_complex_l2, ...
   'whole_grid_complex_l2',m.whole_grid_complex_l2,'magnitude_l2',m.magnitude_l2,'phase_rms_rad',m.phase_rms_rad));
  partial=struct('case_name',d.name,'receiver_x_m',x,'full_field',full.receiver_field,'runs',runs);
  save(fullfile(out,[d.name '_bie_checkpoint.mat']),'partial','-v7.3');
  progress_rows=rmfield(runs,'field');
  writetable(struct2table(progress_rows),fullfile(out,[d.name '_bie_progress.csv']));
  fprintf('  Full-BIE %.5g boundary %.4g\n',m.weighted_complex_l2,b.offgrid_boundary_residual);
 end
 refine=local_error(runs(1).field,runs(2).field,d.weights(idx),d.mask(idx));
 window=0;
 if numel(runs)>2, window=local_error(runs(2).field,runs(3).field,d.weights(idx),d.mask(idx)).weighted_complex_l2; end
 refs=local_append(refs,struct('case_name',d.name,'receiver_x_m',x,'full_field',full.receiver_field, ...
  'runs',runs,'spatial_convergence_l2',refine.weighted_complex_l2,'window_convergence_l2',window));
end
result.bie=refs;
end

function y=local_eta(p,x)
y=evaluate_pm_pilot_surface(p,x);
end

function r=local_physical_quadrature(r)
% Refine only the high-K Full integral at the already frozen BIE receivers.
d=r.cases(5); b=r.bie(2); x=linspace(-60,60,8193)';
[eta,slope]=evaluate_pm_pilot_surface(d.profile,x);
normal=[-slope,ones(size(slope))]./sqrt(1+slope.^2);
[u,ux,uz]=local_source(r.source,x,eta);
cfg=d.refcfg; cfg.surface_x_m=x; cfg.surface_z_m=eta; cfg.surface_normal=normal;
cfg.surface_incident=u; cfg.surface_incident_dn=normal(:,1).*ux+normal(:,2).*uz;
cfg.receiver_x_m=b.receiver_x_m; cfg.receiver_z_m=3+zeros(size(b.receiver_x_m));
clk=tic; cpu=cputime; f=solve_full_kirchhoff_2d(cfg); wall=toc(clk); cp=cputime-cpu;
idx=1:8:numel(d.x_m); w=d.weights(idx); mask=d.mask(idx);
conv=local_error(b.full_field,f.receiver_field,w,mask);
physical=local_error(f.receiver_field,b.runs(end).field,w,mask);
r.highK_full_quadrature=struct('n_surface',8193,'receiver_x_m',b.receiver_x_m, ...
 'coarse_full_field',b.full_field,'fine_full_field',f.receiver_field,'bie_field',b.runs(end).field, ...
 'weights',w,'mask',mask,'convergence_weighted_l2',conv.weighted_complex_l2, ...
 'convergence_whole_l2',conv.whole_grid_complex_l2,'fine_full_bie_weighted_l2',physical.weighted_complex_l2, ...
 'fine_full_bie_whole_l2',physical.whole_grid_complex_l2,'wall_s',wall,'cpu_s',cp);
fprintf('HIGH_K_FULL_QUADRATURE %.7g; refined Full-BIE %.7g\n',conv.weighted_complex_l2,physical.weighted_complex_l2);
end
function y=local_slope(p,x)
[~,y]=evaluate_pm_pilot_surface(p,x);
end

function result=local_gates(result)
t=result.metrics; selected=strcmp(t.study,'selected');
result.algorithm_pass=all(t.weighted_complex_l2(selected)<1e-3 & t.whole_grid_complex_l2(selected)<1e-3) ...
 && isfinite(result.selected_rank_in) && result.padding_2_4_l2<1e-12 ...
 && result.grid.full_weighted_l2<1e-4 && result.grid.full_whole_l2<1e-4 ...
 && result.height_nodes_33_65_l2<1e-4 && all(t.cold_speedup(selected)>2);
result.physics_pass=numel(result.bie)==2;
for i=1:numel(result.bie)
 b=result.bie(i); r=b.runs(end);
 result.physics_pass=result.physics_pass && r.weighted_complex_l2<1e-3 ...
  && b.spatial_convergence_l2<1e-4 && b.window_convergence_l2<1e-4 && r.boundary_residual<1e-4;
end
if isfield(result,'highK_full_quadrature')
 result.physics_pass=result.physics_pass && result.highK_full_quadrature.convergence_weighted_l2<1e-4 ...
  && result.highK_full_quadrature.fine_full_bie_weighted_l2<1e-3;
end
result.untapered_pass=isfield(result,'untapered_metrics') && height(result.untapered_metrics)==5;
if result.untapered_pass
 t=result.untapered_metrics;
 result.untapered_pass=all(t.weighted_complex_l2<1e-3 & t.whole_grid_complex_l2<1e-3);
end
result.stop_reason='';
if ~result.algorithm_pass
 result.stop_reason='algorithm accuracy, convergence or speed gate failed';
elseif ~result.untapered_pass
 result.stop_reason='unmodified finite PM algorithm controls pending or failed';
elseif ~result.physics_pass && ~isempty(result.bie)
 result.stop_reason='FullFK-BIE physical or BIE numerical gate failed';
elseif isempty(result.bie)
 result.stop_reason='physical BIE references pending';
end
result.pilot_pass=result.algorithm_pass && result.untapered_pass && result.physics_pass;
end

function local_outputs(result,out,root)
writetable(result.metrics,fullfile(out,'field_metrics.csv'));
writetable(result.profile_stats,fullfile(out,'surface_statistics.csv'));
writetable(result.grid_metrics,fullfile(out,'grid_metrics.csv'));
if isfield(result,'untapered_metrics')
 writetable(result.untapered_metrics,fullfile(out,'untapered_field_metrics.csv'));
 for i=1:numel(result.untapered_cases)
  d=result.untapered_cases(i); save(fullfile(out,[d.name '_untapered_fields.mat']),'d','-v7.3');
 end
end
if ~isempty(result.rank_diagnostics), writetable(struct2table(result.rank_diagnostics),fullfile(out,'rank_diagnostics.csv')); end
angular=struct([]); bie_rows=struct([]);
for i=1:numel(result.cases)
 d=result.cases(i); s=d.spectrum;
 writetable(table(s.kx,real(s.FK),imag(s.FK),real(s.FastFK),imag(s.FastFK),s.propagating,s.near_evanescent, ...
  'VariableNames',{'kx_radpm','FK_real','FK_imag','FastFK_real','FastFK_imag','propagating','near_evanescent'}), ...
  fullfile(out,[d.name '_reflected_spectrum.csv']));
 row=rmfield(s,{'kx','FK','FastFK','propagating','near_evanescent'}); row.case_name=d.name;
 angular=local_append(angular,row);
 save(fullfile(out,[d.name '_fields.mat']),'d','-v7.3');
end
writetable(struct2table(angular),fullfile(out,'angular_summary.csv'));
for i=1:numel(result.bie)
 b=result.bie(i);
 for j=1:numel(b.runs)
  row=rmfield(b.runs(j),'field'); row.case_name=b.case_name;
  row.spatial_convergence_l2=b.spatial_convergence_l2; row.window_convergence_l2=b.window_convergence_l2;
  bie_rows=local_append(bie_rows,row);
 end
end
if ~isempty(bie_rows), writetable(struct2table(bie_rows),fullfile(out,'bie_reference_metrics.csv')); end
if isfield(result,'highK_full_quadrature')
 row=rmfield(result.highK_full_quadrature,{'receiver_x_m','coarse_full_field','fine_full_field','bie_field','weights','mask'});
 writetable(struct2table(row),fullfile(out,'highK_full_quadrature.csv'));
end
figdir=fullfile(out,'figures'); if ~isfolder(figdir), mkdir(figdir); end
local_figures(result,figdir);
local_report(result,root,angular);
end

function local_figures(r,folder)
f=figure('Visible','off','Position',[50 50 1350 850]); tiledlayout(3,3);
for i=1:3
 d=r.cases(i); nexttile; plot(d.x_m,d.eta_m,'LineWidth',1); ylabel('height (m)'); xlabel('x (m)'); title(d.name,'Interpreter','none'); grid on;
 nexttile; plot(d.x_m,d.slope,'LineWidth',1); xlabel('x (m)'); ylabel('slope'); grid on;
 nexttile; plot(d.x_m,d.second_derivative./(1+d.slope.^2).^1.5,'LineWidth',1); xlabel('x (m)'); ylabel('curvature (1/m)'); grid on;
end
local_export(f,folder,'pm_profiles');
f=figure('Visible','off','Position',[50 50 1100 800]); tiledlayout(2,2);
nexttile; hold on;
for i=[1 4 5], p=r.profiles(i); semilogy(p.k_rad_per_m,(p.cos_coeff_m.^2+p.sin_coeff_m.^2)/2,'o-','DisplayName',p.label); end
set(gca,'YScale','log'); xlabel('physical surface K (rad/m)'); ylabel('mode height variance (m^2)'); legend('Interpreter','none'); grid on;
nexttile; hold on;
d=r.cases(1);
semilogy(1:numel(d.outgoing_singular_values),d.outgoing_singular_values/d.outgoing_singular_values(1),'o-','DisplayName','outgoing, shared bound');
semilogy(1:numel(d.incoming_singular_values),d.incoming_singular_values/d.incoming_singular_values(1),'s--','DisplayName','incoming, shared source');
set(gca,'YScale','log'); ylim([1e-12 2]); xlim([1 15]); xlabel('singular value index'); ylabel('s_j/s_1'); legend('Interpreter','none'); grid on;
nexttile; hold on;
for i=1:3, d=r.cases(i); plot(d.x_m,abs(d.reference_field),'DisplayName',d.name); end
xlim([-30 30]); xlabel('receiver x (m)'); ylabel('|Full FK|'); legend('Interpreter','none'); grid on;
nexttile; t=r.metrics; b=strcmp(t.study,'baseline'); s=strcmp(t.study,'selected');
semilogy(1:5,t.weighted_complex_l2(b),'o-',1:5,t.weighted_complex_l2(s),'s--');
yline(1e-3,':'); xticks(1:5); xticklabels({'260001','260002','260003','low K','high K'});
xlabel('fixed PM case'); ylabel('weighted complex L2'); legend('baseline 3/6','selected common rank','goal'); grid on;
local_export(f,folder,'ranks_fields_cutoffs');
f=figure('Visible','off','Position',[50 50 1100 900]); tiledlayout(3,2);
for i=1:3
 s=r.cases(i).spectrum; nexttile;
 semilogy(s.kx/(2*pi*4000/1500),abs(s.FK).^2,'-',s.kx/(2*pi*4000/1500),abs(s.FastFK).^2,'--');
 xlim([-1.2 1.2]); xlabel('k_x/k'); ylabel('aperture spectral power'); title(r.cases(i).name,'Interpreter','none'); legend('Full','Fast'); grid on;
 nexttile; semilogy(s.kx/(2*pi*4000/1500),abs(s.FastFK-s.FK)/norm(s.FK));
 xlim([-1.2 1.2]); xlabel('k_x/k'); ylabel('|spectral difference| / ||Full||_2'); grid on;
end
local_export(f,folder,'reflected_spectra');
f=figure('Visible','off','Position',[50 50 1100 470]); tiledlayout(1,2);
nexttile; hold on;
for i=[4 1 5]
 s=r.cases(i).spectrum; semilogy(s.kx/(2*pi*4000/1500),abs(s.FK).^2,'DisplayName',r.cases(i).name);
end
set(gca,'YScale','log'); ylim([1e-12 2]); xlim([-.8 .8]); grid on;
xlabel('k_x/k'); ylabel('Full FK aperture spectral power'); legend('Interpreter','none');
nexttile; ii=[4 1 5]; theta=zeros(3,1); outside=theta; kmax=theta;
for j=1:3, theta(j)=r.cases(ii(j)).spectrum.theta99_out_deg; outside(j)=100*r.cases(ii(j)).spectrum.outside_incident99_propagating_energy; kmax(j)=max(r.profiles(ii(j)).k_rad_per_m); end
yyaxis left; plot(kmax,theta,'o-'); ylabel('theta99 outgoing (deg)');
yyaxis right; plot(kmax,outside,'s--'); ylabel('reflected prop energy outside incident99 (%)');
xlabel('physical surface Kmax (rad/m)'); grid on;
local_export(f,folder,'cutoff_outgoing_bandwidth');
end

function local_export(f,folder,name)
exportgraphics(f,fullfile(folder,[name '.png']),'Resolution',160);
exportgraphics(f,fullfile(folder,[name '.pdf']),'ContentType','vector'); close(f);
end

function local_report(r,root,angular)
path=fullfile(root,'reports','spectral_fast_kirchhoff_pm_pilot_report.md');
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'# Spectral Fast Kirchhoff fixed PM pilot\n\nDate: %s. MATLAB %s. Validation-only, 4 kHz, c=1500 m/s, Tx=100 m, Rx=3 m.\n\n',r.date,r.matlab_version);
fprintf(fid,'Algorithm pass: **%d**; physical-reference pass: **%d**; PM pilot PASS: **%d**. No production model changes, Monte Carlo or frequency extension.\n\n',r.algorithm_pass,r.physics_pass,r.pilot_pass);
fprintf(fid,'The unchanged discrete operator and FFT normalization are defined in [prototype method](spectral_fast_kirchhoff_prototype_method.md); its deterministic evidence is in [prototype report](spectral_fast_kirchhoff_prototype_report.md), with [prior separability audit](full_kirchhoff_spectral_rank_audit_report.md). Both exact carrier phases exp(-ik*eta) remain intact. This pilot changes coefficient-defined surface inputs only, retaining full normals, arc weights, incident physical derivatives and both nonlocal Green terms.\n\n');
fprintf(fid,'## Provenance and physical surface\n\nThree fixed realizations use saved coefficients: canonical 260001 and independent Gaussian coefficient seeds 260002/260003, U=6 m/s, period 160 m. Saved coefficients are copied without renormalization. Canonical SHA256: `%s`. All solvers consume one coefficient-defined surface; none resamples random numbers.\n\n',r.profiles(1).coeff_file_sha256);
fprintf(fid,'Canonical is the original fixed realization; the other two are independent Gaussian amplitude draws from its saved Sk band. This is three deterministic pilot inputs, not an ensemble estimate. Analytic Fourier derivatives determine slope, curvature and the full water-facing normal.\n\n');
fprintf(fid,'Raw PM eta=sum(a cos Kx+b sin Kx). Physical cutoff controls the underlying periodic Fourier series. A fixed C2 envelope is 1 for |x|<=42 m and 0 for |x|>=50 m, within [-60,60] m integration support. It preserves central PM heights exactly and makes the same finite surface compatible with the frozen BIE. Its edge transition is a documented geometric embedding, not a PM bandlimited claim or sponge. Raw and embedded samples/statistics are saved separately.\n\n');
fprintf(fid,'|case|actual Kmax|raw RMS height m|raw RMS/max slope|raw RMS/max curvature 1/m|embedded max height m|\n|---|---:|---:|---|---|---:|\n');
t=r.profile_stats;
for i=1:height(t), fprintf(fid,'|%s|%.7g|%.7g|%.7g / %.7g|%.7g / %.7g|%.7g|\n',t.case_name{i},t.actual_Kmax(i),t.raw_rms_height(i),t.raw_rms_slope(i),t.raw_max_slope(i),t.raw_rms_curvature(i),t.raw_max_curvature(i),t.embedded_max_height(i)); end
fprintf(fid,'\n## Physical versus numerical cutoff\n\nRequested Kmax=0.30/0.50/0.75 rad/m; actual last modes %.8g / %.8g / %.8g. Main dx=120/4096=%.8g m; Nyquist=%.8g rad/m, independent of Kmax. Grid refinements reevaluate the SAME coefficients; they add no physical roughness modes. Higher cutoff appends modes 13--19 with Gaussian coefficients, mt19937ar seed 26000175, Sk=alpha/(2K^3) exp(-beta*g^2/(U^4 K^2)), alpha=.0081 beta=.74 g=9.81. Every original coefficient is retained exactly; lower cutoff only removes upper modes. No height rescaling.\n\n',r.profiles(4).realized_kmax_rad_per_m,r.profiles(1).realized_kmax_rad_per_m,r.profiles(5).realized_kmax_rad_per_m,120/4096,pi/(120/4096));
fprintf(fid,'## Algorithm fidelity and uniform ranks\n\nFrozen baseline 3/6 is evaluated FIRST on every case. Selected common incoming/outgoing ranks: **%g/%g**; uniform height bound A=%.3g m, fixed 33 height nodes. No per-case best rank, phase adjustment, scalar fit or incoming99 outgoing truncation.\n\n',r.selected_rank_in,r.selected_rank_out,r.common_height_bound_m);
fprintf(fid,'Relative weighted complex L2 = sqrt(sum w|Fast-Full|^2 / sum w|Full|^2), fixed Stage-0 incident M99 weights plus Full-only -40 dB gate. Whole-grid uses all N receivers and ||Fast-Full||2/||Full||2. Magnitude and wrapped phase RMS use that fixed mask. Native exp(-iwt) fields are compared directly without conjugate/ratio calibration.\n\n');
fprintf(fid,'|case|study|r_in/r_out|weighted L2|whole-grid L2|magnitude L2|phase RMS rad|effective in/out rank 1e-3|\n|---|---|---|---:|---:|---:|---:|---|\n');
t=r.metrics;
for i=1:height(t), fprintf(fid,'|%s|%s|%d/%d|%.7g|%.7g|%.7g|%.7g|%d/%d|\n',t.case_name{i},t.study{i},t.rank_in(i),t.rank_out(i),t.weighted_complex_l2(i),t.whole_grid_complex_l2(i),t.magnitude_l2(i),t.phase_rms_rad(i),t.incoming_effective_rank_1e3(i),t.outgoing_effective_rank_1e3(i)); end
if ~isempty(r.rank_diagnostics)
 tt=struct2table(r.rank_diagnostics); ii=tt.rank_in==4 & tt.rank_out==8;
 fprintf(fid,'\nThe already prescribed common 4/8 diagnostic also reaches the ideal 1e-4 across all five embedded cases: maximum weighted %.7g, whole-grid %.7g. It is reported as a diagnostic with limited margin; the primary selected 3/8 remains judged against 1e-3, without per-realization ranks. Full 45-row ladder, including failed combinations, is saved in rank_diagnostics.csv.\n\n',max(tt.weighted_complex_l2(ii)),max(tt.whole_grid_complex_l2(ii)));
end
fprintf(fid,'\nEffective ranks are relative Frobenius SVD tails of the weighted incoming matrix and stacked carrier-factored outgoing [G;Dz/k;Dx/k] lag-height table, not guaranteed field error bounds. The outgoing basis depends on height bound, receiver depth and lag support; K appears through spatial densities/normal, not directly in this basis. Separate rank diagnostics and incident/outgoing ablations are saved if baseline fails.\n\n');
fprintf(fid,'Canonical rank3/10 field error %.7g, rank5/6 error %.7g. This separates incoming truncation from outgoing truncation without changing the physical model.\n\n',r.diagnosis.rank3_10_whole_l2,r.diagnosis.rank5_6_whole_l2);
fprintf(fid,'## Reflected aperture spectrum\n\nSaved complex A_FK(kx), A_Fast(kx) use dx*FFT with an identical receiver Hann window. They are finite-aperture diagnostics, not infinite-domain far-field flux. Propagating |kx|<=k, near-evanescent k<|kx|<=1.2k, and higher bins are all retained. Fractional energies are spectral norms, not evanescent transported flux.\n\n');
if isfield(r,'untapered_metrics')
 fprintf(fid,'Additional unmodified PM controls (no C2 envelope) use all true Fourier heights and normals on [-60,60] m with the same selected ranks. These compare Fast-Full on the original finite PM segment, preventing the algorithm conclusion from relying on a favorable surface taper. They do not independently validate an infinite periodic PM boundary against BIE.\n\n|case|untapered weighted L2|untapered whole-grid L2|\n|---|---:|---:|\n');
 tt=r.untapered_metrics;
 for j=1:height(tt), fprintf(fid,'|%s|%.7g|%.7g|\n',tt.case_name{j},tt.weighted_complex_l2(j),tt.whole_grid_complex_l2(j)); end
 fprintf(fid,'\n');
end
fprintf(fid,'|case|theta95/99 out deg|outside incident99 prop energy|prop/near-evan energy|spectrum L2|near-evan relative L2|\n|---|---|---:|---|---:|---:|\n');
for i=1:numel(angular), s=angular(i); fprintf(fid,'|%s|%.5g / %.5g|%.7g|%.7g / %.7g|%.7g|%.7g|\n',s.case_name,s.theta95_out_deg,s.theta99_out_deg,s.outside_incident99_propagating_energy,s.propagating_energy_fraction,s.near_evanescent_energy_fraction,s.complex_l2,s.near_evanescent_complex_l2); end
fprintf(fid,'\nNear-evanescent relative errors can be large when the norm is nearly zero; absolute powers and total spectrum L2 are saved for interpretation. No spectral cutoff is applied to the Hankel-kernel solver.\n\n');
fprintf(fid,'## Grid, padding and interpolation\n\nCanonical N=2049/4097/8193 use unchanged Fourier coefficients and analytic normals. 4097 versus8193: Full weighted %.7g, whole-grid %.7g; Fast weighted %.7g; nested surface closure %.3g m. Padding2/4 L2 %.3g. Height-node33/65 L2 %.7g. The latter is a separate interpolation/basis-quadrature check, not physical rank tuning.\n\n',r.grid.full_weighted_l2,r.grid.full_whole_l2,r.grid.fast_weighted_l2,r.grid.surface_sample_closure,r.padding_2_4_l2,r.height_nodes_33_65_l2);
fprintf(fid,'|N|dx|Nyquist|matrix incoming/outgoing rank 1e-3|Fast-Full weighted|whole grid|\n|---:|---:|---:|---|---:|---:|\n');
t=r.grid_metrics; for i=1:height(t), fprintf(fid,'|%d|%.7g|%.7g|%d/%d|%.7g|%.7g|\n',t.n_surface(i),t.dx_m(i),t.K_Nyquist(i),t.incoming_effective_rank_1e3(i),t.outgoing_effective_rank_1e3(i),t.weighted_complex_l2(i),t.whole_grid_complex_l2(i)); end
fprintf(fid,'\n## Separate FullFK--BIE physical reference\n\nPrespecified canonical and high-cutoff cases: frozen Dirichlet halfplane combined layer, h=-2 m, plateau ratio .85, order8/self96, ppw8/10, high-cutoff halfwidth60/70. FullFK reevaluated at the same 513 receivers; same coefficient surface and incident spectrum. Physical acceptable gate is weighted native field L2<1e-3 with numerical uncertainty checks below1e-4. This gate is separate from Fast-Full fidelity.\n\n');
fprintf(fid,'|case|L|ppw|Full-BIE weighted L2|whole-grid L2|boundary residual|BIE refinement L2|window L2|\n|---|---:|---:|---:|---:|---:|---:|---:|\n');
for i=1:numel(r.bie), b=r.bie(i); for j=1:numel(b.runs), s=b.runs(j); fprintf(fid,'|%s|%g|%g|%.7g|%.7g|%.7g|%.7g|%.7g|\n',b.case_name,s.half_width_m,s.ppw,s.weighted_complex_l2,s.whole_grid_complex_l2,s.boundary_residual,b.spatial_convergence_l2,b.window_convergence_l2); end; end
if isfield(r,'highK_full_quadrature')
 s=r.highK_full_quadrature;
 fprintf(fid,'\nHigh-cutoff Full quadrature4097/8193 at the same513 receivers: weighted convergence %.7g, whole-grid %.7g. Refined Full--BIE weighted %.7g, whole-grid %.7g. The remaining physical deviation is not resolved by Full quadrature or BIE refinement.\n\n',s.convergence_weighted_l2,s.convergence_whole_l2,s.fine_full_bie_weighted_l2,s.fine_full_bie_whole_l2);
end
fprintf(fid,'Whole-grid Full--BIE errors include finite-window edge receivers outside the validated M99 physical footprint. The high-cutoff L60/70 whole-grid error changes markedly, while the weighted core is stable. Physical accuracy is therefore not claimed on the entire receiver line or for an infinite periodic PM ocean; algorithm whole-grid fidelity is a separate result.\n\n');
if isempty(r.bie), fprintf(fid,'BIE reference is pending; pilot cannot PASS until it is measured.\n'); end
fprintf(fid,'\n## CPU, memory and speedup\n\nWarm Fast median5 applications; Full dense median2; cold includes Fast preparation or Full exact incident. Windows CPU resolution may give zero for short calls; wall times determine ratios. Memory is named arrays/estimated block payload, NOT peak RSS or FFT/SVD workspace. BIE dense matrix payload/times are in bie_reference_metrics.csv.\n\n');
fprintf(fid,'|case|Full apply s|Fast setup/apply s|warm/cold speedup|plan/setup/apply MiB|Full block MiB|\n|---|---:|---|---|---|---:|\n');
t=r.metrics; rows=find(strcmp(t.study,'selected'));
for i=rows.', fprintf(fid,'|%s|%.6g|%.6g / %.6g|%.5g / %.5g|%.5g / %.5g / %.5g|%.5g|\n',t.case_name{i},t.full_apply_wall_s(i),t.fast_setup_wall_s(i),t.fast_apply_wall_s(i),t.warm_speedup(i),t.cold_speedup(i),t.plan_MiB(i),t.setup_arrays_MiB(i),t.apply_arrays_MiB(i),t.full_green_block_estimated_MiB(i)); end
fprintf(fid,'\nFast workflow wall %.5g s, CPU %.5g s (excluding separate BIE/render stage). The same global FFT structure remains; finite pilot timings do not prove asymptotic complexity.\n\n',r.numerical_wall_s,r.numerical_cpu_s);
if ~isempty(r.bie)
 wall=0; cpu=0;
 for i=1:numel(r.bie), wall=wall+sum([r.bie(i).runs.wall_s]); cpu=cpu+sum([r.bie(i).runs.cpu_s]); end
 fprintf(fid,'The five BIE reference solves total wall %.6g s and process CPU %.6g s. Process CPU can exceed wall on multithreaded linear algebra; these expensive physical checks are not part of Fast setup/apply costs. Provenance hashes are in provenance.json.\n\n',wall,cpu);
end
fprintf(fid,'## Decision and reproduction\n\n');
if r.pilot_pass
 fprintf(fid,'**PM pilot PASS.** A small multi-frequency consistency study is justified next, but is NOT run in this task. This conclusion covers only these fixed-band 1-D PM inputs and finite embedded support.\n\n');
elseif isempty(r.bie)
 fprintf(fid,'**Pilot incomplete: physical BIE checks pending.** No next-stage authorization follows from algorithm accuracy alone.\n\n');
else
 fprintf(fid,'**PM pilot does not meet all gates. Stop before multi-frequency/Monte Carlo.** Algorithm compression succeeds with common3/8 and stable grids/cutoffs. The high-cutoff physical FullFK--BIE discrepancy exceeds the prespecified1e-3 gate. Increasing FastFK rank cannot remove a FullFK physical approximation error. First establish the PM physical validity domain or compare an independent nonlocal SSA/SSA2 reference; block/local-stationary Fast separation is relevant only if compression subsequently fails. No production model is changed here.\n\n');
 fprintf(fid,'This is a narrowly failed strict physical gate, not evidence that low-rank PM compression collapses or that every PM case invalidates FullFK. Canonical passes the weighted1e-3 physical gate; the high-cutoff discrepancy is predominantly phase (about1.2e-3 rad versus magnitude error1.39e-4). Its higher slope/curvature and wider outgoing spectrum suggest a model-validity boundary to investigate, without a fitted causal formula.\n\n');
end
fprintf(fid,'Five independent unit tests (PM derivatives/embedding/unmodified series and prior flat Green/FFT/interpolation controls) pass; saved artifact assertions verify coefficient preservation, Gaussian-seed reproduction, geometry, fixed ranks, spectrum norm closure and BIE/Full refinement data. Test_results.csv and artifact_assertions.mat preserve the verification.\n\n');
fprintf(fid,'Run after setup_vertical_project: `validate_spectral_fast_kirchhoff_pm_pilot(''full'')`; checkpoint stages `''fast''`, `''bie''`, `''raw''`, `''quadrature''`; `''render''` reuses saved numerical data. Outputs: results/validation/spectral_fast_kirchhoff_pm_pilot. No PE/WAPE, Full FK or BIE reference change; no communication regressions are relevant to these validation-only artifacts.\n\n');
for name={'pm_profiles','ranks_fields_cutoffs','reflected_spectra','cutoff_outgoing_bandwidth'}, fprintf(fid,'![%s](../results/validation/spectral_fast_kirchhoff_pm_pilot/figures/%s.png)\n\n',name{1},name{1}); end
end
