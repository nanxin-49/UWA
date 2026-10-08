function s=validate_spectral_fast_kirchhoff_multifrequency_consistency(mode)
%VALIDATE_SPECTRAL_FAST_KIRCHHOFF_MULTIFREQUENCY_CONSISTENCY Frozen PM audit.
% init; next (one sparse frequency); rank (common rank ladder only if failed);
% anchor (one BIE); guard (one frequency); dense (one frequency after gates);
% render/status. References are immutable; each expensive operation checkpoints.
if nargin<1, mode='status'; end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root); setup_vertical_project();
out=fullfile(root,'results','validation','spectral_fast_kirchhoff_multifrequency_consistency');
if ~isfolder(out), mkdir(out); end
file=fullfile(out,'consistency.mat');
if strcmp(mode,'init')
 assert(~isfile(file),'Existing audit cannot be overwritten.');
 p=load(fullfile(root,'results','validation','spectral_fast_kirchhoff_pm_pilot','seed_260001_fields.mat'),'d');
 v=load(fullfile(root,'results','validation','full_kirchhoff_pm_validity_boundary','full_kirchhoff_pm_validity_boundary.mat'),'state');
 s=struct('schema_version','1.0','date','2026-10-08','matlab_version',version, ...
  'profile',p.d.profile,'source',p.d.cfg.spectral_source,'height_bound_m',p.d.cfg.height_bound_m, ...
  'x_m',p.d.x_m,'eta_m',p.d.eta_m,'slope',p.d.slope,'second',p.d.second_derivative, ...
  'weights',p.d.weights,'mask',p.d.mask,'sparse_hz',(4000:1000:8000)', ...
  'dense_hz',(4000:62.5:8000)','ranks',[3 8],'padding',2, ...
  'sparse',struct([]),'dense',struct([]),'anchors',struct([]),'guards',struct([]), ...
  'rank_trials',struct([]),'reference_hashes',local_hashes(root), ...
  'geometry',struct('z_tx',100,'z_rx',3,'z_surface',0,'c0',1500));
 assert(isequal(s.x_m,v.state.surface_x_m));
 assert(isequal(s.profile.k_rad_per_m,v.state.canonical_coefficients.k_rad_per_m));
 assert(isequal(s.profile.cos_coeff_m,v.state.canonical_coefficients.cos_coeff_m));
 assert(isequal(s.profile.sin_coeff_m,v.state.canonical_coefficients.sin_coeff_m));
 [eta,sl,sec]=evaluate_pm_pilot_surface(s.profile,s.x_m);
 assert(isequal(eta,s.eta_m) && isequal(sl,s.slope) && isequal(sec,s.second));
 copyfile(s.profile.coeff_file,fullfile(out,'canonical_coefficients.csv'));
 writetable(table(s.x_m,eta,sl,sec),fullfile(out,'fixed_surface.csv'));
else
 q=load(file,'s'); s=q.s;
 assert(isequal(s.reference_hashes,local_hashes(root)),'Frozen input/reference hash changed.');
end
for family={'sparse','dense','frozen_rank_sparse'}
 if ~isfield(s,family{1}), continue; end
 for j=1:numel(s.(family{1}))
  k=2*pi*s.(family{1})(j).frequency_hz/1500;
  s.(family{1})(j).H_full_native_dsp=conj(s.(family{1})(j).H_full_physical)*exp(1i*k*97);
  s.(family{1})(j).H_fast_native_dsp=conj(s.(family{1})(j).H_fast_physical)*exp(1i*k*97);
  s.(family{1})(j).H_direct_native_dsp=conj(s.(family{1})(j).H_direct_dsp);
 end
end
switch mode
 case 'next'
  i=numel(s.sparse)+1;
  if i<=numel(s.sparse_hz), s.sparse=local_append(s.sparse,local_case(s,s.sparse_hz(i),4097,33,s.ranks)); end
 case 'rank'
  assert(numel(s.sparse)==5 && ~local_algorithm(s.sparse),'Rank changes require an explicit fidelity failure.');
  if ~isfield(s,'frozen_rank_sparse'), s.frozen_rank_sparse=s.sparse; end
  candidates=[3 10;4 10;5 10;3 12;3 14;3 16];
  for pair=candidates.'
   if ~isempty(s.rank_trials) && any([s.rank_trials.rank_in]==pair(1) & [s.rank_trials.rank_out]==pair(2)), continue; end
   trials=struct([]);
   for i=1:5
    d=s.sparse(i); cfg=local_cfg(s,d.frequency_hz,4097,33,pair');
    plan=prepare_spectral_fast_kirchhoff_2d(cfg); f=solve_spectral_fast_kirchhoff_2d(cfg,plan);
    t=local_metric(f.receiver_field,d.full,s.weights,s.mask);
    s.rank_trials=local_append(s.rank_trials,struct('frequency_hz',d.frequency_hz,'rank_in',pair(1), ...
     'rank_out',pair(2),'whole',t.whole,'weighted',t.weighted));
    d=local_fast_record(d,f,plan,s); trials=local_append(trials,d);
   end
   if local_algorithm(trials), s.ranks=pair'; s.sparse=trials; break; end
  end
 case 'anchor'
  assert(numel(s.sparse)==5,'Complete sparse sweep first.');
  freqs=[4000 6000 8000]; i=numel(s.anchors)+1;
  if i<=3, s.anchors=local_append(s.anchors,local_anchor(s,freqs(i))); end
 case 'timing'
  assert(numel(s.sparse)==5);
  clk=tic; cpu=cputime;
  [eta,sl,sec]=evaluate_pm_pilot_surface(s.profile,s.x_m);
  normal=[-sl,ones(size(sl))]./sqrt(1+sl.^2);
  ds=hypot(diff(s.x_m),diff(eta)); qw=[ds(1)/2;(ds(1:end-1)+ds(2:end))/2;ds(end)/2];
  s.geometry_setup_wall=toc(clk); s.geometry_setup_cpu=cputime-cpu;
  s.geometry_arrays_MiB=8*(numel(eta)+numel(sl)+numel(sec)+numel(normal)+numel(qw))/2^20;
  for i=1:5
   cfg=local_cfg(s,s.sparse(i).frequency_hz,4097,33,s.ranks);
   plan=prepare_spectral_fast_kirchhoff_2d(cfg); solve_spectral_fast_kirchhoff_2d(cfg,plan);
   wt=zeros(3,1); ct=wt;
   for j=1:3
    fast=solve_spectral_fast_kirchhoff_2d(cfg,plan); wt(j)=fast.wall_seconds; ct(j)=fast.cpu_seconds;
   end
   fast.wall_seconds=median(wt); fast.cpu_seconds=median(ct);
   s.sparse(i)=local_fast_record(s.sparse(i),fast,plan,s);
  end
 case 'guard'
  assert(numel(s.sparse)==5,'Complete sparse sweep first.');
  freqs=[4000 6000 8000]; i=numel(s.guards)+1;
  if i<=3, s.guards=local_append(s.guards,local_guard(s,freqs(i))); end
 case 'dense'
  s=local_summary(s); assert(s.sparse_pass,'Sparse algorithm/physics/numerics/speed gates must pass before dense sweep.');
  i=numel(s.dense)+1;
  if i<=65
   f=s.dense_hz(i); old=find([s.sparse.frequency_hz]==f,1);
   if isempty(old), d=local_case(s,f,4097,33,s.ranks); else, d=s.sparse(old); end
   s.dense=local_append(s.dense,d);
  end
 case {'init','status','render'}
 otherwise, error('Unknown mode %s.',mode);
end
s=local_summary(s); save(file,'s','-v7.3');
local_outputs(s,out,root,strcmp(mode,'render'));
fprintf('MULTIFREQUENCY sparse=%d/5 BIE=%d/3 guards=%d/3 dense=%d/65 rank=%d/%d sparse_pass=%d PASS=%d\n', ...
 numel(s.sparse),numel(s.anchors),numel(s.guards),numel(s.dense),s.ranks,s.sparse_pass,s.pass);
end

function cfg=local_cfg(s,f,n,q,ranks)
if n==numel(s.x_m)
 x=s.x_m; eta=s.eta_m; sl=s.slope;
else
 x=linspace(-60,60,n)'; [eta,sl]=evaluate_pm_pilot_surface(s.profile,x);
end
src=s.source; k=2*pi*f/1500; src.kz=sqrt(k^2-src.kx.^2);
assert(all(isreal(src.kz)) && all(abs(src.kx)<k));
cfg=struct('frequency_hz',f,'c0_mps',1500,'surface_x_m',x,'surface_z_m',eta, ...
 'surface_normal',[-sl,ones(n,1)]./sqrt(1+sl.^2),'height_bound_m',s.height_bound_m, ...
 'receiver_depth_m',3,'spectral_source',src,'max_rank',max([10 ranks]), ...
 'height_nodes',q,'geometry','B','rank',ranks(1),'outgoing_rank',ranks(2),'padding_factor',2);
end

function [u,ux,uz]=local_source(src,x,z)
x=x(:); z=z(:)+zeros(size(x)); u=complex(zeros(size(x))); ux=u; uz=u;
for first=1:128:numel(x)
 jj=first:min(first+127,numel(x));
 phase=exp(1i*(x(jj)-src.x0_m)*src.kx.').*exp(-1i*(z(jj)-src.z_tx_m).*src.kz.');
 u(jj)=phase*src.coeff/src.n_original;
 ux(jj)=phase*(1i*src.kx.*src.coeff)/src.n_original;
 uz(jj)=phase*(-1i*src.kz.*src.coeff)/src.n_original;
end
end

function [ref,cost]=local_full(cfg,rx)
clk=tic; cpu=cputime; [u,ux,uz]=local_source(cfg.spectral_source,cfg.surface_x_m,cfg.surface_z_m);
cost.setup_wall=toc(clk); cost.setup_cpu=cputime-cpu;
r=struct('frequency_hz',cfg.frequency_hz,'c0_mps',1500,'surface_x_m',cfg.surface_x_m, ...
 'surface_z_m',cfg.surface_z_m,'surface_normal',cfg.surface_normal,'surface_incident',u, ...
 'surface_incident_dn',cfg.surface_normal(:,1).*ux+cfg.surface_normal(:,2).*uz, ...
 'receiver_x_m',rx,'receiver_z_m',3+zeros(size(rx)));
clk=tic; cpu=cputime; ref=solve_full_kirchhoff_2d(r);
cost.apply_wall=toc(clk); cost.apply_cpu=cputime-cpu;
end

function d=local_case(s,f,n,q,ranks)
fprintf('FREQUENCY %g Hz grid%d rank%d/%d\n',f,n,ranks);
cfg=local_cfg(s,f,n,q,ranks); plan=prepare_spectral_fast_kirchhoff_2d(cfg);
[ref,cost]=local_full(cfg,cfg.surface_x_m);
times=zeros(3,1); cpus=times;
solve_spectral_fast_kirchhoff_2d(cfg,plan);
for j=1:3
 fast=solve_spectral_fast_kirchhoff_2d(cfg,plan); times(j)=fast.wall_seconds; cpus(j)=fast.cpu_seconds;
end
fast.wall_seconds=median(times); fast.cpu_seconds=median(cpus);
x=cfg.surface_x_m; w=interp1(s.x_m,s.weights,x,'linear',0);
mask=x>=min(s.x_m(s.mask)) & x<=max(s.x_m(s.mask));
if n==4097, assert(isequal(cfg.surface_z_m,s.eta_m)); mask=s.mask; end
curv=interp1(s.x_m,s.second./(1+s.slope.^2).^1.5,x,'linear');
[inc,~,~]=local_source(cfg.spectral_source,x,zeros(size(x)));
iw=abs(inc).^2; iw=iw/sum(iw);
d=struct('frequency_hz',f,'n_surface',n,'full',ref.receiver_field,'fast',fast.receiver_field, ...
 'cost',cost,'Kmax_over_k',max(s.profile.k_rad_per_m)/plan.k_radpm, ...
 'curvature_rms_over_k',sqrt(mean(curv.^2))/plan.k_radpm, ...
 'illuminated_curvature_rms_over_k',sqrt(sum(iw.*curv.^2))/plan.k_radpm);
d=local_fast_record(d,fast,plan,s,w,mask);
sp=local_spectrum(ref.receiver_field,x,cfg.spectral_source,plan.k_radpm);
d.theta99_out_deg=sp.theta99; d.outside_incident99_energy=sp.outside;
d.spectrum=sp;
fprintf('  Full/Fast weighted %.7g whole %.7g cold speed %.3g\n',d.metrics.weighted,d.metrics.whole,d.cold_speedup);
end

function d=local_fast_record(d,f,p,s,w,mask)
if nargin<5, w=s.weights; mask=s.mask; end
d.fast=f.receiver_field; d.metrics=local_metric(d.fast,d.full,w,mask);
d.rank_in=f.rank; d.rank_out=f.outgoing_rank;
d.effective_rank_in=local_rank(p.incoming_singular_values); d.effective_rank_out=local_rank(p.outgoing_singular_values);
d.incoming_singular_values=p.incoming_singular_values; d.outgoing_singular_values=p.outgoing_singular_values;
d.fast_setup_wall=p.setup_wall_seconds; d.fast_setup_cpu=p.setup_cpu_seconds;
d.fast_apply_wall=f.wall_seconds; d.fast_apply_cpu=f.cpu_seconds;
d.warm_speedup=d.cost.apply_wall/f.wall_seconds;
d.cold_speedup=(d.cost.apply_wall+d.cost.setup_wall)/(p.setup_wall_seconds+f.wall_seconds);
d.plan_MiB=p.plan_bytes/2^20; d.setup_arrays_MiB=p.setup_named_arrays_bytes/2^20;
d.apply_arrays_MiB=f.named_apply_arrays_bytes/2^20;
d.full_block_estimated_MiB=256*d.n_surface*104/2^20;
[~,ic]=min(abs(p.x_m)); d.H_full_physical=d.full(ic); d.H_fast_physical=d.fast(ic);
k=p.k_radpm; d.H_full_reduced=d.H_full_physical*exp(-1i*k*103);
d.H_fast_reduced=d.H_fast_physical*exp(-1i*k*103);
src=s.source; src.kz=sqrt(k^2-src.kx.^2);
direct=local_source(src,0,3); direct_red=direct*exp(-1i*k*97);
a=apply_pe_channel_phase_reference_vertical(d.frequency_hz, ...
 struct('direct_f',direct_red,'reflect_fm',d.H_full_reduced),s.geometry,'direct_dsp');
b=apply_pe_channel_phase_reference_vertical(d.frequency_hz, ...
 struct('direct_f',direct_red,'reflect_fm',d.H_fast_reduced),s.geometry,'direct_dsp');
d.H_full_dsp=a.reflect_fm; d.H_fast_dsp=b.reflect_fm;
d.H_full_native_dsp=conj(d.H_full_physical)*exp(1i*k*97);
d.H_fast_native_dsp=conj(d.H_fast_physical)*exp(1i*k*97);
d.H_direct_native_dsp=conj(direct)*exp(1i*k*97);
adapted=apply_pe_channel_phase_reference_vertical(d.frequency_hz, ...
 struct('direct_f',conj(direct_red),'reflect_fm',conj(d.H_full_reduced)),s.geometry,'direct_dsp');
assert(abs(adapted.reflect_fm-d.H_full_native_dsp)<1e-12);
assert(abs(adapted.direct_f-d.H_direct_native_dsp)<1e-12);
d.H_direct_dsp=a.direct_f; d.H_total_full=NaN; d.H_total_fast=NaN;
assert(abs(a.total_fm-a.direct_f-a.reflect_fm)<1e-14);
assert(abs(d.H_full_reduced*exp(1i*k*103)-d.H_full_physical)<1e-14);
end

function r=local_rank(s)
tail=[sqrt(flipud(cumsum(flipud(abs(s(2:end)).^2)))/sum(abs(s).^2));0];
r=find(tail<=1e-3,1);
end

function m=local_metric(a,b,w,mask)
m.whole=norm(a-b)/norm(b); a=a(mask); b=b(mask); w=w(mask); w=w/sum(w);
m.weighted=sqrt(sum(w.*abs(a-b).^2)/sum(w.*abs(b).^2));
m.magnitude=sqrt(sum(w.*(abs(a)-abs(b)).^2)/sum(w.*abs(b).^2));
m.phase_rms=sqrt(sum(w.*angle(a.*conj(b)).^2));
end

function sp=local_spectrum(u,x,src,k)
n=numel(x); dx=mean(diff(x)); kx=2*pi*(-floor(n/2):ceil(n/2)-1)'/(n*dx);
win=.5-.5*cos(2*pi*(0:n-1)'/(n-1)); uk=dx*fftshift(fft(u.*win)); en=abs(uk).^2;
[ab,ii]=sort(abs(src.kx)); p=cumsum(abs(src.coeff(ii)).^2)/sum(abs(src.coeff).^2);
cut=ab(find(p>=.99,1)); prop=abs(kx)<=k;
[ab,ii]=sort(abs(kx(prop))); ep=en(prop); c=cumsum(ep(ii))/sum(ep);
sp=struct('kx',kx,'field_k',uk,'theta99',asind(ab(find(c>=.99,1))/k), ...
 'incident99_kx',cut,'outside',sum(en(prop & abs(kx)>cut))/sum(en(prop)), ...
 'scope','finite Hann-aperture receiver spectrum; not infinite-ocean flux');
end

function a=local_anchor(s,f)
d=s.sparse([s.sparse.frequency_hz]==f); idx=1:8:4097; x=s.x_m(idx);
cfg=local_cfg(s,f,4097,33,s.ranks); src=cfg.spectral_source;
bcfg=struct('frequency_hz',f,'c0_mps',1500,'half_width_m',60,'points_per_wavelength',8, ...
 'panel_order',8,'self_quadrature_order',96,'halfplane_h_m',-2,'window_plateau_ratio',.85, ...
 'eta_fn',@(xx)local_eta(s.profile,xx),'eta_prime_fn',@(xx)local_slope(s.profile,xx), ...
 'incident_fn',@(xx,zz)local_source(src,xx,zz),'receiver_x_m',x,'receiver_z_m',3);
clk=tic; cpu=cputime; b=solve_helmholtz_bie_halfplane_vertical(bcfg);
a=struct('frequency_hz',f,'field',b.receiver_field,'receiver_x_m',x, ...
 'metrics',local_metric(d.full(idx),b.receiver_field,s.weights(idx),s.mask(idx)), ...
 'boundary_residual',b.offgrid_boundary_residual,'linear_residual',b.linear_residual, ...
 'wall_s',toc(clk),'cpu_s',cputime-cpu,'nodes',b.node_count,'matrix_MiB',16*b.node_count^2/2^20);
fprintf('BIE_ANCHOR %g Hz weighted %.7g boundary %.4g\n',f,a.metrics.weighted,a.boundary_residual);
end

function g=local_guard(s,f)
d=s.sparse([s.sparse.frequency_hz]==f); cfg=local_cfg(s,f,4097,33,s.ranks);
p=prepare_spectral_fast_kirchhoff_2d(cfg); cfg.padding_factor=4;
a=solve_spectral_fast_kirchhoff_2d(cfg,p);
g=struct('frequency_hz',f,'padding_2_4',norm(a.receiver_field-d.fast)/norm(d.fast));
cfg.padding_factor=2; cfg.height_nodes=65; p=prepare_spectral_fast_kirchhoff_2d(cfg);
a=solve_spectral_fast_kirchhoff_2d(cfg,p);
g.interpolation_33_65=norm(a.receiver_field-d.fast)/norm(d.fast);
fine=local_case(s,f,8193,65,s.ranks); idx=1:2:8193;
g.full_grid=local_metric(d.full,fine.full(idx),s.weights,s.mask);
g.fast_grid=local_metric(d.fast,fine.fast(idx),s.weights,s.mask);
g.fine_full=fine.full; g.fine_fast=fine.fast;
[eta,~,~]=evaluate_pm_pilot_surface(s.profile,linspace(-60,60,8193)');
g.surface_closure=max(abs(eta(idx)-s.eta_m));
fprintf('GUARD %g Hz pad %.3g interp %.3g grid %.3g\n',f,g.padding_2_4,g.interpolation_33_65,g.full_grid.weighted);
end

function y=local_eta(p,x)
y=evaluate_pm_pilot_surface(p,x);
end
function y=local_slope(p,x)
[~,y]=evaluate_pm_pilot_surface(p,x);
end
function a=local_append(a,b)
if isempty(a), a=b; else, a(end+1)=orderfields(b,a); end
end
function yes=local_algorithm(d)
yes=~isempty(d);
for i=1:numel(d), yes=yes && d(i).metrics.whole<1e-3 && d(i).metrics.weighted<1e-3; end
end

function s=local_summary(s)
s.algorithm_pass=numel(s.sparse)==5 && local_algorithm(s.sparse);
s.physics_pass=numel(s.anchors)==3;
for a=s.anchors, s.physics_pass=s.physics_pass && a.metrics.weighted<1e-3 && a.boundary_residual<1e-4; end
s.numerical_pass=numel(s.guards)==3;
for g=s.guards
 s.numerical_pass=s.numerical_pass && g.padding_2_4<1e-12 && g.interpolation_33_65<1e-4 ...
  && g.full_grid.weighted<1e-4 && g.full_grid.whole<1e-4 && g.fast_grid.weighted<1e-4 ...
  && g.surface_closure==0;
end
s.speed_pass=numel(s.sparse)==5 && all([s.sparse.cold_speedup]>2);
s.sparse_phase_pass=numel(s.sparse)==5 && ...
 max(abs(angle([s.sparse.H_fast_dsp].*conj([s.sparse.H_full_dsp]))))<1e-3;
s.sparse_pass=s.algorithm_pass && s.physics_pass && s.numerical_pass && s.speed_pass && s.sparse_phase_pass;
s.pass=false; s.phase=struct();
s.reference_control=local_reference_control();
s.reference_bridge_pass=s.reference_control.adapter_complex_l2<1e-12;
s.consistency_pass=false;
if numel(s.dense)==65
 s.phase=local_phase(s.dense);
 nd=s.dense;
 for i=1:65, nd(i).H_full_dsp=nd(i).H_full_native_dsp; nd(i).H_fast_dsp=nd(i).H_fast_native_dsp; end
 s.native_phase=local_phase(nd);
 s.consistency_pass=s.sparse_pass && local_algorithm(s.dense) && s.phase.relative_phase_max<1e-3 ...
  && s.phase.group_delay_difference_max_s<1e-6 && s.phase.cir_l2<1e-3 ...
  && s.phase.reduced_adjacent_phase_max<pi/2 && s.native_phase.cir_l2<1e-3 ...
  && s.native_phase.relative_phase_max<1e-3 && s.native_phase.group_delay_difference_max_s<1e-6 ...
  && all([s.dense.cold_speedup]>2);
 s.pass=s.consistency_pass && s.reference_bridge_pass;
end
if s.consistency_pass
 for name={'sparse','dense'}
  for i=1:numel(s.(name{1}))
   s.(name{1})(i).H_total_full=s.(name{1})(i).H_direct_native_dsp+s.(name{1})(i).H_full_native_dsp;
   s.(name{1})(i).H_total_fast=s.(name{1})(i).H_direct_native_dsp+s.(name{1})(i).H_fast_native_dsp;
  end
 end
end
s.status='INCOMPLETE';
 if s.consistency_pass
 s.status='PASS';
elseif numel(s.sparse)==5 && (~s.algorithm_pass || (numel(s.anchors)==3 && ~s.physics_pass) ...
  || (numel(s.guards)==3 && ~s.numerical_pass))
 s.status='FAIL';
elseif numel(s.dense)==65
 s.status='FAIL';
end
if s.consistency_pass && ~s.reference_bridge_pass, s.status='CONSISTENCY_PASS_REFERENCE_BRIDGE_FAIL'; end
s.phase_reference_meta=struct('source_representation','native exp(-i*omega*t) physical phasor', ...
 'surface_chart','frozen Full/BIE/Fast physical z=eta(x); water-facing normal', ...
 'native_reduced_note','P*exp(-i*omega*tau_reflect0); not assumed equal to public PE reduced fields', ...
 'primary_transfer','H_native_dsp=conj(P)*exp(+i*omega*tau_direct0)', ...
 'primary_H_csv','dense_native_physical_H.csv','primary_CIR_csv','reflected_cir.csv', ...
 'time_origin_s',97/1500,'source_to_helper_bridge','conj(native reduced) supplied only to direct_dsp target', ...
 'public_absolute_physical_helper_note','not used to synthesize native FK physical fields from DSP-adapted envelopes', ...
 'candidate_H_dsp_note','H_full_dsp/H_fast_dsp are the unadapted carrier-only import diagnostic', ...
 'total_reference','H_total_full/H_total_fast sum native physical direct-DSP components after gates', ...
 'production_integration',false);
end

function c=local_reference_control()
% Analytic fixed-kx plane wave: exact physical path phase, no solver fitting.
f=(4000:62.5:8000)'; k=2*pi*f/1500; kx=4; kz=sqrt(k.^2-kx^2);
physical=-exp(1i*kz*103); red=physical.*exp(-1i*k*103);
g=struct('z_tx',100,'z_rx',3,'c0',1500,'z_surface',0);
a=apply_pe_channel_phase_reference_vertical(f,struct('direct_f',ones(65,1),'reflect_fm',red),g,'direct_dsp');
native=conj(physical).*exp(1i*k*97);
project_delay=-gradient(unwrap(angle(a.reflect_fm)),f)/(2*pi);
native_delay=-gradient(unwrap(angle(native)),f)/(2*pi);
theory_native=103/1500*k./kz-97/1500;
theory_project=6/1500-103/1500*(k./kz-1);
c=struct('f_hz',f,'kx_radpm',kx,'physical',physical,'project_dsp',a.reflect_fm, ...
 'native_dsp',native,'project_delay_s',project_delay,'native_delay_s',native_delay, ...
 'theory_native_s',theory_native,'theory_project_s',theory_project, ...
 'delay_difference_max_s',max(abs(native_delay-project_delay)), ...
 'physical_to_ifft','H_native_dsp=conj(P_native)*exp(+i*omega*tau_direct0)', ...
 'meaning','Global Fourier-sign conversion; no phase fitting or peak alignment');
adapted=apply_pe_channel_phase_reference_vertical(f, ...
 struct('direct_f',ones(65,1),'reflect_fm',conj(red)),g,'direct_dsp');
c.adapted_dsp=adapted.reflect_fm;
c.adapter_complex_l2=norm(adapted.reflect_fm-native)/norm(native);
c.adapted_delay_s=-gradient(unwrap(angle(adapted.reflect_fm)),f)/(2*pi);
end

function p=local_phase(d)
p=audit_fk_multifrequency_phase(d);
end

function h=local_hashes(root)
paths={'scripts/validation/support/solve_full_kirchhoff_2d.m', ...
 'scripts/validation/support/solve_helmholtz_bie_halfplane_vertical.m', ...
 'scripts/validation/support/prepare_spectral_fast_kirchhoff_2d.m', ...
 'scripts/validation/support/solve_spectral_fast_kirchhoff_2d.m', ...
 'scripts/validation/support/evaluate_pm_pilot_surface.m', ...
 'src/channel/apply_pe_channel_phase_reference_vertical.m', ...
 'results/validation/spectral_fast_kirchhoff_pm_pilot/seed_260001_fields.mat', ...
 'results/validation/bellhop_internal_pm_fixed_realization/fixed_pm_fourier_coefficients.csv'};
h=struct([]);
for i=1:numel(paths)
 fid=fopen(fullfile(root,paths{i}),'rb'); assert(fid>=0); bytes=fread(fid,Inf,'*uint8'); fclose(fid);
 md=java.security.MessageDigest.getInstance('SHA-256'); md.update(bytes);
 digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
 h=local_append(h,struct('path',paths{i},'sha256',digest));
end
end

function local_outputs(s,out,root,render)
c=s.reference_control;
writetable(table(c.f_hz,c.project_delay_s,c.native_delay_s,c.adapted_delay_s,c.theory_project_s,c.theory_native_s, ...
 'VariableNames',{'f_hz','unadapted_FD_delay_s','native_FD_delay_s','adapted_FD_delay_s','unadapted_exact_delay_s','native_exact_delay_s'}), ...
 fullfile(out,'phase_reference_plane_wave_control.csv'));
for family={'frozen_rank_sparse','sparse','dense'}
 if ~isfield(s,family{1}), continue; end
 d=s.(family{1}); rows=struct([]);
 for a=d
  row=rmfield(a,{'full','fast','cost','metrics','spectrum','incoming_singular_values','outgoing_singular_values'});
  row=rmfield(row,{'H_full_physical','H_fast_physical','H_full_reduced','H_fast_reduced','H_full_dsp', ...
   'H_fast_dsp','H_direct_dsp','H_total_full','H_total_fast','H_full_native_dsp','H_fast_native_dsp','H_direct_native_dsp'});
  row.whole=a.metrics.whole; row.weighted=a.metrics.weighted; row.magnitude=a.metrics.magnitude;
  row.phase_rms=a.metrics.phase_rms; row.full_setup_wall=a.cost.setup_wall; row.full_apply_wall=a.cost.apply_wall;
  row.full_setup_cpu=a.cost.setup_cpu; row.full_apply_cpu=a.cost.apply_cpu; rows=local_append(rows,row);
 end
 if ~isempty(rows), writetable(struct2table(rows),fullfile(out,[family{1} '_metrics.csv'])); end
 if ~isempty(d)
  f=vertcat(d.frequency_hz); a=vertcat(d.H_full_physical); b=vertcat(d.H_fast_physical);
  c=vertcat(d.H_full_dsp); e=vertcat(d.H_fast_dsp);
  writetable(table(f,real(a),imag(a),real(b),imag(b),real(c),imag(c),real(e),imag(e), ...
   angle(a),angle(b),unwrap(angle(a)),unwrap(angle(b)),angle(e.*conj(c)), ...
   'VariableNames',{'f_hz','Full_phys_real','Full_phys_imag','Fast_phys_real','Fast_phys_imag', ...
   'Full_dsp_real','Full_dsp_imag','Fast_dsp_real','Fast_dsp_imag','Full_phase','Fast_phase', ...
   'Full_blind_unwrap','Fast_blind_unwrap','relative_phase'}),fullfile(out,[family{1} '_H.csv']));
  a=vertcat(d.H_full_native_dsp); b=vertcat(d.H_fast_native_dsp);
  writetable(table(f,real(a),imag(a),real(b),imag(b),angle(a),angle(b), ...
   'VariableNames',{'f_hz','Full_native_dsp_real','Full_native_dsp_imag','Fast_native_dsp_real', ...
   'Fast_native_dsp_imag','Full_native_dsp_phase','Fast_native_dsp_phase'}), ...
   fullfile(out,[family{1} '_native_physical_H.csv']));
 end
end
rows=struct([]);
for a=s.anchors
 row=rmfield(a,{'field','receiver_x_m','metrics'}); row.weighted=a.metrics.weighted;
 row.whole=a.metrics.whole; row.magnitude=a.metrics.magnitude; row.phase_rms=a.metrics.phase_rms;
 rows=local_append(rows,row);
end
if ~isempty(rows), writetable(struct2table(rows),fullfile(out,'bie_anchors.csv')); end
rows=struct([]);
for g=s.guards
 row=rmfield(g,{'fine_full','fine_fast','full_grid','fast_grid'});
 row.full_grid_weighted=g.full_grid.weighted; row.full_grid_whole=g.full_grid.whole;
 row.fast_grid_weighted=g.fast_grid.weighted; rows=local_append(rows,row);
end
if ~isempty(rows), writetable(struct2table(rows),fullfile(out,'numerical_guards.csv')); end
if ~isempty(s.rank_trials), writetable(struct2table(s.rank_trials),fullfile(out,'uniform_rank_trials.csv')); end
fid=fopen(fullfile(out,'provenance.json'),'w'); fprintf(fid,'%s',jsonencode(s.reference_hashes,PrettyPrint=true)); fclose(fid);
fid=fopen(fullfile(out,'phase_reference_meta.json'),'w'); fprintf(fid,'%s',jsonencode(s.phase_reference_meta,PrettyPrint=true)); fclose(fid);
if render, local_render(s,out,root); end
end

function local_render(s,out,root)
figdir=fullfile(out,'figures'); if ~isfolder(figdir), mkdir(figdir); end
d=s.sparse; f=[d.frequency_hz]';
fig=figure('Visible','off','Position',[40 40 1200 850]); tiledlayout(2,2);
nexttile; m=[d.metrics]; semilogy(f,[m.weighted],'o-',f,[m.whole],'s--'); yline(1e-3,':');
xlabel('frequency (Hz)'); ylabel('Fast--Full complex L2'); legend('fixed M99 weighted','whole grid','gate'); grid on;
if isfield(s,'frozen_rank_sparse')
 hold on; fm=[s.frozen_rank_sparse.metrics]; semilogy(f,[fm.weighted],'v:');
 legend('selected rank weighted','selected rank whole','gate','frozen3/8 weighted');
end
nexttile; plot(f,angle([d.H_full_physical]),'o-',f,angle([d.H_fast_physical]),'x--');
xlabel('frequency (Hz)'); ylabel('wrapped physical phase (rad)'); legend('Full','Fast'); grid on;
nexttile; plot(f,[d.effective_rank_in],'o-',f,[d.effective_rank_out],'s--');
xlabel('frequency (Hz)'); ylabel('effective rank at 1e-3'); legend('incoming','outgoing'); grid on;
nexttile; plot(f,[d.cold_speedup],'o-',f,[d.warm_speedup],'s--');
xlabel('frequency (Hz)'); ylabel('speedup'); legend('cold','warm'); grid on;
exportgraphics(fig,fullfile(figdir,'sparse_audit.png'),'Resolution',160); close(fig);
if numel(s.dense)==65
 p=s.native_phase; f=p.f_hz;
 writetable(table(f,p.full_physical_unwrapped,p.fast_physical_unwrapped,p.full_dsp_unwrapped, ...
  p.fast_dsp_unwrapped,p.full_physical_group_delay_s,p.fast_physical_group_delay_s, ...
  p.full_dsp_group_delay_s,p.fast_dsp_group_delay_s,p.relative_phase_rad),fullfile(out,'phase_group_delay.csv'));
 writetable(table(p.delay_s,real(p.cir_full),imag(p.cir_full),real(p.cir_fast),imag(p.cir_fast)),fullfile(out,'reflected_cir.csv'));
 fig=figure('Visible','off','Position',[40 40 1200 1000]); tiledlayout(3,2);
 nexttile; plot(f,p.full_physical_unwrapped,'-',f,p.fast_physical_unwrapped,'--'); ylabel('carrier-restored physical phase (rad)');
 nexttile; plot(f,p.full_dsp_unwrapped,'-',f,p.fast_dsp_unwrapped,'--'); ylabel('direct-DSP unwrapped phase (rad)');
 nexttile; plot(f,1e3*p.full_dsp_group_delay_s,'-',f,1e3*p.fast_dsp_group_delay_s,'--'); ylabel('direct-DSP group delay (ms)');
 nexttile; plot(f,p.relative_phase_rad,'o-'); ylabel('Fast--Full phase (rad)');
 nexttile; plot(f,abs([s.dense.H_full_native_dsp]),'-',f,abs([s.dense.H_fast_native_dsp]),'--'); ylabel('|H reflect|');
 nexttile; plot(1e3*p.delay_s,abs(p.cir_full),'-',1e3*p.delay_s,abs(p.cir_fast),'--'); xlabel('relative delay (ms)'); ylabel('|h|, diagnostic');
 for i=1:5, ax=nexttile(i); xlabel(ax,'frequency (Hz)'); grid(ax,'on'); end
 for i=[1 2 3 5 6], ax=nexttile(i); legend(ax,'Full FK','Fast FK'); grid(ax,'on'); end
 exportgraphics(fig,fullfile(figdir,'dense_phase_cir.png'),'Resolution',160); close(fig);
 p=s.native_phase;
 writetable(table(p.delay_s,real(p.cir_full),imag(p.cir_full),real(p.cir_fast),imag(p.cir_fast)), ...
  fullfile(out,'native_physical_reflected_cir.csv'));
 writetable(table(p.f_hz,p.full_dsp_unwrapped,p.fast_dsp_unwrapped,p.full_dsp_group_delay_s, ...
  p.fast_dsp_group_delay_s),fullfile(out,'native_physical_phase_group_delay.csv'));
 fig=figure('Visible','off','Position',[40 40 1100 700]); tiledlayout(2,2);
 nexttile; plot(f,1e3*s.phase.full_dsp_group_delay_s,'-',f,1e3*p.full_dsp_group_delay_s,'--');
 xlabel('frequency (Hz)'); ylabel('group delay (ms)'); legend('unadapted carrier-only import','native physical to IFFT'); grid on;
 nexttile; plot(1e3*p.delay_s,abs(p.cir_full),'-',1e3*p.delay_s,abs(p.cir_fast),'--');
 xlabel('relative delay (ms)'); ylabel('|h|, native physical diagnostic'); legend('Full','Fast'); grid on;
 c=s.reference_control;
 nexttile; plot(c.f_hz,1e3*c.project_delay_s,'-',c.f_hz,1e3*c.native_delay_s,'--',c.f_hz,1e3*c.adapted_delay_s,'x:');
 xlabel('frequency (Hz)'); ylabel('analytic plane-wave delay (ms)'); legend('unadapted import','physical IFFT mapping','adapted helper'); grid on;
 nexttile; plot(c.f_hz,1e3*c.theory_project_s,'-',c.f_hz,1e3*c.theory_native_s,'--');
 xlabel('frequency (Hz)'); ylabel('exact derivative (ms)'); legend('unadapted import formula','physical path formula'); grid on;
 exportgraphics(fig,fullfile(figdir,'native_ifft_reference_control.png'),'Resolution',160); close(fig);
end
fid=fopen(fullfile(root,'reports','spectral_fast_kirchhoff_multifrequency_consistency_report.md'),'w','n','UTF-8');
cl=onCleanup(@()fclose(fid));
fprintf(fid,'# Spectral Fast FK multi-frequency consistency\n\nStatus: **%s**. Date: %s. MATLAB %s.\n\n',s.status,s.date,s.matlab_version);
if s.pass
 fprintf(fid,'本轮 multi-frequency consistency **PASS**：固定 seed260001/C2 海面在4–8 kHz的65个频点上，统一3/12 rank满足FastFK–FullFK复误差门槛，三个BIE物理锚点、数值护栏、相位/群时延/CIR一致性及加速门均通过。冻结3/8在6–8 kHz失败，3/10在8 kHz仍略超标，因此采用整个频带统一3/12，未做逐频tuning。\n\n相位桥接不能只补载波。原生Full/BIE/Fast物理相量到MATLAB IFFT传递函数采用一次明确的全包络Fourier符号转换：conj(P)*exp(+iωτdir0)，三条路径共享同一个直达时刻参考。验证层经该数学转换后再使用现有direct-DSP helper；这不是逐频相位拟合，也没有改动reference或生产API。未转换的carrier-only导入仍保留为失败的独立控制。后续固定PM宽带及小规模realization数据可以沿这个明确记录的表示继续；本轮不把原生FK reduced字段直接当作公共PE reduced字段或宣称生产集成已完成。\n\n');
end
if s.consistency_pass && ~s.reference_bridge_pass
 fprintf(fid,'本轮的 Full/Fast FK 跨频算法、物理锚点、数值稳定性、相位/群时延与 CIR 一致性检查通过，但不能给出无条件 multi-frequency consistency PASS。原因是原生 exp(-iωt) 物理场到项目 direct-DSP 的桥接并非仅补载波：固定非零横向波数控制揭示 reduced envelope 的相位符号也必须明确转换。\n\n统一 rank 已从冻结 3/8 调整为整频带 3/12，未做逐频 tuning。原生物理对应的 IFFT 传递函数和项目当前 carrier-only 表示分别保留，Full/Fast 在两种表示下均一致；只有前者可直接解释为该原生物理场的延迟响应。本轮没有改动任何冻结 reference 或生产接口。正式使用项目统一相位 schema 的宽带 PM 数据生成暂不批准，先解决表示桥接；这不是发现新的散射理论问题。\n\n');
end
fprintf(fid,'Canonical seed260001; 12 immutable Fourier modes; actual Kmax %.9g rad/m; coefficient SHA256 `%s`. Surface x=[-60,60] m, N4097; C2 plateau42/support50 m, common height bound %.3g m. All frequencies use bit-identical eta, analytic slopes and coefficients, no height scaling, complex scalar fitting or phase alignment. Source: frozen984-point Gaussian sigma=.3 m, period192.1875 m, x0=-96.09375 m; all saved bins retained, frequency-dependent exact kz. Uniform c1500, Tx100/Rx3, 2-D line-source validation chain; production PE is unchanged.\n\n', ...
 max(s.profile.k_rad_per_m),s.profile.coeff_file_sha256,s.height_bound_m);
fprintf(fid,'Sparse axis4/5/6/7/8 kHz; BIE anchors4/6/8 kHz; selected common rank%d/%d, padding2. Dense follow-up65 samples, delta-f62.5 Hz, window16 ms, bin spacing16/65=.246154 ms; bandwidth4 kHz gives nominal resolution about.25 ms. Sparse delta-f1 kHz has only1 ms window and aliases nominal4 ms relative delay: no sparse CIR or physical group-delay certification. Dense execution requires sparse algorithm, physics, numerical and cold speed gates.\n\n',s.ranks);
fprintf(fid,'The16 ms window contains the nominal4 ms reflection excess delay, but does not certify all scattering tails as unaliased. A point-source geometric path through x=60 m at the flat integration edge permits excess delay(sqrt(100^2+60^2)+sqrt(3^2+60^2)-97)/1500=53.129 ms. Such contributions may be weak/nonstationary; their absence has not been proved. CIR remains a finite-band periodic diagnostic, not a support-converged final response. Samples are IFFT-ordered as frequency offsets from4 kHz (no ifftshift); a different baseband carrier changes complex modulation, not the compared magnitude/timing.\n\n');
fprintf(fid,'## Phase reference\n\nPhysical exp(-i omega t) Full/Fast outputs carry exp(+ik*103). Native reduced reflection is P*exp(-ik*103); it is not automatically treated as the public PE reduced field. The common-origin MATLAB IFFT transfer is **H_native_dsp=conj(P)*exp(+ik*97)**. Equivalently, supply conj(native reduced) to the unchanged `apply_pe_channel_phase_reference_vertical` helper for its direct_dsp target, which multiplies by exp(-ik*6). The direct source uses the same saved Gaussian coefficients propagated97 m; its DSP envelope is conj(Pdirect*exp(-ik*97)). This is one deterministic Fourier-sign conversion, not a fitted complex scalar, per-frequency phase zero or peak alignment. Native physical fields remain untouched and are not reconstructed by calling the helper absolute_physical target on a DSP-adapted envelope. Total is formed only after reflected-only consistency gates. Public PE/interface semantics are unchanged.\n\nRaw physical phase is saved; blind unwrap is aliased even at62.5 Hz because nominal absolute flight time68.667 ms exceeds16 ms. Carrier-restored physical unwrap is unwrap(native reduced phase)+2*pi*f*103/1500. Its absolute branch is only defined modulo2*pi. Negative -dphi/(2*pi*df) for a positive physical carrier is expected; the IFFT-ready response has the opposite Fourier sign.\n\n');
fprintf(fid,'### Independent native-physical to IFFT bridge control\n\nThe unadapted candidate exp(-ik*6)*(P*exp(-ik*103)) is retained only as a carrier-only import diagnostic (H_full_dsp/H_fast_dsp in MAT and dense_H.csv). It is not the primary physical transfer. A fixed-kx4 rad/m plane wave P=-exp(i*kz*103) gives native relative delay=(103*k/kz-97)/1500; the unadapted import gives6/1500-103*(k/kz-1)/1500. They agree at kx0 and differ for kx nonzero: maximum sampled delay difference %.7g ms. After the whole-envelope conversion, helper output matches the independently constructed native physical IFFT transfer at relative complex L2 %.7g; bridge pass=%d. This validates the external native FK import, not a production API migration. Primary tables are dense_native_physical_H.csv, phase_group_delay.csv and reflected_cir.csv; phase_reference_meta documents the source/target representation and single time origin.\n\n', ...
 1e3*s.reference_control.delay_difference_max_s,s.reference_control.adapter_complex_l2,s.reference_bridge_pass);
fprintf(fid,'## Per-frequency algorithm fidelity and descriptors\n\n|Hz|weighted L2|whole L2|magnitude|phase RMS rad|effective in/out|Kmax/k|curvature RMS/k|illuminated curvature/k|out theta99 deg|outside incident99|\n|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|\n');
for a=s.sparse
 fprintf(fid,'|%g|%.8g|%.8g|%.7g|%.7g|%d/%d|%.7g|%.7g|%.7g|%.7g|%.7g|\n',a.frequency_hz,a.metrics.weighted,a.metrics.whole,a.metrics.magnitude,a.metrics.phase_rms,a.effective_rank_in,a.effective_rank_out,a.Kmax_over_k,a.curvature_rms_over_k,a.illuminated_curvature_rms_over_k,a.theta99_out_deg,a.outside_incident99_energy);
end
fprintf(fid,'\nFixed canonical Stage0 M99 weights/mask are held across frequencies; whole-grid error is separate. All-native complex fields are retained. Curvature is eta_second/(1+slope^2)^(3/2), illumination weights are incident reference-plane energy at each frequency. Spectrum is a finite Hann-aperture diagnostic, not far-field ocean flux; no new universal validity law.\n\n## Physical anchors\n\n|Hz|Full--BIE weighted|whole line diagnostic|phase RMS|boundary residual|nodes|wall s|matrix MiB|\n|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for a=s.anchors, fprintf(fid,'|%g|%.8g|%.7g|%.7g|%.7g|%d|%.5g|%.5g|\n',a.frequency_hz,a.metrics.weighted,a.metrics.whole,a.metrics.phase_rms,a.boundary_residual,a.nodes,a.wall_s,a.matrix_MiB); end
fprintf(fid,'\nBIE uses L60, ppw8/order8/self96, h=-2, plateau.85, the same analytical C2 surface and incident angular spectrum. Physical gate is weighted error<1e-3 and offgrid residual<1e-4, restricted to the frozen M99 core; edge errors are diagnostic. Numerical BIE refinement is required if a residual or near-threshold ambiguity appears.\n\n## Numerical guards\n\n|Hz|padding2/4|height nodes33/65|Full4097/8193 weighted|whole|Fast grid weighted|eta closure|\n|---:|---:|---:|---:|---:|---:|---:|\n');
for g=s.guards, fprintf(fid,'|%g|%.7g|%.7g|%.7g|%.7g|%.7g|%.3g|\n',g.frequency_hz,g.padding_2_4,g.interpolation_33_65,g.full_grid.weighted,g.full_grid.whole,g.fast_grid.weighted,g.surface_closure); end
fprintf(fid,'\nExecution note: the8 kHz BIE exceeded the MCP300-second response wait; the original shared-session computation continued and saved its checkpoint after about481 seconds. A subsequent harmless MCP probe confirmed recovery. No duplicate BIE run or shell MATLAB launch was used.\n\n');
fprintf(fid,'\n## CPU and memory\n\n|Hz|Full setup/apply wall s|Fast setup/apply wall s|warm/cold speed|plan MiB|setup/apply arrays MiB|Full block estimate MiB|\n|---:|---|---|---|---:|---|---:|\n');
for a=s.sparse, fprintf(fid,'|%g|%.5g / %.5g|%.5g / %.5g|%.4g / %.4g|%.5g|%.5g / %.5g|%.5g|\n',a.frequency_hz,a.cost.setup_wall,a.cost.apply_wall,a.fast_setup_wall,a.fast_apply_wall,a.warm_speedup,a.cold_speedup,a.plan_MiB,a.setup_arrays_MiB,a.apply_arrays_MiB,a.full_block_estimated_MiB); end
fprintf(fid,'\nFast apply median of3 after warmup (selected sparse rank remeasured with timing mode); Full one apply per frequency, timing is hardware/session-dependent. CPU columns are in metrics CSV. Memory is named arrays/analytical block estimate, not measured process peak. Geometry eta/slopes/normal/quadrature and Gaussian coefficients are per-realization reusable; kz, both SVD bases, height carrier, Hankel lag kernels and propagation phases are per-frequency. Existing frozen implementation rebuilds them; no claimed cross-frequency setup acceleration. Geometry reuse is already represented by fixed surface state; its negligible load/setup time is not included in cold speed.\n\n');
if isfield(s,'geometry_setup_wall')
 fprintf(fid,'Measured one-time fixed Fourier geometry evaluation/normal/quadrature setup: %.7g wall s, %.7g CPU s, named geometry arrays %.7g MiB. Saved-source/input file IO is outside this kernel benchmark.\n\n',s.geometry_setup_wall,s.geometry_setup_cpu,s.geometry_arrays_MiB);
end
if numel(s.dense)==65
 p=s.native_phase; fprintf(fid,'## Dense phase and CIR diagnostic\n\nAll65-point algorithm gate: %d. Max native relative phase %.8g rad; max DSP group-delay difference %.8g s; CIR relative complex L2 %.8g; max adjacent reduced phase %.7g rad. Sparse/dense boundary is governed by fixed gates, never total-channel cancellation. IFFT uses the rectangular65-sample native-physical direct-DSP reflected-only response with no taper, zero padding or peak alignment; diagnostic only, not a converged final CIR.\n\n',local_algorithm(s.dense),p.relative_phase_max,p.group_delay_difference_max_s,p.cir_l2,p.reduced_adjacent_phase_max);
 m=[s.dense.metrics]; [~,ip]=max(abs(p.cir_full));
 fprintf(fid,'Across65 samples: max weighted/whole-grid complex errors %.9g / %.9g; minimum cold/warm speedups %.7g / %.7g. Native direct-DSP group delay range %.9g--%.9g ms; sampled CIR maximum at %.9g ms. The grid bin size explains why the sampled peak differs from the approximately3.595 ms phase delay; no peak alignment is performed.\n\n', ...
  max([m.weighted]),max([m.whole]),min([s.dense.cold_speedup]),min([s.dense.warm_speedup]), ...
  1000*min(p.full_dsp_group_delay_s),1000*max(p.full_dsp_group_delay_s),1000*p.delay_s(ip));
 fprintf(fid,'The unadapted carrier-only diagnostic also preserves Full/Fast fidelity, but that alone does not establish its native physical delay interpretation. All primary H/CIR results use the validated native-to-DSP bridge and one shared direct time origin.\n\n');
else
 fprintf(fid,'## Dense phase and CIR diagnostic\n\nNot executed: dense axis is conditional on all sparse gates. No final CIR or broadband PASS can be inferred from the sparse response.\n\n');
end
if isfield(s,'frozen_rank_sparse')
 fprintf(fid,'\n## Frozen-rank failure and common-rank decision\n\n|Hz|3/8 weighted error|3/8 whole error|\n|---:|---:|---:|\n');
 for a=s.frozen_rank_sparse, fprintf(fid,'|%g|%.8g|%.8g|\n',a.frequency_hz,a.metrics.weighted,a.metrics.whole); end
 fprintf(fid,'\nThe frozen3/8 fails from6 kHz. Common3/10 narrowly fails at8 kHz; incoming4/5 at outgoing10 does not resolve it. The next outgoing ladder3/12 passes all sparse frequencies. All trials are saved in uniform_rank_trials.csv; no frequency-specific tuning. Chosen rank%d/%d applies to the complete sparse/dense band.\n\n',s.ranks);
end
fprintf(fid,'## Decision and evidence boundary\n\nAlgorithm=%d; physics=%d; numerical=%d; speed=%d; sparse=%d; multi-frequency PASS=%d. ',s.algorithm_pass,s.physics_pass,s.numerical_pass,s.speed_pass,s.sparse_pass,s.pass);
fprintf(fid,'Full/Fast frequency consistency=%d; native-to-project phase bridge=%d. ',s.consistency_pass,s.reference_bridge_pass);
if s.pass, fprintf(fid,'Allows the next separately scoped fixed-PM broadband and small realization data stage under these tested conditions only.\n\n');
else, fprintf(fid,'Does not authorize formal broadband PM data generation. Investigate failed or incomplete gates first; no new scattering theory or production integration performed.\n\n'); end
fprintf(fid,'Verification: five independent tests cover nominal4 ms delay, sparse aliasing, fixed surface nesting, nonzero-kx native/DSP conversion, and the actual phase assembler with a known4.3 ms response. Artifact assertions recompute complex errors from raw fields, compare CSV imaginary signs to scalar native fields, verify exp(i*restored physical phase), check the physical/DSP derivative sign identity, and confirm total-component closure. Code Analyzer is clean for the four new MATLAB files. Complex frequency collection uses vertcat (no implicit conjugate transpose). Reference-phase algebra checks use1e-12 to cover alternative evaluation orders of large carrier arguments; mathematical/scientific fidelity gates remain1e-3 and the delay-difference gate1 microsecond. No production or communication regression run was needed because this task only adds validation workflows and explicitly excludes communication work.\n\n');
fprintf(fid,'Artifacts: `results/validation/spectral_fast_kirchhoff_multifrequency_consistency/consistency.mat`, CSV tables, raw complex fields, fixed surface, input/reference SHA256 provenance and figures. Reference FullFK/BIE/FastFK implementations and public paramsV/output remain unchanged.\n\n![Sparse audit](../results/validation/spectral_fast_kirchhoff_multifrequency_consistency/figures/sparse_audit.png)\n');
if numel(s.dense)==65, fprintf(fid,'\n![Dense phase and CIR](../results/validation/spectral_fast_kirchhoff_multifrequency_consistency/figures/dense_phase_cir.png)\n'); end
if numel(s.dense)==65, fprintf(fid,'\n![Native physical IFFT and phase bridge control](../results/validation/spectral_fast_kirchhoff_multifrequency_consistency/figures/native_ifft_reference_control.png)\n'); end
end
