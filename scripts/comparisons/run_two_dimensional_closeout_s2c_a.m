function out=run_two_dimensional_closeout_s2c_a(root,run_id)
% Saved-field diagnostics and analytic geometry ONLY; no acoustic evaluation.
t0=tic;
base=fullfile(root,'results/comparisons/two_dimensional_channel_closeout');
old=fullfile(base,'20261023_s2b_1305/formal_reference_v2');
out=fullfile(base,run_id); assert(~isfolder(out),'Use a new run_id');
mkdir(out); mkdir(fullfile(out,'mechanism_baseline')); mkdir(fullfile(out,'protocol')); mkdir(fullfile(out,'provenance'));
data=fullfile(out,'mechanism_baseline'); protocol=fullfile(out,'protocol');
S=jsondecode(fileread(fullfile(old,'summary.json')));
assert(strcmp(S.status,'S2B_COMPLETE_REFERENCE_QUALIFIED'));
H=readtable(fullfile(old,'provenance/final_output_hashes.csv'),TextType='string');
for j=1:height(H), assert(strcmp(closeout_s2b_hash(H.path(j)),H.sha256(j))); end
writetable(H,fullfile(out,'provenance/S2b_preserved_hashes.csv'));
C=readtable(fullfile(old,'provenance/code_hashes.csv'),TextType='string');
for j=1:height(C), assert(strcmp(closeout_s2b_hash(C.path(j)),C.sha256(j))); end
writetable(C,fullfile(out,'provenance/S2b_code_hashes.csv'));
D=load(fullfile(old,'reference_data.mat')); P=load(fullfile(old,'plot_data.mat'));
L=readtable(fullfile(old,'fields_line.csv'),TextType='string');
N=readtable(fullfile(old,'fields_core.csv'),TextType='string');
assert(height(L)==4932 && height(N)==108 && all(L.validity));
assert(all(L.frequency_hz==4000) && isscalar(unique(L.source_id)));
assert(all(L.field_version(L.method=="FullFK")=="FullFK_complete_analytic_flat_exterior_GL8_v1"));
assert(all(L.representation=="native exp(-i omega t); outgoing exp(+ikR)"));
copyfile(fullfile(old,'fields_line.csv'),fullfile(data,'fields_line.csv'));
copyfile(fullfile(old,'fields_core.csv'),fullfile(data,'fields_core.csv'));
copyfile(fullfile(old,'reference.csv'),fullfile(data,'reference.csv'));
copyfile(fullfile(old,'resolutions.csv'),fullfile(data,'resolutions.csv'));
names=["flat","weak_low_K","region_II_low_K","strong_height_low_K"];
methods=["PE","BIE","FullFK"]; amps=[0 .01 .05 .20];
x=D.x(:); z=D.z(:); U=D.U; grid=D.grid; ix9=D.ix9;
assert(isequal(size(U),[417 4 3]) && all(z(1:411)==3));
assert(isequal(grid.receiver_id(grid.z_m==3),["P21";"P22";"P23"]));
for c=1:4
 for m=1:3
  a=L(L.case_name==names(c)&L.method==methods(m),:);
  assert(max(abs(complex(a.u_real,a.u_imag)-U(1:411,c,m)))<1e-14);
 end
end
zz=find(abs(P.plot_data.strong.z_m-3)<1e-12); assert(isscalar(zz));
xzcheck=[max(abs(P.plot_data.strong.PE(:,zz)-U(1:411,4,1))),max(abs(P.plot_data.strong.BIE(:,zz)-U(1:411,4,2)))];
assert(max(xzcheck)<1e-12);
src=D.src; save(fullfile(protocol,'frozen_source.mat'),'src');
writetable(grid,fullfile(protocol,'receivers_nine.csv'));
writetable(table(x(1:411),z(1:411),ones(411,1)/411,VariableNames={'x_m','z_m','q'}),fullfile(protocol,'receivers_line.csv'));
% Analytic geometry is evaluated on geometry-only supports, never new fields.
gx=linspace(-60,60,12001)'; geometry=table; geometry_line=table; illumination=table;
for c=1:4
 [eta,slope,curv]=geom(gx,amps(c),.1,0);
 geometry=[geometry;table(repmat(names(c),numel(gx),1),gx,eta,slope,curv,VariableNames={'case_name','x_m','eta_m','slope','curvature_per_m'})]; %#ok<AGROW>
 [eta,slope,curv]=geom(x,amps(c),.1,0);
 geometry_line=[geometry_line;table(repmat(names(c),417,1),(1:417)',x,z,eta,slope,curv,VariableNames={'case_name','sample_index','x_m','z_m','eta_m','slope','curvature_per_m'})]; %#ok<AGROW>
 f=D.fk_runs{c,2}; sx=f.surface_x_m(:); su=f.surface_incident(:); qw=f.quadrature_dx_weights(:).*f.arc_jacobian(:);
 [se,~,~]=geom(sx,amps(c),.1,0); assert(max(abs(se-f.surface_z_m(:)))<1e-12);
 illumination=[illumination;table(repmat(names(c),numel(sx),1),sx,se,real(su),imag(su),qw,VariableNames={'case_name','x_m','eta_m','incident_real','incident_imag','arc_quadrature_weight_m'})]; %#ok<AGROW>
end
writetable(geometry,fullfile(data,'geometry.csv')); writetable(geometry_line,fullfile(data,'geometry_at_receivers.csv'));
writetable(illumination,fullfile(data,'saved_boundary_illumination.csv'));
point=table; physical=table; summary=table; associations=table;
local=find(abs(x(1:411))<=1); edge=find(abs(x(1:411))>=30);
five=ix9([2 4 5 6 8]); supports={(1:411)',local,ix9(:),five(:),edge};
sid=["M99_z3_411","center_abs_x_le_1","P11-P33","P12_P21_P22_P23_P32","edge_abs_x_ge_30"];
for c=1:4
 r=U(:,c,2); f=U(:,1,2); du=r-f;
 ref=zeros(417,1);
 if c>1, rr=D.reference(D.reference.case_name==names(c),:); ref=rr.resolution_estimate; end
 near=abs(r)<=max(max(abs(r))*1e-12,10*ref);
 physical=[physical;table(repmat(names(c),417,1),(1:417)',x,z,abs(r),angle(r),abs(f),angle(f),real(du),imag(du),abs(du),20*log10(abs(r)./abs(f)),angle(r.*conj(f)),(abs(f)+abs(du))./abs(r),ref,near,VariableNames={'case_name','sample_index','x_m','z_m','BIE_rough_abs','BIE_rough_phase_rad','BIE_flat_abs','BIE_flat_phase_rad','delta_real','delta_imag','delta_abs','rough_flat_dB','rough_flat_phase_rad','cancellation_index','BIE_resolution_abs','reference_near_zero'})]; %#ok<AGROW>
 for m=[1 3]
  u=U(:,c,m); e=u-r; ph=angle(u.*conj(r)); a2=(abs(u)-abs(r)).^2; p2=4*abs(u).*abs(r).*sin(ph/2).^2;
  inc=e-(U(:,1,m)-f);
  point=[point;table(repmat(names(c),417,1),repmat(methods(m),417,1),(1:417)',x,z,real(e),imag(e),abs(e),abs(e)./abs(r),20*log10(abs(u)./abs(r)),ph,a2,p2,abs(e).^2,a2+p2-abs(e).^2,real(inc),imag(inc),abs(r).^2,near,VariableNames={'case_name','method','sample_index','x_m','z_m','residual_real','residual_imag','residual_abs','relative_abs','amplitude_dB','wrapped_phase_rad','amplitude_square','phase_square','complex_square','identity_remainder','increment_error_real','increment_error_imag','denominator_abs_square','reference_near_zero'})]; %#ok<AGROW>
  for s=1:numel(supports)
   ii=supports{s}; q=ones(numel(ii),1)/numel(ii); den=sum(q.*abs(r(ii)).^2);
   [worst,wi]=max(abs(e(ii))./abs(r(ii))); wi=ii(wi);
   E=sqrt(sum(q.*abs(e(ii)).^2)/den); Er=sqrt(sum(q.*abs(du(ii)).^2)/sum(q.*abs(f(ii)).^2));
   res=NaN; resolved=false;
   if s==1 || s==3
    a=D.resolutions(D.resolutions.case_name==names(c)&D.resolutions.candidate==methods(m)&D.resolutions.reference=="BIE"&D.resolutions.support_id==sid(s),:);
    res=a.numerical_resolution_estimate; resolved=a.resolved_above_5x_estimate;
   end
   summary=[summary;table(names(c),methods(m),sid(s),numel(ii),"spatial_equal",den,E,sqrt(sum(q.*a2(ii))/den),sqrt(sum(q.*(20*log10(abs(u(ii))./abs(r(ii)))).^2)),sqrt(sum(q.*ph(ii).^2)),sum(q.*a2(ii))/den,sum(q.*p2(ii))/den,Er,worst,x(wi),z(wi),sum(near(ii)),res,resolved,VariableNames={'case_name','method','support_id','n','weight_kind','denominator','complex_L2','amplitude_L2','amplitude_dB_RMS','phase_RMS_rad','amplitude_square_fraction','phase_square_fraction','BIE_physical_increment_L2_flat_denominator','worst_relative_point','worst_x_m','worst_z_m','near_zero_count','empirical_resolution','resolved_above_5x'})]; %#ok<AGROW>
  end
  if c>1
   gg=geometry_line(geometry_line.case_name==names(c),:); descriptors=[abs(gg.eta_m(1:411)),abs(gg.slope(1:411)),abs(gg.curvature_per_m(1:411)),abs(x(1:411)),abs(f(1:411))];
   dn=["abs_projected_eta","abs_projected_slope","abs_projected_curvature","abs_receiver_x","abs_BIE_flat"];
   for k=1:5
    cc=corrcoef(descriptors(:,k),abs(e(1:411)));
    associations=[associations;table(names(c),methods(m),dn(k),cc(1,2),"descriptive_spatial_correlation_no_independent_replicates_no_causality",VariableNames={'case_name','method','descriptor','Pearson_r','interpretation'})]; %#ok<AGROW>
   end
  end
 end
end
writetable(point,fullfile(data,'model_point_diagnostics.csv')); writetable(physical,fullfile(data,'BIE_physical_response.csv'));
writetable(summary,fullfile(data,'support_metrics.csv')); writetable(associations,fullfile(data,'spatial_associations.csv'));
% Freeze reconstructible configurations before any future field observations.
ids=["BASE_M","K_HIGH","PHI_PLUS90","K_LOW","PHI_MINUS90","MATCH_SLOPE","MATCH_CURV"];
AA=[.05 .05 .05 .05 .05 .025 .0125]; KK=[.1 .2 .1 .05 .1 .2 .2]; PP=[0 0 pi/2 0 -pi/2 0 0];
candidate_geometry=table; candidate_stats=table; config_hashes=table;
for j=1:7
 [eta,slope,curv]=geom(gx,AA(j),KK(j),PP(j));
 candidate_geometry=[candidate_geometry;table(repmat(ids(j),numel(gx),1),gx,eta,slope,curv,VariableNames={'case_id','x_m','eta_m','slope','curvature_per_m'})]; %#ok<AGROW>
 mid=find(gx==0); inside=abs(gx)<=50;
 candidate_stats=[candidate_stats;table(ids(j),AA(j),KK(j),PP(j),2*pi/KK(j),100/(2*pi/KK(j)),max(abs(eta)),max(abs(slope)),max(abs(curv)),sqrt(mean(eta(inside).^2)),eta(mid),slope(mid),curv(mid),VariableNames={'case_id','A_m','K_rad_per_m','phi_rad','carrier_period_m','cycles_in_100m_support','actual_max_abs_eta_m','actual_max_abs_slope','actual_max_abs_curvature_per_m','support_RMS_eta_m','center_eta_m','center_slope','center_curvature_per_m'})]; %#ok<AGROW>
 cfg=struct('case_id',ids(j),'A_m',AA(j),'K_rad_per_m',KK(j),'phi_rad',PP(j),'C2_plateau_m',42,'C2_support_m',50,'window_polynomial','1-10t^3+15t^4-6t^5; t=(abs(x)-42)/8','geometry_grid_m',[-60 .01 60],'frequency_hz',4000,'c_m_per_s',1500,'sigma_m',.3,'Tx_m',[0 100],'source_sha256',closeout_s2b_hash(fullfile(protocol,'frozen_source.mat')),'line_sha256',closeout_s2b_hash(fullfile(protocol,'receivers_line.csv')),'nine_sha256',closeout_s2b_hash(fullfile(protocol,'receivers_nine.csv')),'phase','native exp(-i omega t); outgoing exp(+ikR)','quantity','u_reflected_native; unit-source response, not Pa','acoustic_execution_authorized',false);
 fn=fullfile(protocol,ids(j)+"_configuration.json"); closeout_s2b_json(fn,cfg);
 config_hashes=[config_hashes;table(ids(j),string(closeout_s2b_hash(fn)),VariableNames={'case_id','sha256'})]; %#ok<AGROW>
end
writetable(candidate_geometry,fullfile(protocol,'candidate_geometry.csv')); writetable(candidate_stats,fullfile(protocol,'candidate_geometry_stats.csv')); writetable(config_hashes,fullfile(protocol,'configuration_hashes.csv'));
save(fullfile(protocol,'candidate_geometry.mat'),'candidate_geometry','candidate_stats');
% Reuse ONLY saved flat-incident illumination weights on their original GL nodes.
f=D.fk_runs{1,2}; sx=f.surface_x_m(:); sw=abs(f.surface_incident(:)).^2.*f.quadrature_dx_weights(:); sw=sw/sum(sw);
weighted=table;
for j=1:7
 [a,d1,d2]=geom(sx,AA(j),KK(j),PP(j));
 weighted=[weighted;table(ids(j),sqrt(sum(sw.*a.^2)),sqrt(sum(sw.*d1.^2)),sqrt(sum(sw.*d2.^2)),VariableNames={'case_id','saved_flat_incident_weighted_eta_RMS_m','saved_flat_incident_weighted_slope_RMS','saved_flat_incident_weighted_curvature_RMS_per_m'})]; %#ok<AGROW>
end
writetable(weighted,fullfile(protocol,'illumination_weighted_geometry.csv'));
save(fullfile(protocol,'saved_flat_illumination_weights.mat'),'sx','sw');
freeze=struct('stage','S2C_A_DIAGNOSTIC_PROTOCOL_COMPLETE_after_audit','first_candidates',ids(2:3),'baseline',ids(1),'conditional_only',ids(4:7),'new_acoustic_fields',0,'configurations',table2struct(config_hashes),'field_observations_for_candidates','NONE; geometry only','execution_permission','Separate explicit S2c-B authorization required','hypotheses',{{'K_HIGH: compare paired raw BIE response and PE/FK residual versus BASE_M on identical supports; K changes slope, curvature, period and interference together','PHI_PLUS90: counterfactual of zero-height center; center height, slope, curvature and whole illuminated shape change together; no unique height causality'}},'preregistered_metrics',{{'E_line q=1/411; E_array q=1/9; center abs(x)<=1 and five-point cross separate denominators','raw u_rough primary; delta_BIE physical response and model increments secondary','exact amplitude+phase squared decomposition; wrapped phase, near-zero flags, worst point','no p-values; deterministic spatial samples are not independent sea replicates'}},'decision',{{'Retain, weaken or reject candidate mechanisms using paired spatial patterns; no required model PASS','Use 5x sum of matched BIE spatial/window/self, PE probe and FullFK quadrature variations; empirical only','Near-zero or unresolved differences remain explicit; no threshold relaxation, phase fit, new positions or interpolation','If complete trustworthy reference cannot fit cost/RAM gates, stop before new batch; report partial and request a separate budget decision'}},'reference_protocol',{{'For each new surface: BIE full boundary solve on L60 ppw12 qa96 and qa64, same source; new density evaluations count as physics','Matched spatial: L60 ppw14 qa96 against ppw12; matched window: L60 and L70 ppw10 qa96 (two new solves per candidate)','All reference variants evaluated at identical 417 unique positions; estimate sum of absolute spatial, window and self changes at each point','PE original source/grid and corrected negative height screen; half-step and doubled-grid at frozen probes; no wrong-sign reuse','FullFK complete analytic flat exterior GL8 panels400/800, true normals/arclength local_radius=Inf; strict zero outside50; L60/L90 equivalence check'}},'cost',struct('charged_baseline_s',1408.660105,'remaining_s',5791.339895,'suggested_S2c_B_cap_s',1800,'new_budget_s',0,'S3A_reserve_floor_s',3000,'single_job_cap_s',900,'RAM_reserve_fraction',.30,'new_measurement_required','Windows Memory Available Bytes before every physics job; incremental estimate <=0.7*available; not virtual memory'),'order','K_HIGH representative first, all reference/cost gates before PHI_PLUS90; no automatic execution');
closeout_s2b_json(fullfile(protocol,'design_freeze.json'),freeze);
% History is evidence of scale, not a promised cost for new geometry.
ram=table;
for j=1:numel(D.bcost)
 v=D.bcost{j};
 if isfield(v,'available_physical_bytes')
  ram=[ram;table(string(v.job_id),v.node_count,v.self_quadrature_order,v.wall_seconds,v.available_physical_bytes,v.estimated_peak_increment_bytes,VariableNames={'job_id','nodes','self_order','wall_seconds','available_physical_bytes','estimated_increment_bytes'})]; %#ok<AGROW>
 end
end
writetable(ram,fullfile(protocol,'historical_RAM_cost_receipts.csv'));
% Per candidate five new BIE solves, so two fully qualified cases may not fit.
jobs=table(["K_HIGH";"PHI_PLUS90"],[5;5],[845;845],[1450;1450],[false;false],VariableNames={'case_id','new_BIE_solves_for_all_layers','BIE_order_estimate_low_s','BIE_order_estimate_high_s','authorized'});
writetable(jobs,fullfile(protocol,'pending_jobs.csv'));
source_index=table(string(old),string(closeout_s2b_hash(fullfile(old,'reference_data.mat'))),string(S.status),"FullFK_complete_analytic_flat_exterior_GL8_v1",xzcheck(1),xzcheck(2),VariableNames={'formal_source_directory','reference_data_sha256','S2b_status','FullFK_version','strong_xz_PE_z3_max_difference','strong_xz_BIE_z3_max_difference'});
writetable(source_index,fullfile(data,'source_index.csv'));
plot_data=struct('U',U,'x',x,'z',z,'ix9',ix9,'grid',grid,'names',names,'methods',methods,'geometry',geometry,'illumination',illumination,'point',point,'physical',physical,'summary',summary,'local',local,'reference',D.reference,'resolutions',D.resolutions);
save(fullfile(data,'diagnostic_data.mat'),'plot_data','src','candidate_stats','source_index','-v7');
closeout_s2b_json(fullfile(out,'analysis_receipt.json'),struct('status','PENDING_INDEPENDENT_AUDIT','new_acoustic_fields',0,'analysis_IO_seconds',toc(t0),'old_output_hashes_verified',height(H),'old_source_code_hashes_verified',height(C),'S2b_charged_seconds_unchanged',1408.660105,'remaining_seconds_unchanged',5791.339895,'new_physical_permission',false));
render_two_dimensional_closeout_s2c_a(fullfile(data,'diagnostic_data.mat'),fullfile(out,'figures'));
end

function [eta,d1,d2]=geom(x,A,K,p)
r=abs(x); c=double(r<=42); cp=zeros(size(x)); cpp=cp;
ii=r>42&r<50; t=(r(ii)-42)/8;
c(ii)=1-10*t.^3+15*t.^4-6*t.^5;
cp(ii)=(-30*t.^2+60*t.^3-30*t.^4).*sign(x(ii))/8;
cpp(ii)=(-60*t+180*t.^2-120*t.^3)/64;
s=sin(K*x+p); co=cos(K*x+p);
eta=A*s.*c; d1=A*(K*co.*c+s.*cp);
d2=A*(-K^2*s.*c+2*K*co.*cp+s.*cpp);
end
