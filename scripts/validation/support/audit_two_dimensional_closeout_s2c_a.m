function a=audit_two_dimensional_closeout_s2c_a(out)
% Independent Cartesian CSV arithmetic; no main metric/geometry functions.
dest=fullfile(out,'audit'); assert(~isfolder(dest)); mkdir(dest);
data=fullfile(out,'mechanism_baseline');
L=readtable(fullfile(data,'fields_line.csv'),TextType='string');
N=readtable(fullfile(data,'fields_core.csv'),TextType='string');
Q=readtable(fullfile(data,'model_point_diagnostics.csv'),TextType='string');
B=readtable(fullfile(data,'BIE_physical_response.csv'),TextType='string');
T=readtable(fullfile(data,'support_metrics.csv'),TextType='string');
R=load(fullfile(data,'diagnostic_data.mat'),'plot_data'); d=R.plot_data;
fieldmax=0; pointmax=0; physicalmax=0; metricmax=0; identitymax=0;
for c=1:4
 for m=1:3
  ll=L(L.case_name==d.names(c)&L.method==d.methods(m),:);
  nn=N(N.case_name==d.names(c)&N.method==d.methods(m),:);
  fieldmax=max([fieldmax;abs(complex(ll.u_real,ll.u_imag)-d.U(1:411,c,m));abs(complex(nn.u_real,nn.u_imag)-d.U(d.ix9,c,m))]);
 end
 % Assemble the 417 unique positions using true CSV values, no interpolation.
 ul=zeros(417,3); uf=ul;
 for m=1:3
  p=L(L.case_name==d.names(c)&L.method==d.methods(m),:); n=N(N.case_name==d.names(c)&N.method==d.methods(m),:);
  p0=L(L.case_name=="flat"&L.method==d.methods(m),:); n0=N(N.case_name=="flat"&N.method==d.methods(m),:);
  off=n.z_m~=3;
  ul(:,m)=[complex(p.u_real,p.u_imag);complex(n.u_real(off),n.u_imag(off))];
  uf(:,m)=[complex(p0.u_real,p0.u_imag);complex(n0.u_real(off),n0.u_imag(off))];
 end
 rr=ul(:,2); ff=uf(:,2); pp=B(B.case_name==d.names(c),:);
 delta=rr-ff; rot=atan2(imag(rr).*real(ff)-real(rr).*imag(ff),real(rr).*real(ff)+imag(rr).*imag(ff));
 physicalmax=max([physicalmax;abs(delta-complex(pp.delta_real,pp.delta_imag));abs(abs(rr)-pp.BIE_rough_abs);abs(rot-pp.rough_flat_phase_rad);abs((abs(ff)+abs(delta))./abs(rr)-pp.cancellation_index)]);
 for m=[1 3]
  u=ul(:,m); er=real(u)-real(rr); ei=imag(u)-imag(rr); cs=er.^2+ei.^2;
  ab=sqrt(real(u).^2+imag(u).^2); rb=sqrt(real(rr).^2+imag(rr).^2);
  phase=atan2(imag(u).*real(rr)-real(u).*imag(rr),real(u).*real(rr)+imag(u).*imag(rr));
  amp=(ab-rb).^2; ps=2*ab.*rb.*(1-cos(phase));
  qq=Q(Q.case_name==d.names(c)&Q.method==d.methods(m),:);
  pointmax=max([pointmax;abs(er-qq.residual_real);abs(ei-qq.residual_imag);abs(phase-qq.wrapped_phase_rad);abs(20/log(10)*log(ab./rb)-qq.amplitude_dB);abs(cs-qq.complex_square);abs(amp-qq.amplitude_square);abs(ps-qq.phase_square)]);
  identitymax=max([identitymax;abs(amp+ps-cs)]);
  st=T(T.case_name==d.names(c)&T.method==d.methods(m),:);
  for j=1:height(st)
   switch st.support_id(j)
    case "M99_z3_411", ii=(1:411)';
    case "center_abs_x_le_1", ii=find(abs(d.x(1:411))<=1);
    case "P11-P33", ii=d.ix9;
    case "P12_P21_P22_P23_P32", ii=d.ix9([2 4 5 6 8]);
    case "edge_abs_x_ge_30", ii=find(abs(d.x(1:411))>=30);
    otherwise, error('Unfrozen support');
   end
   assert(numel(ii)==st.n(j)); den=mean(rb(ii).^2);
   vals=[sqrt(mean(cs(ii))/den),sqrt(mean(amp(ii))/den),sqrt(mean((20/log(10)*log(ab(ii)./rb(ii))).^2)),sqrt(mean(phase(ii).^2)),mean(amp(ii))/den,mean(ps(ii))/den,sqrt(mean(abs(delta(ii)).^2)/mean(abs(ff(ii)).^2))];
   saved=[st.complex_L2(j),st.amplitude_L2(j),st.amplitude_dB_RMS(j),st.phase_RMS_rad(j),st.amplitude_square_fraction(j),st.phase_square_fraction(j),st.BIE_physical_increment_L2_flat_denominator(j)];
   metricmax=max([metricmax,abs(vals-saved),abs(den-st.denominator(j))]);
  end
 end
end
% Independent analytic C2 formula using factored derivatives.
G=readtable(fullfile(out,'protocol/candidate_geometry.csv'),TextType='string');
J=readtable(fullfile(out,'protocol/candidate_geometry_stats.csv'),TextType='string');
CH=readtable(fullfile(out,'protocol/configuration_hashes.csv'),TextType='string');
geommax=0; derivcheck=0;
for j=1:height(J)
 v=G(G.case_id==J.case_id(j),:); x=v.x_m; r=abs(x); w=zeros(size(x)); w(r<=42)=1; dw=zeros(size(x)); ddw=dw;
 idx=r>42&r<50; t=(r(idx)-42)/8;
 w(idx)=(1-t).^3.*(1+3*t+6*t.^2);
 dw(idx)=-30*t.^2.*(1-t).^2.*sign(x(idx))/8;
 ddw(idx)=-60*t.*(1-t).*(1-2*t)/64;
 aa=J.A_m(j); kk=J.K_rad_per_m(j); ph=J.phi_rad(j); s=sin(kk*x+ph); co=cos(kk*x+ph);
 y=aa*s.*w; dy=aa*kk*co.*w+aa*s.*dw;
 ddy=-aa*kk^2*s.*w+2*aa*kk*co.*dw+aa*s.*ddw;
 geommax=max([geommax;abs(y-v.eta_m);abs(dy-v.slope);abs(ddy-v.curvature_per_m)]);
 % Finite differences test analytic derivatives only, away from C2 junctions.
 dd=diff(v.eta_m)/.01; avg=(v.slope(1:end-1)+v.slope(2:end))/2;
 derivcheck=max(derivcheck,max(abs(dd-avg)));
 cfg=fullfile(out,'protocol',J.case_id(j)+"_configuration.json");
 assert(strcmp(closeout_s2b_hash(cfg),CH.sha256(CH.case_id==J.case_id(j))));
 cf=jsondecode(fileread(cfg)); assert(max(abs([cf.A_m-aa,cf.K_rad_per_m-kk,cf.phi_rad-ph]))<1e-14);
 assert(abs(max(abs(y))-J.actual_max_abs_eta_m(j))<1e-13);
 assert(abs(max(abs(dy))-J.actual_max_abs_slope(j))<1e-13);
 assert(abs(max(abs(ddy))-J.actual_max_abs_curvature_per_m(j))<1e-13);
end
assert(fieldmax<2e-14 && pointmax<1e-11 && physicalmax<1e-11 && metricmax<1e-11 && geommax<1e-13 && derivcheck<1e-6);
assert(identitymax<1e-16);
% Check the actual M00/M01 baseline geometry against a second formula too.
BG=readtable(fullfile(data,'geometry.csv'),TextType='string');
for c=1:4
 v=BG(BG.case_name==d.names(c),:); r=abs(v.x_m); w=double(r<=42); w1=zeros(size(r)); w2=w1;
 ix=r>42&r<50; t=(r(ix)-42)/8;
 w(ix)=(1-t).^3.*(1+3*t+6*t.^2); w1(ix)=-30*t.^2.*(1-t).^2.*sign(v.x_m(ix))/8; w2(ix)=-60*t.*(1-t).*(1-2*t)/64;
 amps=[0 .01 .05 .2]; aa=amps(c); s=sin(.1*v.x_m); co=cos(.1*v.x_m);
 check=[aa*s.*w,aa*.1*co.*w+aa*s.*w1,-aa*.01*s.*w+aa*.2*co.*w1+aa*s.*w2];
 actual=[v.eta_m,v.slope,v.curvature_per_m]; assert(max(abs(check-actual),[],'all')<1e-13);
 saved=d.geometry(d.geometry.case_name==d.names(c),:); assert(max(abs(actual-[saved.eta_m,saved.slope,saved.curvature_per_m]),[],'all')<1e-13);
end
old=load(fullfile(formal_source_directory(out),'reference_data.mat'),'fk_runs');
ill=readtable(fullfile(data,'saved_boundary_illumination.csv'),TextType='string');
for c=1:4
 v=ill(ill.case_name==d.names(c),:); f=old.fk_runs{c,2};
 assert(max(abs(complex(v.incident_real,v.incident_imag)-f.surface_incident(:)))<1e-13);
 assert(max(abs(v.x_m-f.surface_x_m(:)))<1e-12);
end
% Preserve every old result again, after all analysis/rendering.
H=readtable(fullfile(out,'provenance/S2b_preserved_hashes.csv'),TextType='string');
for j=1:height(H), assert(strcmp(closeout_s2b_hash(H.path(j)),H.sha256(j))); end
redraw=fullfile(dest,'redraw'); render_two_dimensional_closeout_s2c_a(fullfile(data,'diagnostic_data.mat'),redraw);
F=dir(fullfile(out,'figures/*.png')); assert(numel(F)==11); pixelmax=0; phasebreaks=true; reopen=0;
for j=1:numel(F)
 p1=imread(fullfile(F(j).folder,F(j).name)); p2=imread(fullfile(redraw,F(j).name));
 assert(isequal(size(p1),size(p2))); pixelmax=max(pixelmax,max(abs(double(p1(:))-double(p2(:)))));
 [~,stem]=fileparts(F(j).name); f=openfig(fullfile(out,'figures',[stem '.fig']),'new','invisible'); reopen=reopen+1;
 ax=findall(f,'Type','axes');
 for k=1:numel(ax)
  label=string(ax(k).YLabel.String);
  if contains(label,'wrapped') || contains(label,'arg(u')
   ls=findall(ax(k),'Type','line');
   for n=1:numel(ls), yy=ls(n).YData; phasebreaks=phasebreaks&&all(abs(diff(yy))<=pi|isnan(diff(yy))); end
  end
 end
 close(f);
end
assert(pixelmax==0&&phasebreaks);
atlas=fileread(fullfile(out,'mechanism_atlas.md'));
links=regexp(atlas,'\]\(([^)]+\.(?:png|fig))\)','tokens');
for j=1:numel(links), assert(isfile(fullfile(out,links{j}{1}))); end
assert(numel(links)==22);
a=struct('pass',true,'line_rows',height(L),'core_rows',height(N),'model_point_rows',height(Q),'physical_response_rows',height(B),'support_metric_rows',height(T),'max_CSV_MAT_field_difference',fieldmax,'max_independent_point_difference',pointmax,'max_independent_physical_response_difference',physicalmax,'max_independent_summary_difference',metricmax,'max_amplitude_phase_identity',identitymax,'max_independent_geometry_difference',geommax,'geometry_finite_difference_slope_error',derivcheck,'old_output_hashes_preserved',height(H),'configuration_hashes_verified',height(CH),'PNG_count',numel(F),'FIG_reopened',reopen,'redraw_max_pixel_difference',pixelmax,'wrapped_phase_breaks_verified',phasebreaks,'atlas_relative_PNG_FIG_links',numel(links),'new_acoustic_field_computations',0,'no_candidate_field_observations',true,'support_rule','Only P21/P22/P23 are z=3; original six off-depth rows remain separate; 411/local/9/five/edge denominators explicit','numerical_scope','empirical resolution; unresolved weak/medium FullFK nine-point differences retained');
closeout_s2b_json(fullfile(dest,'independent_audit.json'),a);
end

function old=formal_source_directory(out)
% Read a numeric provenance path, never a text-specified executable.
s=readtable(fullfile(out,'mechanism_baseline/source_index.csv'),TextType='string',Delimiter=',');
old=char(s.formal_source_directory);
end
