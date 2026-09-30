function audit = validate_incident_reflection_field_covariance()
%VALIDATE_INCIDENT_REFLECTION_FIELD_COVARIANCE
% Validation-only audit of PE incident-field covariance and reflected-field
% differences against the saved Full Kirchhoff/BIE receiver-depth results.
% No production PE, Full Kirchhoff, BIE, or surface model is modified.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root); setup_vertical_project();
c = local_cfg(root);
if ~exist(c.out,'dir'), mkdir(c.out); end
if ~exist(c.fig,'dir'), mkdir(c.fig); end

q = load(c.input,'result'); base = q.result;
if exist(c.audit_input,'file')
    a = load(c.audit_input,'audit'); fixed = a.audit;
else
    fixed = [];
end

% One flat PE run supplies the incident field at the reference plane. The
% field is independent of the rough profile because reflection is applied
% only after this one-way propagation in the validation helper.
flat_cfg = struct('frequency_hz',c.f,'c0_mps',c.c,'xw_m',c.xw,'nx',c.nx, ...
    'z_tx_m',c.ztx,'z_rx_m',3,'sigma_src_m',c.sigma, ...
    'surface_elevation_x_m',zeros(1,c.nx),'surface_reflect_coeff',-1, ...
    'step_m',c.step,'x_rx_m',0,'reflection_model','model1_kz_aware');
pe0 = run_pe_1d_surface_reflection_validation(flat_cfg);
x = base.cases(1).x_m(:); pe_ref = pe0.surface_incident_field(:);
if numel(x) ~= numel(pe_ref), error('Saved and PE incident grids differ.'); end

out_cases = repmat(local_case(),numel(base.cases),1);
for ic = 1:numel(base.cases)
    out_cases(ic) = local_case_audit(c,base.cases(ic),pe_ref,fixed);
end

audit = struct('schema_version','1.0.0','config',c,'cases',out_cases, ...
    'incident_reference',struct('x_m',x,'pe_reference_plane',pe_ref), ...
    'source_note','FK/BIE use the same saved angular-spectrum source field; BIE does not store a separate incident trace.');
audit.table = local_table(out_cases);
mat = fullfile(c.out,'incident_reflection_covariance.mat');
csv = fullfile(c.out,'incident_reflection_covariance.csv');
save(mat,'audit','-v7.3'); writetable(audit.table,csv);
local_write_spectra(c,out_cases);
local_plots(c,out_cases);
local_report(c,audit,root,mat,csv);
disp(audit.table);
end

function c = local_cfg(root)
c = struct('f',4000,'c',1500,'xw',192.1875,'nx',984,'ztx',100,'sigma',1, ...
    'step',.05,'floor_db',-40,'z_select',[2 3], ...
    'case_names',{{'weak_low_K','strong_height_low_K','weak_high_K'}}, ...
    'input',fullfile(root,'results','validation','pe_bie_full_kirchhoff_receiver_depth_sensitivity','receiver_depth_sensitivity.mat'), ...
    'audit_input',fullfile(root,'results','validation','pe_bie_receiver_depth_numerical_audit','receiver_depth_numerical_audit.mat'), ...
    'out',fullfile(root,'results','validation','pe_bie_incident_reflection_covariance'), ...
    'fig',fullfile(root,'results','validation','pe_bie_incident_reflection_covariance','figures'), ...
    'report',fullfile(root,'reports','pe_bie_incident_reflection_covariance_report.md'));
end

function s = local_case()
s = struct('case_name','','A_m',NaN,'K_radpm',NaN,'x_m',[],'surface_x_m',[], ...
    'surface_z_m',[],'z_select',[],'incident',[],'reflection',[],'spectra',[],'mask_source','');
end

function s = local_case_audit(c,old,pe_ref,fixed)
x = old.x_m(:); z = old.z_rx_m(:); idx = zeros(1,numel(c.z_select));
for ii=1:numel(idx), [~,idx(ii)] = min(abs(z-c.z_select(ii))); end
fp = old.footprint; mask = fp.m99.mask(:); w = fp.energy_weights(:); w = w/max(sum(w(mask)),realmin);
src = old.source; sx = old.surface_x_m(:); sz = old.surface_z_m(:);
sz_at_x = interp1(sx,sz,x,'pchip',0);
carrier = exp(1i*(2*pi*c.f/c.c)*c.ztx);
% PE field at z=0 and its continuation to the actual surface. FK source
% evaluation is exactly the angular-spectrum field used by the FK solver.
pe_actual = carrier*local_variable_z(pe_ref,x,sz_at_x,c.f,c.c);
pe_ref_aligned = carrier*pe_ref;
fk_actual = src.eval(x,sz_at_x,c.ztx,'coordinate'); fk_actual = fk_actual(:);
fk_ref = src.eval(x,zeros(size(x)),c.ztx,'coordinate'); fk_ref = fk_ref(:);
incident = struct('reference_plane',local_pair(pe_ref_aligned,fk_ref,mask,w,c.floor_db), ...
    'actual_surface',local_pair(pe_actual,fk_actual,mask,w,c.floor_db), ...
    'pe_reference_plane_raw',pe_ref,'pe_reference_plane',pe_ref_aligned,'carrier_factor',carrier,'fk_reference_plane',fk_ref, ...
    'pe_actual_surface',pe_actual,'fk_actual_surface',fk_actual, ...
    'kx_reference',local_spectrum_pair(pe_ref,fk_ref,x,c,mask));

Gpe = old.G_Model1; Gfk = old.G_FK; Gbie = old.G_BIE;
reflection = repmat(struct('z_rx_m',NaN,'PE_FK',[],'PE_BIE',[],'FK_BIE',[], ...
    'fields',[],'spectrum',[]),numel(idx),1);
for jj=1:numel(idx)
    iz=idx(jj); pp=Gpe(:,iz); ff=Gfk(:,iz); bb=Gbie(:,iz);
    reflection(jj).z_rx_m=z(iz);
    reflection(jj).PE_FK=local_pair(pp,ff,mask,w,c.floor_db);
    reflection(jj).PE_BIE=local_pair(pp,bb,mask,w,c.floor_db);
    reflection(jj).FK_BIE=local_pair(ff,bb,mask,w,c.floor_db);
    reflection(jj).fields=struct('PE',pp,'FK',ff,'BIE',bb);
    reflection(jj).spectrum=local_spectrum_triplet(pp,ff,bb,x,c,mask);
end

s = local_case(); s.case_name=old.case_name; s.A_m=old.A_m; s.K_radpm=old.K_radpm;
s.x_m=x; s.surface_x_m=sx; s.surface_z_m=sz; s.z_select=z(idx);
s.incident=incident; s.reflection=reflection; s.spectra=reflection;
if ~isempty(fixed) && isfield(fixed,'cases')
    % Keep the prior fixed-mask audit identifier for traceability only.
    s.mask_source='receiver_depth_numerical_audit.mat';
else
    s.mask_source='receiver_depth_sensitivity.mat footprint.m99';
end
end

function u = local_variable_z(plane,x,zq,f,c)
k=2*pi*f/c; n=numel(x); dx=median(diff(x)); x0=x(1); kx=2*pi*[0:n/2-1,-n/2:-1]/(n*dx); kz=sqrt(complex(k^2-kx.^2,0)); ck=fft(plane(:)).';
u=zeros(n,1); block=128;
for first=1:block:n
    rows=first:min(first+block-1,n); phx=exp(1i*(x(rows)-x0)*kx); phz=exp(-1i*zq(rows).*kz); u(rows)=sum(phx.*(ck.*phz),2)/n;
end
end

function p=local_pair(a,b,mask,w,floor_db)
a=a(:); b=b(:); base=mask&isfinite(a)&isfinite(b); w0=w(base); w0=w0/max(sum(w0),realmin); a=a(base); b=b(base); pk=max([abs(a);abs(b)]); ok=abs(a)>=10^(floor_db/20)*pk & abs(b)>=10^(floor_db/20)*pk; a=a(ok); b=b(ok); w0=w0(ok); w0=w0/max(sum(w0),realmin); ph=angle(a.*conj(b)); ma=abs(a); mb=abs(b); p=struct('complex_l2',sqrt(sum(w0.*abs(a-b).^2)/max(sum(w0.*mb.^2),realmin)),'magnitude_l2',sqrt(sum(w0.*(ma-mb).^2)/max(sum(w0.*mb.^2),realmin)),'phase_rms',sqrt(sum(w0.*ph.^2)),'magnitude_corr',local_corr(ma,mb,w0),'phase_corr',abs(sum(w0.*exp(1i*ph))),'n',numel(a));
end
function r=local_corr(a,b,w), a=a-sum(w.*a); b=b-sum(w.*b); r=abs(sum(w.*a.*b))/max(sqrt(sum(w.*a.^2)*sum(w.*b.^2)),realmin); end

function S=local_spectrum_pair(a,b,x,c,mask)
[kx,wa,za]=local_fft(a,x,mask); [~,wb,zb]=local_fft(b,x,mask); S=struct('kx_radpm',kx,'kx_over_k',kx/(2*pi*c.f/c.c),'PE_amp',wa,'FK_amp',wb,'relative_amp_error',abs(wa-wb),'phase_error',abs(angle(za.*conj(zb))));
end
function S=local_spectrum_triplet(a,b,d,x,c,mask)
[kx,wa,za]=local_fft(a,x,mask); [~,wb,zb]=local_fft(b,x,mask); [~,wd,zd]=local_fft(d,x,mask); S=struct('kx_radpm',kx,'kx_over_k',kx/(2*pi*c.f/c.c),'PE_amp',wa,'FK_amp',wb,'BIE_amp',wd,'PE_FK_amp_error',abs(wa-wb),'PE_BIE_amp_error',abs(wa-wd),'FK_BIE_amp_error',abs(wb-wd),'PE_FK_phase_error',abs(angle(za.*conj(zb))),'PE_BIE_phase_error',abs(angle(za.*conj(zd))),'FK_BIE_phase_error',abs(angle(zb.*conj(zd))));
end
function [kx,a,z]=local_fft(v,x,mask)
ii=find(mask); i1=ii(1); i2=ii(end); win=zeros(size(x)); nn=i2-i1+1; win(i1:i2)=0.5-0.5*cos(2*pi*(0:nn-1)'/max(nn-1,1)); y=v(:).*win; n=numel(x); dx=median(diff(x)); kx=2*pi*((-floor(n/2)):(ceil(n/2)-1))'/(n*dx); a=fftshift(fft(y)); a=a/max(norm(a),realmin); z=a; a=abs(z);
end

function T=local_table(C)
rows=repmat(struct('case_name','','A_m',NaN,'K_radpm',NaN,'incident_ref_Ec',NaN,'incident_surface_Ec',NaN,'z_rx_m',NaN,'PE_FK_Ec',NaN,'PE_BIE_Ec',NaN,'FK_BIE_Ec',NaN,'PE_FK_mag',NaN,'PE_FK_phase',NaN,'PE_BIE_mag',NaN,'PE_BIE_phase',NaN,'FK_BIE_mag',NaN,'FK_BIE_phase',NaN),0,1);
for i=1:numel(C)
  for j=1:numel(C(i).reflection)
    r=local_row(); r.case_name=C(i).case_name; r.A_m=C(i).A_m; r.K_radpm=C(i).K_radpm; r.incident_ref_Ec=C(i).incident.reference_plane.complex_l2; r.incident_surface_Ec=C(i).incident.actual_surface.complex_l2; r.z_rx_m=C(i).reflection(j).z_rx_m; p=C(i).reflection(j); r.PE_FK_Ec=p.PE_FK.complex_l2; r.PE_BIE_Ec=p.PE_BIE.complex_l2; r.FK_BIE_Ec=p.FK_BIE.complex_l2; r.PE_FK_mag=p.PE_FK.magnitude_l2; r.PE_FK_phase=p.PE_FK.phase_rms; r.PE_BIE_mag=p.PE_BIE.magnitude_l2; r.PE_BIE_phase=p.PE_BIE.phase_rms; r.FK_BIE_mag=p.FK_BIE.magnitude_l2; r.FK_BIE_phase=p.FK_BIE.phase_rms; rows(end+1)=r; %#ok<AGROW>
  end
end
T=struct2table(rows);
end
function r=local_row(), r=struct('case_name','','A_m',NaN,'K_radpm',NaN,'incident_ref_Ec',NaN,'incident_surface_Ec',NaN,'z_rx_m',NaN,'PE_FK_Ec',NaN,'PE_BIE_Ec',NaN,'FK_BIE_Ec',NaN,'PE_FK_mag',NaN,'PE_FK_phase',NaN,'PE_BIE_mag',NaN,'PE_BIE_phase',NaN,'FK_BIE_mag',NaN,'FK_BIE_phase',NaN); end

function local_write_spectra(c,C)
for i=1:numel(C)
 s=C(i).incident.kx_reference; T=table(s.kx_radpm,s.kx_over_k,s.PE_amp,s.FK_amp,s.relative_amp_error,s.phase_error,'VariableNames',{'kx_radpm','kx_over_k','PE_amp','FK_amp','relative_amp_error','phase_error'}); writetable(T,fullfile(c.out,[C(i).case_name '_incident_spectrum.csv']));
 for j=1:numel(C(i).reflection)
  s=C(i).reflection(j).spectrum; T=table(s.kx_radpm,s.kx_over_k,s.PE_amp,s.FK_amp,s.BIE_amp,s.PE_FK_amp_error,s.PE_BIE_amp_error,s.FK_BIE_amp_error,s.PE_FK_phase_error,s.PE_BIE_phase_error,s.FK_BIE_phase_error,'VariableNames',{'kx_radpm','kx_over_k','PE_amp','FK_amp','BIE_amp','PE_FK_amp_error','PE_BIE_amp_error','FK_BIE_amp_error','PE_FK_phase_error','PE_BIE_phase_error','FK_BIE_phase_error'}); writetable(T,fullfile(c.out,sprintf('%s_reflection_z%g_spectrum.csv',C(i).case_name,C(i).reflection(j).z_rx_m)));
 end
end
end

function local_plots(c,C)
for i=1:numel(C)
 s=C(i).incident; f=figure('Visible','off'); tiledlayout(2,2); nexttile; plot(C(i).x_m,abs(s.pe_actual_surface),'--',C(i).x_m,abs(s.fk_actual_surface),'-'); grid on; legend('PE continuation','FK/BIE angular spectrum'); title([C(i).case_name ' incident magnitude']); xlabel('x (m)'); nexttile; plot(C(i).x_m,unwrap(angle(s.pe_actual_surface)),'--',C(i).x_m,unwrap(angle(s.fk_actual_surface)),'-'); grid on; title('incident phase'); xlabel('x (m)'); nexttile([1 2]); plot(s.kx_reference.kx_over_k,s.kx_reference.PE_amp,'--',s.kx_reference.kx_over_k,s.kx_reference.FK_amp,'-'); xlim([-1 1]); grid on; xlabel('k_x/k'); ylabel('normalized |P(k_x)|'); title('incident angular spectrum'); saveas(f,fullfile(c.fig,[C(i).case_name '_incident.png'])); saveas(f,fullfile(c.fig,[C(i).case_name '_incident.pdf'])); close(f);
 for j=1:numel(C(i).reflection)
  rr=C(i).reflection(j); z=rr.z_rx_m; f=figure('Visible','off'); tiledlayout(2,2); nexttile; plot(C(i).x_m,abs(rr.fields.PE),'--',C(i).x_m,abs(rr.fields.FK),'-',C(i).x_m,abs(rr.fields.BIE),':'); grid on; legend('PE Model-1','Full Kirchhoff','BIE'); title(sprintf('%s reflected magnitude z=%g m',C(i).case_name,z),'Interpreter','none'); nexttile; plot(C(i).x_m,angle(rr.fields.PE),'--',C(i).x_m,angle(rr.fields.FK),'-',C(i).x_m,angle(rr.fields.BIE),':'); grid on; title('reflected phase'); nexttile; plot(C(i).x_m,abs(rr.fields.PE-rr.fields.FK)); grid on; title('|PE-FK|'); nexttile; plot(C(i).x_m,angle(rr.fields.PE.*conj(rr.fields.FK))); grid on; title('angle(PE conj(FK))'); saveas(f,fullfile(c.fig,sprintf('%s_reflection_z%g.png',C(i).case_name,z))); saveas(f,fullfile(c.fig,sprintf('%s_reflection_z%g.pdf',C(i).case_name,z))); close(f);
  f=figure('Visible','off'); plot(rr.spectrum.kx_over_k,rr.spectrum.PE_amp,'--',rr.spectrum.kx_over_k,rr.spectrum.FK_amp,'-',rr.spectrum.kx_over_k,rr.spectrum.BIE_amp,':'); xlim([-1 1]); grid on; xlabel('k_x/k'); ylabel('normalized |P(k_x)|'); title(sprintf('%s reflected spectrum z=%g m',C(i).case_name,z),'Interpreter','none'); legend('PE','FK','BIE'); saveas(f,fullfile(c.fig,sprintf('%s_reflection_z%g_spectrum.png',C(i).case_name,z))); saveas(f,fullfile(c.fig,sprintf('%s_reflection_z%g_spectrum.pdf',C(i).case_name,z))); close(f);
  f=figure('Visible','off'); plot(rr.spectrum.kx_over_k,rr.spectrum.PE_FK_phase_error,'-',rr.spectrum.kx_over_k,rr.spectrum.PE_BIE_phase_error,'--',rr.spectrum.kx_over_k,rr.spectrum.FK_BIE_phase_error,':'); xlim([-1 1]); grid on; xlabel('k_x/k'); ylabel('|angle(P conj(Q))| (rad)'); title(sprintf('%s reflected kx phase differences z=%g m',C(i).case_name,z),'Interpreter','none'); legend('PE-FK','PE-BIE','FK-BIE'); saveas(f,fullfile(c.fig,sprintf('%s_reflection_z%g_spectrum_phase.png',C(i).case_name,z))); saveas(f,fullfile(c.fig,sprintf('%s_reflection_z%g_spectrum_phase.pdf',C(i).case_name,z))); close(f);
 end
end
end

function local_report(c,a,root,mat,csv)
fid=fopen(c.report,'w','n','UTF-8'); cl=onCleanup(@()fclose(fid)); fprintf(fid,'# PE / Full Kirchhoff / BIE incident and reflected field covariance audit\n\n'); fprintf(fid,'Validation-only audit at 4 kHz, c=1500 m/s, Tx=100 m, sigma=1 m. Cases are weak_low_K (A=%.3g, K=%.3g), strong_height_low_K (A=%.3g, K=%.3g), and weak_high_K (A=%.3g, K=%.3g). Existing receiver-depth fields are reused at z=2 and 3 m. No production physics was changed.\n\n',a.cases(1).A_m,a.cases(1).K_radpm,a.cases(2).A_m,a.cases(2).K_radpm,a.cases(3).A_m,a.cases(3).K_radpm);
 fprintf(fid,'## Incident field check\n\nAt the reference plane z=0, PE `surface_incident_field` is compared with the same angular-spectrum source evaluation used by Full Kirchhoff. At the true rough surface, PE is continued from z=0 with the same one-way kz spectrum and compared with the saved FK incident field. The BIE solver does not store an independent incident trace; its source term is the same angular-spectrum construction, so FK/BIE incident agreement is a formulation identity rather than a separately measured receiver field.\n\n');
 for i=1:numel(a.cases), fprintf(fid,'- `%s`: reference-plane Ec=%.6g, true-surface Ec=%.6g; reference-plane and true-surface plots are in `%s`.\n',a.cases(i).case_name,a.cases(i).incident.reference_plane.complex_l2,a.cases(i).incident.actual_surface.complex_l2,rel(root,c.fig)); end
 fprintf(fid,'\n## Reflected field checks\n\n'); for i=1:numel(a.cases), for j=1:numel(a.cases(i).reflection), r=a.cases(i).reflection(j); fprintf(fid,'- `%s`, z=%g m: PE--FK Ec=%.6g (magnitude %.6g, phase %.6g), PE--BIE Ec=%.6g, FK--BIE Ec=%.6g.\n',a.cases(i).case_name,r.z_rx_m,r.PE_FK.complex_l2,r.PE_FK.magnitude_l2,r.PE_FK.phase_rms,r.PE_BIE.complex_l2,r.FK_BIE.complex_l2); end; end
 fprintf(fid,'\n### kx-region diagnostics\n\nThe following unweighted spectral summaries use the common M99/Hann support. They are diagnostics because bins with very small spectral amplitude can have unstable phase. Each tuple is `(mean amplitude difference, RMS phase difference)` for `|kx|/k` bands `[0,.25)`, `[.25,.5)`, `[.5,.8)`, `[.8,1)`.\n\n'); regs=[0 .25;.25 .5;.5 .8;.8 1]; for i=1:numel(a.cases), for j=1:numel(a.cases(i).reflection), s=a.cases(i).reflection(j).spectrum; kk=abs(s.kx_over_k); fprintf(fid,'- `%s`, z=%g m PE--FK: ',a.cases(i).case_name,a.cases(i).reflection(j).z_rx_m); for ir=1:size(regs,1), ix=kk>=regs(ir,1)&kk<regs(ir,2); fprintf(fid,'[%.2f,%.2f): %.3g / %.3g; ',regs(ir,1),regs(ir,2),mean(s.PE_FK_amp_error(ix)),sqrt(mean(s.PE_FK_phase_error(ix).^2))); end; fprintf(fid,'\n'); end; end
 fprintf(fid,'\n## Interpretation\n\nThe incident-field test checks the boundary input before reflection. If its fixed-mask complex error is small while the reflected PE--FK/BIE error is larger, the discrepancy begins in the reflection operator or in its spatial/spectral reconstruction. The saved spectra report normalized |P(k_x)| and PE--FK/FK--BIE spectral amplitude differences for the same M99 support; inspect the high-|k_x| tail first. Strong height is expected to emphasize path-phase and nonlocal coherence, while high K can additionally redistribute energy across kx.\n\nFor a reduced or fast Kirchhoff operator, retain the true surface path phase and the nonlocal coherent kernel. The incident angular-spectrum representation and its obliquity/normal derivative should remain explicit; a purely local height phase is insufficient when the reflected spectrum changes.\n\n## Artifacts\n\n- MAT: `%s`\n- CSV: `%s`\n- Figures: `%s`\n\nThe receiver samples are the saved discrete x grid; no depth interpolation was introduced. The audit does not claim that z=1 m is validated.\n',rel(root,mat),rel(root,csv),rel(root,c.fig));
end
function p=rel(root,p), p=strrep(p,[root filesep],''); p=strrep(p,filesep,'/'); end
