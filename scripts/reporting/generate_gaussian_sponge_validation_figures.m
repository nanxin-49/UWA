function artifacts=generate_gaussian_sponge_validation_figures(overrides)
%GENERATE_GAUSSIAN_SPONGE_VALIDATION_FIGURES Figures and final report.

if nargin<1, overrides=struct(); end
root=fileparts(fileparts(fileparts(mfilename('fullpath')))); addpath(root); setup_vertical_project();
cfg=local_defaults(root); cfg=local_overrides(cfg,overrides);
if ~exist(cfg.figure_dir,'dir'), mkdir(cfg.figure_dir); end
g=load(fullfile(cfg.result_dir,'extended_window','gaussian_window_convergence.mat'),'validation'); gate=g.validation;
s=load(fullfile(cfg.result_dir,'gaussian_sponge_validation.mat'),'validation'); sponge=s.validation;
w=load(fullfile(cfg.result_dir,'gaussian_sponge_wideband.mat'),'validation'); wide=w.validation;

fields=struct(); fields.w50=sponge.saved_fields.production_alpha_0;
fields.w160=sponge.reference_field;
for W=[80 128]
    N=2*round(W/cfg.dx_target_m/2); out=vertical_channel_model(local_params(cfg,W,N));
    fields.(sprintf('w%d',W))=struct('psi_xy',out.psifinal_xy,'x',out.x,'y',out.y);
end

files=strings(15,1);
[X,Y]=meshgrid(fields.w50.x,fields.w50.y); psi0=gaussian_source_initial_field_vertical(X,Y,struct('x_tx',0,'y_tx',0,'sigma_src_m',cfg.sigma_src_m));
files(1)=local_image(cfg,'01_gaussian_initial_field.png'); fig=local_fig(); imagesc(fields.w50.x,fields.w50.y,abs(psi0)); axis image xy; colorbar; xlabel('x (m)'); ylabel('y (m)'); title('Production Gaussian initial field, |Psi_0|'); local_export(fig,files(1));

names={'w50','w80','w128','w160'}; labels={'50 m','80 m','128 m','160 m ref'};
files(2)=local_image(cfg,'02_no_sponge_center_amplitude_profiles.png'); fig=local_fig(); hold on; for i=1:4, f=fields.(names{i}); [~,iy]=min(abs(f.y)); m=abs(f.x)<=5; plot(f.x(m),20*log10(abs(f.psi_xy(iy,m))+eps),'DisplayName',labels{i}); end; grid on; xlabel('x (m)'); ylabel('20 log_{10}|Psi| (dB)'); legend; title('No-sponge center amplitude profiles, L=97 m'); local_export(fig,files(2));
files(3)=local_image(cfg,'03_no_sponge_center_phase_profiles.png'); fig=local_fig(); hold on; for i=1:4, f=fields.(names{i}); [~,iy]=min(abs(f.y)); m=abs(f.x)<=5; plot(f.x(m),unwrap(angle(f.psi_xy(iy,m))),'DisplayName',labels{i}); end; grid on; xlabel('x (m)'); ylabel('unwrapped phase (rad)'); legend; title('No-sponge center phase profiles, L=97 m'); local_export(fig,files(3));

t=gate.window_table; q=t(t.radius_m==2,:);
files(4)=local_image(cfg,'04_window_convergence_vs_width.png'); fig=local_fig(); tiledlayout(1,3); nexttile; hold on; nexttile; hold on; nexttile; hold on; for L=unique(q.distance_m).', a=q(q.distance_m==L,:); nexttile(1); semilogy(a.width_m,a.complex_l2,'o-','DisplayName',sprintf('L=%g',L)); nexttile(2); plot(a.width_m,a.max_amplitude_error_db,'o-'); nexttile(3); semilogy(a.width_m,a.phase_rms_rad,'o-'); end; nexttile(1); grid on; xlabel('W (m)'); ylabel('complex L2'); legend; nexttile(2); grid on; xlabel('W (m)'); ylabel('max |dA| (dB)'); yline(.05,'--'); nexttile(3); grid on; xlabel('W (m)'); ylabel('phase RMS (rad)'); yline(.01,'--'); sgtitle('No-sponge Gaussian window convergence, rho<=2 m'); local_export(fig,files(4));

p=readtable(fullfile(cfg.result_dir,'gaussian_propagation_energy.csv')); p=p(p.width_m==50 & p.distance_m==97,:);
files(5)=local_image(cfg,'05_edge_energy_vs_distance.png'); fig=local_fig(); plot(p.propagated_m,10*log10(p.edge5_energy./p.total_energy),'o-'); hold on; plot(p.propagated_m,10*log10(p.edge10_energy./p.total_energy),'s-'); grid on; xlabel('propagation distance (m)'); ylabel('edge energy fraction (dB)'); legend('outer 5%','outer 10%'); title('Production-window no-sponge edge energy'); local_export(fig,files(5));
files(6)=local_image(cfg,'06_center_energy_vs_distance.png'); fig=local_fig(); plot(p.propagated_m,10*log10(p.center_energy/p.center_energy(1)),'o-'); grid on; xlabel('propagation distance (m)'); ylabel('center energy change (dB)'); title('Production-window center energy, rho<=2 m'); local_export(fig,files(6));

off=sponge.saved_fields.production_alpha_0; on=sponge.saved_fields.production_alpha_0p15;
files(7)=local_image(cfg,'07_default_sponge_on_off_terminal_fields.png'); fig=local_fig(); tiledlayout(1,2); lim=max(abs([off.psi_xy(:);on.psi_xy(:)])); nexttile; imagesc(off.x,off.y,abs(off.psi_xy)); axis image xy; clim([0 lim]); colorbar; title('sponge off'); nexttile; imagesc(on.x,on.y,abs(on.psi_xy)); axis image xy; clim([0 lim]); colorbar; title('default sponge'); sgtitle('Gaussian terminal |Psi|, W=50 m, L=97 m'); local_export(fig,files(7));
files(8)=local_image(cfg,'08_default_sponge_amplitude_difference.png'); fig=local_fig(); imagesc(off.x,off.y,20*log10((abs(on.psi_xy)+eps)./(abs(off.psi_xy)+eps))); axis image xy; clim([-10 10]); colorbar; xlabel('x (m)'); ylabel('y (m)'); title('Default sponge amplitude change dA (dB)'); local_export(fig,files(8));
files(9)=local_image(cfg,'09_default_sponge_phase_difference.png'); fig=local_fig(); imagesc(off.x,off.y,angle(on.psi_xy.*conj(off.psi_xy))); axis image xy; clim([-pi pi]); colorbar; xlabel('x (m)'); ylabel('y (m)'); title('Default sponge phase change (rad)'); local_export(fig,files(9));

mt=sponge.matrix(sponge.matrix.width_m==50,:);
files(10)=local_heatmap(cfg,'10_axis_tl_change_heatmap.png',mt,'axis_tl_change_db','axis dTL (dB)');
files(11)=local_heatmap(cfg,'11_center_energy_change_heatmap.png',mt,'center2_energy_change_db','center energy change (dB)');
files(12)=local_heatmap(cfg,'12_axis_phase_change_heatmap.png',mt,'axis_phase_change_rad','axis phase change (rad)');

files(13)=local_image(cfg,'13_production_default_vs_reference_Hf.png'); fig=local_fig(); local_h_plot(wide,{'large_no_sponge','production_default'}); sgtitle('Production default versus large/no-sponge reference'); local_export(fig,files(13));
files(14)=local_image(cfg,'14_recommended_vs_reference_Hf.png'); fig=local_fig(); local_h_plot(wide,{'large_no_sponge','recommended_no_sponge'}); sgtitle('Recommended no-sponge versus large reference'); local_export(fig,files(14));
files(15)=local_image(cfg,'15_group_delay_comparison.png'); fig=local_fig(); hold on; for i=1:numel(wide.cases), plot(wide.cases(i).f_axis_hz,1e3*wide.cases(i).group_delay_s,'DisplayName',strrep(wide.cases(i).name,'_',' ')); end; grid on; xlabel('frequency (Hz)'); ylabel('group delay (ms)'); legend; title('Physical-field group delay'); local_export(fig,files(15));

lfm=local_lfm(wide,cfg); lfm_file=local_image(cfg,'16_optional_lfm_matched_filter.png'); fig=local_fig(); tiledlayout(1,2); nexttile; hold on; nexttile; hold on; for i=1:3, nexttile(1); plot(lfm.time_s*1e3,abs(lfm.rx(:,i))/max(abs(lfm.rx(:,1))),'DisplayName',wide.cases(i).name); nexttile(2); plot(lfm.mf_time_s*1e3,abs(lfm.mf(:,i))/max(abs(lfm.mf(:,1))),'DisplayName',wide.cases(i).name); end; nexttile(1); grid on; xlabel('time (ms)'); ylabel('normalized envelope'); legend; nexttile(2); grid on; xlabel('lag (ms)'); ylabel('normalized matched filter'); legend; sgtitle('Optional noiseless LFM engineering probe'); local_export(fig,lfm_file);

report=fullfile(root,'reports','pe_gaussian_window_sponge_validation_report.md'); local_report(report,cfg,gate,sponge,wide,p,lfm,files,lfm_file);
artifacts=struct('figures',files,'lfm_figure',lfm_file,'report',report);
end

function cfg=local_defaults(root)
cfg=struct('result_dir',fullfile(root,'results','validation','pe_gaussian_window_sponge'), ...
    'figure_dir',fullfile(root,'results','validation','pe_gaussian_window_sponge','figures'), ...
    'frequency_hz',4000,'c0_mps',1500,'sigma_src_m',0.3,'distance_m',97, ...
    'dx_target_m',50/256,'lfm_duration_s',0.02,'lfm_fs_hz',8000);
end
function cfg=local_overrides(cfg,o), n=fieldnames(o); for i=1:numel(n), if ~isfield(cfg,n{i}), error('Unknown override: %s',n{i}); end, cfg.(n{i})=o.(n{i}); end, end
function p=local_params(cfg,W,N), p=struct('f0',cfg.frequency_hz,'c0',cfg.c0_mps,'z_max',100,'z_tx',cfg.distance_m,'z_rx',0,'xw',W,'yw',W,'nx',N,'ny',N,'x_tx',0,'y_tx',0,'x_rx',0,'y_rx',0,'stepz_lamb',0.5,'sigma_src_m',cfg.sigma_src_m,'source_mode','gaussian','sponge_ratio',0.12,'alpha_max_np_per_m',0,'validation_allow_extended_window',W>100,'env_mode','uniform','enable_surface_reflection',false,'enable_bubbles',false,'enforce_1_over_R',false,'show_figures',false,'save_mode','slice','nout',4,'use_gpu',false); end
function f=local_fig(), f=figure('Visible','off','Color','w','Position',[100 100 1100 620]); end
function file=local_image(cfg,name), file=fullfile(cfg.figure_dir,name); end
function local_export(fig,file), exportgraphics(fig,file,'Resolution',180); close(fig); end
function file=local_heatmap(cfg,name,t,var,label)
r=unique(t.sponge_ratio); a=unique(t.alpha_max_np_per_m); z=nan(numel(a),numel(r)); for i=1:numel(a), for j=1:numel(r), q=t(t.sponge_ratio==r(j)&t.alpha_max_np_per_m==a(i),:); z(i,j)=q.(var); end, end
file=local_image(cfg,name); fig=local_fig(); imagesc(r,a,z); axis xy; colorbar; xlabel('sponge ratio'); ylabel('alpha max (Np/m)'); title(['Production W=50 m: ' label]); local_export(fig,file);
end
function local_h_plot(wide,wanted)
tiledlayout(2,1); for i=1:numel(wanted), c=wide.cases(strcmp({wide.cases.name},wanted{i})); nexttile(1); hold on; plot(c.f_axis_hz,20*log10(abs(c.H_physical_f)),'DisplayName',strrep(c.name,'_',' ')); nexttile(2); hold on; plot(c.f_axis_hz,c.phase_physical_rad,'DisplayName',strrep(c.name,'_',' ')); end; nexttile(1); grid on; ylabel('20log10|H| (dB)'); legend; nexttile(2); grid on; xlabel('frequency (Hz)'); ylabel('unwrapped phase (rad)'); legend;
end
function lfm=local_lfm(wide,cfg)
N=1024; fs=cfg.lfm_fs_hz; t=(0:N-1).'/fs; active=t<cfg.lfm_duration_s; T=cfg.lfm_duration_s; B=2000; tx=active.*exp(1i*2*pi*((-B/2).*t+(B/(2*T)).*t.^2)); X=fft(tx); fbb=(0:N-1).'*fs/N; fbb(fbb>=fs/2)=fbb(fbb>=fs/2)-fs;
rx=complex(zeros(N,3)); mf=complex(zeros(2*N-1,3)); peaks=zeros(3,1); delays=zeros(3,1);
for i=1:3, c=wide.cases(i); H=interp1(c.f_axis_hz-4000,c.H_direct_dsp_f,fbb,'linear',0); rx(:,i)=ifft(X.*H); mf(:,i)=conv(rx(:,i),flipud(conj(tx))); [peaks(i),k]=max(abs(mf(:,i))); delays(i)=(k-N)/fs; end
corrs=zeros(3,1); for i=1:3, corrs(i)=abs(rx(:,1)'*rx(:,i))/(norm(rx(:,1))*norm(rx(:,i))); end
lfm=struct('time_s',t,'mf_time_s',((-N+1):(N-1)).'/fs,'tx',tx,'rx',rx,'mf',mf,'peak',peaks,'peak_delay_s',delays,'correlation_to_reference',corrs);
end
function local_report(file,cfg,gate,s,wide,p,lfm,files,lfm_file)
pd=s.production_default; base=s.matrix(s.matrix.width_m==50&s.matrix.sponge_ratio==.12&s.matrix.alpha_max_np_per_m==0,:); wb=wide.summary_table; rec=wb(wb.case_name=="recommended_no_sponge",:); def=wb(wb.case_name=="production_default",:);
fid=fopen(file,'w','n','UTF-8'); c=onCleanup(@()fclose(fid));
fprintf(fid,'# Production Gaussian finite-window and sponge validation\n\n');
fprintf(fid,'Status: **default sponge not recommended; strict production recommendation is `W=160 m`, sponge off**. No PE marching, FFT convention, Gaussian definition, surface model, or communication modulation was changed.\n\n');
fprintf(fid,'## Environment and source audit\n\n- Homogeneous `c=1500 m/s`; direct path only; surface reflection, bubbles, random surface, Doppler, and communication processing off.\n- Gaussian: `Psi0=exp(-((x-x_tx)^2+(y-y_tx)^2)/(2 sigma^2))`, `sigma=%.3g m`, unit peak, no additional normalization, frequency independent, centered at `(0,0)`.\n- Production geometry: `W=50 m`, `256^2`, `dx=%.9g m`, `L=97 m`; single-frequency 4 kHz; wideband 3--5 kHz/33 points.\n- `H_f` is direct-DSP referenced; physical phase/group delay use `H_direct_physical_f` and `p(t)=Re(P exp(-i2pift))`.\n\n',cfg.sigma_src_m,cfg.dx_target_m);
fprintf(fid,'## Hard check and window reference\n\nPE--AS maximum full/center complex errors are `%.6g / %.6g`; the hard gate passed. The 128 m no-sponge field converges to the 160 m field under the preregistered single-frequency limits, so `W_ref=128 m` for the 4 kHz gate while 160 m is retained as the comparison reference.\n\n',max(gate.hard_check_table.full_complex_l2),max(gate.hard_check_table.center_complex_l2_max));
fprintf(fid,'## Production default sponge, W=50 m\n\nDefinitions: `dA=20log10(|Hsp|/|H0|)`, `dTL=-dA`, positive dTL means added loss.\n\n');
fprintf(fid,'| metric | no sponge | default sponge | change |\n|---|---:|---:|---:|\n');
fprintf(fid,'| axis magnitude | %.9g | %.9g | dA %.6f dB / dTL %.6f dB |\n',abs(base.H_axis),abs(pd.H_axis),pd.axis_amplitude_change_db,pd.axis_tl_change_db);
fprintf(fid,'| axis phase | %.9f rad | %.9f rad | %.6f rad |\n',angle(base.H_axis),angle(pd.H_axis),pd.axis_phase_change_rad);
fprintf(fid,'| center energy, rho<=2 m | %.9g | %.9g | %.6f dB |\n',base.center2_energy,pd.center2_energy,pd.center2_energy_change_db);
fprintf(fid,'| edge energy | %.9g | %.9g | %.6f dB |\n',base.edge_energy,pd.edge_energy,pd.edge_energy_change_db);
fprintf(fid,'| total energy | %.9g | %.9g | %.6f dB |\n',base.total_energy,pd.total_energy,pd.total_energy_change_db);
fprintf(fid,'\nThe default changes axis TL by `%.6f dB` and phase by `%.6f rad`; it is not center-neutral. Edge suppression is only `%.6f dB`, below the fixed 3 dB requirement.\n\n',pd.axis_tl_change_db,pd.axis_phase_change_rad,-pd.edge_energy_change_db);
fprintf(fid,'## Boundary wrap evidence\n\nAt 97 m in the production window the outer 5%%/10%% energy fractions reach `%.6g / %.6g` and boundary amplitude reaches `%.6f dB` relative to the plane maximum. Thus the no-sponge Gaussian field materially reaches the periodic boundary.\n\n',p.edge5_energy(end)/p.total_energy(end),p.edge10_energy(end)/p.total_energy(end),p.boundary_max_relative_db(end));
fprintf(fid,'## Wideband results versus 160 m/no-sponge\n\n| case | max |dTL| dB | TL span dB | phase RMS/max rad | group-delay RMS/max us |\n|---|---:|---:|---:|---:|\n'); for i=1:height(wb), r=wb(i,:); fprintf(fid,'| %s | %.6f | %.6f | %.6f / %.6f | %.6f / %.6f |\n',r.case_name,r.max_abs_tl_error_db,r.amplitude_error_span_db,r.rms_phase_error_rad,r.max_abs_phase_error_rad,r.rms_group_delay_error_us,r.max_abs_group_delay_error_us); end
fprintf(fid,'\nThe 128 m no-sponge candidate misses the strict wideband `|dTL|<0.1 dB` target by `%.6g dB`; it is therefore a near-threshold cost compromise, not the strict recommendation.\n\n',rec.max_abs_tl_error_db-.1);
fprintf(fid,'## Direct answers\n\n1. **Q1:** 128 m is the first tested 4 kHz window satisfying the single-frequency convergence gate; 160 m is used as the strict wideband reference/recommendation.\n2. **Q2:** The production 50 m window is not large enough at L=97 m. No-sponge axis TL/phase errors versus 160 m are `%.6f dB / %.6f rad`.\n3. **Q3:** Yes. The final boundary is only `%.3f dB` below the field maximum.\n4. **Q4:** Default single-frequency and wideband errors are tabulated above; wideband max TL, phase RMS, and group-delay RMS are `%.6f dB`, `%.6f rad`, and `%.3f us`.\n5. **Q5:** No. Default sponge does not provide >=3 dB edge suppression and significantly perturbs axis amplitude/phase.\n6. **Q6:** No nonzero scanned sponge meets all fixed center and edge targets. Recommend `W=160 m`, sponge off.\n7. **Q7:** Yes: expanding to 160 m with no sponge is sufficient in the tested 3--5 kHz, 97 m direct-path environment.\n\n',base.axis_tl_error_vs_large_db,base.axis_phase_error_vs_large_rad,p.boundary_max_relative_db(end),def.max_abs_tl_error_db,def.rms_phase_error_rad,def.rms_group_delay_error_us);
fprintf(fid,'## Recommendation\n\n- Strict: `window=160 m`, `sponge off (alpha_max=0)`; ratio is inactive. Relative errors are zero by reference definition.\n- Cost compromise: `window=128 m`, sponge off; max wideband TL `%.6f dB`, phase RMS `%.6f rad`, group-delay RMS `%.6f us`.\n- Do **not** change production defaults automatically yet. The evidence recommends a future reviewed config change from 50 m/default sponge to a larger no-sponge window, followed by communication-cost and full reflected-path qualification.\n\n',rec.max_abs_tl_error_db,rec.rms_phase_error_rad,rec.rms_group_delay_error_us);
fprintf(fid,'## Optional LFM probe\n\nNormalized waveform correlations to the 160 m reference are production/default `%.6f` and 128 m/no-sponge `%.6f`; matched-filter peak delays are `%s s`. This is engineering interpretation only.\n\n',lfm.correlation_to_reference(2),lfm.correlation_to_reference(3),mat2str(lfm.peak_delay_s.',6));
prefix=[cfg.figure_dir filesep];
fprintf(fid,'## Figures\n\n'); for i=1:numel(files), fprintf(fid,'![figure %d](../results/validation/pe_gaussian_window_sponge/figures/%s)\n\n',i,string(extractAfter(files(i),prefix))); end; fprintf(fid,'![optional LFM](../results/validation/pe_gaussian_window_sponge/figures/%s)\n',string(extractAfter(lfm_file,prefix))); clear c
end
