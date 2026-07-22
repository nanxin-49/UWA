function plot_pe_propagation_atlas_vertical(a,out_dir)
%PLOT_PE_PROPAGATION_ATLAS_VERTICAL Render the numbered PE visual atlas.

arguments
    a (1,1) struct
    out_dir (1,:) char
end
if ~isfolder(out_dir), mkdir(out_dir); end
local_geometry(a,out_dir);
local_xz(a,out_dir);
local_surface_fields(a,out_dir);
local_spectra(a,out_dir);
local_branches(a,out_dir);
local_receiver_planes(a,out_dir);
local_profiles(a,out_dir);
local_adjoint(a,out_dir);
local_frequency_response(a,out_dir);
local_cp(a,out_dir);
local_pdp_distribution(a,out_dir);
local_wind_contrast(a,out_dir);
local_lfm(a,out_dir);
local_animation(a,out_dir);
end

function local_geometry(a,out_dir)
fig = local_figure([80 80 1250 650]); ax = axes(fig); hold(ax,'on');
plot(ax,[-30 30],[0 0],'Color',[0.1 0.45 0.75],'LineWidth',3);
patch(ax,[-25 25 25 -25],[0 0 a.cfg.z_tx a.cfg.z_tx],[0.85 0.93 0.98], ...
    'FaceAlpha',0.35,'EdgeColor','none');
plot(ax,[-a.cfg.xw/2 a.cfg.xw/2 a.cfg.xw/2 -a.cfg.xw/2 -a.cfg.xw/2], ...
    [0 0 a.cfg.z_tx a.cfg.z_tx 0],'k-','LineWidth',1.2);
s = a.cfg.sponge_ratio*a.cfg.xw;
patch(ax,[-a.cfg.xw/2 -a.cfg.xw/2+s -a.cfg.xw/2+s -a.cfg.xw/2], ...
    [0 0 a.cfg.z_tx a.cfg.z_tx],[0.85 0.55 0.25],'FaceAlpha',0.22,'EdgeColor','none');
patch(ax,[a.cfg.xw/2-s a.cfg.xw/2 a.cfg.xw/2 a.cfg.xw/2-s], ...
    [0 0 a.cfg.z_tx a.cfg.z_tx],[0.85 0.55 0.25],'FaceAlpha',0.22,'EdgeColor','none');
plot(ax,a.cfg.x_tx,a.cfg.z_tx,'p','MarkerSize',14,'MarkerFaceColor',[0.85 0.2 0.15]);
plot(ax,a.cfg.x_rx,a.cfg.z_rx,'v','MarkerSize',11,'MarkerFaceColor',[0.15 0.55 0.2]);
quiver(ax,-1,a.cfg.z_tx-2,0,-a.cfg.z_tx+4,0,'LineWidth',2,'Color',[0.15 0.3 0.75], ...
    'MaxHeadSize',0.08);
quiver(ax,1,1,0,a.cfg.z_rx-1,0,'LineWidth',2,'Color',[0.75 0.25 0.15], ...
    'MaxHeadSize',0.6);
plot(ax,[0 0],[a.cfg.z_tx a.cfg.z_rx],'--','Color',[0.2 0.55 0.25],'LineWidth',1.5);
text(ax,-29,-2,'Sea surface z=0','FontWeight','bold');
text(ax,2,a.cfg.z_tx-3,'Tx','FontWeight','bold'); text(ax,2,a.cfg.z_rx+4,'Rx','FontWeight','bold');
text(ax,-12,a.cfg.z_tx/2,'Tx → surface incident PE','Rotation',90,'HorizontalAlignment','center');
text(ax,5,a.cfg.z_tx/2,'Direct PE','Rotation',90,'HorizontalAlignment','center');
text(ax,10,a.cfg.z_rx/2,'Surface → Rx reflected PE','HorizontalAlignment','left');
text(ax,-a.cfg.xw/2+0.5,a.cfg.z_tx*0.8,'sponge','Rotation',90);
text(ax,a.cfg.xw/2-0.5,a.cfg.z_tx*0.8,'sponge','Rotation',90,'HorizontalAlignment','right');
set(ax,'YDir','reverse'); axis(ax,'equal'); xlim(ax,[-32 32]); ylim(ax,[-5 a.cfg.z_tx+5]); grid(ax,'on');
xlabel(ax,'transverse coordinate (m)'); ylabel(ax,'z (m, positive downward)');
title(ax,sprintf('Uniform PE geometry: PE %d^2 / PM %d^2, nearest-grid receiver',a.cfg.nx,a.pm.nx));
local_export(fig,out_dir,'01_geometry_and_paths.png');
end

function local_xz(a,out_dir)
m = a.public.wavefield_meta;
ref = max([abs(a.public.direct_xz(:));abs(m.incident_field_slice(:));abs(m.reflected_field_slice(:));eps]);
fig = local_figure([50 80 1550 500]); tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
ax=nexttile; local_db_image(ax,a.public.x_m,a.public.z_track_m,a.public.direct_xz,ref); title(ax,'Direct: Tx → Rx');
ax=nexttile; local_db_image(ax,m.transverse_coordinate_m,m.incident_z_m,m.incident_field_slice,ref); title(ax,'Incident: Tx → surface');
ax=nexttile; local_db_image(ax,m.transverse_coordinate_m,m.reflected_z_m,m.reflected_field_slice,ref); title(ax,'Reflected: surface → Rx');
sgtitle(fig,sprintf('PE complex-envelope propagation, %.0f Hz, common reference',a.reference_frequency_hz));
local_export(fig,out_dir,'02_pe_xz_propagation.png');
end

function local_surface_fields(a,out_dir)
s = a.surface; ref = max(abs(s.incident_xy(:)));
fig = local_figure([40 40 1500 1050]); tiledlayout(fig,3,3,'TileSpacing','compact','Padding','compact');
ax=nexttile; imagesc(ax,a.x_m,a.y_m,s.eta_pe_m); axis(ax,'image'); set(ax,'YDir','normal'); colorbar(ax); title(ax,'Raw-PM elevation eta (m)');
ax=nexttile; local_xy_db(ax,a,s.incident_xy,ref); title(ax,'Incident magnitude');
ax=nexttile; local_xy_phase(ax,a,s.incident_xy,ref); title(ax,'Incident phase');
ax=nexttile; local_xy_db(ax,a,s.coherent_xy,ref); title(ax,'Coherent magnitude');
ax=nexttile; local_xy_phase(ax,a,s.coherent_xy,ref); title(ax,'Coherent phase');
ax=nexttile; local_xy_db(ax,a,s.joint_xy,ref); title(ax,'Joint reflected magnitude');
ax=nexttile; local_xy_phase(ax,a,s.joint_xy,ref); title(ax,'Joint reflected phase');
ax=nexttile; local_xy_db(ax,a,s.joint_scatter_xy,ref); title(ax,'Joint scatter magnitude');
ax=nexttile; local_xy_phase(ax,a,s.joint_scatter_xy,ref); title(ax,'Joint scatter phase');
sgtitle(fig,sprintf('Sea-surface fields, U=5 m/s, %.0f Hz; phase masked below -40 dB',a.reference_frequency_hz));
local_export(fig,out_dir,'03_surface_boundary_fields.png');
end

function local_spectra(a,out_dir)
fields = {a.surface.incident_xy,a.surface.joint_xy,a.surface.joint_scatter_xy};
names = {'Incident','Joint reflected','Joint scatter'};
[KX,KY] = local_kgrid(a.x_m,a.y_m); powers=cell(1,3); common=eps;
for ii=1:3, powers{ii}=abs(fftshift(fft2(fields{ii}))).^2; common=max(common,max(powers{ii}(:))); end
fig=local_figure([60 80 1500 500]); tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
for ii=1:3
    ax=nexttile; imagesc(ax,KX(1,:),KY(:,1),10*log10(max(powers{ii},eps)/common));
    axis(ax,'image'); set(ax,'YDir','normal'); clim(ax,[-60 0]); colorbar(ax); colormap(ax,turbo(256));
    st=local_spectrum_stats(powers{ii},KX,KY);
    title(ax,sprintf('%s\nRMS %.3g, r90 %.3g, high-k %.3g',names{ii},st.rms_k,st.r90,st.high_fraction));
    xlabel(ax,'K_x (rad/m)'); ylabel(ax,'K_y (rad/m)');
end
sgtitle(fig,'Angular-spectrum redistribution, common power reference (dB)');
local_export(fig,out_dir,'04_surface_angular_spectra.png');
end

function local_branches(a,out_dir)
fields={a.surface.flat_xy,a.surface.explicit_xy,a.surface.joint_xy,a.surface.ssa1_xy};
names={'Flat pressure-release','Explicit Kirchhoff','Joint k-stat (illustrative)','SSA1 geometry reference'};
ref=max(abs(a.surface.incident_xy(:)));
fig=local_figure([40 40 1600 800]); tiledlayout(fig,2,4,'TileSpacing','compact','Padding','compact');
for ii=1:4, ax=nexttile; local_xy_db(ax,a,fields{ii},ref); title(ax,[names{ii} ' magnitude']); end
for ii=1:4, ax=nexttile; local_xy_phase(ax,a,fields{ii},ref); title(ax,[names{ii} ' phase']); end
sgtitle(fig,'Physical surface-boundary branches; explicit/joint single fields are not pointwise validation');
local_export(fig,out_dir,'05_physical_boundary_branches.png');
end

function local_receiver_planes(a,out_dir)
r=a.receiver; fields={r.direct_xy,r.coherent_xy,r.joint_scatter_xy,r.joint_xy,r.total_joint_xy};
names={'Direct','Coherent reflected','Scattered reflected','Total reflected','Direct + reflected'};
ref=max(cellfun(@(x)max(abs(x(:))),fields));
fig=local_figure([20 20 1750 760]); tiledlayout(fig,2,5,'TileSpacing','compact','Padding','compact');
for ii=1:5
    ax=nexttile; local_xy_db(ax,a,fields{ii},ref); hold(ax,'on');
    plot(ax,a.cfg.x_rx,a.cfg.y_rx,'wo','MarkerFaceColor','k','MarkerSize',5); title(ax,[names{ii} ' magnitude']);
end
for ii=1:5
    ax=nexttile; local_xy_phase(ax,a,fields{ii},ref); hold(ax,'on');
    plot(ax,a.cfg.x_rx,a.cfg.y_rx,'wo','MarkerFaceColor','k','MarkerSize',5); title(ax,[names{ii} ' phase']);
end
sgtitle(fig,sprintf('Receiver-depth planes at z=%.3g m; marker is nearest-grid sample',a.cfg.z_rx));
local_export(fig,out_dir,'06_receiver_plane_components.png');
end

function local_profiles(a,out_dir)
r=a.receiver; fields={r.direct_xy,r.coherent_xy,r.joint_scatter_xy,r.joint_xy,r.explicit_xy,r.ssa1_xy};
names={'direct','coherent','joint scatter','joint reflected','explicit reflected','SSA1 reflected'};
iy=floor(numel(a.y_m)/2)+1; ix=floor(numel(a.x_m)/2)+1;
ref=max(cellfun(@(x)max(abs(x(:))),fields));
fig=local_figure([80 80 1300 760]); tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
ax=nexttile; hold(ax,'on'); for ii=1:numel(fields), plot(ax,a.x_m,20*log10(max(abs(fields{ii}(iy,:)),eps)/ref),'LineWidth',1.1); end
grid(ax,'on'); ylim(ax,[-60 1]); ylabel(ax,'magnitude (dB)'); xlabel(ax,'x (m)'); title(ax,'Receiver-plane x profile'); legend(ax,names,'Location','eastoutside');
ax=nexttile; hold(ax,'on'); for ii=1:numel(fields), plot(ax,a.y_m,20*log10(max(abs(fields{ii}(:,ix)),eps)/ref),'LineWidth',1.1); end
grid(ax,'on'); ylim(ax,[-60 1]); ylabel(ax,'magnitude (dB)'); xlabel(ax,'y (m)'); title(ax,'Receiver-plane y profile'); legend(ax,names,'Location','eastoutside');
local_export(fig,out_dir,'07_receiver_transverse_profiles.png');
end

function local_adjoint(a,out_dir)
f=a.cfg.f_axis_hz(:)/1000; hc=a.projection.H_cached_scatter_f(:); hp=a.projection.H_projected_scatter_f(:);
err=abs(hc-hp)./max(abs(hc),eps); q=a.projection.q_xy; w=a.projection.a_xy;
hrc=a.projection.H_cached_scatter_reduced_f(:); hrp=a.projection.H_projected_scatter_reduced_f(:);
err_reduced=abs(hrc-hrp)./max(abs(hrc),eps);
qref=max(abs(q(:))); wref=max(abs(w(:)));
fig=local_figure([30 30 1500 850]); tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
ax=nexttile; plot(ax,f,abs(hc),'o-','LineWidth',1.2); hold(ax,'on'); plot(ax,f,abs(hp),'--','LineWidth',1.2); grid(ax,'on'); xlabel(ax,'frequency (kHz)'); ylabel(ax,'|H_{sca}|'); legend(ax,'cached PE','adjoint projection'); title(ax,'Same-input receiver response');
ax=nexttile; semilogy(ax,f,max(err_reduced,eps),'o-'); hold(ax,'on'); semilogy(ax,f,max(err,eps),'--'); grid(ax,'on'); xlabel(ax,'frequency (kHz)'); ylabel(ax,'relative error'); legend(ax,'reduced','direct-DSP'); title(ax,sprintf('max %.3g / %.3g',max(err_reduced),max(err)));
ax=nexttile; local_xy_db(ax,a,q,qref); title(ax,'|q| sensitivity kernel (dB)');
ax=nexttile; local_xy_phase(ax,a,q,qref); title(ax,'phase(q), not physical back-propagated pressure');
ax=nexttile; local_xy_db(ax,a,w,wref); title(ax,'|a| = |conj(psi_{inc})q| (dB)');
ax=nexttile; local_xy_phase(ax,a,w,wref); title(ax,'phase(a)');
sgtitle(fig,'Cached forward PE versus exact discrete-adjoint receiver projection');
local_export(fig,out_dir,'08_cached_vs_adjoint.png');
end

function local_frequency_response(a,out_dir)
f=a.cfg.f_axis_hz(:); he=a.projection.H_explicit_reflect_f(:); hj=a.projection.H_joint_reflect_f(:);
he_red=a.projection.H_explicit_reflect_reduced_f(:); hj_red=a.projection.H_joint_reflect_reduced_f(:);
coh=a.projection.H_ref_coh_f(:); rms_sca=sqrt(max(real(diag(a.stats.analytic_stats.C_H)),0));
tau_e=local_group_delay(f,he); tau_j=local_group_delay(f,hj);
db_floor=-120;
he_db=max(20*log10(max(abs(he),eps)),db_floor);
hj_db=max(20*log10(max(abs(hj),eps)),db_floor);
coh_db=max(20*log10(max(abs(coh),eps)),db_floor);
rms_db=max(20*log10(max(rms_sca,eps)),db_floor);
fig=local_figure([70 70 1400 780]); tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
ax=nexttile; plot(ax,f/1000,he_db,'LineWidth',1.2); hold(ax,'on'); plot(ax,f/1000,hj_db,'--','LineWidth',1.2); plot(ax,f/1000,coh_db,':','LineWidth',1.2); plot(ax,f/1000,rms_db,'-.','LineWidth',1.2); grid(ax,'on'); ylim(ax,[db_floor 0]); ylabel(ax,'magnitude (dB)'); legend(ax,'explicit realization','joint realization','coherent mean','analytic scatter RMS','Location','eastoutside'); title(ax,sprintf('Reflected-only response (display floor %g dB)',db_floor));
theory=-2*pi*f*a.phase_reference_meta.relative_delay_s;
ax=nexttile; plot(ax,f/1000,unwrap(angle(he_red)),'LineWidth',1.0); hold(ax,'on'); plot(ax,f/1000,unwrap(angle(hj_red)),'--','LineWidth',1.0); plot(ax,f/1000,unwrap(angle(he)),'LineWidth',1.4); plot(ax,f/1000,unwrap(angle(hj)),'--','LineWidth',1.4); plot(ax,f/1000,theory-theory(1)+unwrap(angle(he_red(1))),':k','LineWidth',1.1); grid(ax,'on'); ylabel(ax,'unwrapped phase (rad)'); legend(ax,'explicit reduced','joint reduced','explicit direct-DSP','joint direct-DSP','nominal -2\pi f\Delta\tau_0','Location','eastoutside');
ax=nexttile; plot(ax,f(2:end)/1000,tau_e*1e3,'LineWidth',1.2); hold(ax,'on'); plot(ax,f(2:end)/1000,tau_j*1e3,'--','LineWidth',1.2); grid(ax,'on'); ylabel(ax,'group delay (ms)'); xlabel(ax,'frequency (kHz)'); legend(ax,'explicit','joint','Location','eastoutside');
local_export(fig,out_dir,'09_receiver_frequency_response.png');
end

function local_cp(a,out_dir)
s=a.stats; C=s.analytic_stats.C_H; P=s.analytic_stats.P_H; Cs=s.comparison.sample_stats.C; Ps=s.comparison.sample_stats.P;
scale=max(abs(C(:))); mats={C,Cs,Cs-C,P,Ps,Ps-P}; names={'analytic C','sample C','sample-analytic C','analytic P','sample P','sample-analytic P'};
fig=local_figure([30 30 1450 880]); tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
for ii=1:6
    ax=nexttile; imagesc(ax,20*log10(max(abs(mats{ii}),eps)/scale)); axis(ax,'image'); set(ax,'YDir','normal'); clim(ax,[-100 0]); colorbar(ax); colormap(ax,turbo(256)); title(ax,names{ii}); xlabel(ax,'f_j index'); ylabel(ax,'f_i index');
end
sgtitle(fig,sprintf('Receiver C/P on common |C| scale; C/floor %.3f, P(C-norm)/floor %.3f', ...
    s.comparison.C_to_split_floor_ratio,s.comparison.P_Cnorm_to_split_floor_ratio));
local_export(fig,out_dir,'10_receiver_cp_heatmaps.png');
end

function local_pdp_distribution(a,out_dir)
m=a.model5; H=a.distribution.H_scatter_fm; [~,idx]=min(abs(m.frequency_axis-6000)); h=H(idx,:);
fig=local_figure([40 40 1550 900]); tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
pdp_a=a.stats.comparison.pdp_analytic(:); pdp_s=a.stats.comparison.pdp_sample(:);
delay_bin=(0:numel(pdp_a)-1).';
ax=nexttile; plot(ax,delay_bin,pdp_a/max(pdp_a),'LineWidth',1.3); hold(ax,'on'); plot(ax,delay_bin,pdp_s/max(pdp_s),'--','LineWidth',1.2); grid(ax,'on'); xlabel(ax,'delay-bin index'); ylabel(ax,'normalized power'); legend(ax,'analytic FFT C/P','projected realizations'); title(ax,sprintf('Reflected-scatter PDP, correlation %.4f',a.stats.comparison.pdp_correlation));
ax=nexttile; ev=max(real(m.eigenvalues(:)),0); yyaxis(ax,'left'); semilogy(ax,max(ev/ max(ev),eps),'o-'); ylabel(ax,'normalized eigenvalue'); yyaxis(ax,'right'); plot(ax,cumsum(ev)/sum(ev),'LineWidth',1.4); ylabel(ax,'cumulative variance'); grid(ax,'on'); xlabel(ax,'mode'); title(ax,sprintf('Covariance modes: selected %d',m.rank_selected));
ax=nexttile; scatter(ax,real(h),imag(h),8,'filled','MarkerFaceAlpha',0.25); axis(ax,'equal'); grid(ax,'on'); xlabel(ax,'Re H_{sca}'); ylabel(ax,'Im H_{sca}'); title(ax,'6 kHz scatter IQ, conditional full rank');
ax=nexttile; histogram(ax,abs(h),40,'Normalization','pdf'); grid(ax,'on'); xlabel(ax,'|H_{sca}|'); ylabel(ax,'density'); title(ax,sprintf('Amplitude distribution, P/C %.3g',m.properness_ratio));
q=a.phase_audit; ax=nexttile; plot(ax,q.delay_axis_s*1e3,abs(q.direct_cir)/max(abs(q.direct_cir)),'LineWidth',1.2); hold(ax,'on'); plot(ax,q.delay_axis_s*1e3,abs(q.reflect_cir)/max(abs(q.reflect_cir)),'--','LineWidth',1.2); grid(ax,'on'); xlabel(ax,'delay (ms)'); ylabel(ax,'normalized |h|'); legend(ax,'direct','reflected'); title(ax,sprintf('F=65 two-ray: \Delta\tau_0=%.3f ms',1e3*q.relative_delay_s));
ax=nexttile; axis(ax,'off'); text(ax,0.02,0.78,sprintf('F=9 spacing: 0.5 kHz\nUnambiguous window: %.1f ms\n4 ms aliases to zero phase\n6 kHz: %.0f carrier cycles',1e3*q.f9_unambiguous_window_s,q.six_khz_cycle_count),'FontSize',12,'VerticalAlignment','top'); title(ax,'Phase-audit interpretation');
local_export(fig,out_dir,'11_pdp_eigenspectrum_distribution.png');
end

function local_wind_contrast(a,out_dir)
m5=a.model5; m8=a.model8; f=m5.frequency_axis(:)/1000;
pow5=real(diag(m5.C_scatter_f)); pow8=real(diag(m8.C_scatter_f));
R5=abs(local_corr_matrix(m5.C_scatter_f)); R8=abs(local_corr_matrix(m8.C_scatter_f));
fig=local_figure([30 30 1450 880]); tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
ax=nexttile; plot(ax,f,pow5,'LineWidth',1.3); hold(ax,'on'); plot(ax,f,pow8,'--','LineWidth',1.3); grid(ax,'on'); xlabel(ax,'frequency (kHz)'); ylabel(ax,'E|H_{sca}|^2'); legend(ax,'U=5','U=8'); title(ax,'Reflected-scatter power');
ax=nexttile; plot(ax,m5.delay_axis_s*1e3,m5.pdp/max(m5.pdp),'LineWidth',1.3); hold(ax,'on'); plot(ax,m8.delay_axis_s*1e3,m8.pdp/max(m8.pdp),'--','LineWidth',1.3); grid(ax,'on'); xlabel(ax,'delay (ms)'); ylabel(ax,'normalized PDP'); legend(ax,'U=5','U=8'); title(ax,'Validated receiver PDP');
ax=nexttile; imagesc(ax,R5); axis(ax,'image'); set(ax,'YDir','normal'); clim(ax,[0 1]); colorbar(ax); title(ax,sprintf('U=5 |correlation|, rank %d',m5.rank_selected));
ax=nexttile; imagesc(ax,R8); axis(ax,'image'); set(ax,'YDir','normal'); clim(ax,[0 1]); colorbar(ax); title(ax,sprintf('U=8 |correlation|, rank %d',m8.rank_selected));
sgtitle(fig,'Exact-node sea-state contrast; no wind-speed interpolation');
local_export(fig,out_dir,'12_u5_u8_statistical_contrast.png');
end

function local_lfm(a,out_dir)
l=a.lfm_saved; t=l.t_s(:)*1e3; tx=l.tx_bb(:);
fig=local_figure([30 30 1450 880]); tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
ax=nexttile; plot(ax,t,abs(tx)/max(abs(tx)),'k','LineWidth',1.1); hold(ax,'on'); plot(ax,t,abs(l.rx_reflect_kdomain)/max(abs(l.rx_reflect_kdomain)),'LineWidth',1.1); plot(ax,t,abs(l.rx_reflect_joint)/max(abs(l.rx_reflect_joint)),'--','LineWidth',1.1); grid(ax,'on'); xlabel(ax,'time (ms)'); ylabel(ax,'normalized envelope'); legend(ax,'TX LFM','kdomain RX','joint RX'); title(ax,'Current-run noiseless reflected-only LFM');
ax=nexttile; plot(ax,t,abs(l.rx_total_kdomain)/max(abs(l.rx_total_kdomain)),'LineWidth',1.1); hold(ax,'on'); plot(ax,t,abs(l.rx_total_joint)/max(abs(l.rx_total_joint)),'--','LineWidth',1.1); grid(ax,'on'); xlabel(ax,'time (ms)'); ylabel(ax,'normalized envelope'); legend(ax,'kdomain total','joint total'); title(ax,'Total channel shown separately');
ax=nexttile; plot(ax,a.stats.comparison.lfm_analytic,'LineWidth',1.2); hold(ax,'on'); plot(ax,a.stats.comparison.lfm_sample,'--','LineWidth',1.2); grid(ax,'on'); xlabel(ax,'sample'); ylabel(ax,'power'); legend(ax,'analytic','sample'); title(ax,sprintf('Current joint LFM power, corr %.5f',a.stats.comparison.lfm_correlation));
ax=nexttile; plot(ax,a.stats.comparison.matched_filter_analytic,'LineWidth',1.2); hold(ax,'on'); plot(ax,a.stats.comparison.matched_filter_sample,'--','LineWidth',1.2); grid(ax,'on'); xlabel(ax,'sample'); ylabel(ax,'power'); legend(ax,'analytic','sample'); title(ax,sprintf('Matched-filter power, corr %.5f',a.stats.comparison.matched_filter_correlation));
sgtitle(fig,'LFM is a channel-level linear probe: no noise, modulation, synchronization, or equalization');
local_export(fig,out_dir,'13_lfm_and_matched_filter.png');
end

function local_animation(a,out_dir)
m=a.public.wavefield_meta; file=fullfile(out_dir,'14_pe_carrier_reconstruction.mp4');
writer=VideoWriter(file,'MPEG-4'); writer.FrameRate=12; writer.Quality=95; open(writer);
inc_carrier=exp(1i*m.k0_rad_per_m*(a.cfg.z_tx-m.incident_z_m(:).'));
ref_carrier=exp(1i*m.k0_rad_per_m*m.reflected_z_m(:).');
inc=m.incident_field_slice.*inc_carrier; ref=m.reflected_field_slice.*ref_carrier;
scale=max([abs(inc(:));abs(ref(:));eps]); fig=local_figure([50 50 1400 620]);
cleanup=onCleanup(@()local_close_video(writer,fig)); %#ok<NASGU>
for frame=1:48
    phase=2*pi*(frame-1)/48; clf(fig); tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
    ax=nexttile; imagesc(ax,m.transverse_coordinate_m,m.incident_z_m,real(inc.*exp(-1i*phase)).'/scale); set(ax,'YDir','reverse'); clim(ax,[-1 1]); colorbar(ax); colormap(ax,local_redblue(256)); xlabel(ax,'x (m)'); ylabel(ax,'z (m)'); title(ax,'Incident carrier reconstruction');
    ax=nexttile; imagesc(ax,m.transverse_coordinate_m,m.reflected_z_m,real(ref.*exp(-1i*phase)).'/scale); set(ax,'YDir','reverse'); clim(ax,[-1 1]); colorbar(ax); colormap(ax,local_redblue(256)); xlabel(ax,'x (m)'); ylabel(ax,'z (m)'); title(ax,'Reflected carrier reconstruction');
    sgtitle(fig,sprintf('%.0f Hz, U=5 m/s, t/T=%.3f — visualization only, not time-domain PE',a.reference_frequency_hz,(frame-1)/48)); drawnow; writeVideo(writer,getframe(fig));
end
close(writer); close(fig); clear cleanup
end

function local_db_image(ax,x,z,field,ref)
imagesc(ax,x,z,20*log10(max(abs(field),eps)/ref).'); set(ax,'YDir','reverse'); axis(ax,'tight'); clim(ax,[-50 0]); colorbar(ax); colormap(ax,turbo(256)); xlabel(ax,'x (m)'); ylabel(ax,'z (m, positive downward)');
end

function local_xy_db(ax,a,field,ref)
imagesc(ax,a.x_m,a.y_m,20*log10(max(abs(field),eps)/max(ref,eps))); axis(ax,'image'); set(ax,'YDir','normal'); clim(ax,[-50 0]); colorbar(ax); colormap(ax,turbo(256)); xlabel(ax,'x (m)'); ylabel(ax,'y (m)');
end

function local_xy_phase(ax,a,field,ref)
db=20*log10(max(abs(field),eps)/max(ref,eps)); h=imagesc(ax,a.x_m,a.y_m,angle(field)); set(h,'AlphaData',db>=-40); axis(ax,'image'); set(ax,'YDir','normal','Color',[0.75 0.75 0.75]); clim(ax,[-pi pi]); colorbar(ax); colormap(ax,hsv(256)); xlabel(ax,'x (m)'); ylabel(ax,'y (m)');
end

function [KX,KY]=local_kgrid(x,y)
nx=numel(x); ny=numel(y); dx=mean(diff(x)); dy=mean(diff(y));
kx=(2*pi/(nx*dx))*(-nx/2:nx/2-1); ky=(2*pi/(ny*dy))*(-ny/2:ny/2-1); [KX,KY]=meshgrid(kx,ky);
end

function st=local_spectrum_stats(P,KX,KY)
w=P/max(sum(P(:)),eps); kr=sqrt(KX.^2+KY.^2); st.rms_k=sqrt(sum(w(:).*kr(:).^2));
[r,idx]=sort(kr(:)); c=cumsum(w(idx)); j=find(c>=0.9,1); st.r90=r(j); st.high_fraction=sum(w(kr>=0.5*max(kr(:))));
end

function tau=local_group_delay(f,H)
tau=-diff(unwrap(angle(H)))./(2*pi*diff(f));
end

function R=local_corr_matrix(C)
d=sqrt(max(real(diag(C)),eps)); R=C./(d*d.');
end

function fig=local_figure(position)
fig=figure('Visible','off','Color','w','Position',position);
end

function local_export(fig,out_dir,name)
exportgraphics(fig,fullfile(out_dir,name),'Resolution',180); close(fig);
end

function cmap=local_redblue(n)
x=linspace(0,1,n).'; cmap=[min(1,2*x),1-abs(2*x-1),min(1,2*(1-x))];
end

function local_close_video(writer,fig)
try
    close(writer);
catch
end
if isgraphics(fig), close(fig); end
end
