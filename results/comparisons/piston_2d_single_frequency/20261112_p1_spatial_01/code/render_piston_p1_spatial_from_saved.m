function receipt = render_piston_p1_spatial_from_saved(prior_run,output_dir)
%RENDER_PISTON_P1_SPATIAL_FROM_SAVED Evaluate saved weak-C2 solutions in 2-D.
% All three methods reuse the accepted 4 kHz source, surface and BIE density.
arguments
    prior_run (1,:) char
    output_dir (1,:) char
end
assert(~isfolder(output_dir),'Use a new spatial result directory.');
mkdir(output_dir);
diary(fullfile(output_dir,'matlab_diary.log')); cleanup=onCleanup(@()diary('off'));
t0=tic;
q=load(fullfile(prior_run,'common.mat'));
pe=load(fullfile(prior_run,'jobs','weak_pe_cont_l6','data.mat'),'result');
bie=load(fullfile(prior_run,'jobs','weak_bie_l8','data.mat'),'result');
fk=load(fullfile(prior_run,'jobs','weak_fullfk_l3','data.mat'),'result');
assert(isequal(q.cfg.widths_m,[.2 .1]) && q.cfg.frequencies_hz(1)==4000);
assert(strcmp(pe.result.stage,'weak') && strcmp(bie.result.stage,'weak') && strcmp(fk.result.stage,'weak'));
x=(-40:.5:40).'; z=(.5:.25:10).';
[X,Z]=meshgrid(x,z); xv=X(:); zv=Z(:); n=numel(xv);
fields=complex(zeros(numel(z),numel(x),3));
source=struct('frequency_hz',4000,'c_mps',q.cfg.c_mps,'rho_kgpm3',q.cfg.rho_kgpm3, ...
    'v0_mps',q.cfg.v0_mps,'total_width_m',.2,'z_source_m',100);
flat=-baffled_strip_piston_field_vertical(xv,-zv,source);
flat_xy=reshape(flat,size(X));
% The saved compact surface screen is propagated using its frozen FFT grid.
surface=pe.result.surface_inputs{1}; delta=surface.delta_screen_pressure_Pa;
W=pe.result.numeric_config.window_m; dx=pe.result.numeric_config.dx_m;
N=numel(delta); k=2*pi*4000/1500;
kx=2*pi/W*[0:N/2-1 -N/2:-1]'; kz=sqrt(complex(k^2-kx.^2,0));
spec=fft(delta); ix=round((x+W/2)/dx)+1;
for j=1:numel(z)
    value=ifft(spec.*exp(1i*kz*z(j)));
    fields(j,:,1)=flat_xy(j,:)+value(ix).';
end
% Evaluate the existing BIE density with its original half-plane Green kernel.
b=bie.result.boundary; nn=b.panel_order; jj=(1:nn-1)'; beta=jj./sqrt(4*jj.^2-1);
[V,D]=eig(diag(beta,1)+diag(beta,-1)); [~,ord]=sort(diag(D)); gw=2*(V(1,ord).^2)';
weight=(b.half_width_m/b.panel_count)*repmat(gw,b.panel_count,1)./b.node_normal(:,2).*b.window;
xs=b.node_x_m.'; zs=b.node_z_m.'; nx=b.node_normal(:,1).'; nz=b.node_normal(:,2).';
bie_correction=complex(zeros(n,1));
for first=1:24:n
    ii=first:min(first+23,n); ddx=xv(ii)-xs; ddz=zv(ii)-zs; R=hypot(ddx,ddz);
    ddzm=zv(ii)-(2*b.halfplane_h_m-zs); Rm=hypot(ddx,ddzm);
    S=1i/4*(besselh(0,1,k*R)-besselh(0,1,k*Rm));
    K=1i*k/4*(besselh(1,1,k*R).*(ddx.*nx+ddz.*nz)./R ...
        -besselh(1,1,k*Rm).*(ddx.*nx-ddzm.*nz)./Rm);
    bie_correction(ii)=(K-1i*k*S)*(weight.*b.density(:,1));
end
fields(:,:,2)=reshape(flat+bie_correction,size(X));
% FullFK uses the saved source and frozen GL8/800-panel surface prescription.
full=closeout_s2b_fullfk_analytic_exterior(fk.result.fullfk{1}.source,.01,xv,zv,800,flat);
fields(:,:,3)=reshape(full.receiver_field,size(X));
% The receiver row is present on the grid and must reproduce the P1 trace.
iz=find(abs(z-3)<1e-12); ir=round((x+40)/.025)+1;
saved=[pe.result.reflect(ir,1),bie.result.reflect(ir,1),fk.result.reflect(ir,1)];
calculated=squeeze(fields(iz,:,:));
receiver_relative=vecnorm(calculated-saved)./vecnorm(saved);
assert(all(receiver_relative<1e-5),'2-D evaluation must reproduce accepted P1 receiver fields.');
assert(all(isfinite(fields(:))));
eta=closeout_s2b_surface(x,.01);
save(fullfile(output_dir,'spatial_complex_fields.mat'),'x','z','fields','eta','receiver_relative', ...
    'source','W','dx','-v7.3');
amplitude=abs(fields); difference=abs(fields(:,:,1)-fields(:,:,2));
fk_difference=abs(fields(:,:,3)-fields(:,:,2));
phase_difference=angle(fields(:,:,1).*conj(fields(:,:,2)))*180/pi;
fig=figure('Visible','off','Color','w','Position',[20 20 1500 850]);
tl=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
title(tl,'4 kHz weak C2 sinusoidal surface | 20 cm strip | reflected complex pressure','FontSize',14);
titles={'PE Model-0 |p_R|','BIE |p_R|','FullFK |p_R|', ...
    '|PE - BIE|','|FullFK - BIE|','arg(PE / BIE)'};
maps={amplitude(:,:,1),amplitude(:,:,2),amplitude(:,:,3),difference,fk_difference,phase_difference};
common_max=max(amplitude(:));
for j=1:6
    ax=nexttile(tl); imagesc(ax,x,z,maps{j}); axis(ax,'xy');
    if j<=3, clim(ax,[0 common_max]); elseif j==4, clim(ax,[0 max(difference(:))]); ...
    elseif j==5, clim(ax,[0 max(fk_difference(:))]); else, clim(ax,[-max(abs(phase_difference(:))) max(abs(phase_difference(:)))]); end
    title(ax,titles{j}); xlabel(ax,'x (m)'); ylabel(ax,'depth z (m)');
    cb=colorbar(ax); if j==6, ylabel(cb,'deg'); else, ylabel(cb,'peak Pa'); end
    hold(ax,'on'); yline(ax,3,'w--','LineWidth',.75);
end
savefig(fig,fullfile(output_dir,'weak_reflected_spatial_20cm.fig'));
exportgraphics(fig,fullfile(output_dir,'weak_reflected_spatial_20cm.png'),'Resolution',170);
close(fig);
receipt=struct('status','SPATIAL_FIELDS_EVALUATED_FROM_ACCEPTED_P1_DATA', ...
    'remote_matlab_version',version,'x_points',numel(x),'z_points',numel(z), ...
    'x_range_m',[min(x) max(x)],'z_range_m',[min(z) max(z)], ...
    'receiver_relative_error',receiver_relative,'max_peak_Pa',common_max, ...
    'max_PE_BIE_complex_difference_Pa',max(difference(:)), ...
    'max_FullFK_BIE_complex_difference_Pa',max(fk_difference(:)), ...
    'elapsed_s',toc(t0),'failed_checks',0,'source_run',prior_run);
fid=fopen(fullfile(output_dir,'receipt.json'),'w'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid);
fprintf('%s\n',jsonencode(receipt));
end
