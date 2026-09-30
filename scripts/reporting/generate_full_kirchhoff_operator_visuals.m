function result = generate_full_kirchhoff_operator_visuals
%GENERATE_FULL_KIRCHHOFF_OPERATOR_VISUALS
% Regenerate validation-only Full-Kirchhoff operator visuals from saved data.
% This is post-processing only: it does not rerun PE, BIE, or a parameter sweep.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root); setup_vertical_project();
validation_file = fullfile(root,'results','validation','pe_bie_full_kirchhoff', ...
    'pe_bie_full_kirchhoff_validation.mat');
sweep_file = fullfile(root,'results','validation','pe_bie_full_kirchhoff_parameter_sweep', ...
    'pe_bie_full_kirchhoff_parameter_sweep.mat');
adaptive_file = fullfile(root,'results','validation', ...
    'pe_bie_full_kirchhoff_adaptive_aperture_stage1', ...
    'pe_bie_full_kirchhoff_adaptive_aperture_stage1.mat');
stage0_file = fullfile(root,'results','validation','pe_bellhop_controlled_comparison', ...
    'stage0','stage0_validation.mat');
assert(exist(validation_file,'file')==2,'Missing saved Full-Kirchhoff validation MAT.');
assert(exist(sweep_file,'file')==2,'Missing saved parameter-sweep MAT.');
assert(exist(adaptive_file,'file')==2,'Missing saved adaptive-aperture MAT.');
assert(exist(stage0_file,'file')==2,'Missing saved Stage-0 footprint MAT.');

q = load(validation_file,'validation'); validation = q.validation;
q = load(sweep_file,'result'); sweep = q.result;
q = load(adaptive_file,'result'); adaptive = q.result;
q = load(stage0_file,'validation'); stage0 = q.validation;
out_dir = fullfile(root,'results','visualization','pe_bie_full_kirchhoff_operator_atlas');
fig_dir = fullfile(out_dir,'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

tags = ["weak_low_K","strong_height_low_K","weak_high_K"];
summaries = repmat(local_summary(),numel(tags),1);
for ii = 1:numel(tags)
    d = validation.details(strcmp(string({validation.details.case_name}),tags(ii)));
    assert(isscalar(d),'Missing or ambiguous case %s.',tags(ii));
    [mask,weights] = local_mask(d,stage0);
    [reconstructed,relerr] = local_reconstruct(d.nominal);
    assert(relerr < 1e-12,'Saved Full-Kirchhoff reconstruction mismatch for %s.',tags(ii));
    [xr,rx_index] = local_select_receiver(d.nominal.receiver_x_m,mask,20);
    row = local_integrand_row(d.nominal,rx_index);
    local_plot_geometry(fig_dir,d,xr,tags(ii));
    local_plot_fields(fig_dir,d,mask,tags(ii));
    local_plot_fk_bie_formal(fig_dir,d,mask,weights,tags(ii));
    local_plot_errors(fig_dir,d,mask,tags(ii));
    local_plot_kernel(fig_dir,d.nominal,mask,tags(ii));
    local_plot_contribution(fig_dir,d.nominal,row,xr,tags(ii));
    summaries(ii) = local_make_summary(d,tags(ii),xr,rx_index,mask,weights, ...
        relerr,row,reconstructed);
end
local_plot_aperture(fig_dir,adaptive);
local_plot_sweeps(fig_dir,sweep);

summary = struct2table(summaries);
writetable(summary,fullfile(out_dir,'operator_visualization_summary.csv'));
result = struct('schema_version','1.0.0','scope', ...
    'validation-only post-processing of saved Full-Kirchhoff/BIE/PE results', ...
    'integrand_definition',['I(x_r,x_s)=-[G(x_r,x_s) d_n p_inc(x_s)+', ...
    'd_n_s G(x_r,x_s) p_inc(x_s)] w_s'], ...
    'comparison_convention',validation.comparison_convention, ...
    'summary',summary,'figure_dir',fig_dir);
save(fullfile(out_dir,'full_kirchhoff_operator_visualization.mat'),'result','-v7.3');
local_write_report(fullfile(root,'reports','pe_bie_full_kirchhoff_operator_visualization_report.md'), ...
    result,root,out_dir);
disp(summary);
end

function [mask,w] = local_mask(d,stage0)
mask = stage0.footprint.m99.mask(:) & isfinite(d.receiver_full(:)) & ...
    isfinite(d.receiver_bie(:));
floor_amp = 10^(-40/20) * max([abs(d.receiver_full(:)); abs(d.receiver_bie(:))]);
mask = mask & abs(d.receiver_full(:)) >= floor_amp & abs(d.receiver_bie(:)) >= floor_amp;
w = stage0.footprint.energy_weights(:); w = w/max(sum(w(mask)),realmin);
end

function [u,err] = local_reconstruct(nominal)
xr = nominal.receiver_x_m(:); zr = nominal.receiver_z_m(:);
u = complex(zeros(size(xr))); block = 128;
for first = 1:block:numel(xr)
    rows = first:min(first+block-1,numel(xr));
    u(rows) = -sum(local_integrand(nominal,xr(rows),zr(rows)),2);
end
err = norm(u-nominal.receiver_field(:))/max(norm(nominal.receiver_field(:)),realmin);
end

function values = local_integrand(nominal,xr,zr)
xs = nominal.surface_x_m(:); zs = nominal.surface_z_m(:); n = nominal.surface_normal;
dx = xr(:)-xs.'; dz = zr(:)-zs.'; rr = hypot(dx,dz); k = nominal.k_radpm;
G = 1i/4*besselh(0,1,k*rr);
dGdn = 1i*k/4*besselh(1,1,k*rr).*(dx.*n(:,1).'+dz.*n(:,2).')./rr;
values = (G.*nominal.surface_incident_dn(:).' + ...
    dGdn.*nominal.surface_incident(:).') .* nominal.quadrature_weights_m(:).';
end

function row = local_integrand_row(nominal,index)
row = local_integrand(nominal,nominal.receiver_x_m(index),nominal.receiver_z_m(index));
end

function [x,index] = local_select_receiver(grid,mask,target)
eligible = find(mask(:)); [~,jj] = min(abs(grid(eligible)-target)); index = eligible(jj); x = grid(index);
end

function local_plot_geometry(fd,d,xr,tag)
sx = d.surface_x_m(:); sz = d.surface_z_m(:); cfg = d.config;
f = figure('Visible','off'); tiledlayout(1,2,'TileSpacing','compact');
nexttile; hold on;
patch([sx(1);sx(end);sx(end);sx(1)],[min(sz)-0.02;min(sz)-0.02;max(sz)+0.02;max(sz)+0.02], ...
    [0.92 0.92 0.92],'EdgeColor','none','DisplayName','Full surface support');
plot(sx,sz,'k','LineWidth',1.4,'DisplayName','rough surface');
plot(0,cfg.z_tx_m,'rp','MarkerSize',12,'MarkerFaceColor','r','DisplayName','Tx');
plot([min(d.nominal.receiver_x_m) max(d.nominal.receiver_x_m)], ...
    cfg.z_rx_m*[1 1],'b-','LineWidth',1.1,'DisplayName','receiver line');
plot(xr,cfg.z_rx_m,'bo','MarkerFaceColor','b','DisplayName',sprintf('representative Rx x=%.3g m',xr));
for L = [4 16 64 128]
    plot([xr-L xr+L],[min(sz)-0.015 min(sz)-0.015],'LineWidth',2, ...
        'DisplayName',sprintf('aperture L=%g m',L));
end
set(gca,'YDir','reverse'); axis tight; grid on; xlabel('x (m)','Interpreter','latex'); ylabel('z (m)','Interpreter','latex');
title(sprintf('%s: full geometry',tag),'Interpreter','none'); legend('Location','eastoutside');
nexttile; hold on;
plot(sx,sz,'k','LineWidth',1.5,'DisplayName','rough surface');
zbar = min(sz)-0.03 - (0:3)*0.04;
for jj=1:4
    L=[4 16 64 128]; plot([xr-L(jj) xr+L(jj)],[zbar(jj) zbar(jj)],'LineWidth',2, ...
        'DisplayName',sprintf('L=%g m',L(jj)));
end
xlim([min(sx) max(sx)]); ylim([min(sz)-0.20 max(sz)+0.20]); set(gca,'YDir','reverse');
grid on; xlabel('$x_s$ (m)','Interpreter','latex'); ylabel('$z_s$ (m)','Interpreter','latex');
title('surface support and aperture zoom','Interpreter','none');
sgtitle(sprintf('%s: geometry and aperture support',tag),'Interpreter','none');
saveas(f,fullfile(fd,[char(tag) '_geometry_support_aperture.png']));
saveas(f,fullfile(fd,[char(tag) '_geometry_support_aperture.pdf'])); close(f);
end

function local_plot_fields(fd,d,mask,tag)
x=d.nominal.receiver_x_m(:); f=figure('Visible','off'); tiledlayout(2,1);
nexttile; plot(x(mask),abs(d.receiver_model0(mask)),'--',x(mask),abs(d.receiver_model1(mask)),'-', ...
    x(mask),abs(d.receiver_full(mask)),'-.',x(mask),abs(d.receiver_bie(mask)),':','LineWidth',1.1); grid on;
legend('Model-0','Model-1','Full Kirchhoff','BIE','Location','best'); ylabel('normalized magnitude','Interpreter','latex');
title([char(tag) ' receiver line on formal M99 / -40 dB mask'],'Interpreter','none');
nexttile; plot(x(mask),angle(d.receiver_model0(mask)),'--',x(mask),angle(d.receiver_model1(mask)),'-', ...
    x(mask),angle(d.receiver_full(mask)),'-.',x(mask),angle(d.receiver_bie(mask)),':','LineWidth',1.1); grid on;
ylabel('wrapped phase (rad)','Interpreter','latex'); xlabel('receiver x (m)','Interpreter','latex');
saveas(f,fullfile(fd,[char(tag) '_receiver_fields_masked.png']));
saveas(f,fullfile(fd,[char(tag) '_receiver_fields_masked.pdf'])); close(f);
end

function local_plot_fk_bie_formal(fd,d,mask,w,tag)
x=d.nominal.receiver_x_m(:); fk=d.receiver_full(:); bie=d.receiver_bie(:);
peak=max([abs(fk(mask));abs(bie(mask))]); mask=mask(:)&isfinite(fk)&isfinite(bie) & ...
    abs(fk)>=10^(-40/20)*peak & abs(bie)>=10^(-40/20)*peak;
w=w(:); w=w(mask)/max(sum(w(mask)),realmin); a=fk(mask); b=bie(mask);
ec=sqrt(sum(w.*abs(a-b).^2)/max(sum(w.*abs(b).^2),realmin));
f=figure('Visible','off'); tiledlayout(2,1,'TileSpacing','compact');
nexttile; plot(x(mask),abs(a),'k-',x(mask),abs(b),'r--','LineWidth',1.2); grid on;
legend('Full Kirchhoff','BIE','Location','best'); ylabel('normalized magnitude','Interpreter','latex');
title(sprintf('FK/BIE on formal mask: $E_c=%.6g$',ec),'Interpreter','latex');
nexttile; plot(x(mask),angle(a),'k-',x(mask),angle(b),'r--','LineWidth',1.2); grid on;
ylabel('wrapped phase (rad)','Interpreter','latex'); xlabel('receiver x (m)','Interpreter','latex');
saveas(f,fullfile(fd,[char(tag) '_FK_BIE_formal_mask.png']));
saveas(f,fullfile(fd,[char(tag) '_FK_BIE_formal_mask.pdf'])); close(f);
end

function local_plot_errors(fd,d,mask,tag)
x=d.nominal.receiver_x_m(:); pairs = {d.receiver_model1,d.receiver_bie,'PE Model-1 - BIE'; ...
    d.receiver_model1,d.receiver_full,'PE Model-1 - Full Kirchhoff'; ...
    d.receiver_full,d.receiver_bie,'Full Kirchhoff - BIE'};
f=figure('Visible','off'); tiledlayout(3,3,'TileSpacing','compact');
for ii=1:3
    a=pairs{ii,1}; b=pairs{ii,2};
    nexttile(ii); plot(x(mask),20*log10(max(abs(a(mask)),realmin)./max(abs(b(mask)),realmin)));grid on;title(pairs{ii,3},'Interpreter','none');ylabel('magnitude diff (dB)','Interpreter','latex');
    nexttile(ii+3); plot(x(mask),angle(a(mask).*conj(b(mask))));grid on;ylabel('phase diff (rad)','Interpreter','latex');
    nexttile(ii+6); plot(x(mask),abs(a(mask)-b(mask)));grid on;ylabel('|complex error|','Interpreter','latex');xlabel('x (m)','Interpreter','latex');
end
sgtitle([char(tag) ' error decomposition on formal common mask'],'Interpreter','none');
saveas(f,fullfile(fd,[char(tag) '_receiver_error_decomposition.png']));
saveas(f,fullfile(fd,[char(tag) '_receiver_error_decomposition.pdf'])); close(f);
end

function local_plot_kernel(fd,nominal,mask,tag)
rxidx=find(mask); rxidx=rxidx(round(linspace(1,numel(rxidx),min(401,numel(rxidx)))));
sxidx=round(linspace(1,numel(nominal.surface_x_m),min(801,numel(nominal.surface_x_m))));
sub=local_integrand_subset(nominal,rxidx,sxidx);
f=figure('Visible','off'); tiledlayout(2,1);
nexttile; imagesc(nominal.surface_x_m(sxidx),nominal.receiver_x_m(rxidx),log10(max(abs(sub),realmin)));axis xy;colorbar;ylabel('$x_r$ (m)','Interpreter','latex');title([char(tag) ': integrand magnitude'],'Interpreter','none');
nexttile; imagesc(nominal.surface_x_m(sxidx),nominal.receiver_x_m(rxidx),angle(sub));axis xy;colorbar;clim([-pi pi]);xlabel('$x_s$ (m)','Interpreter','latex');ylabel('$x_r$ (m)','Interpreter','latex');title('integrand phase (rad)','Interpreter','none');
saveas(f,fullfile(fd,[char(tag) '_integrand_kernel.png']));
saveas(f,fullfile(fd,[char(tag) '_integrand_kernel.pdf']));close(f);
end

function out=local_integrand_subset(nominal,rxidx,sxidx)
xs=nominal.surface_x_m(sxidx);zs=nominal.surface_z_m(sxidx);n=nominal.surface_normal(sxidx,:);
xr=nominal.receiver_x_m(rxidx);zr=nominal.receiver_z_m(rxidx);dx=xr(:)-xs(:).';dz=zr(:)-zs(:).';rr=hypot(dx,dz);k=nominal.k_radpm;
G=1i/4*besselh(0,1,k*rr);dGdn=1i*k/4*besselh(1,1,k*rr).*(dx.*n(:,1).'+dz.*n(:,2).')./rr;
out=(G.*nominal.surface_incident_dn(sxidx).' + dGdn.*nominal.surface_incident(sxidx).') .* nominal.quadrature_weights_m(sxidx).';
end

function local_plot_contribution(fd,nominal,row,xr,tag)
xs=nominal.surface_x_m(:);row=row(:);f=figure('Visible','off');tiledlayout(4,1,'TileSpacing','compact');
nexttile;plot(xs,abs(row));grid on;ylabel('$|I|$','Interpreter','latex');title(sprintf('%s: surface contribution to receiver x_r=%.3g m',tag,xr),'Interpreter','none');
nexttile;plot(xs,angle(row));grid on;ylabel('$\arg I$ (rad)','Interpreter','latex');
nexttile;plot(xs,real(row),xs,imag(row));grid on;legend('real','imag','Location','best');ylabel('$I$','Interpreter','latex');
nexttile;cs=-cumsum(row);plot(xs,abs(cs));grid on;ylabel('$| -\sum_{j\leq i} I_j |$','Interpreter','latex');xlabel('$x_s$ (m)','Interpreter','latex');
saveas(f,fullfile(fd,[char(tag) '_single_receiver_contribution.png']));
saveas(f,fullfile(fd,[char(tag) '_single_receiver_contribution.pdf']));close(f);
end

function local_plot_aperture(fd,r)
t=r.aperture_rows; names=unique(t.case_name,'stable'); f=figure('Visible','off');hold on;
for ii=1:numel(names)
    q=t(strcmp(t.case_name,names{ii}),:); L=q.L_m; L(isinf(L))=256;
    semilogy(L,q.complex_l2,'o-','DisplayName',names{ii},'LineWidth',1.1);
end
yline(r.complex_l2_gate,'k--','$E_c=10^{-3}$','Interpreter','latex');grid on;xlabel('$L$ (m); 256 denotes Inf','Interpreter','latex');ylabel('Reduced Kirchhoff-BIE complex L2','Interpreter','latex');title('Finite-aperture convergence for saved adaptive database','Interpreter','none');legend('Location','eastoutside');
saveas(f,fullfile(fd,'all_cases_aperture_convergence.png'));saveas(f,fullfile(fd,'all_cases_aperture_convergence.pdf'));close(f);
end

function local_plot_sweeps(fd,r)
c=r.cases; h=c(strcmp(string({c.family}),"height")); k=c(strcmp(string({c.family}),"wavenumber"));
local_plot_one_sweep(fd,[h.A_m],h,'A (m)','K=0.10 rad/m','height_sweep_all_metrics');
local_plot_one_sweep(fd,[k.K_radpm],k,'K (rad/m)','A=0.02 m','wavenumber_sweep_all_metrics');
end

function local_plot_one_sweep(fd,x,c,xlabel_text,title_text,stem)
f=figure('Visible','off');tiledlayout(1,3); labels={'complex L2','magnitude L2','phase RMS (rad)'}; fields={'complex_l2','magnitude_relative_l2','phase_rms_rad'};
for ii=1:3,nexttile;semilogy(x,arrayfun(@(q)q.metrics_PE_BIE.(fields{ii}),c),'o--',x,arrayfun(@(q)q.metrics_FK_BIE.(fields{ii}),c),'s-','LineWidth',1.1);grid on;xlabel(xlabel_text,'Interpreter','latex');ylabel(labels{ii},'Interpreter','latex');title(title_text,'Interpreter','none');if ii==1,legend('PE-BIE','Full Kirchhoff-BIE','Location','best');end,end
sgtitle(['Saved deterministic sweep: ' title_text],'Interpreter','none');saveas(f,fullfile(fd,[stem '.png']));saveas(f,fullfile(fd,[stem '.pdf']));close(f);
end

function s=local_summary(),s=struct('case_name','','receiver_x_m',NaN,'mask_samples',NaN,'reconstruction_relative_l2',NaN,'selected_receiver_native_relative_l2',NaN,'selected_receiver_contribution_peak',NaN,'selected_receiver_cancellation_ratio',NaN);end
function s=local_make_summary(d,tag,xr,index,mask,~,err,row,reconstructed)
s=local_summary();s.case_name=char(tag);s.receiver_x_m=xr;s.mask_samples=sum(mask);s.reconstruction_relative_l2=err;s.selected_receiver_native_relative_l2=abs(reconstructed(index)-d.nominal.receiver_field(index))/max(abs(d.nominal.receiver_field(index)),realmin);s.selected_receiver_contribution_peak=max(abs(row));s.selected_receiver_cancellation_ratio=sum(abs(row))/max(abs(sum(row)),realmin);
end

function local_write_report(path,r,root,out_dir)
fid=fopen(path,'w','n','UTF-8');assert(fid>=0);clean=onCleanup(@()fclose(fid));
fprintf(fid,'# Full Kirchhoff operator visualization atlas\n\n');
fprintf(fid,'Status: **post-processing only**. No PE, BIE, Full-Kirchhoff, or sweep computation was rerun.\n\n');
fprintf(fid,'The plotted scalar density is `%s`. It contains the actual surface position, unit normal, source-normal Green derivative, incident normal derivative, and saved arc-length quadrature weight. It is therefore the correct visual object for this implementation; `G` alone is not called the Full-Kirchhoff kernel.\n\n',r.integrand_definition);
fprintf(fid,'## Reconstruction gate\n\n| case | representative receiver x m | mask samples | field reconstruction relative L2 | selected-row relative error | coherent cancellation ratio |\n|---|---:|---:|---:|---:|---:|\n');
t=r.summary;for ii=1:height(t),fprintf(fid,'| %s | %.6g | %d | %.3g | %.3g | %.6g |\n',t.case_name{ii},t.receiver_x_m(ii),t.mask_samples(ii),t.reconstruction_relative_l2(ii),t.selected_receiver_native_relative_l2(ii),t.selected_receiver_cancellation_ratio(ii));end
fprintf(fid,'\nThe reconstruction gate sums the plotted density over the saved surface samples and compares it with the saved native Full-Kirchhoff receiver field before applying the frozen comparison conjugation/normalization.\n\n');
fprintf(fid,'Figures and summary data: `%s`.\n',local_relative(root,out_dir));
end
function p=local_relative(root,p),p=strrep(p,[root filesep],'');p=strrep(p,filesep,'/');end
