function regenerate_full_kirchhoff_figure_conventions
%REGENERATE_FULL_KIRCHHOFF_FIGURE_CONVENTIONS
% Post-process saved Full-Kirchhoff/sigma-sweep MAT files only.
% No PE, BIE, or Full-Kirchhoff solver is called here.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root); setup_vertical_project();
local_parameter_sweep_figures(root);
local_sigma_figures(root);
fprintf('Visualization-only regeneration completed.\n');
end

function local_parameter_sweep_figures(root)
q=load(fullfile(root,'results','validation','pe_bie_full_kirchhoff_parameter_sweep', ...
    'pe_bie_full_kirchhoff_parameter_sweep.mat'),'result'); r=q.result; out=fullfile(root,'results','validation', ...
    'pe_bie_full_kirchhoff_parameter_sweep','figures'); if ~exist(out,'dir'),mkdir(out);end
names={'height_A0.01_K0.1','height_A0.2_K0.1','wavenumber_A0.02_K0.47'}; idx=[2 6 11];
for ii=1:numel(idx)
    c=r.cases(idx(ii)); x=linspace(-60,60, numel(c.receiver_PE)).';
    f=figure('Visible','off'); tiledlayout(3,1);
    nexttile; plot(c.surface_x_m,c.surface_z_m,'k','LineWidth',1.2); grid on;
    xlabel('$x$ (m)','Interpreter','latex'); ylabel('surface $z$ (m)','Interpreter','latex');
    title(sprintf('%s: A=%.3g m, K=%.3g rad/m; Tx z=%.1f, Rx z=%.1f', ...
        char(c.case_name),c.A_m,c.K_radpm,r.config.z_tx_m,r.config.z_rx_m),'Interpreter','none');
    nexttile; plot(x,abs(c.receiver_PE),'--',x,abs(c.receiver_FK),'-',x,abs(c.receiver_BIE),':','LineWidth',1);
    legend('PE','Full Kirchhoff','BIE','Location','best'); grid on; ylabel('magnitude','Interpreter','latex');
    nexttile; plot(x,angle(c.receiver_PE),'--',x,angle(c.receiver_FK),'-',x,angle(c.receiver_BIE),':','LineWidth',1);
    grid on; ylabel('wrapped phase','Interpreter','latex'); xlabel('receiver $x$ (m)','Interpreter','latex');
    saveas(f,fullfile(out,[names{ii} '_geometry_magnitude.png'])); saveas(f,fullfile(out,[names{ii} '_geometry_magnitude.pdf'])); close(f);
    f=figure('Visible','off'); tiledlayout(3,1);
    nexttile; plot(x,angle(c.receiver_PE.*conj(c.receiver_BIE)),'--',x,angle(c.receiver_FK.*conj(c.receiver_BIE)),'-');
    legend('PE-BIE','FK-BIE','Location','best'); grid on; ylabel('phase diff','Interpreter','latex');
    nexttile; plot(x,abs(c.receiver_PE-c.receiver_BIE),'--',x,abs(c.receiver_FK-c.receiver_BIE),'-');
    legend('PE-BIE','FK-BIE','Location','best'); grid on; ylabel('$|$complex error$|$','Interpreter','latex');
    nexttile; plot(x,20*log10(max(abs(c.receiver_PE),realmin)./max(abs(c.receiver_BIE),realmin)),'--', ...
        x,20*log10(max(abs(c.receiver_FK),realmin)./max(abs(c.receiver_BIE),realmin)),'-');
    grid on; ylabel('magnitude diff (dB)','Interpreter','latex'); xlabel('receiver $x$ (m)','Interpreter','latex');
    saveas(f,fullfile(out,[names{ii} '_errors.png'])); saveas(f,fullfile(out,[names{ii} '_errors.pdf'])); close(f);
end

% Correct the authoritative incident x-z figure: z=0 is at the top.
c=r.cases(2); s=local_source(r.config); x=linspace(-60,60,401); z=linspace(0,100,301);
[xx,zz]=meshgrid(x,z); u=local_source_eval(s,xx(:),zz(:),r.config.z_tx_m); u=reshape(u,size(xx));
f=figure('Visible','off'); imagesc(x,z,abs(u)); set(gca,'YDir','reverse'); colorbar;
xlabel('$x$ (m)','Interpreter','latex'); ylabel('depth $z$ (m)','Interpreter','latex');
title(['Incident magnitude: ' char(c.case_name)],'Interpreter','none');
saveas(f,fullfile(out,'weak_lowK_incident_xz.png')); saveas(f,fullfile(out,'weak_lowK_incident_xz.pdf')); close(f);
f=figure('Visible','off'); tiledlayout(2,1); nexttile; plot(c.surface_x_m,abs(c.surface_incident)); grid on;
ylabel('$|p_{inc}|$','Interpreter','latex'); title(['Surface incident field: ' char(c.case_name)],'Interpreter','none');
nexttile; plot(c.surface_x_m,angle(c.surface_incident)); grid on; ylabel('wrapped phase','Interpreter','latex');
xlabel('surface $x$ (m)','Interpreter','latex'); saveas(f,fullfile(out,'weak_lowK_incident_surface.png')); saveas(f,fullfile(out,'weak_lowK_incident_surface.pdf')); close(f);

h=r.cases([r.cases.family]=="height"); k=r.cases([r.cases.family]=="wavenumber");
f=figure('Visible','off'); semilogy([h.A_m],arrayfun(@(a)a.metrics_PE_BIE.complex_l2,h),'o--', ...
    [h.A_m],arrayfun(@(a)a.metrics_FK_BIE.complex_l2,h),'s-','LineWidth',1.2); grid on;
legend('PE-BIE','Full Kirchhoff-BIE','Location','best'); xlabel('$A$ (m)','Interpreter','latex'); ylabel('normalized complex L2','Interpreter','latex');
title('Height sweep, $K=0.10$ rad/m','Interpreter','latex'); saveas(f,fullfile(out,'height_sweep_complex_error.png')); saveas(f,fullfile(out,'height_sweep_complex_error.pdf')); close(f);
f=figure('Visible','off'); semilogy([k.K_radpm],arrayfun(@(a)a.metrics_PE_BIE.complex_l2,k),'o--', ...
    [k.K_radpm],arrayfun(@(a)a.metrics_FK_BIE.complex_l2,k),'s-','LineWidth',1.2); grid on;
legend('PE-BIE','Full Kirchhoff-BIE','Location','best'); xlabel('$K$ (rad/m)','Interpreter','latex'); ylabel('normalized complex L2','Interpreter','latex');
title('Wavenumber sweep, $A=0.02$ m','Interpreter','latex'); saveas(f,fullfile(out,'wavenumber_sweep_complex_error.png')); saveas(f,fullfile(out,'wavenumber_sweep_complex_error.pdf')); close(f);
end

function local_sigma_figures(root)
q=load(fullfile(root,'results','validation','pe_bie_full_kirchhoff_gaussian_sigma_sensitivity', ...
    'gaussian_sigma_sensitivity.mat'),'result'); r=q.result; c=r.config; d=r.cases; T=r.table; out=c.fig;
f=figure('Visible','off'); tiledlayout(3,2);
for i=1:numel(c.sigmas)
    v=d(i); ok=v.spectrum.propagating_mask; nexttile; plot(v.src.kx(ok)/(2*pi*c.f/c.c),sqrt(v.spectrum.power(ok)/max(v.spectrum.power(ok)))); grid on; xlim([-1 1]);
    xlabel('$k_x/k$','Interpreter','latex'); ylabel('$|A|$','Interpreter','latex');
    title(sprintf('$\\sigma=%.3g$ m, $\\theta_{99}=%.2f^\\circ$',v.sigma_m,v.spectrum.theta99_deg),'Interpreter','latex');
end
saveas(f,fullfile(out,'angular_spectrum_sigma.png')); saveas(f,fullfile(out,'angular_spectrum_sigma.pdf')); close(f);
f=figure('Visible','off'); tiledlayout(2,3); x=linspace(-60,60,241); z=linspace(0,100,241); [xx,zz]=meshgrid(x,z);
for i=1:numel(c.sigmas)
    s=local_source_sigma(c,c.sigmas(i)); u=local_source_eval(s,xx(:),zz(:),c.ztx); nexttile;
    imagesc(x,z,abs(reshape(u,size(xx)))); set(gca,'YDir','reverse'); colorbar;
    title(sprintf('$\\sigma=%.3g$ m',c.sigmas(i)),'Interpreter','latex'); xlabel('$x$ (m)','Interpreter','latex'); ylabel('$z$ (m)','Interpreter','latex');
end
saveas(f,fullfile(out,'source_beam_magnitude_sigma.png')); saveas(f,fullfile(out,'source_beam_magnitude_sigma.pdf')); close(f);
for ic=1:numel(c.cases)
    v=d(strcmp({d.case_name},c.cases(ic).name)); f=figure('Visible','off'); tiledlayout(2,1); colors=lines(numel(c.sigmas));
    for i=1:numel(v)
        nexttile(1); hold on; plot(v(i).x_m,abs(v(i).receiver_PE),'--','Color',colors(i,:), ...
            'DisplayName',sprintf('$\\sigma=%.3g$ m',v(i).sigma_m));
        plot(v(i).x_m,abs(v(i).receiver_FK),'-','Color',colors(i,:),'HandleVisibility','off');
        plot(v(i).x_m,abs(v(i).receiver_BIE),':','Color',colors(i,:),'HandleVisibility','off');
        nexttile(2); hold on; plot(v(i).x_m,angle(v(i).receiver_PE),'--','Color',colors(i,:), ...
            'DisplayName',sprintf('$\\sigma=%.3g$ m',v(i).sigma_m));
        plot(v(i).x_m,angle(v(i).receiver_FK),'-','Color',colors(i,:),'HandleVisibility','off');
        plot(v(i).x_m,angle(v(i).receiver_BIE),':','Color',colors(i,:),'HandleVisibility','off');
    end
    nexttile(1); grid on; ylabel('magnitude','Interpreter','latex'); title(c.cases(ic).name,'Interpreter','none'); legend('Location','best','Interpreter','latex');
    nexttile(2); grid on; ylabel('phase','Interpreter','latex'); xlabel('receiver $x$ (m)','Interpreter','latex');
    saveas(f,fullfile(out,[c.cases(ic).name '_receiver_fields_sigma.png'])); saveas(f,fullfile(out,[c.cases(ic).name '_receiver_fields_sigma.pdf'])); close(f);
end
f=figure('Visible','off'); tiledlayout(1,3);
for ic=1:numel(c.cases), v=T(strcmp(T.case_name,c.cases(ic).name),:); nexttile; plot(v.sigma_m,v.PE_BIE_Ec,'o--',v.sigma_m,v.PE_FK_Ec,'s-',v.sigma_m,v.FK_BIE_Ec,'d:'); set(gca,'XScale','log','YScale','log'); grid on; title(c.cases(ic).name,'Interpreter','none'); xlabel('$\sigma$ (m)','Interpreter','latex'); ylabel('$E_c$','Interpreter','latex'); end
saveas(f,fullfile(out,'error_vs_sigma.png')); saveas(f,fullfile(out,'error_vs_sigma.pdf')); close(f);
f=figure('Visible','off'); hold on; for ic=1:numel(c.cases), v=T(strcmp(T.case_name,c.cases(ic).name),:); plot(v.theta99_deg,v.PE_BIE_Ec,'o-','DisplayName',c.cases(ic).name); end; set(gca,'XScale','log','YScale','log'); grid on; legend('Interpreter','none','Location','best'); title('$\theta_{99}$ versus PE-BIE error','Interpreter','latex'); xlabel('$\theta_{99}$ (deg)','Interpreter','latex'); ylabel('$E_c$','Interpreter','latex'); saveas(f,fullfile(out,'theta99_vs_pe_error.png')); saveas(f,fullfile(out,'theta99_vs_pe_error.pdf')); close(f);
end

function s=local_source(c)
s=local_source_sigma(c,c.sigma_src_m);
end
function s=local_source_sigma(c,sig)
if isfield(c,'xw'), w=c.xw; n=c.nx; f=c.f; c0=c.c;
else, w=c.xw_m; n=c.nx; f=c.frequency_hz; c0=c.c0_mps; end
dx=w/n; x=(-n/2:n/2-1).'*dx; s.x=x; s.kx=(2*pi/w)*[0:n/2-1,-n/2:-1]; s.k=2*pi*f/c0; s.kz=sqrt(complex(s.k^2-s.kx.^2,0)); s.co=fft(exp(-.5*(x/sig).^2)).'; s.x0=x(1); s.n=n;
end
function u=local_source_eval(s,x,z,zref)
x=x(:); z=z(:); px=exp(1i*(x-s.x0)*s.kx); pz=exp(-1i*(z-zref)*s.kz); u=sum(px.*(s.co.*pz),2)/s.n;
end
