function result = validate_reduced_kirchhoff_adaptive_aperture_stage1
%VALIDATE_REDUCED_KIRCHHOFF_ADAPTIVE_APERTURE_STAGE1
% Validation-only first database for adaptive finite-aperture Kirchhoff.
% Uses the accepted deterministic A and K sweeps; production PE, BIE, the
% Full-Kirchhoff formulation, and the surface model are not modified.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root); setup_vertical_project();
in_file = fullfile(root,'results','validation','pe_bie_full_kirchhoff_parameter_sweep', ...
    'pe_bie_full_kirchhoff_parameter_sweep.mat');
stage0_file = fullfile(root,'results','validation','pe_bellhop_controlled_comparison', ...
    'stage0','stage0_validation.mat');
assert(exist(in_file,'file')==2,'Missing accepted Full-Kirchhoff parameter sweep.');
q = load(in_file,'result'); sweep = q.result;
q = load(stage0_file,'validation'); stage0 = q.validation;

out_dir = fullfile(root,'results','validation','pe_bie_full_kirchhoff_adaptive_aperture_stage1');
fig_dir = fullfile(out_dir,'figures');
if ~exist(out_dir,'dir'), mkdir(out_dir); end
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

cfg = sweep.config;
source = local_source(cfg);
x_rx = stage0.x_m(:); z_rx = cfg.z_rx_m + zeros(size(x_rx));
flat_ref = -source.evaluate(x_rx,z_rx,cfg.z_tx_m+cfg.z_rx_m,'path');
L_grid = [4 8 16 32 64 128 Inf];
gate = 1e-3;

% The prior A sweep and K sweep contain one duplicate A=0.02, K=0.10 case.
% Keep the first occurrence: eleven distinct, authoritative geometries.
all_cases = sweep.cases;
keep = true(numel(all_cases),1);
for ii=1:numel(all_cases)
    for jj=1:ii-1
        if abs(all_cases(ii).A_m-all_cases(jj).A_m)<eps && ...
                abs(all_cases(ii).K_radpm-all_cases(jj).K_radpm)<eps
            keep(ii)=false;
            break
        end
    end
end
cases = all_cases(keep);
rows = repmat(local_row(),0,1);
database_rows = repmat(local_database_row(),numel(cases),1);
details = repmat(local_detail(),numel(cases),1);

for cc=1:numel(cases)
    c = cases(cc);
    sx = c.surface_x_m(:); sz = c.surface_z_m(:);
    ep = local_surface_prime(sx,c.A_m,c.K_radpm,cfg.surface_inner_m,cfg.surface_support_m);
    jac = sqrt(1+ep.^2); normal = [-ep./jac,1./jac];
    [ui,ux,uz] = source.evaluate_with_derivatives(sx,sz,'coordinate');
    dni = normal(:,1).*ux + normal(:,2).*uz;
    common = struct('frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps, ...
        'surface_x_m',sx,'surface_z_m',sz,'surface_normal',normal, ...
        'surface_incident',ui,'surface_incident_dn',dni, ...
        'receiver_x_m',x_rx,'receiver_z_m',z_rx);
    [mask,w] = local_fixed_mask(c.receiver_FK,stage0,cfg.phase_floor_db);
    m_pe = local_metrics(c.receiver_PE,c.receiver_BIE,mask,w);
    m_fk = local_metrics(c.receiver_FK,c.receiver_BIE,mask,w);
    rk = complex(zeros(numel(x_rx),numel(L_grid)));
    m_rk = repmat(local_metrics_struct(),numel(L_grid),1);
    full_equiv = NaN;
    for ll=1:numel(L_grid)
        run_cfg = common; run_cfg.local_radius_m = L_grid(ll);
        out = solve_full_kirchhoff_2d(run_cfg);
        rk(:,ll) = conj(out.receiver_field(:)./flat_ref);
        m_rk(ll) = local_metrics(rk(:,ll),c.receiver_BIE,mask,w);
        if isinf(L_grid(ll))
            full_equiv = local_relative_l2(rk(:,ll),c.receiver_FK,mask,w);
        end
        rows(end+1) = local_make_row(local_case_name(c),c.A_m,c.K_radpm, ...
            c.geometry.slope_AK,c.geometry.curvature_AK2,L_grid(ll), ...
            m_rk(ll),m_pe.complex_l2/max(m_rk(ll).complex_l2,realmin)); %#ok<AGROW>
    end
    assert(full_equiv<1e-12,'L=Inf did not reproduce accepted FK for %s.',local_case_name(c));
    finite = find(isfinite(L_grid));
    jj = finite(find([m_rk(finite).complex_l2] < gate,1,'first'));
    requires_full = isempty(jj);
    if requires_full
        chosen = numel(L_grid); % Inf; reports that no finite L passes.
    else
        chosen = jj;
    end
    database_rows(cc) = local_make_database_row(local_case_name(c),c.A_m,c.K_radpm, ...
        c.geometry.slope_AK,c.geometry.curvature_AK2,L_grid(chosen), ...
        m_rk(chosen),requires_full,m_pe,m_fk,full_equiv);
    details(cc) = struct('case_name',local_case_name(c),'family',char(c.family), ...
        'A_m',c.A_m,'K_radpm',c.K_radpm,'slope_AK',c.geometry.slope_AK, ...
        'curvature_AK2',c.geometry.curvature_AK2,'x_m',x_rx, ...
        'receiver_PE',c.receiver_PE,'receiver_BIE',c.receiver_BIE, ...
        'receiver_FK',c.receiver_FK,'receiver_RK',rk,'L_m',L_grid, ...
        'metrics_PE_BIE',m_pe,'metrics_FK_BIE',m_fk,'metrics_RK_BIE',m_rk, ...
        'selected_index',chosen,'requires_full_K',requires_full, ...
        'full_equivalence_l2',full_equiv,'mask_count',sum(mask));
end

aperture_table = struct2table(rows);
database = struct2table(database_rows);
local_plot_representatives(fig_dir,details);
local_plot_convergence(fig_dir,details);
local_plot_adaptive_vs_pe(fig_dir,database);
local_plot_relationships(fig_dir,database);
local_plot_sparse_map(fig_dir,database);

result = struct('schema_version','1.0.0','stage','adaptive_reduced_kirchhoff_stage1', ...
    'scope','validation-only deterministic A/K sweeps; production models unchanged', ...
    'config',cfg,'L_grid_m',L_grid,'complex_l2_gate',gate, ...
    'operator_definition',['Accepted Full-Kirchhoff line integral retaining surface ', ...
        'position, Green and normal derivatives, incident field and ds; truncate ', ...
        'only source points with |xr-xs|>=L.'], ...
    'mask_definition','Fixed Stage-0 M99 plus accepted-FK -40 dB threshold and saved incident-energy weights', ...
    'database',database,'aperture_rows',aperture_table,'details',details);
mat_file = fullfile(out_dir,'pe_bie_full_kirchhoff_adaptive_aperture_stage1.mat');
db_csv = fullfile(out_dir,'adaptive_aperture_database.csv');
scan_csv = fullfile(out_dir,'adaptive_aperture_scan.csv');
save(mat_file,'result','-v7.3'); writetable(database,db_csv); writetable(aperture_table,scan_csv);
report = fullfile(root,'reports','pe_bie_full_kirchhoff_adaptive_aperture_stage1_report.md');
local_write_report(report,result,root,mat_file,db_csv,scan_csv);
disp(database); fprintf('Adaptive-aperture Stage-1 report: %s\n',report);
end

function source = local_source(cfg)
dx=cfg.xw_m/cfg.nx; x=(-cfg.nx/2:cfg.nx/2-1).'*dx;
kx=(2*pi/cfg.xw_m)*[0:(cfg.nx/2-1),-cfg.nx/2:-1];
k=2*pi*cfg.frequency_hz/cfg.c0_mps; kz=sqrt(complex(k^2-kx.^2,0));
source=struct('kx',kx,'kz',kz,'coeff',fft(exp(-0.5*(x/cfg.sigma_src_m).^2)).', ...
    'x0',x(1),'n',cfg.nx,'z_tx',cfg.z_tx_m);
source.evaluate=@(xq,zq,path,mode)local_eval(source,xq,zq,path,mode);
source.evaluate_with_derivatives=@(xq,zq,mode)local_derivatives(source,xq,zq,mode);
end
function u=local_eval(s,xq,zq,path,mode)
xq=xq(:); zq=zq(:); if isscalar(zq),zq=zq+zeros(size(xq));end
px=exp(1i*(xq-s.x0)*s.kx);
if strcmp(mode,'path')
    pz=exp(1i*path*s.kz);
elseif strcmp(mode,'coordinate')
    pz=exp(-1i*(zq-path)*s.kz);
else
    error('Unknown source mode %s.',mode);
end
u=sum(px.*(s.coeff.*pz),2)/s.n;
end
function [u,ux,uz]=local_derivatives(s,xq,zq,mode)
assert(strcmp(mode,'coordinate'),'Derivative requires coordinate source mode.');
xq=xq(:);zq=zq(:);px=exp(1i*(xq-s.x0)*s.kx);pz=exp(-1i*(zq-s.z_tx)*s.kz);
u=sum(px.*(s.coeff.*pz),2)/s.n;
ux=sum(px.*((1i*s.kx).*s.coeff.*pz),2)/s.n;
uz=sum(px.*((-1i*s.kz).*s.coeff.*pz),2)/s.n;
end
function ep=local_surface_prime(x,A,K,inner,outer)
r=abs(x);chi=zeros(size(x));dchi=zeros(size(x));chi(r<=inner)=1;
mid=r>inner&r<outer;t=(r(mid)-inner)/(outer-inner);
chi(mid)=1-10*t.^3+15*t.^4-6*t.^5;
dchi(mid)=(-30*t.^2+60*t.^3-30*t.^4)/(outer-inner);
ep=A*(K*cos(K*x).*chi+sin(K*x).*dchi.*sign(x));ep(x==0)=A*K;
end
function [mask,w]=local_fixed_mask(fk,stage0,floor_db)
fk=fk(:);mask=stage0.footprint.m99.mask(:)&isfinite(fk);
mask=mask&abs(fk)>=10^(floor_db/20)*max(abs(fk));
w=stage0.footprint.energy_weights(:);w=w/max(sum(w(mask)),realmin);
end
function m=local_metrics(a,b,mask,w)
a=a(:);b=b(:);ok=mask&isfinite(a)&isfinite(b);ww=w(ok);ww=ww/max(sum(ww),realmin);
a=a(ok);b=b(ok);phase=angle(a.*conj(b));ma=abs(a);mb=abs(b);
m=struct('complex_l2',sqrt(sum(ww.*abs(a-b).^2)/max(sum(ww.*abs(b).^2),realmin)), ...
 'magnitude_relative_l2',sqrt(sum(ww.*(ma-mb).^2)/max(sum(ww.*mb.^2),realmin)), ...
 'phase_rms_rad',sqrt(sum(ww.*phase.^2)), ...
 'magnitude_correlation',local_corr(ma,mb,ww), ...
 'complex_phase_correlation',abs(sum(ww.*exp(1i*phase))), ...
 'sample_count',numel(a));
end
function e=local_relative_l2(a,b,mask,w)
a=a(:);b=b(:);ok=mask&isfinite(a)&isfinite(b);ww=w(ok);ww=ww/max(sum(ww),realmin);
e=sqrt(sum(ww.*abs(a(ok)-b(ok)).^2)/max(sum(ww.*abs(b(ok)).^2),realmin));
end
function r=local_corr(a,b,w)
a=a-sum(w.*a);b=b-sum(w.*b);
r=abs(sum(w.*a.*b))/max(sqrt(sum(w.*a.^2)*sum(w.*b.^2)),realmin);
end
function name=local_case_name(c),name=sprintf('A%.4g_K%.4g',c.A_m,c.K_radpm);end

function r=local_metrics_struct(),r=struct('complex_l2',NaN,'magnitude_relative_l2',NaN,'phase_rms_rad',NaN,'magnitude_correlation',NaN,'complex_phase_correlation',NaN,'sample_count',NaN);end
function r=local_row(),r=struct('case_name','','A_m',NaN,'K_radpm',NaN,'slope_AK',NaN,'curvature_AK2',NaN,'L_m',NaN,'complex_l2',NaN,'magnitude_relative_l2',NaN,'phase_rms_rad',NaN,'magnitude_correlation',NaN,'complex_phase_correlation',NaN,'reduction_factor_vs_PE',NaN,'sample_count',NaN);end
function r=local_make_row(name,A,K,slope,curvature,L,m,reduction)
r=local_row();r.case_name=name;r.A_m=A;r.K_radpm=K;r.slope_AK=slope;r.curvature_AK2=curvature;r.L_m=L;
r.complex_l2=m.complex_l2;r.magnitude_relative_l2=m.magnitude_relative_l2;r.phase_rms_rad=m.phase_rms_rad;r.magnitude_correlation=m.magnitude_correlation;r.complex_phase_correlation=m.complex_phase_correlation;r.reduction_factor_vs_PE=reduction;r.sample_count=m.sample_count;
end
function r=local_database_row(),r=struct('case_name','','A_m',NaN,'K_radpm',NaN,'slope_AK',NaN,'curvature_AK2',NaN,'L_min_m',NaN,'Ec_at_Lmin',NaN,'phase_error_rad_at_Lmin',NaN,'magnitude_error_at_Lmin',NaN,'requires_full_K',false,'PE_BIE_complex_l2',NaN,'FK_BIE_complex_l2',NaN,'full_equivalence_l2',NaN);end
function r=local_make_database_row(name,A,K,slope,curvature,L,m,requires_full,mpe,mfk,full_equiv)
r=local_database_row();r.case_name=name;r.A_m=A;r.K_radpm=K;r.slope_AK=slope;r.curvature_AK2=curvature;r.L_min_m=L;r.Ec_at_Lmin=m.complex_l2;r.phase_error_rad_at_Lmin=m.phase_rms_rad;r.magnitude_error_at_Lmin=m.magnitude_relative_l2;r.requires_full_K=requires_full;r.PE_BIE_complex_l2=mpe.complex_l2;r.FK_BIE_complex_l2=mfk.complex_l2;r.full_equivalence_l2=full_equiv;
end
function r=local_detail(),r=struct('case_name','','family','','A_m',NaN,'K_radpm',NaN,'slope_AK',NaN,'curvature_AK2',NaN,'x_m',[],'receiver_PE',[],'receiver_BIE',[],'receiver_FK',[],'receiver_RK',[],'L_m',[],'metrics_PE_BIE',[],'metrics_FK_BIE',[],'metrics_RK_BIE',[],'selected_index',NaN,'requires_full_K',false,'full_equivalence_l2',NaN,'mask_count',NaN);end

function local_plot_representatives(fig_dir,details)
tags={'A0.01_K0.1','A0.2_K0.1','A0.02_K0.47'};
for ii=1:numel(tags)
    d=details(strcmp({details.case_name},tags{ii}));if isempty(d),continue;end
    x=d.x_m;rk=d.receiver_RK(:,d.selected_index);ltxt=local_L_text(d.L_m(d.selected_index));
    f=figure('Visible','off');tiledlayout(2,1);
    nexttile;plot(x,abs(d.receiver_PE),'--',x,abs(rk),'-',x,abs(d.receiver_FK),':',x,abs(d.receiver_BIE),'k-.','LineWidth',1);grid on;legend('PE phase screen',['adaptive RK, L=' ltxt],'Full Kirchhoff','BIE','Location','best');ylabel('magnitude');title(sprintf('%s receiver fields',d.case_name));
    nexttile;plot(x,angle(d.receiver_PE),'--',x,angle(rk),'-',x,angle(d.receiver_FK),':',x,angle(d.receiver_BIE),'k-.','LineWidth',1);grid on;ylabel('wrapped phase');xlabel('receiver x (m)');
    saveas(f,fullfile(fig_dir,[d.case_name '_adaptive_fields.png']));saveas(f,fullfile(fig_dir,[d.case_name '_adaptive_fields.pdf']));close(f);
end
end
function local_plot_convergence(fig_dir,details)
tags={'A0.01_K0.1','A0.2_K0.1','A0.02_K0.47'};f=figure('Visible','off');hold on;
for ii=1:numel(tags)
    d=details(strcmp({details.case_name},tags{ii}));if isempty(d),continue;end
    ec=[d.metrics_RK_BIE.complex_l2];L=d.L_m;Lplot=L;Lplot(isinf(Lplot))=256;
    semilogy(Lplot,ec,'o-','LineWidth',1.2,'DisplayName',d.case_name);
    yline(d.metrics_FK_BIE.complex_l2,'--','HandleVisibility','off');
    yline(d.metrics_PE_BIE.complex_l2,':','HandleVisibility','off');
end
yline(1e-3,'k-.','E_c=10^{-3}','HandleVisibility','off');grid on;legend('Location','best');xlabel('L (m); 256 denotes Inf');ylabel('complex L2 relative to BIE');title('Adaptive-aperture convergence; horizontal hidden lines are PE/FK levels');
saveas(f,fullfile(fig_dir,'representative_aperture_convergence.png'));saveas(f,fullfile(fig_dir,'representative_aperture_convergence.pdf'));close(f);
end
function local_plot_adaptive_vs_pe(fig_dir,database)
f=figure('Visible','off');idx=1:height(database);bar(idx,[database.PE_BIE_complex_l2,database.Ec_at_Lmin]);grid on;legend('PE-BIE','adaptive RK-BIE','Location','best');xticks(idx);xticklabels(database.case_name);xtickangle(45);ylabel('complex L2');title('Adaptive-aperture improvement relative to PE');saveas(f,fullfile(fig_dir,'adaptive_vs_PE.png'));saveas(f,fullfile(fig_dir,'adaptive_vs_PE.pdf'));close(f);
end
function local_plot_relationships(fig_dir,database)
finite=isfinite(database.L_min_m);f=figure('Visible','off');tiledlayout(2,2);vars={'A_m','K_radpm','slope_AK','curvature_AK2'};labels={'A (m)','K (rad/m)','A K','A K^2 (1/m)'};
for ii=1:4,nexttile;scatter(database.(vars{ii})(finite),database.L_min_m(finite),60,'filled');grid on;xlabel(labels{ii});ylabel('L_{min} (m)');end
sgtitle('Measured finite L_{min} relationships; censored full-K cases omitted');saveas(f,fullfile(fig_dir,'Lmin_relationships.png'));saveas(f,fullfile(fig_dir,'Lmin_relationships.pdf'));close(f);
end
function local_plot_sparse_map(fig_dir,database)
A=unique(database.A_m);K=unique(database.K_radpm);M=nan(numel(A),numel(K));
for ii=1:height(database),ia=find(A==database.A_m(ii),1);ik=find(K==database.K_radpm(ii),1);M(ia,ik)=database.L_min_m(ii);end
f=figure('Visible','off');imagesc(K,A,M);set(gca,'YDir','normal');colorbar;hold on;[rr,cc]=find(isnan(M));plot(K(cc),A(rr),'wx','MarkerSize',8,'LineWidth',1.2);xlabel('K (rad/m)');ylabel('A (m)');title('Sparse A-K L_{min} map; white x = not sampled in Stage 1');saveas(f,fullfile(fig_dir,'Lmin_sparse_AK_map.png'));saveas(f,fullfile(fig_dir,'Lmin_sparse_AK_map.pdf'));close(f);
end
function s=local_L_text(L),if isinf(L),s='Inf';else,s=sprintf('%.0f m',L);end,end

function local_write_report(path,r,root,mat_file,db_csv,scan_csv)
fid=fopen(path,'w','n','UTF-8');assert(fid>=0);cl=onCleanup(@()fclose(fid));
fprintf(fid,'# Adaptive Reduced Kirchhoff aperture: Stage 1 database\n\nStatus: **validation-only; production PE, BIE, Full Kirchhoff, and surface model unchanged**.\n\n');
fprintf(fid,'This first stage reuses the accepted deterministic A sweep (`K=0.10`) and K sweep (`A=0.02`), omitting their duplicate `A=0.02, K=0.10` point. It is not yet a complete 6-by-6 A-K matrix. The operator retains every accepted Full-Kirchhoff term and truncates only source points outside `|x_r-x_s|<L`.\n\n');
fprintf(fid,'L grid: `4, 8, 16, 32, 64, 128, Inf m`. Gate: `Reduced Kirchhoff-BIE complex L2 < %.6g`. `requires_full_K=true` means no finite tested L passed; `L_min=Inf` is then the complete Full-Kirchhoff reference.\n\n',r.complex_l2_gate);
fprintf(fid,'## Adaptive-aperture database\n\n| case | A m | K rad/m | A K | A K^2 1/m | L_min m | Ec at L_min | phase error rad | magnitude error | requires full K | PE-BIE Ec | FK-BIE Ec |\n|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|\n');
t=r.database;for ii=1:height(t),fprintf(fid,'| %s | %.6g | %.6g | %.6g | %.6g | %.6g | %.6g | %.6g | %.6g | %s | %.6g | %.6g |\n',t.case_name{ii},t.A_m(ii),t.K_radpm(ii),t.slope_AK(ii),t.curvature_AK2(ii),t.L_min_m(ii),t.Ec_at_Lmin(ii),t.phase_error_rad_at_Lmin(ii),t.magnitude_error_at_Lmin(ii),string(t.requires_full_K(ii)),t.PE_BIE_complex_l2(ii),t.FK_BIE_complex_l2(ii));end
fprintf(fid,'\n## Interpretation\n\n');
finite=~t.requires_full_K;fprintf(fid,'- Finite gate passes: `%d/%d`.\n',sum(finite),height(t));
if any(finite),fprintf(fid,'- Finite selected L range: `%.6g` to `%.6g m`.\n',min(t.L_min_m(finite)),max(t.L_min_m(finite)));else,fprintf(fid,'- No sampled case has a finite selected L under this absolute `1e-3` gate.\n');end
fprintf(fid,'- Every `L=Inf` reconstruction was checked directly against its accepted Full-Kirchhoff field; max relative discrepancy is `%.3g`.\n',max(t.full_equivalence_l2));
fprintf(fid,'- Because this Stage-1 database follows two one-dimensional sweeps rather than a complete matrix, it can diagnose monotonic/censored trends but cannot support an empirical law `L=f(A,K)` or a validated adaptive formula.\n\n');
fprintf(fid,'The next step should be chosen from the observed gate pattern. If finite L values vary across the one-dimensional sweeps, run only discriminating cells of the A-K matrix to test whether height, slope, or curvature is the better organizer. If the criterion is uniformly censored at 128 m, the current aperture definition is effectively global at this accuracy target and an adaptive-aperture PE boundary operator is not yet justified; low-rank/local-stationary kernel work is then more informative than fitting L.\n\n');
fprintf(fid,'Artifacts: `%s`, `%s`, `%s`; figures: `results/validation/pe_bie_full_kirchhoff_adaptive_aperture_stage1/figures/`.\n',local_rel(root,mat_file),local_rel(root,db_csv),local_rel(root,scan_csv));
end
function p=local_rel(root,p),p=strrep(p,[root filesep],'');p=strrep(p,filesep,'/');end
