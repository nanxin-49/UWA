function validation=validate_pe_point_source_error_budget_vertical(overrides)
%VALIDATE_PE_POINT_SOURCE_ERROR_BUDGET_VERTICAL Isolate aperture and sponge.
% Validation-only. Does not modify or substitute the production PE operator.

if nargin<1 || isempty(overrides), overrides=struct(); end
root=fileparts(fileparts(fileparts(mfilename('fullpath')))); addpath(root);
cfg=local_defaults(root); cfg=local_overrides(cfg,overrides); local_validate(cfg);
if ~exist(cfg.output_dir,'dir'), mkdir(cfg.output_dir); end

window_rows=repmat(local_window_row(),numel(cfg.widths_m)*2,1); q=0;
for ww=1:numel(cfg.widths_m)
    W=cfg.widths_m(ww); N=round(W/cfg.dx_target_m/2)*2; dx=W/N;
    x=(-W/2):dx:(W/2-dx); [X,Y]=meshgrid(x,x); [~,i0]=min(abs(x));
    source_fns={@virtual_point_source_initial_field_vertical,@weyl_point_source_initial_field_vertical};
    names=["spatial_truncated","weyl_subcell"];
    for mm=1:2
        source_cfg=struct('virtual_source_distance_m',cfg.s0_m,'source_green_amplitude',1/(4*pi), ...
            'c0',cfg.c0_mps,'xw',W,'yw',W,'nx',N,'ny',N,'x_tx',0,'y_tx',0, ...
            'weyl_regularization_np_per_m',cfg.weyl_regularization_np_per_m, ...
            'weyl_subcell_order',cfg.weyl_subcell_order);
        [psi0,source_meta]=source_fns{mm}(X,Y,cfg.frequency_hz,source_cfg);
        [psi_as,~]=exact_angular_spectrum_one_step_vertical(psi0,W,W, ...
            cfg.frequency_hz,cfg.c0_mps,cfg.L_m);
        k=2*pi*cfg.frequency_hz/cfg.c0_mps; R=cfg.s0_m+cfg.L_m;
        h_exact=exp(1i*k*R)/(4*pi*R);
        h_as=psi_as(i0,i0)*exp(1i*k*R);
        p=local_pe_params(cfg,W,N,source_fns{mm},0,cfg.sponge_ratio_values(1));
        pe=vertical_channel_model(p); h_pe=pe.H_direct_physical_f*exp(1i*k*cfg.s0_m);
        q=q+1; row=local_window_row(); row.source_mode=names(mm); row.width_m=W;
        row.nx=N; row.dx_m=dx; row.H_as=h_as; row.H_pe=h_pe; row.H_exact=h_exact;
        row.as_amplitude_error_db=20*log10(abs(h_as/h_exact));
        row.as_tl_error_db=-row.as_amplitude_error_db;
        row.as_phase_error_rad=angle(h_as*conj(h_exact));
        row.pe_amplitude_error_db=20*log10(abs(h_pe/h_exact));
        row.pe_tl_error_db=-row.pe_amplitude_error_db;
        row.pe_phase_error_rad=angle(h_pe*conj(h_exact));
        row.pe_as_relative_error=abs(h_pe-h_as)/abs(h_as);
        row.source_meta={source_meta}; window_rows(q)=row;
    end
end
window_table=struct2table(window_rows);

% Independent continuous Weyl reference. Variable substitutions remove the
% grazing 1/kz singularity, so this is the reliable free-space benchmark.
[~,continuous_weyl_meta]=weyl_point_source_reference_vertical( ...
    cfg.reference_offsets_m,cfg.s0_m+cfg.L_m,cfg.frequency_hz,cfg.c0_mps);
continuous_weyl_reliable=continuous_weyl_meta.relative_error_to_green<= ...
    cfg.continuous_weyl_relative_tolerance;

% A reliable baseline must be selected before sponge scanning. Weyl is used
% only if its largest-window analytic errors satisfy the declared gate.
largest=window_table(window_table.width_m==max(cfg.widths_m),:);
weyl=largest(largest.source_mode=="weyl_subcell",:);
spatial=largest(largest.source_mode=="spatial_truncated",:);
discrete_weyl_reliable=abs(weyl.as_tl_error_db)<=cfg.baseline_tl_tolerance_db && ...
    abs(weyl.as_phase_error_rad)<=cfg.baseline_phase_tolerance_rad;
baseline_source="spatial_truncated";

sponge_rows=repmat(local_sponge_row(),numel(cfg.widths_m)*numel(cfg.sponge_ratio_values)*numel(cfg.sponge_alpha_values),1);
q=0; source_fn=@virtual_point_source_initial_field_vertical;
if baseline_source=="weyl_subcell", source_fn=@weyl_point_source_initial_field_vertical; end
for ww=1:numel(cfg.widths_m)
    W=cfg.widths_m(ww); N=round(W/cfg.dx_target_m/2)*2;
    p0=local_pe_params(cfg,W,N,source_fn,0,cfg.sponge_ratio_values(1),'slice');
    out0=vertical_channel_model(p0);
    k=2*pi*cfg.frequency_hz/cfg.c0_mps;
    h0=out0.H_direct_physical_f*exp(1i*k*cfg.s0_m);
    for rr=1:numel(cfg.sponge_ratio_values)
        ratio=cfg.sponge_ratio_values(rr);
        e0=local_field_energy(out0.psifinal_xy,out0.x,out0.y,W,ratio, ...
            cfg.center_energy_radius_m);
        for aa=1:numel(cfg.sponge_alpha_values)
            q=q+1; alpha=cfg.sponge_alpha_values(aa);
            if alpha==0
                h=h0; e=e0;
            else
                p=local_pe_params(cfg,W,N,source_fn,alpha,ratio,'slice');
                out=vertical_channel_model(p);
                h=out.H_direct_physical_f*exp(1i*k*cfg.s0_m);
                e=local_field_energy(out.psifinal_xy,out.x,out.y,W,ratio, ...
                    cfg.center_energy_radius_m);
            end
            row=local_sponge_row(); row.source_mode=baseline_source;
            row.width_m=W; row.nx=N; row.dx_m=W/N; row.sponge_ratio=ratio;
            row.alpha_max_np_per_m=alpha; row.H_axis=h;
            row.axis_amplitude_change_db=20*log10(abs(h/h0));
            row.axis_tl_change_db=-row.axis_amplitude_change_db;
            row.axis_phase_change_rad=angle(h*conj(h0));
            row.center_energy=e.center; row.edge_energy=e.edge;
            row.total_energy=e.total; row.center_energy_fraction=e.center/e.total;
            row.edge_energy_fraction=e.edge/e.total;
            row.center_energy_change_db=10*log10(e.center/e0.center);
            row.edge_energy_change_db=10*log10(e.edge/e0.edge);
            row.total_energy_change_db=10*log10(e.total/e0.total);
            sponge_rows(q)=row;
        end
    end
end
sponge_table=struct2table(sponge_rows);

finite_window_error_db=spatial.as_tl_error_db;
sponge_max_increment_db=max(abs(sponge_table.axis_tl_change_db));
weyl_rows=window_table(window_table.source_mode=="weyl_subcell",:);
if height(weyl_rows)>=2
    window_delta_db=abs(weyl_rows.as_tl_error_db(end)-weyl_rows.as_tl_error_db(end-1));
else
    window_delta_db=Inf;
end
checks=table(["pe_as_frozen_operator";"continuous_weyl_reference"; ...
    "discrete_weyl_initialization";"window_converged"], ...
    [max(window_table.pe_as_relative_error);continuous_weyl_meta.relative_error_to_green; ...
    double(~discrete_weyl_reliable);window_delta_db], ...
    [cfg.pe_as_tolerance;cfg.continuous_weyl_relative_tolerance;0; ...
    cfg.window_convergence_tolerance_db], ...
    'VariableNames',{'check_name','value','limit'});
checks.passed=checks.value<=checks.limit;
validation=struct('schema_version','1.0.0','config',cfg,'window_table',window_table, ...
    'sponge_table',sponge_table,'baseline_source',baseline_source, ...
    'continuous_weyl_meta',continuous_weyl_meta, ...
    'continuous_weyl_reliable',continuous_weyl_reliable, ...
    'discrete_weyl_reliable',discrete_weyl_reliable, ...
    'finite_window_error_db',finite_window_error_db, ...
    'sponge_max_increment_db',sponge_max_increment_db,'checks',checks, ...
    'ready_for_bellhop',continuous_weyl_reliable);
validation.files=local_outputs(validation,root);
end

function cfg=local_defaults(root)
cfg=struct('output_dir',fullfile(root,'results','validation','pe_point_source_error_budget'), ...
    'frequency_hz',4000,'c0_mps',1500,'s0_m',10,'L_m',70, ...
    'widths_m',[24 32 48 64 80 100],'dx_target_m',0.25,'stepz_lamb',0.5, ...
    'weyl_subcell_order',8,'weyl_regularization_np_per_m',0, ...
    'reference_offsets_m',[0 0.5 1 2], ...
    'sponge_ratio_values',[0.10 0.12 0.15], ...
    'sponge_alpha_values',[0 0.025 0.05 0.1 0.15 0.3], ...
    'baseline_tl_tolerance_db',0.5,'baseline_phase_tolerance_rad',0.1, ...
    'window_convergence_tolerance_db',0.25,'pe_as_tolerance',1e-10, ...
    'continuous_weyl_relative_tolerance',1e-10,'center_energy_radius_m',2);
end

function cfg=local_overrides(cfg,o)
names=fieldnames(o); for ii=1:numel(names), if ~isfield(cfg,names{ii}), error('Unknown override: %s',names{ii}); end, cfg.(names{ii})=o.(names{ii}); end
end

function local_validate(cfg)
if max(cfg.widths_m)>100, error('Public PE validation domain is limited to 100 m.'); end
if max(round(cfg.widths_m/cfg.dx_target_m/2)*2)>2048, error('Grid exceeds the public PE limit.'); end
if ~any(cfg.sponge_alpha_values==0), error('sponge_alpha_values must include zero.'); end
end

function p=local_pe_params(cfg,W,N,source_fn,alpha,ratio,save_mode)
if nargin<7, save_mode='rx_only'; end
p=struct('f0',cfg.frequency_hz,'c0',cfg.c0_mps,'z_max',cfg.L_m,'z_tx',cfg.L_m, ...
    'z_rx',0,'xw',W,'yw',W,'nx',N,'ny',N,'x_tx',0,'y_tx',0,'x_rx',0,'y_rx',0, ...
    'stepz_lamb',cfg.stepz_lamb,'source_mode','custom_field_fn','source_field_fn',source_fn, ...
    'virtual_source_distance_m',cfg.s0_m,'source_green_amplitude',1/(4*pi), ...
    'weyl_regularization_np_per_m',cfg.weyl_regularization_np_per_m, ...
    'weyl_subcell_order',cfg.weyl_subcell_order,'sponge_ratio',ratio, ...
    'alpha_max_np_per_m',alpha,'env_mode','uniform','enable_surface_reflection',false, ...
    'enable_bubbles',false,'enforce_1_over_R',false,'show_figures',false, ...
    'save_mode',save_mode,'use_gpu',false);
end

function r=local_window_row()
r=struct('source_mode',"",'width_m',NaN,'nx',NaN,'dx_m',NaN,'H_as',complex(NaN), ...
    'H_pe',complex(NaN),'H_exact',complex(NaN),'as_amplitude_error_db',NaN, ...
    'as_tl_error_db',NaN,'as_phase_error_rad',NaN,'pe_amplitude_error_db',NaN, ...
    'pe_tl_error_db',NaN,'pe_phase_error_rad',NaN, ...
    'pe_as_relative_error',NaN,'source_meta',{{}});
end

function r=local_sponge_row()
r=struct('source_mode',"",'width_m',NaN,'nx',NaN,'dx_m',NaN, ...
    'sponge_ratio',NaN,'alpha_max_np_per_m',NaN,'H_axis',complex(NaN), ...
    'axis_amplitude_change_db',NaN,'axis_tl_change_db',NaN, ...
    'axis_phase_change_rad',NaN,'center_energy',NaN,'edge_energy',NaN, ...
    'total_energy',NaN,'center_energy_fraction',NaN,'edge_energy_fraction',NaN, ...
    'center_energy_change_db',NaN,'edge_energy_change_db',NaN, ...
    'total_energy_change_db',NaN);
end

function e=local_field_energy(psi_xy,x,y,W,ratio,center_radius)
[X,Y]=meshgrid(x,y); rho=hypot(X,Y);
center_mask=rho<=center_radius;
edge_start=W*(0.5-ratio);
edge_mask=abs(X)>=edge_start | abs(Y)>=edge_start;
dA=abs(x(2)-x(1))*abs(y(2)-y(1)); power_xy=abs(psi_xy).^2;
e=struct('center',sum(power_xy(center_mask))*dA, ...
    'edge',sum(power_xy(edge_mask))*dA,'total',sum(power_xy,'all')*dA);
end

function files=local_outputs(v,root)
out=v.config.output_dir; if ~exist(out,'dir'), mkdir(out); end
w=fullfile(out,'window_convergence.csv'); s=fullfile(out,'sponge_error_budget.csv');
c=fullfile(out,'checks.csv'); m=fullfile(out,'point_source_error_budget.mat');
f=fullfile(out,'point_source_error_budget.png'); report=fullfile(root,'reports','pe_point_source_error_budget_report.md');
writetable(v.window_table,w); writetable(v.sponge_table,s); writetable(v.checks,c);
fig=figure('Visible','off','Color','w','Position',[100 100 1100 460]); cleanup=onCleanup(@()close(fig));
subplot(1,2,1); hold on; for mode=unique(v.window_table.source_mode).', q=v.window_table(v.window_table.source_mode==mode,:); plot(q.width_m,q.as_tl_error_db,'o-','DisplayName',mode); end
yline(0); grid on; xlabel('window width (m)'); ylabel('AS-analytic TL (dB)'); legend('Location','best');
subplot(1,2,2); hold on; W=max(v.sponge_table.width_m); for ratio=unique(v.sponge_table.sponge_ratio).', q=v.sponge_table(v.sponge_table.width_m==W & v.sponge_table.sponge_ratio==ratio,:); plot(q.alpha_max_np_per_m,q.axis_tl_change_db,'o-','DisplayName',sprintf('ratio %.2f',ratio)); end
yline(0); grid on; xlabel('alpha max (Np/m)'); ylabel('sponge TL increment (dB)'); legend('Location','best');
exportgraphics(fig,f,'Resolution',180); clear cleanup
validation=v; schema_version=v.schema_version; save(m,'validation','schema_version','-v7.3');
fid=fopen(report,'w','n','UTF-8'); cl=onCleanup(@()fclose(fid));
fprintf(fid,'# PE point-source finite-window and sponge error budget\n\n');
fprintf(fid,'Status: `ready_for_bellhop=%s`; PE marching and default Gaussian source were frozen.\n\n',string(v.ready_for_bellhop));
fprintf(fid,'- Continuous Weyl reference relative error: `%.6g`; reliable: `%s`.\n',v.continuous_weyl_meta.relative_error_to_green,string(v.continuous_weyl_reliable));
fprintf(fid,'- Discrete FFT-Weyl initialization reliable at largest window: `%s`.\n',string(v.discrete_weyl_reliable));
spatial=v.window_table(v.window_table.source_mode=="spatial_truncated",:);
fprintf(fid,'- Largest-window spatial truncation: `dA=%.6g dB`, `dTL=%.6g dB`.\n', ...
    spatial.as_amplitude_error_db(end),spatial.as_tl_error_db(end));
fprintf(fid,'- Maximum sponge-only TL change across the full matrix: `%.6g dB`.\n',v.sponge_max_increment_db);
fprintf(fid,'- Bellhop rerun is permitted only when `ready_for_bellhop=true`.\n\n');
fprintf(fid,'## Definitions\n\n');
fprintf(fid,'- `axis_amplitude_change_db = 20 log10(|H_sponge|/|H_no_sponge|)`.\n');
fprintf(fid,'- `axis_tl_change_db = TL_sponge - TL_no_sponge = -axis_amplitude_change_db`.\n');
fprintf(fid,'- Center energy integrates `|psi|^2` over `rho <= %.3g m`; edge energy integrates over the ratio-defined sponge band.\n\n',v.config.center_energy_radius_m);
fprintf(fid,'## Full window x sponge ratio x alpha_max table\n\n');
fprintf(fid,'The complete numeric table is saved in [`sponge_error_budget.csv`](../results/validation/pe_point_source_error_budget/sponge_error_budget.csv).\n\n');
fprintf(fid,'| W (m) | ratio | alpha (Np/m) | dA axis (dB) | dTL axis (dB) | dphase (rad) | dE center (dB) | dE edge (dB) | dE total (dB) |\n');
fprintf(fid,'|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for ii=1:height(v.sponge_table)
    r=v.sponge_table(ii,:);
    fprintf(fid,'| %.0f | %.2f | %.3f | %.6f | %.6f | %.6f | %.6f | %.6f | %.6f |\n', ...
        r.width_m,r.sponge_ratio,r.alpha_max_np_per_m,r.axis_amplitude_change_db, ...
        r.axis_tl_change_db,r.axis_phase_change_rad,r.center_energy_change_db, ...
        r.edge_energy_change_db,r.total_energy_change_db);
end
fprintf(fid,'\n');
fprintf(fid,'![error budget](../results/validation/pe_point_source_error_budget/point_source_error_budget.png)\n'); clear cl
files=struct('mat',m,'window_csv',w,'sponge_csv',s,'checks_csv',c,'figure',f,'report',report);
end
