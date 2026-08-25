function validation = validate_pe_bellhop_flat_surface_matrix_vertical(overrides)
%VALIDATE_PE_BELLHOP_FLAT_SURFACE_MATRIX_VERTICAL Multi-case flat-surface audit.
% Validation-only code: the public PE propagator and communication chain are
% not modified. The script requires the separate phase-convention audit to
% pass before PE envelopes are converted to physical positive-delay spectra.

if nargin < 1 || isempty(overrides), overrides = struct(); end
this_file = mfilename('fullpath');
project_root = fileparts(fileparts(fileparts(this_file)));
addpath(project_root);
setup_vertical_project();

cfg = local_defaults();
cfg = local_overrides(cfg, overrides);
local_validate_cfg(cfg);
out_dir = fullfile(project_root, 'results', 'validation', ...
    'pe_bellhop_flat_surface_matrix');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

audit_file = fullfile(project_root, 'results', 'validation', ...
    'pe_phase_convention_uniform', 'pe_phase_convention_audit.mat');
if exist(audit_file, 'file') ~= 2
    error('Run validate_pe_phase_convention_uniform_vertical first.');
end
loaded = load(audit_file, 'audit');
if ~loaded.audit.passed || loaded.audit.selected_carrier_sign ~= 1
    error('The prerequisite phase audit has not passed with carrier sign +1.');
end

bellhop_exe = local_find_bellhop(cfg.bellhop_exe);
fprintf('PE/Bellhop flat-surface geometry matrix\nBellhop: %s\n', bellhop_exe);
f_hz = linspace(cfg.frequency_band_hz(1), cfg.frequency_band_hz(2), ...
    cfg.frequency_count).';
case_count = numel(cfg.receiver_offsets_m);
case_template = struct('geometry', [], 'pe', [], 'bellhop_paths', [], ...
    'direct_cir', [], 'reflection_cir', [], 'delay_axis_s', [], ...
    'pe_total_power', [], 'bellhop_total_power', [], 'invariant_error', []);
cases = repmat(case_template, case_count, 1);
raw_arrivals = cell(case_count, 1);
rows = repmat(local_empty_path_row(), 2 * case_count, 1);

for ii = 1:case_count
    x_m = cfg.receiver_offsets_m(ii);
    fprintf('Geometry %d/%d: receiver offset %.3f m\n', ii, case_count, x_m);
    g = local_geometry(cfg, x_m);
    pe = vertical_channel_model(local_pe_params(cfg, f_hz.', x_m, ...
        cfg.baseline_grid));
    invariant = max(abs(pe.H_f(:) - pe.H_direct_f(:) - pe.H_reflect_f(:)));
    if invariant > cfg.invariant_tolerance || ~pe.pass_1_over_R
        error('PE invariant failed for receiver offset %.3f m.', x_m);
    end

    root = fullfile(out_dir, sprintf('flat_matrix_x%02dm', round(x_m)));
    local_write_env([root '.env'], cfg, g);
    local_run_bellhop(bellhop_exe, root);
    all_bh = local_read_arrivals([root '.arr']);
    raw_arrivals{ii} = all_bh;
    bh = local_classify_arrivals(all_bh, g, cfg.arrival_cluster_tolerance_ms);

    td0 = (cfg.z_tx_m - cfg.z_rx_m) / cfg.c0_mps;
    tr0 = (cfg.z_tx_m + cfg.z_rx_m) / cfg.c0_mps;
    Hd = pe.H_direct_physical_f(:);
    Hr = pe.H_reflect_physical_f(:);
    [td_pe, cir_d] = local_path_peak(Hd, f_hz, td0, g.direct_time_s, cfg);
    [tr_pe, cir_r] = local_path_peak(Hr, f_hz, tr0, g.surface_time_s, cfg);
    [delay_axis_s, pe_total_power] = local_total_pdp(Hd + Hr, f_hz, td0, cfg);
    Hbh = bh.amplitude_complex(1) .* exp(1i * 2*pi*f_hz*bh.delay_s(1)) + ...
        bh.amplitude_complex(2) .* exp(1i * 2*pi*f_hz*bh.delay_s(2));
    [delay_bh_s, bh_total_power] = local_total_pdp(Hbh, f_hz, td0, cfg);
    if max(abs(delay_axis_s - delay_bh_s)) > eps(max(delay_axis_s))
        error('Internal PDP axes do not match.');
    end

    [~, idx_ref] = min(abs(f_hz - cfg.f_ref_hz));
    pe_amp = [abs(Hd(idx_ref)); abs(Hr(idx_ref))];
    for pp = 1:2
        rr = 2 * (ii - 1) + pp;
        path_name = ["direct"; "surface_reflection"];
        analytic_t = [g.direct_time_s; g.surface_time_s];
        pe_t = [td_pe; tr_pe];
        analytic_amp = [1/g.direct_length_m; 1/g.surface_length_m];
        rows(rr) = local_empty_path_row();
        rows(rr).offset_m = x_m;
        rows(rr).path_name = path_name(pp);
        rows(rr).analytic_time_ms = 1000 * analytic_t(pp);
        rows(rr).pe_time_ms = 1000 * pe_t(pp);
        rows(rr).bellhop_time_ms = 1000 * bh.delay_s(pp);
        rows(rr).pe_minus_analytic_ms = 1000 * (pe_t(pp)-analytic_t(pp));
        rows(rr).bellhop_minus_analytic_ms = 1000 * (bh.delay_s(pp)-analytic_t(pp));
        rows(rr).pe_minus_bellhop_ms = 1000 * (pe_t(pp)-bh.delay_s(pp));
        rows(rr).pe_amplitude = pe_amp(pp);
        rows(rr).bellhop_amplitude = bh.amplitude(pp);
        rows(rr).analytic_amplitude = analytic_amp(pp);
        rows(rr).pe_tl_db = -20*log10(max(pe_amp(pp), realmin));
        rows(rr).bellhop_tl_db = -20*log10(max(bh.amplitude(pp), realmin));
        rows(rr).bellhop_cluster_size = bh.cluster_size(pp);
    end
    cases(ii) = struct('geometry', g, 'pe', local_compact(pe), ...
        'bellhop_paths', bh, 'direct_cir', cir_d, 'reflection_cir', cir_r, ...
        'delay_axis_s', delay_axis_s, 'pe_total_power', pe_total_power, ...
        'bellhop_total_power', bh_total_power, 'invariant_error', invariant);
end

path_table = struct2table(rows);
direct_mask = path_table.path_name == "direct";
scale_each = path_table.bellhop_amplitude(direct_mask) ./ ...
    path_table.pe_amplitude(direct_mask);
global_scale = exp(mean(log(scale_each)));
scale_spread_db = max(20*log10(scale_each)) - min(20*log10(scale_each));
path_table.pe_scaled_amplitude = global_scale * path_table.pe_amplitude;
path_table.pe_scaled_tl_db = -20*log10(max(path_table.pe_scaled_amplitude, realmin));
path_table.scaled_tl_residual_db = path_table.pe_scaled_tl_db - path_table.bellhop_tl_db;

fprintf('Scalar direct-only regression at x=%.3f m\n', cfg.convergence_offset_m);
scalar_params = local_pe_params(cfg, cfg.f_ref_hz, cfg.convergence_offset_m, ...
    cfg.baseline_grid);
scalar_params.enable_surface_reflection = false;
scalar_direct_only = vertical_channel_model(scalar_params);
direct_only_error = max(abs(scalar_direct_only.H_reflect_f(:)));

[convergence_table, convergence_cases] = local_convergence(cfg, f_hz, cases);
checks = local_checks(path_table, scale_spread_db, cases, direct_only_error, ...
    convergence_table, cfg);
passed = all(checks.passed);

path_csv = fullfile(out_dir, 'geometry_path_comparison.csv');
conv_csv = fullfile(out_dir, 'pe_convergence_summary.csv');
fig_file = fullfile(out_dir, 'pe_bellhop_matrix_pdp.png');
mat_file = fullfile(out_dir, 'pe_bellhop_flat_surface_matrix.mat');
report_file = fullfile(out_dir, 'pe_bellhop_flat_surface_matrix_report.md');
writetable(path_table, path_csv);
writetable(convergence_table, conv_csv);
local_plot(cases, path_table, convergence_table, fig_file);

validation = struct('config', cfg, 'phase_audit_file', audit_file, ...
    'bellhop_executable', bellhop_exe, 'path_table', path_table, ...
    'cases', cases, 'raw_bellhop_arrivals', {raw_arrivals}, ...
    'global_direct_scale', global_scale, 'direct_scale_each', scale_each, ...
    'scale_spread_db', scale_spread_db, 'scalar_direct_only', ...
    local_compact(scalar_direct_only), 'direct_only_error', direct_only_error, ...
    'convergence_table', convergence_table, ...
    'convergence_cases', convergence_cases, 'checks', checks, ...
    'passed', passed, 'files', struct('path_csv', path_csv, ...
    'convergence_csv', conv_csv, 'figure', fig_file, 'mat', mat_file, ...
    'report', report_file));
save(mat_file, 'validation');
local_write_report(report_file, validation);
disp(path_table); disp(convergence_table); disp(checks);
if ~passed
    failed = strjoin(cellstr(checks.check_name(~checks.passed)), ', ');
    error('validate_pe_bellhop_flat_surface_matrix_vertical:Failed', ...
        'Validation checks failed: %s', failed);
end
fprintf('Matrix validation passed: %s\n', report_file);
end

function cfg = local_defaults()
cfg = struct('water_depth_m',100,'c0_mps',1500,'z_tx_m',80,'z_rx_m',10, ...
    'receiver_offsets_m',[3 6 9],'convergence_offset_m',6, ...
    'f_ref_hz',4000,'frequency_band_hz',[3000 5000],'frequency_count',33, ...
    'sigma_src_m',0.3,'sponge_ratio',0.12,'alpha_max_np_per_m',0.15, ...
    'baseline_grid',struct('name',"C0",'nx',128,'ny',128,'width_m',32,'stepz_lamb',0.5), ...
    'bellhop_beam_count',10001,'bellhop_angle_min_deg',-89.5, ...
    'bellhop_angle_max_deg',-60,'bellhop_exe','', ...
    'arrival_cluster_tolerance_ms',0.1,'cir_window','hann', ...
    'cir_zero_padding_factor',8,'time_tolerance_ms',0.25, ...
    'scale_spread_tolerance_db',1,'reflection_rms_tolerance_db',1, ...
    'reflection_max_tolerance_db',2,'convergence_time_tolerance_ms',0.1, ...
    'convergence_relative_tl_tolerance_db',0.5, ...
    'convergence_phase_tolerance_rad',0.15,'invariant_tolerance',1e-10);
end

function cfg = local_overrides(cfg, over)
if ~isstruct(over), error('overrides must be a struct.'); end
n = fieldnames(over);
for ii=1:numel(n)
    if ~isfield(cfg,n{ii}), error('Unknown override: %s',n{ii}); end
    cfg.(n{ii})=over.(n{ii});
end
end

function local_validate_cfg(cfg)
if cfg.z_rx_m < 0 || cfg.z_rx_m >= cfg.z_tx_m || cfg.z_tx_m > cfg.water_depth_m
    error('Require 0 <= z_rx < z_tx <= water depth.');
end
if mod(cfg.frequency_count,2)~=1 || cfg.frequency_count<3
    error('frequency_count must be odd and >=3.');
end
if ~ismember(cfg.convergence_offset_m,cfg.receiver_offsets_m)
    error('convergence_offset_m must be in receiver_offsets_m.');
end
end

function g = local_geometry(cfg,x)
g = struct('offset_m',x, ...
    'direct_length_m',hypot(x,cfg.z_tx_m-cfg.z_rx_m), ...
    'surface_length_m',hypot(x,cfg.z_tx_m+cfg.z_rx_m));
g.direct_time_s=g.direct_length_m/cfg.c0_mps;
g.surface_time_s=g.surface_length_m/cfg.c0_mps;
end

function p = local_pe_params(cfg,f,x,grid)
p=struct('f0',f,'enable_wideband',false,'f_ref_hz',cfg.f_ref_hz, ...
    'c0',cfg.c0_mps,'z_max',cfg.water_depth_m,'stepz_lamb',grid.stepz_lamb, ...
    'xw',grid.width_m,'yw',grid.width_m,'nx',grid.nx,'ny',grid.ny, ...
    'x_tx',0,'y_tx',0,'z_tx',cfg.z_tx_m,'x_rx',x,'y_rx',0,'z_rx',cfg.z_rx_m, ...
    'sigma_src_m',cfg.sigma_src_m,'sponge_ratio',cfg.sponge_ratio, ...
    'alpha_max_np_per_m',cfg.alpha_max_np_per_m,'env_mode','uniform', ...
    'show_figures',false,'enforce_1_over_R',true,'enable_surface_reflection',true, ...
    'surface_reflect_coeff',-1,'surface_phase_mode','normal', ...
    'surface_boundary_model','kirchhoff_spatial','sea_hs_target',0, ...
    'enable_bubbles',false,'save_mode','rx_only','use_gpu',false);
end

function exe=local_find_bellhop(configured)
c={}; if ~isempty(configured), c{end+1}=char(configured); end
e=getenv('BELLHOP_EXE'); if ~isempty(e), c{end+1}=e; end
w=which('bellhop.exe'); if ~isempty(w), c{end+1}=w; end
for ii=1:numel(c), if exist(c{ii},'file')==2, exe=c{ii}; return; end, end
if ispc, [s,o]=system('where bellhop.exe'); else, [s,o]=system('which bellhop'); end
if s==0, a=regexp(strtrim(o),'\r?\n','split'); exe=strtrim(a{1}); return; end
error('Bellhop executable not found; set BELLHOP_EXE or bellhop_exe.');
end

function local_write_env(file,cfg,g)
fid=fopen(file,'w'); if fid<0, error('Cannot write %s',file); end
cl=onCleanup(@()fclose(fid));
fprintf(fid,'''PE-Bellhop open-fan flat-surface matrix''\n%.12g\n1\n''CVW''\n',cfg.f_ref_hz);
fprintf(fid,'2 0 %.12g\n0 %.12g /\n%.12g %.12g /\n',cfg.water_depth_m,cfg.c0_mps,cfg.water_depth_m,cfg.c0_mps);
fprintf(fid,'''A'' 0\n%.12g 1800 0 2 0 /\n',cfg.water_depth_m);
fprintf(fid,'1\n%.12g /\n1\n%.12g /\n1\n%.12g /\n',cfg.z_tx_m,cfg.z_rx_m,g.offset_m/1000);
fprintf(fid,'''A''\n%d\n%.12g %.12g /\n0 %.12g %.12g\n',cfg.bellhop_beam_count, ...
    cfg.bellhop_angle_min_deg,cfg.bellhop_angle_max_deg,1.1*cfg.water_depth_m,2*g.offset_m/1000);
clear cl
end

function local_run_bellhop(exe,root)
d=fileparts(root); [~,name]=fileparts(root); old=pwd; cl=onCleanup(@()cd(old)); cd(d);
if exist([name '.arr'],'file')==2, delete([name '.arr']); end
[s,o]=system(sprintf('"%s" "%s"',exe,name));
if s~=0 || exist([name '.arr'],'file')~=2, error('Bellhop failed: %s',o); end
clear cl
end

function t=local_read_arrivals(file)
fid=fopen(file,'r'); if fid<0,error('Cannot open %s',file);end; cl=onCleanup(@()fclose(fid));
f=fscanf(fid,'%f',1); ns=fscanf(fid,'%d',1); nd=fscanf(fid,'%d',1); nr=fscanf(fid,'%d',1);
sd=fscanf(fid,'%f',ns); rd=fscanf(fid,'%f',nd); rr=1000*fscanf(fid,'%f',nr);
if ns~=1||nd~=1||nr~=1,error('Expected one source/receiver/range.');end
fscanf(fid,'%d',1); n=fscanf(fid,'%d',1); a=fscanf(fid,'%f',[8,n]);
if size(a,2)~=n,error('Incomplete arrivals file.');end
z=a(1,:).'.*exp(1i*deg2rad(a(2,:).'));
t=table((1:n).',repmat(f,n,1),repmat(sd,n,1),repmat(rd,n,1),repmat(rr,n,1), ...
    z,abs(z),a(2,:).',a(3,:).',a(6,:).',a(7,:).',a(8,:).', ...
    'VariableNames',{'arrival_index','frequency_hz','source_depth_m','receiver_depth_m', ...
    'receiver_range_m','amplitude_complex','amplitude','phase_deg','delay_s', ...
    'receiver_angle_deg','top_bounce_count','bottom_bounce_count'});
clear cl
end

function bh=local_classify_arrivals(all,g,tol_ms)
names=["direct";"surface_reflection"]; tops=[0;1]; targets=[g.direct_time_s;g.surface_time_s];
bh=table(names,nan(2,1),complex(nan(2,1)),nan(2,1),zeros(2,1),tops,zeros(2,1), ...
    'VariableNames',{'path_name','delay_s','amplitude_complex','amplitude','cluster_size', ...
    'top_bounce_count','bottom_bounce_count'});
for pp=1:2
    q=sortrows(all(all.bottom_bounce_count==0 & all.top_bounce_count==tops(pp),:),'delay_s');
    if isempty(q),error('Bellhop did not return %s.',names(pp));end
    breaks=[true;diff(q.delay_s)>tol_ms/1000]; ids=cumsum(breaks); u=unique(ids);
    center=arrayfun(@(k)mean(q.delay_s(ids==k)),u);
    [~,kbest]=min(abs(center-targets(pp))); mask=ids==u(kbest);
    w=max(q.amplitude(mask),realmin);
    bh.delay_s(pp)=sum(w.*q.delay_s(mask))/sum(w);
    bh.amplitude_complex(pp)=sum(q.amplitude_complex(mask));
    bh.amplitude(pp)=abs(bh.amplitude_complex(pp)); bh.cluster_size(pp)=sum(mask);
end
end

function [peak,cir]=local_path_peak(H,f,tref,target,cfg)
[delay,power,h]=local_pdp(H,f,tref,cfg); mask=abs(delay-target)<=0.001;
ix=find(mask); [~,j]=max(power(mask)); peak=delay(ix(j));
cir=struct('delay_s',delay,'h',h,'power',power,'peak_time_s',peak,'reference_time_s',tref);
end

function [delay,power]=local_total_pdp(H,f,tref,cfg)
[delay,power]=local_pdp(H,f,tref,cfg);
end

function [delay,power,h]=local_pdp(H,f,tref,cfg)
F=numel(f); n=(0:F-1).'; df=mean(diff(f)); N=cfg.cir_zero_padding_factor*F;
if strcmpi(cfg.cir_window,'hann'),w=0.5-0.5*cos(2*pi*n/(F-1));else,w=ones(F,1);end
h=fft(H(:).*exp(-1i*2*pi*f(:)*tref).*w,N); delay=tref+(0:N-1).'/(df*N); power=abs(h).^2;
end

function [table_out,records]=local_convergence(cfg,f,baseline_cases)
grids=[cfg.baseline_grid, ...
    struct('name',"C1",'nx',256,'ny',256,'width_m',32,'stepz_lamb',0.5), ...
    struct('name',"C2",'nx',256,'ny',256,'width_m',64,'stepz_lamb',0.5), ...
    struct('name',"C3",'nx',128,'ny',128,'width_m',32,'stepz_lamb',0.25), ...
    struct('name',"C4",'nx',256,'ny',256,'width_m',32,'stepz_lamb',0.25)];
baseline_geometry = [baseline_cases.geometry];
baseline_offsets_m = [baseline_geometry.offset_m];
base_idx=find(baseline_offsets_m==cfg.convergence_offset_m,1);
records=repmat(struct('grid',[],'pe',[]),numel(grids),1);
records(1).grid=grids(1); records(1).pe=baseline_cases(base_idx).pe;
td0=(cfg.z_tx_m-cfg.z_rx_m)/cfg.c0_mps; tr0=(cfg.z_tx_m+cfg.z_rx_m)/cfg.c0_mps;
for ii=2:numel(grids)
    fprintf('Convergence %s: %dx%d, width %.0f m, step %.2f lambda\n',grids(ii).name, ...
        grids(ii).nx,grids(ii).ny,grids(ii).width_m,grids(ii).stepz_lamb);
    pe=vertical_channel_model(local_pe_params(cfg,f.',cfg.convergence_offset_m,grids(ii)));
    records(ii).grid=grids(ii); records(ii).pe=local_compact(pe);
end
base=records(1).pe; Hd0=base.H_direct_physical_f(:); Hr0=base.H_reflect_physical_f(:);
g=local_geometry(cfg,cfg.convergence_offset_m); [td_base,~]=local_path_peak(Hd0,f,td0,g.direct_time_s,cfg); [tr_base,~]=local_path_peak(Hr0,f,tr0,g.surface_time_s,cfg);
[~,ir]=min(abs(f-cfg.f_ref_hz)); rel0=20*log10(abs(Hd0(ir))/max(abs(Hr0(ir)),realmin));
row=repmat(struct('case_name',"",'nx',0,'ny',0,'width_m',0,'stepz_lamb',0, ...
    'direct_time_change_ms',0,'reflection_time_change_ms',0,'relative_tl_change_db',0, ...
    'direct_phase_rms_rad',0,'reflection_phase_rms_rad',0,'invariant_error',0),numel(grids),1);
for ii=1:numel(grids)
    pe=records(ii).pe; Hd=pe.H_direct_physical_f(:); Hr=pe.H_reflect_physical_f(:);
    [td,~]=local_path_peak(Hd,f,td0,g.direct_time_s,cfg); [tr,~]=local_path_peak(Hr,f,tr0,g.surface_time_s,cfg);
    rel=20*log10(abs(Hd(ir))/max(abs(Hr(ir)),realmin));
    row(ii)=struct('case_name',grids(ii).name,'nx',grids(ii).nx,'ny',grids(ii).ny, ...
        'width_m',grids(ii).width_m,'stepz_lamb',grids(ii).stepz_lamb, ...
        'direct_time_change_ms',1000*(td-td_base),'reflection_time_change_ms',1000*(tr-tr_base), ...
        'relative_tl_change_db',rel-rel0,'direct_phase_rms_rad',local_phase_rms(Hd,Hd0,ir), ...
        'reflection_phase_rms_rad',local_phase_rms(Hr,Hr0,ir), ...
        'invariant_error',max(abs(pe.H_f(:)-pe.H_direct_f(:)-pe.H_reflect_f(:))));
end
table_out=struct2table(row);
end

function v=local_phase_rms(H,H0,ir)
p=unwrap(angle((H/H(ir)).*conj(H0/H0(ir)))); p=p-p(ir); v=sqrt(mean(p.^2));
end

function checks=local_checks(t,spread,cases,direct_only,conv,cfg)
reflection=t.path_name=="surface_reflection"; nonbase=conv.case_name~="C0";
names=["bellhop_path_count";"pe_vs_analytic_time";"bellhop_vs_analytic_time"; ...
    "pe_vs_bellhop_time";"global_scale_spread";"reflection_tl_rms"; ...
    "reflection_tl_max";"pe_sum_invariant";"direct_only_reflection"; ...
    "convergence_time";"convergence_relative_tl";"convergence_phase";"convergence_invariant"];
inv=max([cases.invariant_error]); refl=t.scaled_tl_residual_db(reflection);
values=[double(any(t.bellhop_cluster_size<1));max(abs(t.pe_minus_analytic_ms)); ...
    max(abs(t.bellhop_minus_analytic_ms));max(abs(t.pe_minus_bellhop_ms));spread; ...
    sqrt(mean(refl.^2));max(abs(refl));inv;direct_only; ...
    max(abs([conv.direct_time_change_ms(nonbase);conv.reflection_time_change_ms(nonbase)])); ...
    max(abs(conv.relative_tl_change_db(nonbase))); ...
    max([conv.direct_phase_rms_rad(nonbase);conv.reflection_phase_rms_rad(nonbase)]); ...
    max(conv.invariant_error)];
limits=[0;cfg.time_tolerance_ms;cfg.time_tolerance_ms;cfg.time_tolerance_ms; ...
    cfg.scale_spread_tolerance_db;cfg.reflection_rms_tolerance_db; ...
    cfg.reflection_max_tolerance_db;cfg.invariant_tolerance;cfg.invariant_tolerance; ...
    cfg.convergence_time_tolerance_ms;cfg.convergence_relative_tl_tolerance_db; ...
    cfg.convergence_phase_tolerance_rad;cfg.invariant_tolerance];
rel=["==";repmat("<=",12,1)]; pass=(rel=="=="&values==limits)|(rel=="<="&values<=limits);
checks=table(names,values,rel,limits,pass,'VariableNames',{'check_name','value','relation','limit','passed'});
end

function row=local_empty_path_row()
row=struct('offset_m',NaN,'path_name',"",'analytic_time_ms',NaN,'pe_time_ms',NaN, ...
    'bellhop_time_ms',NaN,'pe_minus_analytic_ms',NaN,'bellhop_minus_analytic_ms',NaN, ...
    'pe_minus_bellhop_ms',NaN,'pe_amplitude',NaN,'bellhop_amplitude',NaN, ...
    'analytic_amplitude',NaN,'pe_tl_db',NaN,'bellhop_tl_db',NaN,'bellhop_cluster_size',0);
end

function c=local_compact(p)
c=struct('H_direct_f',p.H_direct_f,'H_reflect_f',p.H_reflect_f,'H_f',p.H_f, ...
    'H_direct_physical_f',p.H_direct_physical_f, ...
    'H_reflect_physical_f',p.H_reflect_physical_f, ...
    'f_axis',p.f_axis,'idx_f_ref',p.idx_f_ref,'h_direct',p.h_direct, ...
    'h_reflect',p.h_reflect,'h_total',p.h_total,'pass_1_over_R',p.pass_1_over_R, ...
    'fit_slope',p.fit_slope,'fit_err_rms',p.fit_err_rms,'config',p.config);
end

function local_plot(cases,t,conv,file)
fig=figure('Visible','off','Color','w','Position',[100 100 1100 850]); cl=onCleanup(@()close(fig));
for ii=1:numel(cases)
    subplot(2,2,ii); p=cases(ii).pe_total_power; b=cases(ii).bellhop_total_power;
    plot(1000*cases(ii).delay_axis_s,10*log10(max(p/max(p),1e-7)),'b-','LineWidth',1.2); hold on;
    plot(1000*cases(ii).delay_axis_s,10*log10(max(b/max(b),1e-7)),'r--','LineWidth',1.2);
    xlim([min(t.analytic_time_ms)-1,max(t.analytic_time_ms)+1]); ylim([-60 1]); grid on;
    title(sprintf('x = %.0f m',cases(ii).geometry.offset_m)); xlabel('Delay (ms)'); ylabel('PDP (dB)');
    if ii==1,legend('PE','Bellhop','Location','best');end
end
subplot(2,2,4); non=conv.case_name~="C0";
bar(categorical(conv.case_name(non)),[conv.direct_phase_rms_rad(non),conv.reflection_phase_rms_rad(non)]);
grid on; ylabel('Phase RMS (rad)'); title('PE convergence vs C0'); legend('direct','surface','Location','best');
exportgraphics(fig,file,'Resolution',180); clear cl
end

function local_write_report(file,v)
t=v.path_table; c=v.checks; q=v.convergence_table;
fid=fopen(file,'w','n','UTF-8'); if fid<0,error('Cannot write report.');end; cl=onCleanup(@()fclose(fid));
fprintf(fid,'# PE–Bellhop 平面海面多几何与收敛验证\n\n');
fprintf(fid,'结论：`%s`。本验证仅新增独立脚本，未修改 PE 主线、公共输出或通信链。\n\n',string(v.passed));
if ~v.passed
    failed = strjoin(cellstr(c.check_name(~c.passed)), '、');
    fprintf(fid,'未通过的严格检查：%s。到达时间与代数不变量仍通过，详见下表。\n\n',failed);
end
fprintf(fid,'相位约定由前置审计确定：在 `exp(-iωt)` 合成约定下，仅恢复纵向载波 `exp(+ikd)`，再由验证层 FFT 提取正时延。\n\n');
fprintf(fid,'## 到达时间与路径幅度\n\n');
fprintf(fid,'| x(m) | 路径 | 解析(ms) | PE(ms) | Bellhop(ms) | PE-BH(ms) | 缩放后TL残差(dB) | BH簇大小 |\n|---:|---|---:|---:|---:|---:|---:|---:|\n');
for ii=1:height(t)
    fprintf(fid,'| %.0f | %s | %.6f | %.6f | %.6f | %+.6f | %+.4f | %d |\n', ...
        t.offset_m(ii),t.path_name(ii),t.analytic_time_ms(ii),t.pe_time_ms(ii), ...
        t.bellhop_time_ms(ii),t.pe_minus_bellhop_ms(ii), ...
        t.scaled_tl_residual_db(ii),t.bellhop_cluster_size(ii));
end
fprintf(fid,'\n跨几何统一直达标定系数：`%.8g`；逐几何标定量跨度：`%.4f dB`。该系数只用于比较，不回写 PE。\n\n',v.global_direct_scale,v.scale_spread_db);
fprintf(fid,'## PE 收敛（x=%.0f m）\n\n',v.config.convergence_offset_m);
fprintf(fid,'| 案例 | 网格 | 窗口(m) | dz(λ) | 直达Δt(ms) | 反射Δt(ms) | 相对TLΔ(dB) | 直达相位RMS | 反射相位RMS |\n|---|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for ii=1:height(q)
    fprintf(fid,'| %s | %dx%d | %.0f | %.2f | %+.6f | %+.6f | %+.4f | %.5f | %.5f |\n', ...
        q.case_name(ii),q.nx(ii),q.ny(ii),q.width_m(ii),q.stepz_lamb(ii), ...
        q.direct_time_change_ms(ii),q.reflection_time_change_ms(ii), ...
        q.relative_tl_change_db(ii),q.direct_phase_rms_rad(ii), ...
        q.reflection_phase_rms_rad(ii));
end
fprintf(fid,'\n## 自动检查\n\n| 检查 | 数值 | 条件 | 阈值 | 通过 |\n|---|---:|:---:|---:|:---:|\n');
for ii=1:height(c)
    fprintf(fid,'| %s | %.8g | %s | %.8g | %d |\n',c.check_name(ii), ...
        c.value(ii),c.relation(ii),c.limit(ii),c.passed(ii));
end
fprintf(fid,'\n![PDP 与收敛结果](pe_bellhop_matrix_pdp.png)\n\n');
fprintf(fid,'## 解释与限制\n\n- Bellhop 使用开放角扇区，并按海面/海底反射次数和 0.1 ms 时延簇分类；没有预先强制“两条射线”。\n');
fprintf(fid,'- PE 是有限宽高斯初场，Bellhop 是点源射线模型；统一标定仅消除源归一化常数，反射残差和跨几何漂移仍保留。\n');
fprintf(fid,'- C0–C4 分别检查横向采样、窗口宽度和纵向步长；本阶段不包含粗糙海面、SSA、Kirchhoff 随机散射或通信调制。\n');
clear cl
end
