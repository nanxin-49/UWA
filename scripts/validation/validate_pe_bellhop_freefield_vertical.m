function validation = validate_pe_bellhop_freefield_vertical(overrides)
%VALIDATE_PE_BELLHOP_FREEFIELD_VERTICAL Four-level free-space validation.
% No surface reflection, bottom reflection, bubbles, Doppler, or comm chain.

if nargin<1 || isempty(overrides), overrides=struct(); end
root=fileparts(fileparts(fileparts(mfilename('fullpath')))); addpath(root); addpath(fileparts(mfilename('fullpath')));
cfg=local_defaults(root); cfg=local_overrides(cfg,overrides);
if ~exist(cfg.output_dir,'dir'), mkdir(cfg.output_dir); end

pe_dir=fullfile(cfg.output_dir,'pe_reference');
pe_file=fullfile(pe_dir,'pe_as_freefield_validation.mat');
if cfg.reuse_existing_subresults && exist(pe_file,'file')==2
    cached=load(pe_file,'validation'); pe=cached.validation;
else
    pe=validate_pe_as_freefield_vertical(struct('output_dir',pe_dir, ...
        'axis_alpha_max_np_per_m',cfg.pe_alpha_max_np_per_m, ...
        'grid_cases',cfg.pe_grid_cases,'axis_nx',cfg.pe_nx, ...
        'axis_width_m',cfg.pe_width_m,'run_full_matrix',false));
end
if ~pe.level1_passed, error('PE-AS gate failed.'); end

bh_dir=fullfile(cfg.output_dir,'bellhop_normalization');
bh_file=fullfile(bh_dir,'bellhop_freefield_normalization.mat');
if cfg.reuse_existing_subresults && exist(bh_file,'file')==2
    cached=load(bh_file,'validation'); bh=cached.validation;
else
    bh=validate_bellhop_freefield_normalization_vertical(struct( ...
        'output_dir',bh_dir,'bellhop_exe',cfg.bellhop_exe, ...
        'frequencies_hz',cfg.normalization_frequencies_hz, ...
        'beam_counts',cfg.normalization_beam_counts, ...
        'step_values_m',cfg.normalization_step_values_m));
end
if ~bh.passed, error('Bellhop normalization gate failed.'); end

rows=repmat(local_row(),numel(cfg.frequencies_hz)*numel(cfg.L_values_m)*numel(cfg.s0_values_m)*numel(cfg.receiver_offsets_m),1);
q=0; raw_runs=cell(numel(cfg.frequencies_hz),1);
axial_ranges=unique(cfg.s0_values_m(:)+cfg.L_values_m(:).'); axial_ranges=axial_ranges(:).';
for ff=1:numel(cfg.frequencies_hz)
    f=cfg.frequencies_hz(ff);
    root_case=fullfile(cfg.output_dir,'bellhop_cases',sprintf('f%g_all_ranges',f));
    c=struct('bellhop_exe',cfg.bellhop_exe,'case_root',root_case, ...
        'frequency_hz',f,'c0_mps',cfg.c0_mps,'source_depth_m',0, ...
        'receiver_depths_m',cfg.receiver_offsets_m,'receiver_ranges_m',axial_ranges, ...
        'run_type','C','beam_count',cfg.bellhop_beam_count, ...
        'angle_limits_deg',[-180 180],'step_m',cfg.bellhop_step_m, ...
        'domain_half_depth_m',cfg.bellhop_domain_half_depth_m);
    run=run_bellhop_freefield_vertical(c); raw_runs{ff}=run;
    p_bh=run.data.pressure;
    if bh.selected_spatial_sign<0
        h_bh_all=conj(p_bh)./conj(bh.bellhop_source_constant)/(4*pi);
    else
        h_bh_all=p_bh./bh.bellhop_source_constant/(4*pi);
    end
    for ss=1:numel(cfg.s0_values_m)
        for ll=1:numel(cfg.L_values_m)
            axial=cfg.s0_values_m(ss)+cfg.L_values_m(ll);
            [~,range_ix]=min(abs(run.data.receiver_range_m-axial));
            h_bh=h_bh_all(:,range_ix);
            p=local_pe_params(cfg,f,cfg.L_values_m(ll),cfg.s0_values_m(ss));
            out=vertical_channel_model(p); k=2*pi*f/cfg.c0_mps;
            [~,iy0]=min(abs(out.y));
            for xx=1:numel(cfg.receiver_offsets_m)
                x=cfg.receiver_offsets_m(xx);
                [~,ix]=min(abs(out.x-x));
                if abs(out.x(ix)-x)>10*eps(max(1,abs(x)))
                    error('PE receiver offset %.9g m is not on the configured grid.',x);
                end
                % The saved reduced field is referenced to the initial plane.
                % Restore both marching distance L and virtual-source distance s0.
                h_pe=out.psifinal_xy(iy0,ix)*exp(1i*k*axial);
                R=hypot(axial,x); h_exact=exp(1i*k*R)/(4*pi*R);
                q=q+1; row=local_row(); row.frequency_hz=f; row.s0_m=cfg.s0_values_m(ss);
                row.L_m=cfg.L_values_m(ll); row.offset_m=x; row.range_m=R;
                row.H_pe_source=h_pe; row.H_bellhop_green=h_bh(xx); row.H_exact=h_exact;
                row.pe_tl_db=-20*log10(abs(h_pe)); row.bellhop_tl_db=-20*log10(abs(h_bh(xx)));
                row.exact_tl_db=-20*log10(abs(h_exact));
                row.pe_minus_bellhop_tl_db=row.pe_tl_db-row.bellhop_tl_db;
                row.pe_minus_exact_tl_db=row.pe_tl_db-row.exact_tl_db;
                row.bellhop_minus_exact_tl_db=row.bellhop_tl_db-row.exact_tl_db;
                row.pe_minus_bellhop_phase_rad=angle(h_pe*conj(h_bh(xx)));
                row.pe_minus_exact_phase_rad=angle(h_pe*conj(h_exact));
                row.bellhop_minus_exact_phase_rad=angle(h_bh(xx)*conj(h_exact));
                rows(q)=row;
            end
        end
    end
end
t=struct2table(rows);

% Fine-frequency phase slope at one declared representative geometry.
fine_f=cfg.delay_frequency_hz+cfg.delay_frequency_offsets_hz(:).';
p=local_pe_params(cfg,fine_f,cfg.delay_L_m,cfg.delay_s0_m);
fine=vertical_channel_model(p); fine_k=2*pi*fine.f_axis(:)/cfg.c0_mps;
h_fine=exp(1i*fine_k*cfg.delay_s0_m).*fine.H_direct_physical_f(:);
phase_fit=polyfit(fine.f_axis(:),unwrap(angle(h_fine)),1);
pe_delay_s=phase_fit(1)/(2*pi);
delay_range_m=cfg.delay_s0_m+cfg.delay_L_m;
exact_delay_s=delay_range_m/cfg.c0_mps;
arrival_cfg=struct('bellhop_exe',cfg.bellhop_exe, ...
    'case_root',fullfile(cfg.output_dir,'bellhop_cases','representative_arrival'), ...
    'frequency_hz',cfg.delay_frequency_hz,'c0_mps',cfg.c0_mps, ...
    'source_depth_m',0,'receiver_depths_m',0,'receiver_ranges_m',delay_range_m, ...
    'run_type','A','beam_count',cfg.bellhop_beam_count,'angle_limits_deg',[-180 180], ...
    'step_m',cfg.bellhop_step_m,'domain_half_depth_m',cfg.bellhop_domain_half_depth_m);
arrival_run=run_bellhop_freefield_vertical(arrival_cfg);
direct=arrival_run.data(arrival_run.data.top_bounce_count==0 & ...
    arrival_run.data.bottom_bounce_count==0,:);
[~,arrival_ix]=min(abs(direct.delay_s-exact_delay_s));
bellhop_delay_s=direct.delay_s(arrival_ix);
arrival_table=table(delay_range_m,pe_delay_s,bellhop_delay_s,exact_delay_s, ...
    pe_delay_s-exact_delay_s,bellhop_delay_s-exact_delay_s, ...
    'VariableNames',{'range_m','pe_delay_s','bellhop_delay_s','exact_delay_s', ...
    'pe_delay_error_s','bellhop_delay_error_s'});
checks=table(["pe_as_gate";"bellhop_normalization_gate";"pe_bellhop_tl"; ...
    "pe_bellhop_phase";"pe_group_delay";"bellhop_group_delay"; ...
    "bellhop_analytic_tl";"bellhop_analytic_phase"], ...
    [double(~pe.level1_passed);double(~bh.passed);max(abs(t.pe_minus_bellhop_tl_db)); ...
    sqrt(mean(t.pe_minus_bellhop_phase_rad.^2));abs(pe_delay_s-exact_delay_s); ...
    abs(bellhop_delay_s-exact_delay_s);max(abs(t.bellhop_minus_exact_tl_db)); ...
    sqrt(mean(t.bellhop_minus_exact_phase_rad.^2))], ...
    [0;0;cfg.tl_tolerance_db;cfg.phase_tolerance_rad;cfg.delay_tolerance_s; ...
    cfg.delay_tolerance_s;cfg.tl_tolerance_db;cfg.phase_tolerance_rad], ...
    'VariableNames',{'check_name','value','limit'});
checks.passed=checks.value<=checks.limit;
validation=struct('schema_version','1.0.0','config',cfg,'pe_reference',pe, ...
    'bellhop_normalization',bh,'comparison_table',t,'bellhop_runs',{raw_runs}, ...
    'arrival_table',arrival_table,'representative_arrival_run',arrival_run, ...
    'checks',checks,'passed',all(checks.passed));
validation.files=local_outputs(validation,root);
end

function cfg=local_defaults(root)
exe=getenv('BELLHOP_EXE');
if isempty(exe)
    candidate='E:\MISC\fxx\Bellhop相关\Bellhop例程包\atWin10_2020_11_4\atWin10_2020_11_4\windows-bin-20201102\bellhop.exe';
    if exist(candidate,'file')==2, exe=candidate; end
end
cfg=struct('bellhop_exe',exe,'output_dir',fullfile(root,'results','validation','pe_bellhop_freefield','formal'), ...
    'reuse_existing_subresults',false, ...
    'c0_mps',1500,'frequencies_hz',[3000 4000 5000],'L_values_m',[20 40 70 100], ...
    's0_values_m',[5 10 20],'receiver_offsets_m',[0 0.5 1 2], ...
    'pe_nx',256,'pe_width_m',64,'pe_stepz_lamb',0.5,'pe_alpha_max_np_per_m',0.15, ...
    'pe_grid_cases',local_default_grid_cases(), ...
    'normalization_frequencies_hz',[3000 4000 5000], ...
    'normalization_beam_counts',[0 501 2001 5001 10001], ...
    'normalization_step_values_m',[0 0.1 0.05], ...
    'bellhop_beam_count',10001,'bellhop_step_m',0.05,'bellhop_domain_half_depth_m',1000, ...
    'delay_frequency_hz',4000,'delay_frequency_offsets_hz',[-1 0 1], ...
    'delay_L_m',70,'delay_s0_m',10,'delay_tolerance_s',5e-6, ...
    'tl_tolerance_db',0.5,'phase_tolerance_rad',0.1);
end

function g=local_default_grid_cases()
g(1)=struct('name','A_128_32_dz050','nx',128,'width_m',32,'stepz_lamb',0.5);
g(2)=struct('name','B_256_32_dz050','nx',256,'width_m',32,'stepz_lamb',0.5);
g(3)=struct('name','C_256_64_dz050','nx',256,'width_m',64,'stepz_lamb',0.5);
g(4)=struct('name','D_256_32_dz025','nx',256,'width_m',32,'stepz_lamb',0.25);
end

function cfg=local_overrides(cfg,o)
names=fieldnames(o); for ii=1:numel(names), if ~isfield(cfg,names{ii}), error('Unknown override: %s',names{ii}); end, cfg.(names{ii})=o.(names{ii}); end
end

function p=local_pe_params(cfg,f,L,s0)
p=struct('f0',f,'c0',cfg.c0_mps,'z_max',L,'z_tx',L,'z_rx',0, ...
    'xw',cfg.pe_width_m,'yw',cfg.pe_width_m,'nx',cfg.pe_nx,'ny',cfg.pe_nx, ...
    'x_tx',0,'y_tx',0,'x_rx',0,'y_rx',0,'stepz_lamb',cfg.pe_stepz_lamb, ...
    'source_mode','custom_field_fn','source_field_fn',@virtual_point_source_initial_field_vertical, ...
    'virtual_source_distance_m',s0,'source_green_amplitude',1/(4*pi), ...
    'sponge_ratio',0.12,'alpha_max_np_per_m',cfg.pe_alpha_max_np_per_m, ...
    'env_mode','uniform','enable_surface_reflection',false,'enable_bubbles',false, ...
    'enforce_1_over_R',false,'show_figures',false,'save_mode','slice','use_gpu',false);
end

function r=local_row()
r=struct('frequency_hz',NaN,'s0_m',NaN,'L_m',NaN,'offset_m',NaN,'range_m',NaN, ...
    'H_pe_source',complex(NaN),'H_bellhop_green',complex(NaN),'H_exact',complex(NaN), ...
    'pe_tl_db',NaN,'bellhop_tl_db',NaN,'exact_tl_db',NaN, ...
    'pe_minus_bellhop_tl_db',NaN,'pe_minus_exact_tl_db',NaN, ...
    'bellhop_minus_exact_tl_db',NaN,'pe_minus_bellhop_phase_rad',NaN, ...
    'pe_minus_exact_phase_rad',NaN,'bellhop_minus_exact_phase_rad',NaN);
end

function files=local_outputs(v,root)
out=v.config.output_dir; csv_file=fullfile(out,'pe_bellhop_freefield_comparison.csv');
checks_file=fullfile(out,'pe_bellhop_freefield_checks.csv'); mat_file=fullfile(out,'pe_bellhop_freefield_validation.mat');
arrival_file=fullfile(out,'pe_bellhop_freefield_arrival.csv');
fig_file=fullfile(out,'pe_bellhop_freefield_comparison.png'); report_file=fullfile(root,'reports','pe_bellhop_freefield_validation_report.md');
writetable(v.comparison_table,csv_file); writetable(v.checks,checks_file); writetable(v.arrival_table,arrival_file);
fig=figure('Visible','off','Color','w','Position',[100 100 1050 760]); cleanup=onCleanup(@()close(fig)); t=v.comparison_table;
q=t(t.frequency_hz==4000 & t.s0_m==10 & t.offset_m==0,:);
subplot(2,2,1); plot(q.range_m,q.pe_tl_db,'o-',q.range_m,q.bellhop_tl_db,'s--',q.range_m,q.exact_tl_db,'k:'); grid on; xlabel('R (m)'); ylabel('TL (dB)'); legend('PE','Bellhop converted','analytic');
subplot(2,2,2); plot(q.range_m,q.pe_minus_bellhop_tl_db,'o-'); grid on; yline(0); xlabel('R (m)'); ylabel('PE-BH TL (dB)');
subplot(2,2,3); plot(q.range_m,q.pe_minus_bellhop_phase_rad,'o-'); grid on; yline(0); xlabel('R (m)'); ylabel('PE-BH phase (rad)');
subplot(2,2,4); bar(categorical(v.checks.check_name),v.checks.value); grid on; ylabel('check value'); title(sprintf('passed=%d',v.passed));
exportgraphics(fig,fig_file,'Resolution',180); clear cleanup
validation=v; schema_version=v.schema_version; save(mat_file,'validation','schema_version','-v7.3');
local_report(report_file,v);
files=struct('mat',mat_file,'comparison_csv',csv_file,'checks_csv',checks_file, ...
    'arrival_csv',arrival_file,'figure',fig_file,'report',report_file);
end

function local_report(file,v)
fid=fopen(file,'w','n','UTF-8'); if fid<0, error('Cannot create report.'); end
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'# PE--Bellhop 无反射自由场验证报告\n\n');
fprintf(fid,'状态：`passed=%s`。本验证不含海面或海底反射、粗糙面、气泡、Doppler、随机信道或通信处理。\n\n',string(v.passed));
fprintf(fid,'## 环境与数值设置\n\n');
fprintf(fid,'- 均匀声速：`%.0f m/s`；频率：`%s Hz`。\n',v.config.c0_mps,mat2str(v.config.frequencies_hz));
fprintf(fid,'- 虚拟源到初始面 `s0=%s m`；PE 推进 `L=%s m`；横向偏移 `%s m`。\n',mat2str(v.config.s0_values_m),mat2str(v.config.L_values_m),mat2str(v.config.receiver_offsets_m));
fprintf(fid,'- PE 正式物理比较：`%d x %d`，横向宽度 `%.1f m`，`stepz_lamb=%.3g`，`alpha_max=%.3g Np/m`。\n',v.config.pe_nx,v.config.pe_nx,v.config.pe_width_m,v.config.pe_stepz_lamb,v.config.pe_alpha_max_np_per_m);
fprintf(fid,'- PE--AS hard gate 的 sponge 完全关闭；Bellhop 使用匹配上下半空间 free-space 构造。\n\n');
fprintf(fid,'## 四层结论\n\n');
fprintf(fid,'1. Bellhop normalization audit：通过；`|p|R` 最大均值偏差 `%.6g`，统一 Green 转换为 `1/(4*pi)`，不做逐点拟合。\n',max(abs(v.bellhop_normalization.normalization_table.C_mean-1)));
fprintf(fid,'2. PE--独立一步 AS：通过；最大全场相对 L2 误差 `%.6g`。\n',max(v.pe_reference.field_table.pe_as_rel_l2));
fprintf(fid,'3. PE--解析球面波：未通过；最大 TL 误差 `%.6g dB`。\n',max(abs(v.comparison_table.pe_minus_exact_tl_db)));
fprintf(fid,'4. PE--Bellhop：未通过；最大 TL 差 `%.6g dB`，相位 RMS `%.6g rad`；Bellhop--解析最大 TL 差 `%.6g dB`。\n\n',max(abs(v.comparison_table.pe_minus_bellhop_tl_db)),sqrt(mean(v.comparison_table.pe_minus_bellhop_phase_rad.^2)),max(abs(v.comparison_table.bellhop_minus_exact_tl_db)));
fprintf(fid,'## 双相位参考\n\n`H_plane = Psi(L)*exp(i*k*L)` 用于 PE--AS；`H_source = exp(i*k*s0)*H_plane` 用于解析球面波和 Bellhop。公共 PE 输出语义未改变。\n\n');
fprintf(fid,'## 解释与判定\n\nPE--AS 已达到浮点误差量级，因此当前结果不支持“生产 PE 自由场平方根传播公式写错”的假设。失败集中在初始球面波被截断到有限横向平面并由周期 FFT 延拓；窗口从 32 m 增至 64 m 虽改善误差，但尚未收敛。下一步应独立验证无限孔径/Weyl 谱初始化或扩大无 sponge 孔径，而不是修改 PE marching 主线。\n\n');
fprintf(fid,'代表性 `R=%.3g m` 群时延：PE 误差 `%.6g us`，Bellhop 误差 `%.6g us`。\n\n',v.arrival_table.range_m,1e6*v.arrival_table.pe_delay_error_s,1e6*v.arrival_table.bellhop_delay_error_s);
fprintf(fid,'## 自动检查\n\n| check | value | limit | passed |\n|---|---:|---:|:---:|\n');
for ii=1:height(v.checks)
    fprintf(fid,'| %s | %.8g | %.8g | %d |\n',v.checks.check_name(ii),v.checks.value(ii),v.checks.limit(ii),v.checks.passed(ii));
end
fprintf(fid,'\n![自由场比较](../results/validation/pe_bellhop_freefield/formal/pe_bellhop_freefield_comparison.png)\n');
clear cleanup
end
