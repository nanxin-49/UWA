function audit = validate_pe_as_bellhop_transverse_vertical(overrides)
%VALIDATE_PE_AS_BELLHOP_TRANSVERSE_VERTICAL
% Compare saved formal PE/Bellhop transverse samples with an independent AS.
% This validation-only audit does not run a PE march or Bellhop.
if nargin<1 || isempty(overrides), overrides=struct(); end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
run(fullfile(root,'scripts','bootstrap_project.m'));
cfg=local_config(root,overrides);
if exist(cfg.formal_mat,'file')~=2, error('Formal validation MAT is missing.'); end
loaded=load(cfg.formal_mat,'validation');
formal=loaded.validation;
if ~isfield(formal,'spatial_profiles') || ~isfield(formal.spatial_profiles,'rows')
    error('Formal MAT has no saved transverse samples.');
end
rows=formal.spatial_profiles.rows;
freq=cfg.frequencies_hz(:).'; offsets=cfg.receiver_offsets_m(:).';
x=formal.pe.chain.x(:).'; y=formal.pe.chain.y(:).';
[X,Y]=meshgrid(x,y);
source_cfg=formal.pe.chain.config;
source_cfg.x_tx=0; source_cfg.y_tx=0; source_cfg.sigma_src_m=cfg.sigma_src_m;
[psi0,source_meta]=gaussian_source_initial_field_vertical(X,Y,source_cfg);
[~,ix0]=min(abs(x)); [~,iy0]=min(abs(y));
template=struct('frequency_hz',NaN,'offset_m',NaN,'grid_offset_error_m',NaN, ...
    'pe_as_direct_tl_db',NaN,'pe_as_direct_phase_rad',NaN,'pe_as_direct_complex_error',NaN, ...
    'pe_as_reflect_tl_db',NaN,'pe_as_reflect_phase_rad',NaN,'pe_as_reflect_complex_error',NaN, ...
    'bh_as_direct_tl_db',NaN,'bh_as_direct_phase_rad',NaN,'bh_as_direct_complex_error',NaN, ...
    'bh_as_reflect_tl_db',NaN,'bh_as_reflect_phase_rad',NaN,'bh_as_reflect_complex_error',NaN, ...
    'pe_bh_direct_tl_db',NaN,'pe_bh_direct_phase_rad',NaN,'pe_bh_direct_complex_error',NaN, ...
    'pe_bh_reflect_tl_db',NaN,'pe_bh_reflect_phase_rad',NaN,'pe_bh_reflect_complex_error',NaN);
out_rows=repmat(template,1,numel(freq)*numel(offsets)); nrow=0;
for ff=1:numel(freq)
    [psi_d,as_meta_d]=exact_angular_spectrum_one_step_vertical(psi0,cfg.pe_width_m,cfg.pe_width_m, ...
        freq(ff),cfg.c0_mps,cfg.direct_range_m);
    [psi_r,as_meta_r]=exact_angular_spectrum_one_step_vertical(psi0,cfg.pe_width_m,cfg.pe_width_m, ...
        freq(ff),cfg.c0_mps,cfg.reflect_range_m);
    psi_r=-psi_r;
    idx=find(abs([rows.frequency_hz]-freq(ff))<1e-9);
    if numel(idx)~=numel(offsets), error('Saved row count mismatch.'); end
    [~,ord]=sort([rows(idx).offset_m]); idx=idx(ord);
    pe_d_axis=rows(idx(1)).pe_direct; pe_r_axis=rows(idx(1)).pe_reflect;
    bh_d_axis=rows(idx(1)).bh_direct; bh_r_axis=rows(idx(1)).bh_reflect;
    as_d_axis=psi_d(iy0,ix0); as_r_axis=psi_r(iy0,ix0);
    for jj=1:numel(offsets)
        rr=rows(idx(jj)); [~,ix]=min(abs(x-offsets(jj))); [~,iy]=min(abs(y));
        q=template; q.frequency_hz=freq(ff); q.offset_m=offsets(jj);
        q.grid_offset_error_m=x(ix)-offsets(jj);
        pe_d=rr.pe_direct/pe_d_axis; pe_r=rr.pe_reflect/pe_r_axis;
        bh_d=rr.bh_direct/bh_d_axis; bh_r=rr.bh_reflect/bh_r_axis;
        as_d=psi_d(iy,ix)/as_d_axis; as_r=psi_r(iy,ix)/as_r_axis;
        q=local_pair_metrics(q,'pe_as_direct',pe_d,as_d);
        q=local_pair_metrics(q,'pe_as_reflect',pe_r,as_r);
        q=local_pair_metrics(q,'bh_as_direct',bh_d,as_d);
        q=local_pair_metrics(q,'bh_as_reflect',bh_r,as_r);
        q=local_pair_metrics(q,'pe_bh_direct',pe_d,bh_d);
        q=local_pair_metrics(q,'pe_bh_reflect',pe_r,bh_r);
        nrow=nrow+1; out_rows(nrow)=q;
    end
end
out_rows=out_rows(1:nrow);
summary=local_summary(out_rows,cfg); checks=local_checks(summary,cfg);
audit=struct('schema_version','1.0.0','config',cfg,'source_meta',source_meta, ...
    'as_meta_direct',as_meta_d,'as_meta_reflect',as_meta_r,'rows',out_rows, ...
    'table',struct2table(out_rows),'summary',summary,'checks',checks, ...
    'conclusion',local_conclusion(summary,checks),'formal_status',formal.passed);
audit.files=local_outputs(audit);
if cfg.fail_on_check && ~all(checks.passed)
    error('PE-AS transverse audit failed; see %s.',audit.files.report);
end
end

function cfg=local_config(root,o)
cfg=struct('formal_mat',fullfile(root,'results','validation','pe_bellhop_unfolded_flat_gaussian', ...
    'pe_bellhop_unfolded_flat_gaussian_validation.mat'), ...
    'output_dir',fullfile(root,'results','validation','pe_as_bellhop_transverse'), ...
    'report_path',fullfile(root,'reports','pe_as_bellhop_transverse_audit_report.md'), ...
    'c0_mps',1500,'sigma_src_m',0.3,'pe_width_m',192.1875, ...
    'frequencies_hz',[4000 6000 8000], ...
    'receiver_offsets_m',[0 1.953125 4.8828125 9.765625 14.6484375 19.53125], ...
    'direct_range_m',97,'reflect_range_m',103, ...
    'pe_as_tl_limit_db',1e-10,'pe_as_phase_limit_rad',1e-10,'pe_as_complex_limit',1e-10, ...
    'bellhop_tl_limit_db',0.25,'bellhop_phase_limit_rad',0.05, ...
    'bellhop_complex_limit',0.02,'fail_on_check',false);
names=fieldnames(o);
for ii=1:numel(names)
    if ~isfield(cfg,names{ii}), error('Unknown override: %s',names{ii}); end
    cfg.(names{ii})=o.(names{ii});
end
end

function q=local_pair_metrics(q,prefix,a,b)
q.([prefix '_tl_db'])=20*log10(max(abs(a),realmin)/max(abs(b),realmin));
q.([prefix '_phase_rad'])=angle(a*conj(b));
q.([prefix '_complex_error'])=abs(a-b)/max(abs(b),realmin);
end

function s=local_summary(rows,cfg)
fields={'pe_as_direct','pe_as_reflect','bh_as_direct','bh_as_reflect','pe_bh_direct','pe_bh_reflect'};
for ii=1:numel(fields)
    p=fields{ii};
    s.(p)=struct('max_tl_db',max(abs([rows.([p '_tl_db'])])), ...
        'max_phase_rad',max(abs([rows.([p '_phase_rad'])])), ...
        'max_complex_error',max([rows.([p '_complex_error'])]));
end
s.max_grid_offset_error_m=max(abs([rows.grid_offset_error_m]));
s.pe_as_pass=max(s.pe_as_direct.max_complex_error,s.pe_as_reflect.max_complex_error)<=cfg.pe_as_complex_limit;
s.bellhop_as_pass=max(s.bh_as_direct.max_tl_db,s.bh_as_reflect.max_tl_db)<=cfg.bellhop_tl_limit_db && ...
    max(s.bh_as_direct.max_phase_rad,s.bh_as_reflect.max_phase_rad)<=cfg.bellhop_phase_limit_rad && ...
    max(s.bh_as_direct.max_complex_error,s.bh_as_reflect.max_complex_error)<=cfg.bellhop_complex_limit;
s.pe_bh_pass=max(s.pe_bh_direct.max_tl_db,s.pe_bh_reflect.max_tl_db)<=cfg.bellhop_tl_limit_db && ...
    max(s.pe_bh_direct.max_phase_rad,s.pe_bh_reflect.max_phase_rad)<=cfg.bellhop_phase_limit_rad && ...
    max(s.pe_bh_direct.max_complex_error,s.pe_bh_reflect.max_complex_error)<=cfg.bellhop_complex_limit;
end

function checks=local_checks(s,cfg)
names=["pe_as_direct";"pe_as_reflect";"bellhop_as";"pe_bellhop";"grid_alignment"];
values=[s.pe_as_direct.max_complex_error;s.pe_as_reflect.max_complex_error; ...
    max(s.bh_as_direct.max_complex_error,s.bh_as_reflect.max_complex_error); ...
    max(s.pe_bh_direct.max_complex_error,s.pe_bh_reflect.max_complex_error); ...
    s.max_grid_offset_error_m];
limits=[cfg.pe_as_complex_limit;cfg.pe_as_complex_limit;cfg.bellhop_complex_limit; ...
    cfg.bellhop_complex_limit;1e-12];
passed=[s.pe_as_pass;s.pe_as_pass;s.bellhop_as_pass;s.pe_bh_pass; ...
    s.max_grid_offset_error_m<=limits(5)];
checks=table(names,values,limits,passed,'VariableNames',{'check_name','value','limit','passed'});
end

function text=local_conclusion(s,checks)
if s.pe_as_pass && ~s.bellhop_as_pass
    text='PE agrees with the independent angular-spectrum transverse reference; the remaining PE--Bellhop transverse discrepancy is attributed to Bellhop/source mapping or the 2-D ray-beam representation, not the PE marching operator.';
elseif ~s.pe_as_pass
    text='PE does not agree with the independent transverse angular-spectrum reference; inspect field extraction, coordinates, and phase convention before scanning Bellhop.';
else
    text='PE, independent angular spectrum, and Bellhop satisfy the configured transverse checks.';
end
if ~all(checks.passed), text=[text ' The audit is diagnostic and does not relax formal thresholds.']; end
end

function files=local_outputs(audit)
out=audit.config.output_dir; if ~exist(out,'dir'), mkdir(out); end
writetable(audit.table,fullfile(out,'pe_as_bellhop_transverse.csv'));
writetable(audit.checks,fullfile(out,'pe_as_bellhop_transverse_checks.csv'));
mat_file=fullfile(out,'pe_as_bellhop_transverse_audit.mat'); save(mat_file,'audit','-v7.3');
fig_file=fullfile(out,'pe_as_bellhop_transverse.png'); local_plot(audit,fig_file);
local_report(audit);
files=struct('mat',mat_file,'table',fullfile(out,'pe_as_bellhop_transverse.csv'), ...
    'checks',fullfile(out,'pe_as_bellhop_transverse_checks.csv'),'figure',fig_file, ...
    'report',audit.config.report_path);
end

function local_plot(audit,file)
r=audit.rows; f=audit.config.frequencies_hz;
fig=figure('Visible','off','Color','w','Position',[100 100 1200 800]); cleanup=onCleanup(@()close(fig));
subplot(2,2,1); hold on;
for ii=1:numel(f)
    q=r([r.frequency_hz]==f(ii)); [~,ord]=sort([q.offset_m]); q=q(ord);
    plot([q.offset_m],[q.pe_as_direct_complex_error],'-o','DisplayName',sprintf('PE-AS %g kHz',f(ii)/1000));
end
grid on; xlabel('offset (m)'); ylabel('direct complex error'); legend('Location','northwest'); title('PE versus AS');
subplot(2,2,2); hold on;
for ii=1:numel(f)
    q=r([r.frequency_hz]==f(ii)); [~,ord]=sort([q.offset_m]); q=q(ord);
    plot([q.offset_m],[q.bh_as_direct_phase_rad],'-s','DisplayName',sprintf('BH-AS %g kHz',f(ii)/1000));
end
grid on; xlabel('offset (m)'); ylabel('direct phase error (rad)'); legend('Location','northwest'); title('Bellhop versus AS');
subplot(2,2,3); hold on;
for ii=1:numel(f)
    q=r([r.frequency_hz]==f(ii)); [~,ord]=sort([q.offset_m]); q=q(ord);
    plot([q.offset_m],[q.pe_bh_direct_complex_error],'-^','DisplayName',sprintf('PE-BH %g kHz',f(ii)/1000));
end
grid on; xlabel('offset (m)'); ylabel('direct complex error'); legend('Location','northwest'); title('Original PE versus Bellhop');
subplot(2,2,4); bar(categorical(audit.checks.check_name),audit.checks.value); grid on; ylabel('maximum error');
title(sprintf('PE-AS pass=%d',audit.summary.pe_as_pass));
exportgraphics(fig,file,'Resolution',180); clear cleanup
end

function local_report(audit)
fid=fopen(audit.config.report_path,'w','n','UTF-8'); if fid<0, error('Cannot create audit report.'); end
cleanup=onCleanup(@()fclose(fid)); s=audit.summary; c=audit.config;
fprintf(fid,'# PE--AS--Bellhop 横向剖面审计\n\n');
fprintf(fid,'本审计读取正式 MAT，不重新运行 PE 或 Bellhop。配置为 c=%.0f m/s、Gaussian sigma=%.3g m、PE 窗口 %.6g m，频率为 4/6/8 kHz，距离为 %.3g/%.3g m。\n\n',c.c0_mps,c.sigma_src_m,c.pe_width_m,c.direct_range_m,c.reflect_range_m);
fprintf(fid,'## 最大误差\n\n| 比较 | TL (dB) | phase (rad) | complex |\n|---|---:|---:|---:|\n');
fprintf(fid,'| PE-AS direct | %.8g | %.8g | %.8g |\n',s.pe_as_direct.max_tl_db,s.pe_as_direct.max_phase_rad,s.pe_as_direct.max_complex_error);
fprintf(fid,'| PE-AS reflect | %.8g | %.8g | %.8g |\n',s.pe_as_reflect.max_tl_db,s.pe_as_reflect.max_phase_rad,s.pe_as_reflect.max_complex_error);
fprintf(fid,'| Bellhop-AS direct | %.8g | %.8g | %.8g |\n',s.bh_as_direct.max_tl_db,s.bh_as_direct.max_phase_rad,s.bh_as_direct.max_complex_error);
fprintf(fid,'| Bellhop-AS reflect | %.8g | %.8g | %.8g |\n',s.bh_as_reflect.max_tl_db,s.bh_as_reflect.max_phase_rad,s.bh_as_reflect.max_complex_error);
fprintf(fid,'| PE-Bellhop direct | %.8g | %.8g | %.8g |\n',s.pe_bh_direct.max_tl_db,s.pe_bh_direct.max_phase_rad,s.pe_bh_direct.max_complex_error);
fprintf(fid,'| PE-Bellhop reflect | %.8g | %.8g | %.8g |\n\n',s.pe_bh_reflect.max_tl_db,s.pe_bh_reflect.max_phase_rad,s.pe_bh_reflect.max_complex_error);
fprintf(fid,'PE-AS pass: %s; Bellhop-AS pass: %s; PE-Bellhop pass: %s.\n\n',string(s.pe_as_pass),string(s.bellhop_as_pass),string(s.pe_bh_pass));
fprintf(fid,'%s\n\n',audit.conclusion);
fprintf(fid,'下一步：PE-AS 通过时继续 Bellhop 步长、角扇区和 SBP 映射审计；PE-AS 失败时先检查场提取、坐标和相位约定。\n\n');
fprintf(fid,'![audit](../results/validation/pe_as_bellhop_transverse/pe_as_bellhop_transverse.png)\n');
clear cleanup
end
