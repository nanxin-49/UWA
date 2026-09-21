function scan = validate_bellhop_transverse_parameter_scan_vertical(overrides)
%VALIDATE_BELLHOP_TRANSVERSE_PARAMETER_SCAN_VERTICAL
% Bellhop-only transverse phase scan against an independent AS reference.
% The saved formal PE result is read; no PE march is executed here.
if nargin<1 || isempty(overrides), overrides=struct(); end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
run(fullfile(root,'scripts','bootstrap_project.m'));
cfg=local_config(root,overrides);
if exist(cfg.formal_mat,'file')~=2, error('Formal validation MAT is missing.'); end
formal=load(cfg.formal_mat,'validation'); formal=formal.validation;
if ~isfield(formal,'spatial_profiles') || ~isfield(formal.spatial_profiles,'rows')
    error('Formal validation MAT has no spatial profiles.');
end
if exist(cfg.bellhop_exe,'file')~=2, error('Bellhop executable is missing.'); end
freq=cfg.frequency_hz; offsets=cfg.receiver_offsets_m(:).';
x=formal.pe.chain.x(:).'; y=formal.pe.chain.y(:).'; [X,Y]=meshgrid(x,y);
source_cfg=formal.pe.chain.config; source_cfg.x_tx=0; source_cfg.y_tx=0; source_cfg.sigma_src_m=cfg.sigma_src_m;
[psi0,source_meta]=gaussian_source_initial_field_vertical(X,Y,source_cfg);
[~,ix0]=min(abs(x)); [~,iy0]=min(abs(y));
[psi_d,as_meta_d]=exact_angular_spectrum_one_step_vertical(psi0,cfg.pe_width_m,cfg.pe_width_m,freq,cfg.c0_mps,cfg.direct_range_m);
[psi_r,as_meta_r]=exact_angular_spectrum_one_step_vertical(psi0,cfg.pe_width_m,cfg.pe_width_m,freq,cfg.c0_mps,cfg.reflect_range_m);
psi_r=-psi_r; as_d=psi_d(iy0,ix0); as_r=psi_r(iy0,ix0);
as_direct=zeros(size(offsets)); as_reflect=zeros(size(offsets));
for jj=1:numel(offsets)
    [~,ix]=min(abs(x-offsets(jj))); [~,iy]=min(abs(y));
    as_direct(jj)=psi_d(iy,ix)/as_d; as_reflect(jj)=psi_r(iy,ix)/as_r;
end
template=struct('step_m',NaN,'angle_half_deg',NaN,'sbp_samples',NaN,'max_direct_tl_db',NaN, ...
    'max_reflect_tl_db',NaN,'max_direct_phase_rad',NaN,'max_reflect_phase_rad',NaN, ...
    'max_direct_complex_error',NaN,'max_reflect_complex_error',NaN,'max_complex_error',NaN, ...
    'max_phase_error_rad',NaN,'passed',false,'case_root','');
rows=repmat(template,1,numel(cfg.step_values_m)*numel(cfg.angle_half_values_deg)*numel(cfg.sbp_sample_values)); nrow=0;
for is=1:numel(cfg.step_values_m)
    for ia=1:numel(cfg.angle_half_values_deg)
        for inn=1:numel(cfg.sbp_sample_values)
            half=cfg.angle_half_values_deg(ia); ns=cfg.sbp_sample_values(inn);
            tag=sprintf('step%03d_ang%02d_n%d',round(cfg.step_values_m(is)*1000),round(half),ns);
            case_root=fullfile(cfg.output_dir,'cases',tag); pat=local_pattern(freq,cfg,half,ns);
            c=struct('bellhop_exe',cfg.bellhop_exe,'case_root',case_root,'frequency_hz',freq, ...
                'c0_mps',cfg.c0_mps,'source_depth_m',0,'receiver_depths_m',offsets, ...
                'receiver_ranges_m',[cfg.direct_range_m cfg.reflect_range_m],'run_type','C', ...
                'beam_count',cfg.beam_count,'angle_limits_deg',[-half half], ...
                'step_m',cfg.step_values_m(is),'domain_half_depth_m',cfg.domain_half_depth_m, ...
                'source_pattern_angles_deg',pat.angles_deg,'source_pattern_level_db',pat.level_db);
            result=run_bellhop_unfolded_gaussian_vertical(c);
            % The installed Bellhop build uses the opposite transverse
            % phasor sign for the unfolded profile.  Match the formal
            % validation convention explicitly; this is not an amplitude fit.
            raw=result.data.pressure; bh_d=conj(raw(:,1).'/raw(1,1)); bh_r=conj(raw(:,2).'/raw(1,2));
            q=template; q.step_m=c.step_m; q.angle_half_deg=half; q.sbp_samples=ns; q.case_root=case_root;
            q.max_direct_tl_db=max(abs(20*log10(max(abs(bh_d),realmin)./max(abs(as_direct),realmin))));
            q.max_reflect_tl_db=max(abs(20*log10(max(abs(bh_r),realmin)./max(abs(as_reflect),realmin))));
            q.max_direct_phase_rad=max(abs(angle(bh_d.*conj(as_direct))));
            q.max_reflect_phase_rad=max(abs(angle(bh_r.*conj(as_reflect))));
            q.max_direct_complex_error=max(abs(bh_d-as_direct)./max(abs(as_direct),realmin));
            q.max_reflect_complex_error=max(abs(bh_r-as_reflect)./max(abs(as_reflect),realmin));
            q.max_complex_error=max(q.max_direct_complex_error,q.max_reflect_complex_error);
            q.max_phase_error_rad=max(q.max_direct_phase_rad,q.max_reflect_phase_rad);
            q.passed=q.max_direct_tl_db<=cfg.tl_limit_db && q.max_reflect_tl_db<=cfg.tl_limit_db && ...
                q.max_phase_error_rad<=cfg.phase_limit_rad && q.max_complex_error<=cfg.complex_limit;
            nrow=nrow+1; rows(nrow)=q;
        end
    end
end
rows=rows(1:nrow); table_out=struct2table(rows); [~,best]=min([rows.max_complex_error]);
summary=struct('best_index',best,'best_case',rows(best),'pass_count',sum([rows.passed]), ...
    'case_count',numel(rows),'max_complex_error',max([rows.max_complex_error]), ...
    'min_complex_error',min([rows.max_complex_error]),'source_meta',source_meta, ...
    'as_meta_direct',as_meta_d,'as_meta_reflect',as_meta_r,'as_direct',as_direct,'as_reflect',as_reflect);
checks=table(["all_cases_complete";"any_case_passes"], ...
    [numel(rows)==numel(cfg.step_values_m)*numel(cfg.angle_half_values_deg)*numel(cfg.sbp_sample_values);summary.pass_count>0], ...
    [true;true],[true;summary.pass_count>0],'VariableNames',{'check_name','value','limit','passed'});
scan=struct('schema_version','1.0.0','config',cfg,'rows',rows,'table',table_out, ...
    'summary',summary,'checks',checks,'conclusion',local_conclusion(summary,cfg));
scan.files=local_outputs(scan);
end

function cfg=local_config(root,o)
exe=getenv('BELLHOP_EXE'); if isempty(exe), exe='AcousticsToolbox_2020/windows-bin-20201102/bellhop.exe'; end
cfg=struct('bellhop_exe',exe,'formal_mat',fullfile(root,'results','validation','pe_bellhop_unfolded_flat_gaussian', ...
    'pe_bellhop_unfolded_flat_gaussian_validation.mat'),'output_dir',fullfile(root,'results','validation','pe_as_bellhop_transverse','bellhop_scan'), ...
    'report_path',fullfile(root,'reports','bellhop_transverse_parameter_scan_report.md'),'c0_mps',1500,'sigma_src_m',0.3, ...
    'pe_width_m',192.1875,'frequency_hz',8000,'direct_range_m',97,'reflect_range_m',103, ...
    'receiver_offsets_m',[0 1.953125 4.8828125 9.765625 14.6484375 19.53125], ...
    'step_values_m',[0.1 0.05 0.025],'angle_half_values_deg',[20 30 45],'sbp_sample_values',[1201 2401 4801], ...
    'beam_count',10001,'domain_half_depth_m',1000,'source_pattern_clip_db',-120, ...
    'tl_limit_db',0.25,'phase_limit_rad',0.05,'complex_limit',0.02);
names=fieldnames(o);
for ii=1:numel(names)
    if ~isfield(cfg,names{ii}), error('Unknown override: %s',names{ii}); end
    cfg.(names{ii})=o.(names{ii});
end
end

function pat=local_pattern(f,cfg,half_deg,n)
angles=linspace(-half_deg,half_deg,n).'; k=2*pi*f/cfg.c0_mps; th=deg2rad(angles);
d=cos(th).*exp(-0.5*(k*cfg.sigma_src_m*sin(th)).^2); d=abs(d)/max(abs(d));
pat=struct('angles_deg',angles,'level_db',20*log10(max(d,10^(cfg.source_pattern_clip_db/20))));
end

function text=local_conclusion(s,cfg)
if s.pass_count>0
    text=sprintf('At least one Bellhop-only configuration satisfies %.3g dB, %.3g rad, and %.3g complex-error limits; inspect the best case before changing formal settings.',cfg.tl_limit_db,cfg.phase_limit_rad,cfg.complex_limit);
else
    text='No scanned Bellhop-only configuration satisfies all transverse limits. The remaining discrepancy is not explained by the tested step, angle, or SBP sampling ranges; continue with the 2-D/3-D source-mapping audit.';
end
end

function files=local_outputs(scan)
out=scan.config.output_dir; if ~exist(out,'dir'), mkdir(out); end
writetable(scan.table,fullfile(out,'bellhop_transverse_parameter_scan.csv'));
writetable(scan.checks,fullfile(out,'bellhop_transverse_parameter_scan_checks.csv'));
mat_file=fullfile(out,'bellhop_transverse_parameter_scan.mat'); save(mat_file,'scan','-v7.3');
fig_file=fullfile(out,'bellhop_transverse_parameter_scan.png'); local_plot(scan,fig_file);
local_report(scan);
files=struct('mat',mat_file,'table',fullfile(out,'bellhop_transverse_parameter_scan.csv'), ...
    'checks',fullfile(out,'bellhop_transverse_parameter_scan_checks.csv'),'figure',fig_file,'report',scan.config.report_path);
end

function local_plot(scan,file)
r=scan.rows; fig=figure('Visible','off','Color','w','Position',[100 100 1200 700]); cleanup=onCleanup(@()close(fig));
subplot(1,2,1); scatter([r.angle_half_deg],[r.max_complex_error],50,[r.step_m],'filled'); colorbar; grid on;
xlabel('half angle (deg)'); ylabel('max complex error'); title('Bellhop transverse scan'); hold on; yline(scan.config.complex_limit,'r--');
subplot(1,2,2); scatter([r.angle_half_deg],[r.max_phase_error_rad],50,[r.step_m],'filled'); colorbar; grid on;
xlabel('half angle (deg)'); ylabel('max phase error (rad)'); title('Phase convergence'); hold on; yline(scan.config.phase_limit_rad,'r--');
exportgraphics(fig,file,'Resolution',180); clear cleanup
end

function local_report(scan)
fid=fopen(scan.config.report_path,'w','n','UTF-8'); if fid<0, error('Cannot create scan report.'); end
cleanup=onCleanup(@()fclose(fid)); s=scan.summary; c=scan.config; b=s.best_case;
fprintf(fid,'# Bellhop 横向参数扫描\n\n');
fprintf(fid,'仅运行 Bellhop，固定 %.0f Hz、10001 beams，并与独立角谱比较；没有重新运行 PE。扫描 %d cases：step=%s m，half-angle=%s deg，SBP samples=%s。Bellhop 横向输出按正式验证的相位约定取共轭，此为固定 convention 转换而非逐点拟合。\n\n',c.frequency_hz,s.case_count,mat2str(c.step_values_m),mat2str(c.angle_half_values_deg),mat2str(c.sbp_sample_values));
fprintf(fid,'## 最优 case\n\n- step=%.6g m，half-angle=%.6g deg，SBP samples=%d。\n- max TL=%.8g dB，max phase=%.8g rad，max complex error=%.8g。\n- 通过 case 数：%d/%d。\n\n',b.step_m,b.angle_half_deg,b.sbp_samples,max(b.max_direct_tl_db,b.max_reflect_tl_db),b.max_phase_error_rad,b.max_complex_error,s.pass_count,s.case_count);
fprintf(fid,'%s\n\n',scan.conclusion);
fprintf(fid,'详细数据：results/validation/pe_as_bellhop_transverse/bellhop_scan/bellhop_transverse_parameter_scan.csv\n');
clear cleanup
end
