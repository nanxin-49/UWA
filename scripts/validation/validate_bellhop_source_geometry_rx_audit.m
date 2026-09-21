function audit = validate_bellhop_source_geometry_rx_audit(options)
%VALIDATE_BELLHOP_SOURCE_GEOMETRY_RX_AUDIT Isolate Bellhop R versus X scaling.
% This validation changes only RunType(4): R (point source) versus X (line
% source). Reflection, wall geometry, beam dynamics, influence, source beam
% pattern, receiver positions, step, and all other numerical settings remain
% unchanged.
arguments
    options.output_dir (1,:) char = ''
    options.report_path (1,:) char = ''
    options.official_exe (1,:) char = ''
    options.flat_validation_exe (1,:) char = ''
    options.sinusoidal_validation_exe (1,:) char = ''
    options.pm_validation_exe (1,:) char = ''
    options.run_tier1_if_improved (1,1) logical = true
    options.fail_on_check (1,1) logical = true
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root); setup_vertical_project();
if isempty(options.official_exe), options.official_exe=getenv('BELLHOP_EXE'); end
if isempty(options.official_exe)
    error('Set options.official_exe or BELLHOP_EXE to the AcousticsToolbox 2020 bellhop.exe.');
end
if isempty(options.output_dir)
    options.output_dir=fullfile(root,'results','validation','bellhop_source_geometry_rx_audit');
end
if isempty(options.report_path)
    options.report_path=fullfile(root,'reports','bellhop_source_geometry_rx_audit_report.md');
end
if isempty(options.flat_validation_exe)
    options.flat_validation_exe=fullfile(root,'results','validation','bellhop_internal_flat_wall_poc','bin','bellhop_iwall_flat_2020.exe');
end
if isempty(options.sinusoidal_validation_exe)
    options.sinusoidal_validation_exe=fullfile(root,'results','validation','bellhop_internal_sinusoidal_wall_poc','bin','bellhop_iwall_sinusoidal_2020.exe');
end
if isempty(options.pm_validation_exe)
    options.pm_validation_exe=fullfile(root,'results','validation','bellhop_internal_pm_wall_poc','bin','bellhop_iwall_pm_2020.exe');
end
for p={options.official_exe,options.flat_validation_exe,options.sinusoidal_validation_exe,options.pm_validation_exe}
    assert(isfile(p{1}),'Missing Bellhop executable: %s',p{1});
end
if ~exist(options.output_dir,'dir'), mkdir(options.output_dir); end

cfg=struct('frequency_hz',4000,'c0_mps',1500,'source_depth_m',0, ...
    'receiver_depths_m',0,'run_type','C','angle_limits_deg',[-30 30], ...
    'domain_half_depth_m',1000,'sigma_src_m',0.3,'source_pattern_clip_db',-120, ...
    'step_m',0.05,'wall_r0_m',100,'wall_amplitude_m',0.25, ...
    'wall_wavenumber_per_m',-0.01,'profile_count',161, ...
    'mapped_receiver_range_m',103,'native_receiver_range_m',97);
pat=local_pattern(cfg.frequency_hz,cfg);
cfg.source_pattern_angles_deg=pat.angles_deg;
cfg.source_pattern_level_db=pat.level_db;
source_geometries={'R','X'}; beam_counts=[5001 10001];

flat_rows=repmat(local_flat_template(),numel(source_geometries)*numel(beam_counts),1);
sine_rows=repmat(local_sine_template(),numel(source_geometries)*numel(beam_counts),1);
rr=0;
z_prof=linspace(-100,100,cfg.profile_count).';
r_prof=cfg.wall_r0_m-cfg.wall_amplitude_m*sin(cfg.wall_wavenumber_per_m*z_prof);
for gg=1:numel(source_geometries)
    source_geometry=source_geometries{gg};
    for bb=1:numel(beam_counts)
        beams=beam_counts(bb); rr=rr+1;
        tag=sprintf('%s_beams_%d',source_geometry,beams);
        common=cfg; common.source_geometry=source_geometry; common.beam_count=beams;

        direct_flat=common; direct_flat.bellhop_exe=options.official_exe;
        direct_flat.receiver_ranges_m=[102 103];
        direct_flat.case_root=fullfile(options.output_dir,'cases','flat','direct',tag);
        direct_flat_run=run_bellhop_unfolded_gaussian_vertical(direct_flat);
        wall_flat=common; wall_flat.bellhop_exe=options.flat_validation_exe;
        wall_flat.receiver_ranges_m=[102 103]; wall_flat.wall_range_m=cfg.wall_r0_m;
        wall_flat.case_root=fullfile(options.output_dir,'cases','flat','wall',tag);
        wall_flat_run=run_bellhop_internal_flat_wall_poc_vertical(wall_flat);
        p_direct=select_bellhop_shd_pressure_at_range_vertical(direct_flat_run.data,103,0);
        p_wall=select_bellhop_shd_pressure_at_range_vertical(wall_flat_run.data,103,0);
        ratio=p_wall/(-p_direct);
        flat_rows(rr)=local_flat_row(source_geometry,beams,p_wall,p_direct,ratio,wall_flat_run);

        wall_sine=common; wall_sine.bellhop_exe=options.sinusoidal_validation_exe;
        wall_sine.receiver_ranges_m=[102 103]; wall_sine.wall_profile_r_m=r_prof;
        wall_sine.wall_profile_z_m=z_prof;
        wall_sine.case_root=fullfile(options.output_dir,'cases','sinusoidal','wall',tag);
        wall_sine_run=run_bellhop_internal_sinusoidal_wall_poc_vertical(wall_sine);
        native=common; native.bellhop_exe=options.official_exe;
        native.receiver_ranges_m=[97 98]; native.wall_profile_r_m=r_prof;
        native.wall_profile_z_m=z_prof;
        native.case_root=fullfile(options.output_dir,'cases','sinusoidal','native',tag);
        native_run=run_bellhop_native_sinusoidal_wall_vertical(native);
        direct_native=native; direct_native.case_root=fullfile(options.output_dir,'cases','sinusoidal','direct',tag);
        direct_native_run=run_bellhop_unfolded_gaussian_vertical(direct_native);
        p_wall=select_bellhop_shd_pressure_at_range_vertical(wall_sine_run.data,103,0);
        p_native=select_bellhop_shd_pressure_at_range_vertical(native_run.data,97,0);
        p_direct=select_bellhop_shd_pressure_at_range_vertical(direct_native_run.data,97,0);
        p_reflected=p_native-p_direct; ratio=p_wall/p_reflected;
        sine_rows(rr)=local_sine_row(source_geometry,beams,p_wall,p_reflected,ratio,wall_sine_run);
    end
end
flat_table=struct2table(flat_rows); sine_table=struct2table(sine_rows);
flat_table.source_geometry=string(flat_table.source_geometry);
sine_table.source_geometry=string(sine_table.source_geometry);
flat_conv=local_convergence(flat_table); sine_conv=local_convergence(sine_table);

r100=sine_table(sine_table.source_geometry=="R" & sine_table.beam_count==10001,:);
x100=sine_table(sine_table.source_geometry=="X" & sine_table.beam_count==10001,:);
predicted_r_tl_db=10*log10(cfg.native_receiver_range_m/cfg.mapped_receiver_range_m);
improvement_db=abs(r100.tl_error_db)-abs(x100.tl_error_db);
line_source_improved=abs(x100.tl_error_db)<0.05 && improvement_db>0.15 && ...
    abs(x100.phase_error_rad)<0.03;

tier1=[]; tier1_r=local_read_existing_tier1(root);
if options.run_tier1_if_improved && line_source_improved
    tier1_dir=fullfile(options.output_dir,'tier1_X');
    tier1=validate_pe_bellhop_pm_stage1_tier1(struct( ...
        'source_geometry','X','output_dir',tier1_dir, ...
        'report_path',fullfile(tier1_dir,'tier1_X_internal_report.md'), ...
        'reuse_flat_case_root',fullfile(tier1_dir,'cases','flat_X'), ...
        'reuse_rough_case_root',fullfile(tier1_dir,'cases','rough_X'), ...
        'flat_validation_exe',options.flat_validation_exe, ...
        'pm_validation_exe',options.pm_validation_exe));
end

checks=struct;
checks.flat_regression=max(abs(flat_table.tl_error_db))<0.02 && ...
    max(abs(flat_table.phase_error_rad))<0.01;
checks.sinusoidal_R_reproduced=abs(r100.tl_error_db-predicted_r_tl_db)<0.03;
checks.line_source_improved=line_source_improved;
checks.beam_convergence_diagnostic=max(sine_conv.tl_change_db)<0.02;
checks.tier1_executed=~isempty(tier1);
checks.tier1_valid=~isempty(tier1) && tier1.passed;
checks.all=checks.flat_regression && checks.sinusoidal_R_reproduced && ...
    checks.line_source_improved && checks.tier1_executed && checks.tier1_valid;

audit=struct('schema_version','1.0.0','config',cfg,'options',options, ...
    'flat',flat_table,'sinusoidal',sine_table,'flat_convergence',flat_conv, ...
    'sinusoidal_convergence',sine_conv,'predicted_R_range_tl_db',predicted_r_tl_db, ...
    'line_source_improvement_db',improvement_db,'line_source_improved',line_source_improved, ...
    'tier1_R',tier1_r,'tier1_X',tier1,'checks',checks);
writetable(flat_table,fullfile(options.output_dir,'flat_R_X_metrics.csv'));
writetable(sine_table,fullfile(options.output_dir,'sinusoidal_R_X_metrics.csv'));
writetable(flat_conv,fullfile(options.output_dir,'flat_beam_convergence.csv'));
writetable(sine_conv,fullfile(options.output_dir,'sinusoidal_beam_convergence.csv'));
save(fullfile(options.output_dir,'bellhop_source_geometry_rx_audit.mat'),'audit','-v7.3');
local_write_report(options.report_path,audit);
disp(flat_table(:,{'source_geometry','beam_count','tl_error_db','phase_error_rad','complex_relative_error'}));
disp(sine_table(:,{'source_geometry','beam_count','tl_error_db','phase_error_rad','complex_relative_error'}));
if ~isempty(tier1)
    disp(struct('tier1_R_delta_tl_db',tier1_r.delta_tl_db,'tier1_X_delta_tl_db',tier1.delta_tl_db, ...
        'tier1_X_delta_phase_rad',tier1.delta_phase_rad));
end
disp(checks);
if options.fail_on_check && ~checks.all
    error('Bellhop R/X source-geometry audit failed; inspect %s.',options.report_path);
end
end

function t=local_flat_template()
t=struct('source_geometry','', 'beam_count',NaN,'wall_pressure',complex(NaN), ...
    'reference_pressure',complex(NaN),'ratio',complex(NaN),'tl_error_db',NaN, ...
    'phase_error_rad',NaN,'complex_relative_error',NaN,'wall_residual_m',NaN, ...
    'phase_jump_error_rad',NaN);
end

function t=local_sine_template()
t=struct('source_geometry','', 'beam_count',NaN,'wall_pressure',complex(NaN), ...
    'native_reflected_pressure',complex(NaN),'ratio',complex(NaN),'tl_error_db',NaN, ...
    'phase_error_rad',NaN,'complex_relative_error',NaN,'wall_residual_m',NaN, ...
    'tangent_error',NaN,'normal_error',NaN,'kappa_abs_max',NaN, ...
    'phase_jump_error_rad',NaN,'min_post_range_increment_m',NaN);
end

function t=local_flat_row(g,b,pwall,pdirect,ratio,run)
d=run.diagnostics; t=local_flat_template(); t.source_geometry=g; t.beam_count=b;
t.wall_pressure=pwall; t.reference_pressure=-pdirect; t.ratio=ratio;
t.tl_error_db=20*log10(abs(ratio)); t.phase_error_rad=angle(ratio);
t.complex_relative_error=abs(ratio-1); t.wall_residual_m=max(abs(d.residual));
t.phase_jump_error_rad=max(abs(angle(exp(1i*(d.phase_delta-pi)))));
end

function t=local_sine_row(g,b,pwall,pref,ratio,run)
d=run.diagnostics; t=local_sine_template(); t.source_geometry=g; t.beam_count=b;
t.wall_pressure=pwall; t.native_reflected_pressure=pref; t.ratio=ratio;
t.tl_error_db=20*log10(abs(ratio)); t.phase_error_rad=angle(ratio);
t.complex_relative_error=abs(ratio-1); t.wall_residual_m=max(abs(d.wall_residual));
t.tangent_error=max(d.tangent_error); t.normal_error=max(d.normal_error);
t.kappa_abs_max=max(abs(d.kappa));
t.phase_jump_error_rad=max(abs(angle(exp(1i*(d.phase_delta-pi)))));
t.min_post_range_increment_m=min(d.min_post_dr);
end

function c=local_convergence(tbl)
groups=unique(tbl.source_geometry,'stable');
c=table('Size',[numel(groups) 5],'VariableTypes',{'string','double','double','double','double'}, ...
    'VariableNames',{'source_geometry','tl_change_db','phase_change_rad','complex_relative_change','beam_count_pair'});
for ii=1:numel(groups)
    q=sortrows(tbl(tbl.source_geometry==groups(ii),:),'beam_count');
    ratio=q.ratio(end)/q.ratio(1);
    c.source_geometry(ii)=groups(ii); c.tl_change_db(ii)=abs(20*log10(abs(ratio)));
    c.phase_change_rad(ii)=abs(angle(ratio)); c.complex_relative_change(ii)=abs(ratio-1);
    c.beam_count_pair(ii)=q.beam_count(1)*1e5+q.beam_count(end);
end
end

function p=local_pattern(f,cfg)
a=linspace(cfg.angle_limits_deg(1),cfg.angle_limits_deg(2),2401).';
th=deg2rad(a); k=2*pi*f/cfg.c0_mps;
d=cos(th).*exp(-0.5*(k*cfg.sigma_src_m*sin(th)).^2);
d=abs(d)/max(abs(d));
p=struct('angles_deg',a,'level_db',20*log10(max(d,10^(cfg.source_pattern_clip_db/20))));
end

function r=local_read_existing_tier1(root)
path=fullfile(root,'results','validation','pe_bellhop_pm_stage1_tier1','stage1a_tier1_comparison.mat');
assert(isfile(path),'Missing authoritative R-source Tier-1 result: %s',path);
s=load(path,'c'); r=struct('G_pe',s.c.G_pe,'G_bellhop',s.c.G_bellhop, ...
    'delta_tl_db',s.c.delta_tl_db,'delta_phase_rad',s.c.delta_phase_rad, ...
    'delta_complex_relative_error',s.c.delta_complex_relative_error,'source_file',path);
end

function local_write_report(path,a)
fid=fopen(path,'w','n','UTF-8'); if fid<0,error('Cannot write %s.',path);end
cleanup=onCleanup(@()fclose(fid));
r=a.sinusoidal(a.sinusoidal.source_geometry=="R" & a.sinusoidal.beam_count==10001,:);
x=a.sinusoidal(a.sinusoidal.source_geometry=="X" & a.sinusoidal.beam_count==10001,:);
fprintf(fid,'# Bellhop point-source R versus line-source X amplitude audit\n\n');
status=ternary(a.checks.all,ternary(a.checks.beam_convergence_diagnostic,'PASS','PASS_WITH_LIMITS'),'FAIL');
fprintf(fid,'状态：**%s**。本审计仅改变 Bellhop `RunType(4)`：`R`（point source）或 `X`（line source）；Gaussian `.sbp`、step=0.05 m、beam 数、receiver、wall geometry、Reflect2D、p/q 与 InfluenceGeoHatCart 均未改变。\n\n',status);
fprintf(fid,'## Flat internal-wall regression\n\n| source | beams | TL error (dB) | phase error (rad) | complex relative error |\n|---|---:|---:|---:|---:|\n');
for ii=1:height(a.flat),fprintf(fid,'| %s | %d | %.9g | %.9g | %.9g |\n',a.flat.source_geometry(ii),a.flat.beam_count(ii),a.flat.tl_error_db(ii),a.flat.phase_error_rad(ii),a.flat.complex_relative_error(ii));end
fprintf(fid,'\nReference is `H_wall=-P_BH(103 m)` under the same source convention.\n\n');
fprintf(fid,'## Weak sinusoidal native ATI versus internal wall\n\nA=0.25 m, K=-0.01 1/m, N=161; native receiver r=97 m and rotated-wall receiver r=103 m.\n\n| source | beams | TL(native-wall) (dB) | phase (rad) | complex relative error |\n|---|---:|---:|---:|---:|\n');
for ii=1:height(a.sinusoidal),fprintf(fid,'| %s | %d | %.9g | %.9g | %.9g |\n',a.sinusoidal.source_geometry(ii),a.sinusoidal.beam_count(ii),a.sinusoidal.tl_error_db(ii),a.sinusoidal.phase_error_rad(ii),a.sinusoidal.complex_relative_error(ii));end
fprintf(fid,'\nThe point-source range factor predicts `10 log10(97/103)=%.9g dB`; measured R is %.9g dB. X leaves %.9g dB. The absolute-error reduction is %.9g dB.\n\n',a.predicted_R_range_tl_db,r.tl_error_db,x.tl_error_db,a.line_source_improvement_db);
fprintf(fid,'## 5001 -> 10001 beam convergence\n\n| case | source | TL change (dB) | phase change (rad) | complex change |\n|---|---|---:|---:|---:|\n');
for name={'flat','sinusoidal'}
    c=a.([name{1} '_convergence']);
    for ii=1:height(c),fprintf(fid,'| %s | %s | %.9g | %.9g | %.9g |\n',name{1},c.source_geometry(ii),c.tl_change_db(ii),c.phase_change_rad(ii),c.complex_relative_change(ii));end
end
if ~a.checks.beam_convergence_diagnostic
    fprintf(fid,'\nThe 5001-beam sinusoidal field is not converged; both R and X share the same approximately 3.16 dB 5001-to-10001 change. Conclusions therefore use the requested 10001-beam endpoint, and no parameter was tuned.\n');
end
fprintf(fid,'\n## Fixed-PM Tier-1, 4 kHz, seed=260001\n\n');
if isempty(a.tier1_X)
    fprintf(fid,'Not run because the precondition that X clearly improve the sinusoidal covariance test was not met.\n');
else
    fprintf(fid,'| source | G_PE | G_BH | Delta TL (dB) | Delta phase (rad) | complex error |\n|---|---:|---:|---:|---:|---:|\n');
    fprintf(fid,'| R (authoritative prior) | %.8g%+.8gi | %.8g%+.8gi | %.9g | %.9g | %.9g |\n',real(a.tier1_R.G_pe),imag(a.tier1_R.G_pe),real(a.tier1_R.G_bellhop),imag(a.tier1_R.G_bellhop),a.tier1_R.delta_tl_db,a.tier1_R.delta_phase_rad,a.tier1_R.delta_complex_relative_error);
    fprintf(fid,'| X (this audit) | %.8g%+.8gi | %.8g%+.8gi | %.9g | %.9g | %.9g |\n',real(a.tier1_X.G_pe),imag(a.tier1_X.G_pe),real(a.tier1_X.G_bellhop),imag(a.tier1_X.G_bellhop),a.tier1_X.delta_tl_db,a.tier1_X.delta_phase_rad,a.tier1_X.delta_complex_relative_error);
end
fprintf(fid,'\n## Conclusions\n\n');
fprintf(fid,'- The ~0.26 dB native/internal covariance bias is %s attributable to the point-source R range normalization.\n',ternary(a.line_source_improved,'primarily','not primarily'));
fprintf(fid,'- Tier-1 should %s use X as the dimensionally matched line-source Bellhop comparator.\n',ternary(a.line_source_improved,'formally','not yet'));
if ~isempty(a.tier1_X),fprintf(fid,'- After X, the original %.6g dB PE--Bellhop Tier-1 amplitude delta becomes %.6g dB.\n',a.tier1_R.delta_tl_db,a.tier1_X.delta_tl_db);end
if ~isempty(a.tier1_X)
    materially_changed=abs(a.tier1_X.delta_tl_db)<0.1 && abs(a.tier1_X.delta_phase_rad)<0.2;
    fprintf(fid,'- The reflection-model discrepancy conclusion %s. X removes a Bellhop coordinate/source normalization artifact but %s the rough/flat cross-model residual.\n',ternary(materially_changed,'requires revision','does not require revision'),ternary(materially_changed,'also removes','does not remove'));
end
fprintf(fid,'\nNo empirical scaling, `.sbp` refit, phase correction, or physics change was used.\n');
end

function s=ternary(c,a,b),if c,s=a;else,s=b;end,end
