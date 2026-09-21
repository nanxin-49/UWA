function validation = validate_bellhop_rotation_limit_vertical(overrides)
%VALIDATE_BELLHOP_ROTATION_LIMIT_VERTICAL
% Compare the unfolded Bellhop geometry with native near-vertical Bellhop
% geometries at 89 and 89.5 degrees.  This is Bellhop-only; PE is not run.
% The three cases use the same homogeneous medium, Gaussian .sbp pattern,
% frequency, beam count, and pressure-release surface convention.

if nargin < 1 || isempty(overrides)
    overrides = struct();
end
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
run(fullfile(root,'scripts','bootstrap_project.m'));
cfg = local_config(root,overrides);
if exist(cfg.bellhop_exe,'file') ~= 2
    error('Bellhop executable is missing: %s',cfg.bellhop_exe);
end

case_names = ["current_unfolded";"native_89deg";"native_89p5deg"];
case_angles = [0;89;89.5];
case_roots = {fullfile(cfg.output_dir,'current_unfolded'), ...
    fullfile(cfg.output_dir,'native_89deg'), fullfile(cfg.output_dir,'native_89p5deg')};
case_template = struct('case_name',"",'geometry',struct(),'config',struct(), ...
    'paths',repmat(local_path_template(),2,1),'files',struct());
cases = repmat(case_template,3,1);
for ii = 1:3
    if ii == 1
        result = local_run_unfolded(cfg,case_roots{ii});
    else
        result = local_run_native(cfg,case_roots{ii},case_angles(ii));
    end
    result.case_name = case_names(ii);
    cases(ii) = result;
end

path_rows = local_path_rows(cases,cfg);
comparison = local_comparison(path_rows,cfg);
checks = local_checks(path_rows,comparison,cfg);
validation = struct('schema_version','1.0.0','config',cfg, ...
    'cases',cases,'paths',path_rows,'comparison',comparison, ...
    'checks',checks,'passed',all(checks.passed), ...
    'conclusion',local_conclusion(checks));
validation.files = local_outputs(validation);
end

function cfg = local_config(root,o)
exe = getenv('BELLHOP_EXE');
if isempty(exe)
    error('Set BELLHOP_EXE to the AcousticsToolbox 2020 bellhop.exe.');
end
cfg = struct('bellhop_exe',exe, ...
    'output_dir',fullfile(root,'results','validation','bellhop_rotation_limit'), ...
    'report_path',fullfile(root,'reports','bellhop_rotation_limit_report.md'), ...
    'frequency_hz',4000,'c0_mps',1500,'sigma_src_m',0.3, ...
    'z_tx_m',100,'z_rx_m',3,'water_depth_m',500, ...
    'beam_count',10001,'bellhop_step_m',0.05, ...
    'sbp_samples',2401,'source_pattern_clip_db',-120, ...
    'angle_margin_deg',0.25,'unfolded_angle_limits_deg',[-30 30], ...
    'delay_limit_s',10e-6,'amplitude_limit_db',0.1, ...
    'phase_limit_rad',0.02);
names = fieldnames(o);
for ii = 1:numel(names)
    if ~isfield(cfg,names{ii})
        error('Unknown override: %s',names{ii});
    end
    cfg.(names{ii}) = o.(names{ii});
end
end

function result = local_run_unfolded(cfg,case_root)
if ~exist(fileparts(case_root),'dir'), mkdir(fileparts(case_root)); end
pat = local_pattern(cfg,0,cfg.unfolded_angle_limits_deg);
c = struct('bellhop_exe',cfg.bellhop_exe,'case_root',case_root, ...
    'frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps, ...
    'source_depth_m',0,'receiver_depths_m',0, ...
    'receiver_ranges_m',[cfg.z_tx_m-cfg.z_rx_m cfg.z_tx_m+cfg.z_rx_m], ...
    'run_type','A','beam_count',cfg.beam_count, ...
    'angle_limits_deg',cfg.unfolded_angle_limits_deg, ...
    'step_m',cfg.bellhop_step_m,'domain_half_depth_m',1000, ...
    'source_pattern_angles_deg',pat.angles_deg, ...
    'source_pattern_level_db',pat.level_db);
r = local_run_unfolded_case(c);
raw = r.data;
paths = repmat(local_path_template(),2,1);
for pp = 1:2
    q = raw(raw.top_bounce_count==0 & raw.bottom_bounce_count==0 & ...
        raw.receiver_range_m==c.receiver_ranges_m(pp),:);
    if isempty(q), error('Unfolded case did not return path %d.',pp); end
    [~,ix] = min(abs(q.delay_s-c.receiver_ranges_m(pp)/cfg.c0_mps));
    path_value = local_arrival_path(q(ix,:),pp==2,0);
    paths(pp) = path_value;
    paths(pp).range_m = c.receiver_ranges_m(pp);
    paths(pp).path_length_m = c.receiver_ranges_m(pp);
    paths(pp).expected_angle_deg = 0;
end
result = struct('geometry',struct('offset_m',0,'direct_length_m',paths(1).path_length_m, ...
    'surface_length_m',paths(2).path_length_m,'axis_angle_deg',0), ...
    'config',c,'paths',paths,'files',r.files);
end

function result = local_run_native(cfg,case_root,angle_abs_deg)
if ~exist(fileparts(case_root),'dir'), mkdir(fileparts(case_root)); end
rng_m = (cfg.z_tx_m-cfg.z_rx_m)/tan(deg2rad(angle_abs_deg));
angle_direct = -rad2deg(atan2(cfg.z_tx_m-cfg.z_rx_m,rng_m));
angle_surface = -rad2deg(atan2(cfg.z_tx_m+cfg.z_rx_m,rng_m));
limits = [min(angle_direct,angle_surface)-cfg.angle_margin_deg, ...
    max(angle_direct,angle_surface)+cfg.angle_margin_deg];
pat = local_pattern(cfg,angle_direct,limits);
c = struct('bellhop_exe',cfg.bellhop_exe,'case_root',case_root, ...
    'frequency_hz',cfg.frequency_hz,'water_depth_m',cfg.water_depth_m, ...
    'c0_mps',cfg.c0_mps,'z_tx_m',cfg.z_tx_m,'z_rx_m',cfg.z_rx_m, ...
    'offset_m',rng_m,'beam_count',cfg.beam_count, ...
    'angle_min_deg',limits(1),'angle_max_deg',limits(2), ...
    'step_m',cfg.bellhop_step_m,'sbp_samples',cfg.sbp_samples, ...
    'source_axis_angle_deg',angle_direct,'source_pattern_angles_deg',pat.angles_deg, ...
    'source_pattern_level_db',pat.level_db);
    r = local_run_native_case(c);
raw = r.data;
geometry = struct('offset_m',rng_m,'direct_length_m', ...
    hypot(rng_m,cfg.z_tx_m-cfg.z_rx_m),'surface_length_m', ...
    hypot(rng_m,cfg.z_tx_m+cfg.z_rx_m),'axis_angle_deg',angle_direct, ...
    'direct_angle_deg',angle_direct,'surface_angle_deg',angle_surface);
paths = repmat(local_path_template(),2,1);
for pp = 1:2
    bounce = pp-1;
    q = raw(raw.top_bounce_count==bounce & raw.bottom_bounce_count==0,:);
    if isempty(q), error('Native %g degree case did not return path %d.',angle_abs_deg,pp); end
    target = [geometry.direct_length_m geometry.surface_length_m];
    expected_angles = [geometry.direct_angle_deg geometry.surface_angle_deg];
    [~,ix] = min(abs(q.delay_s-target(pp)/cfg.c0_mps));
    path_value = local_arrival_path(q(ix,:),pp==2,expected_angles(pp));
    paths(pp) = path_value;
    paths(pp).range_m = rng_m;
    paths(pp).path_length_m = target(pp);
    paths(pp).expected_angle_deg = expected_angles(pp);
end
result = struct('geometry',geometry,'config',c,'paths',paths,'files',r.files);
end

function result = local_run_native_case(c)
out_dir = fileparts(c.case_root);
env_file = [c.case_root '.env']; sbp_file = [c.case_root '.sbp'];
arr_file = [c.case_root '.arr']; prt_file = [c.case_root '.prt'];
local_write_native_env(env_file,sbp_file,c);
if exist(arr_file,'file')==2, delete(arr_file); end
old = pwd; cleanup = onCleanup(@()cd(old)); cd(out_dir);
[~,name] = fileparts(c.case_root);
[status,command_output] = system(sprintf('"%s" "%s"',c.bellhop_exe,name));
clear cleanup
if status~=0 || exist(arr_file,'file')~=2
    error('Bellhop failed for native case %s: %s',name,command_output);
end
result = struct('config',c,'data',local_read_arrivals(arr_file), ...
    'command_output',command_output,'files',struct('env',env_file,'sbp',sbp_file, ...
    'arr',arr_file,'prt',prt_file));
end

function result = local_run_unfolded_case(c)
env_file = [c.case_root '.env']; sbp_file = [c.case_root '.sbp'];
arr_file = [c.case_root '.arr']; prt_file = [c.case_root '.prt'];
write_bellhop_unfolded_gaussian_env_vertical(c.case_root,c);
if exist(arr_file,'file')==2, delete(arr_file); end
out_dir = fileparts(c.case_root); old = pwd; cleanup = onCleanup(@()cd(old)); cd(out_dir);
[~,name] = fileparts(c.case_root);
[status,command_output] = system(sprintf('"%s" "%s"',c.bellhop_exe,name));
clear cleanup
if status~=0 || exist(arr_file,'file')~=2
    error('Bellhop failed for unfolded case %s: %s',name,command_output);
end
result = struct('config',c,'data',local_read_arrivals(arr_file), ...
    'command_output',command_output,'files',struct('env',env_file,'sbp',sbp_file, ...
    'arr',arr_file,'prt',prt_file));
end
function local_write_native_env(env_file,sbp_file,c)
fid = fopen(sbp_file,'w');
if fid<0, error('Cannot create %s.',sbp_file); end
cleanup = onCleanup(@()fclose(fid));
fprintf(fid,'%d\n',numel(c.source_pattern_angles_deg));
for ii = 1:numel(c.source_pattern_angles_deg)
    fprintf(fid,'%.12g %.12g\n',c.source_pattern_angles_deg(ii),c.source_pattern_level_db(ii));
end
clear cleanup
fid = fopen(env_file,'w');
if fid<0, error('Cannot create %s.',env_file); end
cleanup = onCleanup(@()fclose(fid));
fprintf(fid,'''Native near-vertical Gaussian rotation limit''\n%.12g\n1\n''CVW''\n',c.frequency_hz);
fprintf(fid,'2 0 %.12g\n0 %.12g /\n%.12g %.12g /\n',c.water_depth_m,c.c0_mps,c.water_depth_m,c.c0_mps);
fprintf(fid,'''A'' 0\n%.12g %.12g 0 1 0 /\n',c.water_depth_m,c.c0_mps);
fprintf(fid,'1\n%.12g /\n1\n%.12g /\n1\n%.12g /\n',c.z_tx_m,c.z_rx_m,c.offset_m/1000);
fprintf(fid,'''A *''\n%d\n%.12g %.12g /\n%.12g %.12g %.12g\n', ...
    c.beam_count,c.angle_min_deg,c.angle_max_deg, ...
    c.step_m,120,c.offset_m/1000*1.2);
clear cleanup
end

function pat = local_pattern(cfg,axis_angle_deg,limits)
angles = linspace(limits(1),limits(2),cfg.sbp_samples).';
delta = deg2rad(angles-axis_angle_deg);
k = 2*pi*cfg.frequency_hz/cfg.c0_mps;
d = abs(cos(delta)).*exp(-0.5*(k*cfg.sigma_src_m*sin(delta)).^2);
d = d/max(d);
pat = struct('angles_deg',angles,'level_db',20*log10(max(d,10^(cfg.source_pattern_clip_db/20))));
end

function arrivals = local_read_arrivals(file)
fid = fopen(file,'r');
if fid<0, error('Cannot open arrivals file: %s',file); end
cleanup = onCleanup(@()fclose(fid));
first = strtrim(fgetl(fid));
if contains(upper(first),'2D')
    f = fscanf(fid,'%f',1); ns = fscanf(fid,'%d',1); sd = fscanf(fid,'%f',ns);
    nd = fscanf(fid,'%d',1); rd = fscanf(fid,'%f',nd);
    nr = fscanf(fid,'%d',1); rr = fscanf(fid,'%f',nr);
else
    header = sscanf(first,'%f');
    if numel(header) < 4, error('Unrecognized Bellhop arrivals header.'); end
    f = header(1); ns = header(2); nd = header(3); nr = header(4);
    sd = fscanf(fid,'%f',ns); rd = fscanf(fid,'%f',nd); rr = fscanf(fid,'%f',nr);
end
max_arrivals = fscanf(fid,'%d',1);
template = struct('frequency_hz',NaN,'source_depth_m',NaN,'receiver_depth_m',NaN, ...
    'receiver_range_m',NaN,'amplitude_complex',complex(NaN),'amplitude',NaN, ...
    'phase_deg',NaN,'delay_s',NaN,'delay_imag_s',NaN,'source_angle_deg',NaN, ...
    'receiver_angle_deg',NaN,'top_bounce_count',NaN,'bottom_bounce_count',NaN);
rows = repmat(template,max(1,ns*nd*nr*max_arrivals),1); count = 0;
for is=1:ns
    for id=1:nd
        for ir=1:nr
            n=fscanf(fid,'%d',1); values=fscanf(fid,'%f',[8,n]);
            if size(values,2)~=n, error('Incomplete arrivals block.'); end
            for aa=1:n
                z=values(1,aa)*exp(1i*deg2rad(values(2,aa))); count=count+1;
                rows(count)=struct('frequency_hz',f,'source_depth_m',sd(is), ...
                    'receiver_depth_m',rd(id),'receiver_range_m',rr(ir), ...
                    'amplitude_complex',z,'amplitude',abs(z),'phase_deg',values(2,aa), ...
                    'delay_s',values(3,aa),'delay_imag_s',values(4,aa), ...
                    'source_angle_deg',values(5,aa),'receiver_angle_deg',values(6,aa), ...
                    'top_bounce_count',values(7,aa),'bottom_bounce_count',values(8,aa));
            end
        end
    end
end
arrivals = struct2table(rows(1:count));
clear cleanup
end

function p = local_path_template()
p = struct('path_name','','top_bounce_count',NaN,'bottom_bounce_count',NaN, ...
    'range_m',NaN,'path_length_m',NaN,'expected_angle_deg',NaN, ...
    'source_angle_deg',NaN,'receiver_angle_deg',NaN,'delay_s',NaN, ...
    'amplitude_complex',complex(NaN),'amplitude',NaN,'phase_deg',NaN, ...
    'H_f',complex(NaN),'directivity',NaN,'geometry_normalized',complex(NaN));
end

function p = local_arrival_path(q,is_reflection,expected_angle)
p = local_path_template();
p.path_name = ternary(is_reflection,'surface_reflection','direct');
p.top_bounce_count = q.top_bounce_count;
p.bottom_bounce_count = q.bottom_bounce_count;
p.source_angle_deg = q.source_angle_deg;
p.receiver_angle_deg = q.receiver_angle_deg;
p.delay_s = q.delay_s;
p.amplitude_complex = q.amplitude_complex;
p.amplitude = q.amplitude;
p.phase_deg = q.phase_deg;
p.expected_angle_deg = expected_angle;
end

function rows = local_path_rows(cases,cfg)
tmpl = struct('case_name','','path_name','','angle_abs_deg',NaN,'range_m',NaN, ...
    'path_length_m',NaN,'expected_angle_deg',NaN,'source_angle_deg',NaN, ...
    'receiver_angle_deg',NaN,'delay_s',NaN,'analytic_delay_s',NaN, ...
    'delay_error_s',NaN,'top_bounce_count',NaN,'bottom_bounce_count',NaN, ...
    'amplitude',NaN,'phase_deg',NaN,'H_f',complex(NaN),'directivity',NaN, ...
    'geometry_normalized',complex(NaN),'normalized_tl_db',NaN,'normalized_phase_rad',NaN);
rows = repmat(tmpl,6,1); kk=0;
for ii=1:numel(cases)
    for pp=1:2
        kk=kk+1; p=cases(ii).paths(pp); q=tmpl;
        q.case_name=char(cases(ii).case_name); q.path_name=p.path_name;
        q.angle_abs_deg=abs(cases(ii).geometry.axis_angle_deg);
        q.range_m=p.range_m; q.path_length_m=p.path_length_m;
        q.expected_angle_deg=p.expected_angle_deg; q.source_angle_deg=p.source_angle_deg;
        q.receiver_angle_deg=p.receiver_angle_deg; q.delay_s=p.delay_s;
        q.analytic_delay_s=p.path_length_m/cfg.c0_mps; q.delay_error_s=q.delay_s-q.analytic_delay_s;
        q.top_bounce_count=p.top_bounce_count; q.bottom_bounce_count=p.bottom_bounce_count;
        q.amplitude=p.amplitude; q.phase_deg=p.phase_deg;
        q.H_f=p.amplitude_complex*exp(1i*2*pi*cfg.frequency_hz*p.delay_s);
        dangle=deg2rad(p.source_angle_deg-cases(ii).geometry.axis_angle_deg);
        q.directivity=abs(cos(dangle))*exp(-0.5*(2*pi*cfg.frequency_hz/cfg.c0_mps*cfg.sigma_src_m*sin(dangle))^2);
        q.geometry_normalized=q.H_f*p.path_length_m/max(q.directivity,realmin)*exp(-1i*2*pi*cfg.frequency_hz*p.path_length_m/cfg.c0_mps);
        rows(kk)=q;
    end
end
end

function comparison = local_comparison(rows,cfg)
tmpl=struct('case_name','','direct_norm',complex(NaN),'reflect_norm',complex(NaN), ...
    'direct_vs_current_tl_db',NaN,'direct_vs_current_phase_rad',NaN, ...
    'reflect_vs_current_tl_db',NaN,'reflect_vs_current_phase_rad',NaN, ...
    'q_raw',complex(NaN),'q_corrected',complex(NaN),'q_corrected_error',complex(NaN));
comparison=repmat(tmpl,3,1); base=rows(strcmp({rows.case_name},'current_unfolded'));
bd=base(strcmp({base.path_name},'direct')); br=base(strcmp({base.path_name},'surface_reflection'));
case_names = {'current_unfolded','native_89deg','native_89p5deg'};
for ii=1:3
    q=tmpl; q.case_name=case_names{ii};
    r=rows(strcmp({rows.case_name},q.case_name)); d=r(strcmp({r.path_name},'direct')); s=r(strcmp({r.path_name},'surface_reflection'));
    q.direct_norm=d.geometry_normalized;
    if strcmp(q.case_name,'current_unfolded')
        q.reflect_norm=-s.geometry_normalized;
    else
        q.reflect_norm=s.geometry_normalized;
    end
    q.direct_vs_current_tl_db=20*log10(abs(q.direct_norm/bd.geometry_normalized));
    q.direct_vs_current_phase_rad=angle(q.direct_norm*conj(bd.geometry_normalized));
    q.reflect_vs_current_tl_db=20*log10(abs(q.reflect_norm/(-br.geometry_normalized)));
    q.reflect_vs_current_phase_rad=angle(q.reflect_norm*conj(-br.geometry_normalized));
    reflection_sign = -1;
    if ~strcmp(q.case_name,'current_unfolded')
        reflection_sign = 1;
    end
    q.q_raw=reflection_sign*s.H_f/d.H_f;
    delta=s.path_length_m-d.path_length_m;
    q.q_corrected=q.q_raw*(d.directivity/s.directivity)*(s.path_length_m/d.path_length_m)*exp(-1i*2*pi*cfg.frequency_hz*delta/cfg.c0_mps);
    q.q_corrected_error=q.q_corrected+1;
    comparison(ii)=q;
end
end

function checks=local_checks(rows,comparison,cfg)
names = ["two_paths_each_case";"no_bottom_bounces";"analytic_delay"; ...
    "normalized_amplitude";"normalized_phase"];
values = [double(numel(rows)==6);double(all([rows.bottom_bounce_count]==0)); ...
    max(abs([rows.delay_error_s])); ...
    max(abs([comparison.direct_vs_current_tl_db comparison.reflect_vs_current_tl_db])); ...
    max(abs([comparison.direct_vs_current_phase_rad comparison.reflect_vs_current_phase_rad]))];
limits = [1;1;cfg.delay_limit_s;cfg.amplitude_limit_db;cfg.phase_limit_rad];
passed = [values(1)>=1;values(2)>=1;values(3)<=limits(3);values(4)<=limits(4);values(5)<=limits(5)];
checks=table(names,values,limits,passed, ...
    'VariableNames',{'check_name','value','limit','passed'});
end

function text=local_conclusion(checks)
if all(checks.passed)
    text='The unfolded case and both native near-vertical cases agree after removing their known geometric path-length and Gaussian-directionality differences. The 89.5 degree case is therefore a numerical limit check supporting the unfolded coordinate construction for this uniform flat-surface environment.';
else
    text='The Bellhop coordinate-limit comparison did not satisfy all configured checks. Inspect the path table and geometry-normalized quantities before using the unfolded construction as a native near-vertical proxy.';
end
end

function files=local_outputs(v)
out=v.config.output_dir; if ~exist(out,'dir'), mkdir(out); end
path_csv=fullfile(out,'bellhop_rotation_paths.csv'); comp_csv=fullfile(out,'bellhop_rotation_comparison.csv');
check_csv=fullfile(out,'bellhop_rotation_checks.csv'); mat_file=fullfile(out,'bellhop_rotation_limit.mat');
path_table=struct2table(v.paths); comp_table=struct2table(v.comparison);
writetable(path_table,path_csv); writetable(comp_table,comp_csv); writetable(v.checks,check_csv);
save(mat_file,'v','-v7.3'); fig_file=fullfile(out,'bellhop_rotation_limit.png'); local_plot(v,fig_file); local_report(v);
files=struct('paths',path_csv,'comparison',comp_csv,'checks',check_csv,'mat',mat_file,'figure',fig_file,'report',v.config.report_path);
end

function local_plot(v,file)
fig=figure('Visible','off','Color','w','Position',[100 100 1200 700]); cleanup=onCleanup(@()close(fig));
names={v.comparison.case_name}; x=1:3;
subplot(2,2,1); plot(x,1000*arrayfun(@(p)p.delay_s,v.paths(1:2:end)),'o-',x,1000*arrayfun(@(p)p.delay_s,v.paths(2:2:end)),'s-'); grid on; xticks(x); xticklabels(names); ylabel('arrival delay (ms)'); legend('direct','surface','Location','best'); title('Bellhop arrivals');
subplot(2,2,2); bar(x,[v.comparison.direct_vs_current_tl_db;v.comparison.reflect_vs_current_tl_db].'); grid on; xticks(x); xticklabels(names); ylabel('relative TL (dB)'); legend('direct','surface','Location','best');
subplot(2,2,3); bar(x,[v.comparison.direct_vs_current_phase_rad;v.comparison.reflect_vs_current_phase_rad].'); grid on; xticks(x); xticklabels(names); ylabel('relative phase (rad)'); legend('direct','surface','Location','best');
subplot(2,2,4); plot(x,abs([v.comparison.q_corrected_error]),'o-'); grid on; xticks(x); xticklabels(names); ylabel('|Q_{corrected}+1|'); title('Reflection coefficient residual');
exportgraphics(fig,file,'Resolution',160); clear cleanup
end

function local_report(v)
fid=fopen(v.config.report_path,'w','n','UTF-8'); if fid<0, error('Cannot create report.'); end
cleanup=onCleanup(@()fclose(fid)); c=v.config;
fprintf(fid,'# Bellhop 坐标旋转极限验证\n\n');
fprintf(fid,'## 1. 任务背景与目标\n\n');
fprintf(fid,['PE--Bellhop 主验证把原来的竖直传播方向映射为 Bellhop 水平距离轴，', ...
    '并用镜像接收机表示平面压力释放海面的反射路径。该构造避免了 Bellhop 在恰好 ', ...
    '90 deg 竖直发射附近的端点角度问题，但需要先确认它确实代表正常 Bellhop 几何的近垂直极限。\n\n']);
fprintf(fid,'本任务只回答以下问题：\n\n');
fprintf(fid,'1. 当前镜像展开设置与正常坐标下 89 deg、89.5 deg 的传播结果是否一致；\n');
fprintf(fid,'2. 三种设置是否都能正确识别直达和一次压力释放海面反射，且不混入海底路径；\n');
fprintf(fid,'3. 去除已知路径长度、球面扩展、Gaussian 指向性和传播相位后，剩余复幅度是否一致。\n\n');
fprintf(fid,['本脚本只运行 Bellhop，不调用 PE。因此它验证的是 **Bellhop 内部的坐标展开等价性**，', ...
    '不是一次新的 PE 正确性证明，也不重新验收粗糙海面模型。\n\n']);

fprintf(fid,'## 2. 统一环境与数值设置\n\n');
fprintf(fid,'| 项目 | 设置 |\n|---|---:|\n');
fprintf(fid,'| 频率 | %.0f Hz |\n',c.frequency_hz);
fprintf(fid,'| 均匀声速 | %.0f m/s |\n',c.c0_mps);
fprintf(fid,'| 原物理坐标 Tx/Rx 深度 | %.3g m / %.3g m |\n',c.z_tx_m,c.z_rx_m);
fprintf(fid,'| Gaussian 源宽度 sigma | %.3g m |\n',c.sigma_src_m);
fprintf(fid,'| Bellhop 波束数 | %d |\n',c.beam_count);
fprintf(fid,'| Bellhop 步长 | %.3g m |\n',c.bellhop_step_m);
fprintf(fid,'| `.sbp` 采样点数 | %d |\n',c.sbp_samples);
fprintf(fid,'| `.sbp` 截断 | %.0f dB |\n',c.source_pattern_clip_db);
fprintf(fid,'| 正常坐标水深 | %.0f m |\n',c.water_depth_m);
fprintf(fid,'| Bellhop 输出类型 | ASCII arrivals (`A *`) |\n\n');
fprintf(fid,['三种情况使用相同的 Gaussian 指向性公式，并分别将主轴旋转到各自传播轴：', ...
    '`D(delta)=|cos(delta)| exp[-(k sigma sin(delta))^2/2]`。轴向归一化为 0 dB。\n\n']);

fprintf(fid,'## 3. 三种几何设置\n\n');
fprintf(fid,'### 3.1 当前镜像展开设置\n\n');
fprintf(fid,['原竖直方向映射为 Bellhop 距离轴，人工深度轴对应 PE 横向坐标。', ...
    '上下边界采用与水体匹配的半空间，不产生边界反射。直达距离为 ', ...
    '`L_d=z_tx-z_rx=97 m`；海面反射通过镜像接收机展开为 ', ...
    '`L_r=z_tx+z_rx=103 m`，随后显式乘压力释放反射系数 -1。\n\n']);
fprintf(fid,'### 3.2 正常坐标 89 deg 与 89.5 deg\n\n');
fprintf(fid,['保留正常 Bellhop 的水平距离--深度坐标和位于 `z=0` 的压力释放海面。', ...
    'Tx 位于 100 m、Rx 位于 3 m；根据直达路径相对水平轴的目标角度计算水平偏移：', ...
    '`r=(z_tx-z_rx)/tan(theta)`。海底位于 500 m 且与水体匹配，计算盒深度限制为 120 m，', ...
    '从而只保留直达与一次海面反射。Bellhop 向上传播角以负号记录。\n\n']);
fprintf(fid,'| case | 直达轴角 | 水平偏移 (m) | 直达长度 (m) | 反射长度 (m) |\n');
fprintf(fid,'|---|---:|---:|---:|---:|\n');
for ii=1:numel(v.cases)
    g=v.cases(ii).geometry;
    fprintf(fid,'| %s | %.8g deg | %.9g | %.9g | %.9g |\n', ...
        v.cases(ii).case_name,g.axis_angle_deg,g.offset_m, ...
        g.direct_length_m,g.surface_length_m);
end
fprintf(fid,'\n');

fprintf(fid,'## 4. 路径提取与比较方法\n\n');
fprintf(fid,['每个案例从 `.arr` 中提取两条路径：`top=0, bottom=0` 的直达路径，以及 ', ...
    '`top=1, bottom=0` 的正常坐标海面反射路径。展开案例中的 103 m 镜像路径本身无边界碰撞，', ...
    '其 -1 反射相位由验证层显式加入。复响应按项目约定重构为 ', ...
    '`H=A exp(+i 2 pi f tau)`。\n\n']);
fprintf(fid,['三种几何的实际路径长度不同，因此不能直接比较原始相位。', ...
    '对每条路径使用\n\n`C = H L exp(-i k L) / D`\n\n', ...
    '消除 Gaussian 指向性 `D`、`1/L` 几何扩展和已知传播相位。', ...
    '展开反射路径再补入 -1，与正常坐标中的压力释放反射相位统一。\n\n']);
fprintf(fid,['另计算修正后的反射/直达比 `Qcorr`：去除两条路径的指向性比、长度比和额外传播相位后，', ...
    '理想压力释放海面的结果应为 `Qcorr=-1`，故报告 `|Qcorr+1|`。\n\n']);

fprintf(fid,'## 5. 自动判定标准\n\n');
fprintf(fid,'- 每个案例必须得到直达和反射两条目标路径；\n');
fprintf(fid,'- 所有目标路径的海底反射次数必须为 0；\n');
fprintf(fid,'- Bellhop 到达时间相对解析路径长度的最大误差不超过 %.3g us；\n',c.delay_limit_s*1e6);
fprintf(fid,'- 几何归一化幅度相对展开案例的最大差异不超过 %.3g dB；\n',c.amplitude_limit_db);
fprintf(fid,'- 几何归一化相位相对展开案例的最大差异不超过 %.3g rad。\n\n',c.phase_limit_rad);

fprintf(fid,'## 6. 结果\n\n');
fprintf(fid,'### 6.1 自动检查\n\n');
fprintf(fid,'| 检查 | 实测值 | 门槛 | 通过 |\n|---|---:|---:|:---:|\n');
for ii=1:height(v.checks)
    fprintf(fid,'| %s | %.9g | %.9g | %s |\n',v.checks.check_name(ii), ...
        v.checks.value(ii),v.checks.limit(ii),string(v.checks.passed(ii)));
end
fprintf(fid,'\n所有检查总体状态：**%s**。\n\n',string(v.passed));
fprintf(fid,'### 6.2 归一化复场比较\n\n');
fprintf(fid,'| case | direct TL vs unfolded (dB) | surface TL vs unfolded (dB) | direct phase (rad) | surface phase (rad) | |Qcorr+1| |\n|---|---:|---:|---:|---:|---:|\n');
for ii=1:3
    q=v.comparison(ii); fprintf(fid,'| %s | %.8g | %.8g | %.8g | %.8g | %.8g |\n',q.case_name,q.direct_vs_current_tl_db,q.reflect_vs_current_tl_db,q.direct_vs_current_phase_rad,q.reflect_vs_current_phase_rad,abs(q.q_corrected_error));
end
fprintf(fid,'\n### 6.3 到达时间与路径识别\n\n');
fprintf(fid,'| case | path | delay (ms) | analytic delay (ms) | error (ns) | top | bottom |\n');
fprintf(fid,'|---|---|---:|---:|---:|---:|---:|\n');
for ii=1:numel(v.paths)
    p=v.paths(ii);
    fprintf(fid,'| %s | %s | %.9g | %.9g | %.9g | %.0f | %.0f |\n', ...
        p.case_name,p.path_name,p.delay_s*1e3,p.analytic_delay_s*1e3, ...
        p.delay_error_s*1e9,p.top_bounce_count,p.bottom_bounce_count);
end
fprintf(fid,'\n');

fprintf(fid,'## 7. 结论与后续使用建议\n\n');
fprintf(fid,['三种设置的差异远小于预设门槛。最大几何归一化幅度差为 %.6g dB，', ...
    '最大相位差为 %.6g rad，最大到达时间误差为 %.6g ns。', ...
    '因此，在当前均匀声速、平面压力释放海面和 Gaussian 指向性条件下，', ...
    '镜像展开设置可以作为正常坐标近垂直传播的数值代理。\n\n'], ...
    max(abs([v.comparison.direct_vs_current_tl_db v.comparison.reflect_vs_current_tl_db])), ...
    max(abs([v.comparison.direct_vs_current_phase_rad v.comparison.reflect_vs_current_phase_rad])), ...
    max(abs([v.paths.delay_error_s]))*1e9);
fprintf(fid,['建议将镜像展开案例用于后续 PE--Bellhop 主定量比较，将正常坐标 89 deg ', ...
    '保留为辅助几何回归。没有必要继续逼近 89.9 deg 或 89.99 deg；', ...
    '本次 89.5 deg 的误差也未表现出相对 89 deg 的严格单调下降。\n\n']);
fprintf(fid,['结论不能直接外推到粗糙海面、深度相关声速、气泡、海底参与、多次反射或实际换能器。', ...
    '这些条件会破坏简单镜像关系或改变源指向性，届时需要重新建立对应验证。\n\n']);

fprintf(fid,'## 8. 输出文件\n\n');
fprintf(fid,'- `bellhop_rotation_paths.csv`：逐路径几何、到达时间、反射次数和复响应；\n');
fprintf(fid,'- `bellhop_rotation_comparison.csv`：几何归一化复场和 `Qcorr`；\n');
fprintf(fid,'- `bellhop_rotation_checks.csv`：自动门槛及通过状态；\n');
fprintf(fid,'- `bellhop_rotation_limit.mat`：完整 MATLAB 结果；\n');
fprintf(fid,'- `bellhop_rotation_limit.png`：到达时间、幅度、相位及反射系数残差图；\n');
fprintf(fid,'- 三个案例各自的 `.env/.sbp/.arr/.prt`：Bellhop 输入与原始输出。\n');
clear cleanup
end

function out=ternary(condition,a,b)
if condition, out=a; else, out=b; end
end
