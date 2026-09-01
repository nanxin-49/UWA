function audit = validate_bellhop_curved_wall_beam_frame_audit(options)
%VALIDATE_BELLHOP_CURVED_WALL_BEAM_FRAME_AUDIT Locate curved-wall phase bias.
% This validation-only audit compares one weak sinusoid in native ATI and in
% the internal-wall overlay. It logs Reflect2D state, the proper-rotation
% frame, and each accepted InfluenceGeoHatCart contribution at one receiver.
arguments
    options.output_dir (1,:) char = ''
    options.audit_exe (1,:) char = ''
    options.frequency_hz (1,1) double {mustBePositive} = 4000
    options.wall_amplitude_m (1,1) double {mustBePositive} = 0.25
    options.wall_wavenumber_per_m (1,1) double = -0.01
    options.profile_count (1,1) double {mustBeInteger,mustBePositive} = 161
    options.beam_count (1,1) double {mustBeInteger,mustBePositive} = 5001
    options.step_m (1,1) double {mustBePositive} = 0.1
    options.write_report (1,1) logical = true
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup_vertical_project;
if isempty(options.output_dir)
    options.output_dir=fullfile(root,'results','validation','bellhop_curved_wall_beam_frame_audit');
end
if isempty(options.audit_exe)
    options.audit_exe=fullfile(options.output_dir,'bin','bellhop_curved_wall_audit_2020.exe');
end
assert(isfile(options.audit_exe),'Missing curved-wall audit binary: %s',options.audit_exe);
if ~exist(options.output_dir,'dir'), mkdir(options.output_dir); end

R0=100; zc=0; K=options.wall_wavenumber_per_m; A=options.wall_amplitude_m;
zProf=linspace(zc-100,zc+100,options.profile_count).';
rProf=R0-A*sin(K*zProf);
cfg=struct('c0_mps',1500,'wall_r0_m',R0,'wall_amplitude_m',A, ...
    'wall_wavenumber_per_m',K,'wall_profile_r_m',rProf,'wall_profile_z_m',zProf, ...
    'mapped_receiver_range_m',103,'native_receiver_range_m',97, ...
    'source_depth_m',zc,'receiver_depths_m',zc,'receiver_ranges_m',[102 103], ...
    'run_type','C','angle_limits_deg',[-30 30],'domain_half_depth_m',1000, ...
    'sigma_src_m',0.3,'source_pattern_clip_db',-120,'frequency_hz',options.frequency_hz, ...
    'step_m',options.step_m,'beam_count',options.beam_count);
pat=local_pattern(cfg.frequency_hz,cfg);
cfg.source_pattern_angles_deg=pat.angles_deg;
cfg.source_pattern_level_db=pat.level_db;
cfg.bellhop_exe=options.audit_exe;
rot=cfg; rot.audit_mode='rotated'; rot.audit_receiver_range_m=cfg.mapped_receiver_range_m; rot.audit_receiver_depth_m=zc;
rot.case_root=fullfile(options.output_dir,'cases','rotated',sprintf('A_%g_N_%d_S_%g_B_%d',A,options.profile_count,options.step_m,options.beam_count));
nat=cfg; nat.audit_mode='native'; nat.audit_receiver_range_m=cfg.native_receiver_range_m; nat.audit_receiver_depth_m=zc;
nat.receiver_ranges_m=[cfg.native_receiver_range_m 98];
nat.case_root=fullfile(options.output_dir,'cases','native',sprintf('A_%g_N_%d_S_%g_B_%d',A,options.profile_count,options.step_m,options.beam_count));
rotated=run_bellhop_curved_wall_beam_frame_audit_vertical(rot);
native=run_bellhop_curved_wall_beam_frame_audit_vertical(nat);
% A direct-only native-range run supplies the reference needed to isolate
% the reflected field from the coherent total.  It has no audit sidecar and
% therefore does not alter the validation binary or any Bellhop physics.
direct=nat;
direct.audit_mode='disabled';
direct.case_root=fullfile(options.output_dir,'cases','direct',sprintf('A_%g_N_%d_S_%g_B_%d',A,options.profile_count,options.step_m,options.beam_count));
direct=run_bellhop_unfolded_gaussian_vertical(direct);

R=rotated.audit.reflection; N=native.audit.reflection;
[common,reflection_rows]=local_pair_reflections(R,N,[-0.03 0 0.03]);
if isempty(common), error('No common native/rotated TOP reflection records were found.'); end

RC=rotated.audit.contribution; NC=native.audit.contribution;
fan=local_match_contributions(RC,NC);
if isempty(fan), error('No common native/rotated receiver contributions were found.'); end
fan_rows=fan;
fan_rows.phase_difference_rad=angle(fan_rows.rotated_contribution.*conj(fan_rows.native_contribution));
fan_rows.amplitude_ratio=abs(fan_rows.rotated_contribution)./max(abs(fan_rows.native_contribution),realmin);
fan_rows.is_selected=ismember(round(fan_rows.alpha_deg,10),round([-0.03 0 0.03],10));

sum_rot=sum(fan_rows.rotated_contribution); sum_nat=sum(fan_rows.native_contribution);
summary=struct;
summary.sum_phase_difference_rad=angle(sum_rot*conj(sum_nat));
summary.sum_tl_difference_db=20*log10(abs(sum_rot)/max(abs(sum_nat),realmin));
summary.fan_phase_min=min(fan_rows.phase_difference_rad);
summary.fan_phase_max=max(fan_rows.phase_difference_rad);
summary.fan_phase_median=median(fan_rows.phase_difference_rad);
summary.fan_amp_ratio_min=min(fan_rows.amplitude_ratio);
summary.fan_amp_ratio_max=max(fan_rows.amplitude_ratio);
summary.rotation=local_rotation_summary(rotated.audit.rotation);
% Compare the correctly indexed coherent reflected field.  The rotated
% receiver is at range 103 m; native total and direct runs are both sampled
% at range 97 m, so subtraction does not mix receiver columns.
[rot_field,rot_range_meta]=select_bellhop_shd_pressure_at_range_vertical( ...
    rotated.data,cfg.mapped_receiver_range_m,zc);
[nat_total,native_range_meta]=select_bellhop_shd_pressure_at_range_vertical( ...
    native.data,cfg.native_receiver_range_m,zc);
[nat_direct,direct_range_meta]=select_bellhop_shd_pressure_at_range_vertical( ...
    direct.data,cfg.native_receiver_range_m,zc);
nat_reflect=nat_total-nat_direct;
summary.rotated_reflected_field=rot_field;
summary.native_total_field=nat_total;
summary.native_direct_field=nat_direct;
summary.native_reflected_field=nat_reflect;
summary.field_phase_difference_rad=angle(rot_field*conj(nat_reflect));
summary.field_tl_difference_db=20*log10(abs(rot_field)/max(abs(nat_reflect),realmin));
summary.field_complex_abs_error=abs(rot_field-nat_reflect);
summary.field_complex_relative_error=abs(rot_field-nat_reflect)/max(abs(nat_reflect),realmin);
summary.range_selection=struct('rotated',rot_range_meta,'native_total',native_range_meta, ...
    'native_direct',direct_range_meta);
summary.classification=local_classify(common,fan_rows,summary.sum_phase_difference_rad,summary.field_phase_difference_rad);

reflection_csv=fullfile(options.output_dir,'reflection_frame_comparison.csv'); writetable(reflection_rows,reflection_csv);
fan_csv=fullfile(options.output_dir,'receiver_contribution_fan.csv'); writetable(fan_rows,fan_csv);
mat_file=fullfile(options.output_dir,'curved_wall_beam_frame_audit.mat');
audit=struct('config',cfg,'options',options,'rotated',rotated,'native',native,'direct',direct, ...
    'common_reflection',common,'reflection_rows',reflection_rows,'fan_rows',fan_rows, ...
    'summary',summary,'files',struct('reflection_csv',reflection_csv,'fan_csv',fan_csv,'mat',mat_file));
save(mat_file,'audit','-v7.3');
if options.write_report
    local_write_report(audit,fullfile(root,'reports','bellhop_curved_wall_beam_frame_audit_report.md'));
end
disp(summary);
end

function pat=local_pattern(f,cfg)
angles=linspace(cfg.angle_limits_deg(1),cfg.angle_limits_deg(2),2401).';
k=2*pi*f/cfg.c0_mps; th=deg2rad(angles);
D=cos(th).*exp(-0.5*(k*cfg.sigma_src_m*sin(th)).^2); D=abs(D)/max(abs(D));
pat=struct('angles_deg',angles,'level_db',20*log10(max(D,10^(cfg.source_pattern_clip_db/20))));
end

function [common,rows]=local_pair_reflections(R,N,desired)
common=table; rows=table; matched=[];
% Native ATI tracing contains later top/bottom events as well.  For this
% audit retain the first TOP reflection per take-off angle, which is the
% one-reflection branch corresponding to the internal-wall ray.
seen=[]; first_idx=[];
for jj=1:height(N)
    if ~any(abs(seen-N.alpha_deg(jj))<1e-8)
        seen(end+1)=N.alpha_deg(jj); %#ok<AGROW>
        first_idx(end+1)=jj; %#ok<AGROW>
    end
end
Nfirst=N(first_idx,:);
for ii=1:height(R)
    [delta,jj]=min(abs(Nfirst.alpha_deg-R.alpha_deg(ii)));
    if delta<1e-8
        matched(end+1,:)=[ii first_idx(jj)]; %#ok<AGROW>
    end
end
if isempty(matched), return; end
vars={'alpha_deg','x_r','x_z','rayt_r','rayt_z','rayn_r','rayn_z','wallt_r','wallt_z','walln_r','walln_z','kappa','Tg','Th','RN_raw','RN_curv','RN_grad','RN_final','RM','p0_1','p0_2','q0_1','q0_2','p1_1','p1_2','q1_1','q1_2','amp0','amp1','phase0','phase1','tau0_r','tau0_i','tau1_r','tau1_i','ref_t_r','ref_t_z','ref_n_r','ref_n_z','handedness'};
alpha=R.alpha_deg(matched(:,1)); common=table(alpha,'VariableNames',{'alpha_deg'});
for kk=2:numel(vars)
    v=vars{kk}; d=R.(v)(matched(:,1))-N.(v)(matched(:,2)); common.([v '_error'])=d;
end
sel=zeros(numel(desired),1);
for ii=1:numel(desired), [~,jj]=min(abs(alpha-desired(ii))); sel(ii)=jj; end
rows=common(sel,:);
end

function fan=local_match_contributions(R,N)
uR=unique(R.alpha_deg); uN=unique(N.alpha_deg); alpha=[]; cr=[]; cn=[]; nr=[]; nn=[];
for ii=1:numel(uR)
    [delta,jj]=min(abs(uN-uR(ii)));
    if delta<1e-8
        mr=abs(R.alpha_deg-uR(ii))<1e-8; mn=abs(N.alpha_deg-uN(jj))<1e-8;
        alpha(end+1,1)=uR(ii); cr(end+1,1)=sum(R.contribution_r(mr)+1i*R.contribution_i(mr)); %#ok<AGROW>
        cn(end+1,1)=sum(N.contribution_r(mn)+1i*N.contribution_i(mn)); %#ok<AGROW>
        nr(end+1,1)=sum(mr); nn(end+1,1)=sum(mn); %#ok<AGROW>
    end
end
fan=table(alpha,cr,cn,nr,nn,'VariableNames',{'alpha_deg','rotated_contribution','native_contribution','rotated_segments','native_segments'});
end

function summary=local_rotation_summary(T)
summary=struct('record_count',height(T),'max_handedness_error',NaN,'q_sign_values',[]);
if ~isempty(T)
    summary.max_handedness_error=max(abs(T.handedness-1));
    summary.q_sign_values=unique(T.q_sign);
end
end

function class=local_classify(reflection,fan,~,field_phase)
% Use tolerances in the native state scales.  ATI input quantization gives
% micrometre-level x and O(1e-3) absolute q differences on q~1.5e5; those
% are not a frame mismatch.  Geometry/frame and phase errors remain tight.
geo_names={'x_r_error','x_z_error','wallt_r_error','wallt_z_error','walln_r_error','walln_z_error','kappa_error','Tg_error','Th_error','RN_raw_error','RN_curv_error','RN_grad_error','RN_final_error','RM_error','ref_t_r_error','ref_t_z_error','ref_n_r_error','ref_n_z_error'};
p_names={'p0_1_error','p0_2_error','p1_1_error','p1_2_error'};
q_names={'q0_1_error','q0_2_error','q1_1_error','q1_2_error'};
phase_names={'phase0_error','phase1_error'};
tau_names={'tau0_r_error','tau0_i_error','tau1_r_error','tau1_i_error'};
mx_geo=local_max_error(reflection,geo_names);
mx_p=local_max_error(reflection,p_names);
mx_q=local_max_error(reflection,q_names);
mx_phase=local_max_error(reflection,phase_names);
mx_tau=local_max_error(reflection,tau_names);
if mx_geo>2e-6 || mx_p>1e-3 || mx_q>1e-2 || mx_phase>1e-8 || mx_tau>1e-8
    class='A: reflection-frame or p/q mismatch at Reflect2D exit';
elseif any(abs(fan.phase_difference_rad)>0.1)
    class='B: Reflect2D exit agrees, but InfluenceGeoHatCart contribution phase differs';
elseif abs(field_phase)>0.1
    class='C: single-ray contributions agree; discrepancy appears only after coherent accumulation';
else
    class='C (bookkeeping/indexing-only): no physical A/B/C mismatch after correct receiver pairing';
end
end

function mx=local_max_error(T,names)
mx=0;
for ii=1:numel(names)
    if ismember(names{ii},T.Properties.VariableNames)
        x=T.(names{ii});
        if ~isempty(x), mx=max(mx,max(abs(x))); end
    end
end
end

function local_write_report(audit,path)
s=audit.summary; f=fopen(path,'w'); if f<0, error('Cannot create %s.',path); end
cleanup=onCleanup(@()fclose(f));
fprintf(f,'# Bellhop curved-wall beam/frame covariance audit\n\n');
fprintf(f,'## Result\n\nClassification: **%s**\n\n',s.classification);
fprintf(f,'Fixed case: `A=%.6g m`, `K=%.6g 1/m`, profile count `%d`, step `%.6g m`, `%d` beams at `%.6g Hz`. No PE or PM-wall code was changed.\n\n',audit.config.wall_amplitude_m,audit.config.wall_wavenumber_per_m,audit.config.profile_count,audit.config.step_m,audit.config.beam_count,audit.config.frequency_hz);
fprintf(f,'## Reflect2D出口\n\n');
fprintf(f,'Native and rotated reflection records are paired by take-off angle. The three representative rays are stored in `reflection_frame_comparison.csv`, including wall frame, Tg/Th, RN/RM, p/q before and after reflection, reflected ray frame, Amp, Phase, and tau.\n\n');
fprintf(f,'## Proper rotation\n\n');
fprintf(f,'The rotated logger records transformed tangent/ray-normal, p/q, Amp, Phase, tau, handedness, and q sign. No beam state was reset or fitted.\n\n');
fprintf(f,'## InfluenceGeoHatCart\n\n');
fprintf(f,'Accepted target-receiver contributions are aggregated by ray in `receiver_contribution_fan.csv`. Single-ray phase-difference range is `%.9g` to `%.9g` rad (median `%.9g`); amplitude-ratio range is `%.9g` to `%.9g`. The summed logged contribution phase difference is `%.9g` rad and its pre-ScalePressure TL difference is `%.9g dB`.\n\n',s.fan_phase_min,s.fan_phase_max,s.fan_phase_median,s.fan_amp_ratio_min,s.fan_amp_ratio_max,s.sum_phase_difference_rad,s.sum_tl_difference_db);
fprintf(f,'Correctly indexed SHD reflected-field comparison (rotated range 103 m versus native total-minus-direct at 97 m) gives phase difference `%.9g` rad and TL difference `%.9g dB`.\n\n',s.field_phase_difference_rad,s.field_tl_difference_db);
fprintf(f,'Its absolute complex difference is `%.9g` (relative `%.9g`).\n\n',s.field_complex_abs_error,s.field_complex_relative_error);
fprintf(f,'No 2.095-rad correction, amplitude calibration, curvature modification, p/q reset, or receiver normalization change was applied.\n\n');
fprintf(f,'## Decision\n\n');
if startsWith(s.classification,'A')
    fprintf(f,'The mismatch is already present at the Reflect2D exit.\n');
elseif startsWith(s.classification,'B')
    fprintf(f,'Reflect2D exit state is covariant; the first mismatch appears in receiver-side InfluenceGeoHatCart contribution. The remaining suspect range is segment interpolation, ray-frame orientation, KMAH/caustic bookkeeping, or receiver-side Cartesian geometry.\n');
else
    fprintf(f,'Single-ray and corrected coherent reflected-field contributions are covariant. The earlier ~2.095 rad observation is classified as receiver-column/indexing bookkeeping, not a Bellhop frame or accumulation error.\n');
end
fprintf(f,'This audit localizes the issue only; it does not implement a final repair or PM wall.\n');
end
