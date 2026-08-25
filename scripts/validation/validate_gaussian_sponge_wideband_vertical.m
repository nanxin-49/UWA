function validation=validate_gaussian_sponge_wideband_vertical(overrides)
%VALIDATE_GAUSSIAN_SPONGE_WIDEBAND_VERTICAL Direct-path Gaussian H(f) audit.

if nargin<1, overrides=struct(); end
root=fileparts(fileparts(fileparts(mfilename('fullpath')))); addpath(root); setup_vertical_project();
cfg=local_defaults(root); cfg=local_overrides(cfg,overrides);
sponge_file=fullfile(cfg.output_dir,'gaussian_sponge_validation.mat');
if exist(sponge_file,'file')~=2, error('Run validate_gaussian_sponge_vertical first.'); end
d=load(sponge_file,'validation'); sponge=d.validation;

cases=struct( ...
    'name',{'large_no_sponge','production_default','recommended_no_sponge'}, ...
    'width_m',{cfg.reference_width_m,cfg.production_width_m,sponge.recommendation.width_m}, ...
    'ratio',{0.12,cfg.default_ratio,0.12}, ...
    'alpha',{0,cfg.default_alpha,0});
for ii=1:numel(cases)
    N=local_n(cases(ii).width_m,cfg.dx_target_m);
    out=vertical_channel_model(local_params(cfg,cases(ii),N));
    cases(ii).f_axis_hz=out.f_axis(:);
    cases(ii).H_direct_dsp_f=out.H_direct_f(:);
    cases(ii).H_physical_f=out.H_direct_physical_f(:);
    cases(ii).phase_physical_rad=unwrap(angle(cases(ii).H_physical_f));
    cases(ii).group_delay_s=gradient(cases(ii).phase_physical_rad,cases(ii).f_axis_hz)/(2*pi);
end
ref=cases(1); rows=repmat(local_row(),numel(cases)*numel(ref.f_axis_hz),1); q=0;
summary_rows=repmat(local_summary_row(),numel(cases),1);
for ii=1:numel(cases)
    phase_error=unwrap(angle(cases(ii).H_physical_f.*conj(ref.H_physical_f)));
    tl_error=-20*log10(abs(cases(ii).H_physical_f./ref.H_physical_f));
    gd_error=cases(ii).group_delay_s-ref.group_delay_s;
    for ff=1:numel(ref.f_axis_hz)
        q=q+1; rows(q)=local_row(); rows(q).case_name=string(cases(ii).name);
        rows(q).width_m=cases(ii).width_m; rows(q).sponge_ratio=cases(ii).ratio;
        rows(q).alpha_max_np_per_m=cases(ii).alpha; rows(q).frequency_hz=ref.f_axis_hz(ff);
        rows(q).H_physical=cases(ii).H_physical_f(ff); rows(q).magnitude_db=20*log10(abs(cases(ii).H_physical_f(ff)));
        rows(q).unwrapped_phase_rad=cases(ii).phase_physical_rad(ff); rows(q).group_delay_s=cases(ii).group_delay_s(ff);
        rows(q).tl_error_vs_reference_db=tl_error(ff); rows(q).phase_error_vs_reference_rad=phase_error(ff);
        rows(q).group_delay_error_vs_reference_s=gd_error(ff);
    end
    s=local_summary_row(); s.case_name=string(cases(ii).name); s.width_m=cases(ii).width_m;
    s.sponge_ratio=cases(ii).ratio; s.alpha_max_np_per_m=cases(ii).alpha;
    s.max_abs_tl_error_db=max(abs(tl_error)); s.rms_phase_error_rad=sqrt(mean(phase_error.^2));
    s.max_abs_phase_error_rad=max(abs(phase_error)); s.rms_group_delay_error_us=1e6*sqrt(mean(gd_error.^2));
    s.max_abs_group_delay_error_us=1e6*max(abs(gd_error));
    s.amplitude_error_span_db=max(tl_error)-min(tl_error); summary_rows(ii)=s;
end
frequency_table=struct2table(rows); summary_table=struct2table(summary_rows);
validation=struct('schema_version','1.0.0','config',cfg,'sponge_validation',sponge, ...
    'cases',cases,'frequency_table',frequency_table,'summary_table',summary_table, ...
    'group_delay_definition','tau=(1/(2*pi))*d(unwrap(angle(H_physical)))/df under p(t)=real(P exp(-i2pift))');
validation.files=local_save(validation);
end

function cfg=local_defaults(root)
cfg=struct('output_dir',fullfile(root,'results','validation','pe_gaussian_window_sponge'), ...
    'f_axis_hz',linspace(3000,5000,33),'f_ref_hz',4000,'c0_mps',1500, ...
    'sigma_src_m',0.3,'distance_m',97,'dx_target_m',50/256,'stepz_lamb',0.5, ...
    'reference_width_m',160,'production_width_m',50,'default_ratio',0.12,'default_alpha',0.15);
end
function cfg=local_overrides(cfg,o), n=fieldnames(o); for i=1:numel(n), if ~isfield(cfg,n{i}), error('Unknown override: %s',n{i}); end, cfg.(n{i})=o.(n{i}); end, end
function N=local_n(W,dx), N=2*round(W/dx/2); end
function p=local_params(cfg,c,N)
p=struct('f0',cfg.f_axis_hz,'f_ref_hz',cfg.f_ref_hz,'c0',cfg.c0_mps,'z_max',100, ...
    'z_tx',cfg.distance_m,'z_rx',0,'xw',c.width_m,'yw',c.width_m,'nx',N,'ny',N, ...
    'x_tx',0,'y_tx',0,'x_rx',0,'y_rx',0,'stepz_lamb',cfg.stepz_lamb, ...
    'sigma_src_m',cfg.sigma_src_m,'source_mode','gaussian','sponge_ratio',c.ratio, ...
    'alpha_max_np_per_m',c.alpha,'validation_allow_extended_sponge_ratio',true, ...
    'validation_allow_extended_window',c.width_m>100,'env_mode','uniform', ...
    'enable_surface_reflection',false,'enable_bubbles',false,'doppler_fn',[], ...
    'enforce_1_over_R',false,'show_figures',false,'save_mode','rx_only','use_gpu',false);
end
function r=local_row(), r=struct('case_name',"",'width_m',NaN,'sponge_ratio',NaN,'alpha_max_np_per_m',NaN,'frequency_hz',NaN,'H_physical',complex(NaN),'magnitude_db',NaN,'unwrapped_phase_rad',NaN,'group_delay_s',NaN,'tl_error_vs_reference_db',NaN,'phase_error_vs_reference_rad',NaN,'group_delay_error_vs_reference_s',NaN); end
function r=local_summary_row(), r=struct('case_name',"",'width_m',NaN,'sponge_ratio',NaN,'alpha_max_np_per_m',NaN,'max_abs_tl_error_db',NaN,'amplitude_error_span_db',NaN,'rms_phase_error_rad',NaN,'max_abs_phase_error_rad',NaN,'rms_group_delay_error_us',NaN,'max_abs_group_delay_error_us',NaN); end
function files=local_save(v)
out=v.config.output_dir; csv=fullfile(out,'gaussian_sponge_wideband.csv'); sumcsv=fullfile(out,'gaussian_sponge_wideband_summary.csv'); mat=fullfile(out,'gaussian_sponge_wideband.mat');
writetable(v.frequency_table,csv); writetable(v.summary_table,sumcsv); validation=v; save(mat,'validation','-v7.3'); files=struct('csv',csv,'summary_csv',sumcsv,'mat',mat);
end
