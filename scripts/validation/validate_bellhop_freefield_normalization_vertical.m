function validation = validate_bellhop_freefield_normalization_vertical(overrides)
%VALIDATE_BELLHOP_FREEFIELD_NORMALIZATION_VERTICAL Audit Bellhop source scale.
% Uses the matched-halfspace construction from the official omni.env case.

if nargin<1 || isempty(overrides), overrides=struct(); end
root=fileparts(fileparts(fileparts(mfilename('fullpath')))); addpath(root); setup_vertical_project();
cfg=local_defaults(root); cfg=local_overrides(cfg,overrides); local_validate(cfg);
if ~exist(cfg.output_dir,'dir'), mkdir(cfg.output_dir); end
rows=repmat(local_row(),numel(cfg.frequencies_hz)*numel(cfg.beam_counts)*numel(cfg.step_values_m),1);
q=0; runs=cell(size(rows));
for ff=1:numel(cfg.frequencies_hz)
    for bb=1:numel(cfg.beam_counts)
        for st=1:numel(cfg.step_values_m)
            q=q+1; f=cfg.frequencies_hz(ff); beams=cfg.beam_counts(bb); step=cfg.step_values_m(st);
            tag=sprintf('f%g_b%d_s%d',f,beams,round(1000*step));
            c=struct('bellhop_exe',cfg.bellhop_exe,'case_root',fullfile(cfg.output_dir,tag), ...
                'frequency_hz',f,'c0_mps',cfg.c0_mps,'source_depth_m',0, ...
                'receiver_depths_m',0,'receiver_ranges_m',cfg.ranges_m,'run_type','C', ...
                'beam_count',beams,'angle_limits_deg',cfg.angle_limits_deg, ...
                'step_m',step,'domain_half_depth_m',cfg.domain_half_depth_m);
            run=run_bellhop_freefield_vertical(c); p=run.data.pressure(1,:).';
            k=2*pi*f/cfg.c0_mps; C=abs(p).*cfg.ranges_m(:);
            Qplus=p.*cfg.ranges_m(:).*exp(-1i*k*cfg.ranges_m(:));
            Qminus=p.*cfg.ranges_m(:).*exp(1i*k*cfg.ranges_m(:));
            phase_rms_plus=local_constant_phase_rms(Qplus);
            phase_rms_minus=local_constant_phase_rms(Qminus);
            row=local_row(); row.frequency_hz=f; row.beam_count=beams; row.step_m=step;
            row.C_mean=mean(C); row.C_std=std(C); row.C_rel_span=(max(C)-min(C))/mean(C);
            row.Q_plus_phase_rms_rad=phase_rms_plus;
            row.Q_minus_phase_rms_rad=phase_rms_minus;
            row.selected_spatial_sign=double(phase_rms_plus<=phase_rms_minus)-double(phase_rms_minus<phase_rms_plus);
            if row.selected_spatial_sign>0, Q=Qplus; else, Q=Qminus; end
            row.source_constant=mean(Q./abs(Q));
            row.one_over_r_slope=polyfit(log10(cfg.ranges_m(:)),log10(abs(p)),1); row.one_over_r_slope=row.one_over_r_slope(1);
            row.normalization_is_one=abs(row.C_mean-1)<=cfg.normalization_mean_tolerance && row.C_rel_span<=cfg.normalization_span_tolerance;
            rows(q)=row; runs{q}=run;
        end
    end
end
t=struct2table(rows);

% Arrival audit is separate because coherent and arrivals exercise distinct writers/readers.
c=runs{end}.config; c.case_root=fullfile(cfg.output_dir,'arrival_audit'); c.run_type='A';
arrival_run=run_bellhop_freefield_vertical(c); a=arrival_run.data;
direct=a(a.top_bounce_count==0 & a.bottom_bounce_count==0,:);
arrival_rows=repmat(struct('range_m',NaN,'arrival_count',0,'delay_s',NaN, ...
    'exact_delay_s',NaN,'delay_error_s',NaN),numel(cfg.ranges_m),1);
for rr=1:numel(cfg.ranges_m)
    target=cfg.ranges_m(rr); z=direct(abs(direct.receiver_range_m-target)<1e-6,:);
    arrival_rows(rr).range_m=target; arrival_rows(rr).arrival_count=height(z);
    if ~isempty(z)
        [~,ix]=min(abs(z.delay_s-target/cfg.c0_mps)); arrival_rows(rr).delay_s=z.delay_s(ix);
    end
    arrival_rows(rr).exact_delay_s=target/cfg.c0_mps;
    arrival_rows(rr).delay_error_s=arrival_rows(rr).delay_s-arrival_rows(rr).exact_delay_s;
end
arrival_table=struct2table(arrival_rows);

high=t(t.beam_count==cfg.beam_counts(end) & t.step_m==cfg.step_values_m(end),:);
normalization_passed=all(high.normalization_is_one);
if numel(unique(high.selected_spatial_sign))~=1
    error('Bellhop spatial phase sign is not stable across frequency.');
end
selected_sign=high.selected_spatial_sign(1);
source_constant=mean(high.source_constant./abs(high.source_constant));
source_constant=source_constant/abs(source_constant);
[beam_tl,beam_phase]=local_convergence_metrics(t,runs,cfg.beam_counts(end-1), ...
    cfg.beam_counts(end),cfg.step_values_m(end),cfg.step_values_m(end));
[step_tl,step_phase]=local_convergence_metrics(t,runs,cfg.beam_counts(end), ...
    cfg.beam_counts(end),cfg.step_values_m(end-1),cfg.step_values_m(end));
checks=table(["normalization_mean";"normalization_range_invariance";"one_over_r_slope"; ...
    "spatial_phase_convention";"source_phase_constant";"beam_tl_convergence"; ...
    "beam_phase_convergence";"step_tl_convergence";"step_phase_convergence"; ...
    "arrival_time";"zero_bounce_direct_exists"], ...
    [max(abs(high.C_mean-1));max(high.C_rel_span);max(abs(high.one_over_r_slope+1)); ...
    max(min(high.Q_plus_phase_rms_rad,high.Q_minus_phase_rms_rad)); ...
    max(abs(angle(high.source_constant.*conj(source_constant))));beam_tl;beam_phase;step_tl;step_phase; ...
    max(abs(arrival_table.delay_error_s));double(any(arrival_table.arrival_count<1))], ...
    [cfg.normalization_mean_tolerance;cfg.normalization_span_tolerance;cfg.slope_tolerance; ...
    cfg.phase_convention_tolerance_rad;cfg.source_phase_tolerance_rad; ...
    cfg.convergence_tl_tolerance_db;cfg.convergence_phase_tolerance_rad; ...
    cfg.convergence_tl_tolerance_db;cfg.convergence_phase_tolerance_rad; ...
    cfg.delay_tolerance_s;0],'VariableNames',{'check_name','value','limit'});
checks.passed=checks.value<=checks.limit;
validation=struct('schema_version','1.0.0','config',cfg,'normalization_table',t, ...
    'arrival_table',arrival_table,'runs',{runs},'arrival_run',arrival_run, ...
    'normalization_constant',1,'selected_spatial_sign',selected_sign, ...
    'bellhop_source_constant',source_constant,'green_conversion_factor',1/(4*pi), ...
    'conversion_formula','if sign=-1: conj(p_bh)/conj(q0)/(4*pi); if sign=+1: p_bh/q0/(4*pi)', ...
    'green_conversion_tl_db',20*log10(4*pi),'normalization_passed',normalization_passed, ...
    'checks',checks,'passed',normalization_passed && all(checks.passed));
validation.files=local_outputs(validation);
if ~validation.passed
    error('validate_bellhop_freefield_normalization_vertical:Failed', ...
        'Bellhop normalization audit failed; PE-Bellhop absolute comparison is blocked.');
end
end

function cfg=local_defaults(root)
cfg=struct('bellhop_exe',getenv('BELLHOP_EXE'),'output_dir',fullfile(root,'results','validation','pe_bellhop_freefield','bellhop_normalization'), ...
    'c0_mps',1500,'frequencies_hz',[3000 4000 5000],'ranges_m',[20 40 70 100], ...
    'beam_counts',[0 501 2001 5001 10001],'step_values_m',[0 0.1 0.05], ...
    'angle_limits_deg',[-180 180],'domain_half_depth_m',1000, ...
    'normalization_mean_tolerance',0.02,'normalization_span_tolerance',0.02, ...
    'slope_tolerance',0.02,'phase_convention_tolerance_rad',0.02, ...
    'source_phase_tolerance_rad',0.02,'delay_tolerance_s',5e-6, ...
    'convergence_tl_tolerance_db',0.05,'convergence_phase_tolerance_rad',0.02);
end

function cfg=local_overrides(cfg,o)
names=fieldnames(o); for ii=1:numel(names), if ~isfield(cfg,names{ii}), error('Unknown override: %s',names{ii}); end, cfg.(names{ii})=o.(names{ii}); end
end

function local_validate(cfg)
if isempty(cfg.bellhop_exe) || exist(cfg.bellhop_exe,'file')~=2, error('Set overrides.bellhop_exe or BELLHOP_EXE.'); end
if cfg.beam_counts(end)<2 || any(cfg.ranges_m<=0), error('Invalid Bellhop audit grid.'); end
end

function r=local_row()
r=struct('frequency_hz',NaN,'beam_count',NaN,'step_m',NaN,'C_mean',NaN, ...
    'C_std',NaN,'C_rel_span',NaN,'Q_plus_phase_rms_rad',NaN, ...
    'Q_minus_phase_rms_rad',NaN,'selected_spatial_sign',NaN, ...
    'source_constant',complex(NaN),'one_over_r_slope',NaN,'normalization_is_one',false);
end

function value=local_constant_phase_rms(q)
q0=mean(q./abs(q)); q0=q0/abs(q0); value=sqrt(mean(angle(q.*conj(q0)).^2));
end

function [max_tl_db,phase_rms_rad]=local_convergence_metrics(t,runs,beam_a,beam_b,step_a,step_b)
tl=[]; phase=[];
for f=unique(t.frequency_hz).'
    ia=find(t.frequency_hz==f & t.beam_count==beam_a & t.step_m==step_a,1);
    ib=find(t.frequency_hz==f & t.beam_count==beam_b & t.step_m==step_b,1);
    pa=runs{ia}.data.pressure(:); pb=runs{ib}.data.pressure(:);
    tl=[tl;20*log10(abs(pa)./abs(pb))]; %#ok<AGROW>
    phase=[phase;angle(pa.*conj(pb))]; %#ok<AGROW>
end
max_tl_db=max(abs(tl)); phase_rms_rad=sqrt(mean(phase.^2));
end

function files=local_outputs(v)
out=v.config.output_dir; n_csv=fullfile(out,'bellhop_normalization.csv');
a_csv=fullfile(out,'bellhop_arrival_audit.csv'); c_csv=fullfile(out,'bellhop_normalization_checks.csv');
mat_file=fullfile(out,'bellhop_freefield_normalization.mat'); fig_file=fullfile(out,'bellhop_normalization.png');
writetable(v.normalization_table,n_csv); writetable(v.arrival_table,a_csv); writetable(v.checks,c_csv);
fig=figure('Visible','off','Color','w','Position',[100 100 1320 420]); cleanup=onCleanup(@()close(fig));
subplot(1,3,1); hold on; high=v.normalization_table(v.normalization_table.beam_count==v.config.beam_counts(end),:);
for f=unique(high.frequency_hz).', q=high(high.frequency_hz==f,:); plot(q.step_m,q.C_mean,'o-','DisplayName',sprintf('%g Hz',f)); end
grid on; xlabel('Bellhop step (m)'); ylabel('mean |p|R'); legend('Location','best'); title('Normalization audit');
subplot(1,3,2); hold on; fine=v.normalization_table(v.normalization_table.step_m==v.config.step_values_m(end),:);
for f=unique(fine.frequency_hz).', q=fine(fine.frequency_hz==f,:); plot(q.beam_count,q.C_mean,'o-','DisplayName',sprintf('%g Hz',f)); end
grid on; xlabel('beam count (0=auto)'); ylabel('mean |p|R'); title('Beam convergence');
subplot(1,3,3); plot(v.arrival_table.range_m,1e6*v.arrival_table.delay_error_s,'o-'); grid on; xlabel('R (m)'); ylabel('delay error (us)'); title('Zero-bounce arrival');
exportgraphics(fig,fig_file,'Resolution',180); clear cleanup
validation=v; schema_version=v.schema_version; save(mat_file,'validation','schema_version','-v7.3');
files=struct('mat',mat_file,'normalization_csv',n_csv,'arrival_csv',a_csv,'checks_csv',c_csv,'figure',fig_file);
end
