function validation = validate_random_surface_window_robustness_vertical(overrides)
%VALIDATE_RANDOM_SURFACE_WINDOW_ROBUSTNESS_VERTICAL Random-surface window audit.
% The formal run uses three Hs levels, five seeds, and 160/192/256 m windows.

if nargin<1, overrides=struct(); end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
run(fullfile(root,'scripts','bootstrap_project.m'));
cfg=local_config(root,overrides);
if ~exist(cfg.output_dir,'dir'), mkdir(cfg.output_dir); end
if ~exist(cfg.checkpoint_dir,'dir'), mkdir(cfg.checkpoint_dir); end

cases=repmat(local_empty_case(),0,1);
master_rows=repmat(local_empty_master_row(),0,1);
for ss=1:numel(cfg.sea_states)
    sea=cfg.sea_states(ss);
    for kk=1:numel(cfg.seeds)
        seed=cfg.seeds(kk);
        [master_eta,master_meta,master_path]=local_master_surface(cfg,sea,seed);
        master_rows(end+1,1)=local_master_row(sea,seed,master_meta,master_path); %#ok<AGROW>
        order=[cfg.reference_window_m,setdiff(cfg.windows_m,cfg.reference_window_m,'stable')];
        for ww=1:numel(order)
            nominal=order(ww); n=local_n(nominal,cfg.dx_m); actual=n*cfg.dx_m;
            eta=local_crop(master_eta,n);
            checkpoint=local_case_path(cfg,sea.name,seed,nominal);
            fprintf('\nRandom-surface case: %s, Hs=%.3f m, seed=%d, W=%.6f m, N=%d.\n', ...
                sea.name,sea.hs_m,seed,actual,n);
            if cfg.resume_cases && exist(checkpoint,'file')==2
                d=load(checkpoint,'case_result'); case_result=d.case_result;
            else
                out=vertical_channel_model(local_params(cfg,sea,seed,actual,n,eta,true,cfg.frequency_hz));
                case_result=local_extract_case(out,sea,seed,nominal,eta,cfg);
                save(checkpoint,'case_result','-v7.3');
            end
            cases(end+1,1)=case_result; %#ok<AGROW>
            clear out eta case_result
        end
        clear master_eta
    end
end

comparison_rows=repmat(local_empty_comparison(),0,1);
for ss=1:numel(cfg.sea_states)
    sea=cfg.sea_states(ss);
    for seed=cfg.seeds
        idx=find(strcmp({cases.sea_state},sea.name) & [cases.seed]==seed);
        ref_idx=idx([cases(idx).window_nominal_m]==cfg.reference_window_m);
        if numel(ref_idx)~=1, error('Missing unique large-window reference.'); end
        for jj=1:numel(idx)
            comparison_rows(end+1,1)=local_compare(cases(idx(jj)),cases(ref_idx),cfg); %#ok<AGROW>
        end
    end
end

case_table=struct2table(local_case_rows(cases));
comparison_table=struct2table(comparison_rows);
master_table=struct2table(master_rows);
energy_table=local_energy_table(cases);
pass_rate_table=local_pass_rates(comparison_table,cfg);
selection_table=local_select_wideband(comparison_table,cfg);
writetable(master_table,fullfile(cfg.output_dir,'master_surface_summary.csv'));
writetable(case_table,fullfile(cfg.output_dir,'single_frequency_cases.csv'));
writetable(comparison_table,fullfile(cfg.output_dir,'single_frequency_reference_comparison.csv'));
writetable(energy_table,fullfile(cfg.output_dir,'single_frequency_stage_energy.csv'));
writetable(pass_rate_table,fullfile(cfg.output_dir,'single_frequency_pass_rates.csv'));
writetable(selection_table,fullfile(cfg.output_dir,'wideband_case_selection.csv'));

reference_rows=comparison_table(comparison_table.window_nominal_m==cfg.reference_window_m,:);
all_references_valid=all(reference_rows.edge_valid);
validation=struct('schema_version','1.0.0','config',cfg,'cases',cases, ...
    'master_table',master_table,'case_table',case_table, ...
    'comparison_table',comparison_table,'energy_table',energy_table, ...
    'pass_rate_table',pass_rate_table,'selection_table',selection_table, ...
    'all_references_valid',all_references_valid, ...
    'recommendation',local_preliminary_recommendation(comparison_table,cfg));
save(fullfile(cfg.output_dir,'random_surface_window_robustness_4khz.mat'),'validation','-v7.3');
generate_random_surface_window_robustness_figures(validation);
local_write_report(validation);
if ~all_references_valid && cfg.fail_on_invalid_reference
    error('At least one 256 m reference failed the edge-validity gates.');
end
end

function cfg=local_config(root,o)
sea_states=struct('name',{'weak','current','strong'}, ...
    'wind_speed_mps',{5,5,5},'hs_m',{0.05,0.5,1.0});
cfg=struct('output_dir',fullfile(root,'results','validation','pe_random_surface_window_robustness'), ...
    'checkpoint_dir',fullfile(root,'results','validation','pe_random_surface_window_robustness','checkpoints_4khz'), ...
    'report_path',fullfile(root,'reports','pe_random_surface_window_robustness_report.md'), ...
    'frequency_hz',4000,'c0_mps',1500,'z_tx_m',100,'z_rx_m',3, ...
    'sigma_src_m',0.3,'stepz_lamb',0.5,'dx_m',50/256, ...
    'windows_m',[160 192 256],'reference_window_m',256,'seeds',12345:12349, ...
    'sea_states',sea_states,'reflect_coeff',-1,'max_z_samples',33, ...
    'center_radius_m',2,'tl_limit_db',0.1,'phase_limit_rad',0.02, ...
    'center_l2_limit',0.01,'boundary_limit_db',-40,'edge5_fraction_limit',1e-4, ...
    'resume_cases',true,'use_gpu',true,'fail_on_invalid_reference',false);
names=fieldnames(o); for ii=1:numel(names), cfg.(names{ii})=o.(names{ii}); end
if isfield(o,'output_dir') && ~isfield(o,'checkpoint_dir')
    cfg.checkpoint_dir=fullfile(cfg.output_dir,'checkpoints_4khz');
end
if isfield(o,'output_dir') && ~isfield(o,'report_path')
    cfg.report_path=fullfile(cfg.output_dir,'reduced_random_surface_report.md');
end
cfg.windows_m=unique(cfg.windows_m,'stable');
if ~ismember(cfg.reference_window_m,cfg.windows_m)
    error('reference_window_m must be included in windows_m.');
end
end

function [eta,meta,path]=local_master_surface(cfg,sea,seed)
safe=matlab.lang.makeValidName(sea.name);
path=fullfile(cfg.output_dir,sprintf('master_surface_%s_seed%d.mat',safe,seed));
if cfg.resume_cases && exist(path,'file')==2
    d=load(path,'master_eta','master_meta'); eta=d.master_eta; meta=d.master_meta; return
end
n=local_n(cfg.reference_window_m,cfg.dx_m); w=n*cfg.dx_m;
spec=raw_pm_spectrum_grid_vertical(sea.wind_speed_mps,n,n,w,w);
[raw,sample_meta]=sample_raw_pm_surface_vertical(spec,seed);
raw_hs=4*std(raw(:));
if raw_hs>0, eta=raw*(sea.hs_m/raw_hs); else, eta=zeros(size(raw)); end
meta=struct('sea_state',sea.name,'wind_speed_mps',sea.wind_speed_mps, ...
    'target_hs_m',sea.hs_m,'seed',seed,'nominal_width_m',cfg.reference_window_m, ...
    'actual_width_m',w,'nx',n,'ny',n,'dx_m',cfg.dx_m,'raw_hs_m',raw_hs, ...
    'scaled_hs_m',4*std(eta(:)),'mean_m',mean(eta(:)),'std_m',std(eta(:)), ...
    'min_m',min(eta(:)),'max_m',max(eta(:)),'sum_m',sum(eta(:)), ...
    'sum_squares_m2',sum(eta(:).^2),'normalization_count',1, ...
    'crop_rule','exact central crop; no crop is regenerated or renormalized', ...
    'generation_meta',sample_meta);
master_eta=eta; master_meta=meta;
save(path,'master_eta','master_meta','-v7.3');
end

function p=local_params(cfg,sea,seed,w,n,eta,diagnostics,freq)
p=struct('f0',freq,'f_ref_hz',cfg.frequency_hz,'c0',cfg.c0_mps, ...
    'z_max',cfg.z_tx_m,'z_tx',cfg.z_tx_m,'z_rx',cfg.z_rx_m, ...
    'x_tx',0,'y_tx',0,'x_rx',0,'y_rx',0,'xw',w,'yw',w,'nx',n,'ny',n, ...
    'source_mode','gaussian','sigma_src_m',cfg.sigma_src_m,'stepz_lamb',cfg.stepz_lamb, ...
    'env_mode','uniform','show_figures',false,'save_mode','rx_only', ...
    'enforce_1_over_R',false,'enable_surface_reflection',true, ...
    'surface_boundary_model','kirchhoff_spatial','surface_reflect_coeff',cfg.reflect_coeff, ...
    'surface_phase_mode','normal','sea_wind_speed',sea.wind_speed_mps, ...
    'sea_hs_target',sea.hs_m,'surface_roughness_scale_mode','target_hs', ...
    'sea_seed',seed,'surface_elevation_override_xy',eta, ...
    'surface_wavefield_diagnostics',diagnostics, ...
    'surface_wavefield_max_z_samples',cfg.max_z_samples,'surface_wavefield_slice_axis','x', ...
    'validation_allow_extended_window',w>100, ...
    'validation_allow_extended_sponge_ratio',true,'sponge_ratio',0, ...
    'alpha_max_np_per_m',0,'enable_bubbles',false,'use_gpu',cfg.use_gpu);
end

function c=local_extract_case(out,sea,seed,nominal,eta,cfg)
m=out.surface_wavefield_meta;
phase_error=norm(abs(m.surface_reflected_xy(:))-abs(cfg.reflect_coeff)* ...
    abs(m.surface_incident_xy(:)))/max(norm(abs(m.surface_incident_xy(:))),eps);
fields={'surface_incident_xy','surface_reflected_xy','receiver_reflected_xy'};
[X,Y]=meshgrid(m.x_m,m.y_m); mask=hypot(X,Y)<=cfg.center_radius_m;
c=local_empty_case(); c.sea_state=sea.name; c.wind_speed_mps=sea.wind_speed_mps;
c.target_hs_m=sea.hs_m; c.seed=seed; c.window_nominal_m=nominal;
c.window_actual_m=out.config.xw; c.nx=out.config.nx; c.dx_m=out.config.dx;
c.crop_hs_m=4*std(eta(:)); c.crop_mean_m=mean(eta(:));
c.H_reflect_physical=out.H_reflect_physical_f(out.idx_f_ref);
c.H_total_physical=out.H_physical_f(out.idx_f_ref);
c.H_direct=out.H_direct_f(out.idx_f_ref); c.H_reflect=out.H_reflect_f(out.idx_f_ref);
c.H_total=out.H_f(out.idx_f_ref); c.closure_error=abs(c.H_total-c.H_direct-c.H_reflect);
c.phase_screen_magnitude_rel_l2=phase_error;
c.endpoint_error=max(struct2array(m.endpoint_consistency));
c.center_x_m=X(mask); c.center_y_m=Y(mask);
for ii=1:numel(fields), c.(erase(fields{ii},'_xy'))=m.(fields{ii})(mask); end
c.incident_energy_trace=m.incident_energy_trace; c.reflected_energy_trace=m.reflected_energy_trace;
c=local_stage_metrics(c);
end

function c=local_stage_metrics(c)
up=c.incident_energy_trace; down=c.reflected_energy_trace;
c.upward_max_edge5_fraction=max(up.edge5_fraction);
c.upward_max_boundary_db=max(up.boundary_max_relative_db);
c.surface_incident_edge5_fraction=up.edge5_fraction(end);
c.surface_incident_boundary_db=up.boundary_max_relative_db(end);
c.surface_reflected_edge5_fraction=down.edge5_fraction(1);
c.surface_reflected_boundary_db=down.boundary_max_relative_db(1);
c.downward_max_edge5_fraction=max(down.edge5_fraction);
c.downward_max_boundary_db=max(down.boundary_max_relative_db);
c.receiver_edge5_fraction=down.edge5_fraction(end);
c.receiver_boundary_db=down.boundary_max_relative_db(end);
c.max_edge5_fraction=max(c.upward_max_edge5_fraction,c.downward_max_edge5_fraction);
c.max_boundary_db=max(c.upward_max_boundary_db,c.downward_max_boundary_db);
end

function r=local_compare(c,ref,cfg)
r=local_empty_comparison();
copy={'sea_state','wind_speed_mps','target_hs_m','seed','window_nominal_m', ...
    'window_actual_m','crop_hs_m','max_edge5_fraction','max_boundary_db', ...
    'upward_max_edge5_fraction','surface_incident_edge5_fraction', ...
    'surface_reflected_edge5_fraction','downward_max_edge5_fraction', ...
    'receiver_edge5_fraction','closure_error','phase_screen_magnitude_rel_l2','endpoint_error'};
for ii=1:numel(copy), r.(copy{ii})=c.(copy{ii}); end
r.reference_window_actual_m=ref.window_actual_m;
r.reference_edge_valid=ref.max_edge5_fraction<=cfg.edge5_fraction_limit && ...
    ref.max_boundary_db<=cfg.boundary_limit_db;
r.reflect_tl_db=-20*log10(abs(c.H_reflect_physical));
r.reflect_phase_rad=angle(c.H_reflect_physical);
r.reflect_tl_error_db=-20*log10(abs(c.H_reflect_physical/ref.H_reflect_physical));
r.reflect_phase_error_rad=angle(c.H_reflect_physical*conj(ref.H_reflect_physical));
assert(isequal(c.center_x_m,ref.center_x_m)&&isequal(c.center_y_m,ref.center_y_m), ...
    'Center samples do not align on the common-dx grids.');
for name={'surface_incident','surface_reflected','receiver_reflected'}
    key=name{1}; r.([key '_center_l2'])=norm(c.(key)-ref.(key))/max(norm(ref.(key)),eps);
end
r.response_valid=abs(r.reflect_tl_error_db)<=cfg.tl_limit_db && ...
    abs(r.reflect_phase_error_rad)<=cfg.phase_limit_rad;
r.center_field_valid=r.surface_incident_center_l2<=cfg.center_l2_limit && ...
    r.surface_reflected_center_l2<=cfg.center_l2_limit && ...
    r.receiver_reflected_center_l2<=cfg.center_l2_limit;
r.edge_valid=r.max_edge5_fraction<=cfg.edge5_fraction_limit && ...
    r.max_boundary_db<=cfg.boundary_limit_db;
r.passed=r.reference_edge_valid && r.response_valid && r.center_field_valid && r.edge_valid;
r.severity=max([abs(r.reflect_tl_error_db)/cfg.tl_limit_db, ...
    abs(r.reflect_phase_error_rad)/cfg.phase_limit_rad, ...
    r.surface_incident_center_l2/cfg.center_l2_limit, ...
    r.surface_reflected_center_l2/cfg.center_l2_limit, ...
    r.receiver_reflected_center_l2/cfg.center_l2_limit, ...
    r.max_edge5_fraction/cfg.edge5_fraction_limit,10^(r.max_boundary_db/20)/10^(cfg.boundary_limit_db/20)]);
end

function t=local_pass_rates(rows,cfg)
out=repmat(struct('sea_state',"",'target_hs_m',NaN,'window_nominal_m',NaN, ...
    'passed_count',0,'total_count',0,'pass_rate',NaN,'worst_tl_error_db',NaN, ...
    'worst_phase_error_rad',NaN,'worst_receiver_center_l2',NaN, ...
    'worst_edge5_fraction',NaN,'worst_boundary_db',NaN,'all_references_valid',false),0,1);
for ss=1:numel(cfg.sea_states)
    for w=cfg.windows_m
        q=rows(strcmp(rows.sea_state,cfg.sea_states(ss).name) & rows.window_nominal_m==w,:);
        r=local_empty_rate(); r.sea_state=string(cfg.sea_states(ss).name); r.target_hs_m=cfg.sea_states(ss).hs_m;
        r.window_nominal_m=w; r.passed_count=sum(q.passed); r.total_count=height(q); r.pass_rate=r.passed_count/r.total_count;
        r.worst_tl_error_db=max(abs(q.reflect_tl_error_db)); r.worst_phase_error_rad=max(abs(q.reflect_phase_error_rad));
        r.worst_receiver_center_l2=max(q.receiver_reflected_center_l2); r.worst_edge5_fraction=max(q.max_edge5_fraction);
        r.worst_boundary_db=max(q.max_boundary_db); r.all_references_valid=all(q.reference_edge_valid);
        out(end+1,1)=r; %#ok<AGROW>
    end
end
t=struct2table(out);
end

function t=local_select_wideband(rows,cfg)
selected=repmat(struct('sea_state',"",'target_hs_m',NaN,'wind_speed_mps',NaN, ...
    'seed',NaN,'reason',""),0,1);
nonref=sort(setdiff(cfg.windows_m,cfg.reference_window_m));
if isempty(nonref), error('At least one non-reference window is required.'); end
if ismember(192,nonref), primary=192; else, primary=max(nonref); end
if ismember(160,nonref), lowcost=160; else, lowcost=min(nonref); end
q192=rows(rows.window_nominal_m==primary,:); q160=rows(rows.window_nominal_m==lowcost,:);
for ss=1:numel(cfg.sea_states)
    q=q192(strcmp(q192.sea_state,cfg.sea_states(ss).name),:);
    [~,ii]=max(q.severity); selected=local_add_selection(selected,q(ii,:),'sea-state worst 192 m');
end
[~,ii]=max(q192.max_edge5_fraction); selected=local_add_selection(selected,q192(ii,:),'global worst 192 m edge');
[~,ii]=max(q160.severity); selected=local_add_selection(selected,q160(ii,:),'global worst 160 m reference error');
t=struct2table(selected);
end

function a=local_add_selection(a,row,reason)
if isempty(a)
    idx=[];
else
    idx=find(strcmp([a.sea_state],string(row.sea_state)) & [a.seed]==row.seed,1);
end
if isempty(idx)
    r=struct('sea_state',string(row.sea_state),'target_hs_m',row.target_hs_m, ...
        'wind_speed_mps',row.wind_speed_mps,'seed',row.seed,'reason',string(reason));
    a(end+1,1)=r;
else
    a(idx).reason=a(idx).reason+"; "+string(reason);
end
end

function r=local_preliminary_recommendation(rows,~)
q192=rows(rows.window_nominal_m==192,:); q160=rows(rows.window_nominal_m==160,:);
if all(q192.passed)
    choice='A'; text='192 m + no sponge passed every 4 kHz realization.';
elseif all(q192.response_valid & q192.center_field_valid) && any(~q192.edge_valid)
    choice='D'; text='192 m receiver response passed, but at least one strict edge gate failed.';
else
    choice='D'; text='192 m failed at least one response/field gate relative to 256 m.';
end
r=struct('choice',choice,'text',text,'pass_192',sum(q192.passed), ...
    'total_192',height(q192),'pass_160',sum(q160.passed),'total_160',height(q160));
end

function local_write_report(v)
path=v.config.report_path;
fid=fopen(path,'w'); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'# Random-surface reflected-chain window robustness\n\n');
fprintf(fid,'Status after 4 kHz stage: **%s**. PE marching and production defaults were not changed.\n\n',v.recommendation.choice);
fprintf(fid,'## Environment\n\n- Uniform c=1500 m/s; Tx=(0,0,100 m); Rx=(0,0,3 m); production Gaussian sigma=0.3 m.\n');
fprintf(fid,'- Kirchhoff spatial, U=5 m/s, Hs=[0.05,0.5,1.0] m, seeds 12345--12349.\n');
fprintf(fid,'- Fixed dx=50/256 m; 160/192 m are exact central crops of each 256 m master and are never renormalized; sponge is off.\n\n');
fprintf(fid,'## 4 kHz pass rates\n\n');
fprintf(fid,'| sea | Hs (m) | nominal W (m) | pass | rate | worst TL (dB) | worst phase (rad) | worst receiver L2 | worst edge5 | worst boundary (dB) |\n');
fprintf(fid,'|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for ii=1:height(v.pass_rate_table)
    r=v.pass_rate_table(ii,:);
    fprintf(fid,'| %s | %.3f | %.0f | %d/%d | %.3f | %.6g | %.6g | %.6g | %.6g | %.6g |\n', ...
        string(r.sea_state),r.target_hs_m,r.window_nominal_m,r.passed_count,r.total_count, ...
        r.pass_rate,r.worst_tl_error_db,r.worst_phase_error_rad, ...
        r.worst_receiver_center_l2,r.worst_edge5_fraction,r.worst_boundary_db);
end
fprintf(fid,'\nPreliminary recommendation: **%s** -- %s\n',v.recommendation.choice,v.recommendation.text);
fprintf(fid,'\nThe final Q1--Q6 decision is appended after the selected 3--5 kHz qualification.\n');
end

function path=local_case_path(cfg,sea,seed,w)
path=fullfile(cfg.checkpoint_dir,sprintf('%s_seed%d_W%d.mat',matlab.lang.makeValidName(sea),seed,round(w)));
end
function n=local_n(w,dx), n=2*round(w/dx/2); end
function c=local_crop(a,n), s=(size(a,1)-n)/2+1; c=a(s:s+n-1,s:s+n-1); end
function r=local_master_row(sea,seed,m,path)
r=local_empty_master_row(); r.sea_state=string(sea.name); r.wind_speed_mps=sea.wind_speed_mps;
r.target_hs_m=sea.hs_m; r.seed=seed; r.actual_width_m=m.actual_width_m; r.nx=m.nx;
r.raw_hs_m=m.raw_hs_m; r.scaled_hs_m=m.scaled_hs_m; r.mean_m=m.mean_m;
r.std_m=m.std_m; r.sum_m=m.sum_m; r.sum_squares_m2=m.sum_squares_m2; r.file=string(path);
end
function rows=local_case_rows(cases)
template=local_empty_case_row(); rows=repmat(template,numel(cases),1);
names=fieldnames(template); for ii=1:numel(cases), for jj=1:numel(names), rows(ii).(names{jj})=cases(ii).(names{jj}); end, end
end
function t=local_energy_table(cases)
r=repmat(struct('sea_state',"",'target_hs_m',NaN,'seed',NaN,'window_nominal_m',NaN, ...
    'stage',"",'z_m',NaN,'total_energy',NaN,'center_energy_rho_le_2m',NaN, ...
    'edge5_energy',NaN,'edge10_energy',NaN,'edge5_fraction',NaN, ...
    'edge10_fraction',NaN,'boundary_max_relative_db',NaN),0,1);
for ii=1:numel(cases)
    for stage={'upward','downward'}
        if strcmp(stage{1},'upward'), tr=cases(ii).incident_energy_trace; else, tr=cases(ii).reflected_energy_trace; end
        for jj=1:numel(tr.z_m)
            q=numel(r)+1; r(q).sea_state=string(cases(ii).sea_state); r(q).target_hs_m=cases(ii).target_hs_m;
            r(q).seed=cases(ii).seed; r(q).window_nominal_m=cases(ii).window_nominal_m;
            r(q).stage=string(stage{1}); r(q).z_m=tr.z_m(jj);
            for f={'total_energy','center_energy_rho_le_2m','edge5_energy','edge10_energy','edge5_fraction','edge10_fraction','boundary_max_relative_db'}
                r(q).(f{1})=tr.(f{1})(jj);
            end
        end
    end
end
t=struct2table(r);
end
function c=local_empty_case()
c=struct('sea_state','','wind_speed_mps',NaN,'target_hs_m',NaN,'seed',NaN, ...
    'window_nominal_m',NaN,'window_actual_m',NaN,'nx',NaN,'dx_m',NaN, ...
    'crop_hs_m',NaN,'crop_mean_m',NaN,'H_reflect_physical',complex(NaN), ...
    'H_total_physical',complex(NaN),'H_direct',complex(NaN),'H_reflect',complex(NaN), ...
    'H_total',complex(NaN),'closure_error',NaN,'phase_screen_magnitude_rel_l2',NaN, ...
    'endpoint_error',NaN,'center_x_m',[],'center_y_m',[],'surface_incident',[], ...
    'surface_reflected',[],'receiver_reflected',[],'incident_energy_trace',struct(), ...
    'reflected_energy_trace',struct(),'upward_max_edge5_fraction',NaN, ...
    'upward_max_boundary_db',NaN,'surface_incident_edge5_fraction',NaN, ...
    'surface_incident_boundary_db',NaN,'surface_reflected_edge5_fraction',NaN, ...
    'surface_reflected_boundary_db',NaN,'downward_max_edge5_fraction',NaN, ...
    'downward_max_boundary_db',NaN,'receiver_edge5_fraction',NaN, ...
    'receiver_boundary_db',NaN,'max_edge5_fraction',NaN,'max_boundary_db',NaN);
end
function r=local_empty_case_row()
r=rmfield(local_empty_case(),{'center_x_m','center_y_m','surface_incident','surface_reflected', ...
    'receiver_reflected','incident_energy_trace','reflected_energy_trace'});
end
function r=local_empty_comparison()
r=struct('sea_state',"",'wind_speed_mps',NaN,'target_hs_m',NaN,'seed',NaN, ...
    'window_nominal_m',NaN,'window_actual_m',NaN,'reference_window_actual_m',NaN, ...
    'crop_hs_m',NaN,'reflect_tl_db',NaN,'reflect_phase_rad',NaN, ...
    'reflect_tl_error_db',NaN,'reflect_phase_error_rad',NaN, ...
    'surface_incident_center_l2',NaN,'surface_reflected_center_l2',NaN, ...
    'receiver_reflected_center_l2',NaN,'max_edge5_fraction',NaN,'max_boundary_db',NaN, ...
    'upward_max_edge5_fraction',NaN,'surface_incident_edge5_fraction',NaN, ...
    'surface_reflected_edge5_fraction',NaN,'downward_max_edge5_fraction',NaN, ...
    'receiver_edge5_fraction',NaN,'closure_error',NaN, ...
    'phase_screen_magnitude_rel_l2',NaN,'endpoint_error',NaN, ...
    'reference_edge_valid',false,'response_valid',false,'center_field_valid',false, ...
    'edge_valid',false,'passed',false,'severity',NaN);
end
function r=local_empty_master_row()
r=struct('sea_state',"",'wind_speed_mps',NaN,'target_hs_m',NaN,'seed',NaN, ...
    'actual_width_m',NaN,'nx',NaN,'raw_hs_m',NaN,'scaled_hs_m',NaN, ...
    'mean_m',NaN,'std_m',NaN,'sum_m',NaN,'sum_squares_m2',NaN,'file',"");
end
function r=local_empty_rate()
r=struct('sea_state',"",'target_hs_m',NaN,'window_nominal_m',NaN, ...
    'passed_count',0,'total_count',0,'pass_rate',NaN,'worst_tl_error_db',NaN, ...
    'worst_phase_error_rad',NaN,'worst_receiver_center_l2',NaN, ...
    'worst_edge5_fraction',NaN,'worst_boundary_db',NaN,'all_references_valid',false);
end
