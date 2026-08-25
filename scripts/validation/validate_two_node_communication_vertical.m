%% U=5/U=8 physical, joint, and conditional-generator communication validation
% Set TWO_NODE_COMM_MODE=smoke or full. Public communication defaults stay unchanged.
clear; close all; clc;
root=fileparts(fileparts(fileparts(mfilename('fullpath')))); addpath(root); setup_vertical_project();
mode=lower(strtrim(getenv('TWO_NODE_COMM_MODE'))); if isempty(mode), mode='smoke'; end
if ~ismember(mode,{'smoke','full'}), error('TWO_NODE_COMM_MODE must be smoke or full.'); end
out_dir=fullfile(root,'results','validation','two_node_communication'); if ~exist(out_dir,'dir'), mkdir(out_dir); end
if strcmp(mode,'smoke')
    L=4; n_sym=1000; eb=[0,8,16]; bootstrap_count=200;
else
    L=32; n_sym=8000; eb=0:2:20; bootstrap_count=2000;
end
winds=[5,8]; source_names={'kdomain_pe','joint_cached_pe','stats_full','stats_99_9'};
comm_cfg=struct('M',4,'n_sym',n_sym,'EbN0_dB_list',eb,'symbol_rate_hz',1000, ...
    'tap_fft_size',2048,'f_ref_hz',6000,'tap_energy_ratio',0.999, ...
    'circular_peak_align',true, ...
    'mmse_reg_eps',1e-6,'bits_seed',3100001,'noise_seed_base',3200000, ...
    'noise_model','awgn','ebn0_reference','rx_clean', ...
    'bootstrap_count',bootstrap_count,'bootstrap_seed',3300001);
validation_run_meta=pe_phase_release_run_meta_vertical(root,struct( ...
    'frequency_axis_hz',linspace(4000,8000,64).', ...
    'seed_definition',struct('bits',comm_cfg.bits_seed, ...
    'noise_base',comm_cfg.noise_seed_base,'bootstrap',comm_cfg.bootstrap_seed)));

ensemble_file=fullfile(out_dir,['two_node_comm_ensembles_' mode '.mat']);
reuse_ensembles=false;
if exist(ensemble_file,'file')
    load(ensemble_file,'ensembles','generation');
    reuse_ensembles=local_ensembles_are_direct_dsp(ensembles,winds,source_names) && ...
        isfield(ensembles,'validation_run_meta') && ...
        strcmp(ensembles.validation_run_meta.run_id,validation_run_meta.run_id);
end
if ~reuse_ensembles
    ensembles=struct(); generation=struct();
    for ww=1:numel(winds)
        U=winds(ww); key=sprintf('u%d',U); fprintf('Generating fresh U=%d communication channels.\n',U);
        if U==5
            cache_file=fullfile(root,'results','validation','u5_conditional_channel_f64','f64_joint_pe_cache.mat');
            model_file=fullfile(root,'results','validation','u5_conditional_channel_f64','u5_conditional_channel_model_f64_full.mat');
        else
            cache_file=fullfile(root,'results','validation','u8_conditional_channel_f64','f64_streaming_joint_pe_cache.mat');
            model_file=fullfile(root,'results','validation','u8_conditional_channel_f64','u8_conditional_channel_model_f64_full.mat');
        end
        d=load(cache_file,'cache','pm_spec','joint_model'); m=load(model_file,'model');
        assert_pe_phase_release_artifact_vertical(m.model,validation_run_meta.run_id, ...
            sprintf('U=%d conditional model',U));
        geometry=struct('z_tx',d.cache.cfg.z_tx,'z_rx',d.cache.cfg.z_rx, ...
            'z_surface',0,'c0',d.cache.cfg.c0);
        [deterministic,phase_meta]=apply_pe_channel_phase_reference_vertical( ...
            d.cache.f_axis_hz(:),struct('direct_f',d.cache.H_direct_f, ...
            'reflect_coh_f',d.cache.H_ref_coh_f),geometry,'direct_dsp');
        m.model=upgrade_conditional_channel_phase_vertical(m.model,geometry);
        kd_seeds=3400000+U*10000+(1:L);
        [Hkd,tkd,map,ekd]=generate_cached_kdomain_ensemble_vertical(d.cache,d.pm_spec,kd_seeds,4);
        [Hj,tj,ej]=generate_cached_kstat_ensemble_vertical(d.cache,d.joint_model,L,3500000+U*10000,'joint',4);
        Hkd_total=deterministic.direct_f+deterministic.reflect_coh_f+Hkd;
        Hj_total=deterministic.direct_f+deterministic.reflect_coh_f+Hj;
        [sf,tf]=sample_conditional_channel_vertical(m.model,L,3600000+U*10000,struct('path','auto','rank','full'));
        [sl,tl]=sample_conditional_channel_vertical(m.model,L,3700000+U*10000,struct('path','auto','rank','99.9'));
        ensembles.(key)=struct( ...
            'kdomain_pe',local_set('kdomain_pe',U,d.cache.f_axis_hz,Hkd_total,phase_meta), ...
            'joint_cached_pe',local_set('joint_cached_pe',U,d.cache.f_axis_hz,Hj_total,phase_meta), ...
            'stats_full',local_set('stats_full',U,m.model.frequency_axis,sf.H_total_f,m.model.phase_reference_meta), ...
            'stats_99_9',local_set('stats_99_9',U,m.model.frequency_axis,sl.H_total_f,m.model.phase_reference_meta));
        ratio=map.crop_variance_m2./map.full_variance_m2;
        generation.(key)=struct('kdomain_timing',tkd,'joint_timing',tj, ...
            'stats_full_timing',tf,'stats_low_timing',tl,'mapping_mean',mean(ratio), ...
            'mapping_ci95',mean(ratio)+[-1,1]*1.96*std(ratio)/sqrt(L), ...
            'component_error',max(ekd,ej),'fresh_seed_families',true, ...
            'noise_in_channel',false);
        clear d m Hkd Hj Hkd_total Hj_total sf sl
    end
    ensembles.validation_run_meta=validation_run_meta;
    ensembles.schema_version='2.0.0';
    ensembles.phase_reference_meta=phase_meta;
    schema_version=ensembles.schema_version; phase_reference_meta=ensembles.phase_reference_meta;
    save(ensemble_file,'ensembles','generation','schema_version', ...
        'phase_reference_meta','validation_run_meta','-v7.3');
end

results=struct(); total_timer=tic;
for ww=1:numel(winds)
    U=winds(ww); key=sprintf('u%d',U);
    for ss=1:numel(source_names)
        name=source_names{ss}; fprintf('Communication U=%d source=%s.\n',U,name);
        cfg=comm_cfg; cfg.bootstrap_seed=comm_cfg.bootstrap_seed+1000*U+ss;
        results.(key).(name)=evaluate_mpsk_channel_ensemble_vertical(ensembles.(key).(name),cfg);
    end
end
elapsed=toc(total_timer);

% Higher-power paired compression test: full and low rank share latent z.
if strcmp(mode,'full'), L_rank_pair=128; else, L_rank_pair=8; end
rank_pair=struct();
for ww=1:numel(winds)
    U=winds(ww); key=sprintf('u%d',U);
    if U==5, mf=fullfile(root,'results','validation','u5_conditional_channel_f64','u5_conditional_channel_model_f64_full.mat');
    else, mf=fullfile(root,'results','validation','u8_conditional_channel_f64','u8_conditional_channel_model_f64_full.mat'); end
    md=load(mf,'model'); md.model=upgrade_conditional_channel_phase_vertical(md.model, ...
        struct('z_tx',100,'z_rx',3,'z_surface',0,'c0',1500));
    p=sample_conditional_channel_rank_pair_vertical(md.model,L_rank_pair,3800000+U*10000,'99.9');
    cfg=comm_cfg; cfg.bootstrap_seed=3900000+U*10000;
    rf=evaluate_mpsk_channel_ensemble_vertical(local_set('paired_full',U,md.model.frequency_axis,p.H_full_f,md.model.phase_reference_meta),cfg);
    rl=evaluate_mpsk_channel_ensemble_vertical(local_set('paired_99_9',U,md.model.frequency_axis,p.H_low_f,md.model.phase_reference_meta),cfg);
    rank_pair.(key)=struct('full',rf,'low',rl,'rank_full',p.rank_full,'rank_low',p.rank_low, ...
        'ber_difference_low_minus_full',rl.BER-rf.BER, ...
        'ser_difference_low_minus_full',rl.SER-rf.SER, ...
        'ber_difference_cluster95',local_paired_ci(rl.per_channel_BER-rf.per_channel_BER,bootstrap_count,4000000+U), ...
        'ser_difference_cluster95',local_paired_ci(rl.per_channel_SER-rf.per_channel_SER,bootstrap_count,4100000+U));
end

comparison=struct();
for ww=1:numel(winds)
    key=sprintf('u%d',winds(ww)); q=results.(key);
    comparison.(key)=struct( ...
        'joint_vs_kdomain',local_compare(q.joint_cached_pe,q.kdomain_pe), ...
        'stats_full_vs_kdomain',local_compare(q.stats_full,q.kdomain_pe), ...
        'stats_full_vs_joint',local_compare(q.stats_full,q.joint_cached_pe), ...
        'stats_low_vs_full',local_compare(q.stats_99_9,q.stats_full));
end
comparison.wind_trend=local_wind_trend(results.u5,results.u8);

conversion=struct();
for ww=1:numel(winds)
    key=sprintf('u%d',winds(ww));
    for ss=1:numel(source_names)
        name=source_names{ss}; r=results.(key).(name); H=ensembles.(key).(name).H_fm(:,1);
        cir=build_channel_cir_vertical(H,ensembles.(key).(name).f_axis_hz, ...
            struct('input_reference','direct_dsp'));
        [~,legacy]=build_communication_taps_vertical(struct('H_f',H, ...
            'f_axis_hz',ensembles.(key).(name).f_axis_hz, ...
            'phase_reference_meta',ensembles.(key).(name).phase_reference_meta), ...
            struct('symbol_rate_hz',comm_cfg.symbol_rate_hz,'n_fft',comm_cfg.tap_fft_size, ...
            'f_ref_hz',comm_cfg.f_ref_hz,'tap_energy_ratio',comm_cfg.tap_energy_ratio,'circular_peak_align',false));
        conversion.(key).(name)=struct( ...
            'max_ifft_closure_error',max(r.tap_meta.ifft_closure_error), ...
            'median_tap_count',median(r.tap_meta.tap_count), ...
            'max_prepeak_energy_fraction',max(r.tap_meta.prepeak_energy_fraction), ...
            'max_actual_discard_fraction',max(r.tap_meta.actual_discard_fraction), ...
            'min_energy_kept',min(r.tap_meta.energy_kept), ...
            'legacy_first_sample_tap_count',legacy.tap_count, ...
            'legacy_first_sample_prepeak_discard_fraction',legacy.discarded_pre_peak_energy_fraction_if_peak_sync, ...
            'physical_resolution_s',cir.physical_delay_resolution_s, ...
            'symbol_tap_spacing_s',1/comm_cfg.symbol_rate_hz, ...
            'maximum_unambiguous_delay_s',cir.maximum_unambiguous_delay_s);
    end
end

rows=local_rows(results,winds,source_names); summary_table=struct2table(rows);
writetable(summary_table,fullfile(out_dir,['two_node_comm_summary_' mode '.csv']));
validation=struct('mode',mode,'winds_mps',winds,'channel_count_per_source',L, ...
    'comm_cfg',comm_cfg,'generation',generation,'results',results,'comparison',comparison, ...
    'rank_pair',rank_pair,'conversion',conversion,'elapsed_s',elapsed,'summary_table',summary_table, ...
    'noise_in_channel',false,'public_surface_defaults_changed',false, ...
    'channel_phase_reference','direct_dsp');
validation.schema_version='2.0.0'; validation.phase_reference_meta=phase_meta;
validation.validation_run_meta=validation_run_meta;
schema_version=validation.schema_version; phase_reference_meta=validation.phase_reference_meta;
save(fullfile(out_dir,['two_node_communication_validation_' mode '.mat']),'validation', ...
    'schema_version','phase_reference_meta','validation_run_meta','-v7.3');
local_plots(validation,out_dir); local_text(validation,fullfile(out_dir,['summary_' mode '.txt']));
fprintf('Two-node communication %s complete in %.3f s.\n',mode,elapsed);

function s=local_set(name,U,f,H,phase_meta)
s=struct('name',name,'wind_speed_mps',U,'f_axis_hz',f(:),'H_fm',H, ...
    'phase_reference_meta',phase_meta);
end
function passed=local_ensembles_are_direct_dsp(ensembles,winds,names)
passed=true;
for ww=1:numel(winds)
    key=sprintf('u%d',winds(ww));
    if ~isfield(ensembles,key), passed=false; return; end
    for ss=1:numel(names)
        if ~isfield(ensembles.(key),names{ss}), passed=false; return; end
        q=ensembles.(key).(names{ss});
        if ~isfield(q,'phase_reference_meta') || ...
                ~strcmp(q.phase_reference_meta.target_reference,'direct_dsp')
            passed=false; return
        end
    end
end
end
function c=local_compare(a,b)
floorv=1/(2*a.channel_count*a.cfg.n_sym*log2(a.cfg.M));
c=struct('ber_log10_rmse',sqrt(mean((log10(max(a.BER,floorv))-log10(max(b.BER,floorv))).^2)), ...
    'ser_log10_rmse',sqrt(mean((log10(max(a.SER,floorv))-log10(max(b.SER,floorv))).^2)), ...
    'ber_max_abs_difference',max(abs(a.BER-b.BER)), ...
    'ser_max_abs_difference',max(abs(a.SER-b.SER)), ...
    'ber_cluster_ci_overlap_fraction',mean(local_overlap(a.BER_cluster95,b.BER_cluster95)), ...
    'ser_cluster_ci_overlap_fraction',mean(local_overlap(a.SER_cluster95,b.SER_cluster95)));
end
function tf=local_overlap(a,b), tf=max(a(1,:),b(1,:))<=min(a(2,:),b(2,:)); end
function ci=local_paired_ci(x,B,seed)
rng(seed,'twister'); [L,E]=size(x); m=zeros(B,E);
for bb=1:B, idx=randi(L,L,1); m(bb,:)=mean(x(idx,:),1); end
m=sort(m,1); lo=max(1,round(.025*(B-1)+1)); hi=min(B,round(.975*(B-1)+1)); ci=[m(lo,:);m(hi,:)];
end
function t=local_wind_trend(a,b)
names=fieldnames(a); t=struct();
for ii=1:numel(names)
    n=names{ii}; floorv=1/(2*a.(n).channel_count*a.(n).cfg.n_sym*log2(a.(n).cfg.M));
    t.(n)=struct('mean_log10_ber_u8_minus_u5',mean(log10(max(b.(n).BER,floorv))-log10(max(a.(n).BER,floorv))), ...
        'mean_ber_u5',mean(a.(n).BER),'mean_ber_u8',mean(b.(n).BER));
end
end
function rows=local_rows(r,winds,names)
cells=cell(numel(winds)*numel(names)*numel(r.u5.(names{1}).cfg.EbN0_dB_list),1); p=0;
for ww=1:numel(winds), key=sprintf('u%d',winds(ww));
    for ss=1:numel(names), q=r.(key).(names{ss});
        for ee=1:numel(q.cfg.EbN0_dB_list), p=p+1; cells{p}=struct('U_mps',winds(ww),'source',names{ss}, ...
            'EbN0_dB',q.cfg.EbN0_dB_list(ee),'BER',q.BER(ee),'BER_ci_low',q.BER_cluster95(1,ee), ...
            'BER_ci_high',q.BER_cluster95(2,ee),'SER',q.SER(ee),'SER_ci_low',q.SER_cluster95(1,ee), ...
            'SER_ci_high',q.SER_cluster95(2,ee),'BER_pre_eq',q.BER_pre_equalizer(ee), ...
            'SER_pre_eq',q.SER_pre_equalizer(ee)); end
    end
end
rows=vertcat(cells{:});
end
function local_plots(v,out)
colors=lines(4); marks={'o','s','^','d'};
for ww=1:numel(v.winds_mps), U=v.winds_mps(ww); key=sprintf('u%d',U); q=v.results.(key);
    fig=figure('Visible','off','Color','w'); tiledlayout(1,2,'Padding','compact');
    nexttile; hold on; nexttile; hold on;
    names={'kdomain_pe','joint_cached_pe','stats_full','stats_99_9'};
    for ss=1:numel(names), r=q.(names{ss});
        nexttile(1); semilogy(r.cfg.EbN0_dB_list,max(r.BER,1e-7),[marks{ss} '-'],'Color',colors(ss,:),'LineWidth',1.2);
        nexttile(2); semilogy(r.cfg.EbN0_dB_list,max(r.SER,1e-7),[marks{ss} '-'],'Color',colors(ss,:),'LineWidth',1.2); end
    nexttile(1); grid on; xlabel('Eb/N0 (dB)'); ylabel('BER'); legend(names,'Interpreter','none'); title(sprintf('U=%d equalized BER',U));
    nexttile(2); grid on; xlabel('Eb/N0 (dB)'); ylabel('SER'); legend(names,'Interpreter','none'); title(sprintf('U=%d equalized SER',U));
    exportgraphics(fig,fullfile(out,sprintf('u%d_ber_ser_%s.png',U,v.mode)),'Resolution',180); close(fig);
    fig=figure('Visible','off','Color','w'); r=q.stats_full;
    semilogy(r.cfg.EbN0_dB_list,max(r.BER_pre_equalizer,1e-7),'o--',r.cfg.EbN0_dB_list,max(r.BER,1e-7),'s-','LineWidth',1.2);
    grid on; xlabel('Eb/N0 (dB)'); ylabel('BER'); legend('before MMSE','after MMSE'); title(sprintf('U=%d equalizer effect',U));
    exportgraphics(fig,fullfile(out,sprintf('u%d_equalizer_%s.png',U,v.mode)),'Resolution',180); close(fig);
end
end
function local_text(v,path)
fid=fopen(path,'w'); c=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'Two-node communication validation %s L=%d nSym=%d\n',v.mode,v.channel_count_per_source,v.comm_cfg.n_sym);
for ww=1:numel(v.winds_mps), key=sprintf('u%d',v.winds_mps(ww)); fprintf(fid,'U=%d\n',v.winds_mps(ww));
    names=fieldnames(v.comparison.(key)); for ii=1:numel(names), q=v.comparison.(key).(names{ii});
        fprintf(fid,' %s BERlogRMSE=%.6g SERlogRMSE=%.6g BERCIoverlap=%.6g SERCIoverlap=%.6g\n', ...
            names{ii},q.ber_log10_rmse,q.ser_log10_rmse,q.ber_cluster_ci_overlap_fraction,q.ser_cluster_ci_overlap_fraction); end
end
end
