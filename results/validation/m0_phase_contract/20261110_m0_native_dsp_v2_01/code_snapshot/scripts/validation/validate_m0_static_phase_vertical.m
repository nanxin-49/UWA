function result=validate_m0_static_phase_vertical(stage,run_dir)
%VALIDATE_M0_STATIC_PHASE_VERTICAL M0 controls; no M1 scans or external solve.
% Each stage is called through run_m0_phase_job_vertical into a new directory.
checks=struct('name',{},'value',{},'limit',{},'pass',{});
data=struct();
switch stage
    case 'C_receipt_recovery'
        saved=load(fullfile(run_dir,'C_gaussian','payload.mat'),'payload');
        checks=saved.payload.checks; data=saved.payload.data;
        data.recovery_note='Numerics retained from C_gaussian; empty-log SHA receipt failed after numerical PASS. No solve repeated.';
    case 'A_algebra'
        f=[4000;4062.5;6000;8000]; g=local_geometry();
        d=[1+2i;-0.2+0.7i;3-1i;0.1-0.3i];
        r=[-.3+.8i;1.2-.7i;.4+1i;-.6-.2i];
        reduced=struct('direct_f',d,'reflect_coh_f',.4*r,'reflect_sca_fm',r*[.6,1.6]);
        [a,m]=apply_pe_channel_phase_reference_vertical(f,reduced,g,'direct_dsp');
        b=apply_pe_channel_phase_reference_vertical(f,reduced,g,'absolute_physical');
        c=apply_pe_channel_phase_reference_vertical(f,reduced,g,'legacy_reduced');
        D=exp(-1i*2*pi*f*(6/1500));
        checks(end+1)=gate('DSP_direct',max(abs(a.direct_f-conj(d))),1e-12);
        checks(end+1)=gate('DSP_reflect',max(abs(a.reflect_fm-conj(r*[1,2]).*D),[],'all'),1e-12);
        checks(end+1)=gate('legacy',max(abs(c.total_fm-(d+r*[1,2])),[],'all'),1e-12);
        checks(end+1)=gate('physical_to_DSP',max(abs(a.total_fm- ...
            conj(b.total_fm).*exp(1i*2*pi*f*(97/1500))),[],'all'),1e-12);
        checks(end+1)=gate('4062_factor_minus_i',abs(m.reflect_factor_f(2)+1i),1e-12);
        adapted=reduced; names=fieldnames(adapted);
        for j=1:numel(names), adapted.(names{j})=conj(adapted.(names{j})); end
        e=apply_pe_channel_phase_reference_vertical(f,adapted,g,'direct_dsp', ...
            struct('input_representation','dsp_adapted_reduced'));
        checks(end+1)=gate('external_adapter_no_double_conjugation',max(abs(e.total_fm-a.total_fm),[],'all'),1e-12);
        checks(end+1)=gate('component_closure',max(abs(a.total_fm-a.direct_f-a.reflect_fm),[],'all'),1e-12);
        checks(end+1)=gate('reject_DSP_ready_reconversion',~rejected(@() ...
            apply_pe_channel_phase_reference_vertical(f,reduced,g,'direct_dsp', ...
            struct('input_representation','direct_dsp'))),0);
        checks(end+1)=gate('reject_adapted_to_physical',~rejected(@() ...
            apply_pe_channel_phase_reference_vertical(f,adapted,g,'absolute_physical', ...
            struct('input_representation','dsp_adapted_reduced'))),0);
        data=struct('f',f,'native',reduced,'dsp',a,'physical',b,'meta',m);
    case {'B_zero','B_positive','B_negative','C_gaussian','C_direct_only'}
        p=local_public(4062.5,256,50); eta=0;
        if strcmp(stage,'B_positive'), eta=.01; end
        if strcmp(stage,'B_negative'), eta=-.01; end
        if startsWith(stage,'B_')
            p.source_mode='custom_field_fn'; p.source_field_fn=@local_plane_source;
        end
        p.surface_elevation_override_xy=eta*ones(p.ny,p.nx);
        if strcmp(stage,'C_direct_only'), p.enable_surface_reflection=false; end
        out=vertical_channel_model(p); f=out.f_axis(:); k=2*pi*f/1500;
        if startsWith(stage,'B_')
            % Independent physical geometry: z_surface=-eta, Lr=103+2eta.
            Pd=exp(1i*k*97); Pr=-exp(1i*k*(103+2*eta));
            Hd=ones(size(f)); Hr=-exp(-1i*k*(6+2*eta));
            expected_relative_delay_s=(6+2*eta)/1500;
            checks(end+1)=gate('height_phase_sign', ...
                abs(out.delta_phi(1)-2*k*eta),1e-12);
        else
            [Pd,Pr]=local_gaussian_reference(out,p);
            Hd=conj(Pd).*exp(1i*k*97); Hr=conj(Pr).*exp(1i*k*97);
            expected_relative_delay_s=.004; % nominal carrier, not beam group delay
        end
        if ~p.enable_surface_reflection, Pr=0*Pr; Hr=0*Hr; end
        checks(end+1)=gate('direct_physical_relative',rel(out.H_direct_physical_f,Pd),1e-10);
        checks(end+1)=gate('reflect_physical_relative',rel(out.H_reflect_physical_f,Pr),1e-10);
        checks(end+1)=gate('direct_DSP_relative',rel(out.H_direct_f,Hd),1e-10);
        checks(end+1)=gate('reflect_DSP_relative',rel(out.H_reflect_f,Hr),1e-10);
        checks(end+1)=gate('total_direct_scaled',max(abs(out.H_f-Hd-Hr))/max(max(abs(Hd)),eps),1e-10);
        checks(end+1)=gate('component_closure',max(abs(out.H_f-out.H_direct_f-out.H_reflect_f)),1e-12);
        checks(end+1)=gate('scalar_reference',max(abs([out.h_direct-out.H_direct_f(out.idx_f_ref), ...
            out.h_reflect-out.H_reflect_f(out.idx_f_ref),out.h_total-out.H_f(out.idx_f_ref)])),1e-12);
        if startsWith(stage,'B_')
            % Derive delay perturbation from a ratio with |phase|<pi, no fit/unwrap.
            observed_extra_s=-angle(out.H_reflect_f/(-exp(-1i*k*6)))/(2*pi*f);
            checks(end+1)=gate('height_delay_seconds',abs(observed_extra_s-2*eta/1500),1e-12);
        end
        data=struct('params',rmfield_if_present(p,'source_field_fn'),'eta_up_m',eta, ...
            'physical_surface_depth_m',-eta,'expected_relative_delay_s',expected_relative_delay_s, ...
            'out',out,'reference_direct_physical',Pd,'reference_reflect_physical',Pr, ...
            'reference_direct_dsp',Hd,'reference_reflect_dsp',Hr);
        % Function handles are not retained as numerical evidence.
        data.out.config=rmfield_if_present(data.out.config,'source_field_fn');
    case 'D_statistics'
        [checks,data]=local_statistics();
    case 'G_contract_edges'
        q=load(fullfile(run_dir,'D_statistics','payload.mat'),'payload'); s=q.payload.data;
        checks(end+1)=gate('old_DSP_augmented_covariance',rel(s.old_migrated.stats.augmented_covariance, ...
            s.migrated.stats.augmented_covariance),1e-12);
        checks(end+1)=gate('old_DSP_deterministic_direct',rel(s.old_migrated.H_direct_f,s.migrated.H_direct_f),1e-12);
        checks(end+1)=gate('old_DSP_deterministic_coherent',rel(s.old_migrated.H_ref_coh_f,s.migrated.H_ref_coh_f),1e-12);
        U=s.old_migrated.stats.eigenvectors; lambda=s.old_migrated.stats.eigenvalues;
        checks(end+1)=gate('old_DSP_EVD',rel(s.migrated.stats.C_scatter_f*U,U.*lambda.'),1e-10);
        cfg=struct('x_tx',0,'y_tx',0,'sigma_src_m',.3,'f0',4000);
        x=(-128:127)*(50/256); [X,Y]=meshgrid(x,x);
        fields=cell(3,1); freq=[4000,6000,8000];
        for k=1:3, cfg.f0=freq(k); fields{k}=gaussian_source_initial_field_vertical(X,Y,cfg); end
        checks(end+1)=gate('Gaussian_unit_peak',abs(max(abs(fields{1}),[],'all')-1),0);
        checks(end+1)=gate('Gaussian_frequency_independence',max(abs(fields{1}-fields{3}),[],'all'),0);
        checks(end+1)=gate('Gaussian_independent_expression',max(abs(fields{2}-exp(-(X.^2+Y.^2)/.18)),[],'all'),1e-12);
        data=struct('frequencies_hz',freq,'peak',max(abs(fields{1}),[],'all'), ...
            'old_migration_meta',s.old_migrated.phase_reference_meta);
    case 'H_metadata_guards'
        f=4062.5; g=local_geometry(); z=struct('direct_f',1+2i,'reflect_fm',-.3+.4i);
        [dsp,m]=apply_pe_channel_phase_reference_vertical(f,z,g,'direct_dsp');
        checks(end+1)=gate('native_four_argument_contract',abs(dsp.direct_f-(1-2i)),1e-12);
        z.phase_reference_meta=m;
        checks(end+1)=gate('reject_metadata_DSP_reconversion',~rejected(@() ...
            apply_pe_channel_phase_reference_vertical(f,z,g,'direct_dsp')),0);
        m.target_reference='absolute_physical'; z.phase_reference_meta=m;
        checks(end+1)=gate('reject_physical_as_native',~rejected(@() ...
            apply_pe_channel_phase_reference_vertical(f,z,g,'direct_dsp')),0);
        data=struct('meta',m,'dsp',dsp);
    case 'F_external_empirical'
        root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
        suite=testsuite(fullfile(root,'scripts','validation','test_spectral_fast_kirchhoff_multifrequency_phase.m'));
        % Fixed PM grid identity is unrelated to M0 Fourier conversion and
        % requires a separate historical fixture; the four phase tests run.
        suite=suite(~contains({suite.Name},'FixedSurfaceNestedGrid'));
        tests=run(suite);
        checks(end+1)=gate('external_phase_test_count',abs(numel(tests)-4),0);
        checks(end+1)=gate('external_phase_failed_tests',sum(~[tests.Passed]),0);
        data.test_names={tests.Name}; data.test_passed=[tests.Passed];
        e=load(fullfile(run_dir,'E_cir_comm_boundary_fixed','payload.mat'),'payload'); e=e.payload.data;
        base_params=struct('sea_hs_target',0,'sea_wind_speed',5,'f0',6000);
        summary_table=table([1;1],[0;0],[1;1],[0;0],[0;0],[0;0],[0;0], ...
            'VariableNames',{'abs_h_ref','phase_h_ref_rad','reflect_direct_abs','reflect_direct_phase_rad', ...
            'redistribution_reflect_rms_delta_k_rad_per_m','redistribution_reflect_high_k_fraction','tap_rms_delay_symbols'});
        mc_stats=struct('phase_reference_meta',e.meta, ...
            'H_f',struct('samples',[e.H,e.H],'f_axis',e.f,'idx_f_ref',33));
        fixture=fullfile(run_dir,stage,'empirical_fixture.mat');
        save(fixture,'base_params','summary_table','mc_stats');
        empirical=build_surface_empirical_channel_model_vertical(fixture);
        sample=sample_surface_empirical_channel_vertical(empirical,struct(),3,910106,'wideband_hf');
        tapsample=sample_surface_empirical_channel_vertical(empirical, ...
            struct('symbol_rate_hz',1000,'tap_fft_len',64),2,910106,'tap_level');
        checks(end+1)=gate('empirical_DSP_pass_through',max(abs(sample.H_f_samples-e.H),[],'all'),0);
        checks(end+1)=gate('empirical_taps_finite',~all(isfinite(tapsample.tap_h_taps{1})),0);
        old=empirical; old.conditions=rmfield(old.conditions,'phase_reference_meta');
        checks(end+1)=gate('empirical_reject_unmarked',~rejected(@() ...
            sample_surface_empirical_channel_vertical(old,struct(),1,1,'wideband_hf')),0);
        data.empirical=empirical; data.sample=sample; data.tap_sample=tapsample;
    case {'E_cir_comm','E_cir_comm_boundary_fixed'}
        [checks,data]=local_cir_comm();
    case {'static_analysis','static_analysis_final'}
        manifest_path=fullfile(run_dir,'code_manifest_review.json');
        if ~isfile(manifest_path), manifest_path=fullfile(run_dir,'code_manifest.json'); end
        manifest=jsondecode(fileread(manifest_path));
        issues=struct('file',{},'id',{},'message',{},'line',{});
        for j=1:numel(manifest)
            if ~endsWith(manifest(j).path,'.m'), continue; end
            messages=checkcode(manifest(j).path,'-id');
            for k=1:numel(messages)
                issues(end+1)=struct('file',manifest(j).path,'id',messages(k).id, ...
                    'message',messages(k).message,'line',messages(k).line); %#ok<AGROW>
            end
        end
        parse_count=sum(contains({issues.id},'PARSE')|contains({issues.id},'SYNER'));
        checks(end+1)=gate('parse_errors',parse_count,0); data.issues=issues;
    case {'independent_audit','independent_audit_lossless','independent_audit_lossless_v2','independent_audit_final'}
        [checks,data]=audit_m0_static_phase_vertical(run_dir);
    case {'dashboard','dashboard_final'}
        output=render_project_progress_dashboard(); data.output=output;
        checks(end+1)=gate('render_failed_checks',output.failed_checks,0);
    otherwise
        error('Unknown M0 stage %s.',stage);
end
result=struct('stage',stage,'checks',checks,'pass',all([checks.pass]),'data',data);
writetable(struct2table(checks),fullfile(run_dir,stage,'metrics.csv'));
disp(struct2table(checks));
end

function [checks,data]=local_statistics()
checks=struct('name',{},'value',{},'limit',{},'pass',{});
f=[4000;4062.5;6000]; p=local_public(f.',64,25); p.sigma_src_m=.4;
p.surface_elevation_override_xy=zeros(64); p.surface_roughness_scale_mode='raw_pm';
cfg=struct('f_axis_hz',f,'c0',1500,'nx',64,'ny',64,'xw',25,'yw',25, ...
    'x_tx',0,'y_tx',0,'z_tx',100,'x_rx',0,'y_rx',0,'z_rx',3, ...
    'sigma_src_m',.4,'stepz_lamb',.5,'sponge_ratio',.12,'alpha_max_np_per_m',0, ...
    'reflect_coeff',-1,'env_mode','uniform','enable_bubbles',false);
spec=raw_pm_spectrum_grid_vertical(5,64,64,25,25);
cache=build_cached_joint_kstat_pe_executor_vertical(cfg,spec);
[projection,~]=build_adjoint_receiver_projection_vertical(cache,spec);
rng(910101,'twister'); delta=.1*randn(64,64,3,4)+1i*.04*randn(64,64,3,4);
forward=run_cached_joint_kstat_pe_executor_vertical(cache,delta);
adjoint=run_adjoint_receiver_projection_vertical(projection,delta);
checks(end+1)=gate('cached_adjoint',rel(forward.H_total_fm,adjoint.H_total_fm),1e-10);
D=exp(-1i*2*pi*f*.004);
checks(end+1)=gate('cached_independent_DSP',rel(forward.H_ref_sca_fm,D.*conj(forward.H_ref_sca_reduced_fm)),1e-12);
flat_delta=complex(zeros(64,64,3));
for j=1:3, flat_delta(:,:,j)=-1-cache.R_coh_f(j); end
flat=run_cached_joint_kstat_pe_executor_vertical(cache,flat_delta);
public=vertical_channel_model(p);
checks(end+1)=gate('public_cached_total',rel(public.H_f,flat.H_total_fm),1e-10);
checks(end+1)=gate('public_cached_direct',rel(public.H_direct_f,flat.H_direct_f),1e-10);
checks(end+1)=gate('cached_component',max(abs(forward.H_total_fm-forward.H_direct_f-forward.H_ref_coh_f-forward.H_ref_sca_fm),[],'all'),1e-12);
% Dense/FFT contraction on an independently small PM/PE grid, same formulas.
small=cfg; small.nx=8; small.ny=8; small.xw=3.125; small.yw=3.125; small.z_tx=10; small.z_rx=1;
ps=raw_pm_spectrum_grid_vertical(5,16,16,6.25,6.25);
cc=build_cached_joint_kstat_pe_executor_vertical(small,ps);
pp=build_adjoint_receiver_projection_vertical(cc,ps);
dense=contract_kstat_receiver_stats_vertical(ps,pp,struct('method','dense'));
fast=contract_kstat_receiver_stats_vertical(ps,pp,struct('method','fft'));
Ds=diag(exp(-1i*2*pi*f*(2/1500)));
checks(end+1)=gate('dense_fft_C',rel(dense.C_H,fast.C_H),1e-10);
checks(end+1)=gate('dense_fft_P',rel(dense.P_H,fast.P_H),1e-10);
checks(end+1)=gate('analytic_C_transform',rel(dense.C_H,Ds*conj(dense.C_H_reduced)*Ds'),1e-12);
checks(end+1)=gate('analytic_P_transform',rel(dense.P_H,Ds*conj(dense.P_H_reduced)*Ds.'),1e-12);
% Nonreal cross-frequency covariance and nonzero pseudo-covariance.
rng(910102,'twister'); x=randn(3,160); y=randn(3,160);
native=[1+.2i;.3-.4i;-.1+.8i]+(eye(3)+.2i*ones(3))*x+1i*.25*y;
opts=struct('frequency_axis_hz',f,'H_direct_f',cache.H_direct_reduced_f, ...
    'H_ref_coh_f',cache.H_ref_coh_reduced_f,'condition',struct('wind_speed_mps',5));
legacy=estimate_conditional_channel_stats_vertical(native,opts); legacy.schema_version='1.0.0';
legacy=rmfield(legacy,'phase_reference_meta'); g=local_geometry(); g.source_representation='native_pe_reduced';
migrated=upgrade_conditional_channel_phase_vertical(legacy,g);
target=D.*conj(native); centered=target-mean(target,2);
C=centered*centered'/159; P=centered*centered.'/159;
checks(end+1)=gate('sample_C',rel(migrated.stats.C_scatter_f,C),1e-12);
checks(end+1)=gate('sample_P',rel(migrated.stats.P_scatter_f,P),1e-12);
checks(end+1)=gate('nonzero_P_control',double(norm(P,'fro')<.1 || norm(imag(C),'fro')<.01),0);
checks(end+1)=gate('sample_mean',rel(migrated.stats.mu_scatter_f,mean(target,2)),1e-12);
Xa=[real(target);imag(target)]; Xa=Xa-mean(Xa,2); Sigma=Xa*Xa.'/159;
checks(end+1)=gate('augmented_covariance',rel(migrated.stats.augmented_covariance,Sigma),1e-12);
checks(end+1)=gate('complex_EVD',rel(C*migrated.stats.eigenvectors,migrated.stats.eigenvectors.*migrated.stats.eigenvalues.'),1e-10);
checks(end+1)=gate('augmented_EVD',rel(Sigma*migrated.stats.augmented_eigenvectors, ...
    migrated.stats.augmented_eigenvectors.*migrated.stats.augmented_eigenvalues.'),1e-10);
old=estimate_conditional_channel_stats_vertical(D.*native,opts); old.schema_version='2.0.0';
old.H_ref_coh_f=D.*cache.H_ref_coh_reduced_f;
g.source_representation='legacy_carrier_only_direct_dsp'; old_migrated=upgrade_conditional_channel_phase_vertical(old,g);
checks(end+1)=gate('carrier_only_migration_C',rel(old_migrated.stats.C_scatter_f,C),1e-12);
checks(end+1)=gate('carrier_only_migration_P',rel(old_migrated.stats.P_scatter_f,P),1e-12);
checks(end+1)=gate('carrier_only_migration_mean',rel(old_migrated.stats.mu_scatter_f,mean(target,2)),1e-12);
checks(end+1)=gate('reject_unmarked_old_schema',~rejected(@()upgrade_conditional_channel_phase_vertical(old,local_geometry())),0);
draw=sample_conditional_channel_vertical(migrated,12,910103,struct('path','improper','rank','full'));
checks(end+1)=gate('sample_component',max(abs(draw.H_total_f-migrated.H_direct_f-migrated.H_ref_coh_f-draw.H_ref_sca_f),[],'all'),1e-12);
rng(910103,'twister'); s=migrated.stats; rr=s.rank_candidates.augmented_full;
manual=s.augmented_mean+s.augmented_eigenvectors(:,1:rr)*(sqrt(s.augmented_eigenvalues(1:rr)).*randn(rr,12));
checks(end+1)=gate('improper_sampler_exact',rel(draw.H_ref_sca_f,manual(1:3,:)+1i*manual(4:6,:)),1e-12);
other=migrated; other.condition.wind_speed_mps=8;
library=build_conditional_channel_library_vertical({migrated,other},struct('fixed_geometry',local_geometry()));
libdraw=sample_conditional_channel_library_vertical(library,8,4,910104,struct('path','proper','rank','full'));
checks(end+1)=gate('library_metadata',~strcmp(libdraw.phase_reference_meta.phase_contract_version,'2.0.0'),0);
pair=sample_conditional_channel_rank_pair_vertical(migrated,4,910105,'99.9');
checks(end+1)=gate('rank_pair_finite',~all(isfinite(pair.H_full_f),'all'),0);
% Exercise production SSA split conversion with a nonzero roughness control.
sp=p; sp.f0=4062.5; sp.surface_elevation_override_xy=[]; sp.sea_hs_target=.01;
sp.surface_boundary_model='ssa_stat_kernel'; sp.surface_ssa_component_diagnostics=true;
ssa=vertical_channel_model(sp); sm=ssa.surface_ssa_component_meta;
checks(end+1)=gate('SSA_split',max(abs(sm.h_reflect_total_f-sm.h_reflect_coh_f-sm.h_reflect_sca_f)),1e-12);
checks(end+1)=gate('SSA_public_match',rel(sm.h_reflect_total_f,ssa.H_reflect_f),1e-10);
checks(end+1)=gate('SSA_native_to_DSP',rel(sm.h_reflect_sca_f, ...
    ssa.phase_reference_meta.reflect_factor_f.*conj(sm.h_reflect_sca_reduced_f)),1e-12);
data=struct('forward',forward,'adjoint',adjoint,'flat',flat,'public',public, ...
    'dense',dense,'fast',fast,'native_samples',native,'target_samples',target, ...
    'legacy',legacy,'migrated',migrated,'old_carrier_only',old,'old_migrated',old_migrated, ...
    'draw',draw,'library',library,'ssa',ssa,'cfg',cfg,'spec',spec);
end

function [checks,data]=local_cir_comm()
checks=struct('name',{},'value',{},'limit',{},'pass',{});
f=(4000:62.5:7937.5)'; td=97/1500; tr=103/1500; g=.3+.2i;
Hd=(1+.4i)*ones(64,1); Hr=g*exp(-1i*2*pi*f*.004); H=Hd+Hr;
P=conj(H).*exp(1i*2*pi*f*td);
[~,meta]=apply_pe_channel_phase_reference_vertical(f,struct('direct_f',ones(64,1)),local_geometry(),'direct_dsp');
c=build_channel_cir_vertical(H,f,struct('input_reference','direct_dsp','phase_reference_meta',meta));
physical=build_channel_cir_vertical(P,f,struct('input_reference','absolute_physical','time_origin_shift_s',-td));
checks(end+1)=gate('CIR_physical_DSP_relation',max(abs(c.h_tau-conj(physical.h_tau))),1e-12);
checks(end+1)=gate('CIR_fft_roundtrip',max(abs(fft(c.h_tau)-H)),1e-12);
cr=build_channel_cir_vertical(Hr,f,struct('input_reference','direct_dsp','phase_reference_meta',meta));
[~,ir]=max(abs(cr.h_tau)); checks(end+1)=gate('CIR_arrival_seconds',abs(cr.delay_axis_s(ir)-.004),1e-12);
shifted=build_channel_cir_vertical(H,f,struct('time_origin_shift_s',.001,'phase_reference_meta',meta));
checks(end+1)=gate('CIR_common_shift',max(abs(shifted.h_tau-circshift(c.h_tau,4))),1e-12);
old_meta=rmfield(meta,'phase_contract_version');
checks(end+1)=gate('reject_old_CIR_metadata',~rejected(@()build_channel_cir_vertical(H,f, ...
    struct('phase_reference_meta',old_meta))),0);
fs=1000; fc=6000; nu=(-32:31)'*fs/64; taps=[1;.2+.1i;0;.05-.03i];
analytic=exp(-1i*2*pi*nu*(0:3)/fs)*taps;
channel=struct('H_f',analytic,'f_axis_hz',fc+nu,'phase_reference_meta',meta);
[ht,tm]=build_communication_taps_vertical(channel,struct('symbol_rate_hz',fs,'n_fft',64, ...
    'f_ref_hz',fc,'tap_energy_ratio',1,'circular_peak_align',false));
checks(end+1)=gate('tap_fft_closure',tm.frequency_ifft_closure_relative_error,1e-12);
checks(end+1)=gate('analytic_taps',max(abs(tm.h_full-[taps;zeros(60,1)])),1e-12);
set=struct('H_fm',analytic,'f_axis_hz',fc+nu,'phase_reference_meta',meta);
comm=evaluate_mpsk_channel_ensemble_vertical(set,struct('M',4,'n_sym',128,'symbol_rate_hz',fs, ...
    'tap_fft_size',64,'f_ref_hz',fc,'tap_energy_ratio',1,'circular_peak_align',false, ...
    'EbN0_dB_list',[],'bootstrap_count',10,'mmse_reg_eps',1e-10));
checks(end+1)=gate('known_CSI_noiseless_BER',comm.noiseless_BER,0);
checks(end+1)=gate('known_CSI_noiseless_SER',comm.noiseless_SER,0);
channel.phase_reference_meta=old_meta;
checks(end+1)=gate('reject_old_tap_metadata',~rejected(@()build_communication_taps_vertical(channel)),0);
% A newly computed small production wideband output feeds the same consumer.
pp=local_public([5000,6000,7000],64,25); pp.z_tx=10; pp.z_max=10; pp.z_rx=1; pp.sigma_src_m=.4;
pp.surface_elevation_override_xy=zeros(64); produced=vertical_channel_model(pp);
consumer=evaluate_mpsk_channel_ensemble_vertical(struct('H_fm',produced.H_f,'f_axis_hz',produced.f_axis(:), ...
    'phase_reference_meta',produced.phase_reference_meta), ...
    struct('M',4,'n_sym',128,'tap_fft_size',64,'symbol_rate_hz',1000,'f_ref_hz',6000, ...
    'EbN0_dB_list',[],'bootstrap_count',10,'mmse_reg_eps',1e-12));
checks(end+1)=gate('production_wideband_noiseless_BER',consumer.noiseless_BER,0);
data=struct('f',f,'H',H,'P',P,'cir',c,'physical_cir',physical,'reflect_cir',cr, ...
    'meta',meta,'tap_meta',tm,'taps',ht,'comm',comm,'production',produced,'production_comm',consumer);
end

function [Pd,Pr]=local_gaussian_reference(out,p)
% Independent Helmholtz angular spectrum, full carrier (not helper metadata).
[X,Y]=meshgrid(out.x,out.y); initial=exp(-(X.^2+Y.^2)/(2*.3^2));
A=fft2(initial); n=p.nx; K=2*pi/p.xw*[0:n/2-1,-n/2:-1];
[KX,KY]=meshgrid(K,K); k=2*pi*out.f_axis(:)/1500;
Pd=complex(zeros(size(k))); Pr=Pd;
[~,ix]=min(abs(out.x)); [~,iy]=min(abs(out.y));
for j=1:numel(k)
    kz=sqrt(complex(k(j)^2-KX.^2-KY.^2,0));
    direct=ifft2(A.*exp(1i*97*kz)); reflected=-ifft2(A.*exp(1i*103*kz));
    Pd(j)=direct(iy,ix); Pr(j)=reflected(iy,ix);
end
end
function p=local_public(f,n,w)
p=struct('f0',f,'enable_wideband',false,'c0',1500,'z_max',100,'z_tx',100,'z_rx',3, ...
    'x_tx',0,'y_tx',0,'x_rx',0,'y_rx',0,'nx',n,'ny',n,'xw',w,'yw',w, ...
    'sigma_src_m',.3,'stepz_lamb',.5,'show_figures',false,'enforce_1_over_R',false, ...
    'save_mode','rx_only','use_gpu',false,'enable_bubbles',false,'env_mode','uniform', ...
    'enable_surface_reflection',true,'surface_boundary_model','kirchhoff_spatial', ...
    'surface_reflect_coeff',-1,'surface_phase_mode','normal','sea_hs_target',0, ...
    'alpha_max_np_per_m',0,'sea_seed',910100,'channel_phase_reference','direct_dsp');
end
function [psi,meta]=local_plane_source(X,~,~,~)
psi=complex(ones(size(X))); meta=struct('expression','1','normalization','unit amplitude');
end
function g=local_geometry()
g=struct('z_tx',100,'z_rx',3,'z_surface',0,'c0',1500);
end
function s=rmfield_if_present(s,name)
if isfield(s,name), s=rmfield(s,name); end
end
function c=gate(name,value,limit)
c=struct('name',name,'value',double(value),'limit',limit,'pass',isfinite(value)&&value<=limit);
end
function value=rel(a,b)
value=norm(a(:)-b(:))/max(norm(b(:)),eps);
end
function tf=rejected(action)
tf=false; try, action(); catch, tf=true; end
end
