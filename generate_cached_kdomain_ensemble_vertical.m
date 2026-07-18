function [H_ref_sca_fm,timing,mapping,max_component_error]=generate_cached_kdomain_ensemble_vertical(cache,spec,seeds,batch_size)
%GENERATE_CACHED_KDOMAIN_ENSEMBLE_VERTICAL Explicit shared-surface receiver ensemble.
arguments
    cache (1,1) struct
    spec (1,1) struct
    seeds (1,:) double
    batch_size (1,1) double {mustBeInteger,mustBePositive} = 4
end
F=numel(cache.f_axis_hz); L=numel(seeds); H_ref_sca_fm=complex(zeros(F,L));
mapping=struct('full_variance_m2',zeros(L,1),'crop_variance_m2',zeros(L,1));
nb=ceil(L/batch_size); generator_s=zeros(nb,1); pe_s=zeros(nb,1); memory_bytes=nan(nb,1);
max_component_error=0;
for bb=1:nb
    idx=(bb-1)*batch_size+1:min(bb*batch_size,L); B=numel(idx); timer=tic;
    delta=complex(zeros(cache.cfg.ny,cache.cfg.nx,F,B,'single'));
    for mm=1:B
        [eta,meta]=sample_raw_pm_surface_vertical(spec,seeds(idx(mm)));
        eta_pe=eta(cache.pm_mapping.iy,cache.pm_mapping.ix);
        mapping.full_variance_m2(idx(mm))=meta.variance_about_zero_m2;
        mapping.crop_variance_m2(idx(mm))=mean(eta_pe(:).^2);
        for ff=1:F
            alpha=4*pi*cache.f_axis_hz(ff)/cache.cfg.c0;
            delta(:,:,ff,mm)=single(cache.cfg.reflect_coeff*exp(1i*alpha*eta_pe)-cache.R_coh_f(ff));
        end
    end
    generator_s(bb)=toc(timer);
    [r,m]=run_cached_joint_kstat_pe_executor_vertical(cache,delta);
    H_ref_sca_fm(:,idx)=r.H_ref_sca_fm; pe_s(bb)=m.elapsed_s;
    memory_bytes(bb)=m.memory_snapshot_bytes;
    max_component_error=max(max_component_error,r.component_sum_max_abs_error);
end
timing=struct('L',L,'batch_size',batch_size,'generator_s',generator_s,'pe_s',pe_s, ...
    'total_s',sum(generator_s+pe_s),'per_realization_s',sum(generator_s+pe_s)/L, ...
    'peak_memory_snapshot_bytes',max(memory_bytes,[],'omitnan'));
end
