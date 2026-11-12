function [H_ref_sca_fm,timing,max_component_error]=generate_cached_kstat_ensemble_vertical(cache,joint_model,L,seed_base,mode,batch_size)
%GENERATE_CACHED_KSTAT_ENSEMBLE_VERTICAL Joint/independent factor-model ensemble.
arguments
    cache (1,1) struct
    joint_model (1,1) struct
    L (1,1) double {mustBeInteger,mustBePositive}
    seed_base (1,1) double {mustBeFinite}
    mode char {mustBeMember(mode,{'joint','independent'})}
    batch_size (1,1) double {mustBeInteger,mustBePositive} = 4
end
F=joint_model.F; H_ref_sca_fm=complex(zeros(F,L)); nb=ceil(L/batch_size);
generator_s=zeros(nb,1); pe_s=zeros(nb,1); memory_bytes=nan(nb,1); max_component_error=0;
for bb=1:nb
    idx=(bb-1)*batch_size+1:min(bb*batch_size,L);
    [delta,g]=sample_kirchhoff_kstat_factor_model_vertical(joint_model,numel(idx),seed_base+bb,mode);
    [r,p]=run_cached_joint_kstat_pe_executor_vertical(cache,delta);
    H_ref_sca_fm(:,idx)=r.H_ref_sca_fm; generator_s(bb)=g.elapsed_s; pe_s(bb)=p.elapsed_s;
    memory_bytes(bb)=max(g.peak_memory_snapshot_bytes,p.memory_snapshot_bytes);
    max_component_error=max(max_component_error,r.component_sum_max_abs_error);
end
timing=struct('L',L,'batch_size',batch_size,'generator_s',generator_s,'pe_s',pe_s, ...
    'total_s',sum(generator_s+pe_s),'per_realization_s',sum(generator_s+pe_s)/L, ...
    'peak_memory_snapshot_bytes',max(memory_bytes,[],'omitnan'));
end
