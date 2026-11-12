function [model,meta]=expand_kirchhoff_kstat_factor_basis_vertical(model)
%EXPAND_KIRCHHOFF_KSTAT_FACTOR_BASIS_VERTICAL Expand compressed factors in RAM.
%   The saved model remains compact. Call once after loading when repeated
%   sampling throughput is more important than the additional RAM.
arguments
    model (1,1) struct
end
if ~isfield(model,'factor_basis') || isempty(model.factor_basis)
    meta=struct('expanded',false,'elapsed_s',0,'additional_bytes',0); return
end
timer=tic; Q=model.factor_basis; before=sum(cellfun(@numel,model.factor_cells))*4;
for pp=1:numel(model.factor_cells)
    if ~model.pair_is_self(pp)
        expanded=single(Q*model.factor_cells{pp});
        model.factor_cells{pp}=complex(expanded,zeros(size(expanded),'single'));
    end
end
after=sum(cellfun(@numel,model.factor_cells))*4;
model.factor_basis=[]; model.factor_basis_mode='expanded full augmented-frequency factors in RAM';
model.runtime_factor_storage_bytes=after;
meta=struct('expanded',true,'elapsed_s',toc(timer), ...
    'additional_bytes',max(after-before,0),'runtime_factor_bytes',after, ...
    'memory_snapshot_bytes',local_memory_used_bytes(), ...
    'note','Runtime-only expansion; do not overwrite the compact on-disk model unless explicitly desired.');
end

function value=local_memory_used_bytes()
value=NaN; try, info=memory; value=info.MemUsedMATLAB; catch, end
end
