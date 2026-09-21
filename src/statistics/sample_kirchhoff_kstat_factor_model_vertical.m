function [deltaG_xy_fm,meta] = sample_kirchhoff_kstat_factor_model_vertical( ...
    model,n_realizations,rng_seed,mode)
%SAMPLE_KIRCHHOFF_KSTAT_FACTOR_MODEL_VERTICAL Draw a reusable model batch.

arguments
    model (1,1) struct
    n_realizations (1,1) double {mustBeInteger,mustBePositive}
    rng_seed (1,1) double {mustBeFinite}
    mode (1,:) char {mustBeMember(mode,{'joint','independent'})}
end
rng(mod(round(rng_seed),2^32),'twister');
ny=model.ny; nx=model.nx; F=model.F;
timer=tic;
X=complex(zeros(ny,nx,F,n_realizations,'single'));
if strcmp(mode,'independent')
    for ii=1:F
        z=(complex(randn(ny,nx,n_realizations,'single'), ...
            randn(ny,nx,n_realizations,'single')))/sqrt(single(2));
        X(:,:,ii,:)=reshape(model.sqrt_q_independent(:,:,ii),ny,nx,1,1).* ...
            reshape(z,ny,nx,1,n_realizations);
    end
else
    for pp=1:numel(model.factor_cells)
        factor=model.factor_cells{pp};
        rank=size(factor,2);
        if model.pair_is_self(pp)
            z=randn(rank,n_realizations,'single');
            draw=factor*z;
            values=complex(draw(1:F,:),draw(F+1:end,:));
            X(model.pair_iy(pp),model.pair_ix(pp),:,:)= ...
                reshape(values,1,1,F,n_realizations);
        else
            z=complex(randn(rank,n_realizations,'single'), ...
                randn(rank,n_realizations,'single'))/sqrt(single(2));
            if isfield(model,'factor_basis') && ~isempty(model.factor_basis)
                draw=model.factor_basis*(factor*z);
            else
                draw=factor*z;
            end
            X(model.pair_iy(pp),model.pair_ix(pp),:,:)= ...
                reshape(draw(1:F,:),1,1,F,n_realizations);
            X(model.pair_jy(pp),model.pair_jx(pp),:,:)= ...
                reshape(conj(draw(F+1:end,:)),1,1,F,n_realizations);
        end
    end
end

deltaG_xy_fm=complex(zeros(ny,nx,F,n_realizations,'single'));
for mm=1:n_realizations
    for ii=1:F
        deltaG_xy_fm(:,:,ii,mm)=ifft2(X(:,:,ii,mm));
    end
end
peak_memory_snapshot_bytes=local_memory_used_bytes();
meta=struct('mode',mode,'n_realizations',n_realizations, ...
    'rng_seed',mod(round(rng_seed),2^32), ...
    'elapsed_s',toc(timer), ...
    'output_bytes',numel(deltaG_xy_fm)*8, ...
    'ifft_count',F*n_realizations, ...
    'peak_memory_snapshot_bytes',peak_memory_snapshot_bytes);
end

function value=local_memory_used_bytes()
value=NaN;
try
    info=memory;
    value=info.MemUsedMATLAB;
catch
end
end
