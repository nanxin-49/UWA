function [out,snapshots]=solve_deepsea3d_tfe_galerkin(cfg,opts)
%SOLVE_DEEPSEA3D_TFE_GALERKIN Paired collocation versus cubic-exact projection.
% Frozen solvers remain unchanged. Source/geometry interpolants and retained
% space are identical in both branches. Residuals use independent geometry.
timer=tic;cpu0=cputime;p=deepsea3d_fg_plan(cfg,opts);n=p.n;M=numel(p.s);N=opts.order;
estimate=12*16*n^2*M/2^20+18*16*(4*n)^2/2^20;
assert(estimate<=opts.max_workspace_MiB,'FG:MemoryBudget','Named-array estimate exceeds budget.');
forcing=deepsea3d_tfe_forcing(cfg,p.xx,p.yy,p.original_geo.eta,N);
W1=complex(zeros(n^2,M));W2=W1;total=W1;terms=complex(zeros(N,1));
norms=zeros(N,1);conditions=zeros(N,1);linear=zeros(N,1);snaps=opts.checkpoint_orders;
snapshots=repmat(struct('N',0,'W',[]),1,numel(snaps));next_snapshot=0;
out=struct('checkpoints',struct([]));
offset=cfg.rx_m(1:2)-cfg.center_m;
kernel=exp(1i*p.gamma*(cfg.rx_m(3)-p.b)).*exp(1i*p.kv*(offset(1)+p.L)).*exp(1i*p.kv.'*(offset(2)+p.L));
for j=1:N
    [A,~]=deepsea3d_fg_operators(W1,p,opts.product_kind);
    [~,B]=deepsea3d_fg_operators(W2,p,opts.product_kind);
    F=-A-B;conditions(j)=(norm(A(:))+norm(B(:)))/max(norm(F(:)),realmin);clear A B;
    F_f=reshape(fft2(reshape(F,n,n,M)),n^2,M).';clear F;
    bottom=fft2(deepsea3d_fg_resize(forcing{j+1},n));
    top=fft2(deepsea3d_fg_top_product(reshape(W1(:,end),n,n),p,opts.product_kind));
    W_f=deepsea3d_tfe_modal_solve(p.modal,F_f,bottom(:),top(:));
    I=p.modal.interior;R=p.D2(I,:)*W_f+p.gamma(:).'.^2.*W_f(I,:)-F_f(I,:);
    linear(j)=norm(R(:))/max(norm(F_f(I,:),'fro')+p.k^2*norm(W_f,'fro'),realmin);clear R F_f;
    W=reshape(ifft2(reshape(W_f.',n,n,M)),n^2,M);
    % The symmetric Galerkin subspace is explicit even for the original
    % collocation branch; check the solver never populates excluded lines.
    terms(j)=sum(reshape(W_f(end,:),n,n).*kernel,'all')/n^2;clear W_f;
    norms(j)=norm(W(:))*p.h/sqrt(n^2*M);total=total+W;W2=W1;W1=W;clear W;
    if any(j==snaps)
        r3=observe_deepsea3d_fg_residual(total,p,cfg,3);
        r4=observe_deepsea3d_fg_residual(total,p,cfg,4);
        row=struct('N',j,'rough_correction',sum(terms(1:j)), ...
            'terms',terms(1:j),'volume_coefficient_norms',norms(1:j), ...
            'recursion_cancellation',conditions(1:j),'linear_residual',linear(1:j), ...
            'independent3',r3,'independent4',r4);
        if isempty(out.checkpoints),out.checkpoints=row;else,out.checkpoints(end+1)=row;end
        next_snapshot=next_snapshot+1;snapshots(next_snapshot)=struct('N',j,'W',total);
        fprintf('FG %s N%d independent PDE %.6g retained %.6g boundary %.6g DtN %.6g\n', ...
            opts.product_kind,j,r4.physical_PDE,r4.retained_true_mapped_PDE,r4.true_boundary,r4.true_DtN);
    end
end
out.product_kind=opts.product_kind;out.options=opts;out.spacing_m=p.h;
out.geometry_Nyquist_relative=p.geometry_Nyquist_relative;
out.wall_s=toc(timer);out.cpu_s=cputime-cpu0;out.estimated_workspace_MiB=estimate;
out.snapshot_MiB=numel(snaps)*16*n^2*M/2^20;
out.time_convention='exp(-i omega t)';out.support_m=cfg.profile.support_m;
out.dealiasing='no intermediate projection; 2n padding for three primitive factors; top product included';
end
