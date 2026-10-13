function out=deepsea3d_mapped_discrete_residual(W,g,p)
%DEEPSEA3D_MAPPED_DISCRETE_RESIDUAL Unpreconditioned PDE and both boundaries.
n=p.n;M=numel(p.s);[A,B]=deepsea3d_fg_operators(W,p,"galerkin");
F=fft2(reshape(W,n,n,M));
H=reshape(ifft2((p.k^2-p.kv.^2-p.kv.'.^2).*F),n^2,M);
R=H+W*p.D2.'+A+B;I=2:M-1;
bottom=reshape(W(:,1),n,n)-g;
T=ifft2(1i*p.gamma.*fft2(reshape(W(:,end),n,n)));
top=W*p.D(end,:).'-T(:)-reshape(deepsea3d_fg_top_product(reshape(W(:,end),n,n),p,"galerkin"),[],1);
out=struct('mapped_PDE',norm(R(:,I),'fro')/max(p.k^2*norm(W(:,I),'fro'),realmin), ...
    'Dirichlet',norm(bottom,'fro')/max(norm(g,'fro'),realmin), ...
    'DtN',norm(top)/max(norm(T(:)),realmin));
end
