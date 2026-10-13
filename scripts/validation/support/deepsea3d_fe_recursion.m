function [B,terms,bnorm]=deepsea3d_fe_recursion(eta,dz,forcing,rx_kernel)
%DEEPSEA3D_FE_RECURSION Dirichlet recurrence, independent of incident choice.
N=numel(forcing)-1; n2=numel(eta);
B=cell(N,1); terms=complex(zeros(N,1)); bnorm=zeros(N,1);
for m=1:N
    b=-forcing{m+1}; power=ones(size(eta)); factorial=1;
    for j=1:m-1
        power=power.*eta; factorial=factorial*j;
        b=b-power/factorial.*ifft2(dz.^j.*B{m-j});
    end
    bnorm(m)=norm(b(:)); B{m}=fft2(b);
    terms(m)=sum(B{m}.*rx_kernel,'all')/n2;
end
end
