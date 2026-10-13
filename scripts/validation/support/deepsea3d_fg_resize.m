function V=deepsea3d_fg_resize(U,m)
%DEEPSEA3D_FG_RESIZE Fourier interpolation/projection with explicit scaling.
% Even arrays retain symmetric modes |j|<n/2; both Nyquist lines are zero.
n=size(U,1);assert(size(U,2)==n && mod(n,2)==0 && mod(m,2)==0);
F=fftshift(fft2(U));F(1,:)=0;F(:,1)=0;
if m>=n
    G=complex(zeros(m));j=(m-n)/2+(1:n);G(j,j)=F;
else
    j=(n-m)/2+(1:m);G=F(j,j);G(1,:)=0;G(:,1)=0;
end
V=ifft2(ifftshift(G))*(m/n)^2;
end
