function [L1,L2]=deepsea3d_fg_operators(W,p,kind)
%DEEPSEA3D_FG_OPERATORS Unprojected three-factor products, then one projection.
% 2n padding is exact for retained products of three symmetric bandlimited
% factors. Multiplying projected quadratic subproducts would change L2.
n=p.n;M=numel(p.s);Ws=W*p.D.';Wss=W*p.D2.';
L1=complex(zeros(size(W)));L2=L1;
if string(kind)=="galerkin",m=2*n;g=p.pad_geo;else,m=n;g=p.geo;end
e=g.eta/p.b;
for j=1:M
    F=fft2(reshape(W(:,j),n,n));Fs=fft2(reshape(Ws(:,j),n,n));
    H=ifft2((p.k^2-p.kv.^2-p.kv.'.^2).*F);
    X=ifft2(1i*p.kv.*Fs);Y=ifft2(1i*p.kv.'.*Fs);
    S=reshape(Ws(:,j),n,n);SS=reshape(Wss(:,j),n,n);
    if m>n
        H=deepsea3d_fg_resize(H,m);X=deepsea3d_fg_resize(X,m);Y=deepsea3d_fg_resize(Y,m);
        S=deepsea3d_fg_resize(S,m);SS=deepsea3d_fg_resize(SS,m);
    end
    r=1-p.s(j)/p.b;mix=g.ex.*X+g.ey.*Y;
    A=-2*e.*H-2*r*mix-g.lap*r.*S;
    B=e.^2.*H+g.slope2*r^2.*SS+2*e*r.*mix+(e.*g.lap-2*g.slope2/p.b)*r.*S;
    % Both branches must inhabit the same symmetric subspace. Even on the
    % base grid, discard newly generated Nyquist lines before the ODE solve.
    A=deepsea3d_fg_resize(A,n);B=deepsea3d_fg_resize(B,n);
    L1(:,j)=A(:);L2(:,j)=B(:);
end
end
