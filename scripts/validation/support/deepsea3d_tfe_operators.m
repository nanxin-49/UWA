function [L0,L1,L2]=deepsea3d_tfe_operators(W,geo,s,D,D2,kv,k,b)
%DEEPSEA3D_TFE_OPERATORS Polynomial a^2 times mapped physical Helmholtz.
% z=s+t*(1-s/b)*eta, a=1-t*eta/b; L(t)=L0+t L1+t^2 L2.
n=numel(kv);M=numel(s);F=fft2(reshape(W,n,n,M));
horizontal=reshape(ifft2((k^2-kv.^2-kv.'.^2).*F),n^2,M);
Ws=W*D.';Wss=W*D2.';Fs=fft2(reshape(Ws,n,n,M));
mixed=geo.ex(:).*reshape(ifft2(1i*kv.*Fs),n^2,M) ...
    +geo.ey(:).*reshape(ifft2(1i*kv.'.*Fs),n^2,M);
e=geo.eta(:)/b;r=(1-s.'/b);
L0=horizontal+Wss;
L1=-2*e.*horizontal-2*r.*mixed-geo.lap(:).*r.*Ws;
L2=e.^2.*horizontal+geo.slope2(:).*r.^2.*Wss ...
    +2*e.*r.*mixed+(e.*geo.lap(:)-2*geo.slope2(:)/b).*r.*Ws;
end
