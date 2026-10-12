function [x,w] = deepsea3d_quadrature(L,h,kind)
%DEEPSEA3D_QUADRATURE Tensor trapezoid or independent four-point Gauss panels.
assert(L>0 && h>0,'deepsea3d:Grid','Positive L/h required.');
if string(kind)=="trap"
    n=ceil(2*L/h); x=linspace(-L,L,n+1); w=ones(size(x))*(2*L/n);
    w([1 end])=w([1 end])/2;
elseif string(kind)=="gauss4"
    n=ceil(2*L/(4*h)); edges=linspace(-L,L,n+1);
    g=[-.8611363115940526 -.3399810435848563 .3399810435848563 .8611363115940526];
    gw=[.3478548451374538 .6521451548625461 .6521451548625461 .3478548451374538];
    half=L/n; centers=(edges(1:end-1)+edges(2:end))/2;
    x=reshape((centers(:)+half*g).',1,[]); w=repmat(half*gw,1,n);
else
    error('deepsea3d:Quadrature','Unknown quadrature.');
end
end
