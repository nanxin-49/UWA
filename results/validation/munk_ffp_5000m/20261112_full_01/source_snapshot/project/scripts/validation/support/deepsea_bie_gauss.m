function [x,w]=deepsea_bie_gauss(n)
%DEEPSEA_BIE_GAUSS Base-MATLAB Gauss-Legendre on [-1,1].
b=(1:n-1)./sqrt(4*(1:n-1).^2-1);[V,D]=eig(diag(b,1)+diag(b,-1));
[x,i]=sort(diag(D));w=2*V(1,i).'.^2;
end
