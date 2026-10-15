function U=deepsea_bie_gaussian_reference(points,k,aux,sigma,degree)
%DEEPSEA_BIE_GAUSSIAN_REFERENCE Independent radial Fourier-Bessel integral.
% phi=exp(-r^2/(2*sigma^2)), flat trace density; combined-layer field.
[v,w]=deepsea_bie_gauss(degree);rmax=max(k+8/sigma,12/sigma);
edges=[0 k rmax];U=complex(zeros(size(points,1),1));
for j=1:2
    r=(v+1)*(edges(j+1)-edges(j))/2+edges(j);wr=w*(edges(j+1)-edges(j))/2;
    gamma=sqrt(complex(k^2-r.^2));E=exp(2i*gamma*(-aux));
    quotient=-expm1(2i*gamma*(-aux))./gamma;
    symbol=.5*(1+E)+k/2*quotient;
    for first=1:64:size(points,1)
        ids=first:min(first+63,size(points,1));rho=hypot(points(ids,1),points(ids,2));
        B=besselj(0,rho*r.').*exp(1i*points(ids,3)*gamma.');
        U(ids)=U(ids)+sigma^2*B*(wr.*r.*exp(-sigma^2*r.^2/2).*symbol);
    end
end
end
