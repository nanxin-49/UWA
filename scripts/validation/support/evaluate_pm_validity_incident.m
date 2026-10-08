function [u,ux,uz]=evaluate_pm_validity_incident(source,x,z)
%EVALUATE_PM_VALIDITY_INCIDENT Exact frozen discrete Gaussian source.
% Physical partial derivatives, not derivatives along the rough curve.
x=x(:); z=z(:)+zeros(size(x)); n=numel(x);
u=complex(zeros(n,1)); ux=u; uz=u;
for first=1:128:n
 jj=first:min(first+127,n);
 p=exp(1i*(x(jj)-source.x0_m)*source.kx.').*exp(-1i*(z(jj)-source.z_tx_m).*source.kz.');
 u(jj)=p*source.coeff/source.n_original;
 ux(jj)=p*(1i*source.kx.*source.coeff)/source.n_original;
 uz(jj)=p*(-1i*source.kz.*source.coeff)/source.n_original;
end
end
