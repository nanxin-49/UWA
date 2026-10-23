function [u,ux,uz] = fk_single_receiver_incident(src,x,z)
%FK_SINGLE_RECEIVER_INCIDENT Exact frozen angular spectrum, physical gradients.
% Same block summation and normalization as the existing FullFK timing audit.
x=x(:); z=z(:)+zeros(size(x)); u=complex(zeros(size(x))); ux=u; uz=u;
for first=1:128:numel(x)
    jj=first:min(first+127,numel(x));
    phase=exp(1i*(x(jj)-src.x0_m)*src.kx.').*exp(-1i*(z(jj)-src.z_tx_m).*src.kz.');
    u(jj)=phase*src.coeff/src.n_original;
    ux(jj)=phase*(1i*src.kx.*src.coeff)/src.n_original;
    uz(jj)=phase*(-1i*src.kz.*src.coeff)/src.n_original;
end
end
