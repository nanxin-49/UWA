function [U19p5_mps, ustar_mps, z0_m] = li2009_u10_to_u195(U10_mps)
%LI2009_U10_TO_U195 Convert 10-m wind to the PM wind used by Li et al. (2009).
%   The conversion solves Eqs. (12)-(13) of Li et al. (2009):
%     U(z) = u_star/0.4 * log(z/z0), z0*g/u_star^2 = 0.015.

arguments
    U10_mps {mustBeNumeric,mustBeFinite,mustBePositive}
end

g = 9.81;
kappa = 0.4;
charnock = 0.015;
U19p5_mps = zeros(size(U10_mps));
ustar_mps = zeros(size(U10_mps));
z0_m = zeros(size(U10_mps));
for ii = 1:numel(U10_mps)
    U10 = U10_mps(ii);
    residual = @(u) (u/kappa) .* log(10*g/(charnock*u.^2)) - U10;
    ustar = fzero(residual, [1e-4, max(5, U10)]);
    z0 = charnock * ustar^2 / g;
    ustar_mps(ii) = ustar;
    z0_m(ii) = z0;
    U19p5_mps(ii) = (ustar/kappa) * log(19.5/z0);
end
end
