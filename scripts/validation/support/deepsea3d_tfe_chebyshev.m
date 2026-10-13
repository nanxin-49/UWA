function [s,D,D2]=deepsea3d_tfe_chebyshev(M,b)
%DEEPSEA3D_TFE_CHEBYSHEV Ascending physical strip coordinates, base MATLAB.
j=(0:M)';x=cos(pi*j/M);c=[2;ones(M-1,1);2].*(-1).^j;
d=x-x.';D=(c*(1./c).')./(d+eye(M+1));D=D-diag(sum(D,2));
D=-2/b*D;D2=D*D;s=b/2*(1-x);
end
