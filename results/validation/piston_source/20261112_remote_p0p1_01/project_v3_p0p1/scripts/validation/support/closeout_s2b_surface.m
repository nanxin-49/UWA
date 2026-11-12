function [eta,ep]=closeout_s2b_surface(x,A)
% Frozen physical depth C2 sinusoid; no realization or support change.
r=abs(x); chi=zeros(size(x)); dc=chi; chi(r<=42)=1; mid=r>42&r<50;
t=(r(mid)-42)/8; chi(mid)=1-10*t.^3+15*t.^4-6*t.^5;
dc(mid)=(-30*t.^2+60*t.^3-30*t.^4)/8;
eta=A*sin(.1*x).*chi;
ep=A*(.1*cos(.1*x).*chi+sin(.1*x).*dc.*sign(x));
end
