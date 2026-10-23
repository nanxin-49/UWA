function [eta,d1,d2] = closeout_s2c_b_surface(x,A,K,phi)
%CLOSEOUT_S2C_B_SURFACE Frozen C2 physical-depth sinusoid for K_HIGH.
arguments
    x double
    A (1,1) double
    K (1,1) double
    phi (1,1) double
end
r=abs(x); chi=double(r<=42); cp=zeros(size(x)); cpp=cp;
mid=r>42 & r<50; t=(r(mid)-42)/8;
chi(mid)=1-10*t.^3+15*t.^4-6*t.^5;
cp(mid)=(-30*t.^2+60*t.^3-30*t.^4).*sign(x(mid))/8;
cpp(mid)=(-60*t+180*t.^2-120*t.^3)/64;
s=sin(K*x+phi); c=cos(K*x+phi);
eta=A*s.*chi;
d1=A*(K*c.*chi+s.*cp);
d2=A*(-K^2*s.*chi+2*K*c.*cp+s.*cpp);
end
