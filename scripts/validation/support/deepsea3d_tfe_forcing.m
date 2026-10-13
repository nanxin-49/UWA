function forcing=deepsea3d_tfe_forcing(cfg,x,y,eta,N)
%DEEPSEA3D_TFE_FORCING Analytic t coefficients of negative flat total field.
% Point source is outside the strip; no singular volume forcing is expanded.
k=2*pi*cfg.frequency_hz/cfg.c0_mps;
R0=sqrt((x+cfg.center_m(1)-cfg.tx_m(1)).^2+(y+cfg.center_m(2)-cfg.tx_m(2)).^2+cfg.tx_m(3)^2);
R=cell(N+1,1);I=R;E=R;forcing=R;
R{1}=R0;I{1}=1./R0;E{1}=exp(1i*k*R0);forcing{1}=zeros(size(x));
for m=1:N
    a=zeros(size(x));
    if m==1,a=-2*cfg.tx_m(3)*eta;elseif m==2,a=eta.^2;end
    for j=1:m-1,a=a-R{j+1}.*R{m-j+1};end
    R{m+1}=a./(2*R0);a=zeros(size(x));b=complex(a);
    for j=1:m
        a=a+R{j+1}.*I{m-j+1};b=b+j*R{j+1}.*E{m-j+1};
    end
    I{m+1}=-a./R0;E{m+1}=1i*k*b/m;g=complex(zeros(size(x)));
    for j=0:m,g=g+E{j+1}.*I{m-j+1};end
    if mod(m,2)==1,forcing{m+1}=-2*g/(4*pi);else,forcing{m+1}=zeros(size(x));end
end
end
