function out=observe_deepsea3d_fe_cancellation(cfg,opts,reference)
%OBSERVE_DEEPSEA3D_FE_CANCELLATION Diagnostic replay, not a new oracle.
% Reconstruct the frozen recurrence to measure contribution norms. Assert
% its pressure coefficients and b_n norms match the actual frozen solve.
timer=tic;cpu0=cputime;
L=opts.domain_halfwidth_m; n=2*ceil(L/opts.spacing_m); h=2*L/n;
x=-L+(0:n-1)*h;[xx,yy]=meshgrid(x,x);eta=deepsea3d_surface(xx,yy,cfg.profile);
k=2*pi*cfg.frequency_hz/cfg.c0_mps;
kv=2*pi/(2*L)*[0:n/2-1 -n/2:-1];kz=sqrt(complex(k^2-kv.^2-kv.'.^2));dz=1i*kz;
points=[xx(:)+cfg.center_m(1),yy(:)+cfg.center_m(2),zeros(n^2,1)];
R0=reshape(vecnorm(points-cfg.tx_m,2,2),n,n);N=opts.order;
R=cell(N+1,1);invR=R;E=R;forcing=R;
R{1}=R0;invR{1}=1./R0;E{1}=exp(1i*k*R0);forcing{1}=zeros(n,n);
for m=1:N
    a=zeros(n,n);
    if m==1,a=-2*cfg.tx_m(3)*eta;elseif m==2,a=eta.^2;end
    for j=1:m-1,a=a-R{j+1}.*R{m-j+1};end
    R{m+1}=a./(2*R0);a=zeros(n,n);b=complex(a);
    for j=1:m
        a=a+R{j+1}.*invR{m-j+1};b=b+j*R{j+1}.*E{m-j+1};
    end
    invR{m+1}=-a./R0;E{m+1}=1i*k*b/m;g=complex(zeros(n,n));
    for j=0:m,g=g+E{j+1}.*invR{m-j+1};end
    if mod(m,2)==1,forcing{m+1}=2*g/(4*pi);else,forcing{m+1}=zeros(n,n);end
end
clear R invR E R0 points xx yy;
offset=cfg.rx_m(1:2)-cfg.center_m;
kernel=exp(1i*kz*cfg.rx_m(3)).*exp(1i*kv*(offset(1)+L)).*exp(1i*kv.'*(offset(2)+L));
B=cell(N,1);terms=complex(zeros(N,1));bnorm=zeros(N,1);conditions=zeros(N,1);numerators=zeros(N,1);
for m=1:N
    b=-forcing{m+1};amount=norm(b(:));power=ones(n,n);fac=1;
    for j=1:m-1
        power=power.*eta;fac=fac*j; contribution=power/fac.*ifft2(dz.^j.*B{m-j});
        amount=amount+norm(contribution(:)); b=b-contribution;
    end
    bnorm(m)=norm(b(:));numerators(m)=amount;conditions(m)=amount/max(bnorm(m),realmin);
    B{m}=fft2(b);terms(m)=sum(B{m}.*kernel,'all')/n^2;
end
term_error=norm(terms-reference.terms)/norm(reference.terms);
norm_error=norm(bnorm-reference.boundary_coefficient_norms)/norm(reference.boundary_coefficient_norms);
assert(term_error<1e-12 && norm_error<1e-12,'alpha25:ReplayMismatch','Diagnostic differs from frozen FE.');
out=struct('terms',terms,'bnorm',bnorm,'conditions',conditions,'numerators',numerators, ...
    'term_match_error',term_error,'norm_match_error',norm_error,'wall_s',toc(timer),'cpu_s',cputime-cpu0);
end
