function out=observe_deepsea3d_pm_trace(cfg,opts,reference)
%OBSERVE_DEEPSEA3D_PM_TRACE Replay frozen FE, exact finite-mode trace samples.
% Diagnostic only: direct exp evaluation removes M40 Taylor-evaluation bias.
timer=tic;
L=opts.domain_halfwidth_m;n=2*ceil(L/opts.spacing_m);h=2*L/n;N=opts.order;
x=-L+(0:n-1)*h;[xx,yy]=meshgrid(x,x);
eta=deepsea3d_surface(xx,yy,cfg.profile);k=2*pi*cfg.frequency_hz/cfg.c0_mps;
kv=2*pi/(2*L)*[0:n/2-1 -n/2:-1];
kz=sqrt(complex(k^2-kv.^2-kv.'.^2));dz=1i*kz;
pts=[xx(:)+cfg.center_m(1),yy(:)+cfg.center_m(2),zeros(n^2,1)];
R0=reshape(vecnorm(pts-cfg.tx_m,2,2),n,n);R=cell(N+1,1);invR=R;E=R;forcing=R;
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
clear R invR E pts R0;
offset=cfg.rx_m(1:2)-cfg.center_m;
kernel=exp(1i*kz*cfg.rx_m(3)).*exp(1i*kv*(offset(1)+L)).*exp(1i*kv.'*(offset(2)+L));
[B,terms,bnorm]=deepsea3d_fe_recursion(eta,dz,forcing,kernel);
term_error=norm(terms-reference.terms)/norm(reference.terms);
norm_error=norm(bnorm-reference.boundary_coefficient_norms)/norm(reference.boundary_coefficient_norms);
assert(term_error<1e-12 && norm_error<1e-12,'pmTrace:Replay','Frozen replay mismatch.');
totalB=complex(zeros(n,n));for m=1:N,totalB=totalB+B{m};end
clear B forcing;
mask=hypot(xx,yy)<cfg.profile.support_m;
ind=find(mask);[~,lo]=min(eta(ind));[~,hi]=max(eta(ind));
extreme=[xx(ind(lo)) yy(ind(lo));xx(ind(hi)) yy(ind(hi))];
% Explicit fixed stencil, plus same-grid extrema: sparse diagnostics, not a
% whole-boundary convergence certificate. Shift stencil both directions.
xy=[0 0;2 0;-2 0;0 2;0 -2;4 4;-4 4;4 -4;-4 -4;8 0;-8 0;0 8;0 -8]+h/2;
xy=[xy;extreme]; z=deepsea3d_surface(xy(:,1),xy(:,2),cfg.profile);
physical=[xy+cfg.center_m,z];image=cfg.tx_m;image(3)=-image(3);
ui=deepsea3d_point_field(physical,cfg.tx_m,k);
flat=ui-deepsea3d_point_field(physical,image,k);
p=complex(zeros(size(z)));p32=p;p40=p;condition=zeros(size(z));
for j=1:numel(z)
    phase=exp(1i*kv*(xy(j,1)+L)).*exp(1i*kv.'*(xy(j,2)+L));
    modal=totalB.*phase/n^2; argument=dz*z(j);
    contribution=modal.*exp(argument);p(j)=sum(contribution,'all');
    condition(j)=sum(abs(contribution),'all')/max(abs(p(j)),realmin);
    poly=ones(n);power=ones(n);
    for m=1:40
        power=power.*argument/m;poly=poly+power;
        if m==32,p32(j)=sum(modal.*poly,'all');end
    end
    p40(j)=sum(modal.*poly,'all');
end
den=norm(flat);
out=struct('xy_m',xy,'eta_m',z,'flat_native',flat,'incident_native',ui, ...
    'exact_trace_native',p,'M32_trace_native',p32,'M40_trace_native',p40, ...
    'exact_trace_residual_relative',norm(flat+p)/den, ...
    'M32_trace_residual_relative',norm(flat+p32)/den, ...
    'M40_trace_residual_relative',norm(flat+p40)/den, ...
    'M40_to_exact_relative',norm(p40-p)/den, ...
    'modal_cancellation_indicators',condition,'term_match_error',term_error, ...
    'norm_match_error',norm_error,'wall_s',toc(timer), ...
    'scope','15 true-boundary points, finite modal exact exponential; no PDE remainder/roundoff bound');
end
