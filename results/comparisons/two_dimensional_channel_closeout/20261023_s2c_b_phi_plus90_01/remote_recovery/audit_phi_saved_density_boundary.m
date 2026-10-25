% Linux R2023b, read-only boundary evaluation of saved Windows BIE densities.
% Reconstructs the off-node Dirichlet residual; never solves for a density.
root='/home/misc/cwx/UWA_phi_plus90';
addpath(fullfile(root,'scripts','validation','support'));
run_dir=fullfile(root,'results','comparisons','two_dimensional_channel_closeout','20261023_s2c_b_phi_plus90_01');
out=fullfile(run_dir,'audit','remote_saved_density_boundary_20261025.json');
assert(~isfile(out));
src_file=fullfile(root,'results','comparisons','two_dimensional_channel_closeout','20261023_s2c_a_1642','protocol','frozen_source.mat');
t=load(src_file,'src'); src=t.src;
ids={'main96','self64'}; q=cell(2,1);
for j=1:2
    p=fullfile(run_dir,'jobs',['bie_' ids{j} '.mat']);
    d=load(p,'boundary','receipt'); b=d.boundary;
    assert(strcmp(d.receipt.status,'COMPLETE') && b.node_count==3840);
    assert(strcmp(d.receipt.source_sha256,closeout_s2b_hash(src_file)));
    tic_eval=tic;
    y=boundary_residual_again(b,src);
    seconds=toc(tic_eval);
    saved=b.boundary_total_residual(:);
    assert(numel(y)==numel(saved));
    rel=norm(y-saved)/max(norm(saved),realmin);
    maxerr=max(abs(y-saved));
    rho=norm(y)/norm(fk_single_receiver_incident(src,b.boundary_check_x_m(:), ...
        closeout_s2c_b_surface(b.boundary_check_x_m(:),.05,.1,pi/2)));
    q{j}=struct('job_id',ids{j},'points',numel(y), ...
        'linux_residual_relative_to_incident',rho, ...
        'windows_saved_residual_relative_to_incident',b.offgrid_boundary_residual, ...
        'relative_difference_to_saved_residual_vector',rel, ...
        'max_abs_difference_to_saved_residual_vector',maxerr, ...
        'evaluation_seconds',seconds);
end
result=struct('status','SAVED_WINDOWS_DENSITY_BOUNDARY_EVALUATED_ON_LINUX', ...
    'release',version,'levels',{q},'new_boundary_solves',0, ...
    'source_sha256',closeout_s2b_hash(src_file), ...
    'solver_sha256',closeout_s2b_hash(fullfile(root,'scripts','validation','support','solve_helmholtz_bie_halfplane_vertical.m')), ...
    'note','This checks the frozen Windows densities under an independently reconstructed Linux off-node operator. It is not a direct same-grid cross-platform density solve.');
fid=fopen(out,'w','n','UTF-8'); assert(fid>=0);
fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true)); fclose(fid);
disp(result);

function residual=boundary_residual_again(b,src)
    po=b.panel_order; np=b.panel_count; qa=b.self_quadrature_order;
    edges=linspace(-b.half_width_m,b.half_width_m,np+1)';
    [gx,~]=gauss(po); bary=barycentric(gx);
    target=b.boundary_check_x_m(:); z=closeout_s2c_b_surface(target,.05,.1,pi/2);
    u=fk_single_receiver_incident(src,target,z);
    residual=complex(zeros(size(u)));
    psi=b.density(:); k=b.k_radpm; h=b.halfplane_h_m;
    B0=lagrange(gx,bary,0);
    for tt=1:numel(target)
        pp=find(target(tt)>=edges(1:end-1) & target(tt)<=edges(2:end),1,'first');
        assert(~isempty(pp));
        layer=0;
        for ps=1:np
            cols=(ps-1)*po+(1:po);
            w=panel_weights(target(tt),z(tt),edges(ps),edges(ps+1), ...
                gx,bary,k,h,b.half_width_m,b.window_plateau_ratio,qa,ps==pp);
            layer=layer+w*psi(cols);
        end
        cols=(pp-1)*po+(1:po);
        residual(tt)=u(tt)+.5*(B0*psi(cols))+layer;
    end
end

function w=panel_weights(xt,zt,a,b,gx,bary,k,h,L,c,qa,singular)
    [ug,uw]=gauss(qa);
    if singular
        v=.5*(ug+1); vw=.5*uw;
        xi=[-v.^2;v.^2];
        qw=[vw.*(2*v);vw.*(2*v)];
    else
        xi=ug; qw=uw;
    end
    half=.5*(b-a); xs=.5*(a+b)+half*xi;
    [zs,slope]=closeout_s2c_b_surface(xs,.05,.1,pi/2);
    jac=hypot(1,slope); nx=-slope./jac; nz=1./jac;
    win=slow_window(xs,L,c);
    dx=xt-xs; dz=zt-zs; r=hypot(dx,dz);
    rs=r; rs(rs==0)=1;
    phi=1i/4*besselh(0,1,k*rs);
    dphi=1i*k/4*besselh(1,1,k*rs).*((dx./rs).*nx+(dz./rs).*nz);
    dzi=zt-(2*h-zs); ri=hypot(dx,dzi);
    phii=1i/4*besselh(0,1,k*ri);
    dphii=1i*k/4*besselh(1,1,k*ri).*((dx./ri).*nx-(dzi./ri).*nz);
    S=phi-phii; K=dphi-dphii;
    S(r==0)=0; K(r==0)=0;
    kernel=K-1i*k*S;
    B=lagrange(gx,bary,xi);
    physical=half*qw.*jac.*win;
    w=(kernel(:).*physical(:)).'*B;
end

function w=slow_window(x,L,c)
    r=abs(x)/L; w=ones(size(r)); w(r>=1)=0;
    mid=r>c & r<1; t=(r(mid)-c)/(1-c);
    a=exp(-1./t); b=exp(-1./(1-t)); w(mid)=b./(a+b);
end

function B=lagrange(nodes,bary,x)
    nodes=nodes(:).'; bary=bary(:).'; x=x(:);
    B=zeros(numel(x),numel(nodes));
    for ii=1:numel(x)
        d=x(ii)-nodes; [dm,jj]=min(abs(d));
        if dm<=32*eps(max(1,abs(x(ii))))
            B(ii,jj)=1;
        else
            row=bary./d; B(ii,:)=row/sum(row);
        end
    end
end

function bary=barycentric(nodes)
    nodes=nodes(:); n=numel(nodes); bary=ones(n,1);
    for jj=1:n
        bary(jj)=1/prod(nodes(jj)-nodes([1:jj-1,jj+1:n]));
    end
end

function [x,w]=gauss(n)
    j=(1:n-1)'; beta=j./sqrt(4*j.^2-1);
    [V,D]=eig(diag(beta,1)+diag(beta,-1));
    [x,ix]=sort(diag(D)); V=V(:,ix); w=2*(V(1,:).^2)';
end
