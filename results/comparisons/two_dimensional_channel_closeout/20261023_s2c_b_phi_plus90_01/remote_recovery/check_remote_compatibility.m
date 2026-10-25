% Remote-only read/evaluation gate. Saved Windows densities are never changed.
root='/home/misc/cwx/UWA_phi_plus90';
addpath(fullfile(root,'scripts','validation','support'));
run_dir=fullfile(root,'results','comparisons','two_dimensional_channel_closeout','20261023_s2c_b_phi_plus90_01');
assert(~isfile(fullfile(run_dir,'jobs','bie_space14.mat')));
report=struct;
report.release=version;
report.matlab_release=version('-release');
report.started_utc=char(datetime('now','TimeZone','UTC'));
report.solver_sha256=closeout_s2b_hash(fullfile(root,'scripts','validation','support','solve_helmholtz_bie_halfplane_vertical.m'));
assert(strcmp(report.solver_sha256,'76f1a9d52c898d62507ad7d971813c79b1f09cd9ced4f6dbce556ed11c349fa9'));
src_path=fullfile(root,'results','comparisons','two_dimensional_channel_closeout','20261023_s2c_a_1642','protocol','frozen_source.mat');
assert(strcmp(closeout_s2b_hash(src_path),'8d60d3d9d6aa3aa16110ac33d81e0d193a61526aef2c0c2e9edd531cae531f55'));
t=load(src_path,'src'); assert(numel(t.src.coeff)==984);
x0=(-192.1875/2)+(0:983)'*(192.1875/984);
assert(max(abs(ifft(t.src.coeff(:))-exp(-.5*(x0/.3).^2)))<1e-12);
[eta,deta,d2eta]=closeout_s2c_b_surface((-60:.01:60)',.05,.1,pi/2);
assert(max(abs(eta(abs((-60:.01:60)')>=50)))==0);
assert(max(abs(deta(abs((-60:.01:60)')>=50)))==0);
assert(max(abs(d2eta(abs((-60:.01:60)')>=50)))==0);
assert(abs(closeout_s2c_b_surface(0,.05,.1,pi/2)-.05)<1e-14);
assert(abs(besselh(0,1,1)-(0.765197686557967+0.088256964215677i))<1e-13);
ids={'main96','self64'};
report.levels=[];
for j=1:2
    path=fullfile(run_dir,'jobs',['bie_' ids{j} '.mat']);
    q=load(path,'boundary','receipt'); b=q.boundary; v=q.receipt;
    assert(strcmp(v.status,'COMPLETE') && b.finite && b.node_count==3840);
    assert(numel(b.density)==3840 && numel(b.receiver_field)==417);
    assert(numel(b.receiver_x_m)==417 && numel(b.receiver_z_m)==417);
    assert(strcmp(v.source_sha256,closeout_s2b_hash(src_path)));
    assert(strcmp(v.solver_sha256,report.solver_sha256));
    [z,sl]=closeout_s2c_b_surface(b.node_x_m(:),.05,.1,pi/2);
    assert(max(abs(z-b.node_z_m(:)))<1e-12);
    n=[-sl./hypot(1,sl),1./hypot(1,sl)];
    assert(max(abs(n-b.node_normal),[],'all')<1e-12);
    tic_eval=tic;
    u=remote_saved_density_evaluation(b);
    seconds=toc(tic_eval);
    rel=norm(u-b.receiver_field(:))/norm(b.receiver_field(:));
    max_abs=max(abs(u-b.receiver_field(:)));
    level=struct('id',ids{j},'receiver_count',417,'density_count',3840,...
        'relative_cross_platform_field_difference',rel,'max_abs_cross_platform_field_difference',max_abs,...
        'evaluation_seconds',seconds,'saved_wall_seconds',v.wall_seconds);
    if j==1, report.levels=level; else, report.levels(j)=level; end
end
report.cross_platform_pass=all([report.levels.relative_cross_platform_field_difference]<1e-10);
report.completed_utc=char(datetime('now','TimeZone','UTC'));
f=fullfile(run_dir,'remote_recovery_compatibility.json');
assert(~isfile(f)); fid=fopen(f,'w','n','UTF-8'); assert(fid>0);
fprintf(fid,'%s',jsonencode(report,PrettyPrint=true)); fclose(fid);
disp(report);
assert(report.cross_platform_pass,'Saved-density field differs across platforms');

function u=remote_saved_density_evaluation(b)
    n=b.node_count; np=b.panel_count; po=b.panel_order;
    edges=linspace(-b.half_width_m,b.half_width_m,np+1)';
    j=(1:po-1)'; beta=j./sqrt(4*j.^2-1);
    [V,D]=eig(diag(beta,1)+diag(beta,-1));
    [gx,ord]=sort(diag(D)); V=V(:,ord); gw=2*(V(1,:).^2)';
    xp=reshape((.5*(edges(1:end-1)+edges(2:end))+.5*diff(edges).*gx.').',[],1);
    assert(numel(xp)==n && max(abs(xp-b.node_x_m(:)))<1e-12);
    w=repmat(.5*diff(edges),1,po).'; w=w(:).*repmat(gw,np,1);
    jac=1./b.node_normal(:,2);
    weights=w.*jac.*b.window(:);
    src=weights.*b.density(:);
    u=complex(zeros(417,1)); k=b.k_radpm; h=b.halfplane_h_m;
    xs=b.node_x_m(:).'; zs=b.node_z_m(:).';
    nxs=b.node_normal(:,1).'; nzs=b.node_normal(:,2).';
    for first=1:64:417
        ii=first:min(first+63,417);
        dx=b.receiver_x_m(ii)-xs; dz=b.receiver_z_m(ii)-zs;
        rr=hypot(dx,dz); assert(all(rr>0,'all'));
        phi=1i/4*besselh(0,1,k*rr);
        dphi=1i*k/4*besselh(1,1,k*rr).*((dx./rr).*nxs+(dz./rr).*nzs);
        dzi=b.receiver_z_m(ii)-(2*h-zs); ri=hypot(dx,dzi);
        phii=1i/4*besselh(0,1,k*ri);
        dphii=1i*k/4*besselh(1,1,k*ri).*((dx./ri).*nxs-(dzi./ri).*nzs);
        u(ii)=((dphi-dphii)-1i*k*(phi-phii))*src;
    end
end
