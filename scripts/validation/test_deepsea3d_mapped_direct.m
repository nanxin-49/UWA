function tests=test_deepsea3d_mapped_direct
tests=functiontests(localfunctions);
end
function setupOnce(~)
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root);setup_vertical_project;
end
function testPlanarActionIsIdentity(t)
p=local_plan(0);x=p.xx;y=p.yy;
v=exp(1i*pi/4*(x+y));W=v(:)*(1+p.s.');
Y=deepsea3d_mapped_action(W,p);
verifyLessThan(t,norm(Y-W,'fro')/norm(W,'fro'),1e-12);
end
function testPreconditionedEqualsRawResidualInverse(t)
p=local_plan(.01);x=p.xx;y=p.yy;
v=exp(1i*pi/4*(x+y));W=v(:)*(1+p.s.'.^2);
[A,B]=deepsea3d_fg_operators(W,p,"galerkin");
F=fft2(reshape(W,p.n,p.n,[]));
H=reshape(ifft2((p.k^2-p.kv.^2-p.kv.'.^2).*F),size(W));
R=H+W*p.D2.'+A+B;
T=ifft2(1i*p.gamma.*fft2(reshape(W(:,end),p.n,p.n)));
top=reshape(W*p.D(end,:).',p.n,p.n)-T ...
    -deepsea3d_fg_top_product(reshape(W(:,end),p.n,p.n),p,"galerkin");
Z=deepsea3d_mapped_l0_solve(R,reshape(W(:,1),p.n,p.n),top,p);
Y=deepsea3d_mapped_action(W,p);
verifyLessThan(t,norm(Z-Y,'fro')/norm(Y,'fro'),1e-9);
end
function testMatrixFreeGMRESAgainstKnownVector(t)
p=local_plan(.01);v=exp(1i*pi/4*(p.xx+p.yy));W=v(:)*sin(p.s.'+1);
rhs=deepsea3d_mapped_action(W,p);shape=size(W);
action=@(x)reshape(deepsea3d_mapped_action(reshape(x,shape),p),[],1);
[x,flag,relres]=gmres(action,rhs(:),8,1e-11,8);
verifyEqual(t,flag,0);verifyLessThan(t,relres,1e-10);
verifyLessThan(t,norm(x-W(:))/norm(W(:)),1e-9);
end
function testCallBudgetPreservesValidIterate(t)
cfg=struct('frequency_hz',600,'c0_mps',1500,'tx_m',[0 150.5 1500], ...
    'rx_m',[0 0 5],'center_m',[0 0],'profile',struct('height_m',0, ...
    'plateau_m',1,'support_m',2,'amplitudes_m',.01,'wavevectors_radpm',[pi/4 pi/4],'phases_rad',.2));
opts=struct('domain_halfwidth_m',4,'spacing_m',1,'vertical_degree',12,'top_depth_m',.8, ...
    'restart',2,'max_cycles',1,'krylov_tolerance',1e-10,'krylov_wall_limit_s',10, ...
    'max_operator_calls',0,'max_workspace_MiB',900);
[o,W]=solve_deepsea3d_mapped_direct(cfg,opts);
verifyEqual(t,o.status,"Mapped:CallBudget");
verifyTrue(t,all(isfinite(W(:))));verifyFalse(t,o.linear_pass);
end
function p=local_plan(A)
cfg=struct('frequency_hz',600,'c0_mps',1500,'tx_m',[0 150.5 1500], ...
    'rx_m',[0 0 5],'center_m',[0 0], ...
    'profile',struct('height_m',0,'plateau_m',100,'support_m',101, ...
    'amplitudes_m',A,'wavevectors_radpm',[pi/4 pi/4],'phases_rad',.2));
opts=struct('domain_halfwidth_m',4,'spacing_m',1,'vertical_degree',12,'top_depth_m',.8);
p=deepsea3d_fg_plan(cfg,opts);
end
