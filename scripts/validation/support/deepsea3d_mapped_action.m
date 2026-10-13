function Y=deepsea3d_mapped_action(W,p)
%DEEPSEA3D_MAPPED_ACTION L0-preconditioned t=1 system, matrix-free.
% W=S0(-[L1+L2]W, g, -eta/b*T W) iff the discrete mapped PDE/BC hold.
[A,B]=deepsea3d_fg_operators(W,p,"galerkin");
top=deepsea3d_fg_top_product(reshape(W(:,end),p.n,p.n),p,"galerkin");
Y=W-deepsea3d_mapped_l0_solve(-A-B,zeros(p.n),top,p);
end
