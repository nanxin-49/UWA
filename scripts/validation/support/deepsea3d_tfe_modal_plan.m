function p=deepsea3d_tfe_modal_plan(D,D2,gamma)
%DEEPSEA3D_TFE_MODAL_PLAN Diagonalized interior plus exact outgoing Robin.
I=2:size(D,1)-1;[V,Lambda]=eig(D2(I,I));lam=diag(Lambda);
g=gamma(:).';den=lam+g.^2;
p=struct('D',D,'D2',D2,'gamma',g,'V',V,'Vi',V\eye(numel(I)), ...
    'den',den,'interior',I,'basis_rcond',rcond(V));
p.q=V*((V\D2(I,end))./den);
p.robin_den=D(end,end)-1i*g-D(end,I)*p.q;
p.fallback=find(min(abs(den),[],1)./(1+abs(g).^2)<1e-6);
assert(numel(p.fallback)<257,'TFE:ResonanceBudget','Too many near Dirichlet poles.');
end
