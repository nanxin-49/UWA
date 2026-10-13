function W=deepsea3d_tfe_modal_solve(p,F,bottom,top)
%DEEPSEA3D_TFE_MODAL_SOLVE Wss+gamma^2 W=F; W(0)=g; Ws(b)-i gamma W(b)=h.
% Columns are horizontal modes, rows are vertical nodes. A few near-pole
% columns use the full coupled system, avoiding Dirichlet Schur cancellation.
I=p.interior;g=bottom(:).';h=top(:).';
v=p.V*((p.Vi*(F(I,:)-p.D2(I,1)*g))./p.den);
t=(h-p.D(end,I)*v-p.D(end,1)*g)./p.robin_den;
W=[g;v-p.q.*t;t];
for j=p.fallback
    A=p.D2+p.gamma(j)^2*eye(size(p.D));A(1,:)=0;A(1,1)=1;
    A(end,:)=p.D(end,:);A(end,end)=A(end,end)-1i*p.gamma(j);
    rhs=F(:,j);rhs(1)=g(j);rhs(end)=h(j);W(:,j)=A\rhs;
end
end
