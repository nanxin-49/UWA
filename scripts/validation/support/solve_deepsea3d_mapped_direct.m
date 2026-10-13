function [output,W]=solve_deepsea3d_mapped_direct(cfg,opts)
%SOLVE_DEEPSEA3D_MAPPED_DIRECT Restarted GMRES, exact source Dirichlet at t=1.
% Same Galerkin coefficients and outgoing DtN as the verified TFE operators.
timer=tic;cpu0=cputime;
n=2*ceil(opts.domain_halfwidth_m/opts.spacing_m);M=opts.vertical_degree+1;
estimate=(opts.restart+15)*16*n^2*M/2^20+18*16*(4*n)^2/2^20;
[before,sys]=memory;
assert(estimate<=opts.max_workspace_MiB && estimate<.8*sys.PhysicalMemory.Available/2^20, ...
    'Mapped:MemoryBudget','Insufficient bounded physical memory; do not rely on paging.');
p=deepsea3d_fg_plan(cfg,opts);
assert(max(p.original_geo.eta(:))<p.b && p.L>cfg.profile.support_m);
assert(cfg.profile.height_m==0 && p.b>0 && p.b<cfg.rx_m(3) && p.b<cfg.tx_m(3));
points=[p.xx(:)+cfg.center_m(1),p.yy(:)+cfg.center_m(2),p.original_geo.eta(:)];
image=cfg.tx_m;image(3)=-image(3);
flat=deepsea3d_point_field(points,cfg.tx_m,p.k)-deepsea3d_point_field(points,image,p.k);
g=deepsea3d_fg_resize(-reshape(flat,n,n),n);
rhs=deepsea3d_mapped_l0_solve(complex(zeros(n^2,M)),g,zeros(n),p);
rhs=rhs(:);bnorm=norm(rhs);x=rhs;calls=0;limit_active=true;krylov_clock=tic;
cycles=struct([]);allres=[];status="MAX_CYCLES";
if bnorm==0
    status="ZERO_CORRECTION";x=zeros(size(rhs),'like',rhs);
else
    for cycle=1:opts.max_cycles
        previous=x;
        try
            [candidate,flag,relres,iter,resvec]=gmres(@apply, rhs,opts.restart,opts.krylov_tolerance,1,[],[],x);
            x=candidate;
        catch err
            budget_id=local_budget_identifier(err);
            if budget_id~=""
                x=previous;status=budget_id;break
            else
                rethrow(err);
            end
        end
        row=struct('cycle',cycle,'flag',flag,'relative_residual',relres,'outer',iter(1), ...
            'inner',iter(2),'wall_s',toc(krylov_clock),'operator_calls',calls);
        if isempty(cycles),cycles=row;else,cycles(end+1)=row;end %#ok<AGROW>
        allres=[allres;resvec(:)]; %#ok<AGROW>
        fprintf('Mapped cycle %d relres %.6g flag%d calls%d wall %.2fs\n',cycle,relres,flag,calls,toc(krylov_clock));
        if relres<opts.krylov_tolerance && flag==0,status="KRYLOV_CONVERGED";break;end
        if toc(krylov_clock)>=opts.krylov_wall_limit_s,status="Mapped:TimeBudget";break;end
    end
end
function id=local_budget_identifier(err)
% MATLAB iterapp wraps callback failures; inspect the structured cause chain.
id="";
if any(string(err.identifier)==["Mapped:TimeBudget","Mapped:CallBudget"])
    id=string(err.identifier);return
end
for j=1:numel(err.cause)
    id=local_budget_identifier(err.cause{j});
    if id~="",return;end
end
end
krylov_wall=toc(krylov_clock);limit_active=false;
fresh_residual=norm(apply(x)-rhs)/max(bnorm,realmin);
W=reshape(x,n^2,M);discrete=deepsea3d_mapped_discrete_residual(W,g,p);
if bnorm==0
    r3=struct('physical_PDE',0,'true_boundary',0,'true_DtN',0,'wall_s',0);
    r4=r3;
else
    r3=observe_deepsea3d_fg_residual(W,p,cfg,3);
    r4=observe_deepsea3d_fg_residual(W,p,cfg,4);
end
offset=cfg.rx_m(1:2)-cfg.center_m;
kernel=exp(1i*p.gamma*(cfg.rx_m(3)-p.b)).*exp(1i*p.kv*(offset(1)+p.L)).*exp(1i*p.kv.'*(offset(2)+p.L));
delta=sum(fft2(reshape(W(:,end),n,n)).*kernel,'all')/n^2;bg=deepsea3d_flat(cfg,0);
[after,~]=memory;
linear_pass=fresh_residual<1e-8 && discrete.Dirichlet<1e-8 ...
    && discrete.DtN<1e-8 && discrete.mapped_PDE<1e-8;
physical_pass=all([r3.physical_PDE,r3.true_boundary,r3.true_DtN,r4.physical_PDE,r4.true_boundary,r4.true_DtN]<1e-4);
output=struct('status',status,'cycles',cycles,'residual_history',allres, ...
    'fresh_preconditioned_residual',fresh_residual,'discrete',discrete,'independent3',r3,'independent4',r4, ...
    'linear_pass',linear_pass,'initial_physical_pass',physical_pass, ...
    'rough_correction',delta,'reflect_native',bg.reflect_native+delta,'total_native',bg.total_native+delta, ...
    'direct_native',bg.direct_native,'options',opts,'operator_calls',calls, ...
    'krylov_wall_s',krylov_wall,'wall_s',toc(timer),'cpu_s',cputime-cpu0, ...
    'estimated_workspace_MiB',estimate,'vector_MiB',16*n^2*M/2^20, ...
    'process_before_MiB',before.MemUsedMATLAB/2^20,'process_after_MiB',after.MemUsedMATLAB/2^20, ...
    'time_convention','exp(-i omega t)','source','exact point/image on true graph; no Taylor coefficients', ...
    'PM_status','INCONCLUSIVE','global_G1','OPEN','WP2_allowed',false);
fprintf('Mapped %s fresh %.3g discrete PDE %.3g true PDE %.6g boundary %.6g DtN %.6g wall %.3fs\n', ...
    status,fresh_residual,discrete.mapped_PDE,r4.physical_PDE,r4.true_boundary,r4.true_DtN,output.wall_s);
    function y=apply(v)
        if limit_active
            if toc(krylov_clock)>opts.krylov_wall_limit_s,error('Mapped:TimeBudget','Krylov wall limit.');end
            if calls>=opts.max_operator_calls,error('Mapped:CallBudget','Matrix-vector call limit.');end
        end
        calls=calls+1;V=reshape(v,n^2,M);Y=deepsea3d_mapped_action(V,p);y=Y(:);
    end
end
