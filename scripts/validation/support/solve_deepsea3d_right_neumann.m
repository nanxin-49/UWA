function [output,W]=solve_deepsea3d_right_neumann(cfg,opts)
%SOLVE_DEEPSEA3D_RIGHT_NEUMANN One fixed right polynomial, frozen B unchanged.
timer=tic;cpu0=cputime;n=2*ceil(opts.domain_halfwidth_m/opts.spacing_m);M=opts.vertical_degree+1;
estimate=(opts.restart+19)*16*n^2*M/2^20+18*16*(4*n)^2/2^20;
[before,sys]=memory;
assert(estimate<=opts.max_workspace_MiB && estimate<.7*sys.PhysicalMemory.Available/2^20, ...
    'Neumann:MemoryBudget','Require 30 percent physical RAM reserve.');
p=deepsea3d_fg_plan(cfg,opts);
assert(max(p.original_geo.eta(:))<p.b && p.L>cfg.profile.support_m);
points=[p.xx(:)+cfg.center_m(1),p.yy(:)+cfg.center_m(2),p.original_geo.eta(:)];
image=cfg.tx_m;image(3)=-image(3);
flat=deepsea3d_point_field(points,cfg.tx_m,p.k)-deepsea3d_point_field(points,image,p.k);
g=deepsea3d_fg_resize(-reshape(flat,n,n),n);
rhs=deepsea3d_mapped_l0_solve(complex(zeros(n^2,M)),g,zeros(n),p);rhs=rhs(:);
bnorm=norm(rhs);calls=0;active=true;y=rhs;x=rhs;cycles=struct([]);hist=[];
status="MAX_CYCLES";pre_wall=0;pre_calls=0;
if bnorm==0
    x=zeros(size(rhs));status="ZERO_CORRECTION";
else
    for cycle=1:opts.max_cycles
        try
            [candidate,flag,relres,iter,resvec]=gmres(@right,rhs,opts.restart,opts.krylov_tolerance,1,[],[],y);
            candidate_x=pre(candidate);fresh=norm(B(candidate_x)-rhs)/bnorm;
            y=candidate;x=candidate_x;
        catch err
            id=budget_id(err);
            if id~="",status=id;break;else,rethrow(err);end
        end
        row=struct('cycle',cycle,'flag',flag,'right_system_relres',relres, ...
            'fresh_B_residual',fresh,'inner',iter(2),'operator_calls',calls,'wall_s',toc(timer));
        if isempty(cycles),cycles=row;else,cycles(end+1)=row;end %#ok<AGROW>
        hist=[hist;resvec(:)]; %#ok<AGROW>
        fprintf('Right cycle%d B %.6g flag%d calls%d wall%.2fs\n',cycle,fresh,flag,calls,toc(timer));
        if fresh<opts.krylov_tolerance && flag==0,status="KRYLOV_CONVERGED";break;end
    end
end
active=false;fresh=norm(B(x)-rhs)/max(bnorm,realmin);W=reshape(x,n^2,M);
discrete=deepsea3d_mapped_discrete_residual(W,g,p);
if bnorm==0
    r3=struct('physical_PDE',0,'true_boundary',0,'true_DtN',0,'wall_s',0);r4=r3;
else
    r3=observe_deepsea3d_fg_residual(W,p,cfg,3);r4=observe_deepsea3d_fg_residual(W,p,cfg,4);
end
offset=cfg.rx_m(1:2)-cfg.center_m;
kernel=exp(1i*p.gamma*(cfg.rx_m(3)-p.b)).*exp(1i*p.kv*(offset(1)+p.L)).*exp(1i*p.kv.'*(offset(2)+p.L));
delta=sum(fft2(reshape(W(:,end),n,n)).*kernel,'all')/n^2;bg=deepsea3d_flat(cfg,0);
[after,~]=memory;
linear=fresh<1e-8 && max([discrete.Dirichlet,discrete.DtN,discrete.mapped_PDE])<1e-8;
physical=max([r3.physical_PDE,r3.true_boundary,r3.true_DtN,r4.physical_PDE,r4.true_boundary,r4.true_DtN])<1e-4;
output=struct('status',status,'route','right Neumann P=2I-B', ...
    'fresh_preconditioned_residual',fresh,'discrete',discrete,'independent3',r3,'independent4',r4, ...
    'linear_pass',linear,'initial_physical_pass',physical,'cycles',cycles,'residual_history',hist, ...
    'rough_correction',delta,'direct_native',bg.direct_native, ...
    'reflect_native',bg.reflect_native+delta,'total_native',bg.total_native+delta, ...
    'operator_calls',calls,'preconditioner_calls',pre_calls,'preconditioner_wall_s',pre_wall, ...
    'wall_s',toc(timer),'cpu_s',cputime-cpu0,'estimated_workspace_MiB',estimate, ...
    'process_before_MiB',before.MemUsedMATLAB/2^20,'process_after_MiB',after.MemUsedMATLAB/2^20, ...
    'peak_memory_MiB',NaN,'peak_kind','unknown; snapshots only','options',opts, ...
    'G1_target_PM','OPEN','WP2_formal',false,'reference_ready',false);
assert(calls<=opts.max_operator_calls && output.wall_s<=600,'Neumann:FinalBudget','Case cap exceeded.');
    function z=B(v)
        if calls>=opts.max_operator_calls,error('Neumann:CallBudget','B call hard cap.');end
        if active && (calls>=opts.max_operator_calls-1 || toc(timer)>opts.krylov_wall_limit_s)
            error('Neumann:Budget','Reserve final fresh residual and audit.');
        end
        calls=calls+1;z=reshape(deepsea3d_mapped_action(reshape(v,n^2,M),p),[],1);
    end
    function z=pre(v)
        clock=tic;z=deepsea3d_right_neumann_action(v,@B);
        pre_wall=pre_wall+toc(clock);pre_calls=pre_calls+1;
    end
    function z=right(v)
        z=B(pre(v));
    end
end
function id=budget_id(err)
id="";
if startsWith(string(err.identifier),"Neumann:") && contains(string(err.identifier),"Budget")
    id=string(err.identifier);return
end
for j=1:numel(err.cause)
    id=budget_id(err.cause{j});if id~="",return;end
end
end
