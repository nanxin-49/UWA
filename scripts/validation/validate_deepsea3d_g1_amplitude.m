function state=validate_deepsea3d_g1_amplitude(stage)
%VALIDATE_DEEPSEA3D_G1_AMPLITUDE One-variable extension of frozen 3D G1.
% init; repeated next; audit. At most one FE solve per next call.
arguments
    stage (1,1) string = "next"
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
out=fullfile(root,'results','validation','deepsea3d_g1_amplitude');
file=fullfile(out,'amplitude.mat');
if stage=="init"
    assert(~isfile(file),'deepsea3dAmp:Exists','Use next/audit on existing state.');
    if ~isfolder(out), mkdir(out); end
    b=load(fullfile(root,'results','validation','deepsea3d_g1_finite','finite.mat'),'state');
    assert(b.state.case_G1 && b.state.numeric_ready && ~b.state.WP2_allowed);
    state=struct('baseline',b.state,'cfg',b.state.cfg,'cases',struct([]), ...
        'optional_decided',false,'optional_allowed',false,'complete',false,'WP2_allowed',false);
    state.cases=[local_case(.5),local_case(1.5)];
else
    a=load(file,'state'); state=a.state;
end
switch stage
    case "init"
    case "next"
        ix=find(~[state.cases.done],1);
        if isempty(ix) && ~state.optional_decided
            state.optional_decided=true;
            state.optional_allowed=all([state.cases.reference_ready]) && all([state.cases.max_wall_s]<=60);
            if state.optional_allowed, state.cases(end+1)=local_case(2); ix=numel(state.cases); end
        end
        if isempty(ix)
            state.complete=true;
        else
            c=state.cases(ix); cfg=state.cfg; cfg.profile.amplitudes_m=state.cfg.profile.amplitudes_m*c.alpha;
            if isempty(fieldnames(c.geometry)), c.geometry=local_geometry(cfg); end
            if c.phase=="reference"
                job=numel(c.runs)+1; s=c.settings(job,:);
                opts=struct('domain_halfwidth_m',s(1),'spacing_m',s(2),'order',s(3), ...
                    'residual_order',s(4),'offgrid',logical(s(5)));
                r=deepsea3d_dirichlet_fe(cfg,opts);
                c.runs=local_append(c.runs,struct('settings',s,'result',r));
                c.max_wall_s=max(c.max_wall_s,r.wall_s);
                fprintf('AMP %.2f job %d L%g h%g N%d M%d off%d: res %.3g wall %.3f\n', ...
                    c.alpha,job,s(1),s(2),s(3),s(4),s(5),r.boundary_residual_relative,r.wall_s);
                c=local_schedule(c);
            elseif c.phase=="fk"
                cfg.spacing_m=.05; coarse=solve_deepsea3d_fullfk(cfg);
                cfg.spacing_m=.025; fine=solve_deepsea3d_fullfk(cfg);
                cfg.spacing_m=.05; cfg.quadrature='gauss4'; gauss=solve_deepsea3d_fullfk(cfg);
                c.fk=struct('coarse',coarse,'fine',fine,'gauss',gauss);
                c.fk_grid=local_delta(coarse,fine); c.fk_gauss=local_delta(gauss,fine);
                ref=c.runs(c.best_index).result; c.metrics=local_metrics(fine,ref);
                [c.status,c.reason]=classify_deepsea3d_amplitude(c.reference_ready,c.fk_grid,c.fk_gauss,c.metrics);
                c.done=true; c.phase="done";
                fprintf('AMP %.2f %s Eref %.6g Edelta %.6g\n',c.alpha,c.status,c.metrics.reflected_complex,c.metrics.increment_complex);
            end
            state.cases(ix)=c;
        end
    case "audit"
        assert(state.complete && state.optional_decided,'deepsea3dAmp:Incomplete','Finish next jobs.');
        eval_checks=struct([]);
        for c=state.cases
            assert(c.done && c.alpha~=1);
            if c.reference_ready
                verified=local_certify(c);
                assert(verified.reference_ready && isequal(verified.guards,c.guards));
                assert(all(structfun(@(v)v<1e-4,c.guards)));
                m=local_metrics(c.fk.fine,c.runs(c.best_index).result);
                assert(abs(m.reflected_complex-c.metrics.reflected_complex)<1e-14);
                assert(abs(m.increment_complex-c.metrics.increment_complex)<1e-14);
                assert(classify_deepsea3d_amplitude(true,c.fk_grid,c.fk_gauss,m)==c.status);
                a=c.runs(c.best_index).result;
                assert(a.forcing_series_relative<1e-8);
                T=a.max_eta_kz*c.alpha*sum(abs(state.cfg.profile.amplitudes_m))/a.max_height_m;
                bound=sum(a.boundary_coefficient_norms)*exp(T)*T^(a.residual_order+1)/factorial(a.residual_order+1);
                eval_checks=local_append(eval_checks,struct('alpha',c.alpha, ...
                    'forcing_series_relative',a.forcing_series_relative,'evaluation_tail_abs_bound',bound));
            else
                assert(c.status=="INCONCLUSIVE" && isempty(fieldnames(c.metrics)));
            end
        end
        % Exact equality of every non-amplitude input to frozen baseline.
        assert(isequal(state.cfg,state.baseline.cfg));
        tests=readtable(fullfile(out,'classification_tests.csv')); assert(all(tests.Passed) && ~any(tests.Failed));
        checks=struct('complete',true,'new_case_count',numel(state.cases),'baseline_reused',true, ...
            'optional_two_allowed',state.optional_allowed,'WP2_allowed',false,'target_domain','OPEN', ...
            'no_phase_fit',true,'reference_failures_are_inconclusive',true);
        fid=fopen(fullfile(out,'artifact_assertions.json'),'w');fprintf(fid,'%s',jsonencode(checks,PrettyPrint=true));fclose(fid);
        if ~isempty(eval_checks),writetable(struct2table(eval_checks),fullfile(out,'evaluation_checks.csv'));end
    otherwise
        error('deepsea3dAmp:Stage','Unknown stage.');
end
save(file,'state','-v7'); local_outputs(state,out);
fprintf('AMP checkpoint cases=%d done=%d complete=%d\n',numel(state.cases),sum([state.cases.done]),state.complete);
end
function c=local_case(alpha)
c=struct('alpha',alpha,'settings',[24 .125 12 40 0;24 .125 16 40 0], ...
    'runs',struct([]),'phase',"reference",'selected_order',0,'best_index',0, ...
    'reference_ready',false,'done',false,'status',"PENDING",'reason',"", ...
    'max_wall_s',0,'guards',struct(),'metrics',struct(),'geometry',struct(), ...
    'fk',struct(),'fk_grid',NaN,'fk_gauss',NaN,'extra_order',false);
end
function c=local_schedule(c)
n=numel(c.runs);
if c.selected_order==0 && n>=2
    a=c.runs(n-1).result; b=c.runs(n).result;
    if local_delta(a,b)<1e-4 && b.boundary_residual_relative<1e-4
        c.selected_order=b.order;
        N=b.order;
        c.settings=[c.settings;24 .125 N 32 0;24 .0625 N 40 1;48 .125 N 40 0;48 .0625 N 40 1];
    elseif b.order<20
        c.settings=[c.settings;24 .125 20 40 0];
    else
        c.done=true; c.phase="done"; c.status="INCONCLUSIVE";
        c.reason="Coarse order/residual failed at frozen maximum N20";
    end
elseif c.selected_order>0 && n==size(c.settings,1)
    c=local_certify(c);
end
end
function c=local_certify(c)
N=c.selected_order; r=c.runs;
findjob=@(L,h,M,off)find(arrayfun(@(a)isequal(a.settings,[L h N M off]),r),1,'last');
coarse=findjob(24,.125,40,0); eval32=findjob(24,.125,32,0);
fine24=findjob(24,.0625,40,1); coarse48=findjob(48,.125,40,0); best=findjob(48,.0625,40,1);
b=r(best).result;
lower_index=find(arrayfun(@(a)a.settings(1)==24 && a.settings(2)==.125 && ...
    a.settings(3)<N && a.settings(4)==40 && a.settings(5)==0,r),1,'last');
early=r(coarse).result; lower=r(lower_index).result;
g=struct('coarse_order',local_delta(lower,early), ...
    'fine_order',abs(b.partials(end)-b.partials(end-2))/abs(b.rough_correction), ...
    'grid',local_delta(r(coarse48).result,b),'domain',local_delta(r(fine24).result,b), ...
    'evaluation',abs(r(eval32).result.boundary_residual_relative-early.boundary_residual_relative), ...
    'boundary',b.boundary_residual_relative,'offgrid',b.offgrid_boundary_residual_relative);
c.guards=g; c.best_index=best;
c.reference_ready=all(structfun(@(v)v<1e-4,g));
if c.reference_ready
    c.phase="fk";
elseif N<20 && ~c.extra_order && c.max_wall_s<=60
    c.extra_order=true; c.selected_order=20;
    c.settings=[c.settings;24 .125 20 40 0;24 .125 20 32 0; ...
        24 .0625 20 40 1;48 .125 20 40 0;48 .0625 20 40 1];
    c.phase="reference";
else
    c.status="INCONCLUSIVE"; c.done=true; c.phase="done";
    c.reason="Independent order/grid/domain/true-boundary reference guard failed";
end
end
function g=local_geometry(cfg)
x=-12:.025:12; [xx,yy]=meshgrid(x,x); [eta,ex,ey]=deepsea3d_surface(xx,yy,cfg.profile);
mask=hypot(xx,yy)<=12; slope=hypot(ex,ey);
g=struct('height_rms',sqrt(mean(eta(mask).^2)),'height_max',max(abs(eta(mask))), ...
    'slope_rms',sqrt(mean(slope(mask).^2)),'slope_max',max(slope(mask)));
end
function d=local_delta(a,b)
d=abs(a.rough_correction-b.rough_correction)/abs(b.rough_correction);
end
function m=local_metrics(fk,ref)
m=struct('reflected_complex',abs(fk.reflect_native-ref.reflect_native)/abs(ref.reflect_native), ...
    'increment_complex',abs(fk.rough_correction-ref.rough_correction)/abs(ref.rough_correction), ...
    'reflected_phase',angle(fk.reflect_native*conj(ref.reflect_native)), ...
    'increment_phase',angle(fk.rough_correction*conj(ref.rough_correction)), ...
    'reflected_magnitude',abs(abs(fk.reflect_native)-abs(ref.reflect_native))/abs(ref.reflect_native), ...
    'increment_magnitude',abs(abs(fk.rough_correction)-abs(ref.rough_correction))/abs(ref.rough_correction));
end
function a=local_append(a,b)
if isempty(a), a=b; else, a(end+1)=b; end
end
function local_outputs(state,out)
rows=struct([]); terms=struct([]);
for c=state.cases
    for j=1:numel(c.runs)
        a=c.runs(j).result; s=c.runs(j).settings;
        off=NaN; if isfield(a,'offgrid_boundary_residual_relative'),off=a.offgrid_boundary_residual_relative;end
        row=struct('alpha',c.alpha,'job',j,'L',s(1),'h',s(2),'N',s(3),'M',s(4),'offgrid',off, ...
            'boundary',a.boundary_residual_relative,'delta_real',real(a.rough_correction), ...
            'delta_imag',imag(a.rough_correction),'reflect_real',real(a.reflect_native), ...
            'reflect_imag',imag(a.reflect_native),'wall_s',a.wall_s,'cpu_s',a.cpu_s,'spectral_MiB',a.stored_spectral_MiB);
        rows=local_append(rows,row);
        for n=1:a.order
            row=struct('alpha',c.alpha,'job',j,'order',n,'real',real(a.terms(n)),'imag',imag(a.terms(n)), ...
                'boundary_coefficient_norm',a.boundary_coefficient_norms(n)); terms=local_append(terms,row);
        end
    end
end
if ~isempty(rows),writetable(struct2table(rows),fullfile(out,'reference_runs.csv'));end
if ~isempty(terms),writetable(struct2table(terms),fullfile(out,'reference_terms.csv'));end
summary=struct([]); guards=struct([]);
b=state.baseline; geom=local_geometry(state.cfg);
m=struct('reflected_complex',b.comparison.reflected_complex,'increment_complex',b.comparison.increment_complex, ...
    'reflected_phase',b.comparison.reflected_phase_rad,'increment_phase',b.comparison.increment_phase_rad, ...
    'reflected_magnitude',b.comparison.reflected_magnitude,'increment_magnitude',b.comparison.increment_magnitude);
row=local_summary(1,"CONDITIONAL_G1_PASS_REUSED",true,16,geom,m,0,0);
summary=local_append(summary,row);
for c=state.cases
    if ~c.done, continue; end
    m=c.metrics;
    if isempty(fieldnames(m))
        m=struct('reflected_complex',NaN,'increment_complex',NaN,'reflected_phase',NaN, ...
            'increment_phase',NaN,'reflected_magnitude',NaN,'increment_magnitude',NaN);
    end
    wall=sum(arrayfun(@(r)r.result.wall_s,c.runs)); cpu=sum(arrayfun(@(r)r.result.cpu_s,c.runs));
    if ~isempty(fieldnames(c.fk))
        wall=wall+c.fk.coarse.wall_s+c.fk.fine.wall_s+c.fk.gauss.wall_s;
        cpu=cpu+c.fk.coarse.cpu_s+c.fk.fine.cpu_s+c.fk.gauss.cpu_s;
    end
    row=local_summary(c.alpha,c.status,false,c.selected_order,c.geometry,m,wall,cpu);
    summary=local_append(summary,row);
    if ~isempty(fieldnames(c.guards))
        g=c.guards; row=struct('alpha',c.alpha,'coarse_order',g.coarse_order,'fine_order',g.fine_order, ...
            'grid',g.grid,'domain',g.domain,'evaluation',g.evaluation,'boundary',g.boundary,'offgrid',g.offgrid, ...
            'reference_ready',c.reference_ready,'fk_grid',c.fk_grid,'fk_gauss',c.fk_gauss);
        guards=local_append(guards,row);
    end
end
table_summary=sortrows(struct2table(summary),'alpha');
writetable(table_summary,fullfile(out,'summary.csv'));
if ~isempty(guards),writetable(struct2table(guards),fullfile(out,'guards.csv'));end
cost=struct('new_wall_s',sum(table_summary.wall_s),'new_cpu_s',sum(table_summary.cpu_s), ...
    'baseline_new_solves',0,'scope','saved FE/FK kernels; excludes IO, geometry metrics, tests and reporting');
fid=fopen(fullfile(out,'cost_summary.json'),'w');fprintf(fid,'%s',jsonencode(cost,PrettyPrint=true));fclose(fid);
end
function row=local_summary(alpha,status,reused,N,g,m,wall,cpu)
row=struct('alpha',alpha,'status',char(status),'reused',reused,'selected_N',N, ...
    'height_rms_m',g.height_rms,'height_max_m',g.height_max,'slope_rms',g.slope_rms,'slope_max',g.slope_max, ...
    'reflected_complex',m.reflected_complex,'increment_complex',m.increment_complex, ...
    'reflected_phase_rad',m.reflected_phase,'increment_phase_rad',m.increment_phase, ...
    'reflected_magnitude',m.reflected_magnitude,'increment_magnitude',m.increment_magnitude, ...
    'wall_s',wall,'cpu_s',cpu);
end
