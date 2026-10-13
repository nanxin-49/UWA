function state=validate_deepsea3d_g1_finite(stage)
%VALIDATE_DEEPSEA3D_G1_FINITE Bounded one-case finite-roughness G1 campaign.
% preflight first; then init / next (one reference solve) / compare /
% cancellation / audit. Original 3D FK and 2D/production code are untouched.
arguments
    stage (1,1) string = "next"
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
out=fullfile(root,'results','validation','deepsea3d_g1_finite');
file=fullfile(out,'finite.mat');
if stage=="init"
    a=load(fullfile(out,'preflight.mat'),'results'); pre=a.results;
    z=[pre.runs.rough_correction]; change=abs(z(1)-z(2))/abs(z(2));
    assert(change<1e-3 && max([pre.runs.wall_s])<10,'deepsea3dG1:Precheck','Precheck unsuitable.');
    state=struct('cfg',pre.cfg,'precheck_change',change, ...
        'decision','Proceed bounded FE; order6 insufficient, require independent refinement', ...
        'reference',struct([]),'comparison',struct([]),'cancellation',struct([]), ...
        'numeric_ready',false,'case_G1',false,'WP2_allowed',false);
    state.settings=[24 .125 8 32;24 .125 10 32;24 .125 12 32;24 .125 14 32; ...
        24 .0625 14 32;32 .125 14 32;48 .125 14 32;48 .0625 14 32; ...
        48 .0625 16 32;48 .0625 16 40];
else
    a=load(file,'state'); state=a.state;
end
switch stage
    case "init"
    case "next"
        ii=numel(state.reference)+1;
        assert(ii<=size(state.settings,1),'deepsea3dG1:Finished','No remaining reference jobs.');
        s=state.settings(ii,:); opts=local_opts(s);
        a=deepsea3d_dirichlet_fe(state.cfg,opts);
        row=struct('settings',s,'result',a);
        state.reference=local_append(state.reference,row);
        fprintf('G1 reference job %d/%d L=%g h=%g N=%d M=%d res=%.3g wall=%.3f\n', ...
            ii,size(state.settings,1),s(1),a.spacing_m,a.order,a.residual_order,a.boundary_residual_relative,a.wall_s);
        if ii==size(state.settings,1)
            r=state.reference; best=r(end).result;
            state.uncertainty=struct( ...
                'order_change',local_delta(r(8).result,r(9).result), ...
                'mesh_change',local_delta(r(7).result,r(8).result), ...
                'domain_change',local_delta(r(6).result,r(7).result), ...
                'fine_domain_change',local_delta(r(5).result,r(8).result), ...
                'evaluation_change',local_delta(r(9).result,r(10).result), ...
                'residual',best.boundary_residual_relative);
            u=state.uncertainty;
            state.numeric_ready=all(structfun(@(v)v<1e-4,u));
        end
    case "offgrid"
        assert(state.numeric_ready,'deepsea3dG1:Reference','Reference guards first.');
        opts=local_opts(state.settings(end,:)); opts.offgrid=true;
        state.offgrid=deepsea3d_dirichlet_fe(state.cfg,opts);
        fprintf('G1 half-cell true-boundary residual %.6g wall %.3f\n', ...
            state.offgrid.offgrid_boundary_residual_relative,state.offgrid.wall_s);
    case "compare"
        assert(state.numeric_ready,'deepsea3dG1:Reference','Reference guards first.');
        ref=state.reference(end).result; c=state.cfg;
        c.spacing_m=.05; coarse=solve_deepsea3d_fullfk(c);
        c.spacing_m=.025; fine=solve_deepsea3d_fullfk(c);
        c.spacing_m=.05; c.quadrature='gauss4'; gauss=solve_deepsea3d_fullfk(c);
        state.fk_runs=struct('coarse',coarse,'fine',fine,'gauss',gauss);
        state.fk_grid=abs(coarse.rough_correction-fine.rough_correction)/abs(fine.rough_correction);
        state.fk_gauss=abs(gauss.rough_correction-fine.rough_correction)/abs(fine.rough_correction);
        state.comparison=local_metrics(fine,ref);
        state.nonlinear_relative=abs(ref.rough_correction-ref.terms(1))/abs(ref.rough_correction);
        state.case_G1=state.comparison.reflected_complex<.01 && state.comparison.increment_complex<.01 ...
            && state.fk_grid<1e-3 && state.fk_gauss<1e-3;
        disp(state.comparison);
    case "cancellation"
        assert(~isempty(state.comparison),'deepsea3dG1:Compare','Main case comparison first.');
        % Same surface, one nearby frequency diagnosis; no band claims.
        c=state.cfg; c.spacing_m=.05;
        timer=tic; cpu0=cputime;
        [fc,~,~,search]=fminbnd(@(f)local_null(c,f),6000,6075,optimset('Display','off','TolX',.001));
        state.search=struct('frequency_hz',fc,'wall_s',toc(timer),'cpu_s',cputime-cpu0,'calls',search.funcCount);
        c.frequency_hz=fc; c.spacing_m=.025; fk=solve_deepsea3d_fullfk(c);
        low=deepsea3d_dirichlet_fe(c,local_opts([48 .125 14 32]));
        high=deepsea3d_dirichlet_fe(c,local_opts([48 .0625 16 40]));
        state.cancellation=struct('frequency_hz',fc,'fk',fk,'ref_low',low,'ref_high',high, ...
            'reference_increment_change',local_delta(low,high), ...
            'reference_total_change',abs(low.total_native-high.total_native)/abs(high.total_native), ...
            'metrics',local_metrics(fk,high), ...
            'total_condition',(abs(high.reflect_native)+abs(high.total_native-high.reflect_native))/abs(high.total_native));
        disp(state.cancellation.metrics);
    case "audit"
        assert(state.numeric_ready && ~isempty(state.comparison) && ~isempty(state.cancellation));
        assert(numel(state.reference)==size(state.settings,1));
        assert(state.fk_grid<1e-3 && state.fk_gauss<1e-3);
        assert(state.cancellation.reference_increment_change<1e-4 ...
            && state.cancellation.reference_total_change<1e-3 ...
            && state.cancellation.ref_high.boundary_residual_relative<1e-4);
        tests=readtable(fullfile(out,'tests.csv')); assert(all(tests.Passed) && ~any(tests.Failed));
        % Recompute classification from raw native complex fields.
        actual=local_metrics(state.fk_runs.fine,state.reference(end).result);
        assert(abs(actual.reflected_complex-state.comparison.reflected_complex)<1e-14);
        assert(abs(actual.increment_complex-state.comparison.increment_complex)<1e-14);
        assert(all(structfun(@(v)v<1e-4,state.uncertainty)));
        assert(state.reference(end).result.forcing_series_relative<1e-8);
        assert(isfield(state,'offgrid') && state.offgrid.offgrid_boundary_residual_relative<1e-4);
        % Conservative bound for the exponential EVALUATION Taylor tail of
        % the finite Fourier field. Not a bound on FE-order truncation or PDE.
        T=state.offgrid.max_eta_kz*sum(abs(state.cfg.profile.amplitudes_m))/state.offgrid.max_height_m;
        M=state.offgrid.residual_order;
        state.evaluation_tail_abs_bound=sum(state.offgrid.boundary_coefficient_norms)*exp(T)*T^(M+1)/factorial(M+1);
        state.case_G1=actual.reflected_complex<.01 && actual.increment_complex<.01;
        state.status='CONDITIONAL_SINGLE_CASE_PASS';
        if ~state.case_G1, state.status='PHYSICAL_SINGLE_CASE_FAIL'; end
        state.WP2_allowed=false;
        checks=struct('reference_numeric_ready',state.numeric_ready,'finite_case_G1',state.case_G1, ...
            'status',state.status,'WP2_allowed',false, ...
            'scope','one compact finite surface, q0, 6kHz; local cancellation diagnostic separate', ...
            'uncertainty_is_rigorous_bound',false);
        checks.offgrid_boundary_residual=state.offgrid.offgrid_boundary_residual_relative;
        checks.finite_fourier_evaluation_tail_abs_bound=state.evaluation_tail_abs_bound;
        fid=fopen(fullfile(out,'artifact_assertions.json'),'w'); fprintf(fid,'%s',jsonencode(checks,PrettyPrint=true)); fclose(fid);
        writetable(struct2table(state.uncertainty),fullfile(out,'reference_guards.csv'));
        offgrid=struct('residual_relative',state.offgrid.offgrid_boundary_residual_relative, ...
            'residual_abs_max',state.offgrid.offgrid_boundary_abs_max, ...
            'finite_fourier_evaluation_tail_abs_bound',state.evaluation_tail_abs_bound, ...
            'wall_s',state.offgrid.wall_s,'cpu_s',state.offgrid.cpu_s);
        writetable(struct2table(offgrid),fullfile(out,'offgrid.csv'));
        local_cost(state,out);
    otherwise
        error('deepsea3dG1:Stage','Unknown stage.');
end
save(file,'state','-v7'); local_csv(state,out);
end
function opts=local_opts(s)
opts=struct('domain_halfwidth_m',s(1),'spacing_m',s(2),'order',s(3),'residual_order',s(4));
end
function d=local_delta(a,b)
d=abs(a.rough_correction-b.rough_correction)/abs(b.rough_correction);
end
function m=local_metrics(fk,ref)
m=struct('reflected_complex',abs(fk.reflect_native-ref.reflect_native)/abs(ref.reflect_native), ...
    'reflected_magnitude',abs(abs(fk.reflect_native)-abs(ref.reflect_native))/abs(ref.reflect_native), ...
    'reflected_phase_rad',angle(fk.reflect_native*conj(ref.reflect_native)), ...
    'increment_complex',abs(fk.rough_correction-ref.rough_correction)/abs(ref.rough_correction), ...
    'increment_magnitude',abs(abs(fk.rough_correction)-abs(ref.rough_correction))/abs(ref.rough_correction), ...
    'increment_phase_rad',angle(fk.rough_correction*conj(ref.rough_correction)), ...
    'total_complex',abs(fk.total_native-ref.total_native)/abs(ref.total_native), ...
    'absolute_difference',abs(fk.reflect_native-ref.reflect_native));
end
function v=local_null(c,f)
c.frequency_hz=f; a=solve_deepsea3d_fullfk(c); v=abs(a.total_native)/abs(a.direct_native);
end
function a=local_append(a,b)
if isempty(a), a=b; else, a(end+1)=b; end
end
function local_csv(state,out)
if isempty(state.reference), return; end
rows=struct([]); terms=struct([]);
for j=1:numel(state.reference)
    r=state.reference(j); a=r.result;
    row=struct('job',j,'L_m',r.settings(1),'h_m',a.spacing_m,'N',a.order,'M',a.residual_order, ...
        'nodes',a.nodes,'delta_real',real(a.rough_correction),'delta_imag',imag(a.rough_correction), ...
        'reflect_real',real(a.reflect_native),'reflect_imag',imag(a.reflect_native), ...
        'boundary_residual',a.boundary_residual_relative,'forcing_error',a.forcing_series_relative, ...
        'wall_s',a.wall_s,'cpu_s',a.cpu_s,'spectral_MiB',a.stored_spectral_MiB);
    rows=local_append(rows,row);
    for n=1:a.order
        row=struct('job',j,'n',n,'term_real',real(a.terms(n)),'term_imag',imag(a.terms(n)), ...
            'partial_real',real(a.partials(n)),'partial_imag',imag(a.partials(n)), ...
            'boundary_coefficient_norm',a.boundary_coefficient_norms(n));
        terms=local_append(terms,row);
    end
end
writetable(struct2table(rows),fullfile(out,'reference_runs.csv'));
writetable(struct2table(terms),fullfile(out,'reference_terms.csv'));
if ~isempty(state.comparison), writetable(struct2table(state.comparison),fullfile(out,'comparison.csv')); end
if ~isempty(state.cancellation)
    c=state.cancellation; a=c.metrics;
    row=struct('frequency_hz',c.frequency_hz,'reference_increment_change',c.reference_increment_change, ...
        'reference_total_change',c.reference_total_change,'total_condition',c.total_condition, ...
        'reflected_complex',a.reflected_complex,'increment_complex',a.increment_complex, ...
        'reflected_phase_rad',a.reflected_phase_rad,'increment_phase_rad',a.increment_phase_rad, ...
        'total_complex',a.total_complex,'reference_residual',c.ref_high.boundary_residual_relative);
    writetable(struct2table(row),fullfile(out,'cancellation.csv'));
end
end
function local_cost(state,out)
cost=struct('final_reference_wall_s',sum(arrayfun(@(x)x.result.wall_s,state.reference)), ...
    'final_reference_cpu_s',sum(arrayfun(@(x)x.result.cpu_s,state.reference)), ...
    'fk_wall_s',state.fk_runs.coarse.wall_s+state.fk_runs.fine.wall_s+state.fk_runs.gauss.wall_s, ...
    'fk_cpu_s',state.fk_runs.coarse.cpu_s+state.fk_runs.fine.cpu_s+state.fk_runs.gauss.cpu_s, ...
    'cancellation_reference_wall_s',state.cancellation.ref_low.wall_s+state.cancellation.ref_high.wall_s, ...
    'cancellation_reference_cpu_s',state.cancellation.ref_low.cpu_s+state.cancellation.ref_high.cpu_s, ...
    'cancellation_fk_wall_s',state.cancellation.fk.wall_s,'cancellation_fk_cpu_s',state.cancellation.fk.cpu_s, ...
    'search_wall_s',state.search.wall_s,'search_cpu_s',state.search.cpu_s);
cost.offgrid_wall_s=state.offgrid.wall_s; cost.offgrid_cpu_s=state.offgrid.cpu_s;
pre=load(fullfile(out,'preflight.mat'),'results');
cost.preflight_wall_s=sum([pre.results.runs.wall_s]); cost.preflight_cpu_s=sum([pre.results.runs.cpu_s]);
% Optional exploratory probes retained in this case contribute to actual cost,
% never to the authoritative reference classification.
if isfile(fullfile(out,'order_probe.mat'))
    a=load(fullfile(out,'order_probe.mat'),'g1probe');
    cost.order_probe_wall_s=sum([a.g1probe.wall_s]); cost.order_probe_cpu_s=sum([a.g1probe.cpu_s]);
end
if isfile(fullfile(out,'reference_guards_probe.mat'))
    a=load(fullfile(out,'reference_guards_probe.mat'),'g1guard');
    cost.guard_probe_wall_s=sum(arrayfun(@(x)x.result.wall_s,a.g1guard));
    cost.guard_probe_cpu_s=sum(arrayfun(@(x)x.result.cpu_s,a.g1guard));
end
names=string(fieldnames(cost));
cost.total_wall_s=sum(cellfun(@(n)cost.(n),cellstr(names(endsWith(names,'wall_s')))));
cost.total_cpu_s=sum(cellfun(@(n)cost.(n),cellstr(names(endsWith(names,'cpu_s')))));
cost.exclusions='IO, tests, static analysis, documents, unsaved development; memory arrays not process peak';
fid=fopen(fullfile(out,'cost_summary.json'),'w');fprintf(fid,'%s',jsonencode(cost,PrettyPrint=true));fclose(fid);
end
