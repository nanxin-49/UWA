function state=validate_deepsea3d_alpha25_domain(stage)
%VALIDATE_DEEPSEA3D_ALPHA25_DOMAIN Frozen alpha2.5 domain then gated FullFK.
% Each next/observe/fk call checkpoints one computation; status is read-only.
arguments
    stage (1,1) string = "status"
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
out=fullfile(root,'results','validation','deepsea3d_alpha25_domain');
file=fullfile(out,'domain.mat');
if stage=="init"
    assert(~isfile(file),'alpha25domain:Exists','Use existing checkpoint or a separate checkout.');
    if ~isfolder(out),mkdir(out);end
    a=load(fullfile(root,'results','validation','deepsea3d_alpha25_fine','fine.mat'),'state');
    assert(a.state.reference_status=="FINE_REFERENCE_CHECK_PASS");
    state=struct('cfg',a.state.cfg,'baseline',a.state,'domains',32,'runs',struct([]), ...
        'snapshots',struct([]),'observers',struct([]),'observer_snapshots',struct([]), ...
        'reference_ready',false,'reference_decided',false,'reference_guards',struct(), ...
        'fk_settings',[16 .05 0;16 .025 0;16 .05 1;12 .025 0;24 .025 0], ...
        'fk',struct([]),'fk_snapshots',struct([]),'fk_guards',struct(), ...
        'metrics',struct(),'G1_status',"INCONCLUSIVE",'reason',"Reference pending", ...
        'complete',false,'WP2_allowed',false,'global_G1',"OPEN", ...
        'matlab_version',version,'created_local',char(datetime('now')));
    assert(isequal(state.cfg.tx_m,[0 0 1500]) && isequal(state.cfg.rx_m,[0 0 5]));
    assert(state.cfg.frequency_hz==6000 && state.cfg.c0_mps==1500);
    assert(isequal(state.cfg.profile.amplitudes_m,[.05 .025]));
    assert(state.cfg.profile.plateau_m==8 && state.cfg.profile.support_m==12);
else
    a=load(file,'state');state=a.state;
end
switch stage
    case "init"
    case "next"
        assert(~state.reference_decided,'alpha25domain:Decided','Reference decision already made.');
        j=numel(state.runs)+1;
        assert(j<=numel(state.domains),'alpha25domain:Pending','Run observe/reference first.');
        opts=local_opts(state.domains(j));
        [mb,~]=memory;r=deepsea3d_dirichlet_fe(state.cfg,opts);[ma,~]=memory;
        assert(local_finite(r),'alpha25domain:Nonfinite','Nonfinite reference.');
        state.runs=local_append(state.runs,r);
        state.snapshots=local_append(state.snapshots,local_memory(mb,ma));
        fprintf('DOMAIN L%g N20 delta=(%.12g,%.12g), boundary %.3g half %.3g wall %.3fs\n', ...
            state.domains(j),real(r.rough_correction),imag(r.rough_correction), ...
            r.boundary_residual_relative,r.offgrid_boundary_residual_relative,r.wall_s);
    case "observe"
        j=numel(state.observers)+1;
        assert(j<=numel(state.runs),'alpha25domain:Pending','Run next first.');
        opts=local_opts(state.domains(j));
        [mb,~]=memory;r=observe_deepsea3d_fe_cancellation(state.cfg,opts,state.runs(j));[ma,~]=memory;
        state.observers=local_append(state.observers,r);
        state.observer_snapshots=local_append(state.observer_snapshots,local_memory(mb,ma));
        fprintf('OBSERVER L%g max %.6g match %.3g/%.3g wall %.3fs\n', ...
            state.domains(j),max(r.conditions),r.term_match_error,r.norm_match_error,r.wall_s);
    case "reference"
        assert(~state.reference_decided,'alpha25domain:Decided','Reference decision already made.');
        assert(numel(state.runs)==numel(state.domains) && numel(state.observers)==numel(state.runs));
        g=local_reference_guards(state);state.reference_guards=g;
        other=rmfield(g,'domain');
        good=all(structfun(@(v)isfinite(v) && v<1e-4,other)) && ...
            all(arrayfun(@local_tails,state.runs)) && all([state.runs.forcing_series_relative]<1e-8) && ...
            all([state.observers.term_match_error]<1e-12) && all([state.observers.norm_match_error]<1e-12);
        if good && g.domain<1e-4
            state.reference_ready=true;state.reference_decided=true;state.reason="FullFK pending";
        elseif good && state.domains(end)==32
            state.domains(end+1)=48;state.reason="L24 to L32 domain guard requires L48";
        else
            state.reference_decided=true;state.complete=true;
            state.reason="Independent reference domain/order/boundary/tail unresolved";
        end
        fprintf('REFERENCE ready=%d decided=%d next domain=%g\n', ...
            state.reference_ready,state.reference_decided,state.domains(end));disp(g);
    case "fk"
        assert(state.reference_decided && state.reference_ready,'alpha25domain:ReferenceGate','FullFK forbidden before reference PASS.');
        j=numel(state.fk)+1;assert(j<=size(state.fk_settings,1),'alpha25domain:Finished','FullFK jobs finished.');
        s=state.fk_settings(j,:);cfg=state.cfg;cfg.half_width_m=s(1);cfg.spacing_m=s(2);
        if s(3),cfg.quadrature='gauss4';else,cfg.quadrature='trap';end
        [mb,~]=memory;r=solve_deepsea3d_fullfk(cfg);[ma,~]=memory;
        assert(local_finite(r),'alpha25domain:Nonfinite','Nonfinite FullFK.');
        state.fk=local_append(state.fk,r);state.fk_snapshots=local_append(state.fk_snapshots,local_memory(mb,ma));
        fprintf('FK A%g h%g %s delta=(%.12g,%.12g) wall %.3fs\n', ...
            s(1),s(2),cfg.quadrature,real(r.rough_correction),imag(r.rough_correction),r.wall_s);
        if j==size(state.fk_settings,1)
            state.fk_guards=local_fk_guards(state.fk);
            state.metrics=local_metrics(state.fk(2),state.runs(end));
            g=state.fk_guards;
            [state.G1_status,state.reason]=classify_deepsea3d_amplitude(true, ...
                max([g.grid,g.aperture_12_16,g.aperture_16_24]),g.gauss,state.metrics);
            state.complete=true;
        end
    case "audit"
        assert(state.complete,'alpha25domain:Incomplete','Finish the scheduled computations.');
        assert(isequal(state.cfg,state.baseline.cfg) && ~state.WP2_allowed && state.global_G1=="OPEN");
        assert(isequaln(local_reference_guards(state),state.reference_guards));
        for j=1:numel(state.runs)
            r=state.runs(j);o=state.observers(j);
            assert(local_finite(r) && local_finite(o));
            assert(r.forcing_series_relative<1e-8 && o.term_match_error<1e-12 && o.norm_match_error<1e-12);
            assert(norm(o.terms-r.terms)/norm(r.terms)<1e-12);
        end
        if state.reference_ready
            assert(all(structfun(@(v)v<1e-4,state.reference_guards)) && local_tails(state.runs(end)));
            assert(numel(state.fk)==5 && isequaln(local_fk_guards(state.fk),state.fk_guards));
            m=local_metrics(state.fk(2),state.runs(end));assert(isequaln(m,state.metrics));
            g=state.fk_guards;
            assert(classify_deepsea3d_amplitude(true,max([g.grid,g.aperture_12_16,g.aperture_16_24]),g.gauss,m)==state.G1_status);
        else
            assert(isempty(state.fk) && state.G1_status=="INCONCLUSIVE");
        end
        checks=struct('complete',true,'reference_ready',state.reference_ready, ...
            'alpha25_G1',state.G1_status,'global_G1',"OPEN",'WP2_allowed',false, ...
            'FE_new_solves',numel(state.runs),'FK_solves',numel(state.fk), ...
            'baseline_reused',true,'physical_cfg_unchanged',true,'no_phase_fit',true, ...
            'finite_results',true,'diagnostic_replay_matches',true, ...
            'reference_threshold',1e-4,'FK_threshold',1e-3,'physical_threshold',.01);
        local_json(fullfile(out,'artifact_assertions.json'),checks);
        fprintf('AUDIT %s; global G1 OPEN; WP2 false\n',state.G1_status);
    case "status"
        fprintf('reference=%d decided=%d FE=%d observers=%d FK=%d complete=%d G1=%s\n', ...
            state.reference_ready,state.reference_decided,numel(state.runs),numel(state.observers),numel(state.fk),state.complete,state.G1_status);
        return
    otherwise
        error('alpha25domain:Stage','Unknown stage.');
end
save(file,'state','-v7');local_outputs(state,out);
end
function opts=local_opts(L)
opts=struct('domain_halfwidth_m',L,'spacing_m',.0625,'order',20,'residual_order',40,'offgrid',true);
end
function g=local_reference_guards(s)
b=s.runs(end);
if isscalar(s.runs),a=s.baseline.runs(2);else,a=s.runs(end-1);end
g=struct('inherited_grid_N20',s.baseline.guards.grid_N20, ...
    'domain',local_delta(a,b), ...
    'order_N16_N20',abs(b.partials(16)-b.partials(20))/abs(b.rough_correction), ...
    'last_two',abs(b.partials(18)-b.partials(20))/abs(b.rough_correction), ...
    'boundary',b.boundary_residual_relative,'halfgrid',b.offgrid_boundary_residual_relative);
end
function g=local_fk_guards(f)
g=struct('grid',local_delta(f(1),f(2)),'gauss',local_delta(f(3),f(2)), ...
    'aperture_12_16',local_delta(f(4),f(2)),'aperture_16_24',local_delta(f(2),f(5)));
end
function yes=local_tails(r)
yes=all(diff(r.boundary_coefficient_norms(end-4:end))<0) && all(diff(abs(r.terms(end-4:end)))<0);
end
function yes=local_finite(s)
yes=true;names=fieldnames(s);
for j=1:numel(names)
    v=s.(names{j});if isnumeric(v),yes=yes && all(isfinite(v(:)));end
end
end
function d=local_delta(a,b)
d=abs(a.rough_correction-b.rough_correction)/abs(b.rough_correction);
end
function m=local_metrics(f,r)
m=struct('reflected_complex',abs(f.reflect_native-r.reflect_native)/abs(r.reflect_native), ...
    'increment_complex',abs(f.rough_correction-r.rough_correction)/abs(r.rough_correction), ...
    'reflected_magnitude',abs(abs(f.reflect_native)-abs(r.reflect_native))/abs(r.reflect_native), ...
    'increment_magnitude',abs(abs(f.rough_correction)-abs(r.rough_correction))/abs(r.rough_correction), ...
    'reflected_phase_rad',angle(f.reflect_native*conj(r.reflect_native)), ...
    'increment_phase_rad',angle(f.rough_correction*conj(r.rough_correction)));
end
function x=local_append(x,y)
if isempty(x),x=y;else,x(end+1)=y;end
end
function s=local_memory(a,b)
s=struct('before_MiB',a.MemUsedMATLAB/2^20,'after_MiB',b.MemUsedMATLAB/2^20);
end
function local_json(file,s)
fid=fopen(file,'w');assert(fid>=0);clean=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(s,PrettyPrint=true));
end
function local_outputs(s,out)
rows=struct([]);terms=struct([]);cancels=struct([]);obs=struct([]);evaluation=struct([]);
for j=0:numel(s.runs)
    if j==0,r=s.baseline.runs(2);L=24;mem=s.baseline.snapshots(2);else,r=s.runs(j);L=s.domains(j);mem=s.snapshots(j);end
    domain_change=NaN;
    if j==1,domain_change=local_delta(s.baseline.runs(2),r);elseif j>1,domain_change=local_delta(s.runs(j-1),r);end
    row=struct('L',L,'h',r.spacing_m,'N',r.order,'M',r.residual_order,'reused',j==0, ...
        'nodes',r.nodes,'delta_real',real(r.rough_correction),'delta_imag',imag(r.rough_correction), ...
        'reflect_real',real(r.reflect_native),'reflect_imag',imag(r.reflect_native), ...
        'N16_real',real(r.partials(16)),'N16_imag',imag(r.partials(16)), ...
        'domain_step_relative',domain_change, ...
        'N16_N20_relative',abs(r.partials(16)-r.partials(20))/abs(r.rough_correction), ...
        'N18_N20_relative',abs(r.partials(18)-r.partials(20))/abs(r.rough_correction), ...
        'boundary',r.boundary_residual_relative,'halfgrid',r.offgrid_boundary_residual_relative, ...
        'coefficient_tail_decreases',all(diff(r.boundary_coefficient_norms(end-4:end))<0), ...
        'pressure_tail_decreases',all(diff(abs(r.terms(end-4:end)))<0), ...
        'pressure_sum_indicator',sum(abs(r.terms))/abs(r.rough_correction), ...
        'wall_s',r.wall_s,'cpu_s',r.cpu_s,'spectral_MiB',r.stored_spectral_MiB, ...
        'process_before_MiB',mem.before_MiB,'process_after_MiB',mem.after_MiB);
    rows=local_append(rows,row);
    for n=1:r.order
        terms=local_append(terms,struct('L',L,'n',n,'term_real',real(r.terms(n)), ...
            'term_imag',imag(r.terms(n)),'partial_real',real(r.partials(n)), ...
            'partial_imag',imag(r.partials(n)),'boundary_coefficient_norm',r.boundary_coefficient_norms(n), ...
            'h_normalized_L2',r.spacing_m*r.boundary_coefficient_norms(n)));
    end
    T=r.max_eta_kz*sum(abs(s.cfg.profile.amplitudes_m))/r.max_height_m;
    bound=sum(r.boundary_coefficient_norms)*exp(T)*T^(r.residual_order+1)/factorial(r.residual_order+1);
    evaluation=local_append(evaluation,struct('L',L,'forcing_error',r.forcing_series_relative, ...
        'finite_Fourier_exponential_tail_abs_bound',bound));
end
writetable(struct2table(rows),fullfile(out,'reference_runs.csv'));
writetable(struct2table(terms),fullfile(out,'reference_terms.csv'));
writetable(struct2table(evaluation),fullfile(out,'evaluation_checks.csv'));
for j=1:numel(s.observers)
    r=s.observers(j);mem=s.observer_snapshots(j);L=s.domains(j);
    obs=local_append(obs,struct('L',L,'max_indicator',max(r.conditions), ...
        'term_match_error',r.term_match_error,'norm_match_error',r.norm_match_error, ...
        'wall_s',r.wall_s,'cpu_s',r.cpu_s,'process_before_MiB',mem.before_MiB,'process_after_MiB',mem.after_MiB));
    for n=1:numel(r.conditions)
        cancels=local_append(cancels,struct('L',L,'n',n,'coefficient_norm',r.bnorm(n), ...
            'contribution_norm_sum',r.numerators(n),'cancellation_indicator',r.conditions(n)));
    end
end
if ~isempty(obs),writetable(struct2table(obs),fullfile(out,'observers.csv'));end
if ~isempty(cancels),writetable(struct2table(cancels),fullfile(out,'internal_cancellation.csv'));end
fk=struct([]);
for j=1:numel(s.fk)
    r=s.fk(j);setting=s.fk_settings(j,:);mem=s.fk_snapshots(j);
    fk=local_append(fk,struct('A',setting(1),'h',setting(2),'gauss4',logical(setting(3)), ...
        'nodes',r.node_count,'reflect_real',real(r.reflect_native),'reflect_imag',imag(r.reflect_native), ...
        'delta_real',real(r.rough_correction),'delta_imag',imag(r.rough_correction), ...
        'raw_truncated_real',real(r.raw_truncated_reflect),'raw_truncated_imag',imag(r.raw_truncated_reflect), ...
        'identity_residual',r.quadrature_identity_residual,'cancellation_indicator',r.cancellation_condition, ...
        'wall_s',r.wall_s,'cpu_s',r.cpu_s,'block_payload_MiB',r.estimated_block_payload_MiB, ...
        'process_before_MiB',mem.before_MiB,'process_after_MiB',mem.after_MiB));
end
if ~isempty(fk),writetable(struct2table(fk),fullfile(out,'fk_runs.csv'));end
if ~isempty(fieldnames(s.reference_guards)),writetable(struct2table(s.reference_guards),fullfile(out,'reference_guards.csv'));end
if ~isempty(fieldnames(s.fk_guards)),writetable(struct2table(s.fk_guards),fullfile(out,'fk_guards.csv'));end
if ~isempty(fieldnames(s.metrics)),writetable(struct2table(s.metrics),fullfile(out,'physical_metrics.csv'));end
cost=struct('FE_wall_s',local_sum(s.runs,'wall_s'),'FE_cpu_s',local_sum(s.runs,'cpu_s'), ...
    'observer_wall_s',local_sum(s.observers,'wall_s'),'observer_cpu_s',local_sum(s.observers,'cpu_s'), ...
    'FK_wall_s',local_sum(s.fk,'wall_s'),'FK_cpu_s',local_sum(s.fk,'cpu_s'),'baseline_new_solves',0, ...
    'scope','new saved FE, replay and FK kernels; excludes IO, tests, reporting; memory snapshots are not peaks');
local_json(fullfile(out,'cost_summary.json'),cost);
end
function value=local_sum(rows,name)
value=0;for j=1:numel(rows),value=value+rows(j).(name);end
end
