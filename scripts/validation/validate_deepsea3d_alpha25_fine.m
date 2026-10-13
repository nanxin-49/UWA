function state=validate_deepsea3d_alpha25_fine(stage)
%VALIDATE_DEEPSEA3D_ALPHA25_FINE Frozen alpha2.5, L24 h=.0625, N16/N20 only.
arguments
    stage (1,1) string = "next"
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
out=fullfile(root,'results','validation','deepsea3d_alpha25_fine');
file=fullfile(out,'fine.mat');
if stage=="init"
    assert(~isfile(file),'alpha25fine:Exists','Use next/observe/audit on existing checkpoint.');
    if ~isfolder(out),mkdir(out);end
    a=load(fullfile(root,'results','validation','deepsea3d_alpha25_preflight','preflight.mat'),'state');
    assert(a.state.complete && a.state.precheck_status=="PROCEED_FORMAL_G1");
    state=struct('cfg',a.state.cfg,'coarse',a.state,'orders',[16 20],'runs',struct([]), ...
        'snapshots',struct([]),'observers',struct([]),'observer_snapshots',struct([]), ...
        'reference_status',"PENDING",'G1_status',"INCONCLUSIVE",'WP2_allowed',false);
else
    a=load(file,'state');state=a.state;
end
switch stage
    case "init"
    case "next"
        j=numel(state.runs)+1;assert(j<=2,'alpha25fine:Finished','Reference jobs finished.');
        opts=struct('domain_halfwidth_m',24,'spacing_m',.0625,'order',state.orders(j),'residual_order',40,'offgrid',true);
        [mb,~]=memory;r=deepsea3d_dirichlet_fe(state.cfg,opts);[ma,~]=memory;
        state.runs=local_append(state.runs,r);
        state.snapshots=local_append(state.snapshots,struct('before_MiB',mb.MemUsedMATLAB/2^20,'after_MiB',ma.MemUsedMATLAB/2^20));
        fprintf('ALPHA2.5 FINE N%d delta=(%.11g,%.11g) boundary %.4g half %.4g wall %.3f\n', ...
            r.order,real(r.rough_correction),imag(r.rough_correction),r.boundary_residual_relative,r.offgrid_boundary_residual_relative,r.wall_s);
    case "observe"
        assert(numel(state.runs)==2,'alpha25fine:Incomplete','Finish two reference jobs first.');
        j=numel(state.observers)+1;assert(j<=2,'alpha25fine:Observed','Observers finished.');
        opts=struct('domain_halfwidth_m',24,'spacing_m',.0625,'order',state.orders(j));
        [mb,~]=memory;r=observe_deepsea3d_fe_cancellation(state.cfg,opts,state.runs(j));[ma,~]=memory;
        state.observers=local_append(state.observers,r);
        state.observer_snapshots=local_append(state.observer_snapshots,struct('before_MiB',mb.MemUsedMATLAB/2^20,'after_MiB',ma.MemUsedMATLAB/2^20));
        fprintf('FINE observer N%d max indicator %.6g matches %.3g/%.3g wall %.3f\n', ...
            state.orders(j),max(r.conditions),r.term_match_error,r.norm_match_error,r.wall_s);
    case "audit"
        assert(numel(state.runs)==2 && numel(state.observers)==2,'alpha25fine:Incomplete','Finish reference and observers.');
        assert(isequal(state.cfg,state.coarse.cfg));
        a=state.runs(1);b=state.runs(2);
        c16=state.coarse.runs([state.coarse.runs.order]==16);
        c20=state.coarse.runs([state.coarse.runs.order]==20);
        state.guards=struct('fine_order',abs(a.rough_correction-b.rough_correction)/abs(b.rough_correction), ...
            'last_two',abs(b.partials(end)-b.partials(end-2))/abs(b.rough_correction), ...
            'grid_N16',abs(c16.rough_correction-a.rough_correction)/abs(a.rough_correction), ...
            'grid_N20',abs(c20.rough_correction-b.rough_correction)/abs(b.rough_correction), ...
            'boundary_N20',b.boundary_residual_relative,'halfgrid_N20',b.offgrid_boundary_residual_relative);
        state.coefficient_tail_decreases=all(diff(b.boundary_coefficient_norms(end-4:end))<0);
        state.pressure_tail_decreases=all(diff(abs(b.terms(end-4:end)))<0);
        state.coefficient_tail_decreases_N16=all(diff(a.boundary_coefficient_norms(end-4:end))<0);
        state.pressure_tail_decreases_N16=all(diff(abs(a.terms(end-4:end)))<0);
        state.reference_status="INCONCLUSIVE";
        if all(structfun(@(v)v<1e-4,state.guards)) && state.coefficient_tail_decreases && state.pressure_tail_decreases
            state.reference_status="FINE_REFERENCE_CHECK_PASS";
        end
        assert(all([state.runs.forcing_series_relative]<1e-8));
        assert(all([state.observers.term_match_error]<1e-12) && all([state.observers.norm_match_error]<1e-12));
        evaluation=struct([]);
        for r=state.runs
            T=r.max_eta_kz*sum(abs(state.cfg.profile.amplitudes_m))/r.max_height_m;
            bound=sum(r.boundary_coefficient_norms)*exp(T)*T^(r.residual_order+1)/factorial(r.residual_order+1);
            evaluation=local_append(evaluation,struct('N',r.order,'forcing_error',r.forcing_series_relative, ...
                'finite_Fourier_exponential_tail_abs_bound',bound));
        end
        writetable(struct2table(evaluation),fullfile(out,'evaluation_checks.csv'));
        checks=struct('complete',true,'reference_status',state.reference_status,'alpha25_G1',"INCONCLUSIVE", ...
            'global_G1',"OPEN",'WP2_allowed',false,'N_max',20,'FK_solves',0,'domain_extensions',0, ...
            'coefficient_tail_decreases',state.coefficient_tail_decreases,'pressure_tail_decreases',state.pressure_tail_decreases);
        fid=fopen(fullfile(out,'artifact_assertions.json'),'w');fprintf(fid,'%s',jsonencode(checks,PrettyPrint=true));fclose(fid);
        writetable(struct2table(state.guards),fullfile(out,'guards.csv'));
        cost=struct('FE_wall_s',sum([state.runs.wall_s]),'FE_cpu_s',sum([state.runs.cpu_s]), ...
            'observer_wall_s',sum([state.observers.wall_s]),'observer_cpu_s',sum([state.observers.cpu_s]), ...
            'scope','saved FE and diagnostic replay; excludes IO, static analysis and reporting');
        fid=fopen(fullfile(out,'cost_summary.json'),'w');fprintf(fid,'%s',jsonencode(cost,PrettyPrint=true));fclose(fid);
        fprintf('FINE STATUS %s\n',state.reference_status);disp(state.guards);
    otherwise
        error('alpha25fine:Stage','Unknown stage.');
end
save(file,'state','-v7');local_csv(state,out);
end
function x=local_append(x,y)
if isempty(x),x=y;else,x(end+1)=y;end
end
function local_csv(state,out)
rows=struct([]);terms=struct([]);cancels=struct([]);
for i=1:numel(state.runs)
    a=state.runs(i);mem=state.snapshots(i);
    row=struct('N',a.order,'L',24,'h',.0625,'M',40,'delta_real',real(a.rough_correction), ...
        'delta_imag',imag(a.rough_correction),'boundary',a.boundary_residual_relative, ...
        'halfgrid',a.offgrid_boundary_residual_relative,'wall_s',a.wall_s,'cpu_s',a.cpu_s, ...
        'spectral_MiB',a.stored_spectral_MiB,'process_before_MiB',mem.before_MiB,'process_after_MiB',mem.after_MiB, ...
        'pressure_sum_condition',sum(abs(a.terms))/abs(a.rough_correction));
    rows=local_append(rows,row);
    c=state.coarse.runs([state.coarse.runs.order]==a.order);
    for n=1:a.order
        row=struct('N',a.order,'n',n,'term_real',real(a.terms(n)),'term_imag',imag(a.terms(n)), ...
            'partial_real',real(a.partials(n)),'partial_imag',imag(a.partials(n)), ...
            'boundary_coefficient_norm',a.boundary_coefficient_norms(n), ...
            'coarse_L2_approx',.125*c.boundary_coefficient_norms(n),'fine_L2_approx',.0625*a.boundary_coefficient_norms(n));
        terms=local_append(terms,row);
    end
end
for i=1:numel(state.observers)
    a=state.observers(i);
    for n=1:numel(a.conditions)
        row=struct('N',state.orders(i),'n',n,'coefficient_norm',a.bnorm(n), ...
            'contribution_norm_sum',a.numerators(n),'cancellation_indicator',a.conditions(n));
        cancels=local_append(cancels,row);
    end
end
if ~isempty(rows),writetable(struct2table(rows),fullfile(out,'runs.csv'));end
if ~isempty(terms),writetable(struct2table(terms),fullfile(out,'terms.csv'));end
if ~isempty(cancels),writetable(struct2table(cancels),fullfile(out,'internal_cancellation.csv'));end
obs=struct([]);
for i=1:numel(state.observers)
    a=state.observers(i);mem=state.observer_snapshots(i);
    row=struct('N',state.orders(i),'max_indicator',max(a.conditions), ...
        'term_match_error',a.term_match_error,'norm_match_error',a.norm_match_error, ...
        'wall_s',a.wall_s,'cpu_s',a.cpu_s,'process_before_MiB',mem.before_MiB,'process_after_MiB',mem.after_MiB);
    obs=local_append(obs,row);
end
if ~isempty(obs),writetable(struct2table(obs),fullfile(out,'observers.csv'));end
end
