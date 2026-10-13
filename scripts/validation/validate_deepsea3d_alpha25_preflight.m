function state=validate_deepsea3d_alpha25_preflight(stage)
%VALIDATE_DEEPSEA3D_ALPHA25_PREFLIGHT Low-cost only; no FK or domain sweep.
arguments
    stage (1,1) string = "next"
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
out=fullfile(root,'results','validation','deepsea3d_alpha25_preflight');
file=fullfile(out,'preflight.mat');
if stage=="init"
    assert(~isfile(file),'alpha25:Exists','Use next/audit on existing checkpoint.');
    if ~isfolder(out),mkdir(out);end
    a=load(fullfile(root,'results','validation','deepsea3d_g1_amplitude','amplitude.mat'),'state');
    c=a.state.cases([a.state.cases.alpha]==2);
    assert(c.done && c.reference_ready && c.status=="CONDITIONAL_G1_PASS");
    state=struct('cfg',a.state.cfg,'alpha2',c,'base_cfg',a.state.cfg, ...
        'orders',[12 16],'runs',struct([]),'snapshots',struct([]),'complete',false, ...
        'precheck_status',"PENDING",'G1_status',"INCONCLUSIVE",'WP2_allowed',false);
    state.cfg.profile.amplitudes_m=state.base_cfg.profile.amplitudes_m*2.5;
else
    a=load(file,'state');state=a.state;
end
switch stage
    case "init"
    case "next"
        j=numel(state.runs)+1;
        assert(j<=numel(state.orders),'alpha25:Finished','Use audit.');
        N=state.orders(j); opts=struct('domain_halfwidth_m',24,'spacing_m',.125, ...
            'order',N,'residual_order',40,'offgrid',true);
        [mb,~]=memory; r=deepsea3d_dirichlet_fe(state.cfg,opts);[ma,~]=memory;
        state.runs=local_append(state.runs,r);
        state.snapshots=local_append(state.snapshots,struct('before_MiB',mb.MemUsedMATLAB/2^20, ...
            'after_MiB',ma.MemUsedMATLAB/2^20));
        fprintf('ALPHA2.5 N%d delta=(%.9g,%.9g) boundary %.4g offgrid %.4g wall %.3f\n', ...
            N,real(r.rough_correction),imag(r.rough_correction),r.boundary_residual_relative, ...
            r.offgrid_boundary_residual_relative,r.wall_s);
        if j==2
            delta=abs(state.runs(1).rough_correction-r.rough_correction)/abs(r.rough_correction);
            tail=abs(r.partials(end)-r.partials(end-2))/abs(r.rough_correction);
            if max([delta tail r.boundary_residual_relative r.offgrid_boundary_residual_relative])>=1e-4
                state.orders=[12 16 20];
            else
                state.complete=true;
            end
        elseif j==3
            state.complete=true;
        end
    case "cancellation"
        assert(state.complete,'alpha25:Incomplete','Finish planned reference jobs.');
        opts=struct('domain_halfwidth_m',24,'spacing_m',.125,'order',state.runs(end).order);
        state.cancellation=observe_deepsea3d_fe_cancellation(state.cfg,opts,state.runs(end));
        row=table((1:opts.order)',state.cancellation.bnorm,state.cancellation.numerators, ...
            state.cancellation.conditions,'VariableNames',{'order','coefficient_norm','contribution_norm_sum','cancellation_indicator'});
        writetable(row,fullfile(out,'internal_cancellation.csv'));
    case "audit"
        assert(state.complete,'alpha25:Incomplete','Finish planned reference jobs.');
        assert(isfield(state,'cancellation'),'alpha25:MissingObservation','Run cancellation stage.');
        r=state.runs(end); previous=state.runs(end-1);
        state.guards=struct('order_change',abs(previous.rough_correction-r.rough_correction)/abs(r.rough_correction), ...
            'last_two_change',abs(r.partials(end)-r.partials(end-2))/abs(r.rough_correction), ...
            'boundary',r.boundary_residual_relative,'offgrid',r.offgrid_boundary_residual_relative);
        state.coefficient_tail_decreases=all(diff(r.boundary_coefficient_norms(end-4:end))<0);
        state.pressure_tail_decreases=all(diff(abs(r.terms(end-4:end)))<0);
        assert(r.forcing_series_relative<1e-8);
        state.precheck_status="INCONCLUSIVE";
        if all(structfun(@(x)x<1e-4,state.guards)) && state.coefficient_tail_decreases && state.pressure_tail_decreases
            state.precheck_status="PROCEED_FORMAL_G1";
        end
        assert(isequal(state.cfg.profile.amplitudes_m,state.base_cfg.profile.amplitudes_m*2.5));
        unscaled=state.cfg;unscaled.profile.amplitudes_m=state.base_cfg.profile.amplitudes_m;
        assert(isequal(unscaled,state.base_cfg));
        checks=struct('complete',true,'precheck_status',state.precheck_status,'G1_status',"INCONCLUSIVE", ...
            'full_G1_domain',"OPEN",'WP2_allowed',false,'FK_solves',0,'new_reference_runs',numel(state.runs));
        fid=fopen(fullfile(out,'artifact_assertions.json'),'w');fprintf(fid,'%s',jsonencode(checks,PrettyPrint=true));fclose(fid);
        writetable(struct2table(state.guards),fullfile(out,'guards.csv'));
        compared=struct([]);
        for j=1:2
            x=state.alpha2.runs(j).result;y=state.runs(j);
            assert(x.order==y.order && state.alpha2.runs(j).settings(1)==24 ...
                && state.alpha2.runs(j).settings(2)==.125);
            row=struct('N',x.order,'alpha2_boundary',x.boundary_residual_relative, ...
                'alpha25_boundary',y.boundary_residual_relative, ...
                'alpha2_pressure_condition',sum(abs(x.terms))/abs(x.rough_correction), ...
                'alpha25_pressure_condition',sum(abs(y.terms))/abs(y.rough_correction), ...
                'alpha2_wall_s',x.wall_s,'alpha25_wall_s',y.wall_s);
            compared=local_append(compared,row);
        end
        writetable(struct2table(compared),fullfile(out,'alpha2_comparison.csv'));
        cost=struct('FE_wall_s',sum([state.runs.wall_s]),'FE_cpu_s',sum([state.runs.cpu_s]), ...
            'observer_wall_s',state.cancellation.wall_s,'observer_cpu_s',state.cancellation.cpu_s, ...
            'scope','saved kernels and diagnostic replay; excludes IO, tests, static analysis and reports');
        fid=fopen(fullfile(out,'cost_summary.json'),'w');fprintf(fid,'%s',jsonencode(cost,PrettyPrint=true));fclose(fid);
    otherwise
        error('alpha25:Stage','Unknown stage.');
end
save(file,'state','-v7');local_csv(state,out);
end
function x=local_append(x,y)
if isempty(x),x=y;else,x(end+1)=y;end
end
function local_csv(state,out)
rows=struct([]); terms=struct([]);
for i=1:numel(state.runs)
    a=state.runs(i); mem=state.snapshots(i);
    row=struct('alpha',2.5,'N',a.order,'L',24,'h',.125,'M',40, ...
        'delta_real',real(a.rough_correction),'delta_imag',imag(a.rough_correction), ...
        'boundary',a.boundary_residual_relative,'offgrid',a.offgrid_boundary_residual_relative, ...
        'forcing',a.forcing_series_relative,'wall_s',a.wall_s,'cpu_s',a.cpu_s, ...
        'spectral_MiB',a.stored_spectral_MiB,'process_before_MiB',mem.before_MiB, ...
        'process_after_MiB',mem.after_MiB,'pressure_sum_condition',sum(abs(a.terms))/abs(a.rough_correction));
    rows=local_append(rows,row);
    for n=1:a.order
        row=struct('N',a.order,'n',n,'term_real',real(a.terms(n)),'term_imag',imag(a.terms(n)), ...
            'partial_real',real(a.partials(n)),'partial_imag',imag(a.partials(n)), ...
            'boundary_coefficient_norm',a.boundary_coefficient_norms(n));
        terms=local_append(terms,row);
    end
end
if ~isempty(rows),writetable(struct2table(rows),fullfile(out,'runs.csv'));end
if ~isempty(terms),writetable(struct2table(terms),fullfile(out,'terms.csv'));end
end
