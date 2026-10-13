function state=validate_deepsea3d_g1_q10(stage)
%VALIDATE_DEEPSEA3D_G1_Q10 Single horizontal-geometry change, gated physics.
arguments
    stage (1,1) string = "status"
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
out=fullfile(root,'results','validation','deepsea3d_g1_q10');
file=fullfile(out,'q10.mat');
if stage=="init"
    assert(~isfile(file),'deepsea3dQ10:Exists','Use saved checkpoint or a separate checkout.');
    if ~isfolder(out),mkdir(out);end
    a=load(fullfile(root,'results','validation','deepsea3d_alpha25_domain','domain.mat'),'state');
    assert(a.state.complete && a.state.G1_status=="CONDITIONAL_G1_PASS");
    cfg=a.state.cfg;cfg.tx_m=[150.5 0 1500];
    state=struct('cfg',cfg,'baseline',a.state,'phase',"preflight",'jobs',[24 .125 12;24 .125 16], ...
        'runs',struct([]),'preflight_guards',struct(),'selected_N',0,'coarse_index',0,'fine_index',0, ...
        'reference_index',0,'reference_ready',false,'reference_guards',struct(), ...
        'observer',struct(),'observer_memory',struct(), ...
        'fk_settings',[16 .05 0;16 .025 0;16 .05 1;12 .025 0;24 .025 0], ...
        'fk',struct([]),'fk_memory',struct([]),'fk_guards',struct(),'comparison',struct(), ...
        'complete',false,'G1_status',"INCONCLUSIVE",'reason',"Reference pending", ...
        'global_G1',"OPEN",'WP2_allowed',false,'matlab_version',version);
    state.flat_control=local_flat_control(cfg,root);
    restored=cfg;restored.tx_m=a.state.cfg.tx_m;assert(isequal(restored,a.state.cfg));
    assert(isequal(cfg.center_m,[0 0]) && isequal(cfg.rx_m,[0 0 5]));
    assert(isequal(cfg.profile.amplitudes_m,[.05 .025]) && cfg.frequency_hz==6000);
else
    a=load(file,'state');state=a.state;
end
switch stage
    case "init"
    case "next"
        assert(ismember(state.phase,["preflight","upgrade","fine","domain"]), ...
            'deepsea3dQ10:Phase','Use observe/fk/audit after reference jobs.');
        j=numel(state.runs)+1;setting=state.jobs(j,:);
        opts=struct('domain_halfwidth_m',setting(1),'spacing_m',setting(2), ...
            'order',setting(3),'residual_order',40,'offgrid',true);
        [mb,~]=memory;r=deepsea3d_dirichlet_fe(state.cfg,opts);[ma,~]=memory;
        entry=struct('settings',setting,'phase',state.phase,'result',r,'memory',local_memory(mb,ma));
        state.runs=local_append(state.runs,entry);
        fprintf('Q10 %s L%g h%g N%d delta=(%.12g,%.12g) boundary %.3g half %.3g wall %.3fs\n', ...
            state.phase,setting(1),setting(2),setting(3),real(r.rough_correction),imag(r.rough_correction), ...
            r.boundary_residual_relative,r.offgrid_boundary_residual_relative,r.wall_s);
        state=local_schedule(state);
    case "observe"
        assert(state.phase=="observe",'deepsea3dQ10:Phase','Finish reference jobs first.');
        a=state.runs(state.reference_index);v=a.settings;
        opts=struct('domain_halfwidth_m',v(1),'spacing_m',v(2),'order',v(3));
        [mb,~]=memory;r=observe_deepsea3d_fe_cancellation(state.cfg,opts,a.result);[ma,~]=memory;
        state.observer=r;state.observer_memory=local_memory(mb,ma);
        if local_finite(r) && r.term_match_error<1e-12 && r.norm_match_error<1e-12
            state.reference_ready=true;state.phase="fk";
        else
            state=local_stop(state,"Diagnostic replay unresolved");
        end
        fprintf('Q10 observer max %.6g match %.3g/%.3g; ready=%d wall %.3fs\n', ...
            max(r.conditions),r.term_match_error,r.norm_match_error,state.reference_ready,r.wall_s);
    case "fk"
        assert(state.reference_ready && state.phase=="fk",'deepsea3dQ10:ReferenceGate','FullFK requires stable reference.');
        j=numel(state.fk)+1;v=state.fk_settings(j,:);cfg=state.cfg;
        cfg.half_width_m=v(1);cfg.spacing_m=v(2);
        if v(3),cfg.quadrature='gauss4';else,cfg.quadrature='trap';end
        [mb,~]=memory;r=solve_deepsea3d_fullfk(cfg);[ma,~]=memory;
        state.fk=local_append(state.fk,r);state.fk_memory=local_append(state.fk_memory,local_memory(mb,ma));
        fprintf('Q10 FK A%g h%g %s delta=(%.12g,%.12g) wall %.3fs\n', ...
            v(1),v(2),cfg.quadrature,real(r.rough_correction),imag(r.rough_correction),r.wall_s);
        if j==5
            state.fk_guards=local_fk_guards(state.fk);
            state.comparison=summarize_deepsea3d_q10_comparison(state);
            g=state.fk_guards;m=state.comparison.metrics_q10;
            [state.G1_status,state.reason]=classify_deepsea3d_amplitude(true, ...
                max([g.grid,g.aperture_12_16,g.aperture_16_24]),g.gauss,m);
            if ~all(arrayfun(@local_finite,state.fk)),state.G1_status="INCONCLUSIVE";state.reason="Nonfinite FullFK";end
            state.complete=true;state.phase="done";
        end
    case "audit"
        assert(state.complete,'deepsea3dQ10:Incomplete','Finish scoped computations.');
        unchanged=state.cfg;unchanged.tx_m=state.baseline.cfg.tx_m;
        assert(isequal(unchanged,state.baseline.cfg) && ~state.WP2_allowed && state.global_G1=="OPEN");
        assert(state.flat_control.passed);
        if state.reference_ready
            assert(all(structfun(@(v)isfinite(v) && v<1e-4,state.reference_guards)));
            assert(local_good(state.runs(state.reference_index).result));
            assert(numel(state.fk)==5);
            numerical=local_fk_guards(state.fk);names=fieldnames(state.fk_guards);
            for j=1:numel(names),assert(isequaln(numerical.(names{j}),state.fk_guards.(names{j})));end
            state.fk_guards=numerical; % Add saved-pressure diagnostics without rerunning a kernel.
            comparison=summarize_deepsea3d_q10_comparison(state);assert(isequaln(comparison,state.comparison));
            g=state.fk_guards;m=comparison.metrics_q10;
            assert(classify_deepsea3d_amplitude(true,max([g.grid,g.aperture_12_16,g.aperture_16_24]),g.gauss,m)==state.G1_status);
        else
            assert(isempty(state.fk) && state.G1_status=="INCONCLUSIVE");
        end
        checks=struct('complete',true,'flat_control_pass',true,'reference_ready',state.reference_ready, ...
            'q10_G1',state.G1_status,'global_G1',"OPEN",'WP2_allowed',false, ...
            'FE_new_solves',numel(state.runs),'FK_solves',numel(state.fk),'q0_new_solves',0, ...
            'geometry_only_change',true,'surface_center_unchanged',true,'no_phase_fit',true);
        local_json(fullfile(out,'artifact_assertions.json'),checks);
        fprintf('Q10 AUDIT %s\n',state.G1_status);
    case "status"
        fprintf('Q10 phase=%s selected N%d FE=%d FK=%d ready=%d complete=%d G1=%s\n', ...
            state.phase,state.selected_N,numel(state.runs),numel(state.fk),state.reference_ready,state.complete,state.G1_status);return
    otherwise
        error('deepsea3dQ10:Stage','Unknown stage.');
end
save(file,'state','-v7');local_outputs(state,out);
end
function s=local_schedule(s)
j=numel(s.runs);a=s.runs(j);r=a.result;
if s.phase=="preflight" || s.phase=="upgrade"
    if j==1,return;end
    if s.phase=="upgrade",previous=s.runs(s.coarse_index).result;else,previous=s.runs(j-1).result;end
    g=struct('order_change',local_delta(previous,r),'last_two',local_order(r,2), ...
        'boundary',r.boundary_residual_relative,'halfgrid',r.offgrid_boundary_residual_relative);
    s.preflight_guards=g;
    if local_good(r) && all(structfun(@(v)isfinite(v) && v<1e-4,g))
        s.selected_N=r.order;s.coarse_index=j;s.phase="fine";
        s.jobs(end+1,:)=[24 .0625 r.order];
    elseif r.order<20
        s.jobs(end+1,:)=[24 .125 20];
    else
        s=local_stop(s,"Coarse reference unresolved at maximum N20");
    end
elseif s.phase=="fine"
    g=local_order_guards(r);g.grid=local_delta(s.runs(s.coarse_index).result,r);
    if local_good(r) && all(structfun(@(v)isfinite(v) && v<1e-4,g))
        s.fine_index=j;s.phase="domain";s.jobs(end+1,:)=[32 .0625 r.order];
    elseif r.order<20
        s.phase="upgrade";s.jobs(end+1,:)=[24 .125 20];
    else
        s=local_stop(s,"Fine reference order/grid/true boundary unresolved at N20");
    end
elseif s.phase=="domain"
    g=local_order_guards(r);g.coarse_order=s.preflight_guards.order_change;
    g.grid=local_delta(s.runs(s.coarse_index).result,s.runs(s.fine_index).result);
    g.domain=local_delta(s.runs(j-1).result,r);s.reference_guards=g;
    other=rmfield(g,'domain');
    if local_good(r) && all(structfun(@(v)isfinite(v) && v<1e-4,other))
        if g.domain<1e-4
            s.reference_index=j;s.phase="observe";
        elseif a.settings(1)==32
            s.jobs(end+1,:)=[48 .0625 r.order];
        else
            s=local_stop(s,"Domain unresolved at L48");
        end
    else
        s=local_stop(s,"Expanded-domain reference order/true boundary/tail unresolved");
    end
end
end
function g=local_order_guards(r)
g=struct('order_Nminus4',local_order(r,4),'last_two',local_order(r,2), ...
    'boundary',r.boundary_residual_relative,'halfgrid',r.offgrid_boundary_residual_relative);
end
function v=local_order(r,d)
v=abs(r.partials(end)-r.partials(end-d))/abs(r.rough_correction);
end
function yes=local_good(r)
yes=local_finite(r) && r.forcing_series_relative<1e-8 && ...
    all(diff(r.boundary_coefficient_norms(end-4:end))<0) && all(diff(abs(r.terms(end-4:end)))<0);
end
function yes=local_finite(s)
yes=true;names=fieldnames(s);
for j=1:numel(names),v=s.(names{j});if isnumeric(v),yes=yes && all(isfinite(v(:)));end,end
end
function s=local_stop(s,reason)
s.phase="done";s.complete=true;s.G1_status="INCONCLUSIVE";s.reason=reason;
end
function d=local_delta(a,b)
d=abs(a.rough_correction-b.rough_correction)/abs(b.rough_correction);
end
function g=local_fk_guards(f)
g=struct('grid',local_delta(f(1),f(2)),'gauss',local_delta(f(3),f(2)), ...
    'aperture_12_16',local_delta(f(4),f(2)),'aperture_16_24',local_delta(f(2),f(5)), ...
    'reflect_grid',abs(f(1).reflect_native-f(2).reflect_native)/abs(f(2).reflect_native), ...
    'reflect_gauss',abs(f(3).reflect_native-f(2).reflect_native)/abs(f(2).reflect_native), ...
    'total_grid',abs(f(1).total_native-f(2).total_native)/abs(f(2).total_native), ...
    'total_gauss',abs(f(3).total_native-f(2).total_native)/abs(f(2).total_native));
end
function a=local_append(a,b)
if isempty(a),a=b;else,a(end+1)=b;end
end
function s=local_memory(a,b)
s=struct('before_MiB',a.MemUsedMATLAB/2^20,'after_MiB',b.MemUsedMATLAB/2^20);
end
function local_json(file,s)
fid=fopen(file,'w');assert(fid>=0);clean=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(s,PrettyPrint=true));
end
function f=local_flat_control(cfg,root)
saved=readtable(fullfile(root,'results','validation','deepsea3d_wp0_wp1','wp0.csv'));
w=saved(saved.q==.1 & saved.frequency_hz==6000 & saved.height_m==0,:);assert(height(w)==1);
Ldir=sqrt(150.5^2+1495^2);Lref=sqrt(150.5^2+1505^2);k=2*pi*6000/1500;
pd=exp(1i*k*Ldir)/(4*pi*Ldir);pr=-exp(1i*k*Lref)/(4*pi*Lref);
b=deepsea3d_flat(cfg,0);image=cfg.tx_m;image(3)=-image(3);
direct=deepsea3d_point_field(cfg.rx_m,cfg.tx_m,k);reflected=-deepsea3d_point_field(cfg.rx_m,image,k);
f=struct('q',150.5/1505,'Ld_m',Ldir,'Lr_m',Lref,'delay_s',(Lref-Ldir)/1500, ...
    'direct_real',real(pd),'direct_imag',imag(pd),'reflect_real',real(pr),'reflect_imag',imag(pr), ...
    'analytic_direct_error',abs(pd-b.direct_native)/abs(pd),'analytic_reflect_error',abs(pr-b.reflect_native)/abs(pr), ...
    'point_direct_error',abs(direct-pd)/abs(pd),'point_reflect_error',abs(reflected-pr)/abs(pr), ...
    'saved_WP0_direct_error',abs(pd-complex(w.direct_real,w.direct_imag))/abs(pd), ...
    'saved_WP0_reflect_error',abs(pr-complex(w.reflect_real,w.reflect_imag))/abs(pr), ...
    'DSP_direct_error',abs(b.direct_dsp-1/(4*pi*Ldir))/abs(pd), ...
    'DSP_reflect_error',abs(b.reflect_dsp+exp(-1i*k*(Lref-Ldir))/(4*pi*Lref))/abs(pr), ...
    'specular_x_m',150.5*5/1505,'surface_center_x_m',cfg.center_m(1),'passed',false);
names=fieldnames(f);errors=names(contains(names,'error'));f.passed=all(cellfun(@(name)f.(name)<1e-10,errors));
assert(f.passed,'deepsea3dQ10:FlatControl','New geometry flat control failed.');
end
function local_outputs(s,out)
writetable(struct2table(s.flat_control),fullfile(out,'flat_controls.csv'));
rows=struct([]);terms=struct([]);evaluation=struct([]);
for j=1:numel(s.runs)
    a=s.runs(j);r=a.result;v=a.settings;mem=a.memory;
    rows=local_append(rows,struct('job',j,'phase',a.phase,'L',v(1),'h',v(2),'N',v(3),'M',40, ...
        'delta_real',real(r.rough_correction),'delta_imag',imag(r.rough_correction), ...
        'reflect_real',real(r.reflect_native),'reflect_imag',imag(r.reflect_native), ...
        'total_real',real(r.total_native),'total_imag',imag(r.total_native), ...
        'Nminus4_change',local_order(r,4),'last_two_change',local_order(r,2), ...
        'boundary',r.boundary_residual_relative,'halfgrid',r.offgrid_boundary_residual_relative, ...
        'coefficient_tail_decreases',all(diff(r.boundary_coefficient_norms(end-4:end))<0), ...
        'pressure_tail_decreases',all(diff(abs(r.terms(end-4:end)))<0), ...
        'pressure_sum_indicator',sum(abs(r.terms))/abs(r.rough_correction), ...
        'wall_s',r.wall_s,'cpu_s',r.cpu_s,'spectral_MiB',r.stored_spectral_MiB, ...
        'process_before_MiB',mem.before_MiB,'process_after_MiB',mem.after_MiB));
    for n=1:r.order
        terms=local_append(terms,struct('job',j,'n',n,'term_real',real(r.terms(n)), ...
            'term_imag',imag(r.terms(n)),'partial_real',real(r.partials(n)), ...
            'partial_imag',imag(r.partials(n)),'boundary_coefficient_norm',r.boundary_coefficient_norms(n), ...
            'h_normalized_L2',r.spacing_m*r.boundary_coefficient_norms(n)));
    end
    T=r.max_eta_kz*sum(abs(s.cfg.profile.amplitudes_m))/r.max_height_m;
    bound=sum(r.boundary_coefficient_norms)*exp(T)*T^(41)/factorial(41);
    evaluation=local_append(evaluation,struct('job',j,'forcing_error',r.forcing_series_relative, ...
        'finite_Fourier_exponential_tail_abs_bound',bound));
end
if ~isempty(rows),writetable(struct2table(rows),fullfile(out,'reference_runs.csv'));end
if ~isempty(terms),writetable(struct2table(terms),fullfile(out,'reference_terms.csv'));end
if ~isempty(evaluation),writetable(struct2table(evaluation),fullfile(out,'evaluation_checks.csv'));end
if ~isempty(fieldnames(s.preflight_guards)),writetable(struct2table(s.preflight_guards),fullfile(out,'preflight_guards.csv'));end
if ~isempty(fieldnames(s.reference_guards)),writetable(struct2table(s.reference_guards),fullfile(out,'reference_guards.csv'));end
if ~isempty(fieldnames(s.observer))
    r=s.observer;mem=s.observer_memory;
    writetable(table((1:numel(r.conditions))',r.bnorm,r.numerators,r.conditions, ...
        'VariableNames',{'n','coefficient_norm','contribution_norm_sum','cancellation_indicator'}),fullfile(out,'internal_cancellation.csv'));
    row=struct('max_indicator',max(r.conditions),'term_match_error',r.term_match_error,'norm_match_error',r.norm_match_error, ...
        'wall_s',r.wall_s,'cpu_s',r.cpu_s,'process_before_MiB',mem.before_MiB,'process_after_MiB',mem.after_MiB);
    writetable(struct2table(row),fullfile(out,'observer.csv'));
end
rows=struct([]);
for j=1:numel(s.fk)
    r=s.fk(j);v=s.fk_settings(j,:);mem=s.fk_memory(j);
    rows=local_append(rows,struct('A',v(1),'h',v(2),'gauss4',logical(v(3)), ...
        'reflect_real',real(r.reflect_native),'reflect_imag',imag(r.reflect_native), ...
        'delta_real',real(r.rough_correction),'delta_imag',imag(r.rough_correction), ...
        'total_real',real(r.total_native),'total_imag',imag(r.total_native), ...
        'identity_residual',r.quadrature_identity_residual,'density_cancellation_indicator',r.cancellation_condition, ...
        'nodes',r.node_count,'wall_s',r.wall_s,'cpu_s',r.cpu_s,'block_payload_MiB',r.estimated_block_payload_MiB, ...
        'process_before_MiB',mem.before_MiB,'process_after_MiB',mem.after_MiB));
end
if ~isempty(rows),writetable(struct2table(rows),fullfile(out,'fk_runs.csv'));end
if ~isempty(fieldnames(s.fk_guards)),writetable(struct2table(s.fk_guards),fullfile(out,'fk_guards.csv'));end
if ~isempty(fieldnames(s.comparison))
    c=s.comparison;
    writetable(struct2table([c.row_q0,c.row_q10]),fullfile(out,'geometry_comparison.csv'));
    writetable(struct2table(c.envelopes),fullfile(out,'error_envelopes.csv'));
    local_json(fullfile(out,'comparison_summary.json'),c);
end
FEwall=0;FEcpu=0;for j=1:numel(s.runs),FEwall=FEwall+s.runs(j).result.wall_s;FEcpu=FEcpu+s.runs(j).result.cpu_s;end
Owall=0;Ocpu=0;if ~isempty(fieldnames(s.observer)),Owall=s.observer.wall_s;Ocpu=s.observer.cpu_s;end
FKwall=0;FKcpu=0;for j=1:numel(s.fk),FKwall=FKwall+s.fk(j).wall_s;FKcpu=FKcpu+s.fk(j).cpu_s;end
cost=struct('FE_wall_s',FEwall,'FE_cpu_s',FEcpu,'observer_wall_s',Owall,'observer_cpu_s',Ocpu, ...
    'FK_wall_s',FKwall,'FK_cpu_s',FKcpu,'new_wall_s',FEwall+Owall+FKwall,'q0_new_solves',0, ...
    'scope','saved new FE/replay/FK kernels; excludes IO/tests/reporting; memory snapshots are not peaks');
local_json(fullfile(out,'cost_summary.json'),cost);
end
