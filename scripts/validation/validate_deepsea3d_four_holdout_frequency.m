function state=validate_deepsea3d_four_holdout_frequency(stage,frequency_hz)
%VALIDATE_DEEPSEA3D_FOUR_HOLDOUT_FREQUENCY Four prelocked half-kHz FE references.
% Research workflow copy of the frozen 4/8 kHz gate; no solver changes.
arguments
    stage (1,1) string = "status"
    frequency_hz (1,1) double = 0
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
out=fullfile(root,'results','validation','deepsea3d_four_holdout_20261026');file=fullfile(out,'four_holdout.mat');
if stage=="init"
    assert(~isfile(file),'deepsea3dFrequency:Exists','Use saved checkpoint or a separate checkout.');
    a=load(fullfile(root,'results','validation','deepsea3d_g1_azimuth','azimuth.mat'),'state');
    assert(a.state.complete && a.state.G1_status=="CONDITIONAL_G1_PASS");
    assert(isequal(a.state.cfg.tx_m,[0 150.5 1500]) && a.state.cfg.frequency_hz==6000);
    state=struct('baseline6',a.state,'cases',struct([]),'global_G1',"OPEN",'WP2_allowed',false,'matlab_version',version);
    for f=[4500 5500 6500 7500]
        cfg=a.state.cfg;cfg.frequency_hz=f;c=local_case(cfg);c.flat_control=local_flat_control(cfg,root);
        state.cases=local_append(state.cases,c);
    end
else
    a=load(file,'state');state=a.state;
    if ~isfield(state.cases,'fk_refinement')
        for j=1:numel(state.cases),state.cases(j).fk_refinement=struct();state.cases(j).fk_refinement_memory=struct();end
    end
end
if stage=="status"
    for c=state.cases,fprintf('%g Hz phase %s N%d FE%d FK%d ready%d complete%d %s\n', ...
        c.frequency_hz,c.phase,c.selected_N,numel(c.runs),numel(c.fk),c.reference_ready,c.complete,c.G1_status);end
    return
elseif stage=="audit"
    assert(all([state.cases.complete]),'deepsea3dFrequency:Incomplete','Finish four bounded cases.');
    for c=state.cases
        restored=c.cfg;restored.frequency_hz=6000;assert(isequal(restored,state.baseline6.cfg));
        assert(c.flat_control.passed);
        if c.reference_ready
            assert(all(structfun(@(v)isfinite(v) && v<1e-4,c.reference_guards)));
            assert(local_good(c.runs(c.reference_index).result) && numel(c.fk)==5);
            assert(isequaln(local_case_fk_guards(c),c.fk_guards));
            m=local_metrics(c.fk(2),c.runs(c.reference_index).result);assert(isequaln(m,c.metrics));
            g=c.fk_guards;assert(classify_deepsea3d_amplitude(true,max([g.grid,g.aperture_12_16,g.aperture_16_24]),g.gauss,m)==c.G1_status);
        else
            assert(isempty(c.fk) && c.G1_status=="INCONCLUSIVE");
        end
    end
    local_json(fullfile(out,'artifact_assertions.json'),struct('complete',true, ...
        'frequency_hz',[4500 5500 6500 7500],'status',[state.cases.G1_status], ...
        'six_kHz_new_solves',0,'frequency_only_change',true,'no_phase_fit',true,'global_G1',"OPEN",'WP2_allowed',false));
elseif stage~="init"
    if frequency_hz==0,ix=find(~[state.cases.complete],1);else,ix=find([state.cases.frequency_hz]==frequency_hz,1);end
    assert(~isempty(ix),'deepsea3dFrequency:Finished','Use audit/status.');c=state.cases(ix);
    switch stage
        case "next"
            assert(ismember(c.phase,["preflight","upgrade","fine","domain"]),'deepsea3dFrequency:Phase','Use evaluate/observe/fk.');
            j=numel(c.runs)+1;v=c.jobs(j,:);
            opts=struct('domain_halfwidth_m',v(1),'spacing_m',v(2),'order',v(3),'residual_order',40,'offgrid',true);
            mb=struct('MemUsedMATLAB',NaN);r=deepsea3d_dirichlet_fe(c.cfg,opts);ma=struct('MemUsedMATLAB',NaN);
            c.runs=local_append(c.runs,struct('settings',v,'phase',c.phase,'result',r,'memory',local_memory(mb,ma)));
            fprintf('%g Hz %s L%g h%g N%d res %.3g half %.3g wall %.3fs\n',c.frequency_hz,c.phase,v(1),v(2),v(3),r.boundary_residual_relative,r.offgrid_boundary_residual_relative,r.wall_s);
            c=local_schedule(c);
        case "evaluate"
            assert(c.phase=="evaluation",'deepsea3dFrequency:Phase','Finish bounded reference candidate.');
            a=c.runs(c.reference_index);v=a.settings;
            opts=struct('domain_halfwidth_m',v(1),'spacing_m',v(2),'order',v(3),'residual_order',32,'offgrid',true);
            mb=struct('MemUsedMATLAB',NaN);r=deepsea3d_dirichlet_fe(c.cfg,opts);ma=struct('MemUsedMATLAB',NaN);
            c.evaluation=struct('settings',v,'result',r,'memory',local_memory(mb,ma));
            c.evaluation_controls=local_evaluation(c.cfg,a.result,r,v);
            q=c.evaluation_controls;g=c.reference_guards;
            g.M_boundary_change=q.boundary_change;g.M_halfgrid_change=q.halfgrid_change;
            g.M_pressure_change=q.pressure_change;g.M32_boundary=r.boundary_residual_relative;
            g.M32_halfgrid=r.offgrid_boundary_residual_relative;
            g.M32_tail=q.tail32_relative;g.M40_tail=q.tail40_relative;c.reference_guards=g;
            if local_finite(r) && q.term_match<1e-12 && q.norm_match<1e-12 && all(structfun(@(z)isfinite(z) && z<1e-4,g))
                c.phase="observe";
            else
                c=local_stop(c,"M32/M40 boundary evaluation or finite-Fourier tail unresolved");
            end
            fprintf('%g Hz M32/M40 res changes %.3g/%.3g tails %.3g/%.3g phase %s\n', ...
                c.frequency_hz,q.boundary_change,q.halfgrid_change,q.tail32_relative,q.tail40_relative,c.phase);
        case "observe"
            assert(c.phase=="observe",'deepsea3dFrequency:Phase','Finish evaluation controls.');
            a=c.runs(c.reference_index);v=a.settings;opts=struct('domain_halfwidth_m',v(1),'spacing_m',v(2),'order',v(3));
            mb=struct('MemUsedMATLAB',NaN);r=observe_deepsea3d_fe_cancellation(c.cfg,opts,a.result);ma=struct('MemUsedMATLAB',NaN);
            c.observer=r;c.observer_memory=local_memory(mb,ma);
            if local_finite(r) && r.term_match_error<1e-12 && r.norm_match_error<1e-12
                c.reference_ready=true;c.phase="fk";
            else
                c=local_stop(c,"Diagnostic replay unresolved");
            end
            fprintf('%g Hz replay max %.6g match %.3g/%.3g ready%d wall %.3fs\n',c.frequency_hz,max(r.conditions),r.term_match_error,r.norm_match_error,c.reference_ready,r.wall_s);
        case "fk"
            assert(c.reference_ready && c.phase=="fk",'deepsea3dFrequency:ReferenceGate','FullFK forbidden before independent reference.');
            j=numel(c.fk)+1;v=c.fk_settings(j,:);cfg=c.cfg;cfg.half_width_m=v(1);cfg.spacing_m=v(2);
            if v(3),cfg.quadrature='gauss4';else,cfg.quadrature='trap';end
            mb=struct('MemUsedMATLAB',NaN);r=solve_deepsea3d_fullfk(cfg);ma=struct('MemUsedMATLAB',NaN);
            c.fk=local_append(c.fk,r);c.fk_memory=local_append(c.fk_memory,local_memory(mb,ma));
            fprintf('%g Hz FK A%g h%g %s wall %.3fs\n',c.frequency_hz,v(1),v(2),cfg.quadrature,r.wall_s);
            if j==5
                c.fk_guards=local_fk_guards(c.fk);c.metrics=local_metrics(c.fk(2),c.runs(c.reference_index).result);g=c.fk_guards;
                [c.G1_status,c.reason]=classify_deepsea3d_amplitude(true,max([g.grid,g.aperture_12_16,g.aperture_16_24]),g.gauss,c.metrics);
                if ~all(arrayfun(@local_finite,c.fk)),c.G1_status="INCONCLUSIVE";c.reason="Nonfinite FullFK";end
                c.complete=true;c.phase="done";fprintf('%g Hz %s\n',c.frequency_hz,c.G1_status);
            end
        case "quadrature"
            assert(c.reference_ready && numel(c.fk)==5 && c.fk_guards.gauss>=1e-3 && isempty(fieldnames(c.fk_refinement)), ...
                'deepsea3dFrequency:Quadrature','Only one bounded independent Gauss refinement after a failed coarse guard.');
            cfg=c.cfg;cfg.half_width_m=16;cfg.spacing_m=.025;cfg.quadrature='gauss4';
            mb=struct('MemUsedMATLAB',NaN);r=solve_deepsea3d_fullfk(cfg);ma=struct('MemUsedMATLAB',NaN);
            c.fk_refinement=r;c.fk_refinement_memory=local_memory(mb,ma);c.fk_guards=local_case_fk_guards(c);
            g=c.fk_guards;[c.G1_status,c.reason]=classify_deepsea3d_amplitude(true,max([g.grid,g.aperture_12_16,g.aperture_16_24]),g.gauss,c.metrics);
            if ~local_finite(r),c.G1_status="INCONCLUSIVE";c.reason="Nonfinite independent quadrature refinement";end
            fprintf('%g Hz Gauss4 h.025 relative %.6g wall %.3fs %s\n',c.frequency_hz,g.gauss,r.wall_s,c.G1_status);
        otherwise,error('deepsea3dFrequency:Stage','Unknown stage.');
    end
    state.cases(ix)=c;
end
save(file,'state','-v7');
for c=state.cases
    folder=fullfile(out,sprintf('f%d',c.frequency_hz));if ~isfolder(folder),mkdir(folder);end
    local_outputs(c,folder);
    if ~isempty(fieldnames(c.metrics)),writetable(struct2table(c.metrics),fullfile(folder,'physical_metrics.csv'));end
    if ~isempty(fieldnames(c.evaluation_controls)),writetable(struct2table(c.evaluation_controls),fullfile(folder,'M_evaluation_controls.csv'));end
    if ~isempty(fieldnames(c.fk_refinement))
        r=c.fk_refinement;mem=c.fk_refinement_memory;
        row=struct('A',16,'h',.025,'gauss4',true,'delta_real',real(r.rough_correction),'delta_imag',imag(r.rough_correction), ...
            'reflect_real',real(r.reflect_native),'reflect_imag',imag(r.reflect_native),'wall_s',r.wall_s,'cpu_s',r.cpu_s, ...
            'nodes',r.node_count,'block_payload_MiB',r.estimated_block_payload_MiB,'process_before_MiB',mem.before_MiB,'process_after_MiB',mem.after_MiB);
        writetable(struct2table(row),fullfile(folder,'gauss_refinement.csv'));
    end
    local_json(fullfile(folder,'cost_summary.json'),local_cost(c));
end
end
function c=local_case(cfg)
c=struct('frequency_hz',cfg.frequency_hz,'cfg',cfg,'flat_control',struct(), ...
    'phase',"preflight",'jobs',[24 .125 12;24 .125 16],'runs',struct([]),'preflight_guards',struct(), ...
    'selected_N',0,'coarse_index',0,'fine_index',0,'reference_index',0,'reference_ready',false,'reference_guards',struct(), ...
    'evaluation',struct(),'evaluation_controls',struct(),'observer',struct(),'observer_memory',struct(), ...
    'fk_settings',[16 .05 0;16 .025 0;16 .05 1;12 .025 0;24 .025 0], ...
    'fk',struct([]),'fk_memory',struct([]),'fk_guards',struct(),'fk_refinement',struct(),'fk_refinement_memory',struct(),'metrics',struct(),'comparison',struct(), ...
    'complete',false,'G1_status',"INCONCLUSIVE",'reason',"Reference pending");
end
function m=local_metrics(f,r)
m=struct('reflected_complex',abs(f.reflect_native-r.reflect_native)/abs(r.reflect_native), ...
    'increment_complex',abs(f.rough_correction-r.rough_correction)/abs(r.rough_correction), ...
    'reflected_magnitude',abs(abs(f.reflect_native)-abs(r.reflect_native))/abs(r.reflect_native), ...
    'increment_magnitude',abs(abs(f.rough_correction)-abs(r.rough_correction))/abs(r.rough_correction), ...
    'reflected_phase_rad',angle(f.reflect_native*conj(r.reflect_native)), ...
    'increment_phase_rad',angle(f.rough_correction*conj(r.rough_correction)));
end
function q=local_evaluation(cfg,r40,r32,v)
timer=tic;cpu0=cputime;L=v(1);h=r40.spacing_m;n=round(sqrt(r40.nodes));x=-L+(0:n-1)*h;
[xx,yy]=meshgrid(x,x);eta=deepsea3d_surface(xx,yy,cfg.profile);mask=hypot(xx,yy)<cfg.profile.support_m;
points=[xx(mask)+cfg.center_m(1),yy(mask)+cfg.center_m(2),eta(mask)];image=cfg.tx_m;image(3)=-image(3);
k=2*pi*cfg.frequency_hz/cfg.c0_mps;
den=norm(deepsea3d_point_field(points,cfg.tx_m,k)-deepsea3d_point_field(points,image,k));
T=r40.max_eta_kz*sum(abs(cfg.profile.amplitudes_m))/r40.max_height_m;
amount=sum(r40.boundary_coefficient_norms);b32=amount*exp(T)*T^33/factorial(33);b40=amount*exp(T)*T^41/factorial(41);
q=struct('M_low',32,'M_high',40,'T_bound',T,'forcing_norm',den, ...
    'tail32_abs',b32,'tail40_abs',b40,'tail32_relative',b32/den,'tail40_relative',b40/den, ...
    'boundary_change',abs(r32.boundary_residual_relative-r40.boundary_residual_relative), ...
    'halfgrid_change',abs(r32.offgrid_boundary_residual_relative-r40.offgrid_boundary_residual_relative), ...
    'pressure_change',local_delta(r32,r40),'term_match',norm(r32.terms-r40.terms)/norm(r40.terms), ...
    'norm_match',norm(r32.boundary_coefficient_norms-r40.boundary_coefficient_norms)/norm(r40.boundary_coefficient_norms), ...
    'wall_s',toc(timer),'cpu_s',cputime-cpu0);
end
function f=local_flat_control(cfg,root)

Ld=sqrt(150.5^2+1495^2);Lr=sqrt(150.5^2+1505^2);k=2*pi*cfg.frequency_hz/1500;
pd=exp(1i*k*Ld)/(4*pi*Ld);pr=-exp(1i*k*Lr)/(4*pi*Lr);b=deepsea3d_flat(cfg,0);im=cfg.tx_m;im(3)=-im(3);
f=struct('frequency_hz',cfg.frequency_hz,'Ld_m',Ld,'Lr_m',Lr,'delay_s',(Lr-Ld)/1500, ...
    'direct_real',real(pd),'direct_imag',imag(pd),'reflect_real',real(pr),'reflect_imag',imag(pr), ...
    'direct_error',abs(b.direct_native-pd)/abs(pd),'reflect_error',abs(b.reflect_native-pr)/abs(pr), ...
    'point_direct_error',abs(deepsea3d_point_field(cfg.rx_m,cfg.tx_m,k)-pd)/abs(pd), ...
    'point_reflect_error',abs(-deepsea3d_point_field(cfg.rx_m,im,k)-pr)/abs(pr), ...
    'normalization_direct_error',abs(abs(pd)*4*pi*Ld-1),'normalization_reflect_error',abs(abs(pr)*4*pi*Lr-1), ...
    'DSP_direct_error',abs(b.direct_dsp-1/(4*pi*Ld))/abs(pd), ...
    'DSP_reflect_error',abs(b.reflect_dsp+exp(-1i*k*(Lr-Ld))/(4*pi*Lr))/abs(pr),'passed',false);
names=fieldnames(f);errors=names(contains(names,'error'));f.passed=all(cellfun(@(name)f.(name)<1e-10,errors));assert(f.passed);
end
function cost=local_cost(c)
wall=0;cpu=0;for a=c.runs,wall=wall+a.result.wall_s;cpu=cpu+a.result.cpu_s;end
ev=0;evcpu=0;diagnostic=0;dcpu=0;
if ~isempty(fieldnames(c.evaluation)),ev=c.evaluation.result.wall_s;evcpu=c.evaluation.result.cpu_s;diagnostic=c.evaluation_controls.wall_s;dcpu=c.evaluation_controls.cpu_s;end
ob=0;obcpu=0;if ~isempty(fieldnames(c.observer)),ob=c.observer.wall_s;obcpu=c.observer.cpu_s;end
fk=0;fkcpu=0;for a=c.fk,fk=fk+a.wall_s;fkcpu=fkcpu+a.cpu_s;end
if ~isempty(fieldnames(c.fk_refinement)),fk=fk+c.fk_refinement.wall_s;fkcpu=fkcpu+c.fk_refinement.cpu_s;end
cost=struct('FE_wall_s',wall,'FE_cpu_s',cpu,'M32_wall_s',ev,'M32_cpu_s',evcpu, ...
    'evaluation_diagnostic_wall_s',diagnostic,'evaluation_diagnostic_cpu_s',dcpu, ...
    'observer_wall_s',ob,'observer_cpu_s',obcpu,'FK_wall_s',fk,'FK_cpu_s',fkcpu, ...
    'new_wall_s',wall+ev+ob+fk,'new_cpu_s',cpu+evcpu+obcpu+fkcpu, ...
    'scope','saved FE/M32/replay/FK kernels; diagnostic separate; excludes IO/tests/reporting; not peak memory');
end
function g=local_case_fk_guards(c)
g=local_fk_guards(c.fk);
if ~isempty(fieldnames(c.fk_refinement))
    f=c.fk_refinement;g.gauss_coarse=g.gauss;g.gauss=local_delta(f,c.fk(2));
    g.reflect_gauss_coarse=g.reflect_gauss;g.total_gauss_coarse=g.total_gauss;
    g.reflect_gauss=abs(f.reflect_native-c.fk(2).reflect_native)/abs(c.fk(2).reflect_native);
    g.total_gauss=abs(f.total_native-c.fk(2).total_native)/abs(c.fk(2).total_native);
end
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
            s.reference_index=j;s.phase="evaluation";
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
FEwall=0;FEcpu=0;for j=1:numel(s.runs),FEwall=FEwall+s.runs(j).result.wall_s;FEcpu=FEcpu+s.runs(j).result.cpu_s;end
Owall=0;Ocpu=0;if ~isempty(fieldnames(s.observer)),Owall=s.observer.wall_s;Ocpu=s.observer.cpu_s;end
FKwall=0;FKcpu=0;for j=1:numel(s.fk),FKwall=FKwall+s.fk(j).wall_s;FKcpu=FKcpu+s.fk(j).cpu_s;end
cost=struct('FE_wall_s',FEwall,'FE_cpu_s',FEcpu,'observer_wall_s',Owall,'observer_cpu_s',Ocpu, ...
    'FK_wall_s',FKwall,'FK_cpu_s',FKcpu,'new_wall_s',FEwall+Owall+FKwall,'six_kHz_new_solves',0, ...
    'scope','saved new FE/replay/FK kernels; excludes IO/tests/reporting; memory snapshots are not peaks');
local_json(fullfile(out,'cost_summary.json'),cost);
end


