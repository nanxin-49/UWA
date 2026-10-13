function report=audit_deepsea3d_g1_frequency_artifacts()
%AUDIT_DEEPSEA3D_G1_FREQUENCY_ARTIFACTS Native three-point evidence and trends.
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
out=fullfile(root,'results','validation','deepsea3d_g1_frequency');
a=load(fullfile(out,'frequency.mat'),'state');s=a.state;
assert(all([s.cases.complete]) && ~s.WP2_allowed && s.global_G1=="OPEN");
cases={s.cases(1),s.baseline6,s.cases(2)};rows=struct([]);envelopes=struct([]);
for j=1:3
    c=cases{j};f=c.cfg.frequency_hz;restored=c.cfg;restored.frequency_hz=6000;assert(isequal(restored,s.baseline6.cfg));
    assert(c.flat_control.passed);bg=deepsea3d_flat(c.cfg,0);
    if c.reference_ready && c.G1_status~="INCONCLUSIVE"
        r=c.runs(c.reference_index).result;fk=c.fk(2);m=local_metrics(fk,r);
        if j~=2
            folder=fullfile(out,sprintf('f%d',f));guard=readtable(fullfile(folder,'reference_guards.csv'));
            assert(all(guard{1,:}<1e-4) && c.observer.term_match_error<1e-12 && c.observer.norm_match_error<1e-12);
            assert(r.forcing_series_relative<1e-8 && all(diff(r.boundary_coefficient_norms(end-4:end))<0));
            assert(all(diff(abs(r.terms(end-4:end)))<0));
            ev=c.evaluation_controls;assert(ev.M_low==32 && ev.M_high==40 && ev.term_match<1e-12 && ev.norm_match<1e-12);
            assert(ev.tail32_relative<1e-4 && ev.tail40_relative<1e-4 && ev.boundary_change<1e-4 && ev.halfgrid_change<1e-4);
            raw=readtable(fullfile(folder,'physical_metrics.csv'));
            fields=fieldnames(m);
            for n=1:numel(fields),assert(abs(raw.(fields{n})-m.(fields{n}))<1e-13);end
            refcsv=readtable(fullfile(folder,'reference_runs.csv'));fkcsv=readtable(fullfile(folder,'fk_runs.csv'));
            R=complex(refcsv.reflect_real(c.reference_index),refcsv.reflect_imag(c.reference_index));
            F=complex(fkcsv.reflect_real(2),fkcsv.reflect_imag(2));
            D=complex(fkcsv.delta_real(2),fkcsv.delta_imag(2));I=complex(refcsv.delta_real(c.reference_index),refcsv.delta_imag(c.reference_index));
            ratio=[F/R,D/I];assert(max(abs(abs(ratio-1)-[m.reflected_complex,m.increment_complex]))<1e-13);
            assert(max(abs(angle(ratio)-[m.reflected_phase_rad,m.increment_phase_rad]))<1e-13);
            assert(max(abs(abs(ratio-1).^2-((abs(ratio)-1).^2+4*abs(ratio).*sin(angle(ratio)/2).^2)))<1e-16);
        end
        usable_fk=c.fk;
        if isfield(c,'fk_refinement') && ~isempty(fieldnames(c.fk_refinement))
            usable_fk=[c.fk([1 2 4 5]),c.fk_refinement];
            assert(c.fk_guards.gauss_coarse>=1e-3 && c.fk_guards.gauss<1e-3);
        end
        if j~=2
            g=c.fk_guards;assert(classify_deepsea3d_amplitude(true,max([g.grid,g.aperture_12_16,g.aperture_16_24]),g.gauss,m)==c.G1_status);
        else
            assert(abs(m.reflected_complex-c.comparison.metrics_y.reflected_complex)<1e-14);
            assert(abs(m.increment_complex-c.comparison.metrics_y.increment_complex)<1e-14);
        end
        variants=local_variants(c,bg);row=local_row(c,m,r,fk,bg,j==2,variants);
        names=fieldnames(m);
        for n=1:numel(names)
            name=names{n};values=local_errors(usable_fk,variants,name);
            e=struct('frequency_hz',f,'metric',string(name),'minimum',min(values),'maximum',max(values), ...
                'reference_variants',numel(variants),'FK_variants',numel(usable_fk),'combinations',numel(values));
            if isempty(envelopes),envelopes=e;else,envelopes(end+1)=e;end %#ok<AGROW>
        end
    else
        assert(c.G1_status=="INCONCLUSIVE");
        if ~c.reference_ready,assert(isempty(c.fk));end
        row=local_missing(c,j==2);
    end
    if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
end
trends=struct([]);
for name=["reflected_complex","increment_complex","reflected_magnitude","increment_magnitude","reflected_phase_rad","increment_phase_rad"]
    for pair=[1 2]
        lo=rows(pair).frequency_hz;hi=rows(pair+1).frequency_hz;
        e0=envelopes([envelopes.frequency_hz]==lo & string({envelopes.metric})==name);
        e1=envelopes([envelopes.frequency_hz]==hi & string({envelopes.metric})==name);
        direction="INCONCLUSIVE";
        minimum0=NaN;maximum0=NaN;minimum1=NaN;maximum1=NaN;basis="native_relative_error";
        if ~isempty(e0) && ~isempty(e1)
            minimum0=e0.minimum;maximum0=e0.maximum;minimum1=e1.minimum;maximum1=e1.maximum;
            if endsWith(name,"phase_rad")
                basis="absolute_phase_rad";[minimum0,maximum0]=local_abs_range(e0);[minimum1,maximum1]=local_abs_range(e1);
            end
            direction="UNRESOLVED_NUMERICAL_RANGES_OVERLAP";
            if maximum1<minimum0,direction="DISCRETE_DECREASE";elseif minimum1>maximum0,direction="DISCRETE_INCREASE";end
        end
        row=struct('metric',name,'comparison_basis',basis,'frequency_low_hz',lo,'frequency_high_hz',hi, ...
            'minimum_low',minimum0,'maximum_low',maximum0,'minimum_high',minimum1,'maximum_high',maximum1,'decision',direction);
        if isempty(trends),trends=row;else,trends(end+1)=row;end %#ok<AGROW>
    end
end
costs=struct([]);
for c=s.cases
    folder=fullfile(out,sprintf('f%d',c.frequency_hz));k=jsondecode(fileread(fullfile(folder,'cost_summary.json')));
    row=struct('frequency_hz',c.frequency_hz,'new_wall_s',k.new_wall_s,'new_cpu_s',k.new_cpu_s, ...
        'evaluation_diagnostic_wall_s',k.evaluation_diagnostic_wall_s,'FE_solves',numel(c.runs)+double(~isempty(fieldnames(c.evaluation))), ...
        'FK_solves',numel(c.fk)+double(~isempty(fieldnames(c.fk_refinement))));
    if isempty(costs),costs=row;else,costs(end+1)=row;end %#ok<AGROW>
end
report=struct('complete',true,'frequency_only_change',true,'six_kHz_new_solves',0, ...
    'summary',rows,'trends',trends,'global_G1',"OPEN",'WP2_allowed',false, ...
    'new_kernel_wall_s',sum([costs.new_wall_s]),'new_kernel_cpu_s',sum([costs.new_cpu_s]), ...
    'evaluation_diagnostic_wall_s',sum([costs.evaluation_diagnostic_wall_s]), ...
    'scope','three discrete frequencies; empirical ranges including order pressure sensitivity, not a confidence interval or rigorous reference error bound');
writetable(struct2table(rows),fullfile(out,'frequency_summary.csv'));
writetable(struct2table(envelopes),fullfile(out,'numerical_error_envelopes.csv'));
writetable(struct2table(trends),fullfile(out,'frequency_trends.csv'));
writetable(struct2table(costs),fullfile(out,'new_costs.csv'));
local_json(fullfile(out,'saved_data_checks.json'),report);
disp(struct2table(rows));disp(struct2table(trends));
end
function m=local_metrics(f,r)
m=struct('reflected_complex',abs(f.reflect_native-r.reflect_native)/abs(r.reflect_native), ...
    'increment_complex',abs(f.rough_correction-r.rough_correction)/abs(r.rough_correction), ...
    'reflected_magnitude',abs(abs(f.reflect_native)-abs(r.reflect_native))/abs(r.reflect_native), ...
    'increment_magnitude',abs(abs(f.rough_correction)-abs(r.rough_correction))/abs(r.rough_correction), ...
    'reflected_phase_rad',angle(f.reflect_native*conj(r.reflect_native)), ...
    'increment_phase_rad',angle(f.rough_correction*conj(r.rough_correction)));
end
function v=local_variants(c,bg)
v=struct([]);
for index=[c.coarse_index,c.fine_index,c.reference_index]
    r=c.runs(index).result;
    for N=[r.order-4,r.order-2,r.order]
        d=r.partials(N);row=struct('reflect_native',bg.reflect_native+d,'rough_correction',d,'total_native',bg.total_native+d);
        if isempty(v),v=row;else,v(end+1)=row;end %#ok<AGROW>
    end
end
end
function values=local_errors(fk,refs,name)
values=zeros(numel(fk)*numel(refs),1);i=0;
for f=fk
    for r=refs
        m=local_metrics(f,r);i=i+1;values(i)=m.(name);
    end
end
end
function row=local_row(c,m,r,fk,bg,reused,variants)
v=c.runs(c.reference_index).settings;
row=struct('frequency_hz',c.cfg.frequency_hz,'status',c.G1_status,'reused',reused,'N',r.order,'L',v(1),'h',v(2), ...
    'reflected_complex',m.reflected_complex,'increment_complex',m.increment_complex, ...
    'reflected_magnitude',m.reflected_magnitude,'increment_magnitude',m.increment_magnitude, ...
    'reflected_phase_rad',m.reflected_phase_rad,'increment_phase_rad',m.increment_phase_rad, ...
    'reference_reflect_real',real(r.reflect_native),'reference_reflect_imag',imag(r.reflect_native), ...
    'reference_delta_real',real(r.rough_correction),'reference_delta_imag',imag(r.rough_correction), ...
    'reference_reflect_abs',abs(r.reflect_native),'reference_delta_abs',abs(r.rough_correction),'reference_total_abs',abs(r.total_native), ...
    'flat_total_condition',(abs(bg.direct_native)+abs(bg.reflect_native))/abs(bg.total_native), ...
    'rough_total_condition',(abs(bg.direct_native)+abs(r.reflect_native))/abs(r.total_native), ...
    'reflected_background_increment_condition',(abs(bg.reflect_native)+abs(r.rough_correction))/abs(r.reflect_native), ...
    'pressure_order_sum_indicator',sum(abs(r.terms))/abs(r.rough_correction),'FK_density_indicator',fk.cancellation_condition, ...
    'total_complex_error',abs(fk.total_native-r.total_native)/abs(r.total_native), ...
    'reference_total_variant_change',max(arrayfun(@(a)abs(a.total_native-r.total_native)/abs(r.total_native),variants)));
end
function row=local_missing(c,reused)
r=struct('order',0,'reflect_native',NaN,'rough_correction',NaN,'total_native',NaN,'terms',NaN);
m=local_metrics(r,r);fk=r;fk.cancellation_condition=NaN;bg=deepsea3d_flat(c.cfg,0);
dummy=c;dummy.runs=struct('settings',[NaN NaN 0]);dummy.reference_index=1;row=local_row(dummy,m,r,fk,bg,reused,r);
end
function local_json(file,s)
fid=fopen(file,'w');assert(fid>=0);clean=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(s,PrettyPrint=true));
end
function [lo,hi]=local_abs_range(e)
hi=max(abs([e.minimum,e.maximum]));
if e.minimum<=0 && e.maximum>=0,lo=0;else,lo=min(abs([e.minimum,e.maximum]));end
end
