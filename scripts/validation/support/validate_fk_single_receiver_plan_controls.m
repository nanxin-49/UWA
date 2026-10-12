function rows = validate_fk_single_receiver_plan_controls(cfg,plan_in,plan_composite)
%VALIDATE_FK_SINGLE_RECEIVER_PLAN_CONTROLS Actual dependencies and typed rejects.
rows=struct([]);
for method={'fastin','composite'}
    action=method{1}; p=plan_in;
    if strcmp(action,'composite'), p=plan_composite; end
    for variant={'height','source','source_kx','normalization','period','origin','frequency','sound_speed','dispersion','grid','Tx','bound','nodes','rank','Rx_x','Rx_z','metadata','surface'}
        name=variant{1}; c=cfg; expected='fk_single_receiver:PlanMismatch';
        switch name
            case 'height'
                c.surface_z_m(1)=.501; expected='fk_single_receiver:HeightDomain';
            case 'source'
                [~,index]=max(abs(c.spectral_source.coeff));
                c.spectral_source.coeff(index)=c.spectral_source.coeff(index)*.99;
            case 'source_kx'
                [~,index]=max(abs(c.spectral_source.coeff));
                c.spectral_source.kx(index)=c.spectral_source.kx(index)+2*pi/c.spectral_source.period_m;
                c.spectral_source.kz=sqrt((2*pi*c.frequency_hz/c.c0_mps)^2-c.spectral_source.kx.^2);
            case 'normalization'
                c.spectral_source.n_original=c.spectral_source.n_original+1;
            case 'period'
                c.spectral_source.period_m=c.spectral_source.period_m+mean(diff(c.surface_x_m));
            case 'origin'
                c.spectral_source.x0_m=c.spectral_source.x0_m+mean(diff(c.surface_x_m));
            case 'frequency'
                c.frequency_hz=c.frequency_hz+100;
                c.spectral_source.kz=sqrt((2*pi*c.frequency_hz/c.c0_mps)^2-c.spectral_source.kx.^2);
            case 'dispersion'
                c.spectral_source.kz(1)=c.spectral_source.kz(1)+.01;
                expected='fk_single_receiver:InvalidSpectrum';
            case 'sound_speed'
                c.c0_mps=c.c0_mps+1;
                c.spectral_source.kz=sqrt((2*pi*c.frequency_hz/c.c0_mps)^2-c.spectral_source.kx.^2);
            case 'grid'
                c.surface_x_m=c.surface_x_m+mean(diff(c.surface_x_m));
            case 'Tx'
                c.spectral_source.z_tx_m=c.spectral_source.z_tx_m+1;
            case 'bound'
                c.height_bound_m=.6;
            case 'nodes'
                c.height_nodes=65;
            case 'rank'
                c.rank=4;
                if strcmp(action,'composite'), expected=''; end
            case 'Rx_x'
                c.receiver_x_m=1;
                if strcmp(action,'fastin'), expected=''; end
            case 'Rx_z'
                c.receiver_z_m=3.25;
                if strcmp(action,'fastin'), expected=''; end
            case 'metadata'
                c.spectral_source.annotation='inert source metadata'; expected='';
            case 'surface'
                x=c.surface_x_m; c.surface_z_m=.02*cos(.1*x); slope=-.002*sin(.1*x);
                c.surface_normal=[-slope ones(size(slope))]./sqrt(1+slope.^2); expected='';
        end
        value=0; actual=''; tolerance=0;
        if ~isempty(expected)
            e=fk_single_receiver_expect_error(@()fk_single_receiver_prototype(action,c,p),expected);
            actual=e.actual_id;
        else
            reused=fk_single_receiver_prototype(action,c,p);
            fresh=fk_single_receiver_prototype(['prepare_' action],c);
            rebuilt=fk_single_receiver_prototype(action,c,fresh);
            value=abs(reused-rebuilt)/max(abs(rebuilt),1e-14); tolerance=1e-12;
            assert(isfinite(value) && value<=tolerance,'fk_validation:ReuseMismatch','Reuse differs from rebuilding.');
        end
        rows=local_add(rows,cfg.frequency_hz,action,name,expected,actual,value,tolerance);
    end
    other=plan_composite; if strcmp(action,'composite'), other=plan_in; end
    for variant={'kind','missing'}
        name=variant{1}; wrong=other; expected='fk_single_receiver:PlanKind';
        if strcmp(name,'missing'), wrong=[]; expected='fk_single_receiver:MissingPlan'; end
        e=fk_single_receiver_expect_error(@()fk_single_receiver_prototype(action,cfg,wrong),expected);
        rows=local_add(rows,cfg.frequency_hz,action,name,expected,e.actual_id,0,0);
    end
end
% Prove the matcher itself fails for wrong identifiers and absent errors.
fk_single_receiver_expect_error(@()fk_single_receiver_expect_error( ...
    @()error('fk_validation:Unrelated','Injected unrelated exception.'),'fk_single_receiver:PlanMismatch'), ...
    'fk_validation:UnexpectedError');
fk_single_receiver_expect_error(@()fk_single_receiver_expect_error(@()1,'fk_single_receiver:PlanMismatch'), ...
    'fk_validation:UnexpectedError');
for name={'matcher_wrong_identifier','matcher_no_exception'}
    rows=local_add(rows,cfg.frequency_hz,'verifier',name{1}, ...
        'fk_validation:UnexpectedError','fk_validation:UnexpectedError',0,0);
end
for input={NaN,Inf,complex(1,Inf)}
    value=input{1};
    fk_single_receiver_expect_error(@()fk_single_receiver_error_metrics(1,value), ...
        'fk_validation:NonfinitePressure');
    rows=local_add(rows,cfg.frequency_hz,'verifier',sprintf('nonfinite_pressure_%d',numel(rows)), ...
        'fk_validation:NonfinitePressure','fk_validation:NonfinitePressure',0,0);
end
z=fk_single_receiver_error_metrics(0,0,1);
assert(z.finite_metrics && z.near_zero && ~z.phase_defined && ~z.relative_defined && isnan(z.phase_error_rad));
t=fk_single_receiver_error_metrics(1e-16,2e-16,1);
assert(t.finite_metrics && t.near_zero && t.raw_relative_error==1 && t.relative_error<1);
rows=local_add(rows,cfg.frequency_hz,'verifier','exact_zero_undefined_phase','', '',0,0);
rows=local_add(rows,cfg.frequency_hz,'verifier','near_zero_floor_and_raw','', '',0,0);
end
function rows=local_add(rows,f,method,name,expected,actual,value,tolerance)
row=struct('frequency_hz',f,'method',method,'case_name',name,'expected_id',expected, ...
    'actual_id',actual,'value',value,'tolerance',tolerance,'pass',true);
if isempty(rows), rows=row; else, rows(end+1)=row; end
end
