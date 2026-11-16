function gates=r1_bie_continuation_gates(fields,x,z,policy)
%R1_BIE_CONTINUATION_GATES Frozen numerical gates; monotonicity is diagnostic.
policy=r1_bie_continuation_io('policy',policy);
assert(iscell(fields) && ismember(numel(fields),[4 5]));
assert(isequal(size(x),[411 1]) && numel(unique(x))==411 && all(isfinite(x)) && isequal(z,3*ones(411,1)));
for j=1:numel(fields)
    assert(isequal(size(fields{j}),[411 1]) && all(isfinite(fields{j})),'r1v2:Field','Missing/nonfinite/wrong-shaped field');
end
pairing=[1 2;2 4;3 2;2 5]; count=numel(fields)-1; edges=round(linspace(1,411,5));
pairs=struct('name',{},'relative',{},'regions',{},'peak_x_m',{},'peak_abs_Pa',{}, ...
    'point_residual_real_Pa',{},'point_residual_imag_Pa',{},'point_residual_abs_Pa',{});
for j=1:count
    a=fields{pairing(j,1)}; b=fields{pairing(j,2)}; e=a-b;
    supports=[1 411;edges(1:4).' edges(2:5).']; regions=struct('name',{},'relative_L2',{},'left_x_m',{},'right_x_m',{});
    for k=1:5
        ii=supports(k,1):supports(k,2); label=sprintf('Q%d',k-1); if k==1, label='line'; end
        regions(k)=struct('name',label,'relative_L2',norm(e(ii))/max(norm(b(ii)),realmin), ...
            'left_x_m',x(ii(1)),'right_x_m',x(ii(end)));
    end
    [peak,ix]=max(abs(e));
    pairs(j)=struct('name',char(string(policy.pair_names(j))),'relative',max([regions.relative_L2]), ...
        'regions',regions,'peak_x_m',x(ix),'peak_abs_Pa',peak, ...
        'point_residual_real_Pa',real(e),'point_residual_imag_Pa',imag(e),'point_residual_abs_Pa',abs(e));
end
first=pairs(1).relative; second=pairs(2).relative; ratio=NaN; pct=NaN; reason='first is zero; ratio undefined';
if first>0, ratio=second/first; pct=100*(ratio-1); reason=''; end
diagnostic=struct('first',first,'second',second,'increase',second-first,'ratio',ratio, ...
    'percent_increase',pct,'ratio_note',reason,'grid_nonmonotonic',second>first, ...
    'mode','diagnostic','legacy_failure_identifier','r1manual:Nonmonotonic');
status='PASS'; if ~all([pairs.relative]<1e-4), status='FAIL'; end
legacy='PASS'; if second>first, legacy='FAIL'; end
gates=struct('schema_version',2,'status',status,'numeric_gates_pass',strcmp(status,'PASS'), ...
    'policy_version',policy.version,'sensitivity_gate',1e-4,'pairs',pairs, ...
    'grid_diagnostic',diagnostic,'legacy_policy_status',legacy,'x_m',x,'z_m',z, ...
    'final_reference_stage',3,'C4_window_check_complete',numel(fields)==5);
end
