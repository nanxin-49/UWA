function receipt=audit_r1_bie_manual_final(root)
%AUDIT_R1_BIE_MANUAL_FINAL Recompute reference/model resolution from saved data.
q=load(fullfile(root,'common.mat')); a=load(fullfile(root,'comparison.mat'));
ids={'bie_coarse','bie_main','bie_self','bie_space','bie_window'}; B=cell(5,1);
for j=1:5, b=load(fullfile(root,'jobs',ids{j},'data.mat'),'field'); B{j}=b.field; end
values=[local_worst(B{1},B{2}),local_worst(B{2},B{4}),local_worst(B{3},B{2}),local_worst(B{2},B{5})];
assert(all(values<1e-4) && values(2)<=values(1),'Reference hierarchy gate failed');
assert(isequal(a.U(:,2),B{4}) && a.reference_qualified && numel(q.x)==411);
bunc=abs(B{1}-B{2})+abs(B{2}-B{4})+abs(B{3}-B{2})+abs(B{2}-B{5});
punc=abs(q.PS{3}.field-q.PS{4}.field)+abs(q.PS{4}.field-q.PS{5}.field)+abs(q.PS{5}.field-q.PS{6}.field);
func=abs(q.FK{1}.receiver_field-q.FK{2}.receiver_field)+abs(q.FK{2}.receiver_field-q.FK{3}.receiver_field)+ ...
    abs(q.FK{3}.receiver_field-q.FK{4}.receiver_field);
assert(max(abs(bunc-a.bunc))<1e-11 && max(abs(punc-a.punc))<1e-11 && max(abs(func-a.func))<1e-11);
assert(max(abs(values-a.bie_values))<1e-12);
M=readtable(fullfile(root,'model_metrics.csv')); resolution=table;
for m=[1 3]
    if m==1, field=q.PS{6}.field; uncertainty=punc; name="PS";
    else, field=q.FK{2}.receiver_field; uncertainty=func; name="FullFK"; end
    E=norm(field-B{4})/norm(B{4}); limit=5*norm(bunc+uncertainty)/norm(B{4});
    row=find(string(M.method)==name); assert(isscalar(row));
    assert(abs(M.complex_relative_L2(row)-E)<1e-12 && ...
        abs(M.five_times_empirical_resolution(row)-limit)<1e-12 && logical(M.resolved(row))==(E>limit));
    resolution=[resolution;table(name,E,limit,E>limit,'VariableNames', ...
        {'method','complex_relative_L2','five_times_empirical_resolution','resolved'})]; %#ok<AGROW>
end
scale=norm(B{4})/sqrt(411); mask=abs(B{4})>=.01*scale;
assert(abs(scale-a.scale)<1e-12 && isequal(mask,a.mask));
writetable(resolution,fullfile(root,'independent_reference_resolution.csv'));
receipt=struct('status','PASS','failed_checks',0,'reference_sensitivities',values, ...
    'phase_valid_points',nnz(mask),'reference_scope','empirical convergence for fixed 10cm 4kHz seed260001 411-point line');
r1_bie_manual_io('write',fullfile(root,'independent_final_metrics_audit.json'),receipt);
end
function v=local_worst(a,b)
v=norm(a-b)/max(norm(b),realmin); edges=round(linspace(1,numel(a),5));
for j=1:4, ii=edges(j):edges(j+1); v=max(v,norm(a(ii)-b(ii))/max(norm(b(ii)),realmin)); end
end
