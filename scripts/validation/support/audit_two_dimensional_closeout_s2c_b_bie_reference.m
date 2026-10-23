function out = audit_two_dimensional_closeout_s2c_b_bie_reference(run_dir)
%AUDIT_TWO_DIMENSIONAL_CLOSEOUT_S2C_B_BIE_REFERENCE Saved five-level audit.
% Reads completed K_HIGH solves; performs no new field or density evaluation.
arguments
    run_dir (1,:) char
end
dest=fullfile(run_dir,'P1'); out_file=fullfile(dest,'BIE_reference_data.mat');
assert(~isfile(out_file),'Refusing to overwrite a reference audit');
ids={'main96','self64','space14','window60','window70'};
U=complex(zeros(417,5)); receipts=cell(5,1); wall=zeros(1,5);
for j=1:5
    f=fullfile(run_dir,'jobs',['bie_' ids{j} '.mat']); assert(isfile(f));
    q=load(f,'boundary','receipt'); b=q.boundary; v=q.receipt;
    assert(strcmp(v.status,'COMPLETE') && strcmp(v.job_id,ids{j}));
    assert(b.finite && numel(b.receiver_field)==417);
    assert(all(isfinite(real(b.receiver_field))) && all(isfinite(imag(b.receiver_field))));
    assert(b.linear_residual<1e-10 && b.offgrid_boundary_residual<1e-6);
    if j==1
        x=b.receiver_x_m(:); z=b.receiver_z_m(:);
        source_sha=v.source_sha256; cfg_sha=v.configuration_sha256;
        solver_sha=v.solver_sha256;
    else
        assert(isequal(x,b.receiver_x_m(:)) && isequal(z,b.receiver_z_m(:)));
        assert(strcmp(source_sha,v.source_sha256) && strcmp(cfg_sha,v.configuration_sha256));
        assert(strcmp(solver_sha,v.solver_sha256));
        assert(v.authorized_subbudget_s==2400);
    end
    U(:,j)=b.receiver_field(:); receipts{j}=v; wall(j)=v.wall_seconds;
end
assert(all(z(1:411)==3) && nnz(z~=3)==6);
root=fileparts(fileparts(fileparts(fileparts(run_dir))));
f=load(fullfile(root,'results/comparisons/two_dimensional_channel_closeout', ...
    '20261023_s2b_1305/formal_reference_v2/reference_data.mat'),'ix9','x','z');
assert(isequal(x,f.x(:)) && isequal(z,f.z(:)));
g=closeout_s2c_b_effective_gate(run_dir);
total=sum(wall)+g.other_physics_reserve_s+g.failure_contingency_s;
assert(total<g.subbudget_s && g.core_remaining_s-total>=g.S3A_reserve_s);
self_abs=abs(U(:,1)-U(:,2)); space_abs=abs(U(:,1)-U(:,3));
window_abs=abs(U(:,4)-U(:,5)); resolution_abs=self_abs+space_abs+window_abs;
point=table((1:417)',x,z,self_abs,space_abs,window_abs,resolution_abs, ...
    5*resolution_abs,abs(U(:,1)), ...
    'VariableNames',{'sample_index','x_m','z_m','self_abs','space_abs', ...
    'window_abs','empirical_resolution_abs','five_times_resolution_abs','BIE_main_abs'});
version=table;
for j=1:5
    version=[version;table(repmat(string(ids{j}),417,1),(1:417)',x,z, ...
        real(U(:,j)),imag(U(:,j)),repmat(wall(j),417,1), ...
        'VariableNames',{'level','sample_index','x_m','z_m','u_real','u_imag','job_wall_seconds'})]; %#ok<AGROW>
end
center=find(abs(x(1:411))<=1); edge=find(abs(x(1:411))>=30);
five=f.ix9([2 4 5 6 8]); supports={(1:411)',center,f.ix9(:),five(:),edge};
names=["line_411","center_11","nine_9","cross_5","edge_104"];
summary=table;
for j=1:5
    ii=supports{j}; n=numel(ii); q=ones(n,1)/n;
    den=sum(q.*abs(U(ii,1)).^2);
    e=sqrt(sum(q.*resolution_abs(ii).^2)/den);
    summary=[summary;table(names(j),n,den,e,5*e, ...
        sqrt(sum(q.*self_abs(ii).^2)/den),sqrt(sum(q.*space_abs(ii).^2)/den), ...
        sqrt(sum(q.*window_abs(ii).^2)/den),max(resolution_abs(ii)), ...
        'VariableNames',{'support_id','n','BIE_main_denominator', ...
        'empirical_resolution_L2','five_times_resolution_L2', ...
        'self_component_L2','space_component_L2','window_component_L2', ...
        'max_point_resolution_abs'})]; %#ok<AGROW>
end
assert(isequal(summary.n,[411;11;9;5;104]));
if ~isfolder(dest), mkdir(dest); end
writetable(version,fullfile(dest,'BIE_reference_versions.csv'));
writetable(point,fullfile(dest,'BIE_resolution_point.csv'));
writetable(summary,fullfile(dest,'BIE_resolution_support.csv'));
save(out_file,'U','x','z','resolution_abs','self_abs','space_abs','window_abs', ...
    'summary','receipts','wall','-v7.3');
out=struct('status','FIVE_BIE_LEVELS_COMPLETE_EMPIRICAL_REFERENCE', ...
    'new_acoustic_calculation_seconds',0,'actual_BIE_seconds',sum(wall), ...
    'actual_BIE_job_seconds',wall,'K_HIGH_subbudget_s',g.subbudget_s, ...
    'conservative_cost_with_other_and_failure_reserve_s',total, ...
    'S3A_balance_after_conservative_plan_s',g.core_remaining_s-total, ...
    'source_sha256',source_sha,'configuration_sha256',cfg_sha, ...
    'solver_sha256',solver_sha,'BIE_point_count',417, ...
    'line_empirical_resolution_L2',summary.empirical_resolution_L2(1), ...
    'nine_empirical_resolution_L2',summary.empirical_resolution_L2(3), ...
    'empirical_not_strict_bound',true,'model_resolution_not_yet_classified',true, ...
    'budget_amendment_sha256',g.budget_amendment_sha256);
closeout_s2b_json(fullfile(run_dir,'audit','P1_five_BIE_reference_audit.json'),out);
disp(out);
end
