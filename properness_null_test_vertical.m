function result=properness_null_test_vertical(C_train,H_test,n_mc,seed)
%PROPERNESS_NULL_TEST_VERTICAL Finite-sample null test for sample ||P||/||C||.
arguments
    C_train {mustBeNumeric}
    H_test {mustBeNumeric}
    n_mc (1,1) double {mustBeInteger,mustBePositive} = 2000
    seed (1,1) double {mustBeFinite} = 900001
end
[F,L]=size(H_test);
if ~isequal(size(C_train),[F,F]) || L<2, error('C_train must be F-by-F and H_test F-by-L, L>=2.'); end
C_train=0.5*(double(C_train)+double(C_train)');
Xt=double(H_test)-mean(double(H_test),2);
Ct=(Xt*Xt')/(L-1); Pt=(Xt*Xt.')/(L-1);
observed=norm(Pt,'fro')/max(norm(Ct,'fro'),eps);
[Q,d]=eig(C_train,'vector'); d=real(d); keep=d>max(max(d),eps)*1e-12;
A=Q(:,keep).*sqrt(max(d(keep).',0));
rng(mod(round(seed),2^32),'twister'); null_ratio=zeros(n_mc,1);
for bb=1:n_mc
    Z=(randn(sum(keep),L)+1i*randn(sum(keep),L))/sqrt(2);
    X=A*Z; X=X-mean(X,2); C=(X*X')/(L-1); P=(X*X.')/(L-1);
    null_ratio(bb)=norm(P,'fro')/max(norm(C,'fro'),eps);
end
q=sort(null_ratio);
central90=local_q(q,[0.05,0.95]); central95=local_q(q,[0.025,0.975]);
central99=local_q(q,[0.005,0.995]); upper=local_q(q,[0.90,0.95,0.99]);
p_upper=(1+sum(q>=observed))/(n_mc+1);
inside95=observed>=central95(1)&&observed<=central95(2);
result=struct('observed_ratio',observed,'null_mean',mean(q),'null_median',median(q), ...
    'null_interval_90',central90,'null_interval_95',central95, ...
    'null_interval_99',central99,'one_sided_upper_90_95_99',upper, ...
    'upper_tail_p_value',p_upper,'observed_inside_central_95',inside95, ...
    'reject_proper_at_5pct',observed>central95(2), ...
    'reject_proper_at_5pct_one_sided',observed>upper(2), ...
    'reject_proper_at_1pct_one_sided',observed>upper(3), ...
    'decision',local_central_decision(observed,central95,central99), ...
    'one_sided_decision',local_one_sided_decision(observed,upper), ...
    'decision_rule','prescribed central interval rule; one-sided result retained separately', ...
    'n_mc',n_mc,'seed',seed, ...
    'test_sample_count',L,'frequency_count',F,'null_ratio_samples',q);
end

function v=local_q(q,p)
n=numel(q); idx=max(1,min(n,round(1+(n-1).*p))); v=q(idx).';
end
function value=local_one_sided_decision(observed,upper)
if observed>upper(3), value='proper hypothesis rejected above one-sided 99% upper bound';
elseif observed>upper(2), value='proper hypothesis rejected at one-sided 5% level';
else, value='proper hypothesis not rejected'; end
end
function value=local_central_decision(observed,ci95,ci99)
if observed>ci99(2), value='proper hypothesis rejected above central 99% upper bound';
elseif observed>ci95(2), value='proper hypothesis rejected above central 95% upper bound';
else, value='proper hypothesis not rejected'; end
end
