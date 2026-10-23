function m = closeout_s2b_metrics(a,b,q)
%CLOSEOUT_S2B_METRICS Fixed-support statistics, retaining zero/missing cases.
a=a(:); b=b(:); q=q(:); assert(isequal(size(a),size(b),size(q)));
valid=isfinite(a)&isfinite(b)&isfinite(q)&q>=0;
m=struct('n_total',numel(a),'n_valid',nnz(valid),'n_missing',nnz(~valid));
a=a(valid); b=b(valid); q=q(valid); q=q/sum(q);
aa=abs(a); bb=abs(b); d=abs(a-b); den=sum(q.*bb.^2);
phase=angle(a.*conj(b)); nz=aa>0&bb>0; phase(~nz)=NaN;
db=20*log10(aa./bb); db(bb==0)=NaN;
rel=d./bb; rel(bb==0)=NaN;
phase_term=zeros(size(a)); phase_term(nz)=4*aa(nz).*bb(nz).*sin(phase(nz)/2).^2;
m.denominator=den; m.n_zero_candidate=nnz(aa==0); m.n_zero_reference=nnz(bb==0);
m.near_null_threshold=1e-12*max([bb;0]);
m.n_near_null=nnz(bb<=m.near_null_threshold); m.near_null_kind='scale_only_not_reference_error';
m.n_phase_valid=nnz(nz); m.n_db_finite=nnz(isfinite(db)); m.n_db_nonfinite=nnz(~isfinite(db));
m.complex_l2=sqrt(sum(q.*d.^2)/den); m.magnitude_l2=sqrt(sum(q.*(aa-bb).^2)/den);
m.magnitude_squared_term=sum(q.*(aa-bb).^2)/den;
m.phase_squared_term=sum(q.*phase_term)/den;
m.decomposition_residual=m.complex_l2^2-m.magnitude_squared_term-m.phase_squared_term;
m.phase_rms_rad=sqrt(sum(q(nz).*phase(nz).^2)/sum(q(nz)));
fin=isfinite(db); m.db_rms=sqrt(sum(q(fin).*db(fin).^2)/sum(q(fin)));
m.db_mean=sum(q(fin).*db(fin))/sum(q(fin)); m.db_max_abs=max(abs(db(fin)),[],'omitnan');
m.phase_max_abs=max(abs(phase),[],'omitnan');
m.relative_mean=mean(rel,'omitnan'); m.relative_median=median(rel,'omitnan');
[m.relative_max,m.worst_index]=max(rel,[],'omitnan');
m.absolute_mean=mean(d); m.absolute_median=median(d); m.absolute_max=max(d);
if den==0
    m.complex_l2=NaN; m.magnitude_l2=NaN; m.magnitude_squared_term=NaN;
    m.phase_squared_term=NaN; m.decomposition_residual=NaN;
end
% Empty extrema are represented explicitly in tabular exports.
fn=fieldnames(m); for j=1:numel(fn), if isempty(m.(fn{j})), m.(fn{j})=NaN; end, end
end
