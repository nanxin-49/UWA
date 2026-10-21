function m = two_dimensional_closeout_metrics(A, B, weights, requested)
%TWO_DIMENSIONAL_CLOSEOUT_METRICS Saved-complex-sample diagnostics only.
% A and B are in the same declared representation. No field is evaluated.
A=A(:); B=B(:); weights=weights(:); requested=logical(requested(:));
assert(isequal(size(A),size(B),size(weights),size(requested)));
assert(all(isfinite(weights)) && all(weights>=0));
n=numel(A); valid=requested & isfinite(real(A)) & isfinite(imag(A)) ...
    & isfinite(real(B)) & isfinite(imag(B));
m=struct; m.valid=valid; m.requested=requested;
m.reference_nonzero=valid & abs(B)>0;
m.phase_valid=m.reference_nonzero & abs(A)>0;
m.abs_A=nan(n,1); m.abs_B=nan(n,1); m.abs_diff=nan(n,1);
m.rel_diff=nan(n,1); m.mag_db=nan(n,1); m.phase_rad=nan(n,1);
m.mag_sq=nan(n,1); m.phase_sq=nan(n,1); m.closure=nan(n,1);
m.abs_A(valid)=abs(A(valid)); m.abs_B(valid)=abs(B(valid));
m.abs_diff(valid)=abs(A(valid)-B(valid));
ix=m.reference_nonzero;
m.rel_diff(ix)=m.abs_diff(ix)./m.abs_B(ix);
m.mag_db(ix)=20*log10(m.abs_A(ix)./m.abs_B(ix));
ix=m.phase_valid;
m.phase_rad(ix)=angle(A(ix).*conj(B(ix)));
m.mag_sq(valid)=(m.abs_A(valid)-m.abs_B(valid)).^2;
m.phase_sq(valid)=0;
m.phase_sq(ix)=4*m.abs_A(ix).*m.abs_B(ix).*sin(m.phase_rad(ix)/2).^2;
m.closure(valid)=m.abs_diff(valid).^2-m.mag_sq(valid)-m.phase_sq(valid);
w=weights(valid); aa=A(valid); bb=B(valid);
m.n_requested=nnz(requested); m.n_valid=nnz(valid);
m.n_rel=nnz(isfinite(m.rel_diff)); m.n_db=nnz(isfinite(m.mag_db));
m.n_phase=nnz(m.phase_valid); m.n_db_neg_inf=nnz(isinf(m.mag_db)&m.mag_db<0);
m.reference_energy=sum(w.*abs(bb).^2);
m.E=NaN; m.mag_fraction=NaN; m.phase_fraction=NaN;
m.phase_rms_rad=NaN; m.mag_rms_db=NaN;
if m.reference_energy>0
    m.E=sqrt(sum(w.*abs(aa-bb).^2)/m.reference_energy);
end
delta_sq=sum(w.*m.abs_diff(valid).^2);
if delta_sq>0
    m.mag_fraction=sum(w.*m.mag_sq(valid))/delta_sq;
    m.phase_fraction=sum(w.*m.phase_sq(valid))/delta_sq;
end
wp=weights(m.phase_valid);
if sum(wp)>0
    m.phase_rms_rad=sqrt(sum(wp.*m.phase_rad(m.phase_valid).^2)/sum(wp));
end
idb=isfinite(m.mag_db)&requested; wd=weights(idb);
if sum(wd)>0
    m.mag_rms_db=sqrt(sum(wd.*m.mag_db(idb).^2)/sum(wd));
end
end
