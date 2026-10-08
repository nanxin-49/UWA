function m=pm_validity_field_metrics(fk,bie,weights,mask)
%PM_VALIDITY_FIELD_METRICS Native complex, magnitude and phase decomposition.
% No phase alignment or complex scalar fitting is applied.
m.whole_grid_l2=norm(fk-bie)/norm(bie);
a=fk(mask); b=bie(mask); w=weights(mask); w=w/sum(w);
den=sum(w.*abs(b).^2); phi=angle(a.*conj(b));
m.complex_l2=sqrt(sum(w.*abs(a-b).^2)/den);
m.magnitude_l2=sqrt(sum(w.*(abs(a)-abs(b)).^2)/den);
m.phase_rms_rad=sqrt(sum(w.*phi.^2));
m.phase_mean_rad=sum(w.*phi);
m.phase_std_rad=sqrt(sum(w.*(phi-m.phase_mean_rad).^2));
m.phase_contribution_l2=sqrt(sum(w.*4.*abs(a).*abs(b).*sin(phi/2).^2)/den);
m.phase_error_energy_fraction=m.phase_contribution_l2^2/max(m.complex_l2^2,realmin);
m.decomposition_closure=abs(m.complex_l2^2-m.magnitude_l2^2-m.phase_contribution_l2^2);
end
