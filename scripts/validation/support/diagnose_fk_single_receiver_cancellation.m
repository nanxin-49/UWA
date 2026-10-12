function [rows,scan] = diagnose_fk_single_receiver_cancellation(cfg,plan_in)
%DIAGNOSE_FK_SINGLE_RECEIVER_CANCELLATION Fixed existing strong-sine fixture.
% Receiver search is a diagnostic, never a source/phase/amplitude fit.
x=cfg.surface_x_m; c=cfg; c.surface_z_m=.45*cos(.47*x); slope=-.45*.47*sin(.47*x);
c.surface_normal=[-slope ones(size(slope))]./sqrt(1+slope.^2);
[u,ux,uz]=fk_single_receiver_incident(c.spectral_source,x,c.surface_z_m);
fc=struct('frequency_hz',c.frequency_hz,'c0_mps',c.c0_mps,'surface_x_m',x, ...
    'surface_z_m',c.surface_z_m,'surface_normal',c.surface_normal, ...
    'surface_incident',u,'surface_incident_dn',c.surface_normal(:,1).*ux+c.surface_normal(:,2).*uz, ...
    'receiver_x_m',linspace(-60,60,481)','receiver_z_m',3*ones(481,1));
refline=solve_full_kirchhoff_2d(fc); mag=abs(refline.receiver_field); rx=fc.receiver_x_m;
scale=max(mag); ix=find(mag(2:end-1)<mag(1:end-2) & mag(2:end-1)<mag(3:end))+1;
ix=ix(abs(rx(ix))<=50); assert(~isempty(ix),'fk_validation:NoCancellationMinimum','No interior minimum found.');
[~,best]=min(mag(ix)); index=ix(best);
options=optimset('TolX',1e-8,'Display','off');
location=fminbnd(@(r)local_magnitude(fc,r),rx(index-1),rx(index+1),options);
c.receiver_x_m=location; fc.receiver_x_m=location; fc.receiver_z_m=c.receiver_z_m;
ref=fk_single_receiver_prototype('full',c);
check=solve_full_kirchhoff_2d(fc); assert(abs(ref-check.receiver_field)/abs(ref)<1e-12);
contribution=local_density(fc,refline.quadrature_weights_m);
assert(abs(-sum(contribution)-ref)/max(abs(ref),1e-14)<1e-10);
cancellation=abs(ref)/sum(abs(contribution));
fast=fk_single_receiver_prototype('fastin',c,plan_in);
pc=fk_single_receiver_prototype('prepare_composite',c);
composite=fk_single_receiver_prototype('composite',c,pc);
rank5=c; rank5.rank=5; p5=fk_single_receiver_prototype('prepare_fastin',rank5);
fast5=fk_single_receiver_prototype('fastin',rank5,p5);
nodes65=c; nodes65.height_nodes=65; p65=fk_single_receiver_prototype('prepare_composite',nodes65);
composite65=fk_single_receiver_prototype('composite',nodes65,p65);
fine=c; fine.surface_x_m=linspace(-60,60,8193)';
fine.surface_z_m=.45*cos(.47*fine.surface_x_m); sl=-.45*.47*sin(.47*fine.surface_x_m);
fine.surface_normal=[-sl ones(size(sl))]./sqrt(1+sl.^2);
fine_ref=fk_single_receiver_prototype('full',fine);
rows=struct([]); methods={'fastin','composite'}; values=[fast composite];
for j=1:2
    m=fk_single_receiver_error_metrics(ref,values(j),scale);
    row=struct('frequency_hz',c.frequency_hz,'method',methods{j},'receiver_x_m',location, ...
        'full_real',real(ref),'full_imag',imag(ref),'approx_real',real(values(j)),'approx_imag',imag(values(j)), ...
        'cancellation_ratio',cancellation,'peak_relative_magnitude',abs(ref)/scale, ...
        'original_relative_gate_pass',m.relative_error<1e-3, ...
        'rank5_relative_error',abs(fast5-ref)/abs(ref), ...
        'nodes65_relative_error',abs(composite65-ref)/abs(ref), ...
        'quadrature4097_8193_relative',abs(fine_ref-ref)/abs(fine_ref));
    names=fieldnames(m); for t=1:numel(names), row.(names{t})=m.(names{t}); end
    if isempty(rows), rows=row; else, rows(end+1)=row; end %#ok<AGROW>
end
scan=table(c.frequency_hz+zeros(size(rx)),rx,real(refline.receiver_field),imag(refline.receiver_field),mag, ...
    'VariableNames',{'frequency_hz','receiver_x_m','full_real','full_imag','magnitude'});
end
function a=local_magnitude(fc,r)
fc.receiver_x_m=r; fc.receiver_z_m=3; result=solve_full_kirchhoff_2d(fc); a=abs(result.receiver_field);
end
function density=local_density(fc,q)
k=2*pi*fc.frequency_hz/fc.c0_mps;
dx=fc.receiver_x_m-fc.surface_x_m; dz=fc.receiver_z_m-fc.surface_z_m; R=hypot(dx,dz);
G=1i/4*besselh(0,1,k*R); dnG=1i*k/4*besselh(1,1,k*R).* ...
    (dx.*fc.surface_normal(:,1)+dz.*fc.surface_normal(:,2))./R;
density=(G.*fc.surface_incident_dn+dnG.*fc.surface_incident).*q;
end
