function checks=audit_deepsea3d_alpha25_domain_artifacts()
%AUDIT_DEEPSEA3D_ALPHA25_DOMAIN_ARTIFACTS Verify saved CSV against native MAT.
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
out=fullfile(root,'results','validation','deepsea3d_alpha25_domain');
a=load(fullfile(out,'domain.mat'),'state');s=a.state;
assert(s.complete && s.reference_ready && s.G1_status=="CONDITIONAL_G1_PASS");
assert(isequal(s.cfg,s.baseline.cfg) && ~s.WP2_allowed && s.global_G1=="OPEN");
ref=readtable(fullfile(out,'reference_runs.csv'));
fk=readtable(fullfile(out,'fk_runs.csv'));
metrics=readtable(fullfile(out,'physical_metrics.csv'));
reference=ref.reflect_real(end)+1i*ref.reflect_imag(end);
increment=ref.delta_real(end)+1i*ref.delta_imag(end);
full=fk.reflect_real(2)+1i*fk.reflect_imag(2);
delta=fk.delta_real(2)+1i*fk.delta_imag(2);
assert(abs(reference-s.runs(end).reflect_native)/abs(reference)<1e-13);
assert(abs(delta-s.fk(2).rough_correction)/abs(delta)<1e-13);
errors=[abs(full/reference-1),abs(delta/increment-1)];
rho=[abs(full/reference),abs(delta/increment)];
phase=[angle(full/reference),angle(delta/increment)];
stored=[metrics.reflected_complex,metrics.increment_complex];
stored_phase=[metrics.reflected_phase_rad,metrics.increment_phase_rad];
assert(max(abs(errors-stored))<1e-13 && max(abs(phase-stored_phase))<1e-13);
identity_error=max(abs(errors.^2-((rho-1).^2+4*rho.*sin(phase/2).^2)));
assert(identity_error<1e-16);
% Common physical coordinates extracted from each domain, with no rescaling.
h=.0625;x24=-24:h:24-h;x32=-32:h:32-h;
common24=x24(abs(x24)<=12);common32=x32(abs(x32)<=12);
assert(isequal(common24,common32));[xx,yy]=meshgrid(common24,common24);
[e24,ex24,ey24]=deepsea3d_surface(xx,yy,s.baseline.cfg.profile);
[e32,ex32,ey32]=deepsea3d_surface(xx,yy,s.cfg.profile);
geometry_error=max(abs([e24(:)-e32(:);ex24(:)-ex32(:);ey24(:)-ey32(:)]));
assert(geometry_error==0);
[outside,dx,dy]=deepsea3d_surface([12 16 24 32],[0 0 0 0],s.cfg.profile);
assert(all([outside dx dy]==0));
checks=struct('CSV_native_complex_matches',true,'raw_phase_matches',true, ...
    'amplitude_phase_identity_abs_error',identity_error,'common_surface_and_slopes_max_difference',geometry_error, ...
    'outside_support_flat',true,'physical_cfg_unchanged',true,'WP2_allowed',false,'global_G1',"OPEN");
fid=fopen(fullfile(out,'saved_data_checks.json'),'w');assert(fid>=0);
clean=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(checks,PrettyPrint=true));
disp(checks);
end
