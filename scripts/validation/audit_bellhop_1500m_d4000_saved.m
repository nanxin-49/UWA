function receipt=audit_bellhop_1500m_d4000_saved(project_root,run_dir)
%AUDIT_BELLHOP_1500M_D4000_SAVED Independent SHD/CSV/geometry audit; no solver.
timer=tic; diary(fullfile(run_dir,'saved_audit_diary.log')); cleaner=onCleanup(@()diary('off'));
addpath(fullfile(project_root,'scripts','validation','support'));
r=fullfile(run_dir,'returned'); root=fullfile(r,'jobs','01_D_4000_C','case');
s=load(fullfile(r,'jobs','01_D_4000_C','component.mat')); plan=load(fullfile(run_dir,'plan.mat'));
assert(isequaln(s.cfg,plan.frozen.cfg)); assert(s.proc.exit_code==0&&isempty(s.proc.stderr));
assert(~contains(lower(fileread([root '.prt'])),'fatal'));
raw=read_bellhop_rotated_field(root); fid=fopen([root '.shd'],'rb'); assert(fid>=0);
p_native=complex(zeros(5,1));
for k=1:5
 fseek(fid,(raw.range_record+k)*raw.record_bytes,'bof'); a=fread(fid,2,'float32'); p_native(k)=complex(a(1),a(2));
end
fclose(fid); assert(isequal(p_native,raw.pressure));
x=[0 10 50 100 200]'; L=hypot(x,1484); tau=L/1500; ref=exp(1i*2*pi*4000*tau)./L;
p=-conj(p_native); assert(isequal(p(:),s.p(:))); E=1./L;
support=abs(p)>=.01*E; residual=p-ref; phase=angle(p.*conj(ref));
metrics=struct('complex_L2',norm(residual)/norm(ref),'max_complex_difference_over_E',max(abs(residual)./E),'max_phase_rad',max(abs(phase(support))),'phase_support',nnz(support));
t=readtable(fullfile(r,'jobs','01_D_4000_C','field.csv'),'Delimiter',','); assert(isequal(t.x_m,x)&&max(abs(t.p_real+1i*t.p_imag-p))<1e-17);
fields=table(x,L,tau,real(p),imag(p),real(ref),imag(ref),abs(residual),abs(residual)./E,phase,'VariableNames',{'x_m','analytic_L_m','analytic_delay_s','Bellhop_p_real','Bellhop_p_imag','reference_p_real','reference_p_imag','absolute_complex_difference','complex_difference_over_E','phase_difference_rad'});
writetable(fields,fullfile(run_dir,'direct_field_independent.csv'));
d=readmatrix([root '.r25diag'],'FileType','text','CommentStyle','#'); assert(size(d,1)==201&&all(d(:,9:10)==0,'all'));
native_L=hypot(d(:,2),d(:,3)); ray_time_error=d(:,4)-native_L/1500; ray_direction_error=atan2(d(:,3),d(:,2))-deg2rad(d(:,1));
assert(max(abs(ray_time_error))<=1e-7&&max(abs(ray_direction_error))<=1e-6);
axis_index=find(abs(d(:,1))<1e-12); assert(isscalar(axis_index));
axis_crossing_time=d(axis_index,4)*1484/d(axis_index,2); assert(abs(axis_crossing_time-1484/1500)<=1e-7);
ray_end=table(d(:,1),d(:,2),d(:,3),native_L,d(:,4),ray_time_error,ray_direction_error,'VariableNames',{'launch_angle_deg','end_rotated_r_m','end_rotated_q_m','geometric_end_length_m','native_end_tau_s','end_tau_error_s','direction_error_rad'});
writetable(ray_end,fullfile(run_dir,'native_ray_terminal_geometry.csv'));
checks=readtable(fullfile(r,'C_checks.csv'),'Delimiter',',','TextType','string'); assert(nnz(~checks.pass)==2);
expected_failed=["slot1_analytic_complex_L2";"slot1_analytic_max_complex_over_E"]; assert(isequal(checks.check(~checks.pass),expected_failed));
for ext={'.env','.ssp','.r25d'}
 assert(strcmpi(black_file_hash([root ext{1}]),black_file_hash(fullfile(plan.frozen.local_dir,['case' ext{1}]))));
end
protected=readtable(fullfile(project_root,'results','validation','bellhop_1500m_applicability','20261117_design_01','protected_5000m_artifacts_sha256.csv'),'Delimiter',',','TextType','string');
for k=1:height(protected), assert(strcmpi(black_file_hash(char(fullfile(project_root,protected.path(k)))),protected.sha256(k))); end
artifact=table(["SHD";"ARR";"RAY"],[isfile([root '.shd']);isfile([root '.arr']);isfile([root '.ray'])],["C field: required";"A mode: not run";"R/E mode: not run"],'VariableNames',{'artifact','present','interpretation'});
writetable(artifact,fullfile(run_dir,'native_output_inventory.csv'));
fig=figure('Visible','off','Color','w','Position',[100 100 1100 680]); layout=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
a=nexttile(layout); plot(a,x,abs(ref),'k-o',x,abs(p),'b--s'); xlabel(a,'Independent link offset x (m)'); ylabel(a,'Pressure / 1 m reference'); legend(a,{'Analytic direct','Bellhop C'},'Location','best'); grid(a,'on');
a=nexttile(layout); plot(a,x,abs(residual)./E,'b-o'); yline(a,.002,'r--','Frozen limit'); xlabel(a,'x (m)'); ylabel(a,'|complex difference| / E'); grid(a,'on');
a=nexttile(layout); plot(a,x,phase,'b-o'); yline(a,.01,'r--'); yline(a,-.01,'r--'); xlabel(a,'x (m)'); ylabel(a,'Phase difference (rad)'); ylim(a,[-.012 .012]); grid(a,'on');
a=nexttile(layout); plot(a,d(:,1),ray_time_error*1e6,'k-'); xlabel(a,'Native launch angle (deg)'); ylabel(a,'Ray terminal time error (microseconds)'); grid(a,'on');
title(layout,'U0 direct only: 4000 Hz, 201 beams; complex-field gates FAILED'); drawnow;
savefig(fig,fullfile(run_dir,'D4000_C_diagnostic.fig')); exportgraphics(fig,fullfile(run_dir,'D4000_C_diagnostic.png'),'Resolution',150); close(fig);
fig= openfig(fullfile(run_dir,'D4000_C_diagnostic.fig'),'invisible'); objs=findall(fig,'Type','line'); data_ok=false;
for k=1:numel(objs), if isequal(objs(k).XData,x')&&isequal(objs(k).YData,(abs(residual)./E)'), data_ok=true; end; end
close(fig); assert(data_ok,'Editable FIG numerical curve mismatch');
receipt=struct('status','D4000_C_NUMERICAL_FAIL_STOPPED','started_slots',1,'native_exit_code',s.proc.exit_code,'native_elapsed_s',s.proc.elapsed_s,'failed_numerical_gates',2,'geometry_coordinate_checks','PASS limited to native straight-ray end states and actual field coordinates','metrics',metrics,'axis_native_crossing_time_ms',1000*axis_crossing_time,'max_native_ray_length_error_m',max(abs(native_L-1500*d(:,4))),'max_native_ray_time_error_s',max(abs(ray_time_error)),'max_native_ray_direction_error_rad',max(abs(ray_direction_error)),'SHD_present',true,'ARR_present',false,'RAY_present',false,'native_artifact_scope','C writes SHD; A/R/E not authorized or run','frozen_ENV_SSP_R25D_SHA_equal',true,'protected_historical_SHA_checks',height(protected),'protected_historical_SHA_failures',0,'independent_SHD_CSV_reconstruction','PASS','FIG_reopen_data','PASS','elapsed_s',toc(timer),'next_allowed_step','manual review; no retry, beam refinement or P2','matlab_version',version);
save(fullfile(run_dir,'independent_audit.mat'),'receipt','fields','ray_end','p_native','p','ref','metrics','checks','artifact');
fid=fopen(fullfile(run_dir,'independent_audit.json'),'w'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid);
fprintf('SAVED_AUDIT FAIL2 relL2%.10g maxDE%.10g phase%.10g rays201 axis_ms%.12f protected5PASS FIGPASS elapsed%.3fs\n',metrics.complex_L2,metrics.max_complex_difference_over_E,metrics.max_phase_rad,1000*axis_crossing_time,toc(timer));
end