function receipt=prepare_bellhop_C_ffp_comparison(project_root,ffp_root,output_dir)
% Reuse frozen fields only. No solver, interpolation, fit or coordinate changes.
assert(~isfolder(output_dir)); mkdir(output_dir); start=tic; diary(fullfile(output_dir,'diary.log')); guard=onCleanup(@()diary('off'));
[~,ram]=memory; assert(ram.PhysicalMemory.Available>=1.5*2^30,'Offline postprocess launch RAM gate');
Cdir=fullfile(project_root,'results','validation','bellhop_rotated_fig12b','20261116_linux_C_01'); Cfile=fullfile(Cdir,'remote_evidence','20261116_linux_C_01','C_resume','assembly','result.mat'); Ffile=fullfile(ffp_root,'results','validation','munk_ffp_5000m','20261112_full_01','J07.mat');
csha=black_file_hash(Cfile); fsha=black_file_hash(Ffile); assert(strcmpi(csha,'d1931b853a037cdf593c5bfee67a946f051281fe8a426bac2d0dfc8aa9a275b9')&&strcmpi(fsha,'8e123af0beff9904289a5419fdb82f74473beda075792b7b2f9819bd133ce3be'),'Frozen field SHA drift');
C=load(Cfile,'total','envelope','plan'); f=load(Ffile,'output'); F=f.output; clear f
checks=table(strings(0,1),false(0,1),'VariableNames',{'check','pass'}); cfg=C.plan.jobs(1).cfg;
add('exact_x_coordinates',isequal(C.plan.x_m(:),F.range_m(:))); add('exact_z_coordinates',isequal(C.plan.depths_m(:),F.depth_m(:))); add('same_source_mask',isequal(C.plan.source_singularity_mask,F.source_singularity_mask));
add('frequency',cfg.frequency_hz==F.cfg.frequency_hz); add('water_depth',cfg.water_depth_m==F.cfg.water_depth_m); add('source_depth',cfg.source_z_m==F.cfg.source_depth_m); add('Munk_parameters',isequal(cfg.munk_parameters,F.cfg.munk_parameters));
for name={'bottom_cp_mps','bottom_rho_gcm3','bottom_atten_db_wavelength'}, add(name{1},cfg.(name{1})==F.cfg.(name{1})); end
add('FFP_top_pressure_release',strcmp(F.cfg.top_type,'pressure_release')); add('FFP_bottom_fluid',strcmp(F.cfg.bottom_type,'fluid_halfspace')); add('FFP_water_rho1',F.cfg.water_rho_gcm3==1); valid=~F.source_singularity_mask;
add('FFP_pressure_normalization',max(abs(F.pressure_1m(valid)-4*pi*F.green(valid)))==0); add('finite_common_support',all(isfinite(C.total(valid)))&&all(isfinite(F.pressure_1m(valid))));
qual=jsondecode(fileread(fullfile(ffp_root,'results','validation','munk_ffp_5000m','20261112_full_01','audit_03','receipt.json'))); add('FFP_final_audit_pass',strcmp(qual.status,'REQUIRED_NUMERICAL_GATES_PASS')&&qual.failed_checks==0&&qual.final_job==7);
writetable(checks,fullfile(output_dir,'comparability_checks.csv')); assert(all(checks.pass),'Comparability gate failed: no alignment or refit allowed');
p_bellhop=C.total; p_ffp=F.pressure_1m; envelope=C.envelope; x=F.range_m(:).'; z=F.depth_m(:); source_mask=F.source_singularity_mask; clear F C
[zz,xx]=ndgrid(z,x); rr=hypot(xx,zz-4950); regions={valid,valid&zz<=4900,valid&zz>=4900&rr>=1,valid&xx==0&rr>=1,valid&rr<1,valid&abs(xx)<=1}; names={'all','bulk_z0_4900','near_bottom_r_ge1','axis_r_ge1','source_vicinity_r_lt1','central_abs_x_le1'};
metrics=table(strings(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),false(0,1),'VariableNames',{'region','points','complex_relative_L2_FFP','max_absolute_complex','max_complex_over_C_envelope','max_effective_phase_rad','phase_support','TL_difference_p95_dB','within_same_solver_gate_diagnostic_only'});
for j=1:numel(regions)
 use=regions{j}; delta=p_bellhop(use)-p_ffp(use); rel=norm(delta)/norm(p_ffp(use)); normalized=max(abs(delta)./envelope(use)); phase_use=use & abs(p_bellhop)>=.05*envelope & abs(p_ffp)>=.05*envelope; ph=abs(angle(p_bellhop(phase_use).*conj(p_ffp(phase_use)))); if isempty(ph), phase=NaN; tl95=NaN; else, phase=max(ph); db=sort(abs(20*log10(abs(p_bellhop(phase_use))./abs(p_ffp(phase_use))))); tl95=db(max(1,ceil(.95*numel(db)))); end
 metrics(end+1,:)={names{j},nnz(use),rel,max(abs(delta)),normalized,phase,nnz(phase_use),tl95,rel<=.02&&normalized<=.02&&phase<=.05}; %#ok<AGROW>
end
ix=find(x==0); [distance,ir]=min(abs(x-.1)); assert(distance<1e-12,'No original r~.1m grid column'); probe=[1000 3000 4000]; brightness=table(zeros(3,1),zeros(3,1),zeros(3,1),'VariableNames',{'depth_m','Bellhop_axis_minus_r01_TL_dB','FFP_axis_minus_r01_TL_dB'});
for j=1:3, iz=find(z==probe(j)); brightness(j,:)={probe(j),-20*log10(abs(p_bellhop(iz,ix))/abs(p_bellhop(iz,ir))),-20*log10(abs(p_ffp(iz,ix))/abs(p_ffp(iz,ir)))}; end
writetable(metrics,fullfile(output_dir,'common_point_metrics.csv')); writetable(brightness,fullfile(output_dir,'axis_brightness.csv'));
meta=struct('C_file',Cfile,'FFP_file',Ffile,'C_sha256',csha,'FFP_sha256',fsha,'source_reference','Fixed1m: Bellhop -conj(native), FFP4pi Green; no gain/phase fit','differences',{{'Bellhop Q-linear SSP versus FFP midpoint constant layers; both .625m sampling','Bellhop finite N16 Gaussian/ray field versus full layered-wave spectrum including evanescent content','Bellhop beam-count D convergence and full fluid-DtN comparison remain unqualified'}},'native_grid',[numel(z),numel(x)]);
save(fullfile(output_dir,'comparison.mat'),'p_bellhop','p_ffp','envelope','x','z','source_mask','meta','-v7.3');
receipt=struct('status','COMPARABLE_SAVED_FIELDS_COMPARED_NO_FIT','comparability_checks',height(checks),'comparability_failed',nnz(~checks.pass),'same_solver_window_exceeding_regions_diagnostic_only',nnz(~metrics.within_same_solver_gate_diagnostic_only),'metrics',table2struct(metrics),'brightness',table2struct(brightness),'meta',meta,'elapsed_s',toc(start),'new_FFP',0,'ABC_recomputed',0,'prior_remaining_s',135.83157788353276); receipt.remaining_budget_s=receipt.prior_remaining_s-receipt.elapsed_s;
fid=fopen(fullfile(output_dir,'checkpoint.json'),'w','n','UTF-8'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid); fprintf('C_FFP_COMPARE compatible%d failed_regions%d elapsed%.3fs remaining%.3fs\n',height(checks),receipt.same_solver_window_exceeding_regions_diagnostic_only,receipt.elapsed_s,receipt.remaining_budget_s); disp(metrics);
 function add(name,pass), checks(end+1,:)={string(name),logical(pass)}; end
end
