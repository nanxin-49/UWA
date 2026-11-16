function receipt=summarize_bellhop_rotated_C_saved(project_root,run_dir,output_dir)
% Descriptive absolute/envelope differences from saved fields; no new gates/solve.
assert(~isfolder(output_dir)); mkdir(output_dir); start=tic; diary(fullfile(output_dir,'diary.log')); guard=onCleanup(@()diary('off'));
t=fileread('/proc/meminfo'); ram=1024*str2double(regexp(t,'MemAvailable:\s*(\d+)','tokens','once')); assert(ram>=1.5*2^30,'Post-audit RAM gate');
previous=jsondecode(fileread(fullfile(run_dir,'C_resume','checkpoint.json'))); initial=previous.remaining_budget_s; assert(initial>30);
fine=fullfile(run_dir,'C_resume','assembly','result.mat'); coarse={fullfile(project_root,'results','validation','bellhop_rotated_fig12b','20261115_B_run_02','assembly','result.mat'),fullfile(run_dir,'C_resume','assembly','sum_through_order12.mat')}; labels={'B_to_C','C_N12_to_N16'};
metrics=table(strings(0,1),strings(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),'VariableNames',{'comparison','region','relative_L2','max_complex_absolute','max_amplitude_absolute','max_complex_over_fine_envelope','max_envelope_relative','max_effective_phase_rad'});
for k=1:2
 a=load(coarse{k},'total','envelope','plan'); b=load(fine,'total','envelope','plan'); assert(isequal(a.plan.depths_m,b.plan.depths_m)&&isequal(a.plan.x_m,b.plan.x_m));
 [z,x]=ndgrid(b.plan.depths_m(:),b.plan.x_m(:)); r=hypot(x,z-4950); valid=~b.plan.source_singularity_mask; regions={valid,valid&z<=4900,valid&z>=4900&r>=1,valid&x==0&r>=1,valid&r<1}; names={'all','bulk_z0_4900','near_bottom_r_ge1','axis_r_ge1','source_vicinity_r_lt1'};
 for j=1:5
  use=regions{j}; d=a.total(use)-b.total(use); ph=use&abs(b.total)>=.05*b.envelope; phase=max(abs(angle(a.total(ph).*conj(b.total(ph))))); if isempty(phase), phase=0; end
  metrics(end+1,:)={labels{k},names{j},norm(d)/norm(b.total(use)),max(abs(d)),max(abs(abs(a.total(use))-abs(b.total(use)))),max(abs(d)./b.envelope(use)),max(abs(a.envelope(use)-b.envelope(use))./b.envelope(use)),phase}; %#ok<AGROW>
 end
 clear a b z x r regions
end
writetable(metrics,fullfile(output_dir,'absolute_envelope_metrics.csv')); job_wall=posixtime(datetime('now','TimeZone','UTC'))-str2double(getenv('BL_JOB_START_EPOCH')); receipt=struct('status','SAVED_C_DESCRIPTIVE_AUDIT_COMPLETE','new_native_processes',0,'elapsed_s',toc(start),'job_wall_s',job_wall,'previous_remaining_budget_s',initial,'remaining_original_budget_s',initial-job_wall,'accounted_cumulative_s',3600-initial+job_wall,'new_gate_or_tolerance',false,'envelope_definition','E=sum(abs(component)); envelope_relative=abs(Ecoarse-Efine)/Efine; absolute pressure uses fixed1m source reference');
fid=fopen(fullfile(output_dir,'checkpoint.json'),'w','n','UTF-8'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid); disp(jsonencode(receipt)); disp(metrics);
end
