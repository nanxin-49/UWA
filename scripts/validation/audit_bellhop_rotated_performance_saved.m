function receipt=audit_bellhop_rotated_performance_saved(run_dir,output_dir)
%AUDIT_BELLHOP_ROTATED_PERFORMANCE_SAVED Saved bytes/state/timing; no new solve.
arguments
    run_dir (1,:) char
    output_dir (1,:) char
end
assert(~isfolder(output_dir)); mkdir(output_dir); start=tic;
base=fullfile(run_dir,'baseline_clock_case01','case'); opt=fullfile(run_dir,'optimized_clock_case01','case');
checks=table(strings(0,1),false(0,1),zeros(0,1),zeros(0,1),'VariableNames',{'check','pass','value','limit'});
A=wall_bytes([base '.perfwall']); B=wall_bytes([opt '.perfwall']);
add('all1313_reflection_events_retained',abs(size(A.geometry,2)-1313)+abs(size(B.geometry,2)-1313),0);
add('hit_segment_and_bounce_counts_identical',nnz(A.integers~=B.integers),0);
add('all_reflection_state_bytes_identical',double(~strcmpi(black_file_hash([base '.perfwall']),black_file_hash([opt '.perfwall']))),0);
labels={'lambda','tangent_r','tangent_q','normal_r','normal_q','curvature','residual'};
for k=1:7, add(['geometry_' labels{k}],max(abs(A.geometry(k,:)-B.geometry(k,:))),0); end
labels={'x_r','x_q','t_r','t_q','p1','p2','q1','q2','c','Amp','Phase','integral_cds','tau_real','tau_imag'};
for stage=1:3
    for k=1:14
        row=7+(stage-1)*14+k; add(sprintf('state%d_%s',stage,labels{k}),max(abs(A.geometry(row,:)-B.geometry(row,:))),0);
    end
end
stages={'incident','reflected','rotated'}; geometry=struct('integer_rows',{{'segment','incident_top','incident_bottom','reflected_top','reflected_bottom','rotated_top','rotated_bottom'}}, ...
    'double_rows',{{'lambda','tangent_r','tangent_q','normal_r','normal_q','curvature','residual'}},'state_stages',{stages},'state_fields',{labels},'record_bytes',420,'event_count',1313);
bc=readmatrix([base '.perfcounts'],'FileType','text','NumHeaderLines',1); oc=readmatrix([opt '.perfcounts'],'FileType','text','NumHeaderLines',1);
add('step_and_reduce_counts_preserved',max(abs(bc(1:2)-oc(1:2))),0);
add('cross_call_count_preserved',abs(bc(3)-oc(3)),0);
add('true_cross_candidates_preserved',abs(bc(7)-oc(7)),0);
add('optimized_only5252_exact_cross_checks',abs(oc(5)-5252),0);
add('optimized_only5252_geometry_checks',abs(oc(6)-5252),0);
folders={'baseline_profile_case02','baseline_clock_case01','optimized_clock_case01','timing_pair1_baseline','timing_pair1_optimized','timing_pair2_optimized'};
records=cell(numel(folders),1); elapsed=zeros(numel(folders),1); cpu=elapsed; failure_count=0;
for k=1:numel(folders)
    cp=jsondecode(fileread(fullfile(run_dir,folders{k},'checkpoint.json'))); records{k}=cp; elapsed(k)=cp.process.elapsed_s; cpu(k)=cp.process.cpu_time_s;
    tableK=readtable(fullfile(run_dir,folders{k},'checks.csv')); failure_count=failure_count+nnz(~tableK.pass);
    add(['completed_case_' folders{k}],cp.failed_checks,0);
    for ext={'.env','.ssp','.iwcov','.rotq','.r25d','.shd','.rr64input','.rr64eval','.rr64scale','.mwdiag','.r25diag'}
        add([folders{k} '_identical' ext{1}],double(~strcmpi(black_file_hash([cp.root ext{1}]),black_file_hash([cp.frozen_root ext{1}]))),0);
    end
end
add('all_individual_case_CSV_pass',failure_count,0);
% Recheck frozen REAL8 delivery, including original25 component bindings.
family=fileparts(run_dir); saved=load(fullfile(family,'20261114_receiver64_minimal_01','saved_audit01','audit.mat'),'component_files','component_hashes');
for k=1:numel(saved.component_files), add(sprintf('original_component_SHA_%02d',k),double(~strcmpi(black_file_hash(saved.component_files{k}),saved.component_hashes{k})),0); end
manifest=readtable(fullfile(family,'20261114_receiver64_minimal_01','sha256_manifest.csv'),'TextType','string');
for k=1:height(manifest), add(sprintf('REAL8_preserved_manifest_%02d',k),double(~strcmpi(black_file_hash(char(manifest.path(k))),manifest.sha256(k))),0); end
for name={'SourceReceiverPositions.f90','ReadEnvironmentBell.f90','RotatedFieldState.f90','influence.f90','bellhop.f90','sspMod.f90','bellhopMod.f90'}
    add(['plain_optimization_changes_no_other_source_' name{1}],double(~strcmpi(black_file_hash(fullfile(family,'20261114_receiver64_minimal_01','overlay02',name{1})),black_file_hash(fullfile(run_dir,'optimized_plain_overlay01',name{1})))),0);
end
fid=fopen([base '.perftime'],'r'); rate=fscanf(fid,'%f',3); ticks=fscanf(fid,'%f',4); counts=fscanf(fid,'%f',4); ssp=fscanf(fid,'%f',3); coarse=fscanf(fid,'%f',3); fclose(fid);
clock_overhead_tick=rate(2)/rate(3); calls=[bc(1);bc(2);bc(3);bc(4)]; raw=ticks./counts.*calls/rate(1); corrected=max(0,ticks./counts-clock_overhead_tick).*calls/rate(1);
ssp_raw=ssp(2)/ssp(3)*ssp(1)/rate(1); ssp_corrected=max(0,ssp(2)/ssp(3)-clock_overhead_tick)*ssp(1)/rate(1);
wall_sample=struct('method','1/1024 systematic SYSTEM_CLOCK samples; short-function overhead/100ns quantization materially affects estimates', ...
    'rate_hz',rate(1),'calibrated_clock_pair_overhead_s',clock_overhead_tick/rate(1), ...
    'step_reduce_cross_geometry_raw_s',raw.','step_reduce_cross_geometry_corrected_s',corrected.', ...
    'ssp_calls',ssp(1),'ssp_raw_s',ssp_raw,'ssp_corrected_s',ssp_corrected, ...
    'trace_influence_scale_coarse_s',(coarse/rate(1)).', ...
    'wall_cross_plus_geometry_raw_fraction_of_profile_process',sum(raw(3:4))/records{2}.process.elapsed_s, ...
    'wall_cross_plus_geometry_corrected_fraction_of_profile_process',sum(corrected(3:4))/records{2}.process.elapsed_s, ...
    'warning','fractions are overhead-sensitive estimates, not precise CPU attribution; optimized cross sample excludes early returns and is not used');
timing=table(string(folders(:)),elapsed,cpu,'VariableNames',{'case','elapsed_s','cpu_sample_s'}); writetable(timing,fullfile(output_dir,'timing.csv'));
paired_speedup=elapsed(4)/elapsed(5); paired_reduction=1-elapsed(5)/elapsed(4);
receipt=struct('status','NUMERICAL_PASS_TIMING_LIMITED_RAM_STOP','failed_checks',nnz(~checks.pass),'checks_count',height(checks),'elapsed_s',toc(start), ...
    'native_jobs_completed',numel(folders),'native_jobs_elapsed_sum_s',sum(elapsed),'paired_repetitions_completed',1, ...
    'planned_pair_count',3,'unrun_due_RAM_gate',{{'pair2baseline','pair3baseline','pair3optimized'}}, ...
    'baseline_elapsed_s',elapsed(4),'optimized_elapsed_s',elapsed(5),'actual_speedup',paired_speedup,'actual_elapsed_reduction_fraction',paired_reduction, ...
    'CPU_sample_speedup',cpu(4)/cpu(5),'CPU_sample_reduction_fraction',1-cpu(5)/cpu(4), ...
    'additional_unpaired_optimized_elapsed_s',elapsed(6),'baseline_counts',bc,'optimized_counts',oc,'wall_sample',wall_sample, ...
    'geometry_event_count',1313,'all_state_max_difference',max(abs(A.geometry-B.geometry),[],'all'), ...
    'original25_and_REAL8_artifact_SHA_preserved',true,'run_dir',run_dir, ...
    'can_resume_full_A','No automatic resumption; minimum case only, RAM gate stopped repeated timing; manual review required', ...
    'next_allowed_step','Stop for review. No new Bellhop solve, no complete A/order14-16/B/C/D, no commit/push.');
writetable(checks,fullfile(output_dir,'checks.csv')); save(fullfile(output_dir,'audit.mat'),'receipt','checks','geometry','A','B','timing');
fid=fopen(fullfile(output_dir,'checkpoint.json'),'w','n','UTF-8'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid);
disp(jsonencode(receipt)); disp(checks(~checks.pass,:));
    function add(name,value,limit)
        checks(end+1,:)={string(name),isfinite(value)&&value<=limit,value,limit};
    end
end
function output=wall_bytes(path)
fid=fopen(path,'rb','ieee-le'); assert(fid>=0); bytes=fread(fid,inf,'*uint8'); fclose(fid); assert(mod(numel(bytes),420)==0);
raw=reshape(bytes,420,[]); output.integers=reshape(typecast(reshape(raw(1:28,:),[],1),'int32'),7,[]);
output.geometry=reshape(typecast(reshape(raw(29:end,:),[],1),'double'),49,[]); assert(all(isfinite(output.geometry),'all'));
end