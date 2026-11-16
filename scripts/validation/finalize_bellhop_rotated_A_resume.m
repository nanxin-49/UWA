function receipt=finalize_bellhop_rotated_A_resume(run_dir,figure_dir)
%FINALIZE_BELLHOP_ROTATED_A_RESUME StageA receipt; no solver/next group dispatch.
start=tic; pre=load(fullfile(run_dir,'preflight.mat'),'receipt','plan'); p=pre.receipt; plan=pre.plan;
state=jsondecode(fileread(fullfile(run_dir,'checkpoint.json')));
assert(strcmp(state.status,'ASSEMBLED_A_RAW_AUDIT_PENDING'),'Stage already finalized; do not duplicate accounting.');
raw=jsondecode(fileread(fullfile(run_dir,'raw_audit01','checkpoint.json')));
tail=jsondecode(fileread(fullfile(run_dir,'N12_to_N16_audit01','checkpoint.json')));
boundary=jsondecode(fileread(fullfile(run_dir,'boundary_audit01','checkpoint.json')));
fig=jsondecode(fileread(fullfile(figure_dir,'checkpoint.json')));
figaudit=jsondecode(fileread(fullfile(run_dir,'figure_audit02','checkpoint.json')));
assert(state.completed_count==39 && state.new_native_count==8 && raw.failed_checks==0 && tail.failed_checks==0 && boundary.failed_checks==0 && figaudit.failed_checks==0);
timing=table(strings(39,1),false(39,1),zeros(39,1),zeros(39,1),zeros(39,1),'VariableNames',{'tag','reused','original_native_s','native_peak_MiB','component_wall_s'});
for j=1:39
    c=state.components(j); assert(strcmpi(black_file_hash(c.path),c.component_sha256)); cp=jsondecode(fileread(c.checkpoint_path));
    if j<=31
        file=fullfile(fileparts(run_dir),'20261114_full_01','A_Q125_h1_B401',c.tag,'checkpoint.json');
        if isfile(file), cp=jsondecode(fileread(file)); end
    end
    timing.tag(j)=string(c.tag); timing.reused(j)=j<=31; timing.original_native_s(j)=cp.process.elapsed_s;
    timing.native_peak_MiB(j)=cp.process.peak_working_set_bytes/2^20; timing.component_wall_s(j)=cp.elapsed_s;
end
original25=dir(fullfile(fileparts(run_dir),'20261114_full_01','A_Q125_h1_B401','*','component.mat')); assert(numel(original25)==25);
for k=1:25
    path=fullfile(original25(k).folder,original25(k).name); cp=jsondecode(fileread(fullfile(original25(k).folder,'checkpoint.json'))); match=find(strcmp({cp.bindings.path},path)); assert(isscalar(match)); assert(strcmpi(black_file_hash(path),cp.bindings(match).sha256));
end
writetable(timing,fullfile(run_dir,'formal_component_timing.csv'));
extra=raw.elapsed_s+tail.elapsed_s+boundary.elapsed_s+fig.elapsed_s+figaudit.elapsed_s;
state.resume_elapsed_s=state.resume_elapsed_s+extra+toc(start); state.cumulative_accounted_s=state.prior_used_and_reserved_s+state.resume_elapsed_s;
state.status='A_COMPLETE_STAGE_GATES_PASS_STOPPED_FOR_REVIEW'; state.next_allowed_step='Review only. No B/C/D, no commit/push/merge; same original3600s budget if later authorized.';
state.raw_audit=raw; state.order_convergence=tail; state.boundary_audit=boundary; state.figure_audit=figaudit;
state.remaining_budget_s=3600-state.cumulative_accounted_s; assert(state.remaining_budget_s>0);
bytes=0; roots=p.disk_roots; roots{end+1}=figure_dir;
for k=1:numel(roots), d=dir(fullfile(roots{k},'**','*')); bytes=bytes+sum([d(~[d.isdir]).bytes]); end
state.output_bytes=bytes; assert(bytes<=15*2^30); [~,mem]=memory; state.available_RAM_after_bytes=mem.PhysicalMemory.Available;
fid=fopen(fullfile(run_dir,'checkpoint.json'),'w','n','UTF-8'); fprintf(fid,'%s',jsonencode(state,PrettyPrint=true)); fclose(fid);
A_native=sum(timing.original_native_s);
receipt=struct('status',state.status,'completed_formal',39,'reused_formal',31,'new_formal',8, ...
    'original25_SHA_unchanged',true,'new_native_elapsed_s',state.new_native_elapsed_s,'whole_A_saved_native_elapsed_s',A_native, ...
    'new_native_peak_MiB',max(timing.native_peak_MiB(~timing.reused)),'raw_audit_checks',raw.checks_count, ...
    'tail_failed_regions',tail.failed_checks,'figure_audit_checks',figaudit.checks_count,'top_residual',boundary.top_pressure_release_max_over_envelope, ...
    'accounted_elapsed_s',state.cumulative_accounted_s,'remaining_original_budget_s',state.remaining_budget_s,'output_GiB',bytes/2^30, ...
    'available_RAM_after_GiB',mem.PhysicalMemory.Available/2^30,'figure_dir',figure_dir,'field_file',fullfile(run_dir,'assembly','result.mat'), ...
    'field_sha256',black_file_hash(fullfile(run_dir,'assembly','result.mat')),'unrun',{{'B','C','D','full-grid independentFFP comparison','paper claim'}}, ...
    'next_group_plan',struct('jobs_per_group',39,'B_Q_m',.625,'B_step_m',1,'B_beams',401,'C_Q_m',.625,'C_step_m',.5,'C_beams',401,'D_Q_m',.625,'D_step_m',.5,'D_beams',801, ...
        'rough_native_seconds_B',[1.5 2]*A_native,'rough_native_seconds_C',[1.8 2.6]*A_native,'rough_native_seconds_D',[3.3 5]*A_native, ...
        'native_peak_MiB_estimate',[80 160],'MATLAB_extra_peak_MiB_estimate',[300 600],'extra_output_GiB_per_group',[1.5 2.5], ...
        'caution','Planning ranges, not measuredB/C/D. Add audit/I-O costs; upper total may exceed original remaining budget. Review stagedB only, not automatic dispatch.'));
fid=fopen(fullfile(run_dir,'checkpoint_stage.json'),'w','n','UTF-8'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid); disp(jsonencode(receipt));
end