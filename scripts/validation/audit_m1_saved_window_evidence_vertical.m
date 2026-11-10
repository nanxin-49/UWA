function receipt=audit_m1_saved_window_evidence_vertical(output_dir)
%AUDIT_M1_SAVED_WINDOW_EVIDENCE_VERTICAL Read saved evidence; never propagate.
% Only named current result files are loaded. No archive discovery or restore.
arguments
    output_dir (1,:) char
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
assert(java.io.File(output_dir).isAbsolute(),'Absolute output directory required.');
assert(startsWith(char(java.io.File(output_dir).getCanonicalPath()), ...
    fullfile(root,'results','validation','m1_preflight')),'M1 preflight directory required.');
assert(~isfolder(output_dir),'Never overwrite an audit.'); mkdir(output_dir);
diary(fullfile(output_dir,'stdout.log')); cleanup=onCleanup(@()diary('off')); %#ok<NASGU>
timer=tic; receipt=struct('status','FAILED','new_propagation_calls',0, ...
    'matlab_release',version,'exit_code_kind','MATLAB exception status; existing shared session');
try
    [~,ram]=memory; assert(500*2^20<.7*ram.PhysicalMemory.Available,'RAM gate');
    receipt.available_physical_bytes=ram.PhysicalMemory.Available;
    files={ ...
        'results/validation/pe_reflected_chain_window/reflected_chain_window_validation.mat', ...
        'results/validation/pe_reflected_chain_window/reflected_chain_wideband_validation.mat', ...
        'results/validation/pe_random_surface_window_robustness/random_surface_window_robustness_4khz.mat', ...
        'results/validation/pe_gaussian_window_sponge/gaussian_window_convergence.mat', ...
        'results/validation/pe_gaussian_window_sponge/gaussian_sponge_wideband.mat', ...
        'results/validation/pe_gaussian_window_sponge/extended_window/gaussian_window_convergence.mat'};
    evidence=struct(); checks=struct('name',{},'value',{},'limit',{},'pass',{});
    s=load(fullfile(root,files{1}),'validation'); v=s.validation;
    evidence.fixed=struct('config',v.config,'windows_m',[v.cases.width_actual_m], ...
        'nx',[v.cases.nx],'edge_pass',[v.cases.reference_edge_valid], ...
        'max_edge5',[v.cases.reference_max_edge5_fraction], ...
        'max_boundary_db',[v.cases.reference_max_boundary_db], ...
        'upward_saved_planes',numel(v.cases(end).incident_energy_trace.z_m), ...
        'upward_total_steps',ceil(100/(1500/4000*.5)), ...
        'has_code_revision',isfield(v,'code_revision'),'has_run_sha',isfield(v,'code_versions_sha256'));
    t=readtable(fullfile(root,'results','validation','pe_reflected_chain_window','single_frequency_cases.csv'));
    checks(end+1)=gate('fixed_MAT_CSV_edge',max(abs(t.max_edge5_fraction-[v.cases.reference_max_edge5_fraction]')),1e-15);
    checks(end+1)=gate('fixed_192_strict',~v.cases(end).reference_edge_valid,0);
    checks(end+1)=gate('fixed_160_failure_retained',v.cases(end-1).reference_edge_valid,0);
    master=load(fullfile(root,'results','validation','pe_reflected_chain_window','fixed_master_surface.mat'), ...
        'master_eta','master_meta');
    for j=1:numel(v.cases)
        n=v.cases(j).nx; first=(size(master.master_eta,1)-n)/2+1;
        crop=master.master_eta(first:first+n-1,first:first+n-1);
        checks(end+1)=gate(sprintf('fixed_master_exact_crop_%d',n), ...
            ~isequal(crop,v.cases(j).surface_elevation_xy),0); %#ok<AGROW>
    end
    evidence.fixed.master_surface_meta=master.master_meta;
    clear master crop
    % Reconstruct native direct from saved physical components, independently
    % of the current helper. The old direct-DSP array should equal native d.
    c=v.cases(end); k=2*pi*4000/1500;
    pd=c.H_total_physical-c.H_reflect_physical;
    checks(end+1)=gate('saved_native_direct_identity',abs(c.H_direct-pd*exp(-1i*k*97)),1e-12);
    checks(end+1)=gate('saved_native_reflect_identity',abs(c.H_reflect_reduced-c.H_reflect_physical*exp(-1i*k*103)),1e-12);
    evidence.fixed.old_vs_v2_direct_difference=abs(c.H_direct-conj(pd)*exp(1i*k*97));
    clear s v c t
    s=load(fullfile(root,files{2}),'validation'); v=s.validation;
    evidence.fixed_wideband=struct('config',v.config,'case_names',{{v.cases.name}}, ...
        'summary',table2struct(v.summary_table),'has_code_revision',isfield(v,'code_revision'));
    reference=v.cases(v.reference_index); axis=reference.f_axis_hz(:);
    % The correct v2 total is globally conjugated physical total, not a
    % path-independent phase rotation of the old carrier-only total.
    current=conj(reference.H_total_physical_f(:)).*exp(1i*2*pi*axis*97/1500);
    evidence.fixed_wideband.old_vs_v2_total_relative=norm(current-reference.H_f(:))/norm(current);
    checks(end+1)=gate('saved_wideband_axis',max(abs(axis-(3000:62.5:5000)')),0);
    clear s v reference axis current
    s=load(fullfile(root,files{3}),'validation'); v=s.validation;
    evidence.random=struct('config',v.config,'case_count',numel(v.cases), ...
        'pass_rates',table2struct(v.pass_rate_table),'all_references_valid',v.all_references_valid, ...
        'upward_saved_planes',numel(v.cases(end).incident_energy_trace.z_m), ...
        'has_code_revision',isfield(v,'code_revision'),'has_run_sha',isfield(v,'code_versions_sha256'));
    rows=v.comparison_table; mask=rows.window_nominal_m==192;
    checks(end+1)=gate('random_192_all_pass',sum(~rows.passed(mask)),0);
    checks(end+1)=gate('random_160_all_fail',sum(rows.passed(rows.window_nominal_m==160)),0);
    checks(end+1)=gate('random_reference_validity',~v.all_references_valid,0);
    checks(end+1)=gate('random_case_count',abs(numel(v.cases)-45),0);
    clear s v rows mask
    s=load(fullfile(root,files{4}),'validation'); v=s.validation;
    evidence.direct=struct('config',v.config,'hard_max_complex',max(v.hard_check_table.full_complex_l2), ...
        'hard_max_center',max(v.hard_check_table.center_complex_l2_max), ...
        'reference_width_m',v.reference_width_m,'hard_pass',v.hard_check_passed, ...
        'window_converged',v.window_converged,'has_code_revision',isfield(v,'code_revision'), ...
        'actual_dx_range_m',[min(v.window_table.dx_m),max(v.window_table.dx_m)]);
    t=readtable(fullfile(root,'results','validation','pe_gaussian_window_sponge','gaussian_pe_as_hard_check.csv'));
    checks(end+1)=gate('direct_MAT_CSV_hard',max(abs(t.full_complex_l2-v.hard_check_table.full_complex_l2)),1e-25);
    clear s v t
    evidence.direct_initial_scan=evidence.direct;
    s=load(fullfile(root,files{6}),'validation'); v=s.validation;
    evidence.direct=struct('config',v.config,'hard_max_complex',max(v.hard_check_table.full_complex_l2), ...
        'hard_max_center',max(v.hard_check_table.center_complex_l2_max), ...
        'reference_width_m',v.reference_width_m,'hard_pass',v.hard_check_passed, ...
        'window_converged',v.window_converged,'has_code_revision',isfield(v,'code_revision'));
    t=readtable(fullfile(root,'results','validation','pe_gaussian_window_sponge','extended_window','gaussian_pe_as_hard_check.csv'));
    checks(end+1)=gate('extended_direct_MAT_CSV',max(abs(t.full_complex_l2-v.hard_check_table.full_complex_l2)),1e-25);
    checks(end+1)=gate('extended_direct_window_pass',~v.window_converged,0);
    clear s v t
    s=load(fullfile(root,files{5}),'validation'); v=s.validation;
    evidence.direct_wideband=struct('config',v.config,'summary',table2struct(v.summary_table), ...
        'has_code_revision',isfield(v,'code_revision'));
    clear s v
    evidence.provenance_limitations={ ...
        'Historical MAT structs have configuration but no bound run-time Git revision/source SHA.', ...
        'Git history identifies recoverable script versions, not exact dirty run source.', ...
        'Native physical fields remain interpretable; old carrier-only DSP/timing claims require reconstruction.', ...
        'Edge traces cover saved z planes (65/33), not every marching step.', ...
        'The root Gaussian MAT is the initial failed <=100m scan; the accepted gate is under extended_window, as bound by the sponge MAT.', ...
        'No historical CPU wall-time fields: M1 runtime forecasts must be labelled estimates.'};
    receipt.checks=checks; receipt.failed_checks=sum(~[checks.pass]);
    assert(receipt.failed_checks==0,'Saved evidence identity checks failed');
    receipt.status='READONLY_AUDIT_PASS_WITH_PROVENANCE_LIMITS'; receipt.exit_code=0;
    receipt.evidence=evidence;
    writetable(struct2table(checks),fullfile(output_dir,'checks.csv'));
    fid=fopen(fullfile(output_dir,'stderr.log'),'w'); fprintf(fid,'No MATLAB exception.\n'); fclose(fid);
catch exception
    receipt.exit_code=1; receipt.error_identifier=exception.identifier; receipt.error_message=exception.message;
    fid=fopen(fullfile(output_dir,'stderr.log'),'w'); fprintf(fid,'%s',getReport(exception,'extended','hyperlinks','off')); fclose(fid);
end
receipt.elapsed_s=toc(timer); save(fullfile(output_dir,'audit.mat'),'receipt');
fid=fopen(fullfile(output_dir,'audit.json'),'w','n','UTF-8'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid);
short=rmfield_if_present(receipt,'evidence'); short=rmfield_if_present(short,'checks'); disp(short);
end
function c=gate(name,value,limit)
c=struct('name',name,'value',double(value),'limit',limit,'pass',isfinite(value)&&value<=limit);
end
function s=rmfield_if_present(s,name)
if isfield(s,name),s=rmfield(s,name);end
end
