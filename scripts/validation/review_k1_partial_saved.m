function receipt=review_k1_partial_saved(source,output_dir,case_id,mode)
%REVIEW_K1_PARTIAL_SAVED Audit existing N7 arrays and compare saved R1 N12.
% No PS/BIE/FullFK solve, no remote calls; originals are never overwritten.
arguments
    source (1,:) char
    output_dir (1,:) char
    case_id (1,:) char = 'U6_N07'
    mode (1,:) char {mustBeMember(mode,{'full','numeric'})} = 'full'
end
assert(isfolder(source) && ~isfolder(output_dir)); mkdir(output_dir);
timer=tic; oldpath=path; pathguard=onCleanup(@()path(oldpath));
oldthreads=maxNumCompThreads(1); threadguard=onCleanup(@()maxNumCompThreads(oldthreads));
diary(fullfile(output_dir,'diary.log')); dg=onCleanup(@()diary('off'));
receipt=struct('status','RUNNING','new_acoustic_solves',0,'remote_connections',0,'failed_checks',0);
try
assert(ismember(case_id,{'U6_N07','U6_N40','U5_N12','U8_N12','U5_N40','U8_N40'}));
source_case=fullfile(source,'cases',case_id);
source_work=fullfile(source_case,'finalization','attempts','attempt_0001','work');
deps={'audit_r1_bie_manual_stage','audit_piston_pm_r1','audit_r2_pm_final', ...
    'summarize_r2_pm_fields','audit_r2_final_figures','render_k1_pm_stable_saved', ...
    'baffled_strip_piston_field_vertical','piston_independent_angular_spectrum_vertical', ...
    'closeout_s3a_pm_surface','r1_bie_continuation_gates'};
for j=1:numel(deps)
    file=which(deps{j}); frozen=fullfile(source_case,'code',[deps{j} '.m']);
    assert(strcmp(k1_pm_io('sha',file),k1_pm_io('sha',frozen)),'k1review:DependencyDrift','Audit dependency changed: %s',deps{j});
end
sentinels=fullfile(output_dir,'solve_sentinels'); mkdir(sentinels);
names={'solve_helmholtz_bie_halfplane_vertical','solve_helmholtz_bie_halfplane_cached_vertical', ...
    'closeout_s3a_pm_fullfk','run_k1_piston_pm_job','run_k1_pm_case', ...
    'run_r2_piston_pm_job','run_pe_1d_surface_reflection_validation', ...
    'run_k1_recovery_model_job','run_k1_pm_recovery_case'};
for j=1:numel(names)
    fid=fopen(fullfile(sentinels,[names{j} '.m']),'w'); assert(fid>=0);
    fprintf(fid,'function varargout=%s(varargin)\nerror(''k1review:NoSolve'',''Saved audit forbids new solves'');\nend\n',names{j}); fclose(fid);
end
addpath(sentinels,'-begin'); rehash;
for j=1:numel(names), assert(startsWith(which(names{j}),sentinels)); end
[~,sys]=memory; baseline=sys.PhysicalMemory.Available;
runtime=struct('MATLAB_version',version,'pid',feature('getpid'),'mode','existing local MATLAB MCP', ...
    'RAM_baseline_bytes',baseline,'reserve_fraction',.3,'single_thread',true,'wall_limit_s',180);
k1_pm_io('write',fullfile(output_dir,'runtime.json'),runtime);
work=fullfile(output_dir,'work'); mkdir(work);
files={'common.mat','comparison.mat','model_metrics.csv','receiver_fields_and_errors.csv','policy.json'};
versions=struct('path',{},'sha256',{});
for j=1:numel(files)
    p=fullfile(source_work,files{j}); copyfile(p,fullfile(work,files{j}));
    assert(strcmp(k1_pm_io('sha',p),k1_pm_io('sha',fullfile(work,files{j}))));
    versions(end+1)=struct('path',p,'sha256',k1_pm_io('sha',p)); %#ok<AGROW>
end
mkdir(fullfile(work,'input')); copyfile(fullfile(source_case,'input','*'),fullfile(work,'input'));
inventory=struct('common',{whos('-file',fullfile(work,'common.mat'))}, ...
    'comparison',{whos('-file',fullfile(work,'comparison.mat'))});
k1_pm_io('write',fullfile(output_dir,'MAT_native_inventory.json'),inventory);
q=load(fullfile(work,'common.mat')); a=load(fullfile(work,'comparison.mat'));
frozen=k1_pm_io('read',fullfile(source_case,'input','frozen_config.json'));
assert(q.cfg.modes==frozen.modes && q.cfg.U_PM_mps==frozen.U_PM_mps && q.cfg.seed==260001);
assert(q.cfg.source.total_width_m==.1 && q.cfg.source.frequency_hz==4000);
assert(isequal(q.x,(-205:205)'*.1953125) && all(q.z==3));
assert(isequal(q.profile,readtable(fullfile(source_case,'input','coefficients.csv'))));
ids={'bie_coarse','bie_main','bie_self','bie_space','bie_window'};
recipes=[60 10 96;60 12 96;60 12 64;60 14 96;70 12 96]; stages=cell(5,1);
for j=1:5
    local_guard(baseline,timer); dest=fullfile(work,'jobs',ids{j}); mkdir(dest);
    for file={'data.mat','points.csv'}
        p=fullfile(source_case,'jobs',ids{j},file{1}); copyfile(p,fullfile(dest,file{1}));
        assert(strcmp(k1_pm_io('sha',p),k1_pm_io('sha',fullfile(dest,file{1}))));
        versions(end+1)=struct('path',p,'sha256',k1_pm_io('sha',p)); %#ok<AGROW>
    end
    tick=tic; stages{j}=audit_r1_bie_manual_stage(dest,fullfile(work,'common.mat'),recipes(j,:));
    stages{j}.elapsed_s=toc(tick); assert(stages{j}.failed_checks==0);
    k1_pm_io('write',fullfile(output_dir,'numeric_checkpoint.json'),struct('completed_BIE_stage',j-1,'elapsed_s',toc(timer),'new_acoustic_solves',0));
end
local_guard(baseline,timer);
comprehensive=audit_piston_pm_r1(work); diary(fullfile(output_dir,'diary.log'));
final=audit_r2_pm_final(work);
spatial_name='N7_spatial'; if strcmp(mode,'numeric'), spatial_name='spatial'; end
spatial=summarize_r2_pm_fields({work},fullfile(output_dir,spatial_name));
assert(comprehensive.failed_checks==0 && final.failed_checks==0);
if strcmp(mode,'numeric')
    original=audit_r2_final_figures({work},fullfile(source_work,'figures_final'),fullfile(output_dir,'original_FIG_data_audit'),'original-data');
    assert(original.failed_checks==0 && original.checked_error_curves==12);
    k1_pm_io('write',fullfile(output_dir,'input_versions.json'),versions);
    receipt=struct('status','K1_CASE_INDEPENDENT_SAVED_NUMERIC_AND_FIG_DATA_PASS', ...
        'case_id',case_id,'failed_checks',0,'stages',{stages},'comprehensive',comprehensive, ...
        'reference',final,'original_FIG_curves',original.checked_error_curves, ...
        'new_acoustic_solves',0,'elapsed_s',toc(timer),'scope','Saved numeric and FIG data only; PNG/layout strict qualification remains execution-platform audit');
    k1_pm_io('write',fullfile(output_dir,'receipt.json'),receipt);
    fprintf('K1_CASE_SAVED_AUDIT_PASS case=%s elapsed_s=%.6f new_solves=0\n',case_id,receipt.elapsed_s);
    return
end
assert(strcmp(case_id,'U6_N07'),'k1review:Mode','Use numeric mode for remaining cases');

% Paired cutoff comparison uses the exact qualified historical anchor only.
anchor=fullfile(source,'anchor'); proof=k1_pm_io('read',fullfile(anchor,'anchor_proof.json'));
assert(strcmp(proof.status,'STRICT_R1_K0_ANCHOR_MATCH_PASS') && proof.failed_checks==0);
r=load(fullfile(anchor,'common.mat')); b=load(fullfile(anchor,'comparison.mat'));
assert(isequal(q.cfg.source,r.cfg.source) && isequal(q.x,r.x) && isequal(q.z,r.z));
assert(height(r.profile)==12 && isequal(q.profile,r.profile(1:7,:)));
assert(strcmp(q.cfg.phase_contract_version,r.cfg.phase_contract_version));
paired=table; field_change=table;
for index=1:2
    current=a; N=7; if index==2, current=b; N=12; end
    for m=[1 3]
        method="PS"; if m==3, method="FullFK"; end
        model=current.U(:,m); ref=current.U(:,2); e=model-ref;
        phase=angle(model.*conj(ref)); pe=4*abs(model).*abs(ref).*sin(phase/2).^2;
        amplitude=abs(model)-abs(ref);
        assert(max(abs(abs(e).^2-amplitude.^2-pe))<1e-10);
        [peak,ii]=max(abs(e)); row=current.metrics(string(current.metrics.method)==method,:);
        paired=[paired;table(6,N,N*2*pi/160,method,norm(e)/norm(ref), ...
            row.five_times_empirical_resolution,logical(row.resolved),sum(pe)/sum(abs(e).^2), ...
            norm(amplitude)/norm(ref),peak,q.x(ii),nnz(current.mask), ...
            'VariableNames',{'U_PM_mps','N','physical_Kmax_rad_per_m','method','complex_relative_L2', ...
            'five_times_empirical_resolution','resolved','phase_energy_fraction','amplitude_relative_L2', ...
            'peak_residual_Pa','peak_x_m','phase_valid_points'})]; %#ok<AGROW>
    end
end
for m=1:3
    method=["PS" "BIE" "FullFK"];
    field_change=[field_change;table(method(m),norm(a.U(:,m)-b.U(:,m))/norm(b.U(:,m)), ...
        'VariableNames',{'method','N7_vs_N12_field_change_relative_to_N12'})]; %#ok<AGROW>
end
writetable(paired,fullfile(output_dir,'paired_N7_N12_metrics.csv'));
writetable(field_change,fullfile(output_dir,'paired_field_changes.csv'));
save(fullfile(output_dir,'paired_saved_fields.mat'),'q','a','r','b','paired','field_change','-v7');
k1_pm_io('write',fullfile(output_dir,'numeric_receipt.json'),struct( ...
    'status','K1_N7_INDEPENDENT_SAVED_NUMERIC_PASS','failed_checks',0,'stages',{stages}, ...
    'comprehensive',comprehensive,'reference',final,'spatial_status',spatial.status, ...
    'new_acoustic_solves',0,'elapsed_s',toc(timer),'scope','N7 numerical only; separate figure gate follows'));
local_guard(baseline,timer);
original=audit_r2_final_figures({work},fullfile(source_work,'figures_final'),fullfile(output_dir,'original_FIG_data_audit'),'original-data');
render_k1_pm_stable_saved({work},fullfile(output_dir,'figures_N7'));
figures=audit_r2_final_figures({work},fullfile(output_dir,'figures_N7'),fullfile(output_dir,'new_figure_audit'));
assert(original.failed_checks==0 && figures.failed_checks==0);
k1_pm_io('write',fullfile(output_dir,'input_versions.json'),versions);
receipt=struct('status','K1_N7_INDEPENDENT_SAVED_AUDIT_PASS_PARTIAL_MATRIX', ...
    'failed_checks',0,'stages',{stages},'comprehensive',comprehensive,'reference',final, ...
    'spatial_status',spatial.status,'original_figure_curves',original.checked_error_curves, ...
    'new_figure_audit',figures,'new_acoustic_solves',0,'remote_connections',0, ...
    'elapsed_s',toc(timer),'source',source,'work',work,'solve_sentinels',numel(names), ...
    'scope','U6/N7 only independently reaudited; N12 existing qualified anchor; N40 absent; no wind/interaction inference');
k1_pm_io('write',fullfile(output_dir,'receipt.json'),receipt);
fprintf('K1_SAVED_REVIEW_PASS elapsed_s=%.6f stages=5 figures=12 new_solves=0\n',receipt.elapsed_s);
disp(paired);
catch ME
    receipt.status='FAILED_STOPPED'; receipt.failed_checks=1; receipt.error=ME.message; receipt.identifier=ME.identifier; receipt.stack=ME.stack; receipt.elapsed_s=toc(timer);
    k1_pm_io('write',fullfile(output_dir,'failure.json'),receipt);
    fprintf(2,'%s\n',getReport(ME,'extended','hyperlinks','off')); rethrow(ME);
end
end
function local_guard(baseline,timer)
[~,sys]=memory; assert(sys.PhysicalMemory.Available>=.3*baseline,'k1review:RAM','Saved audit RAM reserve');
assert(toc(timer)<180,'k1review:Wall','Saved audit wall limit');
end
