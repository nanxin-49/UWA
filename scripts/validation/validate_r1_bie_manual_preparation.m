function receipt=validate_r1_bie_manual_preparation(package_dir,output_dir)
%VALIDATE_R1_BIE_MANUAL_PREPARATION Local, explicitly solver-free dry-run.
arguments
    package_dir (1,:) char
    output_dir (1,:) char
end
assert(~isfolder(output_dir)); mkdir(output_dir); sentinel=fullfile(output_dir,'acoustic_call_sentinels'); mkdir(sentinel);
diary(fullfile(output_dir,'stdout.log')); cleanup=onCleanup(@()diary('off'));
blocked={'solve_helmholtz_bie_halfplane_vertical','solve_helmholtz_bie_halfplane_cached_vertical', ...
    'baffled_strip_piston_field_vertical','piston_independent_angular_spectrum_vertical', ...
    'closeout_s3a_pm_fullfk','run_pe_1d_surface_reflection_validation','besselh'};
for j=1:numel(blocked)
    fid=fopen(fullfile(sentinel,[blocked{j} '.m']),'w'); assert(fid>=0);
    fprintf(fid,'function varargout=%s(varargin)\nerror(''r1manual:AcousticForbidden'',''Acoustic call forbidden during dry-run'');\nend\n',blocked{j}); fclose(fid);
end
code_dir=fullfile(package_dir,'code'); addpath(code_dir,'-begin'); addpath(sentinel,'-begin');
path_guard=onCleanup(@()local_remove_paths(sentinel,code_dir));
try
    for j=1:numel(blocked), assert(startsWith(which(blocked{j}),sentinel)); end
    files={'run_r1_bie_manual','r1_bie_manual_io','audit_r1_bie_manual_stage','audit_r1_bie_manual_final'};
    analyzer=struct;
    for j=1:numel(files), analyzer.(files{j})=checkcode(fullfile(code_dir,[files{j} '.m']),'-id'); end
    r=run_r1_bie_manual(package_dir,'dry-run');
    assert(strcmp(r.status,'DRY_RUN_INPUT_CHECKS_PASS') && r.new_acoustic_solves==0);
    assert(~isfolder(fullfile(package_dir,'jobs')) && ~isfolder(fullfile(package_dir,'sessions')) && ...
        ~isfolder(fullfile(package_dir,'control')),'Dry-run created execution state');
    save(fullfile(output_dir,'dry_run.mat'),'r','analyzer','blocked','-v7');
    warnings=0; for j=1:numel(files), warnings=warnings+numel(analyzer.(files{j})); end
    receipt=struct('status','SOLVER_FREE_DRY_RUN_PASS','failed_checks',0,'acoustic_calls',0, ...
        'SSH_connections',0,'Code_Analyzer_messages',warnings,'blocked_acoustic_entrypoints',numel(blocked), ...
        'resume_tests',r.resume_selftest,'package_manifest_sha256',r.manifest_sha256);
    r1_bie_manual_io('write',fullfile(output_dir,'receipt.json'),receipt);
    fprintf('%s\n',jsonencode(receipt));
catch ME
    fid=fopen(fullfile(output_dir,'failure.txt'),'w'); fprintf(fid,'%s',getReport(ME,'extended','hyperlinks','off')); fclose(fid);
    rethrow(ME);
end
end
function local_remove_paths(sentinel,code_dir)
rmpath(sentinel); rmpath(code_dir);
end
