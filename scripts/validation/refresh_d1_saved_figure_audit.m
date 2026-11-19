function receipt=refresh_d1_saved_figure_audit(prior,out)
%REFRESH_D1_SAVED_FIGURE_AUDIT Recheck changed graphics; retain valid numeric audits.
assert(~isfolder(out)); mkdir(out); diary(fullfile(out,'diary.log')); guard=onCleanup(@()diary('off')); t0=tic;
old=jsondecode(fileread(fullfile(prior,'checkpoint.json'))); pm_sensitivity_io('verify',prior,old.artifacts_sha256);
receipt=jsondecode(fileread(fullfile(prior,'receipt.json'))); assert(strcmp(receipt.status,'SAVED_TX100_DEPTH_AUDIT_PASS') && receipt.failed_checks==0);
src=fullfile(prior,'audited_work'); work=fullfile(out,'audited_work'); mkdir(work);
listing=dir(fullfile(src,'**','*'));
for j=1:numel(listing)
    if listing(j).isdir, continue; end
    p=fullfile(listing(j).folder,listing(j).name); rel=p(numel(src)+2:end);
    if startsWith(strrep(rel,'\','/'),'figures/'), continue; end
    dest=fullfile(work,rel); parent=fileparts(dest); if ~isfolder(parent), mkdir(parent); end; copyfile(p,dest);
    assert(strcmp(pm_sensitivity_io('sha',p),pm_sensitivity_io('sha',dest)));
end
render_d1_depth_results(work); diary(fullfile(out,'diary.log'));
receipt.figure_audit=jsondecode(fileread(fullfile(work,'figures','figure_audit.json')));
assert(receipt.figure_audit.failed_checks==0 && receipt.figure_audit.differing_pixel_channels==0 && receipt.figure_audit.reopened_curves==24);
pm_sensitivity_io('verify',prior,old.artifacts_sha256);
receipt.numeric_audit_parent_sha256=pm_sensitivity_io('sha',fullfile(prior,'checkpoint.json'));
receipt.figure_refresh_s=toc(t0); receipt.audited_work_dir=work;
pm_sensitivity_io('write',fullfile(out,'receipt.json'),receipt);
repo=fileparts(fileparts(fileparts(mfilename('fullpath')))); files=dir(fullfile(work,'**','*')); paths={'receipt.json'};
for j=1:numel(files)
    if ~files(j).isdir, p=fullfile(files(j).folder,files(j).name); paths{end+1}=strrep(p(numel(out)+2:end),'\','/'); end %#ok<AGROW>
end
pm_sensitivity_io('write',fullfile(out,'checkpoint.json'),struct('status','AUDITED_PASS','numeric_parent_checkpoint_sha256',receipt.numeric_audit_parent_sha256, ...
    'code_versions_sha256',pm_sensitivity_io('manifest',repo,{'scripts/reporting/render_d1_depth_results.m','scripts/validation/refresh_d1_saved_figure_audit.m'}), ...
    'artifacts_sha256',pm_sensitivity_io('manifest',out,paths),'failed_or_unrun_checks',{{'Tx1490 not run'}}));
fprintf('%s\n',jsonencode(receipt));
end
