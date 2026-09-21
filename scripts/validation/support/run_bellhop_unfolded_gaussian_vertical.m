function result = run_bellhop_unfolded_gaussian_vertical(cfg)
%RUN_BELLHOP_UNFOLDED_GAUSSIAN_VERTICAL Run one unfolded Bellhop case.
arguments
    cfg (1,1) struct
end
required={'bellhop_exe','case_root','run_type'};
for ii=1:numel(required)
    if ~isfield(cfg,required{ii}), error('Missing cfg.%s.',required{ii}); end
end
if exist(cfg.bellhop_exe,'file')~=2, error('Bellhop executable is missing.'); end
[env_file,sbp_file]=write_bellhop_unfolded_gaussian_env_vertical(cfg.case_root,cfg);
extensions={'.arr','.shd','.ray','.prt'};
for ii=1:numel(extensions)
    target=[cfg.case_root extensions{ii}]; if exist(target,'file')==2, delete(target); end
end
out_dir=fileparts(cfg.case_root); old=pwd; cleanup=onCleanup(@()cd(old)); cd(out_dir);
[~,name]=fileparts(cfg.case_root);
[status,command_output]=system(sprintf('"%s" "%s"',cfg.bellhop_exe,name));
clear cleanup
if status~=0, error('Bellhop failed for %s: %s',name,command_output); end
run_type=upper(char(cfg.run_type));
if run_type(1)=='A'
    data_file=[cfg.case_root '.arr']; data=read_bellhop_arrivals_ascii_vertical(data_file);
elseif ismember(run_type(1),['C','I'])
    data_file=[cfg.case_root '.shd']; data=read_bellhop_shd_unfolded_vertical(data_file);
else
    data_file=[cfg.case_root '.ray'];
end
if exist(data_file,'file')~=2 || exist([cfg.case_root '.prt'],'file')~=2
    error('Bellhop did not create expected output for %s.',name);
end
result=struct('config',cfg,'data',data,'command_output',command_output, ...
    'files',struct('env',env_file,'sbp',sbp_file,'data',data_file,'prt',[cfg.case_root '.prt']));
end
