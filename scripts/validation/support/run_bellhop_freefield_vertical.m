function result = run_bellhop_freefield_vertical(cfg)
%RUN_BELLHOP_FREEFIELD_VERTICAL Run one matched-halfspace free-space case.
arguments
    cfg (1,1) struct
end
required={'bellhop_exe','case_root','run_type'};
for ii=1:numel(required)
    if ~isfield(cfg,required{ii}), error('Missing cfg.%s.',required{ii}); end
end
if exist(cfg.bellhop_exe,'file')~=2, error('Bellhop executable is missing.'); end
out_dir=fileparts(cfg.case_root); if ~exist(out_dir,'dir'), mkdir(out_dir); end
env_file=write_bellhop_freefield_env_vertical(cfg.case_root,cfg);
extensions={'.arr','.shd','.ray','.prt'};
for ii=1:numel(extensions)
    target=[cfg.case_root extensions{ii}]; if exist(target,'file')==2, delete(target); end
end
old=pwd; cleanup=onCleanup(@()cd(old)); cd(out_dir);
[~,name]=fileparts(cfg.case_root);
[status,command_output]=system(sprintf('"%s" "%s"',cfg.bellhop_exe,name));
clear cleanup
if status~=0, error('Bellhop failed for %s: %s',name,command_output); end
run_type=upper(char(cfg.run_type));
if run_type=='A'
    data_file=[cfg.case_root '.arr']; data=read_bellhop_arrivals_ascii_vertical(data_file);
elseif ismember(run_type,['C','I'])
    data_file=[cfg.case_root '.shd']; data=read_bellhop_shd_2d_vertical(data_file);
else
    data_file=[cfg.case_root '.ray']; data=[];
end
if exist(data_file,'file')~=2 || exist([cfg.case_root '.prt'],'file')~=2
    error('Bellhop did not create the expected output for %s.',name);
end
result=struct('config',cfg,'data',data,'command_output',command_output, ...
    'files',struct('env',env_file,'data',data_file,'prt',[cfg.case_root '.prt']));
end
