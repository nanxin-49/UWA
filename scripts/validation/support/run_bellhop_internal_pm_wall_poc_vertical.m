function result = run_bellhop_internal_pm_wall_poc_vertical(cfg)
%RUN_BELLHOP_INTERNAL_PM_WALL_POC_VERTICAL Run the PM internal-wall binary.
arguments
    cfg (1,1) struct
end
required={'bellhop_exe','case_root','run_type','wall_r0_m','wall_seed', ...
    'wall_profile_r_m','wall_profile_z_m','mapped_receiver_range_m'};
for ii=1:numel(required)
    if ~isfield(cfg,required{ii}), error('Missing cfg.%s.',required{ii}); end
end
if numel(cfg.wall_profile_r_m)~=numel(cfg.wall_profile_z_m) || numel(cfg.wall_profile_r_m)<3 || ...
        any(~isfinite(cfg.wall_profile_r_m(:))) || any(~isfinite(cfg.wall_profile_z_m(:)))
    error('PM wall profile must contain at least three finite paired samples.');
end
if any(diff(cfg.wall_profile_z_m(:))<=0)
    error('PM wall profile parameter samples must increase.');
end
if exist(cfg.bellhop_exe,'file')~=2
    error('Bellhop PM internal-wall validation executable is missing: %s',cfg.bellhop_exe);
end
reuse_existing = isfield(cfg,'reuse_existing') && logical(cfg.reuse_existing);
if reuse_existing && exist([cfg.case_root '.shd'],'file')==2 && ...
        exist([cfg.case_root '.iwdiag'],'file')==2 && ...
        dir([cfg.case_root '.shd']).bytes>0 && dir([cfg.case_root '.iwdiag']).bytes>0
    env_file=[cfg.case_root '.env']; sbp_file=[cfg.case_root '.sbp'];
else
    [env_file,sbp_file]=write_bellhop_unfolded_gaussian_env_vertical(cfg.case_root,cfg);
    reuse_existing = false;
end
iwpm_file=[cfg.case_root '.iwpm'];
if ~reuse_existing
    fid=fopen(iwpm_file,'w');
    if fid<0, error('Cannot create %s.',iwpm_file); end
    cleanup=onCleanup(@()fclose(fid));
    fprintf(fid,'%.17g\n%.17g\n%d\n%d\n',cfg.wall_r0_m,cfg.mapped_receiver_range_m, ...
        cfg.wall_seed,numel(cfg.wall_profile_r_m));
    fprintf(fid,'%.17g %.17g\n',[cfg.wall_profile_r_m(:) cfg.wall_profile_z_m(:)].');
    clear cleanup
end

extensions={'.arr','.shd','.ray','.prt','.iwdiag','.covdiag'};
if ~reuse_existing
    for ii=1:numel(extensions)
        target=[cfg.case_root extensions{ii}];
        if exist(target,'file')==2, delete(target); end
    end
end
out_dir=fileparts(cfg.case_root); old=pwd; cleanup=onCleanup(@()cd(old)); cd(out_dir);
[~,name]=fileparts(cfg.case_root);
data_file=[cfg.case_root '.shd']; diag_file=[cfg.case_root '.iwdiag']; prt_file=[cfg.case_root '.prt'];
[status,command_output]=deal(0,'reused existing raw Bellhop outputs');
if ~reuse_existing
    [status,command_output]=system(sprintf('"%s" "%s"',cfg.bellhop_exe,name));
    clear cleanup
    if status~=0, error('Internal PM-wall Bellhop failed for %s: %s',name,command_output); end
end
for path={data_file,diag_file,prt_file}
    if exist(path{1},'file')~=2, error('Bellhop did not create %s.',path{1}); end
end
data=read_bellhop_shd_unfolded_vertical(data_file);
diag_matrix=readmatrix(diag_file,'FileType','text','CommentStyle','#');
diag_matrix=diag_matrix(~all(isnan(diag_matrix),2),:);
names={'alpha_deg','hit_r','hit_z','wall_residual','wall_t_r','wall_t_z','wall_n_r','wall_n_z', ...
    'tangent_error','normal_error','inc_ur','inc_uz','ref_ur','ref_uz','rot_ur','rot_uz', ...
    'specular_error','rotation_error','phase_in','phase_ref','phase_delta','amp_in','amp_ref','amp_delta', ...
    'p1_in','p2_in','p1_ref','p2_ref','p_ref_error','q1_in','q2_in','q1_ref','q2_ref', ...
    'q_ref_error','p_rot_error','q_rot_error','tau_wall_real','tau_wall_imag', ...
    'tau_receiver_real','tau_receiver_imag','min_post_dr','n_post','kappa', ...
    'wall_seg','wall_lambda','wall_tg','wall_th','wall_rm','wall_rn'};
if size(diag_matrix,2)~=numel(names)
    error('Expected %d internal PM-wall diagnostic columns, found %d.',numel(names),size(diag_matrix,2));
end
diagnostics=array2table(diag_matrix,'VariableNames',names);
result=struct('config',cfg,'data',data,'diagnostics',diagnostics,'command_output',command_output, ...
    'files',struct('env',env_file,'sbp',sbp_file,'iwpm',iwpm_file,'data',data_file, ...
    'diagnostics',diag_file,'prt',prt_file));
end
