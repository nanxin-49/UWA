function probe=probe_deepsea3d_fg_geometry()
%PROBE_DEEPSEA3D_FG_GEOMETRY Saved field only; no new propagation solve.
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
out=fullfile(root,'results','validation','deepsea3d_tfe_dealiasing');
file=fullfile(out,'pm_fg_geometry_probe.mat');assert(~isfile(file),'FG:Exists','Keep saved probe.');
a=load(fullfile(out,'pm_fg.mat'),'output','snapshots','cfg');
p=deepsea3d_fg_plan(a.cfg,a.output.options);W=a.snapshots(end).W;
probe3=observe_deepsea3d_fg_residual(W,p,a.cfg,3);
probe4=observe_deepsea3d_fg_residual(W,p,a.cfg,4);
probe=struct('N',32,'three',probe3,'four',probe4,'new_solves',0, ...
    'wall_s',probe3.wall_s+probe4.wall_s);
save(file,'probe','-v7');disp(probe4);
end
