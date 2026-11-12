task='/home/misc/cwx/piston_p0p1_20261112_01';
root=fullfile(task,'results','comparisons','piston_2d_single_frequency','20261112_remote_p0p1_01');
q=load(fullfile(root,'common.mat'));
rows=[];
for level=1:6
 v=load(fullfile(root,'jobs',sprintf('flat_bie_l%d',level),'data.mat'),'result');
 r=v.result;
 er=vecnorm(r.reflect-q.p_mirror)./vecnorm(q.p_mirror);
 rows=[rows;level r.boundary.node_count er max(r.boundary.offgrid_boundary_residual)];
end
T=array2table(rows,'VariableNames',{'level','nodes','W20_mirror_relative','W10_mirror_relative','boundary_residual'});
writetable(T,fullfile(task,'processes','diagnose_flat_partial_01','mirror_metrics.csv'));
disp(T);
