function receipt=run_bellhop_rays_remote(project_root,run_dir)
% Geometry only, one serial session. No receiver-pressure evaluation.
addpath(fullfile(project_root,'scripts','validation','support')); addpath(fullfile(project_root,'results','validation','bellhop_rotated_fig12b','20261116_linux_C_01','runtime03'),'-begin'); addpath(fullfile(run_dir,'code'),'-begin');
initial=str2double(getenv('BH_RAY_REMAINING_S')); total_limit=str2double(getenv('BH_RAY_TOTAL_LIMIT_S')); assert(isfinite(initial)&&initial>0&&isfinite(total_limit)); epoch=str2double(getenv('BH_RAY_START_EPOCH')); budget=@()initial-(posixtime(datetime('now','TimeZone','UTC'))-epoch);
receipt=struct('status','RAY_ONLY_RUNNING','initial_remaining_s',initial,'native_jobs',0,'groups_complete',0,'new_FFP',0,'ABC_recomputed',0); diary(fullfile(run_dir,'driver_diary.log')); close_diary=onCleanup(@()diary('off'));
try
 assert(license('test','MATLAB')); qual=jsondecode(fileread(fullfile(project_root,'results','validation','bellhop_rotated_fig12b','20261116_linux_C_01','checkpoint_qualification03.json'))); assert(qual.failed_checks==0&&strcmp(qual.status,'QUALIFICATION_PASS_C9_ALLOWED'));
 assert(strcmpi(black_file_hash(fullfile(project_root,'results','validation','bellhop_rotated_fig12b','20261114_full_suite_plan_03','C_Q625_h05_B401.mat')),'cce53c33f09104675eed7ef6131b6b9d6154259b25baaeffa83570d331947e97'));
 assert(budget()>30,'Shared remaining budget gate before R-only');
 rays=run_bellhop_C_representative_rays_linux(project_root,fullfile(run_dir,'rays'),budget()-10); receipt.rays=rays; receipt.native_jobs=rays.native_jobs_launched; receipt.reused_ray_jobs=rays.reused_ray_jobs; receipt.groups_complete=rays.groups_complete;
 assert(strcmp(rays.status,'RAY_ONLY_COMPLETE')&&rays.groups_complete==6,'Ray geometry stop: inspect original case; no retry');
 assert(budget()>12,'Shared budget before graphics'); s=load(fullfile(run_dir,'rays','all_rays.mat'),'groups'); groups=s.groups; figs=fullfile(run_dir,'figures'); mkdir(figs); map_checks=0; no_pressure_checks=0; total_nodes=0; native_s=0;
for k=1:numel(groups)
 g=groups(k); native_s=native_s+g.process.elapsed_s; assert(g.process.exit_code==0&&isempty(g.process.abort_reason)); cfg=g.cfg;
 for j=1:3
  a=g.rays(j); r=a.rq(:,1); q=a.rq(:,2); stage=a.stage; sg=1; if strcmp(cfg.direction,'up'), sg=-1; end
  order=0; if isfield(cfg,'reflection_order'), order=cfg.reflection_order; end
  assert(a.top+a.bottom==order); assert(all(diff(r)>=0));
  direction=sg*(-1).^stage; rr0=zeros(size(r)); zz0=cfg.source_z_m+zeros(size(r)); mask=stage>0;
  if order>0, rr0(mask)=cfg.wall_range_m+(stage(mask)-1)*cfg.water_depth_m; zz0(mask)=cfg.water_depth_m*(1-direction(mask))/2; end
  if order>0; for bounce=1:order; idx=find(diff(stage)==1 & stage(1:end-1)==bounce-1); assert(isscalar(idx)); assert(abs(a.physical_x(idx)-a.physical_x(idx+1))<1e-8 && abs(a.physical_z(idx)-a.physical_z(idx+1))<1e-8,'Physical reflection seam discontinuity'); end; end; rback=rr0+direction.*(a.physical_z-zz0); qback=cfg.transverse_origin_m-direction.*a.physical_x;
  assert(max(abs(rback-r))+max(abs(qback-q))<1e-8); map_checks=map_checks+1; total_nodes=total_nodes+size(a.rq,1);
 end
 ev=dir(fullfile(run_dir,'rays',g.tag,'case.rr64eval')); assert(isempty(ev)||(isscalar(ev)&&ev.bytes==8)); no_pressure_checks=no_pressure_checks+1;
end
checks=table(strings(0,1),false(0,1),'VariableNames',{'check','pass'});
for view=1:2
 assert(budget()>5,'Shared budget during graphics'); fig=figure('Visible','off','Color','w','Position',[20 20 1700 1100]); tiles=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
 for k=1:6
  g=groups(k); ax=nexttile(tiles); hold(ax,'on'); handles=gobjects(1,3); labels=cell(1,3);
  for j=1:3
   a=g.rays(j); keep=a.physical_plot_mask; labels{j}=sprintf('alpha=%.5g deg',a.angle_deg); style={'--','-',':'}; colors={[0 .35 .8],[.1 .1 .1],[.85 .3 0]};
   if view==1
    xx=a.rq(:,1)/1000; yy=a.rq(:,2)-g.cfg.transverse_origin_m; stages=unique(a.stage(keep));
    for st=stages(:).', use=keep&a.stage==st; h=plot(ax,xx(use),yy(use),'LineStyle',style{j},'Color',colors{j},'LineWidth',.8); h.UserData=struct('group',k,'ray',j,'stage',st,'view',view); if ~isgraphics(handles(j)), handles(j)=h; end; end
   else
    h=plot(ax,a.physical_x(keep),a.physical_z(keep),'LineStyle',style{j},'Color',colors{j},'LineWidth',.8); h.UserData=struct('group',k,'ray',j,'stage',-1,'view',view); handles(j)=h;
   end
  end
  set(ax,'YDir','reverse','FontSize',10); grid(ax,'on'); title(ax,strrep(g.tag,'_',' '));
  if view==1, xlabel(ax,'Unfolded range r (km)'); ylabel(ax,'Chart q-q0 (m)'); else, xlabel(ax,'Physical x (m)'); ylabel(ax,'Physical depth z (m)'); end
  legend(ax,handles,labels,'Location','best','Visible','on','FontSize',8);
 end
 if view==1, stem='R01_rotated_rays'; title(tiles,{'Actual Bellhop R-mode nodes | frozen C Q0.625m and original steps','Chart segments are split at pi-rotation seams; no artificial connector'},'FontSize',15); else, stem='R02_physical_rays'; title(tiles,{'Actual rays mapped back to physical x-depth | source4950m, H5000m','3 representative angles per branch, not a pressure/intensity plot'},'FontSize',15); end
 drawnow; savefig(fig,fullfile(figs,[stem '.fig'])); exportgraphics(fig,fullfile(figs,[stem '.png']),'Resolution',130); close(fig);
 fig=openfig(fullfile(figs,[stem '.fig']),'new','invisible'); lines=findall(fig,'Type','line');
 for j=1:numel(lines)
  u=lines(j).UserData; if ~isstruct(u)||~isfield(u,'group'), continue; end; g=groups(u.group); a=g.rays(u.ray); keep=a.physical_plot_mask;
  if view==1, use=keep&a.stage==u.stage; expected_x=a.rq(use,1)/1000; expected_y=a.rq(use,2)-g.cfg.transverse_origin_m; else, expected_x=a.physical_x(keep); expected_y=a.physical_z(keep); end
  checks(end+1,:)={sprintf('view%d_line%d',view,j),isequal(lines(j).XData(:),expected_x)&&isequal(lines(j).YData(:),expected_y)}; %#ok<AGROW>
 end
 close(fig);
end
assert(all(checks.pass)); writetable(checks,fullfile(run_dir,'figure_data_checks.csv'));
html=['<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>Bellhop C代表声线</title><style>body{max-width:1700px;margin:24px auto;padding:0 20px;font:16px/1.65 sans-serif}img{width:100%;height:auto}</style><h1>Bellhop C工况：真实代表声线补齐</h1><p>6类各3束：直达、首海面、首海底及3/4/13阶反射。仅R模式输出，沿用冻结C的介质、原步长、Q0.625m、REAL8和平面优化；没有接收复场求值，也不重算A/B/C。曲线为原生保存节点，不插值或平滑。</p><h2>旋转坐标视图</h2><img src="figures/R01_rotated_rays.png" alt="6类真实声线旋转坐标"><p>每次pi旋转跨chart的接缝单独分段，未用直线连接坐标跳变。<a href="figures/R01_rotated_rays.fig">可编辑FIG</a></p><h2>物理x-depth视图</h2><img src="figures/R02_physical_rays.png" alt="真实声线映射回物理坐标"><p>z向下，分支内使用正确方向和偏移映射。只对超出目标水体的末端guard节点限定显示，原节点完整保存。<a href="figures/R02_physical_rays.fig">可编辑FIG</a></p><p>轴线上垂直射线是alpha=0的几何路线，不代表声场中央亮带。少量代表声线不等同全部声束或独立束数收敛。</p></html>']; fid=fopen(fullfile(run_dir,'index.html'),'w','n','UTF-8'); fprintf(fid,'%s',html); fclose(fid);
receipt.status='RAY_ONLY_VISUALIZATION_COMPLETE'; receipt.native_elapsed_s=native_s; receipt.total_native_nodes=total_nodes; receipt.map_checks=map_checks; receipt.max_map_error_m=rays.max_coordinate_roundtrip_m; receipt.no_receiver_pressure_checks=no_pressure_checks; receipt.figure_data_checks=height(checks); receipt.failed_checks=nnz(~checks.pass); receipt.full_representative_trajectories=true;
catch e
 receipt.status='RAY_ONLY_STOPPED_NO_RETRY'; receipt.error=getReport(e,'extended','hyperlinks','off'); fprintf('%s\n',receipt.error);
end
receipt.matlab_version=version; receipt.job_wall_s=initial-budget(); receipt.remaining_original_budget_s=budget(); receipt.accounted_cumulative_s=total_limit-receipt.remaining_original_budget_s; receipt.approved_total_limit_s=total_limit; receipt.next_allowed_step='Manual review only; no D/newFFP or field recalculation.'; write(fullfile(run_dir,'checkpoint_stage.json'),receipt); fprintf('RAY_REMOTE_FINAL %s groups%d remaining%.3fs\n',receipt.status,receipt.groups_complete,receipt.remaining_original_budget_s);
end
function write(path,value)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); c=onCleanup(@()fclose(fid)); fprintf(fid,'%s',jsonencode(value,PrettyPrint=true));
end
