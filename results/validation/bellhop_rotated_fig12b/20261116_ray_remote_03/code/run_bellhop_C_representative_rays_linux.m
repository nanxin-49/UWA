function receipt=run_bellhop_C_representative_rays_linux(project_root,output_dir,remaining_s)
% R-mode geometry only, original C parameters except3 selected launch angles.
assert(~isfolder(output_dir)); mkdir(output_dir); start=tic; diary(fullfile(output_dir,'diary.log')); guard=onCleanup(@()diary('off'));
family=fullfile(project_root,'results','validation','bellhop_rotated_fig12b'); p=load(fullfile(family,'20261114_full_suite_plan_03','C_Q625_h05_B401.mat'),'plan'); exe=getenv('BH_RAY_EXE'); expected_sha=getenv('BH_RAY_EXE_SHA'); assert(isfile(exe)&&strcmpi(black_file_hash(exe),expected_sha));
tags={'direct_up_part3','up_order1','down_order1','up_order3','down_order4','up_order13'}; groups=struct('tag',{},'cfg',{},'rays',{},'process',{}); max_map=0; launched=0; failure='';
try
for k=1:numel(tags)
 [~,ram]=bellhop_linux_memory(); assert(ram.PhysicalMemory.Available>=1.5*2^30,'Ray-only RAM launch gate'); assert(remaining_s-toc(start)>10,'Ray-only shared budget gate'); j=find(strcmp({p.plan.jobs.tag},tags{k})); assert(isscalar(j)); cfg=p.plan.jobs(j).cfg; cfg.run_type='R'; cfg.beam_count=3;
 d=fullfile(output_dir,tags{k}); mkdir(d); root=fullfile(d,'case'); previous=fullfile(family,'20261116_ray_remote_02','rays','direct_up_part3');
 if strcmp(tags{k},'direct_up_part3') && isfile(fullfile(previous,'rays.mat'))
  cached=load(fullfile(previous,'rays.mat'),'cfg','process'); assert(isequaln(cached.cfg,cfg)); copyfile(fullfile(previous,'*'),d); process=cached.process;
 else
  write_bellhop_rotated_REAL8_C_inputs(root,cfg); launched=launched+1; process=run_bellhop_rotated_local_case(exe,root,min(120,remaining_s-toc(start)-10),allow_task_termination=true,minimum_runtime_available_bytes=256*2^20);
 end
 rays=read_rays([root '.ray']); assert(numel(rays)==3); order=0; if isfield(cfg,'reflection_order'), order=cfg.reflection_order; end
 for n=1:3
  assert(rays(n).top+rays(n).bottom==order); if order>0, assert(rays(n).top==cfg.expected_surface_bounces&&rays(n).bottom==cfg.expected_bottom_bounces); end
  r=rays(n).rq(:,1); q=rays(n).rq(:,2); sg=1; if strcmp(cfg.direction,'up'), sg=-1; end
  stage=zeros(size(r));
  if order>0
   for bounce=1:order
    wall=cfg.wall_range_m+(bounce-1)*cfg.water_depth_m; hit=find(r==wall);
    pairs=hit(ismember(hit+1,hit)); assert(isscalar(pairs),'Expected exactly one incident/reflected native node pair at each wall');
    incident=pairs; reflected=incident+1; assert(abs(q(incident)+q(reflected)-2*cfg.transverse_origin_m)<1e-8,'Pi chart q mapping mismatch'); stage(reflected:end)=bounce;
   end
  end
  s=sg*(-1).^stage; r0=zeros(size(r)); zs=cfg.source_z_m+zeros(size(r)); mask=stage>0;
  if order>0, r0(mask)=cfg.wall_range_m+(stage(mask)-1)*cfg.water_depth_m; zs(mask)=cfg.water_depth_m*(1-s(mask))/2; end
  physical_z=zs+s.*(r-r0); physical_x=-s.*(q-cfg.transverse_origin_m); err=max(abs(r0+s.*(physical_z-zs)-r))+max(abs(cfg.transverse_origin_m-s.*physical_x-q)); max_map=max(max_map,err); assert(err<1e-8,'Ray coordinate map gate');
  keep=r<=cfg.range_limit_m & physical_z>=-1e-8 & physical_z<=cfg.water_depth_m+1e-8; assert(nnz(keep)>2); rays(n).physical_x=physical_x; rays(n).physical_z=physical_z; rays(n).stage=stage; rays(n).physical_plot_mask=keep; rays(n).map_roundtrip_error_m=err;
 end
 ev=[root '.rr64eval']; if isfile(ev), info=dir(ev); assert(info.bytes==8,'R mode unexpectedly evaluated receiver pressure'); end
 groups(k)=struct('tag',tags{k},'cfg',cfg,'rays',rays,'process',process); save(fullfile(d,'rays.mat'),'cfg','rays','process','-v7.3'); fprintf('RAY_ONLY %s native%.3fs rays3 reflections%d\n',tags{k},process.elapsed_s,order);
end
catch e
 failure=getReport(e,'extended','hyperlinks','off'); fprintf('%s\n',failure);
end
save(fullfile(output_dir,'all_rays.mat'),'groups','-v7.3'); status='RAY_ONLY_COMPLETE'; if ~isempty(failure), status='RAY_ONLY_STOPPED_NO_RETRY'; end
receipt=struct('status',status,'groups_complete',numel(groups),'native_jobs_launched',launched,'reused_ray_jobs',1,'receiver_pressure_evaluations',0,'ABC_recomputed',0,'new_FFP',0,'max_coordinate_roundtrip_m',max_map,'elapsed_s',toc(start),'remaining_budget_s',remaining_s-toc(start),'source','Original C frozen configs; only run_type R and beam_count3; no field integration','error',failure);
fid=fopen(fullfile(output_dir,'checkpoint.json'),'w','n','UTF-8'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid); disp(jsonencode(receipt));
end
function rays=read_rays(file)
fid=fopen(file,'r'); assert(fid>=0); close_file=onCleanup(@()fclose(fid)); fgetl(fid); fscanf(fid,'%f',1); fscanf(fid,'%d',3); beams=fscanf(fid,'%d',2); fscanf(fid,'%f',2); fgetl(fid); kind=strtrim(fgetl(fid)); assert(contains(kind,'rz'));
rays=struct('angle_deg',{},'top',{},'bottom',{},'rq',{});
for k=1:beams(1), a=fscanf(fid,'%f',1); info=fscanf(fid,'%d',3); assert(numel(info)==3); rq=fscanf(fid,'%f',[2 info(1)]).'; assert(size(rq,1)==info(1)); rays(k)=struct('angle_deg',a,'top',info(2),'bottom',info(3),'rq',rq); end
end
