function receipt=prepare_bellhop_1500m_p1(project_root,output_dir)
%PREPARE_BELLHOP_1500M_P1 Freeze exactly 18 U0 slots; no acoustic solve.
assert(isfolder(project_root)&&isfolder(output_dir));
addpath(fullfile(project_root,'scripts','validation','support'));
x=[0 10 50 100 200]; tags={'D','S','B','BS'}; v=[1484 1496 1504 1516];
opt='/home/misc/cwx/bellhop_rotated_C_transfer_20261116_01/Explain/results/validation/bellhop_rotated_fig12b/20261116_linux_C_01/bin/bellhop_rotated_receiver64_2020';
base='/home/misc/cwx/bellhop_1500m_P1_20261118_01/bin/bellhop_REAL8_unoptimized';
osh='5d87293a27e7fafc3aedefb569f3f35065294f42426503c4881d6eb55e6769e9';
bsh='61158fd0563332dd2a3ba81d29aac84267e1e469ad09a2918453718ceb2b1827';
remote='/home/misc/cwx/bellhop_1500m_P1_20261118_01';
slots=struct([]); k=0;
for j=1:4
 for f=[4000 4062.5], add(j,f,'C',false); end
end
for mode={'R','A'}, for j=1:4, add(j,4000,mode{1},false); end; end
add(2,4000,'C',true); add(4,4000,'C',true);
assert(k==18); save(fullfile(output_dir,'plan.mat'),'slots','x','v');
writetable(struct2table(rmfield(slots,{'cfg','meta'})),fullfile(output_dir,'task_manifest_18.csv'));
receipt=struct('status','PREPARED_18_NOT_EXECUTED','slots',18,'environment','U0','source_width','N/A ideal omni 3D point; 1m reference','frequencies_hz',[4000 4062.5],'receiver_x_m',x,'water_depth_m',1500,'source_z_m',1490,'receiver_z_m',6,'water_c_mps',1500,'water_rho_gcm3',1,'bottom_cp_mps',1600,'bottom_rho_gcm3',1.8,'bottom_atten_db_wavelength',0,'total_budget_s',300,'unoptimized_source_Step_sha256','cdde8e9bd81b772d7a55d46628a0233d33d96e3e23cab7279adae8d566858cf3','optimized_sha256',osh,'baseline_sha256',bsh,'baseline_cases',{{'S4000','BS4000'}},'later_stages_authorized',false);
writejson(fullfile(output_dir,'checkpoint_prepared.json'),receipt);
 function add(j,f,mode,isbase)
  k=k+1; cfg=struct('frequency_hz',f,'source_z_m',1490,'water_depth_m',1500,'munk_parameters',struct(),'bottom_cp_mps',1600,'bottom_rho_gcm3',1.8,'bottom_atten_db_wavelength',0,'direction','up','transverse_origin_m',0,'artificial_halfwidth_m',1000,'ssp_spacing_m',.625,'step_m',.5,'beam_count',201,'angle_limits_deg',[-1 1]*1.1*atan2d(200,v(j)),'receiver_x_m',x,'run_type',mode,'source_geometry','R','uniform_control',true,'uniform_c_mps',1500,'receiver_range_m',v(j),'range_limit_m',v(j),'nominal_receiver_x_m',x,'nominal_receiver_depths_m',6,'receiver_depths_m',6,'receiver_depth_order',1,'receiver_x_order',1:5);
  % Legacy writer evaluates its named SSP before applying the uniform control.
  cfg.munk_parameters=struct('c0_mps',1500,'axis_depth_m',1300,'scale_m',1300,'epsilon',0);
  if j>1
   if j==2, cfg.reflection_wall='surface'; cfg.wall_range_m=1490; order=1; ns=1; nb=0;
   else, cfg.direction='down'; cfg.reflection_wall='bottom'; cfg.wall_range_m=10; order=j-2; ns=j-3; nb=1; end
   cfg.reflection_order=order; cfg.expected_surface_bounces=ns; cfg.expected_bottom_bounces=nb; cfg.last_wall_range_m=cfg.wall_range_m+(order-1)*1500;
   coeff=-(1-2*strcmp(cfg.direction,'up'))*(-1)^order;
   [~,cfg.receiver_x_order]=sort(coeff*x); cfg.receiver_x_m=x(cfg.receiver_x_order);
  end
  name=sprintf('%02d_%s_%g_%s',k,tags{j},f,mode); name=strrep(name,'.','p'); if isbase, name=[name '_baseline']; end
  dest=fullfile(output_dir,'jobs',name); assert(~isfolder(dest)); mkdir(dest); root=fullfile(dest,'case');
  if j==1, meta=write_bellhop_rotated_munk_q(root,cfg); fid=fopen([root '.iwcov'],'w'); assert(fid>=0); fprintf(fid,'2\n0\n0\n0\n0\n'); fclose(fid); write_bellhop_rotated_field_options(root,cfg.receiver_range_m,meta.q_receiver_m,point3d=true,fast_q=true,native_range_metadata_repair=true);
  else, meta=write_bellhop_rotated_multifield_q(root,cfg); end
  % Give U0 inputs an unambiguous title without changing any numeric line.
  t=fileread([root '.env']); lines=splitlines(string(t)); lines(1)="'U0 1500m P1 "+tags{j}+" "+mode+"; source width N/A'"; fid=fopen([root '.env'],'w'); fprintf(fid,'%s',join(lines,newline)); fclose(fid);
  exe=opt; sha=osh; if isbase, exe=base; sha=bsh; end
  item=struct('slot',k,'path',tags{j},'frequency_hz',f,'mode',mode,'beams',201,'angle_min_deg',cfg.angle_limits_deg(1),'angle_max_deg',cfg.angle_limits_deg(2),'binary_sha256',sha,'binary',exe,'output_dir',[remote '/jobs/' name],'local_dir',dest,'name',name,'cfg',cfg,'meta',meta,'is_baseline',isbase,'status','NOT_EXECUTED');
  if isempty(slots), slots=item; else, slots(k)=item; end
 end
end
function writejson(path,val)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); c=onCleanup(@()fclose(fid)); fprintf(fid,'%s',jsonencode(val,PrettyPrint=true));
end