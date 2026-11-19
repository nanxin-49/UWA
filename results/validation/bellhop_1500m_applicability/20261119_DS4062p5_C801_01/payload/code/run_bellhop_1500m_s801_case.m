function receipt=run_bellhop_1500m_s801_case(project_root,run_dir,remaining_s)
%RERUN_BELLHOP_1500M_D4000_C Exactly one frozen direct slot; no stage advancement.
start=tic; diary(fullfile(run_dir,'C_diary.log')); closer=onCleanup(@()diary('off'));
addpath(fullfile(project_root,'scripts','validation','support'));
addpath(fullfile(fileparts(run_dir),'code'),'-begin'); addpath(fullfile(project_root,'results','validation','bellhop_rotated_fig12b','20261116_linux_C_01','runtime03'),'-end');
assert(~isfile(fullfile(run_dir,'checkpoint_C.json')),'Never reuse an executed run directory');
remaining_s=remaining_s-(posixtime(datetime('now','TimeZone','UTC'))-str2double(getenv('P1_START_EPOCH')));
assert(isfinite(remaining_s)&&remaining_s>30,'Whole-job startup budget');
s=load(fullfile(run_dir,'plan.mat')); slots=s.slots; x=[0 10 50 100 200]; v=[1484 1496]; pressure=complex(nan(8,5)); results=struct([]);
checks=table(strings(0,1),false(0,1),zeros(0,1),zeros(0,1),'VariableNames',{'check','pass','value','limit'});
status='STOPPED'; err=''; completed=0; started=0;
try
 for k=1:1
  job=slots(k); cfg=job.cfg; root=fullfile(job.output_dir,'case');
  assert(numel(slots)==1&&strcmp(job.path,'S')&&cfg.frequency_hz==4062.5&&strcmp(cfg.run_type,'C')&&cfg.beam_count==801);
  assert(isequal(cfg.angle_limits_deg,[-1 1]*1.1*atan2d(200,1496)));
  assert(isfile([root '.rotq'])&&isequal(sscanf(fileread([root '.iwcov']),'%g',1),5),'Required frozen S sidecars');
  assert(toc(start)<remaining_s-15,'TOTAL_BUDGET_GATE');
  assert(strcmpi(black_file_hash(job.binary),job.binary_sha256),'BINARY_IDENTITY_GATE');
  [a,b]=system(['df -B1 --output=avail ' run_dir ' | tail -1']); assert(a==0&&str2double(b)>=20*2^30,'DISK_GATE');
  started=started+1;
  proc=run_bellhop_rotated_local_case(job.binary,root,min(120,remaining_s-toc(start)-15),allow_task_termination=true,minimum_runtime_available_bytes=256*2^20);
  d=readmatrix([root '.r25diag'],'FileType','text','CommentStyle','#');
  add(sprintf('slot%d_all801rays',k),abs(size(d,1)-801),0);
  add(sprintf('slot%d_direction_norm',k),max(d(:,8)),1e-8);
  ns=1; nb=0;
  add(sprintf('slot%d_bounce_counts',k),nnz(d(:,9)~=ns|d(:,10)~=nb),0);
  raw=read_bellhop_rotated_field(root); p=complex(zeros(1,5)); p(cfg.receiver_x_order)=(-conj(raw.pressure)).'; pressure(k,:)=p;
  add(sprintf('slot%d_finite',k),nnz(~isfinite(p)),0);
  fid=fopen([root '.rr64input'],'rb','ieee-le'); assert(fid>=0); h=fread(fid,4,'int32'); env=fread(fid,1,'double'); rr=fread(fid,1,'double'); qq=fread(fid,5,'double'); fclose(fid);
  add(sprintf('slot%d_rr64header',k),double(~isequal(h,[1;5;8;8])),0);
  add(sprintf('slot%d_exact_range_bits',k),bitdiff(rr,cfg.receiver_range_m),0);
  add(sprintf('slot%d_exact_q_bits',k),bitdiff(qq,job.meta.q_receiver_m),0);
  add(sprintf('slot%d_ENV_roundoff_ulp',k),abs(env-rr)/eps(rr),2);
  fid=fopen([root '.rr64eval'],'rb','ieee-le'); assert(fid>=0); eh=fread(fid,2,'int32'); bytes=fread(fid,inf,'*uint8'); fclose(fid); assert(isequal(eh,[1;5])&&mod(numel(bytes),24)==0);
  rec=reshape(bytes,24,[]); ir=double(typecast(reshape(rec(1:4,:),[],1),'int32')); iq=double(typecast(reshape(rec(5:8,:),[],1),'int32')); ar=typecast(reshape(rec(9:16,:),[],1),'double'); aq=typecast(reshape(rec(17:24,:),[],1),'double');
  add(sprintf('slot%d_actual_eval_bits',k),bitdiff(ar,repmat(rr,numel(ar),1))+bitdiff(aq,qq(iq)),0);
  add(sprintf('slot%d_actual_eval_coverage',k),numel(setdiff((1:5)',unique(iq)))+nnz(ir~=1),0);
  kind=2; L=hypot(x,1496);
  add('native_straight_ray_path_length_m',max(abs(hypot(d(:,2),d(:,3))-1500*d(:,4))),1e-4);
  add('native_straight_ray_time_s',max(abs(hypot(d(:,2),d(:,3))/1500-d(:,4))),1e-7);
  add('native_straight_ray_direction_rad',max(abs(atan2(-d(:,3),d(:,2))-deg2rad(d(:,1)))),1e-6);
  phys=[x;6*ones(size(x))]; chart=[1490+phys(2,:);-phys(1,:)]; restored=[-chart(2,:);chart(1,:)-1490];
  add('proper_rotation_roundtrip_m',max(abs(restored-phys),[],'all'),1e-8);
  add('original_order_chart_coordinates',max(abs(chart-[1496*ones(size(x));-x]),[],'all'),0);
  fid=fopen([root '.rr64scale'],'rb','ieee-le'); sh=fread(fid,1,'int32'); ib=fread(fid,1,'int32'); sr=fread(fid,1,'double'); fclose(fid);
  add('scale_exact_range_bits',bitdiff(sr,1496)+abs(sh-1)+abs(ib-1),0);
  ref=exp(1i*2*pi*cfg.frequency_hz/1500*L)./L;
  if kind==2, ref=-ref; end
  metrics=struct();
  if kind<=2
   E=abs(ref); support=abs(p)>=.01*E & abs(ref)>=.01*E;
   metrics=struct('complex_L2',norm(p-ref)/norm(ref),'max_complex_difference_over_E',max(abs(p-ref)./E),'max_phase_rad',max(abs(angle(p(support).*conj(ref(support))))),'phase_support',nnz(support));
   add(sprintf('slot%d_analytic_complex_L2',k),metrics.complex_L2,1e-3);
   add(sprintf('slot%d_analytic_max_complex_over_E',k),metrics.max_complex_difference_over_E,2e-3);
   add(sprintf('slot%d_analytic_phase_rad',k),metrics.max_phase_rad,.01);
  end
  item=struct('slot',k,'path',job.path,'frequency_hz',job.frequency_hz,'process',proc,'metrics',metrics); if isempty(results),results=item;else,results(end+1)=item;end
  save(fullfile(job.output_dir,'component.mat'),'cfg','raw','p','ref','L','proc','metrics');
  writetable(table(x',real(p)',imag(p)',real(ref)',imag(ref)',L',L'/1500,'VariableNames',{'x_m','p_real','p_imag','reference_real','reference_imag','analytic_length_m','analytic_delay_s'}),fullfile(job.output_dir,'field.csv'));
  writetable(checks,fullfile(run_dir,'C_checks.csv'));
  assert(all(checks.pass),'FROZEN_NUMERICAL_GATE_FAILED'); completed=completed+1;
  save_checkpoint('IN_PROGRESS');
 end
 status='S4062p5_C801_NATIVE_AND_FIELD_PASS';
catch e
 err=getReport(e,'extended','hyperlinks','off'); fprintf('%s\n',err); status='P1_NO_GO_STOPPED';
end
save_checkpoint(status);
fprintf('P1_C status=%s started=%d completed=%d failed=%d wall=%.3fs\n',status,started,completed,nnz(~checks.pass),toc(start));
 function add(name,value,lim)
  checks(end+1,:)={string(name),isscalar(value)&&isfinite(value)&&value<=lim,value,lim};
 end
 function save_checkpoint(st)
  receipt=struct('status',st,'started_slots',started,'completed_slots',completed,'elapsed_s',toc(start),'budget_remaining_before_this_call_s',remaining_s,'failed_checks',nnz(~checks.pass)+double(~isempty(err)),'error',err,'results',results,'other_jobs_this_invocation',0,'next_allowed_step','manual review after any failure; no P2-P5','physical_config','U0 point3D 1m, width N/A, D1500 Tx1490 Rx6, x=[0,10,50,100,200], c1500 rho1 bottom1600/rho1.8/loss0','matlab_version',version);
  save(fullfile(run_dir,'C_results.mat'),'pressure','slots','results','checks','receipt','x','v'); writetable(checks,fullfile(run_dir,'C_checks.csv')); fid=fopen(fullfile(run_dir,'checkpoint_C.json'),'w'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid);
 end
end
function n=bitdiff(a,b)
assert(numel(a)==numel(b)); n=nnz(typecast(a(:),'uint64')~=typecast(b(:),'uint64'));
end