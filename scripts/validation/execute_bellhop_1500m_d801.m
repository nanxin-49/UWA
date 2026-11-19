function receipt=execute_bellhop_1500m_d801(project_root,run_dir)
%EXECUTE_BELLHOP_1500M_D801 One new native C801, reuse saved C201 and C401.
timer=tic; addpath(fullfile(project_root,'scripts','validation','support')); addpath(fullfile(run_dir,'code'),'-begin');
base201=load(fullfile(run_dir,'baseline201','independent_audit.mat'),'p','ref','p_native');
base401=load(fullfile(run_dir,'baseline401','comparison.mat'),'p201','p401','ref');
plan=load(fullfile(run_dir,'plan.mat')); cfg=plan.slots.cfg; oldcfg=cfg; oldcfg.beam_count=201;
assert(isequaln(oldcfg,plan.frozen.cfg),'Only beam count may change');
assert(isequal(cfg.nominal_receiver_x_m,[0 10 50 100 200])&&cfg.frequency_hz==4000&&cfg.step_m==.5&&cfg.ssp_spacing_m==.625);
r=run_bellhop_1500m_d801_case(project_root,run_dir,160);
receipt=struct('status','PROCESS_OR_INPUT_FAILED_STOPPED','new_Bellhop_jobs',r.started_slots,'other17_slots_run',0,'P2_run',false,'beams',[201 401 801],'old201_recomputed',false,'old401_recomputed',false,'exception',r.error,'figure_generated',false,'physical_pass',false,'source_width','N/A ideal omni 3D point; pressure_1m','ELF_SHA',plan.slots.binary_sha256);
root=fullfile(run_dir,'jobs','01_D_4000_C801','case');
if ~isfile(fullfile(fileparts(root),'component.mat'))
 writejson(fullfile(run_dir,'comparison_receipt.json'),receipt); writefile(fullfile(run_dir,'chinese_report.md'),['# D4000/C801停止报告' newline newline '输入或执行失败；不重试、不1601。' newline r.error]); return
end
s=load(fullfile(fileparts(root),'component.mat')); checks=readtable(fullfile(run_dir,'C_checks.csv'),'Delimiter',',','TextType','string');
raw=read_bellhop_rotated_field(root); fid=fopen([root '.shd'],'rb'); assert(fid>=0); native=complex(zeros(5,1));
for j=1:5, fseek(fid,(raw.range_record+j)*raw.record_bytes,'bof'); a=fread(fid,2,'float32'); native(j)=complex(a(1),a(2)); end; fclose(fid);
assert(isequal(native,raw.pressure)&&isequal(-conj(native),s.p(:)),'Independent new SHD/MAT mismatch');
x=[0 10 50 100 200]'; L=hypot(1484,x); analytic=exp(1i*2*pi*4000/1500*L)./L; ref=base401.ref(:); E=abs(ref);
assert(max(abs(analytic-ref))<1e-17,'Frozen analytic reference independent audit');
p201=base201.p(:); p401=base401.p401(:); p801=-conj(native); P=[p201 p401 p801]; beams=[201 401 801];
assert(isequal(base401.p201(:),p201)&&max(abs(base201.ref(:)-ref))<1e-14,'Historical reference alignment');
for k=1:2
 nb=beams(k); oldroot=fullfile(run_dir,sprintf('baseline%d',nb),'case'); oldraw=read_bellhop_rotated_field(oldroot); assert(isequal(-conj(oldraw.pressure(:)),P(:,k)),'Historical raw SHD exact check');
end
t=readtable(fullfile(fileparts(root),'field.csv'),'Delimiter',','); assert(max(abs(t.p_real+1i*t.p_imag-p801))<1e-17,'New CSV/MAT mismatch');
M=struct([]); h=zeros(3,1); D=cell(1,3); fields=table(); gridtable=table();
for k=1:3
 p=P(:,k); nb=beams(k); mk=measure(p,ref,E); if k==1,M=mk;else,M(k)=mk;end
 fields=[fields;table(x,nb*ones(5,1),real(p),imag(p),real(ref),imag(ref),100*(abs(p./ref)-1),angle(p.*conj(ref)),abs(p-ref)./E,'VariableNames',{'x_m','beams','p_real','p_imag','ref_real','ref_imag','amplitude_error_percent','phase_error_rad','complex_error_over_E'})]; %#ok<AGROW>
 if k==3, dr=[root '.r25diag']; else, dr=fullfile(run_dir,sprintf('baseline%d',nb),'case.r25diag'); end
 fid=fopen(dr,'r'); fgetl(fid); d=fscanf(fid,'%f',[10 inf]).'; fclose(fid); D{k}=d; assert(size(d,1)==nb);
 alpha=deg2rad(d(:,1)); da=(alpha(end)-alpha(1))/(nb-1); h(k)=da; theta=atan2(x,1484); pos=(theta-alpha(1))/da; nearest=round(pos); mask=abs(pos-nearest)<1e-10; pos(mask)=nearest(mask); frac=pos-floor(pos);
 prediction=-(2*pi*4000/1500)*L/2*da^2.*frac.*(1-frac); phase=angle(p.*conj(ref));
 gridtable=[gridtable;table(x,nb*ones(5,1),rad2deg(theta),frac,repmat(da,5,1),prediction,phase,phase-prediction,'VariableNames',{'x_m','beams','required_angle_deg','angular_cell_fraction','delta_angle_rad','second_order_phase_prediction_rad','measured_phase_error_rad','measured_minus_prediction_rad'})]; %#ok<AGROW>
end
assert(max(abs(D{1}(:,1)-D{3}(1:4:end,1)))<1e-12&&max(abs(D{2}(:,1)-D{3}(1:2:end,1)))<1e-12,'Nested frozen angular fan');
metrics=table(beams',h,rad2deg(h),[M.L2]',[M.maxDE]',[M.maxphase]',[M.support]',[M.L2]'<=.001,[M.maxDE]'<=.002,[M.maxphase]'<=.01,'VariableNames',{'beams','delta_angle_rad','delta_angle_deg','complex_L2','max_complex_over_E','max_phase_rad','phase_support','L2_pass','max_complex_pass','phase_pass'});
convergence=table([201;401],[401;801],log(metrics.complex_L2(1:2)./metrics.complex_L2(2:3))./log(h(1:2)./h(2:3)),log(metrics.max_complex_over_E(1:2)./metrics.max_complex_over_E(2:3))./log(h(1:2)./h(2:3)),'VariableNames',{'coarse_beams','fine_beams','observed_L2_order','observed_max_error_order'});
pair=table();
for k=1:2
 dp=P(:,k+1)-P(:,k); pair=[pair;table(x,repmat(beams(k),5,1),repmat(beams(k+1),5,1),real(dp),imag(dp),abs(dp)./E,angle(P(:,k+1).*conj(P(:,k))),'VariableNames',{'x_m','coarse_beams','fine_beams','delta_p_real','delta_p_imag','abs_delta_over_E','phase_fine_minus_coarse_rad'})]; %#ok<AGROW>
end
writetable(metrics,fullfile(run_dir,'metrics_201_401_801.csv')); writetable(fields,fullfile(run_dir,'fields_201_401_801.csv')); writetable(gridtable,fullfile(run_dir,'angular_grid_comparison.csv')); writetable(convergence,fullfile(run_dir,'observed_orders.csv')); writetable(pair,fullfile(run_dir,'pairwise_complex_differences.csv'));
reload=readtable(fullfile(run_dir,'fields_201_401_801.csv'),'Delimiter',','); assert(height(reload)==15&&max(abs(reload.p_real+1i*reload.p_imag-[p201;p401;p801]))<1e-17,'Independent three-level CSV reload');
pass=all(checks.pass)&&M(3).L2<=.001&&M(3).maxDE<=.002&&M(3).maxphase<=.01;
receipt.status='D801_NUMERICAL_FAIL_STOPPED'; if pass,receipt.status='D801_DIRECT_CONTROL_GO_P1_INCOMPLETE';end
receipt.physical_pass=pass; receipt.native_exit_code=s.proc.exit_code; receipt.native_elapsed_s=s.proc.elapsed_s; receipt.process=s.proc; receipt.old201=M(1); receipt.old401=M(2); receipt.new801=M(3);
receipt.failed_numeric_gate_count=nnz(~[M(3).L2<=.001 M(3).maxDE<=.002 M(3).maxphase<=.01]); receipt.failed_all_saved_checks=nnz(~checks.pass); receipt.frozen_reference_reused_exactly=true;
receipt.L2_order_201_401=convergence.observed_L2_order(1); receipt.L2_order_401_801=convergence.observed_L2_order(2); receipt.second_order_prediction_max_phase_residual_rad=max(abs(gridtable.measured_minus_prediction_rad));
receipt.pair401_801_complex_L2=norm(p801-p401)/norm(p401); receipt.pair401_801_max_difference_over_E=max(abs(p801-p401)./E); receipt.pair_is_diagnostic_not_acceptance_gate=true;
receipt.independent_SHD_MAT_CSV_audit='3 levels PASS'; receipt.actual_field_coordinates='REAL8 bitwise, exact original ordering and proper-rotation roundtrip PASS'; receipt.next_allowed_step='manual review; no 1601, other17, P2, FFP or Munk';
save(fullfile(run_dir,'comparison.mat'),'p201','p401','p801','ref','E','x','fields','metrics','gridtable','convergence','pair','receipt','native','cfg','checks','D');
if pass
 fig=figure('Visible','off','Color','w','Position',[60 60 1120 760]); tl=tiledlayout(fig,2,2,'Padding','compact','TileSpacing','compact');
 vals={real(P-ref)./E,imag(P-ref)./E,abs(P-ref)./E,angle(P.*conj(ref))}; labels={'Re(complex residual)/E','Im(complex residual)/E','|complex residual|/E','Phase error (rad)'};
 for j=1:4
  ax=nexttile(tl); hh=plot(ax,x,vals{j}(:,1),'k--o',x,vals{j}(:,2),'b-.s',x,vals{j}(:,3),'r-^');
  for k=1:3,set(hh(k),'Tag',sprintf('%d_%d',beams(k),j));end
  xlabel(ax,'Independent link offset x (m)'); ylabel(ax,labels{j}); grid(ax,'on'); legend(ax,{'201 saved','401 saved','801 new'},'Location','best'); if j==3,yline(ax,.002,':','Frozen max gate');end
 end
 title(tl,'U0 direct control, 4000 Hz: 201 / 401 / 801 beams'); drawnow; savefig(fig,fullfile(run_dir,'F01_complex_error_201_401_801.fig')); exportgraphics(fig,fullfile(run_dir,'F01_complex_error_201_401_801.png'),'Resolution',140);close(fig);
 fig=openfig(fullfile(run_dir,'F01_complex_error_201_401_801.fig'),'invisible');for j=1:4,for k=1:3,obj=findobj(fig,'Tag',sprintf('%d_%d',beams(k),j));assert(isscalar(obj)&&isequal(obj.XData,x')&&isequal(obj.YData,vals{j}(:,k)'));end;end;close(fig);
 fig=figure('Visible','off','Color','w','Position',[80 80 1050 460]); tl=tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact'); ax=nexttile(tl); ref2=metrics.complex_L2(1)*(h/h(1)).^2;
 hl=loglog(ax,h,metrics.complex_L2,'k-o',h,metrics.max_complex_over_E,'b-s',h,ref2,'r--');set(hl(1),'Tag','L2_scaling');set(hl(2),'Tag','max_scaling');set(hl(3),'Tag','h2_guide'); xlabel(ax,'Angular spacing (rad)');ylabel(ax,'Dimensionless complex error');legend(ax,{'Complex L2','Max difference / E','h^2 guide anchored at saved 201 L2'},'Location','best');grid(ax,'on');
 ax=nexttile(tl); phasevals=[gridtable.measured_phase_error_rad(1:5),gridtable.measured_phase_error_rad(6:10),gridtable.measured_phase_error_rad(11:15)]; pred=[gridtable.second_order_phase_prediction_rad(1:5),gridtable.second_order_phase_prediction_rad(6:10),gridtable.second_order_phase_prediction_rad(11:15)];
 hh=plot(ax,pred(:,1),phasevals(:,1),'ko',pred(:,2),phasevals(:,2),'bs',pred(:,3),phasevals(:,3),'r^'); for k=1:3,set(hh(k),'Tag',sprintf('phase_%d',beams(k)));end; hold(ax,'on');lims=[min(pred,[],'all'),0]; plot(ax,lims,lims,'k:');xlabel(ax,'Second-order projected phase prediction (rad)');ylabel(ax,'Measured phase error (rad)');grid(ax,'on');legend(ax,{'201','401','801','Equality'},'Location','best');title(tl,'Three frozen angular grids; cell position is retained');drawnow;
 savefig(fig,fullfile(run_dir,'F02_angular_error_scaling.fig'));exportgraphics(fig,fullfile(run_dir,'F02_angular_error_scaling.png'),'Resolution',140);close(fig);
 fig=openfig(fullfile(run_dir,'F02_angular_error_scaling.fig'),'invisible'); tags={'L2_scaling','max_scaling','h2_guide'};ys={metrics.complex_L2,metrics.max_complex_over_E,ref2};for z=1:3,o=findobj(fig,'Tag',tags{z});assert(isscalar(o)&&isequal(o.XData,h')&&isequal(o.YData,ys{z}'));end; for k=1:3,obj=findobj(fig,'Tag',sprintf('phase_%d',beams(k)));assert(isscalar(obj)&&isequal(obj.XData,pred(:,k)')&&isequal(obj.YData,phasevals(:,k)'));end;close(fig);
 receipt.figure_generated=true; receipt.FIG_reopen_audit='12 residual + 3 scaling + 3 prediction curves PASS';
end
receipt.MATLAB_workflow_elapsed_s=toc(timer); receipt.matlab_version=version; writejson(fullfile(run_dir,'comparison_receipt.json'),receipt);
report=string('# D/4000Hz/C/801束最后一次最小加密验证')+newline+newline+'日期：2026-11-19。结论：'+receipt.status+'。'+newline+newline+ ...
 '本轮只新增一个801束直达C槽，原201、401保存SHD/MAT/CSV只读复用，不重算。U0等声速1500m/s、水深1500m/Tx1490m/Rx6m、原5个独立链路位置、全向3D点源1m参考/源宽N/A、原扇区±8.443110798927421度、0.5m步长/0.625m Q、REAL8正式C ELF和固定-conj转换全部保持。只改变束数，不修改曲率或其他物理公式。'+newline+newline+ ...
 '## 冻结门与三档结果'+newline+newline+'|束数|角间隔/rad|复L2（门0.001）|max复差/E（门0.002）|max相位/rad（门0.01）|'+newline+'|---:|---:|---:|---:|---:|'+newline;
for k=1:3,report=report+sprintf('|%d|%.12g|%.12g|%.12g|%.12g|\n',beams(k),h(k),M(k).L2,M(k).maxDE,M(k).maxphase);end
report=report+newline+sprintf('801三项复场门失败%d项，全部保存检查失败%d项。原生exit=%d，native=%.6fs；201/401历史失败保持，不追改。\n',receipt.failed_numeric_gate_count,receipt.failed_all_saved_checks,s.proc.exit_code,s.proc.elapsed_s);
report=report+newline+'## 五点801幅相误差'+newline+newline+'|x/m|幅度误差/%|相位误差/rad|复数残差/E|'+newline+'|---:|---:|---:|---:|'+newline;
for j=1:5,idx=10+j;report=report+sprintf('|%g|%.12g|%.12g|%.12g|\n',x(j),fields.amplitude_error_percent(idx),fields.phase_error_rad(idx),fields.complex_error_over_E(idx));end
report=report+newline+'## 角间隔与二阶投影规律'+newline+newline+sprintf('三档角间隔逐次减半，原201和401角度嵌套于801，最大差小于1e-12度。复L2观测阶201→401为%.8g、401→801为%.8g。401/801复声压差L2=%.12g、max差/E=%.12g，仅作诊断，不另设验收门。\n',receipt.L2_order_201_401,receipt.L2_order_401_801,receipt.pair401_801_complex_L2,receipt.pair401_801_max_difference_over_E);
report=report+'采用冻结的二阶投影预测 delta_phi≈−kL h² t(1−t)/2，保留每点在角单元中的分数t；三个等级实测与预测相位的最大差为'+sprintf('%.12g',receipt.second_order_prediction_max_phase_residual_rad)+'rad。h平方引导线锚定201误差，仅作尺度参考；t随加密变化，不能要求每点严格四倍。三档有限控制支持该近似机制的下降规律，不构成任意环境或无限加密的收敛证明。'+newline+newline;
if pass,report=report+'801束通过本直达控制工况三门及坐标/几何检查，可报告这个固定U0/4000Hz/5点工况的数值资格；**不宣布整个P1通过**。两组PNG/可编辑FIG及曲线重开审计已保存，PNG目视检查由总监督器回收后等待完成。'+newline;else,report=report+'801仍失败，立即停止新计算；不自动1601、不修改曲率公式、不运行其他17槽或P2。'+newline;end
report=report+newline+'独立SHD/MAT/CSV复核、实际REAL8坐标逐位检查、终态几何及原始PRT/ENV/SSP/R25D/真实exit/资源记录均保存。C模式仅SHD，不启动额外A/R。完整300秒监督回执为supervisor_receipt.json。所有历史SHA本地执行前后保护。未运行FFP/Munk、未修改官方Bellhop/PE/公共接口，不提交推送；等待人工审查。'+newline;
writefile(fullfile(run_dir,'chinese_report.md'),char(report)); fprintf('D801_RESULT status=%s L2=%.12g maxDE=%.12g phase=%.12g orders=%.6g/%.6g native=%.6fs failed%d\n',receipt.status,M(3).L2,M(3).maxDE,M(3).maxphase,receipt.L2_order_201_401,receipt.L2_order_401_801,s.proc.elapsed_s,receipt.failed_numeric_gate_count);
end
function m=measure(p,ref,E)
support=abs(p)>=.01*E&abs(ref)>=.01*E;m=struct('L2',norm(p-ref)/norm(ref),'maxDE',max(abs(p-ref)./E),'maxphase',max(abs(angle(p(support).*conj(ref(support))))),'support',nnz(support));
end
function writejson(file,val),writefile(file,jsonencode(val,PrettyPrint=true));end
function writefile(file,val),fid=fopen(file,'w','n','UTF-8');assert(fid>=0);c=onCleanup(@()fclose(fid));fprintf(fid,'%s',val);end