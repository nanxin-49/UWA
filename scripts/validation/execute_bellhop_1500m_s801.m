function receipt=execute_bellhop_1500m_s801(project_root,run_dir)
%EXECUTE_BELLHOP_1500M_S801 One new S C801; reuse the qualified saved D801.
timer=tic;addpath(fullfile(project_root,'scripts','validation','support'));addpath(fullfile(run_dir,'code'),'-begin');
assert(license('test','MATLAB')==1,'MATLAB license gate');fprintf('MATLAB_LICENSE_PASS %s\n',version);
f=jsondecode(fileread(fullfile(run_dir,'frozen_task.json'))); original=load(fullfile(run_dir,'original_plan.mat'),'slots','x','v');frozen=original.slots(3);
assert(strcmp(frozen.path,'S')&&frozen.frequency_hz==4000&&strcmp(frozen.mode,'C')&&frozen.beams==201);
slots=frozen;slots.cfg.beam_count=801;slots.beams=801;slots.output_dir=fullfile(run_dir,'jobs','01_S_4000_C801');slots.local_dir=slots.output_dir;slots.name='01_S_4000_C801';cfg=slots.cfg;x=original.x;v=original.v;
assert(f.candidate_C_beams==801&&cfg.source_z_m==1490&&cfg.water_depth_m==1500&&cfg.receiver_range_m==1496&&cfg.frequency_hz==4000&&cfg.step_m==.5&&cfg.ssp_spacing_m==.625);
root=fullfile(slots.output_dir,'case'); oe=fileread(fullfile(run_dir,'originalS201','case.env'));ne=fileread([root '.env']);assert(strcmp(strrep(ne,[newline '801' newline],[newline '201' newline]),oe),'Only beam line may change');
for ext={'.ssp','.r25d','.iwcov','.rotq'},assert(strcmpi(black_file_hash([root ext{1}]),black_file_hash(fullfile(run_dir,'originalS201',['case' ext{1}]))));end
save(fullfile(run_dir,'plan.mat'),'slots','frozen','x','v');
r=run_bellhop_1500m_s801_case(project_root,run_dir,160);
receipt=struct('status','S801_EXECUTION_OR_INPUT_FAIL_STOPPED','physical_pass',false,'figure_generated',false,'new_Bellhop_jobs',r.started_slots,'D_recomputed',false,'B_BS_run',false,'other_P1_jobs_run',0,'P2_P5_run',false,'exception',r.error,'candidate_C_beams',801,'source_width','N/A 3D omni point, pressure_1m');
if ~isfile(fullfile(slots.output_dir,'component.mat')),writejson(fullfile(run_dir,'comparison_receipt.json'),receipt);writefile(fullfile(run_dir,'chinese_report.md'),['# S4000/C801执行停止' newline r.error]);return;end
s=load(fullfile(slots.output_dir,'component.mat'));checks=readtable(fullfile(run_dir,'C_checks.csv'),'Delimiter',',','TextType','string');
iw=readlog([root '.iwdiag']);mw=readlog([root '.mwdiag']);d=readlog([root '.r25diag']);assert(height(iw)==801&&height(mw)==801&&height(d)==801);
geom=table(strings(0,1),false(0,1),zeros(0,1),zeros(0,1),'VariableNames',{'check','pass','value','limit'});
a=deg2rad(iw.alpha_deg);hitExpected=1490*tan(a);inc=[cos(a),sin(a)];refl=[-cos(a),sin(a)];rot=[cos(a),-sin(a)];
gate('all_first_sea_counts',nnz(mw.requested_order~=1|mw.walls_done~=1|mw.surface_bounces~=1|mw.bottom_bounces~=0),0);
gate('terminal_physical_bounces',nnz(d.top_bounces~=1|d.bottom_bounces~=0),0);
gate('reflection_wall_physical_depth_m',max(abs(iw.hit_r-1490)),1e-4);
gate('reflection_point_x_m',max(abs(iw.hit_z-hitExpected)),1e-4);
gate('incident_direction_rad',max(direrr([iw.inc_ur iw.inc_uz],inc)),1e-6);gate('reflected_direction_rad',max(direrr([iw.ref_ur iw.ref_uz],refl)),1e-6);gate('proper_rotated_direction_rad',max(direrr([iw.rot_ur iw.rot_uz],rot)),1e-6);
gate('surface_path_length_m',max(abs(1500*iw.tau_receiver_real-1496./cos(a))),1e-4);gate('surface_arrival_range_time_s',max(abs(iw.tau_receiver_real-1496./cos(a)/1500)),1e-7);gate('time_to_wall_s',max(abs(iw.tau_wall_real-1490./cos(a)/1500)),1e-7);
coeff=(iw.amp_ref./iw.amp_in).*exp(1i*(iw.phase_ref-iw.phase_in));gate('pressure_release_R_minus1',max(abs(coeff+1)),1e-12);gate('vacuum_amplitude_preserved',max(abs(iw.amp_ref-iw.amp_in)),1e-12);gate('pressure_release_phase_pi_rad',max(abs(angle(exp(1i*(iw.phase_ref-iw.phase_in-pi))))),1e-12);
gate('flat_surface_curvature',max(abs(iw.kappa)),0);gate('segment_membership',nnz(iw.wall_seg~=1&iw.wall_seg~=2),0);gate('segment_parameter_range',max([0;-iw.wall_lambda;iw.wall_lambda-1]),1e-10);
seg=round(iw.wall_seg);nodes=[-1000;0;1000];qlambda=nodes(seg)+1000*iw.wall_lambda;gate('segment_parameter_point_m',max(abs(qlambda-iw.hit_z)),1e-4);
writetable(iw,fullfile(run_dir,'all801_first_surface_states.csv'));writetable(geom,fullfile(run_dir,'surface_geometry_checks.csv'));
raw=read_bellhop_rotated_field(root);fid=fopen([root '.shd'],'rb');native=complex(zeros(5,1));for j=1:5,fseek(fid,(raw.range_record+j)*raw.record_bytes,'bof');z=fread(fid,2,'float32');native(j)=complex(z(1),z(2));end;fclose(fid);assert(isequal(native,raw.pressure),'Independent S SHD parser');
pS=zeros(5,1);pS(cfg.receiver_x_order)=-conj(native);assert(isequal(pS,s.p(:)),'S physical receiver permutation');
b=load(fullfile(run_dir,'baselineD801','comparison.mat'),'p801','ref');pD=b.p801(:);refD=b.ref(:);dr=read_bellhop_rotated_field(fullfile(run_dir,'baselineD801','case'));assert(isequal(-conj(dr.pressure(:)),pD),'D saved raw/MAT audit');
x=[0 10 50 100 200]';LD=hypot(1484,x);LS=hypot(1496,x);refS=-exp(1i*2*pi*4000/1500*LS)./LS;refDS=refD+refS;pDS=pD+pS;ES=abs(refS);EDS=abs(refD)+ES;
ms=measure(pS,refS,ES,x);mds=measure(pDS,refDS,EDS,x);metric=table(["S";"D_plus_S"],[ms.L2;mds.L2],[ms.maxDE;mds.maxDE],[ms.maxphase;mds.maxphase],[ms.support;mds.support],[ms.L2<=.001;mds.L2<=.001],[ms.maxDE<=.002;mds.maxDE<=.002],[ms.maxphase<=.01;mds.maxphase<=.01],'VariableNames',{'field','complex_L2','max_complex_over_E','max_phase_rad','phase_support','L2_pass','max_complex_pass','phase_pass'});
writetable(metric,fullfile(run_dir,'S_and_DS_metrics.csv'));
fields=table();for k=1:2,if k==1,p=pS;ref=refS;E=ES;lab="S";else,p=pDS;ref=refDS;E=EDS;lab="D_plus_S";end;support=abs(p)>=.01*E&abs(ref)>=.01*E;fields=[fields;table(x,repmat(lab,5,1),real(p),imag(p),real(ref),imag(ref),E,100*(abs(p)-abs(ref))./E,abs(p-ref)./E,angle(p.*conj(ref)),support,'VariableNames',{'x_m','field','p_real','p_imag','reference_real','reference_imag','independent_envelope','amplitude_error_percent_E','complex_error_over_E','raw_phase_error_rad','phase_support'})];end %#ok<AGROW>
writetable(fields,fullfile(run_dir,'S_and_DS_fields.csv'));
reflectx=1490*x/1496;analytic_geometry=table(x,reflectx,zeros(5,1),LS,LS/1500,LD/1500,(LS-LD)/1500,'VariableNames',{'receiver_x_m','analytic_reflection_x_m','surface_z_m','S_length_m','S_delay_s','D_delay_s','S_minus_D_delay_s'});writetable(analytic_geometry,fullfile(run_dir,'analytic_five_link_geometry.csv'));
reload=readtable(fullfile(run_dir,'S_and_DS_fields.csv'),'Delimiter',',');assert(height(reload)==10&&max(abs(reload.p_real+1i*reload.p_imag-[pS;pDS]))<1e-17,'S/DS CSV reload audit');
pass=all(checks.pass)&&all(geom.pass)&&all(metric.L2_pass)&&all(metric.max_complex_pass)&&all(metric.phase_pass);
receipt.status='S801_OR_DS_NUMERICAL_FAIL_STOPPED';if pass,receipt.status='S801_AND_DS_CONTROL_GO_P1_INCOMPLETE';end
receipt.physical_pass=logical(pass);receipt.native_exit_code=s.proc.exit_code;receipt.native_elapsed_s=s.proc.elapsed_s;receipt.process=s.proc;receipt.S=ms;receipt.DS=mds;receipt.failed_S_field_checks=nnz(~checks.pass);receipt.failed_surface_geometry_checks=nnz(~geom.pass);receipt.failed_DS_and_S_metric_gates=nnz(~[metric.L2_pass;metric.max_complex_pass;metric.phase_pass]);receipt.phase_mask_rule='both amplitudes >=0.01 independent envelope';receipt.combined_field='U0_D_S_pressure_release_halfspace_control';receipt.max_R_plus1=max(abs(coeff+1));receipt.max_reflection_point_error_m=max(abs(iw.hit_z-hitExpected));receipt.max_reflection_direction_error_rad=max(direrr([iw.ref_ur iw.ref_uz],refl));receipt.max_surface_arrival_time_error_s=max(abs(iw.tau_receiver_real-1496./cos(a)/1500));receipt.pq_diagnostics=struct('max_p_ref_error',max(iw.p_ref_error),'max_q_ref_error',max(iw.q_ref_error),'max_p_rot_error',max(iw.p_rot_error),'max_q_rot_error',max(iw.q_rot_error));receipt.independent_SHD_MAT_CSV_audit='PASS';receipt.scope_limit='C first-event diagnostics are actual native rays; five-link reflection geometry CSV is analytic, not ARR/eigenray output';
save(fullfile(run_dir,'comparison.mat'),'pD','pS','pDS','refD','refS','refDS','ES','EDS','x','LD','LS','ms','mds','metric','fields','geom','iw','mw','d','analytic_geometry','receipt','cfg');
if pass
 fig=figure('Visible','off','Color','w','Position',[50 50 1280 560]);ax1=axes('Parent',fig,'Position',[.075 .22 .39 .68]);ax2=axes('Parent',fig,'Position',[.575 .22 .39 .68]);
 h=plot(ax1,x,abs(pS-refS)./ES,'b-o',x,abs(pDS-refDS)./EDS,'r-s');set(h(1),'Tag','S_error');set(h(2),'Tag','DS_error');yline(ax1,.002,'k:','Frozen max gate','HandleVisibility','off');xlabel(ax1,'Independent receiver offset x (m)');ylabel(ax1,'|complex error| / independent envelope');grid(ax1,'on');
 lg=legend(ax1,h,{'S: E=|S analytic|','D+S: E=|D analytic|+|S analytic|'},'Orientation','horizontal','AutoUpdate','off');set(lg,'Units','normalized','Position',[.035 .025 .49 .06],'FontSize',9);
 h2=semilogy(ax2,x,abs(refDS),'k-o',x,abs(pDS),'b--s',x,EDS,'r:',x,.01*EDS,'g-.');tags={'DS_analytic','DS_numeric','envelope','phase_threshold'};for k=1:4,set(h2(k),'Tag',tags{k});end;xlabel(ax2,'Independent receiver offset x (m)');ylabel(ax2,'Pressure / 1 m reference');grid(ax2,'on');lg2=legend(ax2,h2,{'Analytic D+S','Bellhop D+S','Independent E','0.01 E phase threshold'},'NumColumns',2,'AutoUpdate','off');set(lg2,'Units','normalized','Position',[.56 .014 .43 .105],'FontSize',9);
 annotation(fig,'textbox',[0 .95 1 .04],'String','U0 4000 Hz: one S801 run plus saved D801; coherent sum before magnitude','HorizontalAlignment','center','EdgeColor','none','FontSize',13);drawnow;savefig(fig,fullfile(run_dir,'F01_S_and_DS_control.fig'));exportgraphics(fig,fullfile(run_dir,'F01_S_and_DS_control.png'),'Resolution',140);close(fig);
 fig=openfig(fullfile(run_dir,'F01_S_and_DS_control.fig'),'invisible');expected={{'S_error',abs(pS-refS)./ES},{'DS_error',abs(pDS-refDS)./EDS},{'DS_analytic',abs(refDS)},{'DS_numeric',abs(pDS)},{'envelope',EDS},{'phase_threshold',.01*EDS}};for j=1:numel(expected),o=findobj(fig,'Tag',expected{j}{1});assert(isscalar(o)&&isequal(o.XData,x')&&isequal(o.YData,expected{j}{2}'));end;lg=findall(fig,'Type','legend');for j=1:numel(lg),p=lg(j).Position;assert(p(1)>=0&&p(2)>=0&&p(1)+p(3)<=1&&p(2)+p(4)<=1&&strcmp(lg(j).AutoUpdate,'off'));end;close(fig);receipt.figure_generated=true;receipt.FIG_reopen_audit='6 curves and 2 bounded legends PASS';
end
receipt.MATLAB_workflow_elapsed_s=toc(timer);receipt.matlab_version=version;writejson(fullfile(run_dir,'comparison_receipt.json'),receipt);
report=string('# 修订P1最小海面反射验证：S4000/C801及D+S')+newline+newline+'2026-11-19。结论：'+receipt.status+'。'+newline+newline+ ...
 '历史三组PNG/FIG图例已另存修复，24曲线、4图例及6原文件SHA审计/目视完成，无201/401/801重算。该阶段使用本地既有MCP；曾出现许可证警告，未关闭或重启共享会话。'+newline+newline+ ...
 '## 冻结配置与许可范围'+newline+newline+'原201束P1失败保持，修订候选C束数801；本轮只新增S/4000Hz/C801一项，复用通过的D4000/C801。U0水深1500m、Tx1490m、Rx6m、x=0/10/50/100/200m；c1500m/s、水密度1g/cm3/吸收0；流体底1600m/s、1.8g/cm3、损耗0但本两路径不含底反射。理想全向3D点源，1m参考，源宽N/A。原S扇区±1.1atan2d(200,1496)约8.376183136度，步长0.5m/Q0.625m、REAL8正式C ELF与p_abs=-conj(p_native)保持。ENV仅改原S201的束数行，其余原输入SHA冻结。'+newline+newline+ ...
 '## 几何、反射次数与压力释放相位'+newline+newline+sprintf('实际C模式801束全部一次海面、零底反射；独立事件检查失败%d项。反射点x最大误差%.12gm，反射方向最大误差%.12grad，映射接收range处时延最大误差%.12gs；R=-1复系数最大差%.12g。\n',receipt.failed_surface_geometry_checks,receipt.max_reflection_point_error_m,receipt.max_reflection_direction_error_rad,receipt.max_surface_arrival_time_error_s,receipt.max_R_plus1);
report=report+'原生iwdiag/mwdiag检查与逐束状态表均保留；采用L_S=hypot(x,1496)、反射点x_s=1490x/1496独立几何。几何表五个位置为解析参考，不冒充ARR或精确本征声线；原生C事件几何覆盖所有实际801束。source V墙、一次海面计数和相位跃迁pi共同核对边界归属。'+newline+newline+ ...
 '## 独立复压与半空间镜像比较'+newline+newline+'S解析为-exp(ikL_S)/L_S；D+S解析为冻结D解析减exp(ikL_S)/L_S。所有总场先做复数相干相加，未相加dB或幅度包络。S包络为|S解析|，D+S独立包络E=|D解析|+|S解析|；所有五点复值均保留，不用近零总场作逐点相对误差分母。相位仅双方幅度均>=0.01E的点评价。'+newline+newline+'|场|复L2 ≤0.001|max复差/E ≤0.002|max相位/rad ≤0.01|相位支持点|'+newline+'|---|---:|---:|---:|---:|'+newline;
report=report+sprintf('|S|%.12g|%.12g|%.12g|%d|\n|D+S|%.12g|%.12g|%.12g|%d|\n',ms.L2,ms.maxDE,ms.maxphase,ms.support,mds.L2,mds.maxDE,mds.maxphase,mds.support);
report=report+'D+S相位有效x/m：'+join(string(mds.phase_x_m),', ')+'。原始五点复压及支持标记见S_and_DS_fields.csv，状态、几何、MAT与SHA均保留。'+newline+newline;
if pass,report=report+'S及D+S最小海面控制通过；仅本固定U0/4000Hz/五独立位置资格。整个P1仍未通过：B/BS、另一频率、A/R和未执行对照仍缺。'+newline;else,report=report+'存在必需门失败，停止；不加束、不修改公式，不重试。'+newline;end
report=report+newline+'## 资源、审计及停止'+newline+newline+sprintf('唯一新增Bellhop实际exit%d、native%.6fs；原始ENV/SSP/R25D/ROTQ/IWCOV/SHD/PRT、事件日志和真实退出回执保存。统一300s监督包含启动/求解/回收/审计/图件/报告与导航，295s动作截止、远端170s、native120s，启动RAM1.5GiB、运行256MiB及30%%、空闲盘20GiB/输出2GiB门保持，终态见supervisor_receipt.json。\n',s.proc.exit_code,s.proc.elapsed_s);
report=report+'保留原201失败、历史超限、原图件及D801资格。D+S仅标记U0_D_S_pressure_release_halfspace_control，不是含底反射的完整U0波导场、无限阶通信信道或CIR。不启动B/BS、其他P1、P2–P5、FFP/Munk；不改官方Bellhop/PE/公共接口，不提交推送，等待人工审查。'+newline;
writefile(fullfile(run_dir,'chinese_report.md'),char(report));fprintf('S801_RESULT status=%s S_L2%.12g DS_L2%.12g failedGeom%d native%.6fs\n',receipt.status,ms.L2,mds.L2,receipt.failed_surface_geometry_checks,s.proc.elapsed_s);
 function gate(name,value,limit),geom(end+1,:)={string(name),isscalar(value)&&isfinite(value)&&value<=limit,value,limit};end
end
function t=readlog(file)
fid=fopen(file,'r');assert(fid>=0);hdr=fgetl(fid);names=strsplit(strtrim(regexprep(hdr,'^#','')));data=fscanf(fid,'%f',[numel(names) inf]).';fclose(fid);assert(size(data,2)==numel(names));t=array2table(data,'VariableNames',names);
end
function e=direrr(u,v),e=atan2(abs(u(:,1).*v(:,2)-u(:,2).*v(:,1)),sum(u.*v,2));end
function m=measure(p,ref,E,x)
support=abs(p)>=.01*E&abs(ref)>=.01*E;phi=abs(angle(p(support).*conj(ref(support))));maxphi=NaN;if ~isempty(phi),maxphi=max(phi);end;m=struct('L2',norm(p-ref)/norm(ref),'maxDE',max(abs(p-ref)./E),'maxphase',maxphi,'support',nnz(support),'phase_x_m',x(support).');
end
function writejson(file,val),writefile(file,jsonencode(val,PrettyPrint=true));end
function writefile(file,val),fid=fopen(file,'w','n','UTF-8');assert(fid>=0);c=onCleanup(@()fclose(fid));fprintf(fid,'%s',val);end