function receipt=execute_bellhop_1500m_d401(project_root,run_dir)
%EXECUTE_BELLHOP_1500M_D401 Exactly one new C401 slot, then bounded comparison.
timer=tic; addpath(fullfile(project_root,'scripts','validation','support')); addpath(fullfile(run_dir,'code'),'-begin');
base=load(fullfile(run_dir,'baseline201','independent_audit.mat'),'p','ref','p_native');
plan=load(fullfile(run_dir,'plan.mat')); cfg=plan.slots.cfg; oldcfg=cfg; oldcfg.beam_count=201;
assert(isequaln(oldcfg,plan.frozen.cfg),'Only beam count may change');
assert(isequal(cfg.nominal_receiver_x_m,[0 10 50 100 200])&&cfg.frequency_hz==4000&&cfg.step_m==.5&&cfg.ssp_spacing_m==.625);
r=run_bellhop_1500m_d401_case(project_root,run_dir,160);
receipt=struct('status','PROCESS_OR_INPUT_FAILED_STOPPED','new_Bellhop_jobs',r.started_slots,'other17_slots_run',0,'P2_run',false,'beams',[201 401],'old201_recomputed',false,'exception',r.error,'figure_generated',false,'physical_pass',false,'source_width','N/A ideal omni 3D point; pressure_1m','ELF_SHA',plan.slots.binary_sha256);
root=fullfile(run_dir,'jobs','01_D_4000_C401','case');
if ~isfile(fullfile(fileparts(root),'component.mat'))
 writejson(fullfile(run_dir,'comparison_receipt.json'),receipt); report=['# D4000/C401停止报告' newline newline '输入或执行失败；无有效新分量。仅尝试唯一401槽，没有重试。' newline r.error]; writefile(fullfile(run_dir,'chinese_report.md'),report); return
end
s=load(fullfile(fileparts(root),'component.mat')); checks=readtable(fullfile(run_dir,'C_checks.csv'),'Delimiter',',','TextType','string');
raw=read_bellhop_rotated_field(root); fid=fopen([root '.shd'],'rb'); assert(fid>=0); native=complex(zeros(5,1));
for j=1:5, fseek(fid,(raw.range_record+j)*raw.record_bytes,'bof'); a=fread(fid,2,'float32'); native(j)=complex(a(1),a(2)); end; fclose(fid);
assert(isequal(native,raw.pressure)&&isequal(-conj(native),s.p(:)),'Independent SHD/MAT mismatch');
x=[0 10 50 100 200]'; L=hypot(1484,x); ref=exp(1i*2*pi*4000/1500*L)./L; E=abs(ref); p201=base.p(:); p401=-conj(native);
assert(max(abs(base.ref(:)-ref))<1e-14,'Frozen analytic reference mismatch');
oldraw=read_bellhop_rotated_field(fullfile(run_dir,'baseline201','case')); assert(isequal(oldraw.pressure(:),base.p_native(:)),'Frozen 201 SHD/MAT mismatch');
tcsv=readtable(fullfile(run_dir,'baseline201','direct_field_independent.csv'),'Delimiter',','); assert(max(abs(tcsv.Bellhop_p_real+1i*tcsv.Bellhop_p_imag-p201))<1e-17);
m201=measure(p201,ref,E); m401=measure(p401,ref,E);
t=readtable(fullfile(fileparts(root),'field.csv'),'Delimiter',','); assert(max(abs(t.p_real+1i*t.p_imag-p401))<1e-17);
metrics=table([201;401],[m201.L2;m401.L2],[m201.maxDE;m401.maxDE],[m201.maxphase;m401.maxphase],[m201.support;m401.support],[m201.L2<=.001;m401.L2<=.001],[m201.maxDE<=.002;m401.maxDE<=.002],[m201.maxphase<=.01;m401.maxphase<=.01],'VariableNames',{'beams','complex_L2','max_complex_over_E','max_phase_rad','phase_support','L2_pass','max_complex_pass','phase_pass'}); writetable(metrics,fullfile(run_dir,'metrics_201_401.csv'));
fields=table(repmat(x,2,1),[201*ones(5,1);401*ones(5,1)],[real(p201);real(p401)],[imag(p201);imag(p401)],repmat(real(ref),2,1),repmat(imag(ref),2,1),[100*(abs(p201./ref)-1);100*(abs(p401./ref)-1)],[angle(p201.*conj(ref));angle(p401.*conj(ref))],[abs(p201-ref)./E;abs(p401-ref)./E],'VariableNames',{'x_m','beams','p_real','p_imag','ref_real','ref_imag','amplitude_error_percent','phase_error_rad','complex_error_over_E'}); writetable(fields,fullfile(run_dir,'fields_201_401.csv'));
diff=p401-p201; difference=table(x,real(diff),imag(diff),abs(diff),abs(diff)./E,angle(p401.*conj(p201)),'VariableNames',{'x_m','delta_p_real','delta_p_imag','abs_delta_p','abs_delta_over_E','phase401_minus201_rad'}); writetable(difference,fullfile(run_dir,'complex_difference_401_minus201.csv'));
fid=fopen(fullfile(run_dir,'baseline201','case.r25diag'),'r'); fgetl(fid); d201=fscanf(fid,'%f',[10 inf]).'; fclose(fid);
fid=fopen([root '.r25diag'],'r'); fgetl(fid); d401=fscanf(fid,'%f',[10 inf]).'; fclose(fid);
assert(size(d201,1)==201&&size(d401,1)==401&&max(abs(d201(:,1)-d401(1:2:end,1)))<1e-12,'Nested fan identity');
gridtable=table();
for k=1:2
 if k==1, d=d201; p=p201; nb=201; else, d=d401; p=p401; nb=401; end
 alpha=deg2rad(d(:,1)); da=(alpha(end)-alpha(1))/(nb-1); th=atan2(x,1484); pos=(th-alpha(1))/da; nearest=round(pos); mask=abs(pos-nearest)<1e-10; pos(mask)=nearest(mask); frac=pos-floor(pos); prediction=-(2*pi*4000/1500)*L/2*da^2.*frac.*(1-frac);
 gt=table(x,nb*ones(5,1),rad2deg(th),frac,repmat(da,5,1),prediction,angle(p.*conj(ref)),'VariableNames',{'x_m','beams','required_angle_deg','angular_cell_fraction','delta_angle_rad','leading_phase_prediction_rad','measured_phase_error_rad'}); gridtable=[gridtable;gt]; %#ok<AGROW>
end
writetable(gridtable,fullfile(run_dir,'angular_grid_comparison.csv'));
pass=all(checks.pass)&&m401.L2<=.001&&m401.maxDE<=.002&&m401.maxphase<=.01;
receipt.status='D401_NUMERICAL_FAIL_STOPPED'; if pass, receipt.status='D401_SINGLE_SLOT_GO_P1_INCOMPLETE'; end
receipt.physical_pass=pass; receipt.native_exit_code=s.proc.exit_code; receipt.native_elapsed_s=s.proc.elapsed_s; receipt.process=s.proc;
receipt.old201=m201; receipt.new401=m401; receipt.failed_numeric_gate_count=nnz(~[m401.L2<=.001 m401.maxDE<=.002 m401.maxphase<=.01]); receipt.failed_all_saved_checks=nnz(~checks.pass);
receipt.pair_complex_L2=norm(diff)/norm(p201); receipt.pair_max_difference_over_E=max(abs(diff)./E); receipt.pair_is_diagnostic_not_a_frozen_acceptance_gate=true;
receipt.actual_field_coordinates='REAL8 bitwise and roundtrip checks retained'; receipt.independent_SHD_MAT_CSV_audit='PASS'; receipt.next_allowed_step='manual review only; no 801, other P1 slots or P2';
save(fullfile(run_dir,'comparison.mat'),'p201','p401','ref','E','x','fields','difference','gridtable','checks','metrics','receipt','native','cfg');
if pass
 fig=figure('Visible','off','Color','w','Position',[50 50 1100 720]); tl=tiledlayout(fig,2,2,'Padding','compact','TileSpacing','compact');
 vals={real(p201-ref)./E,real(p401-ref)./E;imag(p201-ref)./E,imag(p401-ref)./E;abs(p201-ref)./E,abs(p401-ref)./E;angle(p201.*conj(ref)),angle(p401.*conj(ref))}; labels={'Re(complex residual)/E','Im(complex residual)/E','|complex residual|/E','Phase error (rad)'};
 for j=1:4, ax=nexttile(tl); h=plot(ax,x,vals{j,1},'k--o',x,vals{j,2},'b-s'); set(h(1),'Tag',sprintf('201_%d',j)); set(h(2),'Tag',sprintf('401_%d',j)); xlabel(ax,'Independent link offset x (m)'); ylabel(ax,labels{j}); grid(ax,'on'); legend(ax,{'201 saved','401 new'},'Location','best'); if j==3, yline(ax,.002,'r:','Frozen max gate'); end; end
 title(tl,'U0 direct only, 4000 Hz: saved 201 / new 401 beams'); drawnow; savefig(fig,fullfile(run_dir,'F01_complex_error_201_401.fig')); exportgraphics(fig,fullfile(run_dir,'F01_complex_error_201_401.png'),'Resolution',150); close(fig);
 fig=openfig(fullfile(run_dir,'F01_complex_error_201_401.fig'),'invisible'); for j=1:4, for nb=[201 401], obj=findobj(fig,'Tag',sprintf('%d_%d',nb,j)); col=1+(nb==401); assert(isscalar(obj)&&isequal(obj.XData,x')&&isequal(obj.YData,vals{j,col}'),'FIG reconstruction mismatch'); end; end; close(fig);
 receipt.figure_generated=true; receipt.FIG_reopen_audit='8 curves PASS';
end
receipt.MATLAB_workflow_elapsed_s=toc(timer); receipt.matlab_version=version;
writejson(fullfile(run_dir,'comparison_receipt.json'),receipt);
report=string('# D/4000Hz/C/401束唯一最小复验')+newline+newline+'日期：2026-11-18。结论：'+receipt.status+'。'+newline+newline+ ...
 '本轮只新增401束直达C槽；201束复用原SHD/MAT/CSV，不重算。原物理、5位置、扇区、0.5m步长、0.625m Q间距、正式REAL8 C ELF、点源1m参考及-conj相位桥保持。源宽N/A。其余17槽/P2/801/FFP/Munk未运行。'+newline+newline+ ...
 '## 冻结门与实际结果'+newline+newline+'| 束数 | 复L2（门0.001） | max复差/E（门0.002） | max相位/rad（门0.01） |'+newline+'|---:|---:|---:|---:|'+newline;
report=report+sprintf('|201|%.12g|%.12g|%.12g|\n|401|%.12g|%.12g|%.12g|\n',m201.L2,m201.maxDE,m201.maxphase,m401.L2,m401.maxDE,m401.maxphase);
report=report+newline+sprintf('401三项冻结复场门失败%d项；原生exit=%d，native=%.6fs。201/401配对复差L2=%.12g、max差/E=%.12g，仅作诊断，不另设或改变验收门。\n',receipt.failed_numeric_gate_count,s.proc.exit_code,s.proc.elapsed_s,receipt.pair_complex_L2,receipt.pair_max_difference_over_E);
report=report+newline+'## 五点幅相对照'+newline+newline+'|x/m|201幅度误差/%|401幅度误差/%|201相位/rad|401相位/rad|'+newline+'|---:|---:|---:|---:|---:|'+newline;
for j=1:5, report=report+sprintf('|%g|%.10g|%.10g|%.10g|%.10g|\n',x(j),fields.amplitude_error_percent(j),fields.amplitude_error_percent(j+5),fields.phase_error_rad(j),fields.phase_error_rad(j+5)); end
report=report+newline+'## 角网格与证据范围'+newline+newline+'angular_grid_comparison.csv记录两档每点角单元分数、实测相位和固定hat近似的二阶投影预测。401角网格嵌套原201，原角度逐点复核≤1e-12度。误差包络近似随角间隔平方缩小，但每点单元分数同时变化，不能要求逐点恰好四倍。单个加密槽不能替代完整P1或普适收敛。'+newline+newline+ ...
 '独立SHD读取/MAT/CSV及实际REAL8坐标、401束终态几何已审计。原始SHD/ENV/SSP/R25D/PRT、真实exit、资源采样、comparison.mat、metrics/fields/difference/grid CSV全部保留。本地总监督器统一覆盖启动、回收、SHA、报告和导航封存；最终预算见supervisor_receipt.json。C模式不产生ARR/RAY，没有额外A/R。'+newline+newline;
if pass, report=report+'新401槽三门通过；图F01保留PNG和可编辑FIG，8条曲线重开复核通过。本结果只支持本单槽，不授权P2。'+newline; else, report=report+'数值门失败，立即停止新计算，不生成通过后图件，不修公式、不801、不自动重试。失败原始数据保留。'+newline; end
report=report+newline+'没有修改官方Bellhop、PE或公共接口，不提交推送。原201失败及383.24s历史超限保持。停止等待人工审查。'+newline;
writefile(fullfile(run_dir,'chinese_report.md'),char(report));
fprintf('D401_RESULT status=%s L2=%.12g maxDE=%.12g phase=%.12g native=%.6fs failed%d\n',receipt.status,m401.L2,m401.maxDE,m401.maxphase,s.proc.elapsed_s,receipt.failed_numeric_gate_count);
end
function m=measure(p,ref,E)
support=abs(p)>=.01*E&abs(ref)>=.01*E; m=struct('L2',norm(p-ref)/norm(ref),'maxDE',max(abs(p-ref)./E),'maxphase',max(abs(angle(p(support).*conj(ref(support))))),'support',nnz(support));
end
function writejson(file,val), writefile(file,jsonencode(val,PrettyPrint=true)); end
function writefile(file,val), fid=fopen(file,'w','n','UTF-8'); assert(fid>=0); c=onCleanup(@()fclose(fid)); fprintf(fid,'%s',val); end