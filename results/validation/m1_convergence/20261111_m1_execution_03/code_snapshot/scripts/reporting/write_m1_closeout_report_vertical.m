function output=write_m1_closeout_report_vertical(run_dir)
%WRITE_M1_CLOSEOUT_REPORT_VERTICAL Report saved evidence; never solves fields.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
old_threads=maxNumCompThreads(1); restore_threads=onCleanup(@()maxNumCompThreads(old_threads)); %#ok<NASGU>
assert(startsWith(char(java.io.File(run_dir).getCanonicalPath()), ...
    [fullfile(root,'results','validation','m1_convergence') filesep]));
date=char(datetime('now','TimeZone','Asia/Shanghai','Format','yyyy-MM-dd'));
report_path=fullfile(root,'reports',['m1_static_convergence_qualification_' date '.md']);
assert(~isfile(report_path),'Preserve existing reports; choose a new report version explicitly.');
s1=read_checkpoint('stage1_audit'); s2=read_checkpoint('stage2_audit');
status='M1_INCOMPLETE';
if strcmp(s1.status,'STAGE1_QUALIFIED')
    if strcmp(s2.status,'SPONGE_QUALIFIED'), status='M1_COMPLETE_SPONGE_QUALIFIED';
    elseif strcmp(s2.status,'AUDIT_FAILED') && isfield(s2,'failed_or_unrun_checks') && ...
            ~isempty(s2.failed_or_unrun_checks) && ...
            all(startsWith({s2.failed_or_unrun_checks.name},'sponge_'))
        status='M1_COMPLETE_NO_SPONGE_ONLY';
    end
end
ids={'P4','P6','P8','Z4','Z8','A-C4','A-C6','A-C8','A-W4','A-W6','A-W8','A-G4','A-G6','A-G8','A-R4','SP4','SP6','SP8'};
budget=jsondecode(fileread(fullfile(run_dir,'budget.json')));
partial_summary=struct([]);
partial_file=fullfile(run_dir,'P8','partial.mat');
if isfile(partial_file)
    p=load(partial_file,'trace','cfg'); paths={'direct','incident','reflect'}; lengths=[97 100 3];
    for j=1:3
        tr=p.trace.(paths{j}); expected=ceil(lengths(j)/(p.cfg.stepz_lamb*1500/p.cfg.f_hz))+1;
        if isempty(tr), n=0; max_edge=NaN; max_bd=NaN;
        else
            n=size(tr,1); max_edge=max(tr(:,3)); max_bd=max(tr(:,4));
            tab=array2table(tr,'VariableNames',{'step','z_m','edge5_fraction','boundary_db','peak','total_power'});
            writetable(tab,fullfile(run_dir,'P8',['partial_' paths{j} '_every_step.csv']));
        end
        row=struct('path',paths{j},'observed_planes',n,'expected_planes',expected,'edge5_max_observed',max_edge,'boundary_max_db_observed',max_bd);
        if isempty(partial_summary), partial_summary=row; else, partial_summary(end+1)=row; end %#ok<AGROW>
    end
    jsonfile(fullfile(run_dir,'P8','partial_trace_summary.json'),partial_summary);
    clear p
end
if ~strcmp(status,'M1_INCOMPLETE')
    chosen='P4'; if strcmp(status,'M1_COMPLETE_SPONGE_QUALIFIED'), chosen='SP4'; end
    p=load(fullfile(run_dir,chosen,'payload.mat'),'payload'); paramsV=p.payload.params;
    paramsV.validation_step_observer=[];
    qualification_meta=struct('status',status,'frequencies_hz',[4000 6000 8000], ...
        'receiver_radius_m',2,'phase_contract_version','2.0.0','code_binding','per-job checkpoint source SHA', ...
        'scope','Uniform c=1500, Gaussian sigma=.3 continuous unit peak, Tx100/Rx3, flat pressure release; no broadband/CIR/rough surface claim.');
    save(fullfile(run_dir,'recommended_params.mat'),'paramsV','qualification_meta');
    jsonfile(fullfile(run_dir,'recommended_config.json'),struct('status',status, ...
        'nx',984,'ny',984,'xw_m',192.1875,'yw_m',192.1875,'dx_m',.1953125, ...
        'stepz_lamb',.5,'sponge_ratio',paramsV.sponge_ratio, ...
        'alpha_max_np_per_m',paramsV.alpha_max_np_per_m,'params_mat','recommended_params.mat'));
    clear p paramsV
end
fid=fopen(report_path,'w','n','UTF-8'); cleaner=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# M1：100 m静态信道收敛与原生产吸收层资格收尾\n\n日期：%s（Asia/Shanghai）。状态：**%s**。\n\n',date,status);
fprintf(fid,'本报告读取本轮独立审计、逐步轨迹、MAT/CSV与作业回执，不将历史窗口或随机海面PASS升级为本轮资格。基线compare/c90d0ca，M0生产969127f与phase contract v2保留；未合并main、推送、使用SSH/GPU或外部求解器。\n\n');
fprintf(fid,'证据根目录：`%s`。原连接失败与第一次日志不完整作业分别保留在`20261111_m1_execution_preflight_01`和`20261111_m1_execution_02`，没有覆盖或改写历史检查点。\n\n',strrep(run_dir,'\','/'));
fprintf(fid,'## 固定模型、合同及比较\n\nTx=(0,0,100)m、Rx=(0,0,3)m、c=1500m/s；Gaussian sigma=.3m、连续单位峰值；平面eta=0、R0=-1。直达97m、Tx到海面100m、完整反射103m；不包含海底边界，不与Tx99m Bellhop安装混用。全部新传播为现有本机MATLAB MCP R2025b Update 6、PID28836、CPU double、单线程串行、无bubble/Doppler。\n\n');
fprintf(fid,'原生U(L)=exp(iL(kz-k0))不变。物理分量Pd=d exp(ik0*97)、Pr=r exp(ik0*103)；direct-DSP为conj(d)、exp(-ik0*6)conj(r)，总场H_DSP=conj(Pd+Pr) exp(ik0*97)。native/physical/DSP三表示及direct/reflect/total均留在每作业payload.mat，v2标签保留。\n\n');
fprintf(fid,'共同坐标与rho<=.5/1/2m预先固定；扩窗直接匹配同dx格点，细网格每隔一点取公共位置，不插值、归一化或逐路径移峰。复L2以对应参考分量为分母，总场统一改用参考直达L2；另外记录直达峰值标度最大绝对复差。幅度/相位只排除低于参考区域峰值1e-12的点，排除数逐行报告。总物理场先恢复两条路径载波再相加。\n\n');
fprintf(fid,'验收：PE对同格AS复误差<=1e-10、分量/标量闭合<=1e-12；窗口/网格/步长/吸收层中心复L2<=1e-3、最大幅度差<=.05dB、相位RMS<=.01rad。每传播步外围5%%能量严格<1e-4、离散首末行列相对幅度<=-40dB。AS参考自身在97/100/103m独立检查边缘，不能以参考互相接近授予资格。\n\n');
fprintf(fid,'## 18项冻结作业\n\n|作业|频率Hz|W m / N / dx m|step/lambda|状态|耗时s|最低可用RAM GiB|显式stdout|\n|---|---:|---|---:|---|---:|---:|---|\n');
completed=0; wrapper_extra=0; min_ram=Inf; max_mem=0;
matrix=struct('job_id',{},'frequency_hz',{},'width_m',{},'n',{},'dx_m',{}, ...
    'stepz_lamb',{},'status',{},'elapsed_s',{},'has_complete_payload',{});
for j=1:numel(ids)
    r=read_checkpoint(ids{j});
    if ~isfield(r,'config')
        f=str2double(ids{j}(end))*1000; step=.5; if startsWith(ids{j},'Z'), step=.25; end
        matrix(end+1)=struct('job_id',ids{j},'frequency_hz',f,'width_m',192.1875, ...
            'n',984,'dx_m',.1953125,'stepz_lamb',step,'status','NOT_DISPATCHED', ...
            'elapsed_s',NaN,'has_complete_payload',false); %#ok<AGROW>
        fprintf(fid,'|%s|%.0f|192.1875 / 984 / .1953125|%.2g|NOT_DISPATCHED|NA|NA|无|\n',ids{j},f,step); continue
    end
    c=r.config; captured=isfield(r,'stdout_capture');
    matrix(end+1)=struct('job_id',ids{j},'frequency_hz',c.f_hz,'width_m',c.width_m, ...
        'n',c.n,'dx_m',c.dx_m,'stepz_lamb',c.stepz_lamb,'status',r.status, ...
        'elapsed_s',r.elapsed_s,'has_complete_payload',isfile(fullfile(run_dir,ids{j},'payload.mat'))); %#ok<AGROW>
    completed=completed+1; min_ram=min(min_ram,r.min_available_bytes); max_mem=max(max_mem,r.max_matlab_memory_bytes);
    if isfield(r,'wrapper_elapsed_s'), wrapper_extra=wrapper_extra+max(0,r.wrapper_elapsed_s-r.elapsed_s); end
    fprintf(fid,'|%s|%.0f|%.9g / %d / %.9g|%.2g|%s|%.3f|%.3f|%d|\n',ids{j},c.f_hz,c.width_m,c.n,c.dx_m,c.stepz_lamb,r.status,r.elapsed_s,r.min_available_bytes/2^30,captured);
end
writetable(struct2table(matrix),fullfile(run_dir,'matrix_status.csv'));
fprintf(fid,'\nAS表中的step字段仅为统一配置占位；AS实际调用独立一步参考，不执行该步长的生产march。作业记录不等同资格，必须结合独立审计。每job目录含完整配置/参数、source SHA、payload.mat、fields_ieee754.csv、stdout_captured.log、stderr.log、state.json与checkpoint.json。\n\n');
fprintf(fid,'## 窗口、网格、步长与吸收层结论\n\n');
if strcmp(s1.status,'STAGE1_QUALIFIED')
    fprintf(fid,'192.1875m/984²/dx=.1953125m、step=.5lambda无吸收层候选在冻结4/6/8kHz锚点通过窗口、横向网格与边界资格；4/8kHz半步长实测通过。6kHz .25lambda未直接运行，适用论据为均匀无吸收、平面常系数算子的U(a)U(b)=U(a+b)解析恒等性质及6kHz独立同格AS，不写成直接实测通过。未声称1968²生产PE已运行、联合无限域误差界、连续频带或任意环境资格。\n\n');
else
    fprintf(fid,'第一阶段尚未全通过，窗口、横向网格、步长或边界未形成完整资格。不得推荐生产配置或建议进入M2。缺失/失败以stage1检查点和下表为准。\n\n');
    if isfile(fullfile(run_dir,'partial_closeout_audit','checkpoint.json'))
        fprintf(fid,'已完成的限定证据：三频独立AS扩窗/细网格及更大窗口参考对照通过；4/6kHz完整生产PE对同格AS通过。8kHz生产PE因RAM门中止，Z4/Z8未派发，故这些部分通过不构成一套完整生产资格。\n\n');
    end
end
fprintf(fid,'原生产吸收层真实定义：sponge_ratio=.12表示每侧厚度.12W，起点|x|=.38W；alpha_x/alpha_y在层内按三次函数增长，各方向最大.15Np/m，alpha_xy=alpha_x+alpha_y，角点可接近.30Np/m。W=192.1875m时厚度23.0625m、起点73.03125m；现有分裂步中屏乘法实衰减为exp(-alpha_xy ds)。SP仅恢复.12/.15，源、几何、网格、步长与比较支持不变；未自动调参。\n\n');
if strcmp(status,'M1_COMPLETE_SPONGE_QUALIFIED')
    fprintf(fid,'三频原吸收层对无吸收PE及同格AS的有效圆盘复场均通过门。推荐显式配置W_x=W_y=192.1875m、nx=ny=984、dx=dy=.1953125m、stepz_lamb=.5、sponge_ratio=.12、alpha_max_np_per_m=.15；同网格无吸收方案也在本轮限定范围内合格。扩窗许可validation_allow_extended_window=true需显式设置，不改公共默认值。\n\n');
elseif strcmp(status,'M1_COMPLETE_NO_SPONGE_ONLY')
    fprintf(fid,'原吸收层至少一个有效圆盘影响门不合格，保留无吸收资格。唯一推荐为同窗口/网格/步长，sponge_ratio=0、alpha_max_np_per_m=0，并显式validation_allow_extended_sponge_ratio=true和validation_allow_extended_window=true；不推荐原.12/.15吸收层，不寻找新参数。\n\n');
else
    fprintf(fid,'当前不能发布完整生产推荐配置。吸收层未运行或未完成的项不作PASS/FAIL推断；状态M1_INCOMPLETE。\n\n');
end
fprintf(fid,'## 各频点最差误差与位置\n\n');
tables={};
for folder={'stage1_audit','stage2_audit','partial_closeout_audit'}
    file=fullfile(run_dir,folder{1},'comparisons.csv'); if isfile(file), tables{end+1}=readtable(file,'TextType','string'); end %#ok<AGROW>
end
if isempty(tables)
    fprintf(fid,'完整比较表尚未生成；最差复L2、幅度、相位与坐标均NA。\n\n');
else
    comparisons=vertcat(tables{:}); writetable(comparisons,fullfile(run_dir,'all_comparisons.csv'));
    fprintf(fid,'表中L2位置为该区域内最大直达标度绝对复差的位置；幅度和相位位置分别为对应点误差最大处。每项最差可来自不同对照/路径/半径，不能合并成单一工况。\n\n|f Hz|指标|最差值|对照|路径|rho m|x,y m|\n|---:|---|---:|---|---|---:|---|\n');
    variables={'complex_l2','max_amplitude_db','phase_rms_rad'};
    for f=[4000 6000 8000]
        subset=comparisons(comparisons.frequency_hz==f,:);
        for v=1:3
            name=variables{v}; [val,index]=max(subset.(name)); r=subset(index,:);
            if v==1, xx=r.complex_x_m; yy=r.complex_y_m;
            elseif v==2, xx=r.amplitude_x_m; yy=r.amplitude_y_m;
            else, xx=r.phase_x_m; yy=r.phase_y_m; end
            fprintf(fid,'|%d|%s|%.9g|%s:%s/%s|%s|%g|%g,%g|\n',f,name,val,r.kind,r.candidate,r.reference,r.component,r.radius_m,xx,yy);
        end
    end
    fprintf(fid,'\n直达、反射、总场的每个对照及三个支持的全部指标、排除点与位置见`all_comparisons.csv`，不是只报告上表最大值。\n\n');
    if strcmp(status,'M1_INCOMPLETE')
        fprintf(fid,'上述是已完成对照的最差值；8kHz行仅含AS空间对照，**8kHz生产PE对AS误差为NA**。Z4/Z8步长对照和SP4/SP6/SP8吸收层影响全部NA，不按0或通过处理。\n\n');
    end
end
fprintf(fid,'## 每一步边缘检查与独立参考边界\n\n|作业/路径|实际检查平面数|外围5%%最大比例|所在步/深度m|边界最大dB|所在步/深度m|\n|---|---:|---:|---|---:|---|\n');
for id=ids
    for segment={'direct','incident','reflect'}
        file=fullfile(run_dir,id{1},[segment{1} '_every_step.csv']); if ~isfile(file), continue; end
        tr=readtable(file); [ev,ei]=max(tr.edge5_fraction); [bv,bi]=max(tr.boundary_db);
        fprintf(fid,'|%s/%s|%d|%.9g|%d / %.9g|%.9g|%d / %.9g|\n',id{1},segment{1},height(tr),ev,tr.step(ei),tr.z_m(ei),bv,tr.step(bi),tr.z_m(bi));
    end
end
if ~isempty(partial_summary)
    fprintf(fid,'\nP8以下为资源中止前的**部分轨迹**，不能作为完整边缘资格：\n\n|路径|已观测/要求平面数|已观测外围5%%最大比例|已观测边界最大dB|\n|---|---:|---:|---:|\n');
    for j=1:numel(partial_summary)
        r=partial_summary(j); fprintf(fid,'|%s|%d / %d|%.9g|%.9g|\n',r.path,r.observed_planes,r.expected_planes,r.edge5_max_observed,r.boundary_max_db_observed);
    end
    fprintf(fid,'\nP8的partial.mat与partial_*_every_step.csv保留原始已观测步；未观测传播步明确缺失。\n\n');
end
fprintf(fid,'\n每段含初场和实际每一步；独立审计要求步号严格等于0:n，不以33/65个平面替代。端平面的2D原始数据另重构边缘指标，逐步CSV与MAT轨迹复核。\n\n|AS组|距离m|外围5%%比例|边界dB|\n|---|---:|---:|---:|\n');
for id=ids(startsWith(ids,'A-'))
    file=fullfile(run_dir,id{1},'edge_metrics.csv'); if ~isfile(file), continue; end
    tr=readtable(file); for k=1:height(tr), fprintf(fid,'|%s|%g|%.9g|%.9g|\n',id{1},tr.distance_m(k),tr.edge5_fraction(k),tr.boundary_db(k)); end
end
fprintf(fid,'\n## 执行、资源与保留失败\n\n正式作业已保存%d/18。预算账当前%.3fs/14400s（含原run02控制/P4与保守前检15s）；显式捕获封装尚需封账的额外已测时间%.3fs。单作业门3600s，P4/SP4各900s先导门；各作业启动前以估计PE512MiB、1310² AS512MiB、大AS1GiB检查至少30%%初始可用RAM余量、大AS至少1.5GiB可用，磁盘至少2GiB。逐步每16步检查资源/时间并保存真实PID心跳。最低实测可用RAM %.6gGiB，MATLAB memory接口最大采样占用%.6gGiB；后者是采样的MATLAB分配量，不冒充物理工作集峰值。OS独立工作集/进程生涯峰值见external_resource_samples.csv，进程峰值包含本会话此前用量，不能冒充逐作业独立峰值。\n\n',completed,budget.used_s,wrapper_extra,min_ram/2^30,max_mem/2^30);
fprintf(fid,'最初连接失败未发起求解。run02首次P4正常完成且边缘门通过，用时320.7669698s，但MCP外层300s超时，diary输出为空；数值MAT和状态保留，不能把空日志当完整stdout。先确认原作业终态，再用显式evalc捕获封装在新run03重跑，未重复启动仍运行的任务。run03的stdout_captured.log保留实际控制台输出，checkpoint_solver.json保存封装前回执，最终checkpoint绑定完整日志SHA。所有重跑与控制均计入累计预算，失败历史未抹除。\n\n');
if isfile(fullfile(run_dir,'resource_stop.json'))
    stop=jsondecode(fileread(fullfile(run_dir,'resource_stop.json')));
    fprintf(fid,'**最终停止原因是RAM，不是预算耗尽或物理误差门失败。** P8启动可用RAM=%.0fB（%.9gGiB），要求保留=%.0fB（%.9gGiB），最低观测=%.0fB（%.9gGiB）；在入射段第%d/%d步通过验证回调抛异常停止，用时%.9gs。Z4、Z8与全部SP未派发，没有降低资源门、杀其他进程或自动重跑。本轮停在M1，需新的继续指令和重新核验后的资源条件才能恢复未完成求解。\n\n',stop.available_before_bytes,stop.available_before_bytes/2^30,stop.required_reserve_bytes,stop.required_reserve_bytes/2^30,stop.minimum_observed_available_bytes,stop.minimum_observed_available_bytes/2^30,stop.last_completed_observation.step,stop.last_completed_observation.n_steps,stop.elapsed_s);
end
fprintf(fid,'第一次reference_audit失败也保留：空结构数组追加比较行导致MATLAB:heterogeneousStrucAssignment；独立重构把k0^2-kx^2-ky^2重结合为k0^2-(kx^2+ky^2)，使噪声底边界dB出现最多约.02418dB的重构差，虽然参考自身物理边缘门全部通过。保留failed_audit_code原文后，修正空表初始化并按冻结表达式的从左至右浮点顺序独立重构，在新reference_audit_v2目录重新审计已有MAT，未重跑传播或放宽任何容差。最终资格使用v2审计；原失败不改写。\n\n');
fprintf(fid,'本轮验证性代码只增加默认关闭的CPU单频validation_step_observer（每段初场及每一步只读回调），不改变传播公式、源、边界或相位合同；64²控制比较开启/关闭关键字段逐元素差0并验证direct_only。另以既有M0保存宽带H回放消费者回归，无新宽带频率求解或正式通信实验。AS独立重构不调用生产march或AS helper；MAT与IEEE-754 CSV逐位复核，源码/产物SHA、全步覆盖、标量/分量和v2合同均由独立审计执行。\n\n');
fprintf(fid,'|阶段审计|状态|失败数|耗时s|\n|---|---|---:|---:|\n');
for mode={'reference_audit','partial_closeout_audit','stage1_audit','stage2_audit'}
    r=read_checkpoint(mode{1});
    if isfield(r,'elapsed_s'), fprintf(fid,'|%s|%s|%d|%.3f|\n',mode{1},r.status,r.failed_checks,r.elapsed_s);
    else, fprintf(fid,'|%s|NOT_RUN|NA|NA|\n',mode{1}); end
end
fprintf(fid,'\n失败/未运行项逐条见各审计checkpoint与checks.csv；不能以退出码0或文件存在替代数值PASS。第一阶段若未通过，不执行第二阶段；完整边界配置未合格则不关闭M1。\n\n');
fprintf(fid,'## 复现、版本和后续限制\n\n仅在已有MATLAB MCP中：先setup_vertical_project，使用新的绝对run目录；以run_m1_logged_job_vertical(''P4'',run_dir)先导，确认终态/日志/<=900s后运行run_m1_convergence_stage_vertical(run_dir,''references'')，再audit_m1_convergence_vertical(run_dir,''reference_audit'')。仅参考门通过才运行remaining_pe分发段及stage1_audit；仅STAGE1_QUALIFIED才派发SP4，<=900s且正常后派发sponge_remaining并运行stage2_audit。每个入口拒绝覆盖已有job/审计目录，超时核PID、state、dispatch/checkpoint后等待原作业结束，禁止盲目重派。\n\n');
fprintf(fid,'每作业绑定实际源码SHA，最终源快照/implementation.patch和SHA清单见run目录；输入/输出最终封账见checkpoint_final.json与artifacts_sha256.csv。公共字段h_direct/h_reflect/h_total/H_f/f_axis/idx_f_ref保留；新观察字段默认关闭。PROJECT_CONTEXT、vertical_comm_guide、src/scripts索引和进度JSON同步更新，看板由既有MCP正式渲染并核引用。\n\n');
if ~strcmp(status,'M1_INCOMPLETE')
    fprintf(fid,'完整可复用参数存于recommended_params.mat，摘要见recommended_config.json。新会话连接通过后可load该MAT中的paramsV；f0默认4kHz，改成6000或8000仍仅在本轮锚点资格内。实际调用为output=vertical_channel_model(paramsV)，不依赖当前目录。该配置保留本轮save_mode/nout等输出设置，不静默改变公共默认值。\n\n');
end
fprintf(fid,'本轮最多认证指定Gaussian源、均匀介质、平面一次海面反射、4/6/8kHz离散锚点与rho<=2m。本任务到M1结束；即便M1完成，也只说明未来可另行安排M2，不在本轮启动宽带CIR、随机海面、正式通信、1500m深海或正式S5/S6。\n');
output=struct('status',status,'report_path',report_path,'run_dir',run_dir,'jobs_present',completed, ...
    'budget_used_before_final_packaging_s',budget.used_s,'wrapper_extra_s',wrapper_extra);
jsonfile(fullfile(run_dir,'report_receipt.json'),output);
fprintf('M1_REPORT_RECEIPT %s\n',jsonencode(output));
    function r=read_checkpoint(name)
        file=fullfile(run_dir,name,'checkpoint.json');
        if strcmp(name,'reference_audit') && isfile(fullfile(run_dir,'reference_audit_v2','checkpoint.json'))
            file=fullfile(run_dir,'reference_audit_v2','checkpoint.json');
        end
        if isfile(file), r=jsondecode(fileread(file)); else, r=struct('status','NOT_RUN'); end
    end
end
function jsonfile(file,s)
fid=fopen(file,'w','n','UTF-8'); fprintf(fid,'%s',jsonencode(s,PrettyPrint=true)); fclose(fid);
end
