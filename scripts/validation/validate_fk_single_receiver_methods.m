function result = validate_fk_single_receiver_methods(mode)
%VALIDATE_FK_SINGLE_RECEIVER_METHODS Paired accuracy/cost experiment, no PE edits.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
outdir=fullfile(root,'results','validation','fk_single_receiver_methods');
if ~isfolder(outdir), mkdir(outdir); end
if nargin<1, mode='full'; end
if strcmp(mode,'render')
    loaded=load(fullfile(outdir,'evaluation.mat'),'result'); result=loaded.result;
    local_report(result,root);
    previous=load(fullfile(outdir,'pre_hardening','evaluation.mat'),'result');
    local_hardening_report(result,previous.result,root); return
end
assert(strcmp(mode,'full'),'Use full or render.');
input=fullfile(root,'results','validation','spectral_fast_kirchhoff_multifrequency_consistency','consistency.mat');
loaded=load(input,'s'); s=loaded.s; assert(s.pass);
paths={input,fullfile(root,'scripts','validation','support','solve_full_kirchhoff_2d.m'), ...
    fullfile(root,'scripts','validation','support','solve_spectral_fast_kirchhoff_2d.m'), ...
    fullfile(root,'scripts','validation','support','prepare_spectral_fast_kirchhoff_2d.m'), ...
    fullfile(root,'scripts','validation','support','solve_helmholtz_bie_halfplane_vertical.m')};
for j=1:numel(s.reference_hashes)
    p=fullfile(root,s.reference_hashes(j).path);
    assert(strcmp(local_hash(p),s.reference_hashes(j).sha256),'Frozen reference changed.');
end
before=cellfun(@local_hash,paths,'UniformOutput',false);
code_paths={mfilename('fullpath'),fullfile(root,'scripts','validation','support','fk_single_receiver_prototype.m'), ...
    fullfile(root,'scripts','validation','support','fk_single_receiver_incident.m'), ...
    fullfile(root,'scripts','validation','support','fk_single_receiver_error_metrics.m'), ...
    fullfile(root,'scripts','validation','support','fk_single_receiver_expect_error.m'), ...
    fullfile(root,'scripts','validation','support','validate_fk_single_receiver_plan_controls.m'), ...
    fullfile(root,'scripts','validation','support','diagnose_fk_single_receiver_cancellation.m'), ...
    fullfile(root,'scripts','validation','support','audit_fk_single_receiver_methods_artifacts.m')};
code_paths{1}=[code_paths{1} '.m'];
code_hashes=cellfun(@local_hash,code_paths,'UniformOutput',false);
seeds=260001:260012; profiles=cell(size(seeds));
for j=1:numel(seeds)
    profiles{j}=build_paired_pm_gaussian_realization(s.profile,seeds(j), ...
        fullfile(outdir,sprintf('coefficients_%d.csv',seeds(j))));
end
methods={'full','fastin','composite'}; frequencies=[4000 6000 8000];
rows=struct([]); guards=struct([]); timing=struct([]); setups=struct([]);
plan_controls=struct([]); cancellation=struct([]); cancellation_scan=table;
for f=frequencies
    cfg=local_cfg(s,f); plan_in=fk_single_receiver_prototype('prepare_fastin',cfg);
    plan_composite=fk_single_receiver_prototype('prepare_composite',cfg);
    controls=validate_fk_single_receiver_plan_controls(cfg,plan_in,plan_composite);
    plan_controls=[plan_controls controls]; %#ok<AGROW>
    [diagnostic,scan]=diagnose_fk_single_receiver_cancellation(cfg,plan_in);
    cancellation=[cancellation diagnostic]; %#ok<AGROW>
    cancellation_scan=[cancellation_scan;scan]; %#ok<AGROW>
    % Canonical PM, independently drawn PM, flat/offset/strong sine controls.
    for j=1:numel(seeds)+5
        c=cfg;
        if j<=numel(seeds)
            [eta,slope]=evaluate_pm_pilot_surface(profiles{j},c.surface_x_m);
            label=sprintf('PM_%d',seeds(j));
        else
            jj=j-numel(seeds); x=c.surface_x_m;
            if jj<=3
                offset=[0 -.5 .5]; eta=offset(jj)+zeros(size(x)); slope=zeros(size(x));
                label=sprintf('flat_%g',offset(jj));
            else
                amplitude=[.02 .45]; a=amplitude(jj-3);
                eta=a*cos(.47*x); slope=-a*.47*sin(.47*x); label=sprintf('sine_%g',a);
            end
        end
        c.surface_z_m=eta; c.surface_normal=[-slope ones(size(slope))]./sqrt(1+slope.^2);
        ref=fk_single_receiver_prototype('full',c);
        fast=fk_single_receiver_prototype('fastin',c,plan_in);
        composite=fk_single_receiver_prototype('composite',c,plan_composite);
        exact=fk_single_receiver_prototype('exact_composite',c);
        if j==1
            anchor=s.sparse([s.sparse.frequency_hz]==f);
            assert(abs(ref-anchor.full(s.x_m==0))/abs(ref)<1e-12);
        end
        for method=2:3
            fields=[ref fast composite]; a=fields(method);
            metrics=fk_single_receiver_error_metrics(ref,a);
            row=struct('frequency_hz',f,'case_name',label,'method',methods{method}, ...
                'full_real',real(ref),'full_imag',imag(ref),'approx_real',real(a),'approx_imag',imag(a), ...
                'max_height_m',max(abs(eta)),'max_slope',max(abs(slope)));
            names=fieldnames(metrics);
            for metric=1:numel(names), row.(names{metric})=metrics.(names{metric}); end
            rows=local_append(rows,row);
        end
        guards=local_append(guards,local_guard(f,label,'product_identity',abs(exact-ref)/max(abs(ref),1e-14),1e-12));
        if j==1 || j==numel(seeds)+5
            hi=c; hi.height_nodes=65;
            plan_in65=fk_single_receiver_prototype('prepare_fastin',hi);
            plan_composite65=fk_single_receiver_prototype('prepare_composite',hi);
            guards=local_append(guards,local_guard(f,label,'incoming_nodes33_65', ...
                abs(fk_single_receiver_prototype('fastin',hi,plan_in65)-fast)/abs(ref),1e-5));
            guards=local_append(guards,local_guard(f,label,'composite_nodes33_65', ...
                abs(fk_single_receiver_prototype('composite',hi,plan_composite65)-composite)/abs(ref),1e-5));
            hi=c; hi.rank=5; p5=fk_single_receiver_prototype('prepare_fastin',hi);
            guards=local_append(guards,local_guard(f,label,'incoming_rank3_5', ...
                abs(fk_single_receiver_prototype('fastin',hi,p5)-fast)/abs(ref),1e-3));
            fine=c; fine.surface_x_m=linspace(-60,60,8193)';
            if j==1
                [fine.surface_z_m,sl]=evaluate_pm_pilot_surface(profiles{1},fine.surface_x_m);
            else
                fine.surface_z_m=.45*cos(.47*fine.surface_x_m); sl=-.45*.47*sin(.47*fine.surface_x_m);
            end
            fine.surface_normal=[-sl ones(size(sl))]./sqrt(1+sl.^2);
            finer=fk_single_receiver_prototype('full',fine);
            guards=local_append(guards,local_guard(f,label,'quadrature4097_8193',abs(finer-ref)/abs(finer),1e-3));
        end
    end
    % Five paired timing blocks; three orders cyclically counterbalanced.
    % Six DIFFERENT PM surfaces per apply batch, shared across all methods.
    for block=1:5
        order=circshift(1:3,[0 mod(block-1,3)]);
        for method=order
            clock=tic;
            if method==1
                p=[];
            elseif method==2
                p=fk_single_receiver_prototype('prepare_fastin',cfg);
            else
                p=fk_single_receiver_prototype('prepare_composite',cfg);
            end
            setup=toc(clock);
            clock=tic; fk_single_receiver_prototype(methods{method},cfg,p); first=toc(clock);
            clock=tic;
            for j=1:6
                c=cfg; [eta,sl]=evaluate_pm_pilot_surface(profiles{j},c.surface_x_m);
                c.surface_z_m=eta; c.surface_normal=[-sl ones(size(sl))]./sqrt(1+sl.^2);
                fk_single_receiver_prototype(methods{method},c,p);
            end
            batch=toc(clock)/6;
            integral_s=NaN;
            if method==1
                [u,ux,uz]=fk_single_receiver_incident(cfg.spectral_source,cfg.surface_x_m,cfg.surface_z_m);
                fc=struct('frequency_hz',f,'c0_mps',cfg.c0_mps,'surface_x_m',cfg.surface_x_m, ...
                    'surface_z_m',cfg.surface_z_m,'surface_normal',cfg.surface_normal, ...
                    'surface_incident',u,'surface_incident_dn',cfg.surface_normal(:,1).*ux+cfg.surface_normal(:,2).*uz, ...
                    'receiver_x_m',cfg.receiver_x_m,'receiver_z_m',cfg.receiver_z_m);
                clock=tic;
                for call=1:10, solve_full_kirchhoff_2d(fc); end
                integral_s=toc(clock)/10;
            end
            row=struct('frequency_hz',f,'block',block,'method',methods{method}, ...
                'order_position',find(order==method),'setup_s',setup,'first_apply_s',first, ...
                'cold_s',setup+first,'per_realization_s',batch,'full_integral_only_s',integral_s);
            timing=local_append(timing,row);
        end
    end
    payload=whos('plan_in','plan_composite');
    setups=local_append(setups,struct('frequency_hz',f,'fastin_plan_bytes',payload(strcmp({payload.name},'plan_in')).bytes, ...
        'composite_plan_bytes',payload(strcmp({payload.name},'plan_composite')).bytes));
    fprintf('Finished %g Hz single-Rx accuracy and paired timing.\n',f);
end
after=cellfun(@local_hash,paths,'UniformOutput',false); assert(isequal(before,after));
result=struct('schema_version','1.1','version',version,'seeds',seeds,'frequencies',frequencies, ...
    'cfg',cfg,'profiles',{profiles},'accuracy',struct2table(rows),'guards',struct2table(guards), ...
    'timing',struct2table(timing),'memory',struct2table(setups), ...
    'reference_paths',{paths},'reference_hashes',{before},'references_unchanged',true, ...
    'code_paths',{code_paths},'code_hashes',{code_hashes},'plan_controls',struct2table(plan_controls), ...
    'cancellation',struct2table(cancellation),'cancellation_scan',cancellation_scan);
result.pass=all(result.accuracy.relative_error<1e-3) && all(result.accuracy.finite_metrics) && ...
    all(result.accuracy.phase_defined) && all(result.guards.pass) && all(result.plan_controls.pass) && ...
    all(result.cancellation.finite_metrics);
previous_file=fullfile(outdir,'pre_hardening','evaluation.mat');
previous=load(previous_file,'result'); previous=previous.result;
assert(isequal(result.accuracy.frequency_hz,previous.accuracy.frequency_hz) && ...
    isequal(result.accuracy.case_name,previous.accuracy.case_name) && isequal(result.accuracy.method,previous.accuracy.method));
pressure_delta=max(abs(complex(result.accuracy.approx_real,result.accuracy.approx_imag)- ...
    complex(previous.accuracy.approx_real,previous.accuracy.approx_imag)));
result.hardening=struct('previous_evaluation_sha256',local_hash(previous_file), ...
    'max_pressure_change',pressure_delta,'pressure_unchanged',pressure_delta<=1e-12, ...
    'primary_scope','Tx100 / Rx(0,3), existing 51 paired inputs at 4/6/8 kHz', ...
    'cancellation_relative_gate_pass',all(result.cancellation.original_relative_gate_pass), ...
    'freeze_eligible',result.pass && pressure_delta<=1e-12);
benefit=true; preferred=true;
for f=frequencies
    t=result.timing; base=t(t.frequency_hz==f & strcmp(t.method,'full'),:);
    a=t(t.frequency_hz==f & strcmp(t.method,'fastin'),:); b=t(t.frequency_hz==f & strcmp(t.method,'composite'),:);
    benefit=benefit && median(a.per_realization_s)<median(base.per_realization_s) && ...
        median(b.per_realization_s)<median(base.per_realization_s) && ...
        median(a.cold_s)<median(base.cold_s) && median(b.cold_s)<median(base.cold_s);
    preferred=preferred && median(a.per_realization_s)<median(b.per_realization_s) && median(a.cold_s)<median(b.cold_s);
end
result.hardening.performance_benefit_confirmed=benefit;
result.hardening.fastin_preferred_by_medians=preferred;
result.hardening.freeze_eligible=result.hardening.freeze_eligible && benefit;
writetable(result.accuracy,fullfile(outdir,'accuracy.csv'));
writetable(result.guards,fullfile(outdir,'guards.csv'));
writetable(result.timing,fullfile(outdir,'timings.csv'));
writetable(result.memory,fullfile(outdir,'memory.csv'));
writetable(result.plan_controls,fullfile(outdir,'plan_controls.csv'));
writetable(result.cancellation,fullfile(outdir,'cancellation.csv'));
writetable(result.cancellation_scan,fullfile(outdir,'cancellation_scan.csv'));
artifact_names={'accuracy.csv','guards.csv','timings.csv','memory.csv','plan_controls.csv','cancellation.csv','cancellation_scan.csv'};
result.artifact_names=artifact_names;
result.artifact_hashes=cellfun(@(p)local_hash(fullfile(outdir,p)),artifact_names,'UniformOutput',false);
save(fullfile(outdir,'evaluation.mat'),'result','-v7');
manifest=struct('schema_version','1.1','scope',result.hardening.primary_scope, ...
    'relative_complex_gate',1e-3,'product_identity_gate',1e-12, ...
    'node_refinement_gate',1e-5,'rank_refinement_gate',1e-3,'quadrature_refinement_gate',1e-3, ...
    'reuse_identity_gate',1e-12,'freeze_eligible',result.hardening.freeze_eligible, ...
    'cancellation_relative_gate_pass',result.hardening.cancellation_relative_gate_pass, ...
    'evaluation_sha256',local_hash(fullfile(outdir,'evaluation.mat')), ...
    'previous_evaluation_sha256',result.hardening.previous_evaluation_sha256, ...
    'code_paths',{code_paths},'code_hashes',{code_hashes}, ...
    'artifact_names',{artifact_names},'artifact_hashes',{result.artifact_hashes});
fid=fopen(fullfile(outdir,'freeze_manifest.json'),'w','n','UTF-8'); assert(fid>=0);
fprintf(fid,'%s\n',jsonencode(manifest,PrettyPrint=true)); fclose(fid);
local_report(result,root);
local_hardening_report(result,previous,root);
fprintf('Algorithm gates pass=%d; max error fastin %.4g, composite %.4g\n',result.pass, ...
    max(result.accuracy.relative_error(strcmp(result.accuracy.method,'fastin'))), ...
    max(result.accuracy.relative_error(strcmp(result.accuracy.method,'composite'))));
end

function c=local_cfg(s,f)
source=s.source; source.kz=sqrt((2*pi*f/1500)^2-source.kx.^2);
c=struct('frequency_hz',f,'c0_mps',1500,'surface_x_m',s.x_m(:), ...
    'surface_z_m',s.eta_m(:),'surface_normal',[-s.slope(:) ones(size(s.slope(:)))]./sqrt(1+s.slope(:).^2), ...
    'spectral_source',source,'height_bound_m',.5,'height_nodes',33,'rank',3, ...
    'receiver_x_m',0,'receiver_z_m',3);
end
function g=local_guard(f,label,name,value,tolerance)
g=struct('frequency_hz',f,'case_name',label,'guard',name,'value',value,'tolerance',tolerance, ...
    'pass',isfinite(value) && value<=tolerance);
end
function a=local_append(a,b)
if isempty(a), a=b; else, a(end+1)=orderfields(b,a); end
end
function h=local_hash(path)
fid=fopen(path,'rb'); assert(fid>=0); closer=onCleanup(@()fclose(fid));
bytes=fread(fid,Inf,'*uint8'); md=java.security.MessageDigest.getInstance('SHA-256'); md.update(bytes);
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
function local_report(r,root)
fid=fopen(fullfile(root,'reports','fk_single_receiver_methods_report.md'),'w','n','UTF-8'); assert(fid>=0);
closer=onCleanup(@()fclose(fid));
fprintf(fid,'# 固定 Tx–单 Rx Kirchhoff 原型评价\n\n');
fprintf(fid,'算法门槛通过：%d。完整公式、护栏与预先规定的比较范围见 [方法与协议](fk_single_receiver_methods_protocol.md)。\n\n',r.pass);
fprintf(fid,'此PASS只对应原Rx(0,3)的51组输入。验证补强及相消点的明确失败记录见 [冻结审核](fk_single_receiver_baseline_freeze_report.md)；不以主域PASS覆盖相消诊断失败。\n\n');
fprintf(fid,'Tx100 m、Rx(0,3) m；均匀 c1500 m/s；4/6/8 kHz；既有984谱点 Gaussian 入射；海面[-60,60] m、4097点；高度界±.5 m、33 Chebyshev节点、统一入射rank3。所有方法使用相同真实高度、水侧法向、弧长梯形权重和完整两项积分，无幅相拟合。此为既有2-D研究基准，区别于Tx99 m安装图。\n\n');
fprintf(fid,'12个PM输入包含canonical260001及11个独立Gaussian Fourier系数样本（260002–260012，既有Sk与12模式带宽）；另有flat0/±.5和幅度.02/.45 m、K=.47正弦控制。冻结rank和节点后共51组输入，每种方案均独立对原始FullFK复声压验证。\n\n');
fprintf(fid,'|方案|最大相对复误差|最大相位误差/rad|\n|---|---:|---:|\n');
for m={'fastin','composite'}
    a=r.accuracy(strcmp(r.accuracy.method,m{1}),:);
    fprintf(fid,'|%s|%.6g|%.6g|\n',m{1},max(a.relative_error),max(abs(a.phase_error_rad)));
end
fprintf(fid,'\n数值护栏全部通过：%d；product identity最大%.4g；33→65节点最大%.4g；rank3→5最大%.4g；4097→8193积分最大%.4g。域越界、Rx变动、频率变动均拒绝旧plan。积分离散差单独报告，不归入压缩/插值误差。\n\n',all(r.guards.pass), ...
    max(r.guards.value(strcmp(r.guards.guard,'product_identity'))), ...
    max(r.guards.value(contains(r.guards.guard,'nodes33_65'))), ...
    max(r.guards.value(strcmp(r.guards.guard,'incoming_rank3_5'))), ...
    max(r.guards.value(strcmp(r.guards.guard,'quadrature4097_8193'))));
fprintf(fid,'MATLAB %s，同一共享会话；5配对计时块轮换顺序，每块6个不同海面。冷启动是热JIT下重新建plan加首次apply，非重启MATLAB；预计算未跨频率复用。重复成本包含共同海面求值、法向更新及各方法实际apply。文件I/O、报告、精度护栏不计时。下表为块中位数，原始五轮、顺序和范围见timings.csv。\n\n',r.version);
fprintf(fid,'|Hz|方案|预计算/ms|冷启动/ms|每realization/ms|重复加速Full/方案|冷启动加速|摊销交叉点/海面数|\n|---:|---|---:|---:|---:|---:|---:|---:|\n');
for f=r.frequencies
    base=r.timing(r.timing.frequency_hz==f & strcmp(r.timing.method,'full'),:);
    b=median(base.per_realization_s); bc=median(base.cold_s);
    for m={'full','fastin','composite'}
        t=r.timing(r.timing.frequency_hz==f & strcmp(r.timing.method,m{1}),:);
        prep=median(t.setup_s); a=median(t.per_realization_s); cold=median(t.cold_s);
        cross=Inf; if a<b, cross=max(1,ceil(prep/(b-a))); elseif strcmp(m{1},'full'), cross=0; end
        fprintf(fid,'|%g|%s|%.5g|%.5g|%.5g|%.4g|%.4g|%g|\n',f,m{1},1e3*prep,1e3*cold,1e3*a,b/a,bc/cold,cross);
    end
end
fprintf(fid,'\nplan常驻内存（whos，非峰值进程内存）：\n\n|Hz|FastIn/MiB|复合核/MiB|\n|---:|---:|---:|\n');
for j=1:height(r.memory), a=r.memory(j,:); fprintf(fid,'|%g|%.4g|%.4g|\n',a.frequency_hz,a.fastin_plan_bytes/2^20,a.composite_plan_bytes/2^20); end
fprintf(fid,'\n重复成本的五块min–max范围（不是置信区间）：\n\n|Hz|方案|每realization/ms 范围|\n|---:|---|---:|\n');
for f=r.frequencies
    for m={'full','fastin','composite'}
        t=r.timing(r.timing.frequency_hz==f & strcmp(r.timing.method,m{1}),:);
        fprintf(fid,'|%g|%s|%.5g–%.5g|\n',f,m{1},1e3*min(t.per_realization_s),1e3*max(t.per_realization_s));
    end
end
fprintf(fid,'\n## 结论与限制\n\n');
for f=r.frequencies
    t=r.timing(r.timing.frequency_hz==f,:);
    full=t(strcmp(t.method,'full'),:); ti=t(strcmp(t.method,'fastin'),:); tc=t(strcmp(t.method,'composite'),:);
    fprintf(fid,'%g Hz：FastIn重复加速%.4g倍，复合核%.4g倍；仅Green积分（精确入射已免费提供）%.5g ms。\n\n', ...
        f,median(full.per_realization_s)/median(ti.per_realization_s), ...
        median(full.per_realization_s)/median(tc.per_realization_s),1e3*median(full.full_integral_only_s));
end
fprintf(fid,'**推荐默认保留FastIn + FullOut独立原型**：三频率都通过1e-3复声压门槛，在本机冷启动与重复成本中位数均最低，plan约.63 MiB。其最大误差接近门槛，不能推断到新几何/海况仍通过。**同时保留复合核作为高精度备选**：本轮复误差接近双精度舍入水平、重复计算也明显优于重新构建入射的FullFK，但比FastIn需要更多预计算与约4.19 MiB内存，没有获得额外速度优势；不为追求排名继续强行压缩。两者速度差是当前会话中位数排序，未建立统计显著性。\n\n');
fprintf(fid,'**继续使用原始FullFK的条件**：当前海面精确入射已免费提供时，其标量积分成本见上列独立测量，明显低于两原型包含入射/插值的完整apply成本；越过域界或新场景未通过护栏时，也应回到FullFK验证。单次任务本轮FastIn冷启动仍实测最快，不能把新建plan成本遗漏。FastIn的Rx变化可复用入射plan，但不保证新Rx精度；复合核Rx变化须重建。Tx/频率/声速/源/网格/高度界/节点变化均须重建，rank变化只使FastIn失效。两方案仅作验证原型，不接入生产。\n\n');
fprintf(fid,'FullFK物理模型误差未重新评价；本轮只是对同一离散Kirchhoff的算法保真。既有canonical BIE物理锚点不能当作11个新海面或强正弦的物理验证。未建立PM人口可靠性、3-D、深海传播、统计信道生成、CIR或BER结论。生产PE/公共接口和FullFK/BIE基准均未修改，因此未重复无关生产通信回归。\n\n');
fprintf(fid,'复现：`setup_vertical_project; result=validate_fk_single_receiver_methods();`；只刷新报告：`validate_fk_single_receiver_methods(''render'')`；审计：`audit_fk_single_receiver_methods_artifacts()`。结果：`results/validation/fk_single_receiver_methods/evaluation.mat`、accuracy/guards/timings/memory.csv及12份系数CSV。原始reference SHA256运行前后相同，执行代码SHA256保存在MAT。\n');
end

function local_hardening_report(r,previous,root)
path=fullfile(root,'reports','fk_single_receiver_baseline_freeze_report.md');
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); closer=onCleanup(@()fclose(fid));
fprintf(fid,'# 单接收点 Kirchhoff 研究基线冻结审核（2026-10-12）\n\n');
fprintf(fid,'**冻结结论：原Tx100–Rx(0,3)研究基线满足条件=%d；相消点扩展的1e-3门槛通过=%d。** 原固定点保留FastIn + FullOut为默认候选、复合核为高精度备选。新诊断不支持对任意Rx沿用FastIn精度声明；失败点保留，不放宽门槛、不提高生产rank、不引入新算法。\n\n', ...
    r.hardening.freeze_eligible,r.hardening.cancellation_relative_gate_pass);
fprintf(fid,'## 修改与逐项证据\n\n');
fprintf(fid,'1. 拒绝验证改为精确匹配expected_id/actual_id；无关异常与无异常注入均证明验证器拒绝误通过。%d条plan/验证器控制全部通过。\n',height(r.plan_controls));
fprintf(fid,'2. 原型只修复契约缺陷：高度域、无效谱、旧plan、错误kind、缺失plan有独立标识符；FastIn去掉无依赖Rx绑定，复合核去掉无依赖rank绑定，忽略未参与预计算的源metadata。两种方法的SVD、Green、插值、载波与积分公式均未改变。高度/法向更新、metadata及上述无关参数更新的复用值与重建值按1e-12验证。源coeff/kx/kz、周期、归一化、原点、Tx、频率、声速、网格、高度界、节点及FastIn rank的失效逐项测试。\n');
fprintf(fid,'3. 主域%d条误差记录的复声压、绝对/相对复误差、绝对/相对幅度误差、定义有效的相位和raw相对误差均显式检查有限。NaN/Inf复声压注入被拒绝；exact zero标记相位和raw相对误差undefined（NaN），不以0相位宣称准确。1e-16控制同时保留分母保护值及raw比例，揭示1e-14 floor的局限。\n',height(r.accuracy));
fprintf(fid,'4. 同原51组Gaussian/PM/flat/sine输入重跑4/6/8 kHz；%d条原公式/节点/rank/积分护栏全部通过，原1e-3/1e-12/1e-5门槛不变。与pre_hardening保存证据比较，原近似声压最大绝对变化%.4g（门槛1e-12）。FullFK/BIE/原FastFK SHA256前后一致。\n',height(r.guards),r.hardening.max_pressure_change);
fprintf(fid,'5. timing仍是5配对块、6个不同既有海面/块，同热JIT冷启动定义；未新增海面样本、MC、SSP、3-D或通信计算。PROJECT_CONTEXT、技术指南、运行索引和协议同步本结论。\n\n');
fprintf(fid,'两方案冷/重复成本均优于FullFK=%d；FastIn的冷/重复中位数均最低=%d。不是统计显著性判断。\n\n',r.hardening.performance_benefit_confirmed,r.hardening.fastin_preferred_by_medians);
fprintf(fid,'|方案|主域最大相对复误差|最大幅度相对误差|最大相位误差/rad|\n|---|---:|---:|---:|\n');
for m={'fastin','composite'}
    a=r.accuracy(strcmp(r.accuracy.method,m{1}),:);
    fprintf(fid,'|%s|%.6g|%.6g|%.6g|\n',m{1},max(a.relative_error),max(a.magnitude_relative_error),max(abs(a.phase_error_rad)));
end
fprintf(fid,'\n|Hz|FastIn重复加速：旧→新|复合核重复加速：旧→新|新冷启动加速FastIn/复合核|\n|---:|---:|---:|---:|\n');
for f=r.frequencies
    speeds=zeros(2,3);
    records={previous,r};
    for j=1:2
        t=records{j}.timing; base=t(t.frequency_hz==f & strcmp(t.method,'full'),:);
        for m=1:2
            names={'fastin','composite'}; v=t(t.frequency_hz==f & strcmp(t.method,names{m}),:);
            speeds(j,m)=median(base.per_realization_s)/median(v.per_realization_s);
        end
    end
    t=r.timing; base=t(t.frequency_hz==f & strcmp(t.method,'full'),:);
    a=t(t.frequency_hz==f & strcmp(t.method,'fastin'),:); b=t(t.frequency_hz==f & strcmp(t.method,'composite'),:);
    fprintf(fid,'|%g|%.4g→%.4g|%.4g→%.4g|%.4g / %.4g|\n',f,speeds(1,1),speeds(2,1),speeds(1,2),speeds(2,2), ...
        median(base.cold_s)/median(a.cold_s),median(base.cold_s)/median(b.cold_s));
end
fprintf(fid,'\n## 相消诊断：失败不隐藏\n\n');
fprintf(fid,'使用已有A=.45 m/K=.47 rad/m强正弦面；只以原FullFK在Rx x∈[-60,60]、z3扫描并在|x|≤50的局部低压点细化位置，不拟合或改动源/幅相。另记录`|sum(contribution)|/sum(|contribution|)`量化积分相消。Rx扫描只用于诊断，不计入单点加速计时。near_zero定义为低于扫描峰值1e-3或绝对1e-14，属于解释标签，不替代原误差门槛。\n\n');
fprintf(fid,'|Hz|方案|Rx x/m|峰值归一化压力|相消比|绝对复误差|相对复误差|相位误差/rad|原1e-3通过|\n|---:|---|---:|---:|---:|---:|---:|---:|---|\n');
for j=1:height(r.cancellation)
    d=r.cancellation(j,:);
    fprintf(fid,'|%g|%s|%.7g|%.4g|%.4g|%.4g|%.4g|%.4g|%d|\n',d.frequency_hz,d.method{1},d.receiver_x_m, ...
        d.peak_relative_magnitude,d.cancellation_ratio,d.absolute_error,d.relative_error,d.phase_error_rad,d.original_relative_gate_pass);
end
fprintf(fid,'\n同位置的误差分解控制（仅诊断，不改变基线rank/节点）：\n\n|Hz|FastIn rank5复误差|复合核65节点复误差|FullFK4097→8193差|\n|---:|---:|---:|---:|\n');
for f=r.frequencies
    d=r.cancellation(r.cancellation.frequency_hz==f,:); d=d(1,:);
    fprintf(fid,'|%g|%.4g|%.4g|%.4g|\n',f,d.rank5_relative_error,d.nodes65_relative_error,d.quadrature4097_8193_relative);
end
fprintf(fid,'\n## 所有未关闭问题与冻结边界\n\n');
fprintf(fid,'- FastIn rank3在三个新增低压位置均未满足1e-3；这些Rx不在原x0冻结域内。plan可复用是依赖正确性，不是任意Rx精度保证；不能把主域PASS解释为这些失败点被接受。\n');
fprintf(fid,'- 相消位置的积分离散敏感性按上表单独保留；若超过原1e-3，连FullFK参考的相对收敛也未建立。不能把全部残差归为压缩或插值，更不能宣称物理准确。\n');
fprintf(fid,'- 接近零点时相对误差和相位病态；绝对误差、峰值归一化压力及相消比须同时阅读；exact zero相位无定义。复合核有限节点仍是近似，新Rx不再保证原x0的舍入级误差。\n');
fprintf(fid,'- 无新增BIE物理验证、无限孔径收敛、PM人口可靠性、3-D/深海/SSP/CIR/BER证据；本次冻结是既有2-D离散研究fixture的算法和验收基线。\n');
fprintf(fid,'- 共享会话计时随负载变化，排序按本轮中位数，未证明统计显著性或独立进程冷启动。精确入射已免费提供时继续原标量FullFK。\n\n');
fprintf(fid,'复现：`setup_vertical_project; validate_fk_single_receiver_methods(); audit_fk_single_receiver_methods_artifacts();`。MAT/CSV在`results/validation/fk_single_receiver_methods/`；补强前原始证据保存在该目录`pre_hardening/`，不覆盖删除。新plan_controls/cancellation/cancellation_scan.csv保留原始状态与失败。冻结资格=%d不包含相消点通过声明。\n',r.hardening.freeze_eligible);
end

