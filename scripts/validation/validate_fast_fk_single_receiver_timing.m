function result = validate_fast_fk_single_receiver_timing()
%VALIDATE_FAST_FK_SINGLE_RECEIVER_TIMING Frozen 4 kHz PM single-point audit.
% FullFK receives one scalar coordinate; FastFK keeps its full-grid solver.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
outdir=fullfile(root,'results','validation','fast_fk_single_receiver_timing');
if ~isfolder(outdir), mkdir(outdir); end
input=fullfile(root,'results','validation', ...
    'spectral_fast_kirchhoff_multifrequency_consistency','consistency.mat');
clock=tic; loaded=load(input,'s'); io_seconds=toc(clock); s=loaded.s;
assert(s.pass && isequal(s.ranks,[3 12]) && s.padding==2);
hashes=s.reference_hashes;
for j=1:numel(hashes)
    assert(strcmp(local_hash(fullfile(root,hashes(j).path)),hashes(j).sha256), ...
        'Frozen reference changed: %s',hashes(j).path);
end
frequency=4000; n=numel(s.x_m); assert(n==4097);
idx=find(s.x_m==0); assert(isscalar(idx) && s.mask(idx));
saved=s.sparse([s.sparse.frequency_hz]==frequency);
assert(isscalar(saved) && saved.rank_in==3 && saved.rank_out==12);
clock=tic;
src=s.source; k=2*pi*frequency/1500; src.kz=sqrt(k^2-src.kx.^2);
cfg=struct('frequency_hz',frequency,'c0_mps',1500,'surface_x_m',s.x_m, ...
    'surface_z_m',s.eta_m,'surface_normal',[-s.slope,ones(n,1)]./sqrt(1+s.slope.^2), ...
    'height_bound_m',s.height_bound_m,'receiver_depth_m',3,'spectral_source',src, ...
    'max_rank',12,'height_nodes',33,'geometry','B','rank',3,'outgoing_rank',12,'padding_factor',2);
config_seconds=toc(clock);
% Warmup of both preparation and evaluation paths; no full-grid FullFK call.
fcfg=local_full_input(cfg,s.x_m(idx));
plan=prepare_spectral_fast_kirchhoff_2d(cfg);
full=solve_full_kirchhoff_2d(fcfg);
fast=solve_spectral_fast_kirchhoff_2d(cfg,plan);
assert(isscalar(full.receiver_field) && numel(fast.receiver_field)==n);
assert(abs(full.receiver_field-saved.full(idx))/abs(saved.full(idx))<1e-12);
assert(norm(fast.receiver_field-saved.fast)/norm(saved.fast)<1e-12);
% Nine independent rounds, ten calls per apply batch; alternate method order.
rounds=9; batch=10; prep=zeros(rounds,2); apply=prep;
for j=1:rounds
    order=[1 2]; if mod(j,2)==0, order=[2 1]; end
    for method=order
        clock=tic;
        if method==1
            fcfg=local_full_input(cfg,s.x_m(idx));
        else
            plan=prepare_spectral_fast_kirchhoff_2d(cfg);
        end
        prep(j,method)=toc(clock);
        clock=tic;
        for b=1:batch
            if method==1
                full=solve_full_kirchhoff_2d(fcfg);
            else
                fast=solve_spectral_fast_kirchhoff_2d(cfg,plan);
                fast_point=fast.receiver_field(idx);
            end
        end
        apply(j,method)=toc(clock)/batch;
    end
end
result=struct('schema_version','1.0','matlab_version',version, ...
    'input_path',input,'input_sha256',local_hash(input),'reference_hashes',hashes, ...
    'frequency_hz',frequency,'receiver_index',idx,'receiver_x_m',s.x_m(idx), ...
    'receiver_z_m',3,'cfg',cfg,'rounds',rounds,'calls_per_batch',batch, ...
    'io_wall_seconds',io_seconds,'config_wall_seconds',config_seconds, ...
    'prep_wall_seconds',prep,'apply_wall_seconds',apply, ...
    'full_point',full.receiver_field,'fast_point',fast_point, ...
    'saved_full_point',saved.full(idx),'saved_fast_point',saved.fast(idx), ...
    'historical_full_cost',saved.cost,'historical_fast_setup',saved.fast_setup_wall, ...
    'historical_fast_apply',saved.fast_apply_wall, ...
    'historical_warm_speedup',saved.warm_speedup,'historical_cold_speedup',saved.cold_speedup);
result.absolute_complex_error=abs(result.fast_point-result.full_point);
result.relative_complex_error=result.absolute_complex_error/abs(result.full_point);
result.phase_error_rad=angle(result.fast_point/result.full_point);
result.relative_magnitude_error=abs(abs(result.fast_point)-abs(result.full_point))/abs(result.full_point);
result.warm_speedup=median(apply(:,1))/median(apply(:,2));
result.setup_included_speedup=median(prep(:,1)+apply(:,1))/median(prep(:,2)+apply(:,2));
assert(result.relative_complex_error<1e-3 && isfinite(result.warm_speedup));
save(fullfile(outdir,'single_receiver_timing.mat'),'result','-v7');
t=table((1:rounds)',prep(:,1),prep(:,2),apply(:,1),apply(:,2), ...
    prep(:,1)+apply(:,1),prep(:,2)+apply(:,2),'VariableNames', ...
    {'round','full_prepare_s','fast_prepare_s','full_apply_s','fast_apply_s','full_setup_apply_s','fast_setup_apply_s'});
writetable(t,fullfile(outdir,'timings.csv'));
local_report(result,root);
fprintf('Single-point Full %.7g s vs whole-grid Fast %.7g s; Full/Fast %.7g\n', ...
    median(apply(:,1)),median(apply(:,2)),result.warm_speedup);
fprintf('Setup included Full/Fast %.7g; complex error %.7g\n', ...
    result.setup_included_speedup,result.relative_complex_error);
end

function r=local_full_input(cfg,rx)
% Exact unchanged angular-spectrum incident evaluation from the existing audit.
src=cfg.spectral_source; x=cfg.surface_x_m(:); z=cfg.surface_z_m(:);
u=complex(zeros(size(x))); ux=u; uz=u;
for first=1:128:numel(x)
    jj=first:min(first+127,numel(x));
    phase=exp(1i*(x(jj)-src.x0_m)*src.kx.').*exp(-1i*(z(jj)-src.z_tx_m).*src.kz.');
    u(jj)=phase*src.coeff/src.n_original;
    ux(jj)=phase*(1i*src.kx.*src.coeff)/src.n_original;
    uz(jj)=phase*(-1i*src.kz.*src.coeff)/src.n_original;
end
r=struct('frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps, ...
    'surface_x_m',x,'surface_z_m',z,'surface_normal',cfg.surface_normal, ...
    'surface_incident',u,'surface_incident_dn',cfg.surface_normal(:,1).*ux+cfg.surface_normal(:,2).*uz, ...
    'receiver_x_m',rx,'receiver_z_m',cfg.receiver_depth_m);
assert(isscalar(r.receiver_x_m) && isscalar(r.receiver_z_m));
end

function h=local_hash(path)
fid=fopen(path,'rb'); assert(fid>=0); closer=onCleanup(@()fclose(fid));
bytes=fread(fid,Inf,'*uint8');
md=java.security.MessageDigest.getInstance('SHA-256'); md.update(bytes);
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end

function local_report(r,root)
path=fullfile(root,'reports','fast_fk_single_receiver_timing_report.md');
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); closer=onCleanup(@()fclose(fid));
fprintf(fid,'# 单接收点 FullFK 与当前 FastFK 效率对比\n\n');
fprintf(fid,'日期：2026-10-08。MATLAB %s；共享 MATLAB 会话，CPU double。\n\n',r.matlab_version);
fprintf(fid,'仅新增一个单频、单接收点测试，未修改 FullFK/FastFK 算法，未运行 PE、BIE、新宽带或接收点数量扫描。\n\n');
fprintf(fid,'## 冻结配置\n\n');
fprintf(fid,'复用已通过的 `consistency.mat`：seed260001、12 个固定 PM Fourier 模式、Kmax=0.471238898 rad/m、C2 plateau42/support50 m。频率4000 Hz，c=1500 m/s，Tx=100 m，Rx=3 m；这是原研究夹具。\n\n');
fprintf(fid,'海面积分及 FastFK 接收网格均为 [-60,60] m、4097点，dx=120/4096 m；源为既有984点、sigma=0.3 m Gaussian，保留全部既有传播谱，精确 kz。高度界0.5 m，height nodes33，geometry B，rank_in/out=3/12，max_rank12，padding2（FFT8194）。直接复用保存 eta/slopes/source，所有 reference SHA256 均校验通过。\n\n');
fprintf(fid,'唯一接收点：MATLAB index%d，x=%.9g m，z=3 m，位于既有 M99 核心。FullFK 输入严格为标量 receiver_x/z，仅形成1×4097积分核；FastFK仍求4097点整场后提取此点。\n\n',r.receiver_index,r.receiver_x_m);
fprintf(fid,'## 实测时间与范围\n\n');
fprintf(fid,'两种方法的准备和求值各预热一次；9轮交替执行方法顺序，每轮准备一次、连续求值10次，批次总wall time除以10。下表为9轮中位数及[min,max]。外层tic/toc包含调用、输入检查和输出构造，Fast包含单点提取；不包含磁盘IO、hash检查、报告保存。\n\n');
fprintf(fid,'|阶段|FullFK严格单点 / s|FastFK整网格 / s|\n|---|---:|---:|\n');
arrays={r.prep_wall_seconds,r.apply_wall_seconds,r.prep_wall_seconds+r.apply_wall_seconds};
labels={'准备','实际计算（每次）','准备＋一次计算'};
for j=1:3
    a=arrays{j};
    fprintf(fid,'|%s|%.8g [%.8g, %.8g]|%.8g [%.8g, %.8g]|\n', ...
        labels{j},median(a(:,1)),min(a(:,1)),max(a(:,1)),median(a(:,2)),min(a(:,2)),max(a(:,2)));
end
fprintf(fid,'\nFull准备：在全部海面样点上精确计算入射u/ux/uz并组装dn与标量receiver配置；Fast准备：既有incoming/outgoing SVD、Hankel lag核、插值与plan。Full求值中的弧长求积权重构造仍计入计算。输入准备不是单点Hankel积分本身，两者的准备内容不同。\n\n');
fprintf(fid,'公共输入文件读取 %.8g s，固定配置组装（含kz、normal）%.8g s，单独记录、不计入上述速度比。海面直接复用保存数组；其历史生成时间不在本轮计时内。准备＋一次计算是已预热实现下的plan重建成本，不代表首次启动MATLAB/JIT的耗时。\n\n',r.io_wall_seconds,r.config_wall_seconds);
fprintf(fid,'定义 speedup=T_Full/T_Fast：计算阶段 **%.8g**；含准备 **%.8g**。即Full单点计算比当前Fast整场快 **%.4g倍**，含准备快 **%.4g倍**。计算阶段绝对差 %.8g s；含准备差 %.8g s（各方法中位数之差）。\n\n', ...
    r.warm_speedup,r.setup_included_speedup,1/r.warm_speedup,1/r.setup_included_speedup, ...
    median(r.apply_wall_seconds(:,2))-median(r.apply_wall_seconds(:,1)), ...
    median(sum([r.prep_wall_seconds(:,2),r.apply_wall_seconds(:,2)],2))-median(sum([r.prep_wall_seconds(:,1),r.apply_wall_seconds(:,1)],2)));
fprintf(fid,'## 复数结果与精度\n\n');
fprintf(fid,'原生物理反射相量，exp(-i omega t)，不作相位对齐或复数拟合。\n\n');
fprintf(fid,'- FullFK：`%.16g %+.16gi`\n- FastFK：`%.16g %+.16gi`\n',real(r.full_point),imag(r.full_point),real(r.fast_point),imag(r.fast_point));
fprintf(fid,'- 绝对复误差 |Fast-Full|：%.9g\n- 相对复误差 |Fast-Full|/|Full|：%.9g（%.6g%%）\n- 相对幅度误差：%.9g\n- 相位差 angle(Fast/Full)：%.9g rad\n\n', ...
    r.absolute_complex_error,r.relative_complex_error,100*r.relative_complex_error,r.relative_magnitude_error,r.phase_error_rad);
fprintf(fid,'Full单点与既有保存整场的该点复数值在1e-12相对容差内一致；Fast本轮整场与保存rank3/12整场在1e-12相对L2容差内一致。单点复误差低于1e-3。这里保存整场只用于一致性断言，未用于新Full单点计算或计时。\n\n');
fprintf(fid,'## 与此前整场测试对照\n\n');
fprintf(fid,'同一4 kHz保存记录：Full4097点准备/计算 %.8g / %.8g s；Fast4097点准备/计算 %.8g / %.8g s；warm/cold加速比 %.6g / %.6g。\n\n', ...
    r.historical_full_cost.setup_wall,r.historical_full_cost.apply_wall,r.historical_fast_setup,r.historical_fast_apply,r.historical_warm_speedup,r.historical_cold_speedup);
fprintf(fid,'此前Full整场计时为单次外层调用，Fast计算为预热后3次内部wall_seconds中位数，cold含各自准备、排除共同海面生成和IO。本轮Full只有1点，Fast仍4097点，双方使用相同外层计时且9批×10次。因此旧整场加速比只说明整场任务，不能拿来代表本轮单点任务；本轮并未重新测整场Full，也未混用旧耗时计算新速度比。\n\n');
fprintf(fid,'## 判断\n\n');
fprintf(fid,'当前实现用于此单浮标任务时，FullFK直接单点更快；FastFK整场加速优势不延续到单点。入射场已经可复用时，应看计算阶段；每次重建输入/plan时，应看含准备结果。两种范围均给出同一速度排序。\n\n');
fprintf(fid,'不建议仅为该固定单频单浮标场景继续优化当前整场FastFK，优先使用FullFK单点并复用固定入射数据。若以后另有大量重复查询且可长期摊销plan的需求，才值得另立单点算法研究；本轮未测该需求，不能推算其收益或交叉点。结论仅限本配置与本机共享会话，不推广到其他点数、频率、宽带或新算法。\n\n');
fprintf(fid,'复现实验：`setup_vertical_project; validate_fast_fk_single_receiver_timing();`。原始9轮计时与复数结果位于 `results/validation/fast_fk_single_receiver_timing/timings.csv` 和 `single_receiver_timing.mat`。输入SHA256：`%s`。\n',r.input_sha256);
end
