function [ output ] = CARPE3D( params )
% 主调度函数：
% 根据输入参数 params 完成 3D PE 仿真的数值设置，
% 调用底层宽角传播器 propWAPErev11，
% 画出若干结果图，并将主要结果打包输出。
%
% 输入:
%   params   一个结构体，里面包含频率、声速、网格、环境文件名、
%            声源位置、输出次数等仿真参数
%
% 输出:
%   output   结构体，包含最终场、坐标网格、以及若干后处理量


% ===== 1) 基本物理量与步进参数 =====

lambda0 = params.c0 / params.f0;
% 参考波长 lambda0 = c0 / f0
% 这里用基准声速 c0 和工作频率 f0 计算波长

range_update_interval = 2 * round(.5 * params.range_update_in_meters / lambda0);
% 环境更新间隔（以“步数”或“波长步”为单位）
%
% 解释：
% params.range_update_in_meters 是“每隔多少米更新一次环境”
% 这里先除以 lambda0，把“米”换成“多少个波长”
% 然后通过 2*round(.5*...) 把它取整成一个更规整的整数
%
% 直观理解：
% 不是每一步 x 推进都更新环境，而是隔若干步更新一次


wid = lambda0 * (600/2^11) * params.ny;
% 横向半宽度（half-width）相关参数？
%
% 注释里写的是：horizontal width, domain is TWICE this across
% 即：
% wid 是“半边宽度”的量，整个横向计算域宽度大约是 2*wid
%
% 这里写法说明：
% 原始参考设置可能是 ny=2^11 时半宽约 600 m
% 现在按 ny 做比例缩放

steplength = (params.steplamb) * lambda0;
% x 方向每一步前进的物理距离（米）
% params.steplamb 是“多少个波长为一步”
% 所以实际步长 = 波长数 × 波长

numstep = round(params.rng / steplength)
% 总步数 = 最大传播距离 / 每步长度
% 没有分号，所以会在命令行显示 numstep

range = steplength * numstep;
% 实际总传播距离
% 因为 numstep 是 round 得来的，实际总距离可能和 params.rng 略有差异

disp(['x final range (km) :  ' num2str(range/1000)])
% 在命令行显示最终传播距离（单位 km）


% ===== 2) 频率设置 =====

freqs = [ params.f0 ];
% 当前只算单频
% 注释说明以后可扩展成多频循环（宽带）
% 但目前输出组织和文件保存还没有完全为宽带重构


for freq = freqs
% 频率循环
% 虽然现在只有一个频率，但代码结构已经预留了宽带/多频的外壳

  tic
  % 开始计时

  freq
  % 显示当前频率（因为没有分号）


  % ===== 3) 调用真正的 3D 宽角传播器 =====
  % propagate x direction

  x = 0;
  % 初始化主传播方向坐标起点
  % 注意：真正的 marching 是在 propWAPErev11 里做的
  % 这里的 x=0 更像是逻辑上的起始位置

  [psifinal, x, y, z, M, psi1, Ez, Af] = propWAPErev11( ...
      params.ny, params.nz, ...
      numstep, freq, params.c0, params.Envfile, wid, steplength, ...
      params.aspect, params.depP, params.zs, params.ys, params.filen, ...
      params.nout, params.rdcw_flag, range_update_interval);
  %
  % 这是整个函数最核心的一句：
  % 调用底层宽角 3D PE 传播器
  %
  % 输入的大意：
  %   ny, nz                  y/z 网格点数
  %   numstep                 x 方向总步数
  %   freq                    当前频率
  %   c0                      基准声速
  %   Envfile                 环境脚本名
  %   wid                     横向半宽
  %   steplength              x 步长
  %   aspect                  网格纵横比
  %   depP                    参考水深/绘图用深度参数
  %   zs, ys                  声源深度、横向位置
  %   filen                   环境数据文件名
  %   nout                    输出次数
  %   rdcw_flag               距离相关/无关环境标志
  %   range_update_interval   环境更新间隔
  %
  % 输出的大意（结合变量名推测）：
  %   psifinal   最终/输出的声场
  %   x,y,z      对应坐标
  %   M          某个中间量或传播相关矩阵/变量（此函数后面没用到）
  %   psi1       初始场或起始参考场
  %   Ez         深度平均或某种 y-x 平面量
  %   Af         某个沿 x-z 的截面场（通常在 source y 位置附近）


  % ===== 4) 做归一化与结果可视化 =====

  anorm1 = max(max(abs(psi1)));
  % 初始场 psi1 的最大幅值，用作一个参考归一化量

  anormfin = max(max(abs(psifinal)));
  % 最终场幅值最大值
  % 注意：这个变量后面实际上没有用到

  figure(2); clf
  % 使用图窗 2，并清空旧图

  subplot(211)
  % 图 2 的上半部分：画最终场幅值（dB）

  finaldb = 20 * log10(abs(squeeze(psifinal(end,:,:)) / anorm1));
  % 取 psifinal 的最后一个 x 输出截面（或最后一个保存截面）
  % squeeze 后得到一个 y-z 平面
  % 再相对初始场最大值 anorm1 做 dB 归一化
  %
  % 直观理解：
  % “最终平面上的声压场强弱图（dB）”

  finlin = squeeze(abs(psifinal(end,:,:)));
  % 同一个最终平面的线性幅值
  % 后面虽然算了，但实际绘图用的是 finaldb，不是 finlin

  wat = find(z > -(params.depP+5) & z < .1);
  % 选出要显示的 z 范围索引
  % 大意是选水体附近的深度范围
  %
  % 注意：
  % 这里的 z 坐标正方向带有历史代码约定，
  % 你可以把它粗略理解为“截出接近海面到海底附近的一段来画图”

  % pcolor(y(1:2:end), z(wat), finlin(wat,1:2:end)); shading flat
  % 这是一种备选画法：画线性幅值图
  % 目前被注释掉了

  pcolor(y(1:2:end), z(wat), finaldb(wat,1:2:end)); shading flat
  % 画最终平面上的幅值 dB 图
  % y 方向每隔一个点取样一次，减少绘图量
  % z 只取 wat 对应的范围

  caxis(max(finaldb(:)) + [-20 0])
  % 颜色范围设为“距最大值 20 dB 内”
  % 即只突出最强的那一部分结构

  title('|\psi_{final}| (dB)')
  % 图标题：最终场幅值（dB）

  ylabel('depth (m)')
  xlabel('distance (m)')
  % 轴标签
  % 注意：这里 xlabel 写的是 distance (m)，但横轴实际传入的是 y
  % 这一点从代码上看是有点混乱/沿用了旧图标签
  % 更准确地理解，这里横轴是 y 截面位置

  % caxis([ max(finaldb(:))-40 max(finaldb(:))+2])
  % 另一套颜色范围设置，已注释

  % axis equal
  % 若开启则坐标等比例，当前没用

  colorbar
  % 加颜色条


  subplot(212)
  % 图 2 的下半部分：画最终场相位

  % figure(2); clf
  % 老代码残留，当前不用

  pcolor(y(1:2:end), z(wat), angle(squeeze(psifinal(end,wat,1:2:end)))); shading flat
  % 画最终平面的相位角
  % angle(...) 给出复数场的相位

  title('phase angle \psi_{final}')
  ylabel('depth (m)')
  xlabel('distance (m)')
  % 同样，xlabel 的文字可能沿用了旧命名
  % 但这里实质是在画 y-z 平面上的相位图

  % axis equal
  colorbar

  pause(4)
  % 暂停 4 秒，方便图形刷新或人工查看


  % eval([ 'save ' outdir 'psi' num2str(freq) ' y z x psifinal psi1 freq c0 range steplength -V6' ])
  % 旧版保存方式，当前已注释

  toc
  % 输出该频率计算所花时间

end  % freq loop
% ===== 频率循环结束 =====


% ===== 5) 下列是旧代码残留/备选分析 =====

% anorm2=max(max(abs(psiout2)));
% figure(2)
% imagesc(x,z,20*log10(abs(psiout2)/anorm2),[-80 0]);

% save outputfile z x y psi1 psifinal freq dep zs ...
%    numstep steplength nz range f0 c0 aspect Af Ez


% ===== 6) 构造后处理所需坐标网格 =====

clear finaldb finlin
% 清理临时变量

% eval(['save ' outfile]);
% 旧版保存方式，已注释


xg = (1:numstep) * steplength;
% x 方向每一步对应的物理位置坐标（行向量）

xgg = xg' * ones(1, params.ny);
% 扩展成 numstep × ny 的二维网格
% 每一行是一个 x，每一列对应一个 y 位置
% 便于后面和 y 组成平面坐标

ygg = ones(numstep,1) * y;
% 扩展成 numstep × ny 的二维 y 网格
% 每一行都复制 y 向量

rgg = sqrt(xgg.^2 + ygg.^2);
% 对应 x-y 平面上每个点到原点的斜距（或水平范围半径）

zgg = ones(numstep,1) * (z(1:params.nz/2) - params.zs);
% 生成与 x 对应的 z 网格，并减去声源深度 zs
% 相当于把 z 改写成“相对声源深度”的坐标

rzgg = sqrt(xgg(:,1:params.nz/2).^2 + zgg.^2);
% 构造 x-z 平面上的几何距离量
% 用于后面 Af 图的参考量


% ===== 7) 若 Ez 存在，画 x-y 平面的后处理图 =====

if (exist('Ez'))
% 如果底层传播器输出了 Ez，就画图 3

    figure(3); clf

    refI = max(max(Ez' .* rgg'));
    % 构造一个参考强度量
    % 从写法看，是用 Ez 乘以距离修正因子 rgg 后取最大值作为参考

    pcolor(xg(1:3:end), y(1:3:end), ...
        10*log10((Ez(1:3:end,1:3:end)/refI)' .* rgg(1:3:end,1:3:end)')); shading flat
    % 画 Ez 的 x-y 分布图
    % 标题写的是 “depth average Ez”
    % 因此可以理解为某种深度平均量在水平平面上的分布
    %
    % 这里做了稀疏采样（每隔几个点取一个）来降低绘图量

    caxis([ -40 -10])
    colorbar
    xlabel('x (m)')
    ylabel('y (m)')
    title('depth average Ez','interpreter','none')
end


% ===== 8) 若 Af 存在，画 x-z 平面的后处理图 =====

if (exist('Af'))
% 如果底层传播器输出了 Af，就画图 4

    figure(4); clf

    % pcolor(xg',z(1:nz/2),20*log10(abs(Af).*sqrt(rzgg')));shading flat
    % 旧版本画法，已注释

    pcolor(xg(1:3:end)', z(1:2:params.nz/2), 20*log10(abs(Af(1:2:end,1:3:end)))); shading flat
    % 画 Af 的 x-z 分布图
    % 从标题看，它是在“ys，也就是声源横向位置处”的一个垂向切面

    caxis([ -20*log10(range) -30])
    % 设置 dB 颜色范围

    colorbar
    xlabel('x (m)')
    ylabel('z (m)')
    title('|Af| dB at ys, the source location in the transverse dim.', ...
        'interpreter','none')
    % 标题意思：
    % “在 ys 这个横向位置（也就是声源所在横向位置）的 x-z 切面上，Af 的 dB 图”
end


% ===== 9) 打包输出 =====

output.psifinal = psifinal;
% 最终/输出声场

output.xgg = xgg;
% x 网格

output.ygg = ygg;
% y 网格

output.rgg = rgg;
% x-y 平面上的距离网格

output.zgg = zgg;
% z 网格（相对声源深度）

output.Ez = Ez;
% 深度平均量或水平图量（若有）

output.Af = Af;
% 在 source y 位置附近的 x-z 截面量（若有）

end