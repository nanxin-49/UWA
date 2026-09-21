% Execute the SimCARPE3D code.
% 运行 SimCARPE3D / CARPE3D 的测试入口脚本

% td 28 May 2009
% 原始版本日期

% td 20 July 2023
% 2023 年修订日期

clear
% 清空工作区变量，避免旧变量干扰本次运行

format compact
% 让 MATLAB 命令行输出更紧凑，减少空行


%! rm diary.txt
% 调用操作系统命令删除旧的日志文件 diary.txt
% 注意：这是 Linux / Unix 风格命令；如果你在 Windows 上运行，通常需要改成别的方式

diary diary.txt
% 打开 MATLAB 日志功能，把后续命令行输出保存到 diary.txt


disp(' host:')
% 在命令窗口打印字符串 " host:"

eval('! hostname');
% 调用系统命令 hostname，显示当前计算机主机名
% 用 eval 包一层，本质上还是执行 shell 命令

disp(' ')
% 打印一个空行

%! cat /etc/os-release
% 显示当前操作系统版本信息（Linux 文件）
% 例如 Ubuntu / Debian / CentOS 版本等

disp(' ')
% 打印一个空行

%! lscpu
% 显示 CPU 硬件信息：核心数、架构、线程等

disp(' ')
% 打印一个空行

[ 'script executed : ' date]
% 生成一个字符串，形如 "script executed : 24-Jul-2025"
% 注意：这里没有分号，所以 MATLAB 会把这个字符串显示出来

disp (' ')
% 打印一个空行


%%    specify parameters
% ===== 参数设置区 =====

rng=   20000;  % max x propagation distance
% 最大传播距离，单位米
% 这里的 rng 不是随机数种子，而是 range（传播距离）
% 程序会沿 x 方向把声场推进到 20000 m

f0=200; % Hz   %  single frequency to compute, or initial frequuency for broadband
           % code rev11 needs modification to do broadband
% 计算频率，200 Hz
% 当前 rev11 版本主要按单频运行
% 注释说明：若要真正做宽带，需要进一步修改代码？？

zs= 20;  % source depth
% 声源深度，20 m

depP= 95 ;  % water depth  % used for plotting only, when Enfile='Wedge002' or Wedge003
% 水深 95 m
% 注意：注释说明这个变量在当前楔形环境下主要用于绘图，不一定直接控制实际传播环境

steplamb=1;  % x step length in wavelength units
% x 方向推进步长，用“波长数”来表示
% steplamb = 1 表示每一步推进 1 个波长？？

ny=2^11  % number of y points
% y 方向网格点数 = 2^11 = 2048
% 这里故意没写分号，所以会在命令窗口显示 ny 的值

aspect=4;   % ascpect ratio, delta y over delta z of compute grid
% 计算网格纵横比：Δy / Δz = 4
% 也就是 y 方向网格间距是 z 方向的 4 倍

 % numericalaspect=2;
% 旧参数或备选参数，目前未使用

nz= 2*ny/aspect   % number of z points, including image ocean
% z 方向网格点数
% 按公式 nz = 2*ny/aspect = 2*2048/4 = 1024
% 注释说明“including image ocean”：
% 这通常表示 z 网格中包含镜像区域，用于处理边界/FFT 传播等数值需要
% 同样没写分号，因此会显示 nz 的值

nout=4;  % number of outputs as it steps along, final one is at final range.
% 在传播过程中输出/保存结果的次数
% 最后一个输出对应最终传播距离处

c0=1500; % benchmark sound speed
% 基准声速，1500 m/s
% 常用于构造参考波数、无量纲化、传播算子等


% filen='test267.mat';  % input file for environmental data, this one has sound speed
 % made from enviromental data
% 旧示例：环境数据文件名（当前未使用）

% Envfile='testmoordataZvar' ;  % variable bottom depth
  % uses filen, alternative would be to have this synthesize fields
% 旧示例：环境文件（当前未使用）
% 原意是通过某个环境脚本读取 filen 并生成随距离变化的环境场


% April 2023 wedge case. Uniform in x,y SVP from test301.mat, via make_test302.m .
% 2023 年 4 月楔形算例说明：
% 使用通过 make_test302.m 生成的数据；
% 声速剖面在 x、y 上是统一背景（或基于该文件构造）

rangedep=1;
% 是否采用“距离相关（range-dependent）”环境
% 1 = 距离相关
% 0 = 距离无关

params.range_update_in_meters = 22.5;
% 环境更新间隔，单位米
% 对于距离相关环境，不必每个数值步都更新环境，隔一定距离更新一次即可
% 这里设置为每 22.5 m 更新一次


switch rangedep
% 根据 rangedep 的取值选择不同的算例配置
    
    case 1
        rdcw_flag=1; % zero for range independent
        % 距离相关声速场标志
        % 1 表示 range-dependent sound speed（cw 随 x 变化）

        runcode='layered302RangeDep_Wedge003';
        % 本次运行的名字前缀
        % 后面保存 mat 文件和图片时会用到

        filen='layered302rd' ; %  file with watercolumn sound speed
        % 距离相关水体声速文件名（不带 .mat）
        % rd = range-dependent

        Envfile='Wedge003'  % file that reads filen and makes bathymetry
        % 环境脚本名
        % 它负责读取 filen，并构造/更新 bathymetry（海底地形）
        
    case 0
        rdcw_flag=0; % zero for range independent
        % 距离无关声速场标志
        % 0 表示 range-independent

        runcode='layered302RangeIndep_Wedge003';
        % 距离无关算例的运行名

        filen='layered302ri' ; %    % file with watercolumn sound speed
        % 距离无关水体声速文件名
        % ri = range-independent

        Envfile='Wedge003'  % file that reads filen and makes bathymetry
        % 仍然使用 Wedge003 作为环境脚本
        
end


params.rng=rng;
% 把最大传播距离写入结构体 params

params.f0       =f0 ;
% 频率

params.zs       =zs; 
% 声源深度

params.depP     =depP ; 
% 绘图用水深参数

params.steplamb =steplamb;  
% x 方向步长（以波长为单位）

params.ny       =ny; 
% y 方向网格点数

params.aspect   =aspect;  
% 网格纵横比 Δy/Δz

params.nz       =nz;
% z 方向网格点数

params.nout     =nout;
% 输出次数

params.c0       =c0;   
% 基准声速

params.rdcw_flag=rdcw_flag;
% 距离相关/无关标志

params.filen    =filen;
% 环境数据文件名

params.Envfile  =Envfile;
% 环境脚本名


%% end parameters. Now run 3DPE. Save output and two plots. ........................
% ===== 参数设置结束，开始运行 3D PE，并保存结果与图像 =====
      
      
for ys= [ 500 ] % source y position, meters, can loop on this, for example, as written.
              % or modify to loop on another parameter.
% 遍历声源横向位置 ys ？？什么叫做声源横向距离
% 目前只算一个值：ys = 500 m
% 但这里设计成 for 循环，是为了以后可以方便做参数扫描
% 例如换多个声源横向位置、多频率、多环境参数等

  params.ys=ys;
  % 把当前声源横向位置写入参数结构体
  
  close all
  % 关闭所有已有图窗，避免和上一次运行的图混在一起
  
  outfile=[runcode num2str(round(ys)) ]
  % 生成输出文件名前缀，例如：
  % 'layered302RangeDep_Wedge003500'
  % 没有分号，所以会在命令窗口显示出来
      
  simulata= CARPE3D(params);
  % 调用主程序 CARPE3D
  % 输入：包含所有参数的结构体 params
  % 输出：simulata 结构体，里面应包含这次仿真的主要结果
  
  % save(outfile,'psifinal','xgg','ygg','rgg','zgg','Ez', 'Af')
  % 旧的保存方式：直接把多个独立变量存到 mat 文件
  % 现在改成保存统一的 simulata 结构体，更整洁
  
  save (outfile, 'simulata');
  % 把仿真结果结构体保存成 .mat 文件
  % 文件名就是 outfile 对应的字符串
  
  figure(3)
  % 激活编号为 3 的图窗
  % 一般 CARPE3D 内部应已画好了该图，这里只是在图上加文字并导出

  hold on
  % 保持当前图内容，允许继续往上添加文本
  
  text(0,0,runcode,'interpreter','none');
  % 在图 3 的坐标 (0,0) 处写上 runcode 字符串
  % 'interpreter','none' 表示不要把下划线等字符当作 LaTeX/TeX 解释

  print('-dpng','-r200',[runcode 'Figure3_ys' num2str(round(ys))])
  % 把图 3 导出成 PNG 图片，分辨率 200 dpi
  % 文件名类似：
  % layered302RangeDep_Wedge003Figure3_ys500.png

  pause(2)
  % 暂停 2 秒
  % 常用于确保绘图完成、文件写出稳定，或者给用户看图留出时间

  figure(4)
  % 激活编号为 4 的图窗

  hold on
  % 保持图 4 当前内容

  text(500,-zs,runcode,'interpreter','none');
  % 在图 4 上添加字符串标注
  % 标注位置为 (500, -zs)
  % 这里用 -zs，多半是为了适配该图的深度坐标显示方式

  print('-dpng','-r200',[runcode 'Figure4_ys' num2str(round(ys))])
  % 导出图 4 为 PNG 图片
end
% 一个 ys 值跑完；如果以后 ys 是多个值，会循环执行多次


diary off
% 关闭 diary 日志记录