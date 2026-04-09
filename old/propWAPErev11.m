function [psiout,x,y,z,M,psi1,Ez,Af] = propWAPErev11( ...
    ny,nz,numstep,fr,c0,envfile,wid,steplength,aspect,depP,zs,ys,filen,...
    nout,rdcw_flag,range_update_steps)
% 核心宽角 PE 求解器
%我不理解wid的设置，为什么要设置镜像源，以及psi的计算公式。
% 功能概括：
%   在给定海洋环境下，把初始声场 psi(x=0,y,z) 沿 x 方向一步一步推进，
%   得到三维声传播结果。
%
% 输入参数：
%   ny, nz                y/z 方向网格点数
%   numstep               x 方向推进总步数
%   fr                    频率 (Hz)
%   c0                    参考声速 (m/s)
%   envfile               环境脚本名（通过 eval 调用）
%   wid                   横向半宽相关参数，用于定义 y 网格
%   steplength            x 方向每一步推进距离 dx
%   aspect                dy/dz 比例
%   depP                  绘图/环境参考深度参数
%   zs                    声源深度
%   ys                    声源横向位置
%   filen                 环境数据文件名（供 envfile 使用）
%   nout                  需要保存的输出次数
%   rdcw_flag             是否为距离相关环境（1=相关，0=无关）
%   range_update_steps    环境隔多少步更新一次
%
% 输出参数：
%   psiout   : 保存的输出场（按 nout 次均匀抽取）
%   x        : 最终传播距离
%   y        : y 网格
%   z        : z 网格
%   M        : 旧变量/保留接口，当前基本不用
%   psi1     : 初始场（x=0）
%   Ez       : 每一步的 z 向积分强度，形成 (x,y) 图
%   Af       : 固定在声源横向位置 ys 处的 x-z 截面场


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Parabolic equation simulation of 3-D (x,y,z) acoustic propagation.
% An initial condition psi(x_0,y,z) is marched forward using the split-step
% method to a final range psi (range, y, z). 
%
% T. Duda (WHOI)  originally vertical prop version, march 2006
%
% rev 9 has proper half step
% rev 10 loops in k domain
% rev 11 deletes obsolete variables, is set up for range-independent medium.
%    The medium can have only (y,z) variation, but this code can be modified to 
%    change the medium as it steps along (Read in from a saved enf file, for
%    example.)
%
%  psiout: output psi, nout total number, evenly spaced.
%  x: final range, propagation axis direction
%  y: transverse grid position
%  z: z grid
%  M: used for movles, not typically used, obsolete.
%  psi1: start file at range x=0
%  Ez: Flat (x,y) depth-averaged intensity, integral dz psi^2
%  Af: Complex acoustic field psi(x,z) at the source y position
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
% 开始计时

psiout = [];
% 初始化输出数组，后面按 nout 次保存中间/最终场

nnout = round(numstep*(1:nout)/nout)
% 计算应该在第几步保存输出
% 例如 nout=4，就会在总步数的约 1/4、2/4、3/4、4/4 处保存

% persistent rhosq screen Y Z atten fr
% 老版本可能用过 persistent 变量，这里现在注释掉了

lambda = c0/fr;
% 当前频率对应的参考波长

nss = 60; % plot once every nss steps
% 每 60 步画一次中间图，用来监控传播过程


% ----------- set up grid -----------------------------------------------
dx = steplength;
% x 方向步长

dy = wid/(ny-0.5);
% y 方向网格间距
% 这里 wid 近似是半宽，所以整个 y 方向跨度大约是 2*wid

dz = dy/aspect;
% z 方向网格间距
% 因为 aspect = dy/dz，所以 dz = dy/aspect

ntransz = nz;
ntransy = ny;
nstep   = numstep;
% 简单复制参数，便于后面代码沿用旧变量名

% z     =  dz/2+(0:(nz-1))*dz;
% zt    =  [-z(nz:-1:1) z];
% 老版本网格写法，当前不用

y = (-0.5*ny*dy):dy:(0.5*ny*dy)-dy;
% y 网格：以 0 为中心，左右展开

z = (-0.5*nz*dz):dz:(0.5*nz*dz)-dz;
% z 网格：也以 0 为中心
% 注意：这里的 z 坐标是数值网格坐标，和“深度正方向”的物理习惯可能不同

[Y,Z] = meshgrid(y,z);
% 构造 y-z 平面的二维网格
% 以后 psi、cin、screenm 等都是定义在这个平面上

zso1 = min(find(z > zs));
% 找到 z 网格中“高于/超过 zs”的第一个索引
% 用于确定真实声源在网格上的位置

zso2 = max(find(z < -zs));
% 找到与镜像源相关的另一个 z 索引
% 这个写法说明程序采用“真实源 + 镜像源”的起始场构造

yso = min(find(y > ys))
% 找到声源横向位置 ys 在 y 网格中对应的列索引

rhosq = (Y.^2 + Z.^2);
% 到原点的距离平方（当前版本后面没继续用）

d1 = ((Y - Y(1,yso)).^2 + (Z - Z(zso1,1)).^2);
% 网格上每个点到“真实声源位置”的距离平方

d2 = ((Y - Y(1,yso)).^2 + (Z - Z(zso2,1)).^2);
% 网格上每个点到“镜像声源位置”的距离平方
% 作者注释说它和零点是对称的

Ly = max(y) - min(y);
Lz = max(z) - min(z);
% y/z 方向的总长度，用于后面构造波数网格


% --- Green wide angle starting field from FOR3d ------------------------
k0 = 2*pi*fr/c0;
% 参考波数 k0 = 2*pi / lambda

psi = 0.5 * k0 * (1.4467 - .04201*k0*k0*d1 ) .* exp(-(k0*k0*d1)/3.0512);
% 根据真实源位置生成一个平滑的宽角起始场
% 注意这里 d1 已经是“距离平方”，所以公式里直接用 d1

psi = psi - ( 0.5 * k0 * (1.4467 - .04201*k0*k0*d2 ) .* exp(-(k0*k0*d2)/3.0512) );
% 减去镜像源对应的场
% 通常用于满足海面压力释放边界的要求
%以上psi计算是为了生成一个更适合宽角PE起步的起始声场，而不是简单在一个网格点塞delta函数

% d1 and d2 are distances already squared....

clear d1 d2 rhosq
clear psiz psiy ido dum
% 清掉临时变量和旧变量名


% ---------------------- Free Propagator set up --------------------------
ky = (2*pi/(Ly))*[0:ny/2];
% 先构造 y 方向正频率一半

ky(ntransy:-1:ny/2+1) = -ky(2:ny/2+1);
% 补成 FFT 所需的正负波数排列

ky = ones(ntransz,1)*ky;
% 扩展成与 z 方向同尺寸的二维矩阵，便于和 kz 一起做逐点运算

% kz = aspect*ky';       % make symmetric
% 旧方案：直接由 ky 推 kz，目前不用

kz = (2*pi/(Lz))*[0:nz/2];
% z 方向正频率一半

kz(ntransz:-1:nz/2+1) = -kz(2:nz/2+1);
% 补成 FFT 需要的对称波数序列

kz = ones(ntransy,1)*kz;
kz = kz';
% 扩展成二维矩阵，并转置成与 ky 尺寸匹配


% ----- run environment design code -- produces attbo, cin  -- requires Y and Z
dista = 0;
% 当前传播距离，从 0 开始

eval(envfile);
% 调用环境脚本
% 环境脚本通常会根据当前 dista、Y、Z、filen 等变量，生成：
%   cin   : 当前环境中的声速场 c(y,z)
%   attbo : 海底/介质附加衰减
%   dep   : 海底深度等信息（后面绘图会用到）

hh = max(abs(z));
la400 = 3.75;
% 400 Hz 时波长约 3.75 m lameda=c/f la400=1500/400=3.75m
% 作者注释说测试时用过 400 Hz，这里用来缩放吸收层参数


% ------------------ Attenuation function (Sponge layers all around) -----
alpha = lambda/(50*lambda/la400);%其实是个常数0.075
%也就是海绵层的"过渡宽度"大概是半域宽度的7.5%
beta  = 4*4/dx;
% 这里在 z 边界设置海绵层参数
% 目的是让波到边界时逐渐衰减，避免 FFT 周期延拓带来假反射

attb = exp( -beta*dx*exp( -((z(nz:-1:1)-hh)/(alpha*hh)).^2 ) );
% z 一侧边界的衰减函数 根据经验调好的参数，边界最外层接近e_-16，内部趋于1

attenz = attb .* fliplr(attb);
% 把上下两侧边界拼成对称衰减

attenz = attenz' * ones(1,ny);
% 扩展成 z-y 二维矩阵

hh = max(abs(y));
% 现在转到 y 边界

alpha = lambda/(50*lambda/la400);
beta  = 4*4/dx;
% y 边界的海绵层参数

att = exp( -beta*dx*exp( -((y(ny:-1:1)-hh)/(alpha*hh)).^2 ) );
% y 一侧边界衰减

atteny = att .* fliplr(att);
% y 两侧对称衰减

atteny = atteny' * ones(1,nz);
atteny = atteny';
% 扩展成与场同尺寸的二维矩阵

atten = attenz .* atteny;
% 合并 y 和 z 两个方向的海绵层

clear atteny attenz

atten = atten .* attbo;
% 再乘上环境脚本给出的附加衰减（如海底吸收等）


% --------------------- Phase screen -------------------------------------
U = (cin-c0)./cin;
% 介质相对参考声速 c0 的偏差
% 可理解为“介质引起的相位修正强度”

if max(max(isinf(U)))
    max(max(isinf(U)))
    disp('infs ... won''t run')
    return
end
% 如果 U 中出现无穷大，说明环境数据有问题，直接退出

screenm = exp(-i*k0*dx*U);
% 构造“相位屏”：
% 在实空间中，每一步传播时都乘上这个因子，
% 表示介质变化引起的相位改变

% clear U cin
% 当前不清理，后面更新环境时还要继续用

% standard PE     fr0=exp(-i*0.5*dx*(kz.^2+ky.^2)/k0); 
% 这是窄角/标准 PE 的传播因子，当前注释掉了

% WAPE Thomson Chapman 1983 JASA propagator, 2 half-steps per "range step"
kappasq = (kz.^2 + ky.^2);
% 横向总波数平方 ky^2 + kz^2

fr0 = exp(-i*0.5*dx*kappasq ./ ( (k0^2 - kappasq).^(0.5) + k0 ) );
% 宽角 PE (WAPE) 的自由传播半步算子
% 这是整个求解器最关键的传播因子之一
% 后面每个 x 步都会用它做“半步传播 + 相位屏 + 半步传播”

disp('____________ WAPE rev 10 _________________')
clear kappasq

% dista=0;  % moved to before eval envfile
% 老注释：dista 的初始化已经提前到环境脚本之前了

psi1 = psi;
% 保存初始场，作为输出之一

% [j1,j2]=size(screen);
% midl=find( abs(yin)< 1e-8 );
% nrng=(midl-(ny/2) ):(midl+ (ny/2)-1);
% 旧变量/旧逻辑，当前不用

Ez = ones(nstep,ny);
% 预分配 Ez
% 每一步会保存一个长度为 ny 的向量，表示 z 向积分强度随 y 的分布

nzhalf = nz/2;
% 后面只取 z 的一半用于输出/剖面

Af = ones(nzhalf,nstep);
% 预分配 Af
% 每一步保存一个 z 向剖面（固定在 y=ys）

figure(1); clf
% 用图窗 1 动态显示传播过程

tempout = psi;
% 临时变量，后面会用于保存当前步的实空间场

chat_interval = 2*round(.5*nstep/50);
% 每大约 2% 的进度打印一次进展信息


if(rdcw_flag)
    disp([ 'Environment update interval '   num2str(range_update_steps) ' steps'])
    disp([ 'Environment update interval '   num2str(dx*range_update_steps) ' meters'])
else
    disp('* Range (x) independent environment *')
end
% 如果环境是距离相关的，就告诉你环境每隔多少步/多少米更新一次
% 如果不是，就说明是 x 无关环境


for jj = 1:nstep   % main prop loop --------------------------------------
    % 主传播循环：沿 x 方向推进 numstep 步
    
    % screenm=screen(:,nrng);    
    % 旧代码残留
    
    % medium fixed in march-direction x in this code, need something in
    % this loop to incorporate x-dependent enviromnent. 
    % Can call the  envfile function as one option , i.e.
    %        eval(envfile) 
    % where different things happen for differing jj
    
    % --------------------- Phase screen renewal --------------------------
    % 如果是距离相关环境，则每隔 range_update_steps 步更新一次环境
    if(mod(jj-1,range_update_steps)==0 & jj>1)
        % 第一步已经在 x=0 计算过环境了，所以从 jj>1 开始更新
        
        eval(envfile)
        % 重新读取/生成当前 x 位置的环境
        
        U = (cin-c0)./cin;
        % 新环境对应的新介质偏差
        
        if max(max(isinf(U)))
            max(max(isinf(U)))
            disp('infs ... won''t run')
            return
        end
        
        screenm = exp(-i*k0*dx*U);
        % 更新相位屏
    end
    
    
    % if(mod(jj,25) ==0)
    %     disp('randomize seafloor')
    %     eval(envfile);  %randomize bottom 
    % end  
    % 曾经可能尝试过随机海底，每 25 步刷新一次，这里注释掉了
    
    
    if(jj==1 | mod(jj,nss)==0)  % plot every nss steps
        wat = find( z > -(max(dep(:))+15) & z < .5);
        % 选择一个 z 范围来画图：
        % 大致是从海底以下一点点到接近海面
        
        subplot(211)
        pcolor(y(1:2:end), z(wat), abs(tempout(wat,1:2:end))); shading flat;
        % 上半图：当前步实空间场 |psi(y,z)|
        % y 每隔一个点抽样一次，减少绘图负担
        
        % caxis([ .4 1.2]);
        caxis([ max(max(abs(tempout)))*.1 max(max(abs(tempout))) ])
        % 颜色范围设为当前最大值的 10% 到 100%
        
        colorbar('southoutside')
        title('|field(y,z)|')
        
        subplot(212)
        % pcolor(y(1:2:end),z(wat),angle(screenm(wat,1:2:end)));shading flat
        % 原本也可以画相位屏的相位
        
        pcolor(y(1:2:end), z(wat), (U(wat,1:2:end))); shading flat
        % 下半图：当前介质参数 U(y,z)，用来检测海洋介质是否按要求更新
        
        xlabel('y (m)')
        ylabel('z (m)')
        % axis equal
        
        title([ 'Medium exp(-ik0dxU) @ x= ' num2str(dista) ' m'] );
        title([ 'Medium U @ x= ' num2str(dista) ' m'] );
        % 第二个 title 覆盖第一个，所以最终显示的是 U
        
        colorbar('southoutside')
        caxis([ -.005 0 ])
        pause(0.4)
        % 暂停一下，便于动态看传播过程
    end
   
    
    % rev 10 prop --------------------------------------------------------
    if(jj == 1)
        psi = fft2(psi);
        % 第一步先把起始场变换到二维波数域
    end
    
    % psi=ifft2(fr0.*fft2(atten.*screenm.*ifft2(fr0.*fft2(psi))));
    % rev 9 的写法，当前作者换成更直接的“psi 已在频域”写法
    
    psi = fr0 .* fft2(atten .* screenm .* ifft2(fr0 .* psi));
    % 这是整段代码最核心的一句：
    %
    % 1) 先在频域乘 fr0           -> 半步自由传播
    % 2) ifft2 回实空间           -> 回到 y-z 平面
    % 3) 乘 atten .* screenm      -> 施加边界吸收 + 介质相位屏
    % 4) fft2 回到频域            -> 再回波数空间
    % 5) 再乘一个 fr0             -> 第二个半步自由传播
    %
    % 这就是典型 split-step 宽角 PE 结构：
    %   半步传播 + 介质屏幕 + 半步传播

    tempout = ifft2(psi);
    % 把当前步的场变回实空间，便于计算强度、取剖面、画图
    
    Ez(jj,:) = sum(abs(tempout).^2);
    % 沿 z 方向把 |psi|^2 积分/求和，得到当前 x 处的 (y) 强度分布
    % 这样 Ez 最终就是一个 nstep × ny 的矩阵，可画成 x-y 图

    Af(:,jj) = tempout(1:nzhalf,yso);
    % 取固定在 y=ys 这一列的 x-z 剖面
    % 只取 z 的一半，通常是实际水体部分

    if(jj == nstep)
        psi = tempout;
        % 最后一步结束时，把 psi 保持成实空间场
        % 这样最终返回时更直观
    end
    % -------------------------------------------------------------------
    
    
    dista = dista + dx;
    % 更新当前传播距离
    
    seco = toc;
    % 当前累计用时
    
    if(mod(jj,chat_interval)==0)
        disp(['to ' num2str(dista) ' m ; percent ' num2str(100*jj/nstep) ...
              ';  time: ' num2str(seco) ' sec'])
    end
    % 每隔一定步数打印进度
    
    
    ccc = jj - nnout;
    % 看当前步 jj 是否正好等于需要保存的某个输出步
    
    if(min(abs(ccc)) == 0)
        nuu = find(ccc==0)
        % nuu 表示这是第几个输出
        
        if(nuu == 1)
           psiout = shiftdim(psi(1:(1+nz/2),:), -1);
           % 第一次保存时初始化 psiout
           % 这里只保存 z 的前半部分
        else
           psiout(nuu,:,:) = psi(1:(1+nz/2),:);
           % 后续继续往 psiout 里追加
        end
    end
end  % main prop loop


x = dista;
% 最终传播距离

% psifinal=psi;
% 老变量名注释，当前直接通过输出返回

M = 'off';
% 老接口保留，给做 movie 用；当前基本不重要

% whos
% 调试用，当前不执行