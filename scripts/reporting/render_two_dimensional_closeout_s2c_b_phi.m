function render_two_dimensional_closeout_s2c_b_phi(mat_path,out_dir)
%RENDER_TWO_DIMENSIONAL_CLOSEOUT_S2C_B_PHI Saved-data M05b-e figure family.
% No acoustic evaluation; nine points stay nine pixels/cells, never interpolated.
arguments
    mat_path (1,:) char
    out_dir (1,:) char
end
assert(isfile(mat_path) && ~isfolder(out_dir));
d=load(mat_path,'U','x','z','ix9','physical','point','support','paired');
assert(isequal(size(d.U),[417 2 3]) && numel(d.x)==417 && numel(d.z)==417);
assert(all(d.z(1:411)==3) && numel(d.ix9)==9);
assert(all(isfinite(real(d.U)),'all') && all(isfinite(imag(d.U)),'all'));
assert(all(ismember(["BASE_M","PHI_PLUS90"],string(unique(d.support.case_name)))));
run_dir=fileparts(fileparts(mat_path));
mkdir(out_dir);
copyfile(fullfile(run_dir,'figures_final','M05a_geometry_only.png'),out_dir);
copyfile(fullfile(run_dir,'figures_final','M05a_geometry_only.fig'),out_dir);
x=d.x(1:411); U=d.U; blue=[0 .36 .7]; red=[.84 .27 .12];
% M05b: actual BIE physical response, separate from model errors.
f=newfig(1680,1500); t=tiledlayout(f,3,3,TileSpacing='compact',Padding='compact');
b0=U(1:411,1,2); b1=U(1:411,2,2); dphi=b1-b0;
a=nexttile(t); hold(a,'on'); plot(a,x,abs(b0),Color=blue,LineWidth=1.3);
plot(a,x,abs(b1),Color=red,LineWidth=1.3); ylabel(a,'|u_{BIE}| [unit source]');
title(a,'原生反射场幅度'); legend(a,{'BASE_M','PHI_PLUS90'},Location='best');
a=nexttile(t); hold(a,'on'); phaseplot(a,x,angle(b0),blue);
phaseplot(a,x,angle(b1),red); ylim(a,[-pi pi]); ylabel(a,'wrapped arg(u_{BIE}) [rad]');
title(a,'原生相位；wrap处断线'); legend(a,{'BASE_M','PHI_PLUS90'},Location='best');
a=nexttile(t); hold(a,'on'); plot(a,x,real(dphi),Color=blue,LineWidth=1.1);
plot(a,x,imag(dphi),Color=red,LineWidth=1.1);
ylabel(a,'\Delta_\phi u_{BIE} [unit source]'); title(a,'复变化的实部与虚部');
legend(a,{'real','imag'},Location='best');
a=nexttile(t); plot(a,x,abs(dphi),'k',LineWidth=1.2);
ylabel(a,'|\Delta_\phi u_{BIE}| [unit source]'); title(a,'BIE物理响应；不是模型误差');
a=nexttile(t); plot(a,x,20*log10(abs(b1)./abs(b0)),Color=red,LineWidth=1.2);
ylabel(a,'PHI/BASE amplitude [dB]'); title(a,'同点幅度变化');
a=nexttile(t); phaseplot(a,x,angle(b1.*conj(b0)),red); ylim(a,[-pi pi]);
ylabel(a,'wrapped \Delta_\phi phase [rad]'); title(a,'同点相位变化；wrap处断线');
center=abs(x)<=1;
a7=nexttile(t,7); hold(a7,'on'); plot(a7,x(center),abs(b0(center)),Color=blue,LineWidth=1.3);
plot(a7,x(center),abs(b1(center)),Color=red,LineWidth=1.3);
ylabel(a7,'|u_{BIE}| [unit source]'); title(a7,'中心11点：原生幅度');
a8=nexttile(t,8); hold(a8,'on'); phaseplot(a8,x(center),angle(b0(center)),blue);
phaseplot(a8,x(center),angle(b1(center)),red); ylim(a8,[-pi pi]);
ylabel(a8,'wrapped phase [rad]'); title(a8,'中心11点：原生相位');
a9=nexttile(t,9); plot(a9,x(center),abs(dphi(center)),'k',LineWidth=1.2);
ylabel(a9,'|\Delta_\phi u_{BIE}| [unit source]'); title(a9,'中心11点：物理响应');
title(t,'M05b | 只改初相位；BIE参考场实际响应 | 411个z=3 m真点');
finish(f,t,x); xlim(a7,[-1 1]); xlim(a8,[-1 1]); xlim(a9,[-1 1]);
store(f,out_dir,'M05b_BIE_physical_response');
% M05c: two cases, same per-method axes and exact amplitude/phase squares.
f=newfig(1740,1600); t=tiledlayout(f,4,2,TileSpacing='compact',Padding='compact');
for col=1:2
    m=[1 3]; m=m(col); name={'PE-mode 0','FullFK complete'}; name=name{col};
    ydb=zeros(411,2); yph=ydb; yr=ydb; ya=ydb; yp=ydb;
    for c=1:2
        u=U(1:411,c,m); b=U(1:411,c,2); ph=angle(u.*conj(b));
        ydb(:,c)=20*log10(abs(u)./abs(b)); yph(:,c)=ph;
        yr(:,c)=abs(u-b); ya(:,c)=(abs(u)-abs(b)).^2;
        yp(:,c)=4*abs(u).*abs(b).*sin(ph/2).^2;
    end
    dbmax=1.05*max(abs(ydb),[],'all');
    phmax=1.05*max(abs(yph),[],'all');
    rv=yr(yr>0); sv=[ya(ya>0);yp(yp>0)];
    rlim=[10^floor(log10(min(rv))) 10^ceil(log10(max(rv)))];
    slim=[10^floor(log10(min(sv))) 10^ceil(log10(max(sv)))];
    colors={blue,red};
    a=nexttile(t,col); hold(a,'on'); for c=1:2, plot(a,x,ydb(:,c),Color=colors{c},LineWidth=1.15); end
    ylim(a,[-dbmax dbmax]); ylabel(a,'amplitude diff [dB]');
    if col==2
        % Export can auto-select different ticks for this small signed range.
        % Fix only the display ticks; limits and plotted values are unchanged.
        a.YTick=(-2:2)*1e-5;
    end
    title(a,[name ' | signed amplitude difference']); legend(a,{'BASE_M','PHI_PLUS90'},Location='best');
    a=nexttile(t,2+col); hold(a,'on'); for c=1:2, phaseplot(a,x,yph(:,c),colors{c}); end
    ylim(a,[-phmax phmax]); ylabel(a,'wrapped phase diff [rad]'); title(a,'相对同案BIE相位；wrap处断线');
    a=nexttile(t,4+col); hold(a,'on');
    for c=1:2, semilogy(a,x,max(yr(:,c),realmin),Color=colors{c},LineWidth=1.2); end
    a.YScale='log'; ylim(a,rlim); ylabel(a,'|u_m-u_{BIE}| [unit source]'); title(a,'绝对复残差；显示下限不改原值');
    a=nexttile(t,6+col); hold(a,'on');
    for c=1:2
        semilogy(a,x,max(ya(:,c),realmin),Color=colors{c},LineStyle='-',LineWidth=1.1);
        semilogy(a,x,max(yp(:,c),realmin),Color=colors{c},LineStyle='--',LineWidth=1.1);
    end
    a.YScale='log'; ylim(a,slim); ylabel(a,'square terms [unit source^2]');
    title(a,'|e|^2 = amplitude^2 + phase^2');
    legend(a,{'BASE amplitude','BASE phase','PHI amplitude','PHI phase'},Location='best');
end
title(t,'M05c | 两案同源模型误差；每列两案例同尺度；列间不同纵轴；无拟合');
finish(f,t,x); store(f,out_dir,'M05c_model_error_pair');
% M05d: full line, center, and genuine 3x3 receiver coordinates.
f=newfig(1650,1580); t=tiledlayout(f,4,2,TileSpacing='compact',Padding='compact');
method_col=[0 .38 .70;.15 .55 .28];
yfull=abs(U(1:411,:,[1 3])-U(1:411,:,2));
yfull=yfull(:); yfull=yfull(yfull>0);
yl=[10^floor(log10(min(yfull))),10^ceil(log10(max(yfull)))];
for c=1:2
    ca={'BASE_M','PHI_PLUS90'}; ca=ca{c};
    a=nexttile(t,c); hold(a,'on');
    for j=1:2
        m=[1 3]; m=m(j); e=abs(U(1:411,c,m)-U(1:411,c,2));
        semilogy(a,x,max(e,realmin),Color=method_col(j,:),LineWidth=1.25);
    end
    a.YScale='log'; ylim(a,yl); xlim(a,[x(1) x(end)]); ylabel(a,'|e| [unit source]');
    title(a,[ca ' | 411点整线 q=1/411']); legend(a,{'PE-BIE','FullFK-BIE'},Location='best');
    a=nexttile(t,2+c); hold(a,'on'); ii=find(abs(x)<=1);
    for j=1:2
        m=[1 3]; m=m(j); e=abs(U(ii,c,m)-U(ii,c,2));
        semilogy(a,x(ii),max(e,realmin),Color=method_col(j,:),LineWidth=1.25);
    end
    a.YScale='log'; ylim(a,yl); xlim(a,[-1 1]); ylabel(a,'|e| [unit source]');
    title(a,[ca ' | 中心11点，|x|<=1 m']);
end
for j=1:2
    m=[1 3]; m=m(j);
    v=zeros(9,2);
    for c=1:2
        r=U(d.ix9,c,2); u=U(d.ix9,c,m); v(:,c)=abs(u-r)./abs(r);
    end
    scale=max(v,[],'all');
    for c=1:2
        a=nexttile(t,4+2*(j-1)+c);
        imagesc(a,1:3,1:3,reshape(v(:,c),3,3)');
        a.YDir='reverse'; clim(a,[0 scale]); colormap(a,parula(128));
        a.XTick=1:3; a.XTickLabel=compose('%.6g',unique(d.x(d.ix9)));
        a.YTick=1:3; a.YTickLabel=compose('%.3f',unique(d.z(d.ix9)));
        xlabel(a,'receiver x [m]'); ylabel(a,'depth z [m]');
        cb=colorbar(a); cb.Label.String='|e|/|u_{BIE}| [1]';
        ca={'BASE_M','PHI_PLUS90'}; ca=ca{c};
        md={'PE','FullFK'}; md=md{j};
        row=d.support(d.support.case_name==string(ca)&d.support.method==string(md) ...
            &d.support.support_id=="nine_9",:);
        if row.resolved_above_5x, flag='resolved >5x'; else, flag='UNRESOLVED'; end
        title(a,sprintf('%s %s | E_9=%.3g | %s | worst (%.3g,%.3g) m', ...
            ca,md,row.complex_L2,flag,row.worst_x_m,row.worst_z_m));
        for p=1:9
            text(a,mod(p-1,3)+1,ceil(p/3),sprintf('P%d%d\n%.2g',ceil(p/3),mod(p-1,3)+1,v(p,c)), ...
                HorizontalAlignment='center',FontSize=9,Color='k',BackgroundColor='w',Margin=.3);
        end
    end
end
title(t,'M05d | 共用坐标/色标：整线、中心11点、真实九点；六个离深度点未投影');
style(f); store(f,out_dir,'M05d_line_center_nine');
% M05e: mirror pairs use the actual same-depth receiver coordinates.
[paired,mirror]=ismember([-d.x(:),d.z(:)],[d.x(:),d.z(:)],'rows');
assert(all(paired) && all(mirror>0));
positive=find(x>0); opposite=mirror(positive); assert(all(opposite<=411));
symmetry=table; cases=["BASE_M","PHI_PLUS90"]; methods=["PE","BIE","FullFK"];
for c=1:2
    for m=1:3
        native=U(:,c,m); residual=native-U(:,c,2);
        symmetry=[symmetry;table(repmat(cases(c),417,1),repmat(methods(m),417,1), ...
            (1:417)',d.x(:),d.z(:),mirror,real(native),imag(native), ...
            real(native(mirror)),imag(native(mirror)), ...
            real(native-native(mirror)),imag(native-native(mirror)), ...
            real(residual-residual(mirror)),imag(residual-residual(mirror)), ...
            'VariableNames',{'case_name','method','sample_index','x_m','z_m', ...
            'mirror_index','u_real','u_imag','mirror_u_real','mirror_u_imag', ...
            'native_mirror_difference_real','native_mirror_difference_imag', ...
            'residual_mirror_difference_real','residual_mirror_difference_imag'})]; %#ok<AGROW>
    end
end
assert(height(symmetry)==2502);
writetable(symmetry,fullfile(out_dir,'M05e_mirror_pairs.csv'));
f=newfig(1650,1500); t=tiledlayout(f,3,2,TileSpacing='compact',Padding='compact');
for m=1:3
    a=nexttile(t,2*m-1); hold(a,'on');
    for c=1:2
        y=abs(U(positive,c,m)-U(opposite,c,m));
        col=blue; if c==2, col=red; end
        plot(a,x(positive),y,Color=col,LineWidth=1.15);
    end
    ylabel(a,'|u(+x)-u(-x)| [unit source]');
    title(a,[char(methods(m)) ' 原生复场左右差']);
    legend(a,{'BASE_M','PHI_PLUS90'},Location='best');
    a=nexttile(t,2*m); hold(a,'on');
    if m==2
        delta=U(:,2,2)-U(:,1,2);
        plot(a,x(positive),abs(delta(positive)-delta(opposite)),'k',LineWidth=1.15);
        ylabel(a,'|\Delta_\phi u_{BIE}(+x)-\Delta_\phi u_{BIE}(-x)|');
        title(a,'BIE 物理响应的左右差');
    else
        for c=1:2
            e=U(:,c,m)-U(:,c,2); col=blue; if c==2, col=red; end
            plot(a,x(positive),abs(e(positive)-e(opposite)),Color=col,LineWidth=1.15);
        end
        ylabel(a,'|e(+x)-e(-x)| [unit source]');
        title(a,[char(methods(m)) ' 模型残差左右差']);
        legend(a,{'BASE_M','PHI_PLUS90'},Location='best');
    end
end
title(t,'M05e | 同深度真实\pm x配对；输入近似对称不等于声场严格对称或因果证明');
style(f); aa=findall(f,'Type','axes');
for j=1:numel(aa)
    xlim(aa(j),[0 max(x)]); xlabel(aa(j),'receiver +x [m]'); grid(aa(j),'on');
end
store(f,out_dir,'M05e_mirror_symmetry');
end

function f=newfig(w,h)
f=figure(Visible='off',Color='w',Units='pixels',Position=[40 40 w h]);
end
function phaseplot(a,x,y,col)
jump=abs(diff(y))>pi; pos=(1:numel(x))'+cumsum([0;jump(:)]);
xx=nan(numel(x)+sum(jump),1); yy=xx; xx(pos)=x; yy(pos)=y;
plot(a,xx,yy,Color=col,LineWidth=1.15);
end
function finish(f,t,x)
style(f); aa=findall(f,'Type','axes');
for j=1:numel(aa)
    xlim(aa(j),[x(1) x(end)]); xline(aa(j),0,':',Color=[.55 .55 .55],HandleVisibility='off');
    grid(aa(j),'on'); if isempty(aa(j).XLabel.String), xlabel(aa(j),'receiver x [m]'); end
end
t.TileSpacing='compact';
end
function style(f)
aa=findall(f,'Type','axes');
for j=1:numel(aa), aa(j).FontName='Microsoft YaHei'; aa(j).FontSize=10; aa(j).Box='on'; end
tt=findall(f,'Type','text'); for j=1:numel(tt), tt(j).FontName='Microsoft YaHei'; end
end
function store(f,out_dir,stem)
drawnow; savefig(f,fullfile(out_dir,[stem '.fig']));
exportgraphics(f,fullfile(out_dir,[stem '.png']),Resolution=180); close(f);
end
