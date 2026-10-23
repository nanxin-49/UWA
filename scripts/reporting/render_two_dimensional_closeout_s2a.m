function files = render_two_dimensional_closeout_s2a(plot_file,out_dir)
%RENDER_TWO_DIMENSIONAL_CLOSEOUT_S2A Redraw from retained plot_data only.
arguments
    plot_file (1,:) char
    out_dir (1,:) char
end
assert(isfile(plot_file),'Missing saved plot data.');
assert(~isfolder(out_dir),'Refusing to overwrite figure directory.');
t=load(plot_file,'plot_data'); d=t.plot_data;
assert(numel(d.x_m)==411 && numel(d.PE)==4 && numel(d.BIE)==4 && numel(d.FullFK)==4);
mkdir(out_dir);
x=d.x_m(:); tags={'flat','weak A=.01','medium A=.05','strong A=.20'};
colors=[0 .35 .75; .77 .2 .1; .2 .55 .25];
methods={'PE','BIE','FullFK'};
files=cell(3,1);

f=local_figure(); tiledlayout(f,2,2,'Padding','compact','TileSpacing','compact');
nexttile; hold on
for m=1:3, plot(x,abs(d.(methods{m}){1}),'Color',colors(m,:),'LineWidth',1.2); end
title('F02 · flat reflected amplitude'); ylabel('|u_{flat}| (unit source)');
legend(methods,'Location','best'); grid on
nexttile; hold on
for m=1:3, plot(x,local_phase(d.(methods{m}){1}),'Color',colors(m,:)); end
title('Flat native phase'); ylabel('wrapped phase (rad)'); grid on
nexttile; hold on
for m=[1 3]
    a=d.(methods{m}){1}; b=d.BIE{1};
    plot(x,abs(a-b),'Color',colors(m,:));
end
title('Flat complex residual vs BIE'); xlabel('receiver x (m)'); ylabel('|u - u_{BIE}|'); grid on
nexttile; hold on
for m=[1 3]
    a=d.(methods{m}){1}; b=d.BIE{1};
    plot(x,angle(a.*conj(b)),'Color',colors(m,:));
end
title('Flat phase residual vs BIE'); xlabel('receiver x (m)'); ylabel('Δφ (rad)'); grid on
files{1}=local_save(f,out_dir,'F02_flat_control');

f=local_figure(); tl=tiledlayout(f,3,2,'Padding','compact','TileSpacing','compact');
for j=2:4
    nexttile(tl,2*j-3); hold on
    for m=1:3, plot(x,abs(d.(methods{m}){j}),'Color',colors(m,:),'LineWidth',1); end
    title([tags{j} ' · |u_{rough}|']); ylabel('unit source'); grid on
    if j==4, xlabel('receiver x (m)'); end
    nexttile(tl,2*j-2); hold on
    for m=1:3, plot(x,local_phase(d.(methods{m}){j}),'Color',colors(m,:)); end
    title([tags{j} ' · native phase']); ylabel('rad'); grid on
    if j==4, xlabel('receiver x (m)'); end
end
legend(tl.Children(end),methods,'Location','best');
files{2}=local_save(f,out_dir,'F04_rough_native_lines');

f=local_figure(); tl=tiledlayout(f,3,3,'Padding','compact','TileSpacing','compact');
for j=2:4
    a=d.PE{j}; b=d.BIE{j}; fk=d.FullFK{j};
    nexttile(tl,3*j-5); hold on
    plot(x,abs(a-b),'Color',colors(1,:)); plot(x,abs(fk-b),'Color',colors(3,:));
    title([tags{j} ' · |Δu| vs BIE']); ylabel('unit source'); grid on
    nexttile(tl,3*j-4); hold on
    plot(x,abs(a./d.PE{1}-b./d.BIE{1}),'Color',colors(1,:));
    plot(x,abs(fk./d.FullFK{1}-b./d.BIE{1}),'Color',colors(3,:));
    title('|G - G_{BIE}| (auxiliary)'); ylabel('dimensionless'); grid on
    nexttile(tl,3*j-3); hold on
    plot(x,abs((a-d.PE{1})-(b-d.BIE{1})),'Color',colors(1,:));
    plot(x,abs((fk-d.FullFK{1})-(b-d.BIE{1})),'Color',colors(3,:));
    title('|δu - δu_{BIE}| (auxiliary)'); ylabel('unit source'); grid on
    if j==4
        for k=1:3, nexttile(tl,6+k); xlabel('receiver x (m)'); end
    end
end
nexttile(tl,1);
legend({'PE–BIE','FullFK–BIE'},'Location','northeast');
files{3}=local_save(f,out_dir,'F05_native_residuals');
end

function f=local_figure()
f=figure('Visible','off','Color','w','Units','pixels','Position',[100 100 1500 850]);
end

function p=local_phase(u)
p=angle(u); jump=[false;abs(diff(p))>pi]; p(jump)=NaN;
end

function p=local_save(f,out_dir,name)
drawnow; p=fullfile(out_dir,[name '.png']);
savefig(f,fullfile(out_dir,[name '.fig']));
print(f,p,'-dpng','-r140'); close(f);
end
