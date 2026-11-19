function audit = plot_p1_ds_control(data, stem, caption)
% Render saved control data, then audit the serialized figure independently.
x = data.x(:);
series = {abs(data.pS-data.refS)./data.ES, ...
    abs(data.pDS-data.refDS)./data.EDS, abs(data.refDS), ...
    abs(data.pDS), data.EDS, .01*data.EDS};
labels = {'S complex error / E_S','D+S complex error / E', ...
    'Analytic D+S','Bellhop D+S','Independent E','0.01 E phase threshold'};
styles = {'b-o','r-s','k-o','b--s','r:','g-.'};
fig = figure('Visible','off','Color','w','Position',[50 50 1280 560], ...
    'PaperPositionMode','auto');
closeSource = onCleanup(@() close(fig));
ax(1) = axes('Parent',fig,'Position',[.075 .27 .39 .62]);
ax(2) = axes('Parent',fig,'Position',[.575 .27 .39 .62],'YScale','log');
for j=1:2
    hold(ax(j),'on');
    grid(ax(j),'on');
    xlabel(ax(j),'Independent receiver offset x (m)');
end
ylabel(ax(1),'Normalized complex error');
ylabel(ax(2),'Pressure magnitude / 1 m reference');
for j=1:6
    panel = 1 + (j>2);
    h(j) = plot(ax(panel),x,series{j},styles{j}, ...
        'Tag',sprintf('audit_curve_%d',j),'DisplayName',labels{j}, ...
        'LineWidth',1.4);
end
yline(ax(1),.002,'k:','Frozen gate','HandleVisibility','off');
legend(ax(1),h(1:2),labels(1:2),'Location','southoutside', ...
    'NumColumns',2,'AutoUpdate','off','FontSize',10);
legend(ax(2),h(3:6),labels(3:6),'Location','southoutside', ...
    'NumColumns',2,'AutoUpdate','off','FontSize',10);
annotation(fig,'textbox',[.025 .923 .95 .06],'String',caption, ...
    'HorizontalAlignment','center','EdgeColor','none','FontSize',14);
drawnow;
savefig(fig,[stem '.fig']);
print(fig,[stem '.png'],'-dpng','-r135');
clear closeSource
reopened = openfig([stem '.fig'],'new','invisible');
closeReopened = onCleanup(@() close(reopened)); %#ok<NASGU>
drawnow;
for j=1:6
    line = findobj(reopened,'Tag',sprintf('audit_curve_%d',j));
    assert(isscalar(line) && isequal(line.XData(:),x) && ...
        isequal(line.YData(:),series{j}(:)),'Reopened FIG curve data');
end
legends = findobj(reopened,'Type','legend');
assert(numel(legends)==2,'Reopened FIG legend count');
bounds = zeros(2,4);
counts = zeros(1,2);
for j=1:2
    legends(j).Units='normalized';
    p=legends(j).Position;
    bounds(j,:)=p;
    assert(all(p(1:2)>=.005) && all(p(1:2)+p(3:4)<=.995) && ...
        strcmp(legends(j).AutoUpdate,'off'),'Reopened FIG legend bounds and AutoUpdate');
    strings = cellstr(string(legends(j).String));
    counts(j)=numel(strings);
    assert(isequal(strings(:),labels(1:2)') || ...
        isequal(strings(:),labels(3:6)'),'Reopened FIG legend labels');
end
assert(isequal(sort(counts),[2 4]));
info=imfinfo([stem '.png']);
pixels=imread([stem '.png']);
assert(info.Width>=1200 && info.Height>=500 && std(double(pixels(:)))>1, ...
    'PNG decodes and contains rendered content');
audit=struct('status','FIG_AND_PNG_FILE_CHECKS_PASS','curves',6, ...
    'legend_bounds_normalized',bounds,'legend_auto_update','off', ...
    'png_size',[info.Width,info.Height]);
fid=fopen([stem '_audit.json'],'w','n','UTF-8');
assert(fid>=0);
closeAudit=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(audit,PrettyPrint=true));
end
