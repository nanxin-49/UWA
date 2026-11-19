function receipt = render_d1b_scale_saved(retrieved_dir, figure_dir)
%RENDER_D1B_SCALE_SAVED Render only the closed, fetched resource sample table.
arguments
    retrieved_dir (1,1) string
    figure_dir (1,1) string
end
out=fullfile(retrieved_dir,'output');
exit_info=jsondecode(fileread(fullfile(out,'exit.json')));
log=readtable(fullfile(out,'resources.csv'),TextType='string');
assert(height(log)>1 && all(diff(log.elapsed_s)>0));
assert(all(isfinite([log.elapsed_s log.total_RSS_bytes log.group_RSS_bytes log.main_registered_workers_RSS_bytes ...
    log.available_bytes log.reserve_bytes log.CPU_cores log.main_worker_CPU_cores]),'all'));
assert(all(log.group_RSS_bytes<=log.total_RSS_bytes) && all(log.main_registered_workers_RSS_bytes<=log.total_RSS_bytes));
assert(~isfolder(figure_dir),'Use a new figure output directory'); mkdir(figure_dir);
G=1024^3; t=log.elapsed_s/60;
fig=figure(Visible='off',Color='w',Position=[80 80 1150 850]); cleanup=onCleanup(@() close(fig));
layout=tiledlayout(fig,3,1,TileSpacing='compact',Padding='compact');
blue=[.12 .39 .64]; red=[.77 .26 .13]; gray=[.37 .39 .42];
axes_list=gobjects(3,1);
axes_list(1)=nexttile(layout); hold(axes_list(1),'on');
h=plot(axes_list(1),t,[log.total_RSS_bytes log.main_registered_workers_RSS_bytes log.group_RSS_bytes]/G,LineWidth=1.35);
set(h(1),Color=blue,DisplayName='Owned total'); set(h(2),Color=red,LineStyle='--',DisplayName='Main + registered workers');
set(h(3),Color=gray,LineStyle=':',DisplayName='Process group');
legend(axes_list(1),h,Location='northwest',NumColumns=3); ylabel(axes_list(1),'RSS (GiB)');
axes_list(2)=nexttile(layout); hold(axes_list(2),'on');
h=plot(axes_list(2),t,[log.available_bytes log.reserve_bytes]/G,LineWidth=1.35);
set(h(1),Color=blue,DisplayName='MemAvailable'); set(h(2),Color=red,LineStyle='--',DisplayName='Fixed reserve');
legend(axes_list(2),h,Location='best'); ylabel(axes_list(2),'Host RAM (GiB)');
axes_list(3)=nexttile(layout); hold(axes_list(3),'on');
h=plot(axes_list(3),t,[log.CPU_cores log.main_worker_CPU_cores],LineWidth=1.15);
set(h(1),Color=blue,DisplayName='Owned total'); set(h(2),Color=red,LineStyle='--',DisplayName='Main + registered workers');
legend(axes_list(3),h,Location='best'); ylabel(axes_list(3),'CPU (core equivalents)');
xlabel(axes_list(3),'Elapsed (min)');
for ax=axes_list.'
    grid(ax,'on'); xlim(ax,[t(1) t(end)]);
    for stage=["N4096" "N8192"]
        first=find(startsWith(log.stage,stage+"."),1);
        if ~isempty(first), xline(ax,t(first),':',stage,LabelVerticalAlignment='bottom',HandleVisibility='off'); end
    end
end
title(layout,"D1b scale pilot resources | "+string(exit_info.status),Interpreter='none');
drawnow;
stem=fullfile(figure_dir,'resource_curves'); savefig(fig,stem+'.fig');
exportgraphics(fig,stem+'.png',Resolution=180);
assert(isfile(stem+'.fig') && isfile(stem+'.png'));
reopened=openfig(stem+'.fig','new','invisible'); verify_cleanup=onCleanup(@() close(reopened));
lines=findobj(reopened,'Type','line'); assert(numel(lines)==7);
for line=lines.'
    assert(isequal(line.XData(:),t(:)) && numel(line.YData)==height(log));
end
receipt=struct('status','PASS_SAVED_RESOURCE_CURVES','source_run_status',exit_info.status, ...
    'samples',height(log),'lines_verified',numel(lines),'fig',stem+'.fig','png',stem+'.png');
end
