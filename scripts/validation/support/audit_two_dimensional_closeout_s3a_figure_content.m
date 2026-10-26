function audit = audit_two_dimensional_closeout_s3a_figure_content(root,run_dir)
%AUDIT_TWO_DIMENSIONAL_CLOSEOUT_S3A_FIGURE_CONTENT Verify plotted saved-data values.
arguments
    root (1,:) char
    run_dir (1,:) char
end
assert(isfolder(root) && startsWith(run_dir,root));
d=load(fullfile(run_dir,'analysis_v2','analysis_data.mat'), ...
    'x','z','ix9','flat','rough','bres');
b=jsondecode(fileread(fullfile(run_dir,'protocol','budget.json')));
p=readtable(fullfile(run_dir,'input',sprintf('seed_%d_coefficients.csv',b.seed)));
xx=linspace(-60,60,2401)';
[eta,slope,~,raw]=closeout_s3a_pm_surface(p,xx);
figdir=fullfile(run_dir,'figures');
h=openfig(fullfile(figdir,'F01_PM_geometry_v2.fig'),'invisible');
local_match_lines(h,{raw,eta,slope},{xx,xx,xx}); close(h);
h=openfig(fullfile(figdir,'F02_nine_point_errors_v2.fig'),'invisible');
ii=d.ix9(:); x=d.x(ii); z=d.z(ii);
expected={abs(d.rough(ii,1)-d.rough(ii,2)), ...
    abs(d.rough(ii,3)-d.rough(ii,2)), ...
    angle(d.rough(ii,1).*conj(d.rough(ii,2))), ...
    angle(d.rough(ii,3).*conj(d.rough(ii,2)))};
axes_list=findall(h,'Type','Axes');
matched=false(4,1); labels=compose('P%d%d',ceil((1:9)'/3),mod((1:9)'-1,3)+1);
for k=1:numel(axes_list)
    ax=axes_list(k); sc=findall(ax,'Type','Scatter');
    if isempty(sc), continue; end
    assert(isscalar(sc) && numel(sc.XData)==9);
    assert(max(abs(sc.XData(:)-x))<1e-14 && max(abs(sc.YData(:)-z))<1e-14);
    t=findall(ax,'Type','Text'); strings=string({t.String});
    for j=1:9, assert(sum(strings==labels(j))==1,'Missing or duplicate point label'); end
    assert(sum(ismember(strings,labels))==9);
    for j=1:4
        if max(abs(sc.CData(:)-expected{j}(:)))<1e-14
            assert(~matched(j),'Duplicate scatter field'); matched(j)=true;
        end
    end
end
assert(all(matched),'A nine-point scatter field is missing'); close(h);
h=openfig(fullfile(figdir,'F03_method_line_v2.fig'),'invisible');
line=(1:411)'; x=d.x(line); br=d.rough(line,2);
values={abs(br),abs(d.rough(line,1)),abs(d.rough(line,3)), ...
    abs(d.flat(line,2)),abs(d.rough(line,1)-br), ...
    abs(d.rough(line,3)-br),d.bres(line)};
local_match_lines(h,values,repmat({x},1,numel(values))); close(h);
audit=struct('status',sprintf('S3A_SEED%d_FIGURE_CONTENT_PASS',b.seed), ...
    'F01_saved_geometry_lines',3,'F02_real_points_per_panel',9, ...
    'F02_panels',4,'F03_saved_line_series',7, ...
    'comparison_tolerance_abs',1e-14);
f=fullfile(run_dir,'audit','independent_figure_content.json');
fid=fopen(f,'w'); assert(fid>0); c=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(audit,PrettyPrint=true)); disp(audit);
end

function local_match_lines(h,values,xvalues)
lines=findall(h,'Type','Line'); matched=false(numel(values),1);
for j=1:numel(lines)
    y=lines(j).YData(:); x=lines(j).XData(:);
    for k=1:numel(values)
        if numel(y)==numel(values{k}) && ...
                max(abs(y-values{k}(:)))<1e-14 && ...
                max(abs(x-xvalues{k}(:)))<1e-14
            matched(k)=true;
        end
    end
end
assert(all(matched),'Saved-data line missing from FIG');
end
