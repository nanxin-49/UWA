function result=audit_deepsea3d_5k7k_figure()
%AUDIT_DEEPSEA3D_5K7K_FIGURE Reopen saved editable v2 figure.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
out=fullfile(root,'results','validation','deepsea3d_5k7k_20261026');
file=fullfile(out,'five_frequency_reflection_and_error_v2.fig');
fig=openfig(file,'invisible');cleanup=onCleanup(@()close(fig));
axes=findobj(fig,'Type','axes');lines=findobj(fig,'Type','line');
assert(numel(axes)==2 && numel(lines)>=6);
names=string(get(lines,'DisplayName'));
assert(any(contains(names,'Locked prediction')));
assert(any(contains(names,'5/7 kHz reference')));
assert(any(contains(names,'4/6/8 kHz reference')));
tab=readtable(fullfile(out,'five_frequency_native_summary_v2.csv'));
assert(isequal(tab.frequency_hz,[4000;5000;6000;7000;8000]));
old=lines(contains(names,'4/6/8 kHz reference'));
new=lines(contains(names,'5/7 kHz reference'));
pred=lines(contains(names,'Locked prediction'));
sourceError=max([max(abs(old.XData(:)-tab.R_real([1 3 5])*1e5)), ...
    max(abs(old.YData(:)-tab.R_imag([1 3 5])*1e5)), ...
    max(abs(new.XData(:)-tab.R_real([2 4])*1e5)), ...
    max(abs(new.YData(:)-tab.R_imag([2 4])*1e5)), ...
    max(abs(pred.XData(:)-tab.locked_pred_R_real([2 4])*1e5)), ...
    max(abs(pred.YData(:)-tab.locked_pred_R_imag([2 4])*1e5))]);
errorNames=["|R_{pred}-R|/|R|","|R_{pred}-R|/|\Delta|","|R_{pred}-R|/|H|"];
columns={"reconstruction_R_relative","reconstruction_Delta_relative","reconstruction_H_relative"};
for j=1:3
    curve=lines(names==errorNames(j));assert(isscalar(curve));
    sourceError=max(sourceError,max(abs(curve.XData(:)-[5;7])));
    sourceError=max(sourceError,max(abs(curve.YData(:)-tab.(columns{j})([2 4]))));
end
assert(sourceError<1e-12,'deepsea3d57:FigureSource','FIG coordinates differ from saved CSV.');
result=struct('status','PASS','axes',numel(axes),'line_objects',numel(lines), ...
    'maximum_source_coordinate_error',sourceError,'reopened_fig',file,'release',version);
target=fullfile(out,'figure_object_audit_v3.json');
assert(~isfile(target));
fid=fopen(target,'w');assert(fid>=0);holder=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(result,PrettyPrint=true));
disp(result);
end
