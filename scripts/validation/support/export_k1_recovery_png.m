function receipt=export_k1_recovery_png(f,path_png)
%EXPORT_K1_RECOVERY_PNG Fixed paper-size print, with unchanged canvas gate.
% A new plotting backend only. Saved arrays and curve semantics are untouched.
assert(isgraphics(f,'figure') && ~isfile(path_png));
assert(strcmp(f.Units,'pixels') && isequal(f.Position(3:4),[1800 1400]));
before=f.Position;
if ispc
    % Capture the actual fixed canvas without invoking print's restoration.
    % Visibility is restored; figure position is never reset after capture.
    visibility=f.Visible; guard=onCleanup(@()set(f,'Visible',visibility));
    f.Visible='on'; figure(f); drawnow;
    assert(isequal(f.Position,before),'k1recovery:Canvas','Showing figure changed canvas');
    frame=getframe(f); assert(isequal(size(frame.cdata,1),1400) && isequal(size(frame.cdata,2),1800));
    imwrite(frame.cdata,path_png); clear guard;
    backend='MATLAB native getframe, unscaled 1800x1400 canvas';
else
    % Headless Linux retains the original K1 exporter and all strict gates.
    exportgraphics(f,path_png,'Resolution',160,'Width',1800,'Height',1400,'Units','pixels');
    backend='MATLAB exportgraphics, 1800x1400';
end
assert(isequal(f.Position,before),'k1recovery:Canvas','Exporter changed figure canvas');
info=imfinfo(path_png); assert(info.Width==1800 && info.Height==1400);
receipt=struct('backend',backend,'width_pixels',info.Width,'height_pixels',info.Height,'canvas_unchanged',true);
end
