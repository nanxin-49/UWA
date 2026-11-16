function optimize_bellhop_rotated_flatwall_overlay(source_dir,output_dir)
%OPTIMIZE_BELLHOP_ROTATED_FLATWALL_OVERLAY Conservative cull; exact roots kept.
arguments
    source_dir (1,:) char
    output_dir (1,:) char
end
assert(~isfolder(output_dir)); mkdir(output_dir); files=dir(fullfile(source_dir,'*.f90'));
for k=1:numel(files), copyfile(fullfile(files(k).folder,files(k).name),fullfile(output_dir,files(k).name)); end
step=strrep(fileread(fullfile(output_dir,'Step.f90')),sprintf('\r\n'),newline);
step=once(step,'  INTEGER :: WallNProfile = 0',"  LOGICAL :: WallIsExactVertical=.FALSE."+newline+"  INTEGER :: WallNProfile = 0");
step=once(step,'    WallX( 2, : ) = Z',"    WallX( 2, : ) = Z"+newline+"    WallIsExactVertical=ALL(R==R(1))");
% Bounding slab only rejects segments farther than original closure/debug gates.
% It never changes the accepted root calculation, segment order or tie breaking.
step=once(step,'    tol = 1.0d-10',join(["    tol = 1.0d-10"; ...
"    IF(WallIsExactVertical) THEN"; ...
"       IF(MAX(x0(1),x0(1)+hTrial*u(1))<WallX(1,1)-1d-7.OR. &"; ...
"          MIN(x0(1),x0(1)+hTrial*u(1))>WallX(1,1)+1d-7) RETURN"; ...
"    END IF"],newline));
anchor='    CALL GetWallGeometry( ray2%x, iWallSeg, wallLambda, wallT, wallN, wallKappa, wallResidual )';
step=once(step,anchor,join([ ...
"    iWallSeg=0; wallLambda=0d0; wallT=0d0; wallN=0d0; wallKappa=0d0; wallResidual=HUGE(wallResidual)"; ...
"    IF(WallActive) THEN"; ...
"       IF(WallIsExactVertical) THEN"; ...
"          IF(ABS(ray2%x(1)-WallX(1,1))<=1d-3) THEN";string(anchor); ...
"          END IF";"       ELSE";string(anchor);"       END IF";"    END IF"],newline));
write(fullfile(output_dir,'Step.f90'),step);
if isfile(fullfile(source_dir,'build_performance.ps1')), builder='build_performance.ps1'; else, builder='build_rotated_receiver64.ps1'; end
copyfile(fullfile(source_dir,builder),fullfile(output_dir,'build_optimized.ps1'));
% The REAL8 source modules and field evaluator must be byte-identical.
for name={'SourceReceiverPositions.f90','ReadEnvironmentBell.f90','RotatedFieldState.f90','influence.f90'}
    assert(strcmpi(black_file_hash(fullfile(source_dir,name{1})),black_file_hash(fullfile(output_dir,name{1}))));
end
end
function text=once(text,a,b)
assert(isscalar(strfind(text,char(a))),'Patch anchor drift.'); text=strrep(char(text),char(a),char(b));
end
function write(path,text)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); c=onCleanup(@()fclose(fid)); fprintf(fid,'%s',text);
end