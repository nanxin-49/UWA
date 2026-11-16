function receipt=prepare_bellhop_rotated_performance_overlay(source_dir,output_dir,enable_profile)
%PREPARE_BELLHOP_ROTATED_PERFORMANCE_OVERLAY Isolated unchanged-physics profiler.
arguments
    source_dir (1,:) char
    output_dir (1,:) char
    enable_profile (1,1) logical=true
end
assert(~isfolder(output_dir)); mkdir(output_dir);
files=dir(fullfile(source_dir,'*.f90'));
for k=1:numel(files), copyfile(fullfile(files(k).folder,files(k).name),fullfile(output_dir,files(k).name)); end
step=normtext(fileread(fullfile(output_dir,'Step.f90')));
declarations=join(["  INTEGER :: PerfWallUnit"; ...
"  INTEGER(KIND=8) :: PerfStepCalls=0,PerfReduceCalls=0,PerfCrossCalls=0,PerfGeometryCalls=0"; ...
"  INTEGER(KIND=8) :: PerfCrossChecks=0,PerfGeometryChecks=0,PerfCrossCandidates=0,PerfGeometryCandidates=0"; ...
"  INTEGER(KIND=8) :: PerfInactiveGeometry=0"],newline);
step=once(step,'  LOGICAL :: RotQEnabled = .FALSE.',declarations+newline+"  LOGICAL :: RotQEnabled = .FALSE.");
step=once(step,'    iSeg = 0',"    PerfGeometryCalls=PerfGeometryCalls+1"+newline+"    iSeg = 0");
step=once(step,['    hWall = HUGE( hWall )' newline '    IF ( WallNProfile <= 0 ) RETURN'],"    PerfCrossCalls=PerfCrossCalls+1"+newline+"    hWall = HUGE( hWall )"+newline+"    IF ( WallNProfile <= 0 ) RETURN");
a=strfind(step,'  SUBROUTINE GetWallGeometry'); b=strfind(step,'  END SUBROUTINE GetWallGeometry'); part=step(a:b-1);
part=strrep(part,'    DO i = 1, WallNProfile - 1',"    DO i = 1, WallNProfile - 1"+newline+"       PerfGeometryChecks=PerfGeometryChecks+1");
part=strrep(part,'          dist =',"          PerfGeometryCandidates=PerfGeometryCandidates+1"+newline+"          dist ="); step=[step(1:a-1) char(part) step(b:end)];
a=strfind(step,'  SUBROUTINE WallCrossingStep'); b=strfind(step,'  END SUBROUTINE WallCrossingStep'); part=step(a:b-1);
part=strrep(part,'    DO i = 1, WallNProfile - 1',"    DO i = 1, WallNProfile - 1"+newline+"       PerfCrossChecks=PerfCrossChecks+1");
part=strrep(part,'               hWall = MIN( hWall, hCand )',join(["          THEN";"             PerfCrossCandidates=PerfCrossCandidates+1";"             hWall = MIN( hWall, hCand )";"          END IF"],newline));
part=strrep(part,'lam <= 1.0d0 + tol ) &','lam <= 1.0d0 + tol ) &'); step=[step(1:a-1) char(part) step(b:end)];
step=once(step,"    CALL EvaluateSSP( ray0%x,", "    PerfStepCalls=PerfStepCalls+1"+newline+"    IF(.NOT.WallActive) PerfInactiveGeometry=PerfInactiveGeometry+1"+newline+"    CALL EvaluateSSP( ray0%x,");
step=once(step,'    x = x0 + h * urayt',"    PerfReduceCalls=PerfReduceCalls+1"+newline+"    x = x0 + h * urayt");
walllog=join([ ...
"  SUBROUTINE WritePerformanceWall(seg,lambda,t,n,kappa,residual,a,b,c)"; ...
"    INTEGER,INTENT(IN)::seg"; ...
"    REAL(KIND=8),INTENT(IN)::lambda,t(2),n(2),kappa,residual"; ...
"    TYPE(ray2DPt),INTENT(IN)::a,b,c"; ...
"    WRITE(PerfWallUnit) seg,a%NumTopBnc,a%NumBotBnc,b%NumTopBnc,b%NumBotBnc,c%NumTopBnc,c%NumBotBnc"; ...
"    WRITE(PerfWallUnit) lambda,t,n,kappa,residual"; ...
"    WRITE(PerfWallUnit) a%x,a%t,a%p,a%q,a%c,a%Amp,a%Phase,a%out_of_plane_integral,a%tau"; ...
"    WRITE(PerfWallUnit) b%x,b%t,b%p,b%q,b%c,b%Amp,b%Phase,b%out_of_plane_integral,b%tau"; ...
"    WRITE(PerfWallUnit) c%x,c%t,c%p,c%q,c%c,c%Amp,c%Phase,c%out_of_plane_integral,c%tau"; ...
"  END SUBROUTINE WritePerformanceWall"],newline);
step=once(step,'END MODULE Step',walllog+newline+"END MODULE Step");
write(fullfile(output_dir,'Step.f90'),step);
main=normtext(fileread(fullfile(output_dir,'bellhop.f90')));
main=once(main,'  USE MathConstants',join(["  USE MathConstants"; ...
"  USE Step, ONLY: PerfWallUnit,PerfStepCalls,PerfReduceCalls,PerfCrossCalls,PerfGeometryCalls, &"; ...
"       PerfCrossChecks,PerfGeometryChecks,PerfCrossCandidates,PerfGeometryCandidates,PerfInactiveGeometry"],newline));
main=once(main,'  CALL BellhopCore',"  OPEN(NEWUNIT=PerfWallUnit,FILE=TRIM(FileRoot)//'.perfwall',ACCESS='STREAM',FORM='UNFORMATTED',STATUS='NEW')"+newline+"  CALL BellhopCore");
closecode=join(["  CLOSE(PerfWallUnit)"; ...
"  OPEN(UNIT=93,FILE=TRIM(FileRoot)//'.perfcounts',STATUS='NEW')"; ...
"  WRITE(93,'(A)') '# step reduce cross geometry cross_checks geom_checks cross_candidates geom_candidates inactive_geometry'"; ...
"  WRITE(93,'(9(I0,1X))') PerfStepCalls,PerfReduceCalls,PerfCrossCalls,PerfGeometryCalls, &"; ...
"       PerfCrossChecks,PerfGeometryChecks,PerfCrossCandidates,PerfGeometryCandidates,PerfInactiveGeometry"; ...
"  CLOSE(93)"],newline);
main=once(main,'  CALL CloseExactReceiverAudit()',closecode+newline+"  CALL CloseExactReceiverAudit()");
main=once(main,'        CovRotated        = RotatedState',"        CovRotated        = RotatedState"+newline+"        CALL WritePerformanceWall(WallSeg,WallLambda,WallT,WallN,WallKappa,WallResidual,IncidentState,ReflectedState,RotatedState)");
write(fullfile(output_dir,'bellhop.f90'),main);
build=normtext(fileread(fullfile(source_dir,'build_rotated_receiver64.ps1')));
if enable_profile
    build=once(build,"$compileFlags = @('-O2'","$compileFlags = @('-pg','-O2'");
end
write(fullfile(output_dir,'build_performance.ps1'),build);
receipt=struct('source_dir',source_dir,'output_dir',output_dir,'profile_enabled',enable_profile,'status','PREPARED_NO_SOLVE');
write(fullfile(output_dir,'checkpoint.json'),jsonencode(receipt,PrettyPrint=true)); disp(jsonencode(receipt));
end
function text=normtext(text)
text=strrep(text,sprintf('\r\n'),newline);
end
function text=once(text,a,b)
assert(isscalar(strfind(text,char(a))),'Patch anchor must be unique.'); text=strrep(text,char(a),char(b));
end
function write(path,text)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); close=onCleanup(@()fclose(fid)); fprintf(fid,'%s',text);
end