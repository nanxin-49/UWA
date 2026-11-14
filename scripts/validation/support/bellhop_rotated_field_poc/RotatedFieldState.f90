! Validation-only state for the Bellhop rotated-field experiment.
! Original Bellhop sources/installations are not overwritten.
MODULE RotatedFieldState
  USE bellhopMod, ONLY: Beam
  USE SourceReceiverPositions, ONLY: Pos
  USE sspMod, ONLY: FastQSearch
  USE FatalError
  IMPLICIT NONE
  LOGICAL :: Point3DEnabled = .FALSE., ExactCoordinatesEnabled = .FALSE.
  REAL(KIND=8), ALLOCATABLE :: ReceiverRange64(:), ReceiverDepth64(:)
CONTAINS
  SUBROUTINE ConfigureRotatedField(FileRoot)
    CHARACTER(LEN=*), INTENT(IN) :: FileRoot
    INTEGER :: unit, ios, pointFlag, fastFlag, nr, nz
    LOGICAL :: exists
    IF(ALLOCATED(ReceiverRange64)) DEALLOCATE(ReceiverRange64,ReceiverDepth64)
    ALLOCATE(ReceiverRange64(Pos%NRr),ReceiverDepth64(Pos%NRz))
    ReceiverRange64=REAL(Pos%Rr(1:Pos%NRr),KIND=8)
    ReceiverDepth64=REAL(Pos%Rz(1:Pos%NRz),KIND=8)
    Point3DEnabled=.FALSE.
    ExactCoordinatesEnabled=.FALSE.
    FastQSearch=.FALSE.
    INQUIRE(FILE=TRIM(FileRoot)//'.r25d',EXIST=exists)
    IF(.NOT.exists) RETURN
    OPEN(NEWUNIT=unit,FILE=TRIM(FileRoot)//'.r25d',STATUS='OLD',ACTION='READ',IOSTAT=ios)
    IF(ios/=0) CALL ERROUT('BELLHOP-R25D','Cannot open options/coordinate sidecar')
    READ(unit,*,IOSTAT=ios) pointFlag,fastFlag
    IF(ios/=0) CALL ERROUT('BELLHOP-R25D','Invalid point/fast-search flags')
    IF((pointFlag/=0.AND.pointFlag/=1).OR.(fastFlag/=0.AND.fastFlag/=1)) &
       CALL ERROUT('BELLHOP-R25D','Flags must be zero or one')
    READ(unit,*,IOSTAT=ios) nr,nz
    IF(ios/=0.OR.nr/=Pos%NRr.OR.nz/=Pos%NRz) CALL ERROUT('BELLHOP-R25D','Coordinate-count mismatch')
    READ(unit,*,IOSTAT=ios) ReceiverRange64
    IF(ios/=0) CALL ERROUT('BELLHOP-R25D','Missing exact receiver ranges')
    READ(unit,*,IOSTAT=ios) ReceiverDepth64
    IF(ios/=0) CALL ERROUT('BELLHOP-R25D','Missing exact receiver depths')
    CLOSE(unit)
    IF(ANY(ReceiverRange64<=0d0)) CALL ERROUT('BELLHOP-R25D','Exact range must be positive; use a second chart at r=0')
    IF(ANY(ReceiverRange64(2:nr)<=ReceiverRange64(1:nr-1))) CALL ERROUT('BELLHOP-R25D','Ranges must increase')
    IF(ANY(ReceiverDepth64(2:nz)<=ReceiverDepth64(1:nz-1))) CALL ERROUT('BELLHOP-R25D','Depths must increase')
    ! SHD metadata remains native float32. Exact samples may differ only by its
    ! known coordinate quantization, never by a fitted offset or resampling.
    IF(ANY(ABS(ReceiverRange64-REAL(Pos%Rr(1:nr),KIND=8)) > &
       MAX(1d-9,2d0*REAL(SPACING(Pos%Rr(1:nr)),KIND=8)))) CALL ERROUT('BELLHOP-R25D','Range metadata mismatch')
    IF(ANY(ABS(ReceiverDepth64-REAL(Pos%Rz(1:nz),KIND=8)) > &
       MAX(1d-9,2d0*ABS(REAL(SPACING(Pos%Rz(1:nz)),KIND=8))))) CALL ERROUT('BELLHOP-R25D','Depth metadata mismatch')
    Point3DEnabled=pointFlag==1
    ExactCoordinatesEnabled=.TRUE.
    FastQSearch=fastFlag==1
    IF(Point3DEnabled.AND.Beam%RunType(4:4)/='R') CALL ERROUT('BELLHOP-R25D','3D correction requires point geometry R')
    IF(Point3DEnabled.AND.Beam%Type(1:1)/='G') CALL ERROUT('BELLHOP-R25D','3D pilot requires Cartesian geometric hat beams')
    WRITE(6,*) 'Validation: point3D/exact-coordinates/fast-Q = ',Point3DEnabled,ExactCoordinatesEnabled,FastQSearch
  END SUBROUTINE ConfigureRotatedField
END MODULE RotatedFieldState