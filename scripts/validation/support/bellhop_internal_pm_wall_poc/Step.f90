MODULE Step

  USE bellhopMod
  USE sspMod
  IMPLICIT NONE

  ! Validation-only sampled internal PM wall. The profile points are ordered
  ! exactly as Bellhop's ATI points (in increasing range); two constant-depth
  ! extension points reproduce ReadATI/ComputeBdryTangentNormal's infinite
  ! left/right extensions.
  INTEGER :: WallNProfile = 0
  REAL (KIND=8), ALLOCATABLE :: WallX( :, : ), WallSegT( :, : ), WallNodeT( :, : ), WallNodeN( :, : )
  REAL (KIND=8), ALLOCATABLE :: WallSegLen( : ), WallDx( : ), WallSegKappa( : )

CONTAINS

  SUBROUTINE ConfigureWallProfile( R, Z, N )
    INTEGER, INTENT( IN ) :: N
    REAL (KIND=8), INTENT( IN ) :: R( N ), Z( N )
    INTEGER :: i, nAll
    REAL (KIND=8), ALLOCATABLE :: phi( : )

    IF ( N < 3 ) STOP 'BELLHOP-IWALL: PM wall requires at least 3 profile points'
    WallNProfile = N
    nAll = N + 2
    IF ( ALLOCATED( WallX ) ) DEALLOCATE( WallX, WallSegT, WallNodeT, WallNodeN, WallSegLen, WallDx, WallSegKappa )
    ALLOCATE( WallX( 2, nAll ), WallSegT( 2, nAll - 1 ), WallNodeT( 2, nAll ), WallNodeN( 2, nAll ), &
              WallSegLen( nAll - 1 ), WallDx( nAll ), WallSegKappa( nAll - 1 ), phi( nAll ) )

    WallX( :, 1 ) = [ -SQRT( HUGE( 1.0d0 ) ) / 1.0d5, Z( 1 ) ]
    WallX( 1, 2 : N + 1 ) = R
    WallX( 2, 2 : N + 1 ) = Z
    WallX( :, nAll ) = [ SQRT( HUGE( 1.0d0 ) ) / 1.0d5, Z( N ) ]

    DO i = 1, nAll - 1
       WallSegT( :, i ) = WallX( :, i + 1 ) - WallX( :, i )
       WallSegLen( i ) = NORM2( WallSegT( :, i ) )
       WallSegT( :, i ) = WallSegT( :, i ) / WallSegLen( i )
       WallDx( i ) = ( WallX( 2, i + 1 ) - WallX( 2, i ) ) / ( WallX( 1, i + 1 ) - WallX( 1, i ) )
    END DO
    WallDx( nAll ) = 0.0d0

    WallNodeT( :, 1 ) = [ 1.0d0, 0.0d0 ]
    WallNodeT( :, nAll ) = [ 1.0d0, 0.0d0 ]
    DO i = 2, nAll - 1
       ! Literal averaging used by Bellhop's C ATI path.
       WallNodeT( :, i ) = 0.5d0 * ( WallSegT( :, i - 1 ) + WallSegT( :, i ) )
    END DO
    DO i = 1, nAll
       ! Outward normal for a TOP boundary, copied from bdryMod.f90.
       WallNodeN( 1, i ) = +WallNodeT( 2, i )
       WallNodeN( 2, i ) = -WallNodeT( 1, i )
       phi( i ) = ATAN2( WallNodeT( 2, i ), WallNodeT( 1, i ) )
    END DO

    DO i = 1, nAll - 1
       ! ComputeBdryTangentNormal's C-ATI curvature path, including its
       ! Dss override. Keeping this expression identical is important for
       ! a native-to-internal curvature audit.
       WallSegKappa( i ) = ( phi( i + 1 ) - phi( i ) ) / WallSegLen( i )
       WallSegKappa( i ) = ( WallDx( i + 1 ) - WallDx( i ) ) / &
            ( WallX( 1, i + 1 ) - WallX( 1, i ) ) * WallSegT( 1, i ) ** 3
    END DO
    DEALLOCATE( phi )
  END SUBROUTINE ConfigureWallProfile

  SUBROUTINE GetWallGeometry( x, iSeg, lambda, tWall, nWall, kappa, residual )
    REAL (KIND=8), INTENT( IN ) :: x( 2 )
    INTEGER, INTENT( OUT ) :: iSeg
    REAL (KIND=8), INTENT( OUT ) :: lambda, tWall( 2 ), nWall( 2 ), kappa, residual
    INTEGER :: i
    REAL (KIND=8) :: d( 2 ), q( 2 ), lam, dist, best

    iSeg = 0
    lambda = 0.0d0
    tWall = 0.0d0
    nWall = 0.0d0
    kappa = 0.0d0
    residual = HUGE( residual )
    best = HUGE( best )
    IF ( WallNProfile <= 0 ) RETURN
    DO i = 1, WallNProfile + 1
       d = WallX( :, i + 1 ) - WallX( :, i )
       q = x - WallX( :, i )
       lam = DOT_PRODUCT( q, d ) / DOT_PRODUCT( d, d )
       IF ( lam >= -1.0d-8 .AND. lam <= 1.0d0 + 1.0d-8 ) THEN
          dist = ABS( q( 1 ) * d( 2 ) - q( 2 ) * d( 1 ) ) / WallSegLen( i )
          IF ( dist < best ) THEN
             best = dist
             iSeg = i
             lambda = MIN( 1.0d0, MAX( 0.0d0, lam ) )
          END IF
       END IF
    END DO
    IF ( iSeg > 0 ) THEN
       residual = ( ( x( 1 ) - WallX( 1, iSeg ) ) * ( WallX( 2, iSeg + 1 ) - WallX( 2, iSeg ) ) - &
                    ( x( 2 ) - WallX( 2, iSeg ) ) * ( WallX( 1, iSeg + 1 ) - WallX( 1, iSeg ) ) ) / WallSegLen( iSeg )
       tWall = ( 1.0d0 - lambda ) * WallNodeT( :, iSeg ) + lambda * WallNodeT( :, iSeg + 1 )
       nWall = ( 1.0d0 - lambda ) * WallNodeN( :, iSeg ) + lambda * WallNodeN( :, iSeg + 1 )
       kappa = WallSegKappa( iSeg )
    END IF
  END SUBROUTINE GetWallGeometry

  SUBROUTINE WallCrossingStep( x0, u, hTrial, hWall )
    REAL (KIND=8), INTENT( IN ) :: x0( 2 ), u( 2 ), hTrial
    REAL (KIND=8), INTENT( OUT ) :: hWall
    INTEGER :: i
    REAL (KIND=8) :: d( 2 ), q( 2 ), den, hCand, lam, tol

    hWall = HUGE( hWall )
    IF ( WallNProfile <= 0 ) RETURN
    tol = 1.0d-10
    DO i = 1, WallNProfile + 1
       d = WallX( :, i + 1 ) - WallX( :, i )
       q = WallX( :, i ) - x0
       den = u( 1 ) * d( 2 ) - u( 2 ) * d( 1 )
       IF ( ABS( den ) > 1.0d-14 ) THEN
          hCand = ( q( 1 ) * d( 2 ) - q( 2 ) * d( 1 ) ) / den
          lam   = ( q( 1 ) * u( 2 ) - q( 2 ) * u( 1 ) ) / den
          IF ( hCand > tol .AND. hCand <= hTrial + tol .AND. lam >= -tol .AND. lam <= 1.0d0 + tol ) &
               hWall = MIN( hWall, hCand )
       END IF
    END DO
  END SUBROUTINE WallCrossingStep

  SUBROUTINE Step2D( ray0, ray2, Topx, Topn, Botx, Botn, WallActive, WallHit )

    TYPE( ray2DPt )    :: ray0, ray1, ray2
    REAL (KIND=8 ), INTENT( IN ) :: Topx( 2 ), Topn( 2 ), Botx( 2 ), Botn( 2 )
    LOGICAL,         INTENT( IN ) :: WallActive
    LOGICAL,        INTENT( OUT ) :: WallHit
    INTEGER            :: iSegz0, iSegr0, iWallSeg
    REAL     (KIND=8 ) :: gradc0( 2 ), gradc1( 2 ), gradc2( 2 ), &
         c0, cimag0, crr0, crz0, czz0, csq0, cnn0_csq0, &
         c1, cimag1, crr1, crz1, czz1, csq1, cnn1_csq1, &
         c2, cimag2, crr2, crz2, czz2, csq2, urayt0( 2 ), urayt1( 2 ), &
         h, halfh, hw0, hw1, ray2n( 2 ), RM, RN, gradcjump( 2 ), cnjump, csjump, w0, w1, rho, &
         wallLambda, wallKappa, wallResidual, wallT( 2 ), wallN( 2 )

    CALL EvaluateSSP( ray0%x, c0, cimag0, gradc0, crr0, crz0, czz0, rho, freq, 'TAB' )
    csq0      = c0 * c0
    cnn0_csq0 = crr0 * ray0%t( 2 )**2 - 2.0 * crz0 * ray0%t( 1 ) * ray0%t( 2 ) + czz0 * ray0%t( 1 )**2
    iSegz0    = iSegz
    iSegr0    = iSegr
    h = Beam%deltas
    urayt0 = c0 * ray0%t
    CALL ReduceStep2D( ray0%x, urayt0, iSegz0, iSegr0, Topx, Topn, Botx, Botn, WallActive, h )
    halfh = 0.5 * h

    ray1%x = ray0%x + halfh * urayt0
    ray1%t = ray0%t - halfh * gradc0 / csq0
    ray1%p = ray0%p - halfh * cnn0_csq0 * ray0%q
    ray1%q = ray0%q + halfh * c0        * ray0%p

    CALL EvaluateSSP( ray1%x, c1, cimag1, gradc1, crr1, crz1, czz1, rho, freq, 'TAB' )
    csq1      = c1 * c1
    cnn1_csq1 = crr1 * ray1%t( 2 )**2 - 2.0 * crz1 * ray1%t( 1 ) * ray1%t( 2 ) + czz1 * ray1%t( 1 )**2
    urayt1 = c1 * ray1%t
    CALL ReduceStep2D( ray0%x, urayt1, iSegz0, iSegr0, Topx, Topn, Botx, Botn, WallActive, h )

    w1  = h / ( 2.0d0 * halfh )
    w0  = 1.0d0 - w1
    hw0 = h * w0
    hw1 = h * w1
    ray2%x   = ray0%x   + hw0 * urayt0              + hw1 * urayt1
    ray2%t   = ray0%t   - hw0 * gradc0 / csq0       - hw1 * gradc1 / csq1
    ray2%p   = ray0%p   - hw0 * cnn0_csq0 * ray0%q  - hw1 * cnn1_csq1 * ray1%q
    ray2%q   = ray0%q   + hw0 * c0        * ray0%p  + hw1 * c1        * ray1%p
    ray2%tau = ray0%tau + hw0 / CMPLX( c0, cimag0, KIND=8 ) + hw1 / CMPLX( c1, cimag1, KIND=8 )
    ray2%Amp = ray0%Amp
    ray2%Phase = ray0%Phase
    ray2%NumTopBnc = ray0%NumTopBnc
    ray2%NumBotBnc = ray0%NumBotBnc

    CALL EvaluateSSP( ray2%x, c2, cimag2, gradc2, crr2, crz2, czz2, rho, freq, 'TAB' )
    ray2%c = c2
    CALL GetWallGeometry( ray2%x, iWallSeg, wallLambda, wallT, wallN, wallKappa, wallResidual )
    WallHit = WallActive .AND. iWallSeg > 0 .AND. ABS( wallResidual ) <= 1.0d-8

    IF ( iSegz /= iSegz0 .OR. iSegr /= iSegr0 ) THEN
       gradcjump = gradc2 - gradc0
       ray2n     = [ -ray2%t( 2 ), ray2%t( 1 ) ]
       cnjump    = DOT_PRODUCT( gradcjump, ray2n )
       csjump    = DOT_PRODUCT( gradcjump, ray2%t )
       IF ( iSegz /= iSegz0 ) THEN
          RM = +ray2%t( 1 ) / ray2%t( 2 )
       ELSE
          RM = -ray2%t( 2 ) / ray2%t( 1 )
       END IF
       RN     = RM * ( 2 * cnjump - RM * csjump ) / c2
       ray2%p = ray2%p - ray2%q * RN
    END IF
  END SUBROUTINE Step2D

  SUBROUTINE ReduceStep2D( x0, urayt, iSegz0, iSegr0, Topx, Topn, Botx, Botn, WallActive, h )
    USE BdryMod
    INTEGER,       INTENT( IN    ) :: iSegz0, iSegr0
    REAL (KIND=8), INTENT( IN    ) :: x0( 2 ), urayt( 2 )
    REAL (KIND=8), INTENT( IN    ) :: Topx( 2 ), Topn( 2 ), Botx( 2 ), Botn( 2 )
    LOGICAL,         INTENT( IN    ) :: WallActive
    REAL (KIND=8), INTENT( INOUT ) :: h
    REAL (KIND=8)                  :: x( 2 ), d( 2 ), d0( 2 ), h1, h2, h3, h4, hWall, rSeg( 2 )

    x = x0 + h * urayt
    h1 = HUGE( h1 )
    IF ( ABS( urayt( 2 ) ) > EPSILON( h1 ) ) THEN
       IF      ( SSP%z( iSegz0     ) > x(  2 ) ) THEN
          h1 = ( SSP%z( iSegz0 ) - x0( 2 ) ) / urayt( 2 )
       ELSE IF ( SSP%z( iSegz0 + 1 ) < x(  2 ) ) THEN
          h1 = ( SSP%z( iSegz0 + 1 ) - x0( 2 ) ) / urayt( 2 )
       END IF
    END IF
    h2 = HUGE( h2 )
    d  = x - Topx
    IF ( DOT_PRODUCT( Topn, d ) > EPSILON( h2 ) ) THEN
       d0 = x0 - Topx
       h2 = -DOT_PRODUCT( d0, Topn ) / DOT_PRODUCT( urayt, Topn )
    END IF
    h3 = HUGE( h3 )
    d  = x - Botx
    IF ( DOT_PRODUCT( Botn, d ) > EPSILON( h2 ) ) THEN
       d0 = x0 - Botx
       h3 = -DOT_PRODUCT( d0, Botn ) / DOT_PRODUCT( urayt, Botn )
    END IF
    rSeg( 1 ) = MAX( rTopSeg( 1 ), rBotSeg( 1 ) )
    rSeg( 2 ) = MIN( rTopSeg( 2 ), rBotSeg( 2 ) )
    IF ( SSP%Type == 'Q' ) THEN
       rSeg( 1 ) = MAX( rSeg( 1 ), SSP%Seg%r( iSegr0 ) )
       rSeg( 2 ) = MIN( rSeg( 2 ), SSP%Seg%r( iSegr0 + 1 ) )
    END IF
    h4 = HUGE( h4 )
    IF ( ABS( urayt( 1 ) ) > EPSILON( h4 ) ) THEN
       IF       ( x( 1 ) < rSeg( 1 ) ) THEN
          h4 = -( x0( 1 ) - rSeg( 1 ) ) / urayt( 1 )
       ELSE IF  ( x( 1 ) > rSeg( 2 ) ) THEN
          h4 = -( x0( 1 ) - rSeg( 2 ) ) / urayt( 1 )
       END IF
    END IF
    hWall = HUGE( hWall )
    IF ( WallActive ) CALL WallCrossingStep( x0, urayt, h, hWall )
    h = MIN( h, h1, h2, h3, h4, hWall )
    IF ( hWall < HUGE( hWall ) / 2.0d0 .AND. h == hWall ) THEN
       iSmallStepCtr = 0
    ELSE IF ( h < 1.0d-4 * Beam%deltas ) THEN
       h = 1.0d-5 * Beam%deltas
       iSmallStepCtr = iSmallStepCtr + 1
    ELSE
       iSmallStepCtr = 0
    END IF
  END SUBROUTINE ReduceStep2D

END MODULE Step
