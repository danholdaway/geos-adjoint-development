SUBROUTINE ropp_fm_state2state_gsi_1d(n_lev,state,temp,shum,pres,geop,geop_sfc,use_logq)

  USE ropp_fm_constants

  IMPLICIT NONE

  INTEGER, PARAMETER :: wp = 4
!  REAL(wp), PARAMETER ::  

  INTEGER, INTENT(inout) :: n_lev
  REAL,    INTENT(in)    :: state(3*n_lev+1)
  LOGICAL, INTENT(in)    :: use_logq
  REAL,    INTENT(inout) :: temp(n_lev)
  REAL,    INTENT(inout) :: shum(n_lev)
  REAL,    INTENT(inout) :: pres(n_lev)
  REAL,    INTENT(inout) :: geop(n_lev)
  REAL,    INTENT(inout) :: geop_sfc


  REAL(wp), DIMENSION(:), ALLOCATABLE :: xstatetmp ! Temporary holder
  REAL(wp), DIMENSION(:), ALLOCATABLE :: pedg      ! Pressure at level edges
  REAL(wp), DIMENSION(:), ALLOCATABLE :: Tv        ! Virtual temperauture (K)
  REAL(wp), DIMENSION(:), ALLOCATABLE :: geope     ! Geopot height half level
  REAL(wp), DIMENSION(:), ALLOCATABLE :: lnpm  ! log ratio of pressures
  REAL(wp), DIMENSION(:), ALLOCATABLE :: delp      ! Pressure difference (Pa)
  REAL(wp), DIMENSION(:), ALLOCATABLE :: delz      ! Thickness (m)
  REAL(wp), DIMENSION(:), ALLOCATABLE :: alpha     ! Interpolation coefficient
  REAL(wp), DIMENSION(:), ALLOCATABLE :: logp      ! Log pressure at edges
  REAL(wp), DIMENSION(:), ALLOCATABLE :: pkap      ! p to the kappa at edges
  REAL(wp)                            :: XKappa, XkappaI, dlog, dpkap

  INTEGER :: nlev, elev, k, k1

  nlev = n_lev    !Number of levels
  elev = n_lev+1  !Number of level edges

  ALLOCATE(xstatetmp(SIZE(state)))
  ALLOCATE(Tv       (nlev))
  ALLOCATE(pedg     (elev))
  ALLOCATE(geope    (elev))
  ALLOCATE(lnpm     (nlev))
  ALLOCATE(delp     (nlev))
  ALLOCATE(delz     (nlev))
  ALLOCATE(alpha    (nlev))
  ALLOCATE(logp     (elev))
  ALLOCATE(pkap     (elev))

  !Virtual temperature
  Tv = state(1:nlev)

  !Specific humidity
  shum = state(nlev+1:2*nlev)

  !Convert virtual temperature to temperature
  temp = Tv / (1.0_wp + 0.61_wp * shum)

  !Pressure at the edge of levels from log(p) x10x100 is to hPa to Pa
  pedg = 10.0_wp * 100.0_wp * exp(state(2*nlev+1:3*nlev+1))

  !Pressue at the mid point using the GEOS-5 calculation
  Xkappa = R_dry/C_p
  XkappaI = 1./Xkappa
  DO k = 1,elev
     logp(k)=log(pedg(k))
     pkap(k)=pedg(k)**Xkappa
  ENDDO
  DO k = 1,nlev
     dlog  = logp(k)-logp(k+1)
     dpkap = pkap(k)-pkap(k+1)
     pres(k) = (XkappaI*dpkap/dlog)**XkappaI
  ENDDO  

  !Update specific humidity - circular on tv because of calc on check_qsat?
  IF (use_logq) THEN
     shum = EXP(shum) / 1000.0_wp
     CALL check_qsat(shum, temp, pres, nlev)
  ELSE
     WHERE (shum <= 0.0_wp)
      shum = 1.0e-9_wp
    ENDWHERE
  ENDIF

  !Geopotential height - follows ECMWF version
  delp = pedg(1:elev-1) - pedg(2:elev) !Pressure differences at mid

  lnpm = LOG(pedg(1:elev-1)/pedg(2:elev)) !Log of pressure ratio at mid

  alpha    = 1.0_wp - pedg(2:elev)/delp * lnpm !Interpolation coefficients
  alpha(elev-1) = LOG(2.0_wp)

  delz = R_dry * Tv * lnpm / g_wmo !Function to be integrated

  !Calculate geopotential height integral
  DO k = 1, elev
     geope(k) = geop_sfc + SUM(delz(1:k-1))
  ENDDO

  !Interpolate onto full levels 
  geop = geope(1:elev-1) + alpha * R_dry * Tv / g_wmo
  
  !Clean up
  DEALLOCATE(xstatetmp)
  DEALLOCATE(Tv)
  DEALLOCATE(pedg)
  DEALLOCATE(geope)
  DEALLOCATE(lnpm)
  DEALLOCATE(delp)
  DEALLOCATE(delz)
  DEALLOCATE(alpha)
  DEALLOCATE(logp)
  DEALLOCATE(pkap)

CONTAINS
  
!-------------------------------------------------------------------------------
! 9. Compute Saturation Vapour Pressure - limit log(humidity) within saturation 
!-------------------------------------------------------------------------------

  SUBROUTINE check_qsat(shum, temp, press, nlev)

    IMPLICIT NONE

    INTEGER, PARAMETER :: wp = 4

    INTEGER,                   INTENT(in)    :: nlev
    REAL(wp), DIMENSION(nlev), INTENT(inout) :: shum
    REAL(wp), DIMENSION(nlev), INTENT(in)    :: temp
    REAL(wp), DIMENSION(nlev), INTENT(in)    :: press
    
    REAL(wp), DIMENSION(:), ALLOCATABLE   :: Z
    REAL(wp), DIMENSION(:), ALLOCATABLE   :: X
    REAL(wp), DIMENSION(:), ALLOCATABLE   :: tfrac
    REAL(wp), DIMENSION(:), ALLOCATABLE   :: es
    REAL(wp), DIMENSION(:), ALLOCATABLE   :: qs

    ALLOCATE(es(SIZE(shum)))
    ALLOCATE(qs(SIZE(shum)))
    ALLOCATE(z(SIZE(shum)))
    ALLOCATE(x(SIZE(shum)))
    ALLOCATE(tfrac(SIZE(shum)))

! 9.1 Goff-Gratch equation (saturation vapour pressure over water)
    
    tfrac = 373.16_wp / temp
    Z = (-7.90298_wp * (tfrac - 1.0_wp)) +  &
       (5.02808_wp * LOG10(tfrac)) - (1.3816e-7_wp * (10**(11.344*(1.0_wp - (1.0_wp/tfrac))) - 1.0_wp)) + &
       (8.1328e-3_wp * (10**(-3.49149*(tfrac-1.0_wp)) - 1.0_wp))
    
! 9.2 Goff-Gratch equation (saturation vapour pressure over ice)

    tfrac(:) = 273.16_wp / temp(:)
    X = (-9.09718_wp * (tfrac - 1.0_wp)) + (-3.56654_wp*LOG10(tfrac)) + (0.876793_wp*(1.0_wp - (1.0_wp/tfrac)))

    WHERE(temp >= 273.16)
      es = 1013246.0_wp * (10**Z)
    ELSEWHERE
      es = 610.71_wp * (10**X)
    ENDWHERE
    
! 9.3 Compute saturation humidity

    qs = es * epsilon_water / (MAX(press,es) - ((1.0_wp - epsilon_water)*es))

    WHERE(qs == 1.0_wp)
      qs = 2.0_wp*shum(SIZE(shum))
    ENDWHERE

! 9.4 Limit humidity data below saturation

    WHERE(shum > qs)
      shum = qs
    ENDWHERE

    DEALLOCATE(es)
    DEALLOCATE(qs)
    DEALLOCATE(z)
    DEALLOCATE(x)
    DEALLOCATE(tfrac)

  END SUBROUTINE check_qsat

END SUBROUTINE ropp_fm_state2state_gsi_1d 
