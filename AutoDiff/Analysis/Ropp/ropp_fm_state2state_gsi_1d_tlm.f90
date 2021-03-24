!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.11 (r5903) - 14 Dec 2015 10:32
!
!  Differentiation of ropp_fm_state2state_gsi_1d in forward (tangent) mode:
!   variations   of useful results: temp shum geop pres
!   with respect to varying inputs: state
!   RW status of diff variables: temp:out shum:out state:in geop:out
!                pres:out
SUBROUTINE ROPP_FM_STATE2STATE_GSI_1D_TLM(n_lev, state, state_tl, temp, &
& temp_tl, shum, shum_tl, pres, pres_tl, geop, geop_tl, geop_sfc, &
& use_logq)
  USE ROPP_FM_CONSTANTS
  IMPLICIT NONE
  INTEGER, PARAMETER :: wp=4
!  REAL(wp), PARAMETER ::  
  INTEGER, INTENT(INOUT) :: n_lev
  REAL, INTENT(IN) :: state(3*n_lev+1)
  REAL, INTENT(IN) :: state_tl(3*n_lev+1)
  LOGICAL, INTENT(IN) :: use_logq
  REAL, INTENT(INOUT) :: temp(n_lev)
  REAL, INTENT(INOUT) :: temp_tl(n_lev)
  REAL, INTENT(INOUT) :: shum(n_lev)
  REAL, INTENT(INOUT) :: shum_tl(n_lev)
  REAL, INTENT(INOUT) :: pres(n_lev)
  REAL, INTENT(INOUT) :: pres_tl(n_lev)
  REAL, INTENT(INOUT) :: geop(n_lev)
  REAL, INTENT(INOUT) :: geop_tl(n_lev)
  REAL, INTENT(INOUT) :: geop_sfc
! Temporary holder
  REAL(wp), DIMENSION(:), ALLOCATABLE :: xstatetmp
! Pressure at level edges
  REAL(wp), DIMENSION(:), ALLOCATABLE :: pedg
  REAL(wp), DIMENSION(:), ALLOCATABLE :: pedg_tl
! Virtual temperauture (K)
  REAL(wp), DIMENSION(:), ALLOCATABLE :: tv
  REAL(wp), DIMENSION(:), ALLOCATABLE :: tv_tl
! Geopot height half level
  REAL(wp), DIMENSION(:), ALLOCATABLE :: geope
  REAL(wp), DIMENSION(:), ALLOCATABLE :: geope_tl
! log ratio of pressures
  REAL(wp), DIMENSION(:), ALLOCATABLE :: lnpm
  REAL(wp), DIMENSION(:), ALLOCATABLE :: lnpm_tl
! Pressure difference (Pa)
  REAL(wp), DIMENSION(:), ALLOCATABLE :: delp
  REAL(wp), DIMENSION(:), ALLOCATABLE :: delp_tl
! Thickness (m)
  REAL(wp), DIMENSION(:), ALLOCATABLE :: delz
  REAL(wp), DIMENSION(:), ALLOCATABLE :: delz_tl
! Interpolation coefficient
  REAL(wp), DIMENSION(:), ALLOCATABLE :: alpha
  REAL(wp), DIMENSION(:), ALLOCATABLE :: alpha_tl
! Log pressure at edges
  REAL(wp), DIMENSION(:), ALLOCATABLE :: logp
  REAL(wp), DIMENSION(:), ALLOCATABLE :: logp_tl
! p to the kappa at edges
  REAL(wp), DIMENSION(:), ALLOCATABLE :: pkap
  REAL(wp), DIMENSION(:), ALLOCATABLE :: pkap_tl
  REAL(wp) :: xkappa, xkappai, dlog, dpkap
  REAL(wp) :: dlog_tl, dpkap_tl
  INTEGER :: nlev, elev, k, k1
  INTRINSIC SIZE
  INTRINSIC EXP
  INTRINSIC LOG
  INTRINSIC SUM
  REAL(wp) :: pwx1
  REAL(wp) :: pwx1_tl
!Number of levels
  nlev = n_lev
!Number of level edges
  elev = n_lev + 1
  ALLOCATE(xstatetmp(SIZE(state)))
  ALLOCATE(tv_tl(nlev))
  ALLOCATE(tv(nlev))
  ALLOCATE(pedg_tl(elev))
  ALLOCATE(pedg(elev))
  ALLOCATE(geope_tl(elev))
  ALLOCATE(geope(elev))
  ALLOCATE(lnpm_tl(nlev))
  ALLOCATE(lnpm(nlev))
  ALLOCATE(delp_tl(nlev))
  ALLOCATE(delp(nlev))
  ALLOCATE(delz_tl(nlev))
  ALLOCATE(delz(nlev))
  ALLOCATE(alpha_tl(nlev))
  ALLOCATE(alpha(nlev))
  ALLOCATE(logp_tl(elev))
  ALLOCATE(logp(elev))
  ALLOCATE(pkap_tl(elev))
  ALLOCATE(pkap(elev))
!Virtual temperature
  tv_tl = state_tl(1:nlev)
  tv = state(1:nlev)
!Specific humidity
  shum_tl = state_tl(nlev+1:2*nlev)
  shum = state(nlev+1:2*nlev)
!Convert virtual temperature to temperature
  temp_tl = (tv_tl*(1.0_wp+0.61_wp*shum)-tv*0.61_wp*shum_tl)/(1.0_wp+&
&   0.61_wp*shum)**2
  temp = tv/(1.0_wp+0.61_wp*shum)
!Pressure at the edge of levels from log(p) x10x100 is to hPa to Pa
  pedg_tl = 10.0_wp*100.0_wp*state_tl(2*nlev+1:3*nlev+1)*EXP(state(2*&
&   nlev+1:3*nlev+1))
  pedg = 10.0_wp*100.0_wp*EXP(state(2*nlev+1:3*nlev+1))
!Pressue at the mid point using the GEOS-5 calculation
  xkappa = r_dry/c_p
  xkappai = 1./xkappa
  DO k=1,elev
    logp_tl(k) = pedg_tl(k)/pedg(k)
    logp(k) = LOG(pedg(k))
    IF (pedg(k) .GT. 0.0 .OR. (pedg(k) .LT. 0.0 .AND. xkappa .EQ. INT(&
&       xkappa))) THEN
      pkap_tl(k) = xkappa*pedg(k)**(xkappa-1)*pedg_tl(k)
    ELSE IF (pedg(k) .EQ. 0.0 .AND. xkappa .EQ. 1.0) THEN
      pkap_tl(k) = pedg_tl(k)
    ELSE
      pkap_tl(k) = 0.0
    END IF
    pkap(k) = pedg(k)**xkappa
  END DO
  pres_tl = 0.0
  DO k=1,nlev
    dlog_tl = logp_tl(k) - logp_tl(k+1)
    dlog = logp(k) - logp(k+1)
    dpkap_tl = pkap_tl(k) - pkap_tl(k+1)
    dpkap = pkap(k) - pkap(k+1)
    pwx1_tl = (xkappai*dpkap_tl*dlog-xkappai*dpkap*dlog_tl)/dlog**2
    pwx1 = xkappai*dpkap/dlog
    IF (pwx1 .GT. 0.0 .OR. (pwx1 .LT. 0.0 .AND. xkappai .EQ. INT(xkappai&
&       ))) THEN
      pres_tl(k) = xkappai*pwx1**(xkappai-1)*pwx1_tl
    ELSE IF (pwx1 .EQ. 0.0 .AND. xkappai .EQ. 1.0) THEN
      pres_tl(k) = pwx1_tl
    ELSE
      pres_tl(k) = 0.0
    END IF
    pres(k) = pwx1**xkappai
  END DO
!Update specific humidity - circular on tv because of calc on check_qsat?
  IF (use_logq) THEN
    shum_tl = shum_tl*EXP(shum)/1000.0_wp
    shum = EXP(shum)/1000.0_wp
    CALL CHECK_QSAT_TLM(shum, shum_tl, temp, temp_tl, pres, pres_tl, &
&                 nlev)
  ELSE
    WHERE (shum .LE. 0.0_wp) 
      shum_tl = 0.0
      shum = 1.0e-9_wp
    END WHERE
  END IF
!Geopotential height - follows ECMWF version
!Pressure differences at mid
  delp_tl = pedg_tl(1:elev-1) - pedg_tl(2:elev)
  delp = pedg(1:elev-1) - pedg(2:elev)
!Log of pressure ratio at mid
  lnpm_tl = (pedg_tl(1:elev-1)*pedg(2:elev)-pedg(1:elev-1)*pedg_tl(2:&
&   elev))/(pedg(2:elev)*pedg(1:elev-1))
  lnpm = LOG(pedg(1:elev-1)/pedg(2:elev))
!Interpolation coefficients
  alpha_tl = -((pedg_tl(2:elev)*delp-pedg(2:elev)*delp_tl)*lnpm/delp**2)&
&   - pedg(2:elev)*lnpm_tl/delp
  alpha = 1.0_wp - pedg(2:elev)/delp*lnpm
  alpha_tl(elev-1) = 0.0_4
  alpha(elev-1) = LOG(2.0_wp)
!Function to be integrated
  delz_tl = r_dry*(tv_tl*lnpm+tv*lnpm_tl)/g_wmo
  delz = r_dry*tv*lnpm/g_wmo
!Calculate geopotential height integral
  DO k=1,elev
    geope_tl(k) = SUM(delz_tl(1:k-1))
    geope(k) = geop_sfc + SUM(delz(1:k-1))
  END DO
!Interpolate onto full levels 
  geop_tl = geope_tl(1:elev-1) + r_dry*(alpha_tl*tv+alpha*tv_tl)/g_wmo
  geop = geope(1:elev-1) + alpha*r_dry*tv/g_wmo
!Clean up
  DEALLOCATE(xstatetmp)
  DEALLOCATE(tv_tl)
  DEALLOCATE(tv)
  DEALLOCATE(pedg_tl)
  DEALLOCATE(pedg)
  DEALLOCATE(geope_tl)
  DEALLOCATE(geope)
  DEALLOCATE(lnpm_tl)
  DEALLOCATE(lnpm)
  DEALLOCATE(delp_tl)
  DEALLOCATE(delp)
  DEALLOCATE(delz_tl)
  DEALLOCATE(delz)
  DEALLOCATE(alpha_tl)
  DEALLOCATE(alpha)
  DEALLOCATE(logp_tl)
  DEALLOCATE(logp)
  DEALLOCATE(pkap_tl)
  DEALLOCATE(pkap)

CONTAINS
!  Differentiation of check_qsat in forward (tangent) mode:
!   variations   of useful results: shum
!   with respect to varying inputs: temp shum press
!-------------------------------------------------------------------------------
! 9. Compute Saturation Vapour Pressure - limit log(humidity) within saturation 
!-------------------------------------------------------------------------------
  SUBROUTINE CHECK_QSAT_TLM(shum, shum_tl, temp, temp_tl, press, &
&   press_tl, nlev)
    USE DIFFSIZES
!  Hint: ISIZE1OFpwy1 should be the size of dimension 1 of array pwy1
!  Hint: ISIZE1OFpwy2 should be the size of dimension 1 of array pwy2
    IMPLICIT NONE
    INTEGER, PARAMETER :: wp=4
    INTEGER, INTENT(IN) :: nlev
    REAL(wp), DIMENSION(nlev), INTENT(INOUT) :: shum
    REAL(wp), DIMENSION(nlev), INTENT(INOUT) :: shum_tl
    REAL(wp), DIMENSION(nlev), INTENT(IN) :: temp
    REAL(wp), DIMENSION(nlev), INTENT(IN) :: temp_tl
    REAL(wp), DIMENSION(nlev), INTENT(IN) :: press
    REAL(wp), DIMENSION(nlev), INTENT(IN) :: press_tl
    REAL(wp), DIMENSION(:), ALLOCATABLE :: z
    REAL(wp), DIMENSION(:), ALLOCATABLE :: z_tl
    REAL(wp), DIMENSION(:), ALLOCATABLE :: x
    REAL(wp), DIMENSION(:), ALLOCATABLE :: x_tl
    REAL(wp), DIMENSION(:), ALLOCATABLE :: tfrac
    REAL(wp), DIMENSION(:), ALLOCATABLE :: tfrac_tl
    REAL(wp), DIMENSION(:), ALLOCATABLE :: es
    REAL(wp), DIMENSION(:), ALLOCATABLE :: es_tl
    REAL(wp), DIMENSION(:), ALLOCATABLE :: qs
    REAL(wp), DIMENSION(:), ALLOCATABLE :: qs_tl
    INTRINSIC SIZE
    INTRINSIC LOG10
    INTRINSIC MAX
    REAL(wp), DIMENSION(ISIZE1OFpwy1) :: pwy1
    REAL(wp), DIMENSION(ISIZE1OFpwy1) :: pwy1_tl
    REAL(wp), DIMENSION(ISIZE1OFpwy1) :: pwr1
    REAL(wp), DIMENSION(ISIZE1OFpwy1) :: pwr1_tl
    REAL(wp), DIMENSION(ISIZE1OFpwy2) :: pwy2
    REAL(wp), DIMENSION(ISIZE1OFpwy2) :: pwy2_tl
    REAL(wp), DIMENSION(ISIZE1OFpwy2) :: pwr2
    REAL(wp), DIMENSION(ISIZE1OFpwy2) :: pwr2_tl
    REAL(wp), DIMENSION(ISIZE1OFpwy1) :: pwy10
    REAL(wp), DIMENSION(ISIZE1OFpwy1) :: pwy10_tl
    REAL(wp), DIMENSION(ISIZE1OFpwy1) :: pwy11
    REAL(wp), DIMENSION(ISIZE1OFpwy1) :: pwy11_tl
    REAL*4 :: max1_tl(nlev)
    REAL*4 :: max1(nlev)
    ALLOCATE(es_tl(SIZE(shum)))
    ALLOCATE(es(SIZE(shum)))
    ALLOCATE(qs_tl(SIZE(shum)))
    ALLOCATE(qs(SIZE(shum)))
    ALLOCATE(z_tl(SIZE(shum)))
    ALLOCATE(z(SIZE(shum)))
    ALLOCATE(x_tl(SIZE(shum)))
    ALLOCATE(x(SIZE(shum)))
    ALLOCATE(tfrac_tl(SIZE(shum)))
    ALLOCATE(tfrac(SIZE(shum)))
! 9.1 Goff-Gratch equation (saturation vapour pressure over water)
    tfrac_tl = -(373.16_wp*temp_tl/temp**2)
    tfrac = 373.16_wp/temp
    pwy1_tl = 11.344*tfrac_tl/tfrac**2
    pwy1 = 11.344*(1.0_wp-1.0_wp/tfrac)
    pwr1_tl = 10**pwy1*LOG(10)*pwy1_tl
    pwr1 = 10**pwy1
    pwy2_tl = -(3.49149*tfrac_tl)
    pwy2 = -(3.49149*(tfrac-1.0_wp))
    pwr2_tl = 10**pwy2*LOG(10)*pwy2_tl
    pwr2 = 10**pwy2
    z_tl = 5.02808_wp*tfrac_tl/(tfrac*LOG(10.0)) - 7.90298_wp*tfrac_tl -&
&     1.3816e-7_wp*pwr1_tl + 8.1328e-3_wp*pwr2_tl
    z = -(7.90298_wp*(tfrac-1.0_wp)) + 5.02808_wp*LOG10(tfrac) - &
&     1.3816e-7_wp*(pwr1-1.0_wp) + 8.1328e-3_wp*(pwr2-1.0_wp)
! 9.2 Goff-Gratch equation (saturation vapour pressure over ice)
    tfrac_tl(:) = -(273.16_wp*temp_tl(:)/temp(:)**2)
    tfrac(:) = 273.16_wp/temp(:)
    x_tl = 0.876793_wp*tfrac_tl/tfrac**2 - 3.56654_wp*tfrac_tl/(tfrac*&
&     LOG(10.0)) - 9.09718_wp*tfrac_tl
    x = -(9.09718_wp*(tfrac-1.0_wp)) + (-(3.56654_wp*LOG10(tfrac))) + &
&     0.876793_wp*(1.0_wp-1.0_wp/tfrac)
    pwy10_tl = 0.0_4
    WHERE (temp .GE. 273.16) 
      pwy10_tl = z_tl
      pwy10 = z
      pwr1_tl = 10**pwy10*LOG(10)*pwy10_tl
      pwr1 = 10**pwy10
      es_tl = 1013246.0_wp*pwr1_tl
      es = 1013246.0_wp*pwr1
    END WHERE
    pwy11_tl = 0.0_4
    WHERE (.NOT.temp .GE. 273.16) 
      pwy11_tl = x_tl
      pwy11 = x
      pwr1_tl = 10**pwy11*LOG(10)*pwy11_tl
      pwr1 = 10**pwy11
      es_tl = 610.71_wp*pwr1_tl
      es = 610.71_wp*pwr1
    END WHERE
    max1_tl = 0.0_4
    WHERE (press .LT. es) 
      max1_tl = es_tl
      max1 = es
    ELSEWHERE
      max1_tl = press_tl
      max1 = press
    END WHERE
! 9.3 Compute saturation humidity
    qs_tl = (epsilon_water*es_tl*(max1-(1.0_wp-epsilon_water)*es)-es*&
&     epsilon_water*(max1_tl-(1.0_wp-epsilon_water)*es_tl))/(max1-(&
&     1.0_wp-epsilon_water)*es)**2
    qs = es*epsilon_water/(max1-(1.0_wp-epsilon_water)*es)
    WHERE (qs .EQ. 1.0_wp) 
      qs_tl = 2.0_wp*shum_tl(SIZE(shum))
      qs = 2.0_wp*shum(SIZE(shum))
    END WHERE
    WHERE (shum .GT. qs) 
      shum_tl = qs_tl
      shum = qs
    END WHERE
    DEALLOCATE(es_tl)
    DEALLOCATE(es)
    DEALLOCATE(qs_tl)
    DEALLOCATE(qs)
    DEALLOCATE(z_tl)
    DEALLOCATE(z)
    DEALLOCATE(x_tl)
    DEALLOCATE(x)
    DEALLOCATE(tfrac_tl)
    DEALLOCATE(tfrac)
  END SUBROUTINE CHECK_QSAT_TLM
!-------------------------------------------------------------------------------
! 9. Compute Saturation Vapour Pressure - limit log(humidity) within saturation 
!-------------------------------------------------------------------------------
  SUBROUTINE CHECK_QSAT(shum, temp, press, nlev)
    USE DIFFSIZES
!  Hint: ISIZE1OFpwy1 should be the size of dimension 1 of array pwy1
!  Hint: ISIZE1OFpwy2 should be the size of dimension 1 of array pwy2
    IMPLICIT NONE
    INTEGER, PARAMETER :: wp=4
    INTEGER, INTENT(IN) :: nlev
    REAL(wp), DIMENSION(nlev), INTENT(INOUT) :: shum
    REAL(wp), DIMENSION(nlev), INTENT(IN) :: temp
    REAL(wp), DIMENSION(nlev), INTENT(IN) :: press
    REAL(wp), DIMENSION(:), ALLOCATABLE :: z
    REAL(wp), DIMENSION(:), ALLOCATABLE :: x
    REAL(wp), DIMENSION(:), ALLOCATABLE :: tfrac
    REAL(wp), DIMENSION(:), ALLOCATABLE :: es
    REAL(wp), DIMENSION(:), ALLOCATABLE :: qs
    INTRINSIC SIZE
    INTRINSIC LOG10
    INTRINSIC MAX
    REAL(wp), DIMENSION(ISIZE1OFpwy1) :: pwy1
    REAL(wp), DIMENSION(ISIZE1OFpwy1) :: pwr1
    REAL(wp), DIMENSION(ISIZE1OFpwy2) :: pwy2
    REAL(wp), DIMENSION(ISIZE1OFpwy2) :: pwr2
    REAL(wp), DIMENSION(ISIZE1OFpwy1) :: pwy10
    REAL(wp), DIMENSION(ISIZE1OFpwy1) :: pwy11
    REAL*4 :: max1(nlev)
    ALLOCATE(es(SIZE(shum)))
    ALLOCATE(qs(SIZE(shum)))
    ALLOCATE(z(SIZE(shum)))
    ALLOCATE(x(SIZE(shum)))
    ALLOCATE(tfrac(SIZE(shum)))
! 9.1 Goff-Gratch equation (saturation vapour pressure over water)
    tfrac = 373.16_wp/temp
    pwy1 = 11.344*(1.0_wp-1.0_wp/tfrac)
    pwr1 = 10**pwy1
    pwy2 = -(3.49149*(tfrac-1.0_wp))
    pwr2 = 10**pwy2
    z = -(7.90298_wp*(tfrac-1.0_wp)) + 5.02808_wp*LOG10(tfrac) - &
&     1.3816e-7_wp*(pwr1-1.0_wp) + 8.1328e-3_wp*(pwr2-1.0_wp)
! 9.2 Goff-Gratch equation (saturation vapour pressure over ice)
    tfrac(:) = 273.16_wp/temp(:)
    x = -(9.09718_wp*(tfrac-1.0_wp)) + (-(3.56654_wp*LOG10(tfrac))) + &
&     0.876793_wp*(1.0_wp-1.0_wp/tfrac)
    WHERE (temp .GE. 273.16) 
      pwy10 = z
      pwr1 = 10**pwy10
      es = 1013246.0_wp*pwr1
    ELSEWHERE
      pwy11 = x
      pwr1 = 10**pwy11
      es = 610.71_wp*pwr1
    END WHERE
    WHERE (press .LT. es) 
      max1 = es
    ELSEWHERE
      max1 = press
    END WHERE
! 9.3 Compute saturation humidity
    qs = es*epsilon_water/(max1-(1.0_wp-epsilon_water)*es)
    WHERE (qs .EQ. 1.0_wp) qs = 2.0_wp*shum(SIZE(shum))
    WHERE (shum .GT. qs) shum = qs
    DEALLOCATE(es)
    DEALLOCATE(qs)
    DEALLOCATE(z)
    DEALLOCATE(x)
    DEALLOCATE(tfrac)
  END SUBROUTINE CHECK_QSAT
END SUBROUTINE ROPP_FM_STATE2STATE_GSI_1D_TLM
