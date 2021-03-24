module tp_core_tlm_mod
!BOP
!
! !MODULE: tp_core --- A collection of routines to support FV transport
!
 use fv_arrays_mod,     only: g_precision
 use fv_mp_mod,         only: is,js,ie,je, ng, isd,jsd,ied,jed
 use fv_grid_utils_mod, only: sw_corner, se_corner, ne_corner, nw_corner
 use fv_grid_utils_mod, only: sina_u, sina_v, da_min
 use fv_grid_tools_mod, only: dx, dy, rdxc, rdyc, rarea
 use fv_grid_tools_mod, only: dxa, dya, grid_type

 implicit none

      INTERFACE COPY_CORNERS_TLM
        MODULE PROCEDURE COPY_CORNERS_R8_TLM
      END INTERFACE

      INTERFACE FV_TP_2D_TLM
        MODULE PROCEDURE FV_TP_2D_R8_TLM
      END INTERFACE

 private
 public fv_tp_2d_tlm, copy_corners_tlm

 real, parameter:: r3 = 1./3.

#ifdef WAVE_FORM
! Suresh & Huynh scheme 2.2 (purtabation form)
! The wave-form is more diffusive than scheme 2.1
 real, parameter:: b1 =   0.0375
 real, parameter:: b2 =  -7./30.
 real, parameter:: b3 =  -23./120.
 real, parameter:: b4 =  13./30.
 real, parameter:: b5 = -11./240.
#else
! scheme 2.1: perturbation form
 real, parameter:: b1 =   1./30.
 real, parameter:: b2 = -13./60.
 real, parameter:: b3 = -13./60.
 real, parameter:: b4 =  0.45
 real, parameter:: b5 = -0.05
#endif
 real, parameter:: t11 = 27./28., t12 = -13./28., t13=3./7.
 real, parameter:: s11 = 11./14., s14 = 4./7.,    s15=3./14.
!----------------------------------------------------
! volume-conserving cubic with 2nd drv=0 at end point:
!----------------------------------------------------
! Non-monotonic
  real, parameter:: c1 = -2./14.
  real, parameter:: c2 = 11./14.
  real, parameter:: c3 =  5./14.
!----------------------
! PPM volume mean form:
!----------------------
  real, parameter:: p1 =  7./12.     ! 0.58333333
  real, parameter:: p2 = -1./12.

CONTAINS
!  Differentiation of fv_tp_2d_r8 in forward (tangent) mode (with options r8):
!   variations   of useful results: q fx fy
!   with respect to varying inputs: xfx q mfx mfy ra_x ra_y yfx
!                fx fy crx cry
  SUBROUTINE FV_TP_2D_R8_TLM(q, q_tl, crx, crx_tl, cry, cry_tl, npx, npy&
&   , hord, fx, fx_tl, fy, fy_tl, xfx, xfx_tl, yfx, yfx_tl, area, ra_x, &
&   ra_x_tl, ra_y, ra_y_tl, mfx, mfx_tl, mfy, mfy_tl, ppm_fac, nord, &
&   damp_c)
    IMPLICIT NONE
!#endif
    INTEGER, INTENT(IN) :: npx, npy
    INTEGER, INTENT(IN) :: hord
!
    REAL*8, INTENT(IN) :: crx(is:ie+1, jsd:jed)
    REAL*8, INTENT(IN) :: crx_tl(is:ie+1, jsd:jed)
!
    REAL*8, INTENT(IN) :: xfx(is:ie+1, jsd:jed)
    REAL*8, INTENT(IN) :: xfx_tl(is:ie+1, jsd:jed)
!
    REAL*8, INTENT(IN) :: cry(isd:ied, js:je+1)
    REAL*8, INTENT(IN) :: cry_tl(isd:ied, js:je+1)
!
    REAL*8, INTENT(IN) :: yfx(isd:ied, js:je+1)
    REAL*8, INTENT(IN) :: yfx_tl(isd:ied, js:je+1)
    REAL(g_precision), INTENT(IN) :: area(isd:ied, jsd:jed)
    REAL*8, INTENT(IN) :: ra_x(is:ie, jsd:jed)
    REAL*8, INTENT(IN) :: ra_x_tl(is:ie, jsd:jed)
    REAL*8, INTENT(IN) :: ra_y(isd:ied, js:je)
    REAL*8, INTENT(IN) :: ra_y_tl(isd:ied, js:je)
! transported scalar
    REAL*8, INTENT(INOUT) :: q(isd:ied, jsd:jed)
    REAL*8, INTENT(INOUT) :: q_tl(isd:ied, jsd:jed)
! Flux in x ( E )
    REAL*8, INTENT(OUT) :: fx(is:ie+1, js:je)
    REAL*8, INTENT(OUT) :: fx_tl(is:ie+1, js:je)
! Flux in y ( N )
    REAL*8, INTENT(OUT) :: fy(is:ie, js:je+1)
    REAL*8, INTENT(OUT) :: fy_tl(is:ie, js:je+1)
! optional Arguments:
! Mass Flux X-Dir
    REAL*8, OPTIONAL, INTENT(IN) :: mfx(is:ie+1, js:je)
    REAL*8, OPTIONAL, INTENT(IN) :: mfx_tl(is:ie+1, js:je)
! Mass Flux Y-Dir
    REAL*8, OPTIONAL, INTENT(IN) :: mfy(is:ie, js:je+1)
    REAL*8, OPTIONAL, INTENT(IN) :: mfy_tl(is:ie, js:je+1)
! for ord=4 option
    REAL*8, OPTIONAL, INTENT(IN) :: ppm_fac
    REAL*8, OPTIONAL, INTENT(IN) :: damp_c
    INTEGER, OPTIONAL, INTENT(IN) :: nord
!#ifdef SINGLE_FV
!! Local Single-Precision Copies
!   real ::  crx_r4(is:ie+1,jsd:jed) 
!   real ::  xfx_r4(is:ie+1,jsd:jed) 
!   real ::  cry_r4(isd:ied,js:je+1 ) 
!   real ::  yfx_r4(isd:ied,js:je+1 ) 
!   real :: ra_x_r4(is:ie,jsd:jed)
!   real :: ra_y_r4(isd:ied,js:je)
!   real ::    q_r4(isd:ied,jsd:jed) 
!   real ::  mfx_r4(is:ie+1,js:je  ) 
!   real ::  mfy_r4(is:ie  ,js:je+1) 
!   real ::   fx_r4(is:ie+1 ,js:je) 
!   real ::   fy_r4(is:ie,   js:je+1 )  
!   real :: ppm_fac_r4
!   real :: damp_c_r4
!
!! Double-to-Single Copies of Inputs  
!        ppm_fac_r4 = ppm_fac
!         damp_c_r4 = damp_c
!            crx_r4 = crx
!            xfx_r4 = xfx
!            cry_r4 = cry
!            yfx_r4 = yfx
!           ra_x_r4 = ra_x
!           ra_y_r4 = ra_y
!              q_r4 = q
!   if (present(mfx))  mfx_r4 = mfx
!   if (present(mfy))  mfy_r4 = mfy
!
!   call fv_tp_2d_compute(q_r4, crx_r4, cry_r4, npx, npy, hord, fx_r4, fy_r4, &
!                     xfx_r4, yfx_r4, area, ra_x_r4, ra_y_r4, mfx_r4, mfy_r4, ppm_fac_r4, nord, damp_c_r4)
! 
!! Single-to-Double Copies of Outputs  
!              q = q_r4
!             fx = fx_r4
!             fy = fy_r4
!#else
    CALL FV_TP_2D_COMPUTE_TLM(q, q_tl, crx, crx_tl, cry, cry_tl, npx, &
&                       npy, hord, fx, fx_tl, fy, fy_tl, xfx, xfx_tl, &
&                       yfx, yfx_tl, area, ra_x, ra_x_tl, ra_y, ra_y_tl&
&                       , mfx, mfx_tl, mfy, mfy_tl, ppm_fac, nord, &
&                       damp_c)
  END SUBROUTINE FV_TP_2D_R8_TLM

  SUBROUTINE FV_TP_2D_COMPUTE_TLM(q, q_tl, crx, crx_tl, cry, cry_tl, npx&
&   , npy, hord, fx, fx_tl, fy, fy_tl, xfx, xfx_tl, yfx, yfx_tl, area, &
&   ra_x, ra_x_tl, ra_y, ra_y_tl, mfx, mfx_tl, mfy, mfy_tl, ppm_fac, &
&   nord, damp_c)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx, npy
    INTEGER, INTENT(IN) :: hord
!
    REAL, INTENT(IN) :: crx(is:ie+1, jsd:jed)
    REAL, INTENT(IN) :: crx_tl(is:ie+1, jsd:jed)
!
    REAL, INTENT(IN) :: xfx(is:ie+1, jsd:jed)
    REAL, INTENT(IN) :: xfx_tl(is:ie+1, jsd:jed)
!
    REAL, INTENT(IN) :: cry(isd:ied, js:je+1)
    REAL, INTENT(IN) :: cry_tl(isd:ied, js:je+1)
!
    REAL, INTENT(IN) :: yfx(isd:ied, js:je+1)
    REAL, INTENT(IN) :: yfx_tl(isd:ied, js:je+1)
    REAL(g_precision), INTENT(IN) :: area(isd:ied, jsd:jed)
    REAL, INTENT(IN) :: ra_x(is:ie, jsd:jed)
    REAL, INTENT(IN) :: ra_x_tl(is:ie, jsd:jed)
    REAL, INTENT(IN) :: ra_y(isd:ied, js:je)
    REAL, INTENT(IN) :: ra_y_tl(isd:ied, js:je)
! transported scalar
    REAL, INTENT(INOUT) :: q(isd:ied, jsd:jed)
    REAL, INTENT(INOUT) :: q_tl(isd:ied, jsd:jed)
! Flux in x ( E )
    REAL, INTENT(OUT) :: fx(is:ie+1, js:je)
    REAL, INTENT(OUT) :: fx_tl(is:ie+1, js:je)
! Flux in y ( N )
    REAL, INTENT(OUT) :: fy(is:ie, js:je+1)
    REAL, INTENT(OUT) :: fy_tl(is:ie, js:je+1)
! optional Arguments:
! Mass Flux X-Dir
    REAL, OPTIONAL, INTENT(IN) :: mfx(is:ie+1, js:je)
    REAL, OPTIONAL, INTENT(IN) :: mfx_tl(is:ie+1, js:je)
! Mass Flux Y-Dir
    REAL, OPTIONAL, INTENT(IN) :: mfy(is:ie, js:je+1)
    REAL, OPTIONAL, INTENT(IN) :: mfy_tl(is:ie, js:je+1)
! for ord=4 option
    REAL, OPTIONAL, INTENT(IN) :: ppm_fac
    REAL, OPTIONAL, INTENT(IN) :: damp_c
    INTEGER, OPTIONAL, INTENT(IN) :: nord
! Local:
    INTEGER :: ord, ord_in
    REAL :: q_i(isd:ied, js:je)
    REAL :: q_i_tl(isd:ied, js:je)
    REAL :: q_j(is:ie, jsd:jed)
    REAL :: q_j_tl(is:ie, jsd:jed)
    REAL :: fx2(is:ie+1, jsd:jed)
    REAL :: fx2_tl(is:ie+1, jsd:jed)
    REAL :: fy2(isd:ied, js:je+1)
    REAL :: fy2_tl(isd:ied, js:je+1)
    REAL :: fyy(isd:ied, js:je+1)
    REAL :: fyy_tl(isd:ied, js:je+1)
    REAL :: fx1(is:ie+1)
    REAL :: fx1_tl(is:ie+1)
    REAL :: ppm_limiter, damp
    INTEGER :: i, j
    INTRINSIC ABS
    INTRINSIC PRESENT
    REAL :: pwx1
    INTEGER :: pwy1
    IF (hord .LT. 0) THEN
! more dissipation
      ord_in = 2
      IF (hord .GE. 0.) THEN
        ord = hord
      ELSE
        ord = -hord
      END IF
    ELSE
      ord_in = hord
      ord = hord
    END IF
    IF (PRESENT(ppm_fac)) THEN
      ppm_limiter = ppm_fac
    ELSE
      ppm_limiter = 2.0
    END IF
    CALL COPY_CORNERS_TLM(q, q_tl, npx, npy, 2)
    fy2_tl = 0.0_8
    CALL YTP_TLM(fy2, fy2_tl, q, q_tl, cry, cry_tl, ord_in, isd, ied, js&
&          , je, npx, npy, ppm_limiter)
    fyy_tl = 0.0_8
    DO j=js,je+1
      DO i=isd,ied
        fyy_tl(i, j) = yfx_tl(i, j)*fy2(i, j) + yfx(i, j)*fy2_tl(i, j)
        fyy(i, j) = yfx(i, j)*fy2(i, j)
      END DO
    END DO
    q_i_tl = 0.0_8
    DO j=js,je
      DO i=isd,ied
        q_i_tl(i, j) = ((area(i, j)*q_tl(i, j)+fyy_tl(i, j)-fyy_tl(i, j+&
&         1))*ra_y(i, j)-(q(i, j)*area(i, j)+fyy(i, j)-fyy(i, j+1))*&
&         ra_y_tl(i, j))/ra_y(i, j)**2
        q_i(i, j) = (q(i, j)*area(i, j)+fyy(i, j)-fyy(i, j+1))/ra_y(i, j&
&         )
      END DO
    END DO
    CALL XTP_TLM(fx, fx_tl, q_i, q_i_tl, crx(is, js), crx_tl(is, js), &
&          ord, is, ie, js, je, npx, npy, ppm_limiter)
    CALL COPY_CORNERS_TLM(q, q_tl, npx, npy, 1)
    fx2_tl = 0.0_8
    CALL XTP_TLM(fx2, fx2_tl, q, q_tl, crx, crx_tl, ord_in, is, ie, jsd&
&          , jed, npx, npy, ppm_limiter)
    q_j_tl = 0.0_8
    fx1_tl = 0.0_8
    DO j=jsd,jed
      DO i=is,ie+1
        fx1_tl(i) = xfx_tl(i, j)*fx2(i, j) + xfx(i, j)*fx2_tl(i, j)
        fx1(i) = xfx(i, j)*fx2(i, j)
      END DO
      DO i=is,ie
        q_j_tl(i, j) = ((area(i, j)*q_tl(i, j)+fx1_tl(i)-fx1_tl(i+1))*&
&         ra_x(i, j)-(q(i, j)*area(i, j)+fx1(i)-fx1(i+1))*ra_x_tl(i, j))&
&         /ra_x(i, j)**2
        q_j(i, j) = (q(i, j)*area(i, j)+fx1(i)-fx1(i+1))/ra_x(i, j)
      END DO
    END DO
    CALL YTP_TLM(fy, fy_tl, q_j, q_j_tl, cry, cry_tl, ord, is, ie, js, &
&          je, npx, npy, ppm_limiter)
!----------------
! Flux averaging:
!----------------
    IF (PRESENT(mfx) .AND. PRESENT(mfy)) THEN
!---------------------------------
! For transport of pt and tracers
!---------------------------------
      DO j=js,je
        DO i=is,ie+1
          fx_tl(i, j) = 0.5*((fx_tl(i, j)+fx2_tl(i, j))*mfx(i, j)+(fx(i&
&           , j)+fx2(i, j))*mfx_tl(i, j))
          fx(i, j) = 0.5*(fx(i, j)+fx2(i, j))*mfx(i, j)
        END DO
      END DO
      DO j=js,je+1
        DO i=is,ie
          fy_tl(i, j) = 0.5*((fy_tl(i, j)+fy2_tl(i, j))*mfy(i, j)+(fy(i&
&           , j)+fy2(i, j))*mfy_tl(i, j))
          fy(i, j) = 0.5*(fy(i, j)+fy2(i, j))*mfy(i, j)
        END DO
      END DO
      IF (PRESENT(nord) .AND. PRESENT(damp_c)) THEN
        IF (damp_c .GT. 1.e-4) THEN
          pwx1 = damp_c*da_min
          pwy1 = nord + 1
          damp = pwx1**pwy1
          CALL DELN_FLUX_TLM(nord, npx, npy, damp, q, q_tl, fx, fx_tl, &
&                      fy, fy_tl, mfx, mfy)
        END IF
      END IF
    ELSE
!---------------------------------
! For transport of delp, vorticity
!---------------------------------
      DO j=js,je
        DO i=is,ie+1
          fx_tl(i, j) = 0.5*((fx_tl(i, j)+fx2_tl(i, j))*xfx(i, j)+(fx(i&
&           , j)+fx2(i, j))*xfx_tl(i, j))
          fx(i, j) = 0.5*(fx(i, j)+fx2(i, j))*xfx(i, j)
        END DO
      END DO
      DO j=js,je+1
        DO i=is,ie
          fy_tl(i, j) = 0.5*((fy_tl(i, j)+fy2_tl(i, j))*yfx(i, j)+(fy(i&
&           , j)+fy2(i, j))*yfx_tl(i, j))
          fy(i, j) = 0.5*(fy(i, j)+fy2(i, j))*yfx(i, j)
        END DO
      END DO
      IF (PRESENT(nord) .AND. PRESENT(damp_c)) THEN
        IF (damp_c .GT. 1.e-4) THEN
          pwx1 = damp_c*da_min
          pwy1 = nord + 1
          damp = pwx1**pwy1
          CALL DELN_FLUX_TLM(nord, npx, npy, damp, q, q_tl, fx, fx_tl, &
&                      fy, fy_tl, xfx(is:ie+1, js:je), yfx(is:ie, js:je+&
&                      1))
        END IF
      END IF
    END IF
  END SUBROUTINE FV_TP_2D_COMPUTE_TLM

  SUBROUTINE XTP_TLM(fx, fx_tl, q, q_tl, c, c_tl, iord, ifirst, ilast, &
&   jfirst, jlast, npx, npy, ppm_limiter)
    IMPLICIT NONE
!  X-Dir strip
    INTEGER, INTENT(IN) :: ifirst, ilast
!  Y-Dir strip
    INTEGER, INTENT(IN) :: jfirst, jlast
    INTEGER, INTENT(IN) :: npx, npy
    INTEGER, INTENT(IN) :: iord
! Courant numbers
    REAL, INTENT(IN) :: c(is:ie+1, jfirst:jlast)
    REAL, INTENT(IN) :: c_tl(is:ie+1, jfirst:jlast)
    REAL, INTENT(IN) :: q(isd:ied, jfirst:jlast)
    REAL, INTENT(IN) :: q_tl(isd:ied, jfirst:jlast)
    REAL, INTENT(IN) :: ppm_limiter
    REAL, INTENT(OUT) :: fx(ifirst:ilast+1, jfirst:jlast)
    REAL, INTENT(OUT) :: fx_tl(ifirst:ilast+1, jfirst:jlast)
! Local:
    REAL :: x0, x1
    INTEGER :: i, j
    IF (iord .EQ. 1) THEN
      DO j=jfirst,jlast
        DO i=ifirst,ilast+1
          IF (c(i, j) .GT. 0.) THEN
            fx_tl(i, j) = q_tl(i-1, j)
            fx(i, j) = q(i-1, j)
          ELSE
            fx_tl(i, j) = q_tl(i, j)
            fx(i, j) = q(i, j)
          END IF
        END DO
      END DO
    ELSE IF (iord .EQ. 333) THEN
      DO j=jfirst,jlast
        DO i=ifirst,ilast+1
!At edges keep to first order
          IF (i .EQ. ifirst .OR. i .EQ. ilast + 1) THEN
            IF (c(i, j) .GT. 0.) THEN
              fx_tl(i, j) = q_tl(i-1, j)
              fx(i, j) = q(i-1, j)
            ELSE
              fx_tl(i, j) = q_tl(i, j)
              fx(i, j) = q(i, j)
            END IF
          ELSE IF (c(i, j) .GT. 0.) THEN
!Otherwise use the third order scheme
            fx_tl(i, j) = (2.0*q_tl(i, j)+5.0*q_tl(i-1, j)-q_tl(i-2, j))&
&             /6.0 - 0.5*(c_tl(i, j)*(q(i, j)-q(i-1, j))+c(i, j)*(q_tl(i&
&             , j)-q_tl(i-1, j))) + (c_tl(i, j)*c(i, j)+c(i, j)*c_tl(i, &
&             j))*(q(i, j)-2.0*q(i-1, j)+q(i-2, j))/6.0 + c(i, j)**2*(&
&             q_tl(i, j)-2.0*q_tl(i-1, j)+q_tl(i-2, j))/6.0
            fx(i, j) = (2.0*q(i, j)+5.0*q(i-1, j)-q(i-2, j))/6.0 - 0.5*c&
&             (i, j)*(q(i, j)-q(i-1, j)) + c(i, j)*c(i, j)/6.0*(q(i, j)-&
&             2.0*q(i-1, j)+q(i-2, j))
          ELSE
            fx_tl(i, j) = (2.0*q_tl(i-1, j)+5.0*q_tl(i, j)-q_tl(i+1, j))&
&             /6.0 - 0.5*(c_tl(i, j)*(q(i, j)-q(i-1, j))+c(i, j)*(q_tl(i&
&             , j)-q_tl(i-1, j))) + (c_tl(i, j)*c(i, j)+c(i, j)*c_tl(i, &
&             j))*(q(i+1, j)-2.0*q(i, j)+q(i-1, j))/6.0 + c(i, j)**2*(&
&             q_tl(i+1, j)-2.0*q_tl(i, j)+q_tl(i-1, j))/6.0
            fx(i, j) = (2.0*q(i-1, j)+5.0*q(i, j)-q(i+1, j))/6.0 - 0.5*c&
&             (i, j)*(q(i, j)-q(i-1, j)) + c(i, j)*c(i, j)/6.0*(q(i+1, j&
&             )-2.0*q(i, j)+q(i-1, j))
          END IF
        END DO
      END DO
    ELSE
      CALL FXPPM_TLM(c, c_tl, q, q_tl, fx, fx_tl, iord, ifirst, ilast, &
&              jfirst, jlast, npx, npy, ppm_limiter)
    END IF
  END SUBROUTINE XTP_TLM

  SUBROUTINE YTP_TLM(fy, fy_tl, q, q_tl, c, c_tl, jord, ifirst, ilast, &
&   jfirst, jlast, npx, npy, ppm_limiter)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx, npy
!  X-Dir strip
    INTEGER, INTENT(IN) :: ifirst, ilast
!  Y-Dir strip
    INTEGER, INTENT(IN) :: jfirst, jlast
    INTEGER, INTENT(IN) :: jord
    REAL, INTENT(IN) :: q(ifirst:ilast, jfirst-ng:jlast+ng)
    REAL, INTENT(IN) :: q_tl(ifirst:ilast, jfirst-ng:jlast+ng)
! Courant number
    REAL, INTENT(IN) :: c(isd:ied, js:je+1)
    REAL, INTENT(IN) :: c_tl(isd:ied, js:je+1)
!  Flux
    REAL, INTENT(OUT) :: fy(ifirst:ilast, jfirst:jlast+1)
    REAL, INTENT(OUT) :: fy_tl(ifirst:ilast, jfirst:jlast+1)
    REAL, INTENT(IN) :: ppm_limiter
! !LOCAL VARIABLES:
    REAL :: x0, x1
    INTEGER :: i, j
    IF (jord .EQ. 1) THEN
      DO j=jfirst,jlast+1
        DO i=ifirst,ilast
          IF (c(i, j) .GT. 0.) THEN
            fy_tl(i, j) = q_tl(i, j-1)
            fy(i, j) = q(i, j-1)
          ELSE
            fy_tl(i, j) = q_tl(i, j)
            fy(i, j) = q(i, j)
          END IF
        END DO
      END DO
    ELSE IF (jord .EQ. 333) THEN
      DO j=jfirst,jlast+1
        DO i=ifirst,ilast
!At edges keep to first order
          IF (j .EQ. jfirst .OR. j .EQ. jlast + 1) THEN
            IF (c(i, j) .GT. 0.) THEN
              fy_tl(i, j) = q_tl(i, j-1)
              fy(i, j) = q(i, j-1)
            ELSE
              fy_tl(i, j) = q_tl(i, j)
              fy(i, j) = q(i, j)
            END IF
          ELSE IF (c(i, j) .GT. 0.) THEN
!Otherwise use the third order scheme
            fy_tl(i, j) = (2.0*q_tl(i, j)+5.0*q_tl(i, j-1)-q_tl(i, j-2))&
&             /6.0 - 0.5*(c_tl(i, j)*(q(i, j)-q(i, j-1))+c(i, j)*(q_tl(i&
&             , j)-q_tl(i, j-1))) + (c_tl(i, j)*c(i, j)+c(i, j)*c_tl(i, &
&             j))*(q(i, j)-2.0*q(i, j-1)+q(i, j-2))/6.0 + c(i, j)**2*(&
&             q_tl(i, j)-2.0*q_tl(i, j-1)+q_tl(i, j-2))/6.0
            fy(i, j) = (2.0*q(i, j)+5.0*q(i, j-1)-q(i, j-2))/6.0 - 0.5*c&
&             (i, j)*(q(i, j)-q(i, j-1)) + c(i, j)*c(i, j)/6.0*(q(i, j)-&
&             2.0*q(i, j-1)+q(i, j-2))
          ELSE
            fy_tl(i, j) = (2.0*q_tl(i, j-1)+5.0*q_tl(i, j)-q_tl(i, j+1))&
&             /6.0 - 0.5*(c_tl(i, j)*(q(i, j)-q(i, j-1))+c(i, j)*(q_tl(i&
&             , j)-q_tl(i, j-1))) + (c_tl(i, j)*c(i, j)+c(i, j)*c_tl(i, &
&             j))*(q(i, j+1)-2.0*q(i, j)+q(i, j-1))/6.0 + c(i, j)**2*(&
&             q_tl(i, j+1)-2.0*q_tl(i, j)+q_tl(i, j-1))/6.0
            fy(i, j) = (2.0*q(i, j-1)+5.0*q(i, j)-q(i, j+1))/6.0 - 0.5*c&
&             (i, j)*(q(i, j)-q(i, j-1)) + c(i, j)*c(i, j)/6.0*(q(i, j+1&
&             )-2.0*q(i, j)+q(i, j-1))
          END IF
        END DO
      END DO
    ELSE
      CALL FYPPM_TLM(c, c_tl, q, q_tl, fy, fy_tl, jord, ifirst, ilast, &
&              jfirst, jlast, npx, npy, ppm_limiter)
    END IF
  END SUBROUTINE YTP_TLM

  SUBROUTINE FXPPM_TLM(c, c_tl, q, q_tl, flux, flux_tl, iord, ifirst, &
&   ilast, jfirst, jlast, npx, npy, ppm_limiter)
    IMPLICIT NONE
! !INPUT PARAMETERS:
!  X-Dir strip
    INTEGER, INTENT(IN) :: ifirst, ilast
!  Y-Dir strip
    INTEGER, INTENT(IN) :: jfirst, jlast
    INTEGER, INTENT(IN) :: iord
    INTEGER, INTENT(IN) :: npx, npy
    REAL, INTENT(IN) :: q(ifirst-ng:ilast+ng, jfirst:jlast)
    REAL, INTENT(IN) :: q_tl(ifirst-ng:ilast+ng, jfirst:jlast)
! Courant   N (like FLUX)
    REAL, INTENT(IN) :: c(ifirst:ilast+1, jfirst:jlast)
    REAL, INTENT(IN) :: c_tl(ifirst:ilast+1, jfirst:jlast)
    REAL, INTENT(IN) :: ppm_limiter
! !OUTPUT PARAMETERS:
!  Flux
    REAL, INTENT(OUT) :: flux(ifirst:ilast+1, jfirst:jlast)
    REAL, INTENT(OUT) :: flux_tl(ifirst:ilast+1, jfirst:jlast)
! Local
    REAL :: dm1(ifirst-2:ilast+2)
    REAL :: al(ifirst-1:ilast+2)
    REAL :: bl(ifirst-1:ilast+1)
    REAL :: bl_tl(ifirst-1:ilast+1)
    REAL :: br(ifirst-1:ilast+1)
    REAL :: br_tl(ifirst-1:ilast+1)
    REAL :: dq(ifirst-3:ilast+2)
    REAL :: dl, dr, pmp, lac, ct, qe
    REAL :: xt, x1, x0
    REAL :: xt_tl
    INTEGER :: i, j, is3, ie3, ie2, it
    INTRINSIC MAX
    INTRINSIC MIN
    REAL :: y1_tl
    REAL :: y2_tl
    REAL :: y3_tl
    REAL :: y4_tl
    REAL :: y5_tl
    REAL :: y6_tl
    REAL :: y7_tl
    REAL :: y8_tl
    REAL :: y8
    REAL :: y7
    REAL :: y6
    REAL :: y5
    REAL :: y4
    REAL :: y3
    REAL :: y2
    REAL :: y1
    IF (3 .LT. is - 1) THEN
      is3 = is - 1
    ELSE
      is3 = 3
    END IF
    IF (npx - 3 .GT. ie + 1) THEN
      ie3 = ie + 1
    ELSE
      ie3 = npx - 3
    END IF
    IF (npx - 2 .GT. ie + 2) THEN
      ie2 = ie + 2
    ELSE
      ie2 = npx - 2
    END IF
    IF (iord .EQ. 6) THEN
      bl_tl = 0.0_8
      br_tl = 0.0_8
      DO j=jfirst,jlast
! Non-monotonic "5th order" scheme (not really 5th order)
        DO i=is3,ie3
          bl_tl(i) = b5*q_tl(i-2, j) + b4*q_tl(i-1, j) + b3*q_tl(i, j) +&
&           b2*q_tl(i+1, j) + b1*q_tl(i+2, j)
          bl(i) = b5*q(i-2, j) + b4*q(i-1, j) + b3*q(i, j) + b2*q(i+1, j&
&           ) + b1*q(i+2, j)
          br_tl(i) = b1*q_tl(i-2, j) + b2*q_tl(i-1, j) + b3*q_tl(i, j) +&
&           b4*q_tl(i+1, j) + b5*q_tl(i+2, j)
          br(i) = b1*q(i-2, j) + b2*q(i-1, j) + b3*q(i, j) + b4*q(i+1, j&
&           ) + b5*q(i+2, j)
        END DO
!--------------
! fix the edges
!--------------
        IF (is .EQ. 1) THEN
          br_tl(2) = p1*(q_tl(2, j)+q_tl(3, j)) + p2*(q_tl(1, j)+q_tl(4&
&           , j)) - q_tl(2, j)
          br(2) = p1*(q(2, j)+q(3, j)) + p2*(q(1, j)+q(4, j)) - q(2, j)
          xt_tl = 0.5*((2.*dxa(1, j)+dxa(2, j))*(q_tl(0, j)+q_tl(1, j))-&
&           dxa(1, j)*(q_tl(-1, j)+q_tl(2, j)))/(dxa(1, j)+dxa(2, j))
          xt = 0.5*((2.*dxa(1, j)+dxa(2, j))*(q(0, j)+q(1, j))-dxa(1, j)&
&           *(q(-1, j)+q(2, j)))/(dxa(1, j)+dxa(2, j))
          bl_tl(1) = xt_tl - q_tl(1, j)
          bl(1) = xt - q(1, j)
          br_tl(0) = xt_tl - q_tl(0, j)
          br(0) = xt - q(0, j)
          xt_tl = c1*q_tl(-2, j) + c2*q_tl(-1, j) + c3*q_tl(0, j)
          xt = c1*q(-2, j) + c2*q(-1, j) + c3*q(0, j)
          IF (q(-1, j) .GT. q(0, j)) THEN
            y1_tl = q_tl(0, j)
            y1 = q(0, j)
          ELSE
            y1_tl = q_tl(-1, j)
            y1 = q(-1, j)
          END IF
          IF (xt .LT. y1) THEN
            xt_tl = y1_tl
            xt = y1
          ELSE
            xt = xt
          END IF
          IF (q(-1, j) .LT. q(0, j)) THEN
            y2_tl = q_tl(0, j)
            y2 = q(0, j)
          ELSE
            y2_tl = q_tl(-1, j)
            y2 = q(-1, j)
          END IF
          IF (xt .GT. y2) THEN
            xt_tl = y2_tl
            xt = y2
          ELSE
            xt = xt
          END IF
          bl_tl(0) = xt_tl - q_tl(0, j)
          bl(0) = xt - q(0, j)
          xt_tl = c3*q_tl(1, j) + c2*q_tl(2, j) + c1*q_tl(3, j)
          xt = c3*q(1, j) + c2*q(2, j) + c1*q(3, j)
          IF (q(1, j) .GT. q(2, j)) THEN
            y3_tl = q_tl(2, j)
            y3 = q(2, j)
          ELSE
            y3_tl = q_tl(1, j)
            y3 = q(1, j)
          END IF
          IF (xt .LT. y3) THEN
            xt_tl = y3_tl
            xt = y3
          ELSE
            xt = xt
          END IF
          IF (q(1, j) .LT. q(2, j)) THEN
            y4_tl = q_tl(2, j)
            y4 = q(2, j)
          ELSE
            y4_tl = q_tl(1, j)
            y4 = q(1, j)
          END IF
          IF (xt .GT. y4) THEN
            xt_tl = y4_tl
            xt = y4
          ELSE
            xt = xt
          END IF
          br_tl(1) = xt_tl - q_tl(1, j)
          br(1) = xt - q(1, j)
          bl_tl(2) = xt_tl - q_tl(2, j)
          bl(2) = xt - q(2, j)
        END IF
        IF (ie + 1 .EQ. npx) THEN
          bl_tl(npx-2) = p1*(q_tl(npx-2, j)+q_tl(npx-3, j)) + p2*(q_tl(&
&           npx-4, j)+q_tl(npx-1, j)) - q_tl(npx-2, j)
          bl(npx-2) = p1*(q(npx-2, j)+q(npx-3, j)) + p2*(q(npx-4, j)+q(&
&           npx-1, j)) - q(npx-2, j)
          xt_tl = 0.5*((2.*dxa(npx-1, j)+dxa(npx-2, j))*(q_tl(npx-1, j)+&
&           q_tl(npx, j))-dxa(npx-1, j)*(q_tl(npx-2, j)+q_tl(npx+1, j)))&
&           /(dxa(npx-1, j)+dxa(npx-2, j))
          xt = 0.5*((2.*dxa(npx-1, j)+dxa(npx-2, j))*(q(npx-1, j)+q(npx&
&           , j))-dxa(npx-1, j)*(q(npx-2, j)+q(npx+1, j)))/(dxa(npx-1, j&
&           )+dxa(npx-2, j))
          br_tl(npx-1) = xt_tl - q_tl(npx-1, j)
          br(npx-1) = xt - q(npx-1, j)
          bl_tl(npx) = xt_tl - q_tl(npx, j)
          bl(npx) = xt - q(npx, j)
          xt_tl = c3*q_tl(npx, j) + c2*q_tl(npx+1, j) + c1*q_tl(npx+2, j&
&           )
          xt = c3*q(npx, j) + c2*q(npx+1, j) + c1*q(npx+2, j)
          IF (q(npx, j) .GT. q(npx+1, j)) THEN
            y5_tl = q_tl(npx+1, j)
            y5 = q(npx+1, j)
          ELSE
            y5_tl = q_tl(npx, j)
            y5 = q(npx, j)
          END IF
          IF (xt .LT. y5) THEN
            xt_tl = y5_tl
            xt = y5
          ELSE
            xt = xt
          END IF
          IF (q(npx, j) .LT. q(npx+1, j)) THEN
            y6_tl = q_tl(npx+1, j)
            y6 = q(npx+1, j)
          ELSE
            y6_tl = q_tl(npx, j)
            y6 = q(npx, j)
          END IF
          IF (xt .GT. y6) THEN
            xt_tl = y6_tl
            xt = y6
          ELSE
            xt = xt
          END IF
          br_tl(npx) = xt_tl - q_tl(npx, j)
          br(npx) = xt - q(npx, j)
          xt_tl = c1*q_tl(npx-3, j) + c2*q_tl(npx-2, j) + c3*q_tl(npx-1&
&           , j)
          xt = c1*q(npx-3, j) + c2*q(npx-2, j) + c3*q(npx-1, j)
          IF (q(npx-2, j) .GT. q(npx-1, j)) THEN
            y7_tl = q_tl(npx-1, j)
            y7 = q(npx-1, j)
          ELSE
            y7_tl = q_tl(npx-2, j)
            y7 = q(npx-2, j)
          END IF
          IF (xt .LT. y7) THEN
            xt_tl = y7_tl
            xt = y7
          ELSE
            xt = xt
          END IF
          IF (q(npx-2, j) .LT. q(npx-1, j)) THEN
            y8_tl = q_tl(npx-1, j)
            y8 = q(npx-1, j)
          ELSE
            y8_tl = q_tl(npx-2, j)
            y8 = q(npx-2, j)
          END IF
          IF (xt .GT. y8) THEN
            xt_tl = y8_tl
            xt = y8
          ELSE
            xt = xt
          END IF
          br_tl(npx-2) = xt_tl - q_tl(npx-2, j)
          br(npx-2) = xt - q(npx-2, j)
          bl_tl(npx-1) = xt_tl - q_tl(npx-1, j)
          bl(npx-1) = xt - q(npx-1, j)
        END IF
        DO i=ifirst,ilast+1
          IF (c(i, j) .GT. 0.) THEN
            flux_tl(i, j) = q_tl(i-1, j) + (1.-c(i, j))*(br_tl(i-1)-c_tl&
&             (i, j)*(bl(i-1)+br(i-1))-c(i, j)*(bl_tl(i-1)+br_tl(i-1))) &
&             - c_tl(i, j)*(br(i-1)-c(i, j)*(bl(i-1)+br(i-1)))
            flux(i, j) = q(i-1, j) + (1.-c(i, j))*(br(i-1)-c(i, j)*(bl(i&
&             -1)+br(i-1)))
          ELSE
            flux_tl(i, j) = q_tl(i, j) + c_tl(i, j)*(bl(i)+c(i, j)*(bl(i&
&             )+br(i))) + (1.+c(i, j))*(bl_tl(i)+c_tl(i, j)*(bl(i)+br(i)&
&             )+c(i, j)*(bl_tl(i)+br_tl(i)))
            flux(i, j) = q(i, j) + (1.+c(i, j))*(bl(i)+c(i, j)*(bl(i)+br&
&             (i)))
          END IF
        END DO
      END DO
    END IF
  END SUBROUTINE FXPPM_TLM

  SUBROUTINE FYPPM_TLM(c, c_tl, q, q_tl, flux, flux_tl, jord, ifirst, &
&   ilast, jfirst, jlast, npx, npy, ppm_limiter)
    IMPLICIT NONE
!  X-Dir strip
    INTEGER, INTENT(IN) :: ifirst, ilast
!  Y-Dir strip
    INTEGER, INTENT(IN) :: jfirst, jlast
    INTEGER, INTENT(IN) :: jord
    INTEGER, INTENT(IN) :: npx, npy
    REAL, INTENT(IN) :: q(ifirst:ilast, jfirst-ng:jlast+ng)
    REAL, INTENT(IN) :: q_tl(ifirst:ilast, jfirst-ng:jlast+ng)
! Courant number
    REAL, INTENT(IN) :: c(isd:ied, js:je+1)
    REAL, INTENT(IN) :: c_tl(isd:ied, js:je+1)
    REAL, INTENT(IN) :: ppm_limiter
!  Flux
    REAL, INTENT(OUT) :: flux(ifirst:ilast, jfirst:jlast+1)
    REAL, INTENT(OUT) :: flux_tl(ifirst:ilast, jfirst:jlast+1)
! Local:
    REAL :: al(ifirst:ilast, jfirst-1:jlast+2)
    REAL :: bl(ifirst:ilast, jfirst-1:jlast+1)
    REAL :: bl_tl(ifirst:ilast, jfirst-1:jlast+1)
    REAL :: br(ifirst:ilast, jfirst-1:jlast+1)
    REAL :: br_tl(ifirst:ilast, jfirst-1:jlast+1)
    REAL :: dq(ifirst:ilast, jfirst-3:jlast+2)
    REAL :: dl, dr, pmp, lac, ct, qe
    REAL :: xt, x0, x1
    REAL :: xt_tl
    INTEGER :: i, j, js3, je3, je2, jt
    INTRINSIC MAX
    INTRINSIC MIN
    REAL :: y1_tl
    REAL :: y2_tl
    INTEGER :: min1
    REAL :: y3_tl
    REAL :: y4_tl
    REAL :: y5_tl
    REAL :: y6_tl
    REAL :: y7_tl
    REAL :: y8_tl
    REAL :: y8
    REAL :: y7
    REAL :: y6
    INTEGER :: max1
    REAL :: y5
    REAL :: y4
    REAL :: y3
    REAL :: y2
    REAL :: y1
    IF (3 .LT. js - 1) THEN
      js3 = js - 1
    ELSE
      js3 = 3
    END IF
    IF (npy - 3 .GT. je + 1) THEN
      je3 = je + 1
    ELSE
      je3 = npy - 3
    END IF
    IF (npy - 2 .GT. je + 2) THEN
      je2 = je + 2
    ELSE
      je2 = npy - 2
    END IF
    IF (jord .EQ. 6) THEN
      IF (3 .LT. js - 1) THEN
        max1 = js - 1
      ELSE
        max1 = 3
      END IF
      IF (npy - 3 .GT. je + 1) THEN
        min1 = je + 1
        bl_tl = 0.0_8
        br_tl = 0.0_8
      ELSE
        min1 = npy - 3
        bl_tl = 0.0_8
        br_tl = 0.0_8
      END IF
! Non-monotonic "5th order" scheme (not really 5th order)
      DO j=max1,min1
        DO i=ifirst,ilast
          bl_tl(i, j) = b5*q_tl(i, j-2) + b4*q_tl(i, j-1) + b3*q_tl(i, j&
&           ) + b2*q_tl(i, j+1) + b1*q_tl(i, j+2)
          bl(i, j) = b5*q(i, j-2) + b4*q(i, j-1) + b3*q(i, j) + b2*q(i, &
&           j+1) + b1*q(i, j+2)
          br_tl(i, j) = b1*q_tl(i, j-2) + b2*q_tl(i, j-1) + b3*q_tl(i, j&
&           ) + b4*q_tl(i, j+1) + b5*q_tl(i, j+2)
          br(i, j) = b1*q(i, j-2) + b2*q(i, j-1) + b3*q(i, j) + b4*q(i, &
&           j+1) + b5*q(i, j+2)
        END DO
      END DO
      IF (js .EQ. 1) THEN
        DO i=ifirst,ilast
!           br(i,2) = al(i,3) - q(i,2)
          br_tl(i, 2) = p1*(q_tl(i, 2)+q_tl(i, 3)) + p2*(q_tl(i, 1)+q_tl&
&           (i, 4)) - q_tl(i, 2)
          br(i, 2) = p1*(q(i, 2)+q(i, 3)) + p2*(q(i, 1)+q(i, 4)) - q(i, &
&           2)
          xt_tl = 0.5*((2.*dya(i, 1)+dya(i, 2))*(q_tl(i, 0)+q_tl(i, 1))-&
&           dya(i, 1)*(q_tl(i, -1)+q_tl(i, 2)))/(dya(i, 1)+dya(i, 2))
          xt = 0.5*((2.*dya(i, 1)+dya(i, 2))*(q(i, 0)+q(i, 1))-dya(i, 1)&
&           *(q(i, -1)+q(i, 2)))/(dya(i, 1)+dya(i, 2))
          bl_tl(i, 1) = xt_tl - q_tl(i, 1)
          bl(i, 1) = xt - q(i, 1)
          br_tl(i, 0) = xt_tl - q_tl(i, 0)
          br(i, 0) = xt - q(i, 0)
!           xt = s14*0.25*(q(i,0)-q(i,-2)) - s11*(q(i,0)-q(i,-1)) + q(i,0)
          xt_tl = c1*q_tl(i, -2) + c2*q_tl(i, -1) + c3*q_tl(i, 0)
          xt = c1*q(i, -2) + c2*q(i, -1) + c3*q(i, 0)
          IF (q(i, -1) .LT. q(i, 0)) THEN
            y1_tl = q_tl(i, 0)
            y1 = q(i, 0)
          ELSE
            y1_tl = q_tl(i, -1)
            y1 = q(i, -1)
          END IF
          IF (xt .GT. y1) THEN
            xt_tl = y1_tl
            xt = y1
          ELSE
            xt = xt
          END IF
          IF (q(i, -1) .GT. q(i, 0)) THEN
            y2_tl = q_tl(i, 0)
            y2 = q(i, 0)
          ELSE
            y2_tl = q_tl(i, -1)
            y2 = q(i, -1)
          END IF
          IF (xt .LT. y2) THEN
            xt_tl = y2_tl
            xt = y2
          ELSE
            xt = xt
          END IF
          bl_tl(i, 0) = xt_tl - q_tl(i, 0)
          bl(i, 0) = xt - q(i, 0)
!           xt = s15*q(i,1) + s11*q(i,2) - s14*0.25*(q(i,3)-q(i,1))
          xt_tl = c3*q_tl(i, 1) + c2*q_tl(i, 2) + c1*q_tl(i, 3)
          xt = c3*q(i, 1) + c2*q(i, 2) + c1*q(i, 3)
          IF (q(i, 1) .LT. q(i, 2)) THEN
            y3_tl = q_tl(i, 2)
            y3 = q(i, 2)
          ELSE
            y3_tl = q_tl(i, 1)
            y3 = q(i, 1)
          END IF
          IF (xt .GT. y3) THEN
            xt_tl = y3_tl
            xt = y3
          ELSE
            xt = xt
          END IF
          IF (q(i, 1) .GT. q(i, 2)) THEN
            y4_tl = q_tl(i, 2)
            y4 = q(i, 2)
          ELSE
            y4_tl = q_tl(i, 1)
            y4 = q(i, 1)
          END IF
          IF (xt .LT. y4) THEN
            xt_tl = y4_tl
            xt = y4
          ELSE
            xt = xt
          END IF
          br_tl(i, 1) = xt_tl - q_tl(i, 1)
          br(i, 1) = xt - q(i, 1)
          bl_tl(i, 2) = xt_tl - q_tl(i, 2)
          bl(i, 2) = xt - q(i, 2)
        END DO
      END IF
      IF (je + 1 .EQ. npy) THEN
        DO i=ifirst,ilast
!           bl(i,npy-2) = al(i,npy-2) - q(i,npy-2)
          bl_tl(i, npy-2) = p1*(q_tl(i, npy-3)+q_tl(i, npy-2)) + p2*(&
&           q_tl(i, npy-4)+q_tl(i, npy-1)) - q_tl(i, npy-2)
          bl(i, npy-2) = p1*(q(i, npy-3)+q(i, npy-2)) + p2*(q(i, npy-4)+&
&           q(i, npy-1)) - q(i, npy-2)
          xt_tl = 0.5*((2.*dya(i, npy-1)+dya(i, npy-2))*(q_tl(i, npy-1)+&
&           q_tl(i, npy))-dya(i, npy-1)*(q_tl(i, npy-2)+q_tl(i, npy+1)))&
&           /(dya(i, npy-1)+dya(i, npy-2))
          xt = 0.5*((2.*dya(i, npy-1)+dya(i, npy-2))*(q(i, npy-1)+q(i, &
&           npy))-dya(i, npy-1)*(q(i, npy-2)+q(i, npy+1)))/(dya(i, npy-1&
&           )+dya(i, npy-2))
          br_tl(i, npy-1) = xt_tl - q_tl(i, npy-1)
          br(i, npy-1) = xt - q(i, npy-1)
          bl_tl(i, npy) = xt_tl - q_tl(i, npy)
          bl(i, npy) = xt - q(i, npy)
!           xt = s11*(q(i,npy+1)-q(i,npy)) - s14*0.25*(q(i,npy+2)-q(i,npy)) + q(i,npy)
          xt_tl = c3*q_tl(i, npy) + c2*q_tl(i, npy+1) + c1*q_tl(i, npy+2&
&           )
          xt = c3*q(i, npy) + c2*q(i, npy+1) + c1*q(i, npy+2)
          IF (q(i, npy) .LT. q(i, npy+1)) THEN
            y5_tl = q_tl(i, npy+1)
            y5 = q(i, npy+1)
          ELSE
            y5_tl = q_tl(i, npy)
            y5 = q(i, npy)
          END IF
          IF (xt .GT. y5) THEN
            xt_tl = y5_tl
            xt = y5
          ELSE
            xt = xt
          END IF
          IF (q(i, npy) .GT. q(i, npy+1)) THEN
            y6_tl = q_tl(i, npy+1)
            y6 = q(i, npy+1)
          ELSE
            y6_tl = q_tl(i, npy)
            y6 = q(i, npy)
          END IF
          IF (xt .LT. y6) THEN
            xt_tl = y6_tl
            xt = y6
          ELSE
            xt = xt
          END IF
          br_tl(i, npy) = xt_tl - q_tl(i, npy)
          br(i, npy) = xt - q(i, npy)
!           xt = s15*q(i,npy-1) + s11*q(i,npy-2) + s14*0.25*(q(i,npy-1)-q(i,npy-3))
          xt_tl = c1*q_tl(i, npy-3) + c2*q_tl(i, npy-2) + c3*q_tl(i, npy&
&           -1)
          xt = c1*q(i, npy-3) + c2*q(i, npy-2) + c3*q(i, npy-1)
          IF (q(i, npy-2) .LT. q(i, npy-1)) THEN
            y7_tl = q_tl(i, npy-1)
            y7 = q(i, npy-1)
          ELSE
            y7_tl = q_tl(i, npy-2)
            y7 = q(i, npy-2)
          END IF
          IF (xt .GT. y7) THEN
            xt_tl = y7_tl
            xt = y7
          ELSE
            xt = xt
          END IF
          IF (q(i, npy-2) .GT. q(i, npy-1)) THEN
            y8_tl = q_tl(i, npy-1)
            y8 = q(i, npy-1)
          ELSE
            y8_tl = q_tl(i, npy-2)
            y8 = q(i, npy-2)
          END IF
          IF (xt .LT. y8) THEN
            xt_tl = y8_tl
            xt = y8
          ELSE
            xt = xt
          END IF
          br_tl(i, npy-2) = xt_tl - q_tl(i, npy-2)
          br(i, npy-2) = xt - q(i, npy-2)
          bl_tl(i, npy-1) = xt_tl - q_tl(i, npy-1)
          bl(i, npy-1) = xt - q(i, npy-1)
        END DO
      END IF
      DO j=jfirst,jlast+1
        DO i=ifirst,ilast
          IF (c(i, j) .GT. 0.) THEN
            flux_tl(i, j) = q_tl(i, j-1) + (1.-c(i, j))*(br_tl(i, j-1)-&
&             c_tl(i, j)*(bl(i, j-1)+br(i, j-1))-c(i, j)*(bl_tl(i, j-1)+&
&             br_tl(i, j-1))) - c_tl(i, j)*(br(i, j-1)-c(i, j)*(bl(i, j-&
&             1)+br(i, j-1)))
            flux(i, j) = q(i, j-1) + (1.-c(i, j))*(br(i, j-1)-c(i, j)*(&
&             bl(i, j-1)+br(i, j-1)))
          ELSE
            flux_tl(i, j) = q_tl(i, j) + c_tl(i, j)*(bl(i, j)+c(i, j)*(&
&             bl(i, j)+br(i, j))) + (1.+c(i, j))*(bl_tl(i, j)+c_tl(i, j)&
&             *(bl(i, j)+br(i, j))+c(i, j)*(bl_tl(i, j)+br_tl(i, j)))
            flux(i, j) = q(i, j) + (1.+c(i, j))*(bl(i, j)+c(i, j)*(bl(i&
&             , j)+br(i, j)))
          END IF
        END DO
      END DO
    END IF
  END SUBROUTINE FYPPM_TLM

  SUBROUTINE DELN_FLUX_TLM(nord, npx, npy, damp, q, q_tl, fx, fx_tl, fy&
&   , fy_tl, mfx, mfy)
    IMPLICIT NONE
! Del-n damping for the cell-mean values (A grid)
!------------------
! nord = 0:   del-2
! nord = 1:   del-4
! nord = 2:   del-6
! nord = 3:   del-8 --> requires more ghosting than current
!------------------
! del-n
    INTEGER, INTENT(IN) :: nord
    INTEGER, INTENT(IN) :: npx, npy
    REAL, INTENT(IN) :: damp
! q ghosted on input
    REAL, INTENT(IN) :: q(is-ng:ie+ng, js-ng:je+ng)
    REAL, INTENT(IN) :: q_tl(is-ng:ie+ng, js-ng:je+ng)
! diffusive fluxes: 
    REAL, INTENT(IN) :: mfx(is:ie+1, js:je), mfy(is:ie, js:je+1)
    REAL, INTENT(INOUT) :: fx(is:ie+1, js:je), fy(is:ie, js:je+1)
    REAL, INTENT(INOUT) :: fx_tl(is:ie+1, js:je), fy_tl(is:ie, js:je+1)
! local:
    REAL :: fx2(isd:ied+1, jsd:jed), fy2(isd:ied, jsd:jed+1)
    REAL :: fx2_tl(isd:ied+1, jsd:jed), fy2_tl(isd:ied, jsd:jed+1)
    REAL :: d2(isd:ied, jsd:jed)
    REAL :: d2_tl(isd:ied, jsd:jed)
    INTEGER :: i, j, n, nt
    d2_tl = 0.0_8
    DO j=jsd,jed
      DO i=isd,ied
        d2_tl(i, j) = damp*q_tl(i, j)
        d2(i, j) = damp*q(i, j)
      END DO
    END DO
    IF (nord .GT. 0) THEN
      CALL COPY_CORNERS_TLM(d2, d2_tl, npx, npy, 1)
      fx2_tl = 0.0_8
    ELSE
      fx2_tl = 0.0_8
    END IF
    DO j=js-nord,je+nord
      DO i=is-nord,ie+nord+1
        fx2_tl(i, j) = dy(i, j)*sina_u(i, j)*rdxc(i, j)*(d2_tl(i-1, j)-&
&         d2_tl(i, j))
        fx2(i, j) = dy(i, j)*sina_u(i, j)*(d2(i-1, j)-d2(i, j))*rdxc(i, &
&         j)
      END DO
    END DO
    IF (nord .GT. 0) THEN
      CALL COPY_CORNERS_TLM(d2, d2_tl, npx, npy, 2)
      fy2_tl = 0.0_8
    ELSE
      fy2_tl = 0.0_8
    END IF
    DO j=js-nord,je+nord+1
      DO i=is-nord,ie+nord
        fy2_tl(i, j) = dx(i, j)*sina_v(i, j)*rdyc(i, j)*(d2_tl(i, j-1)-&
&         d2_tl(i, j))
        fy2(i, j) = dx(i, j)*sina_v(i, j)*(d2(i, j-1)-d2(i, j))*rdyc(i, &
&         j)
      END DO
    END DO
    IF (nord .GT. 0) THEN
!----------
! high-order
!----------
      DO n=1,nord
        nt = nord - n
        DO j=js-nt-1,je+nt+1
          DO i=is-nt-1,ie+nt+1
            d2_tl(i, j) = rarea(i, j)*(fx2_tl(i, j)-fx2_tl(i+1, j)+&
&             fy2_tl(i, j)-fy2_tl(i, j+1))
            d2(i, j) = (fx2(i, j)-fx2(i+1, j)+fy2(i, j)-fy2(i, j+1))*&
&             rarea(i, j)
          END DO
        END DO
        CALL COPY_CORNERS_TLM(d2, d2_tl, npx, npy, 1)
        DO j=js-nt,je+nt
          DO i=is-nt,ie+nt+1
            fx2_tl(i, j) = dy(i, j)*sina_u(i, j)*rdxc(i, j)*(d2_tl(i, j)&
&             -d2_tl(i-1, j))
            fx2(i, j) = dy(i, j)*sina_u(i, j)*(d2(i, j)-d2(i-1, j))*rdxc&
&             (i, j)
          END DO
        END DO
        CALL COPY_CORNERS_TLM(d2, d2_tl, npx, npy, 2)
        DO j=js-nt,je+nt+1
          DO i=is-nt,ie+nt
            fy2_tl(i, j) = dx(i, j)*sina_v(i, j)*rdyc(i, j)*(d2_tl(i, j)&
&             -d2_tl(i, j-1))
            fy2(i, j) = dx(i, j)*sina_v(i, j)*(d2(i, j)-d2(i, j-1))*rdyc&
&             (i, j)
          END DO
        END DO
      END DO
    END IF
!---------------------------------------------
! Add the diffusive fluxes to the flux arrays:
!---------------------------------------------
    DO j=js,je
      DO i=is,ie+1
!        fx(i,j) = fx(i,j) + fx2(i,j)*mfx(i,j)
        fx_tl(i, j) = fx_tl(i, j) + fx2_tl(i, j)
        fx(i, j) = fx(i, j) + fx2(i, j)
      END DO
    END DO
    DO j=js,je+1
      DO i=is,ie
!        fy(i,j) = fy(i,j) + fy2(i,j)*mfy(i,j)
        fy_tl(i, j) = fy_tl(i, j) + fy2_tl(i, j)
        fy(i, j) = fy(i, j) + fy2(i, j)
      END DO
    END DO
  END SUBROUTINE DELN_FLUX_TLM
!  Differentiation of copy_corners_r8 in forward (tangent) mode (with options r8):
!   variations   of useful results: q
!   with respect to varying inputs: q
  SUBROUTINE COPY_CORNERS_R8_TLM(q, q_tl, npx, npy, dir)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx, npy, dir
    REAL*8, INTENT(INOUT) :: q(isd:ied, jsd:jed)
    REAL*8, INTENT(INOUT) :: q_tl(isd:ied, jsd:jed)
    INTEGER :: i, j
    IF (dir .EQ. 1) THEN
! XDir:
      IF (sw_corner) THEN
        DO j=1-ng,0
          DO i=1-ng,0
            q_tl(i, j) = q_tl(j, 1-i)
            q(i, j) = q(j, 1-i)
          END DO
        END DO
      END IF
      IF (se_corner) THEN
        DO j=1-ng,0
          DO i=npx,npx+ng-1
            q_tl(i, j) = q_tl(npy-j, i-npx+1)
            q(i, j) = q(npy-j, i-npx+1)
          END DO
        END DO
      END IF
      IF (ne_corner) THEN
        DO j=npy,npy+ng-1
          DO i=npx,npx+ng-1
            q_tl(i, j) = q_tl(j, 2*npx-1-i)
            q(i, j) = q(j, 2*npx-1-i)
          END DO
        END DO
      END IF
      IF (nw_corner) THEN
        DO j=npy,npy+ng-1
          DO i=1-ng,0
            q_tl(i, j) = q_tl(npy-j, i-1+npx)
            q(i, j) = q(npy-j, i-1+npx)
          END DO
        END DO
      END IF
    ELSE IF (dir .EQ. 2) THEN
! YDir:
      IF (sw_corner) THEN
        DO j=1-ng,0
          DO i=1-ng,0
            q_tl(i, j) = q_tl(1-j, i)
            q(i, j) = q(1-j, i)
          END DO
        END DO
      END IF
      IF (se_corner) THEN
        DO j=1-ng,0
          DO i=npx,npx+ng-1
            q_tl(i, j) = q_tl(npy+j-1, npx-i)
            q(i, j) = q(npy+j-1, npx-i)
          END DO
        END DO
      END IF
      IF (ne_corner) THEN
        DO j=npy,npy+ng-1
          DO i=npx,npx+ng-1
            q_tl(i, j) = q_tl(2*npy-1-j, i)
            q(i, j) = q(2*npy-1-j, i)
          END DO
        END DO
      END IF
      IF (nw_corner) THEN
        DO j=npy,npy+ng-1
          DO i=1-ng,0
            q_tl(i, j) = q_tl(j+1-npx, npy-i)
            q(i, j) = q(j+1-npx, npy-i)
          END DO
        END DO
      END IF
    END IF
  END SUBROUTINE COPY_CORNERS_R8_TLM

end module tp_core_tlm_mod
