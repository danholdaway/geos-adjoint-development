module tp_core_adm_mod
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

      INTERFACE copy_corners 
        MODULE PROCEDURE copy_corners_r4
        MODULE PROCEDURE copy_corners_r8
      END INTERFACE

      INTERFACE fv_tp_2d
!        MODULE PROCEDURE fv_tp_2d_r4
        MODULE PROCEDURE fv_tp_2d_r8
      END INTERFACE

      INTERFACE COPY_CORNERS_ADM
        MODULE PROCEDURE COPY_CORNERS_R8_ADM
      END INTERFACE

      INTERFACE FV_TP_2D_ADM
        MODULE PROCEDURE FV_TP_2D_R8_ADM
      END INTERFACE

 private
 public fv_tp_2d, copy_corners, fv_tp_2d_adm, copy_corners_adm

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

  SUBROUTINE FV_TP_2D_R8_ADM(q, q_ad, crx, crx_ad, cry, cry_ad, npx, npy&
&   , hord, fx, fx_ad, fy, fy_ad, xfx, xfx_ad, yfx, yfx_ad, area, ra_x, &
&   ra_x_ad, ra_y, ra_y_ad, mfx, mfx_ad, mfy, mfy_ad, ppm_fac, nord, &
&   damp_c)
    IMPLICIT NONE
!#endif
    INTEGER, INTENT(IN) :: npx, npy
    INTEGER, INTENT(IN) :: hord
!
    REAL*8, INTENT(IN) :: crx(is:ie+1, jsd:jed)
    REAL*8 :: crx_ad(is:ie+1, jsd:jed)
!
    REAL*8, INTENT(IN) :: xfx(is:ie+1, jsd:jed)
    REAL*8 :: xfx_ad(is:ie+1, jsd:jed)
!
    REAL*8, INTENT(IN) :: cry(isd:ied, js:je+1)
    REAL*8 :: cry_ad(isd:ied, js:je+1)
!
    REAL*8, INTENT(IN) :: yfx(isd:ied, js:je+1)
    REAL*8 :: yfx_ad(isd:ied, js:je+1)
    REAL(g_precision), INTENT(IN) :: area(isd:ied, jsd:jed)
    REAL*8, INTENT(IN) :: ra_x(is:ie, jsd:jed)
    REAL*8 :: ra_x_ad(is:ie, jsd:jed)
    REAL*8, INTENT(IN) :: ra_y(isd:ied, js:je)
    REAL*8 :: ra_y_ad(isd:ied, js:je)
! transported scalar
    REAL*8, INTENT(INOUT) :: q(isd:ied, jsd:jed)
    REAL*8 :: q_ad(isd:ied, jsd:jed)
! Flux in x ( E )
    REAL*8 :: fx(is:ie+1, js:je)
    REAL*8 :: fx_ad(is:ie+1, js:je)
! Flux in y ( N )
    REAL*8 :: fy(is:ie, js:je+1)
    REAL*8 :: fy_ad(is:ie, js:je+1)
! optional Arguments:
! Mass Flux X-Dir
    REAL*8, OPTIONAL, INTENT(IN) :: mfx(is:ie+1, js:je)
    REAL*8, OPTIONAL :: mfx_ad(is:ie+1, js:je)
! Mass Flux Y-Dir
    REAL*8, OPTIONAL, INTENT(IN) :: mfy(is:ie, js:je+1)
    REAL*8, OPTIONAL :: mfy_ad(is:ie, js:je+1)
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
    CALL FV_TP_2D_COMPUTE_ADM(q, q_ad, crx, crx_ad, cry, cry_ad, npx, &
&                       npy, hord, fx, fx_ad, fy, fy_ad, xfx, xfx_ad, &
&                       yfx, yfx_ad, area, ra_x, ra_x_ad, ra_y, ra_y_ad&
&                       , mfx, mfx_ad, mfy, mfy_ad, ppm_fac, nord, &
&                       damp_c)
  END SUBROUTINE FV_TP_2D_R8_ADM
  SUBROUTINE FV_TP_2D_R8(q, crx, cry, npx, npy, hord, fx, fy, xfx, yfx, &
&   area, ra_x, ra_y, mfx, mfy, ppm_fac, nord, damp_c)
    IMPLICIT NONE
!#endif
    INTEGER, INTENT(IN) :: npx, npy
    INTEGER, INTENT(IN) :: hord
!
    REAL*8, INTENT(IN) :: crx(is:ie+1, jsd:jed)
!
    REAL*8, INTENT(IN) :: xfx(is:ie+1, jsd:jed)
!
    REAL*8, INTENT(IN) :: cry(isd:ied, js:je+1)
!
    REAL*8, INTENT(IN) :: yfx(isd:ied, js:je+1)
    REAL(g_precision), INTENT(IN) :: area(isd:ied, jsd:jed)
    REAL*8, INTENT(IN) :: ra_x(is:ie, jsd:jed)
    REAL*8, INTENT(IN) :: ra_y(isd:ied, js:je)
! transported scalar
    REAL*8, INTENT(INOUT) :: q(isd:ied, jsd:jed)
! Flux in x ( E )
    REAL*8, INTENT(OUT) :: fx(is:ie+1, js:je)
! Flux in y ( N )
    REAL*8, INTENT(OUT) :: fy(is:ie, js:je+1)
! optional Arguments:
! Mass Flux X-Dir
    REAL*8, OPTIONAL, INTENT(IN) :: mfx(is:ie+1, js:je)
! Mass Flux Y-Dir
    REAL*8, OPTIONAL, INTENT(IN) :: mfy(is:ie, js:je+1)
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
    CALL FV_TP_2D_COMPUTE(q, crx, cry, npx, npy, hord, fx, fy, xfx, yfx&
&                   , area, ra_x, ra_y, mfx, mfy, ppm_fac, nord, damp_c)
  END SUBROUTINE FV_TP_2D_R8

  SUBROUTINE FV_TP_2D_COMPUTE_ADM(q, q_ad, crx, crx_ad, cry, cry_ad, npx&
&   , npy, hord, fx, fx_ad, fy, fy_ad, xfx, xfx_ad, yfx, yfx_ad, area, &
&   ra_x, ra_x_ad, ra_y, ra_y_ad, mfx, mfx_ad, mfy, mfy_ad, ppm_fac, &
&   nord, damp_c)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx, npy
    INTEGER, INTENT(IN) :: hord
!
    REAL, INTENT(IN) :: crx(is:ie+1, jsd:jed)
    REAL :: crx_ad(is:ie+1, jsd:jed)
!
    REAL, INTENT(IN) :: xfx(is:ie+1, jsd:jed)
    REAL :: xfx_ad(is:ie+1, jsd:jed)
!
    REAL, INTENT(IN) :: cry(isd:ied, js:je+1)
    REAL :: cry_ad(isd:ied, js:je+1)
!
    REAL, INTENT(IN) :: yfx(isd:ied, js:je+1)
    REAL :: yfx_ad(isd:ied, js:je+1)
    REAL(g_precision), INTENT(IN) :: area(isd:ied, jsd:jed)
    REAL, INTENT(IN) :: ra_x(is:ie, jsd:jed)
    REAL :: ra_x_ad(is:ie, jsd:jed)
    REAL, INTENT(IN) :: ra_y(isd:ied, js:je)
    REAL :: ra_y_ad(isd:ied, js:je)
! transported scalar
    REAL, INTENT(INOUT) :: q(isd:ied, jsd:jed)
    REAL, INTENT(INOUT) :: q_ad(isd:ied, jsd:jed)
! Flux in x ( E )
    REAL :: fx(is:ie+1, js:je)
    REAL :: fx_ad(is:ie+1, js:je)
! Flux in y ( N )
    REAL :: fy(is:ie, js:je+1)
    REAL :: fy_ad(is:ie, js:je+1)
! optional Arguments:
! Mass Flux X-Dir
    REAL, OPTIONAL, INTENT(IN) :: mfx(is:ie+1, js:je)
    REAL, OPTIONAL :: mfx_ad(is:ie+1, js:je)
! Mass Flux Y-Dir
    REAL, OPTIONAL, INTENT(IN) :: mfy(is:ie, js:je+1)
    REAL, OPTIONAL :: mfy_ad(is:ie, js:je+1)
! for ord=4 option
    REAL, OPTIONAL, INTENT(IN) :: ppm_fac
    REAL, OPTIONAL, INTENT(IN) :: damp_c
    INTEGER, OPTIONAL, INTENT(IN) :: nord
! Local:
    INTEGER :: ord, ord_in
    REAL :: q_i(isd:ied, js:je)
    REAL :: q_i_ad(isd:ied, js:je)
    REAL :: q_j(is:ie, jsd:jed)
    REAL :: q_j_ad(is:ie, jsd:jed)
    REAL :: fx2(is:ie+1, jsd:jed)
    REAL :: fx2_ad(is:ie+1, jsd:jed)
    REAL :: fy2(isd:ied, js:je+1)
    REAL :: fy2_ad(isd:ied, js:je+1)
    REAL :: fyy(isd:ied, js:je+1)
    REAL :: fyy_ad(isd:ied, js:je+1)
    REAL :: fx1(is:ie+1)
    REAL :: fx1_ad(is:ie+1)
    REAL :: ppm_limiter, damp
    INTEGER :: i, j
    INTRINSIC ABS
    INTRINSIC PRESENT
    INTEGER :: branch
    REAL :: temp_ad4
    REAL :: temp_ad3
    REAL :: temp_ad2
    REAL :: temp_ad1
    REAL(g_precision) :: temp_ad0
    REAL(g_precision) :: temp_ad
    IF (hord .LT. 0) THEN
! more dissipation
      ord_in = 2
      IF (hord .GE. 0.) THEN
        CALL PUSHCONTROL1B(1)
        ord = hord
      ELSE
        CALL PUSHCONTROL1B(1)
        ord = -hord
      END IF
    ELSE
      CALL PUSHCONTROL1B(0)
      ord_in = hord
      ord = hord
    END IF
    CALL PUSHREAL8ARRAY(q, (ied-isd+1)*(jed-jsd+1))
    CALL COPY_CORNERS(q, npx, npy, 2)
    CALL YTP(fy2, q, cry, ord_in, isd, ied, js, je, npx, npy, &
&      ppm_limiter)
    DO j=js,je+1
      DO i=isd,ied
        fyy(i, j) = yfx(i, j)*fy2(i, j)
      END DO
    END DO
    DO j=js,je
      DO i=isd,ied
        q_i(i, j) = (q(i, j)*area(i, j)+fyy(i, j)-fyy(i, j+1))/ra_y(i, j&
&         )
      END DO
    END DO
    CALL PUSHREAL8ARRAY(fx, (ie-is+2)*(je-js+1))
    CALL XTP(fx, q_i, crx(is, js), ord, is, ie, js, je, npx, npy, &
&      ppm_limiter)
    CALL PUSHREAL8ARRAY(q, (ied-isd+1)*(jed-jsd+1))
    CALL COPY_CORNERS(q, npx, npy, 1)
    CALL XTP(fx2, q, crx, ord_in, is, ie, jsd, jed, npx, npy, &
&      ppm_limiter)
    DO j=jsd,jed
      DO i=is,ie+1
        CALL PUSHREAL8(fx1(i))
        fx1(i) = xfx(i, j)*fx2(i, j)
      END DO
      DO i=is,ie
        q_j(i, j) = (q(i, j)*area(i, j)+fx1(i)-fx1(i+1))/ra_x(i, j)
      END DO
    END DO
    CALL PUSHREAL8ARRAY(fy, (ie-is+1)*(je-js+2))
    CALL YTP(fy, q_j, cry, ord, is, ie, js, je, npx, npy, ppm_limiter)
!----------------
! Flux averaging:
!----------------
    IF (PRESENT(mfx) .AND. PRESENT(mfy)) THEN
      IF (PRESENT(nord) .AND. PRESENT(damp_c)) THEN
        IF (damp_c .GT. 1.e-4) THEN
          damp = (damp_c*da_min)**(nord+1)
          CALL DELN_FLUX_ADM(nord, npx, npy, damp, q, q_ad, fx, fx_ad, &
&                      fy, fy_ad, mfx, mfy)
        END IF
      END IF
      fy2_ad = 0.0_8
      DO j=je+1,js,-1
        DO i=ie,is,-1
          temp_ad2 = 0.5*mfy(i, j)*fy_ad(i, j)
          fy2_ad(i, j) = fy2_ad(i, j) + temp_ad2
          mfy_ad(i, j) = mfy_ad(i, j) + 0.5*(fy(i, j)+fy2(i, j))*fy_ad(i&
&           , j)
          fy_ad(i, j) = temp_ad2
        END DO
      END DO
      fx2_ad = 0.0_8
      DO j=je,js,-1
        DO i=ie+1,is,-1
          temp_ad1 = 0.5*mfx(i, j)*fx_ad(i, j)
          fx2_ad(i, j) = fx2_ad(i, j) + temp_ad1
          mfx_ad(i, j) = mfx_ad(i, j) + 0.5*(fx(i, j)+fx2(i, j))*fx_ad(i&
&           , j)
          fx_ad(i, j) = temp_ad1
        END DO
      END DO
    ELSE
      IF (PRESENT(nord) .AND. PRESENT(damp_c)) THEN
        IF (damp_c .GT. 1.e-4) THEN
          damp = (damp_c*da_min)**(nord+1)
          CALL DELN_FLUX_ADM(nord, npx, npy, damp, q, q_ad, fx, fx_ad, &
&                      fy, fy_ad, xfx(is:ie+1, js:je), yfx(is:ie, js:je+&
&                      1))
        END IF
      END IF
      fy2_ad = 0.0_8
      DO j=je+1,js,-1
        DO i=ie,is,-1
          temp_ad4 = 0.5*yfx(i, j)*fy_ad(i, j)
          fy2_ad(i, j) = fy2_ad(i, j) + temp_ad4
          yfx_ad(i, j) = yfx_ad(i, j) + 0.5*(fy(i, j)+fy2(i, j))*fy_ad(i&
&           , j)
          fy_ad(i, j) = temp_ad4
        END DO
      END DO
      fx2_ad = 0.0_8
      DO j=je,js,-1
        DO i=ie+1,is,-1
          temp_ad3 = 0.5*xfx(i, j)*fx_ad(i, j)
          fx2_ad(i, j) = fx2_ad(i, j) + temp_ad3
          xfx_ad(i, j) = xfx_ad(i, j) + 0.5*(fx(i, j)+fx2(i, j))*fx_ad(i&
&           , j)
          fx_ad(i, j) = temp_ad3
        END DO
      END DO
    END IF
    CALL POPREAL8ARRAY(fy, (ie-is+1)*(je-js+2))
    q_j_ad = 0.0_8
    CALL YTP_ADM(fy, fy_ad, q_j, q_j_ad, cry, cry_ad, ord, is, ie, js, &
&          je, npx, npy, ppm_limiter)
    fx1_ad = 0.0_8
    DO j=jed,jsd,-1
      DO i=ie,is,-1
        temp_ad0 = q_j_ad(i, j)/ra_x(i, j)
        q_ad(i, j) = q_ad(i, j) + area(i, j)*temp_ad0
        fx1_ad(i) = fx1_ad(i) + temp_ad0
        fx1_ad(i+1) = fx1_ad(i+1) - temp_ad0
        ra_x_ad(i, j) = ra_x_ad(i, j) - (area(i, j)*q(i, j)+fx1(i)-fx1(i&
&         +1))*temp_ad0/ra_x(i, j)
        q_j_ad(i, j) = 0.0_8
      END DO
      DO i=ie+1,is,-1
        CALL POPREAL8(fx1(i))
        xfx_ad(i, j) = xfx_ad(i, j) + fx2(i, j)*fx1_ad(i)
        fx2_ad(i, j) = fx2_ad(i, j) + xfx(i, j)*fx1_ad(i)
        fx1_ad(i) = 0.0_8
      END DO
    END DO
    CALL XTP_ADM(fx2, fx2_ad, q, q_ad, crx, crx_ad, ord_in, is, ie, jsd&
&          , jed, npx, npy, ppm_limiter)
    CALL POPREAL8ARRAY(q, (ied-isd+1)*(jed-jsd+1))
    CALL COPY_CORNERS_ADM(q, q_ad, npx, npy, 1)
    CALL POPREAL8ARRAY(fx, (ie-is+2)*(je-js+1))
    q_i_ad = 0.0_8
    CALL XTP_ADM(fx, fx_ad, q_i, q_i_ad, crx(is, js), crx_ad(is, js), &
&          ord, is, ie, js, je, npx, npy, ppm_limiter)
    fyy_ad = 0.0_8
    DO j=je,js,-1
      DO i=ied,isd,-1
        temp_ad = q_i_ad(i, j)/ra_y(i, j)
        q_ad(i, j) = q_ad(i, j) + area(i, j)*temp_ad
        fyy_ad(i, j) = fyy_ad(i, j) + temp_ad
        fyy_ad(i, j+1) = fyy_ad(i, j+1) - temp_ad
        ra_y_ad(i, j) = ra_y_ad(i, j) - (area(i, j)*q(i, j)+fyy(i, j)-&
&         fyy(i, j+1))*temp_ad/ra_y(i, j)
        q_i_ad(i, j) = 0.0_8
      END DO
    END DO
    DO j=je+1,js,-1
      DO i=ied,isd,-1
        yfx_ad(i, j) = yfx_ad(i, j) + fy2(i, j)*fyy_ad(i, j)
        fy2_ad(i, j) = fy2_ad(i, j) + yfx(i, j)*fyy_ad(i, j)
        fyy_ad(i, j) = 0.0_8
      END DO
    END DO
    CALL YTP_ADM(fy2, fy2_ad, q, q_ad, cry, cry_ad, ord_in, isd, ied, js&
&          , je, npx, npy, ppm_limiter)
    CALL POPREAL8ARRAY(q, (ied-isd+1)*(jed-jsd+1))
    CALL COPY_CORNERS_ADM(q, q_ad, npx, npy, 2)
    CALL POPCONTROL1B(branch)
  END SUBROUTINE FV_TP_2D_COMPUTE_ADM
  SUBROUTINE FV_TP_2D_COMPUTE(q, crx, cry, npx, npy, hord, fx, fy, xfx, &
&   yfx, area, ra_x, ra_y, mfx, mfy, ppm_fac, nord, damp_c)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx, npy
    INTEGER, INTENT(IN) :: hord
!
    REAL, INTENT(IN) :: crx(is:ie+1, jsd:jed)
!
    REAL, INTENT(IN) :: xfx(is:ie+1, jsd:jed)
!
    REAL, INTENT(IN) :: cry(isd:ied, js:je+1)
!
    REAL, INTENT(IN) :: yfx(isd:ied, js:je+1)
    REAL(g_precision), INTENT(IN) :: area(isd:ied, jsd:jed)
    REAL, INTENT(IN) :: ra_x(is:ie, jsd:jed)
    REAL, INTENT(IN) :: ra_y(isd:ied, js:je)
! transported scalar
    REAL, INTENT(INOUT) :: q(isd:ied, jsd:jed)
! Flux in x ( E )
    REAL, INTENT(OUT) :: fx(is:ie+1, js:je)
! Flux in y ( N )
    REAL, INTENT(OUT) :: fy(is:ie, js:je+1)
! optional Arguments:
! Mass Flux X-Dir
    REAL, OPTIONAL, INTENT(IN) :: mfx(is:ie+1, js:je)
! Mass Flux Y-Dir
    REAL, OPTIONAL, INTENT(IN) :: mfy(is:ie, js:je+1)
! for ord=4 option
    REAL, OPTIONAL, INTENT(IN) :: ppm_fac
    REAL, OPTIONAL, INTENT(IN) :: damp_c
    INTEGER, OPTIONAL, INTENT(IN) :: nord
! Local:
    INTEGER :: ord, ord_in
    REAL :: q_i(isd:ied, js:je)
    REAL :: q_j(is:ie, jsd:jed)
    REAL :: fx2(is:ie+1, jsd:jed)
    REAL :: fy2(isd:ied, js:je+1)
    REAL :: fyy(isd:ied, js:je+1)
    REAL :: fx1(is:ie+1)
    REAL :: ppm_limiter, damp
    INTEGER :: i, j
    INTRINSIC ABS
    INTRINSIC PRESENT
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
    CALL COPY_CORNERS(q, npx, npy, 2)
    CALL YTP(fy2, q, cry, ord_in, isd, ied, js, je, npx, npy, &
&      ppm_limiter)
    DO j=js,je+1
      DO i=isd,ied
        fyy(i, j) = yfx(i, j)*fy2(i, j)
      END DO
    END DO
    DO j=js,je
      DO i=isd,ied
        q_i(i, j) = (q(i, j)*area(i, j)+fyy(i, j)-fyy(i, j+1))/ra_y(i, j&
&         )
      END DO
    END DO
    CALL XTP(fx, q_i, crx(is, js), ord, is, ie, js, je, npx, npy, &
&      ppm_limiter)
    CALL COPY_CORNERS(q, npx, npy, 1)
    CALL XTP(fx2, q, crx, ord_in, is, ie, jsd, jed, npx, npy, &
&      ppm_limiter)
    DO j=jsd,jed
      DO i=is,ie+1
        fx1(i) = xfx(i, j)*fx2(i, j)
      END DO
      DO i=is,ie
        q_j(i, j) = (q(i, j)*area(i, j)+fx1(i)-fx1(i+1))/ra_x(i, j)
      END DO
    END DO
    CALL YTP(fy, q_j, cry, ord, is, ie, js, je, npx, npy, ppm_limiter)
!----------------
! Flux averaging:
!----------------
    IF (PRESENT(mfx) .AND. PRESENT(mfy)) THEN
!---------------------------------
! For transport of pt and tracers
!---------------------------------
      DO j=js,je
        DO i=is,ie+1
          fx(i, j) = 0.5*(fx(i, j)+fx2(i, j))*mfx(i, j)
        END DO
      END DO
      DO j=js,je+1
        DO i=is,ie
          fy(i, j) = 0.5*(fy(i, j)+fy2(i, j))*mfy(i, j)
        END DO
      END DO
      IF (PRESENT(nord) .AND. PRESENT(damp_c)) THEN
        IF (damp_c .GT. 1.e-4) THEN
          damp = (damp_c*da_min)**(nord+1)
          CALL DELN_FLUX(nord, npx, npy, damp, q, fx, fy, mfx, mfy)
        END IF
      END IF
    ELSE
!---------------------------------
! For transport of delp, vorticity
!---------------------------------
      DO j=js,je
        DO i=is,ie+1
          fx(i, j) = 0.5*(fx(i, j)+fx2(i, j))*xfx(i, j)
        END DO
      END DO
      DO j=js,je+1
        DO i=is,ie
          fy(i, j) = 0.5*(fy(i, j)+fy2(i, j))*yfx(i, j)
        END DO
      END DO
      IF (PRESENT(nord) .AND. PRESENT(damp_c)) THEN
        IF (damp_c .GT. 1.e-4) THEN
          damp = (damp_c*da_min)**(nord+1)
          CALL DELN_FLUX(nord, npx, npy, damp, q, fx, fy, xfx(is:ie+1, &
&                  js:je), yfx(is:ie, js:je+1))
        END IF
      END IF
    END IF
  END SUBROUTINE FV_TP_2D_COMPUTE
  SUBROUTINE COPY_CORNERS_R4(q, npx, npy, dir)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx, npy, dir
    REAL*4, INTENT(INOUT) :: q(isd:ied, jsd:jed)
    INTEGER :: i, j
    IF (dir .EQ. 1) THEN
! XDir:
      IF (sw_corner) THEN
        DO j=1-ng,0
          DO i=1-ng,0
            q(i, j) = q(j, 1-i)
          END DO
        END DO
      END IF
      IF (se_corner) THEN
        DO j=1-ng,0
          DO i=npx,npx+ng-1
            q(i, j) = q(npy-j, i-npx+1)
          END DO
        END DO
      END IF
      IF (ne_corner) THEN
        DO j=npy,npy+ng-1
          DO i=npx,npx+ng-1
            q(i, j) = q(j, 2*npx-1-i)
          END DO
        END DO
      END IF
      IF (nw_corner) THEN
        DO j=npy,npy+ng-1
          DO i=1-ng,0
            q(i, j) = q(npy-j, i-1+npx)
          END DO
        END DO
      END IF
    ELSE IF (dir .EQ. 2) THEN
! YDir:
      IF (sw_corner) THEN
        DO j=1-ng,0
          DO i=1-ng,0
            q(i, j) = q(1-j, i)
          END DO
        END DO
      END IF
      IF (se_corner) THEN
        DO j=1-ng,0
          DO i=npx,npx+ng-1
            q(i, j) = q(npy+j-1, npx-i)
          END DO
        END DO
      END IF
      IF (ne_corner) THEN
        DO j=npy,npy+ng-1
          DO i=npx,npx+ng-1
            q(i, j) = q(2*npy-1-j, i)
          END DO
        END DO
      END IF
      IF (nw_corner) THEN
        DO j=npy,npy+ng-1
          DO i=1-ng,0
            q(i, j) = q(j+1-npx, npy-i)
          END DO
        END DO
      END IF
    END IF
  END SUBROUTINE COPY_CORNERS_R4
!  Differentiation of xtp in reverse (adjoint) mode (with options r8):
!   gradient     of useful results: q fx c
!   with respect to varying inputs: q fx c
  SUBROUTINE XTP_ADM(fx, fx_ad, q, q_ad, c, c_ad, iord, ifirst, ilast, &
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
    REAL :: c_ad(is:ie+1, jfirst:jlast)
    REAL, INTENT(IN) :: q(isd:ied, jfirst:jlast)
    REAL :: q_ad(isd:ied, jfirst:jlast)
    REAL, INTENT(IN) :: ppm_limiter
    REAL :: fx(ifirst:ilast+1, jfirst:jlast)
    REAL :: fx_ad(ifirst:ilast+1, jfirst:jlast)
! Local:
    REAL :: x0, x1
    INTEGER :: i, j
    INTEGER :: branch
    REAL :: temp_ad4
    REAL :: temp_ad3
    REAL :: temp_ad2
    REAL :: temp_ad1
    REAL :: temp_ad0
    REAL :: temp_ad
    IF (iord .EQ. 1) THEN
      DO j=jfirst,jlast
        DO i=ifirst,ilast+1
          IF (c(i, j) .GT. 0.) THEN
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
      END DO
      DO j=jlast,jfirst,-1
        DO i=ilast+1,ifirst,-1
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            q_ad(i, j) = q_ad(i, j) + fx_ad(i, j)
            fx_ad(i, j) = 0.0_8
          ELSE
            q_ad(i-1, j) = q_ad(i-1, j) + fx_ad(i, j)
            fx_ad(i, j) = 0.0_8
          END IF
        END DO
      END DO
    ELSE IF (iord .EQ. 333) THEN
      DO j=jfirst,jlast
        DO i=ifirst,ilast+1
!At edges keep to first order
          IF (i .EQ. ifirst .OR. i .EQ. ilast + 1) THEN
            IF (c(i, j) .GT. 0.) THEN
              CALL PUSHCONTROL2B(3)
            ELSE
              CALL PUSHCONTROL2B(2)
            END IF
          ELSE IF (c(i, j) .GT. 0.) THEN
!Otherwise use the third order scheme
            CALL PUSHCONTROL2B(1)
          ELSE
            CALL PUSHCONTROL2B(0)
          END IF
        END DO
      END DO
      DO j=jlast,jfirst,-1
        DO i=ilast+1,ifirst,-1
          CALL POPCONTROL2B(branch)
          IF (branch .LT. 2) THEN
            IF (branch .EQ. 0) THEN
              temp_ad2 = fx_ad(i, j)/6.0
              temp_ad3 = -(0.5*c(i, j)*fx_ad(i, j))
              temp_ad4 = c(i, j)**2*fx_ad(i, j)/6.0
              q_ad(i-1, j) = q_ad(i-1, j) + temp_ad4 - temp_ad3 + 2.0*&
&               temp_ad2
              q_ad(i, j) = q_ad(i, j) + temp_ad3 - 2.0*temp_ad4 + 5.0*&
&               temp_ad2
              q_ad(i+1, j) = q_ad(i+1, j) + temp_ad4 - temp_ad2
              c_ad(i, j) = c_ad(i, j) + ((q(i+1, j)-2.0*q(i, j)+q(i-1, j&
&               ))*2*c(i, j)/6.0-0.5*(q(i, j)-q(i-1, j)))*fx_ad(i, j)
              fx_ad(i, j) = 0.0_8
            ELSE
              temp_ad = fx_ad(i, j)/6.0
              temp_ad0 = -(0.5*c(i, j)*fx_ad(i, j))
              temp_ad1 = c(i, j)**2*fx_ad(i, j)/6.0
              q_ad(i, j) = q_ad(i, j) + temp_ad1 + temp_ad0 + 2.0*&
&               temp_ad
              q_ad(i-1, j) = q_ad(i-1, j) + 5.0*temp_ad - temp_ad0 - 2.0&
&               *temp_ad1
              q_ad(i-2, j) = q_ad(i-2, j) + temp_ad1 - temp_ad
              c_ad(i, j) = c_ad(i, j) + ((q(i, j)-2.0*q(i-1, j)+q(i-2, j&
&               ))*2*c(i, j)/6.0-0.5*(q(i, j)-q(i-1, j)))*fx_ad(i, j)
              fx_ad(i, j) = 0.0_8
            END IF
          ELSE IF (branch .EQ. 2) THEN
            q_ad(i, j) = q_ad(i, j) + fx_ad(i, j)
            fx_ad(i, j) = 0.0_8
          ELSE
            q_ad(i-1, j) = q_ad(i-1, j) + fx_ad(i, j)
            fx_ad(i, j) = 0.0_8
          END IF
        END DO
      END DO
    ELSE
      CALL FXPPM_ADM(c, c_ad, q, q_ad, fx, fx_ad, iord, ifirst, ilast, &
&              jfirst, jlast, npx, npy, ppm_limiter)
    END IF
  END SUBROUTINE XTP_ADM
  SUBROUTINE XTP(fx, q, c, iord, ifirst, ilast, jfirst, jlast, npx, npy&
&   , ppm_limiter)
    IMPLICIT NONE
!  X-Dir strip
    INTEGER, INTENT(IN) :: ifirst, ilast
!  Y-Dir strip
    INTEGER, INTENT(IN) :: jfirst, jlast
    INTEGER, INTENT(IN) :: npx, npy
    INTEGER, INTENT(IN) :: iord
! Courant numbers
    REAL, INTENT(IN) :: c(is:ie+1, jfirst:jlast)
    REAL, INTENT(IN) :: q(isd:ied, jfirst:jlast)
    REAL, INTENT(IN) :: ppm_limiter
    REAL, INTENT(OUT) :: fx(ifirst:ilast+1, jfirst:jlast)
! Local:
    REAL :: x0, x1
    INTEGER :: i, j
    IF (iord .EQ. 1) THEN
      DO j=jfirst,jlast
        DO i=ifirst,ilast+1
          IF (c(i, j) .GT. 0.) THEN
            fx(i, j) = q(i-1, j)
          ELSE
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
              fx(i, j) = q(i-1, j)
            ELSE
              fx(i, j) = q(i, j)
            END IF
          ELSE IF (c(i, j) .GT. 0.) THEN
!Otherwise use the third order scheme
            fx(i, j) = (2.0*q(i, j)+5.0*q(i-1, j)-q(i-2, j))/6.0 - 0.5*c&
&             (i, j)*(q(i, j)-q(i-1, j)) + c(i, j)*c(i, j)/6.0*(q(i, j)-&
&             2.0*q(i-1, j)+q(i-2, j))
          ELSE
            fx(i, j) = (2.0*q(i-1, j)+5.0*q(i, j)-q(i+1, j))/6.0 - 0.5*c&
&             (i, j)*(q(i, j)-q(i-1, j)) + c(i, j)*c(i, j)/6.0*(q(i+1, j&
&             )-2.0*q(i, j)+q(i-1, j))
          END IF
        END DO
      END DO
    ELSE
      CALL FXPPM(c, q, fx, iord, ifirst, ilast, jfirst, jlast, npx, npy&
&          , ppm_limiter)
    END IF
  END SUBROUTINE XTP
!  Differentiation of ytp in reverse (adjoint) mode (with options r8):
!   gradient     of useful results: q fy c
!   with respect to varying inputs: q fy c
  SUBROUTINE YTP_ADM(fy, fy_ad, q, q_ad, c, c_ad, jord, ifirst, ilast, &
&   jfirst, jlast, npx, npy, ppm_limiter)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx, npy
!  X-Dir strip
    INTEGER, INTENT(IN) :: ifirst, ilast
!  Y-Dir strip
    INTEGER, INTENT(IN) :: jfirst, jlast
    INTEGER, INTENT(IN) :: jord
    REAL, INTENT(IN) :: q(ifirst:ilast, jfirst-ng:jlast+ng)
    REAL :: q_ad(ifirst:ilast, jfirst-ng:jlast+ng)
! Courant number
    REAL, INTENT(IN) :: c(isd:ied, js:je+1)
    REAL :: c_ad(isd:ied, js:je+1)
!  Flux
    REAL :: fy(ifirst:ilast, jfirst:jlast+1)
    REAL :: fy_ad(ifirst:ilast, jfirst:jlast+1)
    REAL, INTENT(IN) :: ppm_limiter
! !LOCAL VARIABLES:
    REAL :: x0, x1
    INTEGER :: i, j
    INTEGER :: branch
    REAL :: temp_ad4
    REAL :: temp_ad3
    REAL :: temp_ad2
    REAL :: temp_ad1
    REAL :: temp_ad0
    REAL :: temp_ad
    IF (jord .EQ. 1) THEN
      DO j=jfirst,jlast+1
        DO i=ifirst,ilast
          IF (c(i, j) .GT. 0.) THEN
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
      END DO
      DO j=jlast+1,jfirst,-1
        DO i=ilast,ifirst,-1
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            q_ad(i, j) = q_ad(i, j) + fy_ad(i, j)
            fy_ad(i, j) = 0.0_8
          ELSE
            q_ad(i, j-1) = q_ad(i, j-1) + fy_ad(i, j)
            fy_ad(i, j) = 0.0_8
          END IF
        END DO
      END DO
    ELSE IF (jord .EQ. 333) THEN
      DO j=jfirst,jlast+1
        DO i=ifirst,ilast
!At edges keep to first order
          IF (j .EQ. jfirst .OR. j .EQ. jlast + 1) THEN
            IF (c(i, j) .GT. 0.) THEN
              CALL PUSHCONTROL2B(3)
            ELSE
              CALL PUSHCONTROL2B(2)
            END IF
          ELSE IF (c(i, j) .GT. 0.) THEN
!Otherwise use the third order scheme
            CALL PUSHCONTROL2B(1)
          ELSE
            CALL PUSHCONTROL2B(0)
          END IF
        END DO
      END DO
      DO j=jlast+1,jfirst,-1
        DO i=ilast,ifirst,-1
          CALL POPCONTROL2B(branch)
          IF (branch .LT. 2) THEN
            IF (branch .EQ. 0) THEN
              temp_ad2 = fy_ad(i, j)/6.0
              temp_ad3 = -(0.5*c(i, j)*fy_ad(i, j))
              temp_ad4 = c(i, j)**2*fy_ad(i, j)/6.0
              q_ad(i, j-1) = q_ad(i, j-1) + temp_ad4 - temp_ad3 + 2.0*&
&               temp_ad2
              q_ad(i, j) = q_ad(i, j) + temp_ad3 - 2.0*temp_ad4 + 5.0*&
&               temp_ad2
              q_ad(i, j+1) = q_ad(i, j+1) + temp_ad4 - temp_ad2
              c_ad(i, j) = c_ad(i, j) + ((q(i, j+1)-2.0*q(i, j)+q(i, j-1&
&               ))*2*c(i, j)/6.0-0.5*(q(i, j)-q(i, j-1)))*fy_ad(i, j)
              fy_ad(i, j) = 0.0_8
            ELSE
              temp_ad = fy_ad(i, j)/6.0
              temp_ad0 = -(0.5*c(i, j)*fy_ad(i, j))
              temp_ad1 = c(i, j)**2*fy_ad(i, j)/6.0
              q_ad(i, j) = q_ad(i, j) + temp_ad1 + temp_ad0 + 2.0*&
&               temp_ad
              q_ad(i, j-1) = q_ad(i, j-1) + 5.0*temp_ad - temp_ad0 - 2.0&
&               *temp_ad1
              q_ad(i, j-2) = q_ad(i, j-2) + temp_ad1 - temp_ad
              c_ad(i, j) = c_ad(i, j) + ((q(i, j)-2.0*q(i, j-1)+q(i, j-2&
&               ))*2*c(i, j)/6.0-0.5*(q(i, j)-q(i, j-1)))*fy_ad(i, j)
              fy_ad(i, j) = 0.0_8
            END IF
          ELSE IF (branch .EQ. 2) THEN
            q_ad(i, j) = q_ad(i, j) + fy_ad(i, j)
            fy_ad(i, j) = 0.0_8
          ELSE
            q_ad(i, j-1) = q_ad(i, j-1) + fy_ad(i, j)
            fy_ad(i, j) = 0.0_8
          END IF
        END DO
      END DO
    ELSE
      CALL FYPPM_ADM(c, c_ad, q, q_ad, fy, fy_ad, jord, ifirst, ilast, &
&              jfirst, jlast, npx, npy, ppm_limiter)
    END IF
  END SUBROUTINE YTP_ADM
  SUBROUTINE YTP(fy, q, c, jord, ifirst, ilast, jfirst, jlast, npx, npy&
&   , ppm_limiter)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx, npy
!  X-Dir strip
    INTEGER, INTENT(IN) :: ifirst, ilast
!  Y-Dir strip
    INTEGER, INTENT(IN) :: jfirst, jlast
    INTEGER, INTENT(IN) :: jord
    REAL, INTENT(IN) :: q(ifirst:ilast, jfirst-ng:jlast+ng)
! Courant number
    REAL, INTENT(IN) :: c(isd:ied, js:je+1)
!  Flux
    REAL, INTENT(OUT) :: fy(ifirst:ilast, jfirst:jlast+1)
    REAL, INTENT(IN) :: ppm_limiter
! !LOCAL VARIABLES:
    REAL :: x0, x1
    INTEGER :: i, j
    IF (jord .EQ. 1) THEN
      DO j=jfirst,jlast+1
        DO i=ifirst,ilast
          IF (c(i, j) .GT. 0.) THEN
            fy(i, j) = q(i, j-1)
          ELSE
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
              fy(i, j) = q(i, j-1)
            ELSE
              fy(i, j) = q(i, j)
            END IF
          ELSE IF (c(i, j) .GT. 0.) THEN
!Otherwise use the third order scheme
            fy(i, j) = (2.0*q(i, j)+5.0*q(i, j-1)-q(i, j-2))/6.0 - 0.5*c&
&             (i, j)*(q(i, j)-q(i, j-1)) + c(i, j)*c(i, j)/6.0*(q(i, j)-&
&             2.0*q(i, j-1)+q(i, j-2))
          ELSE
            fy(i, j) = (2.0*q(i, j-1)+5.0*q(i, j)-q(i, j+1))/6.0 - 0.5*c&
&             (i, j)*(q(i, j)-q(i, j-1)) + c(i, j)*c(i, j)/6.0*(q(i, j+1&
&             )-2.0*q(i, j)+q(i, j-1))
          END IF
        END DO
      END DO
    ELSE
      CALL FYPPM(c, q, fy, jord, ifirst, ilast, jfirst, jlast, npx, npy&
&          , ppm_limiter)
    END IF
  END SUBROUTINE YTP
!  Differentiation of fxppm in reverse (adjoint) mode (with options r8):
!   gradient     of useful results: q flux c
!   with respect to varying inputs: q flux c
  SUBROUTINE FXPPM_ADM(c, c_ad, q, q_ad, flux, flux_ad, iord, ifirst, &
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
    REAL :: q_ad(ifirst-ng:ilast+ng, jfirst:jlast)
! Courant   N (like FLUX)
    REAL, INTENT(IN) :: c(ifirst:ilast+1, jfirst:jlast)
    REAL :: c_ad(ifirst:ilast+1, jfirst:jlast)
    REAL, INTENT(IN) :: ppm_limiter
! !OUTPUT PARAMETERS:
!  Flux
    REAL :: flux(ifirst:ilast+1, jfirst:jlast)
    REAL :: flux_ad(ifirst:ilast+1, jfirst:jlast)
! Local
    REAL :: dm1(ifirst-2:ilast+2)
    REAL :: al(ifirst-1:ilast+2)
    REAL :: bl(ifirst-1:ilast+1)
    REAL :: bl_ad(ifirst-1:ilast+1)
    REAL :: br(ifirst-1:ilast+1)
    REAL :: br_ad(ifirst-1:ilast+1)
    REAL :: dq(ifirst-3:ilast+2)
    REAL :: dl, dr, pmp, lac, ct, qe
    REAL :: xt, x1, x0
    REAL :: xt_ad
    INTEGER :: i, j, is3, ie3, ie2, it
    INTRINSIC MAX
    INTRINSIC MIN
    INTEGER :: branch
    REAL :: y2_ad
    REAL :: y3_ad
    REAL :: y4_ad
    REAL :: y5_ad
    REAL :: y6_ad
    REAL :: y7_ad
    REAL :: temp_ad8
    REAL :: temp_ad7
    REAL :: temp_ad6
    REAL :: temp_ad5
    REAL :: temp_ad4
    REAL :: temp_ad3
    REAL :: temp_ad2
    REAL :: temp_ad1
    REAL :: y8_ad
    REAL :: temp_ad0
    REAL :: temp_ad
    REAL :: temp
    REAL :: y8
    REAL :: y7
    REAL :: y6
    REAL :: y5
    REAL :: y4
    REAL :: y3
    REAL :: y1_ad
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
    IF (iord .EQ. 6) THEN
      DO j=jfirst,jlast
! Non-monotonic "5th order" scheme (not really 5th order)
        DO i=is3,ie3
          CALL PUSHREAL8(bl(i))
          bl(i) = b5*q(i-2, j) + b4*q(i-1, j) + b3*q(i, j) + b2*q(i+1, j&
&           ) + b1*q(i+2, j)
          CALL PUSHREAL8(br(i))
          br(i) = b1*q(i-2, j) + b2*q(i-1, j) + b3*q(i, j) + b4*q(i+1, j&
&           ) + b5*q(i+2, j)
        END DO
!--------------
! fix the edges
!--------------
        IF (is .EQ. 1) THEN
          CALL PUSHREAL8(br(2))
          br(2) = p1*(q(2, j)+q(3, j)) + p2*(q(1, j)+q(4, j)) - q(2, j)
          xt = 0.5*((2.*dxa(1, j)+dxa(2, j))*(q(0, j)+q(1, j))-dxa(1, j)&
&           *(q(-1, j)+q(2, j)))/(dxa(1, j)+dxa(2, j))
          CALL PUSHREAL8(bl(1))
          bl(1) = xt - q(1, j)
          CALL PUSHREAL8(br(0))
          br(0) = xt - q(0, j)
          xt = c1*q(-2, j) + c2*q(-1, j) + c3*q(0, j)
          IF (q(-1, j) .GT. q(0, j)) THEN
            y1 = q(0, j)
            CALL PUSHCONTROL1B(0)
          ELSE
            y1 = q(-1, j)
            CALL PUSHCONTROL1B(1)
          END IF
          IF (xt .LT. y1) THEN
            xt = y1
            CALL PUSHCONTROL1B(0)
          ELSE
            xt = xt
            CALL PUSHCONTROL1B(1)
          END IF
          IF (q(-1, j) .LT. q(0, j)) THEN
            y2 = q(0, j)
            CALL PUSHCONTROL1B(0)
          ELSE
            y2 = q(-1, j)
            CALL PUSHCONTROL1B(1)
          END IF
          IF (xt .GT. y2) THEN
            xt = y2
            CALL PUSHCONTROL1B(0)
          ELSE
            xt = xt
            CALL PUSHCONTROL1B(1)
          END IF
          CALL PUSHREAL8(bl(0))
          bl(0) = xt - q(0, j)
          xt = c3*q(1, j) + c2*q(2, j) + c1*q(3, j)
          IF (q(1, j) .GT. q(2, j)) THEN
            y3 = q(2, j)
            CALL PUSHCONTROL1B(0)
          ELSE
            y3 = q(1, j)
            CALL PUSHCONTROL1B(1)
          END IF
          IF (xt .LT. y3) THEN
            xt = y3
            CALL PUSHCONTROL1B(0)
          ELSE
            xt = xt
            CALL PUSHCONTROL1B(1)
          END IF
          IF (q(1, j) .LT. q(2, j)) THEN
            y4 = q(2, j)
            CALL PUSHCONTROL1B(0)
          ELSE
            y4 = q(1, j)
            CALL PUSHCONTROL1B(1)
          END IF
          IF (xt .GT. y4) THEN
            xt = y4
            CALL PUSHCONTROL1B(0)
          ELSE
            xt = xt
            CALL PUSHCONTROL1B(1)
          END IF
          CALL PUSHREAL8(br(1))
          br(1) = xt - q(1, j)
          CALL PUSHREAL8(bl(2))
          bl(2) = xt - q(2, j)
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
        IF (ie + 1 .EQ. npx) THEN
          CALL PUSHREAL8(bl(npx-2))
          bl(npx-2) = p1*(q(npx-2, j)+q(npx-3, j)) + p2*(q(npx-4, j)+q(&
&           npx-1, j)) - q(npx-2, j)
          xt = 0.5*((2.*dxa(npx-1, j)+dxa(npx-2, j))*(q(npx-1, j)+q(npx&
&           , j))-dxa(npx-1, j)*(q(npx-2, j)+q(npx+1, j)))/(dxa(npx-1, j&
&           )+dxa(npx-2, j))
          CALL PUSHREAL8(br(npx-1))
          br(npx-1) = xt - q(npx-1, j)
          CALL PUSHREAL8(bl(npx))
          bl(npx) = xt - q(npx, j)
          xt = c3*q(npx, j) + c2*q(npx+1, j) + c1*q(npx+2, j)
          IF (q(npx, j) .GT. q(npx+1, j)) THEN
            y5 = q(npx+1, j)
            CALL PUSHCONTROL1B(0)
          ELSE
            y5 = q(npx, j)
            CALL PUSHCONTROL1B(1)
          END IF
          IF (xt .LT. y5) THEN
            xt = y5
            CALL PUSHCONTROL1B(0)
          ELSE
            xt = xt
            CALL PUSHCONTROL1B(1)
          END IF
          IF (q(npx, j) .LT. q(npx+1, j)) THEN
            y6 = q(npx+1, j)
            CALL PUSHCONTROL1B(0)
          ELSE
            y6 = q(npx, j)
            CALL PUSHCONTROL1B(1)
          END IF
          IF (xt .GT. y6) THEN
            xt = y6
            CALL PUSHCONTROL1B(0)
          ELSE
            xt = xt
            CALL PUSHCONTROL1B(1)
          END IF
          CALL PUSHREAL8(br(npx))
          br(npx) = xt - q(npx, j)
          xt = c1*q(npx-3, j) + c2*q(npx-2, j) + c3*q(npx-1, j)
          IF (q(npx-2, j) .GT. q(npx-1, j)) THEN
            y7 = q(npx-1, j)
            CALL PUSHCONTROL1B(0)
          ELSE
            y7 = q(npx-2, j)
            CALL PUSHCONTROL1B(1)
          END IF
          IF (xt .LT. y7) THEN
            xt = y7
            CALL PUSHCONTROL1B(0)
          ELSE
            xt = xt
            CALL PUSHCONTROL1B(1)
          END IF
          IF (q(npx-2, j) .LT. q(npx-1, j)) THEN
            y8 = q(npx-1, j)
            CALL PUSHCONTROL1B(0)
          ELSE
            y8 = q(npx-2, j)
            CALL PUSHCONTROL1B(1)
          END IF
          IF (xt .GT. y8) THEN
            xt = y8
            CALL PUSHCONTROL1B(0)
          ELSE
            xt = xt
            CALL PUSHCONTROL1B(1)
          END IF
          CALL PUSHREAL8(br(npx-2))
          br(npx-2) = xt - q(npx-2, j)
          CALL PUSHREAL8(bl(npx-1))
          bl(npx-1) = xt - q(npx-1, j)
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
        DO i=ifirst,ilast+1
          IF (c(i, j) .GT. 0.) THEN
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
      END DO
      bl_ad = 0.0_8
      br_ad = 0.0_8
      DO j=jlast,jfirst,-1
        DO i=ilast+1,ifirst,-1
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            temp_ad7 = (c(i, j)+1.)*flux_ad(i, j)
            temp_ad8 = c(i, j)*temp_ad7
            q_ad(i, j) = q_ad(i, j) + flux_ad(i, j)
            c_ad(i, j) = c_ad(i, j) + (bl(i)+br(i))*temp_ad7 + (bl(i)+c(&
&             i, j)*(bl(i)+br(i)))*flux_ad(i, j)
            bl_ad(i) = bl_ad(i) + temp_ad8 + temp_ad7
            br_ad(i) = br_ad(i) + temp_ad8
            flux_ad(i, j) = 0.0_8
          ELSE
            temp = bl(i-1) + br(i-1)
            temp_ad5 = (1.-c(i, j))*flux_ad(i, j)
            temp_ad6 = -(c(i, j)*temp_ad5)
            q_ad(i-1, j) = q_ad(i-1, j) + flux_ad(i, j)
            c_ad(i, j) = c_ad(i, j) - temp*temp_ad5 - (br(i-1)-c(i, j)*&
&             temp)*flux_ad(i, j)
            br_ad(i-1) = br_ad(i-1) + temp_ad6 + temp_ad5
            bl_ad(i-1) = bl_ad(i-1) + temp_ad6
            flux_ad(i, j) = 0.0_8
          END IF
        END DO
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) THEN
          CALL POPREAL8(bl(npx-1))
          xt_ad = br_ad(npx-2) + bl_ad(npx-1)
          q_ad(npx-1, j) = q_ad(npx-1, j) - bl_ad(npx-1)
          bl_ad(npx-1) = 0.0_8
          CALL POPREAL8(br(npx-2))
          q_ad(npx-2, j) = q_ad(npx-2, j) - br_ad(npx-2)
          br_ad(npx-2) = 0.0_8
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            y8_ad = xt_ad
            xt_ad = 0.0_8
          ELSE
            y8_ad = 0.0_8
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            q_ad(npx-1, j) = q_ad(npx-1, j) + y8_ad
          ELSE
            q_ad(npx-2, j) = q_ad(npx-2, j) + y8_ad
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            y7_ad = xt_ad
            xt_ad = 0.0_8
          ELSE
            y7_ad = 0.0_8
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            q_ad(npx-1, j) = q_ad(npx-1, j) + y7_ad
          ELSE
            q_ad(npx-2, j) = q_ad(npx-2, j) + y7_ad
          END IF
          q_ad(npx-3, j) = q_ad(npx-3, j) + c1*xt_ad
          q_ad(npx-2, j) = q_ad(npx-2, j) + c2*xt_ad
          q_ad(npx-1, j) = q_ad(npx-1, j) + c3*xt_ad
          CALL POPREAL8(br(npx))
          xt_ad = br_ad(npx)
          q_ad(npx, j) = q_ad(npx, j) - br_ad(npx)
          br_ad(npx) = 0.0_8
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            y6_ad = xt_ad
            xt_ad = 0.0_8
          ELSE
            y6_ad = 0.0_8
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            q_ad(npx+1, j) = q_ad(npx+1, j) + y6_ad
          ELSE
            q_ad(npx, j) = q_ad(npx, j) + y6_ad
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            y5_ad = xt_ad
            xt_ad = 0.0_8
          ELSE
            y5_ad = 0.0_8
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            q_ad(npx+1, j) = q_ad(npx+1, j) + y5_ad
          ELSE
            q_ad(npx, j) = q_ad(npx, j) + y5_ad
          END IF
          q_ad(npx, j) = q_ad(npx, j) + c3*xt_ad
          q_ad(npx+1, j) = q_ad(npx+1, j) + c2*xt_ad
          q_ad(npx+2, j) = q_ad(npx+2, j) + c1*xt_ad
          CALL POPREAL8(bl(npx))
          xt_ad = br_ad(npx-1) + bl_ad(npx)
          q_ad(npx, j) = q_ad(npx, j) - bl_ad(npx)
          bl_ad(npx) = 0.0_8
          CALL POPREAL8(br(npx-1))
          temp_ad3 = 0.5*xt_ad/(dxa(npx-1, j)+dxa(npx-2, j))
          temp_ad2 = (dxa(npx-1, j)*2.+dxa(npx-2, j))*temp_ad3
          q_ad(npx-1, j) = q_ad(npx-1, j) + temp_ad2 - br_ad(npx-1)
          br_ad(npx-1) = 0.0_8
          temp_ad4 = -(dxa(npx-1, j)*temp_ad3)
          q_ad(npx, j) = q_ad(npx, j) + temp_ad2
          q_ad(npx-2, j) = q_ad(npx-2, j) + temp_ad4
          q_ad(npx+1, j) = q_ad(npx+1, j) + temp_ad4
          CALL POPREAL8(bl(npx-2))
          q_ad(npx-2, j) = q_ad(npx-2, j) + (p1-1.0)*bl_ad(npx-2)
          q_ad(npx-3, j) = q_ad(npx-3, j) + p1*bl_ad(npx-2)
          q_ad(npx-4, j) = q_ad(npx-4, j) + p2*bl_ad(npx-2)
          q_ad(npx-1, j) = q_ad(npx-1, j) + p2*bl_ad(npx-2)
          bl_ad(npx-2) = 0.0_8
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          CALL POPREAL8(bl(2))
          xt_ad = br_ad(1) + bl_ad(2)
          q_ad(2, j) = q_ad(2, j) - bl_ad(2)
          bl_ad(2) = 0.0_8
          CALL POPREAL8(br(1))
          q_ad(1, j) = q_ad(1, j) - br_ad(1)
          br_ad(1) = 0.0_8
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            y4_ad = xt_ad
            xt_ad = 0.0_8
          ELSE
            y4_ad = 0.0_8
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            q_ad(2, j) = q_ad(2, j) + y4_ad
          ELSE
            q_ad(1, j) = q_ad(1, j) + y4_ad
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            y3_ad = xt_ad
            xt_ad = 0.0_8
          ELSE
            y3_ad = 0.0_8
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            q_ad(2, j) = q_ad(2, j) + y3_ad
          ELSE
            q_ad(1, j) = q_ad(1, j) + y3_ad
          END IF
          q_ad(1, j) = q_ad(1, j) + c3*xt_ad
          q_ad(2, j) = q_ad(2, j) + c2*xt_ad
          q_ad(3, j) = q_ad(3, j) + c1*xt_ad
          CALL POPREAL8(bl(0))
          xt_ad = bl_ad(0)
          q_ad(0, j) = q_ad(0, j) - bl_ad(0)
          bl_ad(0) = 0.0_8
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            y2_ad = xt_ad
            xt_ad = 0.0_8
          ELSE
            y2_ad = 0.0_8
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            q_ad(0, j) = q_ad(0, j) + y2_ad
          ELSE
            q_ad(-1, j) = q_ad(-1, j) + y2_ad
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            y1_ad = xt_ad
            xt_ad = 0.0_8
          ELSE
            y1_ad = 0.0_8
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            q_ad(0, j) = q_ad(0, j) + y1_ad
          ELSE
            q_ad(-1, j) = q_ad(-1, j) + y1_ad
          END IF
          q_ad(-2, j) = q_ad(-2, j) + c1*xt_ad
          q_ad(-1, j) = q_ad(-1, j) + c2*xt_ad
          q_ad(0, j) = q_ad(0, j) + c3*xt_ad - br_ad(0)
          CALL POPREAL8(br(0))
          xt_ad = bl_ad(1) + br_ad(0)
          br_ad(0) = 0.0_8
          CALL POPREAL8(bl(1))
          q_ad(1, j) = q_ad(1, j) - bl_ad(1)
          bl_ad(1) = 0.0_8
          temp_ad = 0.5*xt_ad/(dxa(1, j)+dxa(2, j))
          temp_ad0 = (dxa(1, j)*2.+dxa(2, j))*temp_ad
          temp_ad1 = -(dxa(1, j)*temp_ad)
          q_ad(0, j) = q_ad(0, j) + temp_ad0
          q_ad(1, j) = q_ad(1, j) + temp_ad0
          q_ad(-1, j) = q_ad(-1, j) + temp_ad1
          q_ad(2, j) = q_ad(2, j) + (p1-1.0)*br_ad(2) + temp_ad1
          CALL POPREAL8(br(2))
          q_ad(3, j) = q_ad(3, j) + p1*br_ad(2)
          q_ad(1, j) = q_ad(1, j) + p2*br_ad(2)
          q_ad(4, j) = q_ad(4, j) + p2*br_ad(2)
          br_ad(2) = 0.0_8
        END IF
        DO i=ie3,is3,-1
          CALL POPREAL8(br(i))
          q_ad(i-2, j) = q_ad(i-2, j) + b1*br_ad(i)
          q_ad(i-1, j) = q_ad(i-1, j) + b2*br_ad(i)
          q_ad(i, j) = q_ad(i, j) + b3*br_ad(i)
          q_ad(i+1, j) = q_ad(i+1, j) + b4*br_ad(i)
          q_ad(i+2, j) = q_ad(i+2, j) + b5*br_ad(i)
          br_ad(i) = 0.0_8
          CALL POPREAL8(bl(i))
          q_ad(i-2, j) = q_ad(i-2, j) + b5*bl_ad(i)
          q_ad(i-1, j) = q_ad(i-1, j) + b4*bl_ad(i)
          q_ad(i, j) = q_ad(i, j) + b3*bl_ad(i)
          q_ad(i+1, j) = q_ad(i+1, j) + b2*bl_ad(i)
          q_ad(i+2, j) = q_ad(i+2, j) + b1*bl_ad(i)
          bl_ad(i) = 0.0_8
        END DO
      END DO
    END IF
  END SUBROUTINE FXPPM_ADM
  SUBROUTINE FXPPM(c, q, flux, iord, ifirst, ilast, jfirst, jlast, npx, &
&   npy, ppm_limiter)
    IMPLICIT NONE
! !INPUT PARAMETERS:
!  X-Dir strip
    INTEGER, INTENT(IN) :: ifirst, ilast
!  Y-Dir strip
    INTEGER, INTENT(IN) :: jfirst, jlast
    INTEGER, INTENT(IN) :: iord
    INTEGER, INTENT(IN) :: npx, npy
    REAL, INTENT(IN) :: q(ifirst-ng:ilast+ng, jfirst:jlast)
! Courant   N (like FLUX)
    REAL, INTENT(IN) :: c(ifirst:ilast+1, jfirst:jlast)
    REAL, INTENT(IN) :: ppm_limiter
! !OUTPUT PARAMETERS:
!  Flux
    REAL, INTENT(OUT) :: flux(ifirst:ilast+1, jfirst:jlast)
! Local
    REAL :: dm1(ifirst-2:ilast+2)
    REAL :: al(ifirst-1:ilast+2)
    REAL :: bl(ifirst-1:ilast+1)
    REAL :: br(ifirst-1:ilast+1)
    REAL :: dq(ifirst-3:ilast+2)
    REAL :: dl, dr, pmp, lac, ct, qe
    REAL :: xt, x1, x0
    INTEGER :: i, j, is3, ie3, ie2, it
    INTRINSIC MAX
    INTRINSIC MIN
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
      DO j=jfirst,jlast
! Non-monotonic "5th order" scheme (not really 5th order)
        DO i=is3,ie3
          bl(i) = b5*q(i-2, j) + b4*q(i-1, j) + b3*q(i, j) + b2*q(i+1, j&
&           ) + b1*q(i+2, j)
          br(i) = b1*q(i-2, j) + b2*q(i-1, j) + b3*q(i, j) + b4*q(i+1, j&
&           ) + b5*q(i+2, j)
        END DO
!--------------
! fix the edges
!--------------
        IF (is .EQ. 1) THEN
          br(2) = p1*(q(2, j)+q(3, j)) + p2*(q(1, j)+q(4, j)) - q(2, j)
          xt = 0.5*((2.*dxa(1, j)+dxa(2, j))*(q(0, j)+q(1, j))-dxa(1, j)&
&           *(q(-1, j)+q(2, j)))/(dxa(1, j)+dxa(2, j))
          bl(1) = xt - q(1, j)
          br(0) = xt - q(0, j)
          xt = c1*q(-2, j) + c2*q(-1, j) + c3*q(0, j)
          IF (q(-1, j) .GT. q(0, j)) THEN
            y1 = q(0, j)
          ELSE
            y1 = q(-1, j)
          END IF
          IF (xt .LT. y1) THEN
            xt = y1
          ELSE
            xt = xt
          END IF
          IF (q(-1, j) .LT. q(0, j)) THEN
            y2 = q(0, j)
          ELSE
            y2 = q(-1, j)
          END IF
          IF (xt .GT. y2) THEN
            xt = y2
          ELSE
            xt = xt
          END IF
          bl(0) = xt - q(0, j)
          xt = c3*q(1, j) + c2*q(2, j) + c1*q(3, j)
          IF (q(1, j) .GT. q(2, j)) THEN
            y3 = q(2, j)
          ELSE
            y3 = q(1, j)
          END IF
          IF (xt .LT. y3) THEN
            xt = y3
          ELSE
            xt = xt
          END IF
          IF (q(1, j) .LT. q(2, j)) THEN
            y4 = q(2, j)
          ELSE
            y4 = q(1, j)
          END IF
          IF (xt .GT. y4) THEN
            xt = y4
          ELSE
            xt = xt
          END IF
          br(1) = xt - q(1, j)
          bl(2) = xt - q(2, j)
        END IF
        IF (ie + 1 .EQ. npx) THEN
          bl(npx-2) = p1*(q(npx-2, j)+q(npx-3, j)) + p2*(q(npx-4, j)+q(&
&           npx-1, j)) - q(npx-2, j)
          xt = 0.5*((2.*dxa(npx-1, j)+dxa(npx-2, j))*(q(npx-1, j)+q(npx&
&           , j))-dxa(npx-1, j)*(q(npx-2, j)+q(npx+1, j)))/(dxa(npx-1, j&
&           )+dxa(npx-2, j))
          br(npx-1) = xt - q(npx-1, j)
          bl(npx) = xt - q(npx, j)
          xt = c3*q(npx, j) + c2*q(npx+1, j) + c1*q(npx+2, j)
          IF (q(npx, j) .GT. q(npx+1, j)) THEN
            y5 = q(npx+1, j)
          ELSE
            y5 = q(npx, j)
          END IF
          IF (xt .LT. y5) THEN
            xt = y5
          ELSE
            xt = xt
          END IF
          IF (q(npx, j) .LT. q(npx+1, j)) THEN
            y6 = q(npx+1, j)
          ELSE
            y6 = q(npx, j)
          END IF
          IF (xt .GT. y6) THEN
            xt = y6
          ELSE
            xt = xt
          END IF
          br(npx) = xt - q(npx, j)
          xt = c1*q(npx-3, j) + c2*q(npx-2, j) + c3*q(npx-1, j)
          IF (q(npx-2, j) .GT. q(npx-1, j)) THEN
            y7 = q(npx-1, j)
          ELSE
            y7 = q(npx-2, j)
          END IF
          IF (xt .LT. y7) THEN
            xt = y7
          ELSE
            xt = xt
          END IF
          IF (q(npx-2, j) .LT. q(npx-1, j)) THEN
            y8 = q(npx-1, j)
          ELSE
            y8 = q(npx-2, j)
          END IF
          IF (xt .GT. y8) THEN
            xt = y8
          ELSE
            xt = xt
          END IF
          br(npx-2) = xt - q(npx-2, j)
          bl(npx-1) = xt - q(npx-1, j)
        END IF
        DO i=ifirst,ilast+1
          IF (c(i, j) .GT. 0.) THEN
            flux(i, j) = q(i-1, j) + (1.-c(i, j))*(br(i-1)-c(i, j)*(bl(i&
&             -1)+br(i-1)))
          ELSE
            flux(i, j) = q(i, j) + (1.+c(i, j))*(bl(i)+c(i, j)*(bl(i)+br&
&             (i)))
          END IF
        END DO
      END DO
    END IF
  END SUBROUTINE FXPPM
!  Differentiation of fyppm in reverse (adjoint) mode (with options r8):
!   gradient     of useful results: q flux c
!   with respect to varying inputs: q flux c
  SUBROUTINE FYPPM_ADM(c, c_ad, q, q_ad, flux, flux_ad, jord, ifirst, &
&   ilast, jfirst, jlast, npx, npy, ppm_limiter)
    IMPLICIT NONE
!  X-Dir strip
    INTEGER, INTENT(IN) :: ifirst, ilast
!  Y-Dir strip
    INTEGER, INTENT(IN) :: jfirst, jlast
    INTEGER, INTENT(IN) :: jord
    INTEGER, INTENT(IN) :: npx, npy
    REAL, INTENT(IN) :: q(ifirst:ilast, jfirst-ng:jlast+ng)
    REAL :: q_ad(ifirst:ilast, jfirst-ng:jlast+ng)
! Courant number
    REAL, INTENT(IN) :: c(isd:ied, js:je+1)
    REAL :: c_ad(isd:ied, js:je+1)
    REAL, INTENT(IN) :: ppm_limiter
!  Flux
    REAL :: flux(ifirst:ilast, jfirst:jlast+1)
    REAL :: flux_ad(ifirst:ilast, jfirst:jlast+1)
! Local:
    REAL :: al(ifirst:ilast, jfirst-1:jlast+2)
    REAL :: bl(ifirst:ilast, jfirst-1:jlast+1)
    REAL :: bl_ad(ifirst:ilast, jfirst-1:jlast+1)
    REAL :: br(ifirst:ilast, jfirst-1:jlast+1)
    REAL :: br_ad(ifirst:ilast, jfirst-1:jlast+1)
    REAL :: dq(ifirst:ilast, jfirst-3:jlast+2)
    REAL :: dl, dr, pmp, lac, ct, qe
    REAL :: xt, x0, x1
    REAL :: xt_ad
    INTEGER :: i, j, js3, je3, je2, jt
    INTRINSIC MAX
    INTRINSIC MIN
    INTEGER :: branch
    REAL :: temp0
    REAL :: y2_ad
    REAL :: y3_ad
    INTEGER :: min1
    REAL :: y4_ad
    REAL :: y5_ad
    REAL :: y6_ad
    REAL :: y7_ad
    REAL :: temp_ad8
    REAL :: temp_ad7
    REAL :: temp_ad6
    REAL :: temp_ad5
    REAL :: temp_ad4
    REAL :: temp_ad3
    REAL :: temp_ad2
    REAL :: temp_ad1
    REAL :: y8_ad
    REAL :: temp_ad0
    REAL :: temp_ad
    REAL :: temp
    REAL :: y8
    REAL :: y7
    REAL :: y6
    INTEGER :: max1
    REAL :: y5
    REAL :: y4
    REAL :: y3
    REAL :: y1_ad
    REAL :: y2
    REAL :: y1
    IF (jord .EQ. 6) THEN
      IF (3 .LT. js - 1) THEN
        max1 = js - 1
      ELSE
        max1 = 3
      END IF
      IF (npy - 3 .GT. je + 1) THEN
        min1 = je + 1
      ELSE
        min1 = npy - 3
      END IF
! Non-monotonic "5th order" scheme (not really 5th order)
      DO j=max1,min1
        DO i=ifirst,ilast
          bl(i, j) = b5*q(i, j-2) + b4*q(i, j-1) + b3*q(i, j) + b2*q(i, &
&           j+1) + b1*q(i, j+2)
          br(i, j) = b1*q(i, j-2) + b2*q(i, j-1) + b3*q(i, j) + b4*q(i, &
&           j+1) + b5*q(i, j+2)
        END DO
      END DO
      IF (js .EQ. 1) THEN
        DO i=ifirst,ilast
!           br(i,2) = al(i,3) - q(i,2)
          br(i, 2) = p1*(q(i, 2)+q(i, 3)) + p2*(q(i, 1)+q(i, 4)) - q(i, &
&           2)
          xt = 0.5*((2.*dya(i, 1)+dya(i, 2))*(q(i, 0)+q(i, 1))-dya(i, 1)&
&           *(q(i, -1)+q(i, 2)))/(dya(i, 1)+dya(i, 2))
          bl(i, 1) = xt - q(i, 1)
          br(i, 0) = xt - q(i, 0)
!           xt = s14*0.25*(q(i,0)-q(i,-2)) - s11*(q(i,0)-q(i,-1)) + q(i,0)
          xt = c1*q(i, -2) + c2*q(i, -1) + c3*q(i, 0)
          IF (q(i, -1) .LT. q(i, 0)) THEN
            y1 = q(i, 0)
            CALL PUSHCONTROL1B(0)
          ELSE
            y1 = q(i, -1)
            CALL PUSHCONTROL1B(1)
          END IF
          IF (xt .GT. y1) THEN
            xt = y1
            CALL PUSHCONTROL1B(0)
          ELSE
            xt = xt
            CALL PUSHCONTROL1B(1)
          END IF
          IF (q(i, -1) .GT. q(i, 0)) THEN
            y2 = q(i, 0)
            CALL PUSHCONTROL1B(0)
          ELSE
            y2 = q(i, -1)
            CALL PUSHCONTROL1B(1)
          END IF
          IF (xt .LT. y2) THEN
            xt = y2
            CALL PUSHCONTROL1B(0)
          ELSE
            xt = xt
            CALL PUSHCONTROL1B(1)
          END IF
          bl(i, 0) = xt - q(i, 0)
!           xt = s15*q(i,1) + s11*q(i,2) - s14*0.25*(q(i,3)-q(i,1))
          xt = c3*q(i, 1) + c2*q(i, 2) + c1*q(i, 3)
          IF (q(i, 1) .LT. q(i, 2)) THEN
            y3 = q(i, 2)
            CALL PUSHCONTROL1B(0)
          ELSE
            y3 = q(i, 1)
            CALL PUSHCONTROL1B(1)
          END IF
          IF (xt .GT. y3) THEN
            xt = y3
            CALL PUSHCONTROL1B(0)
          ELSE
            xt = xt
            CALL PUSHCONTROL1B(1)
          END IF
          IF (q(i, 1) .GT. q(i, 2)) THEN
            y4 = q(i, 2)
            CALL PUSHCONTROL1B(0)
          ELSE
            y4 = q(i, 1)
            CALL PUSHCONTROL1B(1)
          END IF
          IF (xt .LT. y4) THEN
            xt = y4
            CALL PUSHCONTROL1B(0)
          ELSE
            xt = xt
            CALL PUSHCONTROL1B(1)
          END IF
          br(i, 1) = xt - q(i, 1)
          bl(i, 2) = xt - q(i, 2)
        END DO
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (je + 1 .EQ. npy) THEN
        DO i=ifirst,ilast
!           bl(i,npy-2) = al(i,npy-2) - q(i,npy-2)
          bl(i, npy-2) = p1*(q(i, npy-3)+q(i, npy-2)) + p2*(q(i, npy-4)+&
&           q(i, npy-1)) - q(i, npy-2)
          xt = 0.5*((2.*dya(i, npy-1)+dya(i, npy-2))*(q(i, npy-1)+q(i, &
&           npy))-dya(i, npy-1)*(q(i, npy-2)+q(i, npy+1)))/(dya(i, npy-1&
&           )+dya(i, npy-2))
          br(i, npy-1) = xt - q(i, npy-1)
          bl(i, npy) = xt - q(i, npy)
!           xt = s11*(q(i,npy+1)-q(i,npy)) - s14*0.25*(q(i,npy+2)-q(i,npy)) + q(i,npy)
          xt = c3*q(i, npy) + c2*q(i, npy+1) + c1*q(i, npy+2)
          IF (q(i, npy) .LT. q(i, npy+1)) THEN
            y5 = q(i, npy+1)
            CALL PUSHCONTROL1B(0)
          ELSE
            y5 = q(i, npy)
            CALL PUSHCONTROL1B(1)
          END IF
          IF (xt .GT. y5) THEN
            xt = y5
            CALL PUSHCONTROL1B(0)
          ELSE
            xt = xt
            CALL PUSHCONTROL1B(1)
          END IF
          IF (q(i, npy) .GT. q(i, npy+1)) THEN
            y6 = q(i, npy+1)
            CALL PUSHCONTROL1B(0)
          ELSE
            y6 = q(i, npy)
            CALL PUSHCONTROL1B(1)
          END IF
          IF (xt .LT. y6) THEN
            xt = y6
            CALL PUSHCONTROL1B(0)
          ELSE
            xt = xt
            CALL PUSHCONTROL1B(1)
          END IF
          br(i, npy) = xt - q(i, npy)
!           xt = s15*q(i,npy-1) + s11*q(i,npy-2) + s14*0.25*(q(i,npy-1)-q(i,npy-3))
          xt = c1*q(i, npy-3) + c2*q(i, npy-2) + c3*q(i, npy-1)
          IF (q(i, npy-2) .LT. q(i, npy-1)) THEN
            y7 = q(i, npy-1)
            CALL PUSHCONTROL1B(0)
          ELSE
            y7 = q(i, npy-2)
            CALL PUSHCONTROL1B(1)
          END IF
          IF (xt .GT. y7) THEN
            xt = y7
            CALL PUSHCONTROL1B(0)
          ELSE
            xt = xt
            CALL PUSHCONTROL1B(1)
          END IF
          IF (q(i, npy-2) .GT. q(i, npy-1)) THEN
            y8 = q(i, npy-1)
            CALL PUSHCONTROL1B(0)
          ELSE
            y8 = q(i, npy-2)
            CALL PUSHCONTROL1B(1)
          END IF
          IF (xt .LT. y8) THEN
            xt = y8
            CALL PUSHCONTROL1B(0)
          ELSE
            xt = xt
            CALL PUSHCONTROL1B(1)
          END IF
          br(i, npy-2) = xt - q(i, npy-2)
          bl(i, npy-1) = xt - q(i, npy-1)
        END DO
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL PUSHCONTROL1B(0)
      END IF
      DO j=jfirst,jlast+1
        DO i=ifirst,ilast
          IF (c(i, j) .GT. 0.) THEN
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
      END DO
      bl_ad = 0.0_8
      br_ad = 0.0_8
      DO j=jlast+1,jfirst,-1
        DO i=ilast,ifirst,-1
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            temp0 = bl(i, j) + br(i, j)
            temp_ad7 = (c(i, j)+1.)*flux_ad(i, j)
            temp_ad8 = c(i, j)*temp_ad7
            q_ad(i, j) = q_ad(i, j) + flux_ad(i, j)
            c_ad(i, j) = c_ad(i, j) + temp0*temp_ad7 + (bl(i, j)+c(i, j)&
&             *temp0)*flux_ad(i, j)
            bl_ad(i, j) = bl_ad(i, j) + temp_ad8 + temp_ad7
            br_ad(i, j) = br_ad(i, j) + temp_ad8
            flux_ad(i, j) = 0.0_8
          ELSE
            temp = bl(i, j-1) + br(i, j-1)
            temp_ad5 = (1.-c(i, j))*flux_ad(i, j)
            temp_ad6 = -(c(i, j)*temp_ad5)
            q_ad(i, j-1) = q_ad(i, j-1) + flux_ad(i, j)
            c_ad(i, j) = c_ad(i, j) - temp*temp_ad5 - (br(i, j-1)-c(i, j&
&             )*temp)*flux_ad(i, j)
            br_ad(i, j-1) = br_ad(i, j-1) + temp_ad6 + temp_ad5
            bl_ad(i, j-1) = bl_ad(i, j-1) + temp_ad6
            flux_ad(i, j) = 0.0_8
          END IF
        END DO
      END DO
      CALL POPCONTROL1B(branch)
      IF (branch .NE. 0) THEN
        DO i=ilast,ifirst,-1
          xt_ad = br_ad(i, npy-2) + bl_ad(i, npy-1)
          q_ad(i, npy-1) = q_ad(i, npy-1) - bl_ad(i, npy-1)
          bl_ad(i, npy-1) = 0.0_8
          q_ad(i, npy-2) = q_ad(i, npy-2) - br_ad(i, npy-2)
          br_ad(i, npy-2) = 0.0_8
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            y8_ad = xt_ad
            xt_ad = 0.0_8
          ELSE
            y8_ad = 0.0_8
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            q_ad(i, npy-1) = q_ad(i, npy-1) + y8_ad
          ELSE
            q_ad(i, npy-2) = q_ad(i, npy-2) + y8_ad
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            y7_ad = xt_ad
            xt_ad = 0.0_8
          ELSE
            y7_ad = 0.0_8
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            q_ad(i, npy-1) = q_ad(i, npy-1) + y7_ad
          ELSE
            q_ad(i, npy-2) = q_ad(i, npy-2) + y7_ad
          END IF
          q_ad(i, npy-3) = q_ad(i, npy-3) + c1*xt_ad
          q_ad(i, npy-2) = q_ad(i, npy-2) + c2*xt_ad
          q_ad(i, npy-1) = q_ad(i, npy-1) + c3*xt_ad
          xt_ad = br_ad(i, npy)
          q_ad(i, npy) = q_ad(i, npy) - br_ad(i, npy)
          br_ad(i, npy) = 0.0_8
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            y6_ad = xt_ad
            xt_ad = 0.0_8
          ELSE
            y6_ad = 0.0_8
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            q_ad(i, npy+1) = q_ad(i, npy+1) + y6_ad
          ELSE
            q_ad(i, npy) = q_ad(i, npy) + y6_ad
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            y5_ad = xt_ad
            xt_ad = 0.0_8
          ELSE
            y5_ad = 0.0_8
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            q_ad(i, npy+1) = q_ad(i, npy+1) + y5_ad
          ELSE
            q_ad(i, npy) = q_ad(i, npy) + y5_ad
          END IF
          q_ad(i, npy) = q_ad(i, npy) + c3*xt_ad
          q_ad(i, npy+1) = q_ad(i, npy+1) + c2*xt_ad
          q_ad(i, npy+2) = q_ad(i, npy+2) + c1*xt_ad
          xt_ad = br_ad(i, npy-1) + bl_ad(i, npy)
          q_ad(i, npy) = q_ad(i, npy) - bl_ad(i, npy)
          bl_ad(i, npy) = 0.0_8
          temp_ad3 = 0.5*xt_ad/(dya(i, npy-1)+dya(i, npy-2))
          temp_ad2 = (dya(i, npy-1)*2.+dya(i, npy-2))*temp_ad3
          q_ad(i, npy-1) = q_ad(i, npy-1) + temp_ad2 - br_ad(i, npy-1)
          br_ad(i, npy-1) = 0.0_8
          temp_ad4 = -(dya(i, npy-1)*temp_ad3)
          q_ad(i, npy) = q_ad(i, npy) + temp_ad2
          q_ad(i, npy-2) = q_ad(i, npy-2) + temp_ad4
          q_ad(i, npy+1) = q_ad(i, npy+1) + temp_ad4
          q_ad(i, npy-3) = q_ad(i, npy-3) + p1*bl_ad(i, npy-2)
          q_ad(i, npy-2) = q_ad(i, npy-2) + (p1-1.0)*bl_ad(i, npy-2)
          q_ad(i, npy-4) = q_ad(i, npy-4) + p2*bl_ad(i, npy-2)
          q_ad(i, npy-1) = q_ad(i, npy-1) + p2*bl_ad(i, npy-2)
          bl_ad(i, npy-2) = 0.0_8
        END DO
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO i=ilast,ifirst,-1
          xt_ad = br_ad(i, 1) + bl_ad(i, 2)
          q_ad(i, 2) = q_ad(i, 2) - bl_ad(i, 2)
          bl_ad(i, 2) = 0.0_8
          q_ad(i, 1) = q_ad(i, 1) - br_ad(i, 1)
          br_ad(i, 1) = 0.0_8
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            y4_ad = xt_ad
            xt_ad = 0.0_8
          ELSE
            y4_ad = 0.0_8
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            q_ad(i, 2) = q_ad(i, 2) + y4_ad
          ELSE
            q_ad(i, 1) = q_ad(i, 1) + y4_ad
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            y3_ad = xt_ad
            xt_ad = 0.0_8
          ELSE
            y3_ad = 0.0_8
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            q_ad(i, 2) = q_ad(i, 2) + y3_ad
          ELSE
            q_ad(i, 1) = q_ad(i, 1) + y3_ad
          END IF
          q_ad(i, 1) = q_ad(i, 1) + c3*xt_ad
          q_ad(i, 2) = q_ad(i, 2) + c2*xt_ad
          q_ad(i, 3) = q_ad(i, 3) + c1*xt_ad
          xt_ad = bl_ad(i, 0)
          q_ad(i, 0) = q_ad(i, 0) - bl_ad(i, 0)
          bl_ad(i, 0) = 0.0_8
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            y2_ad = xt_ad
            xt_ad = 0.0_8
          ELSE
            y2_ad = 0.0_8
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            q_ad(i, 0) = q_ad(i, 0) + y2_ad
          ELSE
            q_ad(i, -1) = q_ad(i, -1) + y2_ad
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            y1_ad = xt_ad
            xt_ad = 0.0_8
          ELSE
            y1_ad = 0.0_8
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            q_ad(i, 0) = q_ad(i, 0) + y1_ad
          ELSE
            q_ad(i, -1) = q_ad(i, -1) + y1_ad
          END IF
          q_ad(i, -2) = q_ad(i, -2) + c1*xt_ad
          q_ad(i, -1) = q_ad(i, -1) + c2*xt_ad
          q_ad(i, 0) = q_ad(i, 0) + c3*xt_ad - br_ad(i, 0)
          xt_ad = bl_ad(i, 1) + br_ad(i, 0)
          br_ad(i, 0) = 0.0_8
          q_ad(i, 1) = q_ad(i, 1) - bl_ad(i, 1)
          bl_ad(i, 1) = 0.0_8
          temp_ad = 0.5*xt_ad/(dya(i, 1)+dya(i, 2))
          temp_ad0 = (dya(i, 1)*2.+dya(i, 2))*temp_ad
          temp_ad1 = -(dya(i, 1)*temp_ad)
          q_ad(i, 0) = q_ad(i, 0) + temp_ad0
          q_ad(i, 1) = q_ad(i, 1) + temp_ad0
          q_ad(i, -1) = q_ad(i, -1) + temp_ad1
          q_ad(i, 2) = q_ad(i, 2) + (p1-1.0)*br_ad(i, 2) + temp_ad1
          q_ad(i, 3) = q_ad(i, 3) + p1*br_ad(i, 2)
          q_ad(i, 1) = q_ad(i, 1) + p2*br_ad(i, 2)
          q_ad(i, 4) = q_ad(i, 4) + p2*br_ad(i, 2)
          br_ad(i, 2) = 0.0_8
        END DO
      END IF
      DO j=min1,max1,-1
        DO i=ilast,ifirst,-1
          q_ad(i, j-2) = q_ad(i, j-2) + b1*br_ad(i, j)
          q_ad(i, j-1) = q_ad(i, j-1) + b2*br_ad(i, j)
          q_ad(i, j) = q_ad(i, j) + b3*br_ad(i, j)
          q_ad(i, j+1) = q_ad(i, j+1) + b4*br_ad(i, j)
          q_ad(i, j+2) = q_ad(i, j+2) + b5*br_ad(i, j)
          br_ad(i, j) = 0.0_8
          q_ad(i, j-2) = q_ad(i, j-2) + b5*bl_ad(i, j)
          q_ad(i, j-1) = q_ad(i, j-1) + b4*bl_ad(i, j)
          q_ad(i, j) = q_ad(i, j) + b3*bl_ad(i, j)
          q_ad(i, j+1) = q_ad(i, j+1) + b2*bl_ad(i, j)
          q_ad(i, j+2) = q_ad(i, j+2) + b1*bl_ad(i, j)
          bl_ad(i, j) = 0.0_8
        END DO
      END DO
    END IF
  END SUBROUTINE FYPPM_ADM
  SUBROUTINE FYPPM(c, q, flux, jord, ifirst, ilast, jfirst, jlast, npx, &
&   npy, ppm_limiter)
    IMPLICIT NONE
!  X-Dir strip
    INTEGER, INTENT(IN) :: ifirst, ilast
!  Y-Dir strip
    INTEGER, INTENT(IN) :: jfirst, jlast
    INTEGER, INTENT(IN) :: jord
    INTEGER, INTENT(IN) :: npx, npy
    REAL, INTENT(IN) :: q(ifirst:ilast, jfirst-ng:jlast+ng)
! Courant number
    REAL, INTENT(IN) :: c(isd:ied, js:je+1)
    REAL, INTENT(IN) :: ppm_limiter
!  Flux
    REAL, INTENT(OUT) :: flux(ifirst:ilast, jfirst:jlast+1)
! Local:
    REAL :: al(ifirst:ilast, jfirst-1:jlast+2)
    REAL :: bl(ifirst:ilast, jfirst-1:jlast+1)
    REAL :: br(ifirst:ilast, jfirst-1:jlast+1)
    REAL :: dq(ifirst:ilast, jfirst-3:jlast+2)
    REAL :: dl, dr, pmp, lac, ct, qe
    REAL :: xt, x0, x1
    INTEGER :: i, j, js3, je3, je2, jt
    INTRINSIC MAX
    INTRINSIC MIN
    INTEGER :: min1
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
      ELSE
        min1 = npy - 3
      END IF
! Non-monotonic "5th order" scheme (not really 5th order)
      DO j=max1,min1
        DO i=ifirst,ilast
          bl(i, j) = b5*q(i, j-2) + b4*q(i, j-1) + b3*q(i, j) + b2*q(i, &
&           j+1) + b1*q(i, j+2)
          br(i, j) = b1*q(i, j-2) + b2*q(i, j-1) + b3*q(i, j) + b4*q(i, &
&           j+1) + b5*q(i, j+2)
        END DO
      END DO
      IF (js .EQ. 1) THEN
        DO i=ifirst,ilast
!           br(i,2) = al(i,3) - q(i,2)
          br(i, 2) = p1*(q(i, 2)+q(i, 3)) + p2*(q(i, 1)+q(i, 4)) - q(i, &
&           2)
          xt = 0.5*((2.*dya(i, 1)+dya(i, 2))*(q(i, 0)+q(i, 1))-dya(i, 1)&
&           *(q(i, -1)+q(i, 2)))/(dya(i, 1)+dya(i, 2))
          bl(i, 1) = xt - q(i, 1)
          br(i, 0) = xt - q(i, 0)
!           xt = s14*0.25*(q(i,0)-q(i,-2)) - s11*(q(i,0)-q(i,-1)) + q(i,0)
          xt = c1*q(i, -2) + c2*q(i, -1) + c3*q(i, 0)
          IF (q(i, -1) .LT. q(i, 0)) THEN
            y1 = q(i, 0)
          ELSE
            y1 = q(i, -1)
          END IF
          IF (xt .GT. y1) THEN
            xt = y1
          ELSE
            xt = xt
          END IF
          IF (q(i, -1) .GT. q(i, 0)) THEN
            y2 = q(i, 0)
          ELSE
            y2 = q(i, -1)
          END IF
          IF (xt .LT. y2) THEN
            xt = y2
          ELSE
            xt = xt
          END IF
          bl(i, 0) = xt - q(i, 0)
!           xt = s15*q(i,1) + s11*q(i,2) - s14*0.25*(q(i,3)-q(i,1))
          xt = c3*q(i, 1) + c2*q(i, 2) + c1*q(i, 3)
          IF (q(i, 1) .LT. q(i, 2)) THEN
            y3 = q(i, 2)
          ELSE
            y3 = q(i, 1)
          END IF
          IF (xt .GT. y3) THEN
            xt = y3
          ELSE
            xt = xt
          END IF
          IF (q(i, 1) .GT. q(i, 2)) THEN
            y4 = q(i, 2)
          ELSE
            y4 = q(i, 1)
          END IF
          IF (xt .LT. y4) THEN
            xt = y4
          ELSE
            xt = xt
          END IF
          br(i, 1) = xt - q(i, 1)
          bl(i, 2) = xt - q(i, 2)
        END DO
      END IF
      IF (je + 1 .EQ. npy) THEN
        DO i=ifirst,ilast
!           bl(i,npy-2) = al(i,npy-2) - q(i,npy-2)
          bl(i, npy-2) = p1*(q(i, npy-3)+q(i, npy-2)) + p2*(q(i, npy-4)+&
&           q(i, npy-1)) - q(i, npy-2)
          xt = 0.5*((2.*dya(i, npy-1)+dya(i, npy-2))*(q(i, npy-1)+q(i, &
&           npy))-dya(i, npy-1)*(q(i, npy-2)+q(i, npy+1)))/(dya(i, npy-1&
&           )+dya(i, npy-2))
          br(i, npy-1) = xt - q(i, npy-1)
          bl(i, npy) = xt - q(i, npy)
!           xt = s11*(q(i,npy+1)-q(i,npy)) - s14*0.25*(q(i,npy+2)-q(i,npy)) + q(i,npy)
          xt = c3*q(i, npy) + c2*q(i, npy+1) + c1*q(i, npy+2)
          IF (q(i, npy) .LT. q(i, npy+1)) THEN
            y5 = q(i, npy+1)
          ELSE
            y5 = q(i, npy)
          END IF
          IF (xt .GT. y5) THEN
            xt = y5
          ELSE
            xt = xt
          END IF
          IF (q(i, npy) .GT. q(i, npy+1)) THEN
            y6 = q(i, npy+1)
          ELSE
            y6 = q(i, npy)
          END IF
          IF (xt .LT. y6) THEN
            xt = y6
          ELSE
            xt = xt
          END IF
          br(i, npy) = xt - q(i, npy)
!           xt = s15*q(i,npy-1) + s11*q(i,npy-2) + s14*0.25*(q(i,npy-1)-q(i,npy-3))
          xt = c1*q(i, npy-3) + c2*q(i, npy-2) + c3*q(i, npy-1)
          IF (q(i, npy-2) .LT. q(i, npy-1)) THEN
            y7 = q(i, npy-1)
          ELSE
            y7 = q(i, npy-2)
          END IF
          IF (xt .GT. y7) THEN
            xt = y7
          ELSE
            xt = xt
          END IF
          IF (q(i, npy-2) .GT. q(i, npy-1)) THEN
            y8 = q(i, npy-1)
          ELSE
            y8 = q(i, npy-2)
          END IF
          IF (xt .LT. y8) THEN
            xt = y8
          ELSE
            xt = xt
          END IF
          br(i, npy-2) = xt - q(i, npy-2)
          bl(i, npy-1) = xt - q(i, npy-1)
        END DO
      END IF
      DO j=jfirst,jlast+1
        DO i=ifirst,ilast
          IF (c(i, j) .GT. 0.) THEN
            flux(i, j) = q(i, j-1) + (1.-c(i, j))*(br(i, j-1)-c(i, j)*(&
&             bl(i, j-1)+br(i, j-1)))
          ELSE
            flux(i, j) = q(i, j) + (1.+c(i, j))*(bl(i, j)+c(i, j)*(bl(i&
&             , j)+br(i, j)))
          END IF
        END DO
      END DO
    END IF
  END SUBROUTINE FYPPM
!  Differentiation of deln_flux in reverse (adjoint) mode (with options r8):
!   gradient     of useful results: q fx fy
!   with respect to varying inputs: q fx fy
  SUBROUTINE DELN_FLUX_ADM(nord, npx, npy, damp, q, q_ad, fx, fx_ad, fy&
&   , fy_ad, mfx, mfy)
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
    REAL :: q_ad(is-ng:ie+ng, js-ng:je+ng)
! diffusive fluxes: 
    REAL, INTENT(IN) :: mfx(is:ie+1, js:je), mfy(is:ie, js:je+1)
    REAL, INTENT(INOUT) :: fx(is:ie+1, js:je), fy(is:ie, js:je+1)
    REAL, INTENT(INOUT) :: fx_ad(is:ie+1, js:je)
! local:
    REAL :: fx2(isd:ied+1, jsd:jed), fy2(isd:ied, jsd:jed+1)
    REAL :: fx2_ad(isd:ied+1, jsd:jed), fy2_ad(isd:ied, jsd:jed+1)
    REAL :: d2(isd:ied, jsd:jed)
    REAL :: d2_ad(isd:ied, jsd:jed)
    INTEGER :: i, j, n, nt
    INTEGER :: ad_from
    INTEGER :: ad_to
    INTEGER :: ad_from0
    INTEGER :: ad_to0
    INTEGER :: ad_from1
    INTEGER :: ad_to1
    INTEGER :: ad_from2
    INTEGER :: ad_to2
    INTEGER :: ad_from3
    INTEGER :: ad_to3
    INTEGER :: ad_from4
    INTEGER :: ad_to4
    INTEGER :: branch
    REAL, DIMENSION(is:ie, js:je+1), INTENT(INOUT) :: fy_ad
    REAL :: temp_ad3
    REAL :: temp_ad2
    REAL :: temp_ad1
    REAL :: temp_ad0
    REAL :: temp_ad
    IF (nord .GT. 0) THEN
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHCONTROL1B(0)
    END IF
    IF (nord .GT. 0) THEN
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHCONTROL1B(0)
    END IF
    IF (nord .GT. 0) THEN
!----------
! high-order
!----------
      DO n=1,nord
        nt = nord - n
        ad_from0 = js - nt - 1
        DO j=ad_from0,je+nt+1
          ad_from = is - nt - 1
          i = ie + nt + 2
          CALL PUSHINTEGER4(i - 1)
          CALL PUSHINTEGER4(ad_from)
        END DO
        CALL PUSHINTEGER4(j - 1)
        CALL PUSHINTEGER4(ad_from0)
        ad_from2 = js - nt
        DO j=ad_from2,je+nt
          ad_from1 = is - nt
          i = ie + nt + 2
          CALL PUSHINTEGER4(i - 1)
          CALL PUSHINTEGER4(ad_from1)
        END DO
        CALL PUSHINTEGER4(j - 1)
        CALL PUSHINTEGER4(ad_from2)
        ad_from4 = js - nt
        DO j=ad_from4,je+nt+1
          ad_from3 = is - nt
          i = ie + nt + 1
          CALL PUSHINTEGER4(i - 1)
          CALL PUSHINTEGER4(ad_from3)
        END DO
        CALL PUSHINTEGER4(j - 1)
        CALL PUSHINTEGER4(ad_from4)
      END DO
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHCONTROL1B(0)
    END IF
    fy2_ad = 0.0_8
    DO j=je+1,js,-1
      DO i=ie,is,-1
        fy2_ad(i, j) = fy2_ad(i, j) + fy_ad(i, j)
      END DO
    END DO
    fx2_ad = 0.0_8
    DO j=je,js,-1
      DO i=ie+1,is,-1
        fx2_ad(i, j) = fx2_ad(i, j) + fx_ad(i, j)
      END DO
    END DO
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      d2_ad = 0.0_8
    ELSE
      d2_ad = 0.0_8
      DO n=nord,1,-1
        CALL POPINTEGER4(ad_from4)
        CALL POPINTEGER4(ad_to4)
        DO j=ad_to4,ad_from4,-1
          CALL POPINTEGER4(ad_from3)
          CALL POPINTEGER4(ad_to3)
          DO i=ad_to3,ad_from3,-1
            temp_ad3 = dx(i, j)*sina_v(i, j)*rdyc(i, j)*fy2_ad(i, j)
            d2_ad(i, j) = d2_ad(i, j) + temp_ad3
            d2_ad(i, j-1) = d2_ad(i, j-1) - temp_ad3
            fy2_ad(i, j) = 0.0_8
          END DO
        END DO
        CALL COPY_CORNERS_ADM(d2, d2_ad, npx, npy, 2)
        CALL POPINTEGER4(ad_from2)
        CALL POPINTEGER4(ad_to2)
        DO j=ad_to2,ad_from2,-1
          CALL POPINTEGER4(ad_from1)
          CALL POPINTEGER4(ad_to1)
          DO i=ad_to1,ad_from1,-1
            temp_ad2 = dy(i, j)*sina_u(i, j)*rdxc(i, j)*fx2_ad(i, j)
            d2_ad(i, j) = d2_ad(i, j) + temp_ad2
            d2_ad(i-1, j) = d2_ad(i-1, j) - temp_ad2
            fx2_ad(i, j) = 0.0_8
          END DO
        END DO
        CALL COPY_CORNERS_ADM(d2, d2_ad, npx, npy, 1)
        CALL POPINTEGER4(ad_from0)
        CALL POPINTEGER4(ad_to0)
        DO j=ad_to0,ad_from0,-1
          CALL POPINTEGER4(ad_from)
          CALL POPINTEGER4(ad_to)
          DO i=ad_to,ad_from,-1
            temp_ad1 = rarea(i, j)*d2_ad(i, j)
            fx2_ad(i, j) = fx2_ad(i, j) + temp_ad1
            fx2_ad(i+1, j) = fx2_ad(i+1, j) - temp_ad1
            fy2_ad(i, j) = fy2_ad(i, j) + temp_ad1
            fy2_ad(i, j+1) = fy2_ad(i, j+1) - temp_ad1
            d2_ad(i, j) = 0.0_8
          END DO
        END DO
      END DO
    END IF
    DO j=je+nord+1,js-nord,-1
      DO i=ie+nord,is-nord,-1
        temp_ad0 = dx(i, j)*sina_v(i, j)*rdyc(i, j)*fy2_ad(i, j)
        d2_ad(i, j-1) = d2_ad(i, j-1) + temp_ad0
        d2_ad(i, j) = d2_ad(i, j) - temp_ad0
        fy2_ad(i, j) = 0.0_8
      END DO
    END DO
    CALL POPCONTROL1B(branch)
    IF (branch .NE. 0) CALL COPY_CORNERS_ADM(d2, d2_ad, npx, npy, 2)
    DO j=je+nord,js-nord,-1
      DO i=ie+nord+1,is-nord,-1
        temp_ad = dy(i, j)*sina_u(i, j)*rdxc(i, j)*fx2_ad(i, j)
        d2_ad(i-1, j) = d2_ad(i-1, j) + temp_ad
        d2_ad(i, j) = d2_ad(i, j) - temp_ad
        fx2_ad(i, j) = 0.0_8
      END DO
    END DO
    CALL POPCONTROL1B(branch)
    IF (branch .NE. 0) CALL COPY_CORNERS_ADM(d2, d2_ad, npx, npy, 1)
    DO j=jed,jsd,-1
      DO i=ied,isd,-1
        q_ad(i, j) = q_ad(i, j) + damp*d2_ad(i, j)
        d2_ad(i, j) = 0.0_8
      END DO
    END DO
  END SUBROUTINE DELN_FLUX_ADM
!  Differentiation of copy_corners_r8 in reverse (adjoint) mode (with options r8):
!   gradient     of useful results: q
!   with respect to varying inputs: q
  SUBROUTINE COPY_CORNERS_R8_ADM(q, q_ad, npx, npy, dir)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx, npy, dir
    REAL*8, INTENT(INOUT) :: q(isd:ied, jsd:jed)
    REAL*8 :: q_ad(isd:ied, jsd:jed)
    INTEGER :: i, j
    REAL*8 :: tmp
    REAL*8 :: tmp0
    REAL*8 :: tmp1
    REAL*8 :: tmp2
    REAL*8 :: tmp3
    REAL*8 :: tmp4
    REAL*8 :: tmp5
    REAL*8 :: tmp6
    INTEGER :: branch
    REAL*8 :: tmp_ad6
    REAL*8 :: tmp_ad5
    REAL*8 :: tmp_ad4
    REAL*8 :: tmp_ad3
    REAL*8 :: tmp_ad2
    REAL*8 :: tmp_ad1
    REAL*8 :: tmp_ad0
    REAL*8 :: tmp_ad
    IF (dir .EQ. 1) THEN
! XDir:
      IF (sw_corner) THEN
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (se_corner) THEN
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (ne_corner) THEN
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (nw_corner) THEN
        DO j=npy+ng-1,npy,-1
          DO i=0,1-ng,-1
            tmp_ad2 = q_ad(i, j)
            q_ad(i, j) = 0.0_8
            q_ad(npy-j, i-1+npx) = q_ad(npy-j, i-1+npx) + tmp_ad2
          END DO
        END DO
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO j=npy+ng-1,npy,-1
          DO i=npx+ng-1,npx,-1
            tmp_ad1 = q_ad(i, j)
            q_ad(i, j) = 0.0_8
            q_ad(j, 2*npx-1-i) = q_ad(j, 2*npx-1-i) + tmp_ad1
          END DO
        END DO
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO j=0,1-ng,-1
          DO i=npx+ng-1,npx,-1
            tmp_ad0 = q_ad(i, j)
            q_ad(i, j) = 0.0_8
            q_ad(npy-j, i-npx+1) = q_ad(npy-j, i-npx+1) + tmp_ad0
          END DO
        END DO
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO j=0,1-ng,-1
          DO i=0,1-ng,-1
            tmp_ad = q_ad(i, j)
            q_ad(i, j) = 0.0_8
            q_ad(j, 1-i) = q_ad(j, 1-i) + tmp_ad
          END DO
        END DO
      END IF
    ELSE IF (dir .EQ. 2) THEN
! YDir:
      IF (sw_corner) THEN
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (se_corner) THEN
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (ne_corner) THEN
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (nw_corner) THEN
        DO j=npy+ng-1,npy,-1
          DO i=0,1-ng,-1
            tmp_ad6 = q_ad(i, j)
            q_ad(i, j) = 0.0_8
            q_ad(j+1-npx, npy-i) = q_ad(j+1-npx, npy-i) + tmp_ad6
          END DO
        END DO
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO j=npy+ng-1,npy,-1
          DO i=npx+ng-1,npx,-1
            tmp_ad5 = q_ad(i, j)
            q_ad(i, j) = 0.0_8
            q_ad(2*npy-1-j, i) = q_ad(2*npy-1-j, i) + tmp_ad5
          END DO
        END DO
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO j=0,1-ng,-1
          DO i=npx+ng-1,npx,-1
            tmp_ad4 = q_ad(i, j)
            q_ad(i, j) = 0.0_8
            q_ad(npy+j-1, npx-i) = q_ad(npy+j-1, npx-i) + tmp_ad4
          END DO
        END DO
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO j=0,1-ng,-1
          DO i=0,1-ng,-1
            tmp_ad3 = q_ad(i, j)
            q_ad(i, j) = 0.0_8
            q_ad(1-j, i) = q_ad(1-j, i) + tmp_ad3
          END DO
        END DO
      END IF
    END IF
  END SUBROUTINE COPY_CORNERS_R8_ADM
  SUBROUTINE DELN_FLUX(nord, npx, npy, damp, q, fx, fy, mfx, mfy)
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
! diffusive fluxes: 
    REAL, INTENT(IN) :: mfx(is:ie+1, js:je), mfy(is:ie, js:je+1)
    REAL, INTENT(INOUT) :: fx(is:ie+1, js:je), fy(is:ie, js:je+1)
! local:
    REAL :: fx2(isd:ied+1, jsd:jed), fy2(isd:ied, jsd:jed+1)
    REAL :: d2(isd:ied, jsd:jed)
    INTEGER :: i, j, n, nt
    DO j=jsd,jed
      DO i=isd,ied
        d2(i, j) = damp*q(i, j)
      END DO
    END DO
    IF (nord .GT. 0) CALL COPY_CORNERS(d2, npx, npy, 1)
    DO j=js-nord,je+nord
      DO i=is-nord,ie+nord+1
        fx2(i, j) = dy(i, j)*sina_u(i, j)*(d2(i-1, j)-d2(i, j))*rdxc(i, &
&         j)
      END DO
    END DO
    IF (nord .GT. 0) CALL COPY_CORNERS(d2, npx, npy, 2)
    DO j=js-nord,je+nord+1
      DO i=is-nord,ie+nord
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
            d2(i, j) = (fx2(i, j)-fx2(i+1, j)+fy2(i, j)-fy2(i, j+1))*&
&             rarea(i, j)
          END DO
        END DO
        CALL COPY_CORNERS(d2, npx, npy, 1)
        DO j=js-nt,je+nt
          DO i=is-nt,ie+nt+1
            fx2(i, j) = dy(i, j)*sina_u(i, j)*(d2(i, j)-d2(i-1, j))*rdxc&
&             (i, j)
          END DO
        END DO
        CALL COPY_CORNERS(d2, npx, npy, 2)
        DO j=js-nt,je+nt+1
          DO i=is-nt,ie+nt
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
        fx(i, j) = fx(i, j) + fx2(i, j)
      END DO
    END DO
    DO j=js,je+1
      DO i=is,ie
!        fy(i,j) = fy(i,j) + fy2(i,j)*mfy(i,j)
        fy(i, j) = fy(i, j) + fy2(i, j)
      END DO
    END DO
  END SUBROUTINE DELN_FLUX
  SUBROUTINE COPY_CORNERS_R8(q, npx, npy, dir)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx, npy, dir
    REAL*8, INTENT(INOUT) :: q(isd:ied, jsd:jed)
    INTEGER :: i, j
    IF (dir .EQ. 1) THEN
! XDir:
      IF (sw_corner) THEN
        DO j=1-ng,0
          DO i=1-ng,0
            q(i, j) = q(j, 1-i)
          END DO
        END DO
      END IF
      IF (se_corner) THEN
        DO j=1-ng,0
          DO i=npx,npx+ng-1
            q(i, j) = q(npy-j, i-npx+1)
          END DO
        END DO
      END IF
      IF (ne_corner) THEN
        DO j=npy,npy+ng-1
          DO i=npx,npx+ng-1
            q(i, j) = q(j, 2*npx-1-i)
          END DO
        END DO
      END IF
      IF (nw_corner) THEN
        DO j=npy,npy+ng-1
          DO i=1-ng,0
            q(i, j) = q(npy-j, i-1+npx)
          END DO
        END DO
      END IF
    ELSE IF (dir .EQ. 2) THEN
! YDir:
      IF (sw_corner) THEN
        DO j=1-ng,0
          DO i=1-ng,0
            q(i, j) = q(1-j, i)
          END DO
        END DO
      END IF
      IF (se_corner) THEN
        DO j=1-ng,0
          DO i=npx,npx+ng-1
            q(i, j) = q(npy+j-1, npx-i)
          END DO
        END DO
      END IF
      IF (ne_corner) THEN
        DO j=npy,npy+ng-1
          DO i=npx,npx+ng-1
            q(i, j) = q(2*npy-1-j, i)
          END DO
        END DO
      END IF
      IF (nw_corner) THEN
        DO j=npy,npy+ng-1
          DO i=1-ng,0
            q(i, j) = q(j+1-npx, npy-i)
          END DO
        END DO
      END IF
    END IF
  END SUBROUTINE COPY_CORNERS_R8

end module tp_core_adm_mod
