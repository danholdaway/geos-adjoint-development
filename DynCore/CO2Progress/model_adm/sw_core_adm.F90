 module sw_core_adm_mod

 use fv_arrays_mod,     only: p_precision, i_precision, ke_precision
 use fv_control_mod,    only: shallow_water
 use fv_mp_mod,         only: ng, is,js,ie,je, isd,jsd,ied,jed,  &
                              mp_corner_comm, domain,            &
                              fill_corners, XDir, YDir
 use sw_core_mod,       only: xtp_u, ytp_v
 use fv_grid_tools_mod, only: npx=>npx_g,npy=>npy_g, cosa, sina,  &
                              rdxc, rdyc, dx,dy, dxc,dyc, dxa,dya,  &
                              rdxa, rdya, area, area_c, rarea, rarea_c, rdx, rdy, &
                              grid_type
 use tp_core_mod,       only: fv_tp_2d
 use tp_core_adm_mod,   only: copy_corners, fv_tp_2d_adm, copy_corners_adm
 use fv_grid_utils_mod, only: sw_corner, se_corner, ne_corner, nw_corner,       &
                              cosa_u, cosa_v, cosa_s, sina_s, sina_u, sina_v,   &
                              rsin_u, rsin_v, rsin_v, rsina, ec1, ec2, ew, es,  &
                              big_number, da_min_c, da_min, fC, f0,   &
                              rsin2, divg_u, divg_v, Gnomonic_grid
 use test_cases_mod,    only: test_case

 implicit none

  real, parameter:: r3 =   1./3.
  real, parameter:: t11=27./28., t12=-13./28., t13=3./7., t14=6./7., t15=3./28.
  real, parameter:: s11=11./14., s13=-13./14., s14=4./7., s15=3./14.
!----------------------
! PPM volume mean form:
!----------------------
  real, parameter:: p1 =  7./12.     ! 0.58333333
  real, parameter:: p2 = -1./12.
!----------------------------
! 4-pt Lagrange interpolation
!----------------------------
  real(i_precision), parameter:: a1 =  0.5625
  real(i_precision), parameter:: a2 = -0.0625
!----------------------------------------------
! volume-conserving cubic with 2nd drv=0 at end point:
  real, parameter:: c1 = -2./14.
  real, parameter:: c2 = 11./14.
  real, parameter:: c3 =  5./14.
! 3-pt off-center intp formular:
! real, parameter:: c1 = -0.125
! real, parameter:: c2 =  0.75
! real, parameter:: c3 =  0.375
!----------------------------------------------
! scheme 2.1: perturbation form
  real, parameter:: b1 =   1./30.
  real, parameter:: b2 = -13./60.
  real, parameter:: b3 = -13./60.
  real, parameter:: b4 =  0.45
  real, parameter:: b5 = -0.05

  private
  public :: c_sw_adm, d_sw_adm, d2a2c_adm, d2a2c_vect_adm, divergence_corner_adm, &
            c_sw, d_sw, d2a2c, d2a2c_vect, divergence_corner
  contains

  SUBROUTINE C_SW_ADM(delpc, delpc_ad, delp, delp_ad, ptc, ptc_ad, pt, &
&   pt_ad, u, u_ad, v, v_ad, w, w_ad, uc, uc_ad, vc, vc_ad, ua, ua_ad, &
&   va, va_ad, wc, ut, ut_ad, vt, vt_ad, dt2, hydrostatic, dord4)
    IMPLICIT NONE
!#endif
    REAL, DIMENSION(isd:ied, jsd:jed + 1), INTENT(INOUT) :: u, vc
    REAL, DIMENSION(isd:ied, jsd:jed+1), INTENT(INOUT) :: u_ad
    REAL, DIMENSION(isd:ied + 1, jsd:jed), INTENT(INOUT) :: v, uc
    REAL, DIMENSION(isd:ied+1, jsd:jed), INTENT(INOUT) :: v_ad
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(INOUT) :: delp, pt, ua, va&
&   , w
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(INOUT) :: delp_ad
    REAL, DIMENSION(isd:ied, jsd:jed) :: delpc, ptc, ut, vt, wc
    REAL, DIMENSION(isd:ied, jsd:jed) :: delpc_ad, ptc_ad, ut_ad, vt_ad
    REAL, INTENT(IN) :: dt2
    LOGICAL, INTENT(IN) :: hydrostatic
    LOGICAL, INTENT(IN) :: dord4
! Local:
    REAL, DIMENSION(is - 1:ie + 1, js - 1:je + 1) :: vort, ke
    REAL, DIMENSION(is-1:ie+1, js-1:je+1) :: vort_ad, ke_ad
    REAL, DIMENSION(is - 1:ie + 2, js - 1:je + 1) :: fx, fx1, fx2
    REAL, DIMENSION(is-1:ie+2, js-1:je+1) :: fx_ad, fx1_ad
    REAL, DIMENSION(is - 1:ie + 1, js - 1:je + 2) :: fy, fy1, fy2
    REAL, DIMENSION(is-1:ie+1, js-1:je+2) :: fy_ad, fy1_ad
    REAL :: dt4
    INTEGER :: i, j, is2, ie1
    INTEGER :: iep1, jep1
    INTRINSIC MAX
    INTRINSIC MIN
    INTEGER :: branch
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(INOUT) :: va_ad
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(INOUT) :: pt_ad
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(INOUT) :: w_ad
    REAL, DIMENSION(isd:ied+1, jsd:jed), INTENT(INOUT) :: uc_ad
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(INOUT) :: ua_ad
    REAL, DIMENSION(isd:ied, jsd:jed+1), INTENT(INOUT) :: vc_ad
    REAL :: temp_ad9
    REAL :: temp_ad8
    REAL :: temp_ad7
    REAL :: temp_ad6
    REAL :: temp_ad5
    REAL :: temp_ad4
    REAL :: temp_ad3
    REAL :: temp_ad2
    REAL :: temp_ad1
    REAL :: temp_ad0
    REAL :: temp_ad
    REAL :: temp_ad10
!#ifdef DOUBLE_GRADIENTS_CSW
!      real*8 :: delp_dp, pt_dp, dp_r8
!      real*8, dimension(is-1:ie+1,js-1:je+1):: ke_dp
!      real*8, dimension(is-1:ie+2,js-1:je+1):: fx_dp, fx1_dp
!      real*8, dimension(is-1:ie+1,js-1:je+2):: fy_dp, fy1_dp
!#endif
!#ifdef FIX_C_BOUNDARY 
!      iep1 = ie;   jep1 = je
!#else
    iep1 = ie + 1
    jep1 = je + 1
!#endif
    CALL PUSHREAL8ARRAY(va, (ied-isd+1)*(jed-jsd+1))
    CALL PUSHREAL8ARRAY(ua, (ied-isd+1)*(jed-jsd+1))
    CALL D2A2C_VECT(u, v, ua, va, uc, vc, ut, vt, dord4)
!     call d2a2c_vect_v2(u, v, ua, va, uc, vc, ut, vt)
    DO j=js-1,jep1
      DO i=is-1,iep1+1
        ut(i, j) = dt2*ut(i, j)*dy(i, j)*sina_u(i, j)
      END DO
    END DO
    DO j=js-1,jep1+1
      DO i=is-1,iep1
        vt(i, j) = dt2*vt(i, j)*dx(i, j)*sina_v(i, j)
      END DO
    END DO
!----------------
! Transport delp:
!----------------
! Xdir:
    IF (grid_type .LT. 3) THEN
      CALL PUSHREAL8ARRAY(pt, (ied-isd+1)*(jed-jsd+1))
      CALL FILL2_4CORNERS(delp, pt, 1)
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    IF (hydrostatic) THEN
      IF (shallow_water) THEN
        DO j=js-1,jep1
          DO i=is-1,iep1+1
            IF (ut(i, j) .GT. 0.) THEN
              fx1(i, j) = delp(i-1, j)
              CALL PUSHCONTROL1B(0)
            ELSE
              fx1(i, j) = delp(i, j)
              CALL PUSHCONTROL1B(1)
            END IF
            CALL PUSHREAL8(fx1(i, j))
            fx1(i, j) = ut(i, j)*fx1(i, j)
          END DO
        END DO
        CALL PUSHCONTROL2B(0)
      ELSE
        DO j=js-1,jep1
          DO i=is-1,iep1+1
            IF (ut(i, j) .GT. 0.) THEN
              fx1(i, j) = delp(i-1, j)
              fx(i, j) = pt(i-1, j)
              CALL PUSHCONTROL1B(0)
            ELSE
              fx1(i, j) = delp(i, j)
              fx(i, j) = pt(i, j)
              CALL PUSHCONTROL1B(1)
            END IF
            CALL PUSHREAL8(fx1(i, j))
            fx1(i, j) = ut(i, j)*fx1(i, j)
            CALL PUSHREAL8(fx(i, j))
            fx(i, j) = fx1(i, j)*fx(i, j)
          END DO
        END DO
        CALL PUSHCONTROL2B(1)
      END IF
    ELSE
      IF (grid_type .LT. 3) THEN
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL PUSHCONTROL1B(0)
      END IF
      DO j=js-1,je+1
        DO i=is-1,ie+2
          IF (ut(i, j) .GT. 0.) THEN
            fx1(i, j) = delp(i-1, j)
            fx(i, j) = pt(i-1, j)
            CALL PUSHCONTROL1B(0)
          ELSE
            fx1(i, j) = delp(i, j)
            fx(i, j) = pt(i, j)
            CALL PUSHCONTROL1B(1)
          END IF
          CALL PUSHREAL8(fx1(i, j))
          fx1(i, j) = ut(i, j)*fx1(i, j)
          CALL PUSHREAL8(fx(i, j))
          fx(i, j) = fx1(i, j)*fx(i, j)
        END DO
      END DO
      CALL PUSHCONTROL2B(2)
    END IF
! Ydir:
    IF (grid_type .LT. 3) THEN
      CALL PUSHREAL8ARRAY(pt, (ied-isd+1)*(jed-jsd+1))
      CALL FILL2_4CORNERS(delp, pt, 2)
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    IF (hydrostatic) THEN
      DO j=js-1,jep1+1
        DO i=is-1,iep1
          IF (vt(i, j) .GT. 0.) THEN
            fy1(i, j) = delp(i, j-1)
            fy(i, j) = pt(i, j-1)
            CALL PUSHCONTROL1B(0)
          ELSE
            fy1(i, j) = delp(i, j)
            fy(i, j) = pt(i, j)
            CALL PUSHCONTROL1B(1)
          END IF
          CALL PUSHREAL8(fy1(i, j))
          fy1(i, j) = vt(i, j)*fy1(i, j)
          CALL PUSHREAL8(fy(i, j))
          fy(i, j) = fy1(i, j)*fy(i, j)
        END DO
      END DO
      IF (shallow_water) THEN
        CALL PUSHCONTROL2B(2)
      ELSE
!#ifdef DOUBLE_GRADIENTS_CSW
!           fy1_dp = fy1
!           fy_dp  = fy
!           fx1_dp = fx1
!           fx_dp  = fx
!           do j=js-1,jep1
!              do i=is-1,iep1
!                 delp_dp = delp(i,j)
!                 pt_dp   = pt(i,j)
!                 delpc(i,j) =        delp_dp + (fx1_dp(i,j)-fx1_dp(i+1,j)+fy1_dp(i,j)-fy1_dp(i,j+1))*rarea(i,j)
!                   ptc(i,j) = (pt_dp*delp_dp + ( fx_dp(i,j)- fx_dp(i+1,j)+ fy_dp(i,j)- fy_dp(i,j+1))*rarea(i,j))/delpc(i,j)
!              enddo
!           enddo
!#else
        DO j=js-1,jep1
          DO i=is-1,iep1
            CALL PUSHREAL8(delpc(i, j))
            delpc(i, j) = delp(i, j) + (fx1(i, j)-fx1(i+1, j)+fy1(i, j)-&
&             fy1(i, j+1))*rarea(i, j)
          END DO
        END DO
        CALL PUSHCONTROL2B(1)
      END IF
    ELSE
!#endif
      IF (grid_type .LT. 3) THEN
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL PUSHCONTROL1B(0)
      END IF
      DO j=js-1,je+2
        DO i=is-1,ie+1
          IF (vt(i, j) .GT. 0.) THEN
            fy1(i, j) = delp(i, j-1)
            fy(i, j) = pt(i, j-1)
            CALL PUSHCONTROL1B(0)
          ELSE
            fy1(i, j) = delp(i, j)
            fy(i, j) = pt(i, j)
            CALL PUSHCONTROL1B(1)
          END IF
          CALL PUSHREAL8(fy1(i, j))
          fy1(i, j) = vt(i, j)*fy1(i, j)
          CALL PUSHREAL8(fy(i, j))
          fy(i, j) = fy1(i, j)*fy(i, j)
        END DO
      END DO
      DO j=js-1,je+1
        DO i=is-1,ie+1
          CALL PUSHREAL8(delpc(i, j))
          delpc(i, j) = delp(i, j) + (fx1(i, j)-fx1(i+1, j)+fy1(i, j)-&
&           fy1(i, j+1))*rarea(i, j)
        END DO
      END DO
      CALL PUSHCONTROL2B(0)
    END IF
!------------
! Compute KE:
!------------
    DO j=js-1,jep1
      DO i=is-1,iep1
        IF (ua(i, j) .GT. 0.) THEN
          IF (i .EQ. 1) THEN
            ke(1, j) = uc(1, j)*sina_u(1, j) + v(1, j)*cosa_u(1, j)
            CALL PUSHCONTROL3B(5)
          ELSE IF (i .EQ. npx) THEN
            ke(i, j) = uc(npx, j)*sina_u(npx, j) - v(npx, j)*cosa_u(npx&
&             , j)
            CALL PUSHCONTROL3B(4)
          ELSE
            ke(i, j) = uc(i, j)
            CALL PUSHCONTROL3B(3)
          END IF
        ELSE IF (i .EQ. 0) THEN
          ke(0, j) = uc(1, j)*sina_u(1, j) - v(1, j)*cosa_u(1, j)
          CALL PUSHCONTROL3B(2)
        ELSE IF (i .EQ. npx - 1) THEN
          ke(i, j) = uc(npx, j)*sina_u(npx, j) + v(npx, j)*cosa_u(npx, j&
&           )
          CALL PUSHCONTROL3B(1)
        ELSE
          ke(i, j) = uc(i+1, j)
          CALL PUSHCONTROL3B(0)
        END IF
      END DO
    END DO
    DO j=js-1,jep1
      DO i=is-1,iep1
        IF (va(i, j) .GT. 0.) THEN
          IF (j .EQ. 1) THEN
            vort(i, 1) = vc(i, 1)*sina_v(i, 1) + u(i, 1)*cosa_v(i, 1)
            CALL PUSHCONTROL3B(5)
          ELSE IF (j .EQ. npy) THEN
            vort(i, j) = vc(i, npy)*sina_v(i, npy) - u(i, npy)*cosa_v(i&
&             , npy)
            CALL PUSHCONTROL3B(4)
          ELSE
            vort(i, j) = vc(i, j)
            CALL PUSHCONTROL3B(3)
          END IF
        ELSE IF (j .EQ. 0) THEN
          vort(i, 0) = vc(i, 1)*sina_v(i, 1) - u(i, 1)*cosa_v(i, 1)
          CALL PUSHCONTROL3B(2)
        ELSE IF (j .EQ. npy - 1) THEN
          vort(i, j) = vc(i, npy)*sina_v(i, npy) + u(i, npy)*cosa_v(i, &
&           npy)
          CALL PUSHCONTROL3B(1)
        ELSE
          vort(i, j) = vc(i, j+1)
          CALL PUSHCONTROL3B(0)
        END IF
      END DO
    END DO
    dt4 = 0.5*dt2
    IF (2 .LT. is) THEN
      is2 = is
    ELSE
      is2 = 2
    END IF
    IF (npx - 1 .GT. ie + 1) THEN
      ie1 = ie + 1
    ELSE
      ie1 = npx - 1
    END IF
    DO j=js-1,je+1
      DO i=is2,ie1
        CALL PUSHREAL8(fx(i, j))
        fx(i, j) = uc(i, j)*dxc(i, j)
      END DO
      IF (is .EQ. 1) THEN
        CALL PUSHREAL8(fx(1, j))
        fx(1, j) = uc(1, j)*sina_u(1, j)*dxc(1, j)
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (ie + 1 .EQ. npx) THEN
        CALL PUSHREAL8(fx(npx, j))
        fx(npx, j) = uc(npx, j)*sina_u(npx, j)*dxc(npx, j)
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL PUSHCONTROL1B(0)
      END IF
    END DO
    DO j=js,je+1
      IF (j .EQ. 1 .OR. j .EQ. npy) THEN
        DO i=is-1,ie+1
          CALL PUSHREAL8(fy(i, j))
          fy(i, j) = vc(i, j)*sina_v(i, j)*dyc(i, j)
        END DO
        CALL PUSHCONTROL1B(1)
      ELSE
        DO i=is-1,ie+1
          CALL PUSHREAL8(fy(i, j))
          fy(i, j) = vc(i, j)*dyc(i, j)
        END DO
        CALL PUSHCONTROL1B(0)
      END IF
    END DO
!#endif
    DO j=js,je+1
      DO i=is,ie+1
        CALL PUSHREAL8(vort(i, j))
        vort(i, j) = fx(i, j-1) - fx(i, j) - fy(i-1, j) + fy(i, j)
      END DO
    END DO
! Remove the extra term at the corners:
    IF (sw_corner) THEN
      CALL PUSHREAL8(vort(1, 1))
      vort(1, 1) = vort(1, 1) + fy(0, 1)
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    IF (se_corner) THEN
      CALL PUSHREAL8(vort(npx, 1))
      vort(npx, 1) = vort(npx, 1) - fy(npx, 1)
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    IF (ne_corner) THEN
      CALL PUSHREAL8(vort(npx, npy))
      vort(npx, npy) = vort(npx, npy) - fy(npx, npy)
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    IF (nw_corner) THEN
      CALL PUSHREAL8(vort(1, npy))
      vort(1, npy) = vort(1, npy) + fy(0, npy)
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHCONTROL1B(0)
    END IF
!----------------------------
! Compute absolute vorticity
!----------------------------
    DO j=js,je+1
      DO i=is,ie+1
        CALL PUSHREAL8(vort(i, j))
        vort(i, j) = fc(i, j) + rarea_c(i, j)*vort(i, j)
      END DO
    END DO
!----------------------------------
! Transport absolute vorticity:
!----------------------------------
    DO j=js,je
      DO i=is,iep1
        IF (i .EQ. 1 .OR. i .EQ. npx) THEN
          CALL PUSHREAL8(fy1(i, j))
          fy1(i, j) = dt2*v(i, j)*sina_u(i, j)
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHREAL8(fy1(i, j))
          fy1(i, j) = dt2*(v(i, j)-uc(i, j)*cosa_u(i, j))/sina_u(i, j)
          CALL PUSHCONTROL1B(1)
        END IF
        IF (fy1(i, j) .GT. 0.) THEN
          CALL PUSHREAL8(fy(i, j))
          fy(i, j) = vort(i, j)
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHREAL8(fy(i, j))
          fy(i, j) = vort(i, j+1)
          CALL PUSHCONTROL1B(0)
        END IF
      END DO
    END DO
    DO j=js,jep1
      IF (j .EQ. 1 .OR. j .EQ. npy) THEN
        DO i=is,ie
          CALL PUSHREAL8(fx1(i, j))
          fx1(i, j) = dt2*u(i, j)*sina_v(i, j)
          IF (fx1(i, j) .GT. 0.) THEN
            CALL PUSHREAL8(fx(i, j))
            fx(i, j) = vort(i, j)
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHREAL8(fx(i, j))
            fx(i, j) = vort(i+1, j)
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
        CALL PUSHCONTROL1B(1)
      ELSE
        DO i=is,ie
          CALL PUSHREAL8(fx1(i, j))
          fx1(i, j) = dt2*(u(i, j)-vc(i, j)*cosa_v(i, j))/sina_v(i, j)
          IF (fx1(i, j) .GT. 0.) THEN
            CALL PUSHREAL8(fx(i, j))
            fx(i, j) = vort(i, j)
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHREAL8(fx(i, j))
            fx(i, j) = vort(i+1, j)
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
        CALL PUSHCONTROL1B(0)
      END IF
    END DO
    ke_ad = 0.0_8
    fx_ad = 0.0_8
    fx1_ad = 0.0_8
    DO j=jep1,js,-1
      DO i=ie,is,-1
        temp_ad10 = rdyc(i, j)*vc_ad(i, j)
        fx1_ad(i, j) = fx1_ad(i, j) - fx(i, j)*vc_ad(i, j)
        fx_ad(i, j) = fx_ad(i, j) - fx1(i, j)*vc_ad(i, j)
        ke_ad(i, j-1) = ke_ad(i, j-1) + temp_ad10
        ke_ad(i, j) = ke_ad(i, j) - temp_ad10
      END DO
    END DO
    fy1_ad = 0.0_8
    fy_ad = 0.0_8
    DO j=je,js,-1
      DO i=iep1,is,-1
        temp_ad9 = rdxc(i, j)*uc_ad(i, j)
        fy1_ad(i, j) = fy1_ad(i, j) + fy(i, j)*uc_ad(i, j)
        fy_ad(i, j) = fy_ad(i, j) + fy1(i, j)*uc_ad(i, j)
        ke_ad(i-1, j) = ke_ad(i-1, j) + temp_ad9
        ke_ad(i, j) = ke_ad(i, j) - temp_ad9
      END DO
    END DO
    vort_ad = 0.0_8
    DO j=jep1,js,-1
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO i=ie,is,-1
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL8(fx(i, j))
            vort_ad(i+1, j) = vort_ad(i+1, j) + fx_ad(i, j)
            fx_ad(i, j) = 0.0_8
          ELSE
            CALL POPREAL8(fx(i, j))
            vort_ad(i, j) = vort_ad(i, j) + fx_ad(i, j)
            fx_ad(i, j) = 0.0_8
          END IF
          CALL POPREAL8(fx1(i, j))
          temp_ad8 = dt2*fx1_ad(i, j)/sina_v(i, j)
          u_ad(i, j) = u_ad(i, j) + temp_ad8
          vc_ad(i, j) = vc_ad(i, j) - cosa_v(i, j)*temp_ad8
          fx1_ad(i, j) = 0.0_8
        END DO
      ELSE
        DO i=ie,is,-1
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL8(fx(i, j))
            vort_ad(i+1, j) = vort_ad(i+1, j) + fx_ad(i, j)
            fx_ad(i, j) = 0.0_8
          ELSE
            CALL POPREAL8(fx(i, j))
            vort_ad(i, j) = vort_ad(i, j) + fx_ad(i, j)
            fx_ad(i, j) = 0.0_8
          END IF
          CALL POPREAL8(fx1(i, j))
          u_ad(i, j) = u_ad(i, j) + dt2*sina_v(i, j)*fx1_ad(i, j)
          fx1_ad(i, j) = 0.0_8
        END DO
      END IF
    END DO
    DO j=je,js,-1
      DO i=iep1,is,-1
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          CALL POPREAL8(fy(i, j))
          vort_ad(i, j+1) = vort_ad(i, j+1) + fy_ad(i, j)
          fy_ad(i, j) = 0.0_8
        ELSE
          CALL POPREAL8(fy(i, j))
          vort_ad(i, j) = vort_ad(i, j) + fy_ad(i, j)
          fy_ad(i, j) = 0.0_8
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          CALL POPREAL8(fy1(i, j))
          v_ad(i, j) = v_ad(i, j) + dt2*sina_u(i, j)*fy1_ad(i, j)
          fy1_ad(i, j) = 0.0_8
        ELSE
          CALL POPREAL8(fy1(i, j))
          temp_ad7 = dt2*fy1_ad(i, j)/sina_u(i, j)
          v_ad(i, j) = v_ad(i, j) + temp_ad7
          uc_ad(i, j) = uc_ad(i, j) - cosa_u(i, j)*temp_ad7
          fy1_ad(i, j) = 0.0_8
        END IF
      END DO
    END DO
    DO j=je+1,js,-1
      DO i=ie+1,is,-1
        CALL POPREAL8(vort(i, j))
        vort_ad(i, j) = rarea_c(i, j)*vort_ad(i, j)
      END DO
    END DO
    CALL POPCONTROL1B(branch)
    IF (branch .NE. 0) THEN
      CALL POPREAL8(vort(1, npy))
      fy_ad(0, npy) = fy_ad(0, npy) + vort_ad(1, npy)
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL8(vort(npx, npy))
      fy_ad(npx, npy) = fy_ad(npx, npy) - vort_ad(npx, npy)
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL8(vort(npx, 1))
      fy_ad(npx, 1) = fy_ad(npx, 1) - vort_ad(npx, 1)
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL8(vort(1, 1))
      fy_ad(0, 1) = fy_ad(0, 1) + vort_ad(1, 1)
    END IF
    DO j=je+1,js,-1
      DO i=ie+1,is,-1
        CALL POPREAL8(vort(i, j))
        fx_ad(i, j-1) = fx_ad(i, j-1) + vort_ad(i, j)
        fx_ad(i, j) = fx_ad(i, j) - vort_ad(i, j)
        fy_ad(i, j) = fy_ad(i, j) + vort_ad(i, j)
        fy_ad(i-1, j) = fy_ad(i-1, j) - vort_ad(i, j)
        vort_ad(i, j) = 0.0_8
      END DO
    END DO
    DO j=je+1,js,-1
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO i=ie+1,is-1,-1
          CALL POPREAL8(fy(i, j))
          vc_ad(i, j) = vc_ad(i, j) + dyc(i, j)*fy_ad(i, j)
          fy_ad(i, j) = 0.0_8
        END DO
      ELSE
        DO i=ie+1,is-1,-1
          CALL POPREAL8(fy(i, j))
          vc_ad(i, j) = vc_ad(i, j) + sina_v(i, j)*dyc(i, j)*fy_ad(i, j)
          fy_ad(i, j) = 0.0_8
        END DO
      END IF
    END DO
    DO j=je+1,js-1,-1
      CALL POPCONTROL1B(branch)
      IF (branch .NE. 0) THEN
        CALL POPREAL8(fx(npx, j))
        uc_ad(npx, j) = uc_ad(npx, j) + sina_u(npx, j)*dxc(npx, j)*fx_ad&
&         (npx, j)
        fx_ad(npx, j) = 0.0_8
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        CALL POPREAL8(fx(1, j))
        uc_ad(1, j) = uc_ad(1, j) + sina_u(1, j)*dxc(1, j)*fx_ad(1, j)
        fx_ad(1, j) = 0.0_8
      END IF
      DO i=ie1,is2,-1
        CALL POPREAL8(fx(i, j))
        uc_ad(i, j) = uc_ad(i, j) + dxc(i, j)*fx_ad(i, j)
        fx_ad(i, j) = 0.0_8
      END DO
    END DO
    DO j=jep1,js-1,-1
      DO i=iep1,is-1,-1
        temp_ad6 = dt4*ke_ad(i, j)
        ua_ad(i, j) = ua_ad(i, j) + ke(i, j)*temp_ad6
        va_ad(i, j) = va_ad(i, j) + vort(i, j)*temp_ad6
        vort_ad(i, j) = vort_ad(i, j) + va(i, j)*temp_ad6
        ke_ad(i, j) = ua(i, j)*temp_ad6
      END DO
    END DO
    DO j=jep1,js-1,-1
      DO i=iep1,is-1,-1
        CALL POPCONTROL3B(branch)
        IF (branch .LT. 3) THEN
          IF (branch .EQ. 0) THEN
            vc_ad(i, j+1) = vc_ad(i, j+1) + vort_ad(i, j)
            vort_ad(i, j) = 0.0_8
          ELSE IF (branch .EQ. 1) THEN
            vc_ad(i, npy) = vc_ad(i, npy) + sina_v(i, npy)*vort_ad(i, j)
            u_ad(i, npy) = u_ad(i, npy) + cosa_v(i, npy)*vort_ad(i, j)
            vort_ad(i, j) = 0.0_8
          ELSE
            vc_ad(i, 1) = vc_ad(i, 1) + sina_v(i, 1)*vort_ad(i, 0)
            u_ad(i, 1) = u_ad(i, 1) - cosa_v(i, 1)*vort_ad(i, 0)
            vort_ad(i, 0) = 0.0_8
          END IF
        ELSE IF (branch .EQ. 3) THEN
          vc_ad(i, j) = vc_ad(i, j) + vort_ad(i, j)
          vort_ad(i, j) = 0.0_8
        ELSE IF (branch .EQ. 4) THEN
          vc_ad(i, npy) = vc_ad(i, npy) + sina_v(i, npy)*vort_ad(i, j)
          u_ad(i, npy) = u_ad(i, npy) - cosa_v(i, npy)*vort_ad(i, j)
          vort_ad(i, j) = 0.0_8
        ELSE
          vc_ad(i, 1) = vc_ad(i, 1) + sina_v(i, 1)*vort_ad(i, 1)
          u_ad(i, 1) = u_ad(i, 1) + cosa_v(i, 1)*vort_ad(i, 1)
          vort_ad(i, 1) = 0.0_8
        END IF
      END DO
    END DO
    DO j=jep1,js-1,-1
      DO i=iep1,is-1,-1
        CALL POPCONTROL3B(branch)
        IF (branch .LT. 3) THEN
          IF (branch .EQ. 0) THEN
            uc_ad(i+1, j) = uc_ad(i+1, j) + ke_ad(i, j)
            ke_ad(i, j) = 0.0_8
          ELSE IF (branch .EQ. 1) THEN
            uc_ad(npx, j) = uc_ad(npx, j) + sina_u(npx, j)*ke_ad(i, j)
            v_ad(npx, j) = v_ad(npx, j) + cosa_u(npx, j)*ke_ad(i, j)
            ke_ad(i, j) = 0.0_8
          ELSE
            uc_ad(1, j) = uc_ad(1, j) + sina_u(1, j)*ke_ad(0, j)
            v_ad(1, j) = v_ad(1, j) - cosa_u(1, j)*ke_ad(0, j)
            ke_ad(0, j) = 0.0_8
          END IF
        ELSE IF (branch .EQ. 3) THEN
          uc_ad(i, j) = uc_ad(i, j) + ke_ad(i, j)
          ke_ad(i, j) = 0.0_8
        ELSE IF (branch .EQ. 4) THEN
          uc_ad(npx, j) = uc_ad(npx, j) + sina_u(npx, j)*ke_ad(i, j)
          v_ad(npx, j) = v_ad(npx, j) - cosa_u(npx, j)*ke_ad(i, j)
          ke_ad(i, j) = 0.0_8
        ELSE
          uc_ad(1, j) = uc_ad(1, j) + sina_u(1, j)*ke_ad(1, j)
          v_ad(1, j) = v_ad(1, j) + cosa_u(1, j)*ke_ad(1, j)
          ke_ad(1, j) = 0.0_8
        END IF
      END DO
    END DO
    CALL POPCONTROL2B(branch)
    IF (branch .EQ. 0) THEN
      DO j=je+1,js-1,-1
        DO i=ie+1,is-1,-1
          temp_ad3 = ptc_ad(i, j)/delpc(i, j)
          temp_ad4 = rarea(i, j)*temp_ad3
          pt_ad(i, j) = pt_ad(i, j) + delp(i, j)*temp_ad3
          fx_ad(i, j) = fx_ad(i, j) + temp_ad4
          fx_ad(i+1, j) = fx_ad(i+1, j) - temp_ad4
          fy_ad(i, j) = fy_ad(i, j) + temp_ad4
          fy_ad(i, j+1) = fy_ad(i, j+1) - temp_ad4
          delpc_ad(i, j) = delpc_ad(i, j) - (pt(i, j)*delp(i, j)+rarea(i&
&           , j)*(fx(i, j)-fx(i+1, j)+fy(i, j)-fy(i, j+1)))*temp_ad3/&
&           delpc(i, j)
          delp_ad(i, j) = delp_ad(i, j) + delpc_ad(i, j) + pt(i, j)*&
&           temp_ad3
          ptc_ad(i, j) = 0.0_8
          CALL POPREAL8(delpc(i, j))
          temp_ad5 = rarea(i, j)*delpc_ad(i, j)
          fx1_ad(i, j) = fx1_ad(i, j) + temp_ad5
          fx1_ad(i+1, j) = fx1_ad(i+1, j) - temp_ad5
          fy1_ad(i, j) = fy1_ad(i, j) + temp_ad5
          fy1_ad(i, j+1) = fy1_ad(i, j+1) - temp_ad5
          delpc_ad(i, j) = 0.0_8
        END DO
      END DO
      DO j=je+2,js-1,-1
        DO i=ie+1,is-1,-1
          CALL POPREAL8(fy(i, j))
          fy1_ad(i, j) = fy1_ad(i, j) + fy(i, j)*fy_ad(i, j)
          fy_ad(i, j) = fy1(i, j)*fy_ad(i, j)
          CALL POPREAL8(fy1(i, j))
          vt_ad(i, j) = vt_ad(i, j) + fy1(i, j)*fy1_ad(i, j)
          fy1_ad(i, j) = vt(i, j)*fy1_ad(i, j)
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            pt_ad(i, j-1) = pt_ad(i, j-1) + fy_ad(i, j)
            fy_ad(i, j) = 0.0_8
            delp_ad(i, j-1) = delp_ad(i, j-1) + fy1_ad(i, j)
            fy1_ad(i, j) = 0.0_8
          ELSE
            pt_ad(i, j) = pt_ad(i, j) + fy_ad(i, j)
            fy_ad(i, j) = 0.0_8
            delp_ad(i, j) = delp_ad(i, j) + fy1_ad(i, j)
            fy1_ad(i, j) = 0.0_8
          END IF
        END DO
      END DO
      CALL POPCONTROL1B(branch)
      IF (branch .NE. 0) CALL FILL_4CORNERS_ADM(w, w_ad, 2)
    ELSE
      IF (branch .EQ. 1) THEN
        DO j=jep1,js-1,-1
          DO i=iep1,is-1,-1
            temp_ad0 = ptc_ad(i, j)/delpc(i, j)
            temp_ad1 = rarea(i, j)*temp_ad0
            pt_ad(i, j) = pt_ad(i, j) + delp(i, j)*temp_ad0
            fx_ad(i, j) = fx_ad(i, j) + temp_ad1
            fx_ad(i+1, j) = fx_ad(i+1, j) - temp_ad1
            fy_ad(i, j) = fy_ad(i, j) + temp_ad1
            fy_ad(i, j+1) = fy_ad(i, j+1) - temp_ad1
            delpc_ad(i, j) = delpc_ad(i, j) - (pt(i, j)*delp(i, j)+rarea&
&             (i, j)*(fx(i, j)-fx(i+1, j)+fy(i, j)-fy(i, j+1)))*temp_ad0&
&             /delpc(i, j)
            delp_ad(i, j) = delp_ad(i, j) + delpc_ad(i, j) + pt(i, j)*&
&             temp_ad0
            ptc_ad(i, j) = 0.0_8
            CALL POPREAL8(delpc(i, j))
            temp_ad2 = rarea(i, j)*delpc_ad(i, j)
            fx1_ad(i, j) = fx1_ad(i, j) + temp_ad2
            fx1_ad(i+1, j) = fx1_ad(i+1, j) - temp_ad2
            fy1_ad(i, j) = fy1_ad(i, j) + temp_ad2
            fy1_ad(i, j+1) = fy1_ad(i, j+1) - temp_ad2
            delpc_ad(i, j) = 0.0_8
          END DO
        END DO
      ELSE
        DO j=jep1,js-1,-1
          DO i=iep1,is-1,-1
            pt_ad(i, j) = pt_ad(i, j) + ptc_ad(i, j)
            ptc_ad(i, j) = 0.0_8
            temp_ad = rarea(i, j)*delpc_ad(i, j)
            delp_ad(i, j) = delp_ad(i, j) + delpc_ad(i, j)
            fx1_ad(i, j) = fx1_ad(i, j) + temp_ad
            fx1_ad(i+1, j) = fx1_ad(i+1, j) - temp_ad
            fy1_ad(i, j) = fy1_ad(i, j) + temp_ad
            fy1_ad(i, j+1) = fy1_ad(i, j+1) - temp_ad
            delpc_ad(i, j) = 0.0_8
          END DO
        END DO
      END IF
      DO j=jep1+1,js-1,-1
        DO i=iep1,is-1,-1
          CALL POPREAL8(fy(i, j))
          fy1_ad(i, j) = fy1_ad(i, j) + fy(i, j)*fy_ad(i, j)
          fy_ad(i, j) = fy1(i, j)*fy_ad(i, j)
          CALL POPREAL8(fy1(i, j))
          vt_ad(i, j) = vt_ad(i, j) + fy1(i, j)*fy1_ad(i, j)
          fy1_ad(i, j) = vt(i, j)*fy1_ad(i, j)
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            pt_ad(i, j-1) = pt_ad(i, j-1) + fy_ad(i, j)
            fy_ad(i, j) = 0.0_8
            delp_ad(i, j-1) = delp_ad(i, j-1) + fy1_ad(i, j)
            fy1_ad(i, j) = 0.0_8
          ELSE
            pt_ad(i, j) = pt_ad(i, j) + fy_ad(i, j)
            fy_ad(i, j) = 0.0_8
            delp_ad(i, j) = delp_ad(i, j) + fy1_ad(i, j)
            fy1_ad(i, j) = 0.0_8
          END IF
        END DO
      END DO
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL8ARRAY(pt, (ied-isd+1)*(jed-jsd+1))
      CALL FILL2_4CORNERS_ADM(delp, delp_ad, pt, pt_ad, 2)
    END IF
    CALL POPCONTROL2B(branch)
    IF (branch .EQ. 0) THEN
      DO j=jep1,js-1,-1
        DO i=iep1+1,is-1,-1
          CALL POPREAL8(fx1(i, j))
          ut_ad(i, j) = ut_ad(i, j) + fx1(i, j)*fx1_ad(i, j)
          fx1_ad(i, j) = ut(i, j)*fx1_ad(i, j)
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            delp_ad(i-1, j) = delp_ad(i-1, j) + fx1_ad(i, j)
            fx1_ad(i, j) = 0.0_8
          ELSE
            delp_ad(i, j) = delp_ad(i, j) + fx1_ad(i, j)
            fx1_ad(i, j) = 0.0_8
          END IF
        END DO
      END DO
    ELSE IF (branch .EQ. 1) THEN
      DO j=jep1,js-1,-1
        DO i=iep1+1,is-1,-1
          CALL POPREAL8(fx(i, j))
          fx1_ad(i, j) = fx1_ad(i, j) + fx(i, j)*fx_ad(i, j)
          fx_ad(i, j) = fx1(i, j)*fx_ad(i, j)
          CALL POPREAL8(fx1(i, j))
          ut_ad(i, j) = ut_ad(i, j) + fx1(i, j)*fx1_ad(i, j)
          fx1_ad(i, j) = ut(i, j)*fx1_ad(i, j)
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            pt_ad(i-1, j) = pt_ad(i-1, j) + fx_ad(i, j)
            fx_ad(i, j) = 0.0_8
            delp_ad(i-1, j) = delp_ad(i-1, j) + fx1_ad(i, j)
            fx1_ad(i, j) = 0.0_8
          ELSE
            pt_ad(i, j) = pt_ad(i, j) + fx_ad(i, j)
            fx_ad(i, j) = 0.0_8
            delp_ad(i, j) = delp_ad(i, j) + fx1_ad(i, j)
            fx1_ad(i, j) = 0.0_8
          END IF
        END DO
      END DO
    ELSE
      DO j=je+1,js-1,-1
        DO i=ie+2,is-1,-1
          CALL POPREAL8(fx(i, j))
          fx1_ad(i, j) = fx1_ad(i, j) + fx(i, j)*fx_ad(i, j)
          fx_ad(i, j) = fx1(i, j)*fx_ad(i, j)
          CALL POPREAL8(fx1(i, j))
          ut_ad(i, j) = ut_ad(i, j) + fx1(i, j)*fx1_ad(i, j)
          fx1_ad(i, j) = ut(i, j)*fx1_ad(i, j)
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            pt_ad(i-1, j) = pt_ad(i-1, j) + fx_ad(i, j)
            fx_ad(i, j) = 0.0_8
            delp_ad(i-1, j) = delp_ad(i-1, j) + fx1_ad(i, j)
            fx1_ad(i, j) = 0.0_8
          ELSE
            pt_ad(i, j) = pt_ad(i, j) + fx_ad(i, j)
            fx_ad(i, j) = 0.0_8
            delp_ad(i, j) = delp_ad(i, j) + fx1_ad(i, j)
            fx1_ad(i, j) = 0.0_8
          END IF
        END DO
      END DO
      CALL POPCONTROL1B(branch)
      IF (branch .NE. 0) CALL FILL_4CORNERS_ADM(w, w_ad, 1)
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL8ARRAY(pt, (ied-isd+1)*(jed-jsd+1))
      CALL FILL2_4CORNERS_ADM(delp, delp_ad, pt, pt_ad, 1)
    END IF
    DO j=jep1+1,js-1,-1
      DO i=iep1,is-1,-1
        vt_ad(i, j) = dt2*dx(i, j)*sina_v(i, j)*vt_ad(i, j)
      END DO
    END DO
    DO j=jep1,js-1,-1
      DO i=iep1+1,is-1,-1
        ut_ad(i, j) = dt2*dy(i, j)*sina_u(i, j)*ut_ad(i, j)
      END DO
    END DO
    CALL POPREAL8ARRAY(ua, (ied-isd+1)*(jed-jsd+1))
    CALL POPREAL8ARRAY(va, (ied-isd+1)*(jed-jsd+1))
    CALL D2A2C_VECT_ADM(u, u_ad, v, v_ad, ua, ua_ad, va, va_ad, uc, &
&                 uc_ad, vc, vc_ad, ut, ut_ad, vt, vt_ad, dord4)
  END SUBROUTINE C_SW_ADM
  SUBROUTINE C_SW(delpc, delp, ptc, pt, u, v, w, uc, vc, ua, va, wc, ut&
&   , vt, dt2, hydrostatic, dord4)
    IMPLICIT NONE
!#endif
    REAL, DIMENSION(isd:ied, jsd:jed + 1), INTENT(INOUT) :: u, vc
    REAL, DIMENSION(isd:ied + 1, jsd:jed), INTENT(INOUT) :: v, uc
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(INOUT) :: delp, pt, ua, va&
&   , w
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(OUT) :: delpc, ptc, ut, vt&
&   , wc
    REAL, INTENT(IN) :: dt2
    LOGICAL, INTENT(IN) :: hydrostatic
    LOGICAL, INTENT(IN) :: dord4
! Local:
    REAL, DIMENSION(is - 1:ie + 1, js - 1:je + 1) :: vort, ke
    REAL, DIMENSION(is - 1:ie + 2, js - 1:je + 1) :: fx, fx1, fx2
    REAL, DIMENSION(is - 1:ie + 1, js - 1:je + 2) :: fy, fy1, fy2
    REAL :: dt4
    INTEGER :: i, j, is2, ie1
    INTEGER :: iep1, jep1
    INTRINSIC MAX
    INTRINSIC MIN
!#ifdef DOUBLE_GRADIENTS_CSW
!      real*8 :: delp_dp, pt_dp, dp_r8
!      real*8, dimension(is-1:ie+1,js-1:je+1):: ke_dp
!      real*8, dimension(is-1:ie+2,js-1:je+1):: fx_dp, fx1_dp
!      real*8, dimension(is-1:ie+1,js-1:je+2):: fy_dp, fy1_dp
!#endif
!#ifdef FIX_C_BOUNDARY 
!      iep1 = ie;   jep1 = je
!#else
    iep1 = ie + 1
    jep1 = je + 1
!#endif
    CALL D2A2C_VECT(u, v, ua, va, uc, vc, ut, vt, dord4)
!     call d2a2c_vect_v2(u, v, ua, va, uc, vc, ut, vt)
    DO j=js-1,jep1
      DO i=is-1,iep1+1
        ut(i, j) = dt2*ut(i, j)*dy(i, j)*sina_u(i, j)
      END DO
    END DO
    DO j=js-1,jep1+1
      DO i=is-1,iep1
        vt(i, j) = dt2*vt(i, j)*dx(i, j)*sina_v(i, j)
      END DO
    END DO
!----------------
! Transport delp:
!----------------
! Xdir:
    IF (grid_type .LT. 3) CALL FILL2_4CORNERS(delp, pt, 1)
    IF (hydrostatic) THEN
      IF (shallow_water) THEN
        DO j=js-1,jep1
          DO i=is-1,iep1+1
            IF (ut(i, j) .GT. 0.) THEN
              fx1(i, j) = delp(i-1, j)
            ELSE
              fx1(i, j) = delp(i, j)
            END IF
            fx1(i, j) = ut(i, j)*fx1(i, j)
          END DO
        END DO
      ELSE
        DO j=js-1,jep1
          DO i=is-1,iep1+1
            IF (ut(i, j) .GT. 0.) THEN
              fx1(i, j) = delp(i-1, j)
              fx(i, j) = pt(i-1, j)
            ELSE
              fx1(i, j) = delp(i, j)
              fx(i, j) = pt(i, j)
            END IF
            fx1(i, j) = ut(i, j)*fx1(i, j)
            fx(i, j) = fx1(i, j)*fx(i, j)
          END DO
        END DO
      END IF
    ELSE
      IF (grid_type .LT. 3) CALL FILL_4CORNERS(w, 1)
      DO j=js-1,je+1
        DO i=is-1,ie+2
          IF (ut(i, j) .GT. 0.) THEN
            fx1(i, j) = delp(i-1, j)
            fx(i, j) = pt(i-1, j)
            fx2(i, j) = w(i-1, j)
          ELSE
            fx1(i, j) = delp(i, j)
            fx(i, j) = pt(i, j)
            fx2(i, j) = w(i, j)
          END IF
          fx1(i, j) = ut(i, j)*fx1(i, j)
          fx(i, j) = fx1(i, j)*fx(i, j)
          fx2(i, j) = fx1(i, j)*fx2(i, j)
        END DO
      END DO
    END IF
! Ydir:
    IF (grid_type .LT. 3) CALL FILL2_4CORNERS(delp, pt, 2)
    IF (hydrostatic) THEN
      DO j=js-1,jep1+1
        DO i=is-1,iep1
          IF (vt(i, j) .GT. 0.) THEN
            fy1(i, j) = delp(i, j-1)
            fy(i, j) = pt(i, j-1)
          ELSE
            fy1(i, j) = delp(i, j)
            fy(i, j) = pt(i, j)
          END IF
          fy1(i, j) = vt(i, j)*fy1(i, j)
          fy(i, j) = fy1(i, j)*fy(i, j)
        END DO
      END DO
      IF (shallow_water) THEN
        DO j=js-1,jep1
          DO i=is-1,iep1
            delpc(i, j) = delp(i, j) + (fx1(i, j)-fx1(i+1, j)+fy1(i, j)-&
&             fy1(i, j+1))*rarea(i, j)
            ptc(i, j) = pt(i, j)
          END DO
        END DO
      ELSE
!#ifdef DOUBLE_GRADIENTS_CSW
!           fy1_dp = fy1
!           fy_dp  = fy
!           fx1_dp = fx1
!           fx_dp  = fx
!           do j=js-1,jep1
!              do i=is-1,iep1
!                 delp_dp = delp(i,j)
!                 pt_dp   = pt(i,j)
!                 delpc(i,j) =        delp_dp + (fx1_dp(i,j)-fx1_dp(i+1,j)+fy1_dp(i,j)-fy1_dp(i,j+1))*rarea(i,j)
!                   ptc(i,j) = (pt_dp*delp_dp + ( fx_dp(i,j)- fx_dp(i+1,j)+ fy_dp(i,j)- fy_dp(i,j+1))*rarea(i,j))/delpc(i,j)
!              enddo
!           enddo
!#else
        DO j=js-1,jep1
          DO i=is-1,iep1
            delpc(i, j) = delp(i, j) + (fx1(i, j)-fx1(i+1, j)+fy1(i, j)-&
&             fy1(i, j+1))*rarea(i, j)
            ptc(i, j) = (pt(i, j)*delp(i, j)+(fx(i, j)-fx(i+1, j)+fy(i, &
&             j)-fy(i, j+1))*rarea(i, j))/delpc(i, j)
          END DO
        END DO
      END IF
    ELSE
!#endif
      IF (grid_type .LT. 3) CALL FILL_4CORNERS(w, 2)
      DO j=js-1,je+2
        DO i=is-1,ie+1
          IF (vt(i, j) .GT. 0.) THEN
            fy1(i, j) = delp(i, j-1)
            fy(i, j) = pt(i, j-1)
            fy2(i, j) = w(i, j-1)
          ELSE
            fy1(i, j) = delp(i, j)
            fy(i, j) = pt(i, j)
            fy2(i, j) = w(i, j)
          END IF
          fy1(i, j) = vt(i, j)*fy1(i, j)
          fy(i, j) = fy1(i, j)*fy(i, j)
          fy2(i, j) = fy1(i, j)*fy2(i, j)
        END DO
      END DO
      DO j=js-1,je+1
        DO i=is-1,ie+1
          delpc(i, j) = delp(i, j) + (fx1(i, j)-fx1(i+1, j)+fy1(i, j)-&
&           fy1(i, j+1))*rarea(i, j)
          ptc(i, j) = (pt(i, j)*delp(i, j)+(fx(i, j)-fx(i+1, j)+fy(i, j)&
&           -fy(i, j+1))*rarea(i, j))/delpc(i, j)
          wc(i, j) = (w(i, j)*delp(i, j)+(fx2(i, j)-fx2(i+1, j)+fy2(i, j&
&           )-fy2(i, j+1))*rarea(i, j))/delpc(i, j)
        END DO
      END DO
    END IF
!------------
! Compute KE:
!------------
    DO j=js-1,jep1
      DO i=is-1,iep1
        IF (ua(i, j) .GT. 0.) THEN
          IF (i .EQ. 1) THEN
            ke(1, j) = uc(1, j)*sina_u(1, j) + v(1, j)*cosa_u(1, j)
          ELSE IF (i .EQ. npx) THEN
            ke(i, j) = uc(npx, j)*sina_u(npx, j) - v(npx, j)*cosa_u(npx&
&             , j)
          ELSE
            ke(i, j) = uc(i, j)
          END IF
        ELSE IF (i .EQ. 0) THEN
          ke(0, j) = uc(1, j)*sina_u(1, j) - v(1, j)*cosa_u(1, j)
        ELSE IF (i .EQ. npx - 1) THEN
          ke(i, j) = uc(npx, j)*sina_u(npx, j) + v(npx, j)*cosa_u(npx, j&
&           )
        ELSE
          ke(i, j) = uc(i+1, j)
        END IF
      END DO
    END DO
    DO j=js-1,jep1
      DO i=is-1,iep1
        IF (va(i, j) .GT. 0.) THEN
          IF (j .EQ. 1) THEN
            vort(i, 1) = vc(i, 1)*sina_v(i, 1) + u(i, 1)*cosa_v(i, 1)
          ELSE IF (j .EQ. npy) THEN
            vort(i, j) = vc(i, npy)*sina_v(i, npy) - u(i, npy)*cosa_v(i&
&             , npy)
          ELSE
            vort(i, j) = vc(i, j)
          END IF
        ELSE IF (j .EQ. 0) THEN
          vort(i, 0) = vc(i, 1)*sina_v(i, 1) - u(i, 1)*cosa_v(i, 1)
        ELSE IF (j .EQ. npy - 1) THEN
          vort(i, j) = vc(i, npy)*sina_v(i, npy) + u(i, npy)*cosa_v(i, &
&           npy)
        ELSE
          vort(i, j) = vc(i, j+1)
        END IF
      END DO
    END DO
    dt4 = 0.5*dt2
    DO j=js-1,jep1
      DO i=is-1,iep1
        ke(i, j) = dt4*(ua(i, j)*ke(i, j)+va(i, j)*vort(i, j))
      END DO
    END DO
    IF (2 .LT. is) THEN
      is2 = is
    ELSE
      is2 = 2
    END IF
    IF (npx - 1 .GT. ie + 1) THEN
      ie1 = ie + 1
    ELSE
      ie1 = npx - 1
    END IF
    DO j=js-1,je+1
      DO i=is2,ie1
        fx(i, j) = uc(i, j)*dxc(i, j)
      END DO
      IF (is .EQ. 1) fx(1, j) = uc(1, j)*sina_u(1, j)*dxc(1, j)
      IF (ie + 1 .EQ. npx) fx(npx, j) = uc(npx, j)*sina_u(npx, j)*dxc(&
&         npx, j)
    END DO
    DO j=js,je+1
      IF (j .EQ. 1 .OR. j .EQ. npy) THEN
        DO i=is-1,ie+1
          fy(i, j) = vc(i, j)*sina_v(i, j)*dyc(i, j)
        END DO
      ELSE
        DO i=is-1,ie+1
          fy(i, j) = vc(i, j)*dyc(i, j)
        END DO
      END IF
    END DO
!#endif
    DO j=js,je+1
      DO i=is,ie+1
        vort(i, j) = fx(i, j-1) - fx(i, j) - fy(i-1, j) + fy(i, j)
      END DO
    END DO
! Remove the extra term at the corners:
    IF (sw_corner) vort(1, 1) = vort(1, 1) + fy(0, 1)
    IF (se_corner) vort(npx, 1) = vort(npx, 1) - fy(npx, 1)
    IF (ne_corner) vort(npx, npy) = vort(npx, npy) - fy(npx, npy)
    IF (nw_corner) vort(1, npy) = vort(1, npy) + fy(0, npy)
!----------------------------
! Compute absolute vorticity
!----------------------------
    DO j=js,je+1
      DO i=is,ie+1
        vort(i, j) = fc(i, j) + rarea_c(i, j)*vort(i, j)
      END DO
    END DO
!----------------------------------
! Transport absolute vorticity:
!----------------------------------
    DO j=js,je
      DO i=is,iep1
        IF (i .EQ. 1 .OR. i .EQ. npx) THEN
          fy1(i, j) = dt2*v(i, j)*sina_u(i, j)
        ELSE
          fy1(i, j) = dt2*(v(i, j)-uc(i, j)*cosa_u(i, j))/sina_u(i, j)
        END IF
        IF (fy1(i, j) .GT. 0.) THEN
          fy(i, j) = vort(i, j)
        ELSE
          fy(i, j) = vort(i, j+1)
        END IF
      END DO
    END DO
    DO j=js,jep1
      IF (j .EQ. 1 .OR. j .EQ. npy) THEN
        DO i=is,ie
          fx1(i, j) = dt2*u(i, j)*sina_v(i, j)
          IF (fx1(i, j) .GT. 0.) THEN
            fx(i, j) = vort(i, j)
          ELSE
            fx(i, j) = vort(i+1, j)
          END IF
        END DO
      ELSE
        DO i=is,ie
          fx1(i, j) = dt2*(u(i, j)-vc(i, j)*cosa_v(i, j))/sina_v(i, j)
          IF (fx1(i, j) .GT. 0.) THEN
            fx(i, j) = vort(i, j)
          ELSE
            fx(i, j) = vort(i+1, j)
          END IF
        END DO
      END IF
    END DO
! Update time-centered winds on the C-Grid
!#ifdef DOUBLE_GRADIENTS_CSW
!     ke_dp  = ke
!     fy1_dp = fy1
!     fy_dp  = fy
!     fx1_dp = fx1
!     fx_dp  = fx
!     do j=js,je
!        do i=is,iep1
!           dp_r8 = fy1_dp(i,j)*fy_dp(i,j) + rdxc(i,j)*(ke_dp(i-1,j)-ke_dp(i,j))
!           uc(i,j) = uc(i,j) + dp_r8
!        enddo
!     enddo
!     do j=js,jep1
!        do i=is,ie
!           dp_r8 = fx1_dp(i,j)*fx_dp(i,j) + rdyc(i,j)*(ke_dp(i,j-1)-ke_dp(i,j))
!           vc(i,j) = vc(i,j) - dp_r8
!        enddo
!     enddo
!#else
    DO j=js,je
      DO i=is,iep1
        uc(i, j) = uc(i, j) + fy1(i, j)*fy(i, j) + rdxc(i, j)*(ke(i-1, j&
&         )-ke(i, j))
      END DO
    END DO
    DO j=js,jep1
      DO i=is,ie
        vc(i, j) = vc(i, j) - fx1(i, j)*fx(i, j) + rdyc(i, j)*(ke(i, j-1&
&         )-ke(i, j))
      END DO
    END DO
  END SUBROUTINE C_SW
!  Differentiation of d_sw in reverse (adjoint) mode (with options r8):
!   gradient     of useful results: yfx_adv q crx_adv u v w delp
!                ua xfx_adv uc ptc xflux cry_adv delpc va vc yflux
!                pkz divg_d pt cx cy
!   with respect to varying inputs: yfx_adv q crx_adv u v w delp
!                ua xfx_adv uc ptc xflux cry_adv delpc va vc yflux
!                pkz divg_d pt cx cy
!-------------------------------------------------------------------------------
!
!     d_sw :: D-Grid Shallow Water Routine
!
  SUBROUTINE D_SW_ADM(delpc, delpc_ad, delp, delp_ad, ptc, ptc_ad, pt, &
&   pt_ad, u, u_ad, v, v_ad, w, w_ad, uc, uc_ad, vc, vc_ad, ua, ua_ad, &
&   va, va_ad, divg_d, divg_d_ad, xflux, xflux_ad, yflux, yflux_ad, cx, &
&   cx_ad, cy, cy_ad, crx_adv, crx_adv_ad, cry_adv, cry_adv_ad, xfx_adv&
&   , xfx_adv_ad, yfx_adv, yfx_adv_ad, zvir, sphum, nq, q, q_ad, k, km, &
&   inline_q, pkz, pkz_ad, dt, hord_tr, hord_mt, hord_vt, hord_tm, &
&   hord_dp, hord_tr_pert, hord_mt_pert, hord_vt_pert, hord_tm_pert, &
&   hord_dp_pert, nord, dddmp, d2_bg, d4_bg, vtdm4, d_con, hydrostatic, &
&   ppm_limiter)
    IMPLICIT NONE
! damping applied to relative vorticity:
!   if ( vtdm4>0. ) then
!      damp4 = (vtdm4*da_min)**(nord+1)
!      call del6_flux(nord, npx, npy, damp4, wk, u, v)
!   endif
    INTEGER, INTENT(IN) :: hord_tr, hord_mt, hord_vt, hord_tm, hord_dp
    INTEGER, INTENT(IN) :: hord_tr_pert, hord_mt_pert, hord_vt_pert, hord_tm_pert, hord_dp_pert
! nord=1 (del-4) or 3 (del-8)
    INTEGER, INTENT(IN) :: nord
    INTEGER, INTENT(IN) :: sphum, nq, k, km
    REAL, INTENT(IN) :: dt, dddmp, d2_bg, d4_bg, vtdm4, d_con
    REAL, INTENT(IN) :: zvir
    REAL, INTENT(IN) :: ppm_limiter
    REAL(p_precision), INTENT(IN) :: pkz(is:ie, js:je)
    REAL(p_precision) :: pkz_ad(is:ie, js:je)
! divergence
    REAL, INTENT(INOUT) :: divg_d(isd:ied+1, jsd:jed+1)
    REAL, INTENT(INOUT) :: divg_d_ad(isd:ied+1, jsd:jed+1)
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(INOUT) :: delp, pt, ua, va&
&   , w
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(INOUT) :: delp_ad
    REAL, DIMENSION(isd:ied, jsd:jed + 1), INTENT(INOUT) :: u, vc
    REAL, DIMENSION(isd:ied, jsd:jed+1), INTENT(INOUT) :: u_ad
    REAL, DIMENSION(isd:ied + 1, jsd:jed), INTENT(INOUT) :: v, uc
    REAL, DIMENSION(isd:ied+1, jsd:jed), INTENT(INOUT) :: v_ad
    REAL, INTENT(INOUT) :: q(isd:ied, jsd:jed, km, nq)
    REAL, INTENT(INOUT) :: q_ad(isd:ied, jsd:jed, km, nq)
    REAL, DIMENSION(isd:ied, jsd:jed) :: delpc, ptc
    REAL, DIMENSION(isd:ied, jsd:jed) :: delpc_ad, ptc_ad
! The flux capacitors:
    REAL, INTENT(INOUT) :: xflux(is:ie+1, js:je)
    REAL, INTENT(INOUT) :: xflux_ad(is:ie+1, js:je)
    REAL, INTENT(INOUT) :: yflux(is:ie, js:je+1)
    REAL, INTENT(INOUT) :: yflux_ad(is:ie, js:je+1)
!------------------------
    REAL, INTENT(INOUT) :: cx(is:ie+1, jsd:jed)
    REAL, INTENT(INOUT) :: cx_ad(is:ie+1, jsd:jed)
    REAL, INTENT(INOUT) :: cy(isd:ied, js:je+1)
    REAL, INTENT(INOUT) :: cy_ad(isd:ied, js:je+1)
    LOGICAL, INTENT(IN) :: hydrostatic
    LOGICAL, INTENT(IN) :: inline_q
    REAL, DIMENSION(is:ie + 1, jsd:jed) :: crx_adv, xfx_adv
    REAL, DIMENSION(is:ie+1, jsd:jed) :: crx_adv_ad, xfx_adv_ad
    REAL, DIMENSION(isd:ied, js:je + 1) :: cry_adv, yfx_adv
    REAL, DIMENSION(isd:ied, js:je+1) :: cry_adv_ad, yfx_adv_ad
! Local:
    REAL(i_precision), DIMENSION(isd:ied + 1, jsd:jed) :: ut, uci
    REAL(i_precision), DIMENSION(isd:ied+1, jsd:jed) :: ut_ad, uci_ad
    REAL(i_precision), DIMENSION(isd:ied, jsd:jed + 1) :: vt, vci
    REAL(i_precision), DIMENSION(isd:ied, jsd:jed+1) :: vt_ad, vci_ad
    REAL(i_precision) :: damp
    REAL(i_precision) :: damp_ad
    REAL(i_precision), SAVE :: r1b4=1.0/4.0
    REAL(i_precision), SAVE :: r1b16=1.0/16.0
    REAL(ke_precision), DIMENSION(is:ie + 1, js:je + 1) :: ub, vb
    REAL(ke_precision), DIMENSION(is:ie+1, js:je+1) :: ub_ad, vb_ad
!  needs this for corner_comm
    REAL(ke_precision) :: ke(isd:ied+1, jsd:jed+1)
    REAL(ke_precision) :: ke_ad(isd:ied+1, jsd:jed+1)
    REAL(ke_precision) :: dt4, dt5, dt6
!#ifdef FULL_SMAG
!      real :: vt2(is-1:ie+1,js-1:je+1)     !  work array
!#endif
!  work array
    REAL :: wk(isd:ied, jsd:jed)
    REAL :: wk_ad(isd:ied, jsd:jed)
! Vorticity
    REAL :: vort(isd:ied, jsd:jed)
    REAL :: vort_ad(isd:ied, jsd:jed)
! 1-D X-direction Fluxes
    REAL :: fx(is:ie+1, js:je)
    REAL :: fx_ad(is:ie+1, js:je)
! 1-D Y-direction Fluxes
    REAL :: fy(is:ie, js:je+1)
    REAL :: fy_ad(is:ie, js:je+1)
! 1-D X-direction Fluxes
    REAL :: px(is:ie+1, js:je)
    REAL :: px_ad(is:ie+1, js:je)
! 1-D Y-direction Fluxes
    REAL :: py(is:ie, js:je+1)
    REAL :: py_ad(is:ie, js:je+1)
    REAL :: ra_x(is:ie, jsd:jed)
    REAL :: ra_x_ad(is:ie, jsd:jed)
    REAL :: ra_y(isd:ied, js:je)
    REAL :: ra_y_ad(isd:ied, js:je)
! work x-dir flux
    REAL :: gx(is:ie+1, js:je)
    REAL :: gx_ad(is:ie+1, js:je)
! work Y-dir flux array
    REAL :: gy(is:ie, js:je+1)
    REAL :: gy_ad(is:ie, js:je+1)
    LOGICAL :: fill_c
!#ifdef DOUBLE_GRADIENTS
!      real*8 :: pt_dp, delp_dp0, delp_dp
!      real*8 :: pt_fx_r8(is:ie+1,js:je+1)
!      real*8 :: pt_fy_r8(is:ie  ,js:je+1)
!      real*8 :: dp_fx_r8(is:ie+1,js:je  )
!      real*8 :: dp_fy_r8(is:ie  ,js:je+1)
!#endif
    REAL :: damp2, damp4, dd8, u2, v2, du2, dv2
    REAL :: damp2_ad, u2_ad, v2_ad, du2_ad, dv2_ad
    INTEGER :: i, j, is2, ie1, js2, je1, n, nt, n2, iq
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC ABS
    INTEGER :: arg1
    INTEGER :: arg2
    INTEGER :: branch
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
    INTEGER :: ad_from5
    INTEGER :: ad_to5
    INTEGER :: ad_from6
    INTEGER :: ad_to6
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(INOUT) :: va_ad
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(INOUT) :: pt_ad
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(INOUT) :: w_ad
    REAL, DIMENSION(isd:ied+1, jsd:jed), INTENT(INOUT) :: uc_ad
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(INOUT) :: ua_ad
    REAL, DIMENSION(isd:ied, jsd:jed+1), INTENT(INOUT) :: vc_ad
    REAL :: temp0
    REAL :: y2_ad
    REAL(i_precision) :: temp_ad39
    REAL(i_precision) :: temp_ad38
    INTEGER :: min4
    REAL(i_precision) :: temp_ad37
    INTEGER :: min3
    REAL(i_precision) :: temp_ad36
    INTEGER :: min2
    REAL(i_precision) :: temp_ad35
    REAL :: y3_ad
    INTEGER :: min1
    REAL(i_precision) :: temp_ad34
    REAL(i_precision) :: temp_ad33
    REAL(i_precision) :: temp_ad32
    REAL(i_precision) :: temp_ad69
    REAL(i_precision) :: temp_ad31
    REAL(i_precision) :: temp_ad68
    REAL(i_precision) :: temp_ad30
    REAL*8 :: temp_ad67
    REAL(i_precision) :: temp_ad66
    REAL(i_precision) :: temp_ad65
    REAL(i_precision) :: temp_ad64
    REAL*8 :: temp_ad63
    REAL*8 :: temp_ad62
    REAL(ke_precision) :: temp_ad61
    REAL :: temp_ad60
    REAL :: temp_ad90
    REAL(i_precision) :: temp_ad29
    REAL(i_precision) :: temp_ad28
    REAL(i_precision) :: temp_ad27
    REAL(i_precision) :: temp_ad26
    REAL(i_precision) :: temp_ad25
    REAL(i_precision) :: temp_ad24
    REAL(i_precision) :: temp_ad23
    REAL(i_precision) :: temp_ad22
    REAL :: temp_ad59
    REAL(i_precision) :: temp_ad21
    REAL(i_precision) :: temp_ad58
    REAL(i_precision) :: temp_ad20
    REAL(i_precision) :: temp_ad57
    REAL(i_precision) :: temp_ad9
    REAL(i_precision) :: temp_ad56
    REAL(i_precision) :: temp_ad8
    REAL :: temp_ad55
    REAL(i_precision) :: temp_ad7
    REAL :: temp_ad54
    REAL(i_precision) :: temp_ad6
    REAL :: temp_ad53
    REAL(i_precision) :: temp_ad5
    REAL :: temp_ad52
    REAL :: temp_ad89
    REAL(i_precision) :: temp_ad4
    REAL :: temp_ad51
    REAL*8 :: temp_ad88
    REAL(i_precision) :: temp_ad3
    REAL :: temp_ad50
    REAL :: temp_ad87
    REAL(i_precision) :: temp_ad2
    REAL :: temp_ad86
    REAL(i_precision) :: temp_ad1
    REAL :: temp_ad85
    REAL(i_precision) :: temp_ad0
    REAL :: temp_ad84
    REAL :: temp_ad83
    REAL :: temp_ad82
    REAL :: temp_ad81
    REAL :: temp_ad80
    REAL(i_precision) :: temp_ad
    REAL(i_precision) :: temp_ad19
    REAL(i_precision) :: temp_ad18
    REAL(i_precision) :: temp_ad17
    REAL :: abs0_ad
    REAL(i_precision) :: temp_ad16
    REAL(i_precision) :: temp_ad15
    REAL(i_precision) :: temp_ad14
    REAL(i_precision) :: temp_ad13
    REAL(i_precision) :: temp_ad12
    REAL :: temp_ad49
    REAL :: abs0
    REAL(i_precision) :: temp_ad11
    REAL :: temp_ad48
    REAL(i_precision) :: temp_ad10
    REAL :: temp_ad47
    REAL :: temp_ad46
    REAL :: temp_ad45
    REAL :: temp_ad44
    REAL :: temp_ad43
    REAL :: max5_ad
    REAL(i_precision) :: temp_ad42
    REAL(i_precision) :: temp_ad79
    REAL(i_precision) :: temp_ad41
    REAL(ke_precision) :: temp_ad78
    REAL :: max6
    REAL(i_precision) :: temp_ad40
    REAL(ke_precision) :: temp_ad77
    REAL :: max5
    REAL(ke_precision) :: temp_ad76
    INTEGER :: max4
    REAL :: temp
    REAL(ke_precision) :: temp_ad75
    INTEGER :: max3
    REAL(i_precision) :: temp_ad74
    INTEGER :: max2
    REAL(i_precision) :: temp_ad73
    INTEGER :: max1
    REAL(i_precision) :: temp_ad72
    REAL*8 :: temp_ad71
    REAL :: max6_ad
    REAL(i_precision) :: temp_ad70
    REAL :: y3
    REAL :: y1_ad
    REAL :: y2
    REAL :: y1
! shallow water test case 1
    IF (test_case .EQ. 1) THEN
      DO j=jsd,jed
        DO i=is,ie+1
          xfx_adv(i, j) = dt*uc(i, j)/sina_u(i, j)
          IF (xfx_adv(i, j) .GT. 0.) THEN
            crx_adv(i, j) = xfx_adv(i, j)*rdxa(i-1, j)
            CALL PUSHCONTROL1B(0)
          ELSE
            crx_adv(i, j) = xfx_adv(i, j)*rdxa(i, j)
            CALL PUSHCONTROL1B(1)
          END IF
          xfx_adv(i, j) = dy(i, j)*xfx_adv(i, j)*sina_u(i, j)
        END DO
      END DO
      DO j=js,je+1
        DO i=isd,ied
          yfx_adv(i, j) = dt*vc(i, j)/sina_v(i, j)
          IF (yfx_adv(i, j) .GT. 0.) THEN
            cry_adv(i, j) = yfx_adv(i, j)*rdya(i, j-1)
            CALL PUSHCONTROL1B(0)
          ELSE
            cry_adv(i, j) = yfx_adv(i, j)*rdya(i, j)
            CALL PUSHCONTROL1B(1)
          END IF
          yfx_adv(i, j) = dx(i, j)*yfx_adv(i, j)*sina_v(i, j)
        END DO
      END DO
      CALL PUSHCONTROL1B(1)
    ELSE
! end grid_type choices
      IF (grid_type .LT. 3) THEN
        uci = uc
        vci = vc
! Interior:
        DO j=jsd,jed
          IF (j .NE. 0 .AND. j .NE. 1 .AND. j .NE. npy - 1 .AND. j .NE. &
&             npy) THEN
            DO i=is-1,ie+2
              ut(i, j) = (uci(i, j)-r1b4*cosa_u(i, j)*(vci(i-1, j)+vci(i&
&               , j)+vci(i-1, j+1)+vci(i, j+1)))*rsin_u(i, j)
            END DO
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
        DO j=js-1,je+2
          IF (j .NE. 1 .AND. j .NE. npy) THEN
            DO i=isd,ied
              vt(i, j) = (vci(i, j)-r1b4*cosa_v(i, j)*(uci(i, j-1)+uci(i&
&               +1, j-1)+uci(i, j)+uci(i+1, j)))*rsin_v(i, j)
            END DO
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
! West edge:
        IF (is .EQ. 1) THEN
          DO j=jsd,jed
            ut(1, j) = uci(1, j)*rsin_u(1, j)
          END DO
          IF (3 .LT. js) THEN
            max1 = js
          ELSE
            max1 = 3
          END IF
          IF (npy - 2 .GT. je + 1) THEN
            min1 = je + 1
          ELSE
            min1 = npy - 2
          END IF
          DO j=max1,min1
!             vt(0,j) = vci(0,j) - r1b4*cosa_v(0,j)*   &
            vt(0, j) = vci(0, j) + r1b4*cosa_v(1, j)*(ut(0, j-1)+ut(1, j&
&             -1)+ut(0, j)+ut(1, j))
            vt(1, j) = vci(1, j) - r1b4*cosa_v(1, j)*(ut(1, j-1)+ut(2, j&
&             -1)+ut(1, j)+ut(2, j))
          END DO
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
! East edge:
        IF (ie + 1 .EQ. npx) THEN
          DO j=jsd,jed
            ut(npx, j) = uci(npx, j)*rsin_u(npx, j)
          END DO
          IF (3 .LT. js) THEN
            max2 = js
          ELSE
            max2 = 3
          END IF
          IF (npy - 2 .GT. je + 1) THEN
            min2 = je + 1
          ELSE
            min2 = npy - 2
          END IF
          DO j=max2,min2
            vt(npx-1, j) = vci(npx-1, j) - r1b4*cosa_v(npx-1, j)*(ut(npx&
&             -1, j-1)+ut(npx, j-1)+ut(npx-1, j)+ut(npx, j))
!             vt(npx,j) = vci(npx,j) - r1b4*cosa_v(npx,j)*   &
            vt(npx, j) = vci(npx, j) + r1b4*cosa_v(npx-1, j)*(ut(npx, j-&
&             1)+ut(npx+1, j-1)+ut(npx, j)+ut(npx+1, j))
          END DO
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
! South (Bottom) edge:
        IF (js .EQ. 1) THEN
          DO i=isd,ied
            vt(i, 1) = vci(i, 1)*rsin_v(i, 1)
          END DO
          IF (3 .LT. is) THEN
            max3 = is
          ELSE
            max3 = 3
          END IF
          IF (npx - 2 .GT. ie + 1) THEN
            min3 = ie + 1
          ELSE
            min3 = npx - 2
          END IF
          DO i=max3,min3
!             ut(i,0) = uci(i,0) - r1b4*cosa_u(i,0)*   &
            ut(i, 0) = uci(i, 0) + r1b4*cosa_u(i, 1)*(vt(i-1, 0)+vt(i, 0&
&             )+vt(i-1, 1)+vt(i, 1))
            ut(i, 1) = uci(i, 1) - r1b4*cosa_u(i, 1)*(vt(i-1, 1)+vt(i, 1&
&             )+vt(i-1, 2)+vt(i, 2))
          END DO
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
! North edge:
        IF (je + 1 .EQ. npy) THEN
          DO i=isd,ied
            vt(i, npy) = vci(i, npy)*rsin_v(i, npy)
          END DO
          IF (3 .LT. is) THEN
            max4 = is
          ELSE
            max4 = 3
          END IF
          IF (npx - 2 .GT. ie + 1) THEN
            min4 = ie + 1
          ELSE
            min4 = npx - 2
          END IF
          DO i=max4,min4
            ut(i, npy-1) = uci(i, npy-1) - r1b4*cosa_u(i, npy-1)*(vt(i-1&
&             , npy-1)+vt(i, npy-1)+vt(i-1, npy)+vt(i, npy))
!             ut(i,npy) = uci(i,npy) - r1b4*cosa_u(i,npy)*   &
            ut(i, npy) = uci(i, npy) + r1b4*cosa_u(i, npy-1)*(vt(i-1, &
&             npy)+vt(i, npy)+vt(i-1, npy+1)+vt(i, npy+1))
          END DO
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
        IF (sw_corner) THEN
          damp = 1.0/(1.0-r1b16*cosa_u(2, 1)*cosa_v(1, 2))
          ut(2, 0) = (uci(2, 0)-r1b4*cosa_u(2, 0)*(vt(1, 1)+vt(2, 1)+vt(&
&           2, 0)+vci(1, 0)-r1b4*cosa_v(1, 0)*(ut(1, 0)+ut(1, -1)+ut(2, &
&           -1))))*damp
          ut(2, 1) = (uci(2, 1)-r1b4*cosa_u(2, 1)*(vt(1, 1)+vt(2, 1)+vt(&
&           2, 2)+vci(1, 2)-r1b4*cosa_v(1, 2)*(ut(1, 1)+ut(1, 2)+ut(2, 2&
&           ))))*damp
          vt(1, 2) = (vci(1, 2)-r1b4*cosa_v(1, 2)*(ut(1, 1)+ut(1, 2)+ut(&
&           2, 2)+uci(2, 1)-r1b4*cosa_u(2, 1)*(vt(1, 1)+vt(2, 1)+vt(2, 2&
&           ))))*damp
          vt(0, 2) = (vci(0, 2)-r1b4*cosa_v(0, 2)*(ut(1, 1)+ut(1, 2)+ut(&
&           0, 2)+uci(0, 1)-r1b4*cosa_u(0, 1)*(vt(0, 1)+vt(-1, 1)+vt(-1&
&           , 2))))*damp
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
        IF (se_corner) THEN
          damp = 1.0/(1.0-r1b16*cosa_u(npx-1, 1)*cosa_v(npx-1, 2))
          ut(npx-1, 0) = (uci(npx-1, 0)+r1b4*cosa_u(npx-1, 1)*(vt(npx-1&
&           , 1)+vt(npx-2, 1)+vt(npx-2, 0)+vci(npx-1, 0)+r1b4*cosa_v(npx&
&           -1, 2)*(ut(npx, 0)+ut(npx, -1)+ut(npx-1, -1))))*damp
          ut(npx-1, 1) = (uci(npx-1, 1)-r1b4*cosa_u(npx-1, 1)*(vt(npx-1&
&           , 1)+vt(npx-2, 1)+vt(npx-2, 2)+vci(npx-1, 2)-r1b4*cosa_v(npx&
&           -1, 2)*(ut(npx, 1)+ut(npx, 2)+ut(npx-1, 2))))*damp
          vt(npx-1, 2) = (vci(npx-1, 2)-r1b4*cosa_v(npx-1, 2)*(ut(npx, 1&
&           )+ut(npx, 2)+ut(npx-1, 2)+uci(npx-1, 1)-r1b4*cosa_u(npx-1, 1&
&           )*(vt(npx-1, 1)+vt(npx-2, 1)+vt(npx-2, 2))))*damp
          vt(npx, 2) = (vci(npx, 2)+r1b4*cosa_v(npx-1, 2)*(ut(npx, 1)+ut&
&           (npx, 2)+ut(npx+1, 2)+uci(npx+1, 1)+r1b4*cosa_u(npx-1, 1)*(&
&           vt(npx, 1)+vt(npx+1, 1)+vt(npx+1, 2))))*damp
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
        IF (ne_corner) THEN
          damp = 1.0/(1.0-r1b16*cosa_u(npx-1, npy-1)*cosa_v(npx-1, npy-1&
&           ))
          ut(npx-1, npy) = (uci(npx-1, npy)+r1b4*cosa_u(npx-1, npy-1)*(&
&           vt(npx-1, npy)+vt(npx-2, npy)+vt(npx-2, npy+1)+vci(npx-1, &
&           npy+1)+r1b4*cosa_v(npx-1, npy-1)*(ut(npx, npy)+ut(npx, npy+1&
&           )+ut(npx-1, npy+1))))*damp
          ut(npx-1, npy-1) = (uci(npx-1, npy-1)-r1b4*cosa_u(npx-1, npy-1&
&           )*(vt(npx-1, npy)+vt(npx-2, npy)+vt(npx-2, npy-1)+vci(npx-1&
&           , npy-1)-r1b4*cosa_v(npx-1, npy-1)*(ut(npx, npy-1)+ut(npx, &
&           npy-2)+ut(npx-1, npy-2))))*damp
          vt(npx-1, npy-1) = (vci(npx-1, npy-1)-r1b4*cosa_v(npx-1, npy-1&
&           )*(ut(npx, npy-1)+ut(npx, npy-2)+ut(npx-1, npy-2)+uci(npx-1&
&           , npy-1)-r1b4*cosa_u(npx-1, npy-1)*(vt(npx-1, npy)+vt(npx-2&
&           , npy)+vt(npx-2, npy-1))))*damp
          vt(npx, npy-1) = (vci(npx, npy-1)+r1b4*cosa_v(npx-1, npy-1)*(&
&           ut(npx, npy-1)+ut(npx, npy-2)+ut(npx+1, npy-2)+uci(npx+1, &
&           npy-1)+r1b4*cosa_u(npx-1, npy-1)*(vt(npx, npy)+vt(npx+1, npy&
&           )+vt(npx+1, npy-1))))*damp
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
        IF (nw_corner) THEN
          damp = 1.0/(1.0-r1b16*cosa_u(2, npy-1)*cosa_v(1, npy-1))
          ut(2, npy) = (uci(2, npy)+r1b4*cosa_u(2, npy-1)*(vt(1, npy)+vt&
&           (2, npy)+vt(2, npy+1)+vci(1, npy+1)+r1b4*cosa_v(1, npy-1)*(&
&           ut(1, npy)+ut(1, npy+1)+ut(2, npy+1))))*damp
          ut(2, npy-1) = (uci(2, npy-1)-r1b4*cosa_u(2, npy-1)*(vt(1, npy&
&           )+vt(2, npy)+vt(2, npy-1)+vci(1, npy-1)-r1b4*cosa_v(1, npy-1&
&           )*(ut(1, npy-1)+ut(1, npy-2)+ut(2, npy-2))))*damp
          vt(1, npy-1) = (vci(1, npy-1)-r1b4*cosa_v(1, npy-1)*(ut(1, npy&
&           -1)+ut(1, npy-2)+ut(2, npy-2)+uci(2, npy-1)-r1b4*cosa_u(2, &
&           npy-1)*(vt(1, npy)+vt(2, npy)+vt(2, npy-1))))*damp
          vt(0, npy-1) = (vci(0, npy-1)+r1b4*cosa_v(1, npy-1)*(ut(1, npy&
&           -1)+ut(1, npy-2)+ut(0, npy-2)+uci(0, npy-1)+r1b4*cosa_u(2, &
&           npy-1)*(vt(0, npy)+vt(-1, npy)+vt(-1, npy-1))))*damp
          CALL PUSHCONTROL2B(2)
        ELSE
          CALL PUSHCONTROL2B(1)
        END IF
      ELSE
! grid_type >= 3
        DO j=jsd,jed
          DO i=is-1,ie+2
            ut(i, j) = uc(i, j)
          END DO
        END DO
        DO j=js-1,je+2
          DO i=isd,ied
            vt(i, j) = vc(i, j)
          END DO
        END DO
        CALL PUSHCONTROL2B(0)
      END IF
      DO j=jsd,jed
        DO i=is,ie+1
          xfx_adv(i, j) = dt*ut(i, j)
        END DO
      END DO
      DO j=js,je+1
        DO i=isd,ied
          yfx_adv(i, j) = dt*vt(i, j)
        END DO
      END DO
! Compute E-W CFL number:
      DO j=jsd,jed
        DO i=is,ie+1
          IF (xfx_adv(i, j) .GT. 0.) THEN
            crx_adv(i, j) = xfx_adv(i, j)*rdxa(i-1, j)
            CALL PUSHCONTROL1B(1)
          ELSE
            crx_adv(i, j) = xfx_adv(i, j)*rdxa(i, j)
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
      END DO
      DO j=jsd,jed
        DO i=is,ie+1
          xfx_adv(i, j) = dy(i, j)*xfx_adv(i, j)*sina_u(i, j)
        END DO
      END DO
      DO j=js,je+1
        DO i=isd,ied
          IF (yfx_adv(i, j) .GT. 0.) THEN
            cry_adv(i, j) = yfx_adv(i, j)*rdya(i, j-1)
            CALL PUSHCONTROL1B(1)
          ELSE
            cry_adv(i, j) = yfx_adv(i, j)*rdya(i, j)
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
      END DO
      DO j=js,je+1
        DO i=isd,ied
          yfx_adv(i, j) = dx(i, j)*yfx_adv(i, j)*sina_v(i, j)
        END DO
      END DO
      CALL PUSHCONTROL1B(0)
    END IF
    DO j=jsd,jed
      DO i=is,ie
        ra_x(i, j) = area(i, j) + xfx_adv(i, j) - xfx_adv(i+1, j)
      END DO
    END DO
    DO j=js,je
      DO i=isd,ied
        ra_y(i, j) = area(i, j) + yfx_adv(i, j) - yfx_adv(i, j+1)
      END DO
    END DO
    CALL PUSHREAL8ARRAY(fy, (ie-is+1)*(je-js+2))
    CALL PUSHREAL8ARRAY(fx, (ie-is+2)*(je-js+1))
    CALL PUSHREAL8ARRAY(delp, (ied-isd+1)*(jed-jsd+1))
    CALL FV_TP_2D(delp, crx_adv, cry_adv, npx, npy, hord_dp, fx, fy, &
&           xfx_adv, yfx_adv, area, ra_x, ra_y)
! shallow_water
    IF (shallow_water) THEN
!#ifdef DOUBLE_GRADIENTS
!        dp_fx_r8 = fx
!        dp_fy_r8 = fy
!        do j=js,je
!           do i=is,ie
!              delp(i,j) = delp(i,j) + (dp_fx_r8(i,j)-dp_fx_r8(i+1,j)+dp_fy_r8(i,j)-dp_fy_r8(i,j+1))*rarea(i,j)
!              ptc(i,j) = pt(i,j)
!           enddo
!        enddo
!#else
      DO j=js,je
        DO i=is,ie
          CALL PUSHREAL8(ptc(i, j))
          ptc(i, j) = pt(i, j)
        END DO
      END DO
      CALL PUSHCONTROL2B(0)
    ELSE
      IF (.NOT.hydrostatic) THEN
        CALL PUSHREAL8ARRAY(py, (ie-is+1)*(je-js+2))
        CALL PUSHREAL8ARRAY(px, (ie-is+2)*(je-js+1))
        CALL PUSHREAL8ARRAY(w, (ied-isd+1)*(jed-jsd+1))
        CALL FV_TP_2D(w, crx_adv, cry_adv, npx, npy, hord_vt, px, py, &
&               xfx_adv, yfx_adv, area, ra_x, ra_y, mfx=fx, mfy=fy)
        DO j=js,je
          DO i=is,ie
            CALL PUSHREAL8(w(i, j))
            w(i, j) = w(i, j)*delp(i, j) + (px(i, j)-px(i+1, j)+py(i, j)&
&             -py(i, j+1))*rarea(i, j)
          END DO
        END DO
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (inline_q) THEN
        DO j=jsd,jed
          DO i=isd,ied
            CALL PUSHREAL8(pt(i, j))
            pt(i, j) = pt(i, j)/(1.+zvir*q(i, j, k, sphum))
          END DO
        END DO
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      CALL PUSHREAL8ARRAY(py, (ie-is+1)*(je-js+2))
      CALL PUSHREAL8ARRAY(px, (ie-is+2)*(je-js+1))
      CALL PUSHREAL8ARRAY(pt, (ied-isd+1)*(jed-jsd+1))
      CALL FV_TP_2D(pt, crx_adv, cry_adv, npx, npy, hord_tm, px, py, &
&             xfx_adv, yfx_adv, area, ra_x, ra_y, mfx=fx, mfy=fy)
      IF (inline_q) THEN
        DO j=js,je
          DO i=is,ie
            wk(i, j) = delp(i, j)
            CALL PUSHREAL8(delp(i, j))
            delp(i, j) = wk(i, j) + (fx(i, j)-fx(i+1, j)+fy(i, j)-fy(i, &
&             j+1))*rarea(i, j)
            CALL PUSHREAL8(pt(i, j))
            pt(i, j) = (pt(i, j)*wk(i, j)+(px(i, j)-px(i+1, j)+py(i, j)-&
&             py(i, j+1))*rarea(i, j))/delp(i, j)
          END DO
        END DO
        DO iq=1,nq
          CALL PUSHREAL8ARRAY(py, (ie-is+1)*(je-js+2))
          CALL PUSHREAL8ARRAY(px, (ie-is+2)*(je-js+1))
          CALL PUSHREAL8ARRAY(q(:, :, k, iq), (ied-isd+1)*(jed-jsd+1))
          CALL FV_TP_2D(q(isd:, jsd:, k, iq), crx_adv, cry_adv, npx, npy&
&                 , hord_tr, px, py, xfx_adv, yfx_adv, area, ra_x, ra_y&
&                 , mfx=fx, mfy=fy)
          DO j=js,je
            DO i=is,ie
              CALL PUSHREAL8(q(i, j, k, iq))
              q(i, j, k, iq) = (q(i, j, k, iq)*wk(i, j)+(px(i, j)-px(i+1&
&               , j)+py(i, j)-py(i, j+1))*rarea(i, j))/delp(i, j)
            END DO
          END DO
        END DO
        CALL PUSHCONTROL1B(0)
      ELSE
!#ifdef DOUBLE_GRADIENTS
!        pt_fx_r8 = px
!        pt_fy_r8 = py
!        dp_fx_r8 = fx
!        dp_fy_r8 = fy
!        do j=js,je
!           do i=is,ie
!              pt_dp = pt(i,j)
!              delp_dp = delp(i,j)
!              pt(i,j)   = pt_dp*delp_dp + (pt_fx_r8(i,j)-pt_fx_r8(i+1,j)+pt_fy_r8(i,j)-pt_fy_r8(i,j+1))*rarea(i,j)
!              delp(i,j) =       delp_dp + (dp_fx_r8(i,j)-dp_fx_r8(i+1,j)+dp_fy_r8(i,j)-dp_fy_r8(i,j+1))*rarea(i,j)
!              pt(i,j) = pt(i,j) / delp(i,j)
!           enddo
!        enddo
!#else
        DO j=js,je
          DO i=is,ie
            CALL PUSHREAL8(pt(i, j))
            pt(i, j) = pt(i, j)*delp(i, j) + (px(i, j)-px(i+1, j)+py(i, &
&             j)-py(i, j+1))*rarea(i, j)
            CALL PUSHREAL8(delp(i, j))
            delp(i, j) = delp(i, j) + (fx(i, j)-fx(i+1, j)+fy(i, j)-fy(i&
&             , j+1))*rarea(i, j)
          END DO
        END DO
        CALL PUSHCONTROL1B(1)
      END IF
!#endif
      IF (.NOT.hydrostatic) THEN
        CALL PUSHCONTROL2B(1)
      ELSE
        CALL PUSHCONTROL2B(2)
      END IF
    END IF
! shallow_water test_case 1
    IF (test_case .GT. 1) THEN
!----------------------
! Kinetic Energy Fluxes
!----------------------
! Compute B grid contra-variant components for KE:
      dt5 = 0.5*dt
      dt4 = 0.25*dt
      IF (2 .LT. is) THEN
        is2 = is
      ELSE
        is2 = 2
      END IF
      IF (npx - 1 .GT. ie + 1) THEN
        ie1 = ie + 1
      ELSE
        ie1 = npx - 1
      END IF
      IF (2 .LT. js) THEN
        js2 = js
      ELSE
        js2 = 2
      END IF
      IF (npy - 1 .GT. je + 1) THEN
        je1 = je + 1
      ELSE
        je1 = npy - 1
      END IF
      IF (grid_type .LT. 3) THEN
        IF (js .EQ. 1) THEN
          DO i=is,ie+1
! corner values are incorrect
            vb(i, 1) = dt5*(vt(i-1, 1)+vt(i, 1))
          END DO
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
        DO j=js2,je1
          DO i=is2,ie1
            vb(i, j) = dt5*(vc(i-1, j)+vc(i, j)-(uc(i, j-1)+uc(i, j))*&
&             cosa(i, j))*rsina(i, j)
          END DO
          IF (is .EQ. 1) THEN
!              vb(1,j) = dt5*(vt(0,j)+vt(1,j)) 
! 2-pt extrapolation from both sides:
            vb(1, j) = dt4*(-vt(-1, j)+3.*(vt(0, j)+vt(1, j))-vt(2, j))
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          IF (ie + 1 .EQ. npx) THEN
!              vb(npx,j) = dt5*(vt(npx-1,j)+vt(npx,j))
! 2-pt extrapolation from both sides:
            vb(npx, j) = dt4*(-vt(npx-2, j)+3.*(vt(npx-1, j)+vt(npx, j))&
&             -vt(npx+1, j))
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
        IF (je + 1 .EQ. npy) THEN
          DO i=is,ie+1
! corner values are incorrect
            vb(i, npy) = dt5*(vt(i-1, npy)+vt(i, npy))
          END DO
          CALL PUSHCONTROL2B(0)
        ELSE
          CALL PUSHCONTROL2B(1)
        END IF
      ELSE
        DO j=js,je+1
          DO i=is,ie+1
            vb(i, j) = dt5*(vc(i-1, j)+vc(i, j))
          END DO
        END DO
        CALL PUSHCONTROL2B(2)
      END IF
      CALL YTP_V(vb, u, v, ub, hord_mt)
      IF (grid_type .LT. 3) THEN
        IF (is .EQ. 1) THEN
          DO j=js,je+1
! corner values are incorrect
            CALL PUSHREAL8(ub(1, j))
            ub(1, j) = dt5*(ut(1, j-1)+ut(1, j))
          END DO
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
        DO j=js,je+1
          IF (j .EQ. 1 .OR. j .EQ. npy) THEN
            DO i=is2,ie1
!                 ub(i,j) = dt5*(ut(i,j-1)+ut(i,j))
! 2-pt extrapolation from both sides:
              CALL PUSHREAL8(ub(i, j))
              ub(i, j) = dt4*(-ut(i, j-2)+3.*(ut(i, j-1)+ut(i, j))-ut(i&
&               , j+1))
            END DO
            CALL PUSHCONTROL1B(1)
          ELSE
            DO i=is2,ie1
              CALL PUSHREAL8(ub(i, j))
              ub(i, j) = dt5*(uc(i, j-1)+uc(i, j)-(vc(i-1, j)+vc(i, j))*&
&               cosa(i, j))*rsina(i, j)
            END DO
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
        IF (ie + 1 .EQ. npx) THEN
          DO j=js,je+1
! corner values are incorrect
            CALL PUSHREAL8(ub(npx, j))
            ub(npx, j) = dt5*(ut(npx, j-1)+ut(npx, j))
          END DO
          CALL PUSHCONTROL2B(0)
        ELSE
          CALL PUSHCONTROL2B(1)
        END IF
      ELSE
        DO j=js,je+1
          DO i=is,ie+1
            CALL PUSHREAL8(ub(i, j))
            ub(i, j) = dt5*(uc(i, j-1)+uc(i, j))
          END DO
        END DO
        CALL PUSHCONTROL2B(2)
      END IF
      CALL PUSHREAL8ARRAY(vb, (ie-is+2)*(je-js+2))
      CALL XTP_U(ub, u, v, vb, hord_mt)
!-----------------------------------------
! Fix KE at the 4 corners of the face:
!-----------------------------------------
      IF (gnomonic_grid) THEN
        dt6 = dt/6.
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
          CALL PUSHCONTROL3B(4)
        ELSE
          CALL PUSHCONTROL3B(3)
        END IF
      ELSE IF (grid_type .LT. 3) THEN
!oncall mp_corner_comm(ke, npx, npy) 
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
          CALL PUSHCONTROL3B(2)
        ELSE
          CALL PUSHCONTROL3B(1)
        END IF
      ELSE
        CALL PUSHCONTROL3B(0)
      END IF
      CALL PUSHINTEGER4(j)
! Compute vorticity:
      DO j=jsd,jed+1
        CALL PUSHINTEGER4(i)
        DO i=isd,ied
          CALL PUSHREAL8(vt(i, j))
          vt(i, j) = u(i, j)*dx(i, j)
        END DO
      END DO
      DO j=jsd,jed
        CALL PUSHINTEGER4(i)
        DO i=isd,ied+1
          CALL PUSHREAL8(ut(i, j))
          ut(i, j) = v(i, j)*dy(i, j)
        END DO
      END DO
! wk is "volume-mean" relative vorticity
      DO j=jsd,jed
        CALL PUSHINTEGER4(i)
        DO i=isd,ied
          CALL PUSHREAL8(wk(i, j))
          wk(i, j) = rarea(i, j)*(vt(i, j)-vt(i, j+1)-ut(i, j)+ut(i+1, j&
&           ))
        END DO
      END DO
!-----------------------------
! Compute divergence damping
!-----------------------------
      IF (nord .EQ. 0) THEN
!       area ~ dxb*dyb*sin(alpha)
        DO j=js,je+1
          IF (j .EQ. 1 .OR. j .EQ. npy) THEN
            CALL PUSHINTEGER4(i)
            DO i=is-1,ie+1
              CALL PUSHREAL8(ptc(i, j))
              ptc(i, j) = u(i, j)*dyc(i, j)*sina_v(i, j)
            END DO
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHINTEGER4(i)
            DO i=is-1,ie+1
              CALL PUSHREAL8(ptc(i, j))
              ptc(i, j) = (u(i, j)-0.5*(va(i, j-1)+va(i, j))*cosa_v(i, j&
&               ))*dyc(i, j)*sina_v(i, j)
            END DO
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
        DO j=js-1,je+1
          CALL PUSHINTEGER4(i)
          DO i=is2,ie1
            vort(i, j) = (v(i, j)-0.5*(ua(i-1, j)+ua(i, j))*cosa_u(i, j)&
&             )*dxc(i, j)*sina_u(i, j)
          END DO
          IF (is .EQ. 1) THEN
            vort(1, j) = v(1, j)*dxc(1, j)*sina_u(1, j)
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          IF (ie + 1 .EQ. npx) THEN
            vort(npx, j) = v(npx, j)*dxc(npx, j)*sina_u(npx, j)
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
        DO j=js,je+1
          CALL PUSHINTEGER4(i)
          DO i=is,ie+1
            CALL PUSHREAL8(delpc(i, j))
            delpc(i, j) = vort(i, j-1) - vort(i, j) + ptc(i-1, j) - ptc(&
&             i, j)
          END DO
        END DO
! Remove the extra term at the corners:
        IF (sw_corner) THEN
          CALL PUSHREAL8(delpc(1, 1))
          delpc(1, 1) = delpc(1, 1) - vort(1, 0)
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
        IF (se_corner) THEN
          CALL PUSHREAL8(delpc(npx, 1))
          delpc(npx, 1) = delpc(npx, 1) - vort(npx, 0)
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
        IF (ne_corner) THEN
          CALL PUSHREAL8(delpc(npx, npy))
          delpc(npx, npy) = delpc(npx, npy) + vort(npx, npy)
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
        IF (nw_corner) THEN
          CALL PUSHREAL8(delpc(1, npy))
          delpc(1, npy) = delpc(1, npy) + vort(1, npy)
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
        DO j=js,je+1
          CALL PUSHINTEGER4(i)
          DO i=is,ie+1
            CALL PUSHREAL8(delpc(i, j))
            delpc(i, j) = rarea_c(i, j)*delpc(i, j)
            IF (delpc(i, j)*dt .GE. 0.) THEN
              abs0 = delpc(i, j)*dt
              CALL PUSHCONTROL1B(0)
            ELSE
              abs0 = -(delpc(i, j)*dt)
              CALL PUSHCONTROL1B(1)
            END IF
            y3 = dddmp*abs0
            IF (0.20 .GT. y3) THEN
              y1 = y3
              CALL PUSHCONTROL1B(0)
            ELSE
              y1 = 0.20
              CALL PUSHCONTROL1B(1)
            END IF
            IF (d2_bg .LT. y1) THEN
              max5 = y1
              CALL PUSHCONTROL1B(0)
            ELSE
              max5 = d2_bg
              CALL PUSHCONTROL1B(1)
            END IF
            CALL PUSHREAL8(damp)
            damp = da_min_c*max5
            vort(i, j) = damp*delpc(i, j)
          END DO
        END DO
        CALL PUSHCONTROL1B(0)
      ELSE
!--------------------------
! Higher order divg damping
!--------------------------
        DO j=js,je+1
          CALL PUSHINTEGER4(i)
          DO i=is,ie+1
! Save divergence for external mode filter
            CALL PUSHREAL8(delpc(i, j))
            delpc(i, j) = divg_d(i, j)
          END DO
        END DO
        n2 = nord + 1
        DO n=1,nord
          nt = nord - n
          fill_c = nt .NE. 0 .AND. grid_type .LT. 3 .AND. (((sw_corner &
&           .OR. se_corner) .OR. ne_corner) .OR. nw_corner)
!onif ( fill_c ) call fill_corners(divg_d, npx, npy, FILL=XDir, BGRID=.true.)
          IF (fill_c) THEN
            call fill_corners(divg_d, npx, npy, FILL=XDir, BGRID=.true.)
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
          ad_from0 = js - nt
          DO j=ad_from0,je+1+nt
            ad_from = is - 1 - nt
            CALL PUSHINTEGER4(i)
            DO i=ad_from,ie+1+nt
              vc(i, j) = (divg_d(i+1, j)-divg_d(i, j))*divg_u(i, j)
            END DO
            CALL PUSHINTEGER4(i - 1)
            CALL PUSHINTEGER4(ad_from)
          END DO
          CALL PUSHINTEGER4(j - 1)
          CALL PUSHINTEGER4(ad_from0)
!onif ( fill_c ) call fill_corners(divg_d, npx, npy, FILL=YDir, BGRID=.true.)
          IF (fill_c) THEN
            call fill_corners(divg_d, npx, npy, FILL=YDir, BGRID=.true.)
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
          ad_from2 = js - 1 - nt
          DO j=ad_from2,je+1+nt
            ad_from1 = is - nt
            CALL PUSHINTEGER4(i)
            DO i=ad_from1,ie+1+nt
              uc(i, j) = (divg_d(i, j+1)-divg_d(i, j))*divg_v(i, j)
            END DO
            CALL PUSHINTEGER4(i - 1)
            CALL PUSHINTEGER4(ad_from1)
          END DO
          CALL PUSHINTEGER4(j - 1)
          CALL PUSHINTEGER4(ad_from2)
!onif ( fill_c ) call fill_corners(vc, uc, npx, npy, VECTOR=.true., DGRID=.true.)
          IF (fill_c) THEN
            call fill_corners(vc, uc, npx, npy, VECTOR=.true., DGRID=.true.)
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
          ad_from4 = js - nt
          DO j=ad_from4,je+1+nt
            ad_from3 = is - nt
            CALL PUSHINTEGER4(i)
            DO i=ad_from3,ie+1+nt
              divg_d(i, j) = uc(i, j-1) - uc(i, j) + vc(i-1, j) - vc(i, &
&               j)
            END DO
            CALL PUSHINTEGER4(i - 1)
            CALL PUSHINTEGER4(ad_from3)
          END DO
          CALL PUSHINTEGER4(j - 1)
          CALL PUSHINTEGER4(ad_from4)
! Remove the extra term at the corners:
          IF (sw_corner) THEN
            divg_d(1, 1) = divg_d(1, 1) - uc(1, 0)
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          IF (se_corner) THEN
            divg_d(npx, 1) = divg_d(npx, 1) - uc(npx, 0)
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          IF (ne_corner) THEN
            divg_d(npx, npy) = divg_d(npx, npy) + uc(npx, npy)
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          IF (nw_corner) THEN
            divg_d(1, npy) = divg_d(1, npy) + uc(1, npy)
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
          ad_from6 = js - nt
          DO j=ad_from6,je+1+nt
            ad_from5 = is - nt
            CALL PUSHINTEGER4(i)
            DO i=ad_from5,ie+1+nt
              divg_d(i, j) = divg_d(i, j)*rarea_c(i, j)
            END DO
            CALL PUSHINTEGER4(i - 1)
            CALL PUSHINTEGER4(ad_from5)
          END DO
          CALL PUSHINTEGER4(j - 1)
          CALL PUSHINTEGER4(ad_from6)
        END DO
        IF (dddmp .LT. 1.e-5) THEN
          CALL PUSHCONTROL1B(0)
          vort = 0.
        ELSE
! Compute "time-scale" for del-2 background damping
!#ifdef FULL_SMAG
!        do j=js-1,je+1
!           do i=is-1,ie+1
!              vt2(i,j) = wk(i,j)**2
!           enddo
!        enddo
!        do j=js,je+1
!           do i=is,ie+1
!              vort(i,j) = dt*sqrt(delpc(i,j)**2 + 0.25*(vt2(i-1,j-1) +    &
!                                  vt2(i,j-1) + vt2(i-1,j) + vt2(i,j)))    
!           enddo
!        enddo
!        if (sw_corner) vort(1,1) = dt*sqrt( delpc(1,1)**2 +   &
!                       r3*(vt2(1,0) + vt2(0,1) + vt2(1,1)) )
!
!        if (se_corner) vort(npx,1) = dt*sqrt( delpc(npx,1)**2 +   &
!                       r3*(vt2(npx-1,0) + vt2(npx-1,1) + vt2(npx,1)) )
!
!        if (ne_corner) vort(npx,npy) = dt*sqrt( delpc(npx,npy)**2 +   &
!                       r3*(vt2(npx-1,npy-1) + vt2(npx,npy-1) + vt2(npx-1,npy)) )
!
!        if (nw_corner) vort(1,npy) = dt*sqrt( delpc(1,npy)**2 +   &
!                       r3*(vt2(0,npy-1) + vt2(1,npy-1) + vt2(1,npy)) )
!#else
          DO j=js,je+1
            CALL PUSHINTEGER4(i)
            DO i=is,ie+1
              IF (dt*delpc(i, j) .GE. 0.) THEN
                vort(i, j) = dt*delpc(i, j)
                CALL PUSHCONTROL1B(0)
              ELSE
                vort(i, j) = -(dt*delpc(i, j))
                CALL PUSHCONTROL1B(1)
              END IF
            END DO
          END DO
          CALL PUSHCONTROL1B(1)
        END IF
!#endif
        dd8 = (da_min_c*d4_bg)**n2
        DO j=js,je+1
          CALL PUSHINTEGER4(i)
          DO i=is,ie+1
            IF (0.20 .GT. dddmp*vort(i, j)) THEN
              y2 = dddmp*vort(i, j)
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHCONTROL1B(1)
              y2 = 0.20
            END IF
            IF (d2_bg .LT. y2) THEN
              max6 = y2
              CALL PUSHCONTROL1B(0)
            ELSE
              max6 = d2_bg
              CALL PUSHCONTROL1B(1)
            END IF
! del-2
            CALL PUSHREAL8(damp2)
            damp2 = da_min_c*max6
            vort(i, j) = damp2*delpc(i, j) + dd8*divg_d(i, j)
          END DO
        END DO
        CALL PUSHCONTROL1B(1)
      END IF
!----------------------------------
! Heating due to divergent damping:
!----------------------------------
      IF (d_con .GT. 1.e-5) THEN
!  damp = 0.5*0.25*d_con
        CALL PUSHREAL8(damp)
        damp = 0.25*d_con
        DO j=js,je+1
          CALL PUSHINTEGER4(i)
          DO i=is,ie
! du
            CALL PUSHREAL8(ub(i, j))
            ub(i, j) = (vort(i, j)-vort(i+1, j))*rdx(i, j)
! u*du
            gy(i, j) = u(i, j)*ub(i, j)
          END DO
        END DO
        DO j=js,je
          CALL PUSHINTEGER4(i)
          DO i=is,ie+1
! dv
            CALL PUSHREAL8(vb(i, j))
            vb(i, j) = (vort(i, j)-vort(i, j+1))*rdy(i, j)
! v*dv
            gx(i, j) = v(i, j)*vb(i, j)
          END DO
        END DO
        DO j=js,je
          CALL PUSHINTEGER4(i)
        END DO
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL PUSHCONTROL1B(0)
      END IF
! Vorticity transport
      DO j=jsd,jed
        CALL PUSHINTEGER4(i)
        DO i=isd,ied
          vort(i, j) = wk(i, j) + f0(i, j)
        END DO
      END DO
      DO j=js,je+1
        CALL PUSHINTEGER4(i)
      END DO
      DO j=js,je
        CALL PUSHINTEGER4(i)
      END DO
      ke_ad = 0.0_8
      ut_ad = 0.0_8
      fx_ad = 0.0_8
      DO j=je,js,-1
        DO i=ie+1,is,-1
          ut_ad(i, j) = ut_ad(i, j) + v_ad(i, j)
          ke_ad(i, j) = ke_ad(i, j) + v_ad(i, j)
          ke_ad(i, j+1) = ke_ad(i, j+1) - v_ad(i, j)
          fx_ad(i, j) = fx_ad(i, j) - v_ad(i, j)
          v_ad(i, j) = 0.0_8
        END DO
        CALL POPINTEGER4(i)
      END DO
      vt_ad = 0.0_8
      fy_ad = 0.0_8
      DO j=je+1,js,-1
        DO i=ie,is,-1
          vt_ad(i, j) = vt_ad(i, j) + u_ad(i, j)
          ke_ad(i, j) = ke_ad(i, j) + u_ad(i, j)
          fy_ad(i, j) = fy_ad(i, j) + u_ad(i, j)
          ke_ad(i+1, j) = ke_ad(i+1, j) - u_ad(i, j)
          u_ad(i, j) = 0.0_8
        END DO
        CALL POPINTEGER4(i)
      END DO
      vort_ad = 0.0_8
      ra_x_ad = 0.0_8
      ra_y_ad = 0.0_8
      CALL FV_TP_2D_ADM(vort, vort_ad, crx_adv, crx_adv_ad, cry_adv, &
&                 cry_adv_ad, npx, npy, hord_vt_pert, fx, fx_ad, fy, fy_ad, &
&                 xfx_adv, xfx_adv_ad, yfx_adv, yfx_adv_ad, area, ra_x, &
&                 ra_x_ad, ra_y, ra_y_ad, ppm_fac=ppm_limiter, nord=nord&
&                 , damp_c=vtdm4)
      wk_ad = 0.0_8
      DO j=jed,jsd,-1
        DO i=ied,isd,-1
          wk_ad(i, j) = wk_ad(i, j) + vort_ad(i, j)
          vort_ad(i, j) = 0.0_8
        END DO
        CALL POPINTEGER4(i)
      END DO
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        ub_ad = 0.0_8
        vb_ad = 0.0_8
      ELSE
        gx_ad = 0.0_8
        gy_ad = 0.0_8
        ub_ad = 0.0_8
        vb_ad = 0.0_8
        DO j=je,js,-1
          DO i=ie,is,-1
            dv2 = vb(i, j) + vb(i+1, j)
            v2 = v(i, j) + v(i+1, j)
            du2 = ub(i, j) + ub(i, j+1)
            u2 = u(i, j) + u(i, j+1)
            temp_ad88 = -(damp*rsin2(i, j)*pt_ad(i, j)/pkz(i, j))
            temp_ad89 = 2.*temp_ad88
            temp_ad90 = -(cosa_s(i, j)*temp_ad88)
            ub_ad(i, j) = ub_ad(i, j) + 2*ub(i, j)*temp_ad88
            ub_ad(i, j+1) = ub_ad(i, j+1) + 2*ub(i, j+1)*temp_ad88
            vb_ad(i, j) = vb_ad(i, j) + 2*vb(i, j)*temp_ad88
            vb_ad(i+1, j) = vb_ad(i+1, j) + 2*vb(i+1, j)*temp_ad88
            gy_ad(i, j) = gy_ad(i, j) + temp_ad89
            gy_ad(i, j+1) = gy_ad(i, j+1) + temp_ad89
            gx_ad(i, j) = gx_ad(i, j) + temp_ad89
            gx_ad(i+1, j) = gx_ad(i+1, j) + temp_ad89
            u2_ad = dv2*temp_ad90
            dv2_ad = (du2+u2)*temp_ad90
            v2_ad = du2*temp_ad90
            du2_ad = (dv2+v2)*temp_ad90
            pkz_ad(i, j) = pkz_ad(i, j) - (ub(i, j)**2+ub(i, j+1)**2+vb(&
&             i, j)**2+vb(i+1, j)**2+2.*(gy(i, j)+gy(i, j+1)+gx(i, j)+gx&
&             (i+1, j))-cosa_s(i, j)*(u2*dv2+v2*du2+du2*dv2))*temp_ad88/&
&             pkz(i, j)
            vb_ad(i, j) = vb_ad(i, j) + dv2_ad
            vb_ad(i+1, j) = vb_ad(i+1, j) + dv2_ad
            v_ad(i, j) = v_ad(i, j) + v2_ad
            v_ad(i+1, j) = v_ad(i+1, j) + v2_ad
            ub_ad(i, j) = ub_ad(i, j) + du2_ad
            ub_ad(i, j+1) = ub_ad(i, j+1) + du2_ad
            u_ad(i, j) = u_ad(i, j) + u2_ad
            u_ad(i, j+1) = u_ad(i, j+1) + u2_ad
          END DO
          CALL POPINTEGER4(i)
        END DO
        DO j=je,js,-1
          DO i=ie+1,is,-1
            v_ad(i, j) = v_ad(i, j) + vb(i, j)*gx_ad(i, j)
            vb_ad(i, j) = vb_ad(i, j) + v(i, j)*gx_ad(i, j)
            gx_ad(i, j) = 0.0_8
            CALL POPREAL8(vb(i, j))
            temp_ad87 = rdy(i, j)*vb_ad(i, j)
            vort_ad(i, j) = vort_ad(i, j) + temp_ad87
            vort_ad(i, j+1) = vort_ad(i, j+1) - temp_ad87
            vb_ad(i, j) = 0.0_8
          END DO
          CALL POPINTEGER4(i)
        END DO
        DO j=je+1,js,-1
          DO i=ie,is,-1
            u_ad(i, j) = u_ad(i, j) + ub(i, j)*gy_ad(i, j)
            ub_ad(i, j) = ub_ad(i, j) + u(i, j)*gy_ad(i, j)
            gy_ad(i, j) = 0.0_8
            CALL POPREAL8(ub(i, j))
            temp_ad86 = rdx(i, j)*ub_ad(i, j)
            vort_ad(i, j) = vort_ad(i, j) + temp_ad86
            vort_ad(i+1, j) = vort_ad(i+1, j) - temp_ad86
            ub_ad(i, j) = 0.0_8
          END DO
          CALL POPINTEGER4(i)
        END DO
        CALL POPREAL8(damp)
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO j=je+1,js,-1
          DO i=ie+1,is,-1
            vort_ad(i, j) = vort_ad(i, j) + ke_ad(i, j)
            damp_ad = delpc(i, j)*vort_ad(i, j)
            delpc_ad(i, j) = delpc_ad(i, j) + damp*vort_ad(i, j)
            vort_ad(i, j) = 0.0_8
            CALL POPREAL8(damp)
            max5_ad = da_min_c*damp_ad
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              y1_ad = max5_ad
            ELSE
              y1_ad = 0.0_8
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              y3_ad = y1_ad
            ELSE
              y3_ad = 0.0_8
            END IF
            abs0_ad = dddmp*y3_ad
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              delpc_ad(i, j) = delpc_ad(i, j) + dt*abs0_ad
            ELSE
              delpc_ad(i, j) = delpc_ad(i, j) - dt*abs0_ad
            END IF
            CALL POPREAL8(delpc(i, j))
            delpc_ad(i, j) = rarea_c(i, j)*delpc_ad(i, j)
          END DO
          CALL POPINTEGER4(i)
        END DO
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) THEN
          CALL POPREAL8(delpc(1, npy))
          vort_ad(1, npy) = vort_ad(1, npy) + delpc_ad(1, npy)
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          CALL POPREAL8(delpc(npx, npy))
          vort_ad(npx, npy) = vort_ad(npx, npy) + delpc_ad(npx, npy)
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          CALL POPREAL8(delpc(npx, 1))
          vort_ad(npx, 0) = vort_ad(npx, 0) - delpc_ad(npx, 1)
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          CALL POPREAL8(delpc(1, 1))
          vort_ad(1, 0) = vort_ad(1, 0) - delpc_ad(1, 1)
        END IF
        DO j=je+1,js,-1
          DO i=ie+1,is,-1
            CALL POPREAL8(delpc(i, j))
            vort_ad(i, j-1) = vort_ad(i, j-1) + delpc_ad(i, j)
            vort_ad(i, j) = vort_ad(i, j) - delpc_ad(i, j)
            ptc_ad(i-1, j) = ptc_ad(i-1, j) + delpc_ad(i, j)
            ptc_ad(i, j) = ptc_ad(i, j) - delpc_ad(i, j)
            delpc_ad(i, j) = 0.0_8
          END DO
          CALL POPINTEGER4(i)
        END DO
        DO j=je+1,js-1,-1
          CALL POPCONTROL1B(branch)
          IF (branch .NE. 0) THEN
            v_ad(npx, j) = v_ad(npx, j) + dxc(npx, j)*sina_u(npx, j)*&
&             vort_ad(npx, j)
            vort_ad(npx, j) = 0.0_8
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            v_ad(1, j) = v_ad(1, j) + dxc(1, j)*sina_u(1, j)*vort_ad(1, &
&             j)
            vort_ad(1, j) = 0.0_8
          END IF
          DO i=ie1,is2,-1
            temp_ad82 = dxc(i, j)*sina_u(i, j)*vort_ad(i, j)
            temp_ad83 = -(cosa_u(i, j)*0.5*temp_ad82)
            v_ad(i, j) = v_ad(i, j) + temp_ad82
            ua_ad(i-1, j) = ua_ad(i-1, j) + temp_ad83
            ua_ad(i, j) = ua_ad(i, j) + temp_ad83
            vort_ad(i, j) = 0.0_8
          END DO
          CALL POPINTEGER4(i)
        END DO
        DO j=je+1,js,-1
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            DO i=ie+1,is-1,-1
              CALL POPREAL8(ptc(i, j))
              temp_ad80 = dyc(i, j)*sina_v(i, j)*ptc_ad(i, j)
              temp_ad81 = -(cosa_v(i, j)*0.5*temp_ad80)
              u_ad(i, j) = u_ad(i, j) + temp_ad80
              va_ad(i, j-1) = va_ad(i, j-1) + temp_ad81
              va_ad(i, j) = va_ad(i, j) + temp_ad81
              ptc_ad(i, j) = 0.0_8
            END DO
            CALL POPINTEGER4(i)
          ELSE
            DO i=ie+1,is-1,-1
              CALL POPREAL8(ptc(i, j))
              u_ad(i, j) = u_ad(i, j) + dyc(i, j)*sina_v(i, j)*ptc_ad(i&
&               , j)
              ptc_ad(i, j) = 0.0_8
            END DO
            CALL POPINTEGER4(i)
          END IF
        END DO
      ELSE
        DO j=je+1,js,-1
          DO i=ie+1,is,-1
            vort_ad(i, j) = vort_ad(i, j) + ke_ad(i, j)
            damp2_ad = delpc(i, j)*vort_ad(i, j)
            delpc_ad(i, j) = delpc_ad(i, j) + damp2*vort_ad(i, j)
            divg_d_ad(i, j) = divg_d_ad(i, j) + dd8*vort_ad(i, j)
            vort_ad(i, j) = 0.0_8
            CALL POPREAL8(damp2)
            max6_ad = da_min_c*damp2_ad
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              y2_ad = max6_ad
            ELSE
              y2_ad = 0.0_8
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) vort_ad(i, j) = vort_ad(i, j) + dddmp*&
&               y2_ad
          END DO
          CALL POPINTEGER4(i)
        END DO
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) THEN
          DO j=je+1,js,-1
            DO i=ie+1,is,-1
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                delpc_ad(i, j) = delpc_ad(i, j) + dt*vort_ad(i, j)
                vort_ad(i, j) = 0.0_8
              ELSE
                delpc_ad(i, j) = delpc_ad(i, j) - dt*vort_ad(i, j)
                vort_ad(i, j) = 0.0_8
              END IF
            END DO
            CALL POPINTEGER4(i)
          END DO
        END IF
        DO n=nord,1,-1
          CALL POPINTEGER4(ad_from6)
          CALL POPINTEGER4(ad_to6)
          DO j=ad_to6,ad_from6,-1
            CALL POPINTEGER4(ad_from5)
            CALL POPINTEGER4(ad_to5)
            DO i=ad_to5,ad_from5,-1
              divg_d_ad(i, j) = rarea_c(i, j)*divg_d_ad(i, j)
            END DO
            CALL POPINTEGER4(i)
          END DO
          CALL POPCONTROL1B(branch)
          IF (branch .NE. 0) uc_ad(1, npy) = uc_ad(1, npy) + divg_d_ad(1&
&             , npy)
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) uc_ad(npx, npy) = uc_ad(npx, npy) + &
&             divg_d_ad(npx, npy)
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) uc_ad(npx, 0) = uc_ad(npx, 0) - divg_d_ad(&
&             npx, 1)
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) uc_ad(1, 0) = uc_ad(1, 0) - divg_d_ad(1, 1)
          CALL POPINTEGER4(ad_from4)
          CALL POPINTEGER4(ad_to4)
          DO j=ad_to4,ad_from4,-1
            CALL POPINTEGER4(ad_from3)
            CALL POPINTEGER4(ad_to3)
            DO i=ad_to3,ad_from3,-1
              uc_ad(i, j-1) = uc_ad(i, j-1) + divg_d_ad(i, j)
              uc_ad(i, j) = uc_ad(i, j) - divg_d_ad(i, j)
              vc_ad(i-1, j) = vc_ad(i-1, j) + divg_d_ad(i, j)
              vc_ad(i, j) = vc_ad(i, j) - divg_d_ad(i, j)
              divg_d_ad(i, j) = 0.0_8
            END DO
            CALL POPINTEGER4(i)
          END DO
          CALL POPCONTROL1B(branch)
          IF (branch .NE. 0) THEN
             call fill_corners(vc_ad, uc_ad, npx, npy, VECTOR=.true., DGRID=.true.)
          END IF
          CALL POPINTEGER4(ad_from2)
          CALL POPINTEGER4(ad_to2)
          DO j=ad_to2,ad_from2,-1
            CALL POPINTEGER4(ad_from1)
            CALL POPINTEGER4(ad_to1)
            DO i=ad_to1,ad_from1,-1
              temp_ad85 = divg_v(i, j)*uc_ad(i, j)
              divg_d_ad(i, j+1) = divg_d_ad(i, j+1) + temp_ad85
              divg_d_ad(i, j) = divg_d_ad(i, j) - temp_ad85
              uc_ad(i, j) = 0.0_8
            END DO
            CALL POPINTEGER4(i)
          END DO
          CALL POPCONTROL1B(branch)
          IF (branch .NE. 0) call fill_corners(divg_d_ad, npx, npy, FILL=YDir, BGRID=.true.)
          CALL POPINTEGER4(ad_from0)
          CALL POPINTEGER4(ad_to0)
          DO j=ad_to0,ad_from0,-1
            CALL POPINTEGER4(ad_from)
            CALL POPINTEGER4(ad_to)
            DO i=ad_to,ad_from,-1
              temp_ad84 = divg_u(i, j)*vc_ad(i, j)
              divg_d_ad(i+1, j) = divg_d_ad(i+1, j) + temp_ad84
              divg_d_ad(i, j) = divg_d_ad(i, j) - temp_ad84
              vc_ad(i, j) = 0.0_8
            END DO
            CALL POPINTEGER4(i)
          END DO
          CALL POPCONTROL1B(branch)
          IF (branch .NE. 0) call fill_corners(divg_d_ad, npx, npy, FILL=XDir, BGRID=.true.)
        END DO
        DO j=je+1,js,-1
          DO i=ie+1,is,-1
            CALL POPREAL8(delpc(i, j))
            divg_d_ad(i, j) = divg_d_ad(i, j) + delpc_ad(i, j)
            delpc_ad(i, j) = 0.0_8
          END DO
          CALL POPINTEGER4(i)
        END DO
      END IF
      DO j=jed,jsd,-1
        DO i=ied,isd,-1
          CALL POPREAL8(wk(i, j))
          temp_ad79 = rarea(i, j)*wk_ad(i, j)
          vt_ad(i, j) = vt_ad(i, j) + temp_ad79
          vt_ad(i, j+1) = vt_ad(i, j+1) - temp_ad79
          ut_ad(i+1, j) = ut_ad(i+1, j) + temp_ad79
          ut_ad(i, j) = ut_ad(i, j) - temp_ad79
          wk_ad(i, j) = 0.0_8
        END DO
        CALL POPINTEGER4(i)
      END DO
      DO j=jed,jsd,-1
        DO i=ied+1,isd,-1
          CALL POPREAL8(ut(i, j))
          v_ad(i, j) = v_ad(i, j) + dy(i, j)*ut_ad(i, j)
          ut_ad(i, j) = 0.0_8
        END DO
        CALL POPINTEGER4(i)
      END DO
      DO j=jed+1,jsd,-1
        DO i=ied,isd,-1
          CALL POPREAL8(vt(i, j))
          u_ad(i, j) = u_ad(i, j) + dx(i, j)*vt_ad(i, j)
          vt_ad(i, j) = 0.0_8
        END DO
        CALL POPINTEGER4(i)
      END DO
      CALL POPINTEGER4(j)
      CALL POPCONTROL3B(branch)
      IF (branch .LT. 2) THEN
        IF (branch .EQ. 0) GOTO 100
      ELSE IF (branch .EQ. 2) THEN
        temp_ad78 = r3*ke_ad(1, npy)
        ke_ad(2, npy) = ke_ad(2, npy) + temp_ad78
        ke_ad(1, npy-1) = ke_ad(1, npy-1) + temp_ad78
        ke_ad(0, npy) = ke_ad(0, npy) + temp_ad78
        ke_ad(1, npy) = 0.0_8
      ELSE
        IF (branch .NE. 3) THEN
          j = npy
          temp_ad71 = dt6*ke_ad(1, j)
          temp_ad72 = u(1, j)*temp_ad71
          temp_ad73 = v(1, j-1)*temp_ad71
          temp_ad74 = u(0, j)*temp_ad71
          ut_ad(1, j) = ut_ad(1, j) + temp_ad72
          ut_ad(1, j-1) = ut_ad(1, j-1) + temp_ad74 + temp_ad72
          u_ad(1, j) = u_ad(1, j) + (ut(1, j)+ut(1, j-1))*temp_ad71
          vt_ad(1, j) = vt_ad(1, j) + temp_ad73 - temp_ad74
          vt_ad(0, j) = vt_ad(0, j) + temp_ad73
          v_ad(1, j-1) = v_ad(1, j-1) + (vt(1, j)+vt(0, j))*temp_ad71
          u_ad(0, j) = u_ad(0, j) + (ut(1, j-1)-vt(1, j))*temp_ad71
          ke_ad(1, j) = 0.0_8
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          i = npx
          j = npy
          temp_ad67 = dt6*ke_ad(i, j)
          temp_ad68 = u(i-1, j)*temp_ad67
          temp_ad69 = v(i, j-1)*temp_ad67
          temp_ad70 = u(i, j)*temp_ad67
          ut_ad(i, j) = ut_ad(i, j) + temp_ad68
          ut_ad(i, j-1) = ut_ad(i, j-1) + temp_ad70 + temp_ad68
          u_ad(i-1, j) = u_ad(i-1, j) + (ut(i, j)+ut(i, j-1))*temp_ad67
          vt_ad(i, j) = vt_ad(i, j) + temp_ad69
          vt_ad(i-1, j) = vt_ad(i-1, j) + temp_ad70 + temp_ad69
          v_ad(i, j-1) = v_ad(i, j-1) + (vt(i, j)+vt(i-1, j))*temp_ad67
          u_ad(i, j) = u_ad(i, j) + (ut(i, j-1)+vt(i-1, j))*temp_ad67
          ke_ad(i, j) = 0.0_8
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          i = npx
          temp_ad63 = dt6*ke_ad(i, 1)
          temp_ad64 = u(i-1, 1)*temp_ad63
          temp_ad65 = v(i, 1)*temp_ad63
          temp_ad66 = u(i, 1)*temp_ad63
          ut_ad(i, 1) = ut_ad(i, 1) + temp_ad66 + temp_ad64
          ut_ad(i, 0) = ut_ad(i, 0) + temp_ad64
          u_ad(i-1, 1) = u_ad(i-1, 1) + (ut(i, 1)+ut(i, 0))*temp_ad63
          vt_ad(i, 1) = vt_ad(i, 1) + temp_ad65
          vt_ad(i-1, 1) = vt_ad(i-1, 1) + temp_ad65 - temp_ad66
          v_ad(i, 1) = v_ad(i, 1) + (vt(i, 1)+vt(i-1, 1))*temp_ad63
          u_ad(i, 1) = u_ad(i, 1) + (ut(i, 1)-vt(i-1, 1))*temp_ad63
          ke_ad(i, 1) = 0.0_8
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          temp_ad62 = dt6*ke_ad(1, 1)
          ut_ad(1, 1) = ut_ad(1, 1) + (u(0, 1)+u(1, 1))*temp_ad62
          ut_ad(1, 0) = ut_ad(1, 0) + u(1, 1)*temp_ad62
          u_ad(1, 1) = u_ad(1, 1) + (ut(1, 1)+ut(1, 0))*temp_ad62
          vt_ad(1, 1) = vt_ad(1, 1) + (u(0, 1)+v(1, 1))*temp_ad62
          vt_ad(0, 1) = vt_ad(0, 1) + v(1, 1)*temp_ad62
          v_ad(1, 1) = v_ad(1, 1) + (vt(1, 1)+vt(0, 1))*temp_ad62
          u_ad(0, 1) = u_ad(0, 1) + (ut(1, 1)+vt(1, 1))*temp_ad62
          ke_ad(1, 1) = 0.0_8
        END IF
        GOTO 100
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        temp_ad77 = r3*ke_ad(npx, npy)
        ke_ad(npx+1, npy) = ke_ad(npx+1, npy) + temp_ad77
        ke_ad(npx, npy-1) = ke_ad(npx, npy-1) + temp_ad77
        ke_ad(npx-1, npy) = ke_ad(npx-1, npy) + temp_ad77
        ke_ad(npx, npy) = 0.0_8
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        temp_ad76 = r3*ke_ad(npx, 1)
        ke_ad(npx+1, 1) = ke_ad(npx+1, 1) + temp_ad76
        ke_ad(npx, 2) = ke_ad(npx, 2) + temp_ad76
        ke_ad(npx-1, 1) = ke_ad(npx-1, 1) + temp_ad76
        ke_ad(npx, 1) = 0.0_8
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        temp_ad75 = r3*ke_ad(1, 1)
        ke_ad(2, 1) = ke_ad(2, 1) + temp_ad75
        ke_ad(1, 2) = ke_ad(1, 2) + temp_ad75
        ke_ad(0, 1) = ke_ad(0, 1) + temp_ad75
        ke_ad(1, 1) = 0.0_8
      END IF
      call mp_corner_comm(ke_ad, npx, npy)
 100  DO j=je+1,js,-1
        DO i=ie+1,is,-1
          temp_ad61 = 0.5*ke_ad(i, j)
          ub_ad(i, j) = ub_ad(i, j) + vb(i, j)*temp_ad61
          vb_ad(i, j) = vb_ad(i, j) + ub(i, j)*temp_ad61
          ke_ad(i, j) = temp_ad61
        END DO
      END DO
      CALL POPREAL8ARRAY(vb, (ie-is+2)*(je-js+2))
      CALL XTP_U_ADM(ub, ub_ad, u, u_ad, v, vb, vb_ad, hord_mt_pert)
      CALL POPCONTROL2B(branch)
      IF (branch .EQ. 0) THEN
        DO j=je+1,js,-1
          CALL POPREAL8(ub(npx, j))
          ut_ad(npx, j-1) = ut_ad(npx, j-1) + dt5*ub_ad(npx, j)
          ut_ad(npx, j) = ut_ad(npx, j) + dt5*ub_ad(npx, j)
          ub_ad(npx, j) = 0.0_8
        END DO
      ELSE IF (branch .NE. 1) THEN
        DO j=je+1,js,-1
          DO i=ie+1,is,-1
            CALL POPREAL8(ub(i, j))
            uc_ad(i, j-1) = uc_ad(i, j-1) + dt5*ub_ad(i, j)
            uc_ad(i, j) = uc_ad(i, j) + dt5*ub_ad(i, j)
            ub_ad(i, j) = 0.0_8
          END DO
        END DO
        GOTO 110
      END IF
      DO j=je+1,js,-1
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          DO i=ie1,is2,-1
            CALL POPREAL8(ub(i, j))
            temp_ad59 = dt5*rsina(i, j)*ub_ad(i, j)
            temp_ad60 = -(cosa(i, j)*temp_ad59)
            uc_ad(i, j-1) = uc_ad(i, j-1) + temp_ad59
            uc_ad(i, j) = uc_ad(i, j) + temp_ad59
            vc_ad(i-1, j) = vc_ad(i-1, j) + temp_ad60
            vc_ad(i, j) = vc_ad(i, j) + temp_ad60
            ub_ad(i, j) = 0.0_8
          END DO
        ELSE
          DO i=ie1,is2,-1
            CALL POPREAL8(ub(i, j))
            temp_ad58 = dt4*ub_ad(i, j)
            ut_ad(i, j-1) = ut_ad(i, j-1) + 3.*temp_ad58
            ut_ad(i, j) = ut_ad(i, j) + 3.*temp_ad58
            ut_ad(i, j-2) = ut_ad(i, j-2) - temp_ad58
            ut_ad(i, j+1) = ut_ad(i, j+1) - temp_ad58
            ub_ad(i, j) = 0.0_8
          END DO
        END IF
      END DO
      CALL POPCONTROL1B(branch)
      IF (branch .NE. 0) THEN
        DO j=je+1,js,-1
          CALL POPREAL8(ub(1, j))
          ut_ad(1, j-1) = ut_ad(1, j-1) + dt5*ub_ad(1, j)
          ut_ad(1, j) = ut_ad(1, j) + dt5*ub_ad(1, j)
          ub_ad(1, j) = 0.0_8
        END DO
      END IF
 110  DO j=je+1,js,-1
        DO i=ie+1,is,-1
          vb_ad(i, j) = vb_ad(i, j) + ub(i, j)*ke_ad(i, j)
          ub_ad(i, j) = ub_ad(i, j) + vb(i, j)*ke_ad(i, j)
          ke_ad(i, j) = 0.0_8
        END DO
      END DO
      CALL YTP_V_ADM(vb, vb_ad, u, v, v_ad, ub, ub_ad, hord_mt_pert)
      CALL POPCONTROL2B(branch)
      IF (branch .EQ. 0) THEN
        DO i=ie+1,is,-1
          vt_ad(i-1, npy) = vt_ad(i-1, npy) + dt5*vb_ad(i, npy)
          vt_ad(i, npy) = vt_ad(i, npy) + dt5*vb_ad(i, npy)
          vb_ad(i, npy) = 0.0_8
        END DO
      ELSE IF (branch .NE. 1) THEN
        DO j=je+1,js,-1
          DO i=ie+1,is,-1
            vc_ad(i-1, j) = vc_ad(i-1, j) + dt5*vb_ad(i, j)
            vc_ad(i, j) = vc_ad(i, j) + dt5*vb_ad(i, j)
            vb_ad(i, j) = 0.0_8
          END DO
        END DO
        GOTO 120
      END IF
      DO j=je1,js2,-1
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) THEN
          temp_ad57 = dt4*vb_ad(npx, j)
          vt_ad(npx-1, j) = vt_ad(npx-1, j) + 3.*temp_ad57
          vt_ad(npx, j) = vt_ad(npx, j) + 3.*temp_ad57
          vt_ad(npx-2, j) = vt_ad(npx-2, j) - temp_ad57
          vt_ad(npx+1, j) = vt_ad(npx+1, j) - temp_ad57
          vb_ad(npx, j) = 0.0_8
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          temp_ad56 = dt4*vb_ad(1, j)
          vt_ad(0, j) = vt_ad(0, j) + 3.*temp_ad56
          vt_ad(1, j) = vt_ad(1, j) + 3.*temp_ad56
          vt_ad(-1, j) = vt_ad(-1, j) - temp_ad56
          vt_ad(2, j) = vt_ad(2, j) - temp_ad56
          vb_ad(1, j) = 0.0_8
        END IF
        DO i=ie1,is2,-1
          temp_ad54 = dt5*rsina(i, j)*vb_ad(i, j)
          temp_ad55 = -(cosa(i, j)*temp_ad54)
          vc_ad(i-1, j) = vc_ad(i-1, j) + temp_ad54
          vc_ad(i, j) = vc_ad(i, j) + temp_ad54
          uc_ad(i, j-1) = uc_ad(i, j-1) + temp_ad55
          uc_ad(i, j) = uc_ad(i, j) + temp_ad55
          vb_ad(i, j) = 0.0_8
        END DO
      END DO
      CALL POPCONTROL1B(branch)
      IF (branch .NE. 0) THEN
        DO i=ie+1,is,-1
          vt_ad(i-1, 1) = vt_ad(i-1, 1) + dt5*vb_ad(i, 1)
          vt_ad(i, 1) = vt_ad(i, 1) + dt5*vb_ad(i, 1)
          vb_ad(i, 1) = 0.0_8
        END DO
      END IF
    ELSE
      ut_ad = 0.0_8
      ra_x_ad = 0.0_8
      ra_y_ad = 0.0_8
      vt_ad = 0.0_8
      fx_ad = 0.0_8
      fy_ad = 0.0_8
      wk_ad = 0.0_8
    END IF
 120 CALL POPCONTROL2B(branch)
    IF (branch .EQ. 0) THEN
      DO j=je,js,-1
        DO i=ie,is,-1
          CALL POPREAL8(ptc(i, j))
          pt_ad(i, j) = pt_ad(i, j) + ptc_ad(i, j)
          ptc_ad(i, j) = 0.0_8
          temp_ad43 = rarea(i, j)*delp_ad(i, j)
          fx_ad(i, j) = fx_ad(i, j) + temp_ad43
          fx_ad(i+1, j) = fx_ad(i+1, j) - temp_ad43
          fy_ad(i, j) = fy_ad(i, j) + temp_ad43
          fy_ad(i, j+1) = fy_ad(i, j+1) - temp_ad43
        END DO
      END DO
    ELSE
      IF (branch .EQ. 1) THEN
        DO j=je,js,-1
          DO i=ie,is,-1
            temp_ad53 = w_ad(i, j)/delp(i, j)
            delp_ad(i, j) = delp_ad(i, j) - w(i, j)*temp_ad53/delp(i, j)
            w_ad(i, j) = temp_ad53
          END DO
        END DO
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO j=je,js,-1
          DO i=ie,is,-1
            q_ad(i, j, k, sphum) = q_ad(i, j, k, sphum) + pt(i, j)*zvir*&
&             pt_ad(i, j)
            pt_ad(i, j) = (zvir*q(i, j, k, sphum)+1.)*pt_ad(i, j)
          END DO
        END DO
        px_ad = 0.0_8
        py_ad = 0.0_8
        DO iq=nq,1,-1
          DO j=je,js,-1
            DO i=ie,is,-1
              CALL POPREAL8(q(i, j, k, iq))
              temp_ad48 = q_ad(i, j, k, iq)/delp(i, j)
              temp0 = q(i, j, k, iq)
              temp_ad49 = rarea(i, j)*temp_ad48
              wk_ad(i, j) = wk_ad(i, j) + temp0*temp_ad48
              px_ad(i, j) = px_ad(i, j) + temp_ad49
              px_ad(i+1, j) = px_ad(i+1, j) - temp_ad49
              py_ad(i, j) = py_ad(i, j) + temp_ad49
              py_ad(i, j+1) = py_ad(i, j+1) - temp_ad49
              delp_ad(i, j) = delp_ad(i, j) - (temp0*wk(i, j)+rarea(i, j&
&               )*(px(i, j)-px(i+1, j)+py(i, j)-py(i, j+1)))*temp_ad48/&
&               delp(i, j)
              q_ad(i, j, k, iq) = wk(i, j)*temp_ad48
            END DO
          END DO
          CALL POPREAL8ARRAY(q(:, :, k, iq), (ied-isd+1)*(jed-jsd+1))
          CALL POPREAL8ARRAY(px, (ie-is+2)*(je-js+1))
          CALL POPREAL8ARRAY(py, (ie-is+1)*(je-js+2))
          CALL FV_TP_2D_ADM(q(isd:, jsd:, k, iq), q_ad(isd:, jsd:, k, iq&
&                     ), crx_adv, crx_adv_ad, cry_adv, cry_adv_ad, npx, &
&                     npy, hord_tr_pert, px, px_ad, py, py_ad, xfx_adv, &
&                     xfx_adv_ad, yfx_adv, yfx_adv_ad, area, ra_x, &
&                     ra_x_ad, ra_y, ra_y_ad, mfx=fx, mfx_ad=fx_ad, mfy=&
&                     fy, mfy_ad=fy_ad)
        END DO
        DO j=je,js,-1
          DO i=ie,is,-1
            CALL POPREAL8(pt(i, j))
            temp_ad45 = pt_ad(i, j)/delp(i, j)
            temp_ad46 = rarea(i, j)*temp_ad45
            px_ad(i, j) = px_ad(i, j) + temp_ad46
            px_ad(i+1, j) = px_ad(i+1, j) - temp_ad46
            py_ad(i, j) = py_ad(i, j) + temp_ad46
            py_ad(i, j+1) = py_ad(i, j+1) - temp_ad46
            delp_ad(i, j) = delp_ad(i, j) - (pt(i, j)*wk(i, j)+rarea(i, &
&             j)*(px(i, j)-px(i+1, j)+py(i, j)-py(i, j+1)))*temp_ad45/&
&             delp(i, j)
            wk_ad(i, j) = wk_ad(i, j) + delp_ad(i, j) + pt(i, j)*&
&             temp_ad45
            pt_ad(i, j) = wk(i, j)*temp_ad45
            CALL POPREAL8(delp(i, j))
            temp_ad47 = rarea(i, j)*delp_ad(i, j)
            fx_ad(i, j) = fx_ad(i, j) + temp_ad47
            fx_ad(i+1, j) = fx_ad(i+1, j) - temp_ad47
            fy_ad(i, j) = fy_ad(i, j) + temp_ad47
            fy_ad(i, j+1) = fy_ad(i, j+1) - temp_ad47
            delp_ad(i, j) = wk_ad(i, j)
            wk_ad(i, j) = 0.0_8
          END DO
        END DO
      ELSE
        px_ad = 0.0_8
        py_ad = 0.0_8
        DO j=je,js,-1
          DO i=ie,is,-1
            temp_ad50 = pt_ad(i, j)/delp(i, j)
            delp_ad(i, j) = delp_ad(i, j) - pt(i, j)*temp_ad50/delp(i, j&
&             )
            pt_ad(i, j) = temp_ad50
            CALL POPREAL8(delp(i, j))
            temp_ad51 = rarea(i, j)*delp_ad(i, j)
            fx_ad(i, j) = fx_ad(i, j) + temp_ad51
            fx_ad(i+1, j) = fx_ad(i+1, j) - temp_ad51
            fy_ad(i, j) = fy_ad(i, j) + temp_ad51
            fy_ad(i, j+1) = fy_ad(i, j+1) - temp_ad51
            CALL POPREAL8(pt(i, j))
            temp_ad52 = rarea(i, j)*pt_ad(i, j)
            delp_ad(i, j) = delp_ad(i, j) + pt(i, j)*pt_ad(i, j)
            px_ad(i, j) = px_ad(i, j) + temp_ad52
            px_ad(i+1, j) = px_ad(i+1, j) - temp_ad52
            py_ad(i, j) = py_ad(i, j) + temp_ad52
            py_ad(i, j+1) = py_ad(i, j+1) - temp_ad52
            pt_ad(i, j) = delp(i, j)*pt_ad(i, j)
          END DO
        END DO
      END IF
      CALL POPREAL8ARRAY(pt, (ied-isd+1)*(jed-jsd+1))
      CALL POPREAL8ARRAY(px, (ie-is+2)*(je-js+1))
      CALL POPREAL8ARRAY(py, (ie-is+1)*(je-js+2))
      CALL FV_TP_2D_ADM(pt, pt_ad, crx_adv, crx_adv_ad, cry_adv, &
&                 cry_adv_ad, npx, npy, hord_tm_pert, px, px_ad, py, py_ad, &
&                 xfx_adv, xfx_adv_ad, yfx_adv, yfx_adv_ad, area, ra_x, &
&                 ra_x_ad, ra_y, ra_y_ad, mfx=fx, mfx_ad=fx_ad, mfy=fy, &
&                 mfy_ad=fy_ad)
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO j=jed,jsd,-1
          DO i=ied,isd,-1
            CALL POPREAL8(pt(i, j))
            temp = zvir*q(i, j, k, sphum) + 1.
            q_ad(i, j, k, sphum) = q_ad(i, j, k, sphum) - pt(i, j)*zvir*&
&             pt_ad(i, j)/temp**2
            pt_ad(i, j) = pt_ad(i, j)/temp
          END DO
        END DO
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO j=je,js,-1
          DO i=ie,is,-1
            CALL POPREAL8(w(i, j))
            temp_ad44 = rarea(i, j)*w_ad(i, j)
            delp_ad(i, j) = delp_ad(i, j) + w(i, j)*w_ad(i, j)
            px_ad(i, j) = px_ad(i, j) + temp_ad44
            px_ad(i+1, j) = px_ad(i+1, j) - temp_ad44
            py_ad(i, j) = py_ad(i, j) + temp_ad44
            py_ad(i, j+1) = py_ad(i, j+1) - temp_ad44
            w_ad(i, j) = delp(i, j)*w_ad(i, j)
          END DO
        END DO
        CALL POPREAL8ARRAY(w, (ied-isd+1)*(jed-jsd+1))
        CALL POPREAL8ARRAY(px, (ie-is+2)*(je-js+1))
        CALL POPREAL8ARRAY(py, (ie-is+1)*(je-js+2))
        CALL FV_TP_2D_ADM(w, w_ad, crx_adv, crx_adv_ad, cry_adv, &
&                   cry_adv_ad, npx, npy, hord_vt_pert, px, px_ad, py, py_ad&
&                   , xfx_adv, xfx_adv_ad, yfx_adv, yfx_adv_ad, area, &
&                   ra_x, ra_x_ad, ra_y, ra_y_ad, mfx=fx, mfx_ad=fx_ad, &
&                   mfy=fy, mfy_ad=fy_ad)
      END IF
      DO j=je+1,js,-1
        DO i=ie,is,-1
          fy_ad(i, j) = fy_ad(i, j) + yflux_ad(i, j)
        END DO
        DO i=ied,isd,-1
          cry_adv_ad(i, j) = cry_adv_ad(i, j) + cy_ad(i, j)
        END DO
      END DO
      DO j=je,js,-1
        DO i=ie+1,is,-1
          fx_ad(i, j) = fx_ad(i, j) + xflux_ad(i, j)
        END DO
      END DO
      DO j=jed,jsd,-1
        DO i=ie+1,is,-1
          crx_adv_ad(i, j) = crx_adv_ad(i, j) + cx_ad(i, j)
        END DO
      END DO
    END IF
    CALL POPREAL8ARRAY(delp, (ied-isd+1)*(jed-jsd+1))
    CALL POPREAL8ARRAY(fx, (ie-is+2)*(je-js+1))
    CALL POPREAL8ARRAY(fy, (ie-is+1)*(je-js+2))
    CALL FV_TP_2D_ADM(delp, delp_ad, crx_adv, crx_adv_ad, cry_adv, &
&               cry_adv_ad, npx, npy, hord_dp_pert, fx, fx_ad, fy, fy_ad, &
&               xfx_adv, xfx_adv_ad, yfx_adv, yfx_adv_ad, area, ra_x, &
&               ra_x_ad, ra_y, ra_y_ad)
    DO j=je,js,-1
      DO i=ied,isd,-1
        yfx_adv_ad(i, j) = yfx_adv_ad(i, j) + ra_y_ad(i, j)
        yfx_adv_ad(i, j+1) = yfx_adv_ad(i, j+1) - ra_y_ad(i, j)
        ra_y_ad(i, j) = 0.0_8
      END DO
    END DO
    DO j=jed,jsd,-1
      DO i=ie,is,-1
        xfx_adv_ad(i, j) = xfx_adv_ad(i, j) + ra_x_ad(i, j)
        xfx_adv_ad(i+1, j) = xfx_adv_ad(i+1, j) - ra_x_ad(i, j)
        ra_x_ad(i, j) = 0.0_8
      END DO
    END DO
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      DO j=je+1,js,-1
        DO i=ied,isd,-1
          yfx_adv_ad(i, j) = dx(i, j)*sina_v(i, j)*yfx_adv_ad(i, j)
        END DO
      END DO
      DO j=je+1,js,-1
        DO i=ied,isd,-1
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            yfx_adv_ad(i, j) = yfx_adv_ad(i, j) + rdya(i, j)*cry_adv_ad(&
&             i, j)
            cry_adv_ad(i, j) = 0.0_8
          ELSE
            yfx_adv_ad(i, j) = yfx_adv_ad(i, j) + rdya(i, j-1)*&
&             cry_adv_ad(i, j)
            cry_adv_ad(i, j) = 0.0_8
          END IF
        END DO
      END DO
      DO j=jed,jsd,-1
        DO i=ie+1,is,-1
          xfx_adv_ad(i, j) = dy(i, j)*sina_u(i, j)*xfx_adv_ad(i, j)
        END DO
      END DO
      DO j=jed,jsd,-1
        DO i=ie+1,is,-1
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            xfx_adv_ad(i, j) = xfx_adv_ad(i, j) + rdxa(i, j)*crx_adv_ad(&
&             i, j)
            crx_adv_ad(i, j) = 0.0_8
          ELSE
            xfx_adv_ad(i, j) = xfx_adv_ad(i, j) + rdxa(i-1, j)*&
&             crx_adv_ad(i, j)
            crx_adv_ad(i, j) = 0.0_8
          END IF
        END DO
      END DO
      DO j=je+1,js,-1
        DO i=ied,isd,-1
          vt_ad(i, j) = vt_ad(i, j) + dt*yfx_adv_ad(i, j)
          yfx_adv_ad(i, j) = 0.0_8
        END DO
      END DO
      DO j=jed,jsd,-1
        DO i=ie+1,is,-1
          ut_ad(i, j) = ut_ad(i, j) + dt*xfx_adv_ad(i, j)
          xfx_adv_ad(i, j) = 0.0_8
        END DO
      END DO
      CALL POPCONTROL2B(branch)
      IF (branch .EQ. 0) THEN
        DO j=je+2,js-1,-1
          DO i=ied,isd,-1
            vc_ad(i, j) = vc_ad(i, j) + vt_ad(i, j)
            vt_ad(i, j) = 0.0_8
          END DO
        END DO
        DO j=jed,jsd,-1
          DO i=ie+2,is-1,-1
            uc_ad(i, j) = uc_ad(i, j) + ut_ad(i, j)
            ut_ad(i, j) = 0.0_8
          END DO
        END DO
      ELSE
        IF (branch .EQ. 1) THEN
          uci_ad = 0.0_8
          vci_ad = 0.0_8
        ELSE
          damp = 1.0/(1.0-r1b16*cosa_u(2, npy-1)*cosa_v(1, npy-1))
          uci_ad = 0.0_8
          vci_ad = 0.0_8
          temp_ad35 = damp*r1b4*cosa_v(1, npy-1)*vt_ad(0, npy-1)
          temp_ad36 = r1b4*cosa_u(2, npy-1)*temp_ad35
          vci_ad(0, npy-1) = vci_ad(0, npy-1) + damp*vt_ad(0, npy-1)
          ut_ad(1, npy-1) = ut_ad(1, npy-1) + temp_ad35
          ut_ad(1, npy-2) = ut_ad(1, npy-2) + temp_ad35
          ut_ad(0, npy-2) = ut_ad(0, npy-2) + temp_ad35
          uci_ad(0, npy-1) = uci_ad(0, npy-1) + temp_ad35
          vt_ad(0, npy) = vt_ad(0, npy) + temp_ad36
          vt_ad(-1, npy) = vt_ad(-1, npy) + temp_ad36
          vt_ad(-1, npy-1) = vt_ad(-1, npy-1) + temp_ad36
          vt_ad(0, npy-1) = 0.0_8
          temp_ad37 = -(damp*r1b4*cosa_v(1, npy-1)*vt_ad(1, npy-1))
          temp_ad38 = -(r1b4*cosa_u(2, npy-1)*temp_ad37)
          ut_ad(1, npy-1) = ut_ad(1, npy-1) + temp_ad37
          ut_ad(1, npy-2) = ut_ad(1, npy-2) + temp_ad37
          ut_ad(2, npy-2) = ut_ad(2, npy-2) + temp_ad37
          uci_ad(2, npy-1) = uci_ad(2, npy-1) + damp*ut_ad(2, npy-1) + &
&           temp_ad37
          temp_ad39 = -(damp*r1b4*cosa_u(2, npy-1)*ut_ad(2, npy-1))
          vci_ad(1, npy-1) = vci_ad(1, npy-1) + temp_ad39 + damp*vt_ad(1&
&           , npy-1)
          vt_ad(1, npy) = vt_ad(1, npy) + temp_ad38
          vt_ad(2, npy) = vt_ad(2, npy) + temp_ad38
          vt_ad(2, npy-1) = vt_ad(2, npy-1) + temp_ad38
          vt_ad(1, npy-1) = 0.0_8
          temp_ad40 = -(r1b4*cosa_v(1, npy-1)*temp_ad39)
          vt_ad(1, npy) = vt_ad(1, npy) + temp_ad39
          vt_ad(2, npy) = vt_ad(2, npy) + temp_ad39
          vt_ad(2, npy-1) = vt_ad(2, npy-1) + temp_ad39
          ut_ad(1, npy-1) = ut_ad(1, npy-1) + temp_ad40
          ut_ad(1, npy-2) = ut_ad(1, npy-2) + temp_ad40
          ut_ad(2, npy-2) = ut_ad(2, npy-2) + temp_ad40
          ut_ad(2, npy-1) = 0.0_8
          temp_ad41 = damp*r1b4*cosa_u(2, npy-1)*ut_ad(2, npy)
          temp_ad42 = r1b4*cosa_v(1, npy-1)*temp_ad41
          uci_ad(2, npy) = uci_ad(2, npy) + damp*ut_ad(2, npy)
          vt_ad(1, npy) = vt_ad(1, npy) + temp_ad41
          vt_ad(2, npy) = vt_ad(2, npy) + temp_ad41
          vt_ad(2, npy+1) = vt_ad(2, npy+1) + temp_ad41
          vci_ad(1, npy+1) = vci_ad(1, npy+1) + temp_ad41
          ut_ad(1, npy) = ut_ad(1, npy) + temp_ad42
          ut_ad(1, npy+1) = ut_ad(1, npy+1) + temp_ad42
          ut_ad(2, npy+1) = ut_ad(2, npy+1) + temp_ad42
          ut_ad(2, npy) = 0.0_8
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          damp = 1.0/(1.0-r1b16*cosa_u(npx-1, npy-1)*cosa_v(npx-1, npy-1&
&           ))
          temp_ad27 = damp*r1b4*cosa_v(npx-1, npy-1)*vt_ad(npx, npy-1)
          temp_ad28 = r1b4*cosa_u(npx-1, npy-1)*temp_ad27
          vci_ad(npx, npy-1) = vci_ad(npx, npy-1) + damp*vt_ad(npx, npy-&
&           1)
          ut_ad(npx, npy-1) = ut_ad(npx, npy-1) + temp_ad27
          ut_ad(npx, npy-2) = ut_ad(npx, npy-2) + temp_ad27
          ut_ad(npx+1, npy-2) = ut_ad(npx+1, npy-2) + temp_ad27
          uci_ad(npx+1, npy-1) = uci_ad(npx+1, npy-1) + temp_ad27
          vt_ad(npx, npy) = vt_ad(npx, npy) + temp_ad28
          vt_ad(npx+1, npy) = vt_ad(npx+1, npy) + temp_ad28
          vt_ad(npx+1, npy-1) = vt_ad(npx+1, npy-1) + temp_ad28
          vt_ad(npx, npy-1) = 0.0_8
          temp_ad29 = -(damp*r1b4*cosa_v(npx-1, npy-1)*vt_ad(npx-1, npy-&
&           1))
          temp_ad30 = -(r1b4*cosa_u(npx-1, npy-1)*temp_ad29)
          ut_ad(npx, npy-1) = ut_ad(npx, npy-1) + temp_ad29
          ut_ad(npx, npy-2) = ut_ad(npx, npy-2) + temp_ad29
          ut_ad(npx-1, npy-2) = ut_ad(npx-1, npy-2) + temp_ad29
          uci_ad(npx-1, npy-1) = uci_ad(npx-1, npy-1) + damp*ut_ad(npx-1&
&           , npy-1) + temp_ad29
          temp_ad31 = -(damp*r1b4*cosa_u(npx-1, npy-1)*ut_ad(npx-1, npy-&
&           1))
          vci_ad(npx-1, npy-1) = vci_ad(npx-1, npy-1) + temp_ad31 + damp&
&           *vt_ad(npx-1, npy-1)
          vt_ad(npx-1, npy) = vt_ad(npx-1, npy) + temp_ad30
          vt_ad(npx-2, npy) = vt_ad(npx-2, npy) + temp_ad30
          vt_ad(npx-2, npy-1) = vt_ad(npx-2, npy-1) + temp_ad30
          vt_ad(npx-1, npy-1) = 0.0_8
          temp_ad32 = -(r1b4*cosa_v(npx-1, npy-1)*temp_ad31)
          vt_ad(npx-1, npy) = vt_ad(npx-1, npy) + temp_ad31
          vt_ad(npx-2, npy) = vt_ad(npx-2, npy) + temp_ad31
          vt_ad(npx-2, npy-1) = vt_ad(npx-2, npy-1) + temp_ad31
          ut_ad(npx, npy-1) = ut_ad(npx, npy-1) + temp_ad32
          ut_ad(npx, npy-2) = ut_ad(npx, npy-2) + temp_ad32
          ut_ad(npx-1, npy-2) = ut_ad(npx-1, npy-2) + temp_ad32
          ut_ad(npx-1, npy-1) = 0.0_8
          temp_ad33 = damp*r1b4*cosa_u(npx-1, npy-1)*ut_ad(npx-1, npy)
          temp_ad34 = r1b4*cosa_v(npx-1, npy-1)*temp_ad33
          uci_ad(npx-1, npy) = uci_ad(npx-1, npy) + damp*ut_ad(npx-1, &
&           npy)
          vt_ad(npx-1, npy) = vt_ad(npx-1, npy) + temp_ad33
          vt_ad(npx-2, npy) = vt_ad(npx-2, npy) + temp_ad33
          vt_ad(npx-2, npy+1) = vt_ad(npx-2, npy+1) + temp_ad33
          vci_ad(npx-1, npy+1) = vci_ad(npx-1, npy+1) + temp_ad33
          ut_ad(npx, npy) = ut_ad(npx, npy) + temp_ad34
          ut_ad(npx, npy+1) = ut_ad(npx, npy+1) + temp_ad34
          ut_ad(npx-1, npy+1) = ut_ad(npx-1, npy+1) + temp_ad34
          ut_ad(npx-1, npy) = 0.0_8
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          damp = 1.0/(1.0-r1b16*cosa_u(npx-1, 1)*cosa_v(npx-1, 2))
          temp_ad19 = damp*r1b4*cosa_v(npx-1, 2)*vt_ad(npx, 2)
          temp_ad20 = r1b4*cosa_u(npx-1, 1)*temp_ad19
          vci_ad(npx, 2) = vci_ad(npx, 2) + damp*vt_ad(npx, 2)
          ut_ad(npx, 1) = ut_ad(npx, 1) + temp_ad19
          ut_ad(npx, 2) = ut_ad(npx, 2) + temp_ad19
          ut_ad(npx+1, 2) = ut_ad(npx+1, 2) + temp_ad19
          uci_ad(npx+1, 1) = uci_ad(npx+1, 1) + temp_ad19
          vt_ad(npx, 1) = vt_ad(npx, 1) + temp_ad20
          vt_ad(npx+1, 1) = vt_ad(npx+1, 1) + temp_ad20
          vt_ad(npx+1, 2) = vt_ad(npx+1, 2) + temp_ad20
          vt_ad(npx, 2) = 0.0_8
          temp_ad21 = -(damp*r1b4*cosa_v(npx-1, 2)*vt_ad(npx-1, 2))
          temp_ad22 = -(r1b4*cosa_u(npx-1, 1)*temp_ad21)
          ut_ad(npx, 1) = ut_ad(npx, 1) + temp_ad21
          ut_ad(npx, 2) = ut_ad(npx, 2) + temp_ad21
          ut_ad(npx-1, 2) = ut_ad(npx-1, 2) + temp_ad21
          uci_ad(npx-1, 1) = uci_ad(npx-1, 1) + damp*ut_ad(npx-1, 1) + &
&           temp_ad21
          temp_ad23 = -(damp*r1b4*cosa_u(npx-1, 1)*ut_ad(npx-1, 1))
          vci_ad(npx-1, 2) = vci_ad(npx-1, 2) + temp_ad23 + damp*vt_ad(&
&           npx-1, 2)
          vt_ad(npx-1, 1) = vt_ad(npx-1, 1) + temp_ad22
          vt_ad(npx-2, 1) = vt_ad(npx-2, 1) + temp_ad22
          vt_ad(npx-2, 2) = vt_ad(npx-2, 2) + temp_ad22
          vt_ad(npx-1, 2) = 0.0_8
          temp_ad24 = -(r1b4*cosa_v(npx-1, 2)*temp_ad23)
          vt_ad(npx-1, 1) = vt_ad(npx-1, 1) + temp_ad23
          vt_ad(npx-2, 1) = vt_ad(npx-2, 1) + temp_ad23
          vt_ad(npx-2, 2) = vt_ad(npx-2, 2) + temp_ad23
          ut_ad(npx, 1) = ut_ad(npx, 1) + temp_ad24
          ut_ad(npx, 2) = ut_ad(npx, 2) + temp_ad24
          ut_ad(npx-1, 2) = ut_ad(npx-1, 2) + temp_ad24
          ut_ad(npx-1, 1) = 0.0_8
          temp_ad25 = damp*r1b4*cosa_u(npx-1, 1)*ut_ad(npx-1, 0)
          temp_ad26 = r1b4*cosa_v(npx-1, 2)*temp_ad25
          uci_ad(npx-1, 0) = uci_ad(npx-1, 0) + damp*ut_ad(npx-1, 0)
          vt_ad(npx-1, 1) = vt_ad(npx-1, 1) + temp_ad25
          vt_ad(npx-2, 1) = vt_ad(npx-2, 1) + temp_ad25
          vt_ad(npx-2, 0) = vt_ad(npx-2, 0) + temp_ad25
          vci_ad(npx-1, 0) = vci_ad(npx-1, 0) + temp_ad25
          ut_ad(npx, 0) = ut_ad(npx, 0) + temp_ad26
          ut_ad(npx, -1) = ut_ad(npx, -1) + temp_ad26
          ut_ad(npx-1, -1) = ut_ad(npx-1, -1) + temp_ad26
          ut_ad(npx-1, 0) = 0.0_8
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          damp = 1.0/(1.0-r1b16*cosa_u(2, 1)*cosa_v(1, 2))
          temp_ad11 = -(damp*r1b4*cosa_v(0, 2)*vt_ad(0, 2))
          temp_ad12 = -(r1b4*cosa_u(0, 1)*temp_ad11)
          vci_ad(0, 2) = vci_ad(0, 2) + damp*vt_ad(0, 2)
          ut_ad(1, 1) = ut_ad(1, 1) + temp_ad11
          ut_ad(1, 2) = ut_ad(1, 2) + temp_ad11
          ut_ad(0, 2) = ut_ad(0, 2) + temp_ad11
          uci_ad(0, 1) = uci_ad(0, 1) + temp_ad11
          vt_ad(0, 1) = vt_ad(0, 1) + temp_ad12
          vt_ad(-1, 1) = vt_ad(-1, 1) + temp_ad12
          vt_ad(-1, 2) = vt_ad(-1, 2) + temp_ad12
          vt_ad(0, 2) = 0.0_8
          temp_ad13 = -(damp*r1b4*cosa_v(1, 2)*vt_ad(1, 2))
          temp_ad14 = -(r1b4*cosa_u(2, 1)*temp_ad13)
          ut_ad(1, 1) = ut_ad(1, 1) + temp_ad13
          ut_ad(1, 2) = ut_ad(1, 2) + temp_ad13
          ut_ad(2, 2) = ut_ad(2, 2) + temp_ad13
          uci_ad(2, 1) = uci_ad(2, 1) + damp*ut_ad(2, 1) + temp_ad13
          temp_ad15 = -(damp*r1b4*cosa_u(2, 1)*ut_ad(2, 1))
          vci_ad(1, 2) = vci_ad(1, 2) + temp_ad15 + damp*vt_ad(1, 2)
          vt_ad(1, 1) = vt_ad(1, 1) + temp_ad14
          vt_ad(2, 1) = vt_ad(2, 1) + temp_ad14
          vt_ad(2, 2) = vt_ad(2, 2) + temp_ad14
          vt_ad(1, 2) = 0.0_8
          temp_ad16 = -(r1b4*cosa_v(1, 2)*temp_ad15)
          vt_ad(1, 1) = vt_ad(1, 1) + temp_ad15
          vt_ad(2, 1) = vt_ad(2, 1) + temp_ad15
          vt_ad(2, 2) = vt_ad(2, 2) + temp_ad15
          ut_ad(1, 1) = ut_ad(1, 1) + temp_ad16
          ut_ad(1, 2) = ut_ad(1, 2) + temp_ad16
          ut_ad(2, 2) = ut_ad(2, 2) + temp_ad16
          ut_ad(2, 1) = 0.0_8
          temp_ad17 = -(damp*r1b4*cosa_u(2, 0)*ut_ad(2, 0))
          temp_ad18 = -(r1b4*cosa_v(1, 0)*temp_ad17)
          uci_ad(2, 0) = uci_ad(2, 0) + damp*ut_ad(2, 0)
          vt_ad(1, 1) = vt_ad(1, 1) + temp_ad17
          vt_ad(2, 1) = vt_ad(2, 1) + temp_ad17
          vt_ad(2, 0) = vt_ad(2, 0) + temp_ad17
          vci_ad(1, 0) = vci_ad(1, 0) + temp_ad17
          ut_ad(1, 0) = ut_ad(1, 0) + temp_ad18
          ut_ad(1, -1) = ut_ad(1, -1) + temp_ad18
          ut_ad(2, -1) = ut_ad(2, -1) + temp_ad18
          ut_ad(2, 0) = 0.0_8
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          DO i=min4,max4,-1
            temp_ad9 = r1b4*cosa_u(i, npy-1)*ut_ad(i, npy)
            uci_ad(i, npy) = uci_ad(i, npy) + ut_ad(i, npy)
            vt_ad(i-1, npy) = vt_ad(i-1, npy) + temp_ad9
            vt_ad(i, npy) = vt_ad(i, npy) + temp_ad9
            vt_ad(i-1, npy+1) = vt_ad(i-1, npy+1) + temp_ad9
            vt_ad(i, npy+1) = vt_ad(i, npy+1) + temp_ad9
            ut_ad(i, npy) = 0.0_8
            temp_ad10 = -(r1b4*cosa_u(i, npy-1)*ut_ad(i, npy-1))
            uci_ad(i, npy-1) = uci_ad(i, npy-1) + ut_ad(i, npy-1)
            vt_ad(i-1, npy-1) = vt_ad(i-1, npy-1) + temp_ad10
            vt_ad(i, npy-1) = vt_ad(i, npy-1) + temp_ad10
            vt_ad(i-1, npy) = vt_ad(i-1, npy) + temp_ad10
            vt_ad(i, npy) = vt_ad(i, npy) + temp_ad10
            ut_ad(i, npy-1) = 0.0_8
          END DO
          DO i=ied,isd,-1
            vci_ad(i, npy) = vci_ad(i, npy) + rsin_v(i, npy)*vt_ad(i, &
&             npy)
            vt_ad(i, npy) = 0.0_8
          END DO
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          DO i=min3,max3,-1
            temp_ad7 = -(r1b4*cosa_u(i, 1)*ut_ad(i, 1))
            uci_ad(i, 1) = uci_ad(i, 1) + ut_ad(i, 1)
            vt_ad(i-1, 1) = vt_ad(i-1, 1) + temp_ad7
            vt_ad(i, 1) = vt_ad(i, 1) + temp_ad7
            vt_ad(i-1, 2) = vt_ad(i-1, 2) + temp_ad7
            vt_ad(i, 2) = vt_ad(i, 2) + temp_ad7
            ut_ad(i, 1) = 0.0_8
            temp_ad8 = r1b4*cosa_u(i, 1)*ut_ad(i, 0)
            uci_ad(i, 0) = uci_ad(i, 0) + ut_ad(i, 0)
            vt_ad(i-1, 0) = vt_ad(i-1, 0) + temp_ad8
            vt_ad(i, 0) = vt_ad(i, 0) + temp_ad8
            vt_ad(i-1, 1) = vt_ad(i-1, 1) + temp_ad8
            vt_ad(i, 1) = vt_ad(i, 1) + temp_ad8
            ut_ad(i, 0) = 0.0_8
          END DO
          DO i=ied,isd,-1
            vci_ad(i, 1) = vci_ad(i, 1) + rsin_v(i, 1)*vt_ad(i, 1)
            vt_ad(i, 1) = 0.0_8
          END DO
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          DO j=min2,max2,-1
            temp_ad5 = r1b4*cosa_v(npx-1, j)*vt_ad(npx, j)
            vci_ad(npx, j) = vci_ad(npx, j) + vt_ad(npx, j)
            ut_ad(npx, j-1) = ut_ad(npx, j-1) + temp_ad5
            ut_ad(npx+1, j-1) = ut_ad(npx+1, j-1) + temp_ad5
            ut_ad(npx, j) = ut_ad(npx, j) + temp_ad5
            ut_ad(npx+1, j) = ut_ad(npx+1, j) + temp_ad5
            vt_ad(npx, j) = 0.0_8
            temp_ad6 = -(r1b4*cosa_v(npx-1, j)*vt_ad(npx-1, j))
            vci_ad(npx-1, j) = vci_ad(npx-1, j) + vt_ad(npx-1, j)
            ut_ad(npx-1, j-1) = ut_ad(npx-1, j-1) + temp_ad6
            ut_ad(npx, j-1) = ut_ad(npx, j-1) + temp_ad6
            ut_ad(npx-1, j) = ut_ad(npx-1, j) + temp_ad6
            ut_ad(npx, j) = ut_ad(npx, j) + temp_ad6
            vt_ad(npx-1, j) = 0.0_8
          END DO
          DO j=jed,jsd,-1
            uci_ad(npx, j) = uci_ad(npx, j) + rsin_u(npx, j)*ut_ad(npx, &
&             j)
            ut_ad(npx, j) = 0.0_8
          END DO
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          DO j=min1,max1,-1
            temp_ad3 = -(r1b4*cosa_v(1, j)*vt_ad(1, j))
            vci_ad(1, j) = vci_ad(1, j) + vt_ad(1, j)
            ut_ad(1, j-1) = ut_ad(1, j-1) + temp_ad3
            ut_ad(2, j-1) = ut_ad(2, j-1) + temp_ad3
            ut_ad(1, j) = ut_ad(1, j) + temp_ad3
            ut_ad(2, j) = ut_ad(2, j) + temp_ad3
            vt_ad(1, j) = 0.0_8
            temp_ad4 = r1b4*cosa_v(1, j)*vt_ad(0, j)
            vci_ad(0, j) = vci_ad(0, j) + vt_ad(0, j)
            ut_ad(0, j-1) = ut_ad(0, j-1) + temp_ad4
            ut_ad(1, j-1) = ut_ad(1, j-1) + temp_ad4
            ut_ad(0, j) = ut_ad(0, j) + temp_ad4
            ut_ad(1, j) = ut_ad(1, j) + temp_ad4
            vt_ad(0, j) = 0.0_8
          END DO
          DO j=jed,jsd,-1
            uci_ad(1, j) = uci_ad(1, j) + rsin_u(1, j)*ut_ad(1, j)
            ut_ad(1, j) = 0.0_8
          END DO
        END IF
        DO j=je+2,js-1,-1
          CALL POPCONTROL1B(branch)
          IF (branch .NE. 0) THEN
            DO i=ied,isd,-1
              temp_ad1 = rsin_v(i, j)*vt_ad(i, j)
              temp_ad2 = -(r1b4*cosa_v(i, j)*temp_ad1)
              vci_ad(i, j) = vci_ad(i, j) + temp_ad1
              uci_ad(i, j-1) = uci_ad(i, j-1) + temp_ad2
              uci_ad(i+1, j-1) = uci_ad(i+1, j-1) + temp_ad2
              uci_ad(i, j) = uci_ad(i, j) + temp_ad2
              uci_ad(i+1, j) = uci_ad(i+1, j) + temp_ad2
              vt_ad(i, j) = 0.0_8
            END DO
          END IF
        END DO
        DO j=jed,jsd,-1
          CALL POPCONTROL1B(branch)
          IF (branch .NE. 0) THEN
            DO i=ie+2,is-1,-1
              temp_ad = rsin_u(i, j)*ut_ad(i, j)
              temp_ad0 = -(r1b4*cosa_u(i, j)*temp_ad)
              uci_ad(i, j) = uci_ad(i, j) + temp_ad
              vci_ad(i-1, j) = vci_ad(i-1, j) + temp_ad0
              vci_ad(i, j) = vci_ad(i, j) + temp_ad0
              vci_ad(i-1, j+1) = vci_ad(i-1, j+1) + temp_ad0
              vci_ad(i, j+1) = vci_ad(i, j+1) + temp_ad0
              ut_ad(i, j) = 0.0_8
            END DO
          END IF
        END DO
        vc_ad = vc_ad + vci_ad
        uc_ad = uc_ad + uci_ad
      END IF
    ELSE
      DO j=je+1,js,-1
        DO i=ied,isd,-1
          yfx_adv_ad(i, j) = dx(i, j)*sina_v(i, j)*yfx_adv_ad(i, j)
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            yfx_adv_ad(i, j) = yfx_adv_ad(i, j) + rdya(i, j-1)*&
&             cry_adv_ad(i, j)
            cry_adv_ad(i, j) = 0.0_8
          ELSE
            yfx_adv_ad(i, j) = yfx_adv_ad(i, j) + rdya(i, j)*cry_adv_ad(&
&             i, j)
            cry_adv_ad(i, j) = 0.0_8
          END IF
          vc_ad(i, j) = vc_ad(i, j) + dt*yfx_adv_ad(i, j)/sina_v(i, j)
          yfx_adv_ad(i, j) = 0.0_8
        END DO
      END DO
      DO j=jed,jsd,-1
        DO i=ie+1,is,-1
          xfx_adv_ad(i, j) = dy(i, j)*sina_u(i, j)*xfx_adv_ad(i, j)
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            xfx_adv_ad(i, j) = xfx_adv_ad(i, j) + rdxa(i-1, j)*&
&             crx_adv_ad(i, j)
            crx_adv_ad(i, j) = 0.0_8
          ELSE
            xfx_adv_ad(i, j) = xfx_adv_ad(i, j) + rdxa(i, j)*crx_adv_ad(&
&             i, j)
            crx_adv_ad(i, j) = 0.0_8
          END IF
          uc_ad(i, j) = uc_ad(i, j) + dt*xfx_adv_ad(i, j)/sina_u(i, j)
          xfx_adv_ad(i, j) = 0.0_8
        END DO
      END DO
    END IF
  END SUBROUTINE D_SW_ADM
!-------------------------------------------------------------------------------
!
!     d_sw :: D-Grid Shallow Water Routine
!
  SUBROUTINE D_SW(delpc, delp, ptc, pt, u, v, w, uc, vc, ua, va, divg_d&
&   , xflux, yflux, cx, cy, crx_adv, cry_adv, xfx_adv, yfx_adv, zvir, &
&   sphum, nq, q, k, km, inline_q, pkz, dt, hord_tr, hord_mt, hord_vt, &
&   hord_tm, hord_dp, nord, dddmp, d2_bg, d4_bg, vtdm4, d_con, &
&   hydrostatic, ppm_limiter)
    IMPLICIT NONE
! damping applied to relative vorticity:
!   if ( vtdm4>0. ) then
!      damp4 = (vtdm4*da_min)**(nord+1)
!      call del6_flux(nord, npx, npy, damp4, wk, u, v)
!   endif
    INTEGER, INTENT(IN) :: hord_tr, hord_mt, hord_vt, hord_tm, hord_dp
! nord=1 (del-4) or 3 (del-8)
    INTEGER, INTENT(IN) :: nord
    INTEGER, INTENT(IN) :: sphum, nq, k, km
    REAL, INTENT(IN) :: dt, dddmp, d2_bg, d4_bg, vtdm4, d_con
    REAL, INTENT(IN) :: zvir
    REAL, INTENT(IN) :: ppm_limiter
    REAL(p_precision), INTENT(IN) :: pkz(is:ie, js:je)
! divergence
    REAL, INTENT(INOUT) :: divg_d(isd:ied+1, jsd:jed+1)
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(INOUT) :: delp, pt, ua, va&
&   , w
    REAL, DIMENSION(isd:ied, jsd:jed + 1), INTENT(INOUT) :: u, vc
    REAL, DIMENSION(isd:ied + 1, jsd:jed), INTENT(INOUT) :: v, uc
    REAL, INTENT(INOUT) :: q(isd:ied, jsd:jed, km, nq)
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(OUT) :: delpc, ptc
! The flux capacitors:
    REAL, INTENT(INOUT) :: xflux(is:ie+1, js:je)
    REAL, INTENT(INOUT) :: yflux(is:ie, js:je+1)
!------------------------
    REAL, INTENT(INOUT) :: cx(is:ie+1, jsd:jed)
    REAL, INTENT(INOUT) :: cy(isd:ied, js:je+1)
    LOGICAL, INTENT(IN) :: hydrostatic
    LOGICAL, INTENT(IN) :: inline_q
    REAL, DIMENSION(is:ie + 1, jsd:jed), INTENT(OUT) :: crx_adv, xfx_adv
    REAL, DIMENSION(isd:ied, js:je + 1), INTENT(OUT) :: cry_adv, yfx_adv
! Local:
    REAL(i_precision), DIMENSION(isd:ied + 1, jsd:jed) :: ut, uci
    REAL(i_precision), DIMENSION(isd:ied, jsd:jed + 1) :: vt, vci
    REAL(i_precision) :: damp
    REAL(i_precision), SAVE :: r1b4=1.0/4.0
    REAL(i_precision), SAVE :: r1b16=1.0/16.0
    REAL(ke_precision), DIMENSION(is:ie + 1, js:je + 1) :: ub, vb
!  needs this for corner_comm
    REAL(ke_precision) :: ke(isd:ied+1, jsd:jed+1)
    REAL(ke_precision) :: dt4, dt5, dt6
!#ifdef FULL_SMAG
!      real :: vt2(is-1:ie+1,js-1:je+1)     !  work array
!#endif
!  work array
    REAL :: wk(isd:ied, jsd:jed)
! Vorticity
    REAL :: vort(isd:ied, jsd:jed)
! 1-D X-direction Fluxes
    REAL :: fx(is:ie+1, js:je)
! 1-D Y-direction Fluxes
    REAL :: fy(is:ie, js:je+1)
! 1-D X-direction Fluxes
    REAL :: px(is:ie+1, js:je)
! 1-D Y-direction Fluxes
    REAL :: py(is:ie, js:je+1)
    REAL :: ra_x(is:ie, jsd:jed)
    REAL :: ra_y(isd:ied, js:je)
! work x-dir flux
    REAL :: gx(is:ie+1, js:je)
! work Y-dir flux array
    REAL :: gy(is:ie, js:je+1)
    LOGICAL :: fill_c
!#ifdef DOUBLE_GRADIENTS
!      real*8 :: pt_dp, delp_dp0, delp_dp
!      real*8 :: pt_fx_r8(is:ie+1,js:je+1)
!      real*8 :: pt_fy_r8(is:ie  ,js:je+1)
!      real*8 :: dp_fx_r8(is:ie+1,js:je  )
!      real*8 :: dp_fy_r8(is:ie  ,js:je+1)
!#endif
    REAL :: damp2, damp4, dd8, u2, v2, du2, dv2
    INTEGER :: i, j, is2, ie1, js2, je1, n, nt, n2, iq
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC ABS
    INTEGER :: arg1
    INTEGER :: arg2
    INTEGER :: min4
    INTEGER :: min3
    INTEGER :: min2
    INTEGER :: min1
    REAL :: abs0
    REAL :: max6
    REAL :: max5
    INTEGER :: max4
    INTEGER :: max3
    INTEGER :: max2
    INTEGER :: max1
    REAL :: y3
    REAL :: y2
    REAL :: y1
! shallow water test case 1
    IF (test_case .EQ. 1) THEN
      DO j=jsd,jed
        DO i=is,ie+1
          xfx_adv(i, j) = dt*uc(i, j)/sina_u(i, j)
          IF (xfx_adv(i, j) .GT. 0.) THEN
            crx_adv(i, j) = xfx_adv(i, j)*rdxa(i-1, j)
          ELSE
            crx_adv(i, j) = xfx_adv(i, j)*rdxa(i, j)
          END IF
          xfx_adv(i, j) = dy(i, j)*xfx_adv(i, j)*sina_u(i, j)
        END DO
      END DO
      DO j=js,je+1
        DO i=isd,ied
          yfx_adv(i, j) = dt*vc(i, j)/sina_v(i, j)
          IF (yfx_adv(i, j) .GT. 0.) THEN
            cry_adv(i, j) = yfx_adv(i, j)*rdya(i, j-1)
          ELSE
            cry_adv(i, j) = yfx_adv(i, j)*rdya(i, j)
          END IF
          yfx_adv(i, j) = dx(i, j)*yfx_adv(i, j)*sina_v(i, j)
        END DO
      END DO
    ELSE
! end grid_type choices
      IF (grid_type .LT. 3) THEN
        uci = uc
        vci = vc
! Interior:
        DO j=jsd,jed
          IF (j .NE. 0 .AND. j .NE. 1 .AND. j .NE. npy - 1 .AND. j .NE. &
&             npy) THEN
            DO i=is-1,ie+2
              ut(i, j) = (uci(i, j)-r1b4*cosa_u(i, j)*(vci(i-1, j)+vci(i&
&               , j)+vci(i-1, j+1)+vci(i, j+1)))*rsin_u(i, j)
            END DO
          END IF
        END DO
        DO j=js-1,je+2
          IF (j .NE. 1 .AND. j .NE. npy) THEN
            DO i=isd,ied
              vt(i, j) = (vci(i, j)-r1b4*cosa_v(i, j)*(uci(i, j-1)+uci(i&
&               +1, j-1)+uci(i, j)+uci(i+1, j)))*rsin_v(i, j)
            END DO
          END IF
        END DO
! West edge:
        IF (is .EQ. 1) THEN
          DO j=jsd,jed
            ut(1, j) = uci(1, j)*rsin_u(1, j)
          END DO
          IF (3 .LT. js) THEN
            max1 = js
          ELSE
            max1 = 3
          END IF
          IF (npy - 2 .GT. je + 1) THEN
            min1 = je + 1
          ELSE
            min1 = npy - 2
          END IF
          DO j=max1,min1
!             vt(0,j) = vci(0,j) - r1b4*cosa_v(0,j)*   &
            vt(0, j) = vci(0, j) + r1b4*cosa_v(1, j)*(ut(0, j-1)+ut(1, j&
&             -1)+ut(0, j)+ut(1, j))
            vt(1, j) = vci(1, j) - r1b4*cosa_v(1, j)*(ut(1, j-1)+ut(2, j&
&             -1)+ut(1, j)+ut(2, j))
          END DO
        END IF
! East edge:
        IF (ie + 1 .EQ. npx) THEN
          DO j=jsd,jed
            ut(npx, j) = uci(npx, j)*rsin_u(npx, j)
          END DO
          IF (3 .LT. js) THEN
            max2 = js
          ELSE
            max2 = 3
          END IF
          IF (npy - 2 .GT. je + 1) THEN
            min2 = je + 1
          ELSE
            min2 = npy - 2
          END IF
          DO j=max2,min2
            vt(npx-1, j) = vci(npx-1, j) - r1b4*cosa_v(npx-1, j)*(ut(npx&
&             -1, j-1)+ut(npx, j-1)+ut(npx-1, j)+ut(npx, j))
!             vt(npx,j) = vci(npx,j) - r1b4*cosa_v(npx,j)*   &
            vt(npx, j) = vci(npx, j) + r1b4*cosa_v(npx-1, j)*(ut(npx, j-&
&             1)+ut(npx+1, j-1)+ut(npx, j)+ut(npx+1, j))
          END DO
        END IF
! South (Bottom) edge:
        IF (js .EQ. 1) THEN
          DO i=isd,ied
            vt(i, 1) = vci(i, 1)*rsin_v(i, 1)
          END DO
          IF (3 .LT. is) THEN
            max3 = is
          ELSE
            max3 = 3
          END IF
          IF (npx - 2 .GT. ie + 1) THEN
            min3 = ie + 1
          ELSE
            min3 = npx - 2
          END IF
          DO i=max3,min3
!             ut(i,0) = uci(i,0) - r1b4*cosa_u(i,0)*   &
            ut(i, 0) = uci(i, 0) + r1b4*cosa_u(i, 1)*(vt(i-1, 0)+vt(i, 0&
&             )+vt(i-1, 1)+vt(i, 1))
            ut(i, 1) = uci(i, 1) - r1b4*cosa_u(i, 1)*(vt(i-1, 1)+vt(i, 1&
&             )+vt(i-1, 2)+vt(i, 2))
          END DO
        END IF
! North edge:
        IF (je + 1 .EQ. npy) THEN
          DO i=isd,ied
            vt(i, npy) = vci(i, npy)*rsin_v(i, npy)
          END DO
          IF (3 .LT. is) THEN
            max4 = is
          ELSE
            max4 = 3
          END IF
          IF (npx - 2 .GT. ie + 1) THEN
            min4 = ie + 1
          ELSE
            min4 = npx - 2
          END IF
          DO i=max4,min4
            ut(i, npy-1) = uci(i, npy-1) - r1b4*cosa_u(i, npy-1)*(vt(i-1&
&             , npy-1)+vt(i, npy-1)+vt(i-1, npy)+vt(i, npy))
!             ut(i,npy) = uci(i,npy) - r1b4*cosa_u(i,npy)*   &
            ut(i, npy) = uci(i, npy) + r1b4*cosa_u(i, npy-1)*(vt(i-1, &
&             npy)+vt(i, npy)+vt(i-1, npy+1)+vt(i, npy+1))
          END DO
        END IF
        IF (sw_corner) THEN
          damp = 1.0/(1.0-r1b16*cosa_u(2, 1)*cosa_v(1, 2))
          ut(2, 0) = (uci(2, 0)-r1b4*cosa_u(2, 0)*(vt(1, 1)+vt(2, 1)+vt(&
&           2, 0)+vci(1, 0)-r1b4*cosa_v(1, 0)*(ut(1, 0)+ut(1, -1)+ut(2, &
&           -1))))*damp
          ut(2, 1) = (uci(2, 1)-r1b4*cosa_u(2, 1)*(vt(1, 1)+vt(2, 1)+vt(&
&           2, 2)+vci(1, 2)-r1b4*cosa_v(1, 2)*(ut(1, 1)+ut(1, 2)+ut(2, 2&
&           ))))*damp
          vt(1, 2) = (vci(1, 2)-r1b4*cosa_v(1, 2)*(ut(1, 1)+ut(1, 2)+ut(&
&           2, 2)+uci(2, 1)-r1b4*cosa_u(2, 1)*(vt(1, 1)+vt(2, 1)+vt(2, 2&
&           ))))*damp
          vt(0, 2) = (vci(0, 2)-r1b4*cosa_v(0, 2)*(ut(1, 1)+ut(1, 2)+ut(&
&           0, 2)+uci(0, 1)-r1b4*cosa_u(0, 1)*(vt(0, 1)+vt(-1, 1)+vt(-1&
&           , 2))))*damp
        END IF
        IF (se_corner) THEN
          damp = 1.0/(1.0-r1b16*cosa_u(npx-1, 1)*cosa_v(npx-1, 2))
          ut(npx-1, 0) = (uci(npx-1, 0)+r1b4*cosa_u(npx-1, 1)*(vt(npx-1&
&           , 1)+vt(npx-2, 1)+vt(npx-2, 0)+vci(npx-1, 0)+r1b4*cosa_v(npx&
&           -1, 2)*(ut(npx, 0)+ut(npx, -1)+ut(npx-1, -1))))*damp
          ut(npx-1, 1) = (uci(npx-1, 1)-r1b4*cosa_u(npx-1, 1)*(vt(npx-1&
&           , 1)+vt(npx-2, 1)+vt(npx-2, 2)+vci(npx-1, 2)-r1b4*cosa_v(npx&
&           -1, 2)*(ut(npx, 1)+ut(npx, 2)+ut(npx-1, 2))))*damp
          vt(npx-1, 2) = (vci(npx-1, 2)-r1b4*cosa_v(npx-1, 2)*(ut(npx, 1&
&           )+ut(npx, 2)+ut(npx-1, 2)+uci(npx-1, 1)-r1b4*cosa_u(npx-1, 1&
&           )*(vt(npx-1, 1)+vt(npx-2, 1)+vt(npx-2, 2))))*damp
          vt(npx, 2) = (vci(npx, 2)+r1b4*cosa_v(npx-1, 2)*(ut(npx, 1)+ut&
&           (npx, 2)+ut(npx+1, 2)+uci(npx+1, 1)+r1b4*cosa_u(npx-1, 1)*(&
&           vt(npx, 1)+vt(npx+1, 1)+vt(npx+1, 2))))*damp
        END IF
        IF (ne_corner) THEN
          damp = 1.0/(1.0-r1b16*cosa_u(npx-1, npy-1)*cosa_v(npx-1, npy-1&
&           ))
          ut(npx-1, npy) = (uci(npx-1, npy)+r1b4*cosa_u(npx-1, npy-1)*(&
&           vt(npx-1, npy)+vt(npx-2, npy)+vt(npx-2, npy+1)+vci(npx-1, &
&           npy+1)+r1b4*cosa_v(npx-1, npy-1)*(ut(npx, npy)+ut(npx, npy+1&
&           )+ut(npx-1, npy+1))))*damp
          ut(npx-1, npy-1) = (uci(npx-1, npy-1)-r1b4*cosa_u(npx-1, npy-1&
&           )*(vt(npx-1, npy)+vt(npx-2, npy)+vt(npx-2, npy-1)+vci(npx-1&
&           , npy-1)-r1b4*cosa_v(npx-1, npy-1)*(ut(npx, npy-1)+ut(npx, &
&           npy-2)+ut(npx-1, npy-2))))*damp
          vt(npx-1, npy-1) = (vci(npx-1, npy-1)-r1b4*cosa_v(npx-1, npy-1&
&           )*(ut(npx, npy-1)+ut(npx, npy-2)+ut(npx-1, npy-2)+uci(npx-1&
&           , npy-1)-r1b4*cosa_u(npx-1, npy-1)*(vt(npx-1, npy)+vt(npx-2&
&           , npy)+vt(npx-2, npy-1))))*damp
          vt(npx, npy-1) = (vci(npx, npy-1)+r1b4*cosa_v(npx-1, npy-1)*(&
&           ut(npx, npy-1)+ut(npx, npy-2)+ut(npx+1, npy-2)+uci(npx+1, &
&           npy-1)+r1b4*cosa_u(npx-1, npy-1)*(vt(npx, npy)+vt(npx+1, npy&
&           )+vt(npx+1, npy-1))))*damp
        END IF
        IF (nw_corner) THEN
          damp = 1.0/(1.0-r1b16*cosa_u(2, npy-1)*cosa_v(1, npy-1))
          ut(2, npy) = (uci(2, npy)+r1b4*cosa_u(2, npy-1)*(vt(1, npy)+vt&
&           (2, npy)+vt(2, npy+1)+vci(1, npy+1)+r1b4*cosa_v(1, npy-1)*(&
&           ut(1, npy)+ut(1, npy+1)+ut(2, npy+1))))*damp
          ut(2, npy-1) = (uci(2, npy-1)-r1b4*cosa_u(2, npy-1)*(vt(1, npy&
&           )+vt(2, npy)+vt(2, npy-1)+vci(1, npy-1)-r1b4*cosa_v(1, npy-1&
&           )*(ut(1, npy-1)+ut(1, npy-2)+ut(2, npy-2))))*damp
          vt(1, npy-1) = (vci(1, npy-1)-r1b4*cosa_v(1, npy-1)*(ut(1, npy&
&           -1)+ut(1, npy-2)+ut(2, npy-2)+uci(2, npy-1)-r1b4*cosa_u(2, &
&           npy-1)*(vt(1, npy)+vt(2, npy)+vt(2, npy-1))))*damp
          vt(0, npy-1) = (vci(0, npy-1)+r1b4*cosa_v(1, npy-1)*(ut(1, npy&
&           -1)+ut(1, npy-2)+ut(0, npy-2)+uci(0, npy-1)+r1b4*cosa_u(2, &
&           npy-1)*(vt(0, npy)+vt(-1, npy)+vt(-1, npy-1))))*damp
        END IF
      ELSE
! grid_type >= 3
        DO j=jsd,jed
          DO i=is-1,ie+2
            ut(i, j) = uc(i, j)
          END DO
        END DO
        DO j=js-1,je+2
          DO i=isd,ied
            vt(i, j) = vc(i, j)
          END DO
        END DO
      END IF
      DO j=jsd,jed
        DO i=is,ie+1
          xfx_adv(i, j) = dt*ut(i, j)
        END DO
      END DO
      DO j=js,je+1
        DO i=isd,ied
          yfx_adv(i, j) = dt*vt(i, j)
        END DO
      END DO
! Compute E-W CFL number:
      DO j=jsd,jed
        DO i=is,ie+1
          IF (xfx_adv(i, j) .GT. 0.) THEN
            crx_adv(i, j) = xfx_adv(i, j)*rdxa(i-1, j)
          ELSE
            crx_adv(i, j) = xfx_adv(i, j)*rdxa(i, j)
          END IF
        END DO
      END DO
      DO j=jsd,jed
        DO i=is,ie+1
          xfx_adv(i, j) = dy(i, j)*xfx_adv(i, j)*sina_u(i, j)
        END DO
      END DO
      DO j=js,je+1
        DO i=isd,ied
          IF (yfx_adv(i, j) .GT. 0.) THEN
            cry_adv(i, j) = yfx_adv(i, j)*rdya(i, j-1)
          ELSE
            cry_adv(i, j) = yfx_adv(i, j)*rdya(i, j)
          END IF
        END DO
      END DO
      DO j=js,je+1
        DO i=isd,ied
          yfx_adv(i, j) = dx(i, j)*yfx_adv(i, j)*sina_v(i, j)
        END DO
      END DO
    END IF
    DO j=jsd,jed
      DO i=is,ie
        ra_x(i, j) = area(i, j) + xfx_adv(i, j) - xfx_adv(i+1, j)
      END DO
    END DO
    DO j=js,je
      DO i=isd,ied
        ra_y(i, j) = area(i, j) + yfx_adv(i, j) - yfx_adv(i, j+1)
      END DO
    END DO
    CALL FV_TP_2D(delp, crx_adv, cry_adv, npx, npy, hord_dp, fx, fy, &
&           xfx_adv, yfx_adv, area, ra_x, ra_y)
! shallow_water
    IF (shallow_water) THEN
!#ifdef DOUBLE_GRADIENTS
!        dp_fx_r8 = fx
!        dp_fy_r8 = fy
!        do j=js,je
!           do i=is,ie
!              delp(i,j) = delp(i,j) + (dp_fx_r8(i,j)-dp_fx_r8(i+1,j)+dp_fy_r8(i,j)-dp_fy_r8(i,j+1))*rarea(i,j)
!              ptc(i,j) = pt(i,j)
!           enddo
!        enddo
!#else
      DO j=js,je
        DO i=is,ie
          delp(i, j) = delp(i, j) + (fx(i, j)-fx(i+1, j)+fy(i, j)-fy(i, &
&           j+1))*rarea(i, j)
          ptc(i, j) = pt(i, j)
        END DO
      END DO
    ELSE
!#endif
! <<< Save the mass fluxes to the "Flux Capacitor" for tracer transport >>>
      DO j=jsd,jed
        DO i=is,ie+1
          cx(i, j) = cx(i, j) + crx_adv(i, j)
        END DO
      END DO
      DO j=js,je
        DO i=is,ie+1
          xflux(i, j) = xflux(i, j) + fx(i, j)
        END DO
      END DO
      DO j=js,je+1
        DO i=isd,ied
          cy(i, j) = cy(i, j) + cry_adv(i, j)
        END DO
        DO i=is,ie
          yflux(i, j) = yflux(i, j) + fy(i, j)
        END DO
      END DO
      IF (.NOT.hydrostatic) THEN
        CALL FV_TP_2D(w, crx_adv, cry_adv, npx, npy, hord_vt, px, py, &
&               xfx_adv, yfx_adv, area, ra_x, ra_y, mfx=fx, mfy=fy)
        DO j=js,je
          DO i=is,ie
            w(i, j) = w(i, j)*delp(i, j) + (px(i, j)-px(i+1, j)+py(i, j)&
&             -py(i, j+1))*rarea(i, j)
          END DO
        END DO
      END IF
      IF (inline_q) THEN
        DO j=jsd,jed
          DO i=isd,ied
            pt(i, j) = pt(i, j)/(1.+zvir*q(i, j, k, sphum))
          END DO
        END DO
      END IF
      CALL FV_TP_2D(pt, crx_adv, cry_adv, npx, npy, hord_tm, px, py, &
&             xfx_adv, yfx_adv, area, ra_x, ra_y, mfx=fx, mfy=fy)
      IF (inline_q) THEN
        DO j=js,je
          DO i=is,ie
            wk(i, j) = delp(i, j)
            delp(i, j) = wk(i, j) + (fx(i, j)-fx(i+1, j)+fy(i, j)-fy(i, &
&             j+1))*rarea(i, j)
            pt(i, j) = (pt(i, j)*wk(i, j)+(px(i, j)-px(i+1, j)+py(i, j)-&
&             py(i, j+1))*rarea(i, j))/delp(i, j)
          END DO
        END DO
        DO iq=1,nq
          CALL FV_TP_2D(q(isd:, jsd:, k, iq), crx_adv, cry_adv, npx, npy&
&                 , hord_tr, px, py, xfx_adv, yfx_adv, area, ra_x, ra_y&
&                 , mfx=fx, mfy=fy)
          DO j=js,je
            DO i=is,ie
              q(i, j, k, iq) = (q(i, j, k, iq)*wk(i, j)+(px(i, j)-px(i+1&
&               , j)+py(i, j)-py(i, j+1))*rarea(i, j))/delp(i, j)
            END DO
          END DO
        END DO
        DO j=js,je
          DO i=is,ie
            pt(i, j) = pt(i, j)*(1.+zvir*q(i, j, k, sphum))
          END DO
        END DO
      ELSE
!#ifdef DOUBLE_GRADIENTS
!        pt_fx_r8 = px
!        pt_fy_r8 = py
!        dp_fx_r8 = fx
!        dp_fy_r8 = fy
!        do j=js,je
!           do i=is,ie
!              pt_dp = pt(i,j)
!              delp_dp = delp(i,j)
!              pt(i,j)   = pt_dp*delp_dp + (pt_fx_r8(i,j)-pt_fx_r8(i+1,j)+pt_fy_r8(i,j)-pt_fy_r8(i,j+1))*rarea(i,j)
!              delp(i,j) =       delp_dp + (dp_fx_r8(i,j)-dp_fx_r8(i+1,j)+dp_fy_r8(i,j)-dp_fy_r8(i,j+1))*rarea(i,j)
!              pt(i,j) = pt(i,j) / delp(i,j)
!           enddo
!        enddo
!#else
        DO j=js,je
          DO i=is,ie
            pt(i, j) = pt(i, j)*delp(i, j) + (px(i, j)-px(i+1, j)+py(i, &
&             j)-py(i, j+1))*rarea(i, j)
            delp(i, j) = delp(i, j) + (fx(i, j)-fx(i+1, j)+fy(i, j)-fy(i&
&             , j+1))*rarea(i, j)
            pt(i, j) = pt(i, j)/delp(i, j)
          END DO
        END DO
      END IF
!#endif
      IF (.NOT.hydrostatic) THEN
        DO j=js,je
          DO i=is,ie
            w(i, j) = w(i, j)/delp(i, j)
          END DO
        END DO
      END IF
    END IF
! shallow_water test_case 1
    IF (test_case .GT. 1) THEN
!----------------------
! Kinetic Energy Fluxes
!----------------------
! Compute B grid contra-variant components for KE:
      dt5 = 0.5*dt
      dt4 = 0.25*dt
      IF (2 .LT. is) THEN
        is2 = is
      ELSE
        is2 = 2
      END IF
      IF (npx - 1 .GT. ie + 1) THEN
        ie1 = ie + 1
      ELSE
        ie1 = npx - 1
      END IF
      IF (2 .LT. js) THEN
        js2 = js
      ELSE
        js2 = 2
      END IF
      IF (npy - 1 .GT. je + 1) THEN
        je1 = je + 1
      ELSE
        je1 = npy - 1
      END IF
      IF (grid_type .LT. 3) THEN
        IF (js .EQ. 1) THEN
          DO i=is,ie+1
! corner values are incorrect
            vb(i, 1) = dt5*(vt(i-1, 1)+vt(i, 1))
          END DO
        END IF
        DO j=js2,je1
          DO i=is2,ie1
            vb(i, j) = dt5*(vc(i-1, j)+vc(i, j)-(uc(i, j-1)+uc(i, j))*&
&             cosa(i, j))*rsina(i, j)
          END DO
          IF (is .EQ. 1) vb(1, j) = dt4*(-vt(-1, j)+3.*(vt(0, j)+vt(1, j&
&             ))-vt(2, j))
!              vb(1,j) = dt5*(vt(0,j)+vt(1,j)) 
! 2-pt extrapolation from both sides:
          IF (ie + 1 .EQ. npx) vb(npx, j) = dt4*(-vt(npx-2, j)+3.*(vt(&
&             npx-1, j)+vt(npx, j))-vt(npx+1, j))
!              vb(npx,j) = dt5*(vt(npx-1,j)+vt(npx,j))
! 2-pt extrapolation from both sides:
        END DO
        IF (je + 1 .EQ. npy) THEN
          DO i=is,ie+1
! corner values are incorrect
            vb(i, npy) = dt5*(vt(i-1, npy)+vt(i, npy))
          END DO
        END IF
      ELSE
        DO j=js,je+1
          DO i=is,ie+1
            vb(i, j) = dt5*(vc(i-1, j)+vc(i, j))
          END DO
        END DO
      END IF
      CALL YTP_V(vb, u, v, ub, hord_mt)
      DO j=js,je+1
        DO i=is,ie+1
          ke(i, j) = vb(i, j)*ub(i, j)
        END DO
      END DO
      IF (grid_type .LT. 3) THEN
        IF (is .EQ. 1) THEN
          DO j=js,je+1
! corner values are incorrect
            ub(1, j) = dt5*(ut(1, j-1)+ut(1, j))
          END DO
        END IF
        DO j=js,je+1
          IF (j .EQ. 1 .OR. j .EQ. npy) THEN
            DO i=is2,ie1
!                 ub(i,j) = dt5*(ut(i,j-1)+ut(i,j))
! 2-pt extrapolation from both sides:
              ub(i, j) = dt4*(-ut(i, j-2)+3.*(ut(i, j-1)+ut(i, j))-ut(i&
&               , j+1))
            END DO
          ELSE
            DO i=is2,ie1
              ub(i, j) = dt5*(uc(i, j-1)+uc(i, j)-(vc(i-1, j)+vc(i, j))*&
&               cosa(i, j))*rsina(i, j)
            END DO
          END IF
        END DO
        IF (ie + 1 .EQ. npx) THEN
          DO j=js,je+1
! corner values are incorrect
            ub(npx, j) = dt5*(ut(npx, j-1)+ut(npx, j))
          END DO
        END IF
      ELSE
        DO j=js,je+1
          DO i=is,ie+1
            ub(i, j) = dt5*(uc(i, j-1)+uc(i, j))
          END DO
        END DO
      END IF
      CALL XTP_U(ub, u, v, vb, hord_mt)
      DO j=js,je+1
        DO i=is,ie+1
          ke(i, j) = 0.5*(ke(i, j)+ub(i, j)*vb(i, j))
        END DO
      END DO
!-----------------------------------------
! Fix KE at the 4 corners of the face:
!-----------------------------------------
      IF (gnomonic_grid) THEN
        dt6 = dt/6.
        IF (sw_corner) ke(1, 1) = dt6*((ut(1, 1)+ut(1, 0))*u(1, 1)+(vt(1&
&           , 1)+vt(0, 1))*v(1, 1)+(ut(1, 1)+vt(1, 1))*u(0, 1))
        IF (se_corner) THEN
          i = npx
          ke(i, 1) = dt6*((ut(i, 1)+ut(i, 0))*u(i-1, 1)+(vt(i, 1)+vt(i-1&
&           , 1))*v(i, 1)+(ut(i, 1)-vt(i-1, 1))*u(i, 1))
        END IF
        IF (ne_corner) THEN
          i = npx
          j = npy
          ke(i, j) = dt6*((ut(i, j)+ut(i, j-1))*u(i-1, j)+(vt(i, j)+vt(i&
&           -1, j))*v(i, j-1)+(ut(i, j-1)+vt(i-1, j))*u(i, j))
        END IF
        IF (nw_corner) THEN
          j = npy
          ke(1, j) = dt6*((ut(1, j)+ut(1, j-1))*u(1, j)+(vt(1, j)+vt(0, &
&           j))*v(1, j-1)+(ut(1, j-1)-vt(1, j))*u(0, j))
        END IF
      ELSE IF (grid_type .LT. 3) THEN
!oncall mp_corner_comm(ke, npx, npy) 
        call mp_corner_comm(ke, npx, npy)
        IF (sw_corner) ke(1, 1) = r3*(ke(2, 1)+ke(1, 2)+ke(0, 1))
        IF (se_corner) ke(npx, 1) = r3*(ke(npx+1, 1)+ke(npx, 2)+ke(npx-1&
&           , 1))
        IF (ne_corner) ke(npx, npy) = r3*(ke(npx+1, npy)+ke(npx, npy-1)+&
&           ke(npx-1, npy))
        IF (nw_corner) ke(1, npy) = r3*(ke(2, npy)+ke(1, npy-1)+ke(0, &
&           npy))
      END IF
! Compute vorticity:
      DO j=jsd,jed+1
        DO i=isd,ied
          vt(i, j) = u(i, j)*dx(i, j)
        END DO
      END DO
      DO j=jsd,jed
        DO i=isd,ied+1
          ut(i, j) = v(i, j)*dy(i, j)
        END DO
      END DO
! wk is "volume-mean" relative vorticity
      DO j=jsd,jed
        DO i=isd,ied
          wk(i, j) = rarea(i, j)*(vt(i, j)-vt(i, j+1)-ut(i, j)+ut(i+1, j&
&           ))
        END DO
      END DO
!-----------------------------
! Compute divergence damping
!-----------------------------
      damp = dddmp*da_min_c
      IF (nord .EQ. 0) THEN
!       area ~ dxb*dyb*sin(alpha)
        DO j=js,je+1
          IF (j .EQ. 1 .OR. j .EQ. npy) THEN
            DO i=is-1,ie+1
              ptc(i, j) = u(i, j)*dyc(i, j)*sina_v(i, j)
            END DO
          ELSE
            DO i=is-1,ie+1
              ptc(i, j) = (u(i, j)-0.5*(va(i, j-1)+va(i, j))*cosa_v(i, j&
&               ))*dyc(i, j)*sina_v(i, j)
            END DO
          END IF
        END DO
        DO j=js-1,je+1
          DO i=is2,ie1
            vort(i, j) = (v(i, j)-0.5*(ua(i-1, j)+ua(i, j))*cosa_u(i, j)&
&             )*dxc(i, j)*sina_u(i, j)
          END DO
          IF (is .EQ. 1) vort(1, j) = v(1, j)*dxc(1, j)*sina_u(1, j)
          IF (ie + 1 .EQ. npx) vort(npx, j) = v(npx, j)*dxc(npx, j)*&
&             sina_u(npx, j)
        END DO
        DO j=js,je+1
          DO i=is,ie+1
            delpc(i, j) = vort(i, j-1) - vort(i, j) + ptc(i-1, j) - ptc(&
&             i, j)
          END DO
        END DO
! Remove the extra term at the corners:
        IF (sw_corner) delpc(1, 1) = delpc(1, 1) - vort(1, 0)
        IF (se_corner) delpc(npx, 1) = delpc(npx, 1) - vort(npx, 0)
        IF (ne_corner) delpc(npx, npy) = delpc(npx, npy) + vort(npx, npy&
&           )
        IF (nw_corner) delpc(1, npy) = delpc(1, npy) + vort(1, npy)
        DO j=js,je+1
          DO i=is,ie+1
            delpc(i, j) = rarea_c(i, j)*delpc(i, j)
            IF (delpc(i, j)*dt .GE. 0.) THEN
              abs0 = delpc(i, j)*dt
            ELSE
              abs0 = -(delpc(i, j)*dt)
            END IF
            y3 = dddmp*abs0
            IF (0.20 .GT. y3) THEN
              y1 = y3
            ELSE
              y1 = 0.20
            END IF
            IF (d2_bg .LT. y1) THEN
              max5 = y1
            ELSE
              max5 = d2_bg
            END IF
            damp = da_min_c*max5
            vort(i, j) = damp*delpc(i, j)
            ke(i, j) = ke(i, j) + vort(i, j)
          END DO
        END DO
      ELSE
!--------------------------
! Higher order divg damping
!--------------------------
        DO j=js,je+1
          DO i=is,ie+1
! Save divergence for external mode filter
            delpc(i, j) = divg_d(i, j)
          END DO
        END DO
        n2 = nord + 1
        DO n=1,nord
          nt = nord - n
          fill_c = nt .NE. 0 .AND. grid_type .LT. 3 .AND. (((sw_corner &
&           .OR. se_corner) .OR. ne_corner) .OR. nw_corner)
!onif ( fill_c ) call fill_corners(divg_d, npx, npy, FILL=XDir, BGRID=.true.)
          IF (fill_c) THEN
             call fill_corners(divg_d, npx, npy, FILL=XDir, BGRID=.true.)
          END IF
          DO j=js-nt,je+1+nt
            DO i=is-1-nt,ie+1+nt
              vc(i, j) = (divg_d(i+1, j)-divg_d(i, j))*divg_u(i, j)
            END DO
          END DO
!onif ( fill_c ) call fill_corners(divg_d, npx, npy, FILL=YDir, BGRID=.true.)
          IF (fill_c) THEN
             call fill_corners(divg_d, npx, npy, FILL=YDir, BGRID=.true.)
          END IF
          DO j=js-1-nt,je+1+nt
            DO i=is-nt,ie+1+nt
              uc(i, j) = (divg_d(i, j+1)-divg_d(i, j))*divg_v(i, j)
            END DO
          END DO
!onif ( fill_c ) call fill_corners(vc, uc, npx, npy, VECTOR=.true., DGRID=.true.)
          IF (fill_c) THEN
             call fill_corners(vc, uc, npx, npy, VECTOR=.true., DGRID=.true.)
          END IF
          DO j=js-nt,je+1+nt
            DO i=is-nt,ie+1+nt
              divg_d(i, j) = uc(i, j-1) - uc(i, j) + vc(i-1, j) - vc(i, &
&               j)
            END DO
          END DO
! Remove the extra term at the corners:
          IF (sw_corner) divg_d(1, 1) = divg_d(1, 1) - uc(1, 0)
          IF (se_corner) divg_d(npx, 1) = divg_d(npx, 1) - uc(npx, 0)
          IF (ne_corner) divg_d(npx, npy) = divg_d(npx, npy) + uc(npx, &
&             npy)
          IF (nw_corner) divg_d(1, npy) = divg_d(1, npy) + uc(1, npy)
          DO j=js-nt,je+1+nt
            DO i=is-nt,ie+1+nt
              divg_d(i, j) = divg_d(i, j)*rarea_c(i, j)
            END DO
          END DO
        END DO
        IF (dddmp .LT. 1.e-5) THEN
          vort = 0.
        ELSE
! Compute "time-scale" for del-2 background damping
!#ifdef FULL_SMAG
!        do j=js-1,je+1
!           do i=is-1,ie+1
!              vt2(i,j) = wk(i,j)**2
!           enddo
!        enddo
!        do j=js,je+1
!           do i=is,ie+1
!              vort(i,j) = dt*sqrt(delpc(i,j)**2 + 0.25*(vt2(i-1,j-1) +    &
!                                  vt2(i,j-1) + vt2(i-1,j) + vt2(i,j)))    
!           enddo
!        enddo
!        if (sw_corner) vort(1,1) = dt*sqrt( delpc(1,1)**2 +   &
!                       r3*(vt2(1,0) + vt2(0,1) + vt2(1,1)) )
!
!        if (se_corner) vort(npx,1) = dt*sqrt( delpc(npx,1)**2 +   &
!                       r3*(vt2(npx-1,0) + vt2(npx-1,1) + vt2(npx,1)) )
!
!        if (ne_corner) vort(npx,npy) = dt*sqrt( delpc(npx,npy)**2 +   &
!                       r3*(vt2(npx-1,npy-1) + vt2(npx,npy-1) + vt2(npx-1,npy)) )
!
!        if (nw_corner) vort(1,npy) = dt*sqrt( delpc(1,npy)**2 +   &
!                       r3*(vt2(0,npy-1) + vt2(1,npy-1) + vt2(1,npy)) )
!#else
          DO j=js,je+1
            DO i=is,ie+1
              IF (dt*delpc(i, j) .GE. 0.) THEN
                vort(i, j) = dt*delpc(i, j)
              ELSE
                vort(i, j) = -(dt*delpc(i, j))
              END IF
            END DO
          END DO
        END IF
!#endif
        dd8 = (da_min_c*d4_bg)**n2
        DO j=js,je+1
          DO i=is,ie+1
            IF (0.20 .GT. dddmp*vort(i, j)) THEN
              y2 = dddmp*vort(i, j)
            ELSE
              y2 = 0.20
            END IF
            IF (d2_bg .LT. y2) THEN
              max6 = y2
            ELSE
              max6 = d2_bg
            END IF
! del-2
            damp2 = da_min_c*max6
            vort(i, j) = damp2*delpc(i, j) + dd8*divg_d(i, j)
            ke(i, j) = ke(i, j) + vort(i, j)
          END DO
        END DO
      END IF
!----------------------------------
! Heating due to divergent damping:
!----------------------------------
      IF (d_con .GT. 1.e-5) THEN
!  damp = 0.5*0.25*d_con
        damp = 0.25*d_con
        DO j=js,je+1
          DO i=is,ie
! du
            ub(i, j) = (vort(i, j)-vort(i+1, j))*rdx(i, j)
! u*du
            gy(i, j) = u(i, j)*ub(i, j)
          END DO
        END DO
        DO j=js,je
          DO i=is,ie+1
! dv
            vb(i, j) = (vort(i, j)-vort(i, j+1))*rdy(i, j)
! v*dv
            gx(i, j) = v(i, j)*vb(i, j)
          END DO
        END DO
        DO j=js,je
          DO i=is,ie
            u2 = u(i, j) + u(i, j+1)
            du2 = ub(i, j) + ub(i, j+1)
            v2 = v(i, j) + v(i+1, j)
            dv2 = vb(i, j) + vb(i+1, j)
! Total energy conserving:
! Convert lost KE due to divergence damping to "heat"
            pt(i, j) = pt(i, j) - damp*rsin2(i, j)/pkz(i, j)*(ub(i, j)**&
&             2+ub(i, j+1)**2+vb(i, j)**2+vb(i+1, j)**2+2.*(gy(i, j)+gy(&
&             i, j+1)+gx(i, j)+gx(i+1, j))-cosa_s(i, j)*(u2*dv2+v2*du2+&
&             du2*dv2))
          END DO
        END DO
      END IF
! Vorticity transport
      DO j=jsd,jed
        DO i=isd,ied
          vort(i, j) = wk(i, j) + f0(i, j)
        END DO
      END DO
      CALL FV_TP_2D(vort, crx_adv, cry_adv, npx, npy, hord_vt, fx, fy, &
&             xfx_adv, yfx_adv, area, ra_x, ra_y, ppm_fac=ppm_limiter, &
&             nord=nord, damp_c=vtdm4)
      DO j=js,je+1
        DO i=is,ie
          u(i, j) = vt(i, j) + ke(i, j) - ke(i+1, j) + fy(i, j)
        END DO
      END DO
      DO j=js,je
        DO i=is,ie+1
          v(i, j) = ut(i, j) + ke(i, j) - ke(i, j+1) - fx(i, j)
        END DO
      END DO
    END IF
  END SUBROUTINE D_SW
!  Differentiation of divergence_corner in reverse (adjoint) mode (with options r8):
!   gradient     of useful results: u v ua va divg_d
!   with respect to varying inputs: u v ua va divg_d
  SUBROUTINE DIVERGENCE_CORNER_ADM(u, u_ad, v, v_ad, ua, ua_ad, va, &
&   va_ad, divg_d, divg_d_ad, km)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: km
    REAL, DIMENSION(isd:ied, jsd:jed + 1, km), INTENT(IN) :: u
    REAL, DIMENSION(isd:ied, jsd:jed+1, km) :: u_ad
    REAL, DIMENSION(isd:ied + 1, jsd:jed, km), INTENT(IN) :: v
    REAL, DIMENSION(isd:ied+1, jsd:jed, km) :: v_ad
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(IN) :: ua, va
    REAL, DIMENSION(isd:ied, jsd:jed, km) :: ua_ad, va_ad
    REAL, DIMENSION(isd:ied + 1, jsd:jed + 1, km) :: divg_d
    REAL, DIMENSION(isd:ied+1, jsd:jed+1, km) :: divg_d_ad
! local
    REAL :: uf(is-2:ie+2, js-1:je+2)
    REAL :: uf_ad(is-2:ie+2, js-1:je+2)
    REAL :: vf(is-1:ie+2, js-2:je+2)
    REAL :: vf_ad(is-1:ie+2, js-2:je+2)
    INTEGER :: i, j, k
    INTEGER :: is2, ie1
    INTRINSIC MAX
    INTRINSIC MIN
    INTEGER :: branch
    REAL :: temp_ad3
    REAL :: temp_ad2
    REAL :: temp_ad1
    REAL :: temp_ad0
    REAL :: temp_ad
    IF (2 .LT. is) THEN
      is2 = is
    ELSE
      is2 = 2
    END IF
    IF (npx - 1 .GT. ie + 1) THEN
      ie1 = ie + 1
    ELSE
      ie1 = npx - 1
    END IF
    DO k=1,km
      IF (grid_type .EQ. 4) THEN
        CALL PUSHCONTROL1B(1)
      ELSE
! divg_u(i,j) = sina_v(i,j)*dyc(i,j)/dx(i,j)
        DO j=js,je+1
          IF (j .EQ. 1 .OR. j .EQ. npy) THEN
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
        DO j=js-1,je+1
          IF (is .EQ. 1) THEN
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          IF (ie + 1 .EQ. npx) THEN
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
! Remove the extra term at the corners:
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
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
        CALL PUSHCONTROL1B(0)
      END IF
    END DO
    uf_ad = 0.0_8
    vf_ad = 0.0_8
    DO k=km,1,-1
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO j=je+1,js,-1
          DO i=ie+1,is,-1
            divg_d_ad(i, j, k) = rarea_c(i, j)*divg_d_ad(i, j, k)
          END DO
        END DO
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) vf_ad(1, npy) = vf_ad(1, npy) + divg_d_ad(1, &
&           npy, k)
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) vf_ad(npx, npy) = vf_ad(npx, npy) + divg_d_ad&
&           (npx, npy, k)
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) vf_ad(npx, 0) = vf_ad(npx, 0) - divg_d_ad(npx&
&           , 1, k)
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) vf_ad(1, 0) = vf_ad(1, 0) - divg_d_ad(1, 1, k&
&           )
        DO j=je+1,js,-1
          DO i=ie+1,is,-1
            vf_ad(i, j-1) = vf_ad(i, j-1) + divg_d_ad(i, j, k)
            vf_ad(i, j) = vf_ad(i, j) - divg_d_ad(i, j, k)
            uf_ad(i-1, j) = uf_ad(i-1, j) + divg_d_ad(i, j, k)
            uf_ad(i, j) = uf_ad(i, j) - divg_d_ad(i, j, k)
            divg_d_ad(i, j, k) = 0.0_8
          END DO
        END DO
        DO j=je+1,js-1,-1
          CALL POPCONTROL1B(branch)
          IF (branch .NE. 0) THEN
            v_ad(npx, j, k) = v_ad(npx, j, k) + dxc(npx, j)*sina_u(npx, &
&             j)*vf_ad(npx, j)
            vf_ad(npx, j) = 0.0_8
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            v_ad(1, j, k) = v_ad(1, j, k) + dxc(1, j)*sina_u(1, j)*vf_ad&
&             (1, j)
            vf_ad(1, j) = 0.0_8
          END IF
          DO i=ie1,is2,-1
            temp_ad2 = dxc(i, j)*sina_u(i, j)*vf_ad(i, j)
            temp_ad3 = -(cosa_u(i, j)*0.5*temp_ad2)
            v_ad(i, j, k) = v_ad(i, j, k) + temp_ad2
            ua_ad(i-1, j, k) = ua_ad(i-1, j, k) + temp_ad3
            ua_ad(i, j, k) = ua_ad(i, j, k) + temp_ad3
            vf_ad(i, j) = 0.0_8
          END DO
        END DO
        DO j=je+1,js,-1
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            DO i=ie+1,is-1,-1
              temp_ad0 = dyc(i, j)*sina_v(i, j)*uf_ad(i, j)
              temp_ad1 = -(cosa_v(i, j)*0.5*temp_ad0)
              u_ad(i, j, k) = u_ad(i, j, k) + temp_ad0
              va_ad(i, j-1, k) = va_ad(i, j-1, k) + temp_ad1
              va_ad(i, j, k) = va_ad(i, j, k) + temp_ad1
              uf_ad(i, j) = 0.0_8
            END DO
          ELSE
            DO i=ie+1,is-1,-1
              u_ad(i, j, k) = u_ad(i, j, k) + dyc(i, j)*sina_v(i, j)*&
&               uf_ad(i, j)
              uf_ad(i, j) = 0.0_8
            END DO
          END IF
        END DO
      ELSE
        DO j=je+2,js-1,-1
          DO i=ie+2,is-1,-1
            temp_ad = rarea_c(i, j)*divg_d_ad(i, j, k)
            vf_ad(i, j-1) = vf_ad(i, j-1) + temp_ad
            vf_ad(i, j) = vf_ad(i, j) - temp_ad
            uf_ad(i-1, j) = uf_ad(i-1, j) + temp_ad
            uf_ad(i, j) = uf_ad(i, j) - temp_ad
            divg_d_ad(i, j, k) = 0.0_8
          END DO
        END DO
        DO j=je+2,js-2,-1
          DO i=ie+2,is-1,-1
            v_ad(i, j, k) = v_ad(i, j, k) + dxc(i, j)*vf_ad(i, j)
            vf_ad(i, j) = 0.0_8
          END DO
        END DO
        DO j=je+2,js-1,-1
          DO i=ie+2,is-2,-1
            u_ad(i, j, k) = u_ad(i, j, k) + dyc(i, j)*uf_ad(i, j)
            uf_ad(i, j) = 0.0_8
          END DO
        END DO
      END IF
    END DO
  END SUBROUTINE DIVERGENCE_CORNER_ADM
  SUBROUTINE DIVERGENCE_CORNER(u, v, ua, va, divg_d, km)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: km
    REAL, DIMENSION(isd:ied, jsd:jed + 1, km), INTENT(IN) :: u
    REAL, DIMENSION(isd:ied + 1, jsd:jed, km), INTENT(IN) :: v
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(IN) :: ua, va
    REAL, DIMENSION(isd:ied + 1, jsd:jed + 1, km), INTENT(OUT) :: divg_d
! local
    REAL :: uf(is-2:ie+2, js-1:je+2)
    REAL :: vf(is-1:ie+2, js-2:je+2)
    INTEGER :: i, j, k
    INTEGER :: is2, ie1
    INTRINSIC MAX
    INTRINSIC MIN
    IF (2 .LT. is) THEN
      is2 = is
    ELSE
      is2 = 2
    END IF
    IF (npx - 1 .GT. ie + 1) THEN
      ie1 = ie + 1
    ELSE
      ie1 = npx - 1
    END IF
    DO k=1,km
      IF (grid_type .EQ. 4) THEN
        DO j=js-1,je+2
          DO i=is-2,ie+2
            uf(i, j) = u(i, j, k)*dyc(i, j)
          END DO
        END DO
        DO j=js-2,je+2
          DO i=is-1,ie+2
            vf(i, j) = v(i, j, k)*dxc(i, j)
          END DO
        END DO
        DO j=js-1,je+2
          DO i=is-1,ie+2
            divg_d(i, j, k) = rarea_c(i, j)*(vf(i, j-1)-vf(i, j)+uf(i-1&
&             , j)-uf(i, j))
          END DO
        END DO
      ELSE
! divg_u(i,j) = sina_v(i,j)*dyc(i,j)/dx(i,j)
        DO j=js,je+1
          IF (j .EQ. 1 .OR. j .EQ. npy) THEN
            DO i=is-1,ie+1
              uf(i, j) = u(i, j, k)*dyc(i, j)*sina_v(i, j)
            END DO
          ELSE
            DO i=is-1,ie+1
              uf(i, j) = (u(i, j, k)-0.5*(va(i, j-1, k)+va(i, j, k))*&
&               cosa_v(i, j))*dyc(i, j)*sina_v(i, j)
            END DO
          END IF
        END DO
        DO j=js-1,je+1
          DO i=is2,ie1
            vf(i, j) = (v(i, j, k)-0.5*(ua(i-1, j, k)+ua(i, j, k))*&
&             cosa_u(i, j))*dxc(i, j)*sina_u(i, j)
          END DO
          IF (is .EQ. 1) vf(1, j) = v(1, j, k)*dxc(1, j)*sina_u(1, j)
          IF (ie + 1 .EQ. npx) vf(npx, j) = v(npx, j, k)*dxc(npx, j)*&
&             sina_u(npx, j)
        END DO
        DO j=js,je+1
          DO i=is,ie+1
            divg_d(i, j, k) = vf(i, j-1) - vf(i, j) + uf(i-1, j) - uf(i&
&             , j)
          END DO
        END DO
! Remove the extra term at the corners:
        IF (sw_corner) divg_d(1, 1, k) = divg_d(1, 1, k) - vf(1, 0)
        IF (se_corner) divg_d(npx, 1, k) = divg_d(npx, 1, k) - vf(npx, 0&
&           )
        IF (ne_corner) divg_d(npx, npy, k) = divg_d(npx, npy, k) + vf(&
&           npx, npy)
        IF (nw_corner) divg_d(1, npy, k) = divg_d(1, npy, k) + vf(1, npy&
&           )
        DO j=js,je+1
          DO i=is,ie+1
            divg_d(i, j, k) = rarea_c(i, j)*divg_d(i, j, k)
          END DO
        END DO
      END IF
    END DO
  END SUBROUTINE DIVERGENCE_CORNER
!  Differentiation of xtp_u in reverse (adjoint) mode (with options r8):
!   gradient     of useful results: flux u c
!   with respect to varying inputs: flux u c
  SUBROUTINE XTP_U_ADM(c, c_ad, u, u_ad, v, flux, flux_ad, iord)
    IMPLICIT NONE
    REAL, INTENT(IN) :: u(isd:ied, jsd:jed+1)
    REAL :: u_ad(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: v(isd:ied+1, jsd:jed)
    REAL(ke_precision), INTENT(IN) :: c(is:ie+1, js:je+1)
    REAL(ke_precision) :: c_ad(is:ie+1, js:je+1)
    REAL(ke_precision) :: flux(is:ie+1, js:je+1)
    REAL(ke_precision) :: flux_ad(is:ie+1, js:je+1)
    INTEGER, INTENT(IN) :: iord
! Local
    REAL(ke_precision) :: al(is-1:ie+2), dm(is-2:ie+2)
    REAL(ke_precision) :: bl(is-1:ie+1)
    REAL(ke_precision) :: bl_ad(is-1:ie+1)
    REAL(ke_precision) :: br(is-1:ie+1)
    REAL(ke_precision) :: br_ad(is-1:ie+1)
    REAL(ke_precision) :: dq(is-3:ie+2)
    REAL(ke_precision) :: dl, dr, xt, pmp, lac, dqt, cfl
    REAL(ke_precision) :: xt_ad, cfl_ad
    REAL(ke_precision) :: x0, x1
    INTEGER :: i, j
    INTRINSIC MAX
    INTRINSIC MIN
    INTEGER :: branch
    INTEGER :: ad_from
    INTEGER :: ad_to
    INTEGER :: min1
    REAL :: temp_ad9
    REAL :: temp_ad8
    REAL :: temp_ad7
    REAL :: temp_ad6
    REAL :: temp_ad5
    REAL :: temp_ad4
    REAL :: temp_ad3
    REAL :: temp_ad2
    REAL :: temp_ad1
    REAL :: temp_ad0
    REAL :: temp_ad
    REAL(ke_precision) :: temp_ad12
    REAL(ke_precision) :: temp_ad11
    REAL :: temp_ad10
    REAL(ke_precision) :: temp
    INTEGER :: max1
    SELECT CASE  (iord) 
    CASE (1) 
      DO j=js,je+1
        DO i=is,ie+1
          IF (c(i, j) .GT. 0.) THEN
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
      END DO
      DO j=je+1,js,-1
        DO i=ie+1,is,-1
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            u_ad(i, j) = u_ad(i, j) + flux_ad(i, j)
            flux_ad(i, j) = 0.0_8
          ELSE
            u_ad(i-1, j) = u_ad(i-1, j) + flux_ad(i, j)
            flux_ad(i, j) = 0.0_8
          END IF
        END DO
      END DO
    CASE (333) 
      DO j=js,je+1
        DO i=is,ie+1
!Apply first order at the edges
          IF (i .EQ. is .OR. i .EQ. ie + 1) THEN
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
      DO j=je+1,js,-1
        DO i=ie+1,is,-1
          CALL POPCONTROL2B(branch)
          IF (branch .LT. 2) THEN
            IF (branch .EQ. 0) THEN
              cfl = c(i, j)*rdx(i, j)
              temp_ad2 = flux_ad(i, j)/6.0
              temp_ad3 = -(0.5*cfl*flux_ad(i, j))
              temp_ad4 = cfl**2*flux_ad(i, j)/6.0
              u_ad(i-1, j) = u_ad(i-1, j) + temp_ad4 - temp_ad3 + 2.0*&
&               temp_ad2
              u_ad(i, j) = u_ad(i, j) + temp_ad3 - 2.0*temp_ad4 + 5.0*&
&               temp_ad2
              u_ad(i+1, j) = u_ad(i+1, j) + temp_ad4 - temp_ad2
              cfl_ad = ((u(i+1, j)-2.0*u(i, j)+u(i-1, j))*2*cfl/6.0-0.5*&
&               (u(i, j)-u(i-1, j)))*flux_ad(i, j)
              flux_ad(i, j) = 0.0_8
              c_ad(i, j) = c_ad(i, j) + rdx(i, j)*cfl_ad
            ELSE
              cfl = c(i, j)*rdx(i-1, j)
              temp_ad = flux_ad(i, j)/6.0
              temp_ad0 = -(0.5*cfl*flux_ad(i, j))
              temp_ad1 = cfl**2*flux_ad(i, j)/6.0
              u_ad(i, j) = u_ad(i, j) + temp_ad1 + temp_ad0 + 2.0*&
&               temp_ad
              u_ad(i-1, j) = u_ad(i-1, j) + 5.0*temp_ad - temp_ad0 - 2.0&
&               *temp_ad1
              u_ad(i-2, j) = u_ad(i-2, j) + temp_ad1 - temp_ad
              cfl_ad = ((u(i, j)-2.0*u(i-1, j)+u(i-2, j))*2*cfl/6.0-0.5*&
&               (u(i, j)-u(i-1, j)))*flux_ad(i, j)
              flux_ad(i, j) = 0.0_8
              c_ad(i, j) = c_ad(i, j) + rdx(i-1, j)*cfl_ad
            END IF
          ELSE IF (branch .EQ. 2) THEN
            u_ad(i, j) = u_ad(i, j) + flux_ad(i, j)
            flux_ad(i, j) = 0.0_8
          ELSE
            u_ad(i-1, j) = u_ad(i-1, j) + flux_ad(i, j)
            flux_ad(i, j) = 0.0_8
          END IF
        END DO
      END DO
    CASE (6) 
      DO j=js,je+1
        IF (3 .LT. is - 1) THEN
          max1 = is - 1
        ELSE
          max1 = 3
        END IF
        IF (npx - 3 .GT. ie + 1) THEN
          min1 = ie + 1
        ELSE
          min1 = npx - 3
        END IF
        ad_from = max1
        DO i=ad_from,min1
          CALL PUSHREAL8(bl(i))
          bl(i) = b5*u(i-2, j) + b4*u(i-1, j) + b3*u(i, j) + b2*u(i+1, j&
&           ) + b1*u(i+2, j)
          CALL PUSHREAL8(br(i))
          br(i) = b1*u(i-2, j) + b2*u(i-1, j) + b3*u(i, j) + b4*u(i+1, j&
&           ) + b5*u(i+2, j)
        END DO
        CALL PUSHINTEGER4(i - 1)
        CALL PUSHINTEGER4(ad_from)
        IF (grid_type .LT. 3) THEN
          IF (is .EQ. 1) THEN
            CALL PUSHREAL8(br(2))
            br(2) = p1*(u(2, j)+u(3, j)) + p2*(u(1, j)+u(4, j)) - u(2, j&
&             )
            xt = c3*u(1, j) + c2*u(2, j) + c1*u(3, j)
            CALL PUSHREAL8(bl(2))
            bl(2) = xt - u(2, j)
            IF (j .EQ. 1 .OR. j .EQ. npy) THEN
! out
              CALL PUSHREAL8(bl(0))
              bl(0) = 0.
! edge
              CALL PUSHREAL8(br(0))
              br(0) = 0.
! edge
              CALL PUSHREAL8(bl(1))
              bl(1) = 0.
! in
              CALL PUSHREAL8(br(1))
              br(1) = 0.
              CALL PUSHCONTROL2B(0)
            ELSE
              CALL PUSHREAL8(br(1))
              br(1) = xt - u(1, j)
              xt = 0.5*((2.*dx(1, j)+dx(2, j))*(u(0, j)+u(1, j))-dx(1, j&
&               )*(u(-1, j)+u(2, j)))/(dx(1, j)+dx(2, j))
              CALL PUSHREAL8(bl(1))
              bl(1) = xt - u(1, j)
              CALL PUSHREAL8(br(0))
              br(0) = xt - u(0, j)
              xt = c1*u(-2, j) + c2*u(-1, j) + c3*u(0, j)
              CALL PUSHREAL8(bl(0))
              bl(0) = xt - u(0, j)
              CALL PUSHCONTROL2B(1)
            END IF
          ELSE
            CALL PUSHCONTROL2B(2)
          END IF
          IF (ie + 1 .EQ. npx) THEN
            CALL PUSHREAL8(bl(npx-2))
            bl(npx-2) = p1*(u(npx-2, j)+u(npx-3, j)) + p2*(u(npx-4, j)+u&
&             (npx-1, j)) - u(npx-2, j)
            xt = c1*u(npx-3, j) + c2*u(npx-2, j) + c3*u(npx-1, j)
            CALL PUSHREAL8(br(npx-2))
            br(npx-2) = xt - u(npx-2, j)
            IF (j .EQ. 1 .OR. j .EQ. npy) THEN
! in
              CALL PUSHREAL8(bl(npx-1))
              bl(npx-1) = 0.
! edge
              CALL PUSHREAL8(br(npx-1))
              br(npx-1) = 0.
! edge
              CALL PUSHREAL8(bl(npx))
              bl(npx) = 0.
! out
              CALL PUSHREAL8(br(npx))
              br(npx) = 0.
              CALL PUSHCONTROL2B(3)
            ELSE
              CALL PUSHREAL8(bl(npx-1))
              bl(npx-1) = xt - u(npx-1, j)
              xt = 0.5*((2.*dx(npx-1, j)+dx(npx-2, j))*(u(npx-1, j)+u(&
&               npx, j))-dx(npx-1, j)*(u(npx-2, j)+u(npx+1, j)))/(dx(npx&
&               -1, j)+dx(npx-2, j))
              CALL PUSHREAL8(br(npx-1))
              br(npx-1) = xt - u(npx-1, j)
              CALL PUSHREAL8(bl(npx))
              bl(npx) = xt - u(npx, j)
              xt = c3*u(npx, j) + c2*u(npx+1, j) + c1*u(npx+2, j)
              CALL PUSHREAL8(br(npx))
              br(npx) = xt - u(npx, j)
              CALL PUSHCONTROL2B(2)
            END IF
          ELSE
            CALL PUSHCONTROL2B(1)
          END IF
        ELSE
          CALL PUSHCONTROL2B(0)
        END IF
        DO i=is,ie+1
          IF (c(i, j) .GT. 0.) THEN
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
      END DO
      bl_ad = 0.0_8
      br_ad = 0.0_8
      DO j=je+1,js,-1
        DO i=ie+1,is,-1
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            cfl = c(i, j)*rdx(i, j)
            temp_ad12 = (cfl+1.)*flux_ad(i, j)
            u_ad(i, j) = u_ad(i, j) + flux_ad(i, j)
            cfl_ad = (bl(i)+br(i))*temp_ad12 + (bl(i)+cfl*(bl(i)+br(i)))&
&             *flux_ad(i, j)
            bl_ad(i) = bl_ad(i) + (cfl+1.0_8)*temp_ad12
            br_ad(i) = br_ad(i) + cfl*temp_ad12
            flux_ad(i, j) = 0.0_8
            c_ad(i, j) = c_ad(i, j) + rdx(i, j)*cfl_ad
          ELSE
            cfl = c(i, j)*rdx(i-1, j)
            temp = bl(i-1) + br(i-1)
            temp_ad11 = (1.-cfl)*flux_ad(i, j)
            u_ad(i-1, j) = u_ad(i-1, j) + flux_ad(i, j)
            cfl_ad = -(temp*temp_ad11) - (br(i-1)-cfl*temp)*flux_ad(i, j&
&             )
            br_ad(i-1) = br_ad(i-1) + (1.0_8-cfl)*temp_ad11
            bl_ad(i-1) = bl_ad(i-1) - cfl*temp_ad11
            flux_ad(i, j) = 0.0_8
            c_ad(i, j) = c_ad(i, j) + rdx(i-1, j)*cfl_ad
          END IF
        END DO
        CALL POPCONTROL2B(branch)
        IF (branch .LT. 2) THEN
          IF (branch .EQ. 0) GOTO 100
        ELSE
          IF (branch .EQ. 2) THEN
            CALL POPREAL8(br(npx))
            xt_ad = br_ad(npx)
            u_ad(npx, j) = u_ad(npx, j) + c3*xt_ad - br_ad(npx)
            br_ad(npx) = 0.0_8
            u_ad(npx+1, j) = u_ad(npx+1, j) + c2*xt_ad
            u_ad(npx+2, j) = u_ad(npx+2, j) + c1*xt_ad
            CALL POPREAL8(bl(npx))
            xt_ad = br_ad(npx-1) + bl_ad(npx)
            u_ad(npx, j) = u_ad(npx, j) - bl_ad(npx)
            bl_ad(npx) = 0.0_8
            CALL POPREAL8(br(npx-1))
            temp_ad9 = 0.5*xt_ad/(dx(npx-1, j)+dx(npx-2, j))
            temp_ad8 = (dx(npx-1, j)*2.+dx(npx-2, j))*temp_ad9
            u_ad(npx-1, j) = u_ad(npx-1, j) + temp_ad8 - br_ad(npx-1)
            br_ad(npx-1) = 0.0_8
            temp_ad10 = -(dx(npx-1, j)*temp_ad9)
            u_ad(npx, j) = u_ad(npx, j) + temp_ad8
            u_ad(npx-2, j) = u_ad(npx-2, j) + temp_ad10
            u_ad(npx+1, j) = u_ad(npx+1, j) + temp_ad10
            CALL POPREAL8(bl(npx-1))
            xt_ad = bl_ad(npx-1)
            u_ad(npx-1, j) = u_ad(npx-1, j) - bl_ad(npx-1)
            bl_ad(npx-1) = 0.0_8
          ELSE
            CALL POPREAL8(br(npx))
            br_ad(npx) = 0.0_8
            CALL POPREAL8(bl(npx))
            bl_ad(npx) = 0.0_8
            CALL POPREAL8(br(npx-1))
            br_ad(npx-1) = 0.0_8
            CALL POPREAL8(bl(npx-1))
            bl_ad(npx-1) = 0.0_8
            xt_ad = 0.0_8
          END IF
          CALL POPREAL8(br(npx-2))
          xt_ad = xt_ad + br_ad(npx-2)
          u_ad(npx-2, j) = u_ad(npx-2, j) - br_ad(npx-2)
          br_ad(npx-2) = 0.0_8
          u_ad(npx-3, j) = u_ad(npx-3, j) + c1*xt_ad
          u_ad(npx-2, j) = u_ad(npx-2, j) + c2*xt_ad
          u_ad(npx-1, j) = u_ad(npx-1, j) + c3*xt_ad
          CALL POPREAL8(bl(npx-2))
          u_ad(npx-2, j) = u_ad(npx-2, j) + (p1-1.0)*bl_ad(npx-2)
          u_ad(npx-3, j) = u_ad(npx-3, j) + p1*bl_ad(npx-2)
          u_ad(npx-4, j) = u_ad(npx-4, j) + p2*bl_ad(npx-2)
          u_ad(npx-1, j) = u_ad(npx-1, j) + p2*bl_ad(npx-2)
          bl_ad(npx-2) = 0.0_8
        END IF
        CALL POPCONTROL2B(branch)
        IF (branch .EQ. 0) THEN
          CALL POPREAL8(br(1))
          br_ad(1) = 0.0_8
          CALL POPREAL8(bl(1))
          bl_ad(1) = 0.0_8
          CALL POPREAL8(br(0))
          br_ad(0) = 0.0_8
          CALL POPREAL8(bl(0))
          bl_ad(0) = 0.0_8
          xt_ad = 0.0_8
        ELSE IF (branch .EQ. 1) THEN
          CALL POPREAL8(bl(0))
          xt_ad = bl_ad(0)
          u_ad(0, j) = u_ad(0, j) - bl_ad(0)
          bl_ad(0) = 0.0_8
          u_ad(-2, j) = u_ad(-2, j) + c1*xt_ad
          u_ad(-1, j) = u_ad(-1, j) + c2*xt_ad
          u_ad(0, j) = u_ad(0, j) + c3*xt_ad - br_ad(0)
          CALL POPREAL8(br(0))
          xt_ad = bl_ad(1) + br_ad(0)
          br_ad(0) = 0.0_8
          CALL POPREAL8(bl(1))
          u_ad(1, j) = u_ad(1, j) - bl_ad(1)
          bl_ad(1) = 0.0_8
          temp_ad5 = 0.5*xt_ad/(dx(1, j)+dx(2, j))
          temp_ad6 = (dx(1, j)*2.+dx(2, j))*temp_ad5
          temp_ad7 = -(dx(1, j)*temp_ad5)
          u_ad(0, j) = u_ad(0, j) + temp_ad6
          u_ad(1, j) = u_ad(1, j) + temp_ad6
          u_ad(-1, j) = u_ad(-1, j) + temp_ad7
          u_ad(2, j) = u_ad(2, j) + temp_ad7
          CALL POPREAL8(br(1))
          xt_ad = br_ad(1)
          u_ad(1, j) = u_ad(1, j) - br_ad(1)
          br_ad(1) = 0.0_8
        ELSE
          GOTO 100
        END IF
        CALL POPREAL8(bl(2))
        xt_ad = xt_ad + bl_ad(2)
        u_ad(2, j) = u_ad(2, j) - bl_ad(2)
        bl_ad(2) = 0.0_8
        u_ad(1, j) = u_ad(1, j) + c3*xt_ad
        u_ad(2, j) = u_ad(2, j) + c2*xt_ad
        u_ad(3, j) = u_ad(3, j) + c1*xt_ad
        CALL POPREAL8(br(2))
        u_ad(2, j) = u_ad(2, j) + (p1-1.0)*br_ad(2)
        u_ad(3, j) = u_ad(3, j) + p1*br_ad(2)
        u_ad(1, j) = u_ad(1, j) + p2*br_ad(2)
        u_ad(4, j) = u_ad(4, j) + p2*br_ad(2)
        br_ad(2) = 0.0_8
 100    CALL POPINTEGER4(ad_from)
        CALL POPINTEGER4(ad_to)
        DO i=ad_to,ad_from,-1
          CALL POPREAL8(br(i))
          u_ad(i-2, j) = u_ad(i-2, j) + b1*br_ad(i)
          u_ad(i-1, j) = u_ad(i-1, j) + b2*br_ad(i)
          u_ad(i, j) = u_ad(i, j) + b3*br_ad(i)
          u_ad(i+1, j) = u_ad(i+1, j) + b4*br_ad(i)
          u_ad(i+2, j) = u_ad(i+2, j) + b5*br_ad(i)
          br_ad(i) = 0.0_8
          CALL POPREAL8(bl(i))
          u_ad(i-2, j) = u_ad(i-2, j) + b5*bl_ad(i)
          u_ad(i-1, j) = u_ad(i-1, j) + b4*bl_ad(i)
          u_ad(i, j) = u_ad(i, j) + b3*bl_ad(i)
          u_ad(i+1, j) = u_ad(i+1, j) + b2*bl_ad(i)
          u_ad(i+2, j) = u_ad(i+2, j) + b1*bl_ad(i)
          bl_ad(i) = 0.0_8
        END DO
      END DO
    END SELECT
  END SUBROUTINE XTP_U_ADM
!  SUBROUTINE XTP_U(c, u, v, flux, iord)
!    IMPLICIT NONE
!    REAL, INTENT(IN) :: u(isd:ied, jsd:jed+1)
!    REAL, INTENT(IN) :: v(isd:ied+1, jsd:jed)
!    REAL(ke_precision), INTENT(IN) :: c(is:ie+1, js:je+1)
!    REAL(ke_precision), INTENT(OUT) :: flux(is:ie+1, js:je+1)
!    INTEGER, INTENT(IN) :: iord
!! Local
!    REAL(ke_precision) :: al(is-1:ie+2), dm(is-2:ie+2)
!    REAL(ke_precision) :: bl(is-1:ie+1)
!    REAL(ke_precision) :: br(is-1:ie+1)
!    REAL(ke_precision) :: dq(is-3:ie+2)
!    REAL(ke_precision) :: dl, dr, xt, pmp, lac, dqt, cfl
!    REAL(ke_precision) :: x0, x1
!    INTEGER :: i, j
!    INTRINSIC MAX
!    INTRINSIC MIN
!    INTEGER :: min1
!    INTEGER :: max1
!    SELECT CASE  (iord) 
!    CASE (1) 
!      DO j=js,je+1
!        DO i=is,ie+1
!          IF (c(i, j) .GT. 0.) THEN
!            flux(i, j) = u(i-1, j)
!          ELSE
!            flux(i, j) = u(i, j)
!          END IF
!        END DO
!      END DO
!    CASE (333) 
!      DO j=js,je+1
!        DO i=is,ie+1
!!Apply first order at the edges
!          IF (i .EQ. is .OR. i .EQ. ie + 1) THEN
!            IF (c(i, j) .GT. 0.) THEN
!              flux(i, j) = u(i-1, j)
!            ELSE
!              flux(i, j) = u(i, j)
!            END IF
!          ELSE IF (c(i, j) .GT. 0.) THEN
!!Otherwise use the third order scheme
!            cfl = c(i, j)*rdx(i-1, j)
!            flux(i, j) = (2.0*u(i, j)+5.0*u(i-1, j)-u(i-2, j))/6.0 - 0.5&
!&             *cfl*(u(i, j)-u(i-1, j)) + cfl*cfl/6.0*(u(i, j)-2.0*u(i-1&
!&             , j)+u(i-2, j))
!          ELSE
!            cfl = c(i, j)*rdx(i, j)
!            flux(i, j) = (2.0*u(i-1, j)+5.0*u(i, j)-u(i+1, j))/6.0 - 0.5&
!&             *cfl*(u(i, j)-u(i-1, j)) + cfl*cfl/6.0*(u(i+1, j)-2.0*u(i&
!&             , j)+u(i-1, j))
!          END IF
!        END DO
!      END DO
!    CASE (6) 
!      DO j=js,je+1
!        IF (3 .LT. is - 1) THEN
!          max1 = is - 1
!        ELSE
!          max1 = 3
!        END IF
!        IF (npx - 3 .GT. ie + 1) THEN
!          min1 = ie + 1
!        ELSE
!          min1 = npx - 3
!        END IF
!        DO i=max1,min1
!          bl(i) = b5*u(i-2, j) + b4*u(i-1, j) + b3*u(i, j) + b2*u(i+1, j&
!&           ) + b1*u(i+2, j)
!          br(i) = b1*u(i-2, j) + b2*u(i-1, j) + b3*u(i, j) + b4*u(i+1, j&
!&           ) + b5*u(i+2, j)
!        END DO
!        IF (grid_type .LT. 3) THEN
!          IF (is .EQ. 1) THEN
!            br(2) = p1*(u(2, j)+u(3, j)) + p2*(u(1, j)+u(4, j)) - u(2, j&
!&             )
!            xt = c3*u(1, j) + c2*u(2, j) + c1*u(3, j)
!            bl(2) = xt - u(2, j)
!            IF (j .EQ. 1 .OR. j .EQ. npy) THEN
!! out
!              bl(0) = 0.
!! edge
!              br(0) = 0.
!! edge
!              bl(1) = 0.
!! in
!              br(1) = 0.
!            ELSE
!              br(1) = xt - u(1, j)
!              xt = 0.5*((2.*dx(1, j)+dx(2, j))*(u(0, j)+u(1, j))-dx(1, j&
!&               )*(u(-1, j)+u(2, j)))/(dx(1, j)+dx(2, j))
!              bl(1) = xt - u(1, j)
!              br(0) = xt - u(0, j)
!              xt = c1*u(-2, j) + c2*u(-1, j) + c3*u(0, j)
!              bl(0) = xt - u(0, j)
!            END IF
!          END IF
!          IF (ie + 1 .EQ. npx) THEN
!            bl(npx-2) = p1*(u(npx-2, j)+u(npx-3, j)) + p2*(u(npx-4, j)+u&
!&             (npx-1, j)) - u(npx-2, j)
!            xt = c1*u(npx-3, j) + c2*u(npx-2, j) + c3*u(npx-1, j)
!            br(npx-2) = xt - u(npx-2, j)
!            IF (j .EQ. 1 .OR. j .EQ. npy) THEN
!! in
!              bl(npx-1) = 0.
!! edge
!              br(npx-1) = 0.
!! edge
!              bl(npx) = 0.
!! out
!              br(npx) = 0.
!            ELSE
!              bl(npx-1) = xt - u(npx-1, j)
!              xt = 0.5*((2.*dx(npx-1, j)+dx(npx-2, j))*(u(npx-1, j)+u(&
!&               npx, j))-dx(npx-1, j)*(u(npx-2, j)+u(npx+1, j)))/(dx(npx&
!&               -1, j)+dx(npx-2, j))
!              br(npx-1) = xt - u(npx-1, j)
!              bl(npx) = xt - u(npx, j)
!              xt = c3*u(npx, j) + c2*u(npx+1, j) + c1*u(npx+2, j)
!              br(npx) = xt - u(npx, j)
!            END IF
!          END IF
!        END IF
!        DO i=is,ie+1
!          IF (c(i, j) .GT. 0.) THEN
!            cfl = c(i, j)*rdx(i-1, j)
!            flux(i, j) = u(i-1, j) + (1.-cfl)*(br(i-1)-cfl*(bl(i-1)+br(i&
!&             -1)))
!          ELSE
!            cfl = c(i, j)*rdx(i, j)
!            flux(i, j) = u(i, j) + (1.+cfl)*(bl(i)+cfl*(bl(i)+br(i)))
!          END IF
!        END DO
!      END DO
!    END SELECT
!  END SUBROUTINE XTP_U
!  Differentiation of ytp_v in reverse (adjoint) mode:
!   gradient     of useful results: flux v c
!   with respect to varying inputs: v c
  SUBROUTINE YTP_V_ADM(c, c_ad, u, v, v_ad, flux, flux_ad, jord)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: jord
    REAL, INTENT(IN) :: u(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: v(isd:ied+1, jsd:jed)
    REAL :: v_ad(isd:ied+1, jsd:jed)
!  Courant   N (like FLUX)
    REAL(ke_precision), INTENT(IN) :: c(is:ie+1, js:je+1)
    REAL(ke_precision) :: c_ad(is:ie+1, js:je+1)
    REAL(ke_precision) :: flux(is:ie+1, js:je+1)
    REAL(ke_precision) :: flux_ad(is:ie+1, js:je+1)
! Local:
    REAL(ke_precision) :: dm(is:ie+1, js-2:je+2)
    REAL(ke_precision) :: al(is:ie+1, js-1:je+2)
    REAL(ke_precision) :: bl(is:ie+1, js-1:je+1)
    REAL(ke_precision) :: bl_ad(is:ie+1, js-1:je+1)
    REAL(ke_precision) :: br(is:ie+1, js-1:je+1)
    REAL(ke_precision) :: br_ad(is:ie+1, js-1:je+1)
    REAL(ke_precision) :: dq(is:ie+1, js-3:je+2)
    REAL(ke_precision) :: xt, dl, dr, pmp, lac, dqt, cfl
    REAL(ke_precision) :: xt_ad, cfl_ad
    REAL(ke_precision) :: x0, x1
    INTEGER :: i, j
    INTRINSIC MAX
    INTRINSIC MIN
    INTEGER :: branch
    REAL(ke_precision) :: temp0
    INTEGER :: min1
    REAL :: temp_ad9
    REAL :: temp_ad8
    REAL :: temp_ad7
    REAL :: temp_ad6
    REAL :: temp_ad5
    REAL :: temp_ad4
    REAL :: temp_ad3
    REAL :: temp_ad2
    REAL :: temp_ad1
    REAL :: temp_ad0
    REAL :: temp_ad
    REAL(ke_precision) :: temp_ad12
    REAL(ke_precision) :: temp_ad11
    REAL :: temp_ad10
    REAL(ke_precision) :: temp
    INTEGER :: max1
    SELECT CASE  (jord) 
    CASE (1) 
      DO j=js,je+1
        DO i=is,ie+1
          IF (c(i, j) .GT. 0.) THEN
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
      END DO
      DO j=je+1,js,-1
        DO i=ie+1,is,-1
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            v_ad(i, j) = v_ad(i, j) + flux_ad(i, j)
            flux_ad(i, j) = 0.0_8
          ELSE
            v_ad(i, j-1) = v_ad(i, j-1) + flux_ad(i, j)
            flux_ad(i, j) = 0.0_8
          END IF
        END DO
      END DO
    CASE (333) 
      DO j=js,je+1
        DO i=is,ie+1
!Apply first order at the edges
          IF (j .EQ. js .OR. j .EQ. je + 1) THEN
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
      DO j=je+1,js,-1
        DO i=ie+1,is,-1
          CALL POPCONTROL2B(branch)
          IF (branch .LT. 2) THEN
            IF (branch .EQ. 0) THEN
              cfl = c(i, j)*rdy(i, j)
              temp_ad2 = flux_ad(i, j)/6.0
              temp_ad3 = -(0.5*cfl*flux_ad(i, j))
              temp_ad4 = cfl**2*flux_ad(i, j)/6.0
              v_ad(i, j-1) = v_ad(i, j-1) + temp_ad4 - temp_ad3 + 2.0*&
&               temp_ad2
              v_ad(i, j) = v_ad(i, j) + temp_ad3 - 2.0*temp_ad4 + 5.0*&
&               temp_ad2
              v_ad(i, j+1) = v_ad(i, j+1) + temp_ad4 - temp_ad2
              cfl_ad = ((v(i, j+1)-2.0*v(i, j)+v(i, j-1))*2*cfl/6.0-0.5*&
&               (v(i, j)-v(i, j-1)))*flux_ad(i, j)
              flux_ad(i, j) = 0.0_8
              c_ad(i, j) = c_ad(i, j) + rdy(i, j)*cfl_ad
            ELSE
              cfl = c(i, j)*rdy(i, j-1)
              temp_ad = flux_ad(i, j)/6.0
              temp_ad0 = -(0.5*cfl*flux_ad(i, j))
              temp_ad1 = cfl**2*flux_ad(i, j)/6.0
              v_ad(i, j) = v_ad(i, j) + temp_ad1 + temp_ad0 + 2.0*&
&               temp_ad
              v_ad(i, j-1) = v_ad(i, j-1) + 5.0*temp_ad - temp_ad0 - 2.0&
&               *temp_ad1
              v_ad(i, j-2) = v_ad(i, j-2) + temp_ad1 - temp_ad
              cfl_ad = ((v(i, j)-2.0*v(i, j-1)+v(i, j-2))*2*cfl/6.0-0.5*&
&               (v(i, j)-v(i, j-1)))*flux_ad(i, j)
              flux_ad(i, j) = 0.0_8
              c_ad(i, j) = c_ad(i, j) + rdy(i, j-1)*cfl_ad
            END IF
          ELSE IF (branch .EQ. 2) THEN
            v_ad(i, j) = v_ad(i, j) + flux_ad(i, j)
            flux_ad(i, j) = 0.0_8
          ELSE
            v_ad(i, j-1) = v_ad(i, j-1) + flux_ad(i, j)
            flux_ad(i, j) = 0.0_8
          END IF
        END DO
      END DO
    CASE (6) 
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
      DO j=max1,min1
        DO i=is,ie+1
          bl(i, j) = b5*v(i, j-2) + b4*v(i, j-1) + b3*v(i, j) + b2*v(i, &
&           j+1) + b1*v(i, j+2)
          br(i, j) = b1*v(i, j-2) + b2*v(i, j-1) + b3*v(i, j) + b4*v(i, &
&           j+1) + b5*v(i, j+2)
        END DO
      END DO
      IF (grid_type .LT. 3) THEN
        IF (js .EQ. 1) THEN
          DO i=is,ie+1
            br(i, 2) = p1*(v(i, 2)+v(i, 3)) + p2*(v(i, 1)+v(i, 4)) - v(i&
&             , 2)
            xt = c3*v(i, 1) + c2*v(i, 2) + c1*v(i, 3)
            br(i, 1) = xt - v(i, 1)
            bl(i, 2) = xt - v(i, 2)
            bl(i, 0) = c1*v(i, -2) + c2*v(i, -1) + c3*v(i, 0) - v(i, 0)
            xt = 0.5*((2.*dy(i, 1)+dy(i, 2))*(v(i, 0)+v(i, 1))-dy(i, 1)*&
&             (v(i, -1)+v(i, 2)))/(dy(i, 1)+dy(i, 2))
            bl(i, 1) = xt - v(i, 1)
            br(i, 0) = xt - v(i, 0)
          END DO
          IF (is .EQ. 1) THEN
! out
            bl(1, 0) = 0.
! edge
            br(1, 0) = 0.
! edge
            bl(1, 1) = 0.
! in
            br(1, 1) = 0.
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          IF (ie + 1 .EQ. npx) THEN
! out
            bl(npx, 0) = 0.
! edge
            br(npx, 0) = 0.
! edge
            bl(npx, 1) = 0.
! in
            br(npx, 1) = 0.
            CALL PUSHCONTROL2B(0)
          ELSE
            CALL PUSHCONTROL2B(1)
          END IF
        ELSE
          CALL PUSHCONTROL2B(2)
        END IF
        IF (je + 1 .EQ. npy) THEN
          DO i=is,ie+1
            bl(i, npy-2) = p1*(v(i, npy-3)+v(i, npy-2)) + p2*(v(i, npy-4&
&             )+v(i, npy-1)) - v(i, npy-2)
            xt = c1*v(i, npy-3) + c2*v(i, npy-2) + c3*v(i, npy-1)
            br(i, npy-2) = xt - v(i, npy-2)
            bl(i, npy-1) = xt - v(i, npy-1)
            br(i, npy) = c3*v(i, npy) + c2*v(i, npy+1) + c1*v(i, npy+2) &
&             - v(i, npy)
            xt = 0.5*((2.*dy(i, npy-1)+dy(i, npy-2))*(v(i, npy-1)+v(i, &
&             npy))-dy(i, npy-1)*(v(i, npy-2)+v(i, npy+1)))/(dy(i, npy-1&
&             )+dy(i, npy-2))
            br(i, npy-1) = xt - v(i, npy-1)
            bl(i, npy) = xt - v(i, npy)
          END DO
          IF (is .EQ. 1) THEN
! in
            bl(1, npy-1) = 0.
! edge
            br(1, npy-1) = 0.
! edge
            bl(1, npy) = 0.
! out
            br(1, npy) = 0.
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          IF (ie + 1 .EQ. npx) THEN
! in
            bl(npx, npy-1) = 0.
! edge
            br(npx, npy-1) = 0.
! edge
            bl(npx, npy) = 0.
! out
            br(npx, npy) = 0.
            CALL PUSHCONTROL2B(3)
          ELSE
            CALL PUSHCONTROL2B(2)
          END IF
        ELSE
          CALL PUSHCONTROL2B(1)
        END IF
      ELSE
        CALL PUSHCONTROL2B(0)
      END IF
      DO j=js,je+1
        DO i=is,ie+1
          IF (c(i, j) .GT. 0.) THEN
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
      END DO
      bl_ad = 0.0_8
      br_ad = 0.0_8
      DO j=je+1,js,-1
        DO i=ie+1,is,-1
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            cfl = c(i, j)*rdy(i, j)
            temp0 = bl(i, j) + br(i, j)
            temp_ad12 = (cfl+1.)*flux_ad(i, j)
            v_ad(i, j) = v_ad(i, j) + flux_ad(i, j)
            cfl_ad = temp0*temp_ad12 + (bl(i, j)+cfl*temp0)*flux_ad(i, j&
&             )
            bl_ad(i, j) = bl_ad(i, j) + (cfl+1.0_8)*temp_ad12
            br_ad(i, j) = br_ad(i, j) + cfl*temp_ad12
            flux_ad(i, j) = 0.0_8
            c_ad(i, j) = c_ad(i, j) + rdy(i, j)*cfl_ad
          ELSE
            cfl = c(i, j)*rdy(i, j-1)
            temp = bl(i, j-1) + br(i, j-1)
            temp_ad11 = (1.-cfl)*flux_ad(i, j)
            v_ad(i, j-1) = v_ad(i, j-1) + flux_ad(i, j)
            cfl_ad = -(temp*temp_ad11) - (br(i, j-1)-cfl*temp)*flux_ad(i&
&             , j)
            br_ad(i, j-1) = br_ad(i, j-1) + (1.0_8-cfl)*temp_ad11
            bl_ad(i, j-1) = bl_ad(i, j-1) - cfl*temp_ad11
            flux_ad(i, j) = 0.0_8
            c_ad(i, j) = c_ad(i, j) + rdy(i, j-1)*cfl_ad
          END IF
        END DO
      END DO
      CALL POPCONTROL2B(branch)
      IF (branch .LT. 2) THEN
        IF (branch .EQ. 0) GOTO 100
      ELSE
        IF (branch .NE. 2) THEN
          br_ad(npx, npy) = 0.0_8
          bl_ad(npx, npy) = 0.0_8
          br_ad(npx, npy-1) = 0.0_8
          bl_ad(npx, npy-1) = 0.0_8
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          br_ad(1, npy) = 0.0_8
          bl_ad(1, npy) = 0.0_8
          br_ad(1, npy-1) = 0.0_8
          bl_ad(1, npy-1) = 0.0_8
        END IF
        DO i=ie+1,is,-1
          xt_ad = br_ad(i, npy-1) + bl_ad(i, npy)
          v_ad(i, npy) = v_ad(i, npy) - bl_ad(i, npy)
          bl_ad(i, npy) = 0.0_8
          temp_ad9 = 0.5*xt_ad/(dy(i, npy-1)+dy(i, npy-2))
          temp_ad8 = (dy(i, npy-1)*2.+dy(i, npy-2))*temp_ad9
          v_ad(i, npy-1) = v_ad(i, npy-1) + temp_ad8 - br_ad(i, npy-1)
          br_ad(i, npy-1) = 0.0_8
          temp_ad10 = -(dy(i, npy-1)*temp_ad9)
          v_ad(i, npy) = v_ad(i, npy) + temp_ad8
          v_ad(i, npy-2) = v_ad(i, npy-2) + temp_ad10
          v_ad(i, npy+1) = v_ad(i, npy+1) + temp_ad10
          v_ad(i, npy) = v_ad(i, npy) + (c3-1.0)*br_ad(i, npy)
          v_ad(i, npy+1) = v_ad(i, npy+1) + c2*br_ad(i, npy)
          v_ad(i, npy+2) = v_ad(i, npy+2) + c1*br_ad(i, npy)
          br_ad(i, npy) = 0.0_8
          xt_ad = br_ad(i, npy-2) + bl_ad(i, npy-1)
          v_ad(i, npy-1) = v_ad(i, npy-1) - bl_ad(i, npy-1)
          bl_ad(i, npy-1) = 0.0_8
          v_ad(i, npy-2) = v_ad(i, npy-2) - br_ad(i, npy-2)
          br_ad(i, npy-2) = 0.0_8
          v_ad(i, npy-3) = v_ad(i, npy-3) + c1*xt_ad
          v_ad(i, npy-2) = v_ad(i, npy-2) + c2*xt_ad
          v_ad(i, npy-1) = v_ad(i, npy-1) + c3*xt_ad
          v_ad(i, npy-3) = v_ad(i, npy-3) + p1*bl_ad(i, npy-2)
          v_ad(i, npy-2) = v_ad(i, npy-2) + (p1-1.0)*bl_ad(i, npy-2)
          v_ad(i, npy-4) = v_ad(i, npy-4) + p2*bl_ad(i, npy-2)
          v_ad(i, npy-1) = v_ad(i, npy-1) + p2*bl_ad(i, npy-2)
          bl_ad(i, npy-2) = 0.0_8
        END DO
      END IF
      CALL POPCONTROL2B(branch)
      IF (branch .EQ. 0) THEN
        br_ad(npx, 1) = 0.0_8
        bl_ad(npx, 1) = 0.0_8
        br_ad(npx, 0) = 0.0_8
        bl_ad(npx, 0) = 0.0_8
      ELSE IF (branch .NE. 1) THEN
        GOTO 100
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        br_ad(1, 1) = 0.0_8
        bl_ad(1, 1) = 0.0_8
        br_ad(1, 0) = 0.0_8
        bl_ad(1, 0) = 0.0_8
      END IF
      DO i=ie+1,is,-1
        xt_ad = bl_ad(i, 1) + br_ad(i, 0)
        v_ad(i, 0) = v_ad(i, 0) - br_ad(i, 0)
        br_ad(i, 0) = 0.0_8
        v_ad(i, 1) = v_ad(i, 1) - bl_ad(i, 1)
        bl_ad(i, 1) = 0.0_8
        temp_ad5 = 0.5*xt_ad/(dy(i, 1)+dy(i, 2))
        temp_ad6 = (dy(i, 1)*2.+dy(i, 2))*temp_ad5
        temp_ad7 = -(dy(i, 1)*temp_ad5)
        v_ad(i, 0) = v_ad(i, 0) + temp_ad6
        v_ad(i, 1) = v_ad(i, 1) + temp_ad6
        v_ad(i, -1) = v_ad(i, -1) + temp_ad7
        v_ad(i, 2) = v_ad(i, 2) + temp_ad7
        v_ad(i, -2) = v_ad(i, -2) + c1*bl_ad(i, 0)
        v_ad(i, -1) = v_ad(i, -1) + c2*bl_ad(i, 0)
        v_ad(i, 0) = v_ad(i, 0) + (c3-1.0)*bl_ad(i, 0)
        bl_ad(i, 0) = 0.0_8
        xt_ad = br_ad(i, 1) + bl_ad(i, 2)
        v_ad(i, 2) = v_ad(i, 2) - bl_ad(i, 2)
        bl_ad(i, 2) = 0.0_8
        v_ad(i, 1) = v_ad(i, 1) + c3*xt_ad - br_ad(i, 1)
        br_ad(i, 1) = 0.0_8
        v_ad(i, 2) = v_ad(i, 2) + c2*xt_ad
        v_ad(i, 3) = v_ad(i, 3) + c1*xt_ad
        v_ad(i, 2) = v_ad(i, 2) + (p1-1.0)*br_ad(i, 2)
        v_ad(i, 3) = v_ad(i, 3) + p1*br_ad(i, 2)
        v_ad(i, 1) = v_ad(i, 1) + p2*br_ad(i, 2)
        v_ad(i, 4) = v_ad(i, 4) + p2*br_ad(i, 2)
        br_ad(i, 2) = 0.0_8
      END DO
 100  DO j=min1,max1,-1
        DO i=ie+1,is,-1
          v_ad(i, j-2) = v_ad(i, j-2) + b1*br_ad(i, j)
          v_ad(i, j-1) = v_ad(i, j-1) + b2*br_ad(i, j)
          v_ad(i, j) = v_ad(i, j) + b3*br_ad(i, j)
          v_ad(i, j+1) = v_ad(i, j+1) + b4*br_ad(i, j)
          v_ad(i, j+2) = v_ad(i, j+2) + b5*br_ad(i, j)
          br_ad(i, j) = 0.0_8
          v_ad(i, j-2) = v_ad(i, j-2) + b5*bl_ad(i, j)
          v_ad(i, j-1) = v_ad(i, j-1) + b4*bl_ad(i, j)
          v_ad(i, j) = v_ad(i, j) + b3*bl_ad(i, j)
          v_ad(i, j+1) = v_ad(i, j+1) + b2*bl_ad(i, j)
          v_ad(i, j+2) = v_ad(i, j+2) + b1*bl_ad(i, j)
          bl_ad(i, j) = 0.0_8
        END DO
      END DO
    END SELECT
  END SUBROUTINE YTP_V_ADM
!  SUBROUTINE YTP_V(c, u, v, flux, jord)
!    IMPLICIT NONE
!    INTEGER, INTENT(IN) :: jord
!    REAL, INTENT(IN) :: u(isd:ied, jsd:jed+1)
!    REAL, INTENT(IN) :: v(isd:ied+1, jsd:jed)
!!  Courant   N (like FLUX)
!    REAL(ke_precision), INTENT(IN) :: c(is:ie+1, js:je+1)
!    REAL(ke_precision), INTENT(OUT) :: flux(is:ie+1, js:je+1)
!! Local:
!    REAL(ke_precision) :: dm(is:ie+1, js-2:je+2)
!    REAL(ke_precision) :: al(is:ie+1, js-1:je+2)
!    REAL(ke_precision) :: bl(is:ie+1, js-1:je+1)
!    REAL(ke_precision) :: br(is:ie+1, js-1:je+1)
!    REAL(ke_precision) :: dq(is:ie+1, js-3:je+2)
!    REAL(ke_precision) :: xt, dl, dr, pmp, lac, dqt, cfl
!    REAL(ke_precision) :: x0, x1
!    INTEGER :: i, j
!    INTRINSIC MAX
!    INTRINSIC MIN
!    INTEGER :: min1
!    INTEGER :: max1
!    SELECT CASE  (jord) 
!    CASE (1) 
!      DO j=js,je+1
!        DO i=is,ie+1
!          IF (c(i, j) .GT. 0.) THEN
!            flux(i, j) = v(i, j-1)
!          ELSE
!            flux(i, j) = v(i, j)
!          END IF
!        END DO
!      END DO
!    CASE (333) 
!      DO j=js,je+1
!        DO i=is,ie+1
!!Apply first order at the edges
!          IF (j .EQ. js .OR. j .EQ. je + 1) THEN
!            IF (c(i, j) .GT. 0.) THEN
!              flux(i, j) = v(i, j-1)
!            ELSE
!              flux(i, j) = v(i, j)
!            END IF
!          ELSE IF (c(i, j) .GT. 0.) THEN
!!Otherwise use the third order scheme
!            cfl = c(i, j)*rdy(i, j-1)
!            flux(i, j) = (2.0*v(i, j)+5.0*v(i, j-1)-v(i, j-2))/6.0 - 0.5&
!&             *cfl*(v(i, j)-v(i, j-1)) + cfl*cfl/6.0*(v(i, j)-2.0*v(i, j&
!&             -1)+v(i, j-2))
!          ELSE
!            cfl = c(i, j)*rdy(i, j)
!            flux(i, j) = (2.0*v(i, j-1)+5.0*v(i, j)-v(i, j+1))/6.0 - 0.5&
!&             *cfl*(v(i, j)-v(i, j-1)) + cfl*cfl/6.0*(v(i, j+1)-2.0*v(i&
!&             , j)+v(i, j-1))
!          END IF
!        END DO
!      END DO
!    CASE (6) 
!      IF (3 .LT. js - 1) THEN
!        max1 = js - 1
!      ELSE
!        max1 = 3
!      END IF
!      IF (npy - 3 .GT. je + 1) THEN
!        min1 = je + 1
!      ELSE
!        min1 = npy - 3
!      END IF
!      DO j=max1,min1
!        DO i=is,ie+1
!          bl(i, j) = b5*v(i, j-2) + b4*v(i, j-1) + b3*v(i, j) + b2*v(i, &
!&           j+1) + b1*v(i, j+2)
!          br(i, j) = b1*v(i, j-2) + b2*v(i, j-1) + b3*v(i, j) + b4*v(i, &
!&           j+1) + b5*v(i, j+2)
!        END DO
!      END DO
!      IF (grid_type .LT. 3) THEN
!        IF (js .EQ. 1) THEN
!          DO i=is,ie+1
!            br(i, 2) = p1*(v(i, 2)+v(i, 3)) + p2*(v(i, 1)+v(i, 4)) - v(i&
!&             , 2)
!            xt = c3*v(i, 1) + c2*v(i, 2) + c1*v(i, 3)
!            br(i, 1) = xt - v(i, 1)
!            bl(i, 2) = xt - v(i, 2)
!            bl(i, 0) = c1*v(i, -2) + c2*v(i, -1) + c3*v(i, 0) - v(i, 0)
!            xt = 0.5*((2.*dy(i, 1)+dy(i, 2))*(v(i, 0)+v(i, 1))-dy(i, 1)*&
!&             (v(i, -1)+v(i, 2)))/(dy(i, 1)+dy(i, 2))
!            bl(i, 1) = xt - v(i, 1)
!            br(i, 0) = xt - v(i, 0)
!          END DO
!          IF (is .EQ. 1) THEN
!! out
!            bl(1, 0) = 0.
!! edge
!            br(1, 0) = 0.
!! edge
!            bl(1, 1) = 0.
!! in
!            br(1, 1) = 0.
!          END IF
!          IF (ie + 1 .EQ. npx) THEN
!! out
!            bl(npx, 0) = 0.
!! edge
!            br(npx, 0) = 0.
!! edge
!            bl(npx, 1) = 0.
!! in
!            br(npx, 1) = 0.
!          END IF
!        END IF
!        IF (je + 1 .EQ. npy) THEN
!          DO i=is,ie+1
!            bl(i, npy-2) = p1*(v(i, npy-3)+v(i, npy-2)) + p2*(v(i, npy-4&
!&             )+v(i, npy-1)) - v(i, npy-2)
!            xt = c1*v(i, npy-3) + c2*v(i, npy-2) + c3*v(i, npy-1)
!            br(i, npy-2) = xt - v(i, npy-2)
!            bl(i, npy-1) = xt - v(i, npy-1)
!            br(i, npy) = c3*v(i, npy) + c2*v(i, npy+1) + c1*v(i, npy+2) &
!&             - v(i, npy)
!            xt = 0.5*((2.*dy(i, npy-1)+dy(i, npy-2))*(v(i, npy-1)+v(i, &
!&             npy))-dy(i, npy-1)*(v(i, npy-2)+v(i, npy+1)))/(dy(i, npy-1&
!&             )+dy(i, npy-2))
!            br(i, npy-1) = xt - v(i, npy-1)
!            bl(i, npy) = xt - v(i, npy)
!          END DO
!          IF (is .EQ. 1) THEN
!! in
!            bl(1, npy-1) = 0.
!! edge
!            br(1, npy-1) = 0.
!! edge
!            bl(1, npy) = 0.
!! out
!            br(1, npy) = 0.
!          END IF
!          IF (ie + 1 .EQ. npx) THEN
!! in
!            bl(npx, npy-1) = 0.
!! edge
!            br(npx, npy-1) = 0.
!! edge
!            bl(npx, npy) = 0.
!! out
!            br(npx, npy) = 0.
!          END IF
!        END IF
!      END IF
!      DO j=js,je+1
!        DO i=is,ie+1
!          IF (c(i, j) .GT. 0.) THEN
!            cfl = c(i, j)*rdy(i, j-1)
!            flux(i, j) = v(i, j-1) + (1.-cfl)*(br(i, j-1)-cfl*(bl(i, j-1&
!&             )+br(i, j-1)))
!          ELSE
!            cfl = c(i, j)*rdy(i, j)
!            flux(i, j) = v(i, j) + (1.+cfl)*(bl(i, j)+cfl*(bl(i, j)+br(i&
!&             , j)))
!          END IF
!        END DO
!      END DO
!    END SELECT
!  END SUBROUTINE YTP_V
!  Differentiation of d2a2c_vect in reverse (adjoint) mode:
!   gradient     of useful results: u v ua uc ut va vc vt
!   with respect to varying inputs: u v ua uc ut va vc vt
  SUBROUTINE D2A2C_VECT_ADM(u, u_ad, v, v_ad, ua, ua_ad, va, va_ad, uc, &
&   uc_ad, vc, vc_ad, ut, ut_ad, vt, vt_ad, dord4)
    IMPLICIT NONE
    REAL, INTENT(IN) :: u(isd:ied, jsd:jed+1)
    REAL :: u_ad(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: v(isd:ied+1, jsd:jed)
    REAL :: v_ad(isd:ied+1, jsd:jed)
    LOGICAL, INTENT(IN) :: dord4
    REAL, DIMENSION(isd:ied + 1, jsd:jed) :: uc
    REAL, DIMENSION(isd:ied+1, jsd:jed) :: uc_ad
    REAL, DIMENSION(isd:ied, jsd:jed + 1) :: vc
    REAL, DIMENSION(isd:ied, jsd:jed+1) :: vc_ad
    REAL, DIMENSION(isd:ied, jsd:jed) :: ua, va, ut, vt
    REAL, DIMENSION(isd:ied, jsd:jed) :: ua_ad, va_ad, ut_ad, vt_ad
! Local 
    REAL, DIMENSION(isd:ied, jsd:jed) :: utmp, vtmp
    REAL, DIMENSION(isd:ied, jsd:jed) :: utmp_ad, vtmp_ad
    INTEGER :: npt, i, j, ifirst, ilast, id
    INTRINSIC MAX
    INTRINSIC MIN
    INTEGER :: ad_from
    INTEGER :: ad_to
    INTEGER :: ad_from0
    INTEGER :: ad_to0
    INTEGER :: branch
    INTEGER :: min6
    INTEGER :: min5
    INTEGER :: min4
    INTEGER :: min3
    INTEGER :: min2
    INTEGER :: min1
    REAL :: temp_ad9
    REAL :: temp_ad8
    REAL :: temp_ad7
    REAL :: temp_ad6
    REAL :: temp_ad5
    REAL :: temp_ad4
    REAL :: temp_ad3
    REAL :: temp_ad2
    REAL :: temp_ad1
    REAL :: temp_ad0
    REAL :: temp_ad
    REAL :: temp_ad12
    REAL :: temp_ad11
    REAL :: temp_ad10
    INTEGER :: max6
    INTEGER :: max5
    INTEGER :: max4
    INTEGER :: max3
    INTEGER :: max2
    INTEGER :: max1
    IF (dord4) THEN
      id = 1
    ELSE
      id = 0
    END IF
    IF (grid_type .LT. 3) THEN
      npt = 4
    ELSE
      npt = -2
    END IF
! DSK This routine does not properly fill the halo region of utmp
    IF (npt .LT. js - 1) THEN
      max1 = js - 1
    ELSE
      max1 = npt
    END IF
    IF (npy - npt .GT. je + 1) THEN
      min1 = je + 1
    ELSE
      min1 = npy - npt
    END IF
!----------
! Interior:
!----------
    DO j=max1,min1
      IF (npt .LT. isd) THEN
        max2 = isd
      ELSE
        max2 = npt
      END IF
      IF (npx - npt .GT. ied) THEN
        min2 = ied
      ELSE
        min2 = npx - npt
      END IF
      ad_from = max2
      i = min2 + 1
      CALL PUSHINTEGER4(i - 1)
      CALL PUSHINTEGER4(ad_from)
    END DO
    IF (npt .LT. jsd) THEN
      max3 = jsd
    ELSE
      max3 = npt
    END IF
    IF (npy - npt .GT. jed) THEN
      min3 = jed
    ELSE
      min3 = npy - npt
    END IF
    DO j=max3,min3
      IF (npt .LT. is - 1) THEN
        max4 = is - 1
      ELSE
        max4 = npt
      END IF
      IF (npx - npt .GT. ie + 1) THEN
        min4 = ie + 1
      ELSE
        min4 = npx - npt
      END IF
      ad_from0 = max4
      i = min4 + 1
      CALL PUSHINTEGER4(i - 1)
      CALL PUSHINTEGER4(ad_from0)
    END DO
!----------
! edges:
!----------
    IF (grid_type .LT. 3) THEN
      IF (js .EQ. 1 .OR. jsd .LT. npt) THEN
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
!#endif
      IF (je + 1 .EQ. npy .OR. jed .GE. npy - npt) THEN
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
!#endif
      IF (is .EQ. 1 .OR. isd .LT. npt) THEN
        IF (npt .LT. jsd) THEN
          max5 = jsd
        ELSE
          max5 = npt
        END IF
        IF (npy - npt .GT. jed) THEN
          min5 = jed
        ELSE
          min5 = npy - npt
        END IF
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
!#endif
      IF (ie + 1 .EQ. npx .OR. ied .GE. npx - npt) THEN
        IF (npt .LT. jsd) THEN
          max6 = jsd
        ELSE
          max6 = npt
        END IF
        IF (npy - npt .GT. jed) THEN
          min6 = jed
        ELSE
          min6 = npy - npt
        END IF
        CALL PUSHCONTROL2B(2)
      ELSE
        CALL PUSHCONTROL2B(1)
      END IF
    ELSE
      CALL PUSHCONTROL2B(0)
    END IF
! A -> C
!--------------
! Fix the edges
!--------------
! Xdir:
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
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    IF (grid_type .LT. 3) THEN
      IF (3 .LT. is - 1) THEN
        ifirst = is - 1
      ELSE
        ifirst = 3
      END IF
      IF (npx - 2 .GT. ie + 2) THEN
        CALL PUSHCONTROL1B(1)
        ilast = ie + 2
      ELSE
        CALL PUSHCONTROL1B(1)
        ilast = npx - 2
      END IF
    ELSE
      CALL PUSHCONTROL1B(0)
      ifirst = is - 1
      ilast = ie + 2
    END IF
    IF (grid_type .LT. 3) THEN
      IF (is .EQ. 1) THEN
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (ie + 1 .EQ. npx) THEN
        CALL PUSHCONTROL2B(0)
      ELSE
        CALL PUSHCONTROL2B(1)
      END IF
    ELSE
      CALL PUSHCONTROL2B(2)
    END IF
!------
! Ydir:
!------
    IF (sw_corner) THEN
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    IF (nw_corner) THEN
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
    IF (grid_type .LT. 3) THEN
      DO j=js-1,je+2
        IF (j .EQ. 1) THEN
          CALL PUSHCONTROL3B(4)
        ELSE IF (j .EQ. 0 .OR. j .EQ. npy - 1) THEN
          CALL PUSHCONTROL3B(3)
        ELSE IF (j .EQ. 2 .OR. j .EQ. npy + 1) THEN
          CALL PUSHCONTROL3B(2)
        ELSE IF (j .EQ. npy) THEN
          CALL PUSHCONTROL3B(1)
        ELSE
          CALL PUSHCONTROL3B(0)
        END IF
      END DO
      vtmp_ad = 0.0_8
      DO j=je+2,js-1,-1
        CALL POPCONTROL3B(branch)
        IF (branch .LT. 2) THEN
          IF (branch .EQ. 0) THEN
            DO i=ie+1,is-1,-1
              temp_ad12 = rsin_v(i, j)*vt_ad(i, j)
              vc_ad(i, j) = vc_ad(i, j) + temp_ad12
              u_ad(i, j) = u_ad(i, j) - cosa_v(i, j)*temp_ad12
              vt_ad(i, j) = 0.0_8
              vtmp_ad(i, j-2) = vtmp_ad(i, j-2) + a2*vc_ad(i, j)
              vtmp_ad(i, j+1) = vtmp_ad(i, j+1) + a2*vc_ad(i, j)
              vtmp_ad(i, j-1) = vtmp_ad(i, j-1) + a1*vc_ad(i, j)
              vtmp_ad(i, j) = vtmp_ad(i, j) + a1*vc_ad(i, j)
              vc_ad(i, j) = 0.0_8
            END DO
          ELSE
            DO i=ie+1,is-1,-1
              vc_ad(i, npy) = vc_ad(i, npy) + rsin_v(i, npy)*vt_ad(i, &
&               npy)
              vt_ad(i, npy) = 0.0_8
              temp_ad11 = rsin_v(i, npy)*vc_ad(i, npy)
              vtmp_ad(i, npy-1) = vtmp_ad(i, npy-1) + t14*temp_ad11
              vtmp_ad(i, npy) = vtmp_ad(i, npy) + t14*temp_ad11
              vtmp_ad(i, npy-2) = vtmp_ad(i, npy-2) + t12*temp_ad11
              vtmp_ad(i, npy+1) = vtmp_ad(i, npy+1) + t12*temp_ad11
              vtmp_ad(i, npy-3) = vtmp_ad(i, npy-3) + t15*temp_ad11
              vtmp_ad(i, npy+2) = vtmp_ad(i, npy+2) + t15*temp_ad11
              vc_ad(i, npy) = 0.0_8
            END DO
          END IF
        ELSE IF (branch .EQ. 2) THEN
          DO i=ie+1,is-1,-1
            temp_ad10 = rsin_v(i, j)*vt_ad(i, j)
            vc_ad(i, j) = vc_ad(i, j) + temp_ad10
            u_ad(i, j) = u_ad(i, j) - cosa_v(i, j)*temp_ad10
            vt_ad(i, j) = 0.0_8
            vtmp_ad(i, j+1) = vtmp_ad(i, j+1) + c1*vc_ad(i, j)
            vtmp_ad(i, j) = vtmp_ad(i, j) + c2*vc_ad(i, j)
            vtmp_ad(i, j-1) = vtmp_ad(i, j-1) + c3*vc_ad(i, j)
            vc_ad(i, j) = 0.0_8
          END DO
        ELSE IF (branch .EQ. 3) THEN
          DO i=ie+1,is-1,-1
            temp_ad9 = rsin_v(i, j)*vt_ad(i, j)
            vc_ad(i, j) = vc_ad(i, j) + temp_ad9
            u_ad(i, j) = u_ad(i, j) - cosa_v(i, j)*temp_ad9
            vt_ad(i, j) = 0.0_8
            vtmp_ad(i, j-2) = vtmp_ad(i, j-2) + c1*vc_ad(i, j)
            vtmp_ad(i, j-1) = vtmp_ad(i, j-1) + c2*vc_ad(i, j)
            vtmp_ad(i, j) = vtmp_ad(i, j) + c3*vc_ad(i, j)
            vc_ad(i, j) = 0.0_8
          END DO
        ELSE
          DO i=ie+1,is-1,-1
            vc_ad(i, 1) = vc_ad(i, 1) + rsin_v(i, 1)*vt_ad(i, 1)
            vt_ad(i, 1) = 0.0_8
            temp_ad8 = rsin_v(i, 1)*vc_ad(i, 1)
            vtmp_ad(i, 0) = vtmp_ad(i, 0) + t14*temp_ad8
            vtmp_ad(i, 1) = vtmp_ad(i, 1) + t14*temp_ad8
            vtmp_ad(i, -1) = vtmp_ad(i, -1) + t12*temp_ad8
            vtmp_ad(i, 2) = vtmp_ad(i, 2) + t12*temp_ad8
            vtmp_ad(i, -2) = vtmp_ad(i, -2) + t15*temp_ad8
            vtmp_ad(i, 3) = vtmp_ad(i, 3) + t15*temp_ad8
            vc_ad(i, 1) = 0.0_8
          END DO
        END IF
      END DO
    ELSE
      vtmp_ad = 0.0_8
      DO j=je+2,js-1,-1
        DO i=ie+1,is-1,-1
          vc_ad(i, j) = vc_ad(i, j) + vt_ad(i, j)
          vt_ad(i, j) = 0.0_8
          vtmp_ad(i, j-2) = vtmp_ad(i, j-2) + a2*vc_ad(i, j)
          vtmp_ad(i, j+1) = vtmp_ad(i, j+1) + a2*vc_ad(i, j)
          vtmp_ad(i, j-1) = vtmp_ad(i, j-1) + a1*vc_ad(i, j)
          vtmp_ad(i, j) = vtmp_ad(i, j) + a1*vc_ad(i, j)
          vc_ad(i, j) = 0.0_8
        END DO
      END DO
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      utmp_ad = 0.0_8
      DO j=2,0,-1
        utmp_ad(ie-j, npy) = utmp_ad(ie-j, npy) - vtmp_ad(npx, npy+j)
        vtmp_ad(npx, npy+j) = 0.0_8
      END DO
    ELSE
      utmp_ad = 0.0_8
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      DO j=0,-2,-1
        utmp_ad(ie+j, 0) = utmp_ad(ie+j, 0) + vtmp_ad(npx, j)
        vtmp_ad(npx, j) = 0.0_8
      END DO
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      DO j=2,0,-1
        utmp_ad(j+1, npy) = utmp_ad(j+1, npy) + vtmp_ad(0, npy+j)
        vtmp_ad(0, npy+j) = 0.0_8
      END DO
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      DO j=0,-2,-1
        utmp_ad(1-j, 0) = utmp_ad(1-j, 0) - vtmp_ad(0, j)
        vtmp_ad(0, j) = 0.0_8
      END DO
    END IF
    CALL POPCONTROL2B(branch)
    IF (branch .EQ. 0) THEN
      DO j=je+1,js-1,-1
        temp_ad5 = rsin_u(npx+1, j)*ut_ad(npx+1, j)
        uc_ad(npx+1, j) = uc_ad(npx+1, j) + temp_ad5
        v_ad(npx+1, j) = v_ad(npx+1, j) - cosa_u(npx+1, j)*temp_ad5
        ut_ad(npx+1, j) = 0.0_8
        uc_ad(npx, j) = uc_ad(npx, j) + rsin_u(npx, j)*ut_ad(npx, j)
        ut_ad(npx, j) = 0.0_8
        temp_ad6 = rsin_u(npx-1, j)*ut_ad(npx-1, j)
        uc_ad(npx-1, j) = uc_ad(npx-1, j) + temp_ad6
        v_ad(npx-1, j) = v_ad(npx-1, j) - cosa_u(npx-1, j)*temp_ad6
        ut_ad(npx-1, j) = 0.0_8
        utmp_ad(npx, j) = utmp_ad(npx, j) + c3*uc_ad(npx+1, j)
        utmp_ad(npx+1, j) = utmp_ad(npx+1, j) + c2*uc_ad(npx+1, j)
        utmp_ad(npx+2, j) = utmp_ad(npx+2, j) + c1*uc_ad(npx+1, j)
        uc_ad(npx+1, j) = 0.0_8
        temp_ad7 = rsin_u(npx, j)*uc_ad(npx, j)
        utmp_ad(npx-1, j) = utmp_ad(npx-1, j) + t14*temp_ad7
        utmp_ad(npx, j) = utmp_ad(npx, j) + t14*temp_ad7
        utmp_ad(npx-2, j) = utmp_ad(npx-2, j) + t12*temp_ad7
        utmp_ad(npx+1, j) = utmp_ad(npx+1, j) + t12*temp_ad7
        utmp_ad(npx-3, j) = utmp_ad(npx-3, j) + t15*temp_ad7
        utmp_ad(npx+2, j) = utmp_ad(npx+2, j) + t15*temp_ad7
        uc_ad(npx, j) = 0.0_8
        utmp_ad(npx-3, j) = utmp_ad(npx-3, j) + c1*uc_ad(npx-1, j)
        utmp_ad(npx-2, j) = utmp_ad(npx-2, j) + c2*uc_ad(npx-1, j)
        utmp_ad(npx-1, j) = utmp_ad(npx-1, j) + c3*uc_ad(npx-1, j)
        uc_ad(npx-1, j) = 0.0_8
      END DO
    ELSE IF (branch .NE. 1) THEN
      GOTO 100
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      DO j=je+1,js-1,-1
        temp_ad2 = rsin_u(2, j)*ut_ad(2, j)
        uc_ad(2, j) = uc_ad(2, j) + temp_ad2
        v_ad(2, j) = v_ad(2, j) - cosa_u(2, j)*temp_ad2
        ut_ad(2, j) = 0.0_8
        uc_ad(1, j) = uc_ad(1, j) + rsin_u(1, j)*ut_ad(1, j)
        ut_ad(1, j) = 0.0_8
        temp_ad3 = rsin_u(0, j)*ut_ad(0, j)
        uc_ad(0, j) = uc_ad(0, j) + temp_ad3
        v_ad(0, j) = v_ad(0, j) - cosa_u(0, j)*temp_ad3
        ut_ad(0, j) = 0.0_8
        utmp_ad(3, j) = utmp_ad(3, j) + c1*uc_ad(2, j)
        utmp_ad(2, j) = utmp_ad(2, j) + c2*uc_ad(2, j)
        utmp_ad(1, j) = utmp_ad(1, j) + c3*uc_ad(2, j)
        uc_ad(2, j) = 0.0_8
        temp_ad4 = rsin_u(1, j)*uc_ad(1, j)
        utmp_ad(0, j) = utmp_ad(0, j) + t14*temp_ad4
        utmp_ad(1, j) = utmp_ad(1, j) + t14*temp_ad4
        utmp_ad(-1, j) = utmp_ad(-1, j) + t12*temp_ad4
        utmp_ad(2, j) = utmp_ad(2, j) + t12*temp_ad4
        utmp_ad(-2, j) = utmp_ad(-2, j) + t15*temp_ad4
        utmp_ad(3, j) = utmp_ad(3, j) + t15*temp_ad4
        uc_ad(1, j) = 0.0_8
        utmp_ad(-2, j) = utmp_ad(-2, j) + c1*uc_ad(0, j)
        utmp_ad(-1, j) = utmp_ad(-1, j) + c2*uc_ad(0, j)
        utmp_ad(0, j) = utmp_ad(0, j) + c3*uc_ad(0, j)
        uc_ad(0, j) = 0.0_8
      END DO
    END IF
 100 DO j=je+1,js-1,-1
      DO i=ilast,ifirst,-1
        temp_ad1 = rsin_u(i, j)*ut_ad(i, j)
        uc_ad(i, j) = uc_ad(i, j) + temp_ad1
        v_ad(i, j) = v_ad(i, j) - cosa_u(i, j)*temp_ad1
        ut_ad(i, j) = 0.0_8
        utmp_ad(i-1, j) = utmp_ad(i-1, j) + a1*uc_ad(i, j)
        utmp_ad(i, j) = utmp_ad(i, j) + a1*uc_ad(i, j)
        utmp_ad(i-2, j) = utmp_ad(i-2, j) + a2*uc_ad(i, j)
        utmp_ad(i+1, j) = utmp_ad(i+1, j) + a2*uc_ad(i, j)
        uc_ad(i, j) = 0.0_8
      END DO
    END DO
    CALL POPCONTROL1B(branch)
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      DO i=0,-2,-1
        vtmp_ad(0, je+i) = vtmp_ad(0, je+i) + utmp_ad(i, npy)
        utmp_ad(i, npy) = 0.0_8
      END DO
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      DO i=2,0,-1
        vtmp_ad(npx, je-i) = vtmp_ad(npx, je-i) - utmp_ad(npx+i, npy)
        utmp_ad(npx+i, npy) = 0.0_8
      END DO
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      DO i=2,0,-1
        vtmp_ad(npx, i+1) = vtmp_ad(npx, i+1) + utmp_ad(npx+i, 0)
        utmp_ad(npx+i, 0) = 0.0_8
      END DO
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      DO i=0,-2,-1
        vtmp_ad(0, 1-i) = vtmp_ad(0, 1-i) - utmp_ad(i, 0)
        utmp_ad(i, 0) = 0.0_8
      END DO
    END IF
    DO j=je+id+1,js-1-id,-1
      DO i=ie+id+1,is-1-id,-1
        temp_ad0 = rsin2(i, j)*ua_ad(i, j)
        temp_ad = rsin2(i, j)*va_ad(i, j)
        vtmp_ad(i, j) = vtmp_ad(i, j) + temp_ad - cosa_s(i, j)*temp_ad0
        utmp_ad(i, j) = utmp_ad(i, j) + temp_ad0 - cosa_s(i, j)*temp_ad
        va_ad(i, j) = 0.0_8
        ua_ad(i, j) = 0.0_8
      END DO
    END DO
    CALL POPCONTROL2B(branch)
    IF (branch .NE. 0) THEN
      IF (branch .NE. 1) THEN
        DO j=min6,max6,-1
          DO i=ied,npx-npt+1,-1
            v_ad(i, j) = v_ad(i, j) + 0.5*vtmp_ad(i, j)
            v_ad(i+1, j) = v_ad(i+1, j) + 0.5*vtmp_ad(i, j)
            vtmp_ad(i, j) = 0.0_8
            u_ad(i, j) = u_ad(i, j) + 0.5*utmp_ad(i, j)
            u_ad(i, j+1) = u_ad(i, j+1) + 0.5*utmp_ad(i, j)
            utmp_ad(i, j) = 0.0_8
          END DO
        END DO
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO j=min5,max5,-1
          DO i=npt-1,isd,-1
            v_ad(i, j) = v_ad(i, j) + 0.5*vtmp_ad(i, j)
            v_ad(i+1, j) = v_ad(i+1, j) + 0.5*vtmp_ad(i, j)
            vtmp_ad(i, j) = 0.0_8
            u_ad(i, j) = u_ad(i, j) + 0.5*utmp_ad(i, j)
            u_ad(i, j+1) = u_ad(i, j+1) + 0.5*utmp_ad(i, j)
            utmp_ad(i, j) = 0.0_8
          END DO
        END DO
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO j=jed,npy-npt+1,-1
          DO i=ied,isd,-1
            v_ad(i, j) = v_ad(i, j) + 0.5*vtmp_ad(i, j)
            v_ad(i+1, j) = v_ad(i+1, j) + 0.5*vtmp_ad(i, j)
            vtmp_ad(i, j) = 0.0_8
            u_ad(i, j) = u_ad(i, j) + 0.5*utmp_ad(i, j)
            u_ad(i, j+1) = u_ad(i, j+1) + 0.5*utmp_ad(i, j)
            utmp_ad(i, j) = 0.0_8
          END DO
        END DO
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO j=npt-1,jsd,-1
          DO i=ied,isd,-1
            v_ad(i, j) = v_ad(i, j) + 0.5*vtmp_ad(i, j)
            v_ad(i+1, j) = v_ad(i+1, j) + 0.5*vtmp_ad(i, j)
            vtmp_ad(i, j) = 0.0_8
            u_ad(i, j) = u_ad(i, j) + 0.5*utmp_ad(i, j)
            u_ad(i, j+1) = u_ad(i, j+1) + 0.5*utmp_ad(i, j)
            utmp_ad(i, j) = 0.0_8
          END DO
        END DO
      END IF
    END IF
    DO j=min3,max3,-1
      CALL POPINTEGER4(ad_from0)
      CALL POPINTEGER4(ad_to0)
      DO i=ad_to0,ad_from0,-1
        v_ad(i-1, j) = v_ad(i-1, j) + a2*vtmp_ad(i, j)
        v_ad(i+2, j) = v_ad(i+2, j) + a2*vtmp_ad(i, j)
        v_ad(i, j) = v_ad(i, j) + a1*vtmp_ad(i, j)
        v_ad(i+1, j) = v_ad(i+1, j) + a1*vtmp_ad(i, j)
        vtmp_ad(i, j) = 0.0_8
      END DO
    END DO
    DO j=min1,max1,-1
      CALL POPINTEGER4(ad_from)
      CALL POPINTEGER4(ad_to)
      DO i=ad_to,ad_from,-1
        u_ad(i, j-1) = u_ad(i, j-1) + a2*utmp_ad(i, j)
        u_ad(i, j+2) = u_ad(i, j+2) + a2*utmp_ad(i, j)
        u_ad(i, j) = u_ad(i, j) + a1*utmp_ad(i, j)
        u_ad(i, j+1) = u_ad(i, j+1) + a1*utmp_ad(i, j)
        utmp_ad(i, j) = 0.0_8
      END DO
    END DO
  END SUBROUTINE D2A2C_VECT_ADM
  SUBROUTINE D2A2C_VECT(u, v, ua, va, uc, vc, ut, vt, dord4)
    IMPLICIT NONE
    REAL, INTENT(IN) :: u(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: v(isd:ied+1, jsd:jed)
    LOGICAL, INTENT(IN) :: dord4
    REAL, DIMENSION(isd:ied + 1, jsd:jed), INTENT(OUT) :: uc
    REAL, DIMENSION(isd:ied, jsd:jed + 1), INTENT(OUT) :: vc
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(OUT) :: ua, va, ut, vt
! Local 
    REAL, DIMENSION(isd:ied, jsd:jed) :: utmp, vtmp
    INTEGER :: npt, i, j, ifirst, ilast, id
    INTRINSIC MAX
    INTRINSIC MIN
    INTEGER :: min6
    INTEGER :: min5
    INTEGER :: min4
    INTEGER :: min3
    INTEGER :: min2
    INTEGER :: min1
    INTEGER :: max6
    INTEGER :: max5
    INTEGER :: max4
    INTEGER :: max3
    INTEGER :: max2
    INTEGER :: max1
    IF (dord4) THEN
      id = 1
    ELSE
      id = 0
    END IF
    IF (grid_type .LT. 3) THEN
      npt = 4
    ELSE
      npt = -2
    END IF
! DSK This routine does not properly fill the halo region of utmp
    utmp = 0.0
    vtmp = 0.0
    IF (npt .LT. js - 1) THEN
      max1 = js - 1
    ELSE
      max1 = npt
    END IF
    IF (npy - npt .GT. je + 1) THEN
      min1 = je + 1
    ELSE
      min1 = npy - npt
    END IF
!----------
! Interior:
!----------
    DO j=max1,min1
      IF (npt .LT. isd) THEN
        max2 = isd
      ELSE
        max2 = npt
      END IF
      IF (npx - npt .GT. ied) THEN
        min2 = ied
      ELSE
        min2 = npx - npt
      END IF
      DO i=max2,min2
        utmp(i, j) = a2*(u(i, j-1)+u(i, j+2)) + a1*(u(i, j)+u(i, j+1))
      END DO
    END DO
    IF (npt .LT. jsd) THEN
      max3 = jsd
    ELSE
      max3 = npt
    END IF
    IF (npy - npt .GT. jed) THEN
      min3 = jed
    ELSE
      min3 = npy - npt
    END IF
    DO j=max3,min3
      IF (npt .LT. is - 1) THEN
        max4 = is - 1
      ELSE
        max4 = npt
      END IF
      IF (npx - npt .GT. ie + 1) THEN
        min4 = ie + 1
      ELSE
        min4 = npx - npt
      END IF
      DO i=max4,min4
        vtmp(i, j) = a2*(v(i-1, j)+v(i+2, j)) + a1*(v(i, j)+v(i+1, j))
      END DO
    END DO
!----------
! edges:
!----------
    IF (grid_type .LT. 3) THEN
      IF (js .EQ. 1 .OR. jsd .LT. npt) THEN
!#ifdef CONSV_VT
!      do j=jsd,npt-1
!         do i=isd,ied+1
!            uc(i,j) = v(i,j)*dy(i,j)
!         enddo
!      enddo
!      do j=jsd,npt
!         do i=isd,ied
!            vc(i,j) = u(i,j)*dx(i,j)
!         enddo
!      enddo
!#endif
        DO j=jsd,npt-1
          DO i=isd,ied
!#ifdef CONSV_VT
!            utmp(i,j) = 0.5*(vc(i,j) + vc(i,j+1)) * rdxa(i,j)
!            vtmp(i,j) = 0.5*(uc(i,j) + uc(i+1,j)) * rdya(i,j)
!#else
            utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
            vtmp(i, j) = 0.5*(v(i, j)+v(i+1, j))
          END DO
        END DO
      END IF
!#endif
      IF (je + 1 .EQ. npy .OR. jed .GE. npy - npt) THEN
!#ifdef CONSV_VT
!      do j=npy-npt+1,jed
!         do i=isd,ied+1
!            uc(i,j) = v(i,j)*dy(i,j)
!         enddo
!      enddo
!      do j=npy-npt+1,jed+1
!         do i=isd,ied
!            vc(i,j) = u(i,j)*dx(i,j)
!         enddo
!      enddo
!#endif
        DO j=npy-npt+1,jed
          DO i=isd,ied
!#ifdef CONSV_VT
!            utmp(i,j) = 0.5*(vc(i,j) + vc(i,j+1)) * rdxa(i,j)
!            vtmp(i,j) = 0.5*(uc(i,j) + uc(i+1,j)) * rdya(i,j)
!#else
            utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
            vtmp(i, j) = 0.5*(v(i, j)+v(i+1, j))
          END DO
        END DO
      END IF
!#endif
      IF (is .EQ. 1 .OR. isd .LT. npt) THEN
        IF (npt .LT. jsd) THEN
          max5 = jsd
        ELSE
          max5 = npt
        END IF
        IF (npy - npt .GT. jed) THEN
          min5 = jed
        ELSE
          min5 = npy - npt
        END IF
!#ifdef CONSV_VT
!      do j=max(npt,jsd),min(npy-npt,jed)
!         do i=isd,npt
!            uc(i,j) = v(i,j)*dy(i,j)
!         enddo
!      enddo
!      do j=max(npt,jsd),min(npy-npt+1,jed+1)
!         do i=isd,npt-1
!            vc(i,j) = u(i,j)*dx(i,j)
!         enddo
!      enddo
!#endif
        DO j=max5,min5
          DO i=isd,npt-1
!#ifdef CONSV_VT
!            utmp(i,j) = 0.5*(vc(i,j) + vc(i,j+1)) * rdxa(i,j)
!            vtmp(i,j) = 0.5*(uc(i,j) + uc(i+1,j)) * rdya(i,j)
!#else
            utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
            vtmp(i, j) = 0.5*(v(i, j)+v(i+1, j))
          END DO
        END DO
      END IF
!#endif
      IF (ie + 1 .EQ. npx .OR. ied .GE. npx - npt) THEN
        IF (npt .LT. jsd) THEN
          max6 = jsd
        ELSE
          max6 = npt
        END IF
        IF (npy - npt .GT. jed) THEN
          min6 = jed
        ELSE
          min6 = npy - npt
        END IF
!#ifdef CONSV_VT
!      do j=max(npt,jsd),min(npy-npt,jed)
!         do i=npx-npt+1,ied+1
!            uc(i,j) = v(i,j)*dy(i,j)
!         enddo
!      enddo
!      do j=max(npt,jsd),min(npy-npt+1,jed+1)
!         do i=npx-npt+1,ied
!            vc(i,j) = u(i,j)*dx(i,j)
!         enddo
!      enddo
!#endif
        DO j=max6,min6
          DO i=npx-npt+1,ied
!#ifdef CONSV_VT
!            utmp(i,j) = 0.5*(vc(i,j) + vc(i,j+1)) * rdxa(i,j)
!            vtmp(i,j) = 0.5*(uc(i,j) + uc(i+1,j)) * rdya(i,j)
!#else
            utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
            vtmp(i, j) = 0.5*(v(i, j)+v(i+1, j))
          END DO
        END DO
      END IF
    END IF
!#endif
    DO j=js-1-id,je+1+id
      DO i=is-1-id,ie+1+id
        ua(i, j) = (utmp(i, j)-vtmp(i, j)*cosa_s(i, j))*rsin2(i, j)
        va(i, j) = (vtmp(i, j)-utmp(i, j)*cosa_s(i, j))*rsin2(i, j)
      END DO
    END DO
! A -> C
!--------------
! Fix the edges
!--------------
! Xdir:
    IF (sw_corner) THEN
      DO i=-2,0
        utmp(i, 0) = -vtmp(0, 1-i)
      END DO
    END IF
    IF (se_corner) THEN
      DO i=0,2
        utmp(npx+i, 0) = vtmp(npx, i+1)
      END DO
    END IF
    IF (ne_corner) THEN
      DO i=0,2
        utmp(npx+i, npy) = -vtmp(npx, je-i)
      END DO
    END IF
    IF (nw_corner) THEN
      DO i=-2,0
        utmp(i, npy) = vtmp(0, je+i)
      END DO
    END IF
    IF (grid_type .LT. 3) THEN
      IF (3 .LT. is - 1) THEN
        ifirst = is - 1
      ELSE
        ifirst = 3
      END IF
      IF (npx - 2 .GT. ie + 2) THEN
        ilast = ie + 2
      ELSE
        ilast = npx - 2
      END IF
    ELSE
      ifirst = is - 1
      ilast = ie + 2
    END IF
!---------------------------------------------
! 4th order interpolation for interior points:
!---------------------------------------------
    DO j=js-1,je+1
      DO i=ifirst,ilast
        uc(i, j) = a1*(utmp(i-1, j)+utmp(i, j)) + a2*(utmp(i-2, j)+utmp(&
&         i+1, j))
        ut(i, j) = (uc(i, j)-v(i, j)*cosa_u(i, j))*rsin_u(i, j)
      END DO
    END DO
    IF (grid_type .LT. 3) THEN
      IF (is .EQ. 1) THEN
        DO j=js-1,je+1
          uc(0, j) = c1*utmp(-2, j) + c2*utmp(-1, j) + c3*utmp(0, j)
! 3-pt extrapolation --------------------------------------------------
          uc(1, j) = (t14*(utmp(0, j)+utmp(1, j))+t12*(utmp(-1, j)+utmp(&
&           2, j))+t15*(utmp(-2, j)+utmp(3, j)))*rsin_u(1, j)
          uc(2, j) = c1*utmp(3, j) + c2*utmp(2, j) + c3*utmp(1, j)
          ut(0, j) = (uc(0, j)-v(0, j)*cosa_u(0, j))*rsin_u(0, j)
          ut(1, j) = uc(1, j)*rsin_u(1, j)
          ut(2, j) = (uc(2, j)-v(2, j)*cosa_u(2, j))*rsin_u(2, j)
        END DO
      END IF
      IF (ie + 1 .EQ. npx) THEN
        DO j=js-1,je+1
          uc(npx-1, j) = c1*utmp(npx-3, j) + c2*utmp(npx-2, j) + c3*utmp&
&           (npx-1, j)
! 3-pt extrapolation --------------------------------------------------------
          uc(npx, j) = (t14*(utmp(npx-1, j)+utmp(npx, j))+t12*(utmp(npx-&
&           2, j)+utmp(npx+1, j))+t15*(utmp(npx-3, j)+utmp(npx+2, j)))*&
&           rsin_u(npx, j)
          uc(npx+1, j) = c3*utmp(npx, j) + c2*utmp(npx+1, j) + c1*utmp(&
&           npx+2, j)
          ut(npx-1, j) = (uc(npx-1, j)-v(npx-1, j)*cosa_u(npx-1, j))*&
&           rsin_u(npx-1, j)
          ut(npx, j) = uc(npx, j)*rsin_u(npx, j)
          ut(npx+1, j) = (uc(npx+1, j)-v(npx+1, j)*cosa_u(npx+1, j))*&
&           rsin_u(npx+1, j)
        END DO
      END IF
    END IF
!------
! Ydir:
!------
    IF (sw_corner) THEN
      DO j=-2,0
        vtmp(0, j) = -utmp(1-j, 0)
      END DO
    END IF
    IF (nw_corner) THEN
      DO j=0,2
        vtmp(0, npy+j) = utmp(j+1, npy)
      END DO
    END IF
    IF (se_corner) THEN
      DO j=-2,0
        vtmp(npx, j) = utmp(ie+j, 0)
      END DO
    END IF
    IF (ne_corner) THEN
      DO j=0,2
        vtmp(npx, npy+j) = -utmp(ie-j, npy)
      END DO
    END IF
    IF (grid_type .LT. 3) THEN
      DO j=js-1,je+2
        IF (j .EQ. 1) THEN
          DO i=is-1,ie+1
! 3-pt extrapolation -----------------------------------------
            vc(i, 1) = (t14*(vtmp(i, 0)+vtmp(i, 1))+t12*(vtmp(i, -1)+&
&             vtmp(i, 2))+t15*(vtmp(i, -2)+vtmp(i, 3)))*rsin_v(i, 1)
            vt(i, 1) = vc(i, 1)*rsin_v(i, 1)
          END DO
        ELSE IF (j .EQ. 0 .OR. j .EQ. npy - 1) THEN
          DO i=is-1,ie+1
            vc(i, j) = c1*vtmp(i, j-2) + c2*vtmp(i, j-1) + c3*vtmp(i, j)
            vt(i, j) = (vc(i, j)-u(i, j)*cosa_v(i, j))*rsin_v(i, j)
          END DO
        ELSE IF (j .EQ. 2 .OR. j .EQ. npy + 1) THEN
          DO i=is-1,ie+1
            vc(i, j) = c1*vtmp(i, j+1) + c2*vtmp(i, j) + c3*vtmp(i, j-1)
            vt(i, j) = (vc(i, j)-u(i, j)*cosa_v(i, j))*rsin_v(i, j)
          END DO
        ELSE IF (j .EQ. npy) THEN
          DO i=is-1,ie+1
! 3-pt extrapolation --------------------------------------------------------
            vc(i, npy) = (t14*(vtmp(i, npy-1)+vtmp(i, npy))+t12*(vtmp(i&
&             , npy-2)+vtmp(i, npy+1))+t15*(vtmp(i, npy-3)+vtmp(i, npy+2&
&             )))*rsin_v(i, npy)
            vt(i, npy) = vc(i, npy)*rsin_v(i, npy)
          END DO
        ELSE
! 4th order interpolation for interior points:
          DO i=is-1,ie+1
            vc(i, j) = a2*(vtmp(i, j-2)+vtmp(i, j+1)) + a1*(vtmp(i, j-1)&
&             +vtmp(i, j))
            vt(i, j) = (vc(i, j)-u(i, j)*cosa_v(i, j))*rsin_v(i, j)
          END DO
        END IF
      END DO
    ELSE
! 4th order interpolation:
      DO j=js-1,je+2
        DO i=is-1,ie+1
          vc(i, j) = a2*(vtmp(i, j-2)+vtmp(i, j+1)) + a1*(vtmp(i, j-1)+&
&           vtmp(i, j))
          vt(i, j) = vc(i, j)
        END DO
      END DO
    END IF
  END SUBROUTINE D2A2C_VECT
  SUBROUTINE D2A2C_VECT_V2(u, v, ua, va, uc, vc, ut, vt)
    IMPLICIT NONE
    REAL, INTENT(IN) :: u(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: v(isd:ied+1, jsd:jed)
    REAL, DIMENSION(isd:ied + 1, jsd:jed), INTENT(OUT) :: uc
    REAL, DIMENSION(isd:ied, jsd:jed + 1), INTENT(OUT) :: vc
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(OUT) :: ua, va, ut, vt
! Local 
    REAL, DIMENSION(is - 2:ie + 2, js - 2:je + 2) :: wk
    REAL :: utmp, vtmp
    INTEGER :: i, j
! needs only ut[is-1:ie+2,js-1:je+1], vt[is-1:ie+1,js-1:je+2]
    DO j=js-2,je+2
      DO i=is-2,ie+3
        uc(i, j) = v(i, j)*dy(i, j)
      END DO
    END DO
    DO j=js-2,je+3
      DO i=is-2,ie+2
        vc(i, j) = u(i, j)*dx(i, j)
      END DO
    END DO
! D --> A
! Co-variant to Co-variant "vorticity-conserving" interpolation
    DO j=js-2,je+2
      DO i=is-2,ie+2
        utmp = 0.5*(vc(i, j)+vc(i, j+1))*rdxa(i, j)
        vtmp = 0.5*(uc(i, j)+uc(i+1, j))*rdya(i, j)
        ua(i, j) = (utmp-vtmp*cosa_s(i, j))*rsin2(i, j)
        va(i, j) = (vtmp-utmp*cosa_s(i, j))*rsin2(i, j)
      END DO
    END DO
! Xdir:
    IF (sw_corner) THEN
      ua(-1, 0) = -va(0, 2)
      ua(0, 0) = -va(0, 1)
    END IF
    IF (se_corner) THEN
      ua(npx, 0) = va(npx, 1)
      ua(npx+1, 0) = va(npx, 2)
    END IF
    IF (ne_corner) THEN
      ua(npx, npy) = -va(npx, npy-1)
      ua(npx+1, npy) = -va(npx, npy-2)
    END IF
    IF (nw_corner) THEN
      ua(-1, npy) = va(0, npy-2)
      ua(0, npy) = va(0, npy-1)
    END IF
! A -> C
!--------------------------------------------
! Divergence conserving interp to cell walls
!--------------------------------------------
    DO j=js-1,je+1
      DO i=is-2,ie+2
        wk(i, j) = ua(i, j)*dya(i, j)*sina_s(i, j)
      END DO
    END DO
    DO j=js-1,je+1
      DO i=is-1,ie+2
        ut(i, j) = 0.5*(wk(i-1, j)+wk(i, j))/(dy(i, j)*sina_u(i, j))
        uc(i, j) = ut(i, j) + 0.5*(va(i-1, j)*cosa_s(i-1, j)+va(i, j)*&
&         cosa_s(i, j))
      END DO
    END DO
    IF (grid_type .LT. 3) THEN
      IF (is .EQ. 1) THEN
        i = 1
        DO j=js-1,je+1
!          ut(i,j) = 0.75*(ua(i-1,j)+ua(i,j))-0.25*(ua(i-2,j)+ua(i+1,j))
          ut(i, j) = 0.25*(-ua(-1, j)+3.*(ua(0, j)+ua(1, j))-ua(2, j))
          uc(i, j) = ut(i, j)*sina_u(i, j)
        END DO
      END IF
      IF (ie + 1 .EQ. npx) THEN
        i = npx
        DO j=js-1,je+1
!          ut(i,j) = 0.75*(ua(i-1,j)+ua(i,j))-0.25*(ua(i-2,j)+ua(i+1,j))
          ut(i, j) = 0.25*(-ua(i-2, j)+3.*(ua(i-1, j)+ua(i, j))-ua(i+1, &
&           j))
          uc(i, j) = ut(i, j)*sina_u(i, j)
        END DO
      END IF
    END IF
! Ydir:
    IF (sw_corner) THEN
      va(0, -1) = -ua(2, 0)
      va(0, 0) = -ua(1, 0)
    END IF
    IF (se_corner) THEN
      va(npx, 0) = ua(npx-1, 0)
      va(npx, -1) = ua(npx-2, 0)
    END IF
    IF (ne_corner) THEN
      va(npx, npy) = -ua(npx-1, npy)
      va(npx, npy+1) = -ua(npx-2, npy)
    END IF
    IF (nw_corner) THEN
      va(0, npy) = ua(1, npy)
      va(0, npy+1) = ua(2, npy)
    END IF
    DO j=js-2,je+2
      DO i=is-1,ie+1
        wk(i, j) = va(i, j)*dxa(i, j)*sina_s(i, j)
      END DO
    END DO
    IF (grid_type .LT. 3) THEN
      DO j=js-1,je+2
        IF (j .EQ. 1 .OR. j .EQ. npy) THEN
          DO i=is-1,ie+1
            vt(i, j) = 0.25*(-va(i, j-2)+3.*(va(i, j-1)+va(i, j))-va(i, &
&             j+1))
            vc(i, j) = vt(i, j)*sina_v(i, j)
          END DO
        ELSE
          DO i=is-1,ie+1
            vt(i, j) = 0.5*(wk(i, j-1)+wk(i, j))/(dx(i, j)*sina_v(i, j))
            vc(i, j) = vt(i, j) + 0.5*(ua(i, j-1)*cosa_s(i, j-1)+ua(i, j&
&             )*cosa_s(i, j))
          END DO
        END IF
      END DO
    ELSE
      DO j=js-1,je+2
        DO i=is-1,ie+1
          vt(i, j) = 0.5*(wk(i, j-1)+wk(i, j))/(dx(i, j)*sina_v(i, j))
          vc(i, j) = vt(i, j) + 0.5*(ua(i, j-1)*cosa_s(i, j-1)+ua(i, j)*&
&           cosa_s(i, j))
        END DO
      END DO
    END IF
  END SUBROUTINE D2A2C_VECT_V2
  SUBROUTINE D2A2C_VECT_V1(u, v, ua, va, uc, vc, ut, vt)
    IMPLICIT NONE
    REAL, INTENT(IN) :: u(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: v(isd:ied+1, jsd:jed)
    REAL, DIMENSION(isd:ied + 1, jsd:jed), INTENT(OUT) :: uc
    REAL, DIMENSION(isd:ied, jsd:jed + 1), INTENT(OUT) :: vc
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(OUT) :: ua, va, ut, vt
! Local 
    REAL, DIMENSION(isd:ied, jsd:jed) :: v1, v2, v3
    REAL, DIMENSION(isd:ied, jsd:jed) :: utmp, vtmp
    REAL :: vw1, vw2, vw3
    REAL :: vs1, vs2, vs3
    REAL :: up, vp
    INTEGER :: i, j
    INTRINSIC MAX
    INTRINSIC MIN
    INTEGER :: min1
    INTEGER :: max1
! Needs ut[is-1:ie+2,js-1:je+1], vt[is-1:ie+1,js-1:je+2]
    DO j=jsd,jed
      DO i=isd,ied+1
        uc(i, j) = v(i, j)*dy(i, j)
      END DO
    END DO
    DO j=jsd,jed+1
      DO i=isd,ied
        vc(i, j) = u(i, j)*dx(i, j)
      END DO
    END DO
! D --> A
    DO j=jsd,jed
      DO i=isd,ied
        up = 0.5*(vc(i, j)+vc(i, j+1))*rdxa(i, j)
        vp = 0.5*(uc(i, j)+uc(i+1, j))*rdya(i, j)
        ua(i, j) = (up-vp*cosa_s(i, j))*rsin2(i, j)
        va(i, j) = (vp-up*cosa_s(i, j))*rsin2(i, j)
        v1(i, j) = ua(i, j)*ec1(1, i, j) + va(i, j)*ec2(1, i, j)
        v2(i, j) = ua(i, j)*ec1(2, i, j) + va(i, j)*ec2(2, i, j)
        v3(i, j) = ua(i, j)*ec1(3, i, j) + va(i, j)*ec2(3, i, j)
      END DO
    END DO
! A -> C (across face averaging taking place here):
! Xdir
    CALL FILL3_4CORNERS(v1, v2, v3, 1)
!    call copy_corners(v1, npx, npy, 1)
!    call copy_corners(v2, npx, npy, 1)
!    call copy_corners(v3, npx, npy, 1)
! 4th order interpolation:
    DO j=js-1,je+1
      IF (3 .LT. is - 1) THEN
        max1 = is - 1
      ELSE
        max1 = 3
      END IF
      IF (npx - 2 .GT. ie + 2) THEN
        min1 = ie + 2
      ELSE
        min1 = npx - 2
      END IF
      DO i=max1,min1
        vw1 = a2*(v1(i-2, j)+v1(i+1, j)) + a1*(v1(i-1, j)+v1(i, j))
        vw2 = a2*(v2(i-2, j)+v2(i+1, j)) + a1*(v2(i-1, j)+v2(i, j))
        vw3 = a2*(v3(i-2, j)+v3(i+1, j)) + a1*(v3(i-1, j)+v3(i, j))
        uc(i, j) = vw1*ew(1, i, j, 1) + vw2*ew(2, i, j, 1) + vw3*ew(3, i&
&         , j, 1)
        ut(i, j) = (uc(i, j)-v(i, j)*cosa_u(i, j))*rsin_u(i, j)
      END DO
    END DO
! Fix the edge:
    IF (is .EQ. 1) THEN
      DO j=js-1,je+1
        i = 0
        vw1 = c1*v1(-2, j) + c2*v1(-1, j) + c3*v1(0, j)
        vw2 = c1*v2(-2, j) + c2*v2(-1, j) + c3*v2(0, j)
        vw3 = c1*v3(-2, j) + c2*v3(-1, j) + c3*v3(0, j)
        uc(i, j) = vw1*ew(1, i, j, 1) + vw2*ew(2, i, j, 1) + vw3*ew(3, i&
&         , j, 1)
        ut(i, j) = (uc(i, j)-v(i, j)*cosa_u(i, j))*rsin_u(i, j)
        i = 1
        vw1 = 3.*(v1(0, j)+v1(1, j)) - (v1(-1, j)+v1(2, j))
        vw2 = 3.*(v2(0, j)+v2(1, j)) - (v2(-1, j)+v2(2, j))
        vw3 = 3.*(v3(0, j)+v3(1, j)) - (v3(-1, j)+v3(2, j))
        uc(i, j) = 0.25*(vw1*ew(1, i, j, 1)+vw2*ew(2, i, j, 1)+vw3*ew(3&
&         , i, j, 1))
        ut(i, j) = uc(i, j)*rsin_u(i, j)
        i = 2
        vw1 = c3*v1(1, j) + c2*v1(2, j) + c1*v1(3, j)
        vw2 = c3*v2(1, j) + c2*v2(2, j) + c1*v2(3, j)
        vw3 = c3*v3(1, j) + c2*v3(2, j) + c1*v3(3, j)
        uc(i, j) = vw1*ew(1, i, j, 1) + vw2*ew(2, i, j, 1) + vw3*ew(3, i&
&         , j, 1)
        ut(i, j) = (uc(i, j)-v(i, j)*cosa_u(i, j))*rsin_u(i, j)
      END DO
    END IF
    IF (ie + 1 .EQ. npx) THEN
      DO j=js-1,je+1
        i = npx - 1
        vw1 = c1*v1(npx-3, j) + c2*v1(npx-2, j) + c3*v1(npx-1, j)
        vw2 = c1*v2(npx-3, j) + c2*v2(npx-2, j) + c3*v2(npx-1, j)
        vw3 = c1*v3(npx-3, j) + c2*v3(npx-2, j) + c3*v3(npx-1, j)
        uc(i, j) = vw1*ew(1, i, j, 1) + vw2*ew(2, i, j, 1) + vw3*ew(3, i&
&         , j, 1)
        ut(i, j) = (uc(i, j)-v(i, j)*cosa_u(i, j))*rsin_u(i, j)
        i = npx
        vw1 = 3.*(v1(i-1, j)+v1(i, j)) - (v1(i-2, j)+v1(i+1, j))
        vw2 = 3.*(v2(i-1, j)+v2(i, j)) - (v2(i-2, j)+v2(i+1, j))
        vw3 = 3.*(v3(i-1, j)+v3(i, j)) - (v3(i-2, j)+v3(i+1, j))
        uc(i, j) = 0.25*(vw1*ew(1, i, j, 1)+vw2*ew(2, i, j, 1)+vw3*ew(3&
&         , i, j, 1))
        ut(i, j) = uc(i, j)*rsin_u(i, j)
        i = npx + 1
        vw1 = c3*v1(npx, j) + c2*v1(npx+1, j) + c1*v1(npx+2, j)
        vw2 = c3*v2(npx, j) + c2*v2(npx+1, j) + c1*v2(npx+2, j)
        vw3 = c3*v3(npx, j) + c2*v3(npx+1, j) + c1*v3(npx+2, j)
        uc(i, j) = vw1*ew(1, i, j, 1) + vw2*ew(2, i, j, 1) + vw3*ew(3, i&
&         , j, 1)
        ut(i, j) = (uc(i, j)-v(i, j)*cosa_u(i, j))*rsin_u(i, j)
      END DO
    END IF
! Ydir:
    CALL FILL3_4CORNERS(v1, v2, v3, 2)
!    call copy_corners(v1, npx, npy, 2)
!    call copy_corners(v2, npx, npy, 2)
!    call copy_corners(v3, npx, npy, 2)
    DO j=js-1,je+2
      IF (j .EQ. 0 .OR. j .EQ. npy - 1) THEN
        DO i=is-1,ie+1
          vs1 = c1*v1(i, j-2) + c2*v1(i, j-1) + c3*v1(i, j)
          vs2 = c1*v2(i, j-2) + c2*v2(i, j-1) + c3*v2(i, j)
          vs3 = c1*v3(i, j-2) + c2*v3(i, j-1) + c3*v3(i, j)
          vc(i, j) = vs1*es(1, i, j, 2) + vs2*es(2, i, j, 2) + vs3*es(3&
&           , i, j, 2)
          vt(i, j) = (vc(i, j)-u(i, j)*cosa_v(i, j))*rsin_v(i, j)
        END DO
      ELSE IF (j .EQ. 2 .OR. j .EQ. npy + 1) THEN
        DO i=is-1,ie+1
          vs1 = c3*v1(i, j-1) + c2*v1(i, j) + c1*v1(i, j+1)
          vs2 = c3*v2(i, j-1) + c2*v2(i, j) + c1*v2(i, j+1)
          vs3 = c3*v3(i, j-1) + c2*v3(i, j) + c1*v3(i, j+1)
          vc(i, j) = vs1*es(1, i, j, 2) + vs2*es(2, i, j, 2) + vs3*es(3&
&           , i, j, 2)
          vt(i, j) = (vc(i, j)-u(i, j)*cosa_v(i, j))*rsin_v(i, j)
        END DO
      ELSE IF (j .EQ. 1 .OR. j .EQ. npy) THEN
        DO i=is-1,ie+1
          vs1 = 3.*(v1(i, j-1)+v1(i, j)) - (v1(i, j-2)+v1(i, j+1))
          vs2 = 3.*(v2(i, j-1)+v2(i, j)) - (v2(i, j-2)+v2(i, j+1))
          vs3 = 3.*(v3(i, j-1)+v3(i, j)) - (v3(i, j-2)+v3(i, j+1))
          vc(i, j) = 0.25*(vs1*es(1, i, j, 2)+vs2*es(2, i, j, 2)+vs3*es(&
&           3, i, j, 2))
          vt(i, j) = vc(i, j)*rsin_v(i, j)
        END DO
      ELSE
! Interior: 4th order
        DO i=is-1,ie+1
          vs1 = a2*(v1(i, j-2)+v1(i, j+1)) + a1*(v1(i, j-1)+v1(i, j))
          vs2 = a2*(v2(i, j-2)+v2(i, j+1)) + a1*(v2(i, j-1)+v2(i, j))
          vs3 = a2*(v3(i, j-2)+v3(i, j+1)) + a1*(v3(i, j-1)+v3(i, j))
          vc(i, j) = vs1*es(1, i, j, 2) + vs2*es(2, i, j, 2) + vs3*es(3&
&           , i, j, 2)
          vt(i, j) = (vc(i, j)-u(i, j)*cosa_v(i, j))*rsin_v(i, j)
        END DO
      END IF
    END DO
  END SUBROUTINE D2A2C_VECT_V1
  SUBROUTINE FILL3_4CORNERS(q1, q2, q3, dir)
    IMPLICIT NONE
! This routine fill the 4 corners of the scalar fileds only as needed by c_core
! 1: x-dir; 2: y-dir
    INTEGER, INTENT(IN) :: dir
    REAL, INTENT(INOUT) :: q1(isd:ied, jsd:jed)
    REAL, INTENT(INOUT) :: q2(isd:ied, jsd:jed)
    REAL, INTENT(INOUT) :: q3(isd:ied, jsd:jed)
    INTEGER :: i, j
    SELECT CASE  (dir) 
    CASE (1) 
      IF (sw_corner) THEN
        q1(-1, 0) = q1(0, 2)
        q1(0, 0) = q1(0, 1)
        q1(0, -1) = q1(-1, 1)
        q2(-1, 0) = q2(0, 2)
        q2(0, 0) = q2(0, 1)
        q2(0, -1) = q2(-1, 1)
        q3(-1, 0) = q3(0, 2)
        q3(0, 0) = q3(0, 1)
        q3(0, -1) = q3(-1, 1)
      END IF
      IF (se_corner) THEN
        q1(npx+1, 0) = q1(npx, 2)
        q1(npx, 0) = q1(npx, 1)
        q1(npx, -1) = q1(npx+1, 1)
        q2(npx+1, 0) = q2(npx, 2)
        q2(npx, 0) = q2(npx, 1)
        q2(npx, -1) = q2(npx+1, 1)
        q3(npx+1, 0) = q3(npx, 2)
        q3(npx, 0) = q3(npx, 1)
        q3(npx, -1) = q3(npx+1, 1)
      END IF
      IF (ne_corner) THEN
        q1(npx, npy) = q1(npx, npy-1)
        q1(npx+1, npy) = q1(npx, npy-2)
        q1(npx, npy+1) = q1(npx+1, npy-1)
        q2(npx, npy) = q2(npx, npy-1)
        q2(npx+1, npy) = q2(npx, npy-2)
        q2(npx, npy+1) = q2(npx+1, npy-1)
        q3(npx, npy) = q3(npx, npy-1)
        q3(npx+1, npy) = q3(npx, npy-2)
        q3(npx, npy+1) = q3(npx+1, npy-1)
      END IF
      IF (nw_corner) THEN
        q1(0, npy) = q1(0, npy-1)
        q1(-1, npy) = q1(0, npy-2)
        q1(0, npy+1) = q1(-1, npy-1)
        q2(0, npy) = q2(0, npy-1)
        q2(-1, npy) = q2(0, npy-2)
        q2(0, npy+1) = q2(-1, npy-1)
        q3(0, npy) = q3(0, npy-1)
        q3(-1, npy) = q3(0, npy-2)
        q3(0, npy+1) = q3(-1, npy-1)
      END IF
    CASE (2) 
      IF (sw_corner) THEN
        q1(0, 0) = q1(1, 0)
        q1(0, -1) = q1(2, 0)
        q1(-1, 0) = q1(1, -1)
        q2(0, 0) = q2(1, 0)
        q2(0, -1) = q2(2, 0)
        q2(-1, 0) = q2(1, -1)
        q3(0, 0) = q3(1, 0)
        q3(0, -1) = q3(2, 0)
        q3(-1, 0) = q3(1, -1)
      END IF
      IF (se_corner) THEN
        q1(npx, 0) = q1(npx-1, 0)
        q1(npx, -1) = q1(npx-2, 0)
        q1(npx+1, 0) = q1(npx-1, -1)
        q2(npx, 0) = q2(npx-1, 0)
        q2(npx, -1) = q2(npx-2, 0)
        q2(npx+1, 0) = q2(npx-1, -1)
        q3(npx, 0) = q3(npx-1, 0)
        q3(npx, -1) = q3(npx-2, 0)
        q3(npx+1, 0) = q3(npx-1, -1)
      END IF
      IF (ne_corner) THEN
        q1(npx, npy) = q1(npx-1, npy)
        q1(npx, npy+1) = q1(npx-2, npy)
        q1(npx+1, npy) = q1(npx-1, npy+1)
        q2(npx, npy) = q2(npx-1, npy)
        q2(npx, npy+1) = q2(npx-2, npy)
        q2(npx+1, npy) = q2(npx-1, npy+1)
        q3(npx, npy) = q3(npx-1, npy)
        q3(npx, npy+1) = q3(npx-2, npy)
        q3(npx+1, npy) = q3(npx-1, npy+1)
      END IF
      IF (nw_corner) THEN
        q1(0, npy) = q1(1, npy)
        q1(0, npy+1) = q1(2, npy)
        q1(-1, npy) = q1(1, npy+1)
        q2(0, npy) = q2(1, npy)
        q2(0, npy+1) = q2(2, npy)
        q2(-1, npy) = q2(1, npy+1)
        q3(0, npy) = q3(1, npy)
        q3(0, npy+1) = q3(2, npy)
        q3(-1, npy) = q3(1, npy+1)
      END IF
    END SELECT
  END SUBROUTINE FILL3_4CORNERS
!  Differentiation of fill2_4corners in reverse (adjoint) mode (with options r8):
!   gradient     of useful results: q1 q2
!   with respect to varying inputs: q1 q2
  SUBROUTINE FILL2_4CORNERS_ADM(q1, q1_ad, q2, q2_ad, dir)
    IMPLICIT NONE
! This routine fill the 4 corners of the scalar fileds only as needed by c_core
! 1: x-dir; 2: y-dir
    INTEGER, INTENT(IN) :: dir
    REAL, INTENT(INOUT) :: q1(isd:ied, jsd:jed)
    REAL, INTENT(INOUT) :: q1_ad(isd:ied, jsd:jed)
    REAL, INTENT(INOUT) :: q2(isd:ied, jsd:jed)
    REAL, INTENT(INOUT) :: q2_ad(isd:ied, jsd:jed)
    INTEGER :: branch
    SELECT CASE  (dir) 
    CASE (1) 
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
      IF (nw_corner) THEN
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (ne_corner) THEN
        q2_ad(npx, npy-2) = q2_ad(npx, npy-2) + q2_ad(npx+1, npy)
        q2_ad(npx+1, npy) = 0.0_8
        q2_ad(npx, npy-1) = q2_ad(npx, npy-1) + q2_ad(npx, npy)
        q2_ad(npx, npy) = 0.0_8
        q1_ad(npx, npy-2) = q1_ad(npx, npy-2) + q1_ad(npx+1, npy)
        q1_ad(npx+1, npy) = 0.0_8
        q1_ad(npx, npy-1) = q1_ad(npx, npy-1) + q1_ad(npx, npy)
        q1_ad(npx, npy) = 0.0_8
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        q2_ad(0, npy-2) = q2_ad(0, npy-2) + q2_ad(-1, npy)
        q2_ad(-1, npy) = 0.0_8
        q2_ad(0, npy-1) = q2_ad(0, npy-1) + q2_ad(0, npy)
        q2_ad(0, npy) = 0.0_8
        q1_ad(0, npy-2) = q1_ad(0, npy-2) + q1_ad(-1, npy)
        q1_ad(-1, npy) = 0.0_8
        q1_ad(0, npy-1) = q1_ad(0, npy-1) + q1_ad(0, npy)
        q1_ad(0, npy) = 0.0_8
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        q2_ad(npx, 1) = q2_ad(npx, 1) + q2_ad(npx, 0)
        q2_ad(npx, 0) = 0.0_8
        q2_ad(npx, 2) = q2_ad(npx, 2) + q2_ad(npx+1, 0)
        q2_ad(npx+1, 0) = 0.0_8
        q1_ad(npx, 1) = q1_ad(npx, 1) + q1_ad(npx, 0)
        q1_ad(npx, 0) = 0.0_8
        q1_ad(npx, 2) = q1_ad(npx, 2) + q1_ad(npx+1, 0)
        q1_ad(npx+1, 0) = 0.0_8
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        q2_ad(0, 1) = q2_ad(0, 1) + q2_ad(0, 0)
        q2_ad(0, 0) = 0.0_8
        q2_ad(0, 2) = q2_ad(0, 2) + q2_ad(-1, 0)
        q2_ad(-1, 0) = 0.0_8
        q1_ad(0, 1) = q1_ad(0, 1) + q1_ad(0, 0)
        q1_ad(0, 0) = 0.0_8
        q1_ad(0, 2) = q1_ad(0, 2) + q1_ad(-1, 0)
        q1_ad(-1, 0) = 0.0_8
      END IF
    CASE (2) 
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
      IF (nw_corner) THEN
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (ne_corner) THEN
        q2_ad(npx-2, npy) = q2_ad(npx-2, npy) + q2_ad(npx, npy+1)
        q2_ad(npx, npy+1) = 0.0_8
        q2_ad(npx-1, npy) = q2_ad(npx-1, npy) + q2_ad(npx, npy)
        q2_ad(npx, npy) = 0.0_8
        q1_ad(npx-2, npy) = q1_ad(npx-2, npy) + q1_ad(npx, npy+1)
        q1_ad(npx, npy+1) = 0.0_8
        q1_ad(npx-1, npy) = q1_ad(npx-1, npy) + q1_ad(npx, npy)
        q1_ad(npx, npy) = 0.0_8
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        q2_ad(2, npy) = q2_ad(2, npy) + q2_ad(0, npy+1)
        q2_ad(0, npy+1) = 0.0_8
        q2_ad(1, npy) = q2_ad(1, npy) + q2_ad(0, npy)
        q2_ad(0, npy) = 0.0_8
        q1_ad(2, npy) = q1_ad(2, npy) + q1_ad(0, npy+1)
        q1_ad(0, npy+1) = 0.0_8
        q1_ad(1, npy) = q1_ad(1, npy) + q1_ad(0, npy)
        q1_ad(0, npy) = 0.0_8
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        q2_ad(npx-2, 0) = q2_ad(npx-2, 0) + q2_ad(npx, -1)
        q2_ad(npx, -1) = 0.0_8
        q2_ad(npx-1, 0) = q2_ad(npx-1, 0) + q2_ad(npx, 0)
        q2_ad(npx, 0) = 0.0_8
        q1_ad(npx-2, 0) = q1_ad(npx-2, 0) + q1_ad(npx, -1)
        q1_ad(npx, -1) = 0.0_8
        q1_ad(npx-1, 0) = q1_ad(npx-1, 0) + q1_ad(npx, 0)
        q1_ad(npx, 0) = 0.0_8
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        q2_ad(2, 0) = q2_ad(2, 0) + q2_ad(0, -1)
        q2_ad(0, -1) = 0.0_8
        q2_ad(1, 0) = q2_ad(1, 0) + q2_ad(0, 0)
        q2_ad(0, 0) = 0.0_8
        q1_ad(2, 0) = q1_ad(2, 0) + q1_ad(0, -1)
        q1_ad(0, -1) = 0.0_8
        q1_ad(1, 0) = q1_ad(1, 0) + q1_ad(0, 0)
        q1_ad(0, 0) = 0.0_8
      END IF
    END SELECT
  END SUBROUTINE FILL2_4CORNERS_ADM
  SUBROUTINE FILL2_4CORNERS(q1, q2, dir)
    IMPLICIT NONE
! This routine fill the 4 corners of the scalar fileds only as needed by c_core
! 1: x-dir; 2: y-dir
    INTEGER, INTENT(IN) :: dir
    REAL, INTENT(INOUT) :: q1(isd:ied, jsd:jed)
    REAL, INTENT(INOUT) :: q2(isd:ied, jsd:jed)
    SELECT CASE  (dir) 
    CASE (1) 
      IF (sw_corner) THEN
        q1(-1, 0) = q1(0, 2)
        q1(0, 0) = q1(0, 1)
        q2(-1, 0) = q2(0, 2)
        q2(0, 0) = q2(0, 1)
      END IF
      IF (se_corner) THEN
        q1(npx+1, 0) = q1(npx, 2)
        q1(npx, 0) = q1(npx, 1)
        q2(npx+1, 0) = q2(npx, 2)
        q2(npx, 0) = q2(npx, 1)
      END IF
      IF (nw_corner) THEN
        q1(0, npy) = q1(0, npy-1)
        q1(-1, npy) = q1(0, npy-2)
        q2(0, npy) = q2(0, npy-1)
        q2(-1, npy) = q2(0, npy-2)
      END IF
      IF (ne_corner) THEN
        q1(npx, npy) = q1(npx, npy-1)
        q1(npx+1, npy) = q1(npx, npy-2)
        q2(npx, npy) = q2(npx, npy-1)
        q2(npx+1, npy) = q2(npx, npy-2)
      END IF
    CASE (2) 
      IF (sw_corner) THEN
        q1(0, 0) = q1(1, 0)
        q1(0, -1) = q1(2, 0)
        q2(0, 0) = q2(1, 0)
        q2(0, -1) = q2(2, 0)
      END IF
      IF (se_corner) THEN
        q1(npx, 0) = q1(npx-1, 0)
        q1(npx, -1) = q1(npx-2, 0)
        q2(npx, 0) = q2(npx-1, 0)
        q2(npx, -1) = q2(npx-2, 0)
      END IF
      IF (nw_corner) THEN
        q1(0, npy) = q1(1, npy)
        q1(0, npy+1) = q1(2, npy)
        q2(0, npy) = q2(1, npy)
        q2(0, npy+1) = q2(2, npy)
      END IF
      IF (ne_corner) THEN
        q1(npx, npy) = q1(npx-1, npy)
        q1(npx, npy+1) = q1(npx-2, npy)
        q2(npx, npy) = q2(npx-1, npy)
        q2(npx, npy+1) = q2(npx-2, npy)
      END IF
    END SELECT
  END SUBROUTINE FILL2_4CORNERS
!  Differentiation of fill_4corners in reverse (adjoint) mode (with options r8):
!   gradient     of useful results: q
!   with respect to varying inputs: q
  SUBROUTINE FILL_4CORNERS_ADM(q, q_ad, dir)
    IMPLICIT NONE
! This routine fill the 4 corners of the scalar fileds only as needed by c_core
! 1: x-dir; 2: y-dir
    INTEGER, INTENT(IN) :: dir
    REAL, INTENT(INOUT) :: q(isd:ied, jsd:jed)
    REAL, INTENT(INOUT) :: q_ad(isd:ied, jsd:jed)
    INTEGER :: branch
    SELECT CASE  (dir) 
    CASE (1) 
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
      IF (nw_corner) THEN
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (ne_corner) THEN
        q_ad(npx, npy-2) = q_ad(npx, npy-2) + q_ad(npx+1, npy)
        q_ad(npx+1, npy) = 0.0_8
        q_ad(npx, npy-1) = q_ad(npx, npy-1) + q_ad(npx, npy)
        q_ad(npx, npy) = 0.0_8
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        q_ad(0, npy-2) = q_ad(0, npy-2) + q_ad(-1, npy)
        q_ad(-1, npy) = 0.0_8
        q_ad(0, npy-1) = q_ad(0, npy-1) + q_ad(0, npy)
        q_ad(0, npy) = 0.0_8
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        q_ad(npx, 1) = q_ad(npx, 1) + q_ad(npx, 0)
        q_ad(npx, 0) = 0.0_8
        q_ad(npx, 2) = q_ad(npx, 2) + q_ad(npx+1, 0)
        q_ad(npx+1, 0) = 0.0_8
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        q_ad(0, 1) = q_ad(0, 1) + q_ad(0, 0)
        q_ad(0, 0) = 0.0_8
        q_ad(0, 2) = q_ad(0, 2) + q_ad(-1, 0)
        q_ad(-1, 0) = 0.0_8
      END IF
    CASE (2) 
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
      IF (nw_corner) THEN
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (ne_corner) THEN
        q_ad(npx-2, npy) = q_ad(npx-2, npy) + q_ad(npx, npy+1)
        q_ad(npx, npy+1) = 0.0_8
        q_ad(npx-1, npy) = q_ad(npx-1, npy) + q_ad(npx, npy)
        q_ad(npx, npy) = 0.0_8
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        q_ad(2, npy) = q_ad(2, npy) + q_ad(0, npy+1)
        q_ad(0, npy+1) = 0.0_8
        q_ad(1, npy) = q_ad(1, npy) + q_ad(0, npy)
        q_ad(0, npy) = 0.0_8
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        q_ad(npx-2, 0) = q_ad(npx-2, 0) + q_ad(npx, -1)
        q_ad(npx, -1) = 0.0_8
        q_ad(npx-1, 0) = q_ad(npx-1, 0) + q_ad(npx, 0)
        q_ad(npx, 0) = 0.0_8
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        q_ad(2, 0) = q_ad(2, 0) + q_ad(0, -1)
        q_ad(0, -1) = 0.0_8
        q_ad(1, 0) = q_ad(1, 0) + q_ad(0, 0)
        q_ad(0, 0) = 0.0_8
      END IF
    END SELECT
  END SUBROUTINE FILL_4CORNERS_ADM
  SUBROUTINE FILL_4CORNERS(q, dir)
    IMPLICIT NONE
! This routine fill the 4 corners of the scalar fileds only as needed by c_core
! 1: x-dir; 2: y-dir
    INTEGER, INTENT(IN) :: dir
    REAL, INTENT(INOUT) :: q(isd:ied, jsd:jed)
    SELECT CASE  (dir) 
    CASE (1) 
      IF (sw_corner) THEN
        q(-1, 0) = q(0, 2)
        q(0, 0) = q(0, 1)
      END IF
      IF (se_corner) THEN
        q(npx+1, 0) = q(npx, 2)
        q(npx, 0) = q(npx, 1)
      END IF
      IF (nw_corner) THEN
        q(0, npy) = q(0, npy-1)
        q(-1, npy) = q(0, npy-2)
      END IF
      IF (ne_corner) THEN
        q(npx, npy) = q(npx, npy-1)
        q(npx+1, npy) = q(npx, npy-2)
      END IF
    CASE (2) 
      IF (sw_corner) THEN
        q(0, 0) = q(1, 0)
        q(0, -1) = q(2, 0)
      END IF
      IF (se_corner) THEN
        q(npx, 0) = q(npx-1, 0)
        q(npx, -1) = q(npx-2, 0)
      END IF
      IF (nw_corner) THEN
        q(0, npy) = q(1, npy)
        q(0, npy+1) = q(2, npy)
      END IF
      IF (ne_corner) THEN
        q(npx, npy) = q(npx-1, npy)
        q(npx, npy+1) = q(npx-2, npy)
      END IF
    END SELECT
  END SUBROUTINE FILL_4CORNERS
!  Differentiation of d2a2c in reverse (adjoint) mode (with options r8):
!   gradient     of useful results: u v ua uc um va vc vm
!   with respect to varying inputs: u v ua uc um va vc vm
!#ifdef REL_VOR_DMP
! subroutine del6_flux( nord, npx, npy, damp, q, u, v)
!! Del-nord damping for the relative vorticity
!!------------------
!! nord = 0:   del-2
!! nord = 1:   del-4
!! nord = 2:   del-6
!!------------------
!   integer, intent(in):: nord            ! del-n
!   integer, intent(in):: npx, npy
!   real, intent(in):: damp
!   real, intent(in):: q(isd:ied, jsd:jed)  ! q ghosted on input
!   real, intent(inout),  dimension(isd:ied,  jsd:jed+1):: u
!   real, intent(inout),  dimension(isd:ied+1,jsd:jed  ):: v
!! local:
!   real fx2(isd:ied+1,jsd:jed), fy2(isd:ied,jsd:jed+1)
!   real d2(isd:ied,jsd:jed)
!   integer i,j, n, nt
!
!
!   do j=jsd,jed
!      do i=isd,ied
!         d2(i,j) = damp*q(i,j)
!      enddo
!   enddo
!
!   if( nord>0 ) call copy_corners(d2, npx, npy, 1)
!   do j=js-nord,je+nord
!      do i=is-nord,ie+nord+1
!         fx2(i,j) = dy(i,j)*sina_u(i,j)*(d2(i-1,j)-d2(i,j))*rdxc(i,j)
!      enddo
!   enddo
!
!   if( nord>0 ) call copy_corners(d2, npx, npy, 2)
!   do j=js-nord,je+nord+1
!      do i=is-nord,ie+nord
!         fy2(i,j) = dx(i,j)*sina_v(i,j)*(d2(i,j-1)-d2(i,j))*rdyc(i,j)
!      enddo
!   enddo
!
!   if ( nord>0 ) then
!
!!----------
!! high-order
!!----------
!
!   do n=1, nord
!
!      nt = nord-n
!
!      do j=js-nt-1,je+nt+1
!         do i=is-nt-1,ie+nt+1
!            d2(i,j) = (fx2(i,j)-fx2(i+1,j)+fy2(i,j)-fy2(i,j+1))*rarea(i,j)
!         enddo
!      enddo
!
!      call copy_corners(d2, npx, npy, 1)
!      do j=js-nt,je+nt
!         do i=is-nt,ie+nt+1
!            fx2(i,j) = dy(i,j)*sina_u(i,j)*(d2(i,j)-d2(i-1,j))*rdxc(i,j)
!         enddo
!      enddo
!
!      call copy_corners(d2, npx, npy, 2)
!      do j=js-nt,je+nt+1
!         do i=is-nt,ie+nt
!            fy2(i,j) = dx(i,j)*sina_v(i,j)*(d2(i,j)-d2(i,j-1))*rdyc(i,j)
!         enddo
!      enddo
!   enddo
!
!   endif
!
!   do j=js,je
!      do i=is,ie+1
!         u(i,j) = u(i,j) + fx2(i,j)
!      enddo
!   enddo
!
!   do j=js,je+1
!      do i=is,ie
!         v(i,j) = v(i,j) - fy2(i,j)
!      enddo
!   enddo
!
! end subroutine del6_flux
!#endif
  SUBROUTINE D2A2C_ADM(u, u_ad, v, v_ad, um, um_ad, vm, vm_ad, ua, ua_ad&
&   , va, va_ad, uc, uc_ad, vc, vc_ad, dord4)
    IMPLICIT NONE
    REAL, DIMENSION(isd:ied, jsd:jed + 1), INTENT(IN) :: u, um
    REAL, DIMENSION(isd:ied, jsd:jed+1) :: u_ad, um_ad
    REAL, DIMENSION(isd:ied + 1, jsd:jed), INTENT(IN) :: v, vm
    REAL, DIMENSION(isd:ied+1, jsd:jed) :: v_ad, vm_ad
    LOGICAL, INTENT(IN) :: dord4
    REAL, DIMENSION(isd:ied + 1, jsd:jed) :: uc
    REAL, DIMENSION(isd:ied+1, jsd:jed) :: uc_ad
    REAL, DIMENSION(isd:ied, jsd:jed + 1) :: vc
    REAL, DIMENSION(isd:ied, jsd:jed+1) :: vc_ad
    REAL, DIMENSION(isd:ied, jsd:jed) :: ua, va
    REAL, DIMENSION(isd:ied, jsd:jed) :: ua_ad, va_ad
! Local 
    REAL, DIMENSION(isd:ied, jsd:jed) :: utmp, vtmp
    REAL, DIMENSION(isd:ied, jsd:jed) :: utmp_ad, vtmp_ad
    INTEGER :: npt, i, j, ifirst, ilast, id
    INTRINSIC MAX
    INTRINSIC MIN
    INTEGER :: ad_from
    INTEGER :: ad_to
    INTEGER :: ad_from0
    INTEGER :: ad_to0
    INTEGER :: ad_from1
    INTEGER :: ad_to1
    INTEGER :: ad_from2
    INTEGER :: ad_to2
    INTEGER :: branch
    INTEGER :: min9
    INTEGER :: min8
    INTEGER :: min7
    INTEGER :: min6
    INTEGER :: min5
    INTEGER :: min4
    INTEGER :: min3
    INTEGER :: min2
    INTEGER :: min1
    REAL :: temp_ad4
    REAL :: temp_ad3
    REAL :: temp_ad2
    REAL :: temp_ad1
    REAL :: temp_ad0
    REAL :: temp_ad
    INTEGER :: min12
    INTEGER :: min11
    INTEGER :: min10
    INTEGER :: max9
    INTEGER :: max8
    INTEGER :: max7
    INTEGER :: max6
    INTEGER :: max5
    INTEGER :: max4
    INTEGER :: max3
    INTEGER :: max2
    INTEGER :: max1
    INTEGER :: max12
    INTEGER :: max11
    INTEGER :: max10
    IF (dord4) THEN
      id = 1
    ELSE
      id = 0
    END IF
    IF (grid_type .LT. 3) THEN
      npt = 4
    ELSE
      npt = -2
    END IF
    IF (npt .LT. js - 1) THEN
      max1 = js - 1
    ELSE
      max1 = npt
    END IF
    IF (npy - npt .GT. je + 1) THEN
      min1 = je + 1
    ELSE
      min1 = npy - npt
    END IF
!----------
! Interior:
!----------
    DO j=max1,min1
      IF (npt .LT. isd) THEN
        max2 = isd
      ELSE
        max2 = npt
      END IF
      IF (npx - npt .GT. ied) THEN
        min2 = ied
      ELSE
        min2 = npx - npt
      END IF
      ad_from = max2
      i = min2 + 1
      CALL PUSHINTEGER4(i - 1)
      CALL PUSHINTEGER4(ad_from)
    END DO
    IF (npt .LT. jsd) THEN
      max3 = jsd
    ELSE
      max3 = npt
    END IF
    IF (npy - npt .GT. jed) THEN
      min3 = jed
    ELSE
      min3 = npy - npt
    END IF
    DO j=max3,min3
      IF (npt .LT. is - 1) THEN
        max4 = is - 1
      ELSE
        max4 = npt
      END IF
      IF (npx - npt .GT. ie + 1) THEN
        min4 = ie + 1
      ELSE
        min4 = npx - npt
      END IF
      ad_from0 = max4
      i = min4 + 1
      CALL PUSHINTEGER4(i - 1)
      CALL PUSHINTEGER4(ad_from0)
    END DO
!----------
! edges:
!----------
    IF (grid_type .LT. 3) THEN
      IF (js .EQ. 1 .OR. jsd .LT. npt) THEN
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (je + 1 .EQ. npy .OR. jed .GE. npy - npt) THEN
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (is .EQ. 1 .OR. isd .LT. npt) THEN
        IF (npt .LT. jsd) THEN
          max5 = jsd
        ELSE
          max5 = npt
        END IF
        IF (npy - npt .GT. jed) THEN
          min5 = jed
        ELSE
          min5 = npy - npt
        END IF
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (ie + 1 .EQ. npx .OR. ied .GE. npx - npt) THEN
        IF (npt .LT. jsd) THEN
          max6 = jsd
        ELSE
          max6 = npt
        END IF
        IF (npy - npt .GT. jed) THEN
          min6 = jed
        ELSE
          min6 = npy - npt
        END IF
        CALL PUSHCONTROL2B(2)
      ELSE
        CALL PUSHCONTROL2B(1)
      END IF
    ELSE
      CALL PUSHCONTROL2B(0)
    END IF
    IF (npt .LT. js - 1) THEN
      max7 = js - 1
    ELSE
      max7 = npt
    END IF
    IF (npy - npt .GT. je + 1) THEN
      min7 = je + 1
    ELSE
      min7 = npy - npt
    END IF
! Re-compute (utmp, vtmp) using (um,vm)
!----------
! Interior:
!----------
    DO j=max7,min7
      IF (npt .LT. isd) THEN
        max8 = isd
      ELSE
        max8 = npt
      END IF
      IF (npx - npt .GT. ied) THEN
        min8 = ied
      ELSE
        min8 = npx - npt
      END IF
      ad_from1 = max8
      i = min8 + 1
      CALL PUSHINTEGER4(i - 1)
      CALL PUSHINTEGER4(ad_from1)
    END DO
    IF (npt .LT. jsd) THEN
      max9 = jsd
    ELSE
      max9 = npt
    END IF
    IF (npy - npt .GT. jed) THEN
      min9 = jed
    ELSE
      min9 = npy - npt
    END IF
    DO j=max9,min9
      IF (npt .LT. is - 1) THEN
        max10 = is - 1
      ELSE
        max10 = npt
      END IF
      IF (npx - npt .GT. ie + 1) THEN
        min10 = ie + 1
      ELSE
        min10 = npx - npt
      END IF
      ad_from2 = max10
      i = min10 + 1
      CALL PUSHINTEGER4(i - 1)
      CALL PUSHINTEGER4(ad_from2)
    END DO
!----------
! edges:
!----------
    IF (grid_type .LT. 3) THEN
      IF (js .EQ. 1 .OR. jsd .LT. npt) THEN
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (je + 1 .EQ. npy .OR. jed .GE. npy - npt) THEN
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (is .EQ. 1 .OR. isd .LT. npt) THEN
        IF (npt .LT. jsd) THEN
          max11 = jsd
        ELSE
          max11 = npt
        END IF
        IF (npy - npt .GT. jed) THEN
          min11 = jed
        ELSE
          min11 = npy - npt
        END IF
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (ie + 1 .EQ. npx .OR. ied .GE. npx - npt) THEN
        IF (npt .LT. jsd) THEN
          max12 = jsd
        ELSE
          max12 = npt
        END IF
        IF (npy - npt .GT. jed) THEN
          min12 = jed
        ELSE
          min12 = npy - npt
        END IF
        CALL PUSHCONTROL2B(0)
      ELSE
        CALL PUSHCONTROL2B(1)
      END IF
    ELSE
      CALL PUSHCONTROL2B(2)
    END IF
! A -> C
!--------------
! Fix the edges
!--------------
! Xdir:
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
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    IF (grid_type .LT. 3) THEN
      IF (3 .LT. is) THEN
        ifirst = is
      ELSE
        ifirst = 3
      END IF
      IF (npx - 2 .GT. ie + 1) THEN
        CALL PUSHCONTROL1B(1)
        ilast = ie + 1
      ELSE
        CALL PUSHCONTROL1B(1)
        ilast = npx - 2
      END IF
    ELSE
      CALL PUSHCONTROL1B(0)
      ifirst = is
      ilast = ie + 1
    END IF
    IF (grid_type .LT. 3) THEN
      IF (is .EQ. 1) THEN
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (ie + 1 .EQ. npx) THEN
        CALL PUSHCONTROL2B(0)
      ELSE
        CALL PUSHCONTROL2B(1)
      END IF
    ELSE
      CALL PUSHCONTROL2B(2)
    END IF
!------
! Ydir:
!------
    IF (sw_corner) THEN
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    IF (nw_corner) THEN
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
    IF (grid_type .LT. 3) THEN
      DO j=js,je+1
        IF (j .EQ. 1) THEN
          CALL PUSHCONTROL3B(4)
        ELSE IF (j .EQ. npy - 1) THEN
          CALL PUSHCONTROL3B(3)
        ELSE IF (j .EQ. 2) THEN
          CALL PUSHCONTROL3B(2)
        ELSE IF (j .EQ. npy) THEN
          CALL PUSHCONTROL3B(1)
        ELSE
          CALL PUSHCONTROL3B(0)
        END IF
      END DO
      vtmp_ad = 0.0_8
      DO j=je+1,js,-1
        CALL POPCONTROL3B(branch)
        IF (branch .LT. 2) THEN
          IF (branch .EQ. 0) THEN
            DO i=ie,is,-1
              vtmp_ad(i, j-2) = vtmp_ad(i, j-2) + a2*vc_ad(i, j)
              vtmp_ad(i, j+1) = vtmp_ad(i, j+1) + a2*vc_ad(i, j)
              vtmp_ad(i, j-1) = vtmp_ad(i, j-1) + a1*vc_ad(i, j)
              vtmp_ad(i, j) = vtmp_ad(i, j) + a1*vc_ad(i, j)
              vc_ad(i, j) = 0.0_8
            END DO
          ELSE
            DO i=ie,is,-1
              temp_ad4 = rsin_v(i, npy)*vc_ad(i, npy)
              vtmp_ad(i, npy-1) = vtmp_ad(i, npy-1) + t14*temp_ad4
              vtmp_ad(i, npy) = vtmp_ad(i, npy) + t14*temp_ad4
              vtmp_ad(i, npy-2) = vtmp_ad(i, npy-2) + t12*temp_ad4
              vtmp_ad(i, npy+1) = vtmp_ad(i, npy+1) + t12*temp_ad4
              vtmp_ad(i, npy-3) = vtmp_ad(i, npy-3) + t15*temp_ad4
              vtmp_ad(i, npy+2) = vtmp_ad(i, npy+2) + t15*temp_ad4
              vc_ad(i, npy) = 0.0_8
            END DO
          END IF
        ELSE IF (branch .EQ. 2) THEN
          DO i=ie,is,-1
            vtmp_ad(i, j-1) = vtmp_ad(i, j-1) + c3*vc_ad(i, j)
            vtmp_ad(i, j) = vtmp_ad(i, j) + c2*vc_ad(i, j)
            vtmp_ad(i, j+1) = vtmp_ad(i, j+1) + c1*vc_ad(i, j)
            vc_ad(i, j) = 0.0_8
          END DO
        ELSE IF (branch .EQ. 3) THEN
          DO i=ie,is,-1
            vtmp_ad(i, j-2) = vtmp_ad(i, j-2) + c1*vc_ad(i, j)
            vtmp_ad(i, j-1) = vtmp_ad(i, j-1) + c2*vc_ad(i, j)
            vtmp_ad(i, j) = vtmp_ad(i, j) + c3*vc_ad(i, j)
            vc_ad(i, j) = 0.0_8
          END DO
        ELSE
          DO i=ie,is,-1
            temp_ad3 = rsin_v(i, 1)*vc_ad(i, 1)
            vtmp_ad(i, 0) = vtmp_ad(i, 0) + t14*temp_ad3
            vtmp_ad(i, 1) = vtmp_ad(i, 1) + t14*temp_ad3
            vtmp_ad(i, -1) = vtmp_ad(i, -1) + t12*temp_ad3
            vtmp_ad(i, 2) = vtmp_ad(i, 2) + t12*temp_ad3
            vtmp_ad(i, -2) = vtmp_ad(i, -2) + t15*temp_ad3
            vtmp_ad(i, 3) = vtmp_ad(i, 3) + t15*temp_ad3
            vc_ad(i, 1) = 0.0_8
          END DO
        END IF
      END DO
    ELSE
      vtmp_ad = 0.0_8
      DO j=je+1,js,-1
        DO i=ie,is,-1
          vtmp_ad(i, j-2) = vtmp_ad(i, j-2) + a2*vc_ad(i, j)
          vtmp_ad(i, j+1) = vtmp_ad(i, j+1) + a2*vc_ad(i, j)
          vtmp_ad(i, j-1) = vtmp_ad(i, j-1) + a1*vc_ad(i, j)
          vtmp_ad(i, j) = vtmp_ad(i, j) + a1*vc_ad(i, j)
          vc_ad(i, j) = 0.0_8
        END DO
      END DO
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      utmp_ad = 0.0_8
      DO j=2,0,-1
        utmp_ad(ie-j, npy) = utmp_ad(ie-j, npy) - vtmp_ad(npx, npy+j)
        vtmp_ad(npx, npy+j) = 0.0_8
      END DO
    ELSE
      utmp_ad = 0.0_8
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      DO j=0,-2,-1
        utmp_ad(ie+j, 0) = utmp_ad(ie+j, 0) + vtmp_ad(npx, j)
        vtmp_ad(npx, j) = 0.0_8
      END DO
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      DO j=2,0,-1
        utmp_ad(j+1, npy) = utmp_ad(j+1, npy) + vtmp_ad(0, npy+j)
        vtmp_ad(0, npy+j) = 0.0_8
      END DO
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      DO j=0,-2,-1
        utmp_ad(1-j, 0) = utmp_ad(1-j, 0) - vtmp_ad(0, j)
        vtmp_ad(0, j) = 0.0_8
      END DO
    END IF
    CALL POPCONTROL2B(branch)
    IF (branch .EQ. 0) THEN
      DO j=je,js,-1
        temp_ad2 = rsin_u(npx, j)*uc_ad(npx, j)
        utmp_ad(npx-1, j) = utmp_ad(npx-1, j) + t14*temp_ad2
        utmp_ad(npx, j) = utmp_ad(npx, j) + t14*temp_ad2
        utmp_ad(npx-2, j) = utmp_ad(npx-2, j) + t12*temp_ad2
        utmp_ad(npx+1, j) = utmp_ad(npx+1, j) + t12*temp_ad2
        utmp_ad(npx-3, j) = utmp_ad(npx-3, j) + t15*temp_ad2
        utmp_ad(npx+2, j) = utmp_ad(npx+2, j) + t15*temp_ad2
        uc_ad(npx, j) = 0.0_8
        utmp_ad(npx-3, j) = utmp_ad(npx-3, j) + c1*uc_ad(npx-1, j)
        utmp_ad(npx-2, j) = utmp_ad(npx-2, j) + c2*uc_ad(npx-1, j)
        utmp_ad(npx-1, j) = utmp_ad(npx-1, j) + c3*uc_ad(npx-1, j)
        uc_ad(npx-1, j) = 0.0_8
      END DO
    ELSE IF (branch .NE. 1) THEN
      GOTO 100
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      DO j=je,js,-1
        utmp_ad(3, j) = utmp_ad(3, j) + c1*uc_ad(2, j)
        utmp_ad(2, j) = utmp_ad(2, j) + c2*uc_ad(2, j)
        utmp_ad(1, j) = utmp_ad(1, j) + c3*uc_ad(2, j)
        uc_ad(2, j) = 0.0_8
        temp_ad1 = rsin_u(1, j)*uc_ad(1, j)
        utmp_ad(0, j) = utmp_ad(0, j) + t14*temp_ad1
        utmp_ad(1, j) = utmp_ad(1, j) + t14*temp_ad1
        utmp_ad(-1, j) = utmp_ad(-1, j) + t12*temp_ad1
        utmp_ad(2, j) = utmp_ad(2, j) + t12*temp_ad1
        utmp_ad(-2, j) = utmp_ad(-2, j) + t15*temp_ad1
        utmp_ad(3, j) = utmp_ad(3, j) + t15*temp_ad1
        uc_ad(1, j) = 0.0_8
      END DO
    END IF
 100 DO j=je,js,-1
      DO i=ilast,ifirst,-1
        utmp_ad(i-1, j) = utmp_ad(i-1, j) + a1*uc_ad(i, j)
        utmp_ad(i, j) = utmp_ad(i, j) + a1*uc_ad(i, j)
        utmp_ad(i-2, j) = utmp_ad(i-2, j) + a2*uc_ad(i, j)
        utmp_ad(i+1, j) = utmp_ad(i+1, j) + a2*uc_ad(i, j)
        uc_ad(i, j) = 0.0_8
      END DO
    END DO
    CALL POPCONTROL1B(branch)
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      DO i=0,-2,-1
        vtmp_ad(0, je+i) = vtmp_ad(0, je+i) + utmp_ad(i, npy)
        utmp_ad(i, npy) = 0.0_8
      END DO
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      DO i=2,0,-1
        vtmp_ad(npx, je-i) = vtmp_ad(npx, je-i) - utmp_ad(npx+i, npy)
        utmp_ad(npx+i, npy) = 0.0_8
      END DO
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      DO i=2,0,-1
        vtmp_ad(npx, i+1) = vtmp_ad(npx, i+1) + utmp_ad(npx+i, 0)
        utmp_ad(npx+i, 0) = 0.0_8
      END DO
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      DO i=0,-2,-1
        vtmp_ad(0, 1-i) = vtmp_ad(0, 1-i) - utmp_ad(i, 0)
        utmp_ad(i, 0) = 0.0_8
      END DO
    END IF
    CALL POPCONTROL2B(branch)
    IF (branch .EQ. 0) THEN
      DO j=min12,max12,-1
        DO i=ied,npx-npt+1,-1
          vm_ad(i, j) = vm_ad(i, j) + 0.5*vtmp_ad(i, j)
          vm_ad(i+1, j) = vm_ad(i+1, j) + 0.5*vtmp_ad(i, j)
          vtmp_ad(i, j) = 0.0_8
          um_ad(i, j) = um_ad(i, j) + 0.5*utmp_ad(i, j)
          um_ad(i, j+1) = um_ad(i, j+1) + 0.5*utmp_ad(i, j)
          utmp_ad(i, j) = 0.0_8
        END DO
      END DO
    ELSE IF (branch .NE. 1) THEN
      GOTO 110
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      DO j=min11,max11,-1
        DO i=npt-1,isd,-1
          vm_ad(i, j) = vm_ad(i, j) + 0.5*vtmp_ad(i, j)
          vm_ad(i+1, j) = vm_ad(i+1, j) + 0.5*vtmp_ad(i, j)
          vtmp_ad(i, j) = 0.0_8
          um_ad(i, j) = um_ad(i, j) + 0.5*utmp_ad(i, j)
          um_ad(i, j+1) = um_ad(i, j+1) + 0.5*utmp_ad(i, j)
          utmp_ad(i, j) = 0.0_8
        END DO
      END DO
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      DO j=jed,npy-npt+1,-1
        DO i=ied,isd,-1
          vm_ad(i, j) = vm_ad(i, j) + 0.5*vtmp_ad(i, j)
          vm_ad(i+1, j) = vm_ad(i+1, j) + 0.5*vtmp_ad(i, j)
          vtmp_ad(i, j) = 0.0_8
          um_ad(i, j) = um_ad(i, j) + 0.5*utmp_ad(i, j)
          um_ad(i, j+1) = um_ad(i, j+1) + 0.5*utmp_ad(i, j)
          utmp_ad(i, j) = 0.0_8
        END DO
      END DO
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      DO j=npt-1,jsd,-1
        DO i=ied,isd,-1
          vm_ad(i, j) = vm_ad(i, j) + 0.5*vtmp_ad(i, j)
          vm_ad(i+1, j) = vm_ad(i+1, j) + 0.5*vtmp_ad(i, j)
          vtmp_ad(i, j) = 0.0_8
          um_ad(i, j) = um_ad(i, j) + 0.5*utmp_ad(i, j)
          um_ad(i, j+1) = um_ad(i, j+1) + 0.5*utmp_ad(i, j)
          utmp_ad(i, j) = 0.0_8
        END DO
      END DO
    END IF
 110 DO j=min9,max9,-1
      CALL POPINTEGER4(ad_from2)
      CALL POPINTEGER4(ad_to2)
      DO i=ad_to2,ad_from2,-1
        vm_ad(i-1, j) = vm_ad(i-1, j) + a2*vtmp_ad(i, j)
        vm_ad(i+2, j) = vm_ad(i+2, j) + a2*vtmp_ad(i, j)
        vm_ad(i, j) = vm_ad(i, j) + a1*vtmp_ad(i, j)
        vm_ad(i+1, j) = vm_ad(i+1, j) + a1*vtmp_ad(i, j)
        vtmp_ad(i, j) = 0.0_8
      END DO
    END DO
    DO j=min7,max7,-1
      CALL POPINTEGER4(ad_from1)
      CALL POPINTEGER4(ad_to1)
      DO i=ad_to1,ad_from1,-1
        um_ad(i, j-1) = um_ad(i, j-1) + a2*utmp_ad(i, j)
        um_ad(i, j+2) = um_ad(i, j+2) + a2*utmp_ad(i, j)
        um_ad(i, j) = um_ad(i, j) + a1*utmp_ad(i, j)
        um_ad(i, j+1) = um_ad(i, j+1) + a1*utmp_ad(i, j)
        utmp_ad(i, j) = 0.0_8
      END DO
    END DO
    DO j=je+id+1,js-1-id,-1
      DO i=ie+id+1,is-1-id,-1
        temp_ad0 = rsin2(i, j)*ua_ad(i, j)
        temp_ad = rsin2(i, j)*va_ad(i, j)
        vtmp_ad(i, j) = vtmp_ad(i, j) + temp_ad - cosa_s(i, j)*temp_ad0
        utmp_ad(i, j) = utmp_ad(i, j) + temp_ad0 - cosa_s(i, j)*temp_ad
        va_ad(i, j) = 0.0_8
        ua_ad(i, j) = 0.0_8
      END DO
    END DO
    CALL POPCONTROL2B(branch)
    IF (branch .NE. 0) THEN
      IF (branch .NE. 1) THEN
        DO j=min6,max6,-1
          DO i=ied,npx-npt+1,-1
            v_ad(i, j) = v_ad(i, j) + 0.5*vtmp_ad(i, j)
            v_ad(i+1, j) = v_ad(i+1, j) + 0.5*vtmp_ad(i, j)
            vtmp_ad(i, j) = 0.0_8
            u_ad(i, j) = u_ad(i, j) + 0.5*utmp_ad(i, j)
            u_ad(i, j+1) = u_ad(i, j+1) + 0.5*utmp_ad(i, j)
            utmp_ad(i, j) = 0.0_8
          END DO
        END DO
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO j=min5,max5,-1
          DO i=npt-1,isd,-1
            v_ad(i, j) = v_ad(i, j) + 0.5*vtmp_ad(i, j)
            v_ad(i+1, j) = v_ad(i+1, j) + 0.5*vtmp_ad(i, j)
            vtmp_ad(i, j) = 0.0_8
            u_ad(i, j) = u_ad(i, j) + 0.5*utmp_ad(i, j)
            u_ad(i, j+1) = u_ad(i, j+1) + 0.5*utmp_ad(i, j)
            utmp_ad(i, j) = 0.0_8
          END DO
        END DO
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO j=jed,npy-npt+1,-1
          DO i=ied,isd,-1
            v_ad(i, j) = v_ad(i, j) + 0.5*vtmp_ad(i, j)
            v_ad(i+1, j) = v_ad(i+1, j) + 0.5*vtmp_ad(i, j)
            vtmp_ad(i, j) = 0.0_8
            u_ad(i, j) = u_ad(i, j) + 0.5*utmp_ad(i, j)
            u_ad(i, j+1) = u_ad(i, j+1) + 0.5*utmp_ad(i, j)
            utmp_ad(i, j) = 0.0_8
          END DO
        END DO
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO j=npt-1,jsd,-1
          DO i=ied,isd,-1
            v_ad(i, j) = v_ad(i, j) + 0.5*vtmp_ad(i, j)
            v_ad(i+1, j) = v_ad(i+1, j) + 0.5*vtmp_ad(i, j)
            vtmp_ad(i, j) = 0.0_8
            u_ad(i, j) = u_ad(i, j) + 0.5*utmp_ad(i, j)
            u_ad(i, j+1) = u_ad(i, j+1) + 0.5*utmp_ad(i, j)
            utmp_ad(i, j) = 0.0_8
          END DO
        END DO
      END IF
    END IF
    DO j=min3,max3,-1
      CALL POPINTEGER4(ad_from0)
      CALL POPINTEGER4(ad_to0)
      DO i=ad_to0,ad_from0,-1
        v_ad(i-1, j) = v_ad(i-1, j) + a2*vtmp_ad(i, j)
        v_ad(i+2, j) = v_ad(i+2, j) + a2*vtmp_ad(i, j)
        v_ad(i, j) = v_ad(i, j) + a1*vtmp_ad(i, j)
        v_ad(i+1, j) = v_ad(i+1, j) + a1*vtmp_ad(i, j)
        vtmp_ad(i, j) = 0.0_8
      END DO
    END DO
    DO j=min1,max1,-1
      CALL POPINTEGER4(ad_from)
      CALL POPINTEGER4(ad_to)
      DO i=ad_to,ad_from,-1
        u_ad(i, j-1) = u_ad(i, j-1) + a2*utmp_ad(i, j)
        u_ad(i, j+2) = u_ad(i, j+2) + a2*utmp_ad(i, j)
        u_ad(i, j) = u_ad(i, j) + a1*utmp_ad(i, j)
        u_ad(i, j+1) = u_ad(i, j+1) + a1*utmp_ad(i, j)
        utmp_ad(i, j) = 0.0_8
      END DO
    END DO
  END SUBROUTINE D2A2C_ADM
!#ifdef REL_VOR_DMP
! subroutine del6_flux( nord, npx, npy, damp, q, u, v)
!! Del-nord damping for the relative vorticity
!!------------------
!! nord = 0:   del-2
!! nord = 1:   del-4
!! nord = 2:   del-6
!!------------------
!   integer, intent(in):: nord            ! del-n
!   integer, intent(in):: npx, npy
!   real, intent(in):: damp
!   real, intent(in):: q(isd:ied, jsd:jed)  ! q ghosted on input
!   real, intent(inout),  dimension(isd:ied,  jsd:jed+1):: u
!   real, intent(inout),  dimension(isd:ied+1,jsd:jed  ):: v
!! local:
!   real fx2(isd:ied+1,jsd:jed), fy2(isd:ied,jsd:jed+1)
!   real d2(isd:ied,jsd:jed)
!   integer i,j, n, nt
!
!
!   do j=jsd,jed
!      do i=isd,ied
!         d2(i,j) = damp*q(i,j)
!      enddo
!   enddo
!
!   if( nord>0 ) call copy_corners(d2, npx, npy, 1)
!   do j=js-nord,je+nord
!      do i=is-nord,ie+nord+1
!         fx2(i,j) = dy(i,j)*sina_u(i,j)*(d2(i-1,j)-d2(i,j))*rdxc(i,j)
!      enddo
!   enddo
!
!   if( nord>0 ) call copy_corners(d2, npx, npy, 2)
!   do j=js-nord,je+nord+1
!      do i=is-nord,ie+nord
!         fy2(i,j) = dx(i,j)*sina_v(i,j)*(d2(i,j-1)-d2(i,j))*rdyc(i,j)
!      enddo
!   enddo
!
!   if ( nord>0 ) then
!
!!----------
!! high-order
!!----------
!
!   do n=1, nord
!
!      nt = nord-n
!
!      do j=js-nt-1,je+nt+1
!         do i=is-nt-1,ie+nt+1
!            d2(i,j) = (fx2(i,j)-fx2(i+1,j)+fy2(i,j)-fy2(i,j+1))*rarea(i,j)
!         enddo
!      enddo
!
!      call copy_corners(d2, npx, npy, 1)
!      do j=js-nt,je+nt
!         do i=is-nt,ie+nt+1
!            fx2(i,j) = dy(i,j)*sina_u(i,j)*(d2(i,j)-d2(i-1,j))*rdxc(i,j)
!         enddo
!      enddo
!
!      call copy_corners(d2, npx, npy, 2)
!      do j=js-nt,je+nt+1
!         do i=is-nt,ie+nt
!            fy2(i,j) = dx(i,j)*sina_v(i,j)*(d2(i,j)-d2(i,j-1))*rdyc(i,j)
!         enddo
!      enddo
!   enddo
!
!   endif
!
!   do j=js,je
!      do i=is,ie+1
!         u(i,j) = u(i,j) + fx2(i,j)
!      enddo
!   enddo
!
!   do j=js,je+1
!      do i=is,ie
!         v(i,j) = v(i,j) - fy2(i,j)
!      enddo
!   enddo
!
! end subroutine del6_flux
!#endif
  SUBROUTINE D2A2C(u, v, um, vm, ua, va, uc, vc, dord4)
    IMPLICIT NONE
    REAL, DIMENSION(isd:ied, jsd:jed + 1), INTENT(IN) :: u, um
    REAL, DIMENSION(isd:ied + 1, jsd:jed), INTENT(IN) :: v, vm
    LOGICAL, INTENT(IN) :: dord4
    REAL, DIMENSION(isd:ied + 1, jsd:jed), INTENT(OUT) :: uc
    REAL, DIMENSION(isd:ied, jsd:jed + 1), INTENT(OUT) :: vc
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(OUT) :: ua, va
! Local 
    REAL, DIMENSION(isd:ied, jsd:jed) :: utmp, vtmp
    INTEGER :: npt, i, j, ifirst, ilast, id
    INTRINSIC MAX
    INTRINSIC MIN
    INTEGER :: min9
    INTEGER :: min8
    INTEGER :: min7
    INTEGER :: min6
    INTEGER :: min5
    INTEGER :: min4
    INTEGER :: min3
    INTEGER :: min2
    INTEGER :: min1
    INTEGER :: min12
    INTEGER :: min11
    INTEGER :: min10
    INTEGER :: max9
    INTEGER :: max8
    INTEGER :: max7
    INTEGER :: max6
    INTEGER :: max5
    INTEGER :: max4
    INTEGER :: max3
    INTEGER :: max2
    INTEGER :: max1
    INTEGER :: max12
    INTEGER :: max11
    INTEGER :: max10
    IF (dord4) THEN
      id = 1
    ELSE
      id = 0
    END IF
    IF (grid_type .LT. 3) THEN
      npt = 4
    ELSE
      npt = -2
    END IF
    IF (npt .LT. js - 1) THEN
      max1 = js - 1
    ELSE
      max1 = npt
    END IF
    IF (npy - npt .GT. je + 1) THEN
      min1 = je + 1
    ELSE
      min1 = npy - npt
    END IF
!----------
! Interior:
!----------
    DO j=max1,min1
      IF (npt .LT. isd) THEN
        max2 = isd
      ELSE
        max2 = npt
      END IF
      IF (npx - npt .GT. ied) THEN
        min2 = ied
      ELSE
        min2 = npx - npt
      END IF
      DO i=max2,min2
        utmp(i, j) = a2*(u(i, j-1)+u(i, j+2)) + a1*(u(i, j)+u(i, j+1))
      END DO
    END DO
    IF (npt .LT. jsd) THEN
      max3 = jsd
    ELSE
      max3 = npt
    END IF
    IF (npy - npt .GT. jed) THEN
      min3 = jed
    ELSE
      min3 = npy - npt
    END IF
    DO j=max3,min3
      IF (npt .LT. is - 1) THEN
        max4 = is - 1
      ELSE
        max4 = npt
      END IF
      IF (npx - npt .GT. ie + 1) THEN
        min4 = ie + 1
      ELSE
        min4 = npx - npt
      END IF
      DO i=max4,min4
        vtmp(i, j) = a2*(v(i-1, j)+v(i+2, j)) + a1*(v(i, j)+v(i+1, j))
      END DO
    END DO
!----------
! edges:
!----------
    IF (grid_type .LT. 3) THEN
      IF (js .EQ. 1 .OR. jsd .LT. npt) THEN
        DO j=jsd,npt-1
          DO i=isd,ied
            utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
            vtmp(i, j) = 0.5*(v(i, j)+v(i+1, j))
          END DO
        END DO
      END IF
      IF (je + 1 .EQ. npy .OR. jed .GE. npy - npt) THEN
        DO j=npy-npt+1,jed
          DO i=isd,ied
            utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
            vtmp(i, j) = 0.5*(v(i, j)+v(i+1, j))
          END DO
        END DO
      END IF
      IF (is .EQ. 1 .OR. isd .LT. npt) THEN
        IF (npt .LT. jsd) THEN
          max5 = jsd
        ELSE
          max5 = npt
        END IF
        IF (npy - npt .GT. jed) THEN
          min5 = jed
        ELSE
          min5 = npy - npt
        END IF
        DO j=max5,min5
          DO i=isd,npt-1
            utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
            vtmp(i, j) = 0.5*(v(i, j)+v(i+1, j))
          END DO
        END DO
      END IF
      IF (ie + 1 .EQ. npx .OR. ied .GE. npx - npt) THEN
        IF (npt .LT. jsd) THEN
          max6 = jsd
        ELSE
          max6 = npt
        END IF
        IF (npy - npt .GT. jed) THEN
          min6 = jed
        ELSE
          min6 = npy - npt
        END IF
        DO j=max6,min6
          DO i=npx-npt+1,ied
            utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
            vtmp(i, j) = 0.5*(v(i, j)+v(i+1, j))
          END DO
        END DO
      END IF
    END IF
    DO j=js-1-id,je+1+id
      DO i=is-1-id,ie+1+id
        ua(i, j) = (utmp(i, j)-vtmp(i, j)*cosa_s(i, j))*rsin2(i, j)
        va(i, j) = (vtmp(i, j)-utmp(i, j)*cosa_s(i, j))*rsin2(i, j)
      END DO
    END DO
    IF (npt .LT. js - 1) THEN
      max7 = js - 1
    ELSE
      max7 = npt
    END IF
    IF (npy - npt .GT. je + 1) THEN
      min7 = je + 1
    ELSE
      min7 = npy - npt
    END IF
! Re-compute (utmp, vtmp) using (um,vm)
!----------
! Interior:
!----------
    DO j=max7,min7
      IF (npt .LT. isd) THEN
        max8 = isd
      ELSE
        max8 = npt
      END IF
      IF (npx - npt .GT. ied) THEN
        min8 = ied
      ELSE
        min8 = npx - npt
      END IF
      DO i=max8,min8
        utmp(i, j) = a2*(um(i, j-1)+um(i, j+2)) + a1*(um(i, j)+um(i, j+1&
&         ))
      END DO
    END DO
    IF (npt .LT. jsd) THEN
      max9 = jsd
    ELSE
      max9 = npt
    END IF
    IF (npy - npt .GT. jed) THEN
      min9 = jed
    ELSE
      min9 = npy - npt
    END IF
    DO j=max9,min9
      IF (npt .LT. is - 1) THEN
        max10 = is - 1
      ELSE
        max10 = npt
      END IF
      IF (npx - npt .GT. ie + 1) THEN
        min10 = ie + 1
      ELSE
        min10 = npx - npt
      END IF
      DO i=max10,min10
        vtmp(i, j) = a2*(vm(i-1, j)+vm(i+2, j)) + a1*(vm(i, j)+vm(i+1, j&
&         ))
      END DO
    END DO
!----------
! edges:
!----------
    IF (grid_type .LT. 3) THEN
      IF (js .EQ. 1 .OR. jsd .LT. npt) THEN
        DO j=jsd,npt-1
          DO i=isd,ied
            utmp(i, j) = 0.5*(um(i, j)+um(i, j+1))
            vtmp(i, j) = 0.5*(vm(i, j)+vm(i+1, j))
          END DO
        END DO
      END IF
      IF (je + 1 .EQ. npy .OR. jed .GE. npy - npt) THEN
        DO j=npy-npt+1,jed
          DO i=isd,ied
            utmp(i, j) = 0.5*(um(i, j)+um(i, j+1))
            vtmp(i, j) = 0.5*(vm(i, j)+vm(i+1, j))
          END DO
        END DO
      END IF
      IF (is .EQ. 1 .OR. isd .LT. npt) THEN
        IF (npt .LT. jsd) THEN
          max11 = jsd
        ELSE
          max11 = npt
        END IF
        IF (npy - npt .GT. jed) THEN
          min11 = jed
        ELSE
          min11 = npy - npt
        END IF
        DO j=max11,min11
          DO i=isd,npt-1
            utmp(i, j) = 0.5*(um(i, j)+um(i, j+1))
            vtmp(i, j) = 0.5*(vm(i, j)+vm(i+1, j))
          END DO
        END DO
      END IF
      IF (ie + 1 .EQ. npx .OR. ied .GE. npx - npt) THEN
        IF (npt .LT. jsd) THEN
          max12 = jsd
        ELSE
          max12 = npt
        END IF
        IF (npy - npt .GT. jed) THEN
          min12 = jed
        ELSE
          min12 = npy - npt
        END IF
        DO j=max12,min12
          DO i=npx-npt+1,ied
            utmp(i, j) = 0.5*(um(i, j)+um(i, j+1))
            vtmp(i, j) = 0.5*(vm(i, j)+vm(i+1, j))
          END DO
        END DO
      END IF
    END IF
! A -> C
!--------------
! Fix the edges
!--------------
! Xdir:
    IF (sw_corner) THEN
      DO i=-2,0
        utmp(i, 0) = -vtmp(0, 1-i)
      END DO
    END IF
    IF (se_corner) THEN
      DO i=0,2
        utmp(npx+i, 0) = vtmp(npx, i+1)
      END DO
    END IF
    IF (ne_corner) THEN
      DO i=0,2
        utmp(npx+i, npy) = -vtmp(npx, je-i)
      END DO
    END IF
    IF (nw_corner) THEN
      DO i=-2,0
        utmp(i, npy) = vtmp(0, je+i)
      END DO
    END IF
    IF (grid_type .LT. 3) THEN
      IF (3 .LT. is) THEN
        ifirst = is
      ELSE
        ifirst = 3
      END IF
      IF (npx - 2 .GT. ie + 1) THEN
        ilast = ie + 1
      ELSE
        ilast = npx - 2
      END IF
    ELSE
      ifirst = is
      ilast = ie + 1
    END IF
!---------------------------------------------
! 4th order interpolation for interior points:
!---------------------------------------------
    DO j=js,je
      DO i=ifirst,ilast
        uc(i, j) = a1*(utmp(i-1, j)+utmp(i, j)) + a2*(utmp(i-2, j)+utmp(&
&         i+1, j))
      END DO
    END DO
    IF (grid_type .LT. 3) THEN
      IF (is .EQ. 1) THEN
        DO j=js,je
! 3-pt extrapolation --------------------------------------------------
          uc(1, j) = (t14*(utmp(0, j)+utmp(1, j))+t12*(utmp(-1, j)+utmp(&
&           2, j))+t15*(utmp(-2, j)+utmp(3, j)))*rsin_u(1, j)
          uc(2, j) = c1*utmp(3, j) + c2*utmp(2, j) + c3*utmp(1, j)
        END DO
      END IF
      IF (ie + 1 .EQ. npx) THEN
        DO j=js,je
          uc(npx-1, j) = c1*utmp(npx-3, j) + c2*utmp(npx-2, j) + c3*utmp&
&           (npx-1, j)
! 3-pt extrapolation --------------------------------------------------------
          uc(npx, j) = (t14*(utmp(npx-1, j)+utmp(npx, j))+t12*(utmp(npx-&
&           2, j)+utmp(npx+1, j))+t15*(utmp(npx-3, j)+utmp(npx+2, j)))*&
&           rsin_u(npx, j)
        END DO
      END IF
    END IF
!------
! Ydir:
!------
    IF (sw_corner) THEN
      DO j=-2,0
        vtmp(0, j) = -utmp(1-j, 0)
      END DO
    END IF
    IF (nw_corner) THEN
      DO j=0,2
        vtmp(0, npy+j) = utmp(j+1, npy)
      END DO
    END IF
    IF (se_corner) THEN
      DO j=-2,0
        vtmp(npx, j) = utmp(ie+j, 0)
      END DO
    END IF
    IF (ne_corner) THEN
      DO j=0,2
        vtmp(npx, npy+j) = -utmp(ie-j, npy)
      END DO
    END IF
    IF (grid_type .LT. 3) THEN
      DO j=js,je+1
        IF (j .EQ. 1) THEN
          DO i=is,ie
! 3-pt extrapolation -----------------------------------------
            vc(i, 1) = (t14*(vtmp(i, 0)+vtmp(i, 1))+t12*(vtmp(i, -1)+&
&             vtmp(i, 2))+t15*(vtmp(i, -2)+vtmp(i, 3)))*rsin_v(i, 1)
          END DO
        ELSE IF (j .EQ. npy - 1) THEN
          DO i=is,ie
            vc(i, j) = c1*vtmp(i, j-2) + c2*vtmp(i, j-1) + c3*vtmp(i, j)
          END DO
        ELSE IF (j .EQ. 2) THEN
          DO i=is,ie
!          vc(i,j) = c1*vtmp(i,j+1) + c2*vtmp(i,j) + c3*vtmp(i,j-1)
            vc(i, j) = c3*vtmp(i, j-1) + c2*vtmp(i, j) + c1*vtmp(i, j+1)
          END DO
        ELSE IF (j .EQ. npy) THEN
          DO i=is,ie
! 3-pt extrapolation --------------------------------------------------------
            vc(i, npy) = (t14*(vtmp(i, npy-1)+vtmp(i, npy))+t12*(vtmp(i&
&             , npy-2)+vtmp(i, npy+1))+t15*(vtmp(i, npy-3)+vtmp(i, npy+2&
&             )))*rsin_v(i, npy)
          END DO
        ELSE
! 4th order interpolation for interior points:
          DO i=is,ie
            vc(i, j) = a2*(vtmp(i, j-2)+vtmp(i, j+1)) + a1*(vtmp(i, j-1)&
&             +vtmp(i, j))
          END DO
        END IF
      END DO
    ELSE
! 4th order interpolation:
      DO j=js,je+1
        DO i=is,ie
          vc(i, j) = a2*(vtmp(i, j-2)+vtmp(i, j+1)) + a1*(vtmp(i, j-1)+&
&           vtmp(i, j))
        END DO
      END DO
    END IF
  END SUBROUTINE D2A2C

 end module sw_core_adm_mod
