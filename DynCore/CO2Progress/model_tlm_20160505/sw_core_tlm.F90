 module sw_core_tlm_mod

 use fv_arrays_mod,     only: p_precision, i_precision, ke_precision
 use fv_control_mod,    only: shallow_water
 use fv_mp_mod,         only: ng, is,js,ie,je, isd,jsd,ied,jed,  &
                              mp_corner_comm, domain,            &
                              fill_corners, XDir, YDir
 use fv_grid_tools_mod, only: npx=>npx_g,npy=>npy_g, cosa, sina,  &
                              rdxc, rdyc, dx,dy, dxc,dyc, dxa,dya,  &
                              rdxa, rdya, area, area_c, rarea, rarea_c, rdx, rdy, &
                              grid_type
 use tp_core_mod,       only: fv_tp_2d
 use sw_core_mod,       only: xtp_u, ytp_v
 use tp_core_tlm_mod,   only: fv_tp_2d_tlm, copy_corners_tlm
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
  public :: c_sw_tlm, d_sw_tlm, d2a2c_tlm, d2a2c_vect_tlm, divergence_corner_tlm
  contains

!  Differentiation of c_sw in forward (tangent) mode (with options r8):
!   variations   of useful results: w delp ua uc ptc ut delpc va
!                vc vt pt
!   with respect to varying inputs: u v w delp ua uc ptc ut delpc
!                va vc vt pt
  SUBROUTINE C_SW_TLM(delpc, delpc_tl, delp, delp_tl, ptc, ptc_tl, pt, &
&   pt_tl, u, u_tl, v, v_tl, w, w_tl, uc, uc_tl, vc, vc_tl, ua, ua_tl, &
&   va, va_tl, wc, ut, ut_tl, vt, vt_tl, dt2, hydrostatic, dord4)
    IMPLICIT NONE
!#endif
    REAL, DIMENSION(isd:ied, jsd:jed + 1), INTENT(INOUT) :: u, vc
    REAL, DIMENSION(isd:ied, jsd:jed+1), INTENT(INOUT) :: u_tl, vc_tl
    REAL, DIMENSION(isd:ied + 1, jsd:jed), INTENT(INOUT) :: v, uc
    REAL, DIMENSION(isd:ied+1, jsd:jed), INTENT(INOUT) :: v_tl, uc_tl
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(INOUT) :: delp, pt, ua, va&
&   , w
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(INOUT) :: delp_tl, pt_tl, &
&   ua_tl, va_tl, w_tl
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(OUT) :: delpc, ptc, ut, vt&
&   , wc
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(OUT) :: delpc_tl, ptc_tl, &
&   ut_tl, vt_tl
    REAL, INTENT(IN) :: dt2
    LOGICAL, INTENT(IN) :: hydrostatic
    LOGICAL, INTENT(IN) :: dord4
! Local:
    REAL, DIMENSION(is - 1:ie + 1, js - 1:je + 1) :: vort, ke
    REAL, DIMENSION(is-1:ie+1, js-1:je+1) :: vort_tl, ke_tl
    REAL, DIMENSION(is - 1:ie + 2, js - 1:je + 1) :: fx, fx1, fx2
    REAL, DIMENSION(is-1:ie+2, js-1:je+1) :: fx_tl, fx1_tl
    REAL, DIMENSION(is - 1:ie + 1, js - 1:je + 2) :: fy, fy1, fy2
    REAL, DIMENSION(is-1:ie+1, js-1:je+2) :: fy_tl, fy1_tl
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
    CALL D2A2C_VECT_TLM(u, u_tl, v, v_tl, ua, ua_tl, va, va_tl, uc, &
&                 uc_tl, vc, vc_tl, ut, ut_tl, vt, vt_tl, dord4)
!     call d2a2c_vect_v2(u, v, ua, va, uc, vc, ut, vt)
    DO j=js-1,jep1
      DO i=is-1,iep1+1
        ut_tl(i, j) = dt2*dy(i, j)*sina_u(i, j)*ut_tl(i, j)
        ut(i, j) = dt2*ut(i, j)*dy(i, j)*sina_u(i, j)
      END DO
    END DO
    DO j=js-1,jep1+1
      DO i=is-1,iep1
        vt_tl(i, j) = dt2*dx(i, j)*sina_v(i, j)*vt_tl(i, j)
        vt(i, j) = dt2*vt(i, j)*dx(i, j)*sina_v(i, j)
      END DO
    END DO
!----------------
! Transport delp:
!----------------
! Xdir:
    IF (grid_type .LT. 3) CALL FILL2_4CORNERS_TLM(delp, delp_tl, pt, &
&                                           pt_tl, 1)
    IF (hydrostatic) THEN
      IF (shallow_water) THEN
        fx1_tl = 0.0_8
        DO j=js-1,jep1
          DO i=is-1,iep1+1
            IF (ut(i, j) .GT. 0.) THEN
              fx1_tl(i, j) = delp_tl(i-1, j)
              fx1(i, j) = delp(i-1, j)
            ELSE
              fx1_tl(i, j) = delp_tl(i, j)
              fx1(i, j) = delp(i, j)
            END IF
            fx1_tl(i, j) = ut_tl(i, j)*fx1(i, j) + ut(i, j)*fx1_tl(i, j)
            fx1(i, j) = ut(i, j)*fx1(i, j)
          END DO
        END DO
        fx_tl = 0.0_8
      ELSE
        fx_tl = 0.0_8
        fx1_tl = 0.0_8
        DO j=js-1,jep1
          DO i=is-1,iep1+1
            IF (ut(i, j) .GT. 0.) THEN
              fx1_tl(i, j) = delp_tl(i-1, j)
              fx1(i, j) = delp(i-1, j)
              fx_tl(i, j) = pt_tl(i-1, j)
              fx(i, j) = pt(i-1, j)
            ELSE
              fx1_tl(i, j) = delp_tl(i, j)
              fx1(i, j) = delp(i, j)
              fx_tl(i, j) = pt_tl(i, j)
              fx(i, j) = pt(i, j)
            END IF
            fx1_tl(i, j) = ut_tl(i, j)*fx1(i, j) + ut(i, j)*fx1_tl(i, j)
            fx1(i, j) = ut(i, j)*fx1(i, j)
            fx_tl(i, j) = fx1_tl(i, j)*fx(i, j) + fx1(i, j)*fx_tl(i, j)
            fx(i, j) = fx1(i, j)*fx(i, j)
          END DO
        END DO
      END IF
    ELSE
      IF (grid_type .LT. 3) THEN
        CALL FILL_4CORNERS_TLM(w, w_tl, 1)
        fx_tl = 0.0_8
        fx1_tl = 0.0_8
      ELSE
        fx_tl = 0.0_8
        fx1_tl = 0.0_8
      END IF
      DO j=js-1,je+1
        DO i=is-1,ie+2
          IF (ut(i, j) .GT. 0.) THEN
            fx1_tl(i, j) = delp_tl(i-1, j)
            fx1(i, j) = delp(i-1, j)
            fx_tl(i, j) = pt_tl(i-1, j)
            fx(i, j) = pt(i-1, j)
            fx2(i, j) = w(i-1, j)
          ELSE
            fx1_tl(i, j) = delp_tl(i, j)
            fx1(i, j) = delp(i, j)
            fx_tl(i, j) = pt_tl(i, j)
            fx(i, j) = pt(i, j)
            fx2(i, j) = w(i, j)
          END IF
          fx1_tl(i, j) = ut_tl(i, j)*fx1(i, j) + ut(i, j)*fx1_tl(i, j)
          fx1(i, j) = ut(i, j)*fx1(i, j)
          fx_tl(i, j) = fx1_tl(i, j)*fx(i, j) + fx1(i, j)*fx_tl(i, j)
          fx(i, j) = fx1(i, j)*fx(i, j)
          fx2(i, j) = fx1(i, j)*fx2(i, j)
        END DO
      END DO
    END IF
! Ydir:
    IF (grid_type .LT. 3) CALL FILL2_4CORNERS_TLM(delp, delp_tl, pt, &
&                                           pt_tl, 2)
    IF (hydrostatic) THEN
      fy1_tl = 0.0_8
      fy_tl = 0.0_8
      DO j=js-1,jep1+1
        DO i=is-1,iep1
          IF (vt(i, j) .GT. 0.) THEN
            fy1_tl(i, j) = delp_tl(i, j-1)
            fy1(i, j) = delp(i, j-1)
            fy_tl(i, j) = pt_tl(i, j-1)
            fy(i, j) = pt(i, j-1)
          ELSE
            fy1_tl(i, j) = delp_tl(i, j)
            fy1(i, j) = delp(i, j)
            fy_tl(i, j) = pt_tl(i, j)
            fy(i, j) = pt(i, j)
          END IF
          fy1_tl(i, j) = vt_tl(i, j)*fy1(i, j) + vt(i, j)*fy1_tl(i, j)
          fy1(i, j) = vt(i, j)*fy1(i, j)
          fy_tl(i, j) = fy1_tl(i, j)*fy(i, j) + fy1(i, j)*fy_tl(i, j)
          fy(i, j) = fy1(i, j)*fy(i, j)
        END DO
      END DO
      IF (shallow_water) THEN
        DO j=js-1,jep1
          DO i=is-1,iep1
            delpc_tl(i, j) = delp_tl(i, j) + rarea(i, j)*(fx1_tl(i, j)-&
&             fx1_tl(i+1, j)+fy1_tl(i, j)-fy1_tl(i, j+1))
            delpc(i, j) = delp(i, j) + (fx1(i, j)-fx1(i+1, j)+fy1(i, j)-&
&             fy1(i, j+1))*rarea(i, j)
            ptc_tl(i, j) = pt_tl(i, j)
            ptc(i, j) = pt(i, j)
          END DO
        END DO
        ke_tl = 0.0_8
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
            delpc_tl(i, j) = delp_tl(i, j) + rarea(i, j)*(fx1_tl(i, j)-&
&             fx1_tl(i+1, j)+fy1_tl(i, j)-fy1_tl(i, j+1))
            delpc(i, j) = delp(i, j) + (fx1(i, j)-fx1(i+1, j)+fy1(i, j)-&
&             fy1(i, j+1))*rarea(i, j)
            ptc_tl(i, j) = ((pt_tl(i, j)*delp(i, j)+pt(i, j)*delp_tl(i, &
&             j)+rarea(i, j)*(fx_tl(i, j)-fx_tl(i+1, j)+fy_tl(i, j)-&
&             fy_tl(i, j+1)))*delpc(i, j)-(pt(i, j)*delp(i, j)+(fx(i, j)&
&             -fx(i+1, j)+fy(i, j)-fy(i, j+1))*rarea(i, j))*delpc_tl(i, &
&             j))/delpc(i, j)**2
            ptc(i, j) = (pt(i, j)*delp(i, j)+(fx(i, j)-fx(i+1, j)+fy(i, &
&             j)-fy(i, j+1))*rarea(i, j))/delpc(i, j)
          END DO
        END DO
        ke_tl = 0.0_8
      END IF
    ELSE
!#endif
      IF (grid_type .LT. 3) THEN
        CALL FILL_4CORNERS_TLM(w, w_tl, 2)
        fy1_tl = 0.0_8
        fy_tl = 0.0_8
      ELSE
        fy1_tl = 0.0_8
        fy_tl = 0.0_8
      END IF
      DO j=js-1,je+2
        DO i=is-1,ie+1
          IF (vt(i, j) .GT. 0.) THEN
            fy1_tl(i, j) = delp_tl(i, j-1)
            fy1(i, j) = delp(i, j-1)
            fy_tl(i, j) = pt_tl(i, j-1)
            fy(i, j) = pt(i, j-1)
            fy2(i, j) = w(i, j-1)
          ELSE
            fy1_tl(i, j) = delp_tl(i, j)
            fy1(i, j) = delp(i, j)
            fy_tl(i, j) = pt_tl(i, j)
            fy(i, j) = pt(i, j)
            fy2(i, j) = w(i, j)
          END IF
          fy1_tl(i, j) = vt_tl(i, j)*fy1(i, j) + vt(i, j)*fy1_tl(i, j)
          fy1(i, j) = vt(i, j)*fy1(i, j)
          fy_tl(i, j) = fy1_tl(i, j)*fy(i, j) + fy1(i, j)*fy_tl(i, j)
          fy(i, j) = fy1(i, j)*fy(i, j)
          fy2(i, j) = fy1(i, j)*fy2(i, j)
        END DO
      END DO
      DO j=js-1,je+1
        DO i=is-1,ie+1
          delpc_tl(i, j) = delp_tl(i, j) + rarea(i, j)*(fx1_tl(i, j)-&
&           fx1_tl(i+1, j)+fy1_tl(i, j)-fy1_tl(i, j+1))
          delpc(i, j) = delp(i, j) + (fx1(i, j)-fx1(i+1, j)+fy1(i, j)-&
&           fy1(i, j+1))*rarea(i, j)
          ptc_tl(i, j) = ((pt_tl(i, j)*delp(i, j)+pt(i, j)*delp_tl(i, j)&
&           +rarea(i, j)*(fx_tl(i, j)-fx_tl(i+1, j)+fy_tl(i, j)-fy_tl(i&
&           , j+1)))*delpc(i, j)-(pt(i, j)*delp(i, j)+(fx(i, j)-fx(i+1, &
&           j)+fy(i, j)-fy(i, j+1))*rarea(i, j))*delpc_tl(i, j))/delpc(i&
&           , j)**2
          ptc(i, j) = (pt(i, j)*delp(i, j)+(fx(i, j)-fx(i+1, j)+fy(i, j)&
&           -fy(i, j+1))*rarea(i, j))/delpc(i, j)
          wc(i, j) = (w(i, j)*delp(i, j)+(fx2(i, j)-fx2(i+1, j)+fy2(i, j&
&           )-fy2(i, j+1))*rarea(i, j))/delpc(i, j)
        END DO
      END DO
      ke_tl = 0.0_8
    END IF
!------------
! Compute KE:
!------------
    DO j=js-1,jep1
      DO i=is-1,iep1
        IF (ua(i, j) .GT. 0.) THEN
          IF (i .EQ. 1) THEN
            ke_tl(1, j) = sina_u(1, j)*uc_tl(1, j) + cosa_u(1, j)*v_tl(1&
&             , j)
            ke(1, j) = uc(1, j)*sina_u(1, j) + v(1, j)*cosa_u(1, j)
          ELSE IF (i .EQ. npx) THEN
            ke_tl(i, j) = sina_u(npx, j)*uc_tl(npx, j) - cosa_u(npx, j)*&
&             v_tl(npx, j)
            ke(i, j) = uc(npx, j)*sina_u(npx, j) - v(npx, j)*cosa_u(npx&
&             , j)
          ELSE
            ke_tl(i, j) = uc_tl(i, j)
            ke(i, j) = uc(i, j)
          END IF
        ELSE IF (i .EQ. 0) THEN
          ke_tl(0, j) = sina_u(1, j)*uc_tl(1, j) - cosa_u(1, j)*v_tl(1, &
&           j)
          ke(0, j) = uc(1, j)*sina_u(1, j) - v(1, j)*cosa_u(1, j)
        ELSE IF (i .EQ. npx - 1) THEN
          ke_tl(i, j) = sina_u(npx, j)*uc_tl(npx, j) + cosa_u(npx, j)*&
&           v_tl(npx, j)
          ke(i, j) = uc(npx, j)*sina_u(npx, j) + v(npx, j)*cosa_u(npx, j&
&           )
        ELSE
          ke_tl(i, j) = uc_tl(i+1, j)
          ke(i, j) = uc(i+1, j)
        END IF
      END DO
    END DO
    vort_tl = 0.0_8
    DO j=js-1,jep1
      DO i=is-1,iep1
        IF (va(i, j) .GT. 0.) THEN
          IF (j .EQ. 1) THEN
            vort_tl(i, 1) = sina_v(i, 1)*vc_tl(i, 1) + cosa_v(i, 1)*u_tl&
&             (i, 1)
            vort(i, 1) = vc(i, 1)*sina_v(i, 1) + u(i, 1)*cosa_v(i, 1)
          ELSE IF (j .EQ. npy) THEN
            vort_tl(i, j) = sina_v(i, npy)*vc_tl(i, npy) - cosa_v(i, npy&
&             )*u_tl(i, npy)
            vort(i, j) = vc(i, npy)*sina_v(i, npy) - u(i, npy)*cosa_v(i&
&             , npy)
          ELSE
            vort_tl(i, j) = vc_tl(i, j)
            vort(i, j) = vc(i, j)
          END IF
        ELSE IF (j .EQ. 0) THEN
          vort_tl(i, 0) = sina_v(i, 1)*vc_tl(i, 1) - cosa_v(i, 1)*u_tl(i&
&           , 1)
          vort(i, 0) = vc(i, 1)*sina_v(i, 1) - u(i, 1)*cosa_v(i, 1)
        ELSE IF (j .EQ. npy - 1) THEN
          vort_tl(i, j) = sina_v(i, npy)*vc_tl(i, npy) + cosa_v(i, npy)*&
&           u_tl(i, npy)
          vort(i, j) = vc(i, npy)*sina_v(i, npy) + u(i, npy)*cosa_v(i, &
&           npy)
        ELSE
          vort_tl(i, j) = vc_tl(i, j+1)
          vort(i, j) = vc(i, j+1)
        END IF
      END DO
    END DO
    dt4 = 0.5*dt2
    DO j=js-1,jep1
      DO i=is-1,iep1
        ke_tl(i, j) = dt4*(ua_tl(i, j)*ke(i, j)+ua(i, j)*ke_tl(i, j)+&
&         va_tl(i, j)*vort(i, j)+va(i, j)*vort_tl(i, j))
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
        fx_tl(i, j) = dxc(i, j)*uc_tl(i, j)
        fx(i, j) = uc(i, j)*dxc(i, j)
      END DO
      IF (is .EQ. 1) THEN
        fx_tl(1, j) = sina_u(1, j)*dxc(1, j)*uc_tl(1, j)
        fx(1, j) = uc(1, j)*sina_u(1, j)*dxc(1, j)
      END IF
      IF (ie + 1 .EQ. npx) THEN
        fx_tl(npx, j) = sina_u(npx, j)*dxc(npx, j)*uc_tl(npx, j)
        fx(npx, j) = uc(npx, j)*sina_u(npx, j)*dxc(npx, j)
      END IF
    END DO
    DO j=js,je+1
      IF (j .EQ. 1 .OR. j .EQ. npy) THEN
        DO i=is-1,ie+1
          fy_tl(i, j) = sina_v(i, j)*dyc(i, j)*vc_tl(i, j)
          fy(i, j) = vc(i, j)*sina_v(i, j)*dyc(i, j)
        END DO
      ELSE
        DO i=is-1,ie+1
          fy_tl(i, j) = dyc(i, j)*vc_tl(i, j)
          fy(i, j) = vc(i, j)*dyc(i, j)
        END DO
      END IF
    END DO
!#endif
    DO j=js,je+1
      DO i=is,ie+1
        vort_tl(i, j) = fx_tl(i, j-1) - fx_tl(i, j) - fy_tl(i-1, j) + &
&         fy_tl(i, j)
        vort(i, j) = fx(i, j-1) - fx(i, j) - fy(i-1, j) + fy(i, j)
      END DO
    END DO
! Remove the extra term at the corners:
    IF (sw_corner) THEN
      vort_tl(1, 1) = vort_tl(1, 1) + fy_tl(0, 1)
      vort(1, 1) = vort(1, 1) + fy(0, 1)
    END IF
    IF (se_corner) THEN
      vort_tl(npx, 1) = vort_tl(npx, 1) - fy_tl(npx, 1)
      vort(npx, 1) = vort(npx, 1) - fy(npx, 1)
    END IF
    IF (ne_corner) THEN
      vort_tl(npx, npy) = vort_tl(npx, npy) - fy_tl(npx, npy)
      vort(npx, npy) = vort(npx, npy) - fy(npx, npy)
    END IF
    IF (nw_corner) THEN
      vort_tl(1, npy) = vort_tl(1, npy) + fy_tl(0, npy)
      vort(1, npy) = vort(1, npy) + fy(0, npy)
    END IF
!----------------------------
! Compute absolute vorticity
!----------------------------
    DO j=js,je+1
      DO i=is,ie+1
        vort_tl(i, j) = rarea_c(i, j)*vort_tl(i, j)
        vort(i, j) = fc(i, j) + rarea_c(i, j)*vort(i, j)
      END DO
    END DO
!----------------------------------
! Transport absolute vorticity:
!----------------------------------
    DO j=js,je
      DO i=is,iep1
        IF (i .EQ. 1 .OR. i .EQ. npx) THEN
          fy1_tl(i, j) = dt2*sina_u(i, j)*v_tl(i, j)
          fy1(i, j) = dt2*v(i, j)*sina_u(i, j)
        ELSE
          fy1_tl(i, j) = dt2*(v_tl(i, j)-cosa_u(i, j)*uc_tl(i, j))/&
&           sina_u(i, j)
          fy1(i, j) = dt2*(v(i, j)-uc(i, j)*cosa_u(i, j))/sina_u(i, j)
        END IF
        IF (fy1(i, j) .GT. 0.) THEN
          fy_tl(i, j) = vort_tl(i, j)
          fy(i, j) = vort(i, j)
        ELSE
          fy_tl(i, j) = vort_tl(i, j+1)
          fy(i, j) = vort(i, j+1)
        END IF
      END DO
    END DO
    DO j=js,jep1
      IF (j .EQ. 1 .OR. j .EQ. npy) THEN
        DO i=is,ie
          fx1_tl(i, j) = dt2*sina_v(i, j)*u_tl(i, j)
          fx1(i, j) = dt2*u(i, j)*sina_v(i, j)
          IF (fx1(i, j) .GT. 0.) THEN
            fx_tl(i, j) = vort_tl(i, j)
            fx(i, j) = vort(i, j)
          ELSE
            fx_tl(i, j) = vort_tl(i+1, j)
            fx(i, j) = vort(i+1, j)
          END IF
        END DO
      ELSE
        DO i=is,ie
          fx1_tl(i, j) = dt2*(u_tl(i, j)-cosa_v(i, j)*vc_tl(i, j))/&
&           sina_v(i, j)
          fx1(i, j) = dt2*(u(i, j)-vc(i, j)*cosa_v(i, j))/sina_v(i, j)
          IF (fx1(i, j) .GT. 0.) THEN
            fx_tl(i, j) = vort_tl(i, j)
            fx(i, j) = vort(i, j)
          ELSE
            fx_tl(i, j) = vort_tl(i+1, j)
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
        uc_tl(i, j) = uc_tl(i, j) + fy1_tl(i, j)*fy(i, j) + fy1(i, j)*&
&         fy_tl(i, j) + rdxc(i, j)*(ke_tl(i-1, j)-ke_tl(i, j))
        uc(i, j) = uc(i, j) + fy1(i, j)*fy(i, j) + rdxc(i, j)*(ke(i-1, j&
&         )-ke(i, j))
      END DO
    END DO
    DO j=js,jep1
      DO i=is,ie
        vc_tl(i, j) = vc_tl(i, j) - fx1_tl(i, j)*fx(i, j) - fx1(i, j)*&
&         fx_tl(i, j) + rdyc(i, j)*(ke_tl(i, j-1)-ke_tl(i, j))
        vc(i, j) = vc(i, j) - fx1(i, j)*fx(i, j) + rdyc(i, j)*(ke(i, j-1&
&         )-ke(i, j))
      END DO
    END DO
  END SUBROUTINE C_SW_TLM
!  Differentiation of d_sw in forward (tangent) mode (with options r8):
!   variations   of useful results: yfx_adv q crx_adv u v w delp
!                xfx_adv uc ptc xflux cry_adv delpc vc yflux divg_d
!                pt cx cy
!   with respect to varying inputs: yfx_adv q crx_adv u v w delp
!                ua xfx_adv uc ptc xflux cry_adv delpc va vc yflux
!                pkz divg_d pt cx cy
!-------------------------------------------------------------------------------
!
!     d_sw :: D-Grid Shallow Water Routine
!
  SUBROUTINE D_SW_TLM(delpc, delpc_tl, delp, delp_tl, ptc, ptc_tl, pt, &
&   pt_tl, u, u_tl, v, v_tl, w, w_tl, uc, uc_tl, vc, vc_tl, ua, ua_tl, &
&   va, va_tl, divg_d, divg_d_tl, xflux, xflux_tl, yflux, yflux_tl, cx, &
&   cx_tl, cy, cy_tl, crx_adv, crx_adv_tl, cry_adv, cry_adv_tl, xfx_adv&
&   , xfx_adv_tl, yfx_adv, yfx_adv_tl, zvir, sphum, nq, q, q_tl, k, km, &
&   inline_q, pkz, pkz_tl, dt, hord_tr, hord_mt, hord_vt, hord_tm, &
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
    REAL(p_precision), INTENT(IN) :: pkz_tl(is:ie, js:je)
! divergence
    REAL, INTENT(INOUT) :: divg_d(isd:ied+1, jsd:jed+1)
    REAL, INTENT(INOUT) :: divg_d_tl(isd:ied+1, jsd:jed+1)
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(INOUT) :: delp, pt, ua, va&
&   , w
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(INOUT) :: delp_tl, pt_tl, &
&   ua_tl, va_tl, w_tl
    REAL, DIMENSION(isd:ied, jsd:jed) :: delp_tj, pt_tj, w_tj
    REAL                              :: q_tj(isd:ied, jsd:jed, km, nq)
    REAL, DIMENSION(isd:ied, jsd:jed + 1), INTENT(INOUT) :: u, vc
    REAL, DIMENSION(isd:ied, jsd:jed+1), INTENT(INOUT) :: u_tl, vc_tl
    REAL, DIMENSION(isd:ied + 1, jsd:jed), INTENT(INOUT) :: v, uc
    REAL, DIMENSION(isd:ied+1, jsd:jed), INTENT(INOUT) :: v_tl, uc_tl
    REAL, INTENT(INOUT) :: q(isd:ied, jsd:jed, km, nq)
    REAL, INTENT(INOUT) :: q_tl(isd:ied, jsd:jed, km, nq)
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(OUT) :: delpc, ptc
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(OUT) :: delpc_tl, ptc_tl
! The flux capacitors:
    REAL, INTENT(INOUT) :: xflux(is:ie+1, js:je)
    REAL, INTENT(INOUT) :: xflux_tl(is:ie+1, js:je)
    REAL, INTENT(INOUT) :: yflux(is:ie, js:je+1)
    REAL, INTENT(INOUT) :: yflux_tl(is:ie, js:je+1)
!------------------------
    REAL, INTENT(INOUT) :: cx(is:ie+1, jsd:jed)
    REAL, INTENT(INOUT) :: cx_tl(is:ie+1, jsd:jed)
    REAL, INTENT(INOUT) :: cy(isd:ied, js:je+1)
    REAL, INTENT(INOUT) :: cy_tl(isd:ied, js:je+1)
    LOGICAL, INTENT(IN) :: hydrostatic
    LOGICAL, INTENT(IN) :: inline_q
    REAL, DIMENSION(is:ie + 1, jsd:jed), INTENT(OUT) :: crx_adv, xfx_adv
    REAL, DIMENSION(is:ie+1, jsd:jed), INTENT(OUT) :: crx_adv_tl, &
&   xfx_adv_tl
    REAL, DIMENSION(isd:ied, js:je + 1), INTENT(OUT) :: cry_adv, yfx_adv
    REAL, DIMENSION(isd:ied, js:je+1), INTENT(OUT) :: cry_adv_tl, &
&   yfx_adv_tl
! Local:
    REAL(i_precision), DIMENSION(isd:ied + 1, jsd:jed) :: ut, uci
    REAL(i_precision), DIMENSION(isd:ied+1, jsd:jed) :: ut_tl, uci_tl
    REAL(i_precision), DIMENSION(isd:ied, jsd:jed + 1) :: vt, vci
    REAL(i_precision), DIMENSION(isd:ied, jsd:jed+1) :: vt_tl, vci_tl
    REAL(i_precision) :: damp
    REAL(i_precision) :: damp_tl
    REAL(i_precision), SAVE :: r1b4=1.0/4.0
    REAL(i_precision), SAVE :: r1b16=1.0/16.0
    REAL(ke_precision), DIMENSION(is:ie + 1, js:je + 1) :: ub, vb
    REAL(ke_precision), DIMENSION(is:ie+1, js:je+1) :: ub_tl, vb_tl
    REAL(ke_precision), DIMENSION(is:ie + 1, js:je + 1) :: ub_tj, vb_tj
!  needs this for corner_comm
    REAL(ke_precision) :: ke(isd:ied+1, jsd:jed+1)
    REAL(ke_precision) :: ke_tl(isd:ied+1, jsd:jed+1)
    REAL(ke_precision) :: dt4, dt5, dt6
!#ifdef FULL_SMAG
!      real :: vt2(is-1:ie+1,js-1:je+1)     !  work array
!#endif
!  work array
    REAL :: wk(isd:ied, jsd:jed)
    REAL :: wk_tl(isd:ied, jsd:jed)
! Vorticity
    REAL :: vort(isd:ied, jsd:jed)
    REAL :: vort_tj(isd:ied, jsd:jed)
    REAL :: vort_tl(isd:ied, jsd:jed)
! 1-D X-direction Fluxes
    REAL :: fx(is:ie+1, js:je)
    REAL :: fx_tl(is:ie+1, js:je)
    REAL :: fx_tj(is:ie+1, js:je)
! 1-D Y-direction Fluxes
    REAL :: fy(is:ie, js:je+1)
    REAL :: fy_tl(is:ie, js:je+1)
    REAL :: fy_tj(is:ie, js:je+1)
! 1-D X-direction Fluxes
    REAL :: px(is:ie+1, js:je)
    REAL :: px_tl(is:ie+1, js:je)
    REAL :: px_tj(is:ie+1, js:je)
! 1-D Y-direction Fluxes
    REAL :: py(is:ie, js:je+1)
    REAL :: py_tl(is:ie, js:je+1)
    REAL :: py_tj(is:ie, js:je+1)
    REAL :: ra_x(is:ie, jsd:jed)
    REAL :: ra_x_tl(is:ie, jsd:jed)
    REAL :: ra_y(isd:ied, js:je)
    REAL :: ra_y_tl(isd:ied, js:je)
! work x-dir flux
    REAL :: gx(is:ie+1, js:je)
    REAL :: gx_tl(is:ie+1, js:je)
! work Y-dir flux array
    REAL :: gy(is:ie, js:je+1)
    REAL :: gy_tl(is:ie, js:je+1)
    LOGICAL :: fill_c
!#ifdef DOUBLE_GRADIENTS
!      real*8 :: pt_dp, delp_dp0, delp_dp
!      real*8 :: pt_fx_r8(is:ie+1,js:je+1)
!      real*8 :: pt_fy_r8(is:ie  ,js:je+1)
!      real*8 :: dp_fx_r8(is:ie+1,js:je  )
!      real*8 :: dp_fy_r8(is:ie  ,js:je+1)
!#endif
    REAL :: damp2, damp4, dd8, u2, v2, du2, dv2
    REAL :: damp2_tl, u2_tl, v2_tl, du2_tl, dv2_tl
    INTEGER :: i, j, is2, ie1, js2, je1, n, nt, n2, iq
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC ABS
    REAL :: pwx1
    REAL :: y1_tl
    REAL :: y2_tl
    INTEGER :: min4
    INTEGER :: min3
    INTEGER :: min2
    INTEGER :: min1
    REAL :: y3_tl
    REAL :: abs0_tl
    REAL :: abs0
    REAL :: max6
    REAL :: max5_tl
    REAL :: max5
    INTEGER :: max4
    INTEGER :: max3
    INTEGER :: max2
    INTEGER :: max1
    REAL :: y3
    REAL :: y2
    REAL :: max6_tl
    REAL :: y1
! shallow water test case 1
    IF (test_case .EQ. 1) THEN
      DO j=jsd,jed
        DO i=is,ie+1
          xfx_adv_tl(i, j) = dt*uc_tl(i, j)/sina_u(i, j)
          xfx_adv(i, j) = dt*uc(i, j)/sina_u(i, j)
          IF (xfx_adv(i, j) .GT. 0.) THEN
            crx_adv_tl(i, j) = rdxa(i-1, j)*xfx_adv_tl(i, j)
            crx_adv(i, j) = xfx_adv(i, j)*rdxa(i-1, j)
          ELSE
            crx_adv_tl(i, j) = rdxa(i, j)*xfx_adv_tl(i, j)
            crx_adv(i, j) = xfx_adv(i, j)*rdxa(i, j)
          END IF
          xfx_adv_tl(i, j) = dy(i, j)*sina_u(i, j)*xfx_adv_tl(i, j)
          xfx_adv(i, j) = dy(i, j)*xfx_adv(i, j)*sina_u(i, j)
        END DO
      END DO
      DO j=js,je+1
        DO i=isd,ied
          yfx_adv_tl(i, j) = dt*vc_tl(i, j)/sina_v(i, j)
          yfx_adv(i, j) = dt*vc(i, j)/sina_v(i, j)
          IF (yfx_adv(i, j) .GT. 0.) THEN
            cry_adv_tl(i, j) = rdya(i, j-1)*yfx_adv_tl(i, j)
            cry_adv(i, j) = yfx_adv(i, j)*rdya(i, j-1)
          ELSE
            cry_adv_tl(i, j) = rdya(i, j)*yfx_adv_tl(i, j)
            cry_adv(i, j) = yfx_adv(i, j)*rdya(i, j)
          END IF
          yfx_adv_tl(i, j) = dx(i, j)*sina_v(i, j)*yfx_adv_tl(i, j)
          yfx_adv(i, j) = dx(i, j)*yfx_adv(i, j)*sina_v(i, j)
        END DO
      END DO
      ut_tl = 0.0_8
      ra_x_tl = 0.0_8
      vt_tl = 0.0_8
    ELSE
! end grid_type choices
      IF (grid_type .LT. 3) THEN
        uci_tl = uc_tl
        uci = uc
        vci_tl = vc_tl
        vci = vc
        ut_tl = 0.0_8
! Interior:
        DO j=jsd,jed
          IF (j .NE. 0 .AND. j .NE. 1 .AND. j .NE. npy - 1 .AND. j .NE. &
&             npy) THEN
            DO i=is-1,ie+2
              ut_tl(i, j) = rsin_u(i, j)*(uci_tl(i, j)-r1b4*cosa_u(i, j)&
&               *(vci_tl(i-1, j)+vci_tl(i, j)+vci_tl(i-1, j+1)+vci_tl(i&
&               , j+1)))
              ut(i, j) = (uci(i, j)-r1b4*cosa_u(i, j)*(vci(i-1, j)+vci(i&
&               , j)+vci(i-1, j+1)+vci(i, j+1)))*rsin_u(i, j)
            END DO
          END IF
        END DO
        vt_tl = 0.0_8
        DO j=js-1,je+2
          IF (j .NE. 1 .AND. j .NE. npy) THEN
            DO i=isd,ied
              vt_tl(i, j) = rsin_v(i, j)*(vci_tl(i, j)-r1b4*cosa_v(i, j)&
&               *(uci_tl(i, j-1)+uci_tl(i+1, j-1)+uci_tl(i, j)+uci_tl(i+&
&               1, j)))
              vt(i, j) = (vci(i, j)-r1b4*cosa_v(i, j)*(uci(i, j-1)+uci(i&
&               +1, j-1)+uci(i, j)+uci(i+1, j)))*rsin_v(i, j)
            END DO
          END IF
        END DO
! West edge:
        IF (is .EQ. 1) THEN
          DO j=jsd,jed
            ut_tl(1, j) = rsin_u(1, j)*uci_tl(1, j)
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
            vt_tl(0, j) = vci_tl(0, j) + r1b4*cosa_v(1, j)*(ut_tl(0, j-1&
&             )+ut_tl(1, j-1)+ut_tl(0, j)+ut_tl(1, j))
            vt(0, j) = vci(0, j) + r1b4*cosa_v(1, j)*(ut(0, j-1)+ut(1, j&
&             -1)+ut(0, j)+ut(1, j))
            vt_tl(1, j) = vci_tl(1, j) - r1b4*cosa_v(1, j)*(ut_tl(1, j-1&
&             )+ut_tl(2, j-1)+ut_tl(1, j)+ut_tl(2, j))
            vt(1, j) = vci(1, j) - r1b4*cosa_v(1, j)*(ut(1, j-1)+ut(2, j&
&             -1)+ut(1, j)+ut(2, j))
          END DO
        END IF
! East edge:
        IF (ie + 1 .EQ. npx) THEN
          DO j=jsd,jed
            ut_tl(npx, j) = rsin_u(npx, j)*uci_tl(npx, j)
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
            vt_tl(npx-1, j) = vci_tl(npx-1, j) - r1b4*cosa_v(npx-1, j)*(&
&             ut_tl(npx-1, j-1)+ut_tl(npx, j-1)+ut_tl(npx-1, j)+ut_tl(&
&             npx, j))
            vt(npx-1, j) = vci(npx-1, j) - r1b4*cosa_v(npx-1, j)*(ut(npx&
&             -1, j-1)+ut(npx, j-1)+ut(npx-1, j)+ut(npx, j))
!             vt(npx,j) = vci(npx,j) - r1b4*cosa_v(npx,j)*   &
            vt_tl(npx, j) = vci_tl(npx, j) + r1b4*cosa_v(npx-1, j)*(&
&             ut_tl(npx, j-1)+ut_tl(npx+1, j-1)+ut_tl(npx, j)+ut_tl(npx+&
&             1, j))
            vt(npx, j) = vci(npx, j) + r1b4*cosa_v(npx-1, j)*(ut(npx, j-&
&             1)+ut(npx+1, j-1)+ut(npx, j)+ut(npx+1, j))
          END DO
        END IF
! South (Bottom) edge:
        IF (js .EQ. 1) THEN
          DO i=isd,ied
            vt_tl(i, 1) = rsin_v(i, 1)*vci_tl(i, 1)
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
            ut_tl(i, 0) = uci_tl(i, 0) + r1b4*cosa_u(i, 1)*(vt_tl(i-1, 0&
&             )+vt_tl(i, 0)+vt_tl(i-1, 1)+vt_tl(i, 1))
            ut(i, 0) = uci(i, 0) + r1b4*cosa_u(i, 1)*(vt(i-1, 0)+vt(i, 0&
&             )+vt(i-1, 1)+vt(i, 1))
            ut_tl(i, 1) = uci_tl(i, 1) - r1b4*cosa_u(i, 1)*(vt_tl(i-1, 1&
&             )+vt_tl(i, 1)+vt_tl(i-1, 2)+vt_tl(i, 2))
            ut(i, 1) = uci(i, 1) - r1b4*cosa_u(i, 1)*(vt(i-1, 1)+vt(i, 1&
&             )+vt(i-1, 2)+vt(i, 2))
          END DO
        END IF
! North edge:
        IF (je + 1 .EQ. npy) THEN
          DO i=isd,ied
            vt_tl(i, npy) = rsin_v(i, npy)*vci_tl(i, npy)
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
            ut_tl(i, npy-1) = uci_tl(i, npy-1) - r1b4*cosa_u(i, npy-1)*(&
&             vt_tl(i-1, npy-1)+vt_tl(i, npy-1)+vt_tl(i-1, npy)+vt_tl(i&
&             , npy))
            ut(i, npy-1) = uci(i, npy-1) - r1b4*cosa_u(i, npy-1)*(vt(i-1&
&             , npy-1)+vt(i, npy-1)+vt(i-1, npy)+vt(i, npy))
!             ut(i,npy) = uci(i,npy) - r1b4*cosa_u(i,npy)*   &
            ut_tl(i, npy) = uci_tl(i, npy) + r1b4*cosa_u(i, npy-1)*(&
&             vt_tl(i-1, npy)+vt_tl(i, npy)+vt_tl(i-1, npy+1)+vt_tl(i, &
&             npy+1))
            ut(i, npy) = uci(i, npy) + r1b4*cosa_u(i, npy-1)*(vt(i-1, &
&             npy)+vt(i, npy)+vt(i-1, npy+1)+vt(i, npy+1))
          END DO
        END IF
        IF (sw_corner) THEN
          damp = 1.0/(1.0-r1b16*cosa_u(2, 1)*cosa_v(1, 2))
          ut_tl(2, 0) = damp*(uci_tl(2, 0)-r1b4*cosa_u(2, 0)*(vt_tl(1, 1&
&           )+vt_tl(2, 1)+vt_tl(2, 0)+vci_tl(1, 0)-r1b4*cosa_v(1, 0)*(&
&           ut_tl(1, 0)+ut_tl(1, -1)+ut_tl(2, -1))))
          ut(2, 0) = (uci(2, 0)-r1b4*cosa_u(2, 0)*(vt(1, 1)+vt(2, 1)+vt(&
&           2, 0)+vci(1, 0)-r1b4*cosa_v(1, 0)*(ut(1, 0)+ut(1, -1)+ut(2, &
&           -1))))*damp
          ut_tl(2, 1) = damp*(uci_tl(2, 1)-r1b4*cosa_u(2, 1)*(vt_tl(1, 1&
&           )+vt_tl(2, 1)+vt_tl(2, 2)+vci_tl(1, 2)-r1b4*cosa_v(1, 2)*(&
&           ut_tl(1, 1)+ut_tl(1, 2)+ut_tl(2, 2))))
          ut(2, 1) = (uci(2, 1)-r1b4*cosa_u(2, 1)*(vt(1, 1)+vt(2, 1)+vt(&
&           2, 2)+vci(1, 2)-r1b4*cosa_v(1, 2)*(ut(1, 1)+ut(1, 2)+ut(2, 2&
&           ))))*damp
          vt_tl(1, 2) = damp*(vci_tl(1, 2)-r1b4*cosa_v(1, 2)*(ut_tl(1, 1&
&           )+ut_tl(1, 2)+ut_tl(2, 2)+uci_tl(2, 1)-r1b4*cosa_u(2, 1)*(&
&           vt_tl(1, 1)+vt_tl(2, 1)+vt_tl(2, 2))))
          vt(1, 2) = (vci(1, 2)-r1b4*cosa_v(1, 2)*(ut(1, 1)+ut(1, 2)+ut(&
&           2, 2)+uci(2, 1)-r1b4*cosa_u(2, 1)*(vt(1, 1)+vt(2, 1)+vt(2, 2&
&           ))))*damp
          vt_tl(0, 2) = damp*(vci_tl(0, 2)-r1b4*cosa_v(0, 2)*(ut_tl(1, 1&
&           )+ut_tl(1, 2)+ut_tl(0, 2)+uci_tl(0, 1)-r1b4*cosa_u(0, 1)*(&
&           vt_tl(0, 1)+vt_tl(-1, 1)+vt_tl(-1, 2))))
          vt(0, 2) = (vci(0, 2)-r1b4*cosa_v(0, 2)*(ut(1, 1)+ut(1, 2)+ut(&
&           0, 2)+uci(0, 1)-r1b4*cosa_u(0, 1)*(vt(0, 1)+vt(-1, 1)+vt(-1&
&           , 2))))*damp
        END IF
        IF (se_corner) THEN
          damp = 1.0/(1.0-r1b16*cosa_u(npx-1, 1)*cosa_v(npx-1, 2))
          ut_tl(npx-1, 0) = damp*(uci_tl(npx-1, 0)+r1b4*cosa_u(npx-1, 1)&
&           *(vt_tl(npx-1, 1)+vt_tl(npx-2, 1)+vt_tl(npx-2, 0)+vci_tl(npx&
&           -1, 0)+r1b4*cosa_v(npx-1, 2)*(ut_tl(npx, 0)+ut_tl(npx, -1)+&
&           ut_tl(npx-1, -1))))
          ut(npx-1, 0) = (uci(npx-1, 0)+r1b4*cosa_u(npx-1, 1)*(vt(npx-1&
&           , 1)+vt(npx-2, 1)+vt(npx-2, 0)+vci(npx-1, 0)+r1b4*cosa_v(npx&
&           -1, 2)*(ut(npx, 0)+ut(npx, -1)+ut(npx-1, -1))))*damp
          ut_tl(npx-1, 1) = damp*(uci_tl(npx-1, 1)-r1b4*cosa_u(npx-1, 1)&
&           *(vt_tl(npx-1, 1)+vt_tl(npx-2, 1)+vt_tl(npx-2, 2)+vci_tl(npx&
&           -1, 2)-r1b4*cosa_v(npx-1, 2)*(ut_tl(npx, 1)+ut_tl(npx, 2)+&
&           ut_tl(npx-1, 2))))
          ut(npx-1, 1) = (uci(npx-1, 1)-r1b4*cosa_u(npx-1, 1)*(vt(npx-1&
&           , 1)+vt(npx-2, 1)+vt(npx-2, 2)+vci(npx-1, 2)-r1b4*cosa_v(npx&
&           -1, 2)*(ut(npx, 1)+ut(npx, 2)+ut(npx-1, 2))))*damp
          vt_tl(npx-1, 2) = damp*(vci_tl(npx-1, 2)-r1b4*cosa_v(npx-1, 2)&
&           *(ut_tl(npx, 1)+ut_tl(npx, 2)+ut_tl(npx-1, 2)+uci_tl(npx-1, &
&           1)-r1b4*cosa_u(npx-1, 1)*(vt_tl(npx-1, 1)+vt_tl(npx-2, 1)+&
&           vt_tl(npx-2, 2))))
          vt(npx-1, 2) = (vci(npx-1, 2)-r1b4*cosa_v(npx-1, 2)*(ut(npx, 1&
&           )+ut(npx, 2)+ut(npx-1, 2)+uci(npx-1, 1)-r1b4*cosa_u(npx-1, 1&
&           )*(vt(npx-1, 1)+vt(npx-2, 1)+vt(npx-2, 2))))*damp
          vt_tl(npx, 2) = damp*(vci_tl(npx, 2)+r1b4*cosa_v(npx-1, 2)*(&
&           ut_tl(npx, 1)+ut_tl(npx, 2)+ut_tl(npx+1, 2)+uci_tl(npx+1, 1)&
&           +r1b4*cosa_u(npx-1, 1)*(vt_tl(npx, 1)+vt_tl(npx+1, 1)+vt_tl(&
&           npx+1, 2))))
          vt(npx, 2) = (vci(npx, 2)+r1b4*cosa_v(npx-1, 2)*(ut(npx, 1)+ut&
&           (npx, 2)+ut(npx+1, 2)+uci(npx+1, 1)+r1b4*cosa_u(npx-1, 1)*(&
&           vt(npx, 1)+vt(npx+1, 1)+vt(npx+1, 2))))*damp
        END IF
        IF (ne_corner) THEN
          damp = 1.0/(1.0-r1b16*cosa_u(npx-1, npy-1)*cosa_v(npx-1, npy-1&
&           ))
          ut_tl(npx-1, npy) = damp*(uci_tl(npx-1, npy)+r1b4*cosa_u(npx-1&
&           , npy-1)*(vt_tl(npx-1, npy)+vt_tl(npx-2, npy)+vt_tl(npx-2, &
&           npy+1)+vci_tl(npx-1, npy+1)+r1b4*cosa_v(npx-1, npy-1)*(ut_tl&
&           (npx, npy)+ut_tl(npx, npy+1)+ut_tl(npx-1, npy+1))))
          ut(npx-1, npy) = (uci(npx-1, npy)+r1b4*cosa_u(npx-1, npy-1)*(&
&           vt(npx-1, npy)+vt(npx-2, npy)+vt(npx-2, npy+1)+vci(npx-1, &
&           npy+1)+r1b4*cosa_v(npx-1, npy-1)*(ut(npx, npy)+ut(npx, npy+1&
&           )+ut(npx-1, npy+1))))*damp
          ut_tl(npx-1, npy-1) = damp*(uci_tl(npx-1, npy-1)-r1b4*cosa_u(&
&           npx-1, npy-1)*(vt_tl(npx-1, npy)+vt_tl(npx-2, npy)+vt_tl(npx&
&           -2, npy-1)+vci_tl(npx-1, npy-1)-r1b4*cosa_v(npx-1, npy-1)*(&
&           ut_tl(npx, npy-1)+ut_tl(npx, npy-2)+ut_tl(npx-1, npy-2))))
          ut(npx-1, npy-1) = (uci(npx-1, npy-1)-r1b4*cosa_u(npx-1, npy-1&
&           )*(vt(npx-1, npy)+vt(npx-2, npy)+vt(npx-2, npy-1)+vci(npx-1&
&           , npy-1)-r1b4*cosa_v(npx-1, npy-1)*(ut(npx, npy-1)+ut(npx, &
&           npy-2)+ut(npx-1, npy-2))))*damp
          vt_tl(npx-1, npy-1) = damp*(vci_tl(npx-1, npy-1)-r1b4*cosa_v(&
&           npx-1, npy-1)*(ut_tl(npx, npy-1)+ut_tl(npx, npy-2)+ut_tl(npx&
&           -1, npy-2)+uci_tl(npx-1, npy-1)-r1b4*cosa_u(npx-1, npy-1)*(&
&           vt_tl(npx-1, npy)+vt_tl(npx-2, npy)+vt_tl(npx-2, npy-1))))
          vt(npx-1, npy-1) = (vci(npx-1, npy-1)-r1b4*cosa_v(npx-1, npy-1&
&           )*(ut(npx, npy-1)+ut(npx, npy-2)+ut(npx-1, npy-2)+uci(npx-1&
&           , npy-1)-r1b4*cosa_u(npx-1, npy-1)*(vt(npx-1, npy)+vt(npx-2&
&           , npy)+vt(npx-2, npy-1))))*damp
          vt_tl(npx, npy-1) = damp*(vci_tl(npx, npy-1)+r1b4*cosa_v(npx-1&
&           , npy-1)*(ut_tl(npx, npy-1)+ut_tl(npx, npy-2)+ut_tl(npx+1, &
&           npy-2)+uci_tl(npx+1, npy-1)+r1b4*cosa_u(npx-1, npy-1)*(vt_tl&
&           (npx, npy)+vt_tl(npx+1, npy)+vt_tl(npx+1, npy-1))))
          vt(npx, npy-1) = (vci(npx, npy-1)+r1b4*cosa_v(npx-1, npy-1)*(&
&           ut(npx, npy-1)+ut(npx, npy-2)+ut(npx+1, npy-2)+uci(npx+1, &
&           npy-1)+r1b4*cosa_u(npx-1, npy-1)*(vt(npx, npy)+vt(npx+1, npy&
&           )+vt(npx+1, npy-1))))*damp
        END IF
        IF (nw_corner) THEN
          damp = 1.0/(1.0-r1b16*cosa_u(2, npy-1)*cosa_v(1, npy-1))
          ut_tl(2, npy) = damp*(uci_tl(2, npy)+r1b4*cosa_u(2, npy-1)*(&
&           vt_tl(1, npy)+vt_tl(2, npy)+vt_tl(2, npy+1)+vci_tl(1, npy+1)&
&           +r1b4*cosa_v(1, npy-1)*(ut_tl(1, npy)+ut_tl(1, npy+1)+ut_tl(&
&           2, npy+1))))
          ut(2, npy) = (uci(2, npy)+r1b4*cosa_u(2, npy-1)*(vt(1, npy)+vt&
&           (2, npy)+vt(2, npy+1)+vci(1, npy+1)+r1b4*cosa_v(1, npy-1)*(&
&           ut(1, npy)+ut(1, npy+1)+ut(2, npy+1))))*damp
          ut_tl(2, npy-1) = damp*(uci_tl(2, npy-1)-r1b4*cosa_u(2, npy-1)&
&           *(vt_tl(1, npy)+vt_tl(2, npy)+vt_tl(2, npy-1)+vci_tl(1, npy-&
&           1)-r1b4*cosa_v(1, npy-1)*(ut_tl(1, npy-1)+ut_tl(1, npy-2)+&
&           ut_tl(2, npy-2))))
          ut(2, npy-1) = (uci(2, npy-1)-r1b4*cosa_u(2, npy-1)*(vt(1, npy&
&           )+vt(2, npy)+vt(2, npy-1)+vci(1, npy-1)-r1b4*cosa_v(1, npy-1&
&           )*(ut(1, npy-1)+ut(1, npy-2)+ut(2, npy-2))))*damp
          vt_tl(1, npy-1) = damp*(vci_tl(1, npy-1)-r1b4*cosa_v(1, npy-1)&
&           *(ut_tl(1, npy-1)+ut_tl(1, npy-2)+ut_tl(2, npy-2)+uci_tl(2, &
&           npy-1)-r1b4*cosa_u(2, npy-1)*(vt_tl(1, npy)+vt_tl(2, npy)+&
&           vt_tl(2, npy-1))))
          vt(1, npy-1) = (vci(1, npy-1)-r1b4*cosa_v(1, npy-1)*(ut(1, npy&
&           -1)+ut(1, npy-2)+ut(2, npy-2)+uci(2, npy-1)-r1b4*cosa_u(2, &
&           npy-1)*(vt(1, npy)+vt(2, npy)+vt(2, npy-1))))*damp
          vt_tl(0, npy-1) = damp*(vci_tl(0, npy-1)+r1b4*cosa_v(1, npy-1)&
&           *(ut_tl(1, npy-1)+ut_tl(1, npy-2)+ut_tl(0, npy-2)+uci_tl(0, &
&           npy-1)+r1b4*cosa_u(2, npy-1)*(vt_tl(0, npy)+vt_tl(-1, npy)+&
&           vt_tl(-1, npy-1))))
          vt(0, npy-1) = (vci(0, npy-1)+r1b4*cosa_v(1, npy-1)*(ut(1, npy&
&           -1)+ut(1, npy-2)+ut(0, npy-2)+uci(0, npy-1)+r1b4*cosa_u(2, &
&           npy-1)*(vt(0, npy)+vt(-1, npy)+vt(-1, npy-1))))*damp
        END IF
      ELSE
        ut_tl = 0.0_8
! grid_type >= 3
        DO j=jsd,jed
          DO i=is-1,ie+2
            ut_tl(i, j) = uc_tl(i, j)
            ut(i, j) = uc(i, j)
          END DO
        END DO
        vt_tl = 0.0_8
        DO j=js-1,je+2
          DO i=isd,ied
            vt_tl(i, j) = vc_tl(i, j)
            vt(i, j) = vc(i, j)
          END DO
        END DO
      END IF
      DO j=jsd,jed
        DO i=is,ie+1
          xfx_adv_tl(i, j) = dt*ut_tl(i, j)
          xfx_adv(i, j) = dt*ut(i, j)
        END DO
      END DO
      DO j=js,je+1
        DO i=isd,ied
          yfx_adv_tl(i, j) = dt*vt_tl(i, j)
          yfx_adv(i, j) = dt*vt(i, j)
        END DO
      END DO
! Compute E-W CFL number:
      DO j=jsd,jed
        DO i=is,ie+1
          IF (xfx_adv(i, j) .GT. 0.) THEN
            crx_adv_tl(i, j) = rdxa(i-1, j)*xfx_adv_tl(i, j)
            crx_adv(i, j) = xfx_adv(i, j)*rdxa(i-1, j)
          ELSE
            crx_adv_tl(i, j) = rdxa(i, j)*xfx_adv_tl(i, j)
            crx_adv(i, j) = xfx_adv(i, j)*rdxa(i, j)
          END IF
        END DO
      END DO
      DO j=jsd,jed
        DO i=is,ie+1
          xfx_adv_tl(i, j) = dy(i, j)*sina_u(i, j)*xfx_adv_tl(i, j)
          xfx_adv(i, j) = dy(i, j)*xfx_adv(i, j)*sina_u(i, j)
        END DO
      END DO
      DO j=js,je+1
        DO i=isd,ied
          IF (yfx_adv(i, j) .GT. 0.) THEN
            cry_adv_tl(i, j) = rdya(i, j-1)*yfx_adv_tl(i, j)
            cry_adv(i, j) = yfx_adv(i, j)*rdya(i, j-1)
          ELSE
            cry_adv_tl(i, j) = rdya(i, j)*yfx_adv_tl(i, j)
            cry_adv(i, j) = yfx_adv(i, j)*rdya(i, j)
          END IF
        END DO
      END DO
      DO j=js,je+1
        DO i=isd,ied
          yfx_adv_tl(i, j) = dx(i, j)*sina_v(i, j)*yfx_adv_tl(i, j)
          yfx_adv(i, j) = dx(i, j)*yfx_adv(i, j)*sina_v(i, j)
        END DO
      END DO
      ra_x_tl = 0.0_8
    END IF
    DO j=jsd,jed
      DO i=is,ie
        ra_x_tl(i, j) = xfx_adv_tl(i, j) - xfx_adv_tl(i+1, j)
        ra_x(i, j) = area(i, j) + xfx_adv(i, j) - xfx_adv(i+1, j)
      END DO
    END DO
    ra_y_tl = 0.0_8
    DO j=js,je
      DO i=isd,ied
        ra_y_tl(i, j) = yfx_adv_tl(i, j) - yfx_adv_tl(i, j+1)
        ra_y(i, j) = area(i, j) + yfx_adv(i, j) - yfx_adv(i, j+1)
      END DO
    END DO
    fy_tl = 0.0
    fx_tl = 0.0
    fy_tj = fy
    fx_tj = fx
    delp_tj = delp
    CALL FV_TP_2D(delp, crx_adv, cry_adv, &
&               npx, npy, hord_dp, fx, fy, &
&               xfx_adv, yfx_adv, area, ra_x, &
&               ra_y)
    CALL FV_TP_2D_TLM(delp_tj, delp_tl, crx_adv, crx_adv_tl, cry_adv, &
&               cry_adv_tl, npx, npy, hord_dp_pert, fx_tj, fx_tl, fy_tj, fy_tl, &
&               xfx_adv, xfx_adv_tl, yfx_adv, yfx_adv_tl, area, ra_x, &
&               ra_x_tl, ra_y, ra_y_tl)
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
          delp_tl(i, j) = delp_tl(i, j) + rarea(i, j)*(fx_tl(i, j)-fx_tl&
&           (i+1, j)+fy_tl(i, j)-fy_tl(i, j+1))
          delp(i, j) = delp(i, j) + (fx(i, j)-fx(i+1, j)+fy(i, j)-fy(i, &
&           j+1))*rarea(i, j)
          ptc_tl(i, j) = pt_tl(i, j)
          ptc(i, j) = pt(i, j)
        END DO
      END DO
      wk_tl = 0.0_8
    ELSE
!#endif
! <<< Save the mass fluxes to the "Flux Capacitor" for tracer transport >>>
      DO j=jsd,jed
        DO i=is,ie+1
          cx_tl(i, j) = cx_tl(i, j) + crx_adv_tl(i, j)
          cx(i, j) = cx(i, j) + crx_adv(i, j)
        END DO
      END DO
      DO j=js,je
        DO i=is,ie+1
          xflux_tl(i, j) = xflux_tl(i, j) + fx_tl(i, j)
          xflux(i, j) = xflux(i, j) + fx(i, j)
        END DO
      END DO
      DO j=js,je+1
        DO i=isd,ied
          cy_tl(i, j) = cy_tl(i, j) + cry_adv_tl(i, j)
          cy(i, j) = cy(i, j) + cry_adv(i, j)
        END DO
        DO i=is,ie
          yflux_tl(i, j) = yflux_tl(i, j) + fy_tl(i, j)
          yflux(i, j) = yflux(i, j) + fy(i, j)
        END DO
      END DO
      IF (.NOT.hydrostatic) THEN
        py_tl = 0.0
        px_tl = 0.0
        w_tj = w
        px_tj = px
        py_tj = py
        CALL FV_TP_2D(w, crx_adv, cry_adv, &
&                   npx, npy, hord_vt, px, py &
&                   , xfx_adv, yfx_adv, area, &
&                   ra_x, ra_y, mfx=fx, mfy=fy)
        CALL FV_TP_2D_TLM(w_tj, w_tl, crx_adv, crx_adv_tl, cry_adv, &
&                   cry_adv_tl, npx, npy, hord_vt_pert, px_tj, px_tl, py_tj, py_tl&
&                   , xfx_adv, xfx_adv_tl, yfx_adv, yfx_adv_tl, area, &
&                   ra_x, ra_x_tl, ra_y, ra_y_tl, mfx=fx, mfx_tl=fx_tl, &
&                   mfy=fy, mfy_tl=fy_tl)
        DO j=js,je
          DO i=is,ie
            w_tl(i, j) = w_tl(i, j)*delp(i, j) + w(i, j)*delp_tl(i, j) +&
&             rarea(i, j)*(px_tl(i, j)-px_tl(i+1, j)+py_tl(i, j)-py_tl(i&
&             , j+1))
            w(i, j) = w(i, j)*delp(i, j) + (px(i, j)-px(i+1, j)+py(i, j)&
&             -py(i, j+1))*rarea(i, j)
          END DO
        END DO
      ELSE
        px_tl = 0.0_8
        py_tl = 0.0_8
      END IF
      IF (inline_q) THEN
        DO j=jsd,jed
          DO i=isd,ied
            pt_tl(i, j) = (pt_tl(i, j)*(1.+zvir*q(i, j, k, sphum))-pt(i&
&             , j)*zvir*q_tl(i, j, k, sphum))/(1.+zvir*q(i, j, k, sphum)&
&             )**2
            pt(i, j) = pt(i, j)/(1.+zvir*q(i, j, k, sphum))
          END DO
        END DO
      END IF
      px_tj = px
      py_tj = py
      pt_tj = pt
      CALL FV_TP_2D(pt, crx_adv, cry_adv, &
&                 npx, npy, hord_tm, px, py, &
&                 xfx_adv, yfx_adv, area, ra_x, &
&                 ra_y, mfx=fx, mfy=fy)
      CALL FV_TP_2D_TLM(pt_tj, pt_tl, crx_adv, crx_adv_tl, cry_adv, &
&                 cry_adv_tl, npx, npy, hord_tm_pert, px_tj, px_tl, py_tj, py_tl, &
&                 xfx_adv, xfx_adv_tl, yfx_adv, yfx_adv_tl, area, ra_x, &
&                 ra_x_tl, ra_y, ra_y_tl, mfx=fx, mfx_tl=fx_tl, mfy=fy, &
&                 mfy_tl=fy_tl)
      IF (inline_q) THEN
        wk_tl = 0.0_8
        DO j=js,je
          DO i=is,ie
            wk_tl(i, j) = delp_tl(i, j)
            wk(i, j) = delp(i, j)
            delp_tl(i, j) = wk_tl(i, j) + rarea(i, j)*(fx_tl(i, j)-fx_tl&
&             (i+1, j)+fy_tl(i, j)-fy_tl(i, j+1))
            delp(i, j) = wk(i, j) + (fx(i, j)-fx(i+1, j)+fy(i, j)-fy(i, &
&             j+1))*rarea(i, j)
            pt_tl(i, j) = ((pt_tl(i, j)*wk(i, j)+pt(i, j)*wk_tl(i, j)+&
&             rarea(i, j)*(px_tl(i, j)-px_tl(i+1, j)+py_tl(i, j)-py_tl(i&
&             , j+1)))*delp(i, j)-(pt(i, j)*wk(i, j)+(px(i, j)-px(i+1, j&
&             )+py(i, j)-py(i, j+1))*rarea(i, j))*delp_tl(i, j))/delp(i&
&             , j)**2
            pt(i, j) = (pt(i, j)*wk(i, j)+(px(i, j)-px(i+1, j)+py(i, j)-&
&             py(i, j+1))*rarea(i, j))/delp(i, j)
          END DO
        END DO
        DO iq=1,nq
          px_tj = px
          py_tj = py
          q_tj = q
          CALL FV_TP_2D(q(isd:, jsd:, k, iq), crx_adv, cry_adv, npx, &
&                     npy, hord_tr, px, py, xfx_adv, &
&                     yfx_adv, area, ra_x, &
&                     ra_y, mfx=fx, mfy=fy)
          CALL FV_TP_2D_TLM(q_tj(isd:, jsd:, k, iq), q_tl(isd:, jsd:, k, iq&
&                     ), crx_adv, crx_adv_tl, cry_adv, cry_adv_tl, npx, &
&                     npy, hord_tr_pert, px_tj, px_tl, py_tj, py_tl, xfx_adv, &
&                     xfx_adv_tl, yfx_adv, yfx_adv_tl, area, ra_x, &
&                     ra_x_tl, ra_y, ra_y_tl, mfx=fx, mfx_tl=fx_tl, mfy=&
&                     fy, mfy_tl=fy_tl)
          DO j=js,je
            DO i=is,ie
              q_tl(i, j, k, iq) = ((q_tl(i, j, k, iq)*wk(i, j)+q(i, j, k&
&               , iq)*wk_tl(i, j)+rarea(i, j)*(px_tl(i, j)-px_tl(i+1, j)&
&               +py_tl(i, j)-py_tl(i, j+1)))*delp(i, j)-(q(i, j, k, iq)*&
&               wk(i, j)+(px(i, j)-px(i+1, j)+py(i, j)-py(i, j+1))*rarea&
&               (i, j))*delp_tl(i, j))/delp(i, j)**2
              q(i, j, k, iq) = (q(i, j, k, iq)*wk(i, j)+(px(i, j)-px(i+1&
&               , j)+py(i, j)-py(i, j+1))*rarea(i, j))/delp(i, j)
            END DO
          END DO
        END DO
        DO j=js,je
          DO i=is,ie
            pt_tl(i, j) = pt_tl(i, j)*(1.+zvir*q(i, j, k, sphum)) + pt(i&
&             , j)*zvir*q_tl(i, j, k, sphum)
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
            pt_tl(i, j) = pt_tl(i, j)*delp(i, j) + pt(i, j)*delp_tl(i, j&
&             ) + rarea(i, j)*(px_tl(i, j)-px_tl(i+1, j)+py_tl(i, j)-&
&             py_tl(i, j+1))
            pt(i, j) = pt(i, j)*delp(i, j) + (px(i, j)-px(i+1, j)+py(i, &
&             j)-py(i, j+1))*rarea(i, j)
            delp_tl(i, j) = delp_tl(i, j) + rarea(i, j)*(fx_tl(i, j)-&
&             fx_tl(i+1, j)+fy_tl(i, j)-fy_tl(i, j+1))
            delp(i, j) = delp(i, j) + (fx(i, j)-fx(i+1, j)+fy(i, j)-fy(i&
&             , j+1))*rarea(i, j)
            pt_tl(i, j) = (pt_tl(i, j)*delp(i, j)-pt(i, j)*delp_tl(i, j)&
&             )/delp(i, j)**2
            pt(i, j) = pt(i, j)/delp(i, j)
          END DO
        END DO
        wk_tl = 0.0_8
      END IF
!#endif
      IF (.NOT.hydrostatic) THEN
        DO j=js,je
          DO i=is,ie
            w_tl(i, j) = (w_tl(i, j)*delp(i, j)-w(i, j)*delp_tl(i, j))/&
&             delp(i, j)**2
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
          vb_tl = 0.0_8
          DO i=is,ie+1
! corner values are incorrect
            vb_tl(i, 1) = dt5*(vt_tl(i-1, 1)+vt_tl(i, 1))
            vb(i, 1) = dt5*(vt(i-1, 1)+vt(i, 1))
          END DO
        ELSE
          vb_tl = 0.0_8
        END IF
        DO j=js2,je1
          DO i=is2,ie1
            vb_tl(i, j) = dt5*rsina(i, j)*(vc_tl(i-1, j)+vc_tl(i, j)-&
&             cosa(i, j)*(uc_tl(i, j-1)+uc_tl(i, j)))
            vb(i, j) = dt5*(vc(i-1, j)+vc(i, j)-(uc(i, j-1)+uc(i, j))*&
&             cosa(i, j))*rsina(i, j)
          END DO
          IF (is .EQ. 1) THEN
!              vb(1,j) = dt5*(vt(0,j)+vt(1,j)) 
! 2-pt extrapolation from both sides:
            vb_tl(1, j) = dt4*(3.*(vt_tl(0, j)+vt_tl(1, j))-vt_tl(-1, j)&
&             -vt_tl(2, j))
            vb(1, j) = dt4*(-vt(-1, j)+3.*(vt(0, j)+vt(1, j))-vt(2, j))
          END IF
          IF (ie + 1 .EQ. npx) THEN
!              vb(npx,j) = dt5*(vt(npx-1,j)+vt(npx,j))
! 2-pt extrapolation from both sides:
            vb_tl(npx, j) = dt4*(3.*(vt_tl(npx-1, j)+vt_tl(npx, j))-&
&             vt_tl(npx-2, j)-vt_tl(npx+1, j))
            vb(npx, j) = dt4*(-vt(npx-2, j)+3.*(vt(npx-1, j)+vt(npx, j))&
&             -vt(npx+1, j))
          END IF
        END DO
        IF (je + 1 .EQ. npy) THEN
          DO i=is,ie+1
! corner values are incorrect
            vb_tl(i, npy) = dt5*(vt_tl(i-1, npy)+vt_tl(i, npy))
            vb(i, npy) = dt5*(vt(i-1, npy)+vt(i, npy))
          END DO
        END IF
      ELSE
        vb_tl = 0.0_8
        DO j=js,je+1
          DO i=is,ie+1
            vb_tl(i, j) = dt5*(vc_tl(i-1, j)+vc_tl(i, j))
            vb(i, j) = dt5*(vc(i-1, j)+vc(i, j))
          END DO
        END DO
      END IF
      ub_tj = ub
      CALL YTP_V(vb, u, v, ub, hord_mt)
      CALL YTP_V_TLM(vb, vb_tl, u, v, v_tl, ub_tj, ub_tl, hord_mt_pert)
      ke_tl = 0.0_8
      DO j=js,je+1
        DO i=is,ie+1
          ke_tl(i, j) = vb_tl(i, j)*ub(i, j) + vb(i, j)*ub_tl(i, j)
          ke(i, j) = vb(i, j)*ub(i, j)
        END DO
      END DO
      IF (grid_type .LT. 3) THEN
        IF (is .EQ. 1) THEN
          DO j=js,je+1
! corner values are incorrect
            ub_tl(1, j) = dt5*(ut_tl(1, j-1)+ut_tl(1, j))
            ub(1, j) = dt5*(ut(1, j-1)+ut(1, j))
          END DO
        END IF
        DO j=js,je+1
          IF (j .EQ. 1 .OR. j .EQ. npy) THEN
            DO i=is2,ie1
!                 ub(i,j) = dt5*(ut(i,j-1)+ut(i,j))
! 2-pt extrapolation from both sides:
              ub_tl(i, j) = dt4*(3.*(ut_tl(i, j-1)+ut_tl(i, j))-ut_tl(i&
&               , j-2)-ut_tl(i, j+1))
              ub(i, j) = dt4*(-ut(i, j-2)+3.*(ut(i, j-1)+ut(i, j))-ut(i&
&               , j+1))
            END DO
          ELSE
            DO i=is2,ie1
              ub_tl(i, j) = dt5*rsina(i, j)*(uc_tl(i, j-1)+uc_tl(i, j)-&
&               cosa(i, j)*(vc_tl(i-1, j)+vc_tl(i, j)))
              ub(i, j) = dt5*(uc(i, j-1)+uc(i, j)-(vc(i-1, j)+vc(i, j))*&
&               cosa(i, j))*rsina(i, j)
            END DO
          END IF
        END DO
        IF (ie + 1 .EQ. npx) THEN
          DO j=js,je+1
! corner values are incorrect
            ub_tl(npx, j) = dt5*(ut_tl(npx, j-1)+ut_tl(npx, j))
            ub(npx, j) = dt5*(ut(npx, j-1)+ut(npx, j))
          END DO
        END IF
      ELSE
        DO j=js,je+1
          DO i=is,ie+1
            ub_tl(i, j) = dt5*(uc_tl(i, j-1)+uc_tl(i, j))
            ub(i, j) = dt5*(uc(i, j-1)+uc(i, j))
          END DO
        END DO
      END IF
      vb_tj = vb
      CALL XTP_U(ub, u, v, vb, hord_mt)
      CALL XTP_U_TLM(ub, ub_tl, u, u_tl, v, vb_tj, vb_tl, hord_mt_pert)
      DO j=js,je+1
        DO i=is,ie+1
          ke_tl(i, j) = 0.5*(ke_tl(i, j)+ub_tl(i, j)*vb(i, j)+ub(i, j)*&
&           vb_tl(i, j))
          ke(i, j) = 0.5*(ke(i, j)+ub(i, j)*vb(i, j))
        END DO
      END DO
!-----------------------------------------
! Fix KE at the 4 corners of the face:
!-----------------------------------------
      IF (gnomonic_grid) THEN
        dt6 = dt/6.
        IF (sw_corner) THEN
          ke_tl(1, 1) = dt6*((ut_tl(1, 1)+ut_tl(1, 0))*u(1, 1)+(ut(1, 1)&
&           +ut(1, 0))*u_tl(1, 1)+(vt_tl(1, 1)+vt_tl(0, 1))*v(1, 1)+(vt(&
&           1, 1)+vt(0, 1))*v_tl(1, 1)+(ut_tl(1, 1)+vt_tl(1, 1))*u(0, 1)&
&           +(ut(1, 1)+vt(1, 1))*u_tl(0, 1))
          ke(1, 1) = dt6*((ut(1, 1)+ut(1, 0))*u(1, 1)+(vt(1, 1)+vt(0, 1)&
&           )*v(1, 1)+(ut(1, 1)+vt(1, 1))*u(0, 1))
        END IF
        IF (se_corner) THEN
          i = npx
          ke_tl(i, 1) = dt6*((ut_tl(i, 1)+ut_tl(i, 0))*u(i-1, 1)+(ut(i, &
&           1)+ut(i, 0))*u_tl(i-1, 1)+(vt_tl(i, 1)+vt_tl(i-1, 1))*v(i, 1&
&           )+(vt(i, 1)+vt(i-1, 1))*v_tl(i, 1)+(ut_tl(i, 1)-vt_tl(i-1, 1&
&           ))*u(i, 1)+(ut(i, 1)-vt(i-1, 1))*u_tl(i, 1))
          ke(i, 1) = dt6*((ut(i, 1)+ut(i, 0))*u(i-1, 1)+(vt(i, 1)+vt(i-1&
&           , 1))*v(i, 1)+(ut(i, 1)-vt(i-1, 1))*u(i, 1))
        END IF
        IF (ne_corner) THEN
          i = npx
          j = npy
          ke_tl(i, j) = dt6*((ut_tl(i, j)+ut_tl(i, j-1))*u(i-1, j)+(ut(i&
&           , j)+ut(i, j-1))*u_tl(i-1, j)+(vt_tl(i, j)+vt_tl(i-1, j))*v(&
&           i, j-1)+(vt(i, j)+vt(i-1, j))*v_tl(i, j-1)+(ut_tl(i, j-1)+&
&           vt_tl(i-1, j))*u(i, j)+(ut(i, j-1)+vt(i-1, j))*u_tl(i, j))
          ke(i, j) = dt6*((ut(i, j)+ut(i, j-1))*u(i-1, j)+(vt(i, j)+vt(i&
&           -1, j))*v(i, j-1)+(ut(i, j-1)+vt(i-1, j))*u(i, j))
        END IF
        IF (nw_corner) THEN
          j = npy
          ke_tl(1, j) = dt6*((ut_tl(1, j)+ut_tl(1, j-1))*u(1, j)+(ut(1, &
&           j)+ut(1, j-1))*u_tl(1, j)+(vt_tl(1, j)+vt_tl(0, j))*v(1, j-1&
&           )+(vt(1, j)+vt(0, j))*v_tl(1, j-1)+(ut_tl(1, j-1)-vt_tl(1, j&
&           ))*u(0, j)+(ut(1, j-1)-vt(1, j))*u_tl(0, j))
          ke(1, j) = dt6*((ut(1, j)+ut(1, j-1))*u(1, j)+(vt(1, j)+vt(0, &
&           j))*v(1, j-1)+(ut(1, j-1)-vt(1, j))*u(0, j))
        END IF
      ELSE IF (grid_type .LT. 3) THEN
         call mp_corner_comm(ke, npx, npy) 
         call mp_corner_comm(ke_tl, npx, npy) 
        IF (sw_corner) THEN
          ke_tl(1, 1) = r3*(ke_tl(2, 1)+ke_tl(1, 2)+ke_tl(0, 1))
          ke(1, 1) = r3*(ke(2, 1)+ke(1, 2)+ke(0, 1))
        END IF
        IF (se_corner) THEN
          ke_tl(npx, 1) = r3*(ke_tl(npx+1, 1)+ke_tl(npx, 2)+ke_tl(npx-1&
&           , 1))
          ke(npx, 1) = r3*(ke(npx+1, 1)+ke(npx, 2)+ke(npx-1, 1))
        END IF
        IF (ne_corner) THEN
          ke_tl(npx, npy) = r3*(ke_tl(npx+1, npy)+ke_tl(npx, npy-1)+&
&           ke_tl(npx-1, npy))
          ke(npx, npy) = r3*(ke(npx+1, npy)+ke(npx, npy-1)+ke(npx-1, npy&
&           ))
        END IF
        IF (nw_corner) THEN
          ke_tl(1, npy) = r3*(ke_tl(2, npy)+ke_tl(1, npy-1)+ke_tl(0, npy&
&           ))
          ke(1, npy) = r3*(ke(2, npy)+ke(1, npy-1)+ke(0, npy))
        END IF
      END IF
! Compute vorticity:
      DO j=jsd,jed+1
        DO i=isd,ied
          vt_tl(i, j) = dx(i, j)*u_tl(i, j)
          vt(i, j) = u(i, j)*dx(i, j)
        END DO
      END DO
      DO j=jsd,jed
        DO i=isd,ied+1
          ut_tl(i, j) = dy(i, j)*v_tl(i, j)
          ut(i, j) = v(i, j)*dy(i, j)
        END DO
      END DO
! wk is "volume-mean" relative vorticity
      DO j=jsd,jed
        DO i=isd,ied
          wk_tl(i, j) = rarea(i, j)*(vt_tl(i, j)-vt_tl(i, j+1)-ut_tl(i, &
&           j)+ut_tl(i+1, j))
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
              ptc_tl(i, j) = dyc(i, j)*sina_v(i, j)*u_tl(i, j)
              ptc(i, j) = u(i, j)*dyc(i, j)*sina_v(i, j)
            END DO
          ELSE
            DO i=is-1,ie+1
              ptc_tl(i, j) = dyc(i, j)*sina_v(i, j)*(u_tl(i, j)-0.5*&
&               cosa_v(i, j)*(va_tl(i, j-1)+va_tl(i, j)))
              ptc(i, j) = (u(i, j)-0.5*(va(i, j-1)+va(i, j))*cosa_v(i, j&
&               ))*dyc(i, j)*sina_v(i, j)
            END DO
          END IF
        END DO
        vort_tl = 0.0_8
        DO j=js-1,je+1
          DO i=is2,ie1
            vort_tl(i, j) = dxc(i, j)*sina_u(i, j)*(v_tl(i, j)-0.5*&
&             cosa_u(i, j)*(ua_tl(i-1, j)+ua_tl(i, j)))
            vort(i, j) = (v(i, j)-0.5*(ua(i-1, j)+ua(i, j))*cosa_u(i, j)&
&             )*dxc(i, j)*sina_u(i, j)
          END DO
          IF (is .EQ. 1) THEN
            vort_tl(1, j) = dxc(1, j)*sina_u(1, j)*v_tl(1, j)
            vort(1, j) = v(1, j)*dxc(1, j)*sina_u(1, j)
          END IF
          IF (ie + 1 .EQ. npx) THEN
            vort_tl(npx, j) = dxc(npx, j)*sina_u(npx, j)*v_tl(npx, j)
            vort(npx, j) = v(npx, j)*dxc(npx, j)*sina_u(npx, j)
          END IF
        END DO
        DO j=js,je+1
          DO i=is,ie+1
            delpc_tl(i, j) = vort_tl(i, j-1) - vort_tl(i, j) + ptc_tl(i-&
&             1, j) - ptc_tl(i, j)
            delpc(i, j) = vort(i, j-1) - vort(i, j) + ptc(i-1, j) - ptc(&
&             i, j)
          END DO
        END DO
! Remove the extra term at the corners:
        IF (sw_corner) THEN
          delpc_tl(1, 1) = delpc_tl(1, 1) - vort_tl(1, 0)
          delpc(1, 1) = delpc(1, 1) - vort(1, 0)
        END IF
        IF (se_corner) THEN
          delpc_tl(npx, 1) = delpc_tl(npx, 1) - vort_tl(npx, 0)
          delpc(npx, 1) = delpc(npx, 1) - vort(npx, 0)
        END IF
        IF (ne_corner) THEN
          delpc_tl(npx, npy) = delpc_tl(npx, npy) + vort_tl(npx, npy)
          delpc(npx, npy) = delpc(npx, npy) + vort(npx, npy)
        END IF
        IF (nw_corner) THEN
          delpc_tl(1, npy) = delpc_tl(1, npy) + vort_tl(1, npy)
          delpc(1, npy) = delpc(1, npy) + vort(1, npy)
        END IF
        DO j=js,je+1
          DO i=is,ie+1
            delpc_tl(i, j) = rarea_c(i, j)*delpc_tl(i, j)
            delpc(i, j) = rarea_c(i, j)*delpc(i, j)
            IF (delpc(i, j)*dt .GE. 0.) THEN
              abs0_tl = dt*delpc_tl(i, j)
              abs0 = delpc(i, j)*dt
            ELSE
              abs0_tl = -(dt*delpc_tl(i, j))
              abs0 = -(delpc(i, j)*dt)
            END IF
            y3_tl = dddmp*abs0_tl
            y3 = dddmp*abs0
            IF (0.20 .GT. y3) THEN
              y1_tl = y3_tl
              y1 = y3
            ELSE
              y1 = 0.20
              y1_tl = 0.0_8
            END IF
            IF (d2_bg .LT. y1) THEN
              max5_tl = y1_tl
              max5 = y1
            ELSE
              max5 = d2_bg
              max5_tl = 0.0_8
            END IF
            damp_tl = da_min_c*max5_tl
            damp = da_min_c*max5
            vort_tl(i, j) = damp_tl*delpc(i, j) + damp*delpc_tl(i, j)
            vort(i, j) = damp*delpc(i, j)
            ke_tl(i, j) = ke_tl(i, j) + vort_tl(i, j)
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
            delpc_tl(i, j) = divg_d_tl(i, j)
            delpc(i, j) = divg_d(i, j)
          END DO
        END DO
        n2 = nord + 1
        DO n=1,nord
          nt = nord - n
          fill_c = nt .NE. 0 .AND. grid_type .LT. 3 .AND. (((sw_corner &
&           .OR. se_corner) .OR. ne_corner) .OR. nw_corner)
          if ( fill_c ) call fill_corners(divg_d, npx, npy, FILL=XDir, BGRID=.true.)
          if ( fill_c ) call fill_corners(divg_d_tl, npx, npy, FILL=XDir, BGRID=.true.)
          DO j=js-nt,je+1+nt
            DO i=is-1-nt,ie+1+nt
              vc_tl(i, j) = divg_u(i, j)*(divg_d_tl(i+1, j)-divg_d_tl(i&
&               , j))
              vc(i, j) = (divg_d(i+1, j)-divg_d(i, j))*divg_u(i, j)
            END DO
          END DO
          if ( fill_c ) call fill_corners(divg_d, npx, npy, FILL=YDir, BGRID=.true.)
          if ( fill_c ) call fill_corners(divg_d_tl, npx, npy, FILL=YDir, BGRID=.true.)
          DO j=js-1-nt,je+1+nt
            DO i=is-nt,ie+1+nt
              uc_tl(i, j) = divg_v(i, j)*(divg_d_tl(i, j+1)-divg_d_tl(i&
&               , j))
              uc(i, j) = (divg_d(i, j+1)-divg_d(i, j))*divg_v(i, j)
            END DO
          END DO
          if ( fill_c ) call fill_corners(vc, uc, npx, npy, VECTOR=.true., DGRID=.true.)
          if ( fill_c ) call fill_corners(vc_tl, uc_tl, npx, npy, VECTOR=.true., DGRID=.true.)
          DO j=js-nt,je+1+nt
            DO i=is-nt,ie+1+nt
              divg_d_tl(i, j) = uc_tl(i, j-1) - uc_tl(i, j) + vc_tl(i-1&
&               , j) - vc_tl(i, j)
              divg_d(i, j) = uc(i, j-1) - uc(i, j) + vc(i-1, j) - vc(i, &
&               j)
            END DO
          END DO
! Remove the extra term at the corners:
          IF (sw_corner) THEN
            divg_d_tl(1, 1) = divg_d_tl(1, 1) - uc_tl(1, 0)
            divg_d(1, 1) = divg_d(1, 1) - uc(1, 0)
          END IF
          IF (se_corner) THEN
            divg_d_tl(npx, 1) = divg_d_tl(npx, 1) - uc_tl(npx, 0)
            divg_d(npx, 1) = divg_d(npx, 1) - uc(npx, 0)
          END IF
          IF (ne_corner) THEN
            divg_d_tl(npx, npy) = divg_d_tl(npx, npy) + uc_tl(npx, npy)
            divg_d(npx, npy) = divg_d(npx, npy) + uc(npx, npy)
          END IF
          IF (nw_corner) THEN
            divg_d_tl(1, npy) = divg_d_tl(1, npy) + uc_tl(1, npy)
            divg_d(1, npy) = divg_d(1, npy) + uc(1, npy)
          END IF
          DO j=js-nt,je+1+nt
            DO i=is-nt,ie+1+nt
              divg_d_tl(i, j) = rarea_c(i, j)*divg_d_tl(i, j)
              divg_d(i, j) = divg_d(i, j)*rarea_c(i, j)
            END DO
          END DO
        END DO
        IF (dddmp .LT. 1.e-5) THEN
          vort = 0.
          vort_tl = 0.0_8
        ELSE
          vort_tl = 0.0_8
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
                vort_tl(i, j) = dt*delpc_tl(i, j)
                vort(i, j) = dt*delpc(i, j)
              ELSE
                vort_tl(i, j) = -(dt*delpc_tl(i, j))
                vort(i, j) = -(dt*delpc(i, j))
              END IF
            END DO
          END DO
        END IF
!#endif
        pwx1 = da_min_c*d4_bg
        dd8 = pwx1**n2
        DO j=js,je+1
          DO i=is,ie+1
            IF (0.20 .GT. dddmp*vort(i, j)) THEN
              y2_tl = dddmp*vort_tl(i, j)
              y2 = dddmp*vort(i, j)
            ELSE
              y2 = 0.20
              y2_tl = 0.0_8
            END IF
            IF (d2_bg .LT. y2) THEN
              max6_tl = y2_tl
              max6 = y2
            ELSE
              max6 = d2_bg
              max6_tl = 0.0_8
            END IF
! del-2
            damp2_tl = da_min_c*max6_tl
            damp2 = da_min_c*max6
            vort_tl(i, j) = damp2_tl*delpc(i, j) + damp2*delpc_tl(i, j) &
&             + dd8*divg_d_tl(i, j)
            vort(i, j) = damp2*delpc(i, j) + dd8*divg_d(i, j)
            ke_tl(i, j) = ke_tl(i, j) + vort_tl(i, j)
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
        gy_tl = 0.0_8
        DO j=js,je+1
          DO i=is,ie
! du
            ub_tl(i, j) = rdx(i, j)*(vort_tl(i, j)-vort_tl(i+1, j))
            ub(i, j) = (vort(i, j)-vort(i+1, j))*rdx(i, j)
! u*du
            gy_tl(i, j) = u_tl(i, j)*ub(i, j) + u(i, j)*ub_tl(i, j)
            gy(i, j) = u(i, j)*ub(i, j)
          END DO
        END DO
        gx_tl = 0.0_8
        DO j=js,je
          DO i=is,ie+1
! dv
            vb_tl(i, j) = rdy(i, j)*(vort_tl(i, j)-vort_tl(i, j+1))
            vb(i, j) = (vort(i, j)-vort(i, j+1))*rdy(i, j)
! v*dv
            gx_tl(i, j) = v_tl(i, j)*vb(i, j) + v(i, j)*vb_tl(i, j)
            gx(i, j) = v(i, j)*vb(i, j)
          END DO
        END DO
        DO j=js,je
          DO i=is,ie
            u2_tl = u_tl(i, j) + u_tl(i, j+1)
            u2 = u(i, j) + u(i, j+1)
            du2_tl = ub_tl(i, j) + ub_tl(i, j+1)
            du2 = ub(i, j) + ub(i, j+1)
            v2_tl = v_tl(i, j) + v_tl(i+1, j)
            v2 = v(i, j) + v(i+1, j)
            dv2_tl = vb_tl(i, j) + vb_tl(i+1, j)
            dv2 = vb(i, j) + vb(i+1, j)
! Total energy conserving:
! Convert lost KE due to divergence damping to "heat"
            pt_tl(i, j) = pt_tl(i, j) - damp*rsin2(i, j)*(2*ub(i, j)*&
&             ub_tl(i, j)+2*ub(i, j+1)*ub_tl(i, j+1)+2*vb(i, j)*vb_tl(i&
&             , j)+2*vb(i+1, j)*vb_tl(i+1, j)+2.*(gy_tl(i, j)+gy_tl(i, j&
&             +1)+gx_tl(i, j)+gx_tl(i+1, j))-cosa_s(i, j)*(u2_tl*dv2+u2*&
&             dv2_tl+v2_tl*du2+v2*du2_tl+du2_tl*dv2+du2*dv2_tl))/pkz(i, &
&             j) + damp*rsin2(i, j)*pkz_tl(i, j)*(ub(i, j)**2+ub(i, j+1)&
&             **2+vb(i, j)**2+vb(i+1, j)**2+2.*(gy(i, j)+gy(i, j+1)+gx(i&
&             , j)+gx(i+1, j))-cosa_s(i, j)*(u2*dv2+v2*du2+du2*dv2))/pkz&
&             (i, j)**2
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
          vort_tl(i, j) = wk_tl(i, j)
          vort(i, j) = wk(i, j) + f0(i, j)
        END DO
      END DO
      fx_tj = fx
      fy_tj = fy
      vort_tj = vort
      CALL FV_TP_2D(vort,  crx_adv,  cry_adv, &
&                 npx, npy, hord_vt, fx, fy, &
&                 xfx_adv, yfx_adv, area, ra_x, &
&                 ra_y, ppm_fac=ppm_limiter, nord=nord&
&                 , damp_c=vtdm4)
      CALL FV_TP_2D_TLM(vort_tj, vort_tl, crx_adv, crx_adv_tl, cry_adv, &
&                 cry_adv_tl, npx, npy, hord_vt_pert, fx_tj, fx_tl, fy_tj, fy_tl, &
&                 xfx_adv, xfx_adv_tl, yfx_adv, yfx_adv_tl, area, ra_x, &
&                 ra_x_tl, ra_y, ra_y_tl, ppm_fac=ppm_limiter, nord=nord&
&                 , damp_c=vtdm4)
      DO j=js,je+1
        DO i=is,ie
          u_tl(i, j) = vt_tl(i, j) + ke_tl(i, j) - ke_tl(i+1, j) + fy_tl&
&           (i, j)
          u(i, j) = vt(i, j) + ke(i, j) - ke(i+1, j) + fy(i, j)
        END DO
      END DO
      DO j=js,je
        DO i=is,ie+1
          v_tl(i, j) = ut_tl(i, j) + ke_tl(i, j) - ke_tl(i, j+1) - fx_tl&
&           (i, j)
          v(i, j) = ut(i, j) + ke(i, j) - ke(i, j+1) - fx(i, j)
        END DO
      END DO
    END IF
  END SUBROUTINE D_SW_TLM
!  Differentiation of divergence_corner in forward (tangent) mode (with options r8):
!   variations   of useful results: divg_d
!   with respect to varying inputs: u v ua va divg_d
  SUBROUTINE DIVERGENCE_CORNER_TLM(u, u_tl, v, v_tl, ua, ua_tl, va, &
&   va_tl, divg_d, divg_d_tl, km)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: km
    REAL, DIMENSION(isd:ied, jsd:jed + 1, km), INTENT(IN) :: u
    REAL, DIMENSION(isd:ied, jsd:jed+1, km), INTENT(IN) :: u_tl
    REAL, DIMENSION(isd:ied + 1, jsd:jed, km), INTENT(IN) :: v
    REAL, DIMENSION(isd:ied+1, jsd:jed, km), INTENT(IN) :: v_tl
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(IN) :: ua, va
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(IN) :: ua_tl, va_tl
    REAL, DIMENSION(isd:ied + 1, jsd:jed + 1, km), INTENT(OUT) :: divg_d
    REAL, DIMENSION(isd:ied+1, jsd:jed+1, km), INTENT(OUT) :: divg_d_tl
! local
    REAL :: uf(is-2:ie+2, js-1:je+2)
    REAL :: uf_tl(is-2:ie+2, js-1:je+2)
    REAL :: vf(is-1:ie+2, js-2:je+2)
    REAL :: vf_tl(is-1:ie+2, js-2:je+2)
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
    uf_tl = 0.0_8
    vf_tl = 0.0_8
    DO k=1,km
      IF (grid_type .EQ. 4) THEN
        DO j=js-1,je+2
          DO i=is-2,ie+2
            uf_tl(i, j) = dyc(i, j)*u_tl(i, j, k)
            uf(i, j) = u(i, j, k)*dyc(i, j)
          END DO
        END DO
        DO j=js-2,je+2
          DO i=is-1,ie+2
            vf_tl(i, j) = dxc(i, j)*v_tl(i, j, k)
            vf(i, j) = v(i, j, k)*dxc(i, j)
          END DO
        END DO
        DO j=js-1,je+2
          DO i=is-1,ie+2
            divg_d_tl(i, j, k) = rarea_c(i, j)*(vf_tl(i, j-1)-vf_tl(i, j&
&             )+uf_tl(i-1, j)-uf_tl(i, j))
            divg_d(i, j, k) = rarea_c(i, j)*(vf(i, j-1)-vf(i, j)+uf(i-1&
&             , j)-uf(i, j))
          END DO
        END DO
      ELSE
! divg_u(i,j) = sina_v(i,j)*dyc(i,j)/dx(i,j)
        DO j=js,je+1
          IF (j .EQ. 1 .OR. j .EQ. npy) THEN
            DO i=is-1,ie+1
              uf_tl(i, j) = dyc(i, j)*sina_v(i, j)*u_tl(i, j, k)
              uf(i, j) = u(i, j, k)*dyc(i, j)*sina_v(i, j)
            END DO
          ELSE
            DO i=is-1,ie+1
              uf_tl(i, j) = dyc(i, j)*sina_v(i, j)*(u_tl(i, j, k)-0.5*&
&               cosa_v(i, j)*(va_tl(i, j-1, k)+va_tl(i, j, k)))
              uf(i, j) = (u(i, j, k)-0.5*(va(i, j-1, k)+va(i, j, k))*&
&               cosa_v(i, j))*dyc(i, j)*sina_v(i, j)
            END DO
          END IF
        END DO
        DO j=js-1,je+1
          DO i=is2,ie1
            vf_tl(i, j) = dxc(i, j)*sina_u(i, j)*(v_tl(i, j, k)-0.5*&
&             cosa_u(i, j)*(ua_tl(i-1, j, k)+ua_tl(i, j, k)))
            vf(i, j) = (v(i, j, k)-0.5*(ua(i-1, j, k)+ua(i, j, k))*&
&             cosa_u(i, j))*dxc(i, j)*sina_u(i, j)
          END DO
          IF (is .EQ. 1) THEN
            vf_tl(1, j) = dxc(1, j)*sina_u(1, j)*v_tl(1, j, k)
            vf(1, j) = v(1, j, k)*dxc(1, j)*sina_u(1, j)
          END IF
          IF (ie + 1 .EQ. npx) THEN
            vf_tl(npx, j) = dxc(npx, j)*sina_u(npx, j)*v_tl(npx, j, k)
            vf(npx, j) = v(npx, j, k)*dxc(npx, j)*sina_u(npx, j)
          END IF
        END DO
        DO j=js,je+1
          DO i=is,ie+1
            divg_d_tl(i, j, k) = vf_tl(i, j-1) - vf_tl(i, j) + uf_tl(i-1&
&             , j) - uf_tl(i, j)
            divg_d(i, j, k) = vf(i, j-1) - vf(i, j) + uf(i-1, j) - uf(i&
&             , j)
          END DO
        END DO
! Remove the extra term at the corners:
        IF (sw_corner) THEN
          divg_d_tl(1, 1, k) = divg_d_tl(1, 1, k) - vf_tl(1, 0)
          divg_d(1, 1, k) = divg_d(1, 1, k) - vf(1, 0)
        END IF
        IF (se_corner) THEN
          divg_d_tl(npx, 1, k) = divg_d_tl(npx, 1, k) - vf_tl(npx, 0)
          divg_d(npx, 1, k) = divg_d(npx, 1, k) - vf(npx, 0)
        END IF
        IF (ne_corner) THEN
          divg_d_tl(npx, npy, k) = divg_d_tl(npx, npy, k) + vf_tl(npx, &
&           npy)
          divg_d(npx, npy, k) = divg_d(npx, npy, k) + vf(npx, npy)
        END IF
        IF (nw_corner) THEN
          divg_d_tl(1, npy, k) = divg_d_tl(1, npy, k) + vf_tl(1, npy)
          divg_d(1, npy, k) = divg_d(1, npy, k) + vf(1, npy)
        END IF
        DO j=js,je+1
          DO i=is,ie+1
            divg_d_tl(i, j, k) = rarea_c(i, j)*divg_d_tl(i, j, k)
            divg_d(i, j, k) = rarea_c(i, j)*divg_d(i, j, k)
          END DO
        END DO
      END IF
    END DO
  END SUBROUTINE DIVERGENCE_CORNER_TLM
!  Differentiation of xtp_u in forward (tangent) mode (with options r8):
!   variations   of useful results: flux
!   with respect to varying inputs: flux u c
  SUBROUTINE XTP_U_TLM(c, c_tl, u, u_tl, v, flux, flux_tl, iord)
    IMPLICIT NONE
    REAL, INTENT(IN) :: u(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: u_tl(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: v(isd:ied+1, jsd:jed)
    REAL(ke_precision), INTENT(IN) :: c(is:ie+1, js:je+1)
    REAL(ke_precision), INTENT(IN) :: c_tl(is:ie+1, js:je+1)
    REAL(ke_precision), INTENT(OUT) :: flux(is:ie+1, js:je+1)
    REAL(ke_precision), INTENT(OUT) :: flux_tl(is:ie+1, js:je+1)
    INTEGER, INTENT(IN) :: iord
! Local
    REAL(ke_precision) :: al(is-1:ie+2), dm(is-2:ie+2)
    REAL(ke_precision) :: bl(is-1:ie+1)
    REAL(ke_precision) :: bl_tl(is-1:ie+1)
    REAL(ke_precision) :: br(is-1:ie+1)
    REAL(ke_precision) :: br_tl(is-1:ie+1)
    REAL(ke_precision) :: dq(is-3:ie+2)
    REAL(ke_precision) :: dl, dr, xt, pmp, lac, dqt, cfl
    REAL(ke_precision) :: xt_tl, cfl_tl
    REAL(ke_precision) :: x0, x1
    INTEGER :: i, j
    INTRINSIC MAX
    INTRINSIC MIN
    INTEGER :: min1
    INTEGER :: max1
    SELECT CASE  (iord) 
    CASE (1) 
      DO j=js,je+1
        DO i=is,ie+1
          IF (c(i, j) .GT. 0.) THEN
            flux_tl(i, j) = u_tl(i-1, j)
            flux(i, j) = u(i-1, j)
          ELSE
            flux_tl(i, j) = u_tl(i, j)
            flux(i, j) = u(i, j)
          END IF
        END DO
      END DO
    CASE (333) 
      DO j=js,je+1
        DO i=is,ie+1
!Apply first order at the edges
          IF (i .EQ. is .OR. i .EQ. ie + 1) THEN
            IF (c(i, j) .GT. 0.) THEN
              flux_tl(i, j) = u_tl(i-1, j)
              flux(i, j) = u(i-1, j)
            ELSE
              flux_tl(i, j) = u_tl(i, j)
              flux(i, j) = u(i, j)
            END IF
          ELSE IF (c(i, j) .GT. 0.) THEN
!Otherwise use the third order scheme
            cfl_tl = rdx(i-1, j)*c_tl(i, j)
            cfl = c(i, j)*rdx(i-1, j)
            flux_tl(i, j) = (2.0*u_tl(i, j)+5.0*u_tl(i-1, j)-u_tl(i-2, j&
&             ))/6.0 - 0.5*(cfl_tl*(u(i, j)-u(i-1, j))+cfl*(u_tl(i, j)-&
&             u_tl(i-1, j))) + (cfl_tl*cfl+cfl*cfl_tl)*(u(i, j)-2.0*u(i-&
&             1, j)+u(i-2, j))/6.0 + cfl**2*(u_tl(i, j)-2.0*u_tl(i-1, j)&
&             +u_tl(i-2, j))/6.0
            flux(i, j) = (2.0*u(i, j)+5.0*u(i-1, j)-u(i-2, j))/6.0 - 0.5&
&             *cfl*(u(i, j)-u(i-1, j)) + cfl*cfl/6.0*(u(i, j)-2.0*u(i-1&
&             , j)+u(i-2, j))
          ELSE
            cfl_tl = rdx(i, j)*c_tl(i, j)
            cfl = c(i, j)*rdx(i, j)
            flux_tl(i, j) = (2.0*u_tl(i-1, j)+5.0*u_tl(i, j)-u_tl(i+1, j&
&             ))/6.0 - 0.5*(cfl_tl*(u(i, j)-u(i-1, j))+cfl*(u_tl(i, j)-&
&             u_tl(i-1, j))) + (cfl_tl*cfl+cfl*cfl_tl)*(u(i+1, j)-2.0*u(&
&             i, j)+u(i-1, j))/6.0 + cfl**2*(u_tl(i+1, j)-2.0*u_tl(i, j)&
&             +u_tl(i-1, j))/6.0
            flux(i, j) = (2.0*u(i-1, j)+5.0*u(i, j)-u(i+1, j))/6.0 - 0.5&
&             *cfl*(u(i, j)-u(i-1, j)) + cfl*cfl/6.0*(u(i+1, j)-2.0*u(i&
&             , j)+u(i-1, j))
          END IF
        END DO
      END DO
    CASE (6) 
      bl_tl = 0.0_8
      br_tl = 0.0_8
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
        DO i=max1,min1
          bl_tl(i) = b5*u_tl(i-2, j) + b4*u_tl(i-1, j) + b3*u_tl(i, j) +&
&           b2*u_tl(i+1, j) + b1*u_tl(i+2, j)
          bl(i) = b5*u(i-2, j) + b4*u(i-1, j) + b3*u(i, j) + b2*u(i+1, j&
&           ) + b1*u(i+2, j)
          br_tl(i) = b1*u_tl(i-2, j) + b2*u_tl(i-1, j) + b3*u_tl(i, j) +&
&           b4*u_tl(i+1, j) + b5*u_tl(i+2, j)
          br(i) = b1*u(i-2, j) + b2*u(i-1, j) + b3*u(i, j) + b4*u(i+1, j&
&           ) + b5*u(i+2, j)
        END DO
        IF (grid_type .LT. 3) THEN
          IF (is .EQ. 1) THEN
            br_tl(2) = p1*(u_tl(2, j)+u_tl(3, j)) + p2*(u_tl(1, j)+u_tl(&
&             4, j)) - u_tl(2, j)
            br(2) = p1*(u(2, j)+u(3, j)) + p2*(u(1, j)+u(4, j)) - u(2, j&
&             )
            xt_tl = c3*u_tl(1, j) + c2*u_tl(2, j) + c1*u_tl(3, j)
            xt = c3*u(1, j) + c2*u(2, j) + c1*u(3, j)
            bl_tl(2) = xt_tl - u_tl(2, j)
            bl(2) = xt - u(2, j)
            IF (j .EQ. 1 .OR. j .EQ. npy) THEN
! out
              bl_tl(0) = 0.0_8
              bl(0) = 0.
! edge
              br_tl(0) = 0.0_8
              br(0) = 0.
! edge
              bl_tl(1) = 0.0_8
              bl(1) = 0.
! in
              br_tl(1) = 0.0_8
              br(1) = 0.
            ELSE
              br_tl(1) = xt_tl - u_tl(1, j)
              br(1) = xt - u(1, j)
              xt_tl = 0.5*((2.*dx(1, j)+dx(2, j))*(u_tl(0, j)+u_tl(1, j)&
&               )-dx(1, j)*(u_tl(-1, j)+u_tl(2, j)))/(dx(1, j)+dx(2, j))
              xt = 0.5*((2.*dx(1, j)+dx(2, j))*(u(0, j)+u(1, j))-dx(1, j&
&               )*(u(-1, j)+u(2, j)))/(dx(1, j)+dx(2, j))
              bl_tl(1) = xt_tl - u_tl(1, j)
              bl(1) = xt - u(1, j)
              br_tl(0) = xt_tl - u_tl(0, j)
              br(0) = xt - u(0, j)
              xt_tl = c1*u_tl(-2, j) + c2*u_tl(-1, j) + c3*u_tl(0, j)
              xt = c1*u(-2, j) + c2*u(-1, j) + c3*u(0, j)
              bl_tl(0) = xt_tl - u_tl(0, j)
              bl(0) = xt - u(0, j)
            END IF
          END IF
          IF (ie + 1 .EQ. npx) THEN
            bl_tl(npx-2) = p1*(u_tl(npx-2, j)+u_tl(npx-3, j)) + p2*(u_tl&
&             (npx-4, j)+u_tl(npx-1, j)) - u_tl(npx-2, j)
            bl(npx-2) = p1*(u(npx-2, j)+u(npx-3, j)) + p2*(u(npx-4, j)+u&
&             (npx-1, j)) - u(npx-2, j)
            xt_tl = c1*u_tl(npx-3, j) + c2*u_tl(npx-2, j) + c3*u_tl(npx-&
&             1, j)
            xt = c1*u(npx-3, j) + c2*u(npx-2, j) + c3*u(npx-1, j)
            br_tl(npx-2) = xt_tl - u_tl(npx-2, j)
            br(npx-2) = xt - u(npx-2, j)
            IF (j .EQ. 1 .OR. j .EQ. npy) THEN
! in
              bl_tl(npx-1) = 0.0_8
              bl(npx-1) = 0.
! edge
              br_tl(npx-1) = 0.0_8
              br(npx-1) = 0.
! edge
              bl_tl(npx) = 0.0_8
              bl(npx) = 0.
! out
              br_tl(npx) = 0.0_8
              br(npx) = 0.
            ELSE
              bl_tl(npx-1) = xt_tl - u_tl(npx-1, j)
              bl(npx-1) = xt - u(npx-1, j)
              xt_tl = 0.5*((2.*dx(npx-1, j)+dx(npx-2, j))*(u_tl(npx-1, j&
&               )+u_tl(npx, j))-dx(npx-1, j)*(u_tl(npx-2, j)+u_tl(npx+1&
&               , j)))/(dx(npx-1, j)+dx(npx-2, j))
              xt = 0.5*((2.*dx(npx-1, j)+dx(npx-2, j))*(u(npx-1, j)+u(&
&               npx, j))-dx(npx-1, j)*(u(npx-2, j)+u(npx+1, j)))/(dx(npx&
&               -1, j)+dx(npx-2, j))
              br_tl(npx-1) = xt_tl - u_tl(npx-1, j)
              br(npx-1) = xt - u(npx-1, j)
              bl_tl(npx) = xt_tl - u_tl(npx, j)
              bl(npx) = xt - u(npx, j)
              xt_tl = c3*u_tl(npx, j) + c2*u_tl(npx+1, j) + c1*u_tl(npx+&
&               2, j)
              xt = c3*u(npx, j) + c2*u(npx+1, j) + c1*u(npx+2, j)
              br_tl(npx) = xt_tl - u_tl(npx, j)
              br(npx) = xt - u(npx, j)
            END IF
          END IF
        END IF
        DO i=is,ie+1
          IF (c(i, j) .GT. 0.) THEN
            cfl_tl = rdx(i-1, j)*c_tl(i, j)
            cfl = c(i, j)*rdx(i-1, j)
            flux_tl(i, j) = u_tl(i-1, j) + (1.-cfl)*(br_tl(i-1)-cfl_tl*(&
&             bl(i-1)+br(i-1))-cfl*(bl_tl(i-1)+br_tl(i-1))) - cfl_tl*(br&
&             (i-1)-cfl*(bl(i-1)+br(i-1)))
            flux(i, j) = u(i-1, j) + (1.-cfl)*(br(i-1)-cfl*(bl(i-1)+br(i&
&             -1)))
          ELSE
            cfl_tl = rdx(i, j)*c_tl(i, j)
            cfl = c(i, j)*rdx(i, j)
            flux_tl(i, j) = u_tl(i, j) + cfl_tl*(bl(i)+cfl*(bl(i)+br(i))&
&             ) + (1.+cfl)*(bl_tl(i)+cfl_tl*(bl(i)+br(i))+cfl*(bl_tl(i)+&
&             br_tl(i)))
            flux(i, j) = u(i, j) + (1.+cfl)*(bl(i)+cfl*(bl(i)+br(i)))
          END IF
        END DO
      END DO
    END SELECT
  END SUBROUTINE XTP_U_TLM
!  Differentiation of ytp_v in forward (tangent) mode (with options r8):
!   variations   of useful results: flux
!   with respect to varying inputs: v c
  SUBROUTINE YTP_V_TLM(c, c_tl, u, v, v_tl, flux, flux_tl, jord)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: jord
    REAL, INTENT(IN) :: u(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: v(isd:ied+1, jsd:jed)
    REAL, INTENT(IN) :: v_tl(isd:ied+1, jsd:jed)
!  Courant   N (like FLUX)
    REAL(ke_precision), INTENT(IN) :: c(is:ie+1, js:je+1)
    REAL(ke_precision), INTENT(IN) :: c_tl(is:ie+1, js:je+1)
    REAL(ke_precision), INTENT(OUT) :: flux(is:ie+1, js:je+1)
    REAL(ke_precision), INTENT(OUT) :: flux_tl(is:ie+1, js:je+1)
! Local:
    REAL(ke_precision) :: dm(is:ie+1, js-2:je+2)
    REAL(ke_precision) :: al(is:ie+1, js-1:je+2)
    REAL(ke_precision) :: bl(is:ie+1, js-1:je+1)
    REAL(ke_precision) :: bl_tl(is:ie+1, js-1:je+1)
    REAL(ke_precision) :: br(is:ie+1, js-1:je+1)
    REAL(ke_precision) :: br_tl(is:ie+1, js-1:je+1)
    REAL(ke_precision) :: dq(is:ie+1, js-3:je+2)
    REAL(ke_precision) :: xt, dl, dr, pmp, lac, dqt, cfl
    REAL(ke_precision) :: xt_tl, cfl_tl
    REAL(ke_precision) :: x0, x1
    INTEGER :: i, j
    INTRINSIC MAX
    INTRINSIC MIN
    INTEGER :: min1
    INTEGER :: max1
    SELECT CASE  (jord) 
    CASE (1) 
      flux_tl = 0.0_8
      DO j=js,je+1
        DO i=is,ie+1
          IF (c(i, j) .GT. 0.) THEN
            flux_tl(i, j) = v_tl(i, j-1)
            flux(i, j) = v(i, j-1)
          ELSE
            flux_tl(i, j) = v_tl(i, j)
            flux(i, j) = v(i, j)
          END IF
        END DO
      END DO
    CASE (333) 
      flux_tl = 0.0_8
      DO j=js,je+1
        DO i=is,ie+1
!Apply first order at the edges
          IF (j .EQ. js .OR. j .EQ. je + 1) THEN
            IF (c(i, j) .GT. 0.) THEN
              flux_tl(i, j) = v_tl(i, j-1)
              flux(i, j) = v(i, j-1)
            ELSE
              flux_tl(i, j) = v_tl(i, j)
              flux(i, j) = v(i, j)
            END IF
          ELSE IF (c(i, j) .GT. 0.) THEN
!Otherwise use the third order scheme
            cfl_tl = rdy(i, j-1)*c_tl(i, j)
            cfl = c(i, j)*rdy(i, j-1)
            flux_tl(i, j) = (2.0*v_tl(i, j)+5.0*v_tl(i, j-1)-v_tl(i, j-2&
&             ))/6.0 - 0.5*(cfl_tl*(v(i, j)-v(i, j-1))+cfl*(v_tl(i, j)-&
&             v_tl(i, j-1))) + (cfl_tl*cfl+cfl*cfl_tl)*(v(i, j)-2.0*v(i&
&             , j-1)+v(i, j-2))/6.0 + cfl**2*(v_tl(i, j)-2.0*v_tl(i, j-1&
&             )+v_tl(i, j-2))/6.0
            flux(i, j) = (2.0*v(i, j)+5.0*v(i, j-1)-v(i, j-2))/6.0 - 0.5&
&             *cfl*(v(i, j)-v(i, j-1)) + cfl*cfl/6.0*(v(i, j)-2.0*v(i, j&
&             -1)+v(i, j-2))
          ELSE
            cfl_tl = rdy(i, j)*c_tl(i, j)
            cfl = c(i, j)*rdy(i, j)
            flux_tl(i, j) = (2.0*v_tl(i, j-1)+5.0*v_tl(i, j)-v_tl(i, j+1&
&             ))/6.0 - 0.5*(cfl_tl*(v(i, j)-v(i, j-1))+cfl*(v_tl(i, j)-&
&             v_tl(i, j-1))) + (cfl_tl*cfl+cfl*cfl_tl)*(v(i, j+1)-2.0*v(&
&             i, j)+v(i, j-1))/6.0 + cfl**2*(v_tl(i, j+1)-2.0*v_tl(i, j)&
&             +v_tl(i, j-1))/6.0
            flux(i, j) = (2.0*v(i, j-1)+5.0*v(i, j)-v(i, j+1))/6.0 - 0.5&
&             *cfl*(v(i, j)-v(i, j-1)) + cfl*cfl/6.0*(v(i, j+1)-2.0*v(i&
&             , j)+v(i, j-1))
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
        bl_tl = 0.0_8
        br_tl = 0.0_8
      ELSE
        min1 = npy - 3
        bl_tl = 0.0_8
        br_tl = 0.0_8
      END IF
      DO j=max1,min1
        DO i=is,ie+1
          bl_tl(i, j) = b5*v_tl(i, j-2) + b4*v_tl(i, j-1) + b3*v_tl(i, j&
&           ) + b2*v_tl(i, j+1) + b1*v_tl(i, j+2)
          bl(i, j) = b5*v(i, j-2) + b4*v(i, j-1) + b3*v(i, j) + b2*v(i, &
&           j+1) + b1*v(i, j+2)
          br_tl(i, j) = b1*v_tl(i, j-2) + b2*v_tl(i, j-1) + b3*v_tl(i, j&
&           ) + b4*v_tl(i, j+1) + b5*v_tl(i, j+2)
          br(i, j) = b1*v(i, j-2) + b2*v(i, j-1) + b3*v(i, j) + b4*v(i, &
&           j+1) + b5*v(i, j+2)
        END DO
      END DO
      IF (grid_type .LT. 3) THEN
        IF (js .EQ. 1) THEN
          DO i=is,ie+1
            br_tl(i, 2) = p1*(v_tl(i, 2)+v_tl(i, 3)) + p2*(v_tl(i, 1)+&
&             v_tl(i, 4)) - v_tl(i, 2)
            br(i, 2) = p1*(v(i, 2)+v(i, 3)) + p2*(v(i, 1)+v(i, 4)) - v(i&
&             , 2)
            xt_tl = c3*v_tl(i, 1) + c2*v_tl(i, 2) + c1*v_tl(i, 3)
            xt = c3*v(i, 1) + c2*v(i, 2) + c1*v(i, 3)
            br_tl(i, 1) = xt_tl - v_tl(i, 1)
            br(i, 1) = xt - v(i, 1)
            bl_tl(i, 2) = xt_tl - v_tl(i, 2)
            bl(i, 2) = xt - v(i, 2)
            bl_tl(i, 0) = c1*v_tl(i, -2) + c2*v_tl(i, -1) + c3*v_tl(i, 0&
&             ) - v_tl(i, 0)
            bl(i, 0) = c1*v(i, -2) + c2*v(i, -1) + c3*v(i, 0) - v(i, 0)
            xt_tl = 0.5*((2.*dy(i, 1)+dy(i, 2))*(v_tl(i, 0)+v_tl(i, 1))-&
&             dy(i, 1)*(v_tl(i, -1)+v_tl(i, 2)))/(dy(i, 1)+dy(i, 2))
            xt = 0.5*((2.*dy(i, 1)+dy(i, 2))*(v(i, 0)+v(i, 1))-dy(i, 1)*&
&             (v(i, -1)+v(i, 2)))/(dy(i, 1)+dy(i, 2))
            bl_tl(i, 1) = xt_tl - v_tl(i, 1)
            bl(i, 1) = xt - v(i, 1)
            br_tl(i, 0) = xt_tl - v_tl(i, 0)
            br(i, 0) = xt - v(i, 0)
          END DO
          IF (is .EQ. 1) THEN
! out
            bl_tl(1, 0) = 0.0_8
            bl(1, 0) = 0.
! edge
            br_tl(1, 0) = 0.0_8
            br(1, 0) = 0.
! edge
            bl_tl(1, 1) = 0.0_8
            bl(1, 1) = 0.
! in
            br_tl(1, 1) = 0.0_8
            br(1, 1) = 0.
          END IF
          IF (ie + 1 .EQ. npx) THEN
! out
            bl_tl(npx, 0) = 0.0_8
            bl(npx, 0) = 0.
! edge
            br_tl(npx, 0) = 0.0_8
            br(npx, 0) = 0.
! edge
            bl_tl(npx, 1) = 0.0_8
            bl(npx, 1) = 0.
! in
            br_tl(npx, 1) = 0.0_8
            br(npx, 1) = 0.
          END IF
        END IF
        IF (je + 1 .EQ. npy) THEN
          DO i=is,ie+1
            bl_tl(i, npy-2) = p1*(v_tl(i, npy-3)+v_tl(i, npy-2)) + p2*(&
&             v_tl(i, npy-4)+v_tl(i, npy-1)) - v_tl(i, npy-2)
            bl(i, npy-2) = p1*(v(i, npy-3)+v(i, npy-2)) + p2*(v(i, npy-4&
&             )+v(i, npy-1)) - v(i, npy-2)
            xt_tl = c1*v_tl(i, npy-3) + c2*v_tl(i, npy-2) + c3*v_tl(i, &
&             npy-1)
            xt = c1*v(i, npy-3) + c2*v(i, npy-2) + c3*v(i, npy-1)
            br_tl(i, npy-2) = xt_tl - v_tl(i, npy-2)
            br(i, npy-2) = xt - v(i, npy-2)
            bl_tl(i, npy-1) = xt_tl - v_tl(i, npy-1)
            bl(i, npy-1) = xt - v(i, npy-1)
            br_tl(i, npy) = c3*v_tl(i, npy) + c2*v_tl(i, npy+1) + c1*&
&             v_tl(i, npy+2) - v_tl(i, npy)
            br(i, npy) = c3*v(i, npy) + c2*v(i, npy+1) + c1*v(i, npy+2) &
&             - v(i, npy)
            xt_tl = 0.5*((2.*dy(i, npy-1)+dy(i, npy-2))*(v_tl(i, npy-1)+&
&             v_tl(i, npy))-dy(i, npy-1)*(v_tl(i, npy-2)+v_tl(i, npy+1))&
&             )/(dy(i, npy-1)+dy(i, npy-2))
            xt = 0.5*((2.*dy(i, npy-1)+dy(i, npy-2))*(v(i, npy-1)+v(i, &
&             npy))-dy(i, npy-1)*(v(i, npy-2)+v(i, npy+1)))/(dy(i, npy-1&
&             )+dy(i, npy-2))
            br_tl(i, npy-1) = xt_tl - v_tl(i, npy-1)
            br(i, npy-1) = xt - v(i, npy-1)
            bl_tl(i, npy) = xt_tl - v_tl(i, npy)
            bl(i, npy) = xt - v(i, npy)
          END DO
          IF (is .EQ. 1) THEN
! in
            bl_tl(1, npy-1) = 0.0_8
            bl(1, npy-1) = 0.
! edge
            br_tl(1, npy-1) = 0.0_8
            br(1, npy-1) = 0.
! edge
            bl_tl(1, npy) = 0.0_8
            bl(1, npy) = 0.
! out
            br_tl(1, npy) = 0.0_8
            br(1, npy) = 0.
          END IF
          IF (ie + 1 .EQ. npx) THEN
! in
            bl_tl(npx, npy-1) = 0.0_8
            bl(npx, npy-1) = 0.
! edge
            br_tl(npx, npy-1) = 0.0_8
            br(npx, npy-1) = 0.
! edge
            bl_tl(npx, npy) = 0.0_8
            bl(npx, npy) = 0.
! out
            br_tl(npx, npy) = 0.0_8
            br(npx, npy) = 0.
            flux_tl = 0.0_8
          ELSE
            flux_tl = 0.0_8
          END IF
        ELSE
          flux_tl = 0.0_8
        END IF
      ELSE
        flux_tl = 0.0_8
      END IF
      DO j=js,je+1
        DO i=is,ie+1
          IF (c(i, j) .GT. 0.) THEN
            cfl_tl = rdy(i, j-1)*c_tl(i, j)
            cfl = c(i, j)*rdy(i, j-1)
            flux_tl(i, j) = v_tl(i, j-1) + (1.-cfl)*(br_tl(i, j-1)-&
&             cfl_tl*(bl(i, j-1)+br(i, j-1))-cfl*(bl_tl(i, j-1)+br_tl(i&
&             , j-1))) - cfl_tl*(br(i, j-1)-cfl*(bl(i, j-1)+br(i, j-1)))
            flux(i, j) = v(i, j-1) + (1.-cfl)*(br(i, j-1)-cfl*(bl(i, j-1&
&             )+br(i, j-1)))
          ELSE
            cfl_tl = rdy(i, j)*c_tl(i, j)
            cfl = c(i, j)*rdy(i, j)
            flux_tl(i, j) = v_tl(i, j) + cfl_tl*(bl(i, j)+cfl*(bl(i, j)+&
&             br(i, j))) + (1.+cfl)*(bl_tl(i, j)+cfl_tl*(bl(i, j)+br(i, &
&             j))+cfl*(bl_tl(i, j)+br_tl(i, j)))
            flux(i, j) = v(i, j) + (1.+cfl)*(bl(i, j)+cfl*(bl(i, j)+br(i&
&             , j)))
          END IF
        END DO
      END DO
    CASE DEFAULT
      flux_tl = 0.0_8
    END SELECT
  END SUBROUTINE YTP_V_TLM
!  Differentiation of d2a2c_vect in forward (tangent) mode (with options r8):
!   variations   of useful results: ua uc ut va vc vt
!   with respect to varying inputs: u v ua uc ut va vc vt
  SUBROUTINE D2A2C_VECT_TLM(u, u_tl, v, v_tl, ua, ua_tl, va, va_tl, uc, &
&   uc_tl, vc, vc_tl, ut, ut_tl, vt, vt_tl, dord4)
    IMPLICIT NONE
    REAL, INTENT(IN) :: u(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: u_tl(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: v(isd:ied+1, jsd:jed)
    REAL, INTENT(IN) :: v_tl(isd:ied+1, jsd:jed)
    LOGICAL, INTENT(IN) :: dord4
    REAL, DIMENSION(isd:ied + 1, jsd:jed), INTENT(OUT) :: uc
    REAL, DIMENSION(isd:ied+1, jsd:jed), INTENT(OUT) :: uc_tl
    REAL, DIMENSION(isd:ied, jsd:jed + 1), INTENT(OUT) :: vc
    REAL, DIMENSION(isd:ied, jsd:jed+1), INTENT(OUT) :: vc_tl
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(OUT) :: ua, va, ut, vt
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(OUT) :: ua_tl, va_tl, &
&   ut_tl, vt_tl
! Local 
    REAL, DIMENSION(isd:ied, jsd:jed) :: utmp, vtmp
    REAL, DIMENSION(isd:ied, jsd:jed) :: utmp_tl, vtmp_tl
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
      utmp_tl = 0.0_8
    ELSE
      min1 = npy - npt
      utmp_tl = 0.0_8
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
        utmp_tl(i, j) = a2*(u_tl(i, j-1)+u_tl(i, j+2)) + a1*(u_tl(i, j)+&
&         u_tl(i, j+1))
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
      vtmp_tl = 0.0_8
    ELSE
      min3 = npy - npt
      vtmp_tl = 0.0_8
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
        vtmp_tl(i, j) = a2*(v_tl(i-1, j)+v_tl(i+2, j)) + a1*(v_tl(i, j)+&
&         v_tl(i+1, j))
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
            utmp_tl(i, j) = 0.5*(u_tl(i, j)+u_tl(i, j+1))
            utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
            vtmp_tl(i, j) = 0.5*(v_tl(i, j)+v_tl(i+1, j))
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
            utmp_tl(i, j) = 0.5*(u_tl(i, j)+u_tl(i, j+1))
            utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
            vtmp_tl(i, j) = 0.5*(v_tl(i, j)+v_tl(i+1, j))
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
            utmp_tl(i, j) = 0.5*(u_tl(i, j)+u_tl(i, j+1))
            utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
            vtmp_tl(i, j) = 0.5*(v_tl(i, j)+v_tl(i+1, j))
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
            utmp_tl(i, j) = 0.5*(u_tl(i, j)+u_tl(i, j+1))
            utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
            vtmp_tl(i, j) = 0.5*(v_tl(i, j)+v_tl(i+1, j))
            vtmp(i, j) = 0.5*(v(i, j)+v(i+1, j))
          END DO
        END DO
      END IF
    END IF
!#endif
    DO j=js-1-id,je+1+id
      DO i=is-1-id,ie+1+id
        ua_tl(i, j) = rsin2(i, j)*(utmp_tl(i, j)-cosa_s(i, j)*vtmp_tl(i&
&         , j))
        ua(i, j) = (utmp(i, j)-vtmp(i, j)*cosa_s(i, j))*rsin2(i, j)
        va_tl(i, j) = rsin2(i, j)*(vtmp_tl(i, j)-cosa_s(i, j)*utmp_tl(i&
&         , j))
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
        utmp_tl(i, 0) = -vtmp_tl(0, 1-i)
        utmp(i, 0) = -vtmp(0, 1-i)
      END DO
    END IF
    IF (se_corner) THEN
      DO i=0,2
        utmp_tl(npx+i, 0) = vtmp_tl(npx, i+1)
        utmp(npx+i, 0) = vtmp(npx, i+1)
      END DO
    END IF
    IF (ne_corner) THEN
      DO i=0,2
        utmp_tl(npx+i, npy) = -vtmp_tl(npx, je-i)
        utmp(npx+i, npy) = -vtmp(npx, je-i)
      END DO
    END IF
    IF (nw_corner) THEN
      DO i=-2,0
        utmp_tl(i, npy) = vtmp_tl(0, je+i)
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
        uc_tl(i, j) = a1*(utmp_tl(i-1, j)+utmp_tl(i, j)) + a2*(utmp_tl(i&
&         -2, j)+utmp_tl(i+1, j))
        uc(i, j) = a1*(utmp(i-1, j)+utmp(i, j)) + a2*(utmp(i-2, j)+utmp(&
&         i+1, j))
        ut_tl(i, j) = rsin_u(i, j)*(uc_tl(i, j)-cosa_u(i, j)*v_tl(i, j))
        ut(i, j) = (uc(i, j)-v(i, j)*cosa_u(i, j))*rsin_u(i, j)
      END DO
    END DO
    IF (grid_type .LT. 3) THEN
      IF (is .EQ. 1) THEN
        DO j=js-1,je+1
          uc_tl(0, j) = c1*utmp_tl(-2, j) + c2*utmp_tl(-1, j) + c3*&
&           utmp_tl(0, j)
          uc(0, j) = c1*utmp(-2, j) + c2*utmp(-1, j) + c3*utmp(0, j)
! 3-pt extrapolation --------------------------------------------------
          uc_tl(1, j) = rsin_u(1, j)*(t14*(utmp_tl(0, j)+utmp_tl(1, j))+&
&           t12*(utmp_tl(-1, j)+utmp_tl(2, j))+t15*(utmp_tl(-2, j)+&
&           utmp_tl(3, j)))
          uc(1, j) = (t14*(utmp(0, j)+utmp(1, j))+t12*(utmp(-1, j)+utmp(&
&           2, j))+t15*(utmp(-2, j)+utmp(3, j)))*rsin_u(1, j)
          uc_tl(2, j) = c1*utmp_tl(3, j) + c2*utmp_tl(2, j) + c3*utmp_tl&
&           (1, j)
          uc(2, j) = c1*utmp(3, j) + c2*utmp(2, j) + c3*utmp(1, j)
          ut_tl(0, j) = rsin_u(0, j)*(uc_tl(0, j)-cosa_u(0, j)*v_tl(0, j&
&           ))
          ut(0, j) = (uc(0, j)-v(0, j)*cosa_u(0, j))*rsin_u(0, j)
          ut_tl(1, j) = rsin_u(1, j)*uc_tl(1, j)
          ut(1, j) = uc(1, j)*rsin_u(1, j)
          ut_tl(2, j) = rsin_u(2, j)*(uc_tl(2, j)-cosa_u(2, j)*v_tl(2, j&
&           ))
          ut(2, j) = (uc(2, j)-v(2, j)*cosa_u(2, j))*rsin_u(2, j)
        END DO
      END IF
      IF (ie + 1 .EQ. npx) THEN
        DO j=js-1,je+1
          uc_tl(npx-1, j) = c1*utmp_tl(npx-3, j) + c2*utmp_tl(npx-2, j) &
&           + c3*utmp_tl(npx-1, j)
          uc(npx-1, j) = c1*utmp(npx-3, j) + c2*utmp(npx-2, j) + c3*utmp&
&           (npx-1, j)
! 3-pt extrapolation --------------------------------------------------------
          uc_tl(npx, j) = rsin_u(npx, j)*(t14*(utmp_tl(npx-1, j)+utmp_tl&
&           (npx, j))+t12*(utmp_tl(npx-2, j)+utmp_tl(npx+1, j))+t15*(&
&           utmp_tl(npx-3, j)+utmp_tl(npx+2, j)))
          uc(npx, j) = (t14*(utmp(npx-1, j)+utmp(npx, j))+t12*(utmp(npx-&
&           2, j)+utmp(npx+1, j))+t15*(utmp(npx-3, j)+utmp(npx+2, j)))*&
&           rsin_u(npx, j)
          uc_tl(npx+1, j) = c3*utmp_tl(npx, j) + c2*utmp_tl(npx+1, j) + &
&           c1*utmp_tl(npx+2, j)
          uc(npx+1, j) = c3*utmp(npx, j) + c2*utmp(npx+1, j) + c1*utmp(&
&           npx+2, j)
          ut_tl(npx-1, j) = rsin_u(npx-1, j)*(uc_tl(npx-1, j)-cosa_u(npx&
&           -1, j)*v_tl(npx-1, j))
          ut(npx-1, j) = (uc(npx-1, j)-v(npx-1, j)*cosa_u(npx-1, j))*&
&           rsin_u(npx-1, j)
          ut_tl(npx, j) = rsin_u(npx, j)*uc_tl(npx, j)
          ut(npx, j) = uc(npx, j)*rsin_u(npx, j)
          ut_tl(npx+1, j) = rsin_u(npx+1, j)*(uc_tl(npx+1, j)-cosa_u(npx&
&           +1, j)*v_tl(npx+1, j))
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
        vtmp_tl(0, j) = -utmp_tl(1-j, 0)
        vtmp(0, j) = -utmp(1-j, 0)
      END DO
    END IF
    IF (nw_corner) THEN
      DO j=0,2
        vtmp_tl(0, npy+j) = utmp_tl(j+1, npy)
        vtmp(0, npy+j) = utmp(j+1, npy)
      END DO
    END IF
    IF (se_corner) THEN
      DO j=-2,0
        vtmp_tl(npx, j) = utmp_tl(ie+j, 0)
        vtmp(npx, j) = utmp(ie+j, 0)
      END DO
    END IF
    IF (ne_corner) THEN
      DO j=0,2
        vtmp_tl(npx, npy+j) = -utmp_tl(ie-j, npy)
        vtmp(npx, npy+j) = -utmp(ie-j, npy)
      END DO
    END IF
    IF (grid_type .LT. 3) THEN
      DO j=js-1,je+2
        IF (j .EQ. 1) THEN
          DO i=is-1,ie+1
! 3-pt extrapolation -----------------------------------------
            vc_tl(i, 1) = rsin_v(i, 1)*(t14*(vtmp_tl(i, 0)+vtmp_tl(i, 1)&
&             )+t12*(vtmp_tl(i, -1)+vtmp_tl(i, 2))+t15*(vtmp_tl(i, -2)+&
&             vtmp_tl(i, 3)))
            vc(i, 1) = (t14*(vtmp(i, 0)+vtmp(i, 1))+t12*(vtmp(i, -1)+&
&             vtmp(i, 2))+t15*(vtmp(i, -2)+vtmp(i, 3)))*rsin_v(i, 1)
            vt_tl(i, 1) = rsin_v(i, 1)*vc_tl(i, 1)
            vt(i, 1) = vc(i, 1)*rsin_v(i, 1)
          END DO
        ELSE IF (j .EQ. 0 .OR. j .EQ. npy - 1) THEN
          DO i=is-1,ie+1
            vc_tl(i, j) = c1*vtmp_tl(i, j-2) + c2*vtmp_tl(i, j-1) + c3*&
&             vtmp_tl(i, j)
            vc(i, j) = c1*vtmp(i, j-2) + c2*vtmp(i, j-1) + c3*vtmp(i, j)
            vt_tl(i, j) = rsin_v(i, j)*(vc_tl(i, j)-cosa_v(i, j)*u_tl(i&
&             , j))
            vt(i, j) = (vc(i, j)-u(i, j)*cosa_v(i, j))*rsin_v(i, j)
          END DO
        ELSE IF (j .EQ. 2 .OR. j .EQ. npy + 1) THEN
          DO i=is-1,ie+1
            vc_tl(i, j) = c1*vtmp_tl(i, j+1) + c2*vtmp_tl(i, j) + c3*&
&             vtmp_tl(i, j-1)
            vc(i, j) = c1*vtmp(i, j+1) + c2*vtmp(i, j) + c3*vtmp(i, j-1)
            vt_tl(i, j) = rsin_v(i, j)*(vc_tl(i, j)-cosa_v(i, j)*u_tl(i&
&             , j))
            vt(i, j) = (vc(i, j)-u(i, j)*cosa_v(i, j))*rsin_v(i, j)
          END DO
        ELSE IF (j .EQ. npy) THEN
          DO i=is-1,ie+1
! 3-pt extrapolation --------------------------------------------------------
            vc_tl(i, npy) = rsin_v(i, npy)*(t14*(vtmp_tl(i, npy-1)+&
&             vtmp_tl(i, npy))+t12*(vtmp_tl(i, npy-2)+vtmp_tl(i, npy+1))&
&             +t15*(vtmp_tl(i, npy-3)+vtmp_tl(i, npy+2)))
            vc(i, npy) = (t14*(vtmp(i, npy-1)+vtmp(i, npy))+t12*(vtmp(i&
&             , npy-2)+vtmp(i, npy+1))+t15*(vtmp(i, npy-3)+vtmp(i, npy+2&
&             )))*rsin_v(i, npy)
            vt_tl(i, npy) = rsin_v(i, npy)*vc_tl(i, npy)
            vt(i, npy) = vc(i, npy)*rsin_v(i, npy)
          END DO
        ELSE
! 4th order interpolation for interior points:
          DO i=is-1,ie+1
            vc_tl(i, j) = a2*(vtmp_tl(i, j-2)+vtmp_tl(i, j+1)) + a1*(&
&             vtmp_tl(i, j-1)+vtmp_tl(i, j))
            vc(i, j) = a2*(vtmp(i, j-2)+vtmp(i, j+1)) + a1*(vtmp(i, j-1)&
&             +vtmp(i, j))
            vt_tl(i, j) = rsin_v(i, j)*(vc_tl(i, j)-cosa_v(i, j)*u_tl(i&
&             , j))
            vt(i, j) = (vc(i, j)-u(i, j)*cosa_v(i, j))*rsin_v(i, j)
          END DO
        END IF
      END DO
    ELSE
! 4th order interpolation:
      DO j=js-1,je+2
        DO i=is-1,ie+1
          vc_tl(i, j) = a2*(vtmp_tl(i, j-2)+vtmp_tl(i, j+1)) + a1*(&
&           vtmp_tl(i, j-1)+vtmp_tl(i, j))
          vc(i, j) = a2*(vtmp(i, j-2)+vtmp(i, j+1)) + a1*(vtmp(i, j-1)+&
&           vtmp(i, j))
          vt_tl(i, j) = vc_tl(i, j)
          vt(i, j) = vc(i, j)
        END DO
      END DO
    END IF
  END SUBROUTINE D2A2C_VECT_TLM
!  Differentiation of fill2_4corners in forward (tangent) mode (with options r8):
!   variations   of useful results: q1 q2
!   with respect to varying inputs: q1 q2
  SUBROUTINE FILL2_4CORNERS_TLM(q1, q1_tl, q2, q2_tl, dir)
    IMPLICIT NONE
! This routine fill the 4 corners of the scalar fileds only as needed by c_core
! 1: x-dir; 2: y-dir
    INTEGER, INTENT(IN) :: dir
    REAL, INTENT(INOUT) :: q1(isd:ied, jsd:jed)
    REAL, INTENT(INOUT) :: q1_tl(isd:ied, jsd:jed)
    REAL, INTENT(INOUT) :: q2(isd:ied, jsd:jed)
    REAL, INTENT(INOUT) :: q2_tl(isd:ied, jsd:jed)
    SELECT CASE  (dir) 
    CASE (1) 
      IF (sw_corner) THEN
        q1_tl(-1, 0) = q1_tl(0, 2)
        q1(-1, 0) = q1(0, 2)
        q1_tl(0, 0) = q1_tl(0, 1)
        q1(0, 0) = q1(0, 1)
        q2_tl(-1, 0) = q2_tl(0, 2)
        q2(-1, 0) = q2(0, 2)
        q2_tl(0, 0) = q2_tl(0, 1)
        q2(0, 0) = q2(0, 1)
      END IF
      IF (se_corner) THEN
        q1_tl(npx+1, 0) = q1_tl(npx, 2)
        q1(npx+1, 0) = q1(npx, 2)
        q1_tl(npx, 0) = q1_tl(npx, 1)
        q1(npx, 0) = q1(npx, 1)
        q2_tl(npx+1, 0) = q2_tl(npx, 2)
        q2(npx+1, 0) = q2(npx, 2)
        q2_tl(npx, 0) = q2_tl(npx, 1)
        q2(npx, 0) = q2(npx, 1)
      END IF
      IF (nw_corner) THEN
        q1_tl(0, npy) = q1_tl(0, npy-1)
        q1(0, npy) = q1(0, npy-1)
        q1_tl(-1, npy) = q1_tl(0, npy-2)
        q1(-1, npy) = q1(0, npy-2)
        q2_tl(0, npy) = q2_tl(0, npy-1)
        q2(0, npy) = q2(0, npy-1)
        q2_tl(-1, npy) = q2_tl(0, npy-2)
        q2(-1, npy) = q2(0, npy-2)
      END IF
      IF (ne_corner) THEN
        q1_tl(npx, npy) = q1_tl(npx, npy-1)
        q1(npx, npy) = q1(npx, npy-1)
        q1_tl(npx+1, npy) = q1_tl(npx, npy-2)
        q1(npx+1, npy) = q1(npx, npy-2)
        q2_tl(npx, npy) = q2_tl(npx, npy-1)
        q2(npx, npy) = q2(npx, npy-1)
        q2_tl(npx+1, npy) = q2_tl(npx, npy-2)
        q2(npx+1, npy) = q2(npx, npy-2)
      END IF
    CASE (2) 
      IF (sw_corner) THEN
        q1_tl(0, 0) = q1_tl(1, 0)
        q1(0, 0) = q1(1, 0)
        q1_tl(0, -1) = q1_tl(2, 0)
        q1(0, -1) = q1(2, 0)
        q2_tl(0, 0) = q2_tl(1, 0)
        q2(0, 0) = q2(1, 0)
        q2_tl(0, -1) = q2_tl(2, 0)
        q2(0, -1) = q2(2, 0)
      END IF
      IF (se_corner) THEN
        q1_tl(npx, 0) = q1_tl(npx-1, 0)
        q1(npx, 0) = q1(npx-1, 0)
        q1_tl(npx, -1) = q1_tl(npx-2, 0)
        q1(npx, -1) = q1(npx-2, 0)
        q2_tl(npx, 0) = q2_tl(npx-1, 0)
        q2(npx, 0) = q2(npx-1, 0)
        q2_tl(npx, -1) = q2_tl(npx-2, 0)
        q2(npx, -1) = q2(npx-2, 0)
      END IF
      IF (nw_corner) THEN
        q1_tl(0, npy) = q1_tl(1, npy)
        q1(0, npy) = q1(1, npy)
        q1_tl(0, npy+1) = q1_tl(2, npy)
        q1(0, npy+1) = q1(2, npy)
        q2_tl(0, npy) = q2_tl(1, npy)
        q2(0, npy) = q2(1, npy)
        q2_tl(0, npy+1) = q2_tl(2, npy)
        q2(0, npy+1) = q2(2, npy)
      END IF
      IF (ne_corner) THEN
        q1_tl(npx, npy) = q1_tl(npx-1, npy)
        q1(npx, npy) = q1(npx-1, npy)
        q1_tl(npx, npy+1) = q1_tl(npx-2, npy)
        q1(npx, npy+1) = q1(npx-2, npy)
        q2_tl(npx, npy) = q2_tl(npx-1, npy)
        q2(npx, npy) = q2(npx-1, npy)
        q2_tl(npx, npy+1) = q2_tl(npx-2, npy)
        q2(npx, npy+1) = q2(npx-2, npy)
      END IF
    END SELECT
  END SUBROUTINE FILL2_4CORNERS_TLM
!  Differentiation of fill_4corners in forward (tangent) mode (with options r8):
!   variations   of useful results: q
!   with respect to varying inputs: q
  SUBROUTINE FILL_4CORNERS_TLM(q, q_tl, dir)
    IMPLICIT NONE
! This routine fill the 4 corners of the scalar fileds only as needed by c_core
! 1: x-dir; 2: y-dir
    INTEGER, INTENT(IN) :: dir
    REAL, INTENT(INOUT) :: q(isd:ied, jsd:jed)
    REAL, INTENT(INOUT) :: q_tl(isd:ied, jsd:jed)
    SELECT CASE  (dir) 
    CASE (1) 
      IF (sw_corner) THEN
        q_tl(-1, 0) = q_tl(0, 2)
        q(-1, 0) = q(0, 2)
        q_tl(0, 0) = q_tl(0, 1)
        q(0, 0) = q(0, 1)
      END IF
      IF (se_corner) THEN
        q_tl(npx+1, 0) = q_tl(npx, 2)
        q(npx+1, 0) = q(npx, 2)
        q_tl(npx, 0) = q_tl(npx, 1)
        q(npx, 0) = q(npx, 1)
      END IF
      IF (nw_corner) THEN
        q_tl(0, npy) = q_tl(0, npy-1)
        q(0, npy) = q(0, npy-1)
        q_tl(-1, npy) = q_tl(0, npy-2)
        q(-1, npy) = q(0, npy-2)
      END IF
      IF (ne_corner) THEN
        q_tl(npx, npy) = q_tl(npx, npy-1)
        q(npx, npy) = q(npx, npy-1)
        q_tl(npx+1, npy) = q_tl(npx, npy-2)
        q(npx+1, npy) = q(npx, npy-2)
      END IF
    CASE (2) 
      IF (sw_corner) THEN
        q_tl(0, 0) = q_tl(1, 0)
        q(0, 0) = q(1, 0)
        q_tl(0, -1) = q_tl(2, 0)
        q(0, -1) = q(2, 0)
      END IF
      IF (se_corner) THEN
        q_tl(npx, 0) = q_tl(npx-1, 0)
        q(npx, 0) = q(npx-1, 0)
        q_tl(npx, -1) = q_tl(npx-2, 0)
        q(npx, -1) = q(npx-2, 0)
      END IF
      IF (nw_corner) THEN
        q_tl(0, npy) = q_tl(1, npy)
        q(0, npy) = q(1, npy)
        q_tl(0, npy+1) = q_tl(2, npy)
        q(0, npy+1) = q(2, npy)
      END IF
      IF (ne_corner) THEN
        q_tl(npx, npy) = q_tl(npx-1, npy)
        q(npx, npy) = q(npx-1, npy)
        q_tl(npx, npy+1) = q_tl(npx-2, npy)
        q(npx, npy+1) = q(npx-2, npy)
      END IF
    END SELECT
  END SUBROUTINE FILL_4CORNERS_TLM
  SUBROUTINE D2A2C_TLM(u, u_tl, v, v_tl, um, um_tl, vm, vm_tl, ua, ua_tl&
&   , va, va_tl, uc, uc_tl, vc, vc_tl, dord4)
    IMPLICIT NONE
    REAL, DIMENSION(isd:ied, jsd:jed + 1), INTENT(IN) :: u, um
    REAL, DIMENSION(isd:ied, jsd:jed+1), INTENT(IN) :: u_tl, um_tl
    REAL, DIMENSION(isd:ied + 1, jsd:jed), INTENT(IN) :: v, vm
    REAL, DIMENSION(isd:ied+1, jsd:jed), INTENT(IN) :: v_tl, vm_tl
    LOGICAL, INTENT(IN) :: dord4
    REAL, DIMENSION(isd:ied + 1, jsd:jed), INTENT(OUT) :: uc
    REAL, DIMENSION(isd:ied+1, jsd:jed), INTENT(OUT) :: uc_tl
    REAL, DIMENSION(isd:ied, jsd:jed + 1), INTENT(OUT) :: vc
    REAL, DIMENSION(isd:ied, jsd:jed+1), INTENT(OUT) :: vc_tl
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(OUT) :: ua, va
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(OUT) :: ua_tl, va_tl
! Local 
    REAL, DIMENSION(isd:ied, jsd:jed) :: utmp, vtmp
    REAL, DIMENSION(isd:ied, jsd:jed) :: utmp_tl, vtmp_tl
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
      utmp_tl = 0.0_8
    ELSE
      min1 = npy - npt
      utmp_tl = 0.0_8
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
        utmp_tl(i, j) = a2*(u_tl(i, j-1)+u_tl(i, j+2)) + a1*(u_tl(i, j)+&
&         u_tl(i, j+1))
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
      vtmp_tl = 0.0_8
    ELSE
      min3 = npy - npt
      vtmp_tl = 0.0_8
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
        vtmp_tl(i, j) = a2*(v_tl(i-1, j)+v_tl(i+2, j)) + a1*(v_tl(i, j)+&
&         v_tl(i+1, j))
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
            utmp_tl(i, j) = 0.5*(u_tl(i, j)+u_tl(i, j+1))
            utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
            vtmp_tl(i, j) = 0.5*(v_tl(i, j)+v_tl(i+1, j))
            vtmp(i, j) = 0.5*(v(i, j)+v(i+1, j))
          END DO
        END DO
      END IF
      IF (je + 1 .EQ. npy .OR. jed .GE. npy - npt) THEN
        DO j=npy-npt+1,jed
          DO i=isd,ied
            utmp_tl(i, j) = 0.5*(u_tl(i, j)+u_tl(i, j+1))
            utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
            vtmp_tl(i, j) = 0.5*(v_tl(i, j)+v_tl(i+1, j))
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
            utmp_tl(i, j) = 0.5*(u_tl(i, j)+u_tl(i, j+1))
            utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
            vtmp_tl(i, j) = 0.5*(v_tl(i, j)+v_tl(i+1, j))
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
            utmp_tl(i, j) = 0.5*(u_tl(i, j)+u_tl(i, j+1))
            utmp(i, j) = 0.5*(u(i, j)+u(i, j+1))
            vtmp_tl(i, j) = 0.5*(v_tl(i, j)+v_tl(i+1, j))
            vtmp(i, j) = 0.5*(v(i, j)+v(i+1, j))
          END DO
        END DO
      END IF
    END IF
    DO j=js-1-id,je+1+id
      DO i=is-1-id,ie+1+id
        ua_tl(i, j) = rsin2(i, j)*(utmp_tl(i, j)-cosa_s(i, j)*vtmp_tl(i&
&         , j))
        ua(i, j) = (utmp(i, j)-vtmp(i, j)*cosa_s(i, j))*rsin2(i, j)
        va_tl(i, j) = rsin2(i, j)*(vtmp_tl(i, j)-cosa_s(i, j)*utmp_tl(i&
&         , j))
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
        utmp_tl(i, j) = a2*(um_tl(i, j-1)+um_tl(i, j+2)) + a1*(um_tl(i, &
&         j)+um_tl(i, j+1))
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
        vtmp_tl(i, j) = a2*(vm_tl(i-1, j)+vm_tl(i+2, j)) + a1*(vm_tl(i, &
&         j)+vm_tl(i+1, j))
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
            utmp_tl(i, j) = 0.5*(um_tl(i, j)+um_tl(i, j+1))
            utmp(i, j) = 0.5*(um(i, j)+um(i, j+1))
            vtmp_tl(i, j) = 0.5*(vm_tl(i, j)+vm_tl(i+1, j))
            vtmp(i, j) = 0.5*(vm(i, j)+vm(i+1, j))
          END DO
        END DO
      END IF
      IF (je + 1 .EQ. npy .OR. jed .GE. npy - npt) THEN
        DO j=npy-npt+1,jed
          DO i=isd,ied
            utmp_tl(i, j) = 0.5*(um_tl(i, j)+um_tl(i, j+1))
            utmp(i, j) = 0.5*(um(i, j)+um(i, j+1))
            vtmp_tl(i, j) = 0.5*(vm_tl(i, j)+vm_tl(i+1, j))
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
            utmp_tl(i, j) = 0.5*(um_tl(i, j)+um_tl(i, j+1))
            utmp(i, j) = 0.5*(um(i, j)+um(i, j+1))
            vtmp_tl(i, j) = 0.5*(vm_tl(i, j)+vm_tl(i+1, j))
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
            utmp_tl(i, j) = 0.5*(um_tl(i, j)+um_tl(i, j+1))
            utmp(i, j) = 0.5*(um(i, j)+um(i, j+1))
            vtmp_tl(i, j) = 0.5*(vm_tl(i, j)+vm_tl(i+1, j))
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
        utmp_tl(i, 0) = -vtmp_tl(0, 1-i)
        utmp(i, 0) = -vtmp(0, 1-i)
      END DO
    END IF
    IF (se_corner) THEN
      DO i=0,2
        utmp_tl(npx+i, 0) = vtmp_tl(npx, i+1)
        utmp(npx+i, 0) = vtmp(npx, i+1)
      END DO
    END IF
    IF (ne_corner) THEN
      DO i=0,2
        utmp_tl(npx+i, npy) = -vtmp_tl(npx, je-i)
        utmp(npx+i, npy) = -vtmp(npx, je-i)
      END DO
    END IF
    IF (nw_corner) THEN
      DO i=-2,0
        utmp_tl(i, npy) = vtmp_tl(0, je+i)
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
        uc_tl(i, j) = a1*(utmp_tl(i-1, j)+utmp_tl(i, j)) + a2*(utmp_tl(i&
&         -2, j)+utmp_tl(i+1, j))
        uc(i, j) = a1*(utmp(i-1, j)+utmp(i, j)) + a2*(utmp(i-2, j)+utmp(&
&         i+1, j))
      END DO
    END DO
    IF (grid_type .LT. 3) THEN
      IF (is .EQ. 1) THEN
        DO j=js,je
! 3-pt extrapolation --------------------------------------------------
          uc_tl(1, j) = rsin_u(1, j)*(t14*(utmp_tl(0, j)+utmp_tl(1, j))+&
&           t12*(utmp_tl(-1, j)+utmp_tl(2, j))+t15*(utmp_tl(-2, j)+&
&           utmp_tl(3, j)))
          uc(1, j) = (t14*(utmp(0, j)+utmp(1, j))+t12*(utmp(-1, j)+utmp(&
&           2, j))+t15*(utmp(-2, j)+utmp(3, j)))*rsin_u(1, j)
          uc_tl(2, j) = c1*utmp_tl(3, j) + c2*utmp_tl(2, j) + c3*utmp_tl&
&           (1, j)
          uc(2, j) = c1*utmp(3, j) + c2*utmp(2, j) + c3*utmp(1, j)
        END DO
      END IF
      IF (ie + 1 .EQ. npx) THEN
        DO j=js,je
          uc_tl(npx-1, j) = c1*utmp_tl(npx-3, j) + c2*utmp_tl(npx-2, j) &
&           + c3*utmp_tl(npx-1, j)
          uc(npx-1, j) = c1*utmp(npx-3, j) + c2*utmp(npx-2, j) + c3*utmp&
&           (npx-1, j)
! 3-pt extrapolation --------------------------------------------------------
          uc_tl(npx, j) = rsin_u(npx, j)*(t14*(utmp_tl(npx-1, j)+utmp_tl&
&           (npx, j))+t12*(utmp_tl(npx-2, j)+utmp_tl(npx+1, j))+t15*(&
&           utmp_tl(npx-3, j)+utmp_tl(npx+2, j)))
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
        vtmp_tl(0, j) = -utmp_tl(1-j, 0)
        vtmp(0, j) = -utmp(1-j, 0)
      END DO
    END IF
    IF (nw_corner) THEN
      DO j=0,2
        vtmp_tl(0, npy+j) = utmp_tl(j+1, npy)
        vtmp(0, npy+j) = utmp(j+1, npy)
      END DO
    END IF
    IF (se_corner) THEN
      DO j=-2,0
        vtmp_tl(npx, j) = utmp_tl(ie+j, 0)
        vtmp(npx, j) = utmp(ie+j, 0)
      END DO
    END IF
    IF (ne_corner) THEN
      DO j=0,2
        vtmp_tl(npx, npy+j) = -utmp_tl(ie-j, npy)
        vtmp(npx, npy+j) = -utmp(ie-j, npy)
      END DO
    END IF
    IF (grid_type .LT. 3) THEN
      DO j=js,je+1
        IF (j .EQ. 1) THEN
          DO i=is,ie
! 3-pt extrapolation -----------------------------------------
            vc_tl(i, 1) = rsin_v(i, 1)*(t14*(vtmp_tl(i, 0)+vtmp_tl(i, 1)&
&             )+t12*(vtmp_tl(i, -1)+vtmp_tl(i, 2))+t15*(vtmp_tl(i, -2)+&
&             vtmp_tl(i, 3)))
            vc(i, 1) = (t14*(vtmp(i, 0)+vtmp(i, 1))+t12*(vtmp(i, -1)+&
&             vtmp(i, 2))+t15*(vtmp(i, -2)+vtmp(i, 3)))*rsin_v(i, 1)
          END DO
        ELSE IF (j .EQ. npy - 1) THEN
          DO i=is,ie
            vc_tl(i, j) = c1*vtmp_tl(i, j-2) + c2*vtmp_tl(i, j-1) + c3*&
&             vtmp_tl(i, j)
            vc(i, j) = c1*vtmp(i, j-2) + c2*vtmp(i, j-1) + c3*vtmp(i, j)
          END DO
        ELSE IF (j .EQ. 2) THEN
          DO i=is,ie
!          vc(i,j) = c1*vtmp(i,j+1) + c2*vtmp(i,j) + c3*vtmp(i,j-1)
            vc_tl(i, j) = c3*vtmp_tl(i, j-1) + c2*vtmp_tl(i, j) + c1*&
&             vtmp_tl(i, j+1)
            vc(i, j) = c3*vtmp(i, j-1) + c2*vtmp(i, j) + c1*vtmp(i, j+1)
          END DO
        ELSE IF (j .EQ. npy) THEN
          DO i=is,ie
! 3-pt extrapolation --------------------------------------------------------
            vc_tl(i, npy) = rsin_v(i, npy)*(t14*(vtmp_tl(i, npy-1)+&
&             vtmp_tl(i, npy))+t12*(vtmp_tl(i, npy-2)+vtmp_tl(i, npy+1))&
&             +t15*(vtmp_tl(i, npy-3)+vtmp_tl(i, npy+2)))
            vc(i, npy) = (t14*(vtmp(i, npy-1)+vtmp(i, npy))+t12*(vtmp(i&
&             , npy-2)+vtmp(i, npy+1))+t15*(vtmp(i, npy-3)+vtmp(i, npy+2&
&             )))*rsin_v(i, npy)
          END DO
        ELSE
! 4th order interpolation for interior points:
          DO i=is,ie
            vc_tl(i, j) = a2*(vtmp_tl(i, j-2)+vtmp_tl(i, j+1)) + a1*(&
&             vtmp_tl(i, j-1)+vtmp_tl(i, j))
            vc(i, j) = a2*(vtmp(i, j-2)+vtmp(i, j+1)) + a1*(vtmp(i, j-1)&
&             +vtmp(i, j))
          END DO
        END IF
      END DO
    ELSE
! 4th order interpolation:
      DO j=js,je+1
        DO i=is,ie
          vc_tl(i, j) = a2*(vtmp_tl(i, j-2)+vtmp_tl(i, j+1)) + a1*(&
&           vtmp_tl(i, j-1)+vtmp_tl(i, j))
          vc(i, j) = a2*(vtmp(i, j-2)+vtmp(i, j+1)) + a1*(vtmp(i, j-1)+&
&           vtmp(i, j))
        END DO
      END DO
    END IF
  END SUBROUTINE D2A2C_TLM

END MODULE SW_CORE_TLM_MOD
