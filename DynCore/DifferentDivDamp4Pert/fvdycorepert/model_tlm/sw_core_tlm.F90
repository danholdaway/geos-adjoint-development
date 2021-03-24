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
 use tp_core_tlm_mod,   only: fv_tp_2d_tlm, pert_ppm_tlm, copy_corners_tlm
 use tp_core_mod,       only: fv_tp_2d
 use fv_grid_utils_mod, only: sw_corner, se_corner, ne_corner, nw_corner,       &
                              cosa_u, cosa_v, cosa_s, sina_s, sina_u, sina_v,   &
                              rsin_u, rsin_v, rsin_v, rsina, ec1, ec2, ew, es,  &
                              big_number, da_min_c, da_min, fC, f0,   &
                              rsin2, divg_u, divg_v, Gnomonic_grid
 use test_cases_mod,    only: test_case
 use sw_core_mod,       only: xtp_u, ytp_v


 implicit none

  real, parameter:: r3 =   1./3.
  real, parameter:: t11=27./28., t12=-13./28., t13=3./7., t14=6./7., t15=3./28.
  real, parameter:: s11=11./14., s13=-13./14., s14=4./7., s15=3./14.
  real, parameter:: p1 =  7./12.     ! 0.58333333
  real, parameter:: p2 = -1./12.
  real(i_precision), parameter:: a1 =  0.5625
  real(i_precision), parameter:: a2 = -0.0625
  real, parameter:: c1 = -2./14.
  real, parameter:: c2 = 11./14.
  real, parameter:: c3 =  5./14.
  real, parameter:: b1 =   1./30.
  real, parameter:: b2 = -13./60.
  real, parameter:: b3 = -13./60.
  real, parameter:: b4 =  0.45
  real, parameter:: b5 = -0.05

  private
  public :: c_sw_tlm, d_sw_tlm, d2a2c_tlm, d2a2c_vect_tlm, divergence_corner_tlm, &
            divergence_corner_onlytlm

  contains

   SUBROUTINE C_SW_TLM(delpc, delpc_tl, delp, delp_tl, ptc, ptc_tl, pt, &
&   pt_tl, u, u_tl, v, v_tl, w, w_tl, uc, uc_tl, vc, vc_tl, ua, ua_tl, &
&   va, va_tl, wc, wc_tl, ut, ut_tl, vt, vt_tl, dt2, hydrostatic, dord4, dord4_tj)
    IMPLICIT NONE
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
&   ut_tl, vt_tl, wc_tl
    REAL, INTENT(IN) :: dt2
    LOGICAL, INTENT(IN) :: hydrostatic
    LOGICAL, INTENT(IN) :: dord4, dord4_tj
! Local:
    REAL, DIMENSION(is - 1:ie + 1, js - 1:je + 1) :: vort, ke
    REAL, DIMENSION(is-1:ie+1, js-1:je+1) :: vort_tl, ke_tl
    REAL, DIMENSION(is - 1:ie + 2, js - 1:je + 1) :: fx, fx1, fx2
    REAL, DIMENSION(is-1:ie+2, js-1:je+1) :: fx_tl, fx1_tl, fx2_tl
    REAL, DIMENSION(is - 1:ie + 1, js - 1:je + 2) :: fy, fy1, fy2
    REAL, DIMENSION(is-1:ie+1, js-1:je+2) :: fy_tl, fy1_tl, fy2_tl
    REAL :: dt4
    INTEGER :: i, j, is2, ie1
    INTEGER :: iep1, jep1
    INTRINSIC MAX
    INTRINSIC MIN
    iep1 = ie + 1
    jep1 = je + 1
    CALL D2A2C_VECT_TLM(u, u_tl, v, v_tl, ua, ua_tl, va, va_tl, uc, &
&                 uc_tl, vc, vc_tl, ut, ut_tl, vt, vt_tl, dord4, dord4_tj)
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
        fx1_tl = 0.0
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
        fx_tl = 0.0
        fx2_tl = 0.0
      ELSE
        fx_tl = 0.0
        fx1_tl = 0.0
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
        fx2_tl = 0.0
      END IF
    ELSE
      IF (grid_type .LT. 3) THEN
        CALL FILL_4CORNERS_TLM(w, w_tl, 1)
        fx_tl = 0.0
        fx1_tl = 0.0
        fx2_tl = 0.0
      ELSE
        fx_tl = 0.0
        fx1_tl = 0.0
        fx2_tl = 0.0
      END IF
      DO j=js-1,je+1
        DO i=is-1,ie+2
          IF (ut(i, j) .GT. 0.) THEN
            fx1_tl(i, j) = delp_tl(i-1, j)
            fx1(i, j) = delp(i-1, j)
            fx_tl(i, j) = pt_tl(i-1, j)
            fx(i, j) = pt(i-1, j)
            fx2_tl(i, j) = w_tl(i-1, j)
            fx2(i, j) = w(i-1, j)
          ELSE
            fx1_tl(i, j) = delp_tl(i, j)
            fx1(i, j) = delp(i, j)
            fx_tl(i, j) = pt_tl(i, j)
            fx(i, j) = pt(i, j)
            fx2_tl(i, j) = w_tl(i, j)
            fx2(i, j) = w(i, j)
          END IF
          fx1_tl(i, j) = ut_tl(i, j)*fx1(i, j) + ut(i, j)*fx1_tl(i, j)
          fx1(i, j) = ut(i, j)*fx1(i, j)
          fx_tl(i, j) = fx1_tl(i, j)*fx(i, j) + fx1(i, j)*fx_tl(i, j)
          fx(i, j) = fx1(i, j)*fx(i, j)
          fx2_tl(i, j) = fx1_tl(i, j)*fx2(i, j) + fx1(i, j)*fx2_tl(i, j)
          fx2(i, j) = fx1(i, j)*fx2(i, j)
        END DO
      END DO
    END IF
! Ydir:
    IF (grid_type .LT. 3) CALL FILL2_4CORNERS_TLM(delp, delp_tl, pt, &
&                                           pt_tl, 2)
    IF (hydrostatic) THEN
      fy1_tl = 0.0
      fy_tl = 0.0
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
        ke_tl = 0.0
      ELSE
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
        ke_tl = 0.0
      END IF
    ELSE
      IF (grid_type .LT. 3) THEN
        CALL FILL_4CORNERS_TLM(w, w_tl, 2)
        fy1_tl = 0.0
        fy2_tl = 0.0
        fy_tl = 0.0
      ELSE
        fy1_tl = 0.0
        fy2_tl = 0.0
        fy_tl = 0.0
      END IF
      DO j=js-1,je+2
        DO i=is-1,ie+1
          IF (vt(i, j) .GT. 0.) THEN
            fy1_tl(i, j) = delp_tl(i, j-1)
            fy1(i, j) = delp(i, j-1)
            fy_tl(i, j) = pt_tl(i, j-1)
            fy(i, j) = pt(i, j-1)
            fy2_tl(i, j) = w_tl(i, j-1)
            fy2(i, j) = w(i, j-1)
          ELSE
            fy1_tl(i, j) = delp_tl(i, j)
            fy1(i, j) = delp(i, j)
            fy_tl(i, j) = pt_tl(i, j)
            fy(i, j) = pt(i, j)
            fy2_tl(i, j) = w_tl(i, j)
            fy2(i, j) = w(i, j)
          END IF
          fy1_tl(i, j) = vt_tl(i, j)*fy1(i, j) + vt(i, j)*fy1_tl(i, j)
          fy1(i, j) = vt(i, j)*fy1(i, j)
          fy_tl(i, j) = fy1_tl(i, j)*fy(i, j) + fy1(i, j)*fy_tl(i, j)
          fy(i, j) = fy1(i, j)*fy(i, j)
          fy2_tl(i, j) = fy1_tl(i, j)*fy2(i, j) + fy1(i, j)*fy2_tl(i, j)
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
          wc_tl(i, j) = ((w_tl(i, j)*delp(i, j)+w(i, j)*delp_tl(i, j)+&
&           rarea(i, j)*(fx2_tl(i, j)-fx2_tl(i+1, j)+fy2_tl(i, j)-fy2_tl&
&           (i, j+1)))*delpc(i, j)-(w(i, j)*delp(i, j)+(fx2(i, j)-fx2(i+&
&           1, j)+fy2(i, j)-fy2(i, j+1))*rarea(i, j))*delpc_tl(i, j))/&
&           delpc(i, j)**2
          wc(i, j) = (w(i, j)*delp(i, j)+(fx2(i, j)-fx2(i+1, j)+fy2(i, j&
&           )-fy2(i, j+1))*rarea(i, j))/delpc(i, j)
        END DO
      END DO
      ke_tl = 0.0
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
    vort_tl = 0.0
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

  SUBROUTINE D_SW_TLM(delpc, delpc_tj, delpc_tl, delp, delp_tl, ptc, ptc_tj, ptc_tl, pt, &
&   pt_tl, u, u_tl, v, v_tl, w, w_tl, uc, uc_tl, vc, vc_tl, ua, ua_tl, &
&   va, va_tl, divg_d, divg_d_tl, xflux, xflux_tl, yflux, yflux_tl, cx, &
&   cx_tl, cy, cy_tl, crx_adv, crx_adv_tl, cry_adv, cry_adv_tl, xfx_adv&
&   , xfx_adv_tl, yfx_adv, yfx_adv_tl, zvir, sphum, nq, q, q_tl, k, km, &
&   inline_q, pkz, pkz_tl, dt, hord_tr, hord_mt, hord_vt, hord_tm, &
&   hord_dp, hord_tr_tj, hord_mt_tj, hord_vt_tj, hord_tm_tj, &
&   hord_dp_tj, nord, nord_tj, dddmp, dddmp_tj, d2_bg, d2_bg_tj, d4_bg, d4_bg_tj, &
&   vtdm4, d_con, hydrostatic, ppm_limiter,proc_id,it)
    IMPLICIT NONE
! damping applied to relative vorticity:
!   if ( vtdm4>0. ) then
!      damp4 = (vtdm4*da_min)**(nord+1)
!      call del6_flux(nord, npx, npy, damp4, wk, u, v)
!   endif
    INTEGER, INTENT(IN) :: hord_tr, hord_mt, hord_vt, hord_tm, hord_dp, proc_id,it
    INTEGER, INTENT(IN) :: hord_tr_tj, hord_mt_tj, hord_vt_tj, hord_tm_tj, hord_dp_tj
! nord=1 (del-4) or 3 (del-8)
    INTEGER, INTENT(IN) :: nord, nord_tj
    INTEGER, INTENT(IN) :: sphum, nq, k, km
    REAL, INTENT(IN) :: dt, dddmp, dddmp_tj, d2_bg, d2_bg_tj, d4_bg, d4_bg_tj, vtdm4, d_con
    REAL, INTENT(IN) :: zvir
    REAL, INTENT(IN) :: ppm_limiter
    REAL(p_precision), INTENT(IN) :: pkz(is:ie, js:je)
    REAL(p_precision), INTENT(IN) :: pkz_tl(is:ie, js:je)
! divergence
    REAL, INTENT(INOUT) :: divg_d(isd:ied+1, jsd:jed+1)
    REAL :: divg_d_tj(isd:ied+1, jsd:jed+1)
    REAL, INTENT(INOUT) :: divg_d_tl(isd:ied+1, jsd:jed+1)
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(INOUT) :: delp, pt, ua, va&
&   , w
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(INOUT) :: delp_tl, pt_tl, &
&   ua_tl, va_tl, w_tl
    REAL, DIMENSION(isd:ied, jsd:jed + 1), INTENT(INOUT) :: u, vc
    REAL, DIMENSION(isd:ied, jsd:jed + 1) :: vc_tj
    REAL, DIMENSION(isd:ied, jsd:jed+1), INTENT(INOUT) :: u_tl, vc_tl
    REAL, DIMENSION(isd:ied + 1, jsd:jed), INTENT(INOUT) :: v, uc
    REAL, DIMENSION(isd:ied + 1, jsd:jed) :: uc_tj
    REAL, DIMENSION(isd:ied+1, jsd:jed), INTENT(INOUT) :: v_tl, uc_tl
    REAL, INTENT(INOUT) :: q(isd:ied, jsd:jed, km, nq)
    REAL, INTENT(INOUT) :: q_tl(isd:ied, jsd:jed, km, nq)
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(OUT) :: delpc, ptc
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(OUT) :: delpc_tl, ptc_tl
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(OUT) :: delpc_tj, ptc_tj
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
    REAL(i_precision) :: damp_tj
    REAL(i_precision) :: damp_tl
    REAL(i_precision), SAVE :: r1b4=1.0/4.0
    REAL(i_precision), SAVE :: r1b16=1.0/16.0
    REAL(ke_precision), DIMENSION(is:ie + 1, js:je + 1) :: ub, vb
    REAL(ke_precision), DIMENSION(is:ie + 1, js:je + 1) :: ub_tj, vb_tj
    REAL(ke_precision), DIMENSION(is:ie + 1, js:je + 1) :: ubtmp, vbtmp
    REAL(ke_precision), DIMENSION(is:ie+1, js:je+1) :: ub_tl, vb_tl
!  needs this for corner_comm
    REAL(ke_precision) :: ke(isd:ied+1, jsd:jed+1)
    REAL(ke_precision) :: ke_tj(isd:ied+1, jsd:jed+1)
    REAL(ke_precision) :: ke_tl(isd:ied+1, jsd:jed+1)
    REAL(ke_precision) :: dt4, dt5, dt6
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
! 1-D Y-direction Fluxes
    REAL :: fy(is:ie, js:je+1)
    REAL :: fy_tl(is:ie, js:je+1)
! 1-D temporary X and Y-direction Fluxes
    REAL :: fieldtmp(isd:ied,jsd:jed)
    REAL :: fxtmp(is:ie+1, js:je)
    REAL :: fytmp(is:ie+1, js:je)
! 1-D X-direction Fluxes
    REAL :: px(is:ie+1, js:je)
    REAL :: px_tl(is:ie+1, js:je)
! 1-D Y-direction Fluxes
    REAL :: py(is:ie, js:je+1)
    REAL :: py_tl(is:ie, js:je+1)
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
    REAL :: pwx1_tj, damp2_tj, abs0_tj, y3_tj, y2_tj, y1_tj, max5_tj, max6_tj, dd8_tj, du2_tj, dv2_tj

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
      ra_x_tl = 0.0
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
      ra_x_tl = 0.0
    END IF
    DO j=jsd,jed
      DO i=is,ie
        ra_x_tl(i, j) = xfx_adv_tl(i, j) - xfx_adv_tl(i+1, j)
        ra_x(i, j) = area(i, j) + xfx_adv(i, j) - xfx_adv(i+1, j)
      END DO
    END DO
    ra_y_tl = 0.0
    DO j=js,je
      DO i=isd,ied
        ra_y_tl(i, j) = yfx_adv_tl(i, j) - yfx_adv_tl(i, j+1)
        ra_y(i, j) = area(i, j) + yfx_adv(i, j) - yfx_adv(i, j+1)
      END DO
    END DO
    fy_tl = 0.0
    fx_tl = 0.0
    fxtmp = 0.0
    fytmp = 0.0
    fieldtmp = delp
    CALL FV_TP_2D_TLM(fieldtmp, delp_tl, crx_adv, crx_adv_tl, cry_adv, &
&               cry_adv_tl, npx, npy, hord_dp_tj, fxtmp, fx_tl, fytmp, fy_tl, &
&               xfx_adv, xfx_adv_tl, yfx_adv, yfx_adv_tl, area, ra_x, &
&               ra_x_tl, ra_y, ra_y_tl)
    CALL FV_TP_2D(delp, crx_adv, cry_adv, npx, npy, hord_dp, fx, fy,  &
&                 xfx_adv,  yfx_adv,  area, ra_x, ra_y) 
! shallow_water
    IF (shallow_water) THEN
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
      wk_tl = 0.0
    ELSE
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
    fxtmp = 0.0
    fytmp = 0.0
    fieldtmp = w
        CALL FV_TP_2D_TLM(fieldtmp, w_tl, crx_adv, crx_adv_tl, cry_adv, &
&                   cry_adv_tl, npx, npy, hord_vt_tj, fxtmp, px_tl, fytmp, py_tl&
&                   , xfx_adv, xfx_adv_tl, yfx_adv, yfx_adv_tl, area, &
&                   ra_x, ra_x_tl, ra_y, ra_y_tl, mfx=fx, mfx_tl=fx_tl, &
&                   mfy=fy, mfy_tl=fy_tl)
        CALL FV_TP_2D( w,crx_adv,cry_adv,npx,npy,hord_vt,px,py,xfx_adv,&
                       yfx_adv,area,ra_x,ra_y,mfx=fx,mfy=fy )
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
        px_tl = 0.0
        py_tl = 0.0
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
    fxtmp = 0.0
    fytmp = 0.0
    fieldtmp = pt
      CALL FV_TP_2D_TLM(fieldtmp, pt_tl, crx_adv, crx_adv_tl, cry_adv, &
&                 cry_adv_tl, npx, npy, hord_tm_tj, fxtmp, px_tl, fytmp, py_tl, &
&                 xfx_adv, xfx_adv_tl, yfx_adv, yfx_adv_tl, area, ra_x, &
&                 ra_x_tl, ra_y, ra_y_tl, mfx=fx, mfx_tl=fx_tl, mfy=fy, &
&                 mfy_tl=fy_tl)
      CALL FV_TP_2D(pt, crx_adv, cry_adv, npx, npy, hord_tm, px, py, &
&                 xfx_adv, yfx_adv, area, ra_x, ra_y, mfx=fx, mfy=fy)
      IF (inline_q) THEN
        wk_tl = 0.0
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
    fxtmp = 0.0
    fytmp = 0.0
    fieldtmp = q(isd:, jsd:, k, iq)
          CALL FV_TP_2D_TLM(fieldtmp, q_tl(isd:, jsd:, k, iq&
&                     ), crx_adv, crx_adv_tl, cry_adv, cry_adv_tl, npx, &
&                     npy, hord_tr_tj, fxtmp, px_tl, fytmp, py_tl, xfx_adv, &
&                     xfx_adv_tl, yfx_adv, yfx_adv_tl, area, ra_x, &
&                     ra_x_tl, ra_y, ra_y_tl, mfx=fx, mfx_tl=fx_tl, mfy=&
&                     fy, mfy_tl=fy_tl)
           CALL FV_TP_2D(q(isd:, jsd:, k, iq), crx_adv, cry_adv, npx, &
&                     npy, hord_tr, px, py, xfx_adv, &
&                     yfx_adv, area, ra_x, ra_y, mfx=fx, mfy=fy)
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
        wk_tl = 0.0
      END IF
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
      ubtmp = 0.0
      CALL YTP_V_TLM(vb, vb_tl, u, v, v_tl, ubtmp, ub_tl, hord_mt_tj)
      CALL YTP_V    (vb,        u, v,       ub,           hord_mt)
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
      vbtmp = 0.0
      CALL XTP_U_TLM(ub, ub_tl, u, u_tl, v, vbtmp, vb_tl, hord_mt_tj)
      CALL XTP_U    (ub,        u,       v,        vb,    hord_mt)
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
!call mp_corner_comm(ke, npx, npy) 
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
      
! PERTURBATIONS
! -------------
! Compute divergence damping for lienarized variables. nord_tj allows
! one to choose a different order of damping for the TLM/ADM.
! Because of nonlinear dependence and different routes through the code 
! there will be two trajectories for this part.
! We could then carry this trajectory forward into future TLM calculations.
! This is not done at this point but could be adapted for better performance.
! _tj is used to denote a new trajectory calculation.
! -------------

      damp_tj = dddmp_tj*da_min_c

      divg_d_tj = divg_d
      ptc_tj = ptc
      delpc_tj = delpc

      vort_tj = vort
      ke_tj = ke
      uc_tj = uc
      vc_tj = vc

      IF (nord_tj .EQ. 0) THEN
!       area ~ dxb*dyb*sin(alpha)
        DO j=js,je+1
          IF (j .EQ. 1 .OR. j .EQ. npy) THEN
            DO i=is-1,ie+1
              ptc_tl(i, j) = dyc(i, j)*sina_v(i, j)*u_tl(i, j)
              ptc_tj(i, j) = dyc(i, j)*sina_v(i, j)*u   (i, j)
            END DO
          ELSE
            DO i=is-1,ie+1
              ptc_tl(i, j) = dyc(i, j)*sina_v(i, j)*(u_tl(i, j)-0.5*&
&               cosa_v(i, j)*(va_tl(i, j-1)+va_tl(i, j)))
              ptc_tj(i, j) = dyc(i, j)*sina_v(i, j)*(u   (i, j)-0.5*&
&               cosa_v(i, j)*(va   (i, j-1)+va   (i, j)))
            END DO
          END IF
        END DO
        vort_tl = 0.0
        DO j=js-1,je+1
          DO i=is2,ie1
            vort_tl(i, j) = dxc(i, j)*sina_u(i, j)*(v_tl(i, j)-0.5*cosa_u(i, j)*(ua_tl(i-1, j)+ua_tl(i, j)))
            vort_tj(i, j) = dxc(i, j)*sina_u(i, j)*(v   (i, j)-0.5*cosa_u(i, j)*(ua   (i-1, j)+ua   (i, j)))
          END DO
          IF (is .EQ. 1) THEN
            vort_tl(1, j) = dxc(1, j)*sina_u(1, j)*v_tl(1, j)
            vort_tj(1, j) = dxc(1, j)*sina_u(1, j)*v   (1, j)
          END IF
          IF (ie + 1 .EQ. npx) THEN
            vort_tl(npx, j) = dxc(npx, j)*sina_u(npx, j)*v_tl(npx, j)
            vort_tj(npx, j) = dxc(npx, j)*sina_u(npx, j)*v   (npx, j)
          END IF
        END DO
        DO j=js,je+1
          DO i=is,ie+1
            delpc_tl(i, j) = vort_tl(i, j-1) - vort_tl(i, j) + ptc_tl(i-1, j) - ptc_tl(i, j)
            delpc_tj(i, j) = vort_tj(i, j-1) - vort_tj(i, j) + ptc_tj(i-1, j) - ptc_tj(i, j)
          END DO
        END DO
! Remove the extra term at the corners:
        IF (sw_corner) THEN
          delpc_tl(1, 1) = delpc_tl(1, 1) - vort_tl(1, 0)
          delpc_tj(1, 1) = dyc(i, j)*sina_v(i, j)*delpc_tj(1, 1) - vort_tj(1, 0)
        END IF
        IF (se_corner) THEN
          delpc_tl(npx, 1) = delpc_tl(npx, 1) - vort_tl(npx, 0)
          delpc_tj(npx, 1) = delpc_tj(npx, 1) - vort_tj(npx, 0)
        END IF
        IF (ne_corner) THEN
          delpc_tl(npx, npy) = delpc_tl(npx, npy) + vort_tl(npx, npy)
          delpc_tj(npx, npy) = delpc_tj(npx, npy) + vort_tj(npx, npy)
        END IF
        IF (nw_corner) THEN
          delpc_tl(1, npy) = delpc_tl(1, npy) + vort_tl(1, npy)
          delpc_tj(1, npy) = delpc_tj(1, npy) + vort_tj(1, npy)
        END IF
        DO j=js,je+1
          DO i=is,ie+1
            delpc_tl(i, j) = rarea_c(i, j)*delpc_tl(i, j)
            delpc_tj(i, j) = rarea_c(i, j)*delpc_tj(i, j)
            IF (delpc_tj(i, j)*dt .GE. 0.) THEN
              abs0_tl = dt*delpc_tl(i, j)
              abs0_tj = dt*delpc_tj(i, j)
            ELSE
              abs0_tl = -(dt*delpc_tl(i, j))
              abs0_tj = -(dt*delpc_tj(i, j))
            END IF
            y3_tl = dddmp_tj*abs0_tl
            y3_tj = dddmp_tj*abs0_tj
            IF (0.20 .GT. y3_tj) THEN
              y1_tj = y3_tj
              y1_tl = y3_tl
            ELSE
              y1_tj = 0.20
              y1_tl = 0.0
            END IF
            IF (d2_bg_tj .LT. y1_tj) THEN
              max5_tj = y1_tj
              max5_tl = y1_tl
            ELSE
              max5_tj = d2_bg_tj
              max5_tl = 0.0
            END IF
            damp_tl = da_min_c*max5_tl
            damp_tj = da_min_c*max5_tj
            vort_tl(i, j) = damp_tl*delpc_tj(i, j) + damp_tj*delpc_tl(i, j)
            vort_tj(i, j) = damp_tj*delpc_tj(i, j)
            ke_tl(i, j) = ke_tl(i, j) + vort_tl(i, j)
            ke_tj(i, j) = ke_tj(i, j) + vort_tj(i, j)
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
            delpc_tj(i, j) = divg_d_tj(i, j)
          END DO
        END DO
        n2 = nord_tj + 1
        DO n=1,nord_tj
          nt = nord_tj - n
          fill_c = nt .NE. 0 .AND. grid_type .LT. 3 .AND. (((sw_corner &
&           .OR. se_corner) .OR. ne_corner) .OR. nw_corner)
          DO j=js-nt,je+1+nt
            DO i=is-1-nt,ie+1+nt
              vc_tl(i, j) = divg_u(i, j)*(divg_d_tl(i+1, j)-divg_d_tl(i, j))
              vc_tj(i, j) = divg_u(i, j)*(divg_d_tj(i+1, j)-divg_d_tj(i, j))
            END DO
          END DO
          DO j=js-1-nt,je+1+nt
            DO i=is-nt,ie+1+nt
              uc_tl(i, j) = divg_v(i, j)*(divg_d_tl(i, j+1)-divg_d_tl(i, j))
              uc_tj(i, j) = divg_v(i, j)*(divg_d_tj(i, j+1)-divg_d_tj(i, j))
            END DO
          END DO
          DO j=js-nt,je+1+nt
            DO i=is-nt,ie+1+nt
              divg_d_tl(i, j) = uc_tl(i, j-1) - uc_tl(i, j) + vc_tl(i-1&
&               , j) - vc_tl(i, j)
              divg_d_tj(i, j) = uc_tj(i, j-1) - uc_tj(i, j) + vc_tj(i-1&
&               , j) - vc_tj(i, j)
            END DO
          END DO
! Remove the extra term at the corners:
          IF (sw_corner) THEN
            divg_d_tl(1, 1) = divg_d_tl(1, 1) - uc_tl(1, 0)
            divg_d_tj(1, 1) = divg_d_tj(1, 1) - uc_tj(1, 0)
          END IF
          IF (se_corner) THEN
            divg_d_tl(npx, 1) = divg_d_tl(npx, 1) - uc_tl(npx, 0)
            divg_d_tj(npx, 1) = divg_d_tj(npx, 1) - uc_tj(npx, 0)
          END IF
          IF (ne_corner) THEN
            divg_d_tl(npx, npy) = divg_d_tl(npx, npy) + uc_tl(npx, npy)
            divg_d_tj(npx, npy) = divg_d_tj(npx, npy) + uc_tj(npx, npy)
          END IF
          IF (nw_corner) THEN
            divg_d_tl(1, npy) = divg_d_tl(1, npy) + uc_tl(1, npy)
            divg_d_tj(1, npy) = divg_d_tj(1, npy) + uc_tj(1, npy)
          END IF
          DO j=js-nt,je+1+nt
            DO i=is-nt,ie+1+nt
              divg_d_tl(i, j) = rarea_c(i, j)*divg_d_tl(i, j)
              divg_d_tj(i, j) = rarea_c(i, j)*divg_d_tj(i, j)
            END DO
          END DO
        END DO
        IF (dddmp_tj .LT. 1.e-5) THEN
          vort_tj = 0.
          vort_tl = 0.0
        ELSE
          vort_tl = 0.0
! Compute "time-scale" for del-2 background damping
          DO j=js,je+1
            DO i=is,ie+1
              IF (dt*delpc_tj(i, j) .GE. 0.) THEN
                vort_tj(i, j) = dt*delpc_tj(i, j)
                vort_tl(i, j) = dt*delpc_tl(i, j)
              ELSE
                vort_tj(i, j) = -(dt*delpc_tj(i, j))
                vort_tl(i, j) = -(dt*delpc_tl(i, j))
              END IF
            END DO
          END DO
        END IF
        pwx1_tj = da_min_c*d4_bg_tj
        dd8_tj = pwx1_tj**n2
        DO j=js,je+1
          DO i=is,ie+1
            IF (0.20 .GT. dddmp_tj*vort_tj(i, j)) THEN
              y2_tj = dddmp_tj*vort_tj(i, j)
              y2_tl = dddmp_tj*vort_tl(i, j)
            ELSE
              y2_tj = 0.20
              y2_tl = 0.0
            END IF
            IF (d2_bg_tj .LT. y2) THEN
              max6_tj = y2_tj
              max6_tl = y2_tl
            ELSE
              max6_tj = d2_bg_tj
              max6_tl = 0.0
            END IF
! del-2
            damp2_tl = da_min_c*max6_tl
            damp2_tj = da_min_c*max6_tj
            vort_tl(i, j) = damp2_tl*delpc_tj(i, j) + damp2_tj*delpc_tl(i, j) + dd8_tj*divg_d_tl(i, j)
            vort_tj(i, j) = damp2_tj*delpc_tj(i, j) +                           dd8_tj*divg_d_tj(i, j)
            ke_tl(i, j) = ke_tl(i, j) + vort_tl(i, j)
            ke_tj(i, j) = ke_tj(i, j) + vort_tj(i, j)

          END DO
        END DO

      END IF


! Recompute for TRACETORY so damping is the same as NLM
! -----------------------------------------------------

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
          IF (is .EQ. 1) THEN
            vort(1, j) = v(1, j)*dxc(1, j)*sina_u(1, j)
          END IF
          IF (ie + 1 .EQ. npx) THEN
            vort(npx, j) = v(npx, j)*dxc(npx, j)*sina_u(npx, j)
          END IF
        END DO
        DO j=js,je+1
          DO i=is,ie+1
            delpc(i, j) = vort(i, j-1) - vort(i, j) + ptc(i-1, j) - ptc(&
&             i, j)
          END DO
        END DO
! Remove the extra term at the corners:
        IF (sw_corner) THEN
          delpc(1, 1) = delpc(1, 1) - vort(1, 0)
        END IF
        IF (se_corner) THEN
          delpc(npx, 1) = delpc(npx, 1) - vort(npx, 0)
        END IF
        IF (ne_corner) THEN
          delpc(npx, npy) = delpc(npx, npy) + vort(npx, npy)
        END IF
        IF (nw_corner) THEN
          delpc(1, npy) = delpc(1, npy) + vort(1, npy)
        END IF
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
!        if ( fill_c ) call fill_corners(divg_d, npx, npy, FILL=XDir, BGRID=.true.)
          DO j=js-nt,je+1+nt
            DO i=is-1-nt,ie+1+nt
              vc(i, j) = (divg_d(i+1, j)-divg_d(i, j))*divg_u(i, j)
            END DO
          END DO
!        if ( fill_c ) call fill_corners(divg_d, npx, npy, FILL=YDir, BGRID=.true.)
          DO j=js-1-nt,je+1+nt
            DO i=is-nt,ie+1+nt
              uc(i, j) = (divg_d(i, j+1)-divg_d(i, j))*divg_v(i, j)
            END DO
          END DO
!        if ( fill_c ) call fill_corners(vc, uc, npx, npy, VECTOR=.true., DGRID=.true.)
          DO j=js-nt,je+1+nt
            DO i=is-nt,ie+1+nt
              divg_d(i, j) = uc(i, j-1) - uc(i, j) + vc(i-1, j) - vc(i, &
&               j)
            END DO
          END DO
! Remove the extra term at the corners:
          IF (sw_corner) THEN
            divg_d(1, 1) = divg_d(1, 1) - uc(1, 0)
          END IF
          IF (se_corner) THEN
            divg_d(npx, 1) = divg_d(npx, 1) - uc(npx, 0)
          END IF
          IF (ne_corner) THEN
            divg_d(npx, npy) = divg_d(npx, npy) + uc(npx, npy)
          END IF
          IF (nw_corner) THEN
            divg_d(1, npy) = divg_d(1, npy) + uc(1, npy)
          END IF
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
        pwx1 = da_min_c*d4_bg
        dd8 = pwx1**n2
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
        gy_tl = 0.0
        DO j=js,je+1
          DO i=is,ie
! du
            ub_tl(i, j) = rdx(i, j)*(vort_tl(i, j)-vort_tl(i+1, j))
            ub_tj(i, j) = rdx(i, j)*(vort_tj(i, j)-vort_tj(i+1, j))
            ub(i, j) = (vort(i, j)-vort(i+1, j))*rdx(i, j)
! u*du
            gy_tl(i, j) = u_tl(i, j)*ub_tj(i, j) + u(i, j)*ub_tl(i, j)
            gy(i, j) = u(i, j)*ub(i, j)
          END DO
        END DO
        gx_tl = 0.0
        DO j=js,je
          DO i=is,ie+1
! dv
            vb_tl(i, j) = rdy(i, j)*(vort_tl(i, j)-vort_tl(i, j+1))
            vb_tj(i, j) = (vort_tj(i, j)-vort_tj(i, j+1))*rdy(i, j)
            vb(i, j) = (vort(i, j)-vort(i, j+1))*rdy(i, j)
! v*dv
            gx_tl(i, j) = v_tl(i, j)*vb_tj(i, j) + v(i, j)*vb_tl(i, j)
            gx(i, j) = v(i, j)*vb(i, j)
          END DO
        END DO
        DO j=js,je
          DO i=is,ie
            u2_tl = u_tl(i, j) + u_tl(i, j+1)
            u2 = u(i, j) + u(i, j+1)
            du2_tl = ub_tl(i, j) + ub_tl(i, j+1)
            du2_tj = ub_tj(i, j) + ub_tj(i, j+1)
            du2 = ub(i, j) + ub(i, j+1)
            v2_tl = v_tl(i, j) + v_tl(i+1, j)
            v2 = v(i, j) + v(i+1, j)
            dv2_tl = vb_tl(i, j) + vb_tl(i+1, j)
            dv2_tj = vb_tj(i, j) + vb_tj(i+1, j)
            dv2 = vb(i, j) + vb(i+1, j)
! Total energy conserving:
! Convert lost KE due to divergence damping to "heat"
            pt_tl(i, j) = pt_tl(i, j) - damp_tj*rsin2(i, j)*(2*ub(i, j)*&
&             ub_tl(i, j)+2*ub(i, j+1)*ub_tl(i, j+1)+2*vb(i, j)*vb_tl(i&
&             , j)+2*vb(i+1, j)*vb_tl(i+1, j)+2.*(gy_tl(i, j)+gy_tl(i, j&
&             +1)+gx_tl(i, j)+gx_tl(i+1, j))-cosa_s(i, j)*(u2_tl*dv2_tj+u2*&
&             dv2_tl+v2_tl*du2_tj+v2*du2_tl+du2_tl*dv2_tj+du2_tj*dv2_tl))/pkz(i, &
&             j) + damp_tj*rsin2(i, j)*pkz_tl(i, j)*(ub(i, j)**2+ub(i, j+1)&
&             **2+vb(i, j)**2+vb(i+1, j)**2+2.*(gy(i, j)+gy(i, j+1)+gx(i&
&             , j)+gx(i+1, j))-cosa_s(i, j)*(u2*dv2_tj+v2*du2_tj+du2_tj*dv2_tj))/pkz&
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
    fxtmp = 0.0
    fytmp = 0.0
    fieldtmp = vort
      CALL FV_TP_2D_TLM(fieldtmp, vort_tl, crx_adv, crx_adv_tl, cry_adv, &
&                 cry_adv_tl, npx, npy, hord_vt_tj, fxtmp, fx_tl, fytmp, fy_tl, &
&                 xfx_adv, xfx_adv_tl, yfx_adv, yfx_adv_tl, area, ra_x, &
&                 ra_x_tl, ra_y, ra_y_tl, ppm_fac=ppm_limiter, nord=nord_tj&
&                 , damp_c=vtdm4)
      CALL FV_TP_2D(vort, crx_adv, cry_adv, npx, npy, hord_vt, fx, fy, &
&                 xfx_adv, yfx_adv, area, ra_x, &
&                 ra_y, ppm_fac=ppm_limiter, nord=nord, damp_c=vtdm4)
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
    uf_tl = 0.0
    vf_tl = 0.0
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

  SUBROUTINE DIVERGENCE_CORNER_ONLYTLM(u_tl, v_tl, ua_tl, va_tl, divg_d_tl, km)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: km
    REAL, DIMENSION(isd:ied, jsd:jed+1, km), INTENT(IN) :: u_tl
    REAL, DIMENSION(isd:ied+1, jsd:jed, km), INTENT(IN) :: v_tl
    REAL, DIMENSION(isd:ied, jsd:jed, km), INTENT(IN) :: ua_tl, va_tl
    REAL, DIMENSION(isd:ied+1, jsd:jed+1, km), INTENT(OUT) :: divg_d_tl
! local
    REAL :: uf_tl(is-2:ie+2, js-1:je+2)
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
    uf_tl = 0.0
    vf_tl = 0.0
    DO k=1,km
      IF (grid_type .EQ. 4) THEN
        DO j=js-1,je+2
          DO i=is-2,ie+2
            uf_tl(i, j) = dyc(i, j)*u_tl(i, j, k)
          END DO
        END DO
        DO j=js-2,je+2
          DO i=is-1,ie+2
            vf_tl(i, j) = dxc(i, j)*v_tl(i, j, k)
          END DO
        END DO
        DO j=js-1,je+2
          DO i=is-1,ie+2
            divg_d_tl(i, j, k) = rarea_c(i, j)*(vf_tl(i, j-1)-vf_tl(i, j&
&             )+uf_tl(i-1, j)-uf_tl(i, j))
          END DO
        END DO
      ELSE
! divg_u(i,j) = sina_v(i,j)*dyc(i,j)/dx(i,j)
        DO j=js,je+1
          IF (j .EQ. 1 .OR. j .EQ. npy) THEN
            DO i=is-1,ie+1
              uf_tl(i, j) = dyc(i, j)*sina_v(i, j)*u_tl(i, j, k)
            END DO
          ELSE
            DO i=is-1,ie+1
              uf_tl(i, j) = dyc(i, j)*sina_v(i, j)*(u_tl(i, j, k)-0.5*&
&               cosa_v(i, j)*(va_tl(i, j-1, k)+va_tl(i, j, k)))
            END DO
          END IF
        END DO
        DO j=js-1,je+1
          DO i=is2,ie1
            vf_tl(i, j) = dxc(i, j)*sina_u(i, j)*(v_tl(i, j, k)-0.5*&
&             cosa_u(i, j)*(ua_tl(i-1, j, k)+ua_tl(i, j, k)))
          END DO
          IF (is .EQ. 1) THEN
            vf_tl(1, j) = dxc(1, j)*sina_u(1, j)*v_tl(1, j, k)
          END IF
          IF (ie + 1 .EQ. npx) THEN
            vf_tl(npx, j) = dxc(npx, j)*sina_u(npx, j)*v_tl(npx, j, k)
          END IF
        END DO
        DO j=js,je+1
          DO i=is,ie+1
            divg_d_tl(i, j, k) = vf_tl(i, j-1) - vf_tl(i, j) + uf_tl(i-1&
&             , j) - uf_tl(i, j)
          END DO
        END DO
! Remove the extra term at the corners:
        IF (sw_corner) THEN
          divg_d_tl(1, 1, k) = divg_d_tl(1, 1, k) - vf_tl(1, 0)
        END IF
        IF (se_corner) THEN
          divg_d_tl(npx, 1, k) = divg_d_tl(npx, 1, k) - vf_tl(npx, 0)
        END IF
        IF (ne_corner) THEN
          divg_d_tl(npx, npy, k) = divg_d_tl(npx, npy, k) + vf_tl(npx, &
&           npy)
        END IF
        IF (nw_corner) THEN
          divg_d_tl(1, npy, k) = divg_d_tl(1, npy, k) + vf_tl(1, npy)
        END IF
        DO j=js,je+1
          DO i=is,ie+1
            divg_d_tl(i, j, k) = rarea_c(i, j)*divg_d_tl(i, j, k)
          END DO
        END DO
      END IF
    END DO
  END SUBROUTINE DIVERGENCE_CORNER_ONLYTLM

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
    REAL(ke_precision) :: al_tl(is-1:ie+2), dm_tl(is-2:ie+2)
    REAL(ke_precision) :: bl(is-1:ie+1)
    REAL(ke_precision) :: bl_tl(is-1:ie+1)
    REAL(ke_precision) :: br(is-1:ie+1)
    REAL(ke_precision) :: br_tl(is-1:ie+1)
    REAL(ke_precision) :: dq(is-3:ie+2)
    REAL(ke_precision) :: dq_tl(is-3:ie+2)
    REAL(ke_precision) :: dl, dr, xt, pmp, lac, dqt, cfl
    REAL(ke_precision) :: dl_tl, dr_tl, xt_tl, pmp_tl, lac_tl, cfl_tl
    REAL(ke_precision) :: x0, x1
    REAL(ke_precision) :: x0_tl, x1_tl
    INTEGER :: i, j
    INTRINSIC ABS
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC SIGN
    REAL*8 :: y1_tl
    REAL*8 :: min6_tl
    REAL*8 :: y18_tl
    REAL :: min22_tl
    REAL*8 :: x18
    REAL*8 :: x17
    REAL*8 :: min10_tl
    REAL*8 :: x16
    REAL*8 :: z4_tl
    REAL*8 :: x13_tl
    REAL*8 :: x15
    REAL :: max7_tl
    REAL*8 :: y20_tl
    REAL*8 :: x14
    REAL*8 :: min9
    REAL*8 :: x13
    REAL*8 :: min8
    REAL*8 :: y2_tl
    REAL*8 :: min7_tl
    REAL*8 :: y19_tl
    REAL*8 :: x12
    REAL*8 :: min7
    REAL*8 :: x11
    REAL*8 :: min6
    REAL*8 :: x10
    REAL*8 :: min5
    REAL*8 :: min4
    REAL*8 :: min3
    REAL*8 :: min11_tl
    REAL*8 :: min2
    REAL*8 :: z5_tl
    REAL*8 :: x14_tl
    REAL*8 :: min1
    REAL*8 :: max8_tl
    REAL*8 :: y21_tl
    REAL*8 :: y3_tl
    REAL*8 :: min8_tl
    REAL*8 :: y21
    REAL*8 :: y20
    REAL*8 :: x9
    REAL*8 :: x8
    REAL*8 :: min12_tl
    REAL*8 :: x7
    REAL*8 :: z6_tl
    REAL*8 :: x15_tl
    REAL*8 :: x6
    REAL*8 :: max9_tl
    REAL*8 :: x5
    REAL*8 :: x4
    REAL*8 :: y4_tl
    REAL*8 :: min9_tl
    REAL*8 :: x3
    REAL*8 :: y10_tl
    REAL*8 :: x2
    REAL*8 :: x2_tl
    REAL*8 :: min13_tl
    REAL*8 :: z7_tl
    REAL*8 :: x16_tl
    REAL*8 :: y5_tl
    REAL*8 :: y11_tl
    REAL*8 :: x3_tl
    INTEGER :: min26
    REAL :: min14_tl
    REAL*8 :: y19
    INTEGER :: min25
    REAL*8 :: x17_tl
    REAL*8 :: y18
    INTEGER :: min24
    REAL*8 :: y17
    INTEGER :: min23
    REAL*8 :: y16
    REAL :: min22
    REAL*8 :: y6_tl
    REAL*8 :: y15
    INTEGER :: min21
    REAL*8 :: y12_tl
    REAL*8 :: y14
    INTEGER :: min20
    REAL*8 :: y13
    REAL*8 :: x4_tl
    REAL*8 :: y12
    REAL*8 :: y11
    REAL*8 :: min15_tl
    REAL*8 :: y10
    REAL*8 :: max10_tl
    REAL*8 :: x18_tl
    REAL*8 :: y7_tl
    REAL*8 :: min1_tl
    REAL*8 :: y13_tl
    REAL*8 :: x5_tl
    REAL*8 :: min16_tl
    REAL*8 :: max11_tl
    REAL*8 :: z7
    REAL*8 :: z6
    REAL*8 :: y8_tl
    REAL*8 :: z5
    REAL*8 :: min2_tl
    REAL*8 :: y14_tl
    REAL*8 :: z4
    REAL*8 :: z3
    REAL*8 :: x6_tl
    REAL*8 :: z2
    REAL*8 :: z1
    REAL*8 :: min17_tl
    REAL :: min19
    REAL :: max12_tl
    REAL*8 :: min18
    REAL*8 :: min17
    REAL*8 :: min16
    REAL*8 :: y9_tl
    REAL*8 :: min15
    REAL*8 :: min3_tl
    REAL*8 :: y15_tl
    REAL :: min14
    REAL*8 :: min13
    REAL*8 :: x7_tl
    REAL*8 :: min12
    REAL*8 :: min11
    REAL*8 :: min18_tl
    REAL*8 :: min10
    REAL :: max13_tl
    REAL*8 :: z1_tl
    REAL*8 :: x10_tl
    REAL*8 :: min4_tl
    REAL*8 :: y16_tl
    REAL*8 :: x8_tl
    REAL*8 :: max9
    REAL :: min19_tl
    REAL*8 :: max8
    REAL :: max7
    REAL*8 :: z2_tl
    REAL*8 :: x11_tl
    INTEGER :: max6
    INTEGER :: max5
    REAL*8 :: y9
    INTEGER :: max4
    REAL*8 :: min5_tl
    REAL*8 :: y17_tl
    REAL*8 :: y8
    INTEGER :: max3
    REAL*8 :: y7
    INTEGER :: max2
    REAL :: max13
    REAL*8 :: x9_tl
    REAL*8 :: y6
    INTEGER :: max1
    REAL :: max12
    REAL*8 :: y5
    REAL*8 :: max11
    REAL*8 :: y4
    REAL*8 :: max10
    REAL*8 :: y3
    REAL*8 :: z3_tl
    REAL*8 :: x12_tl
    REAL*8 :: y2
    REAL*8 :: y1
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

      flux_tl = 0.0

      !First order
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

      !Upgrade to third order scheme internally
      DO j=js,je+1
        DO i=is+1,ie
          IF (c(i, j) .GT. 0.) THEN
            cfl = c(i, j)*rdx(i-1, j)
            cfl_tl = c_tl(i, j)*rdx(i-1, j)
            flux_tl(i, j) = (2.0*u_tl(i, j)+5.0*u_tl(i-1, j)-u_tl(i-2, j&
&             ))/6.0 - 0.5*(cfl_tl*(u(i, j)-u(i-1, j))+cfl*(u_tl&
&             (i, j)-u_tl(i-1, j))) + (cfl_tl*cfl+cfl*cfl_tl&
&             )*(u(i, j)-2.0*u(i-1, j)+u(i-2, j))/6.0 + cfl**2*(&
&             u_tl(i, j)-2.0*u_tl(i-1, j)+u_tl(i-2, j))/6.0
            flux(i, j) = (2.0*u(i, j)+5.0*u(i-1, j)-u(i-2, j))/6.0 - 0.5&
&             *cfl*(u(i, j)-u(i-1, j)) + cfl*cfl/6.0*(u(i, j&
&             )-2.0*u(i-1, j)+u(i-2, j))
          ELSE
            cfl = c(i, j)*rdx(i, j)
            cfl_tl = c_tl(i, j)*rdx(i, j)
            flux_tl(i, j) = (2.0*u_tl(i-1, j)+5.0*u_tl(i, j)-u_tl(i+1, j&
&             ))/6.0 - 0.5*(cfl_tl*(u(i, j)-u(i-1, j))+cfl*(u_tl&
&             (i, j)-u_tl(i-1, j))) + (cfl_tl*cfl+cfl*cfl_tl&
&             )*(u(i+1, j)-2.0*u(i, j)+u(i-1, j))/6.0 + cfl**2*(&
&             u_tl(i+1, j)-2.0*u_tl(i, j)+u_tl(i-1, j))/6.0
            flux(i, j) = (2.0*u(i-1, j)+5.0*u(i, j)-u(i+1, j))/6.0 - 0.5&
&             *cfl*(u(i, j)-u(i-1, j)) + cfl*cfl/6.0*(u(i+1&
&             , j)-2.0*u(i, j)+u(i-1, j))
          END IF
        END DO
      END DO

      !Downgrade to first order for the edges
      IF (grid_type .LT. 3) THEN
        IF (is .EQ. 1) THEN
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
        END IF
        IF (ie + 1 .EQ. npx) THEN
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
        END IF
      END IF

    CASE (2) 
      dm_tl = 0.0_8
      DO j=js,je+1
        DO i=is-2,ie+2
          xt_tl = 0.25*(u_tl(i+1, j)-u_tl(i-1, j))
          xt = 0.25*(u(i+1, j)-u(i-1, j))
          IF (xt .GE. 0.) THEN
            x2_tl = xt_tl
            x2 = xt
          ELSE
            x2_tl = -xt_tl
            x2 = -xt
          END IF
          IF (u(i-1, j) .LT. u(i, j)) THEN
            IF (u(i, j) .LT. u(i+1, j)) THEN
              max7_tl = u_tl(i+1, j)
              max7 = u(i+1, j)
            ELSE
              max7_tl = u_tl(i, j)
              max7 = u(i, j)
            END IF
          ELSE IF (u(i-1, j) .LT. u(i+1, j)) THEN
            max7_tl = u_tl(i+1, j)
            max7 = u(i+1, j)
          ELSE
            max7_tl = u_tl(i-1, j)
            max7 = u(i-1, j)
          END IF
          y1_tl = max7_tl - u_tl(i, j)
          y1 = max7 - u(i, j)
          IF (u(i-1, j) .GT. u(i, j)) THEN
            IF (u(i, j) .GT. u(i+1, j)) THEN
              min14_tl = u_tl(i+1, j)
              min14 = u(i+1, j)
            ELSE
              min14_tl = u_tl(i, j)
              min14 = u(i, j)
            END IF
          ELSE IF (u(i-1, j) .GT. u(i+1, j)) THEN
            min14_tl = u_tl(i+1, j)
            min14 = u(i+1, j)
          ELSE
            min14_tl = u_tl(i-1, j)
            min14 = u(i-1, j)
          END IF
          z1_tl = u_tl(i, j) - min14_tl
          z1 = u(i, j) - min14
          IF (x2 .GT. y1) THEN
            IF (y1 .GT. z1) THEN
              min1_tl = z1_tl
              min1 = z1
            ELSE
              min1_tl = y1_tl
              min1 = y1
            END IF
          ELSE IF (x2 .GT. z1) THEN
            min1_tl = z1_tl
            min1 = z1
          ELSE
            min1_tl = x2_tl
            min1 = x2
          END IF
          dm_tl(i) = min1_tl*SIGN(1.d0, min1*xt)
          dm(i) = SIGN(min1, xt)
        END DO
! Fix slopes near edges:
        IF (grid_type .LT. 3) THEN
          IF (is .EQ. 1) THEN
            IF (j .EQ. 1 .OR. j .EQ. npy) THEN
              dm_tl(0) = 0.0_8
              dm(0) = 0.
              dm_tl(1) = 0.0_8
              dm(1) = 0.
            ELSE
              x0_tl = 0.5*((2.*dx(1, j)+dx(2, j))*(u_tl(0, j)+u_tl(1, j)&
&               )-dx(1, j)*(u_tl(-1, j)+u_tl(2, j)))/(dx(1, j)+dx(2, j))
              x0 = 0.5*((2.*dx(1, j)+dx(2, j))*(u(0, j)+u(1, j))-dx(1, j&
&               )*(u(-1, j)+u(2, j)))/(dx(1, j)+dx(2, j))
              x1_tl = s15*u_tl(0, j) + s11*u_tl(-1, j) + s14*dm_tl(-1)
              x1 = s15*u(0, j) + s11*u(-1, j) + s14*dm(-1)
!          dm(0) = u(0,j) - x1
              dm_tl(0) = 0.5*(x0_tl-x1_tl)
              dm(0) = 0.5*(x0-x1)
              IF (dm(0) .GE. 0.) THEN
                x3_tl = dm_tl(0)
                x3 = dm(0)
              ELSE
                x3_tl = -dm_tl(0)
                x3 = -dm(0)
              END IF
              IF (u(0, j) .LT. x0) THEN
                IF (x0 .LT. x1) THEN
                  max8_tl = x1_tl
                  max8 = x1
                ELSE
                  max8_tl = x0_tl
                  max8 = x0
                END IF
              ELSE IF (u(0, j) .LT. x1) THEN
                max8_tl = x1_tl
                max8 = x1
              ELSE
                max8_tl = u_tl(0, j)
                max8 = u(0, j)
              END IF
              y2_tl = max8_tl - u_tl(0, j)
              y2 = max8 - u(0, j)
              IF (u(0, j) .GT. x0) THEN
                IF (x0 .GT. x1) THEN
                  min15_tl = x1_tl
                  min15 = x1
                ELSE
                  min15_tl = x0_tl
                  min15 = x0
                END IF
              ELSE IF (u(0, j) .GT. x1) THEN
                min15_tl = x1_tl
                min15 = x1
              ELSE
                min15_tl = u_tl(0, j)
                min15 = u(0, j)
              END IF
              z2_tl = u_tl(0, j) - min15_tl
              z2 = u(0, j) - min15
              IF (x3 .GT. y2) THEN
                IF (y2 .GT. z2) THEN
                  min2_tl = z2_tl
                  min2 = z2
                ELSE
                  min2_tl = y2_tl
                  min2 = y2
                END IF
              ELSE IF (x3 .GT. z2) THEN
                min2_tl = z2_tl
                min2 = z2
              ELSE
                min2_tl = x3_tl
                min2 = x3
              END IF
              dm_tl(0) = min2_tl*SIGN(1.d0, min2*dm(0))
              dm(0) = SIGN(min2, dm(0))
              x1_tl = s15*u_tl(1, j) + s11*u_tl(2, j) - s14*dm_tl(2)
              x1 = s15*u(1, j) + s11*u(2, j) - s14*dm(2)
!          dm(1) = x1 - u(1,j)
              dm_tl(1) = 0.5*(x1_tl-x0_tl)
              dm(1) = 0.5*(x1-x0)
              IF (dm(1) .GE. 0.) THEN
                x4_tl = dm_tl(1)
                x4 = dm(1)
              ELSE
                x4_tl = -dm_tl(1)
                x4 = -dm(1)
              END IF
              IF (u(1, j) .LT. x0) THEN
                IF (x0 .LT. x1) THEN
                  max9_tl = x1_tl
                  max9 = x1
                ELSE
                  max9_tl = x0_tl
                  max9 = x0
                END IF
              ELSE IF (u(1, j) .LT. x1) THEN
                max9_tl = x1_tl
                max9 = x1
              ELSE
                max9_tl = u_tl(1, j)
                max9 = u(1, j)
              END IF
              y3_tl = max9_tl - u_tl(1, j)
              y3 = max9 - u(1, j)
              IF (u(1, j) .GT. x0) THEN
                IF (x0 .GT. x1) THEN
                  min16_tl = x1_tl
                  min16 = x1
                ELSE
                  min16_tl = x0_tl
                  min16 = x0
                END IF
              ELSE IF (u(1, j) .GT. x1) THEN
                min16_tl = x1_tl
                min16 = x1
              ELSE
                min16_tl = u_tl(1, j)
                min16 = u(1, j)
              END IF
              z3_tl = u_tl(1, j) - min16_tl
              z3 = u(1, j) - min16
              IF (x4 .GT. y3) THEN
                IF (y3 .GT. z3) THEN
                  min3_tl = z3_tl
                  min3 = z3
                ELSE
                  min3_tl = y3_tl
                  min3 = y3
                END IF
              ELSE IF (x4 .GT. z3) THEN
                min3_tl = z3_tl
                min3 = z3
              ELSE
                min3_tl = x4_tl
                min3 = x4
              END IF
              dm_tl(1) = min3_tl*SIGN(1.d0, min3*dm(1))
              dm(1) = SIGN(min3, dm(1))
            END IF
          END IF
          IF (ie + 1 .EQ. npx) THEN
            IF (j .EQ. 1 .OR. j .EQ. npy) THEN
              dm_tl(npx-1) = 0.0_8
              dm(npx-1) = 0.
              dm_tl(npx) = 0.0_8
              dm(npx) = 0.
            ELSE
              x0_tl = 0.5*((2.*dx(npx-1, j)+dx(npx-2, j))*(u_tl(npx-1, j&
&               )+u_tl(npx, j))-dx(npx-1, j)*(u_tl(npx-2, j)+u_tl(npx+1&
&               , j)))/(dx(npx-1, j)+dx(npx-2, j))
              x0 = 0.5*((2.*dx(npx-1, j)+dx(npx-2, j))*(u(npx-1, j)+u(&
&               npx, j))-dx(npx-1, j)*(u(npx-2, j)+u(npx+1, j)))/(dx(npx&
&               -1, j)+dx(npx-2, j))
              x1_tl = s15*u_tl(npx-1, j) + s11*u_tl(npx-2, j) + s14*&
&               dm_tl(npx-2)
              x1 = s15*u(npx-1, j) + s11*u(npx-2, j) + s14*dm(npx-2)
!          dm(npx-1) = u(npx-1,j) - x1
              dm_tl(npx-1) = 0.5*(x0_tl-x1_tl)
              dm(npx-1) = 0.5*(x0-x1)
              IF (dm(npx-1) .GE. 0.) THEN
                x5_tl = dm_tl(npx-1)
                x5 = dm(npx-1)
              ELSE
                x5_tl = -dm_tl(npx-1)
                x5 = -dm(npx-1)
              END IF
              IF (u(npx-1, j) .LT. x0) THEN
                IF (x0 .LT. x1) THEN
                  max10_tl = x1_tl
                  max10 = x1
                ELSE
                  max10_tl = x0_tl
                  max10 = x0
                END IF
              ELSE IF (u(npx-1, j) .LT. x1) THEN
                max10_tl = x1_tl
                max10 = x1
              ELSE
                max10_tl = u_tl(npx-1, j)
                max10 = u(npx-1, j)
              END IF
              y4_tl = max10_tl - u_tl(npx-1, j)
              y4 = max10 - u(npx-1, j)
              IF (u(npx-1, j) .GT. x0) THEN
                IF (x0 .GT. x1) THEN
                  min17_tl = x1_tl
                  min17 = x1
                ELSE
                  min17_tl = x0_tl
                  min17 = x0
                END IF
              ELSE IF (u(npx-1, j) .GT. x1) THEN
                min17_tl = x1_tl
                min17 = x1
              ELSE
                min17_tl = u_tl(npx-1, j)
                min17 = u(npx-1, j)
              END IF
              z4_tl = u_tl(npx-1, j) - min17_tl
              z4 = u(npx-1, j) - min17
              IF (x5 .GT. y4) THEN
                IF (y4 .GT. z4) THEN
                  min4_tl = z4_tl
                  min4 = z4
                ELSE
                  min4_tl = y4_tl
                  min4 = y4
                END IF
              ELSE IF (x5 .GT. z4) THEN
                min4_tl = z4_tl
                min4 = z4
              ELSE
                min4_tl = x5_tl
                min4 = x5
              END IF
              dm_tl(npx-1) = min4_tl*SIGN(1.d0, min4*dm(npx-1))
              dm(npx-1) = SIGN(min4, dm(npx-1))
              x1_tl = s15*u_tl(npx, j) + s11*u_tl(npx+1, j) - s14*dm_tl(&
&               npx+1)
              x1 = s15*u(npx, j) + s11*u(npx+1, j) - s14*dm(npx+1)
!          dm(npx) = x1 - u(npx,j)
              dm_tl(npx) = 0.5*(x1_tl-x0_tl)
              dm(npx) = 0.5*(x1-x0)
              IF (dm(npx) .GE. 0.) THEN
                x6_tl = dm_tl(npx)
                x6 = dm(npx)
              ELSE
                x6_tl = -dm_tl(npx)
                x6 = -dm(npx)
              END IF
              IF (u(npx, j) .LT. x0) THEN
                IF (x0 .LT. x1) THEN
                  max11_tl = x1_tl
                  max11 = x1
                ELSE
                  max11_tl = x0_tl
                  max11 = x0
                END IF
              ELSE IF (u(npx, j) .LT. x1) THEN
                max11_tl = x1_tl
                max11 = x1
              ELSE
                max11_tl = u_tl(npx, j)
                max11 = u(npx, j)
              END IF
              y5_tl = max11_tl - u_tl(npx, j)
              y5 = max11 - u(npx, j)
              IF (u(npx, j) .GT. x0) THEN
                IF (x0 .GT. x1) THEN
                  min18_tl = x1_tl
                  min18 = x1
                ELSE
                  min18_tl = x0_tl
                  min18 = x0
                END IF
              ELSE IF (u(npx, j) .GT. x1) THEN
                min18_tl = x1_tl
                min18 = x1
              ELSE
                min18_tl = u_tl(npx, j)
                min18 = u(npx, j)
              END IF
              z5_tl = u_tl(npx, j) - min18_tl
              z5 = u(npx, j) - min18
              IF (x6 .GT. y5) THEN
                IF (y5 .GT. z5) THEN
                  min5_tl = z5_tl
                  min5 = z5
                ELSE
                  min5_tl = y5_tl
                  min5 = y5
                END IF
              ELSE IF (x6 .GT. z5) THEN
                min5_tl = z5_tl
                min5 = z5
              ELSE
                min5_tl = x6_tl
                min5 = x6
              END IF
              dm_tl(npx) = min5_tl*SIGN(1.d0, min5*dm(npx))
              dm(npx) = SIGN(min5, dm(npx))
            END IF
          END IF
        END IF
        DO i=is,ie+1
          IF (c(i, j) .GT. 0.) THEN
            flux_tl(i, j) = u_tl(i-1, j) + (1.-c(i, j)*rdx(i-1, j))*&
&             dm_tl(i-1) - rdx(i-1, j)*c_tl(i, j)*dm(i-1)
            flux(i, j) = u(i-1, j) + (1.-c(i, j)*rdx(i-1, j))*dm(i-1)
          ELSE
            flux_tl(i, j) = u_tl(i, j) - rdx(i, j)*c_tl(i, j)*dm(i) - (&
&             1.+c(i, j)*rdx(i, j))*dm_tl(i)
            flux(i, j) = u(i, j) - (1.+c(i, j)*rdx(i, j))*dm(i)
          END IF
        END DO
      END DO
    CASE (4) 
      dm_tl = 0.0_8
      al_tl = 0.0_8
      DO j=js,je+1
        DO i=is-2,ie+2
          xt_tl = 0.25*(u_tl(i+1, j)-u_tl(i-1, j))
          xt = 0.25*(u(i+1, j)-u(i-1, j))
          IF (xt .GE. 0.) THEN
            x7_tl = xt_tl
            x7 = xt
          ELSE
            x7_tl = -xt_tl
            x7 = -xt
          END IF
          IF (u(i-1, j) .LT. u(i, j)) THEN
            IF (u(i, j) .LT. u(i+1, j)) THEN
              max12_tl = u_tl(i+1, j)
              max12 = u(i+1, j)
            ELSE
              max12_tl = u_tl(i, j)
              max12 = u(i, j)
            END IF
          ELSE IF (u(i-1, j) .LT. u(i+1, j)) THEN
            max12_tl = u_tl(i+1, j)
            max12 = u(i+1, j)
          ELSE
            max12_tl = u_tl(i-1, j)
            max12 = u(i-1, j)
          END IF
          y6_tl = max12_tl - u_tl(i, j)
          y6 = max12 - u(i, j)
          IF (u(i-1, j) .GT. u(i, j)) THEN
            IF (u(i, j) .GT. u(i+1, j)) THEN
              min19_tl = u_tl(i+1, j)
              min19 = u(i+1, j)
            ELSE
              min19_tl = u_tl(i, j)
              min19 = u(i, j)
            END IF
          ELSE IF (u(i-1, j) .GT. u(i+1, j)) THEN
            min19_tl = u_tl(i+1, j)
            min19 = u(i+1, j)
          ELSE
            min19_tl = u_tl(i-1, j)
            min19 = u(i-1, j)
          END IF
          z6_tl = u_tl(i, j) - min19_tl
          z6 = u(i, j) - min19
          IF (x7 .GT. y6) THEN
            IF (y6 .GT. z6) THEN
              min6_tl = z6_tl
              min6 = z6
            ELSE
              min6_tl = y6_tl
              min6 = y6
            END IF
          ELSE IF (x7 .GT. z6) THEN
            min6_tl = z6_tl
            min6 = z6
          ELSE
            min6_tl = x7_tl
            min6 = x7
          END IF
          dm_tl(i) = min6_tl*SIGN(1.d0, min6*xt)
          dm(i) = SIGN(min6, xt)
        END DO
        IF (3 .LT. is - 1) THEN
          max1 = is - 1
        ELSE
          max1 = 3
        END IF
        IF (npx - 2 .GT. ie + 2) THEN
          min20 = ie + 2
        ELSE
          min20 = npx - 2
        END IF
        DO i=max1,min20
          al_tl(i) = 0.5*(u_tl(i-1, j)+u_tl(i, j)) + r3*(dm_tl(i-1)-&
&           dm_tl(i))
          al(i) = 0.5*(u(i-1, j)+u(i, j)) + r3*(dm(i-1)-dm(i))
        END DO
! Fix slopes near edges:
        IF (grid_type .LT. 3) THEN
          IF (is .EQ. 1) THEN
            IF (j .EQ. 1 .OR. j .EQ. npy) THEN
              dm_tl(0) = 0.0_8
              dm(0) = 0.
              dm_tl(1) = 0.0_8
              dm(1) = 0.
              al_tl(0) = 0.5*(u_tl(-1, j)+u_tl(0, j)) + r3*dm_tl(-1)
              al(0) = 0.5*(u(-1, j)+u(0, j)) + r3*dm(-1)
              al_tl(1) = 0.5*((2.*dx(1, j)+dx(2, j))*(u_tl(0, j)+u_tl(1&
&               , j))-dx(1, j)*(u_tl(-1, j)+u_tl(2, j)))/(dx(1, j)+dx(2&
&               , j))
              al(1) = 0.5*((2.*dx(1, j)+dx(2, j))*(u(0, j)+u(1, j))-dx(1&
&               , j)*(u(-1, j)+u(2, j)))/(dx(1, j)+dx(2, j))
              al_tl(2) = 0.5*(u_tl(1, j)+u_tl(2, j)) - r3*dm_tl(2)
              al(2) = 0.5*(u(1, j)+u(2, j)) - r3*dm(2)
            ELSE
              x0_tl = 0.5*((2.*dx(1, j)+dx(2, j))*(u_tl(0, j)+u_tl(1, j)&
&               )-dx(1, j)*(u_tl(-1, j)+u_tl(2, j)))/(dx(1, j)+dx(2, j))
              x0 = 0.5*((2.*dx(1, j)+dx(2, j))*(u(0, j)+u(1, j))-dx(1, j&
&               )*(u(-1, j)+u(2, j)))/(dx(1, j)+dx(2, j))
              x1_tl = s15*u_tl(1, j) + s11*u_tl(2, j) - s14*dm_tl(2)
              x1 = s15*u(1, j) + s11*u(2, j) - s14*dm(2)
              dm_tl(1) = 0.5*(x1_tl-x0_tl)
              dm(1) = 0.5*(x1-x0)
!          dm(1) = sign(min(abs(dm(1)), max(u(1,j), x0, x1) - u(1,j),   &
!                              u(1,j) - min(u(1,j), x0, x1)), dm(1))
              x1_tl = s15*u_tl(0, j) + s11*u_tl(-1, j) + s14*dm_tl(-1)
              x1 = s15*u(0, j) + s11*u(-1, j) + s14*dm(-1)
              dm_tl(0) = 0.5*(x0_tl-x1_tl)
              dm(0) = 0.5*(x0-x1)
!          dm(0) = sign(min(abs(dm(0)), max(u(0,j), x0, x1) - u(0,j),   &
!                              u(0,j) - min(u(0,j), x0, x1)), dm(0))
              al_tl(0) = 0.5*(u_tl(-1, j)+u_tl(0, j)) + r3*(dm_tl(-1)-&
&               dm_tl(0))
              al(0) = 0.5*(u(-1, j)+u(0, j)) + r3*(dm(-1)-dm(0))
              al_tl(1) = x0_tl
              al(1) = x0
              al_tl(2) = 0.5*(u_tl(1, j)+u_tl(2, j)) + r3*(dm_tl(1)-&
&               dm_tl(2))
              al(2) = 0.5*(u(1, j)+u(2, j)) + r3*(dm(1)-dm(2))
            END IF
          END IF
          IF (ie + 1 .EQ. npx) THEN
            IF (j .EQ. 1 .OR. j .EQ. npy) THEN
              dm_tl(npx-1) = 0.0_8
              dm(npx-1) = 0.
              dm_tl(npx) = 0.0_8
              dm(npx) = 0.
              al_tl(npx-1) = 0.5*(u_tl(npx-2, j)+u_tl(npx-1, j)) + r3*&
&               dm_tl(npx-2)
              al(npx-1) = 0.5*(u(npx-2, j)+u(npx-1, j)) + r3*dm(npx-2)
              al_tl(npx) = 0.5*((2.*dx(npx-1, j)+dx(npx-2, j))*(u_tl(npx&
&               -1, j)+u_tl(npx, j))-dx(npx-1, j)*(u_tl(npx-2, j)+u_tl(&
&               npx+1, j)))/(dx(npx-1, j)+dx(npx-2, j))
              al(npx) = 0.5*((2.*dx(npx-1, j)+dx(npx-2, j))*(u(npx-1, j)&
&               +u(npx, j))-dx(npx-1, j)*(u(npx-2, j)+u(npx+1, j)))/(dx(&
&               npx-1, j)+dx(npx-2, j))
              al_tl(npx+1) = 0.5*(u_tl(npx, j)+u_tl(npx+1, j)) - r3*&
&               dm_tl(npx+1)
              al(npx+1) = 0.5*(u(npx, j)+u(npx+1, j)) - r3*dm(npx+1)
            ELSE
              x0_tl = 0.5*((2.*dx(npx-1, j)+dx(npx-2, j))*(u_tl(npx-1, j&
&               )+u_tl(npx, j))-dx(npx-1, j)*(u_tl(npx-2, j)+u_tl(npx+1&
&               , j)))/(dx(npx-1, j)+dx(npx-2, j))
              x0 = 0.5*((2.*dx(npx-1, j)+dx(npx-2, j))*(u(npx-1, j)+u(&
&               npx, j))-dx(npx-1, j)*(u(npx-2, j)+u(npx+1, j)))/(dx(npx&
&               -1, j)+dx(npx-2, j))
              x1_tl = s15*u_tl(npx-1, j) + s11*u_tl(npx-2, j) + s14*&
&               dm_tl(npx-2)
              x1 = s15*u(npx-1, j) + s11*u(npx-2, j) + s14*dm(npx-2)
              dm_tl(npx-1) = 0.5*(x0_tl-x1_tl)
              dm(npx-1) = 0.5*(x0-x1)
!          dm(npx-1) = sign(min(abs(dm(npx-1)), max(u(npx-1,j), x0, x1) - u(npx-1,j),  &
!                                  u(npx-1,j) - min(u(npx-1,j), x0, x1)), dm(npx-1))
              x1_tl = s15*u_tl(npx, j) + s11*u_tl(npx+1, j) - s14*dm_tl(&
&               npx+1)
              x1 = s15*u(npx, j) + s11*u(npx+1, j) - s14*dm(npx+1)
              dm_tl(npx) = 0.5*(x1_tl-x0_tl)
              dm(npx) = 0.5*(x1-x0)
!          dm(npx) = sign(min(abs(dm(npx)), max(u(npx,j), x0, x1) - u(npx,j),   &
!                                u(npx,j) - min(u(npx,j), x0, x1)), dm(npx))
              al_tl(npx-1) = 0.5*(u_tl(npx-2, j)+u_tl(npx-1, j)) + r3*(&
&               dm_tl(npx-2)-dm_tl(npx-1))
              al(npx-1) = 0.5*(u(npx-2, j)+u(npx-1, j)) + r3*(dm(npx-2)-&
&               dm(npx-1))
              al_tl(npx) = x0_tl
              al(npx) = x0
              al_tl(npx+1) = 0.5*(u_tl(npx, j)+u_tl(npx+1, j)) + r3*(&
&               dm_tl(npx)-dm_tl(npx+1))
              al(npx+1) = 0.5*(u(npx, j)+u(npx+1, j)) + r3*(dm(npx)-dm(&
&               npx+1))
            END IF
          END IF
        END IF
        DO i=is,ie+1
          IF (c(i, j) .GT. 0.) THEN
            xt_tl = 2.*dm_tl(i-1)
            xt = 2.*dm(i-1)
            IF (xt .GE. 0.) THEN
              x8_tl = xt_tl
              x8 = xt
            ELSE
              x8_tl = -xt_tl
              x8 = -xt
            END IF
            IF (al(i-1) - u(i-1, j) .GE. 0.) THEN
              y7_tl = al_tl(i-1) - u_tl(i-1, j)
              y7 = al(i-1) - u(i-1, j)
            ELSE
              y7_tl = -(al_tl(i-1)-u_tl(i-1, j))
              y7 = -(al(i-1)-u(i-1, j))
            END IF
            IF (x8 .GT. y7) THEN
              min7_tl = y7_tl
              min7 = y7
            ELSE
              min7_tl = x8_tl
              min7 = x8
            END IF
            dl_tl = min7_tl*SIGN(1.d0, min7*xt)
            dl = SIGN(min7, xt)
            IF (xt .GE. 0.) THEN
              x9_tl = xt_tl
              x9 = xt
            ELSE
              x9_tl = -xt_tl
              x9 = -xt
            END IF
            IF (al(i) - u(i-1, j) .GE. 0.) THEN
              y8_tl = al_tl(i) - u_tl(i-1, j)
              y8 = al(i) - u(i-1, j)
            ELSE
              y8_tl = -(al_tl(i)-u_tl(i-1, j))
              y8 = -(al(i)-u(i-1, j))
            END IF
            IF (x9 .GT. y8) THEN
              min8_tl = y8_tl
              min8 = y8
            ELSE
              min8_tl = x9_tl
              min8 = x9
            END IF
            dr_tl = min8_tl*SIGN(1.d0, min8*xt)
            dr = SIGN(min8, xt)
            cfl_tl = rdx(i-1, j)*c_tl(i, j)
            cfl = c(i, j)*rdx(i-1, j)
            flux_tl(i, j) = u_tl(i-1, j) + (1.-cfl)*(dr_tl+cfl_tl*(dl-dr&
&             )+cfl*(dl_tl-dr_tl)) - cfl_tl*(dr+cfl*(dl-dr))
            flux(i, j) = u(i-1, j) + (1.-cfl)*(dr+cfl*(dl-dr))
          ELSE
            xt_tl = 2.*dm_tl(i)
            xt = 2.*dm(i)
            IF (xt .GE. 0.) THEN
              x10_tl = xt_tl
              x10 = xt
            ELSE
              x10_tl = -xt_tl
              x10 = -xt
            END IF
            IF (al(i) - u(i, j) .GE. 0.) THEN
              y9_tl = al_tl(i) - u_tl(i, j)
              y9 = al(i) - u(i, j)
            ELSE
              y9_tl = -(al_tl(i)-u_tl(i, j))
              y9 = -(al(i)-u(i, j))
            END IF
            IF (x10 .GT. y9) THEN
              min9_tl = y9_tl
              min9 = y9
            ELSE
              min9_tl = x10_tl
              min9 = x10
            END IF
            dl_tl = min9_tl*SIGN(1.d0, min9*xt)
            dl = SIGN(min9, xt)
            IF (xt .GE. 0.) THEN
              x11_tl = xt_tl
              x11 = xt
            ELSE
              x11_tl = -xt_tl
              x11 = -xt
            END IF
            IF (al(i+1) - u(i, j) .GE. 0.) THEN
              y10_tl = al_tl(i+1) - u_tl(i, j)
              y10 = al(i+1) - u(i, j)
            ELSE
              y10_tl = -(al_tl(i+1)-u_tl(i, j))
              y10 = -(al(i+1)-u(i, j))
            END IF
            IF (x11 .GT. y10) THEN
              min10_tl = y10_tl
              min10 = y10
            ELSE
              min10_tl = x11_tl
              min10 = x11
            END IF
            dr_tl = min10_tl*SIGN(1.d0, min10*xt)
            dr = SIGN(min10, xt)
            cfl_tl = rdx(i, j)*c_tl(i, j)
            cfl = c(i, j)*rdx(i, j)
            flux_tl(i, j) = u_tl(i, j) - cfl_tl*(dl+cfl*(dl-dr)) - (1.+&
&             cfl)*(dl_tl+cfl_tl*(dl-dr)+cfl*(dl_tl-dr_tl))
            flux(i, j) = u(i, j) - (1.+cfl)*(dl+cfl*(dl-dr))
          END IF
        END DO
      END DO
    CASE (6) 
      bl_tl = 0.0_8
      br_tl = 0.0_8
      DO j=js,je+1
        IF (3 .LT. is - 1) THEN
          max2 = is - 1
        ELSE
          max2 = 3
        END IF
        IF (npx - 3 .GT. ie + 1) THEN
          min21 = ie + 1
        ELSE
          min21 = npx - 3
        END IF
        DO i=max2,min21
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
    CASE DEFAULT
      dm_tl = 0.0_8
      dq_tl = 0.0_8
      al_tl = 0.0_8
      bl_tl = 0.0_8
      br_tl = 0.0_8
! iord = 8, 9, 10
      DO j=js,je+1
        DO i=is-2,ie+2
          xt_tl = 0.25*(u_tl(i+1, j)-u_tl(i-1, j))
          xt = 0.25*(u(i+1, j)-u(i-1, j))
          IF (xt .GE. 0.) THEN
            x12_tl = xt_tl
            x12 = xt
          ELSE
            x12_tl = -xt_tl
            x12 = -xt
          END IF
          IF (u(i-1, j) .LT. u(i, j)) THEN
            IF (u(i, j) .LT. u(i+1, j)) THEN
              max13_tl = u_tl(i+1, j)
              max13 = u(i+1, j)
            ELSE
              max13_tl = u_tl(i, j)
              max13 = u(i, j)
            END IF
          ELSE IF (u(i-1, j) .LT. u(i+1, j)) THEN
            max13_tl = u_tl(i+1, j)
            max13 = u(i+1, j)
          ELSE
            max13_tl = u_tl(i-1, j)
            max13 = u(i-1, j)
          END IF
          y11_tl = max13_tl - u_tl(i, j)
          y11 = max13 - u(i, j)
          IF (u(i-1, j) .GT. u(i, j)) THEN
            IF (u(i, j) .GT. u(i+1, j)) THEN
              min22_tl = u_tl(i+1, j)
              min22 = u(i+1, j)
            ELSE
              min22_tl = u_tl(i, j)
              min22 = u(i, j)
            END IF
          ELSE IF (u(i-1, j) .GT. u(i+1, j)) THEN
            min22_tl = u_tl(i+1, j)
            min22 = u(i+1, j)
          ELSE
            min22_tl = u_tl(i-1, j)
            min22 = u(i-1, j)
          END IF
          z7_tl = u_tl(i, j) - min22_tl
          z7 = u(i, j) - min22
          IF (x12 .GT. y11) THEN
            IF (y11 .GT. z7) THEN
              min11_tl = z7_tl
              min11 = z7
            ELSE
              min11_tl = y11_tl
              min11 = y11
            END IF
          ELSE IF (x12 .GT. z7) THEN
            min11_tl = z7_tl
            min11 = z7
          ELSE
            min11_tl = x12_tl
            min11 = x12
          END IF
          dm_tl(i) = min11_tl*SIGN(1.d0, min11*xt)
          dm(i) = SIGN(min11, xt)
        END DO
        DO i=is-3,ie+2
          dq_tl(i) = u_tl(i+1, j) - u_tl(i, j)
          dq(i) = u(i+1, j) - u(i, j)
        END DO
        IF (grid_type .LT. 3) THEN
          IF (3 .LT. is - 1) THEN
            max3 = is - 1
          ELSE
            max3 = 3
          END IF
          IF (npx - 2 .GT. ie + 2) THEN
            min23 = ie + 2
          ELSE
            min23 = npx - 2
          END IF
          DO i=max3,min23
            al_tl(i) = 0.5*(u_tl(i-1, j)+u_tl(i, j)) + r3*(dm_tl(i-1)-&
&             dm_tl(i))
            al(i) = 0.5*(u(i-1, j)+u(i, j)) + r3*(dm(i-1)-dm(i))
          END DO
! Perturbation form:
          IF (iord .EQ. 8) THEN
            IF (3 .LT. is - 1) THEN
              max4 = is - 1
            ELSE
              max4 = 3
            END IF
            IF (npx - 3 .GT. ie + 1) THEN
              min24 = ie + 1
            ELSE
              min24 = npx - 3
            END IF
            DO i=max4,min24
              xt_tl = 2.*dm_tl(i)
              xt = 2.*dm(i)
              IF (xt .GE. 0.) THEN
                x13_tl = xt_tl
                x13 = xt
              ELSE
                x13_tl = -xt_tl
                x13 = -xt
              END IF
              IF (al(i) - u(i, j) .GE. 0.) THEN
                y12_tl = al_tl(i) - u_tl(i, j)
                y12 = al(i) - u(i, j)
              ELSE
                y12_tl = -(al_tl(i)-u_tl(i, j))
                y12 = -(al(i)-u(i, j))
              END IF
              IF (x13 .GT. y12) THEN
                min12_tl = y12_tl
                min12 = y12
              ELSE
                min12_tl = x13_tl
                min12 = x13
              END IF
              bl_tl(i) = -(min12_tl*SIGN(1.d0, min12*xt))
              bl(i) = -SIGN(min12, xt)
              IF (xt .GE. 0.) THEN
                x14_tl = xt_tl
                x14 = xt
              ELSE
                x14_tl = -xt_tl
                x14 = -xt
              END IF
              IF (al(i+1) - u(i, j) .GE. 0.) THEN
                y13_tl = al_tl(i+1) - u_tl(i, j)
                y13 = al(i+1) - u(i, j)
              ELSE
                y13_tl = -(al_tl(i+1)-u_tl(i, j))
                y13 = -(al(i+1)-u(i, j))
              END IF
              IF (x14 .GT. y13) THEN
                min13_tl = y13_tl
                min13 = y13
              ELSE
                min13_tl = x14_tl
                min13 = x14
              END IF
              br_tl(i) = min13_tl*SIGN(1.d0, min13*xt)
              br(i) = SIGN(min13, xt)
            END DO
          ELSE IF (iord .EQ. 9) THEN
            IF (3 .LT. is - 1) THEN
              max5 = is - 1
            ELSE
              max5 = 3
            END IF
            IF (npx - 3 .GT. ie + 1) THEN
              min25 = ie + 1
            ELSE
              min25 = npx - 3
            END IF
            DO i=max5,min25
              pmp_tl = 2.*dq_tl(i-1)
              pmp = 2.*dq(i-1)
              lac_tl = pmp_tl - 1.5*dq_tl(i-2)
              lac = pmp - 1.5*dq(i-2)
              IF (0. .LT. pmp) THEN
                IF (pmp .LT. lac) THEN
                  x15_tl = lac_tl
                  x15 = lac
                ELSE
                  x15_tl = pmp_tl
                  x15 = pmp
                END IF
              ELSE IF (0. .LT. lac) THEN
                x15_tl = lac_tl
                x15 = lac
              ELSE
                x15 = 0.
                x15_tl = 0.0_8
              END IF
              IF (0. .GT. pmp) THEN
                IF (pmp .GT. lac) THEN
                  y18_tl = lac_tl
                  y18 = lac
                ELSE
                  y18_tl = pmp_tl
                  y18 = pmp
                END IF
              ELSE IF (0. .GT. lac) THEN
                y18_tl = lac_tl
                y18 = lac
              ELSE
                y18 = 0.
                y18_tl = 0.0_8
              END IF
              IF (al(i+1) - u(i, j) .LT. y18) THEN
                y14_tl = y18_tl
                y14 = y18
              ELSE
                y14_tl = al_tl(i+1) - u_tl(i, j)
                y14 = al(i+1) - u(i, j)
              END IF
              IF (x15 .GT. y14) THEN
                br_tl(i) = y14_tl
                br(i) = y14
              ELSE
                br_tl(i) = x15_tl
                br(i) = x15
              END IF
              pmp_tl = -(2.*dq_tl(i))
              pmp = -(2.*dq(i))
              lac_tl = pmp_tl + 1.5*dq_tl(i+1)
              lac = pmp + 1.5*dq(i+1)
              IF (0. .LT. pmp) THEN
                IF (pmp .LT. lac) THEN
                  x16_tl = lac_tl
                  x16 = lac
                ELSE
                  x16_tl = pmp_tl
                  x16 = pmp
                END IF
              ELSE IF (0. .LT. lac) THEN
                x16_tl = lac_tl
                x16 = lac
              ELSE
                x16 = 0.
                x16_tl = 0.0_8
              END IF
              IF (0. .GT. pmp) THEN
                IF (pmp .GT. lac) THEN
                  y19_tl = lac_tl
                  y19 = lac
                ELSE
                  y19_tl = pmp_tl
                  y19 = pmp
                END IF
              ELSE IF (0. .GT. lac) THEN
                y19_tl = lac_tl
                y19 = lac
              ELSE
                y19 = 0.
                y19_tl = 0.0_8
              END IF
              IF (al(i) - u(i, j) .LT. y19) THEN
                y15_tl = y19_tl
                y15 = y19
              ELSE
                y15_tl = al_tl(i) - u_tl(i, j)
                y15 = al(i) - u(i, j)
              END IF
              IF (x16 .GT. y15) THEN
                bl_tl(i) = y15_tl
                bl(i) = y15
              ELSE
                bl_tl(i) = x16_tl
                bl(i) = x16
              END IF
            END DO
          ELSE
            IF (3 .LT. is - 1) THEN
              max6 = is - 1
            ELSE
              max6 = 3
            END IF
            IF (npx - 3 .GT. ie + 1) THEN
              min26 = ie + 1
            ELSE
              min26 = npx - 3
            END IF
! un-limited:
            DO i=max6,min26
              bl_tl(i) = al_tl(i) - u_tl(i, j)
              bl(i) = al(i) - u(i, j)
              br_tl(i) = al_tl(i+1) - u_tl(i, j)
              br(i) = al(i+1) - u(i, j)
            END DO
          END IF
!--------------
! fix the edges
!--------------
          IF (is .EQ. 1) THEN
            br_tl(2) = al_tl(3) - u_tl(2, j)
            br(2) = al(3) - u(2, j)
            xt_tl = s15*u_tl(1, j) + s11*u_tl(2, j) - s14*dm_tl(2)
            xt = s15*u(1, j) + s11*u(2, j) - s14*dm(2)
            bl_tl(2) = xt_tl - u_tl(2, j)
            bl(2) = xt - u(2, j)
            br_tl(1) = xt_tl - u_tl(1, j)
            br(1) = xt - u(1, j)
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
              bl_tl(0) = s14*dm_tl(-1) - s11*dq_tl(-1)
              bl(0) = s14*dm(-1) - s11*dq(-1)
!---------------------------------------------------------------
              xt_tl = 0.5*((2.*dx(1, j)+dx(2, j))*(u_tl(0, j)+u_tl(1, j)&
&               )-dx(1, j)*(u_tl(-1, j)+u_tl(2, j)))/(dx(1, j)+dx(2, j))
              xt = 0.5*((2.*dx(1, j)+dx(2, j))*(u(0, j)+u(1, j))-dx(1, j&
&               )*(u(-1, j)+u(2, j)))/(dx(1, j)+dx(2, j))
              br_tl(0) = xt_tl - u_tl(0, j)
              br(0) = xt - u(0, j)
              bl_tl(1) = xt_tl - u_tl(1, j)
              bl(1) = xt - u(1, j)
!                br(0) = xt - 0.5*(v(1,j-1)+v(1,j))*cosa(1,j) - u(0,j)
!                bl(1) = xt + 0.5*(v(1,j-1)+v(1,j))*cosa(1,j) - u(1,j)
!---------------------------------------------------------------
            END IF
            IF (iord .LT. 10) CALL PERT_PPM_TLM(1, u(2:, j), bl(2:), &
&                                         bl_tl(2:), br(2:), br_tl(2:), &
&                                         -1)
          END IF
          IF (ie + 1 .EQ. npx) THEN
            bl_tl(npx-2) = al_tl(npx-2) - u_tl(npx-2, j)
            bl(npx-2) = al(npx-2) - u(npx-2, j)
            xt_tl = s15*u_tl(npx-1, j) + s11*u_tl(npx-2, j) + s14*dm_tl(&
&             npx-2)
            xt = s15*u(npx-1, j) + s11*u(npx-2, j) + s14*dm(npx-2)
            br_tl(npx-2) = xt_tl - u_tl(npx-2, j)
            br(npx-2) = xt - u(npx-2, j)
            bl_tl(npx-1) = xt_tl - u_tl(npx-1, j)
            bl(npx-1) = xt - u(npx-1, j)
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
              br_tl(npx) = s11*dq_tl(npx) - s14*dm_tl(npx+1)
              br(npx) = s11*dq(npx) - s14*dm(npx+1)
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
!               br(npx-1) = xt + 0.5*(v(npx,j-1)+v(npx,j))*cosa(npx,j) - u(npx-1,j)
!               bl(npx  ) = xt - 0.5*(v(npx,j-1)+v(npx,j))*cosa(npx,j) - u(npx  ,j)
            END IF
            IF (iord .LT. 10) CALL PERT_PPM_TLM(1, u(npx-2:, j), bl(npx-&
&                                         2:), bl_tl(npx-2:), br(npx-2:)&
&                                         , br_tl(npx-2:), -1)
          END IF
        ELSE
          DO i=is-1,ie+2
            al_tl(i) = 0.5*(u_tl(i-1, j)+u_tl(i, j)) + r3*(dm_tl(i-1)-&
&             dm_tl(i))
            al(i) = 0.5*(u(i-1, j)+u(i, j)) + r3*(dm(i-1)-dm(i))
          END DO
          DO i=is-1,ie+1
            pmp_tl = -(2.*dq_tl(i))
            pmp = -(2.*dq(i))
            lac_tl = pmp_tl + 1.5*dq_tl(i+1)
            lac = pmp + 1.5*dq(i+1)
            IF (0. .LT. pmp) THEN
              IF (pmp .LT. lac) THEN
                x17_tl = lac_tl
                x17 = lac
              ELSE
                x17_tl = pmp_tl
                x17 = pmp
              END IF
            ELSE IF (0. .LT. lac) THEN
              x17_tl = lac_tl
              x17 = lac
            ELSE
              x17 = 0.
              x17_tl = 0.0_8
            END IF
            IF (0. .GT. pmp) THEN
              IF (pmp .GT. lac) THEN
                y20_tl = lac_tl
                y20 = lac
              ELSE
                y20_tl = pmp_tl
                y20 = pmp
              END IF
            ELSE IF (0. .GT. lac) THEN
              y20_tl = lac_tl
              y20 = lac
            ELSE
              y20 = 0.
              y20_tl = 0.0_8
            END IF
            IF (al(i) - u(i, j) .LT. y20) THEN
              y16_tl = y20_tl
              y16 = y20
            ELSE
              y16_tl = al_tl(i) - u_tl(i, j)
              y16 = al(i) - u(i, j)
            END IF
            IF (x17 .GT. y16) THEN
              bl_tl(i) = y16_tl
              bl(i) = y16
            ELSE
              bl_tl(i) = x17_tl
              bl(i) = x17
            END IF
            pmp_tl = 2.*dq_tl(i-1)
            pmp = 2.*dq(i-1)
            lac_tl = pmp_tl - 1.5*dq_tl(i-2)
            lac = pmp - 1.5*dq(i-2)
            IF (0. .LT. pmp) THEN
              IF (pmp .LT. lac) THEN
                x18_tl = lac_tl
                x18 = lac
              ELSE
                x18_tl = pmp_tl
                x18 = pmp
              END IF
            ELSE IF (0. .LT. lac) THEN
              x18_tl = lac_tl
              x18 = lac
            ELSE
              x18 = 0.
              x18_tl = 0.0_8
            END IF
            IF (0. .GT. pmp) THEN
              IF (pmp .GT. lac) THEN
                y21_tl = lac_tl
                y21 = lac
              ELSE
                y21_tl = pmp_tl
                y21 = pmp
              END IF
            ELSE IF (0. .GT. lac) THEN
              y21_tl = lac_tl
              y21 = lac
            ELSE
              y21 = 0.
              y21_tl = 0.0_8
            END IF
            IF (al(i+1) - u(i, j) .LT. y21) THEN
              y17_tl = y21_tl
              y17 = y21
            ELSE
              y17_tl = al_tl(i+1) - u_tl(i, j)
              y17 = al(i+1) - u(i, j)
            END IF
            IF (x18 .GT. y17) THEN
              br_tl(i) = y17_tl
              br(i) = y17
            ELSE
              br_tl(i) = x18_tl
              br(i) = x18
            END IF
          END DO
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
    REAL(ke_precision) :: dm_tl(is:ie+1, js-2:je+2)
    REAL(ke_precision) :: al(is:ie+1, js-1:je+2)
    REAL(ke_precision) :: al_tl(is:ie+1, js-1:je+2)
    REAL(ke_precision) :: bl(is:ie+1, js-1:je+1)
    REAL(ke_precision) :: bl_tl(is:ie+1, js-1:je+1)
    REAL(ke_precision) :: br(is:ie+1, js-1:je+1)
    REAL(ke_precision) :: br_tl(is:ie+1, js-1:je+1)
    REAL(ke_precision) :: dq(is:ie+1, js-3:je+2)
    REAL(ke_precision) :: dq_tl(is:ie+1, js-3:je+2)
    REAL(ke_precision) :: xt, dl, dr, pmp, lac, dqt, cfl
    REAL(ke_precision) :: xt_tl, dl_tl, dr_tl, pmp_tl, lac_tl, cfl_tl
    REAL(ke_precision) :: x0, x1
    REAL(ke_precision) :: x0_tl, x1_tl
    INTEGER :: i, j
    INTRINSIC ABS
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC SIGN
    REAL*8 :: y1_tl
    REAL*8 :: min6_tl
    REAL*8 :: y18_tl
    REAL*8 :: x18
    REAL*8 :: x17
    REAL*8 :: min10_tl
    REAL*8 :: x16
    REAL*8 :: z4_tl
    REAL*8 :: x13_tl
    REAL*8 :: x15
    REAL*8 :: max7_tl
    REAL*8 :: y20_tl
    REAL*8 :: x14
    REAL*8 :: min9
    REAL*8 :: x13
    REAL*8 :: min8
    REAL*8 :: y2_tl
    REAL*8 :: min7_tl
    REAL*8 :: y19_tl
    REAL*8 :: x12
    REAL*8 :: min7
    REAL*8 :: x11
    REAL*8 :: min6
    REAL*8 :: x10
    REAL*8 :: min5
    REAL*8 :: min4
    REAL*8 :: min3
    REAL*8 :: min11_tl
    REAL*8 :: min2
    REAL*8 :: z5_tl
    REAL*8 :: x14_tl
    REAL*8 :: min1
    REAL*8 :: max8_tl
    REAL*8 :: y21_tl
    REAL*8 :: y3_tl
    REAL*8 :: min8_tl
    REAL*8 :: y21
    REAL*8 :: y20
    REAL*8 :: x9
    REAL*8 :: x8
    REAL*8 :: min12_tl
    REAL*8 :: x7
    REAL*8 :: z6_tl
    REAL*8 :: x15_tl
    REAL*8 :: x6
    REAL*8 :: max9_tl
    REAL*8 :: x5
    REAL*8 :: x4
    REAL*8 :: y4_tl
    REAL*8 :: min9_tl
    REAL*8 :: x3
    REAL*8 :: y10_tl
    REAL*8 :: x2
    REAL*8 :: x2_tl
    REAL*8 :: min13_tl
    REAL*8 :: z7_tl
    REAL*8 :: x16_tl
    REAL*8 :: y5_tl
    REAL*8 :: y11_tl
    REAL*8 :: x3_tl
    REAL :: min14_tl
    REAL*8 :: y19
    INTEGER :: min25
    REAL*8 :: x17_tl
    REAL*8 :: y18
    INTEGER :: min24
    REAL*8 :: y17
    INTEGER :: min23
    REAL*8 :: y16
    INTEGER :: min22
    REAL*8 :: y6_tl
    REAL*8 :: y15
    REAL :: min21
    REAL*8 :: y12_tl
    REAL*8 :: y14
    INTEGER :: min20
    REAL*8 :: y13
    REAL*8 :: x4_tl
    REAL*8 :: y12
    REAL*8 :: y11
    REAL*8 :: min15_tl
    REAL*8 :: y10
    REAL*8 :: max10_tl
    REAL*8 :: x18_tl
    REAL*8 :: y7_tl
    REAL*8 :: min1_tl
    REAL*8 :: y13_tl
    REAL*8 :: x5_tl
    REAL*8 :: min16_tl
    REAL :: max11_tl
    REAL*8 :: z7
    REAL*8 :: z6
    REAL*8 :: y8_tl
    REAL*8 :: z5
    REAL*8 :: min2_tl
    REAL*8 :: y14_tl
    REAL*8 :: z4
    REAL*8 :: z3
    REAL*8 :: x6_tl
    REAL*8 :: z2
    REAL*8 :: z1
    REAL*8 :: min17_tl
    REAL :: min19
    REAL :: max12_tl
    REAL*8 :: min18
    REAL*8 :: min17
    REAL*8 :: min16
    REAL*8 :: y9_tl
    REAL*8 :: min15
    REAL*8 :: min3_tl
    REAL*8 :: y15_tl
    REAL :: min14
    REAL*8 :: min13
    REAL*8 :: x7_tl
    REAL*8 :: min12
    REAL*8 :: min11
    REAL*8 :: min18_tl
    REAL*8 :: min10
    REAL*8 :: z1_tl
    REAL*8 :: x10_tl
    REAL*8 :: min4_tl
    REAL*8 :: y16_tl
    REAL*8 :: x8_tl
    REAL*8 :: max9
    REAL :: min19_tl
    REAL*8 :: max8
    REAL*8 :: max7
    REAL*8 :: z2_tl
    REAL*8 :: x11_tl
    REAL :: max6
    INTEGER :: max5
    REAL*8 :: y9
    INTEGER :: max4
    REAL*8 :: min5_tl
    REAL*8 :: y17_tl
    REAL*8 :: y8
    INTEGER :: max3
    REAL*8 :: y7
    INTEGER :: max2
    REAL*8 :: x9_tl
    REAL :: min21_tl
    REAL*8 :: y6
    INTEGER :: max1
    REAL :: max12
    REAL*8 :: y5
    REAL :: max11
    REAL*8 :: y4
    REAL*8 :: max10
    REAL*8 :: y3
    REAL*8 :: z3_tl
    REAL*8 :: x12_tl
    REAL*8 :: y2
    REAL :: max6_tl
    REAL*8 :: y1
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

      !First order
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

      !Upgrade to third order scheme internally
      DO j=js+1,je
        DO i=is,ie+1
          IF (c(i, j) .GT. 0.) THEN
            cfl = c(i, j)*rdy(i, j-1)
            cfl_tl = c_tl(i, j)*rdy(i, j-1)
            flux_tl(i, j) = (2.0*v_tl(i, j)+5.0*v_tl(i, j-1)-v_tl(i, j-2&
&             ))/6.0 - 0.5*(cfl_tl*(v(i, j)-v(i, j-1))+cfl*(v_tl&
&             (i, j)-v_tl(i, j-1))) + (cfl_tl*cfl+cfl*cfl_tl&
&             )*(v(i, j)-2.0*v(i, j-1)+v(i, j-2))/6.0 + cfl**2*(&
&             v_tl(i, j)-2.0*v_tl(i, j-1)+v_tl(i, j-2))/6.0
            flux(i, j) = (2.0*v(i, j)+5.0*v(i, j-1)-v(i, j-2))/6.0 - 0.5&
&             *cfl*(v(i, j)-v(i, j-1)) + cfl*cfl/6.0*(v(i, j&
&             )-2.0*v(i, j-1)+v(i, j-2))
          ELSE
            cfl = c(i, j)*rdy(i, j)
            cfl_tl = c_tl(i, j)*rdy(i, j)
            flux_tl(i, j) = (2.0*v_tl(i, j-1)+5.0*v_tl(i, j)-v_tl(i, j+1&
&             ))/6.0 - 0.5*(cfl_tl*(v(i, j)-v(i, j-1))+cfl*(v_tl&
&             (i, j)-v_tl(i, j-1))) + (cfl_tl*cfl+cfl*cfl_tl&
&             )*(v(i, j+1)-2.0*v(i, j)+v(i, j-1))/6.0 + cfl**2*(&
&             v_tl(i, j+1)-2.0*v_tl(i, j)+v_tl(i, j-1))/6.0
            flux(i, j) = (2.0*v(i, j-1)+5.0*v(i, j)-v(i, j+1))/6.0 - 0.5&
&             *cfl*(v(i, j)-v(i, j-1)) + cfl*cfl/6.0*(v(i, j&
&             +1)-2.0*v(i, j)+v(i, j-1))
          END IF
        END DO
      END DO

      !Downgrade to third order for the edges
      IF (grid_type .LT. 3) THEN
        IF (js .EQ. 1) THEN
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
        ENDIF
        IF (je + 1 .EQ. npy) THEN
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
        ENDIF
      ENDIF

    CASE (2) 
      dm_tl = 0.0_8
      DO j=js-2,je+2
        DO i=is,ie+1
          xt_tl = 0.25*(v_tl(i, j+1)-v_tl(i, j-1))
          xt = 0.25*(v(i, j+1)-v(i, j-1))
          IF (xt .GE. 0.) THEN
            x2_tl = xt_tl
            x2 = xt
          ELSE
            x2_tl = -xt_tl
            x2 = -xt
          END IF
          IF (v(i, j-1) .LT. v(i, j)) THEN
            IF (v(i, j) .LT. v(i, j+1)) THEN
              max6_tl = v_tl(i, j+1)
              max6 = v(i, j+1)
            ELSE
              max6_tl = v_tl(i, j)
              max6 = v(i, j)
            END IF
          ELSE IF (v(i, j-1) .LT. v(i, j+1)) THEN
            max6_tl = v_tl(i, j+1)
            max6 = v(i, j+1)
          ELSE
            max6_tl = v_tl(i, j-1)
            max6 = v(i, j-1)
          END IF
          y1_tl = max6_tl - v_tl(i, j)
          y1 = max6 - v(i, j)
          IF (v(i, j-1) .GT. v(i, j)) THEN
            IF (v(i, j) .GT. v(i, j+1)) THEN
              min14_tl = v_tl(i, j+1)
              min14 = v(i, j+1)
            ELSE
              min14_tl = v_tl(i, j)
              min14 = v(i, j)
            END IF
          ELSE IF (v(i, j-1) .GT. v(i, j+1)) THEN
            min14_tl = v_tl(i, j+1)
            min14 = v(i, j+1)
          ELSE
            min14_tl = v_tl(i, j-1)
            min14 = v(i, j-1)
          END IF
          z1_tl = v_tl(i, j) - min14_tl
          z1 = v(i, j) - min14
          IF (x2 .GT. y1) THEN
            IF (y1 .GT. z1) THEN
              min1_tl = z1_tl
              min1 = z1
            ELSE
              min1_tl = y1_tl
              min1 = y1
            END IF
          ELSE IF (x2 .GT. z1) THEN
            min1_tl = z1_tl
            min1 = z1
          ELSE
            min1_tl = x2_tl
            min1 = x2
          END IF
          dm_tl(i, j) = min1_tl*SIGN(1.d0, min1*xt)
          dm(i, j) = SIGN(min1, xt)
        END DO
      END DO
      IF (grid_type .LT. 3) THEN
        IF (js .EQ. 1) THEN
          DO i=is,ie+1
            x0_tl = 0.5*((2.*dy(i, 1)+dy(i, 2))*(v_tl(i, 0)+v_tl(i, 1))-&
&             dy(i, 1)*(v_tl(i, -1)+v_tl(i, 2)))/(dy(i, 1)+dy(i, 2))
            x0 = 0.5*((2.*dy(i, 1)+dy(i, 2))*(v(i, 0)+v(i, 1))-dy(i, 1)*&
&             (v(i, -1)+v(i, 2)))/(dy(i, 1)+dy(i, 2))
            x1_tl = s15*v_tl(i, 1) + s11*v_tl(i, 2) - s14*dm_tl(i, 2)
            x1 = s15*v(i, 1) + s11*v(i, 2) - s14*dm(i, 2)
!           dm(i,1) = x1 - v(i,1)
            dm_tl(i, 1) = 0.5*(x1_tl-x0_tl)
            dm(i, 1) = 0.5*(x1-x0)
            IF (dm(i, 1) .GE. 0.) THEN
              x3_tl = dm_tl(i, 1)
              x3 = dm(i, 1)
            ELSE
              x3_tl = -dm_tl(i, 1)
              x3 = -dm(i, 1)
            END IF
            IF (v(i, 1) .LT. x0) THEN
              IF (x0 .LT. x1) THEN
                max7_tl = x1_tl
                max7 = x1
              ELSE
                max7_tl = x0_tl
                max7 = x0
              END IF
            ELSE IF (v(i, 1) .LT. x1) THEN
              max7_tl = x1_tl
              max7 = x1
            ELSE
              max7_tl = v_tl(i, 1)
              max7 = v(i, 1)
            END IF
            y2_tl = max7_tl - v_tl(i, 1)
            y2 = max7 - v(i, 1)
            IF (v(i, 1) .GT. x0) THEN
              IF (x0 .GT. x1) THEN
                min15_tl = x1_tl
                min15 = x1
              ELSE
                min15_tl = x0_tl
                min15 = x0
              END IF
            ELSE IF (v(i, 1) .GT. x1) THEN
              min15_tl = x1_tl
              min15 = x1
            ELSE
              min15_tl = v_tl(i, 1)
              min15 = v(i, 1)
            END IF
            z2_tl = v_tl(i, 1) - min15_tl
            z2 = v(i, 1) - min15
            IF (x3 .GT. y2) THEN
              IF (y2 .GT. z2) THEN
                min2_tl = z2_tl
                min2 = z2
              ELSE
                min2_tl = y2_tl
                min2 = y2
              END IF
            ELSE IF (x3 .GT. z2) THEN
              min2_tl = z2_tl
              min2 = z2
            ELSE
              min2_tl = x3_tl
              min2 = x3
            END IF
            dm_tl(i, 1) = min2_tl*SIGN(1.d0, min2*dm(i, 1))
            dm(i, 1) = SIGN(min2, dm(i, 1))
            x1_tl = s15*v_tl(i, 0) + s11*v_tl(i, -1) + s14*dm_tl(i, -1)
            x1 = s15*v(i, 0) + s11*v(i, -1) + s14*dm(i, -1)
!           dm(i,0) = v(i,0) - x1
            dm_tl(i, 0) = 0.5*(x0_tl-x1_tl)
            dm(i, 0) = 0.5*(x0-x1)
            IF (dm(i, 0) .GE. 0.) THEN
              x4_tl = dm_tl(i, 0)
              x4 = dm(i, 0)
            ELSE
              x4_tl = -dm_tl(i, 0)
              x4 = -dm(i, 0)
            END IF
            IF (v(i, 0) .LT. x0) THEN
              IF (x0 .LT. x1) THEN
                max8_tl = x1_tl
                max8 = x1
              ELSE
                max8_tl = x0_tl
                max8 = x0
              END IF
            ELSE IF (v(i, 0) .LT. x1) THEN
              max8_tl = x1_tl
              max8 = x1
            ELSE
              max8_tl = v_tl(i, 0)
              max8 = v(i, 0)
            END IF
            y3_tl = max8_tl - v_tl(i, 0)
            y3 = max8 - v(i, 0)
            IF (v(i, 0) .GT. x0) THEN
              IF (x0 .GT. x1) THEN
                min16_tl = x1_tl
                min16 = x1
              ELSE
                min16_tl = x0_tl
                min16 = x0
              END IF
            ELSE IF (v(i, 0) .GT. x1) THEN
              min16_tl = x1_tl
              min16 = x1
            ELSE
              min16_tl = v_tl(i, 0)
              min16 = v(i, 0)
            END IF
            z3_tl = v_tl(i, 0) - min16_tl
            z3 = v(i, 0) - min16
            IF (x4 .GT. y3) THEN
              IF (y3 .GT. z3) THEN
                min3_tl = z3_tl
                min3 = z3
              ELSE
                min3_tl = y3_tl
                min3 = y3
              END IF
            ELSE IF (x4 .GT. z3) THEN
              min3_tl = z3_tl
              min3 = z3
            ELSE
              min3_tl = x4_tl
              min3 = x4
            END IF
            dm_tl(i, 0) = min3_tl*SIGN(1.d0, min3*dm(i, 0))
            dm(i, 0) = SIGN(min3, dm(i, 0))
          END DO
          IF (is .EQ. 1) THEN
            dm_tl(1, 0) = 0.0_8
            dm(1, 0) = 0.
            dm_tl(1, 1) = 0.0_8
            dm(1, 1) = 0.
          END IF
          IF (ie + 1 .EQ. npx) THEN
            dm_tl(npx, 0) = 0.0_8
            dm(npx, 0) = 0.
            dm_tl(npx, 1) = 0.0_8
            dm(npx, 1) = 0.
          END IF
        END IF
        IF (je + 1 .EQ. npy) THEN
          DO i=is,ie+1
            x0_tl = 0.5*((2.*dy(i, npy-1)+dy(i, npy-2))*(v_tl(i, npy-1)+&
&             v_tl(i, npy))-dy(i, npy-1)*(v_tl(i, npy-2)+v_tl(i, npy+1))&
&             )/(dy(i, npy-1)+dy(i, npy-2))
            x0 = 0.5*((2.*dy(i, npy-1)+dy(i, npy-2))*(v(i, npy-1)+v(i, &
&             npy))-dy(i, npy-1)*(v(i, npy-2)+v(i, npy+1)))/(dy(i, npy-1&
&             )+dy(i, npy-2))
            x1_tl = s15*v_tl(i, npy-1) + s11*v_tl(i, npy-2) + s14*dm_tl(&
&             i, npy-2)
            x1 = s15*v(i, npy-1) + s11*v(i, npy-2) + s14*dm(i, npy-2)
!           dm(i,npy-1) = v(i,npy-1) - x1
            dm_tl(i, npy-1) = 0.5*(x0_tl-x1_tl)
            dm(i, npy-1) = 0.5*(x0-x1)
            IF (dm(i, npy-1) .GE. 0.) THEN
              x5_tl = dm_tl(i, npy-1)
              x5 = dm(i, npy-1)
            ELSE
              x5_tl = -dm_tl(i, npy-1)
              x5 = -dm(i, npy-1)
            END IF
            IF (v(i, npy-1) .LT. x0) THEN
              IF (x0 .LT. x1) THEN
                max9_tl = x1_tl
                max9 = x1
              ELSE
                max9_tl = x0_tl
                max9 = x0
              END IF
            ELSE IF (v(i, npy-1) .LT. x1) THEN
              max9_tl = x1_tl
              max9 = x1
            ELSE
              max9_tl = v_tl(i, npy-1)
              max9 = v(i, npy-1)
            END IF
            y4_tl = max9_tl - v_tl(i, npy-1)
            y4 = max9 - v(i, npy-1)
            IF (v(i, npy-1) .GT. x0) THEN
              IF (x0 .GT. x1) THEN
                min17_tl = x1_tl
                min17 = x1
              ELSE
                min17_tl = x0_tl
                min17 = x0
              END IF
            ELSE IF (v(i, npy-1) .GT. x1) THEN
              min17_tl = x1_tl
              min17 = x1
            ELSE
              min17_tl = v_tl(i, npy-1)
              min17 = v(i, npy-1)
            END IF
            z4_tl = v_tl(i, npy-1) - min17_tl
            z4 = v(i, npy-1) - min17
            IF (x5 .GT. y4) THEN
              IF (y4 .GT. z4) THEN
                min4_tl = z4_tl
                min4 = z4
              ELSE
                min4_tl = y4_tl
                min4 = y4
              END IF
            ELSE IF (x5 .GT. z4) THEN
              min4_tl = z4_tl
              min4 = z4
            ELSE
              min4_tl = x5_tl
              min4 = x5
            END IF
            dm_tl(i, npy-1) = min4_tl*SIGN(1.d0, min4*dm(i, npy-1))
            dm(i, npy-1) = SIGN(min4, dm(i, npy-1))
            x1_tl = s15*v_tl(i, npy) + s11*v_tl(i, npy+1) - s14*dm_tl(i&
&             , npy+1)
            x1 = s15*v(i, npy) + s11*v(i, npy+1) - s14*dm(i, npy+1)
!           dm(i,npy) = x1 - v(i,npy)
            dm_tl(i, npy) = 0.5*(x1_tl-x0_tl)
            dm(i, npy) = 0.5*(x1-x0)
            IF (dm(i, npy) .GE. 0.) THEN
              x6_tl = dm_tl(i, npy)
              x6 = dm(i, npy)
            ELSE
              x6_tl = -dm_tl(i, npy)
              x6 = -dm(i, npy)
            END IF
            IF (v(i, npy) .LT. x0) THEN
              IF (x0 .LT. x1) THEN
                max10_tl = x1_tl
                max10 = x1
              ELSE
                max10_tl = x0_tl
                max10 = x0
              END IF
            ELSE IF (v(i, npy) .LT. x1) THEN
              max10_tl = x1_tl
              max10 = x1
            ELSE
              max10_tl = v_tl(i, npy)
              max10 = v(i, npy)
            END IF
            y5_tl = max10_tl - v_tl(i, npy)
            y5 = max10 - v(i, npy)
            IF (v(i, npy) .GT. x0) THEN
              IF (x0 .GT. x1) THEN
                min18_tl = x1_tl
                min18 = x1
              ELSE
                min18_tl = x0_tl
                min18 = x0
              END IF
            ELSE IF (v(i, npy) .GT. x1) THEN
              min18_tl = x1_tl
              min18 = x1
            ELSE
              min18_tl = v_tl(i, npy)
              min18 = v(i, npy)
            END IF
            z5_tl = v_tl(i, npy) - min18_tl
            z5 = v(i, npy) - min18
            IF (x6 .GT. y5) THEN
              IF (y5 .GT. z5) THEN
                min5_tl = z5_tl
                min5 = z5
              ELSE
                min5_tl = y5_tl
                min5 = y5
              END IF
            ELSE IF (x6 .GT. z5) THEN
              min5_tl = z5_tl
              min5 = z5
            ELSE
              min5_tl = x6_tl
              min5 = x6
            END IF
            dm_tl(i, npy) = min5_tl*SIGN(1.d0, min5*dm(i, npy))
            dm(i, npy) = SIGN(min5, dm(i, npy))
          END DO
          IF (is .EQ. 1) THEN
            dm_tl(1, npy-1) = 0.0_8
            dm(1, npy-1) = 0.
            dm_tl(1, npy) = 0.0_8
            dm(1, npy) = 0.
          END IF
          IF (ie + 1 .EQ. npx) THEN
            dm_tl(npx, npy-1) = 0.0_8
            dm(npx, npy-1) = 0.
            dm_tl(npx, npy) = 0.0_8
            dm(npx, npy) = 0.
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
            flux_tl(i, j) = v_tl(i, j-1) + (1.-c(i, j)*rdy(i, j-1))*&
&             dm_tl(i, j-1) - rdy(i, j-1)*c_tl(i, j)*dm(i, j-1)
            flux(i, j) = v(i, j-1) + (1.-c(i, j)*rdy(i, j-1))*dm(i, j-1)
          ELSE
            flux_tl(i, j) = v_tl(i, j) - rdy(i, j)*c_tl(i, j)*dm(i, j) -&
&             (1.+c(i, j)*rdy(i, j))*dm_tl(i, j)
            flux(i, j) = v(i, j) - (1.+c(i, j)*rdy(i, j))*dm(i, j)
          END IF
        END DO
      END DO
    CASE (4) 
      dm_tl = 0.0_8
      DO j=js-2,je+2
        DO i=is,ie+1
          xt_tl = 0.25*(v_tl(i, j+1)-v_tl(i, j-1))
          xt = 0.25*(v(i, j+1)-v(i, j-1))
          IF (xt .GE. 0.) THEN
            x7_tl = xt_tl
            x7 = xt
          ELSE
            x7_tl = -xt_tl
            x7 = -xt
          END IF
          IF (v(i, j-1) .LT. v(i, j)) THEN
            IF (v(i, j) .LT. v(i, j+1)) THEN
              max11_tl = v_tl(i, j+1)
              max11 = v(i, j+1)
            ELSE
              max11_tl = v_tl(i, j)
              max11 = v(i, j)
            END IF
          ELSE IF (v(i, j-1) .LT. v(i, j+1)) THEN
            max11_tl = v_tl(i, j+1)
            max11 = v(i, j+1)
          ELSE
            max11_tl = v_tl(i, j-1)
            max11 = v(i, j-1)
          END IF
          y6_tl = max11_tl - v_tl(i, j)
          y6 = max11 - v(i, j)
          IF (v(i, j-1) .GT. v(i, j)) THEN
            IF (v(i, j) .GT. v(i, j+1)) THEN
              min19_tl = v_tl(i, j+1)
              min19 = v(i, j+1)
            ELSE
              min19_tl = v_tl(i, j)
              min19 = v(i, j)
            END IF
          ELSE IF (v(i, j-1) .GT. v(i, j+1)) THEN
            min19_tl = v_tl(i, j+1)
            min19 = v(i, j+1)
          ELSE
            min19_tl = v_tl(i, j-1)
            min19 = v(i, j-1)
          END IF
          z6_tl = v_tl(i, j) - min19_tl
          z6 = v(i, j) - min19
          IF (x7 .GT. y6) THEN
            IF (y6 .GT. z6) THEN
              min6_tl = z6_tl
              min6 = z6
            ELSE
              min6_tl = y6_tl
              min6 = y6
            END IF
          ELSE IF (x7 .GT. z6) THEN
            min6_tl = z6_tl
            min6 = z6
          ELSE
            min6_tl = x7_tl
            min6 = x7
          END IF
          dm_tl(i, j) = min6_tl*SIGN(1.d0, min6*xt)
          dm(i, j) = SIGN(min6, xt)
        END DO
      END DO
      al_tl = 0.0_8
      DO j=js-1,je+2
        DO i=is,ie+1
          al_tl(i, j) = 0.5*(v_tl(i, j-1)+v_tl(i, j)) + r3*(dm_tl(i, j-1&
&           )-dm_tl(i, j))
          al(i, j) = 0.5*(v(i, j-1)+v(i, j)) + r3*(dm(i, j-1)-dm(i, j))
        END DO
      END DO
      IF (grid_type .LT. 3) THEN
        IF (js .EQ. 1) THEN
          DO i=is,ie+1
            x0_tl = 0.5*((2.*dy(i, 1)+dy(i, 2))*(v_tl(i, 0)+v_tl(i, 1))-&
&             dy(i, 1)*(v_tl(i, -1)+v_tl(i, 2)))/(dy(i, 1)+dy(i, 2))
            x0 = 0.5*((2.*dy(i, 1)+dy(i, 2))*(v(i, 0)+v(i, 1))-dy(i, 1)*&
&             (v(i, -1)+v(i, 2)))/(dy(i, 1)+dy(i, 2))
            x1_tl = s15*v_tl(i, 1) + s11*v_tl(i, 2) - s14*dm_tl(i, 2)
            x1 = s15*v(i, 1) + s11*v(i, 2) - s14*dm(i, 2)
            dm_tl(i, 1) = 0.5*(x1_tl-x0_tl)
            dm(i, 1) = 0.5*(x1-x0)
!           dm(i,1) = sign(min(abs(dm(i,1)), max(v(i,1), x0, x1) - v(i,1),   &
!                                   v(i,1) - min(v(i,1), x0, x1)), dm(i,1))
            x1_tl = s15*v_tl(i, 0) + s11*v_tl(i, -1) + s14*dm_tl(i, -1)
            x1 = s15*v(i, 0) + s11*v(i, -1) + s14*dm(i, -1)
            dm_tl(i, 0) = 0.5*(x0_tl-x1_tl)
            dm(i, 0) = 0.5*(x0-x1)
!           dm(i,0) = sign(min(abs(dm(i,0)), max(v(i,0), x0, x1) - v(i,0),   &
!                                   v(i,0) - min(v(i,0), x0, x1)), dm(i,0))
            al_tl(i, 0) = 0.5*(v_tl(i, -1)+v_tl(i, 0)) + r3*(dm_tl(i, -1&
&             )-dm_tl(i, 0))
            al(i, 0) = 0.5*(v(i, -1)+v(i, 0)) + r3*(dm(i, -1)-dm(i, 0))
            al_tl(i, 1) = x0_tl
            al(i, 1) = x0
            al_tl(i, 2) = 0.5*(v_tl(i, 1)+v_tl(i, 2)) + r3*(dm_tl(i, 1)-&
&             dm_tl(i, 2))
            al(i, 2) = 0.5*(v(i, 1)+v(i, 2)) + r3*(dm(i, 1)-dm(i, 2))
          END DO
          IF (is .EQ. 1) THEN
            dm_tl(1, 0) = 0.0_8
            dm(1, 0) = 0.
            dm_tl(1, 1) = 0.0_8
            dm(1, 1) = 0.
            i = 1
            al_tl(i, 0) = 0.5*(v_tl(i, -1)+v_tl(i, 0)) + r3*(dm_tl(i, -1&
&             )-dm_tl(i, 0))
            al(i, 0) = 0.5*(v(i, -1)+v(i, 0)) + r3*(dm(i, -1)-dm(i, 0))
            al_tl(i, 2) = 0.5*(v_tl(i, 1)+v_tl(i, 2)) + r3*(dm_tl(i, 1)-&
&             dm_tl(i, 2))
            al(i, 2) = 0.5*(v(i, 1)+v(i, 2)) + r3*(dm(i, 1)-dm(i, 2))
          END IF
          IF (ie + 1 .EQ. npx) THEN
            dm_tl(npx, 0) = 0.0_8
            dm(npx, 0) = 0.
            dm_tl(npx, 1) = 0.0_8
            dm(npx, 1) = 0.
            i = npx
            al_tl(i, 0) = 0.5*(v_tl(i, -1)+v_tl(i, 0)) + r3*dm_tl(i, -1)
            al(i, 0) = 0.5*(v(i, -1)+v(i, 0)) + r3*dm(i, -1)
            al_tl(i, 2) = 0.5*(v_tl(i, 1)+v_tl(i, 2)) - r3*dm_tl(i, 2)
            al(i, 2) = 0.5*(v(i, 1)+v(i, 2)) - r3*dm(i, 2)
          END IF
        END IF
        IF (je + 1 .EQ. npy) THEN
          DO i=is,ie+1
            x0_tl = 0.5*((2.*dy(i, npy-1)+dy(i, npy-2))*(v_tl(i, npy-1)+&
&             v_tl(i, npy))-dy(i, npy-1)*(v_tl(i, npy-2)+v_tl(i, npy+1))&
&             )/(dy(i, npy-1)+dy(i, npy-2))
            x0 = 0.5*((2.*dy(i, npy-1)+dy(i, npy-2))*(v(i, npy-1)+v(i, &
&             npy))-dy(i, npy-1)*(v(i, npy-2)+v(i, npy+1)))/(dy(i, npy-1&
&             )+dy(i, npy-2))
            x1_tl = s15*v_tl(i, npy-1) + s11*v_tl(i, npy-2) + s14*dm_tl(&
&             i, npy-2)
            x1 = s15*v(i, npy-1) + s11*v(i, npy-2) + s14*dm(i, npy-2)
            dm_tl(i, npy-1) = 0.5*(x0_tl-x1_tl)
            dm(i, npy-1) = 0.5*(x0-x1)
!           dm(i,npy-1) = sign(min(abs(dm(i,npy-1)), max(v(i,npy-1), x0, x1) - v(i,npy-1),  &
!                                       v(i,npy-1) - min(v(i,npy-1), x0, x1)), dm(i,npy-1))
            x1_tl = s15*v_tl(i, npy) + s11*v_tl(i, npy+1) - s14*dm_tl(i&
&             , npy+1)
            x1 = s15*v(i, npy) + s11*v(i, npy+1) - s14*dm(i, npy+1)
            dm_tl(i, npy) = 0.5*(x1_tl-x0_tl)
            dm(i, npy) = 0.5*(x1-x0)
!           dm(i,npy) = sign(min(abs(dm(i,npy)), max(v(i,npy), x0, x1) - v(i,npy),   &
!                                     v(i,npy) - min(v(i,npy), x0, x1)), dm(i,npy))
            al_tl(i, npy-1) = 0.5*(v_tl(i, npy-2)+v_tl(i, npy-1)) + r3*(&
&             dm_tl(i, npy-2)-dm_tl(i, npy-1))
            al(i, npy-1) = 0.5*(v(i, npy-2)+v(i, npy-1)) + r3*(dm(i, npy&
&             -2)-dm(i, npy-1))
            al_tl(i, npy) = x0_tl
            al(i, npy) = x0
            al_tl(i, npy+1) = 0.5*(v_tl(i, npy)+v_tl(i, npy+1)) + r3*(&
&             dm_tl(i, npy)-dm_tl(i, npy+1))
            al(i, npy+1) = 0.5*(v(i, npy)+v(i, npy+1)) + r3*(dm(i, npy)-&
&             dm(i, npy+1))
          END DO
          IF (is .EQ. 1) THEN
            dm_tl(1, npy-1) = 0.0_8
            dm(1, npy-1) = 0.
            dm_tl(1, npy) = 0.0_8
            dm(1, npy) = 0.
            i = 1
            al_tl(i, npy-1) = 0.5*(v_tl(i, npy-2)+v_tl(i, npy-1)) + r3*&
&             dm_tl(i, npy-2)
            al(i, npy-1) = 0.5*(v(i, npy-2)+v(i, npy-1)) + r3*dm(i, npy-&
&             2)
            al_tl(i, npy+1) = 0.5*(v_tl(i, npy)+v_tl(i, npy+1)) - r3*&
&             dm_tl(i, npy+1)
            al(i, npy+1) = 0.5*(v(i, npy)+v(i, npy+1)) - r3*dm(i, npy+1)
          END IF
          IF (ie + 1 .EQ. npx) THEN
            dm_tl(npx, npy-1) = 0.0_8
            dm(npx, npy-1) = 0.
            dm_tl(npx, npy) = 0.0_8
            dm(npx, npy) = 0.
            i = npx
            al_tl(i, npy-1) = 0.5*(v_tl(i, npy-2)+v_tl(i, npy-1)) + r3*&
&             dm_tl(i, npy-2)
            al(i, npy-1) = 0.5*(v(i, npy-2)+v(i, npy-1)) + r3*dm(i, npy-&
&             2)
            al_tl(i, npy+1) = 0.5*(v_tl(i, npy)+v_tl(i, npy+1)) - r3*&
&             dm_tl(i, npy+1)
            al(i, npy+1) = 0.5*(v(i, npy)+v(i, npy+1)) - r3*dm(i, npy+1)
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
            xt_tl = 2.*dm_tl(i, j-1)
            xt = 2.*dm(i, j-1)
            IF (xt .GE. 0.) THEN
              x8_tl = xt_tl
              x8 = xt
            ELSE
              x8_tl = -xt_tl
              x8 = -xt
            END IF
            IF (al(i, j-1) - v(i, j-1) .GE. 0.) THEN
              y7_tl = al_tl(i, j-1) - v_tl(i, j-1)
              y7 = al(i, j-1) - v(i, j-1)
            ELSE
              y7_tl = -(al_tl(i, j-1)-v_tl(i, j-1))
              y7 = -(al(i, j-1)-v(i, j-1))
            END IF
            IF (x8 .GT. y7) THEN
              min7_tl = y7_tl
              min7 = y7
            ELSE
              min7_tl = x8_tl
              min7 = x8
            END IF
            dl_tl = min7_tl*SIGN(1.d0, min7*xt)
            dl = SIGN(min7, xt)
            IF (xt .GE. 0.) THEN
              x9_tl = xt_tl
              x9 = xt
            ELSE
              x9_tl = -xt_tl
              x9 = -xt
            END IF
            IF (al(i, j) - v(i, j-1) .GE. 0.) THEN
              y8_tl = al_tl(i, j) - v_tl(i, j-1)
              y8 = al(i, j) - v(i, j-1)
            ELSE
              y8_tl = -(al_tl(i, j)-v_tl(i, j-1))
              y8 = -(al(i, j)-v(i, j-1))
            END IF
            IF (x9 .GT. y8) THEN
              min8_tl = y8_tl
              min8 = y8
            ELSE
              min8_tl = x9_tl
              min8 = x9
            END IF
            dr_tl = min8_tl*SIGN(1.d0, min8*xt)
            dr = SIGN(min8, xt)
            cfl_tl = rdy(i, j-1)*c_tl(i, j)
            cfl = c(i, j)*rdy(i, j-1)
            flux_tl(i, j) = v_tl(i, j-1) + (1.-cfl)*(dr_tl+cfl_tl*(dl-dr&
&             )+cfl*(dl_tl-dr_tl)) - cfl_tl*(dr+cfl*(dl-dr))
            flux(i, j) = v(i, j-1) + (1.-cfl)*(dr+cfl*(dl-dr))
          ELSE
            xt_tl = 2.*dm_tl(i, j)
            xt = 2.*dm(i, j)
            IF (xt .GE. 0.) THEN
              x10_tl = xt_tl
              x10 = xt
            ELSE
              x10_tl = -xt_tl
              x10 = -xt
            END IF
            IF (al(i, j) - v(i, j) .GE. 0.) THEN
              y9_tl = al_tl(i, j) - v_tl(i, j)
              y9 = al(i, j) - v(i, j)
            ELSE
              y9_tl = -(al_tl(i, j)-v_tl(i, j))
              y9 = -(al(i, j)-v(i, j))
            END IF
            IF (x10 .GT. y9) THEN
              min9_tl = y9_tl
              min9 = y9
            ELSE
              min9_tl = x10_tl
              min9 = x10
            END IF
            dl_tl = min9_tl*SIGN(1.d0, min9*xt)
            dl = SIGN(min9, xt)
            IF (xt .GE. 0.) THEN
              x11_tl = xt_tl
              x11 = xt
            ELSE
              x11_tl = -xt_tl
              x11 = -xt
            END IF
            IF (al(i, j+1) - v(i, j) .GE. 0.) THEN
              y10_tl = al_tl(i, j+1) - v_tl(i, j)
              y10 = al(i, j+1) - v(i, j)
            ELSE
              y10_tl = -(al_tl(i, j+1)-v_tl(i, j))
              y10 = -(al(i, j+1)-v(i, j))
            END IF
            IF (x11 .GT. y10) THEN
              min10_tl = y10_tl
              min10 = y10
            ELSE
              min10_tl = x11_tl
              min10 = x11
            END IF
            dr_tl = min10_tl*SIGN(1.d0, min10*xt)
            dr = SIGN(min10, xt)
            cfl_tl = rdy(i, j)*c_tl(i, j)
            cfl = c(i, j)*rdy(i, j)
            flux_tl(i, j) = v_tl(i, j) - cfl_tl*(dl+cfl*(dl-dr)) - (1.+&
&             cfl)*(dl_tl+cfl_tl*(dl-dr)+cfl*(dl_tl-dr_tl))
            flux(i, j) = v(i, j) - (1.+cfl)*(dl+cfl*(dl-dr))
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
        min20 = je + 1
        bl_tl = 0.0_8
        br_tl = 0.0_8
      ELSE
        min20 = npy - 3
        bl_tl = 0.0_8
        br_tl = 0.0_8
      END IF
      DO j=max1,min20
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
      dm_tl = 0.0_8
! jord= 8, 9, 10
      DO j=js-2,je+2
        DO i=is,ie+1
          xt_tl = 0.25*(v_tl(i, j+1)-v_tl(i, j-1))
          xt = 0.25*(v(i, j+1)-v(i, j-1))
          IF (xt .GE. 0.) THEN
            x12_tl = xt_tl
            x12 = xt
          ELSE
            x12_tl = -xt_tl
            x12 = -xt
          END IF
          IF (v(i, j-1) .LT. v(i, j)) THEN
            IF (v(i, j) .LT. v(i, j+1)) THEN
              max12_tl = v_tl(i, j+1)
              max12 = v(i, j+1)
            ELSE
              max12_tl = v_tl(i, j)
              max12 = v(i, j)
            END IF
          ELSE IF (v(i, j-1) .LT. v(i, j+1)) THEN
            max12_tl = v_tl(i, j+1)
            max12 = v(i, j+1)
          ELSE
            max12_tl = v_tl(i, j-1)
            max12 = v(i, j-1)
          END IF
          y11_tl = max12_tl - v_tl(i, j)
          y11 = max12 - v(i, j)
          IF (v(i, j-1) .GT. v(i, j)) THEN
            IF (v(i, j) .GT. v(i, j+1)) THEN
              min21_tl = v_tl(i, j+1)
              min21 = v(i, j+1)
            ELSE
              min21_tl = v_tl(i, j)
              min21 = v(i, j)
            END IF
          ELSE IF (v(i, j-1) .GT. v(i, j+1)) THEN
            min21_tl = v_tl(i, j+1)
            min21 = v(i, j+1)
          ELSE
            min21_tl = v_tl(i, j-1)
            min21 = v(i, j-1)
          END IF
          z7_tl = v_tl(i, j) - min21_tl
          z7 = v(i, j) - min21
          IF (x12 .GT. y11) THEN
            IF (y11 .GT. z7) THEN
              min11_tl = z7_tl
              min11 = z7
            ELSE
              min11_tl = y11_tl
              min11 = y11
            END IF
          ELSE IF (x12 .GT. z7) THEN
            min11_tl = z7_tl
            min11 = z7
          ELSE
            min11_tl = x12_tl
            min11 = x12
          END IF
          dm_tl(i, j) = min11_tl*SIGN(1.d0, min11*xt)
          dm(i, j) = SIGN(min11, xt)
        END DO
      END DO
      dq_tl = 0.0_8
      DO j=js-3,je+2
        DO i=is,ie+1
          dq_tl(i, j) = v_tl(i, j+1) - v_tl(i, j)
          dq(i, j) = v(i, j+1) - v(i, j)
        END DO
      END DO
      IF (grid_type .LT. 3) THEN
        IF (3 .LT. js - 1) THEN
          max2 = js - 1
        ELSE
          max2 = 3
        END IF
        IF (npy - 2 .GT. je + 2) THEN
          min22 = je + 2
          al_tl = 0.0_8
        ELSE
          min22 = npy - 2
          al_tl = 0.0_8
        END IF
        DO j=max2,min22
          DO i=is,ie+1
            al_tl(i, j) = 0.5*(v_tl(i, j-1)+v_tl(i, j)) + r3*(dm_tl(i, j&
&             -1)-dm_tl(i, j))
            al(i, j) = 0.5*(v(i, j-1)+v(i, j)) + r3*(dm(i, j-1)-dm(i, j)&
&             )
          END DO
        END DO
        IF (jord .EQ. 8) THEN
          IF (3 .LT. js - 1) THEN
            max3 = js - 1
          ELSE
            max3 = 3
          END IF
          IF (npy - 3 .GT. je + 1) THEN
            min23 = je + 1
            bl_tl = 0.0_8
            br_tl = 0.0_8
          ELSE
            min23 = npy - 3
            bl_tl = 0.0_8
            br_tl = 0.0_8
          END IF
          DO j=max3,min23
            DO i=is,ie+1
              xt_tl = 2.*dm_tl(i, j)
              xt = 2.*dm(i, j)
              IF (xt .GE. 0.) THEN
                x13_tl = xt_tl
                x13 = xt
              ELSE
                x13_tl = -xt_tl
                x13 = -xt
              END IF
              IF (al(i, j) - v(i, j) .GE. 0.) THEN
                y12_tl = al_tl(i, j) - v_tl(i, j)
                y12 = al(i, j) - v(i, j)
              ELSE
                y12_tl = -(al_tl(i, j)-v_tl(i, j))
                y12 = -(al(i, j)-v(i, j))
              END IF
              IF (x13 .GT. y12) THEN
                min12_tl = y12_tl
                min12 = y12
              ELSE
                min12_tl = x13_tl
                min12 = x13
              END IF
              bl_tl(i, j) = -(min12_tl*SIGN(1.d0, min12*xt))
              bl(i, j) = -SIGN(min12, xt)
              IF (xt .GE. 0.) THEN
                x14_tl = xt_tl
                x14 = xt
              ELSE
                x14_tl = -xt_tl
                x14 = -xt
              END IF
              IF (al(i, j+1) - v(i, j) .GE. 0.) THEN
                y13_tl = al_tl(i, j+1) - v_tl(i, j)
                y13 = al(i, j+1) - v(i, j)
              ELSE
                y13_tl = -(al_tl(i, j+1)-v_tl(i, j))
                y13 = -(al(i, j+1)-v(i, j))
              END IF
              IF (x14 .GT. y13) THEN
                min13_tl = y13_tl
                min13 = y13
              ELSE
                min13_tl = x14_tl
                min13 = x14
              END IF
              br_tl(i, j) = min13_tl*SIGN(1.d0, min13*xt)
              br(i, j) = SIGN(min13, xt)
            END DO
          END DO
        ELSE IF (jord .EQ. 9) THEN
          IF (3 .LT. js - 1) THEN
            max4 = js - 1
          ELSE
            max4 = 3
          END IF
          IF (npy - 3 .GT. je + 1) THEN
            min24 = je + 1
            bl_tl = 0.0_8
            br_tl = 0.0_8
          ELSE
            min24 = npy - 3
            bl_tl = 0.0_8
            br_tl = 0.0_8
          END IF
          DO j=max4,min24
            DO i=is,ie+1
              pmp_tl = 2.*dq_tl(i, j-1)
              pmp = 2.*dq(i, j-1)
              lac_tl = pmp_tl - 1.5*dq_tl(i, j-2)
              lac = pmp - 1.5*dq(i, j-2)
              IF (0. .LT. pmp) THEN
                IF (pmp .LT. lac) THEN
                  x15_tl = lac_tl
                  x15 = lac
                ELSE
                  x15_tl = pmp_tl
                  x15 = pmp
                END IF
              ELSE IF (0. .LT. lac) THEN
                x15_tl = lac_tl
                x15 = lac
              ELSE
                x15 = 0.
                x15_tl = 0.0_8
              END IF
              IF (0. .GT. pmp) THEN
                IF (pmp .GT. lac) THEN
                  y18_tl = lac_tl
                  y18 = lac
                ELSE
                  y18_tl = pmp_tl
                  y18 = pmp
                END IF
              ELSE IF (0. .GT. lac) THEN
                y18_tl = lac_tl
                y18 = lac
              ELSE
                y18 = 0.
                y18_tl = 0.0_8
              END IF
              IF (al(i, j+1) - v(i, j) .LT. y18) THEN
                y14_tl = y18_tl
                y14 = y18
              ELSE
                y14_tl = al_tl(i, j+1) - v_tl(i, j)
                y14 = al(i, j+1) - v(i, j)
              END IF
              IF (x15 .GT. y14) THEN
                br_tl(i, j) = y14_tl
                br(i, j) = y14
              ELSE
                br_tl(i, j) = x15_tl
                br(i, j) = x15
              END IF
              pmp_tl = -(2.*dq_tl(i, j))
              pmp = -(2.*dq(i, j))
              lac_tl = pmp_tl + 1.5*dq_tl(i, j+1)
              lac = pmp + 1.5*dq(i, j+1)
              IF (0. .LT. pmp) THEN
                IF (pmp .LT. lac) THEN
                  x16_tl = lac_tl
                  x16 = lac
                ELSE
                  x16_tl = pmp_tl
                  x16 = pmp
                END IF
              ELSE IF (0. .LT. lac) THEN
                x16_tl = lac_tl
                x16 = lac
              ELSE
                x16 = 0.
                x16_tl = 0.0_8
              END IF
              IF (0. .GT. pmp) THEN
                IF (pmp .GT. lac) THEN
                  y19_tl = lac_tl
                  y19 = lac
                ELSE
                  y19_tl = pmp_tl
                  y19 = pmp
                END IF
              ELSE IF (0. .GT. lac) THEN
                y19_tl = lac_tl
                y19 = lac
              ELSE
                y19 = 0.
                y19_tl = 0.0_8
              END IF
              IF (al(i, j) - v(i, j) .LT. y19) THEN
                y15_tl = y19_tl
                y15 = y19
              ELSE
                y15_tl = al_tl(i, j) - v_tl(i, j)
                y15 = al(i, j) - v(i, j)
              END IF
              IF (x16 .GT. y15) THEN
                bl_tl(i, j) = y15_tl
                bl(i, j) = y15
              ELSE
                bl_tl(i, j) = x16_tl
                bl(i, j) = x16
              END IF
            END DO
          END DO
        ELSE
          IF (3 .LT. js - 1) THEN
            max5 = js - 1
          ELSE
            max5 = 3
          END IF
          IF (npy - 3 .GT. je + 1) THEN
            min25 = je + 1
            bl_tl = 0.0_8
            br_tl = 0.0_8
          ELSE
            min25 = npy - 3
            bl_tl = 0.0_8
            br_tl = 0.0_8
          END IF
! Unlimited:
          DO j=max5,min25
            DO i=is,ie+1
              bl_tl(i, j) = al_tl(i, j) - v_tl(i, j)
              bl(i, j) = al(i, j) - v(i, j)
              br_tl(i, j) = al_tl(i, j+1) - v_tl(i, j)
              br(i, j) = al(i, j+1) - v(i, j)
            END DO
          END DO
        END IF
!--------------
! fix the edges
!--------------
        IF (js .EQ. 1) THEN
          DO i=is,ie+1
            br_tl(i, 2) = al_tl(i, 3) - v_tl(i, 2)
            br(i, 2) = al(i, 3) - v(i, 2)
            xt_tl = s15*v_tl(i, 1) + s11*v_tl(i, 2) - s14*dm_tl(i, 2)
            xt = s15*v(i, 1) + s11*v(i, 2) - s14*dm(i, 2)
            br_tl(i, 1) = xt_tl - v_tl(i, 1)
            br(i, 1) = xt - v(i, 1)
            bl_tl(i, 2) = xt_tl - v_tl(i, 2)
            bl(i, 2) = xt - v(i, 2)
            bl_tl(i, 0) = s14*dm_tl(i, -1) - s11*dq_tl(i, -1)
            bl(i, 0) = s14*dm(i, -1) - s11*dq(i, -1)
            xt_tl = 0.5*((2.*dy(i, 1)+dy(i, 2))*(v_tl(i, 0)+v_tl(i, 1))-&
&             dy(i, 1)*(v_tl(i, -1)+v_tl(i, 2)))/(dy(i, 1)+dy(i, 2))
            xt = 0.5*((2.*dy(i, 1)+dy(i, 2))*(v(i, 0)+v(i, 1))-dy(i, 1)*&
&             (v(i, -1)+v(i, 2)))/(dy(i, 1)+dy(i, 2))
            bl_tl(i, 1) = xt_tl - v_tl(i, 1)
            bl(i, 1) = xt - v(i, 1)
            br_tl(i, 0) = xt_tl - v_tl(i, 0)
            br(i, 0) = xt - v(i, 0)
          END DO
!           br(i,0) = xt - 0.5*(u(i-1,1)+u(i,1))*cosa(i,1) - v(i,0)
!           bl(i,1) = xt + 0.5*(u(i-1,1)+u(i,1))*cosa(i,1) - v(i,1)
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
          j = 2
          IF (jord .LT. 10) CALL PERT_PPM_TLM(ie - is + 2, v(is:, j), bl&
&                                       (is:, j), bl_tl(is:, j), br(is:&
&                                       , j), br_tl(is:, j), -1)
        END IF
        IF (je + 1 .EQ. npy) THEN
          DO i=is,ie+1
            bl_tl(i, npy-2) = al_tl(i, npy-2) - v_tl(i, npy-2)
            bl(i, npy-2) = al(i, npy-2) - v(i, npy-2)
            xt_tl = s15*v_tl(i, npy-1) + s11*v_tl(i, npy-2) + s14*dm_tl(&
&             i, npy-2)
            xt = s15*v(i, npy-1) + s11*v(i, npy-2) + s14*dm(i, npy-2)
            br_tl(i, npy-2) = xt_tl - v_tl(i, npy-2)
            br(i, npy-2) = xt - v(i, npy-2)
            bl_tl(i, npy-1) = xt_tl - v_tl(i, npy-1)
            bl(i, npy-1) = xt - v(i, npy-1)
            br_tl(i, npy) = s11*dq_tl(i, npy) - s14*dm_tl(i, npy+1)
            br(i, npy) = s11*dq(i, npy) - s14*dm(i, npy+1)
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
!           br(i,npy-1) = xt + 0.5*(u(i-1,npy)+u(i,npy))*cosa(i,npy) - v(i,npy-1)
!           bl(i,npy  ) = xt - 0.5*(u(i-1,npy)+u(i,npy))*cosa(i,npy) - v(i,npy)
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
          END IF
          j = npy - 2
          IF (jord .LT. 10) THEN
            CALL PERT_PPM_TLM(ie - is + 2, v(is:, j), bl(is:, j), bl_tl(&
&                       is:, j), br(is:, j), br_tl(is:, j), -1)
            flux_tl = 0.0_8
          ELSE
            flux_tl = 0.0_8
          END IF
        ELSE
          flux_tl = 0.0_8
        END IF
      ELSE
        al_tl = 0.0_8
        DO j=js-1,je+2
          DO i=is,ie+1
            al_tl(i, j) = 0.5*(v_tl(i, j-1)+v_tl(i, j)) + r3*(dm_tl(i, j&
&             -1)-dm_tl(i, j))
            al(i, j) = 0.5*(v(i, j-1)+v(i, j)) + r3*(dm(i, j-1)-dm(i, j)&
&             )
          END DO
        END DO
        bl_tl = 0.0_8
        br_tl = 0.0_8
        DO j=js-1,je+1
          DO i=is,ie+1
            pmp_tl = 2.*dq_tl(i, j-1)
            pmp = 2.*dq(i, j-1)
            lac_tl = pmp_tl - 1.5*dq_tl(i, j-2)
            lac = pmp - 1.5*dq(i, j-2)
            IF (0. .LT. pmp) THEN
              IF (pmp .LT. lac) THEN
                x17_tl = lac_tl
                x17 = lac
              ELSE
                x17_tl = pmp_tl
                x17 = pmp
              END IF
            ELSE IF (0. .LT. lac) THEN
              x17_tl = lac_tl
              x17 = lac
            ELSE
              x17 = 0.
              x17_tl = 0.0_8
            END IF
            IF (0. .GT. pmp) THEN
              IF (pmp .GT. lac) THEN
                y20_tl = lac_tl
                y20 = lac
              ELSE
                y20_tl = pmp_tl
                y20 = pmp
              END IF
            ELSE IF (0. .GT. lac) THEN
              y20_tl = lac_tl
              y20 = lac
            ELSE
              y20 = 0.
              y20_tl = 0.0_8
            END IF
            IF (al(i, j+1) - v(i, j) .LT. y20) THEN
              y16_tl = y20_tl
              y16 = y20
            ELSE
              y16_tl = al_tl(i, j+1) - v_tl(i, j)
              y16 = al(i, j+1) - v(i, j)
            END IF
            IF (x17 .GT. y16) THEN
              br_tl(i, j) = y16_tl
              br(i, j) = y16
            ELSE
              br_tl(i, j) = x17_tl
              br(i, j) = x17
            END IF
            pmp_tl = -(2.*dq_tl(i, j))
            pmp = -(2.*dq(i, j))
            lac_tl = pmp_tl + 1.5*dq_tl(i, j+1)
            lac = pmp + 1.5*dq(i, j+1)
            IF (0. .LT. pmp) THEN
              IF (pmp .LT. lac) THEN
                x18_tl = lac_tl
                x18 = lac
              ELSE
                x18_tl = pmp_tl
                x18 = pmp
              END IF
            ELSE IF (0. .LT. lac) THEN
              x18_tl = lac_tl
              x18 = lac
            ELSE
              x18 = 0.
              x18_tl = 0.0_8
            END IF
            IF (0. .GT. pmp) THEN
              IF (pmp .GT. lac) THEN
                y21_tl = lac_tl
                y21 = lac
              ELSE
                y21_tl = pmp_tl
                y21 = pmp
              END IF
            ELSE IF (0. .GT. lac) THEN
              y21_tl = lac_tl
              y21 = lac
            ELSE
              y21 = 0.
              y21_tl = 0.0_8
            END IF
            IF (al(i, j) - v(i, j) .LT. y21) THEN
              y17_tl = y21_tl
              y17 = y21
            ELSE
              y17_tl = al_tl(i, j) - v_tl(i, j)
              y17 = al(i, j) - v(i, j)
            END IF
            IF (x18 .GT. y17) THEN
              bl_tl(i, j) = y17_tl
              bl(i, j) = y17
            ELSE
              bl_tl(i, j) = x18_tl
              bl(i, j) = x18
            END IF
          END DO
        END DO
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
    END SELECT
  END SUBROUTINE YTP_V_TLM

  SUBROUTINE D2A2C_VECT_TLM(u, u_tl, v, v_tl, ua, ua_tl, va, va_tl, uc, &
&   uc_tl, vc, vc_tl, ut, ut_tl, vt, vt_tl, dord4, dord4_tj)
    IMPLICIT NONE
    REAL, INTENT(IN) :: u(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: u_tl(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: v(isd:ied+1, jsd:jed)
    REAL, INTENT(IN) :: v_tl(isd:ied+1, jsd:jed)
    LOGICAL, INTENT(IN) :: dord4, dord4_tj
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
    INTEGER :: npt, i, j, ifirst, ilast, id, id_tj
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
    IF (dord4_tj) THEN
      id_tj = 1
    ELSE
      id_tj = 0
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
      utmp_tl = 0.0
    ELSE
      min1 = npy - npt
      utmp_tl = 0.0
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
      vtmp_tl = 0.0
    ELSE
      min3 = npy - npt
      vtmp_tl = 0.0
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
        ua(i, j) = (utmp(i, j)-vtmp(i, j)*cosa_s(i, j))*rsin2(i, j)
        va(i, j) = (vtmp(i, j)-utmp(i, j)*cosa_s(i, j))*rsin2(i, j)
      END DO
    END DO
    DO j=js-1-id_tj,je+1+id_tj
      DO i=is-1-id_tj,ie+1+id_tj
        ua_tl(i, j) = rsin2(i, j)*(utmp_tl(i, j)-cosa_s(i, j)*vtmp_tl(i&
&         , j))
        va_tl(i, j) = rsin2(i, j)*(vtmp_tl(i, j)-cosa_s(i, j)*utmp_tl(i&
&         , j))
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
&   , va, va_tl, uc, uc_tl, vc, vc_tl, dord4, dord4_tj)
    IMPLICIT NONE
    REAL, DIMENSION(isd:ied, jsd:jed + 1), INTENT(IN) :: u, um
    REAL, DIMENSION(isd:ied, jsd:jed+1), INTENT(IN) :: u_tl, um_tl
    REAL, DIMENSION(isd:ied + 1, jsd:jed), INTENT(IN) :: v, vm
    REAL, DIMENSION(isd:ied+1, jsd:jed), INTENT(IN) :: v_tl, vm_tl
    LOGICAL, INTENT(IN) :: dord4, dord4_tj
    REAL, DIMENSION(isd:ied + 1, jsd:jed), INTENT(OUT) :: uc
    REAL, DIMENSION(isd:ied+1, jsd:jed), INTENT(OUT) :: uc_tl
    REAL, DIMENSION(isd:ied, jsd:jed + 1), INTENT(OUT) :: vc
    REAL, DIMENSION(isd:ied, jsd:jed+1), INTENT(OUT) :: vc_tl
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(OUT) :: ua, va
    REAL, DIMENSION(isd:ied, jsd:jed), INTENT(OUT) :: ua_tl, va_tl
! Local 
    REAL, DIMENSION(isd:ied, jsd:jed) :: utmp, vtmp
    REAL, DIMENSION(isd:ied, jsd:jed) :: utmp_tl, vtmp_tl
    INTEGER :: npt, i, j, ifirst, ilast, id, id_tj
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
    IF (dord4_tj) THEN
      id_tj = 1
    ELSE
      id_tj = 0
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
      utmp_tl = 0.0
    ELSE
      min1 = npy - npt
      utmp_tl = 0.0
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
      vtmp_tl = 0.0
    ELSE
      min3 = npy - npt
      vtmp_tl = 0.0
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
        ua(i, j) = (utmp(i, j)-vtmp(i, j)*cosa_s(i, j))*rsin2(i, j)
        va(i, j) = (vtmp(i, j)-utmp(i, j)*cosa_s(i, j))*rsin2(i, j)
      END DO
    END DO
    DO j=js-1-id_tj,je+1+id_tj
      DO i=is-1-id_tj,ie+1+id_tj
        ua_tl(i, j) = rsin2(i, j)*(utmp_tl(i, j)-cosa_s(i, j)*vtmp_tl(i&
&         , j))
        va_tl(i, j) = rsin2(i, j)*(vtmp_tl(i, j)-cosa_s(i, j)*utmp_tl(i&
&         , j))
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

 end module sw_core_tlm_mod

