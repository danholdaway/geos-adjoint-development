SUBROUTINE g_v2_tracer_2d_1l(q, g_q, dp0, g_dp0, mfx, g_mfx, mfy, g_mfy, cx, g_cx, cy, g_cy, npx, npy, npz, nq, hord, q_split, k, &
&q3, g_q3, dt, uniform_ppm, id_divg)

!==============================================
! referencing used modules
!==============================================
use g_v2_tp_core_mod, only : g_v2_fv_tp_2d_trc
use g_fv_my_mpp, only : g_mp_reduce_max_dummy, g_mpp_update_domains_dummy4

!==============================================
! all entries are defined explicitly
!==============================================
implicit none

!==============================================
! declare arguments
!==============================================

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: npx, npy, npz
  INTEGER, INTENT(IN) :: k
! number of tracers to be advected
  INTEGER, INTENT(IN) :: nq
  INTEGER, INTENT(IN) :: hord
  INTEGER, INTENT(IN) :: q_split
  INTEGER, INTENT(IN) :: id_divg
  REAL, INTENT(IN) :: dt
  INTEGER, INTENT(IN) :: is, ie, js, je, isd, ied, jsd, jed
! 2D Tracers
  REAL, INTENT(INOUT) :: q(isd:ied, jsd:jed, nq)
  REAL, INTENT(INOUT) :: g_q(isd:ied, jsd:jed, nq)
! Tracers
  REAL, INTENT(INOUT) :: q3(isd:ied, jsd:jed, npz, nq)
  REAL, INTENT(INOUT) :: g_q3(isd:ied, jsd:jed, npz, nq)
! DELP before dyn_core
  REAL, INTENT(INOUT) :: dp0(is:ie, js:je)
  REAL, INTENT(INOUT) :: g_dp0(is:ie, js:je)
! Mass Flux X-Dir
  REAL, INTENT(IN) :: mfx(is:ie+1, js:je)
  REAL, INTENT(IN) :: g_mfx(is:ie+1, js:je)
! Mass Flux Y-Dir
  REAL, INTENT(IN) :: mfy(is:ie, js:je+1)
  REAL, INTENT(IN) :: g_mfy(is:ie, js:je+1)
! Courant Number X-Dir
  REAL, INTENT(IN) :: cx(is:ie+1, jsd:jed)
  REAL, INTENT(IN) :: g_cx(is:ie+1, jsd:jed)
! Courant Number Y-Dir
  REAL, INTENT(IN) :: cy(isd:ied, js:je+1)
  REAL, INTENT(IN) :: g_cy(isd:ied, js:je+1)
  REAL, INTENT(IN) :: dx(is:ie+1, jsd:jed), dy(isd:ied, js:je+1)
  REAL, INTENT(IN) :: dxa(is:ie+1, jsd:jed), dya(isd:ied, js:je+1)
  REAL, INTENT(IN) :: sina_u(is:ie+1, jsd:jed), sina_v(isd:ied, js:je+1)
  REAL, INTENT(IN) :: area(is:ie, jsd:jed), rarea(is:ie, js:je)
  logical, intent(in) :: uniform_ppm
! Local Arrays
  REAL :: mfx2(is:ie+1, js:je)
  REAL :: g_mfx2(is:ie+1, js:je)
  REAL :: mfy2(is:ie, js:je+1)
  REAL :: g_mfy2(is:ie, js:je+1)
  REAL :: cx2(is:ie+1, jsd:jed)
  REAL :: g_cx2(is:ie+1, jsd:jed)
  REAL :: cy2(isd:ied, js:je+1)
  REAL :: g_cy2(isd:ied, js:je+1)
  REAL :: dp1(is:ie, js:je)
  REAL :: g_dp1(is:ie, js:je)
  REAL :: dp2(is:ie, js:je)
  REAL :: g_dp2(is:ie, js:je)
  REAL :: fx(is:ie+1, js:je)
  REAL :: g_fx(is:ie+1, js:je)
  REAL :: fy(is:ie, js:je+1)
  REAL :: g_fy(is:ie, js:je+1)
  REAL :: ra_x(is:ie, jsd:jed)
  REAL :: g_ra_x(is:ie, jsd:jed)
  REAL :: ra_y(isd:ied, js:je)
  REAL :: g_ra_y(isd:ied, js:je)
  REAL :: xfx(is:ie+1, jsd:jed)
  REAL :: g_xfx(is:ie+1, jsd:jed)
  REAL :: yfx(isd:ied, js:je+1)
  REAL :: g_yfx(isd:ied, js:je+1)
  REAL :: cmax
  REAL :: frac, rdt
  INTEGER :: nsplt
  INTEGER :: i, j, it, iq
  integer :: iq_taf
  INTRINSIC ABS
  INTRINSIC MAX
  INTRINSIC INT
  INTRINSIC REAL
  REAL :: x1
  REAL :: abs1
  REAL :: abs0
  REAL :: y1
  g_xfx = 0.0
  DO j=jsd,jed
    DO i=is,ie+1
      IF (cx(i, j) .GT. 0.) THEN
        g_xfx(i, j) = dxa(i-1, j)*dy(i, j)*sina_u(i, j)*g_cx(i, j)
        xfx(i, j) = cx(i, j)*dxa(i-1, j)*dy(i, j)*sina_u(i, j)
      ELSE
        g_xfx(i, j) = dxa(i, j)*dy(i, j)*sina_u(i, j)*g_cx(i, j)
        xfx(i, j) = cx(i, j)*dxa(i, j)*dy(i, j)*sina_u(i, j)
      END IF
    END DO
  END DO
  g_yfx = 0.0
  DO j=js,je+1
    DO i=isd,ied
      IF (cy(i, j) .GT. 0.) THEN
        g_yfx(i, j) = dya(i, j-1)*dx(i, j)*sina_v(i, j)*g_cy(i, j)
        yfx(i, j) = cy(i, j)*dya(i, j-1)*dx(i, j)*sina_v(i, j)
      ELSE
        g_yfx(i, j) = dya(i, j)*dx(i, j)*sina_v(i, j)*g_cy(i, j)
        yfx(i, j) = cy(i, j)*dya(i, j)*dx(i, j)*sina_v(i, j)
      END IF
    END DO
  END DO
  IF (q_split .EQ. 0) THEN
! Determine nsplt for tracer advection
    cmax = 0.
    DO j=js,je
      DO i=is,ie
        IF (cx(i, j) .GE. 0.) THEN
          abs0 = cx(i, j)
        ELSE
          abs0 = -cx(i, j)
        END IF
        x1 = abs0 + (1.-sina_u(i, j))
        IF (cy(i, j) .GE. 0.) THEN
          abs1 = cy(i, j)
        ELSE
          abs1 = -cy(i, j)
        END IF
        y1 = abs1 + (1.-sina_v(i, j))
        IF (x1 .LT. y1) THEN
          IF (y1 .LT. cmax) THEN
            cmax = cmax
          ELSE
            cmax = y1
          END IF
        ELSE IF (x1 .LT. cmax) THEN
          cmax = cmax
        ELSE
          cmax = x1
        END IF
      END DO
    END DO
    CALL MP_REDUCE_MAX(cmax)
    nsplt = INT(1.01 + cmax)
  ELSE
    nsplt = q_split
  END IF
  frac = 1./REAL(nsplt)
  g_cx2 = 0.0
  DO j=jsd,jed
    DO i=is,ie+1
      g_cx2(i, j) = frac*g_cx(i, j)
      cx2(i, j) = cx(i, j)*frac
      g_xfx(i, j) = frac*g_xfx(i, j)
      xfx(i, j) = xfx(i, j)*frac
    END DO
  END DO
  g_mfx2 = 0.0
  DO j=js,je
    DO i=is,ie+1
      g_mfx2(i, j) = frac*g_mfx(i, j)
      mfx2(i, j) = mfx(i, j)*frac
    END DO
  END DO
  g_cy2 = 0.0
  DO j=js,je+1
    DO i=isd,ied
      g_cy2(i, j) = frac*g_cy(i, j)
      cy2(i, j) = cy(i, j)*frac
      g_yfx(i, j) = frac*g_yfx(i, j)
      yfx(i, j) = yfx(i, j)*frac
    END DO
  END DO
  g_mfy2 = 0.0
  DO j=js,je+1
    DO i=is,ie
      g_mfy2(i, j) = frac*g_mfy(i, j)
      mfy2(i, j) = mfy(i, j)*frac
    END DO
  END DO
  g_ra_x = 0.0
  DO j=jsd,jed
    DO i=is,ie
      g_ra_x(i, j) = g_xfx(i, j) - g_xfx(i+1, j)
      ra_x(i, j) = area(i, j) + xfx(i, j) - xfx(i+1, j)
    END DO
  END DO
  g_ra_y = 0.0
  DO j=js,je
    DO i=isd,ied
      g_ra_y(i, j) = g_yfx(i, j) - g_yfx(i, j+1)
      ra_y(i, j) = area(i, j) + yfx(i, j) - yfx(i, j+1)
    END DO
  END DO
  g_dp1 = 0.0
  DO j=js,je
    DO i=is,ie
      g_dp1(i, j) = g_dp0(i, j)
      dp1(i, j) = dp0(i, j)
    END DO
  END DO
  g_dp2 = 0.0
  DO it=1,nsplt
    DO j=js,je
      DO i=is,ie
        g_dp2(i, j) = g_dp1(i, j) + rarea(i, j)*(g_mfx2(i, j)-g_mfx2(i+1, j)&
&         +g_mfy2(i, j)-g_mfy2(i, j+1))
        dp2(i, j) = dp1(i, j) + (mfx2(i, j)-mfx2(i+1, j)+mfy2(i, j)-mfy2&
&         (i, j+1))*rarea(i, j)
      END DO
    END DO
!call timing_on('COMM_TOTAL')
!call timing_on('COMM_TRAC')
    call g_mpp_update_domains_dummy4( q,g_q,is,ie,js,je,isd,ied,jsd,jed,1,nq )
!call timing_off('COMM_TRAC')
!call timing_off('COMM_TOTAL')
    DO iq=1,nq
      iq_taf = (k-1)*nq*nsplt_max+(it-1)*nq+iq
      CALL g_v2_fv_tp_2d_trc(q(isd:, jsd:, iq), g_q(isd:, jsd:, iq), cx2, g_cx2, &
&               cy2, g_cy2, npx, npy, hord, fx, g_fx, fy, g_fy, xfx, g_xfx, &
&               yfx, g_yfx, area, ra_x, g_ra_x, ra_y, g_ra_y, mfx=mfx2, &
&               g_mfx=g_mfx2, mfy=mfy2, g_mfy=g_mfy2, .true. ,iq_taf)
      IF (it .EQ. nsplt) THEN
        DO j=js,je
          DO i=is,ie
            g_q3(i, j, k, iq) = ((g_q(i, j, iq)*dp1(i, j)+q(i, j, iq)*g_dp1&
&             (i, j)+rarea(i, j)*(g_fx(i, j)-g_fx(i+1, j)+g_fy(i, j)-g_fy(i&
&             , j+1)))*dp2(i, j)-(q(i, j, iq)*dp1(i, j)+(fx(i, j)-fx(i+1&
&             , j)+fy(i, j)-fy(i, j+1))*rarea(i, j))*g_dp2(i, j))/dp2(i, &
&             j)**2
            q3(i, j, k, iq) = (q(i, j, iq)*dp1(i, j)+(fx(i, j)-fx(i+1, j&
&             )+fy(i, j)-fy(i, j+1))*rarea(i, j))/dp2(i, j)
          END DO
        END DO
      ELSE
        DO j=js,je
          DO i=is,ie
            g_q(i, j, iq) = ((g_q(i, j, iq)*dp1(i, j)+q(i, j, iq)*g_dp1(i, &
&             j)+rarea(i, j)*(g_fx(i, j)-g_fx(i+1, j)+g_fy(i, j)-g_fy(i, j+1&
&             )))*dp2(i, j)-(q(i, j, iq)*dp1(i, j)+(fx(i, j)-fx(i+1, j)+&
&             fy(i, j)-fy(i, j+1))*rarea(i, j))*g_dp2(i, j))/dp2(i, j)**2
            q(i, j, iq) = (q(i, j, iq)*dp1(i, j)+(fx(i, j)-fx(i+1, j)+fy&
&             (i, j)-fy(i, j+1))*rarea(i, j))/dp2(i, j)
          END DO
        END DO
      END IF
    END DO
    IF (it .NE. nsplt) THEN
      DO j=js,je
        DO i=is,ie
          g_dp1(i, j) = g_dp2(i, j)
          dp1(i, j) = dp2(i, j)
        END DO
      END DO
    END IF
  END DO
! nsplt
  IF (id_divg .GT. 0) THEN
    rdt = 1./(frac*dt)
    DO j=js,je
      DO i=is,ie
        g_dp0(i, j) = rarea(i, j)*rdt*(g_xfx(i+1, j)-g_xfx(i, j)+g_yfx(i, j+&
&         1)-g_yfx(i, j))
        dp0(i, j) = (xfx(i+1, j)-xfx(i, j)+yfx(i, j+1)-yfx(i, j))*rarea(&
&         i, j)*rdt
      END DO
    END DO
  END IF
END SUBROUTINE g_v2_tracer_2d_1l
