module fv_tracer2d_adm_mod
      use tp_core_adm_mod,      only: fv_tp_2d, fv_tp_2d_adm
      use fv_grid_tools_mod,      only: area, rarea, dxa, dya, dx, dy
      use fv_grid_utils_mod,      only: sina_u, sina_v
      use fv_mp_mod,          only: gid, domain, mp_reduce_max,   &
                                 ng,isd,ied,jsd,jed,is,js,ie,je
      use mpp_domains_mod, only: mpp_update_domains
      use fv_timing_mod,    only: timing_on, timing_off

implicit none
private

public :: tracer_2d, tracer_2d_1L, tracer_2d_adm, tracer_2d_1L_adm


contains

  SUBROUTINE TRACER_2D_1L_ADM(q, q_ad, dp0, dp0_ad, mfx, mfx_ad, mfy, &
&   mfy_ad, cx, cx_ad, cy, cy_ad, npx, npy, npz, nq, hord, q_split, k, &
&   q3, q3_ad, dt, id_divg)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx, npy, npz
    INTEGER, INTENT(IN) :: k
! number of tracers to be advected
    INTEGER, INTENT(IN) :: nq
    INTEGER, INTENT(IN) :: hord
    INTEGER, INTENT(IN) :: q_split
    INTEGER, INTENT(IN) :: id_divg
    REAL, INTENT(IN) :: dt
! 2D Tracers
    REAL, INTENT(INOUT) :: q(isd:ied, jsd:jed, nq)
    REAL, INTENT(INOUT) :: q_ad(isd:ied, jsd:jed, nq)
! Tracers
    REAL, INTENT(INOUT) :: q3(isd:ied, jsd:jed, npz, nq)
    REAL, INTENT(INOUT) :: q3_ad(isd:ied, jsd:jed, npz, nq)
! DELP before dyn_core
    REAL, INTENT(INOUT) :: dp0(is:ie, js:je)
    REAL, INTENT(INOUT) :: dp0_ad(is:ie, js:je)
! Mass Flux X-Dir
    REAL, INTENT(IN) :: mfx(is:ie+1, js:je)
    REAL :: mfx_ad(is:ie+1, js:je)
! Mass Flux Y-Dir
    REAL, INTENT(IN) :: mfy(is:ie, js:je+1)
    REAL :: mfy_ad(is:ie, js:je+1)
! Courant Number X-Dir
    REAL, INTENT(IN) :: cx(is:ie+1, jsd:jed)
    REAL :: cx_ad(is:ie+1, jsd:jed)
! Courant Number Y-Dir
    REAL, INTENT(IN) :: cy(isd:ied, js:je+1)
    REAL :: cy_ad(isd:ied, js:je+1)
! Local Arrays
    REAL :: mfx2(is:ie+1, js:je)
    REAL :: mfx2_ad(is:ie+1, js:je)
    REAL :: mfy2(is:ie, js:je+1)
    REAL :: mfy2_ad(is:ie, js:je+1)
    REAL :: cx2(is:ie+1, jsd:jed)
    REAL :: cx2_ad(is:ie+1, jsd:jed)
    REAL :: cy2(isd:ied, js:je+1)
    REAL :: cy2_ad(isd:ied, js:je+1)
    REAL :: dp1(is:ie, js:je)
    REAL :: dp1_ad(is:ie, js:je)
    REAL :: dp2(is:ie, js:je)
    REAL :: dp2_ad(is:ie, js:je)
    REAL :: fx(is:ie+1, js:je)
    REAL :: fx_ad(is:ie+1, js:je)
    REAL :: fy(is:ie, js:je+1)
    REAL :: fy_ad(is:ie, js:je+1)
    REAL :: ra_x(is:ie, jsd:jed)
    REAL :: ra_x_ad(is:ie, jsd:jed)
    REAL :: ra_y(isd:ied, js:je)
    REAL :: ra_y_ad(isd:ied, js:je)
    REAL :: xfx(is:ie+1, jsd:jed)
    REAL :: xfx_ad(is:ie+1, jsd:jed)
    REAL :: yfx(isd:ied, js:je+1)
    REAL :: yfx_ad(isd:ied, js:je+1)
    REAL :: cmax
    REAL :: frac, rdt
    INTEGER :: nsplt
    INTEGER :: i, j, it, iq
    INTRINSIC ABS
    INTRINSIC MAX
    INTRINSIC INT
    INTRINSIC REAL
    INTEGER :: branch
    REAL :: x1
    REAL :: temp_ad4
    REAL :: temp_ad3
    REAL :: temp_ad2
    REAL :: temp_ad1
    REAL :: temp_ad0
    REAL :: temp_ad
    REAL :: abs1
    REAL :: abs0
    REAL :: y1
    DO j=jsd,jed
      DO i=is,ie+1
        IF (cx(i, j) .GT. 0.) THEN
          xfx(i, j) = cx(i, j)*dxa(i-1, j)*dy(i, j)*sina_u(i, j)
          CALL PUSHCONTROL1B(1)
        ELSE
          xfx(i, j) = cx(i, j)*dxa(i, j)*dy(i, j)*sina_u(i, j)
          CALL PUSHCONTROL1B(0)
        END IF
      END DO
    END DO
    DO j=js,je+1
      DO i=isd,ied
        IF (cy(i, j) .GT. 0.) THEN
          yfx(i, j) = cy(i, j)*dya(i, j-1)*dx(i, j)*sina_v(i, j)
          CALL PUSHCONTROL1B(1)
        ELSE
          yfx(i, j) = cy(i, j)*dya(i, j)*dx(i, j)*sina_v(i, j)
          CALL PUSHCONTROL1B(0)
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
              CALL PUSHCONTROL2B(0)
              cmax = cmax
            ELSE
              CALL PUSHCONTROL2B(1)
              cmax = y1
            END IF
          ELSE IF (x1 .LT. cmax) THEN
            CALL PUSHCONTROL2B(2)
            cmax = cmax
          ELSE
            CALL PUSHCONTROL2B(3)
            cmax = x1
          END IF
        END DO
      END DO
      CALL PUSHCONTROL1B(0)
!oncall mp_reduce_max(cmax)
      call mp_reduce_max(cmax)
      nsplt = INT(1.01 + cmax)
!        if ( gid == 0 .and. nsplt > 5 )  write(6,*) k, 'Tracer_2d_split=', nsplt, cmax
    ELSE
      CALL PUSHCONTROL1B(1)
      nsplt = q_split
    END IF
    frac = 1./REAL(nsplt)
    DO j=jsd,jed
      DO i=is,ie+1
        cx2(i, j) = cx(i, j)*frac
        xfx(i, j) = xfx(i, j)*frac
      END DO
    END DO
    DO j=js,je
      DO i=is,ie+1
        mfx2(i, j) = mfx(i, j)*frac
      END DO
    END DO
    DO j=js,je+1
      DO i=isd,ied
        cy2(i, j) = cy(i, j)*frac
        yfx(i, j) = yfx(i, j)*frac
      END DO
    END DO
    DO j=js,je+1
      DO i=is,ie
        mfy2(i, j) = mfy(i, j)*frac
      END DO
    END DO
    DO j=jsd,jed
      DO i=is,ie
        ra_x(i, j) = area(i, j) + xfx(i, j) - xfx(i+1, j)
      END DO
    END DO
    DO j=js,je
      DO i=isd,ied
        ra_y(i, j) = area(i, j) + yfx(i, j) - yfx(i, j+1)
      END DO
    END DO
    DO j=js,je
      DO i=is,ie
        dp1(i, j) = dp0(i, j)
      END DO
    END DO
    DO it=1,nsplt
      DO j=js,je
        DO i=is,ie
          CALL PUSHREAL8(dp2(i, j))
          dp2(i, j) = dp1(i, j) + (mfx2(i, j)-mfx2(i+1, j)+mfy2(i, j)-&
&           mfy2(i, j+1))*rarea(i, j)
        END DO
      END DO
!call timing_on('COMM_TOTAL')
!     call timing_on('COMM_TRAC')
!oncall mpp_update_domains( q, domain, complete= .true. )
      CALL PUSHREAL8ARRAY(q, (ied-isd+1)*(jed-jsd+1)*nq)
      call mpp_update_domains( q, domain, complete= .true. )
!     call timing_off('COMM_TRAC')
!call timing_off('COMM_TOTAL')
      DO iq=1,nq
        CALL PUSHREAL8ARRAY(fy, (ie-is+1)*(je-js+2))
        CALL PUSHREAL8ARRAY(fx, (ie-is+2)*(je-js+1))
        CALL PUSHREAL8ARRAY(q(:, :, iq), (ied-isd+1)*(jed-jsd+1))
        CALL FV_TP_2D(q(isd:, jsd:, iq), cx2, cy2, npx, npy, hord, fx, &
&               fy, xfx, yfx, area, ra_x, ra_y, mfx=mfx2, mfy=mfy2)
        IF (it .EQ. nsplt) THEN
          CALL PUSHCONTROL1B(1)
        ELSE
          DO j=js,je
            DO i=is,ie
              CALL PUSHREAL8(q(i, j, iq))
              q(i, j, iq) = (q(i, j, iq)*dp1(i, j)+(fx(i, j)-fx(i+1, j)+&
&               fy(i, j)-fy(i, j+1))*rarea(i, j))/dp2(i, j)
            END DO
          END DO
          CALL PUSHCONTROL1B(0)
        END IF
      END DO
      IF (it .NE. nsplt) THEN
        DO j=js,je
          DO i=is,ie
            CALL PUSHREAL8(dp1(i, j))
            dp1(i, j) = dp2(i, j)
          END DO
        END DO
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL PUSHCONTROL1B(0)
      END IF
    END DO
! nsplt
    IF (id_divg .GT. 0) THEN
      rdt = 1./(frac*dt)
      xfx_ad = 0.0_8
      yfx_ad = 0.0_8
      DO j=je,js,-1
        DO i=ie,is,-1
          temp_ad4 = rarea(i, j)*rdt*dp0_ad(i, j)
          xfx_ad(i+1, j) = xfx_ad(i+1, j) + temp_ad4
          xfx_ad(i, j) = xfx_ad(i, j) - temp_ad4
          yfx_ad(i, j+1) = yfx_ad(i, j+1) + temp_ad4
          yfx_ad(i, j) = yfx_ad(i, j) - temp_ad4
          dp0_ad(i, j) = 0.0_8
        END DO
      END DO
    ELSE
      xfx_ad = 0.0_8
      yfx_ad = 0.0_8
    END IF
    cx2_ad = 0.0_8
    mfx2_ad = 0.0_8
    dp1_ad = 0.0_8
    dp2_ad = 0.0_8
    cy2_ad = 0.0_8
    mfy2_ad = 0.0_8
    ra_x_ad = 0.0_8
    ra_y_ad = 0.0_8
    fx_ad = 0.0_8
    fy_ad = 0.0_8
    DO it=nsplt,1,-1
      CALL POPCONTROL1B(branch)
      IF (branch .NE. 0) THEN
        DO j=je,js,-1
          DO i=ie,is,-1
            CALL POPREAL8(dp1(i, j))
            dp2_ad(i, j) = dp2_ad(i, j) + dp1_ad(i, j)
            dp1_ad(i, j) = 0.0_8
          END DO
        END DO
      END IF
      DO iq=nq,1,-1
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          DO j=je,js,-1
            DO i=ie,is,-1
              CALL POPREAL8(q(i, j, iq))
              temp_ad2 = q_ad(i, j, iq)/dp2(i, j)
              temp_ad3 = rarea(i, j)*temp_ad2
              dp1_ad(i, j) = dp1_ad(i, j) + q(i, j, iq)*temp_ad2
              fx_ad(i, j) = fx_ad(i, j) + temp_ad3
              fx_ad(i+1, j) = fx_ad(i+1, j) - temp_ad3
              fy_ad(i, j) = fy_ad(i, j) + temp_ad3
              fy_ad(i, j+1) = fy_ad(i, j+1) - temp_ad3
              dp2_ad(i, j) = dp2_ad(i, j) - (q(i, j, iq)*dp1(i, j)+rarea&
&               (i, j)*(fx(i, j)-fx(i+1, j)+fy(i, j)-fy(i, j+1)))*&
&               temp_ad2/dp2(i, j)
              q_ad(i, j, iq) = dp1(i, j)*temp_ad2
            END DO
          END DO
        ELSE
          DO j=je,js,-1
            DO i=ie,is,-1
              temp_ad0 = q3_ad(i, j, k, iq)/dp2(i, j)
              temp_ad1 = rarea(i, j)*temp_ad0
              q_ad(i, j, iq) = q_ad(i, j, iq) + dp1(i, j)*temp_ad0
              dp1_ad(i, j) = dp1_ad(i, j) + q(i, j, iq)*temp_ad0
              fx_ad(i, j) = fx_ad(i, j) + temp_ad1
              fx_ad(i+1, j) = fx_ad(i+1, j) - temp_ad1
              fy_ad(i, j) = fy_ad(i, j) + temp_ad1
              fy_ad(i, j+1) = fy_ad(i, j+1) - temp_ad1
              dp2_ad(i, j) = dp2_ad(i, j) - (q(i, j, iq)*dp1(i, j)+rarea&
&               (i, j)*(fx(i, j)-fx(i+1, j)+fy(i, j)-fy(i, j+1)))*&
&               temp_ad0/dp2(i, j)
              q3_ad(i, j, k, iq) = 0.0_8
            END DO
          END DO
        END IF
        CALL POPREAL8ARRAY(q(:, :, iq), (ied-isd+1)*(jed-jsd+1))
        CALL POPREAL8ARRAY(fx, (ie-is+2)*(je-js+1))
        CALL POPREAL8ARRAY(fy, (ie-is+1)*(je-js+2))
        CALL FV_TP_2D_ADM(q(isd:, jsd:, iq), q_ad(isd:, jsd:, iq), cx2, &
&                   cx2_ad, cy2, cy2_ad, npx, npy, hord, fx, fx_ad, fy, &
&                   fy_ad, xfx, xfx_ad, yfx, yfx_ad, area, ra_x, ra_x_ad&
&                   , ra_y, ra_y_ad, mfx=mfx2, mfx_ad=mfx2_ad, mfy=mfy2&
&                   , mfy_ad=mfy2_ad)
      END DO
      CALL POPREAL8ARRAY(q, (ied-isd+1)*(jed-jsd+1)*nq)
      call mpp_update_domains( q_ad, domain, complete= .true. )
      DO j=je,js,-1
        DO i=ie,is,-1
          CALL POPREAL8(dp2(i, j))
          temp_ad = rarea(i, j)*dp2_ad(i, j)
          dp1_ad(i, j) = dp1_ad(i, j) + dp2_ad(i, j)
          mfx2_ad(i, j) = mfx2_ad(i, j) + temp_ad
          mfx2_ad(i+1, j) = mfx2_ad(i+1, j) - temp_ad
          mfy2_ad(i, j) = mfy2_ad(i, j) + temp_ad
          mfy2_ad(i, j+1) = mfy2_ad(i, j+1) - temp_ad
          dp2_ad(i, j) = 0.0_8
        END DO
      END DO
    END DO
    DO j=je,js,-1
      DO i=ie,is,-1
        dp0_ad(i, j) = dp0_ad(i, j) + dp1_ad(i, j)
        dp1_ad(i, j) = 0.0_8
      END DO
    END DO
    DO j=je,js,-1
      DO i=ied,isd,-1
        yfx_ad(i, j) = yfx_ad(i, j) + ra_y_ad(i, j)
        yfx_ad(i, j+1) = yfx_ad(i, j+1) - ra_y_ad(i, j)
        ra_y_ad(i, j) = 0.0_8
      END DO
    END DO
    DO j=jed,jsd,-1
      DO i=ie,is,-1
        xfx_ad(i, j) = xfx_ad(i, j) + ra_x_ad(i, j)
        xfx_ad(i+1, j) = xfx_ad(i+1, j) - ra_x_ad(i, j)
        ra_x_ad(i, j) = 0.0_8
      END DO
    END DO
    DO j=je+1,js,-1
      DO i=ie,is,-1
        mfy_ad(i, j) = mfy_ad(i, j) + frac*mfy2_ad(i, j)
        mfy2_ad(i, j) = 0.0_8
      END DO
    END DO
    DO j=je+1,js,-1
      DO i=ied,isd,-1
        yfx_ad(i, j) = frac*yfx_ad(i, j)
        cy_ad(i, j) = cy_ad(i, j) + frac*cy2_ad(i, j)
        cy2_ad(i, j) = 0.0_8
      END DO
    END DO
    DO j=je,js,-1
      DO i=ie+1,is,-1
        mfx_ad(i, j) = mfx_ad(i, j) + frac*mfx2_ad(i, j)
        mfx2_ad(i, j) = 0.0_8
      END DO
    END DO
    DO j=jed,jsd,-1
      DO i=ie+1,is,-1
        xfx_ad(i, j) = frac*xfx_ad(i, j)
        cx_ad(i, j) = cx_ad(i, j) + frac*cx2_ad(i, j)
        cx2_ad(i, j) = 0.0_8
      END DO
    END DO
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      DO j=je,js,-1
        DO i=ie,is,-1
          CALL POPCONTROL2B(branch)
        END DO
      END DO
    END IF
    DO j=je+1,js,-1
      DO i=ied,isd,-1
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          cy_ad(i, j) = cy_ad(i, j) + dya(i, j)*dx(i, j)*sina_v(i, j)*&
&           yfx_ad(i, j)
          yfx_ad(i, j) = 0.0_8
        ELSE
          cy_ad(i, j) = cy_ad(i, j) + dx(i, j)*sina_v(i, j)*dya(i, j-1)*&
&           yfx_ad(i, j)
          yfx_ad(i, j) = 0.0_8
        END IF
      END DO
    END DO
    DO j=jed,jsd,-1
      DO i=ie+1,is,-1
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          cx_ad(i, j) = cx_ad(i, j) + dxa(i, j)*dy(i, j)*sina_u(i, j)*&
&           xfx_ad(i, j)
          xfx_ad(i, j) = 0.0_8
        ELSE
          cx_ad(i, j) = cx_ad(i, j) + dy(i, j)*sina_u(i, j)*dxa(i-1, j)*&
&           xfx_ad(i, j)
          xfx_ad(i, j) = 0.0_8
        END IF
      END DO
    END DO
  END SUBROUTINE TRACER_2D_1L_ADM
!-----------------------------------------------------------------------
! !ROUTINE: Perform 2D horizontal-to-lagrangian transport
!-----------------------------------------------------------------------
  SUBROUTINE TRACER_2D_1L(q, dp0, mfx, mfy, cx, cy, npx, npy, npz, nq, &
&   hord, q_split, k, q3, dt, id_divg)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx, npy, npz
    INTEGER, INTENT(IN) :: k
! number of tracers to be advected
    INTEGER, INTENT(IN) :: nq
    INTEGER, INTENT(IN) :: hord
    INTEGER, INTENT(IN) :: q_split
    INTEGER, INTENT(IN) :: id_divg
    REAL, INTENT(IN) :: dt
! 2D Tracers
    REAL, INTENT(INOUT) :: q(isd:ied, jsd:jed, nq)
! Tracers
    REAL, INTENT(INOUT) :: q3(isd:ied, jsd:jed, npz, nq)
! DELP before dyn_core
    REAL, INTENT(INOUT) :: dp0(is:ie, js:je)
! Mass Flux X-Dir
    REAL, INTENT(IN) :: mfx(is:ie+1, js:je)
! Mass Flux Y-Dir
    REAL, INTENT(IN) :: mfy(is:ie, js:je+1)
! Courant Number X-Dir
    REAL, INTENT(IN) :: cx(is:ie+1, jsd:jed)
! Courant Number Y-Dir
    REAL, INTENT(IN) :: cy(isd:ied, js:je+1)
! Local Arrays
    REAL :: mfx2(is:ie+1, js:je)
    REAL :: mfy2(is:ie, js:je+1)
    REAL :: cx2(is:ie+1, jsd:jed)
    REAL :: cy2(isd:ied, js:je+1)
    REAL :: dp1(is:ie, js:je)
    REAL :: dp2(is:ie, js:je)
    REAL :: fx(is:ie+1, js:je)
    REAL :: fy(is:ie, js:je+1)
    REAL :: ra_x(is:ie, jsd:jed)
    REAL :: ra_y(isd:ied, js:je)
    REAL :: xfx(is:ie+1, jsd:jed)
    REAL :: yfx(isd:ied, js:je+1)
    REAL :: cmax
    REAL :: frac, rdt
    INTEGER :: nsplt
    INTEGER :: i, j, it, iq
    INTRINSIC ABS
    INTRINSIC MAX
    INTRINSIC INT
    INTRINSIC REAL
    REAL :: x1
    REAL :: abs1
    REAL :: abs0
    REAL :: y1
    DO j=jsd,jed
      DO i=is,ie+1
        IF (cx(i, j) .GT. 0.) THEN
          xfx(i, j) = cx(i, j)*dxa(i-1, j)*dy(i, j)*sina_u(i, j)
        ELSE
          xfx(i, j) = cx(i, j)*dxa(i, j)*dy(i, j)*sina_u(i, j)
        END IF
      END DO
    END DO
    DO j=js,je+1
      DO i=isd,ied
        IF (cy(i, j) .GT. 0.) THEN
          yfx(i, j) = cy(i, j)*dya(i, j-1)*dx(i, j)*sina_v(i, j)
        ELSE
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
!oncall mp_reduce_max(cmax)
      call mp_reduce_max(cmax)
      nsplt = INT(1.01 + cmax)
!        if ( gid == 0 .and. nsplt > 5 )  write(6,*) k, 'Tracer_2d_split=', nsplt, cmax
    ELSE
      nsplt = q_split
    END IF
    frac = 1./REAL(nsplt)
    DO j=jsd,jed
      DO i=is,ie+1
        cx2(i, j) = cx(i, j)*frac
        xfx(i, j) = xfx(i, j)*frac
      END DO
    END DO
    DO j=js,je
      DO i=is,ie+1
        mfx2(i, j) = mfx(i, j)*frac
      END DO
    END DO
    DO j=js,je+1
      DO i=isd,ied
        cy2(i, j) = cy(i, j)*frac
        yfx(i, j) = yfx(i, j)*frac
      END DO
    END DO
    DO j=js,je+1
      DO i=is,ie
        mfy2(i, j) = mfy(i, j)*frac
      END DO
    END DO
    DO j=jsd,jed
      DO i=is,ie
        ra_x(i, j) = area(i, j) + xfx(i, j) - xfx(i+1, j)
      END DO
    END DO
    DO j=js,je
      DO i=isd,ied
        ra_y(i, j) = area(i, j) + yfx(i, j) - yfx(i, j+1)
      END DO
    END DO
    DO j=js,je
      DO i=is,ie
        dp1(i, j) = dp0(i, j)
      END DO
    END DO
    DO it=1,nsplt
      DO j=js,je
        DO i=is,ie
          dp2(i, j) = dp1(i, j) + (mfx2(i, j)-mfx2(i+1, j)+mfy2(i, j)-&
&           mfy2(i, j+1))*rarea(i, j)
        END DO
      END DO
!call timing_on('COMM_TOTAL')
!     call timing_on('COMM_TRAC')
!oncall mpp_update_domains( q, domain, complete= .true. )
      call mpp_update_domains( q, domain, complete= .true. )
!     call timing_off('COMM_TRAC')
!call timing_off('COMM_TOTAL')
      DO iq=1,nq
        CALL FV_TP_2D(q(isd:, jsd:, iq), cx2, cy2, npx, npy, hord, fx, &
&               fy, xfx, yfx, area, ra_x, ra_y, mfx=mfx2, mfy=mfy2)
        IF (it .EQ. nsplt) THEN
          DO j=js,je
            DO i=is,ie
              q3(i, j, k, iq) = (q(i, j, iq)*dp1(i, j)+(fx(i, j)-fx(i+1&
&               , j)+fy(i, j)-fy(i, j+1))*rarea(i, j))/dp2(i, j)
            END DO
          END DO
        ELSE
          DO j=js,je
            DO i=is,ie
              q(i, j, iq) = (q(i, j, iq)*dp1(i, j)+(fx(i, j)-fx(i+1, j)+&
&               fy(i, j)-fy(i, j+1))*rarea(i, j))/dp2(i, j)
            END DO
          END DO
        END IF
      END DO
      IF (it .NE. nsplt) THEN
        DO j=js,je
          DO i=is,ie
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
          dp0(i, j) = (xfx(i+1, j)-xfx(i, j)+yfx(i, j+1)-yfx(i, j))*&
&           rarea(i, j)*rdt
        END DO
      END DO
    END IF
  END SUBROUTINE TRACER_2D_1L
!  Differentiation of tracer_2d in reverse (adjoint) mode (with options r8):
!   gradient     of useful results: q dp1 mfx mfy cx cy
!   with respect to varying inputs: q dp1 mfx mfy cx cy
  SUBROUTINE TRACER_2D_ADM(q, q_ad, dp1, dp1_ad, mfx, mfx_ad, mfy, &
&   mfy_ad, cx, cx_ad, cy, cy_ad, npx, npy, npz, nq, hord, q_split, dt, &
&   id_divg)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx
    INTEGER, INTENT(IN) :: npy
    INTEGER, INTENT(IN) :: npz
! number of tracers to be advected
    INTEGER, INTENT(IN) :: nq
    INTEGER, INTENT(IN) :: hord
    INTEGER, INTENT(IN) :: q_split
    INTEGER, INTENT(IN) :: id_divg
    REAL, INTENT(IN) :: dt
! Tracers
    REAL, INTENT(INOUT) :: q(isd:ied, jsd:jed, npz, nq)
    REAL, INTENT(INOUT) :: q_ad(isd:ied, jsd:jed, npz, nq)
! DELP before dyn_core
    REAL, INTENT(INOUT) :: dp1(is:ie, js:je, npz)
    REAL, INTENT(INOUT) :: dp1_ad(is:ie, js:je, npz)
! Mass Flux X-Dir
    REAL, INTENT(INOUT) :: mfx(is:ie+1, js:je, npz)
    REAL, INTENT(INOUT) :: mfx_ad(is:ie+1, js:je, npz)
! Mass Flux Y-Dir
    REAL, INTENT(INOUT) :: mfy(is:ie, js:je+1, npz)
    REAL, INTENT(INOUT) :: mfy_ad(is:ie, js:je+1, npz)
! Courant Number X-Dir
    REAL, INTENT(INOUT) :: cx(is:ie+1, jsd:jed, npz)
    REAL, INTENT(INOUT) :: cx_ad(is:ie+1, jsd:jed, npz)
! Courant Number Y-Dir
    REAL, INTENT(INOUT) :: cy(isd:ied, js:je+1, npz)
    REAL, INTENT(INOUT) :: cy_ad(isd:ied, js:je+1, npz)
! Local Arrays
    REAL :: dp2(is:ie, js:je)
    REAL :: dp2_ad(is:ie, js:je)
    REAL :: fx(is:ie+1, js:je)
    REAL :: fx_ad(is:ie+1, js:je)
    REAL :: fy(is:ie, js:je+1)
    REAL :: fy_ad(is:ie, js:je+1)
    REAL :: ra_x(is:ie, jsd:jed)
    REAL :: ra_x_ad(is:ie, jsd:jed)
    REAL :: ra_y(isd:ied, js:je)
    REAL :: ra_y_ad(isd:ied, js:je)
    REAL :: xfx(is:ie+1, jsd:jed, npz)
    REAL :: xfx_ad(is:ie+1, jsd:jed, npz)
    REAL :: yfx(isd:ied, js:je+1, npz)
    REAL :: yfx_ad(isd:ied, js:je+1, npz)
    REAL :: cmax(npz)
    REAL :: c_global
    REAL :: frac, rdt
    INTEGER :: nsplt
    INTEGER :: i, j, k, it, iq
    INTRINSIC ABS
    INTRINSIC MAX
    INTRINSIC INT
    INTRINSIC REAL
    INTEGER :: branch
    REAL :: x1
    REAL :: temp_ad2
    REAL :: temp_ad1
    REAL :: temp_ad0
    REAL :: temp_ad
    REAL :: abs1
    REAL :: abs0
    REAL :: temp
    REAL :: y1
    DO k=1,npz
      DO j=jsd,jed
        DO i=is,ie+1
          IF (cx(i, j, k) .GT. 0.) THEN
            xfx(i, j, k) = cx(i, j, k)*dxa(i-1, j)*dy(i, j)*sina_u(i, j)
            CALL PUSHCONTROL1B(1)
          ELSE
            xfx(i, j, k) = cx(i, j, k)*dxa(i, j)*dy(i, j)*sina_u(i, j)
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
      END DO
      DO j=js,je+1
        DO i=isd,ied
          IF (cy(i, j, k) .GT. 0.) THEN
            yfx(i, j, k) = cy(i, j, k)*dya(i, j-1)*dx(i, j)*sina_v(i, j)
            CALL PUSHCONTROL1B(1)
          ELSE
            yfx(i, j, k) = cy(i, j, k)*dya(i, j)*dx(i, j)*sina_v(i, j)
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
      END DO
    END DO
!--------------------------------------------------------------------------------
    IF (q_split .EQ. 0) THEN
! Determine nsplt
      DO k=1,npz
        cmax(k) = 0.
        DO j=js,je
          DO i=is,ie
            IF (cx(i, j, k) .GE. 0.) THEN
              abs0 = cx(i, j, k)
            ELSE
              abs0 = -cx(i, j, k)
            END IF
            x1 = abs0 + 1. - sina_u(i, j)
            IF (cy(i, j, k) .GE. 0.) THEN
              abs1 = cy(i, j, k)
            ELSE
              abs1 = -cy(i, j, k)
            END IF
            y1 = abs1 + 1. - sina_v(i, j)
            IF (x1 .LT. y1) THEN
              IF (y1 .LT. cmax(k)) THEN
                CALL PUSHCONTROL2B(0)
                cmax(k) = cmax(k)
              ELSE
                CALL PUSHCONTROL2B(1)
                cmax(k) = y1
              END IF
            ELSE IF (x1 .LT. cmax(k)) THEN
              CALL PUSHCONTROL2B(2)
              cmax(k) = cmax(k)
            ELSE
              CALL PUSHCONTROL2B(3)
              cmax(k) = x1
            END IF
          END DO
        END DO
      END DO
!oncall mp_reduce_max(cmax,npz)
      call mp_reduce_max(cmax,npz)
! find global max courant number and define nsplt to scale cx,cy,mfx,mfy
      c_global = cmax(1)
      IF (npz .NE. 1) THEN
! if NOT shallow water test case
        DO k=2,npz
          IF (cmax(k) .LT. c_global) THEN
            CALL PUSHCONTROL1B(0)
            c_global = c_global
          ELSE
            CALL PUSHCONTROL1B(1)
            c_global = cmax(k)
          END IF
        END DO
        CALL PUSHCONTROL2B(0)
      ELSE
        CALL PUSHCONTROL2B(1)
      END IF
      nsplt = INT(1. + c_global)
!     if ( gid == 0 .and. nsplt > 5 )  write(6,*) 'Tracer_2d_split=', nsplt, c_global
    ELSE
      CALL PUSHCONTROL2B(2)
      nsplt = q_split
    END IF
!--------------------------------------------------------------------------------
    frac = 1./REAL(nsplt)
    IF (nsplt .NE. 1) THEN
      DO k=1,npz
        DO j=jsd,jed
          DO i=is,ie+1
            cx(i, j, k) = cx(i, j, k)*frac
            xfx(i, j, k) = xfx(i, j, k)*frac
          END DO
        END DO
        DO j=js,je
          DO i=is,ie+1
            mfx(i, j, k) = mfx(i, j, k)*frac
          END DO
        END DO
        DO j=js,je+1
          DO i=isd,ied
            cy(i, j, k) = cy(i, j, k)*frac
            yfx(i, j, k) = yfx(i, j, k)*frac
          END DO
        END DO
        DO j=js,je+1
          DO i=is,ie
            mfy(i, j, k) = mfy(i, j, k)*frac
          END DO
        END DO
      END DO
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHCONTROL1B(0)
    END IF
    DO it=1,nsplt
!     call timing_on('COMM_TOTAL')
!       call timing_on('COMM_TRAC')
!call mpp_update_domains( q, domain, complete=.true. )
      call mpp_update_domains( q, domain, complete=.true. )
!              call timing_off('COMM_TRAC')
!     call timing_off('COMM_TOTAL')
      DO k=1,npz
        DO j=jsd,jed
          DO i=is,ie
            CALL PUSHREAL8(ra_x(i, j))
            ra_x(i, j) = area(i, j) + xfx(i, j, k) - xfx(i+1, j, k)
          END DO
        END DO
        DO j=js,je
          DO i=isd,ied
            CALL PUSHREAL8(ra_y(i, j))
            ra_y(i, j) = area(i, j) + yfx(i, j, k) - yfx(i, j+1, k)
          END DO
        END DO
        DO j=js,je
          DO i=is,ie
            CALL PUSHREAL8(dp2(i, j))
            dp2(i, j) = dp1(i, j, k) + (mfx(i, j, k)-mfx(i+1, j, k)+mfy(&
&             i, j, k)-mfy(i, j+1, k))*rarea(i, j)
          END DO
        END DO
        DO iq=1,nq
          CALL PUSHREAL8ARRAY(fy, (ie-is+1)*(je-js+2))
          CALL PUSHREAL8ARRAY(fx, (ie-is+2)*(je-js+1))
          CALL PUSHREAL8ARRAY(q(:, :, k, iq), (ied-isd+1)*(jed-jsd+1))
          CALL FV_TP_2D(q(isd:, jsd:, k, iq), cx(is:, jsd:, k), cy(isd:&
&                 , js:, k), npx, npy, hord, fx, fy, xfx(is:, jsd:, k), &
&                 yfx(isd:, js:, k), area, ra_x, ra_y, mfx=mfx(is:, js:&
&                 , k), mfy=mfy(is:, js:, k))
          DO j=js,je
            DO i=is,ie
              CALL PUSHREAL8(q(i, j, k, iq))
              q(i, j, k, iq) = (q(i, j, k, iq)*dp1(i, j, k)+(fx(i, j)-fx&
&               (i+1, j)+fy(i, j)-fy(i, j+1))*rarea(i, j))/dp2(i, j)
            END DO
          END DO
        END DO
        DO j=js,je
          DO i=is,ie
            CALL PUSHREAL8(dp1(i, j, k))
            dp1(i, j, k) = dp2(i, j)
          END DO
        END DO
      END DO
    END DO
! npz
! nsplt
    IF (id_divg .GT. 0) THEN
      rdt = 1./(frac*dt)
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
! Rescale Mass fluxes for exports to offline tracer transport
    IF (nsplt .NE. 1) THEN
      DO k=npz,1,-1
        DO j=je+1,js,-1
          DO i=ie,is,-1
            mfy_ad(i, j, k) = mfy_ad(i, j, k)/frac
          END DO
        END DO
        DO j=je,js,-1
          DO i=ie+1,is,-1
            mfx_ad(i, j, k) = mfx_ad(i, j, k)/frac
          END DO
        END DO
      END DO
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      xfx_ad = 0.0_8
      yfx_ad = 0.0_8
      DO k=npz,1,-1
        DO j=je,js,-1
          DO i=ie,is,-1
            temp_ad2 = rarea(i, j)*rdt*dp1_ad(i, j, k)
            xfx_ad(i+1, j, k) = xfx_ad(i+1, j, k) + temp_ad2
            xfx_ad(i, j, k) = xfx_ad(i, j, k) - temp_ad2
            yfx_ad(i, j+1, k) = yfx_ad(i, j+1, k) + temp_ad2
            yfx_ad(i, j, k) = yfx_ad(i, j, k) - temp_ad2
            dp1_ad(i, j, k) = 0.0_8
          END DO
        END DO
      END DO
    ELSE
      xfx_ad = 0.0_8
      yfx_ad = 0.0_8
    END IF
    dp2_ad = 0.0_8
    ra_x_ad = 0.0_8
    ra_y_ad = 0.0_8
    fx_ad = 0.0_8
    fy_ad = 0.0_8
    DO it=nsplt,1,-1
      DO k=npz,1,-1
        DO j=je,js,-1
          DO i=ie,is,-1
            CALL POPREAL8(dp1(i, j, k))
            dp2_ad(i, j) = dp2_ad(i, j) + dp1_ad(i, j, k)
            dp1_ad(i, j, k) = 0.0_8
          END DO
        END DO
        DO iq=nq,1,-1
          DO j=je,js,-1
            DO i=ie,is,-1
              CALL POPREAL8(q(i, j, k, iq))
              temp_ad0 = q_ad(i, j, k, iq)/dp2(i, j)
              temp = q(i, j, k, iq)
              temp_ad1 = rarea(i, j)*temp_ad0
              dp1_ad(i, j, k) = dp1_ad(i, j, k) + temp*temp_ad0
              fx_ad(i, j) = fx_ad(i, j) + temp_ad1
              fx_ad(i+1, j) = fx_ad(i+1, j) - temp_ad1
              fy_ad(i, j) = fy_ad(i, j) + temp_ad1
              fy_ad(i, j+1) = fy_ad(i, j+1) - temp_ad1
              dp2_ad(i, j) = dp2_ad(i, j) - (temp*dp1(i, j, k)+rarea(i, &
&               j)*(fx(i, j)-fx(i+1, j)+fy(i, j)-fy(i, j+1)))*temp_ad0/&
&               dp2(i, j)
              q_ad(i, j, k, iq) = dp1(i, j, k)*temp_ad0
            END DO
          END DO
          CALL POPREAL8ARRAY(q(:, :, k, iq), (ied-isd+1)*(jed-jsd+1))
          CALL POPREAL8ARRAY(fx, (ie-is+2)*(je-js+1))
          CALL POPREAL8ARRAY(fy, (ie-is+1)*(je-js+2))
          CALL FV_TP_2D_ADM(q(isd:, jsd:, k, iq), q_ad(isd:, jsd:, k, iq&
&                     ), cx(is:, jsd:, k), cx_ad(is:, jsd:, k), cy(isd:&
&                     , js:, k), cy_ad(isd:, js:, k), npx, npy, hord, fx&
&                     , fx_ad, fy, fy_ad, xfx(is:, jsd:, k), xfx_ad(is:&
&                     , jsd:, k), yfx(isd:, js:, k), yfx_ad(isd:, js:, k&
&                     ), area, ra_x, ra_x_ad, ra_y, ra_y_ad, mfx(is:, js&
&                     :, k), mfx_ad(is:, js:, k), mfy(is:, js:, k), &
&                     mfy_ad(is:, js:, k))
        END DO
        DO j=je,js,-1
          DO i=ie,is,-1
            CALL POPREAL8(dp2(i, j))
            temp_ad = rarea(i, j)*dp2_ad(i, j)
            dp1_ad(i, j, k) = dp1_ad(i, j, k) + dp2_ad(i, j)
            mfx_ad(i, j, k) = mfx_ad(i, j, k) + temp_ad
            mfx_ad(i+1, j, k) = mfx_ad(i+1, j, k) - temp_ad
            mfy_ad(i, j, k) = mfy_ad(i, j, k) + temp_ad
            mfy_ad(i, j+1, k) = mfy_ad(i, j+1, k) - temp_ad
            dp2_ad(i, j) = 0.0_8
          END DO
        END DO
        DO j=je,js,-1
          DO i=ied,isd,-1
            CALL POPREAL8(ra_y(i, j))
            yfx_ad(i, j, k) = yfx_ad(i, j, k) + ra_y_ad(i, j)
            yfx_ad(i, j+1, k) = yfx_ad(i, j+1, k) - ra_y_ad(i, j)
            ra_y_ad(i, j) = 0.0_8
          END DO
        END DO
        DO j=jed,jsd,-1
          DO i=ie,is,-1
            CALL POPREAL8(ra_x(i, j))
            xfx_ad(i, j, k) = xfx_ad(i, j, k) + ra_x_ad(i, j)
            xfx_ad(i+1, j, k) = xfx_ad(i+1, j, k) - ra_x_ad(i, j)
            ra_x_ad(i, j) = 0.0_8
          END DO
        END DO
      END DO
      call mpp_update_domains( q_ad, domain, complete= .true. )
    END DO
    CALL POPCONTROL1B(branch)
    IF (branch .NE. 0) THEN
      DO k=npz,1,-1
        DO j=je+1,js,-1
          DO i=ie,is,-1
            mfy_ad(i, j, k) = frac*mfy_ad(i, j, k)
          END DO
        END DO
        DO j=je+1,js,-1
          DO i=ied,isd,-1
            yfx_ad(i, j, k) = frac*yfx_ad(i, j, k)
            cy_ad(i, j, k) = frac*cy_ad(i, j, k)
          END DO
        END DO
        DO j=je,js,-1
          DO i=ie+1,is,-1
            mfx_ad(i, j, k) = frac*mfx_ad(i, j, k)
          END DO
        END DO
        DO j=jed,jsd,-1
          DO i=ie+1,is,-1
            xfx_ad(i, j, k) = frac*xfx_ad(i, j, k)
            cx_ad(i, j, k) = frac*cx_ad(i, j, k)
          END DO
        END DO
      END DO
    END IF
    CALL POPCONTROL2B(branch)
    IF (branch .EQ. 0) THEN
      DO k=npz,2,-1
        CALL POPCONTROL1B(branch)
      END DO
    ELSE IF (branch .NE. 1) THEN
      GOTO 100
    END IF
    DO k=npz,1,-1
      DO j=je,js,-1
        DO i=ie,is,-1
          CALL POPCONTROL2B(branch)
        END DO
      END DO
    END DO
 100 DO k=npz,1,-1
      DO j=je+1,js,-1
        DO i=ied,isd,-1
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            cy_ad(i, j, k) = cy_ad(i, j, k) + dya(i, j)*dx(i, j)*sina_v(&
&             i, j)*yfx_ad(i, j, k)
            yfx_ad(i, j, k) = 0.0_8
          ELSE
            cy_ad(i, j, k) = cy_ad(i, j, k) + dx(i, j)*sina_v(i, j)*dya(&
&             i, j-1)*yfx_ad(i, j, k)
            yfx_ad(i, j, k) = 0.0_8
          END IF
        END DO
      END DO
      DO j=jed,jsd,-1
        DO i=ie+1,is,-1
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            cx_ad(i, j, k) = cx_ad(i, j, k) + dxa(i, j)*dy(i, j)*sina_u(&
&             i, j)*xfx_ad(i, j, k)
            xfx_ad(i, j, k) = 0.0_8
          ELSE
            cx_ad(i, j, k) = cx_ad(i, j, k) + dy(i, j)*sina_u(i, j)*dxa(&
&             i-1, j)*xfx_ad(i, j, k)
            xfx_ad(i, j, k) = 0.0_8
          END IF
        END DO
      END DO
    END DO
  END SUBROUTINE TRACER_2D_ADM
  SUBROUTINE TRACER_2D(q, dp1, mfx, mfy, cx, cy, npx, npy, npz, nq, hord&
&   , q_split, dt, id_divg)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx
    INTEGER, INTENT(IN) :: npy
    INTEGER, INTENT(IN) :: npz
! number of tracers to be advected
    INTEGER, INTENT(IN) :: nq
    INTEGER, INTENT(IN) :: hord
    INTEGER, INTENT(IN) :: q_split
    INTEGER, INTENT(IN) :: id_divg
    REAL, INTENT(IN) :: dt
! Tracers
    REAL, INTENT(INOUT) :: q(isd:ied, jsd:jed, npz, nq)
! DELP before dyn_core
    REAL, INTENT(INOUT) :: dp1(is:ie, js:je, npz)
! Mass Flux X-Dir
    REAL, INTENT(INOUT) :: mfx(is:ie+1, js:je, npz)
! Mass Flux Y-Dir
    REAL, INTENT(INOUT) :: mfy(is:ie, js:je+1, npz)
! Courant Number X-Dir
    REAL, INTENT(INOUT) :: cx(is:ie+1, jsd:jed, npz)
! Courant Number Y-Dir
    REAL, INTENT(INOUT) :: cy(isd:ied, js:je+1, npz)
! Local Arrays
    REAL :: dp2(is:ie, js:je)
    REAL :: fx(is:ie+1, js:je)
    REAL :: fy(is:ie, js:je+1)
    REAL :: ra_x(is:ie, jsd:jed)
    REAL :: ra_y(isd:ied, js:je)
    REAL :: xfx(is:ie+1, jsd:jed, npz)
    REAL :: yfx(isd:ied, js:je+1, npz)
    REAL :: cmax(npz)
    REAL :: c_global
    REAL :: frac, rdt
    INTEGER :: nsplt
    INTEGER :: i, j, k, it, iq
    INTRINSIC ABS
    INTRINSIC MAX
    INTRINSIC INT
    INTRINSIC REAL
    REAL :: x1
    REAL :: abs1
    REAL :: abs0
    REAL :: y1
    DO k=1,npz
      DO j=jsd,jed
        DO i=is,ie+1
          IF (cx(i, j, k) .GT. 0.) THEN
            xfx(i, j, k) = cx(i, j, k)*dxa(i-1, j)*dy(i, j)*sina_u(i, j)
          ELSE
            xfx(i, j, k) = cx(i, j, k)*dxa(i, j)*dy(i, j)*sina_u(i, j)
          END IF
        END DO
      END DO
      DO j=js,je+1
        DO i=isd,ied
          IF (cy(i, j, k) .GT. 0.) THEN
            yfx(i, j, k) = cy(i, j, k)*dya(i, j-1)*dx(i, j)*sina_v(i, j)
          ELSE
            yfx(i, j, k) = cy(i, j, k)*dya(i, j)*dx(i, j)*sina_v(i, j)
          END IF
        END DO
      END DO
    END DO
!--------------------------------------------------------------------------------
    IF (q_split .EQ. 0) THEN
! Determine nsplt
      DO k=1,npz
        cmax(k) = 0.
        DO j=js,je
          DO i=is,ie
            IF (cx(i, j, k) .GE. 0.) THEN
              abs0 = cx(i, j, k)
            ELSE
              abs0 = -cx(i, j, k)
            END IF
            x1 = abs0 + 1. - sina_u(i, j)
            IF (cy(i, j, k) .GE. 0.) THEN
              abs1 = cy(i, j, k)
            ELSE
              abs1 = -cy(i, j, k)
            END IF
            y1 = abs1 + 1. - sina_v(i, j)
            IF (x1 .LT. y1) THEN
              IF (y1 .LT. cmax(k)) THEN
                cmax(k) = cmax(k)
              ELSE
                cmax(k) = y1
              END IF
            ELSE IF (x1 .LT. cmax(k)) THEN
              cmax(k) = cmax(k)
            ELSE
              cmax(k) = x1
            END IF
          END DO
        END DO
      END DO
!oncall mp_reduce_max(cmax,npz)
      call mp_reduce_max(cmax,npz)
! find global max courant number and define nsplt to scale cx,cy,mfx,mfy
      c_global = cmax(1)
      IF (npz .NE. 1) THEN
! if NOT shallow water test case
        DO k=2,npz
          IF (cmax(k) .LT. c_global) THEN
            c_global = c_global
          ELSE
            c_global = cmax(k)
          END IF
        END DO
      END IF
      nsplt = INT(1. + c_global)
!     if ( gid == 0 .and. nsplt > 5 )  write(6,*) 'Tracer_2d_split=', nsplt, c_global
    ELSE
      nsplt = q_split
    END IF
!--------------------------------------------------------------------------------
    frac = 1./REAL(nsplt)
    IF (nsplt .NE. 1) THEN
      DO k=1,npz
        DO j=jsd,jed
          DO i=is,ie+1
            cx(i, j, k) = cx(i, j, k)*frac
            xfx(i, j, k) = xfx(i, j, k)*frac
          END DO
        END DO
        DO j=js,je
          DO i=is,ie+1
            mfx(i, j, k) = mfx(i, j, k)*frac
          END DO
        END DO
        DO j=js,je+1
          DO i=isd,ied
            cy(i, j, k) = cy(i, j, k)*frac
            yfx(i, j, k) = yfx(i, j, k)*frac
          END DO
        END DO
        DO j=js,je+1
          DO i=is,ie
            mfy(i, j, k) = mfy(i, j, k)*frac
          END DO
        END DO
      END DO
    END IF
    DO it=1,nsplt
!     call timing_on('COMM_TOTAL')
!       call timing_on('COMM_TRAC')
!call mpp_update_domains( q, domain, complete=.true. )
      call mpp_update_domains( q, domain, complete=.true. )
!              call timing_off('COMM_TRAC')
!     call timing_off('COMM_TOTAL')
      DO k=1,npz
        DO j=jsd,jed
          DO i=is,ie
            ra_x(i, j) = area(i, j) + xfx(i, j, k) - xfx(i+1, j, k)
          END DO
        END DO
        DO j=js,je
          DO i=isd,ied
            ra_y(i, j) = area(i, j) + yfx(i, j, k) - yfx(i, j+1, k)
          END DO
        END DO
        DO j=js,je
          DO i=is,ie
            dp2(i, j) = dp1(i, j, k) + (mfx(i, j, k)-mfx(i+1, j, k)+mfy(&
&             i, j, k)-mfy(i, j+1, k))*rarea(i, j)
          END DO
        END DO
        DO iq=1,nq
          CALL FV_TP_2D(q(isd:, jsd:, k, iq), cx(is:, jsd:, k), cy(isd:&
&                 , js:, k), npx, npy, hord, fx, fy, xfx(is:, jsd:, k), &
&                 yfx(isd:, js:, k), area, ra_x, ra_y, mfx=mfx(is:, js:&
&                 , k), mfy=mfy(is:, js:, k))
          DO j=js,je
            DO i=is,ie
              q(i, j, k, iq) = (q(i, j, k, iq)*dp1(i, j, k)+(fx(i, j)-fx&
&               (i+1, j)+fy(i, j)-fy(i, j+1))*rarea(i, j))/dp2(i, j)
            END DO
          END DO
        END DO
        DO j=js,je
          DO i=is,ie
            dp1(i, j, k) = dp2(i, j)
          END DO
        END DO
      END DO
    END DO
! npz
! nsplt
    IF (id_divg .GT. 0) THEN
      rdt = 1./(frac*dt)
      DO k=1,npz
        DO j=js,je
          DO i=is,ie
            dp1(i, j, k) = (xfx(i+1, j, k)-xfx(i, j, k)+yfx(i, j+1, k)-&
&             yfx(i, j, k))*rarea(i, j)*rdt
          END DO
        END DO
      END DO
    END IF
! Rescale Mass fluxes for exports to offline tracer transport
    IF (nsplt .NE. 1) THEN
      DO k=1,npz
        DO j=js,je
          DO i=is,ie+1
            mfx(i, j, k) = mfx(i, j, k)/frac
          END DO
        END DO
        DO j=js,je+1
          DO i=is,ie
            mfy(i, j, k) = mfy(i, j, k)/frac
          END DO
        END DO
      END DO
    END IF
  END SUBROUTINE TRACER_2D

end module fv_tracer2d_adm_mod
