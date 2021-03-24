module fv_tracer2d_tlm_mod

 use tp_core_tlm_mod,     only: fv_tp_2d_tlm
 use fv_grid_tools_mod,   only: area, rarea, dxa, dya, dx, dy
 use fv_grid_utils_mod,   only: sina_u, sina_v
 use fv_mp_mod,           only: domain,mp_reduce_max,isd,ied,jsd,jed,is,js,ie,je
 use mpp_domains_mod,     only: mpp_update_domains

 implicit none
 private

 public :: tracer_2d_tlm, tracer_2d_1L_tlm

 contains

  SUBROUTINE TRACER_2D_1L_TLM(q, q_tl, dp0, dp0_tl, mfx, mfx_tl, mfy, &
&   mfy_tl, cx, cx_tl, cy, cy_tl, npx, npy, npz, nq, hord, q_split, k, &
&   q3, q3_tl, dt, id_divg)
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
    REAL, INTENT(INOUT) :: q_tl(isd:ied, jsd:jed, nq)
! Tracers
    REAL, INTENT(INOUT) :: q3(isd:ied, jsd:jed, npz, nq)
    REAL, INTENT(INOUT) :: q3_tl(isd:ied, jsd:jed, npz, nq)
! DELP before dyn_core
    REAL, INTENT(INOUT) :: dp0(is:ie, js:je)
    REAL, INTENT(INOUT) :: dp0_tl(is:ie, js:je)
! Mass Flux X-Dir
    REAL, INTENT(IN) :: mfx(is:ie+1, js:je)
    REAL, INTENT(IN) :: mfx_tl(is:ie+1, js:je)
! Mass Flux Y-Dir
    REAL, INTENT(IN) :: mfy(is:ie, js:je+1)
    REAL, INTENT(IN) :: mfy_tl(is:ie, js:je+1)
! Courant Number X-Dir
    REAL, INTENT(IN) :: cx(is:ie+1, jsd:jed)
    REAL, INTENT(IN) :: cx_tl(is:ie+1, jsd:jed)
! Courant Number Y-Dir
    REAL, INTENT(IN) :: cy(isd:ied, js:je+1)
    REAL, INTENT(IN) :: cy_tl(isd:ied, js:je+1)
! Local Arrays
    REAL :: mfx2(is:ie+1, js:je)
    REAL :: mfx2_tl(is:ie+1, js:je)
    REAL :: mfy2(is:ie, js:je+1)
    REAL :: mfy2_tl(is:ie, js:je+1)
    REAL :: cx2(is:ie+1, jsd:jed)
    REAL :: cx2_tl(is:ie+1, jsd:jed)
    REAL :: cy2(isd:ied, js:je+1)
    REAL :: cy2_tl(isd:ied, js:je+1)
    REAL :: dp1(is:ie, js:je)
    REAL :: dp1_tl(is:ie, js:je)
    REAL :: dp2(is:ie, js:je)
    REAL :: dp2_tl(is:ie, js:je)
    REAL :: fx(is:ie+1, js:je)
    REAL :: fx_tl(is:ie+1, js:je)
    REAL :: fy(is:ie, js:je+1)
    REAL :: fy_tl(is:ie, js:je+1)
    REAL :: ra_x(is:ie, jsd:jed)
    REAL :: ra_x_tl(is:ie, jsd:jed)
    REAL :: ra_y(isd:ied, js:je)
    REAL :: ra_y_tl(isd:ied, js:je)
    REAL :: xfx(is:ie+1, jsd:jed)
    REAL :: xfx_tl(is:ie+1, jsd:jed)
    REAL :: yfx(isd:ied, js:je+1)
    REAL :: yfx_tl(isd:ied, js:je+1)
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
    xfx_tl = 0.0
    DO j=jsd,jed
      DO i=is,ie+1
        IF (cx(i, j) .GT. 0.) THEN
          xfx_tl(i, j) = dxa(i-1, j)*dy(i, j)*sina_u(i, j)*cx_tl(i, j)
          xfx(i, j) = cx(i, j)*dxa(i-1, j)*dy(i, j)*sina_u(i, j)
        ELSE
          xfx_tl(i, j) = dxa(i, j)*dy(i, j)*sina_u(i, j)*cx_tl(i, j)
          xfx(i, j) = cx(i, j)*dxa(i, j)*dy(i, j)*sina_u(i, j)
        END IF
      END DO
    END DO
    yfx_tl = 0.0
    DO j=js,je+1
      DO i=isd,ied
        IF (cy(i, j) .GT. 0.) THEN
          yfx_tl(i, j) = dya(i, j-1)*dx(i, j)*sina_v(i, j)*cy_tl(i, j)
          yfx(i, j) = cy(i, j)*dya(i, j-1)*dx(i, j)*sina_v(i, j)
        ELSE
          yfx_tl(i, j) = dya(i, j)*dx(i, j)*sina_v(i, j)*cy_tl(i, j)
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
      call mp_reduce_max(cmax)
      !CALL MP_REDUCE_MAX_DUM0(cmax)
      nsplt = INT(1.01 + cmax)
!        if ( gid == 0 .and. nsplt > 5 )  write(6,*) k, 'Tracer_2d_split=', nsplt, cmax
    ELSE
      nsplt = q_split
    END IF
    frac = 1./REAL(nsplt)
    cx2_tl = 0.0
    DO j=jsd,jed
      DO i=is,ie+1
        cx2_tl(i, j) = frac*cx_tl(i, j)
        cx2(i, j) = cx(i, j)*frac
        xfx_tl(i, j) = frac*xfx_tl(i, j)
        xfx(i, j) = xfx(i, j)*frac
      END DO
    END DO
    mfx2_tl = 0.0
    DO j=js,je
      DO i=is,ie+1
        mfx2_tl(i, j) = frac*mfx_tl(i, j)
        mfx2(i, j) = mfx(i, j)*frac
      END DO
    END DO
    cy2_tl = 0.0
    DO j=js,je+1
      DO i=isd,ied
        cy2_tl(i, j) = frac*cy_tl(i, j)
        cy2(i, j) = cy(i, j)*frac
        yfx_tl(i, j) = frac*yfx_tl(i, j)
        yfx(i, j) = yfx(i, j)*frac
      END DO
    END DO
    mfy2_tl = 0.0
    DO j=js,je+1
      DO i=is,ie
        mfy2_tl(i, j) = frac*mfy_tl(i, j)
        mfy2(i, j) = mfy(i, j)*frac
      END DO
    END DO
    ra_x_tl = 0.0
    DO j=jsd,jed
      DO i=is,ie
        ra_x_tl(i, j) = xfx_tl(i, j) - xfx_tl(i+1, j)
        ra_x(i, j) = area(i, j) + xfx(i, j) - xfx(i+1, j)
      END DO
    END DO
    ra_y_tl = 0.0
    DO j=js,je
      DO i=isd,ied
        ra_y_tl(i, j) = yfx_tl(i, j) - yfx_tl(i, j+1)
        ra_y(i, j) = area(i, j) + yfx(i, j) - yfx(i, j+1)
      END DO
    END DO
    dp1_tl = 0.0
    DO j=js,je
      DO i=is,ie
        dp1_tl(i, j) = dp0_tl(i, j)
        dp1(i, j) = dp0(i, j)
      END DO
    END DO
    dp2_tl = 0.0
    fx_tl = 0.0
    fy_tl = 0.0
    DO it=1,nsplt
      DO j=js,je
        DO i=is,ie
          dp2_tl(i, j) = dp1_tl(i, j) + rarea(i, j)*(mfx2_tl(i, j)-&
&           mfx2_tl(i+1, j)+mfy2_tl(i, j)-mfy2_tl(i, j+1))
          dp2(i, j) = dp1(i, j) + (mfx2(i, j)-mfx2(i+1, j)+mfy2(i, j)-&
&           mfy2(i, j+1))*rarea(i, j)
        END DO
      END DO
      call mpp_update_domains( q_tl, domain, complete= .true. )
      call mpp_update_domains( q   , domain, complete= .true. )
      !CALL MPP_UPDATE_DOMAINS_DUM3_TLM(q, q_tl, isd, ied, jsd, jed, nq)
      DO iq=1,nq
        CALL FV_TP_2D_TLM(q(isd:, jsd:, iq), q_tl(isd:, jsd:, iq), cx2, &
&                   cx2_tl, cy2, cy2_tl, npx, npy, hord, fx, fx_tl, fy, &
&                   fy_tl, xfx, xfx_tl, yfx, yfx_tl, area, ra_x, ra_x_tl&
&                   , ra_y, ra_y_tl, mfx=mfx2, mfx_tl=mfx2_tl, mfy=mfy2&
&                   , mfy_tl=mfy2_tl)
        IF (it .EQ. nsplt) THEN
          DO j=js,je
            DO i=is,ie
              q3_tl(i, j, k, iq) = ((q_tl(i, j, iq)*dp1(i, j)+q(i, j, iq&
&               )*dp1_tl(i, j)+rarea(i, j)*(fx_tl(i, j)-fx_tl(i+1, j)+&
&               fy_tl(i, j)-fy_tl(i, j+1)))*dp2(i, j)-(q(i, j, iq)*dp1(i&
&               , j)+(fx(i, j)-fx(i+1, j)+fy(i, j)-fy(i, j+1))*rarea(i, &
&               j))*dp2_tl(i, j))/dp2(i, j)**2
              q3(i, j, k, iq) = (q(i, j, iq)*dp1(i, j)+(fx(i, j)-fx(i+1&
&               , j)+fy(i, j)-fy(i, j+1))*rarea(i, j))/dp2(i, j)
            END DO
          END DO
        ELSE
          DO j=js,je
            DO i=is,ie
              q_tl(i, j, iq) = ((q_tl(i, j, iq)*dp1(i, j)+q(i, j, iq)*&
&               dp1_tl(i, j)+rarea(i, j)*(fx_tl(i, j)-fx_tl(i+1, j)+&
&               fy_tl(i, j)-fy_tl(i, j+1)))*dp2(i, j)-(q(i, j, iq)*dp1(i&
&               , j)+(fx(i, j)-fx(i+1, j)+fy(i, j)-fy(i, j+1))*rarea(i, &
&               j))*dp2_tl(i, j))/dp2(i, j)**2
              q(i, j, iq) = (q(i, j, iq)*dp1(i, j)+(fx(i, j)-fx(i+1, j)+&
&               fy(i, j)-fy(i, j+1))*rarea(i, j))/dp2(i, j)
            END DO
          END DO
        END IF
      END DO
      IF (it .NE. nsplt) THEN
        DO j=js,je
          DO i=is,ie
            dp1_tl(i, j) = dp2_tl(i, j)
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
          dp0_tl(i, j) = rarea(i, j)*rdt*(xfx_tl(i+1, j)-xfx_tl(i, j)+&
&           yfx_tl(i, j+1)-yfx_tl(i, j))
          dp0(i, j) = (xfx(i+1, j)-xfx(i, j)+yfx(i, j+1)-yfx(i, j))*&
&           rarea(i, j)*rdt
        END DO
      END DO
    END IF
  END SUBROUTINE TRACER_2D_1L_TLM

  SUBROUTINE TRACER_2D_TLM(q, q_tl, dp1, dp1_tl, mfx, mfx_tl, mfy, &
&   mfy_tl, cx, cx_tl, cy, cy_tl, npx, npy, npz, nq, hord, q_split, dt, &
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
    REAL, INTENT(INOUT) :: q_tl(isd:ied, jsd:jed, npz, nq)
! DELP before dyn_core
    REAL, INTENT(INOUT) :: dp1(is:ie, js:je, npz)
    REAL, INTENT(INOUT) :: dp1_tl(is:ie, js:je, npz)
! Mass Flux X-Dir
    REAL, INTENT(INOUT) :: mfx(is:ie+1, js:je, npz)
    REAL, INTENT(INOUT) :: mfx_tl(is:ie+1, js:je, npz)
! Mass Flux Y-Dir
    REAL, INTENT(INOUT) :: mfy(is:ie, js:je+1, npz)
    REAL, INTENT(INOUT) :: mfy_tl(is:ie, js:je+1, npz)
! Courant Number X-Dir
    REAL, INTENT(INOUT) :: cx(is:ie+1, jsd:jed, npz)
    REAL, INTENT(INOUT) :: cx_tl(is:ie+1, jsd:jed, npz)
! Courant Number Y-Dir
    REAL, INTENT(INOUT) :: cy(isd:ied, js:je+1, npz)
    REAL, INTENT(INOUT) :: cy_tl(isd:ied, js:je+1, npz)
! Local Arrays
    REAL :: dp2(is:ie, js:je)
    REAL :: dp2_tl(is:ie, js:je)
    REAL :: fx(is:ie+1, js:je)
    REAL :: fx_tl(is:ie+1, js:je)
    REAL :: fy(is:ie, js:je+1)
    REAL :: fy_tl(is:ie, js:je+1)
    REAL :: ra_x(is:ie, jsd:jed)
    REAL :: ra_x_tl(is:ie, jsd:jed)
    REAL :: ra_y(isd:ied, js:je)
    REAL :: ra_y_tl(isd:ied, js:je)
    REAL :: xfx(is:ie+1, jsd:jed, npz)
    REAL :: xfx_tl(is:ie+1, jsd:jed, npz)
    REAL :: yfx(isd:ied, js:je+1, npz)
    REAL :: yfx_tl(isd:ied, js:je+1, npz)
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
    xfx_tl = 0.0
    yfx_tl = 0.0
    DO k=1,npz
      DO j=jsd,jed
        DO i=is,ie+1
          IF (cx(i, j, k) .GT. 0.) THEN
            xfx_tl(i, j, k) = dxa(i-1, j)*dy(i, j)*sina_u(i, j)*cx_tl(i&
&             , j, k)
            xfx(i, j, k) = cx(i, j, k)*dxa(i-1, j)*dy(i, j)*sina_u(i, j)
          ELSE
            xfx_tl(i, j, k) = dxa(i, j)*dy(i, j)*sina_u(i, j)*cx_tl(i, j&
&             , k)
            xfx(i, j, k) = cx(i, j, k)*dxa(i, j)*dy(i, j)*sina_u(i, j)
          END IF
        END DO
      END DO
      DO j=js,je+1
        DO i=isd,ied
          IF (cy(i, j, k) .GT. 0.) THEN
            yfx_tl(i, j, k) = dya(i, j-1)*dx(i, j)*sina_v(i, j)*cy_tl(i&
&             , j, k)
            yfx(i, j, k) = cy(i, j, k)*dya(i, j-1)*dx(i, j)*sina_v(i, j)
          ELSE
            yfx_tl(i, j, k) = dya(i, j)*dx(i, j)*sina_v(i, j)*cy_tl(i, j&
&             , k)
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
      call mp_reduce_max(cmax,npz)
      !CALL MP_REDUCE_MAX_DUM1(cmax, npz)
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
            cx_tl(i, j, k) = frac*cx_tl(i, j, k)
            cx(i, j, k) = cx(i, j, k)*frac
            xfx_tl(i, j, k) = frac*xfx_tl(i, j, k)
            xfx(i, j, k) = xfx(i, j, k)*frac
          END DO
        END DO
        DO j=js,je
          DO i=is,ie+1
            mfx_tl(i, j, k) = frac*mfx_tl(i, j, k)
            mfx(i, j, k) = mfx(i, j, k)*frac
          END DO
        END DO
        DO j=js,je+1
          DO i=isd,ied
            cy_tl(i, j, k) = frac*cy_tl(i, j, k)
            cy(i, j, k) = cy(i, j, k)*frac
            yfx_tl(i, j, k) = frac*yfx_tl(i, j, k)
            yfx(i, j, k) = yfx(i, j, k)*frac
          END DO
        END DO
        DO j=js,je+1
          DO i=is,ie
            mfy_tl(i, j, k) = frac*mfy_tl(i, j, k)
            mfy(i, j, k) = mfy(i, j, k)*frac
          END DO
        END DO
      END DO
      dp2_tl = 0.0
      ra_x_tl = 0.0
      ra_y_tl = 0.0
      fx_tl = 0.0
      fy_tl = 0.0
    ELSE
      dp2_tl = 0.0
      ra_x_tl = 0.0
      ra_y_tl = 0.0
      fx_tl = 0.0
      fy_tl = 0.0
    END IF
    DO it=1,nsplt
      call mpp_update_domains( q_tl, domain, complete=.true. )
      call mpp_update_domains( q   , domain, complete=.true. )
      !CALL MPP_UPDATE_DOMAINS_DUM4_TLM(q, q_tl, isd, ied, jsd, jed, npz&
!&                                , nq)
      DO k=1,npz
        DO j=jsd,jed
          DO i=is,ie
            ra_x_tl(i, j) = xfx_tl(i, j, k) - xfx_tl(i+1, j, k)
            ra_x(i, j) = area(i, j) + xfx(i, j, k) - xfx(i+1, j, k)
          END DO
        END DO
        DO j=js,je
          DO i=isd,ied
            ra_y_tl(i, j) = yfx_tl(i, j, k) - yfx_tl(i, j+1, k)
            ra_y(i, j) = area(i, j) + yfx(i, j, k) - yfx(i, j+1, k)
          END DO
        END DO
        DO j=js,je
          DO i=is,ie
            dp2_tl(i, j) = dp1_tl(i, j, k) + rarea(i, j)*(mfx_tl(i, j, k&
&             )-mfx_tl(i+1, j, k)+mfy_tl(i, j, k)-mfy_tl(i, j+1, k))
            dp2(i, j) = dp1(i, j, k) + (mfx(i, j, k)-mfx(i+1, j, k)+mfy(&
&             i, j, k)-mfy(i, j+1, k))*rarea(i, j)
          END DO
        END DO
        DO iq=1,nq
          CALL FV_TP_2D_TLM(q(isd:, jsd:, k, iq), q_tl(isd:, jsd:, k, iq&
&                     ), cx(is:, jsd:, k), cx_tl(is:, jsd:, k), cy(isd:&
&                     , js:, k), cy_tl(isd:, js:, k), npx, npy, hord, fx&
&                     , fx_tl, fy, fy_tl, xfx(is:, jsd:, k), xfx_tl(is:&
&                     , jsd:, k), yfx(isd:, js:, k), yfx_tl(isd:, js:, k&
&                     ), area, ra_x, ra_x_tl, ra_y, ra_y_tl, mfx=mfx(is:&
&                     , js:, k), mfx_tl=mfx_tl(is:, js:, k), mfy=mfy(is:&
&                     , js:, k), mfy_tl=mfy_tl(is:, js:, k))
          DO j=js,je
            DO i=is,ie
              q_tl(i, j, k, iq) = ((q_tl(i, j, k, iq)*dp1(i, j, k)+q(i, &
&               j, k, iq)*dp1_tl(i, j, k)+rarea(i, j)*(fx_tl(i, j)-fx_tl&
&               (i+1, j)+fy_tl(i, j)-fy_tl(i, j+1)))*dp2(i, j)-(q(i, j, &
&               k, iq)*dp1(i, j, k)+(fx(i, j)-fx(i+1, j)+fy(i, j)-fy(i, &
&               j+1))*rarea(i, j))*dp2_tl(i, j))/dp2(i, j)**2
              q(i, j, k, iq) = (q(i, j, k, iq)*dp1(i, j, k)+(fx(i, j)-fx&
&               (i+1, j)+fy(i, j)-fy(i, j+1))*rarea(i, j))/dp2(i, j)
            END DO
          END DO
        END DO
        DO j=js,je
          DO i=is,ie
            dp1_tl(i, j, k) = dp2_tl(i, j)
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
            dp1_tl(i, j, k) = rarea(i, j)*rdt*(xfx_tl(i+1, j, k)-xfx_tl(&
&             i, j, k)+yfx_tl(i, j+1, k)-yfx_tl(i, j, k))
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
            mfx_tl(i, j, k) = mfx_tl(i, j, k)/frac
            mfx(i, j, k) = mfx(i, j, k)/frac
          END DO
        END DO
        DO j=js,je+1
          DO i=is,ie
            mfy_tl(i, j, k) = mfy_tl(i, j, k)/frac
            mfy(i, j, k) = mfy(i, j, k)/frac
          END DO
        END DO
      END DO
    END IF
  END SUBROUTINE TRACER_2D_TLM

end module fv_tracer2d_tlm_mod
