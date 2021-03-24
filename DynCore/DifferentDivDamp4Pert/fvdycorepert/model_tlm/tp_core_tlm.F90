module tp_core_tlm_mod

 use fv_arrays_mod,     only: g_precision
 use fv_mp_mod,         only: is,js,ie,je, ng, isd,jsd,ied,jed
 use fv_grid_utils_mod, only: sw_corner, se_corner, ne_corner, nw_corner
 use fv_grid_utils_mod, only: sina_u, sina_v, da_min
 use fv_grid_tools_mod, only: dx, dy, rdxc, rdyc, rarea
 use fv_grid_tools_mod, only: dxa, dya, grid_type, rdx, rdy

 implicit none

 private
 public fv_tp_2d_tlm, pert_ppm_tlm, copy_corners_tlm

 real, parameter:: r3 = 1./3.

! scheme 2.1: perturbation form
 real, parameter:: b1 =   1./30.
 real, parameter:: b2 = -13./60.
 real, parameter:: b3 = -13./60.
 real, parameter:: b4 =  0.45
 real, parameter:: b5 = -0.05
 real, parameter:: t11 = 27./28., t12 = -13./28., t13=3./7.
 real, parameter:: s11 = 11./14., s14 = 4./7.,    s15=3./14.
  real, parameter:: c1 = -2./14.
  real, parameter:: c2 = 11./14.
  real, parameter:: c3 =  5./14.
  real, parameter:: p1 =  7./12.     ! 0.58333333
  real, parameter:: p2 = -1./12.

CONTAINS

  SUBROUTINE FV_TP_2D_TLM(q, q_tl, crx, crx_tl, cry, cry_tl, npx, npy, &
&   hord, fx, fx_tl, fy, fy_tl, xfx, xfx_tl, yfx, yfx_tl, area, ra_x, &
&   ra_x_tl, ra_y, ra_y_tl, mfx, mfx_tl, mfy, mfy_tl, ppm_fac, nord, &
&   damp_c)
    IMPLICIT NONE
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
    CALL FV_TP_2D_COMPUTE_TLM(q, q_tl, crx, crx_tl, cry, cry_tl, npx, &
&                       npy, hord, fx, fx_tl, fy, fy_tl, xfx, xfx_tl, &
&                       yfx, yfx_tl, area, ra_x, ra_x_tl, ra_y, ra_y_tl&
&                       , mfx, mfx_tl, mfy, mfy_tl, ppm_fac, nord, &
&                       damp_c)
  END SUBROUTINE FV_TP_2D_TLM

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
    fy2_tl = 0.0
    CALL YTP_TLM(fy2, fy2_tl, q, q_tl, cry, cry_tl, ord_in, isd, ied, js&
&          , je, npx, npy, ppm_limiter)
    fyy_tl = 0.0
    DO j=js,je+1
      DO i=isd,ied
        fyy_tl(i, j) = yfx_tl(i, j)*fy2(i, j) + yfx(i, j)*fy2_tl(i, j)
        fyy(i, j) = yfx(i, j)*fy2(i, j)
      END DO
    END DO
    q_i_tl = 0.0
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
    fx2_tl = 0.0
    CALL XTP_TLM(fx2, fx2_tl, q, q_tl, crx, crx_tl, ord_in, is, ie, jsd&
&          , jed, npx, npy, ppm_limiter)
    q_j_tl = 0.0
    fx1_tl = 0.0
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
    d2_tl = 0.0
    DO j=jsd,jed
      DO i=isd,ied
        d2_tl(i, j) = damp*q_tl(i, j)
        d2(i, j) = damp*q(i, j)
      END DO
    END DO
    IF (nord .GT. 0) THEN
      CALL COPY_CORNERS_TLM(d2, d2_tl, npx, npy, 1)
      fx2_tl = 0.0
    ELSE
      fx2_tl = 0.0
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
      fy2_tl = 0.0
    ELSE
      fy2_tl = 0.0
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

  SUBROUTINE COPY_CORNERS_TLM(q, q_tl, npx, npy, dir)
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
  END SUBROUTINE COPY_CORNERS_TLM

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
    REAL :: dm(ifirst:ilast, jfirst-2:jlast+2)
    REAL :: dm_tl(ifirst:ilast, jfirst-2:jlast+2)
    REAL :: x0, x1
    REAL :: x0_tl, x1_tl
    INTEGER :: i, j
    INTRINSIC ABS
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC SIGN
    REAL :: y1_tl
    REAL :: min6_tl
    REAL :: min10_tl
    REAL :: z4_tl
    REAL :: min9
    REAL :: min8
    REAL :: y2_tl
    REAL :: min7_tl
    REAL :: min7
    REAL :: min6
    REAL :: min5
    REAL :: min4
    REAL :: min3
    REAL :: min2
    REAL :: z5_tl
    REAL :: min1
    REAL :: y3_tl
    REAL :: min8_tl
    REAL :: x6
    REAL :: x5
    REAL :: x4
    REAL :: y4_tl
    REAL :: min9_tl
    REAL :: x3
    REAL :: x2
    REAL :: x2_tl
    REAL :: y5_tl
    REAL :: x3_tl
    REAL :: x4_tl
    REAL :: max1_tl
    REAL :: min1_tl
    REAL :: x5_tl
    REAL :: max2_tl
    REAL :: z5
    REAL :: min2_tl
    REAL :: z4
    REAL :: z3
    REAL :: x6_tl
    REAL :: z2
    REAL :: z1
    REAL :: max3_tl
    REAL :: min3_tl
    REAL :: min10
    REAL :: z1_tl
    REAL :: max4_tl
    REAL :: min4_tl
    REAL :: z2_tl
    REAL :: max5_tl
    REAL :: max5
    REAL :: max4
    REAL :: min5_tl
    REAL :: max3
    REAL :: max2
    REAL :: max1
    REAL :: y5
    REAL :: y4
    REAL :: y3
    REAL :: z3_tl
    REAL :: y2
    REAL :: y1
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
          IF (c(i, j) .GT. 0.) THEN
            fy_tl(i, j) = q_tl(i, j-1)
            fy(i, j) = q(i, j-1)
          ELSE
            fy_tl(i, j) = q_tl(i, j)
            fy(i, j) = q(i, j)
          END IF
        END DO
      END DO

      !Update to third order
      DO j=jfirst+1,jlast
        DO i=ifirst,ilast
          IF (c(i, j) .GT. 0.) THEN
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

      !Edge change to first order
      IF (grid_type .LT. 3) THEN
        IF (js .EQ. 1) THEN
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
        END IF
        IF (je + 1 .EQ. npy) THEN
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
        END IF
      END IF

    ELSE IF (jord .EQ. 2) THEN
      dm_tl = 0.0
      DO j=jfirst-2,jlast+2
        DO i=ifirst,ilast
          dm_tl(i, j) = 0.25*(q_tl(i, j+1)-q_tl(i, j-1))
          dm(i, j) = 0.25*(q(i, j+1)-q(i, j-1))
          IF (dm(i, j) .GE. 0.) THEN
            x2_tl = dm_tl(i, j)
            x2 = dm(i, j)
          ELSE
            x2_tl = -dm_tl(i, j)
            x2 = -dm(i, j)
          END IF
          IF (q(i, j-1) .LT. q(i, j)) THEN
            IF (q(i, j) .LT. q(i, j+1)) THEN
              max1_tl = q_tl(i, j+1)
              max1 = q(i, j+1)
            ELSE
              max1_tl = q_tl(i, j)
              max1 = q(i, j)
            END IF
          ELSE IF (q(i, j-1) .LT. q(i, j+1)) THEN
            max1_tl = q_tl(i, j+1)
            max1 = q(i, j+1)
          ELSE
            max1_tl = q_tl(i, j-1)
            max1 = q(i, j-1)
          END IF
          y1_tl = max1_tl - q_tl(i, j)
          y1 = max1 - q(i, j)
          IF (q(i, j-1) .GT. q(i, j)) THEN
            IF (q(i, j) .GT. q(i, j+1)) THEN
              min6_tl = q_tl(i, j+1)
              min6 = q(i, j+1)
            ELSE
              min6_tl = q_tl(i, j)
              min6 = q(i, j)
            END IF
          ELSE IF (q(i, j-1) .GT. q(i, j+1)) THEN
            min6_tl = q_tl(i, j+1)
            min6 = q(i, j+1)
          ELSE
            min6_tl = q_tl(i, j-1)
            min6 = q(i, j-1)
          END IF
          z1_tl = q_tl(i, j) - min6_tl
          z1 = q(i, j) - min6
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
          dm_tl(i, j) = min1_tl*SIGN(1.d0, min1*dm(i, j))
          dm(i, j) = SIGN(min1, dm(i, j))
        END DO
      END DO
!--------------
! Fix the edges:
!--------------
      IF (grid_type .LT. 3) THEN
        IF (js .EQ. 1) THEN
          DO i=ifirst,ilast
            x0_tl = 0.5*((2.*dya(i, 1)+dya(i, 2))*(q_tl(i, 0)+q_tl(i, 1)&
&             )-dya(i, 1)*(q_tl(i, -1)+q_tl(i, 2)))/(dya(i, 1)+dya(i, 2)&
&             )
            x0 = 0.5*((2.*dya(i, 1)+dya(i, 2))*(q(i, 0)+q(i, 1))-dya(i, &
&             1)*(q(i, -1)+q(i, 2)))/(dya(i, 1)+dya(i, 2))
            x1_tl = s15*q_tl(i, 1) + s11*q_tl(i, 2) - s14*dm_tl(i, 2)
            x1 = s15*q(i, 1) + s11*q(i, 2) - s14*dm(i, 2)
            dm_tl(i, 1) = 0.5*(x1_tl-x0_tl)
            dm(i, 1) = 0.5*(x1-x0)
            IF (dm(i, 1) .GE. 0.) THEN
              x3_tl = dm_tl(i, 1)
              x3 = dm(i, 1)
            ELSE
              x3_tl = -dm_tl(i, 1)
              x3 = -dm(i, 1)
            END IF
            IF (q(i, 1) .LT. x0) THEN
              IF (x0 .LT. x1) THEN
                max2_tl = x1_tl
                max2 = x1
              ELSE
                max2_tl = x0_tl
                max2 = x0
              END IF
            ELSE IF (q(i, 1) .LT. x1) THEN
              max2_tl = x1_tl
              max2 = x1
            ELSE
              max2_tl = q_tl(i, 1)
              max2 = q(i, 1)
            END IF
            y2_tl = max2_tl - q_tl(i, 1)
            y2 = max2 - q(i, 1)
            IF (q(i, 1) .GT. x0) THEN
              IF (x0 .GT. x1) THEN
                min7_tl = x1_tl
                min7 = x1
              ELSE
                min7_tl = x0_tl
                min7 = x0
              END IF
            ELSE IF (q(i, 1) .GT. x1) THEN
              min7_tl = x1_tl
              min7 = x1
            ELSE
              min7_tl = q_tl(i, 1)
              min7 = q(i, 1)
            END IF
            z2_tl = q_tl(i, 1) - min7_tl
            z2 = q(i, 1) - min7
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
!
            x1_tl = s15*q_tl(i, 0) + s11*q_tl(i, -1) + s14*dm_tl(i, -1)
            x1 = s15*q(i, 0) + s11*q(i, -1) + s14*dm(i, -1)
            dm_tl(i, 0) = 0.5*(x0_tl-x1_tl)
            dm(i, 0) = 0.5*(x0-x1)
            IF (dm(i, 0) .GE. 0.) THEN
              x4_tl = dm_tl(i, 0)
              x4 = dm(i, 0)
            ELSE
              x4_tl = -dm_tl(i, 0)
              x4 = -dm(i, 0)
            END IF
            IF (q(i, 0) .LT. x0) THEN
              IF (x0 .LT. x1) THEN
                max3_tl = x1_tl
                max3 = x1
              ELSE
                max3_tl = x0_tl
                max3 = x0
              END IF
            ELSE IF (q(i, 0) .LT. x1) THEN
              max3_tl = x1_tl
              max3 = x1
            ELSE
              max3_tl = q_tl(i, 0)
              max3 = q(i, 0)
            END IF
            y3_tl = max3_tl - q_tl(i, 0)
            y3 = max3 - q(i, 0)
            IF (q(i, 0) .GT. x0) THEN
              IF (x0 .GT. x1) THEN
                min8_tl = x1_tl
                min8 = x1
              ELSE
                min8_tl = x0_tl
                min8 = x0
              END IF
            ELSE IF (q(i, 0) .GT. x1) THEN
              min8_tl = x1_tl
              min8 = x1
            ELSE
              min8_tl = q_tl(i, 0)
              min8 = q(i, 0)
            END IF
            z3_tl = q_tl(i, 0) - min8_tl
            z3 = q(i, 0) - min8
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
        END IF
        IF (je + 1 .EQ. npy) THEN
          DO i=ifirst,ilast
            x0_tl = 0.5*((2.*dya(i, npy-1)+dya(i, npy-2))*(q_tl(i, npy-1&
&             )+q_tl(i, npy))-dya(i, npy-1)*(q_tl(i, npy-2)+q_tl(i, npy+&
&             1)))/(dya(i, npy-1)+dya(i, npy-2))
            x0 = 0.5*((2.*dya(i, npy-1)+dya(i, npy-2))*(q(i, npy-1)+q(i&
&             , npy))-dya(i, npy-1)*(q(i, npy-2)+q(i, npy+1)))/(dya(i, &
&             npy-1)+dya(i, npy-2))
            x1_tl = s15*q_tl(i, npy-1) + s11*q_tl(i, npy-2) + s14*dm_tl(&
&             i, npy-2)
            x1 = s15*q(i, npy-1) + s11*q(i, npy-2) + s14*dm(i, npy-2)
            dm_tl(i, npy-1) = 0.5*(x0_tl-x1_tl)
            dm(i, npy-1) = 0.5*(x0-x1)
            IF (dm(i, npy-1) .GE. 0.) THEN
              x5_tl = dm_tl(i, npy-1)
              x5 = dm(i, npy-1)
            ELSE
              x5_tl = -dm_tl(i, npy-1)
              x5 = -dm(i, npy-1)
            END IF
            IF (q(i, npy-1) .LT. x0) THEN
              IF (x0 .LT. x1) THEN
                max4_tl = x1_tl
                max4 = x1
              ELSE
                max4_tl = x0_tl
                max4 = x0
              END IF
            ELSE IF (q(i, npy-1) .LT. x1) THEN
              max4_tl = x1_tl
              max4 = x1
            ELSE
              max4_tl = q_tl(i, npy-1)
              max4 = q(i, npy-1)
            END IF
            y4_tl = max4_tl - q_tl(i, npy-1)
            y4 = max4 - q(i, npy-1)
            IF (q(i, npy-1) .GT. x0) THEN
              IF (x0 .GT. x1) THEN
                min9_tl = x1_tl
                min9 = x1
              ELSE
                min9_tl = x0_tl
                min9 = x0
              END IF
            ELSE IF (q(i, npy-1) .GT. x1) THEN
              min9_tl = x1_tl
              min9 = x1
            ELSE
              min9_tl = q_tl(i, npy-1)
              min9 = q(i, npy-1)
            END IF
            z4_tl = q_tl(i, npy-1) - min9_tl
            z4 = q(i, npy-1) - min9
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
!
            x1_tl = s15*q_tl(i, npy) + s11*q_tl(i, npy+1) - s14*dm_tl(i&
&             , npy+1)
            x1 = s15*q(i, npy) + s11*q(i, npy+1) - s14*dm(i, npy+1)
            dm_tl(i, npy) = 0.5*(x1_tl-x0_tl)
            dm(i, npy) = 0.5*(x1-x0)
            IF (dm(i, npy) .GE. 0.) THEN
              x6_tl = dm_tl(i, npy)
              x6 = dm(i, npy)
            ELSE
              x6_tl = -dm_tl(i, npy)
              x6 = -dm(i, npy)
            END IF
            IF (q(i, npy) .LT. x0) THEN
              IF (x0 .LT. x1) THEN
                max5_tl = x1_tl
                max5 = x1
              ELSE
                max5_tl = x0_tl
                max5 = x0
              END IF
            ELSE IF (q(i, npy) .LT. x1) THEN
              max5_tl = x1_tl
              max5 = x1
            ELSE
              max5_tl = q_tl(i, npy)
              max5 = q(i, npy)
            END IF
            y5_tl = max5_tl - q_tl(i, npy)
            y5 = max5 - q(i, npy)
            IF (q(i, npy) .GT. x0) THEN
              IF (x0 .GT. x1) THEN
                min10_tl = x1_tl
                min10 = x1
              ELSE
                min10_tl = x0_tl
                min10 = x0
              END IF
            ELSE IF (q(i, npy) .GT. x1) THEN
              min10_tl = x1_tl
              min10 = x1
            ELSE
              min10_tl = q_tl(i, npy)
              min10 = q(i, npy)
            END IF
            z5_tl = q_tl(i, npy) - min10_tl
            z5 = q(i, npy) - min10
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
        END IF
      END IF
      DO j=jfirst,jlast+1
        DO i=ifirst,ilast
          IF (c(i, j) .GT. 0.) THEN
            fy_tl(i, j) = q_tl(i, j-1) + (1.-c(i, j))*dm_tl(i, j-1) - &
&             c_tl(i, j)*dm(i, j-1)
            fy(i, j) = q(i, j-1) + (1.-c(i, j))*dm(i, j-1)
          ELSE
            fy_tl(i, j) = q_tl(i, j) - c_tl(i, j)*dm(i, j) - (1.+c(i, j)&
&             )*dm_tl(i, j)
            fy(i, j) = q(i, j) - (1.+c(i, j))*dm(i, j)
          END IF
        END DO
      END DO
    ELSE
      CALL FYPPM_TLM(c, c_tl, q, q_tl, fy, fy_tl, jord, ifirst, ilast, &
&              jfirst, jlast, npx, npy, dm, dm_tl, ppm_limiter)
    END IF
  END SUBROUTINE YTP_TLM

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
    REAL :: dm(is-2:ie+2)
    REAL :: dm_tl(is-2:ie+2)
    REAL :: x0, x1
    REAL :: x0_tl, x1_tl
    INTEGER :: i, j
    INTRINSIC ABS
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC SIGN
    REAL :: y1_tl
    REAL :: min6_tl
    REAL :: min10_tl
    REAL :: z4_tl
    REAL :: min9
    REAL :: min8
    REAL :: y2_tl
    REAL :: min7_tl
    REAL :: min7
    REAL :: min6
    REAL :: min5
    REAL :: min4
    REAL :: min3
    REAL :: min2
    REAL :: z5_tl
    REAL :: min1
    REAL :: y3_tl
    REAL :: min8_tl
    REAL :: x6
    REAL :: x5
    REAL :: x4
    REAL :: y4_tl
    REAL :: min9_tl
    REAL :: x3
    REAL :: x2
    REAL :: x2_tl
    REAL :: y5_tl
    REAL :: x3_tl
    REAL :: x4_tl
    REAL :: max1_tl
    REAL :: min1_tl
    REAL :: x5_tl
    REAL :: max2_tl
    REAL :: z5
    REAL :: min2_tl
    REAL :: z4
    REAL :: z3
    REAL :: x6_tl
    REAL :: z2
    REAL :: z1
    REAL :: max3_tl
    REAL :: min3_tl
    REAL :: min10
    REAL :: z1_tl
    REAL :: max4_tl
    REAL :: min4_tl
    REAL :: z2_tl
    REAL :: max5_tl
    REAL :: max5
    REAL :: max4
    REAL :: min5_tl
    REAL :: max3
    REAL :: max2
    REAL :: max1
    REAL :: y5
    REAL :: y4
    REAL :: y3
    REAL :: z3_tl
    REAL :: y2
    REAL :: y1
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
          IF (c(i, j) .GT. 0.) THEN
            fx_tl(i, j) = q_tl(i-1, j)
            fx(i, j) = q(i-1, j)
          ELSE
            fx_tl(i, j) = q_tl(i, j)
            fx(i, j) = q(i, j)
          END IF
        END DO
      END DO

      !Upgrade to third order scheme
      DO j=jfirst,jlast
        DO i=ifirst+1,ilast
          IF (c(i, j) .GT. 0.) THEN
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

      !Downgrade to first order for the edges
      IF (grid_type .LT. 3) THEN
        IF (is .EQ. 1) THEN
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
        ENDIF
        IF (ie + 1 .EQ. npx) THEN
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
        ENDIF
      ENDIF

    ELSE IF (iord .EQ. 2) THEN
      dm_tl = 0.0
      DO j=jfirst,jlast
        DO i=is-2,ie+2
          dm_tl(i) = 0.25*(q_tl(i+1, j)-q_tl(i-1, j))
          dm(i) = 0.25*(q(i+1, j)-q(i-1, j))
          IF (dm(i) .GE. 0.) THEN
            x2_tl = dm_tl(i)
            x2 = dm(i)
          ELSE
            x2_tl = -dm_tl(i)
            x2 = -dm(i)
          END IF
          IF (q(i-1, j) .LT. q(i, j)) THEN
            IF (q(i, j) .LT. q(i+1, j)) THEN
              max1_tl = q_tl(i+1, j)
              max1 = q(i+1, j)
            ELSE
              max1_tl = q_tl(i, j)
              max1 = q(i, j)
            END IF
          ELSE IF (q(i-1, j) .LT. q(i+1, j)) THEN
            max1_tl = q_tl(i+1, j)
            max1 = q(i+1, j)
          ELSE
            max1_tl = q_tl(i-1, j)
            max1 = q(i-1, j)
          END IF
          y1_tl = max1_tl - q_tl(i, j)
          y1 = max1 - q(i, j)
          IF (q(i-1, j) .GT. q(i, j)) THEN
            IF (q(i, j) .GT. q(i+1, j)) THEN
              min6_tl = q_tl(i+1, j)
              min6 = q(i+1, j)
            ELSE
              min6_tl = q_tl(i, j)
              min6 = q(i, j)
            END IF
          ELSE IF (q(i-1, j) .GT. q(i+1, j)) THEN
            min6_tl = q_tl(i+1, j)
            min6 = q(i+1, j)
          ELSE
            min6_tl = q_tl(i-1, j)
            min6 = q(i-1, j)
          END IF
          z1_tl = q_tl(i, j) - min6_tl
          z1 = q(i, j) - min6
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
          dm_tl(i) = min1_tl*SIGN(1.d0, min1*dm(i))
          dm(i) = SIGN(min1, dm(i))
        END DO
        IF (grid_type .LT. 3) THEN
!--------------
! fix the edges
!--------------
          IF (is .EQ. 1) THEN
            x0_tl = 0.5*((2.*dxa(1, j)+dxa(2, j))*(q_tl(0, j)+q_tl(1, j)&
&             )-dxa(1, j)*(q_tl(-1, j)+q_tl(2, j)))/(dxa(1, j)+dxa(2, j)&
&             )
            x0 = 0.5*((2.*dxa(1, j)+dxa(2, j))*(q(0, j)+q(1, j))-dxa(1, &
&             j)*(q(-1, j)+q(2, j)))/(dxa(1, j)+dxa(2, j))
            x1_tl = s15*q_tl(1, j) + s11*q_tl(2, j) - s14*dm_tl(2)
            x1 = s15*q(1, j) + s11*q(2, j) - s14*dm(2)
            dm_tl(1) = 0.5*(x1_tl-x0_tl)
            dm(1) = 0.5*(x1-x0)
            IF (dm(1) .GE. 0.) THEN
              x3_tl = dm_tl(1)
              x3 = dm(1)
            ELSE
              x3_tl = -dm_tl(1)
              x3 = -dm(1)
            END IF
            IF (q(1, j) .LT. x0) THEN
              IF (x0 .LT. x1) THEN
                max2_tl = x1_tl
                max2 = x1
              ELSE
                max2_tl = x0_tl
                max2 = x0
              END IF
            ELSE IF (q(1, j) .LT. x1) THEN
              max2_tl = x1_tl
              max2 = x1
            ELSE
              max2_tl = q_tl(1, j)
              max2 = q(1, j)
            END IF
            y2_tl = max2_tl - q_tl(1, j)
            y2 = max2 - q(1, j)
            IF (q(1, j) .GT. x0) THEN
              IF (x0 .GT. x1) THEN
                min7_tl = x1_tl
                min7 = x1
              ELSE
                min7_tl = x0_tl
                min7 = x0
              END IF
            ELSE IF (q(1, j) .GT. x1) THEN
              min7_tl = x1_tl
              min7 = x1
            ELSE
              min7_tl = q_tl(1, j)
              min7 = q(1, j)
            END IF
            z2_tl = q_tl(1, j) - min7_tl
            z2 = q(1, j) - min7
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
            dm_tl(1) = min2_tl*SIGN(1.d0, min2*dm(1))
            dm(1) = SIGN(min2, dm(1))
!
            x1_tl = s15*q_tl(0, j) + s11*q_tl(-1, j) + s14*dm_tl(-1)
            x1 = s15*q(0, j) + s11*q(-1, j) + s14*dm(-1)
            dm_tl(0) = 0.5*(x0_tl-x1_tl)
            dm(0) = 0.5*(x0-x1)
            IF (dm(0) .GE. 0.) THEN
              x4_tl = dm_tl(0)
              x4 = dm(0)
            ELSE
              x4_tl = -dm_tl(0)
              x4 = -dm(0)
            END IF
            IF (q(0, j) .LT. x0) THEN
              IF (x0 .LT. x1) THEN
                max3_tl = x1_tl
                max3 = x1
              ELSE
                max3_tl = x0_tl
                max3 = x0
              END IF
            ELSE IF (q(0, j) .LT. x1) THEN
              max3_tl = x1_tl
              max3 = x1
            ELSE
              max3_tl = q_tl(0, j)
              max3 = q(0, j)
            END IF
            y3_tl = max3_tl - q_tl(0, j)
            y3 = max3 - q(0, j)
            IF (q(0, j) .GT. x0) THEN
              IF (x0 .GT. x1) THEN
                min8_tl = x1_tl
                min8 = x1
              ELSE
                min8_tl = x0_tl
                min8 = x0
              END IF
            ELSE IF (q(0, j) .GT. x1) THEN
              min8_tl = x1_tl
              min8 = x1
            ELSE
              min8_tl = q_tl(0, j)
              min8 = q(0, j)
            END IF
            z3_tl = q_tl(0, j) - min8_tl
            z3 = q(0, j) - min8
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
            dm_tl(0) = min3_tl*SIGN(1.d0, min3*dm(0))
            dm(0) = SIGN(min3, dm(0))
          END IF
          IF (ie + 1 .EQ. npx) THEN
            x0_tl = 0.5*((2.*dxa(npx-1, j)+dxa(npx-2, j))*(q_tl(npx-1, j&
&             )+q_tl(npx, j))-dxa(npx-1, j)*(q_tl(npx-2, j)+q_tl(npx+1, &
&             j)))/(dxa(npx-1, j)+dxa(npx-2, j))
            x0 = 0.5*((2.*dxa(npx-1, j)+dxa(npx-2, j))*(q(npx-1, j)+q(&
&             npx, j))-dxa(npx-1, j)*(q(npx-2, j)+q(npx+1, j)))/(dxa(npx&
&             -1, j)+dxa(npx-2, j))
            x1_tl = s15*q_tl(npx-1, j) + s11*q_tl(npx-2, j) + s14*dm_tl(&
&             npx-2)
            x1 = s15*q(npx-1, j) + s11*q(npx-2, j) + s14*dm(npx-2)
            dm_tl(npx-1) = 0.5*(x0_tl-x1_tl)
            dm(npx-1) = 0.5*(x0-x1)
            IF (dm(npx-1) .GE. 0.) THEN
              x5_tl = dm_tl(npx-1)
              x5 = dm(npx-1)
            ELSE
              x5_tl = -dm_tl(npx-1)
              x5 = -dm(npx-1)
            END IF
            IF (q(npx-1, j) .LT. x0) THEN
              IF (x0 .LT. x1) THEN
                max4_tl = x1_tl
                max4 = x1
              ELSE
                max4_tl = x0_tl
                max4 = x0
              END IF
            ELSE IF (q(npx-1, j) .LT. x1) THEN
              max4_tl = x1_tl
              max4 = x1
            ELSE
              max4_tl = q_tl(npx-1, j)
              max4 = q(npx-1, j)
            END IF
            y4_tl = max4_tl - q_tl(npx-1, j)
            y4 = max4 - q(npx-1, j)
            IF (q(npx-1, j) .GT. x0) THEN
              IF (x0 .GT. x1) THEN
                min9_tl = x1_tl
                min9 = x1
              ELSE
                min9_tl = x0_tl
                min9 = x0
              END IF
            ELSE IF (q(npx-1, j) .GT. x1) THEN
              min9_tl = x1_tl
              min9 = x1
            ELSE
              min9_tl = q_tl(npx-1, j)
              min9 = q(npx-1, j)
            END IF
            z4_tl = q_tl(npx-1, j) - min9_tl
            z4 = q(npx-1, j) - min9
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
!
            x1_tl = s15*q_tl(npx, j) + s11*q_tl(npx+1, j) - s14*dm_tl(&
&             npx+1)
            x1 = s15*q(npx, j) + s11*q(npx+1, j) - s14*dm(npx+1)
            dm_tl(npx) = 0.5*(x1_tl-x0_tl)
            dm(npx) = 0.5*(x1-x0)
            IF (dm(npx) .GE. 0.) THEN
              x6_tl = dm_tl(npx)
              x6 = dm(npx)
            ELSE
              x6_tl = -dm_tl(npx)
              x6 = -dm(npx)
            END IF
            IF (q(npx, j) .LT. x0) THEN
              IF (x0 .LT. x1) THEN
                max5_tl = x1_tl
                max5 = x1
              ELSE
                max5_tl = x0_tl
                max5 = x0
              END IF
            ELSE IF (q(npx, j) .LT. x1) THEN
              max5_tl = x1_tl
              max5 = x1
            ELSE
              max5_tl = q_tl(npx, j)
              max5 = q(npx, j)
            END IF
            y5_tl = max5_tl - q_tl(npx, j)
            y5 = max5 - q(npx, j)
            IF (q(npx, j) .GT. x0) THEN
              IF (x0 .GT. x1) THEN
                min10_tl = x1_tl
                min10 = x1
              ELSE
                min10_tl = x0_tl
                min10 = x0
              END IF
            ELSE IF (q(npx, j) .GT. x1) THEN
              min10_tl = x1_tl
              min10 = x1
            ELSE
              min10_tl = q_tl(npx, j)
              min10 = q(npx, j)
            END IF
            z5_tl = q_tl(npx, j) - min10_tl
            z5 = q(npx, j) - min10
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
        DO i=is,ie+1
          IF (c(i, j) .GT. 0.) THEN
            fx_tl(i, j) = q_tl(i-1, j) + (1.-c(i, j))*dm_tl(i-1) - c_tl(&
&             i, j)*dm(i-1)
            fx(i, j) = q(i-1, j) + (1.-c(i, j))*dm(i-1)
          ELSE
            fx_tl(i, j) = q_tl(i, j) - c_tl(i, j)*dm(i) - (1.+c(i, j))*&
&             dm_tl(i)
            fx(i, j) = q(i, j) - (1.+c(i, j))*dm(i)
          END IF
        END DO
      END DO
    ELSE
      CALL FXPPM_TLM(c, c_tl, q, q_tl, fx, fx_tl, iord, ifirst, ilast, &
&              jfirst, jlast, npx, npy, ppm_limiter)
    END IF
  END SUBROUTINE XTP_TLM

  SUBROUTINE FYPPM_TLM(c, c_tl, q, q_tl, flux, flux_tl, jord, ifirst, &
&   ilast, jfirst, jlast, npx, npy, dm, dm_tl, ppm_limiter)
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
    REAL, INTENT(OUT) :: dm(ifirst:ilast, jfirst-2:jlast+2)
    REAL, INTENT(OUT) :: dm_tl(ifirst:ilast, jfirst-2:jlast+2)
! Local:
    REAL :: al(ifirst:ilast, jfirst-1:jlast+2)
    REAL :: al_tl(ifirst:ilast, jfirst-1:jlast+2)
    REAL :: bl(ifirst:ilast, jfirst-1:jlast+1)
    REAL :: bl_tl(ifirst:ilast, jfirst-1:jlast+1)
    REAL :: br(ifirst:ilast, jfirst-1:jlast+1)
    REAL :: br_tl(ifirst:ilast, jfirst-1:jlast+1)
    REAL :: dq(ifirst:ilast, jfirst-3:jlast+2)
    REAL :: dq_tl(ifirst:ilast, jfirst-3:jlast+2)
    REAL :: dl, dr, pmp, lac, ct, qe
    REAL :: dl_tl, dr_tl, pmp_tl, lac_tl
    REAL :: xt, x0, x1
    REAL :: xt_tl, x0_tl, x1_tl
    INTEGER :: i, j, js3, je3, je2, jt
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC ABS
    INTRINSIC SIGN
    REAL :: y1_tl
    REAL :: min6_tl
    REAL :: y18_tl
    REAL :: y44_tl
    REAL :: x19
    REAL :: x25_tl
    REAL :: x18
    REAL :: y32_tl
    REAL :: x17
    REAL :: min10_tl
    REAL :: x16
    REAL :: z4_tl
    REAL :: x13_tl
    REAL :: x15
    REAL :: y20_tl
    REAL :: x14
    REAL :: min9
    REAL :: min35_tl
    REAL :: x13
    REAL :: min8
    REAL :: y2_tl
    REAL :: min7_tl
    REAL :: y19_tl
    REAL :: x12
    REAL :: min7
    REAL :: y45_tl
    REAL :: x11
    REAL :: y29
    REAL :: min6
    REAL :: min35
    REAL :: min23_tl
    REAL :: x10
    REAL :: y28
    REAL :: min5
    INTEGER :: min34
    REAL :: x26_tl
    REAL :: y27
    REAL :: min4
    INTEGER :: min33
    REAL :: y33_tl
    REAL :: y26
    REAL :: min3
    INTEGER :: min32
    REAL :: min11_tl
    REAL :: y25
    REAL :: min2
    INTEGER :: min31
    REAL :: z5_tl
    REAL :: x14_tl
    REAL :: y24
    REAL :: min1
    REAL :: min30
    REAL :: max8_tl
    REAL :: y21_tl
    REAL :: y23
    REAL :: y22
    REAL :: y3_tl
    REAL :: min8_tl
    REAL :: y21
    REAL :: y46_tl
    REAL :: y20
    REAL :: min24_tl
    REAL :: y56
    REAL :: x27_tl
    REAL :: x9
    REAL :: y55
    REAL :: y34_tl
    REAL :: x8
    REAL :: y54
    REAL :: min12_tl
    REAL :: x7
    REAL :: y53
    REAL :: z6_tl
    REAL :: x15_tl
    REAL :: x6
    REAL :: y52
    REAL :: max9_tl
    REAL :: y22_tl
    REAL :: x5
    REAL :: y51
    REAL :: x4
    REAL :: y50
    REAL :: y4_tl
    REAL :: min9_tl
    REAL :: x3
    REAL :: y10_tl
    REAL :: y47_tl
    REAL :: x2
    REAL :: min25_tl
    REAL :: x2_tl
    REAL :: x28_tl
    REAL :: y35_tl
    REAL :: min13_tl
    REAL :: x16_tl
    REAL :: z7_tl
    REAL :: y23_tl
    REAL :: y5_tl
    REAL :: y11_tl
    REAL :: y48_tl
    REAL :: x30_tl
    INTEGER :: min29
    REAL :: min26_tl
    INTEGER :: min28
    REAL :: x3_tl
    REAL :: x29_tl
    REAL :: min27
    REAL :: y36_tl
    REAL :: min26
    REAL :: min14_tl
    REAL :: y19
    REAL :: min25
    REAL :: x17_tl
    REAL :: z8_tl
    REAL :: y18
    REAL :: min24
    REAL :: y24_tl
    REAL :: y17
    REAL :: min23
    REAL :: y50_tl
    REAL :: y16
    REAL :: x35
    INTEGER :: min22
    REAL :: y6_tl
    REAL :: y15
    REAL :: x34
    REAL :: min21
    REAL :: y12_tl
    REAL :: y49_tl
    REAL :: x31_tl
    REAL :: y14
    REAL :: x33
    REAL :: min20
    REAL :: min27_tl
    REAL :: y13
    REAL :: x32
    REAL :: x4_tl
    REAL :: y12
    REAL :: x31
    REAL :: y49
    REAL :: y37_tl
    REAL :: y11
    REAL :: x30
    REAL :: y48
    REAL :: min15_tl
    REAL :: y10
    REAL :: y47
    REAL :: max10_tl
    REAL :: x18_tl
    REAL :: y46
    REAL :: y25_tl
    REAL :: y45
    REAL :: y51_tl
    REAL :: y44
    REAL :: y7_tl
    REAL :: y43
    REAL :: min1_tl
    REAL :: y13_tl
    REAL :: x32_tl
    REAL :: y42
    REAL :: y41
    REAL :: x5_tl
    REAL :: y40
    REAL :: x20_tl
    REAL :: y38_tl
    REAL :: min16_tl
    REAL :: max11_tl
    REAL :: x19_tl
    REAL :: z8
    REAL :: y26_tl
    REAL :: z7
    REAL :: y52_tl
    REAL :: z6
    REAL :: y8_tl
    REAL :: min30_tl
    REAL :: z5
    REAL :: min2_tl
    REAL :: y14_tl
    REAL :: x33_tl
    REAL :: z4
    REAL :: y40_tl
    REAL :: z3
    REAL :: x6_tl
    REAL :: z2
    REAL :: x21_tl
    REAL :: y39_tl
    REAL :: z1
    REAL :: min17_tl
    REAL :: min19
    REAL :: max12_tl
    REAL :: min18
    REAL :: y27_tl
    REAL :: min17
    REAL :: y53_tl
    REAL :: x29
    REAL :: min16
    REAL :: y9_tl
    REAL :: x28
    REAL :: min15
    REAL :: min3_tl
    REAL :: y15_tl
    REAL :: x34_tl
    REAL :: x27
    REAL :: min14
    REAL :: y41_tl
    REAL :: x26
    REAL :: min13
    REAL :: x7_tl
    REAL :: x25
    REAL :: min12
    REAL :: x22_tl
    REAL :: x24
    REAL :: min11
    REAL :: min18_tl
    REAL :: x23
    REAL :: min10
    REAL :: max13_tl
    REAL :: x22
    REAL :: z1_tl
    REAL :: x10_tl
    REAL :: y28_tl
    REAL :: x21
    REAL :: y39
    REAL :: y54_tl
    REAL :: x20
    REAL :: y38
    REAL :: y37
    REAL :: min4_tl
    REAL :: y16_tl
    REAL :: x35_tl
    REAL :: y36
    REAL :: y42_tl
    REAL :: y35
    REAL :: x8_tl
    REAL :: min20_tl
    REAL :: y34
    REAL :: x23_tl
    REAL :: y33
    REAL :: max9
    REAL :: y30_tl
    REAL :: min19_tl
    REAL :: y32
    REAL :: max8
    REAL :: max14_tl
    REAL :: y31
    INTEGER :: max7
    REAL :: z2_tl
    REAL :: x11_tl
    REAL :: y29_tl
    REAL :: y30
    INTEGER :: max6
    REAL :: y55_tl
    INTEGER :: max5
    REAL :: y9
    INTEGER :: max4
    REAL :: max15
    REAL :: min5_tl
    REAL :: y17_tl
    REAL :: y8
    INTEGER :: max3
    REAL :: max14
    REAL :: y43_tl
    REAL :: y7
    INTEGER :: max2
    REAL :: max13
    REAL :: min21_tl
    REAL :: x9_tl
    REAL :: y6
    INTEGER :: max1
    REAL :: max12
    REAL :: x24_tl
    REAL :: y5
    REAL :: max11
    REAL :: y31_tl
    REAL :: y4
    REAL :: max10
    REAL :: max15_tl
    REAL :: y3
    REAL :: z3_tl
    REAL :: x12_tl
    REAL :: y2
    REAL :: y56_tl
    REAL :: y1
    REAL :: cfl,cfl_tl
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
    IF (jord .LE. 4) THEN
      dm_tl = 0.0
      DO j=js-2,je+2
        DO i=ifirst,ilast
          xt_tl = 0.25*(q_tl(i, j+1)-q_tl(i, j-1))
          xt = 0.25*(q(i, j+1)-q(i, j-1))
          IF (xt .GE. 0.) THEN
            x2_tl = xt_tl
            x2 = xt
          ELSE
            x2_tl = -xt_tl
            x2 = -xt
          END IF
          IF (q(i, j-1) .LT. q(i, j)) THEN
            IF (q(i, j) .LT. q(i, j+1)) THEN
              max8_tl = q_tl(i, j+1)
              max8 = q(i, j+1)
            ELSE
              max8_tl = q_tl(i, j)
              max8 = q(i, j)
            END IF
          ELSE IF (q(i, j-1) .LT. q(i, j+1)) THEN
            max8_tl = q_tl(i, j+1)
            max8 = q(i, j+1)
          ELSE
            max8_tl = q_tl(i, j-1)
            max8 = q(i, j-1)
          END IF
          y1_tl = max8_tl - q_tl(i, j)
          y1 = max8 - q(i, j)
          IF (q(i, j-1) .GT. q(i, j)) THEN
            IF (q(i, j) .GT. q(i, j+1)) THEN
              min21_tl = q_tl(i, j+1)
              min21 = q(i, j+1)
            ELSE
              min21_tl = q_tl(i, j)
              min21 = q(i, j)
            END IF
          ELSE IF (q(i, j-1) .GT. q(i, j+1)) THEN
            min21_tl = q_tl(i, j+1)
            min21 = q(i, j+1)
          ELSE
            min21_tl = q_tl(i, j-1)
            min21 = q(i, j-1)
          END IF
          z1_tl = q_tl(i, j) - min21_tl
          z1 = q(i, j) - min21
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
        IF (3 .LT. js - 1) THEN
          max1 = js - 1
        ELSE
          max1 = 3
        END IF
        IF (npy - 2 .GT. je + 2) THEN
          min22 = je + 2
          al_tl = 0.0
        ELSE
          min22 = npy - 2
          al_tl = 0.0
        END IF
        DO j=max1,min22
          DO i=ifirst,ilast
            al_tl(i, j) = 0.5*(q_tl(i, j-1)+q_tl(i, j)) + r3*(dm_tl(i, j&
&             -1)-dm_tl(i, j))
            al(i, j) = 0.5*(q(i, j-1)+q(i, j)) + r3*(dm(i, j-1)-dm(i, j)&
&             )
          END DO
        END DO
!--------------
! Fix the edges:
!--------------
        IF (js .EQ. 1) THEN
          DO i=ifirst,ilast
            x0_tl = 0.5*((2.*dya(i, 1)+dya(i, 2))*(q_tl(i, 0)+q_tl(i, 1)&
&             )-dya(i, 1)*(q_tl(i, -1)+q_tl(i, 2)))/(dya(i, 1)+dya(i, 2)&
&             )
            x0 = 0.5*((2.*dya(i, 1)+dya(i, 2))*(q(i, 0)+q(i, 1))-dya(i, &
&             1)*(q(i, -1)+q(i, 2)))/(dya(i, 1)+dya(i, 2))
            al_tl(i, 1) = x0_tl
            al(i, 1) = x0
            x1_tl = s15*q_tl(i, 0) + s11*q_tl(i, -1) + s14*dm_tl(i, -1)
            x1 = s15*q(i, 0) + s11*q(i, -1) + s14*dm(i, -1)
            dm_tl(i, 0) = 0.5*(x0_tl-x1_tl)
            dm(i, 0) = 0.5*(x0-x1)
            IF (dm(i, 0) .GE. 0.) THEN
              x3_tl = dm_tl(i, 0)
              x3 = dm(i, 0)
            ELSE
              x3_tl = -dm_tl(i, 0)
              x3 = -dm(i, 0)
            END IF
            IF (q(i, 0) .LT. x0) THEN
              IF (x0 .LT. x1) THEN
                max9_tl = x1_tl
                max9 = x1
              ELSE
                max9_tl = x0_tl
                max9 = x0
              END IF
            ELSE IF (q(i, 0) .LT. x1) THEN
              max9_tl = x1_tl
              max9 = x1
            ELSE
              max9_tl = q_tl(i, 0)
              max9 = q(i, 0)
            END IF
            y2_tl = max9_tl - q_tl(i, 0)
            y2 = max9 - q(i, 0)
            IF (q(i, 0) .GT. x0) THEN
              IF (x0 .GT. x1) THEN
                min23_tl = x1_tl
                min23 = x1
              ELSE
                min23_tl = x0_tl
                min23 = x0
              END IF
            ELSE IF (q(i, 0) .GT. x1) THEN
              min23_tl = x1_tl
              min23 = x1
            ELSE
              min23_tl = q_tl(i, 0)
              min23 = q(i, 0)
            END IF
            z2_tl = q_tl(i, 0) - min23_tl
            z2 = q(i, 0) - min23
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
            dm_tl(i, 0) = min2_tl*SIGN(1.d0, min2*dm(i, 0))
            dm(i, 0) = SIGN(min2, dm(i, 0))
            al_tl(i, 0) = 0.5*(q_tl(i, -1)+q_tl(i, 0)) + r3*(dm_tl(i, -1&
&             )-dm_tl(i, 0))
            al(i, 0) = 0.5*(q(i, -1)+q(i, 0)) + r3*(dm(i, -1)-dm(i, 0))
!
            x1_tl = s15*q_tl(i, 1) + s11*q_tl(i, 2) - s14*dm_tl(i, 2)
            x1 = s15*q(i, 1) + s11*q(i, 2) - s14*dm(i, 2)
            dm_tl(i, 1) = 0.5*(x1_tl-x0_tl)
            dm(i, 1) = 0.5*(x1-x0)
            IF (dm(i, 1) .GE. 0.) THEN
              x4_tl = dm_tl(i, 1)
              x4 = dm(i, 1)
            ELSE
              x4_tl = -dm_tl(i, 1)
              x4 = -dm(i, 1)
            END IF
            IF (q(i, 1) .LT. x0) THEN
              IF (x0 .LT. x1) THEN
                max10_tl = x1_tl
                max10 = x1
              ELSE
                max10_tl = x0_tl
                max10 = x0
              END IF
            ELSE IF (q(i, 1) .LT. x1) THEN
              max10_tl = x1_tl
              max10 = x1
            ELSE
              max10_tl = q_tl(i, 1)
              max10 = q(i, 1)
            END IF
            y3_tl = max10_tl - q_tl(i, 1)
            y3 = max10 - q(i, 1)
            IF (q(i, 1) .GT. x0) THEN
              IF (x0 .GT. x1) THEN
                min24_tl = x1_tl
                min24 = x1
              ELSE
                min24_tl = x0_tl
                min24 = x0
              END IF
            ELSE IF (q(i, 1) .GT. x1) THEN
              min24_tl = x1_tl
              min24 = x1
            ELSE
              min24_tl = q_tl(i, 1)
              min24 = q(i, 1)
            END IF
            z3_tl = q_tl(i, 1) - min24_tl
            z3 = q(i, 1) - min24
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
            dm_tl(i, 1) = min3_tl*SIGN(1.d0, min3*dm(i, 1))
            dm(i, 1) = SIGN(min3, dm(i, 1))
            al_tl(i, 2) = 0.5*(q_tl(i, 1)+q_tl(i, 2)) + r3*(dm_tl(i, 1)-&
&             dm_tl(i, 2))
            al(i, 2) = 0.5*(q(i, 1)+q(i, 2)) + r3*(dm(i, 1)-dm(i, 2))
          END DO
        END IF
        IF (je + 1 .EQ. npy) THEN
          DO i=ifirst,ilast
            x0_tl = 0.5*((2.*dya(i, npy-1)+dya(i, npy-2))*(q_tl(i, npy-1&
&             )+q_tl(i, npy))-dya(i, npy-1)*(q_tl(i, npy-2)+q_tl(i, npy+&
&             1)))/(dya(i, npy-1)+dya(i, npy-2))
            x0 = 0.5*((2.*dya(i, npy-1)+dya(i, npy-2))*(q(i, npy-1)+q(i&
&             , npy))-dya(i, npy-1)*(q(i, npy-2)+q(i, npy+1)))/(dya(i, &
&             npy-1)+dya(i, npy-2))
            al_tl(i, npy) = x0_tl
            al(i, npy) = x0
            x1_tl = s15*q_tl(i, npy-1) + s11*q_tl(i, npy-2) + s14*dm_tl(&
&             i, npy-2)
            x1 = s15*q(i, npy-1) + s11*q(i, npy-2) + s14*dm(i, npy-2)
            dm_tl(i, npy-1) = 0.5*(x0_tl-x1_tl)
            dm(i, npy-1) = 0.5*(x0-x1)
            IF (dm(i, npy-1) .GE. 0.) THEN
              x5_tl = dm_tl(i, npy-1)
              x5 = dm(i, npy-1)
            ELSE
              x5_tl = -dm_tl(i, npy-1)
              x5 = -dm(i, npy-1)
            END IF
            IF (q(i, npy-1) .LT. x0) THEN
              IF (x0 .LT. x1) THEN
                max11_tl = x1_tl
                max11 = x1
              ELSE
                max11_tl = x0_tl
                max11 = x0
              END IF
            ELSE IF (q(i, npy-1) .LT. x1) THEN
              max11_tl = x1_tl
              max11 = x1
            ELSE
              max11_tl = q_tl(i, npy-1)
              max11 = q(i, npy-1)
            END IF
            y4_tl = max11_tl - q_tl(i, npy-1)
            y4 = max11 - q(i, npy-1)
            IF (q(i, npy-1) .GT. x0) THEN
              IF (x0 .GT. x1) THEN
                min25_tl = x1_tl
                min25 = x1
              ELSE
                min25_tl = x0_tl
                min25 = x0
              END IF
            ELSE IF (q(i, npy-1) .GT. x1) THEN
              min25_tl = x1_tl
              min25 = x1
            ELSE
              min25_tl = q_tl(i, npy-1)
              min25 = q(i, npy-1)
            END IF
            z4_tl = q_tl(i, npy-1) - min25_tl
            z4 = q(i, npy-1) - min25
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
            al_tl(i, npy-1) = 0.5*(q_tl(i, npy-2)+q_tl(i, npy-1)) + r3*(&
&             dm_tl(i, npy-2)-dm_tl(i, npy-1))
            al(i, npy-1) = 0.5*(q(i, npy-2)+q(i, npy-1)) + r3*(dm(i, npy&
&             -2)-dm(i, npy-1))
!
            x1_tl = s15*q_tl(i, npy) + s11*q_tl(i, npy+1) - s14*dm_tl(i&
&             , npy+1)
            x1 = s15*q(i, npy) + s11*q(i, npy+1) - s14*dm(i, npy+1)
            dm_tl(i, npy) = 0.5*(x1_tl-x0_tl)
            dm(i, npy) = 0.5*(x1-x0)
            IF (dm(i, npy) .GE. 0.) THEN
              x6_tl = dm_tl(i, npy)
              x6 = dm(i, npy)
            ELSE
              x6_tl = -dm_tl(i, npy)
              x6 = -dm(i, npy)
            END IF
            IF (q(i, npy) .LT. x0) THEN
              IF (x0 .LT. x1) THEN
                max12_tl = x1_tl
                max12 = x1
              ELSE
                max12_tl = x0_tl
                max12 = x0
              END IF
            ELSE IF (q(i, npy) .LT. x1) THEN
              max12_tl = x1_tl
              max12 = x1
            ELSE
              max12_tl = q_tl(i, npy)
              max12 = q(i, npy)
            END IF
            y5_tl = max12_tl - q_tl(i, npy)
            y5 = max12 - q(i, npy)
            IF (q(i, npy) .GT. x0) THEN
              IF (x0 .GT. x1) THEN
                min26_tl = x1_tl
                min26 = x1
              ELSE
                min26_tl = x0_tl
                min26 = x0
              END IF
            ELSE IF (q(i, npy) .GT. x1) THEN
              min26_tl = x1_tl
              min26 = x1
            ELSE
              min26_tl = q_tl(i, npy)
              min26 = q(i, npy)
            END IF
            z5_tl = q_tl(i, npy) - min26_tl
            z5 = q(i, npy) - min26
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
            al_tl(i, npy+1) = 0.5*(q_tl(i, npy)+q_tl(i, npy+1)) + r3*(&
&             dm_tl(i, npy)-dm_tl(i, npy+1))
            al(i, npy+1) = 0.5*(q(i, npy)+q(i, npy+1)) + r3*(dm(i, npy)-&
&             dm(i, npy+1))
          END DO
        END IF
      ELSE
        al_tl = 0.0
! Doubly periodic BC:
        DO j=js-1,je+2
          DO i=ifirst,ilast
            al_tl(i, j) = 0.5*(q_tl(i, j-1)+q_tl(i, j)) + r3*(dm_tl(i, j&
&             -1)-dm_tl(i, j))
            al(i, j) = 0.5*(q(i, j-1)+q(i, j)) + r3*(dm(i, j-1)-dm(i, j)&
&             )
          END DO
        END DO
      END IF
      IF (jord .EQ. 3) THEN
        bl_tl = 0.0
        br_tl = 0.0
        DO j=js-1,je+1
          DO i=ifirst,ilast
            bl_tl(i, j) = al_tl(i, j) - q_tl(i, j)
            bl(i, j) = al(i, j) - q(i, j)
            br_tl(i, j) = al_tl(i, j+1) - q_tl(i, j)
            br(i, j) = al(i, j+1) - q(i, j)
          END DO
          CALL PERT_PPM_TLM(ilast - ifirst + 1, q(ifirst:, j), bl(ifirst&
&                     :, j), bl_tl(ifirst:, j), br(ifirst:, j), br_tl(&
&                     ifirst:, j), 1)
        END DO
        DO j=js,je+1
          DO i=ifirst,ilast
            IF (c(i, j) .GT. 0.) THEN
              flux_tl(i, j) = q_tl(i, j-1) + (1.-c(i, j))*(br_tl(i, j-1)&
&               -c_tl(i, j)*(bl(i, j-1)+br(i, j-1))-c(i, j)*(bl_tl(i, j-&
&               1)+br_tl(i, j-1))) - c_tl(i, j)*(br(i, j-1)-c(i, j)*(bl(&
&               i, j-1)+br(i, j-1)))
              flux(i, j) = q(i, j-1) + (1.-c(i, j))*(br(i, j-1)-c(i, j)*&
&               (bl(i, j-1)+br(i, j-1)))
            ELSE
              flux_tl(i, j) = q_tl(i, j) + c_tl(i, j)*(bl(i, j)+c(i, j)*&
&               (bl(i, j)+br(i, j))) + (1.+c(i, j))*(bl_tl(i, j)+c_tl(i&
&               , j)*(bl(i, j)+br(i, j))+c(i, j)*(bl_tl(i, j)+br_tl(i, j&
&               )))
              flux(i, j) = q(i, j) + (1.+c(i, j))*(bl(i, j)+c(i, j)*(bl(&
&               i, j)+br(i, j)))
            END IF
          END DO
        END DO
      ELSE
! Inlined limiter
        DO j=js,je+1
          DO i=ifirst,ilast
            IF (c(i, j) .GT. 0.) THEN
              xt_tl = ppm_limiter*dm_tl(i, j-1)
              xt = ppm_limiter*dm(i, j-1)
              IF (xt .GE. 0.) THEN
                x7_tl = xt_tl
                x7 = xt
              ELSE
                x7_tl = -xt_tl
                x7 = -xt
              END IF
              IF (al(i, j-1) - q(i, j-1) .GE. 0.) THEN
                y6_tl = al_tl(i, j-1) - q_tl(i, j-1)
                y6 = al(i, j-1) - q(i, j-1)
              ELSE
                y6_tl = -(al_tl(i, j-1)-q_tl(i, j-1))
                y6 = -(al(i, j-1)-q(i, j-1))
              END IF
              IF (x7 .GT. y6) THEN
                min6_tl = y6_tl
                min6 = y6
              ELSE
                min6_tl = x7_tl
                min6 = x7
              END IF
              dl_tl = min6_tl*SIGN(1.d0, min6*xt)
              dl = SIGN(min6, xt)
              IF (xt .GE. 0.) THEN
                x8_tl = xt_tl
                x8 = xt
              ELSE
                x8_tl = -xt_tl
                x8 = -xt
              END IF
              IF (al(i, j) - q(i, j-1) .GE. 0.) THEN
                y7_tl = al_tl(i, j) - q_tl(i, j-1)
                y7 = al(i, j) - q(i, j-1)
              ELSE
                y7_tl = -(al_tl(i, j)-q_tl(i, j-1))
                y7 = -(al(i, j)-q(i, j-1))
              END IF
              IF (x8 .GT. y7) THEN
                min7_tl = y7_tl
                min7 = y7
              ELSE
                min7_tl = x8_tl
                min7 = x8
              END IF
              dr_tl = min7_tl*SIGN(1.d0, min7*xt)
              dr = SIGN(min7, xt)
              flux_tl(i, j) = q_tl(i, j-1) + (1.-c(i, j))*(c_tl(i, j)*(&
&               dl-dr)+c(i, j)*(dl_tl-dr_tl)+dr_tl) - c_tl(i, j)*(c(i, j&
&               )*(dl-dr)+dr)
              flux(i, j) = q(i, j-1) + (1.-c(i, j))*(c(i, j)*(dl-dr)+dr)
            ELSE
              xt_tl = ppm_limiter*dm_tl(i, j)
              xt = ppm_limiter*dm(i, j)
              IF (xt .GE. 0.) THEN
                x9_tl = xt_tl
                x9 = xt
              ELSE
                x9_tl = -xt_tl
                x9 = -xt
              END IF
              IF (al(i, j) - q(i, j) .GE. 0.) THEN
                y8_tl = al_tl(i, j) - q_tl(i, j)
                y8 = al(i, j) - q(i, j)
              ELSE
                y8_tl = -(al_tl(i, j)-q_tl(i, j))
                y8 = -(al(i, j)-q(i, j))
              END IF
              IF (x9 .GT. y8) THEN
                min8_tl = y8_tl
                min8 = y8
              ELSE
                min8_tl = x9_tl
                min8 = x9
              END IF
              dl_tl = min8_tl*SIGN(1.d0, min8*xt)
              dl = SIGN(min8, xt)
              IF (xt .GE. 0.) THEN
                x10_tl = xt_tl
                x10 = xt
              ELSE
                x10_tl = -xt_tl
                x10 = -xt
              END IF
              IF (al(i, j+1) - q(i, j) .GE. 0.) THEN
                y9_tl = al_tl(i, j+1) - q_tl(i, j)
                y9 = al(i, j+1) - q(i, j)
              ELSE
                y9_tl = -(al_tl(i, j+1)-q_tl(i, j))
                y9 = -(al(i, j+1)-q(i, j))
              END IF
              IF (x10 .GT. y9) THEN
                min9_tl = y9_tl
                min9 = y9
              ELSE
                min9_tl = x10_tl
                min9 = x10
              END IF
              dr_tl = min9_tl*SIGN(1.d0, min9*xt)
              dr = SIGN(min9, xt)
              flux_tl(i, j) = q_tl(i, j) - c_tl(i, j)*(c(i, j)*(dl-dr)+&
&               dl) - (1.+c(i, j))*(c_tl(i, j)*(dl-dr)+c(i, j)*(dl_tl-&
&               dr_tl)+dl_tl)
              flux(i, j) = q(i, j) - (1.+c(i, j))*(c(i, j)*(dl-dr)+dl)
            END IF
          END DO
        END DO
      END IF
    ELSE IF (jord .EQ. 5) THEN
      dq_tl = 0.0
! PPM with Hunyh's 2nd constraint
      DO j=jfirst-3,jlast+2
        DO i=ifirst,ilast
          dq_tl(i, j) = q_tl(i, j+1) - q_tl(i, j)
          dq(i, j) = q(i, j+1) - q(i, j)
        END DO
      END DO
      dm_tl = 0.0
      DO j=jfirst-2,jlast+2
        DO i=ifirst,ilast
          xt_tl = 0.25*(q_tl(i, j+1)-q_tl(i, j-1))
          xt = 0.25*(q(i, j+1)-q(i, j-1))
          IF (xt .GE. 0.) THEN
            x11_tl = xt_tl
            x11 = xt
          ELSE
            x11_tl = -xt_tl
            x11 = -xt
          END IF
          IF (q(i, j-1) .LT. q(i, j)) THEN
            IF (q(i, j) .LT. q(i, j+1)) THEN
              max13_tl = q_tl(i, j+1)
              max13 = q(i, j+1)
            ELSE
              max13_tl = q_tl(i, j)
              max13 = q(i, j)
            END IF
          ELSE IF (q(i, j-1) .LT. q(i, j+1)) THEN
            max13_tl = q_tl(i, j+1)
            max13 = q(i, j+1)
          ELSE
            max13_tl = q_tl(i, j-1)
            max13 = q(i, j-1)
          END IF
          y10_tl = max13_tl - q_tl(i, j)
          y10 = max13 - q(i, j)
          IF (q(i, j-1) .GT. q(i, j)) THEN
            IF (q(i, j) .GT. q(i, j+1)) THEN
              min27_tl = q_tl(i, j+1)
              min27 = q(i, j+1)
            ELSE
              min27_tl = q_tl(i, j)
              min27 = q(i, j)
            END IF
          ELSE IF (q(i, j-1) .GT. q(i, j+1)) THEN
            min27_tl = q_tl(i, j+1)
            min27 = q(i, j+1)
          ELSE
            min27_tl = q_tl(i, j-1)
            min27 = q(i, j-1)
          END IF
          z6_tl = q_tl(i, j) - min27_tl
          z6 = q(i, j) - min27
          IF (x11 .GT. y10) THEN
            IF (y10 .GT. z6) THEN
              min10_tl = z6_tl
              min10 = z6
            ELSE
              min10_tl = y10_tl
              min10 = y10
            END IF
          ELSE IF (x11 .GT. z6) THEN
            min10_tl = z6_tl
            min10 = z6
          ELSE
            min10_tl = x11_tl
            min10 = x11
          END IF
          dm_tl(i, j) = min10_tl*SIGN(1.d0, min10*xt)
          dm(i, j) = SIGN(min10, xt)
        END DO
      END DO
      al_tl = 0.0
      DO j=jfirst-1,jlast+2
        DO i=ifirst,ilast
          al_tl(i, j) = 0.5*(q_tl(i, j-1)+q_tl(i, j)) + r3*(dm_tl(i, j-1&
&           )-dm_tl(i, j))
          al(i, j) = 0.5*(q(i, j-1)+q(i, j)) + r3*(dm(i, j-1)-dm(i, j))
        END DO
      END DO
      bl_tl = 0.0
      br_tl = 0.0
      DO j=jfirst-1,jlast+1
        DO i=ifirst,ilast
          pmp_tl = -(2.*dq_tl(i, j))
          pmp = -(2.*dq(i, j))
          lac_tl = pmp_tl + 1.5*dq_tl(i, j+1)
          lac = pmp + 1.5*dq(i, j+1)
          IF (0. .LT. pmp) THEN
            IF (pmp .LT. lac) THEN
              x12_tl = lac_tl
              x12 = lac
            ELSE
              x12_tl = pmp_tl
              x12 = pmp
            END IF
          ELSE IF (0. .LT. lac) THEN
            x12_tl = lac_tl
            x12 = lac
          ELSE
            x12 = 0.
            x12_tl = 0.0
          END IF
          IF (0. .GT. pmp) THEN
            IF (pmp .GT. lac) THEN
              y43_tl = lac_tl
              y43 = lac
            ELSE
              y43_tl = pmp_tl
              y43 = pmp
            END IF
          ELSE IF (0. .GT. lac) THEN
            y43_tl = lac_tl
            y43 = lac
          ELSE
            y43 = 0.
            y43_tl = 0.0
          END IF
          IF (al(i, j) - q(i, j) .LT. y43) THEN
            y11_tl = y43_tl
            y11 = y43
          ELSE
            y11_tl = al_tl(i, j) - q_tl(i, j)
            y11 = al(i, j) - q(i, j)
          END IF
          IF (x12 .GT. y11) THEN
            bl_tl(i, j) = y11_tl
            bl(i, j) = y11
          ELSE
            bl_tl(i, j) = x12_tl
            bl(i, j) = x12
          END IF
          pmp_tl = 2.*dq_tl(i, j-1)
          pmp = 2.*dq(i, j-1)
          lac_tl = pmp_tl - 1.5*dq_tl(i, j-2)
          lac = pmp - 1.5*dq(i, j-2)
          IF (0. .LT. pmp) THEN
            IF (pmp .LT. lac) THEN
              x13_tl = lac_tl
              x13 = lac
            ELSE
              x13_tl = pmp_tl
              x13 = pmp
            END IF
          ELSE IF (0. .LT. lac) THEN
            x13_tl = lac_tl
            x13 = lac
          ELSE
            x13 = 0.
            x13_tl = 0.0
          END IF
          IF (0. .GT. pmp) THEN
            IF (pmp .GT. lac) THEN
              y44_tl = lac_tl
              y44 = lac
            ELSE
              y44_tl = pmp_tl
              y44 = pmp
            END IF
          ELSE IF (0. .GT. lac) THEN
            y44_tl = lac_tl
            y44 = lac
          ELSE
            y44 = 0.
            y44_tl = 0.0
          END IF
          IF (al(i, j+1) - q(i, j) .LT. y44) THEN
            y12_tl = y44_tl
            y12 = y44
          ELSE
            y12_tl = al_tl(i, j+1) - q_tl(i, j)
            y12 = al(i, j+1) - q(i, j)
          END IF
          IF (x13 .GT. y12) THEN
            br_tl(i, j) = y12_tl
            br(i, j) = y12
          ELSE
            br_tl(i, j) = x13_tl
            br(i, j) = x13
          END IF
        END DO
      END DO
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

    ELSE IF (jord .EQ. 6) THEN

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
      IF (grid_type .LT. 3) THEN
        IF (js .EQ. 1) THEN
          DO i=ifirst,ilast
            br_tl(i, 2) = p1*(q_tl(i, 2)+q_tl(i, 3)) + p2*(q_tl(i, 1)+&
&             q_tl(i, 4)) - q_tl(i, 2)
            br(i, 2) = p1*(q(i, 2)+q(i, 3)) + p2*(q(i, 1)+q(i, 4)) - q(i&
&             , 2)
            xt_tl = c3*q_tl(i, 1) + c2*q_tl(i, 2) + c1*q_tl(i, 3)
            xt = c3*q(i, 1) + c2*q(i, 2) + c1*q(i, 3)
            br_tl(i, 1) = xt_tl - q_tl(i, 1)
            br(i, 1) = xt - q(i, 1)
            bl_tl(i, 2) = xt_tl - q_tl(i, 2)
            bl(i, 2) = xt - q(i, 2)
            bl_tl(i, 0) = c1*q_tl(i, -2) + c2*q_tl(i, -1) + c3*q_tl(i, 0&
&             ) - q_tl(i, 0)
            bl(i, 0) = c1*q(i, -2) + c2*q(i, -1) + c3*q(i, 0) - q(i, 0)
            xt_tl = 0.5*((2.*dy(i, 1)+dy(i, 2))*(q_tl(i, 0)+q_tl(i, 1))-&
&             dy(i, 1)*(q_tl(i, -1)+q_tl(i, 2)))/(dy(i, 1)+dy(i, 2))
            xt = 0.5*((2.*dy(i, 1)+dy(i, 2))*(q(i, 0)+q(i, 1))-dy(i, 1)*&
&             (q(i, -1)+q(i, 2)))/(dy(i, 1)+dy(i, 2))
            bl_tl(i, 1) = xt_tl - q_tl(i, 1)
            bl(i, 1) = xt - q(i, 1)
            br_tl(i, 0) = xt_tl - q_tl(i, 0)
            br(i, 0) = xt - q(i, 0)
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
          DO i=ifirst,ilast
            bl_tl(i, npy-2) = p1*(q_tl(i, npy-3)+q_tl(i, npy-2)) + p2*(&
&             q_tl(i, npy-4)+q_tl(i, npy-1)) - q_tl(i, npy-2)
            bl(i, npy-2) = p1*(q(i, npy-3)+q(i, npy-2)) + p2*(q(i, npy-4&
&             )+q(i, npy-1)) - q(i, npy-2)
            xt_tl = c1*q_tl(i, npy-3) + c2*q_tl(i, npy-2) + c3*q_tl(i, &
&             npy-1)
            xt = c1*q(i, npy-3) + c2*q(i, npy-2) + c3*q(i, npy-1)
            br_tl(i, npy-2) = xt_tl - q_tl(i, npy-2)
            br(i, npy-2) = xt - q(i, npy-2)
            bl_tl(i, npy-1) = xt_tl - q_tl(i, npy-1)
            bl(i, npy-1) = xt - q(i, npy-1)
            br_tl(i, npy) = c3*q_tl(i, npy) + c2*q_tl(i, npy+1) + c1*&
&             q_tl(i, npy+2) - q_tl(i, npy)
            br(i, npy) = c3*q(i, npy) + c2*q(i, npy+1) + c1*q(i, npy+2) &
&             - q(i, npy)
            xt_tl = 0.5*((2.*dy(i, npy-1)+dy(i, npy-2))*(q_tl(i, npy-1)+&
&             q_tl(i, npy))-dy(i, npy-1)*(q_tl(i, npy-2)+q_tl(i, npy+1))&
&             )/(dy(i, npy-1)+dy(i, npy-2))
            xt = 0.5*((2.*dy(i, npy-1)+dy(i, npy-2))*(q(i, npy-1)+q(i, &
&             npy))-dy(i, npy-1)*(q(i, npy-2)+q(i, npy+1)))/(dy(i, npy-1&
&             )+dy(i, npy-2))
            br_tl(i, npy-1) = xt_tl - q_tl(i, npy-1)
            br(i, npy-1) = xt - q(i, npy-1)
            bl_tl(i, npy) = xt_tl - q_tl(i, npy)
            bl(i, npy) = xt - q(i, npy)
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
      DO j=jfirst,jlast+1
        DO i=ifirst,ilast
          IF (c(i, j) .GT. 0.) THEN
            cfl_tl = c_tl(i, j)
            cfl = c(i, j)
            flux_tl(i, j) = q_tl(i, j-1) + (1.-cfl)*(br_tl(i, j-1)-&
&             cfl_tl*(bl(i, j-1)+br(i, j-1))-cfl*(bl_tl(i, j-1)+br_tl(i&
&             , j-1))) - cfl_tl*(br(i, j-1)-cfl*(bl(i, j-1)+br(i, j-1)))
            flux(i, j) = q(i, j-1) + (1.-cfl)*(br(i, j-1)-cfl*(bl(i, j-1&
&             )+br(i, j-1)))
          ELSE
            cfl_tl = c_tl(i, j)
            cfl = c(i, j)
            flux_tl(i, j) = q_tl(i, j) + cfl_tl*(bl(i, j)+cfl*(bl(i, j)+&
&             br(i, j))) + (1.+cfl)*(bl_tl(i, j)+cfl_tl*(bl(i, j)+br(i, &
&             j))+cfl*(bl_tl(i, j)+br_tl(i, j)))
            flux(i, j) = q(i, j) + (1.+cfl)*(bl(i, j)+cfl*(bl(i, j)+br(i&
&             , j)))
          END IF
        END DO
      END DO


    ELSE IF (jord .GE. 8 .and. jord .LE. 10) THEN
      dm_tl = 0.0
! jord=8, 9, 10
      DO j=js-2,je+2
        DO i=ifirst,ilast
          xt_tl = 0.25*(q_tl(i, j+1)-q_tl(i, j-1))
          xt = 0.25*(q(i, j+1)-q(i, j-1))
          IF (xt .GE. 0.) THEN
            x18_tl = xt_tl
            x18 = xt
          ELSE
            x18_tl = -xt_tl
            x18 = -xt
          END IF
          IF (q(i, j-1) .LT. q(i, j)) THEN
            IF (q(i, j) .LT. q(i, j+1)) THEN
              max14_tl = q_tl(i, j+1)
              max14 = q(i, j+1)
            ELSE
              max14_tl = q_tl(i, j)
              max14 = q(i, j)
            END IF
          ELSE IF (q(i, j-1) .LT. q(i, j+1)) THEN
            max14_tl = q_tl(i, j+1)
            max14 = q(i, j+1)
          ELSE
            max14_tl = q_tl(i, j-1)
            max14 = q(i, j-1)
          END IF
          y25_tl = max14_tl - q_tl(i, j)
          y25 = max14 - q(i, j)
          IF (q(i, j-1) .GT. q(i, j)) THEN
            IF (q(i, j) .GT. q(i, j+1)) THEN
              min30_tl = q_tl(i, j+1)
              min30 = q(i, j+1)
            ELSE
              min30_tl = q_tl(i, j)
              min30 = q(i, j)
            END IF
          ELSE IF (q(i, j-1) .GT. q(i, j+1)) THEN
            min30_tl = q_tl(i, j+1)
            min30 = q(i, j+1)
          ELSE
            min30_tl = q_tl(i, j-1)
            min30 = q(i, j-1)
          END IF
          z7_tl = q_tl(i, j) - min30_tl
          z7 = q(i, j) - min30
          IF (x18 .GT. y25) THEN
            IF (y25 .GT. z7) THEN
              min13_tl = z7_tl
              min13 = z7
            ELSE
              min13_tl = y25_tl
              min13 = y25
            END IF
          ELSE IF (x18 .GT. z7) THEN
            min13_tl = z7_tl
            min13 = z7
          ELSE
            min13_tl = x18_tl
            min13 = x18
          END IF
          dm_tl(i, j) = min13_tl*SIGN(1.d0, min13*xt)
          dm(i, j) = SIGN(min13, xt)
        END DO
      END DO
      IF (grid_type .LT. 3) THEN
        IF (3 .LT. js - 1) THEN
          max4 = js - 1
        ELSE
          max4 = 3
        END IF
        IF (npy - 2 .GT. je + 2) THEN
          min31 = je + 2
          al_tl = 0.0
        ELSE
          min31 = npy - 2
          al_tl = 0.0
        END IF
        DO j=max4,min31
          DO i=ifirst,ilast
            al_tl(i, j) = 0.5*(q_tl(i, j-1)+q_tl(i, j)) + r3*(dm_tl(i, j&
&             -1)-dm_tl(i, j))
            al(i, j) = 0.5*(q(i, j-1)+q(i, j)) + r3*(dm(i, j-1)-dm(i, j)&
&             )
          END DO
        END DO
        dq_tl = 0.0
        DO j=js-3,je+2
          DO i=ifirst,ilast
            dq_tl(i, j) = q_tl(i, j+1) - q_tl(i, j)
            dq(i, j) = q(i, j+1) - q(i, j)
          END DO
        END DO
        IF (jord .EQ. 8) THEN
          IF (3 .LT. js - 1) THEN
            max5 = js - 1
          ELSE
            max5 = 3
          END IF
          IF (npy - 3 .GT. je + 1) THEN
            min32 = je + 1
            bl_tl = 0.0
            br_tl = 0.0
          ELSE
            min32 = npy - 3
            bl_tl = 0.0
            br_tl = 0.0
          END IF
          DO j=max5,min32
            DO i=ifirst,ilast
              xt_tl = 2.*dm_tl(i, j)
              xt = 2.*dm(i, j)
              IF (xt .GE. 0.) THEN
                x19_tl = xt_tl
                x19 = xt
              ELSE
                x19_tl = -xt_tl
                x19 = -xt
              END IF
              IF (al(i, j) - q(i, j) .GE. 0.) THEN
                y26_tl = al_tl(i, j) - q_tl(i, j)
                y26 = al(i, j) - q(i, j)
              ELSE
                y26_tl = -(al_tl(i, j)-q_tl(i, j))
                y26 = -(al(i, j)-q(i, j))
              END IF
              IF (x19 .GT. y26) THEN
                min14_tl = y26_tl
                min14 = y26
              ELSE
                min14_tl = x19_tl
                min14 = x19
              END IF
              bl_tl(i, j) = -(min14_tl*SIGN(1.d0, min14*xt))
              bl(i, j) = -SIGN(min14, xt)
              IF (xt .GE. 0.) THEN
                x20_tl = xt_tl
                x20 = xt
              ELSE
                x20_tl = -xt_tl
                x20 = -xt
              END IF
              IF (al(i, j+1) - q(i, j) .GE. 0.) THEN
                y27_tl = al_tl(i, j+1) - q_tl(i, j)
                y27 = al(i, j+1) - q(i, j)
              ELSE
                y27_tl = -(al_tl(i, j+1)-q_tl(i, j))
                y27 = -(al(i, j+1)-q(i, j))
              END IF
              IF (x20 .GT. y27) THEN
                min15_tl = y27_tl
                min15 = y27
              ELSE
                min15_tl = x20_tl
                min15 = x20
              END IF
              br_tl(i, j) = min15_tl*SIGN(1.d0, min15*xt)
              br(i, j) = SIGN(min15, xt)
            END DO
          END DO
        ELSE IF (jord .EQ. 9) THEN
          IF (3 .LT. js - 1) THEN
            max6 = js - 1
          ELSE
            max6 = 3
          END IF
          IF (npy - 3 .GT. je + 1) THEN
            min33 = je + 1
            bl_tl = 0.0
            br_tl = 0.0
          ELSE
            min33 = npy - 3
            bl_tl = 0.0
            br_tl = 0.0
          END IF
          DO j=max6,min33
            DO i=ifirst,ilast
              pmp_tl = -(2.*dq_tl(i, j))
              pmp = -(2.*dq(i, j))
              lac_tl = pmp_tl + 1.5*dq_tl(i, j+1)
              lac = pmp + 1.5*dq(i, j+1)
              IF (0. .LT. pmp) THEN
                IF (pmp .LT. lac) THEN
                  x21_tl = lac_tl
                  x21 = lac
                ELSE
                  x21_tl = pmp_tl
                  x21 = pmp
                END IF
              ELSE IF (0. .LT. lac) THEN
                x21_tl = lac_tl
                x21 = lac
              ELSE
                x21 = 0.
                x21_tl = 0.0
              END IF
              IF (0. .GT. pmp) THEN
                IF (pmp .GT. lac) THEN
                  y47_tl = lac_tl
                  y47 = lac
                ELSE
                  y47_tl = pmp_tl
                  y47 = pmp
                END IF
              ELSE IF (0. .GT. lac) THEN
                y47_tl = lac_tl
                y47 = lac
              ELSE
                y47 = 0.
                y47_tl = 0.0
              END IF
              IF (al(i, j) - q(i, j) .LT. y47) THEN
                y28_tl = y47_tl
                y28 = y47
              ELSE
                y28_tl = al_tl(i, j) - q_tl(i, j)
                y28 = al(i, j) - q(i, j)
              END IF
              IF (x21 .GT. y28) THEN
                bl_tl(i, j) = y28_tl
                bl(i, j) = y28
              ELSE
                bl_tl(i, j) = x21_tl
                bl(i, j) = x21
              END IF
              pmp_tl = 2.*dq_tl(i, j-1)
              pmp = 2.*dq(i, j-1)
              lac_tl = pmp_tl - 1.5*dq_tl(i, j-2)
              lac = pmp - 1.5*dq(i, j-2)
              IF (0. .LT. pmp) THEN
                IF (pmp .LT. lac) THEN
                  x22_tl = lac_tl
                  x22 = lac
                ELSE
                  x22_tl = pmp_tl
                  x22 = pmp
                END IF
              ELSE IF (0. .LT. lac) THEN
                x22_tl = lac_tl
                x22 = lac
              ELSE
                x22 = 0.
                x22_tl = 0.0
              END IF
              IF (0. .GT. pmp) THEN
                IF (pmp .GT. lac) THEN
                  y48_tl = lac_tl
                  y48 = lac
                ELSE
                  y48_tl = pmp_tl
                  y48 = pmp
                END IF
              ELSE IF (0. .GT. lac) THEN
                y48_tl = lac_tl
                y48 = lac
              ELSE
                y48 = 0.
                y48_tl = 0.0
              END IF
              IF (al(i, j+1) - q(i, j) .LT. y48) THEN
                y29_tl = y48_tl
                y29 = y48
              ELSE
                y29_tl = al_tl(i, j+1) - q_tl(i, j)
                y29 = al(i, j+1) - q(i, j)
              END IF
              IF (x22 .GT. y29) THEN
                br_tl(i, j) = y29_tl
                br(i, j) = y29
              ELSE
                br_tl(i, j) = x22_tl
                br(i, j) = x22
              END IF
            END DO
          END DO
        ELSE
          IF (3 .LT. js - 1) THEN
            max7 = js - 1
          ELSE
            max7 = 3
          END IF
          IF (npy - 3 .GT. je + 1) THEN
            min34 = je + 1
            bl_tl = 0.0
            br_tl = 0.0
          ELSE
            min34 = npy - 3
            bl_tl = 0.0
            br_tl = 0.0
          END IF
          DO j=max7,min34
            DO i=ifirst,ilast
              bl_tl(i, j) = al_tl(i, j) - q_tl(i, j)
              bl(i, j) = al(i, j) - q(i, j)
              br_tl(i, j) = al_tl(i, j+1) - q_tl(i, j)
              br(i, j) = al(i, j+1) - q(i, j)
              IF (dq(i, j-1)*dq(i, j) .LE. 0.) THEN
                pmp_tl = -(2.*dq_tl(i, j))
                pmp = -(2.*dq(i, j))
                lac_tl = pmp_tl + 1.5*dq_tl(i, j+1)
                lac = pmp + 1.5*dq(i, j+1)
                IF (0. .LT. pmp) THEN
                  IF (pmp .LT. lac) THEN
                    x23_tl = lac_tl
                    x23 = lac
                  ELSE
                    x23_tl = pmp_tl
                    x23 = pmp
                  END IF
                ELSE IF (0. .LT. lac) THEN
                  x23_tl = lac_tl
                  x23 = lac
                ELSE
                  x23 = 0.
                  x23_tl = 0.0
                END IF
                IF (0. .GT. pmp) THEN
                  IF (pmp .GT. lac) THEN
                    y49_tl = lac_tl
                    y49 = lac
                  ELSE
                    y49_tl = pmp_tl
                    y49 = pmp
                  END IF
                ELSE IF (0. .GT. lac) THEN
                  y49_tl = lac_tl
                  y49 = lac
                ELSE
                  y49 = 0.
                  y49_tl = 0.0
                END IF
                IF (bl(i, j) .LT. y49) THEN
                  y30_tl = y49_tl
                  y30 = y49
                ELSE
                  y30_tl = bl_tl(i, j)
                  y30 = bl(i, j)
                END IF
                IF (x23 .GT. y30) THEN
                  bl_tl(i, j) = y30_tl
                  bl(i, j) = y30
                ELSE
                  bl_tl(i, j) = x23_tl
                  bl(i, j) = x23
                END IF
                pmp_tl = 2.*dq_tl(i, j-1)
                pmp = 2.*dq(i, j-1)
                lac_tl = pmp_tl - 1.5*dq_tl(i, j-2)
                lac = pmp - 1.5*dq(i, j-2)
                IF (0. .LT. pmp) THEN
                  IF (pmp .LT. lac) THEN
                    x24_tl = lac_tl
                    x24 = lac
                  ELSE
                    x24_tl = pmp_tl
                    x24 = pmp
                  END IF
                ELSE IF (0. .LT. lac) THEN
                  x24_tl = lac_tl
                  x24 = lac
                ELSE
                  x24 = 0.
                  x24_tl = 0.0
                END IF
                IF (0. .GT. pmp) THEN
                  IF (pmp .GT. lac) THEN
                    y50_tl = lac_tl
                    y50 = lac
                  ELSE
                    y50_tl = pmp_tl
                    y50 = pmp
                  END IF
                ELSE IF (0. .GT. lac) THEN
                  y50_tl = lac_tl
                  y50 = lac
                ELSE
                  y50 = 0.
                  y50_tl = 0.0
                END IF
                IF (br(i, j) .LT. y50) THEN
                  y31_tl = y50_tl
                  y31 = y50
                ELSE
                  y31_tl = br_tl(i, j)
                  y31 = br(i, j)
                END IF
                IF (x24 .GT. y31) THEN
                  br_tl(i, j) = y31_tl
                  br(i, j) = y31
                ELSE
                  br_tl(i, j) = x24_tl
                  br(i, j) = x24
                END IF
              END IF
            END DO
          END DO
        END IF
!--------------
! Fix the edges:
!--------------
        IF (js .EQ. 1) THEN
          DO i=ifirst,ilast
            br_tl(i, 2) = al_tl(i, 3) - q_tl(i, 2)
            br(i, 2) = al(i, 3) - q(i, 2)
!           xt = t11*(q(i,0)+q(i,1)) + t12*(q(i,-1)+q(i,2))   &
!                                  + t13*(dm(i,2)-dm(i,-1))
            xt_tl = 0.5*((2.*dya(i, 1)+dya(i, 2))*(q_tl(i, 0)+q_tl(i, 1)&
&             )-dya(i, 1)*(q_tl(i, -1)+q_tl(i, 2)))/(dya(i, 1)+dya(i, 2)&
&             )
            xt = 0.5*((2.*dya(i, 1)+dya(i, 2))*(q(i, 0)+q(i, 1))-dya(i, &
&             1)*(q(i, -1)+q(i, 2)))/(dya(i, 1)+dya(i, 2))
            bl_tl(i, 1) = xt_tl - q_tl(i, 1)
            bl(i, 1) = xt - q(i, 1)
            br_tl(i, 0) = xt_tl - q_tl(i, 0)
            br(i, 0) = xt - q(i, 0)
            xt_tl = s14*dm_tl(i, -1) - s11*dq_tl(i, -1) + q_tl(i, 0)
            xt = s14*dm(i, -1) - s11*dq(i, -1) + q(i, 0)
!           xt = min( xt, max(q(i,-1), q(i,0)) )
!           xt = max( xt, min(q(i,-1), q(i,0)) )
            bl_tl(i, 0) = xt_tl - q_tl(i, 0)
            bl(i, 0) = xt - q(i, 0)
            xt_tl = s15*q_tl(i, 1) + s11*q_tl(i, 2) - s14*dm_tl(i, 2)
            xt = s15*q(i, 1) + s11*q(i, 2) - s14*dm(i, 2)
!           xt = min( xt, max(q(i,1), q(i,2)) )
!           xt = max( xt, min(q(i,1), q(i,2)) )
            br_tl(i, 1) = xt_tl - q_tl(i, 1)
            br(i, 1) = xt - q(i, 1)
            bl_tl(i, 2) = xt_tl - q_tl(i, 2)
            bl(i, 2) = xt - q(i, 2)
          END DO
!        if ( jord<=9 ) then
          DO j=0,2
            CALL PERT_PPM_TLM(ilast - ifirst + 1, q(ifirst:, j), bl(&
&                       ifirst:, j), bl_tl(ifirst:, j), br(ifirst:, j), &
&                       br_tl(ifirst:, j), 1)
          END DO
        END IF
!        endif
        IF (je + 1 .EQ. npy) THEN
          DO i=ifirst,ilast
            bl_tl(i, npy-2) = al_tl(i, npy-2) - q_tl(i, npy-2)
            bl(i, npy-2) = al(i, npy-2) - q(i, npy-2)
!           xt = t11*( q(i,npy-1)+q(i,npy)) + t12*(q(i,npy-2)+q(i,npy+1))   &
!                                         + t13*(dm(i,npy+1)-dm(i,npy-2))
            xt_tl = 0.5*((2.*dya(i, npy-1)+dya(i, npy-2))*(q_tl(i, npy-1&
&             )+q_tl(i, npy))-dya(i, npy-1)*(q_tl(i, npy-2)+q_tl(i, npy+&
&             1)))/(dya(i, npy-1)+dya(i, npy-2))
            xt = 0.5*((2.*dya(i, npy-1)+dya(i, npy-2))*(q(i, npy-1)+q(i&
&             , npy))-dya(i, npy-1)*(q(i, npy-2)+q(i, npy+1)))/(dya(i, &
&             npy-1)+dya(i, npy-2))
            br_tl(i, npy-1) = xt_tl - q_tl(i, npy-1)
            br(i, npy-1) = xt - q(i, npy-1)
            bl_tl(i, npy) = xt_tl - q_tl(i, npy)
            bl(i, npy) = xt - q(i, npy)
            xt_tl = s11*dq_tl(i, npy) - s14*dm_tl(i, npy+1) + q_tl(i, &
&             npy)
            xt = s11*dq(i, npy) - s14*dm(i, npy+1) + q(i, npy)
!           xt = min( xt, max( q(i,npy), q(i,npy+1)) )
!           xt = max( xt, min( q(i,npy), q(i,npy+1)) )
            br_tl(i, npy) = xt_tl - q_tl(i, npy)
            br(i, npy) = xt - q(i, npy)
            xt_tl = s15*q_tl(i, npy-1) + s11*q_tl(i, npy-2) + s14*dm_tl(&
&             i, npy-2)
            xt = s15*q(i, npy-1) + s11*q(i, npy-2) + s14*dm(i, npy-2)
!           xt = min( xt, max( q(i,npy-2), q(i,npy-1)) )
!           xt = max( xt, min( q(i,npy-2), q(i,npy-1)) )
            br_tl(i, npy-2) = xt_tl - q_tl(i, npy-2)
            br(i, npy-2) = xt - q(i, npy-2)
            bl_tl(i, npy-1) = xt_tl - q_tl(i, npy-1)
            bl(i, npy-1) = xt - q(i, npy-1)
          END DO
!        if ( jord<=9 ) then
          DO j=npy-2,npy
            CALL PERT_PPM_TLM(ilast - ifirst + 1, q(ifirst:, j), bl(&
&                       ifirst:, j), bl_tl(ifirst:, j), br(ifirst:, j), &
&                       br_tl(ifirst:, j), 1)
          END DO
        END IF
      ELSE
        al_tl = 0.0
!        endif
!---------------
! grid_type == 4
!---------------
        DO j=jfirst-1,jlast+2
          DO i=ifirst,ilast
            al_tl(i, j) = 0.5*(q_tl(i, j-1)+q_tl(i, j)) + r3*(dm_tl(i, j&
&             -1)-dm_tl(i, j))
            al(i, j) = 0.5*(q(i, j-1)+q(i, j)) + r3*(dm(i, j-1)-dm(i, j)&
&             )
          END DO
        END DO
        dq_tl = 0.0
        DO j=jfirst-3,jlast+2
          DO i=ifirst,ilast
            dq_tl(i, j) = q_tl(i, j+1) - q_tl(i, j)
            dq(i, j) = q(i, j+1) - q(i, j)
          END DO
        END DO
        bl_tl = 0.0
        br_tl = 0.0
        DO j=jfirst-1,jlast+1
          DO i=ifirst,ilast
            pmp_tl = -(2.*dq_tl(i, j))
            pmp = -(2.*dq(i, j))
            lac_tl = pmp_tl + 1.5*dq_tl(i, j+1)
            lac = pmp + 1.5*dq(i, j+1)
            IF (0. .LT. pmp) THEN
              IF (pmp .LT. lac) THEN
                x25_tl = lac_tl
                x25 = lac
              ELSE
                x25_tl = pmp_tl
                x25 = pmp
              END IF
            ELSE IF (0. .LT. lac) THEN
              x25_tl = lac_tl
              x25 = lac
            ELSE
              x25 = 0.
              x25_tl = 0.0
            END IF
            IF (0. .GT. pmp) THEN
              IF (pmp .GT. lac) THEN
                y51_tl = lac_tl
                y51 = lac
              ELSE
                y51_tl = pmp_tl
                y51 = pmp
              END IF
            ELSE IF (0. .GT. lac) THEN
              y51_tl = lac_tl
              y51 = lac
            ELSE
              y51 = 0.
              y51_tl = 0.0
            END IF
            IF (al(i, j) - q(i, j) .LT. y51) THEN
              y32_tl = y51_tl
              y32 = y51
            ELSE
              y32_tl = al_tl(i, j) - q_tl(i, j)
              y32 = al(i, j) - q(i, j)
            END IF
            IF (x25 .GT. y32) THEN
              bl_tl(i, j) = y32_tl
              bl(i, j) = y32
            ELSE
              bl_tl(i, j) = x25_tl
              bl(i, j) = x25
            END IF
            pmp_tl = 2.*dq_tl(i, j-1)
            pmp = 2.*dq(i, j-1)
            lac_tl = pmp_tl - 1.5*dq_tl(i, j-2)
            lac = pmp - 1.5*dq(i, j-2)
            IF (0. .LT. pmp) THEN
              IF (pmp .LT. lac) THEN
                x26_tl = lac_tl
                x26 = lac
              ELSE
                x26_tl = pmp_tl
                x26 = pmp
              END IF
            ELSE IF (0. .LT. lac) THEN
              x26_tl = lac_tl
              x26 = lac
            ELSE
              x26 = 0.
              x26_tl = 0.0
            END IF
            IF (0. .GT. pmp) THEN
              IF (pmp .GT. lac) THEN
                y52_tl = lac_tl
                y52 = lac
              ELSE
                y52_tl = pmp_tl
                y52 = pmp
              END IF
            ELSE IF (0. .GT. lac) THEN
              y52_tl = lac_tl
              y52 = lac
            ELSE
              y52 = 0.
              y52_tl = 0.0
            END IF
            IF (al(i, j+1) - q(i, j) .LT. y52) THEN
              y33_tl = y52_tl
              y33 = y52
            ELSE
              y33_tl = al_tl(i, j+1) - q_tl(i, j)
              y33 = al(i, j+1) - q(i, j)
            END IF
            IF (x26 .GT. y33) THEN
              br_tl(i, j) = y33_tl
              br(i, j) = y33
            ELSE
              br_tl(i, j) = x26_tl
              br(i, j) = x26
            END IF
          END DO
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
    ELSE
      dm_tl = 0.0
!-------------------------------
! For positive definite tracers:
!-------------------------------
! jord=11: PPM mono constraint (Lin 2004)
! jord=12: Huynh 2nd constraint (Lin 2004) + positive definite (Lin & Rood 1996)
! jord>12: positive definite only (Lin & Rood 1996)
      DO j=js-2,je+2
        DO i=ifirst,ilast
          xt_tl = 0.25*(q_tl(i, j+1)-q_tl(i, j-1))
          xt = 0.25*(q(i, j+1)-q(i, j-1))
          IF (xt .GE. 0.) THEN
            x27_tl = xt_tl
            x27 = xt
          ELSE
            x27_tl = -xt_tl
            x27 = -xt
          END IF
          IF (q(i, j-1) .LT. q(i, j)) THEN
            IF (q(i, j) .LT. q(i, j+1)) THEN
              max15_tl = q_tl(i, j+1)
              max15 = q(i, j+1)
            ELSE
              max15_tl = q_tl(i, j)
              max15 = q(i, j)
            END IF
          ELSE IF (q(i, j-1) .LT. q(i, j+1)) THEN
            max15_tl = q_tl(i, j+1)
            max15 = q(i, j+1)
          ELSE
            max15_tl = q_tl(i, j-1)
            max15 = q(i, j-1)
          END IF
          y34_tl = max15_tl - q_tl(i, j)
          y34 = max15 - q(i, j)
          IF (q(i, j-1) .GT. q(i, j)) THEN
            IF (q(i, j) .GT. q(i, j+1)) THEN
              min35_tl = q_tl(i, j+1)
              min35 = q(i, j+1)
            ELSE
              min35_tl = q_tl(i, j)
              min35 = q(i, j)
            END IF
          ELSE IF (q(i, j-1) .GT. q(i, j+1)) THEN
            min35_tl = q_tl(i, j+1)
            min35 = q(i, j+1)
          ELSE
            min35_tl = q_tl(i, j-1)
            min35 = q(i, j-1)
          END IF
          z8_tl = q_tl(i, j) - min35_tl
          z8 = q(i, j) - min35
          IF (x27 .GT. y34) THEN
            IF (y34 .GT. z8) THEN
              min16_tl = z8_tl
              min16 = z8
            ELSE
              min16_tl = y34_tl
              min16 = y34
            END IF
          ELSE IF (x27 .GT. z8) THEN
            min16_tl = z8_tl
            min16 = z8
          ELSE
            min16_tl = x27_tl
            min16 = x27
          END IF
          dm_tl(i, j) = min16_tl*SIGN(1.d0, min16*xt)
          dm(i, j) = SIGN(min16, xt)
        END DO
      END DO
      IF (grid_type .LT. 3) THEN
        al_tl = 0.0
        DO j=js3,je2
          DO i=ifirst,ilast
            al_tl(i, j) = 0.5*(q_tl(i, j-1)+q_tl(i, j)) + r3*(dm_tl(i, j&
&             -1)-dm_tl(i, j))
            al(i, j) = 0.5*(q(i, j-1)+q(i, j)) + r3*(dm(i, j-1)-dm(i, j)&
&             )
          END DO
        END DO
        IF (jord .EQ. 11) THEN
          bl_tl = 0.0
          br_tl = 0.0
          DO j=js3,je3
            DO i=ifirst,ilast
              xt_tl = 2.*dm_tl(i, j)
              xt = 2.*dm(i, j)
              IF (xt .GE. 0.) THEN
                x28_tl = xt_tl
                x28 = xt
              ELSE
                x28_tl = -xt_tl
                x28 = -xt
              END IF
              IF (al(i, j) - q(i, j) .GE. 0.) THEN
                y35_tl = al_tl(i, j) - q_tl(i, j)
                y35 = al(i, j) - q(i, j)
              ELSE
                y35_tl = -(al_tl(i, j)-q_tl(i, j))
                y35 = -(al(i, j)-q(i, j))
              END IF
              IF (x28 .GT. y35) THEN
                min17_tl = y35_tl
                min17 = y35
              ELSE
                min17_tl = x28_tl
                min17 = x28
              END IF
              bl_tl(i, j) = -(min17_tl*SIGN(1.d0, min17*xt))
              bl(i, j) = -SIGN(min17, xt)
              IF (xt .GE. 0.) THEN
                x29_tl = xt_tl
                x29 = xt
              ELSE
                x29_tl = -xt_tl
                x29 = -xt
              END IF
              IF (al(i, j+1) - q(i, j) .GE. 0.) THEN
                y36_tl = al_tl(i, j+1) - q_tl(i, j)
                y36 = al(i, j+1) - q(i, j)
              ELSE
                y36_tl = -(al_tl(i, j+1)-q_tl(i, j))
                y36 = -(al(i, j+1)-q(i, j))
              END IF
              IF (x29 .GT. y36) THEN
                min18_tl = y36_tl
                min18 = y36
              ELSE
                min18_tl = x29_tl
                min18 = x29
              END IF
              br_tl(i, j) = min18_tl*SIGN(1.d0, min18*xt)
              br(i, j) = SIGN(min18, xt)
            END DO
          END DO
        ELSE IF (jord .EQ. 12) THEN
          dq_tl = 0.0
          DO j=js-3,je+2
            DO i=ifirst,ilast
              dq_tl(i, j) = q_tl(i, j+1) - q_tl(i, j)
              dq(i, j) = q(i, j+1) - q(i, j)
            END DO
          END DO
          bl_tl = 0.0
          br_tl = 0.0
          DO j=js3,je3
            DO i=ifirst,ilast
              pmp_tl = -(2.*dq_tl(i, j))
              pmp = -(2.*dq(i, j))
              lac_tl = pmp_tl + 1.5*dq_tl(i, j+1)
              lac = pmp + 1.5*dq(i, j+1)
              IF (0. .LT. pmp) THEN
                IF (pmp .LT. lac) THEN
                  x30_tl = lac_tl
                  x30 = lac
                ELSE
                  x30_tl = pmp_tl
                  x30 = pmp
                END IF
              ELSE IF (0. .LT. lac) THEN
                x30_tl = lac_tl
                x30 = lac
              ELSE
                x30 = 0.
                x30_tl = 0.0
              END IF
              IF (0. .GT. pmp) THEN
                IF (pmp .GT. lac) THEN
                  y53_tl = lac_tl
                  y53 = lac
                ELSE
                  y53_tl = pmp_tl
                  y53 = pmp
                END IF
              ELSE IF (0. .GT. lac) THEN
                y53_tl = lac_tl
                y53 = lac
              ELSE
                y53 = 0.
                y53_tl = 0.0
              END IF
              IF (al(i, j) - q(i, j) .LT. y53) THEN
                y37_tl = y53_tl
                y37 = y53
              ELSE
                y37_tl = al_tl(i, j) - q_tl(i, j)
                y37 = al(i, j) - q(i, j)
              END IF
              IF (x30 .GT. y37) THEN
                bl_tl(i, j) = y37_tl
                bl(i, j) = y37
              ELSE
                bl_tl(i, j) = x30_tl
                bl(i, j) = x30
              END IF
              pmp_tl = 2.*dq_tl(i, j-1)
              pmp = 2.*dq(i, j-1)
              lac_tl = pmp_tl - 1.5*dq_tl(i, j-2)
              lac = pmp - 1.5*dq(i, j-2)
              IF (0. .LT. pmp) THEN
                IF (pmp .LT. lac) THEN
                  x31_tl = lac_tl
                  x31 = lac
                ELSE
                  x31_tl = pmp_tl
                  x31 = pmp
                END IF
              ELSE IF (0. .LT. lac) THEN
                x31_tl = lac_tl
                x31 = lac
              ELSE
                x31 = 0.
                x31_tl = 0.0
              END IF
              IF (0. .GT. pmp) THEN
                IF (pmp .GT. lac) THEN
                  y54_tl = lac_tl
                  y54 = lac
                ELSE
                  y54_tl = pmp_tl
                  y54 = pmp
                END IF
              ELSE IF (0. .GT. lac) THEN
                y54_tl = lac_tl
                y54 = lac
              ELSE
                y54 = 0.
                y54_tl = 0.0
              END IF
              IF (al(i, j+1) - q(i, j) .LT. y54) THEN
                y38_tl = y54_tl
                y38 = y54
              ELSE
                y38_tl = al_tl(i, j+1) - q_tl(i, j)
                y38 = al(i, j+1) - q(i, j)
              END IF
              IF (x31 .GT. y38) THEN
                br_tl(i, j) = y38_tl
                br(i, j) = y38
              ELSE
                br_tl(i, j) = x31_tl
                br(i, j) = x31
              END IF
            END DO
          END DO
        ELSE
          bl_tl = 0.0
          br_tl = 0.0
          DO j=js3,je3
            DO i=ifirst,ilast
              bl_tl(i, j) = al_tl(i, j) - q_tl(i, j)
              bl(i, j) = al(i, j) - q(i, j)
              br_tl(i, j) = al_tl(i, j+1) - q_tl(i, j)
              br(i, j) = al(i, j+1) - q(i, j)
            END DO
          END DO
        END IF
        IF (jord .NE. 11) THEN
! Positive definite constraint:
          DO j=js3,je3
            CALL PERT_PPM_TLM(ilast - ifirst + 1, q(ifirst:, j), bl(&
&                       ifirst:, j), bl_tl(ifirst:, j), br(ifirst:, j), &
&                       br_tl(ifirst:, j), 0)
          END DO
        END IF
!--------------
! Fix the edges:
!--------------
        IF (js .EQ. 1) THEN
          DO i=ifirst,ilast
            br_tl(i, 2) = al_tl(i, 3) - q_tl(i, 2)
            br(i, 2) = al(i, 3) - q(i, 2)
!           xt = t11*(q(i,0)+q(i,1)) + t12*(q(i,-1)+q(i,2))   &
!              + t13*(dm(i,2)-dm(i,-1))
!!!         xt = 0.75*(q(i,0)+q(i,1)) - 0.25*(q(i,-1)+q(i,2))
            xt_tl = 0.5*((2.*dya(i, 1)+dya(i, 2))*(q_tl(i, 0)+q_tl(i, 1)&
&             )-dya(i, 1)*(q_tl(i, -1)+q_tl(i, 2)))/(dya(i, 1)+dya(i, 2)&
&             )
            xt = 0.5*((2.*dya(i, 1)+dya(i, 2))*(q(i, 0)+q(i, 1))-dya(i, &
&             1)*(q(i, -1)+q(i, 2)))/(dya(i, 1)+dya(i, 2))
            IF (0. .LT. xt) THEN
              xt = xt
            ELSE
              xt = 0.
              xt_tl = 0.0
            END IF
            bl_tl(i, 1) = xt_tl - q_tl(i, 1)
            bl(i, 1) = xt - q(i, 1)
            br_tl(i, 0) = xt_tl - q_tl(i, 0)
            br(i, 0) = xt - q(i, 0)
            xt_tl = 4.*dm_tl(i, -1)/7. + 11.*q_tl(i, -1)/14. + 3.*q_tl(i&
&             , 0)/14.
            xt = 4./7.*dm(i, -1) + 11./14.*q(i, -1) + 3./14.*q(i, 0)
            IF (0. .LT. xt) THEN
              xt = xt
            ELSE
              xt = 0.
              xt_tl = 0.0
            END IF
            bl_tl(i, 0) = xt_tl - q_tl(i, 0)
            bl(i, 0) = xt - q(i, 0)
            xt_tl = 3.*q_tl(i, 1)/14. + 11.*q_tl(i, 2)/14. - 4.*dm_tl(i&
&             , 2)/7.
            xt = 3./14.*q(i, 1) + 11./14.*q(i, 2) - 4./7.*dm(i, 2)
            IF (0. .LT. xt) THEN
              xt = xt
            ELSE
              xt = 0.
              xt_tl = 0.0
            END IF
            br_tl(i, 1) = xt_tl - q_tl(i, 1)
            br(i, 1) = xt - q(i, 1)
            bl_tl(i, 2) = xt_tl - q_tl(i, 2)
            bl(i, 2) = xt - q(i, 2)
          END DO
          DO j=0,2
            CALL PERT_PPM_TLM(ilast - ifirst + 1, q(ifirst:, j), bl(&
&                       ifirst:, j), bl_tl(ifirst:, j), br(ifirst:, j), &
&                       br_tl(ifirst:, j), 1)
          END DO
        END IF
        IF (je + 1 .EQ. npy) THEN
          DO i=ifirst,ilast
            bl_tl(i, npy-2) = al_tl(i, npy-2) - q_tl(i, npy-2)
            bl(i, npy-2) = al(i, npy-2) - q(i, npy-2)
!           xt = t11*(q(i,npy-1)+q(i,npy)) + t12*(q(i,npy-2)+q(i,npy+1))   &
!               + t13*(dm(i,npy+1)-dm(i,npy-2))
!!!         xt = 0.75*(q(i,npy-1)+q(i,npy)) - 0.25*(q(i,npy-2)+q(i,npy+1))
            xt_tl = 0.5*((2.*dya(i, npy-1)+dya(i, npy-2))*(q_tl(i, npy-1&
&             )+q_tl(i, npy))-dya(i, npy-1)*(q_tl(i, npy-2)+q_tl(i, npy+&
&             1)))/(dya(i, npy-1)+dya(i, npy-2))
            xt = 0.5*((2.*dya(i, npy-1)+dya(i, npy-2))*(q(i, npy-1)+q(i&
&             , npy))-dya(i, npy-1)*(q(i, npy-2)+q(i, npy+1)))/(dya(i, &
&             npy-1)+dya(i, npy-2))
            IF (0. .LT. xt) THEN
              xt = xt
            ELSE
              xt = 0.
              xt_tl = 0.0
            END IF
            br_tl(i, npy-1) = xt_tl - q_tl(i, npy-1)
            br(i, npy-1) = xt - q(i, npy-1)
            bl_tl(i, npy) = xt_tl - q_tl(i, npy)
            bl(i, npy) = xt - q(i, npy)
            xt_tl = 3.*q_tl(i, npy)/14. + 11.*q_tl(i, npy+1)/14. - 4.*&
&             dm_tl(i, npy+1)/7.
            xt = 3./14.*q(i, npy) + 11./14.*q(i, npy+1) - 4./7.*dm(i, &
&             npy+1)
            IF (0. .LT. xt) THEN
              xt = xt
            ELSE
              xt = 0.
              xt_tl = 0.0
            END IF
            br_tl(i, npy) = xt_tl - q_tl(i, npy)
            br(i, npy) = xt - q(i, npy)
            xt_tl = 3.*q_tl(i, npy-1)/14. + 11.*q_tl(i, npy-2)/14. + 4.*&
&             dm_tl(i, npy-2)/7.
            xt = 3./14.*q(i, npy-1) + 11./14.*q(i, npy-2) + 4./7.*dm(i, &
&             npy-2)
            IF (0. .LT. xt) THEN
              xt = xt
            ELSE
              xt = 0.
              xt_tl = 0.0
            END IF
            br_tl(i, npy-2) = xt_tl - q_tl(i, npy-2)
            br(i, npy-2) = xt - q(i, npy-2)
            bl_tl(i, npy-1) = xt_tl - q_tl(i, npy-1)
            bl(i, npy-1) = xt - q(i, npy-1)
          END DO
          DO j=npy-2,npy
            CALL PERT_PPM_TLM(ilast - ifirst + 1, q(ifirst:, j), bl(&
&                       ifirst:, j), bl_tl(ifirst:, j), br(ifirst:, j), &
&                       br_tl(ifirst:, j), 1)
          END DO
        END IF
      ELSE
        al_tl = 0.0
        DO j=js-1,je+2
          DO i=ifirst,ilast
            al_tl(i, j) = 0.5*(q_tl(i, j-1)+q_tl(i, j)) + r3*(dm_tl(i, j&
&             -1)-dm_tl(i, j))
            al(i, j) = 0.5*(q(i, j-1)+q(i, j)) + r3*(dm(i, j-1)-dm(i, j)&
&             )
          END DO
        END DO
        IF (jord .EQ. 11) THEN
          bl_tl = 0.0
          br_tl = 0.0
          DO j=js-1,je+1
            DO i=ifirst,ilast
              xt_tl = 2.*dm_tl(i, j)
              xt = 2.*dm(i, j)
              IF (xt .GE. 0.) THEN
                x32_tl = xt_tl
                x32 = xt
              ELSE
                x32_tl = -xt_tl
                x32 = -xt
              END IF
              IF (al(i, j) - q(i, j) .GE. 0.) THEN
                y39_tl = al_tl(i, j) - q_tl(i, j)
                y39 = al(i, j) - q(i, j)
              ELSE
                y39_tl = -(al_tl(i, j)-q_tl(i, j))
                y39 = -(al(i, j)-q(i, j))
              END IF
              IF (x32 .GT. y39) THEN
                min19_tl = y39_tl
                min19 = y39
              ELSE
                min19_tl = x32_tl
                min19 = x32
              END IF
              bl_tl(i, j) = -(min19_tl*SIGN(1.d0, min19*xt))
              bl(i, j) = -SIGN(min19, xt)
              IF (xt .GE. 0.) THEN
                x33_tl = xt_tl
                x33 = xt
              ELSE
                x33_tl = -xt_tl
                x33 = -xt
              END IF
              IF (al(i, j+1) - q(i, j) .GE. 0.) THEN
                y40_tl = al_tl(i, j+1) - q_tl(i, j)
                y40 = al(i, j+1) - q(i, j)
              ELSE
                y40_tl = -(al_tl(i, j+1)-q_tl(i, j))
                y40 = -(al(i, j+1)-q(i, j))
              END IF
              IF (x33 .GT. y40) THEN
                min20_tl = y40_tl
                min20 = y40
              ELSE
                min20_tl = x33_tl
                min20 = x33
              END IF
              br_tl(i, j) = min20_tl*SIGN(1.d0, min20*xt)
              br(i, j) = SIGN(min20, xt)
            END DO
          END DO
        ELSE IF (jord .EQ. 12) THEN
          dq_tl = 0.0
          DO j=js-3,je+2
            DO i=ifirst,ilast
              dq_tl(i, j) = q_tl(i, j+1) - q_tl(i, j)
              dq(i, j) = q(i, j+1) - q(i, j)
            END DO
          END DO
          bl_tl = 0.0
          br_tl = 0.0
          DO j=js-1,je+1
            DO i=ifirst,ilast
              pmp_tl = -(2.*dq_tl(i, j))
              pmp = -(2.*dq(i, j))
              lac_tl = pmp_tl + 1.5*dq_tl(i, j+1)
              lac = pmp + 1.5*dq(i, j+1)
              IF (0. .LT. pmp) THEN
                IF (pmp .LT. lac) THEN
                  x34_tl = lac_tl
                  x34 = lac
                ELSE
                  x34_tl = pmp_tl
                  x34 = pmp
                END IF
              ELSE IF (0. .LT. lac) THEN
                x34_tl = lac_tl
                x34 = lac
              ELSE
                x34 = 0.
                x34_tl = 0.0
              END IF
              IF (0. .GT. pmp) THEN
                IF (pmp .GT. lac) THEN
                  y55_tl = lac_tl
                  y55 = lac
                ELSE
                  y55_tl = pmp_tl
                  y55 = pmp
                END IF
              ELSE IF (0. .GT. lac) THEN
                y55_tl = lac_tl
                y55 = lac
              ELSE
                y55 = 0.
                y55_tl = 0.0
              END IF
              IF (al(i, j) - q(i, j) .LT. y55) THEN
                y41_tl = y55_tl
                y41 = y55
              ELSE
                y41_tl = al_tl(i, j) - q_tl(i, j)
                y41 = al(i, j) - q(i, j)
              END IF
              IF (x34 .GT. y41) THEN
                bl_tl(i, j) = y41_tl
                bl(i, j) = y41
              ELSE
                bl_tl(i, j) = x34_tl
                bl(i, j) = x34
              END IF
              pmp_tl = 2.*dq_tl(i, j-1)
              pmp = 2.*dq(i, j-1)
              lac_tl = pmp_tl - 1.5*dq_tl(i, j-2)
              lac = pmp - 1.5*dq(i, j-2)
              IF (0. .LT. pmp) THEN
                IF (pmp .LT. lac) THEN
                  x35_tl = lac_tl
                  x35 = lac
                ELSE
                  x35_tl = pmp_tl
                  x35 = pmp
                END IF
              ELSE IF (0. .LT. lac) THEN
                x35_tl = lac_tl
                x35 = lac
              ELSE
                x35 = 0.
                x35_tl = 0.0
              END IF
              IF (0. .GT. pmp) THEN
                IF (pmp .GT. lac) THEN
                  y56_tl = lac_tl
                  y56 = lac
                ELSE
                  y56_tl = pmp_tl
                  y56 = pmp
                END IF
              ELSE IF (0. .GT. lac) THEN
                y56_tl = lac_tl
                y56 = lac
              ELSE
                y56 = 0.
                y56_tl = 0.0
              END IF
              IF (al(i, j+1) - q(i, j) .LT. y56) THEN
                y42_tl = y56_tl
                y42 = y56
              ELSE
                y42_tl = al_tl(i, j+1) - q_tl(i, j)
                y42 = al(i, j+1) - q(i, j)
              END IF
              IF (x35 .GT. y42) THEN
                br_tl(i, j) = y42_tl
                br(i, j) = y42
              ELSE
                br_tl(i, j) = x35_tl
                br(i, j) = x35
              END IF
            END DO
          END DO
        ELSE
          bl_tl = 0.0
          br_tl = 0.0
          DO j=js-1,je+1
            DO i=ifirst,ilast
              bl_tl(i, j) = al_tl(i, j) - q_tl(i, j)
              bl(i, j) = al(i, j) - q(i, j)
              br_tl(i, j) = al_tl(i, j+1) - q_tl(i, j)
              br(i, j) = al(i, j+1) - q(i, j)
            END DO
          END DO
        END IF
        IF (jord .NE. 11) THEN
! Positive definite constraint:
          DO j=js-1,je+1
            CALL PERT_PPM_TLM(ilast - ifirst + 1, q(ifirst:, j), bl(&
&                       ifirst:, j), bl_tl(ifirst:, j), br(ifirst:, j), &
&                       br_tl(ifirst:, j), 0)
          END DO
        END IF
      END IF
      DO j=js,je+1
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
    REAL :: dm1_tl(ifirst-2:ilast+2)
    REAL :: al(ifirst-1:ilast+2)
    REAL :: al_tl(ifirst-1:ilast+2)
    REAL :: bl(ifirst-1:ilast+1)
    REAL :: bl_tl(ifirst-1:ilast+1)
    REAL :: br(ifirst-1:ilast+1)
    REAL :: br_tl(ifirst-1:ilast+1)
    REAL :: dq(ifirst-3:ilast+2)
    REAL :: dq_tl(ifirst-3:ilast+2)
    REAL :: dl, dr, pmp, lac, ct, qe
    REAL :: dl_tl, dr_tl, pmp_tl, lac_tl
    REAL :: xt, x1, x0
    REAL :: xt_tl, x1_tl, x0_tl
    INTEGER :: i, j, is3, ie3, ie2, it
    INTRINSIC MAX
    INTRINSIC MIN
    INTRINSIC ABS
    INTRINSIC SIGN
    REAL :: y1_tl
    REAL :: min6_tl
    REAL :: y18_tl
    REAL :: y44_tl
    REAL :: min22_tl
    REAL :: x19
    REAL :: x25_tl
    REAL :: x18
    REAL :: y32_tl
    REAL :: x17
    REAL :: min10_tl
    REAL :: x16
    REAL :: z4_tl
    REAL :: x13_tl
    REAL :: x15
    REAL :: max7_tl
    REAL :: y20_tl
    REAL :: y57_tl
    REAL :: x14
    REAL :: min9
    REAL :: x13
    REAL :: min8
    REAL :: y2_tl
    REAL :: min7_tl
    REAL :: y19_tl
    REAL :: x12
    REAL :: min7
    REAL :: y45_tl
    REAL :: x11
    REAL :: y29
    REAL :: min6
    REAL :: min23_tl
    REAL :: x10
    REAL :: y28
    REAL :: min5
    REAL :: x26_tl
    REAL :: y27
    REAL :: min4
    REAL :: y33_tl
    REAL :: y26
    REAL :: min3
    REAL :: min32
    REAL :: min11_tl
    REAL :: y25
    REAL :: min2
    REAL :: min31
    REAL :: z5_tl
    REAL :: x14_tl
    REAL :: y24
    REAL :: min1
    REAL :: min30
    REAL :: y21_tl
    REAL :: max8_tl
    REAL :: y23
    REAL :: y22
    REAL :: y3_tl
    REAL :: min8_tl
    REAL :: y21
    REAL :: y46_tl
    REAL :: y20
    REAL :: y57
    REAL :: y56
    REAL :: x27_tl
    REAL :: x9
    REAL :: y55
    REAL :: y34_tl
    REAL :: x8
    REAL :: y54
    REAL :: min12_tl
    REAL :: x7
    REAL :: y53
    REAL :: z6_tl
    REAL :: x15_tl
    REAL :: x6
    REAL :: y52
    REAL :: y22_tl
    REAL :: max9_tl
    REAL :: x5
    REAL :: y51
    REAL :: x4
    REAL :: y50
    REAL :: y4_tl
    REAL :: min9_tl
    REAL :: x3
    REAL :: y10_tl
    REAL :: y47_tl
    REAL :: x2
    REAL :: min25_tl
    REAL :: x2_tl
    REAL :: x28_tl
    REAL :: y35_tl
    REAL :: min13_tl
    REAL :: x16_tl
    REAL :: z7_tl
    REAL :: y23_tl
    REAL :: y5_tl
    REAL :: y11_tl
    REAL :: y48_tl
    REAL :: x30_tl
    REAL :: min29
    REAL :: min26_tl
    REAL :: min28
    REAL :: x3_tl
    REAL :: x29_tl
    REAL :: min27
    REAL :: y36_tl
    REAL :: min26
    REAL :: y19
    REAL :: min25
    REAL :: x17_tl
    REAL :: z8_tl
    REAL :: y18
    INTEGER :: min24
    REAL :: y24_tl
    REAL :: y17
    REAL :: x36
    REAL :: min23
    REAL :: y50_tl
    REAL :: y16
    REAL :: x35
    REAL :: min22
    REAL :: y6_tl
    REAL :: y15
    REAL :: x34
    REAL :: min21
    REAL :: y12_tl
    REAL :: y49_tl
    REAL :: x31_tl
    REAL :: y14
    REAL :: x33
    REAL :: min20
    REAL :: min27_tl
    REAL :: y13
    REAL :: x32
    REAL :: x4_tl
    REAL :: y12
    REAL :: x31
    REAL :: y49
    REAL :: y37_tl
    REAL :: y11
    REAL :: x30
    REAL :: y48
    REAL :: min15_tl
    REAL :: y10
    REAL :: y47
    REAL :: x18_tl
    REAL :: max10_tl
    REAL :: z9_tl
    REAL :: y46
    REAL :: y25_tl
    REAL :: y45
    REAL :: y51_tl
    REAL :: y44
    REAL :: y7_tl
    REAL :: y43
    REAL :: min1_tl
    REAL :: y13_tl
    REAL :: x32_tl
    REAL :: y42
    REAL :: min28_tl
    REAL :: y41
    REAL :: x5_tl
    REAL :: y40
    REAL :: x20_tl
    REAL :: y38_tl
    REAL :: min16_tl
    REAL :: z9
    REAL :: x19_tl
    REAL :: z8
    REAL :: y26_tl
    REAL :: z7
    REAL :: max2_tl
    REAL :: y52_tl
    REAL :: z6
    REAL :: y8_tl
    REAL :: min30_tl
    REAL :: z5
    REAL :: min2_tl
    REAL :: y14_tl
    REAL :: x33_tl
    REAL :: z4
    REAL :: min29_tl
    REAL :: y40_tl
    REAL :: z3
    REAL :: x6_tl
    REAL :: z2
    REAL :: x21_tl
    REAL :: y39_tl
    REAL :: z1
    REAL :: min17_tl
    REAL :: min19
    REAL :: min18
    REAL :: y27_tl
    REAL :: min17
    REAL :: max3_tl
    REAL :: y53_tl
    REAL :: x29
    REAL :: min16
    REAL :: y9_tl
    REAL :: min31_tl
    REAL :: x28
    REAL :: min15
    REAL :: min3_tl
    REAL :: y15_tl
    REAL :: x34_tl
    REAL :: x27
    INTEGER :: min14
    REAL :: y41_tl
    REAL :: x26
    REAL :: min13
    REAL :: x7_tl
    REAL :: x25
    REAL :: min12
    REAL :: x22_tl
    REAL :: x24
    REAL :: min11
    REAL :: min18_tl
    REAL :: x23
    REAL :: min10
    REAL :: x22
    REAL :: z1_tl
    REAL :: x10_tl
    REAL :: y28_tl
    REAL :: x21
    REAL :: y39
    REAL :: max4_tl
    REAL :: y54_tl
    REAL :: x20
    REAL :: y38
    REAL :: min32_tl
    REAL :: y37
    REAL :: min4_tl
    REAL :: y16_tl
    REAL :: x35_tl
    REAL :: y36
    REAL :: y42_tl
    REAL :: y35
    REAL :: x8_tl
    REAL :: min20_tl
    REAL :: y34
    REAL :: x23_tl
    REAL :: y33
    REAL :: max9
    REAL :: y30_tl
    REAL :: min19_tl
    REAL :: y32
    REAL :: max8
    REAL :: y31
    REAL :: max7
    REAL :: z2_tl
    REAL :: x11_tl
    REAL :: y29_tl
    REAL :: y30
    REAL :: max6
    REAL :: max5_tl
    REAL :: y55_tl
    REAL :: max5
    REAL :: y9
    REAL :: max4
    REAL :: min5_tl
    REAL :: y17_tl
    REAL :: x36_tl
    REAL :: y8
    REAL :: max3
    REAL :: y43_tl
    REAL :: y7
    REAL :: max2
    REAL :: x9_tl
    REAL :: min21_tl
    REAL :: y6
    INTEGER :: max1
    REAL :: x24_tl
    REAL :: y5
    REAL :: y31_tl
    REAL :: y4
    REAL :: max10
    REAL :: y3
    REAL :: z3_tl
    REAL :: x12_tl
    REAL :: y2
    REAL :: max6_tl
    REAL :: y56_tl
    REAL :: y1
    REAL :: cfl, cfl_tl
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
    IF (iord .LE. 4) THEN
      al_tl = 0.0
      dm1_tl = 0.0
      bl_tl = 0.0
      br_tl = 0.0
      DO j=jfirst,jlast
        DO i=is-2,ie+2
          xt_tl = 0.25*(q_tl(i+1, j)-q_tl(i-1, j))
          xt = 0.25*(q(i+1, j)-q(i-1, j))
          IF (xt .GE. 0.) THEN
            x2_tl = xt_tl
            x2 = xt
          ELSE
            x2_tl = -xt_tl
            x2 = -xt
          END IF
          IF (q(i-1, j) .LT. q(i, j)) THEN
            IF (q(i, j) .LT. q(i+1, j)) THEN
              max2_tl = q_tl(i+1, j)
              max2 = q(i+1, j)
            ELSE
              max2_tl = q_tl(i, j)
              max2 = q(i, j)
            END IF
          ELSE IF (q(i-1, j) .LT. q(i+1, j)) THEN
            max2_tl = q_tl(i+1, j)
            max2 = q(i+1, j)
          ELSE
            max2_tl = q_tl(i-1, j)
            max2 = q(i-1, j)
          END IF
          y1_tl = max2_tl - q_tl(i, j)
          y1 = max2 - q(i, j)
          IF (q(i-1, j) .GT. q(i, j)) THEN
            IF (q(i, j) .GT. q(i+1, j)) THEN
              min23_tl = q_tl(i+1, j)
              min23 = q(i+1, j)
            ELSE
              min23_tl = q_tl(i, j)
              min23 = q(i, j)
            END IF
          ELSE IF (q(i-1, j) .GT. q(i+1, j)) THEN
            min23_tl = q_tl(i+1, j)
            min23 = q(i+1, j)
          ELSE
            min23_tl = q_tl(i-1, j)
            min23 = q(i-1, j)
          END IF
          z1_tl = q_tl(i, j) - min23_tl
          z1 = q(i, j) - min23
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
          dm1_tl(i) = min1_tl*SIGN(1.d0, min1*xt)
          dm1(i) = SIGN(min1, xt)
        END DO
        IF (grid_type .LT. 3) THEN
          IF (3 .LT. is - 1) THEN
            max1 = is - 1
          ELSE
            max1 = 3
          END IF
          IF (npx - 2 .GT. ie + 2) THEN
            min24 = ie + 2
          ELSE
            min24 = npx - 2
          END IF
          DO i=max1,min24
            al_tl(i) = 0.5*(q_tl(i-1, j)+q_tl(i, j)) + r3*(dm1_tl(i-1)-&
&             dm1_tl(i))
            al(i) = 0.5*(q(i-1, j)+q(i, j)) + r3*(dm1(i-1)-dm1(i))
          END DO
! Fix the edges:
          IF (is .EQ. 1) THEN
            x0_tl = 0.5*((2.*dxa(1, j)+dxa(2, j))*(q_tl(0, j)+q_tl(1, j)&
&             )-dxa(1, j)*(q_tl(-1, j)+q_tl(2, j)))/(dxa(1, j)+dxa(2, j)&
&             )
            x0 = 0.5*((2.*dxa(1, j)+dxa(2, j))*(q(0, j)+q(1, j))-dxa(1, &
&             j)*(q(-1, j)+q(2, j)))/(dxa(1, j)+dxa(2, j))
            al_tl(1) = x0_tl
            al(1) = x0
            x1_tl = s15*q_tl(0, j) + s11*q_tl(-1, j) + s14*dm1_tl(-1)
            x1 = s15*q(0, j) + s11*q(-1, j) + s14*dm1(-1)
            dm1_tl(0) = 0.5*(x0_tl-x1_tl)
            dm1(0) = 0.5*(x0-x1)
            IF (dm1(0) .GE. 0.) THEN
              x3_tl = dm1_tl(0)
              x3 = dm1(0)
            ELSE
              x3_tl = -dm1_tl(0)
              x3 = -dm1(0)
            END IF
            IF (q(0, j) .LT. x0) THEN
              IF (x0 .LT. x1) THEN
                max3_tl = x1_tl
                max3 = x1
              ELSE
                max3_tl = x0_tl
                max3 = x0
              END IF
            ELSE IF (q(0, j) .LT. x1) THEN
              max3_tl = x1_tl
              max3 = x1
            ELSE
              max3_tl = q_tl(0, j)
              max3 = q(0, j)
            END IF
            y2_tl = max3_tl - q_tl(0, j)
            y2 = max3 - q(0, j)
            IF (q(0, j) .GT. x0) THEN
              IF (x0 .GT. x1) THEN
                min25_tl = x1_tl
                min25 = x1
              ELSE
                min25_tl = x0_tl
                min25 = x0
              END IF
            ELSE IF (q(0, j) .GT. x1) THEN
              min25_tl = x1_tl
              min25 = x1
            ELSE
              min25_tl = q_tl(0, j)
              min25 = q(0, j)
            END IF
            z2_tl = q_tl(0, j) - min25_tl
            z2 = q(0, j) - min25
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
            dm1_tl(0) = min2_tl*SIGN(1.d0, min2*dm1(0))
            dm1(0) = SIGN(min2, dm1(0))
            al_tl(0) = 0.5*(q_tl(-1, j)+q_tl(0, j)) + r3*(dm1_tl(-1)-&
&             dm1_tl(0))
            al(0) = 0.5*(q(-1, j)+q(0, j)) + r3*(dm1(-1)-dm1(0))
!
            x1_tl = s15*q_tl(1, j) + s11*q_tl(2, j) - s14*dm1_tl(2)
            x1 = s15*q(1, j) + s11*q(2, j) - s14*dm1(2)
            dm1_tl(1) = 0.5*(x1_tl-x0_tl)
            dm1(1) = 0.5*(x1-x0)
            IF (dm1(1) .GE. 0.) THEN
              x4_tl = dm1_tl(1)
              x4 = dm1(1)
            ELSE
              x4_tl = -dm1_tl(1)
              x4 = -dm1(1)
            END IF
            IF (q(1, j) .LT. x0) THEN
              IF (x0 .LT. x1) THEN
                max4_tl = x1_tl
                max4 = x1
              ELSE
                max4_tl = x0_tl
                max4 = x0
              END IF
            ELSE IF (q(1, j) .LT. x1) THEN
              max4_tl = x1_tl
              max4 = x1
            ELSE
              max4_tl = q_tl(1, j)
              max4 = q(1, j)
            END IF
            y3_tl = max4_tl - q_tl(1, j)
            y3 = max4 - q(1, j)
            IF (q(1, j) .GT. x0) THEN
              IF (x0 .GT. x1) THEN
                min26_tl = x1_tl
                min26 = x1
              ELSE
                min26_tl = x0_tl
                min26 = x0
              END IF
            ELSE IF (q(1, j) .GT. x1) THEN
              min26_tl = x1_tl
              min26 = x1
            ELSE
              min26_tl = q_tl(1, j)
              min26 = q(1, j)
            END IF
            z3_tl = q_tl(1, j) - min26_tl
            z3 = q(1, j) - min26
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
            dm1_tl(1) = min3_tl*SIGN(1.d0, min3*dm1(1))
            dm1(1) = SIGN(min3, dm1(1))
            al_tl(2) = 0.5*(q_tl(1, j)+q_tl(2, j)) + r3*(dm1_tl(1)-&
&             dm1_tl(2))
            al(2) = 0.5*(q(1, j)+q(2, j)) + r3*(dm1(1)-dm1(2))
          END IF
          IF (ie + 1 .EQ. npx) THEN
            x0_tl = 0.5*((2.*dxa(npx-1, j)+dxa(npx-2, j))*(q_tl(npx-1, j&
&             )+q_tl(npx, j))-dxa(npx-1, j)*(q_tl(npx-2, j)+q_tl(npx+1, &
&             j)))/(dxa(npx-1, j)+dxa(npx-2, j))
            x0 = 0.5*((2.*dxa(npx-1, j)+dxa(npx-2, j))*(q(npx-1, j)+q(&
&             npx, j))-dxa(npx-1, j)*(q(npx-2, j)+q(npx+1, j)))/(dxa(npx&
&             -1, j)+dxa(npx-2, j))
            al_tl(npx) = x0_tl
            al(npx) = x0
            x1_tl = s15*q_tl(npx-1, j) + s11*q_tl(npx-2, j) + s14*dm1_tl&
&             (npx-2)
            x1 = s15*q(npx-1, j) + s11*q(npx-2, j) + s14*dm1(npx-2)
            dm1_tl(npx-1) = 0.5*(x0_tl-x1_tl)
            dm1(npx-1) = 0.5*(x0-x1)
            IF (dm1(npx-1) .GE. 0.) THEN
              x5_tl = dm1_tl(npx-1)
              x5 = dm1(npx-1)
            ELSE
              x5_tl = -dm1_tl(npx-1)
              x5 = -dm1(npx-1)
            END IF
            IF (q(npx-1, j) .LT. x0) THEN
              IF (x0 .LT. x1) THEN
                max5_tl = x1_tl
                max5 = x1
              ELSE
                max5_tl = x0_tl
                max5 = x0
              END IF
            ELSE IF (q(npx-1, j) .LT. x1) THEN
              max5_tl = x1_tl
              max5 = x1
            ELSE
              max5_tl = q_tl(npx-1, j)
              max5 = q(npx-1, j)
            END IF
            y4_tl = max5_tl - q_tl(npx-1, j)
            y4 = max5 - q(npx-1, j)
            IF (q(npx-1, j) .GT. x0) THEN
              IF (x0 .GT. x1) THEN
                min27_tl = x1_tl
                min27 = x1
              ELSE
                min27_tl = x0_tl
                min27 = x0
              END IF
            ELSE IF (q(npx-1, j) .GT. x1) THEN
              min27_tl = x1_tl
              min27 = x1
            ELSE
              min27_tl = q_tl(npx-1, j)
              min27 = q(npx-1, j)
            END IF
            z4_tl = q_tl(npx-1, j) - min27_tl
            z4 = q(npx-1, j) - min27
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
            dm1_tl(npx-1) = min4_tl*SIGN(1.d0, min4*dm1(npx-1))
            dm1(npx-1) = SIGN(min4, dm1(npx-1))
            al_tl(npx-1) = 0.5*(q_tl(npx-2, j)+q_tl(npx-1, j)) + r3*(&
&             dm1_tl(npx-2)-dm1_tl(npx-1))
            al(npx-1) = 0.5*(q(npx-2, j)+q(npx-1, j)) + r3*(dm1(npx-2)-&
&             dm1(npx-1))
!
            x1_tl = s15*q_tl(npx, j) + s11*q_tl(npx+1, j) - s14*dm1_tl(&
&             npx+1)
            x1 = s15*q(npx, j) + s11*q(npx+1, j) - s14*dm1(npx+1)
            dm1_tl(npx) = 0.5*(x1_tl-x0_tl)
            dm1(npx) = 0.5*(x1-x0)
            IF (dm1(npx) .GE. 0.) THEN
              x6_tl = dm1_tl(npx)
              x6 = dm1(npx)
            ELSE
              x6_tl = -dm1_tl(npx)
              x6 = -dm1(npx)
            END IF
            IF (q(npx, j) .LT. x0) THEN
              IF (x0 .LT. x1) THEN
                max6_tl = x1_tl
                max6 = x1
              ELSE
                max6_tl = x0_tl
                max6 = x0
              END IF
            ELSE IF (q(npx, j) .LT. x1) THEN
              max6_tl = x1_tl
              max6 = x1
            ELSE
              max6_tl = q_tl(npx, j)
              max6 = q(npx, j)
            END IF
            y5_tl = max6_tl - q_tl(npx, j)
            y5 = max6 - q(npx, j)
            IF (q(npx, j) .GT. x0) THEN
              IF (x0 .GT. x1) THEN
                min28_tl = x1_tl
                min28 = x1
              ELSE
                min28_tl = x0_tl
                min28 = x0
              END IF
            ELSE IF (q(npx, j) .GT. x1) THEN
              min28_tl = x1_tl
              min28 = x1
            ELSE
              min28_tl = q_tl(npx, j)
              min28 = q(npx, j)
            END IF
            z5_tl = q_tl(npx, j) - min28_tl
            z5 = q(npx, j) - min28
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
            dm1_tl(npx) = min5_tl*SIGN(1.d0, min5*dm1(npx))
            dm1(npx) = SIGN(min5, dm1(npx))
            al_tl(npx+1) = 0.5*(q_tl(npx, j)+q_tl(npx+1, j)) + r3*(&
&             dm1_tl(npx)-dm1_tl(npx+1))
            al(npx+1) = 0.5*(q(npx, j)+q(npx+1, j)) + r3*(dm1(npx)-dm1(&
&             npx+1))
          END IF
        ELSE
! For doubly periodic BC
          DO i=is-1,ie+2
            al_tl(i) = 0.5*(q_tl(i-1, j)+q_tl(i, j)) + r3*(dm1_tl(i-1)-&
&             dm1_tl(i))
            al(i) = 0.5*(q(i-1, j)+q(i, j)) + r3*(dm1(i-1)-dm1(i))
          END DO
        END IF
        IF (iord .EQ. 3) THEN
          DO i=is-1,ie+1
            bl_tl(i) = al_tl(i) - q_tl(i, j)
            bl(i) = al(i) - q(i, j)
            br_tl(i) = al_tl(i+1) - q_tl(i, j)
            br(i) = al(i+1) - q(i, j)
          END DO
          CALL PERT_PPM_TLM(ie - is + 3, q(is-1:, j), bl(is-1:), bl_tl(&
&                     is-1:), br(is-1:), br_tl(is-1:), 1)
          DO i=is,ie+1
            IF (c(i, j) .GT. 0.) THEN
              flux_tl(i, j) = q_tl(i-1, j) + (1.-c(i, j))*(br_tl(i-1)-&
&               c_tl(i, j)*(bl(i-1)+br(i-1))-c(i, j)*(bl_tl(i-1)+br_tl(i&
&               -1))) - c_tl(i, j)*(br(i-1)-c(i, j)*(bl(i-1)+br(i-1)))
              flux(i, j) = q(i-1, j) + (1.-c(i, j))*(br(i-1)-c(i, j)*(bl&
&               (i-1)+br(i-1)))
            ELSE
              flux_tl(i, j) = q_tl(i, j) + c_tl(i, j)*(bl(i)+c(i, j)*(bl&
&               (i)+br(i))) + (1.+c(i, j))*(bl_tl(i)+c_tl(i, j)*(bl(i)+&
&               br(i))+c(i, j)*(bl_tl(i)+br_tl(i)))
              flux(i, j) = q(i, j) + (1.+c(i, j))*(bl(i)+c(i, j)*(bl(i)+&
&               br(i)))
            END IF
          END DO
        ELSE
          DO i=is,ie+1
            IF (c(i, j) .GT. 0.) THEN
              xt_tl = ppm_limiter*dm1_tl(i-1)
              xt = ppm_limiter*dm1(i-1)
              IF (xt .GE. 0.) THEN
                x7_tl = xt_tl
                x7 = xt
              ELSE
                x7_tl = -xt_tl
                x7 = -xt
              END IF
              IF (al(i-1) - q(i-1, j) .GE. 0.) THEN
                y6_tl = al_tl(i-1) - q_tl(i-1, j)
                y6 = al(i-1) - q(i-1, j)
              ELSE
                y6_tl = -(al_tl(i-1)-q_tl(i-1, j))
                y6 = -(al(i-1)-q(i-1, j))
              END IF
              IF (x7 .GT. y6) THEN
                min6_tl = y6_tl
                min6 = y6
              ELSE
                min6_tl = x7_tl
                min6 = x7
              END IF
              dl_tl = min6_tl*SIGN(1.d0, min6*xt)
              dl = SIGN(min6, xt)
              IF (xt .GE. 0.) THEN
                x8_tl = xt_tl
                x8 = xt
              ELSE
                x8_tl = -xt_tl
                x8 = -xt
              END IF
              IF (al(i) - q(i-1, j) .GE. 0.) THEN
                y7_tl = al_tl(i) - q_tl(i-1, j)
                y7 = al(i) - q(i-1, j)
              ELSE
                y7_tl = -(al_tl(i)-q_tl(i-1, j))
                y7 = -(al(i)-q(i-1, j))
              END IF
              IF (x8 .GT. y7) THEN
                min7_tl = y7_tl
                min7 = y7
              ELSE
                min7_tl = x8_tl
                min7 = x8
              END IF
              dr_tl = min7_tl*SIGN(1.d0, min7*xt)
              dr = SIGN(min7, xt)
              flux_tl(i, j) = q_tl(i-1, j) + (1.-c(i, j))*(c_tl(i, j)*(&
&               dl-dr)+c(i, j)*(dl_tl-dr_tl)+dr_tl) - c_tl(i, j)*(c(i, j&
&               )*(dl-dr)+dr)
              flux(i, j) = q(i-1, j) + (1.-c(i, j))*(c(i, j)*(dl-dr)+dr)
            ELSE
              xt_tl = ppm_limiter*dm1_tl(i)
              xt = ppm_limiter*dm1(i)
              IF (xt .GE. 0.) THEN
                x9_tl = xt_tl
                x9 = xt
              ELSE
                x9_tl = -xt_tl
                x9 = -xt
              END IF
              IF (al(i) - q(i, j) .GE. 0.) THEN
                y8_tl = al_tl(i) - q_tl(i, j)
                y8 = al(i) - q(i, j)
              ELSE
                y8_tl = -(al_tl(i)-q_tl(i, j))
                y8 = -(al(i)-q(i, j))
              END IF
              IF (x9 .GT. y8) THEN
                min8_tl = y8_tl
                min8 = y8
              ELSE
                min8_tl = x9_tl
                min8 = x9
              END IF
              dl_tl = min8_tl*SIGN(1.d0, min8*xt)
              dl = SIGN(min8, xt)
              IF (xt .GE. 0.) THEN
                x10_tl = xt_tl
                x10 = xt
              ELSE
                x10_tl = -xt_tl
                x10 = -xt
              END IF
              IF (al(i+1) - q(i, j) .GE. 0.) THEN
                y9_tl = al_tl(i+1) - q_tl(i, j)
                y9 = al(i+1) - q(i, j)
              ELSE
                y9_tl = -(al_tl(i+1)-q_tl(i, j))
                y9 = -(al(i+1)-q(i, j))
              END IF
              IF (x10 .GT. y9) THEN
                min9_tl = y9_tl
                min9 = y9
              ELSE
                min9_tl = x10_tl
                min9 = x10
              END IF
              dr_tl = min9_tl*SIGN(1.d0, min9*xt)
              dr = SIGN(min9, xt)
              flux_tl(i, j) = q_tl(i, j) - c_tl(i, j)*(c(i, j)*(dl-dr)+&
&               dl) - (1.+c(i, j))*(c_tl(i, j)*(dl-dr)+c(i, j)*(dl_tl-&
&               dr_tl)+dl_tl)
              flux(i, j) = q(i, j) - (1.+c(i, j))*(c(i, j)*(dl-dr)+dl)
            END IF
          END DO
        END IF
      END DO
    ELSE IF (iord .EQ. 5) THEN
      dq_tl = 0.0
      al_tl = 0.0
      dm1_tl = 0.0
      bl_tl = 0.0
      br_tl = 0.0
! PPM with Hunyh's 2nd constraint
      DO j=jfirst,jlast
        DO i=ifirst-3,ilast+2
          dq_tl(i) = q_tl(i+1, j) - q_tl(i, j)
          dq(i) = q(i+1, j) - q(i, j)
        END DO
        DO i=ifirst-2,ilast+2
          xt_tl = 0.25*(q_tl(i+1, j)-q_tl(i-1, j))
          xt = 0.25*(q(i+1, j)-q(i-1, j))
          IF (xt .GE. 0.) THEN
            x11_tl = xt_tl
            x11 = xt
          ELSE
            x11_tl = -xt_tl
            x11 = -xt
          END IF
          IF (q(i-1, j) .LT. q(i, j)) THEN
            IF (q(i, j) .LT. q(i+1, j)) THEN
              max7_tl = q_tl(i+1, j)
              max7 = q(i+1, j)
            ELSE
              max7_tl = q_tl(i, j)
              max7 = q(i, j)
            END IF
          ELSE IF (q(i-1, j) .LT. q(i+1, j)) THEN
            max7_tl = q_tl(i+1, j)
            max7 = q(i+1, j)
          ELSE
            max7_tl = q_tl(i-1, j)
            max7 = q(i-1, j)
          END IF
          y10_tl = max7_tl - q_tl(i, j)
          y10 = max7 - q(i, j)
          IF (q(i-1, j) .GT. q(i, j)) THEN
            IF (q(i, j) .GT. q(i+1, j)) THEN
              min29_tl = q_tl(i+1, j)
              min29 = q(i+1, j)
            ELSE
              min29_tl = q_tl(i, j)
              min29 = q(i, j)
            END IF
          ELSE IF (q(i-1, j) .GT. q(i+1, j)) THEN
            min29_tl = q_tl(i+1, j)
            min29 = q(i+1, j)
          ELSE
            min29_tl = q_tl(i-1, j)
            min29 = q(i-1, j)
          END IF
          z6_tl = q_tl(i, j) - min29_tl
          z6 = q(i, j) - min29
          IF (x11 .GT. y10) THEN
            IF (y10 .GT. z6) THEN
              min10_tl = z6_tl
              min10 = z6
            ELSE
              min10_tl = y10_tl
              min10 = y10
            END IF
          ELSE IF (x11 .GT. z6) THEN
            min10_tl = z6_tl
            min10 = z6
          ELSE
            min10_tl = x11_tl
            min10 = x11
          END IF
          dm1_tl(i) = min10_tl*SIGN(1.d0, min10*xt)
          dm1(i) = SIGN(min10, xt)
        END DO
        DO i=ifirst-1,ilast+2
          al_tl(i) = 0.5*(q_tl(i-1, j)+q_tl(i, j)) + r3*(dm1_tl(i-1)-&
&           dm1_tl(i))
          al(i) = 0.5*(q(i-1, j)+q(i, j)) + r3*(dm1(i-1)-dm1(i))
        END DO
        DO i=ifirst-1,ilast+1
          pmp_tl = -(2.*dq_tl(i))
          pmp = -(2.*dq(i))
          lac_tl = pmp_tl + 1.5*dq_tl(i+1)
          lac = pmp + 1.5*dq(i+1)
          IF (0. .LT. pmp) THEN
            IF (pmp .LT. lac) THEN
              x12_tl = lac_tl
              x12 = lac
            ELSE
              x12_tl = pmp_tl
              x12 = pmp
            END IF
          ELSE IF (0. .LT. lac) THEN
            x12_tl = lac_tl
            x12 = lac
          ELSE
            x12 = 0.
            x12_tl = 0.0
          END IF
          IF (0. .GT. pmp) THEN
            IF (pmp .GT. lac) THEN
              y44_tl = lac_tl
              y44 = lac
            ELSE
              y44_tl = pmp_tl
              y44 = pmp
            END IF
          ELSE IF (0. .GT. lac) THEN
            y44_tl = lac_tl
            y44 = lac
          ELSE
            y44 = 0.
            y44_tl = 0.0
          END IF
          IF (al(i) - q(i, j) .LT. y44) THEN
            y11_tl = y44_tl
            y11 = y44
          ELSE
            y11_tl = al_tl(i) - q_tl(i, j)
            y11 = al(i) - q(i, j)
          END IF
          IF (x12 .GT. y11) THEN
            bl_tl(i) = y11_tl
            bl(i) = y11
          ELSE
            bl_tl(i) = x12_tl
            bl(i) = x12
          END IF
          pmp_tl = 2.*dq_tl(i-1)
          pmp = 2.*dq(i-1)
          lac_tl = pmp_tl - 1.5*dq_tl(i-2)
          lac = pmp - 1.5*dq(i-2)
          IF (0. .LT. pmp) THEN
            IF (pmp .LT. lac) THEN
              x13_tl = lac_tl
              x13 = lac
            ELSE
              x13_tl = pmp_tl
              x13 = pmp
            END IF
          ELSE IF (0. .LT. lac) THEN
            x13_tl = lac_tl
            x13 = lac
          ELSE
            x13 = 0.
            x13_tl = 0.0
          END IF
          IF (0. .GT. pmp) THEN
            IF (pmp .GT. lac) THEN
              y45_tl = lac_tl
              y45 = lac
            ELSE
              y45_tl = pmp_tl
              y45 = pmp
            END IF
          ELSE IF (0. .GT. lac) THEN
            y45_tl = lac_tl
            y45 = lac
          ELSE
            y45 = 0.
            y45_tl = 0.0
          END IF
          IF (al(i+1) - q(i, j) .LT. y45) THEN
            y12_tl = y45_tl
            y12 = y45
          ELSE
            y12_tl = al_tl(i+1) - q_tl(i, j)
            y12 = al(i+1) - q(i, j)
          END IF
          IF (x13 .GT. y12) THEN
            br_tl(i) = y12_tl
            br(i) = y12
          ELSE
            br_tl(i) = x13_tl
            br(i) = x13
          END IF
        END DO
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

    ELSE IF (iord .EQ. 6) THEN

      bl_tl = 0.0_8
      br_tl = 0.0_8
      DO j=jfirst,jlast
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
          bl_tl(i) = b5*q_tl(i-2, j) + b4*q_tl(i-1, j) + b3*q_tl(i, j) +&
&           b2*q_tl(i+1, j) + b1*q_tl(i+2, j)
          bl(i) = b5*q(i-2, j) + b4*q(i-1, j) + b3*q(i, j) + b2*q(i+1, j&
&           ) + b1*q(i+2, j)
          br_tl(i) = b1*q_tl(i-2, j) + b2*q_tl(i-1, j) + b3*q_tl(i, j) +&
&           b4*q_tl(i+1, j) + b5*q_tl(i+2, j)
          br(i) = b1*q(i-2, j) + b2*q(i-1, j) + b3*q(i, j) + b4*q(i+1, j&
&           ) + b5*q(i+2, j)
        END DO
        IF (grid_type .LT. 3) THEN
          IF (is .EQ. 1) THEN
            br_tl(2) = p1*(q_tl(2, j)+q_tl(3, j)) + p2*(q_tl(1, j)+q_tl(&
&             4, j)) - q_tl(2, j)
            br(2) = p1*(q(2, j)+q(3, j)) + p2*(q(1, j)+q(4, j)) - q(2, j&
&             )
            xt_tl = c3*q_tl(1, j) + c2*q_tl(2, j) + c1*q_tl(3, j)
            xt = c3*q(1, j) + c2*q(2, j) + c1*q(3, j)
            bl_tl(2) = xt_tl - q_tl(2, j)
            bl(2) = xt - q(2, j)
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
              br_tl(1) = xt_tl - q_tl(1, j)
              br(1) = xt - q(1, j)
              xt_tl = 0.5*((2.*dx(1, j)+dx(2, j))*(q_tl(0, j)+q_tl(1, j)&
&               )-dx(1, j)*(q_tl(-1, j)+q_tl(2, j)))/(dx(1, j)+dx(2, j))
              xt = 0.5*((2.*dx(1, j)+dx(2, j))*(q(0, j)+q(1, j))-dx(1, j&
&               )*(q(-1, j)+q(2, j)))/(dx(1, j)+dx(2, j))
              bl_tl(1) = xt_tl - q_tl(1, j)
              bl(1) = xt - q(1, j)
              br_tl(0) = xt_tl - q_tl(0, j)
              br(0) = xt - q(0, j)
              xt_tl = c1*q_tl(-2, j) + c2*q_tl(-1, j) + c3*q_tl(0, j)
              xt = c1*q(-2, j) + c2*q(-1, j) + c3*q(0, j)
              bl_tl(0) = xt_tl - q_tl(0, j)
              bl(0) = xt - q(0, j)
            END IF
          END IF
          IF (ie + 1 .EQ. npx) THEN
            bl_tl(npx-2) = p1*(q_tl(npx-2, j)+q_tl(npx-3, j)) + p2*(q_tl&
&             (npx-4, j)+q_tl(npx-1, j)) - q_tl(npx-2, j)
            bl(npx-2) = p1*(q(npx-2, j)+q(npx-3, j)) + p2*(q(npx-4, j)+q&
&             (npx-1, j)) - q(npx-2, j)
            xt_tl = c1*q_tl(npx-3, j) + c2*q_tl(npx-2, j) + c3*q_tl(npx-&
&             1, j)
            xt = c1*q(npx-3, j) + c2*q(npx-2, j) + c3*q(npx-1, j)
            br_tl(npx-2) = xt_tl - q_tl(npx-2, j)
            br(npx-2) = xt - q(npx-2, j)
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
              bl_tl(npx-1) = xt_tl - q_tl(npx-1, j)
              bl(npx-1) = xt - q(npx-1, j)
              xt_tl = 0.5*((2.*dx(npx-1, j)+dx(npx-2, j))*(q_tl(npx-1, j&
&               )+q_tl(npx, j))-dx(npx-1, j)*(q_tl(npx-2, j)+q_tl(npx+1&
&               , j)))/(dx(npx-1, j)+dx(npx-2, j))
              xt = 0.5*((2.*dx(npx-1, j)+dx(npx-2, j))*(q(npx-1, j)+q(&
&               npx, j))-dx(npx-1, j)*(q(npx-2, j)+q(npx+1, j)))/(dx(npx&
&               -1, j)+dx(npx-2, j))
              br_tl(npx-1) = xt_tl - q_tl(npx-1, j)
              br(npx-1) = xt - q(npx-1, j)
              bl_tl(npx) = xt_tl - q_tl(npx, j)
              bl(npx) = xt - q(npx, j)
              xt_tl = c3*q_tl(npx, j) + c2*q_tl(npx+1, j) + c1*q_tl(npx+&
&               2, j)
              xt = c3*q(npx, j) + c2*q(npx+1, j) + c1*q(npx+2, j)
              br_tl(npx) = xt_tl - q_tl(npx, j)
              br(npx) = xt - q(npx, j)
            END IF
          END IF
        END IF
        DO i=ifirst,ilast+1
          IF (c(i, j) .GT. 0.) THEN
            cfl_tl = c_tl(i, j)
            cfl = c(i, j)
            flux_tl(i, j) = q_tl(i-1, j) + (1.-cfl)*(br_tl(i-1)-cfl_tl*(&
&             bl(i-1)+br(i-1))-cfl*(bl_tl(i-1)+br_tl(i-1))) - cfl_tl*(br&
&             (i-1)-cfl*(bl(i-1)+br(i-1)))
            flux(i, j) = q(i-1, j) + (1.-cfl)*(br(i-1)-cfl*(bl(i-1)+br(i&
&             -1)))
          ELSE
            cfl_tl = c_tl(i, j)
            cfl = c(i, j)
            flux_tl(i, j) = q_tl(i, j) + cfl_tl*(bl(i)+cfl*(bl(i)+br(i))&
&             ) + (1.+cfl)*(bl_tl(i)+cfl_tl*(bl(i)+br(i))+cfl*(bl_tl(i)+&
&             br_tl(i)))
            flux(i, j) = q(i, j) + (1.+cfl)*(bl(i)+cfl*(bl(i)+br(i)))
          END IF
        END DO
      END DO

    ELSE IF (iord .GE. 8 .and. iord .LE. 10) THEN
      dq_tl = 0.0
      al_tl = 0.0
      dm1_tl = 0.0
      bl_tl = 0.0
      br_tl = 0.0
! iord=8, 9, 10
      DO j=jfirst,jlast
! grid_type check
        IF (grid_type .LT. 3) THEN
          DO i=is-3,ie+2
            dq_tl(i) = q_tl(i+1, j) - q_tl(i, j)
            dq(i) = q(i+1, j) - q(i, j)
          END DO
          DO i=is-2,ie+2
            xt_tl = 0.25*(q_tl(i+1, j)-q_tl(i-1, j))
            xt = 0.25*(q(i+1, j)-q(i-1, j))
            IF (xt .GE. 0.) THEN
              x18_tl = xt_tl
              x18 = xt
            ELSE
              x18_tl = -xt_tl
              x18 = -xt
            END IF
            IF (q(i-1, j) .LT. q(i, j)) THEN
              IF (q(i, j) .LT. q(i+1, j)) THEN
                max8_tl = q_tl(i+1, j)
                max8 = q(i+1, j)
              ELSE
                max8_tl = q_tl(i, j)
                max8 = q(i, j)
              END IF
            ELSE IF (q(i-1, j) .LT. q(i+1, j)) THEN
              max8_tl = q_tl(i+1, j)
              max8 = q(i+1, j)
            ELSE
              max8_tl = q_tl(i-1, j)
              max8 = q(i-1, j)
            END IF
            y25_tl = max8_tl - q_tl(i, j)
            y25 = max8 - q(i, j)
            IF (q(i-1, j) .GT. q(i, j)) THEN
              IF (q(i, j) .GT. q(i+1, j)) THEN
                min30_tl = q_tl(i+1, j)
                min30 = q(i+1, j)
              ELSE
                min30_tl = q_tl(i, j)
                min30 = q(i, j)
              END IF
            ELSE IF (q(i-1, j) .GT. q(i+1, j)) THEN
              min30_tl = q_tl(i+1, j)
              min30 = q(i+1, j)
            ELSE
              min30_tl = q_tl(i-1, j)
              min30 = q(i-1, j)
            END IF
            z7_tl = q_tl(i, j) - min30_tl
            z7 = q(i, j) - min30
            IF (x18 .GT. y25) THEN
              IF (y25 .GT. z7) THEN
                min13_tl = z7_tl
                min13 = z7
              ELSE
                min13_tl = y25_tl
                min13 = y25
              END IF
            ELSE IF (x18 .GT. z7) THEN
              min13_tl = z7_tl
              min13 = z7
            ELSE
              min13_tl = x18_tl
              min13 = x18
            END IF
            dm1_tl(i) = min13_tl*SIGN(1.d0, min13*xt)
            dm1(i) = SIGN(min13, xt)
          END DO
          IF (npx - 2 .GT. ie + 2) THEN
            min14 = ie + 2
          ELSE
            min14 = npx - 2
          END IF
          DO i=is3,min14
            al_tl(i) = 0.5*(q_tl(i-1, j)+q_tl(i, j)) + r3*(dm1_tl(i-1)-&
&             dm1_tl(i))
            al(i) = 0.5*(q(i-1, j)+q(i, j)) + r3*(dm1(i-1)-dm1(i))
          END DO
          IF (iord .EQ. 8) THEN
            DO i=is3,ie3
              xt_tl = 2.*dm1_tl(i)
              xt = 2.*dm1(i)
              IF (xt .GE. 0.) THEN
                x19_tl = xt_tl
                x19 = xt
              ELSE
                x19_tl = -xt_tl
                x19 = -xt
              END IF
              IF (al(i) - q(i, j) .GE. 0.) THEN
                y26_tl = al_tl(i) - q_tl(i, j)
                y26 = al(i) - q(i, j)
              ELSE
                y26_tl = -(al_tl(i)-q_tl(i, j))
                y26 = -(al(i)-q(i, j))
              END IF
              IF (x19 .GT. y26) THEN
                min15_tl = y26_tl
                min15 = y26
              ELSE
                min15_tl = x19_tl
                min15 = x19
              END IF
              bl_tl(i) = -(min15_tl*SIGN(1.d0, min15*xt))
              bl(i) = -SIGN(min15, xt)
              IF (xt .GE. 0.) THEN
                x20_tl = xt_tl
                x20 = xt
              ELSE
                x20_tl = -xt_tl
                x20 = -xt
              END IF
              IF (al(i+1) - q(i, j) .GE. 0.) THEN
                y27_tl = al_tl(i+1) - q_tl(i, j)
                y27 = al(i+1) - q(i, j)
              ELSE
                y27_tl = -(al_tl(i+1)-q_tl(i, j))
                y27 = -(al(i+1)-q(i, j))
              END IF
              IF (x20 .GT. y27) THEN
                min16_tl = y27_tl
                min16 = y27
              ELSE
                min16_tl = x20_tl
                min16 = x20
              END IF
              br_tl(i) = min16_tl*SIGN(1.d0, min16*xt)
              br(i) = SIGN(min16, xt)
            END DO
          ELSE IF (iord .EQ. 9) THEN
            DO i=is3,ie3
              pmp_tl = -(2.*dq_tl(i))
              pmp = -(2.*dq(i))
              lac_tl = pmp_tl + 1.5*dq_tl(i+1)
              lac = pmp + 1.5*dq(i+1)
              IF (0. .LT. pmp) THEN
                IF (pmp .LT. lac) THEN
                  x21_tl = lac_tl
                  x21 = lac
                ELSE
                  x21_tl = pmp_tl
                  x21 = pmp
                END IF
              ELSE IF (0. .LT. lac) THEN
                x21_tl = lac_tl
                x21 = lac
              ELSE
                x21 = 0.
                x21_tl = 0.0
              END IF
              IF (0. .GT. pmp) THEN
                IF (pmp .GT. lac) THEN
                  y48_tl = lac_tl
                  y48 = lac
                ELSE
                  y48_tl = pmp_tl
                  y48 = pmp
                END IF
              ELSE IF (0. .GT. lac) THEN
                y48_tl = lac_tl
                y48 = lac
              ELSE
                y48 = 0.
                y48_tl = 0.0
              END IF
              IF (al(i) - q(i, j) .LT. y48) THEN
                y28_tl = y48_tl
                y28 = y48
              ELSE
                y28_tl = al_tl(i) - q_tl(i, j)
                y28 = al(i) - q(i, j)
              END IF
              IF (x21 .GT. y28) THEN
                bl_tl(i) = y28_tl
                bl(i) = y28
              ELSE
                bl_tl(i) = x21_tl
                bl(i) = x21
              END IF
              pmp_tl = 2.*dq_tl(i-1)
              pmp = 2.*dq(i-1)
              lac_tl = pmp_tl - 1.5*dq_tl(i-2)
              lac = pmp - 1.5*dq(i-2)
              IF (0. .LT. pmp) THEN
                IF (pmp .LT. lac) THEN
                  x22_tl = lac_tl
                  x22 = lac
                ELSE
                  x22_tl = pmp_tl
                  x22 = pmp
                END IF
              ELSE IF (0. .LT. lac) THEN
                x22_tl = lac_tl
                x22 = lac
              ELSE
                x22 = 0.
                x22_tl = 0.0
              END IF
              IF (0. .GT. pmp) THEN
                IF (pmp .GT. lac) THEN
                  y49_tl = lac_tl
                  y49 = lac
                ELSE
                  y49_tl = pmp_tl
                  y49 = pmp
                END IF
              ELSE IF (0. .GT. lac) THEN
                y49_tl = lac_tl
                y49 = lac
              ELSE
                y49 = 0.
                y49_tl = 0.0
              END IF
              IF (al(i+1) - q(i, j) .LT. y49) THEN
                y29_tl = y49_tl
                y29 = y49
              ELSE
                y29_tl = al_tl(i+1) - q_tl(i, j)
                y29 = al(i+1) - q(i, j)
              END IF
              IF (x22 .GT. y29) THEN
                br_tl(i) = y29_tl
                br(i) = y29
              ELSE
                br_tl(i) = x22_tl
                br(i) = x22
              END IF
            END DO
          ELSE
            DO i=is3,ie3
              bl_tl(i) = al_tl(i) - q_tl(i, j)
              bl(i) = al(i) - q(i, j)
              br_tl(i) = al_tl(i+1) - q_tl(i, j)
              br(i) = al(i+1) - q(i, j)
              IF (dq(i-1)*dq(i) .LE. 0.) THEN
                pmp_tl = -(2.*dq_tl(i))
                pmp = -(2.*dq(i))
                lac_tl = pmp_tl + 1.5*dq_tl(i+1)
                lac = pmp + 1.5*dq(i+1)
                IF (0. .LT. pmp) THEN
                  IF (pmp .LT. lac) THEN
                    x23_tl = lac_tl
                    x23 = lac
                  ELSE
                    x23_tl = pmp_tl
                    x23 = pmp
                  END IF
                ELSE IF (0. .LT. lac) THEN
                  x23_tl = lac_tl
                  x23 = lac
                ELSE
                  x23 = 0.
                  x23_tl = 0.0
                END IF
                IF (0. .GT. pmp) THEN
                  IF (pmp .GT. lac) THEN
                    y50_tl = lac_tl
                    y50 = lac
                  ELSE
                    y50_tl = pmp_tl
                    y50 = pmp
                  END IF
                ELSE IF (0. .GT. lac) THEN
                  y50_tl = lac_tl
                  y50 = lac
                ELSE
                  y50 = 0.
                  y50_tl = 0.0
                END IF
                IF (bl(i) .LT. y50) THEN
                  y30_tl = y50_tl
                  y30 = y50
                ELSE
                  y30_tl = bl_tl(i)
                  y30 = bl(i)
                END IF
                IF (x23 .GT. y30) THEN
                  bl_tl(i) = y30_tl
                  bl(i) = y30
                ELSE
                  bl_tl(i) = x23_tl
                  bl(i) = x23
                END IF
                pmp_tl = 2.*dq_tl(i-1)
                pmp = 2.*dq(i-1)
                lac_tl = pmp_tl - 1.5*dq_tl(i-2)
                lac = pmp - 1.5*dq(i-2)
                IF (0. .LT. pmp) THEN
                  IF (pmp .LT. lac) THEN
                    x24_tl = lac_tl
                    x24 = lac
                  ELSE
                    x24_tl = pmp_tl
                    x24 = pmp
                  END IF
                ELSE IF (0. .LT. lac) THEN
                  x24_tl = lac_tl
                  x24 = lac
                ELSE
                  x24 = 0.
                  x24_tl = 0.0
                END IF
                IF (0. .GT. pmp) THEN
                  IF (pmp .GT. lac) THEN
                    y51_tl = lac_tl
                    y51 = lac
                  ELSE
                    y51_tl = pmp_tl
                    y51 = pmp
                  END IF
                ELSE IF (0. .GT. lac) THEN
                  y51_tl = lac_tl
                  y51 = lac
                ELSE
                  y51 = 0.
                  y51_tl = 0.0
                END IF
                IF (br(i) .LT. y51) THEN
                  y31_tl = y51_tl
                  y31 = y51
                ELSE
                  y31_tl = br_tl(i)
                  y31 = br(i)
                END IF
                IF (x24 .GT. y31) THEN
                  br_tl(i) = y31_tl
                  br(i) = y31
                ELSE
                  br_tl(i) = x24_tl
                  br(i) = x24
                END IF
              END IF
            END DO
          END IF
!--------------
! fix the edges
!--------------
          IF (is .EQ. 1) THEN
            br_tl(2) = al_tl(3) - q_tl(2, j)
            br(2) = al(3) - q(2, j)
!             xt = t11*(q(0,j)+q(1,j)) + t12*(q(-1,j)+q(2,j)) + t13*(dm1(2)-dm1(-1))
            xt_tl = 0.5*((2.*dxa(1, j)+dxa(2, j))*(q_tl(0, j)+q_tl(1, j)&
&             )-dxa(1, j)*(q_tl(-1, j)+q_tl(2, j)))/(dxa(1, j)+dxa(2, j)&
&             )
            xt = 0.5*((2.*dxa(1, j)+dxa(2, j))*(q(0, j)+q(1, j))-dxa(1, &
&             j)*(q(-1, j)+q(2, j)))/(dxa(1, j)+dxa(2, j))
            bl_tl(1) = xt_tl - q_tl(1, j)
            bl(1) = xt - q(1, j)
            br_tl(0) = xt_tl - q_tl(0, j)
            br(0) = xt - q(0, j)
            xt_tl = s14*dm1_tl(-1) - s11*dq_tl(-1) + q_tl(0, j)
            xt = s14*dm1(-1) - s11*dq(-1) + q(0, j)
!             xt = max( xt, min(q(-1,j),q(0,j)) )
!             xt = min( xt, max(q(-1,j),q(0,j)) )
            bl_tl(0) = xt_tl - q_tl(0, j)
            bl(0) = xt - q(0, j)
            xt_tl = s15*q_tl(1, j) + s11*q_tl(2, j) - s14*dm1_tl(2)
            xt = s15*q(1, j) + s11*q(2, j) - s14*dm1(2)
!             xt = max( xt, min(q(1,j),q(2,j)) )
!             xt = min( xt, max(q(1,j),q(2,j)) )
            br_tl(1) = xt_tl - q_tl(1, j)
            br(1) = xt - q(1, j)
            bl_tl(2) = xt_tl - q_tl(2, j)
            bl(2) = xt - q(2, j)
            CALL PERT_PPM_TLM(3, q(0:, j), bl(0:), bl_tl(0:), br(0:), &
&                       br_tl(0:), 1)
          END IF
          IF (ie + 1 .EQ. npx) THEN
            bl_tl(npx-2) = al_tl(npx-2) - q_tl(npx-2, j)
            bl(npx-2) = al(npx-2) - q(npx-2, j)
!             xt = t11*(q(npx-1,j)+q(npx,j)) + t12*(q(npx-2,j)+q(npx+1,j))   &
!                                            + t13*(dm1(npx+1)-dm1(npx-2))
            xt_tl = 0.5*((2.*dxa(npx-1, j)+dxa(npx-2, j))*(q_tl(npx-1, j&
&             )+q_tl(npx, j))-dxa(npx-1, j)*(q_tl(npx-2, j)+q_tl(npx+1, &
&             j)))/(dxa(npx-1, j)+dxa(npx-2, j))
            xt = 0.5*((2.*dxa(npx-1, j)+dxa(npx-2, j))*(q(npx-1, j)+q(&
&             npx, j))-dxa(npx-1, j)*(q(npx-2, j)+q(npx+1, j)))/(dxa(npx&
&             -1, j)+dxa(npx-2, j))
            br_tl(npx-1) = xt_tl - q_tl(npx-1, j)
            br(npx-1) = xt - q(npx-1, j)
            bl_tl(npx) = xt_tl - q_tl(npx, j)
            bl(npx) = xt - q(npx, j)
            xt_tl = s11*dq_tl(npx) - s14*dm1_tl(npx+1) + q_tl(npx, j)
            xt = s11*dq(npx) - s14*dm1(npx+1) + q(npx, j)
!             xt = min( xt, max(q(npx,j), q(npx+1,j)) )
!             xt = max( xt, min(q(npx,j), q(npx+1,j)) )
            br_tl(npx) = xt_tl - q_tl(npx, j)
            br(npx) = xt - q(npx, j)
            xt_tl = s15*q_tl(npx-1, j) + s11*q_tl(npx-2, j) + s14*dm1_tl&
&             (npx-2)
            xt = s15*q(npx-1, j) + s11*q(npx-2, j) + s14*dm1(npx-2)
!             xt = min( xt, max(q(npx-2,j), q(npx-1,j)) )
!             xt = max( xt, min(q(npx-2,j), q(npx-1,j)) )
            br_tl(npx-2) = xt_tl - q_tl(npx-2, j)
            br(npx-2) = xt - q(npx-2, j)
            bl_tl(npx-1) = xt_tl - q_tl(npx-1, j)
            bl(npx-1) = xt - q(npx-1, j)
            CALL PERT_PPM_TLM(3, q(npx-2:, j), bl(npx-2:), bl_tl(npx-2:)&
&                       , br(npx-2:), br_tl(npx-2:), 1)
          END IF
        ELSE
!---------------
! grid_type == 4
!---------------
          DO i=ifirst-2,ilast+2
            xt_tl = 0.25*(q_tl(i+1, j)-q_tl(i-1, j))
            xt = 0.25*(q(i+1, j)-q(i-1, j))
            IF (xt .GE. 0.) THEN
              x25_tl = xt_tl
              x25 = xt
            ELSE
              x25_tl = -xt_tl
              x25 = -xt
            END IF
            IF (q(i-1, j) .LT. q(i, j)) THEN
              IF (q(i, j) .LT. q(i+1, j)) THEN
                max9_tl = q_tl(i+1, j)
                max9 = q(i+1, j)
              ELSE
                max9_tl = q_tl(i, j)
                max9 = q(i, j)
              END IF
            ELSE IF (q(i-1, j) .LT. q(i+1, j)) THEN
              max9_tl = q_tl(i+1, j)
              max9 = q(i+1, j)
            ELSE
              max9_tl = q_tl(i-1, j)
              max9 = q(i-1, j)
            END IF
            y32_tl = max9_tl - q_tl(i, j)
            y32 = max9 - q(i, j)
            IF (q(i-1, j) .GT. q(i, j)) THEN
              IF (q(i, j) .GT. q(i+1, j)) THEN
                min31_tl = q_tl(i+1, j)
                min31 = q(i+1, j)
              ELSE
                min31_tl = q_tl(i, j)
                min31 = q(i, j)
              END IF
            ELSE IF (q(i-1, j) .GT. q(i+1, j)) THEN
              min31_tl = q_tl(i+1, j)
              min31 = q(i+1, j)
            ELSE
              min31_tl = q_tl(i-1, j)
              min31 = q(i-1, j)
            END IF
            z8_tl = q_tl(i, j) - min31_tl
            z8 = q(i, j) - min31
            IF (x25 .GT. y32) THEN
              IF (y32 .GT. z8) THEN
                min17_tl = z8_tl
                min17 = z8
              ELSE
                min17_tl = y32_tl
                min17 = y32
              END IF
            ELSE IF (x25 .GT. z8) THEN
              min17_tl = z8_tl
              min17 = z8
            ELSE
              min17_tl = x25_tl
              min17 = x25
            END IF
            dm1_tl(i) = min17_tl*SIGN(1.d0, min17*xt)
            dm1(i) = SIGN(min17, xt)
          END DO
          DO i=ifirst-1,ilast+2
            al_tl(i) = 0.5*(q_tl(i-1, j)+q_tl(i, j)) + r3*(dm1_tl(i-1)-&
&             dm1_tl(i))
            al(i) = 0.5*(q(i-1, j)+q(i, j)) + r3*(dm1(i-1)-dm1(i))
          END DO
          DO i=ifirst-3,ilast+2
            dq_tl(i) = q_tl(i+1, j) - q_tl(i, j)
            dq(i) = q(i+1, j) - q(i, j)
          END DO
          DO i=ifirst-1,ilast+1
            pmp_tl = -(2.*dq_tl(i))
            pmp = -(2.*dq(i))
            lac_tl = pmp_tl + 1.5*dq_tl(i+1)
            lac = pmp + 1.5*dq(i+1)
            IF (0. .LT. pmp) THEN
              IF (pmp .LT. lac) THEN
                x26_tl = lac_tl
                x26 = lac
              ELSE
                x26_tl = pmp_tl
                x26 = pmp
              END IF
            ELSE IF (0. .LT. lac) THEN
              x26_tl = lac_tl
              x26 = lac
            ELSE
              x26 = 0.
              x26_tl = 0.0
            END IF
            IF (0. .GT. pmp) THEN
              IF (pmp .GT. lac) THEN
                y52_tl = lac_tl
                y52 = lac
              ELSE
                y52_tl = pmp_tl
                y52 = pmp
              END IF
            ELSE IF (0. .GT. lac) THEN
              y52_tl = lac_tl
              y52 = lac
            ELSE
              y52 = 0.
              y52_tl = 0.0
            END IF
            IF (al(i) - q(i, j) .LT. y52) THEN
              y33_tl = y52_tl
              y33 = y52
            ELSE
              y33_tl = al_tl(i) - q_tl(i, j)
              y33 = al(i) - q(i, j)
            END IF
            IF (x26 .GT. y33) THEN
              bl_tl(i) = y33_tl
              bl(i) = y33
            ELSE
              bl_tl(i) = x26_tl
              bl(i) = x26
            END IF
            pmp_tl = 2.*dq_tl(i-1)
            pmp = 2.*dq(i-1)
            lac_tl = pmp_tl - 1.5*dq_tl(i-2)
            lac = pmp - 1.5*dq(i-2)
            IF (0. .LT. pmp) THEN
              IF (pmp .LT. lac) THEN
                x27_tl = lac_tl
                x27 = lac
              ELSE
                x27_tl = pmp_tl
                x27 = pmp
              END IF
            ELSE IF (0. .LT. lac) THEN
              x27_tl = lac_tl
              x27 = lac
            ELSE
              x27 = 0.
              x27_tl = 0.0
            END IF
            IF (0. .GT. pmp) THEN
              IF (pmp .GT. lac) THEN
                y53_tl = lac_tl
                y53 = lac
              ELSE
                y53_tl = pmp_tl
                y53 = pmp
              END IF
            ELSE IF (0. .GT. lac) THEN
              y53_tl = lac_tl
              y53 = lac
            ELSE
              y53 = 0.
              y53_tl = 0.0
            END IF
            IF (al(i+1) - q(i, j) .LT. y53) THEN
              y34_tl = y53_tl
              y34 = y53
            ELSE
              y34_tl = al_tl(i+1) - q_tl(i, j)
              y34 = al(i+1) - q(i, j)
            END IF
            IF (x27 .GT. y34) THEN
              br_tl(i) = y34_tl
              br(i) = y34
            ELSE
              br_tl(i) = x27_tl
              br(i) = x27
            END IF
          END DO
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
    ELSE
      dq_tl = 0.0
      al_tl = 0.0
      dm1_tl = 0.0
      bl_tl = 0.0
      br_tl = 0.0
!------------------------------
! For positive definite tracers:
!------------------------------
! iord=11: PPM mono constraint (Lin 2004)
! iord=12: Huynh 2nd constraint (Lin 2004) + positive definite (Lin & Rood 1996)
! iord>12: positive definite only (Lin & Rood 1996)
      DO j=jfirst,jlast
        DO i=is-2,ie+2
          xt_tl = 0.25*(q_tl(i+1, j)-q_tl(i-1, j))
          xt = 0.25*(q(i+1, j)-q(i-1, j))
          IF (xt .GE. 0.) THEN
            x28_tl = xt_tl
            x28 = xt
          ELSE
            x28_tl = -xt_tl
            x28 = -xt
          END IF
          IF (q(i-1, j) .LT. q(i, j)) THEN
            IF (q(i, j) .LT. q(i+1, j)) THEN
              max10_tl = q_tl(i+1, j)
              max10 = q(i+1, j)
            ELSE
              max10_tl = q_tl(i, j)
              max10 = q(i, j)
            END IF
          ELSE IF (q(i-1, j) .LT. q(i+1, j)) THEN
            max10_tl = q_tl(i+1, j)
            max10 = q(i+1, j)
          ELSE
            max10_tl = q_tl(i-1, j)
            max10 = q(i-1, j)
          END IF
          y35_tl = max10_tl - q_tl(i, j)
          y35 = max10 - q(i, j)
          IF (q(i-1, j) .GT. q(i, j)) THEN
            IF (q(i, j) .GT. q(i+1, j)) THEN
              min32_tl = q_tl(i+1, j)
              min32 = q(i+1, j)
            ELSE
              min32_tl = q_tl(i, j)
              min32 = q(i, j)
            END IF
          ELSE IF (q(i-1, j) .GT. q(i+1, j)) THEN
            min32_tl = q_tl(i+1, j)
            min32 = q(i+1, j)
          ELSE
            min32_tl = q_tl(i-1, j)
            min32 = q(i-1, j)
          END IF
          z9_tl = q_tl(i, j) - min32_tl
          z9 = q(i, j) - min32
          IF (x28 .GT. y35) THEN
            IF (y35 .GT. z9) THEN
              min18_tl = z9_tl
              min18 = z9
            ELSE
              min18_tl = y35_tl
              min18 = y35
            END IF
          ELSE IF (x28 .GT. z9) THEN
            min18_tl = z9_tl
            min18 = z9
          ELSE
            min18_tl = x28_tl
            min18 = x28
          END IF
          dm1_tl(i) = min18_tl*SIGN(1.d0, min18*xt)
          dm1(i) = SIGN(min18, xt)
        END DO
        IF (grid_type .LT. 3) THEN
          DO i=is3,ie2
            al_tl(i) = 0.5*(q_tl(i-1, j)+q_tl(i, j)) + r3*(dm1_tl(i-1)-&
&             dm1_tl(i))
            al(i) = 0.5*(q(i-1, j)+q(i, j)) + r3*(dm1(i-1)-dm1(i))
          END DO
          IF (iord .EQ. 11) THEN
            DO i=is3,ie3
              xt_tl = 2.*dm1_tl(i)
              xt = 2.*dm1(i)
              IF (xt .GE. 0.) THEN
                x29_tl = xt_tl
                x29 = xt
              ELSE
                x29_tl = -xt_tl
                x29 = -xt
              END IF
              IF (al(i) - q(i, j) .GE. 0.) THEN
                y36_tl = al_tl(i) - q_tl(i, j)
                y36 = al(i) - q(i, j)
              ELSE
                y36_tl = -(al_tl(i)-q_tl(i, j))
                y36 = -(al(i)-q(i, j))
              END IF
              IF (x29 .GT. y36) THEN
                min19_tl = y36_tl
                min19 = y36
              ELSE
                min19_tl = x29_tl
                min19 = x29
              END IF
              bl_tl(i) = -(min19_tl*SIGN(1.d0, min19*xt))
              bl(i) = -SIGN(min19, xt)
              IF (xt .GE. 0.) THEN
                x30_tl = xt_tl
                x30 = xt
              ELSE
                x30_tl = -xt_tl
                x30 = -xt
              END IF
              IF (al(i+1) - q(i, j) .GE. 0.) THEN
                y37_tl = al_tl(i+1) - q_tl(i, j)
                y37 = al(i+1) - q(i, j)
              ELSE
                y37_tl = -(al_tl(i+1)-q_tl(i, j))
                y37 = -(al(i+1)-q(i, j))
              END IF
              IF (x30 .GT. y37) THEN
                min20_tl = y37_tl
                min20 = y37
              ELSE
                min20_tl = x30_tl
                min20 = x30
              END IF
              br_tl(i) = min20_tl*SIGN(1.d0, min20*xt)
              br(i) = SIGN(min20, xt)
            END DO
          ELSE IF (iord .EQ. 12) THEN
            DO i=is-3,ie+2
              dq_tl(i) = q_tl(i+1, j) - q_tl(i, j)
              dq(i) = q(i+1, j) - q(i, j)
            END DO
            DO i=is3,ie3
              pmp_tl = -(2.*dq_tl(i))
              pmp = -(2.*dq(i))
              lac_tl = pmp_tl + 1.5*dq_tl(i+1)
              lac = pmp + 1.5*dq(i+1)
              IF (0. .LT. pmp) THEN
                IF (pmp .LT. lac) THEN
                  x31_tl = lac_tl
                  x31 = lac
                ELSE
                  x31_tl = pmp_tl
                  x31 = pmp
                END IF
              ELSE IF (0. .LT. lac) THEN
                x31_tl = lac_tl
                x31 = lac
              ELSE
                x31 = 0.
                x31_tl = 0.0
              END IF
              IF (0. .GT. pmp) THEN
                IF (pmp .GT. lac) THEN
                  y54_tl = lac_tl
                  y54 = lac
                ELSE
                  y54_tl = pmp_tl
                  y54 = pmp
                END IF
              ELSE IF (0. .GT. lac) THEN
                y54_tl = lac_tl
                y54 = lac
              ELSE
                y54 = 0.
                y54_tl = 0.0
              END IF
              IF (al(i) - q(i, j) .LT. y54) THEN
                y38_tl = y54_tl
                y38 = y54
              ELSE
                y38_tl = al_tl(i) - q_tl(i, j)
                y38 = al(i) - q(i, j)
              END IF
              IF (x31 .GT. y38) THEN
                bl_tl(i) = y38_tl
                bl(i) = y38
              ELSE
                bl_tl(i) = x31_tl
                bl(i) = x31
              END IF
              pmp_tl = 2.*dq_tl(i-1)
              pmp = 2.*dq(i-1)
              lac_tl = pmp_tl - 1.5*dq_tl(i-2)
              lac = pmp - 1.5*dq(i-2)
              IF (0. .LT. pmp) THEN
                IF (pmp .LT. lac) THEN
                  x32_tl = lac_tl
                  x32 = lac
                ELSE
                  x32_tl = pmp_tl
                  x32 = pmp
                END IF
              ELSE IF (0. .LT. lac) THEN
                x32_tl = lac_tl
                x32 = lac
              ELSE
                x32 = 0.
                x32_tl = 0.0
              END IF
              IF (0. .GT. pmp) THEN
                IF (pmp .GT. lac) THEN
                  y55_tl = lac_tl
                  y55 = lac
                ELSE
                  y55_tl = pmp_tl
                  y55 = pmp
                END IF
              ELSE IF (0. .GT. lac) THEN
                y55_tl = lac_tl
                y55 = lac
              ELSE
                y55 = 0.
                y55_tl = 0.0
              END IF
              IF (al(i+1) - q(i, j) .LT. y55) THEN
                y39_tl = y55_tl
                y39 = y55
              ELSE
                y39_tl = al_tl(i+1) - q_tl(i, j)
                y39 = al(i+1) - q(i, j)
              END IF
              IF (x32 .GT. y39) THEN
                br_tl(i) = y39_tl
                br(i) = y39
              ELSE
                br_tl(i) = x32_tl
                br(i) = x32
              END IF
            END DO
          ELSE
            DO i=is3,ie3
              bl_tl(i) = al_tl(i) - q_tl(i, j)
              bl(i) = al(i) - q(i, j)
              br_tl(i) = al_tl(i+1) - q_tl(i, j)
              br(i) = al(i+1) - q(i, j)
            END DO
          END IF
! Positive definite constraint:
          IF (iord .NE. 11) CALL PERT_PPM_TLM(ie3 - is3 + 1, q(is3:, j)&
&                                       , bl(is3:), bl_tl(is3:), br(is3:&
&                                       ), br_tl(is3:), 0)
!--------------
! fix the edges
!--------------
          IF (is .EQ. 1) THEN
            br_tl(2) = al_tl(3) - q_tl(2, j)
            br(2) = al(3) - q(2, j)
!               xt = t11*(q(0,j)+q(1,j)) + t12*(q(-1,j)+q(2,j)) + t13*(dm1(2)-dm1(-1))
!!!             xt = 0.75*(q(0,j)+q(1,j)) - 0.25*(q(-1,j)+q(2,j))
            xt_tl = 0.5*((2.*dxa(1, j)+dxa(2, j))*(q_tl(0, j)+q_tl(1, j)&
&             )-dxa(1, j)*(q_tl(-1, j)+q_tl(2, j)))/(dxa(1, j)+dxa(2, j)&
&             )
            xt = 0.5*((2.*dxa(1, j)+dxa(2, j))*(q(0, j)+q(1, j))-dxa(1, &
&             j)*(q(-1, j)+q(2, j)))/(dxa(1, j)+dxa(2, j))
            IF (0. .LT. xt) THEN
              xt = xt
            ELSE
              xt = 0.
              xt_tl = 0.0
            END IF
            bl_tl(1) = xt_tl - q_tl(1, j)
            bl(1) = xt - q(1, j)
            br_tl(0) = xt_tl - q_tl(0, j)
            br(0) = xt - q(0, j)
            xt_tl = 4.*dm1_tl(-1)/7. + 11.*q_tl(-1, j)/14. + 3.*q_tl(0, &
&             j)/14.
            xt = 4./7.*dm1(-1) + 11./14.*q(-1, j) + 3./14.*q(0, j)
            IF (0. .LT. xt) THEN
              xt = xt
            ELSE
              xt = 0.
              xt_tl = 0.0
            END IF
            bl_tl(0) = xt_tl - q_tl(0, j)
            bl(0) = xt - q(0, j)
            xt_tl = 3.*q_tl(1, j)/14. + 11.*q_tl(2, j)/14. - 4.*dm1_tl(2&
&             )/7.
            xt = 3./14.*q(1, j) + 11./14.*q(2, j) - 4./7.*dm1(2)
            IF (0. .LT. xt) THEN
              xt = xt
            ELSE
              xt = 0.
              xt_tl = 0.0
            END IF
            br_tl(1) = xt_tl - q_tl(1, j)
            br(1) = xt - q(1, j)
            bl_tl(2) = xt_tl - q_tl(2, j)
            bl(2) = xt - q(2, j)
            CALL PERT_PPM_TLM(3, q(0:, j), bl(0:), bl_tl(0:), br(0:), &
&                       br_tl(0:), 1)
          END IF
          IF (ie + 1 .EQ. npx) THEN
            bl_tl(npx-2) = al_tl(npx-2) - q_tl(npx-2, j)
            bl(npx-2) = al(npx-2) - q(npx-2, j)
!               xt = t11*(q(npx-1,j)+q(npx,j)) + t12*(q(npx-2,j)+q(npx+1,j))   &
!                  + t13*(dm1(npx+1)-dm1(npx-2))
!!!             xt = 0.75*(q(npx-1,j)+q(npx,j)) - 0.25*(q(npx-2,j)+q(npx+1,j))
            xt_tl = 0.5*((2.*dxa(npx-1, j)+dxa(npx-2, j))*(q_tl(npx-1, j&
&             )+q_tl(npx, j))-dxa(npx-1, j)*(q_tl(npx-2, j)+q_tl(npx+1, &
&             j)))/(dxa(npx-1, j)+dxa(npx-2, j))
            xt = 0.5*((2.*dxa(npx-1, j)+dxa(npx-2, j))*(q(npx-1, j)+q(&
&             npx, j))-dxa(npx-1, j)*(q(npx-2, j)+q(npx+1, j)))/(dxa(npx&
&             -1, j)+dxa(npx-2, j))
            IF (0. .LT. xt) THEN
              xt = xt
            ELSE
              xt = 0.
              xt_tl = 0.0
            END IF
            br_tl(npx-1) = xt_tl - q_tl(npx-1, j)
            br(npx-1) = xt - q(npx-1, j)
            bl_tl(npx) = xt_tl - q_tl(npx, j)
            bl(npx) = xt - q(npx, j)
!               br(npx) = 11./14.*q(npx+1,j) + 3./14.*q(npx,j) - 4./7.*dm1(npx+1)
            xt_tl = 11.*q_tl(npx+1, j)/14. + 3.*q_tl(npx, j)/14. - 4.*&
&             dm1_tl(npx+1)/7.
            xt = 11./14.*q(npx+1, j) + 3./14.*q(npx, j) - 4./7.*dm1(npx+&
&             1)
            IF (0. .LT. xt) THEN
              xt = xt
            ELSE
              xt = 0.
              xt_tl = 0.0
            END IF
            br_tl(npx) = xt_tl - q_tl(npx, j)
            br(npx) = xt - q(npx, j)
            xt_tl = 3.*q_tl(npx-1, j)/14. + 11.*q_tl(npx-2, j)/14. + 4.*&
&             dm1_tl(npx-2)/7.
            xt = 3./14.*q(npx-1, j) + 11./14.*q(npx-2, j) + 4./7.*dm1(&
&             npx-2)
            IF (0. .LT. xt) THEN
              xt = xt
            ELSE
              xt = 0.
              xt_tl = 0.0
            END IF
            br_tl(npx-2) = xt_tl - q_tl(npx-2, j)
            br(npx-2) = xt - q(npx-2, j)
            bl_tl(npx-1) = xt_tl - q_tl(npx-1, j)
            bl(npx-1) = xt - q(npx-1, j)
            CALL PERT_PPM_TLM(3, q(npx-2:, j), bl(npx-2:), bl_tl(npx-2:)&
&                       , br(npx-2:), br_tl(npx-2:), 1)
          END IF
        ELSE
!--------------
! grid_type >=4
!--------------
          DO i=ifirst-1,ilast+2
            al_tl(i) = 0.5*(q_tl(i-1, j)+q_tl(i, j)) + r3*(dm1_tl(i-1)-&
&             dm1_tl(i))
            al(i) = 0.5*(q(i-1, j)+q(i, j)) + r3*(dm1(i-1)-dm1(i))
          END DO
          IF (iord .EQ. 11) THEN
            DO i=ifirst-1,ilast+1
              xt_tl = 2.*dm1_tl(i)
              xt = 2.*dm1(i)
              IF (xt .GE. 0.) THEN
                x33_tl = xt_tl
                x33 = xt
              ELSE
                x33_tl = -xt_tl
                x33 = -xt
              END IF
              IF (al(i) - q(i, j) .GE. 0.) THEN
                y40_tl = al_tl(i) - q_tl(i, j)
                y40 = al(i) - q(i, j)
              ELSE
                y40_tl = -(al_tl(i)-q_tl(i, j))
                y40 = -(al(i)-q(i, j))
              END IF
              IF (x33 .GT. y40) THEN
                min21_tl = y40_tl
                min21 = y40
              ELSE
                min21_tl = x33_tl
                min21 = x33
              END IF
              bl_tl(i) = -(min21_tl*SIGN(1.d0, min21*xt))
              bl(i) = -SIGN(min21, xt)
              IF (xt .GE. 0.) THEN
                x34_tl = xt_tl
                x34 = xt
              ELSE
                x34_tl = -xt_tl
                x34 = -xt
              END IF
              IF (al(i+1) - q(i, j) .GE. 0.) THEN
                y41_tl = al_tl(i+1) - q_tl(i, j)
                y41 = al(i+1) - q(i, j)
              ELSE
                y41_tl = -(al_tl(i+1)-q_tl(i, j))
                y41 = -(al(i+1)-q(i, j))
              END IF
              IF (x34 .GT. y41) THEN
                min22_tl = y41_tl
                min22 = y41
              ELSE
                min22_tl = x34_tl
                min22 = x34
              END IF
              br_tl(i) = min22_tl*SIGN(1.d0, min22*xt)
              br(i) = SIGN(min22, xt)
            END DO
          ELSE IF (iord .EQ. 12) THEN
            DO i=ifirst-3,ilast+2
              dq_tl(i) = q_tl(i+1, j) - q_tl(i, j)
              dq(i) = q(i+1, j) - q(i, j)
            END DO
            DO i=ifirst-1,ilast+1
              pmp_tl = -(2.*dq_tl(i))
              pmp = -(2.*dq(i))
              lac_tl = pmp_tl + 1.5*dq_tl(i+1)
              lac = pmp + 1.5*dq(i+1)
              IF (0. .LT. pmp) THEN
                IF (pmp .LT. lac) THEN
                  x35_tl = lac_tl
                  x35 = lac
                ELSE
                  x35_tl = pmp_tl
                  x35 = pmp
                END IF
              ELSE IF (0. .LT. lac) THEN
                x35_tl = lac_tl
                x35 = lac
              ELSE
                x35 = 0.
                x35_tl = 0.0
              END IF
              IF (0. .GT. pmp) THEN
                IF (pmp .GT. lac) THEN
                  y56_tl = lac_tl
                  y56 = lac
                ELSE
                  y56_tl = pmp_tl
                  y56 = pmp
                END IF
              ELSE IF (0. .GT. lac) THEN
                y56_tl = lac_tl
                y56 = lac
              ELSE
                y56 = 0.
                y56_tl = 0.0
              END IF
              IF (al(i) - q(i, j) .LT. y56) THEN
                y42_tl = y56_tl
                y42 = y56
              ELSE
                y42_tl = al_tl(i) - q_tl(i, j)
                y42 = al(i) - q(i, j)
              END IF
              IF (x35 .GT. y42) THEN
                bl_tl(i) = y42_tl
                bl(i) = y42
              ELSE
                bl_tl(i) = x35_tl
                bl(i) = x35
              END IF
              pmp_tl = 2.*dq_tl(i-1)
              pmp = 2.*dq(i-1)
              lac_tl = pmp_tl - 1.5*dq_tl(i-2)
              lac = pmp - 1.5*dq(i-2)
              IF (0. .LT. pmp) THEN
                IF (pmp .LT. lac) THEN
                  x36_tl = lac_tl
                  x36 = lac
                ELSE
                  x36_tl = pmp_tl
                  x36 = pmp
                END IF
              ELSE IF (0. .LT. lac) THEN
                x36_tl = lac_tl
                x36 = lac
              ELSE
                x36 = 0.
                x36_tl = 0.0
              END IF
              IF (0. .GT. pmp) THEN
                IF (pmp .GT. lac) THEN
                  y57_tl = lac_tl
                  y57 = lac
                ELSE
                  y57_tl = pmp_tl
                  y57 = pmp
                END IF
              ELSE IF (0. .GT. lac) THEN
                y57_tl = lac_tl
                y57 = lac
              ELSE
                y57 = 0.
                y57_tl = 0.0
              END IF
              IF (al(i+1) - q(i, j) .LT. y57) THEN
                y43_tl = y57_tl
                y43 = y57
              ELSE
                y43_tl = al_tl(i+1) - q_tl(i, j)
                y43 = al(i+1) - q(i, j)
              END IF
              IF (x36 .GT. y43) THEN
                br_tl(i) = y43_tl
                br(i) = y43
              ELSE
                br_tl(i) = x36_tl
                br(i) = x36
              END IF
            END DO
          ELSE
            DO i=is-1,ie+1
              bl_tl(i) = al_tl(i) - q_tl(i, j)
              bl(i) = al(i) - q(i, j)
              br_tl(i) = al_tl(i+1) - q_tl(i, j)
              br(i) = al(i+1) - q(i, j)
            END DO
          END IF
! Positive definite constraint:
          IF (iord .NE. 11) CALL PERT_PPM_TLM(ie - is + 3, q(is-1:, j), &
&                                       bl(is-1:), bl_tl(is-1:), br(is-1&
&                                       :), br_tl(is-1:), 0)
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

  SUBROUTINE PERT_PPM_TLM(im, a0, al, al_tl, ar, ar_tl, iv)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: im
    INTEGER, INTENT(IN) :: iv
    REAL*8, INTENT(IN) :: a0(im)
    REAL*8, INTENT(INOUT) :: al(im), ar(im)
    REAL*8, INTENT(INOUT) :: al_tl(im), ar_tl(im)
! Local:
    REAL*8 :: a4, da1, da2, a6da, fmin
    INTEGER :: i
    REAL*8, PARAMETER :: r12=1./12.
    INTRINSIC ABS
    REAL*8 :: abs0
!-----------------------------------
! Optimized PPM in perturbation form:
!-----------------------------------
    IF (iv .EQ. 0) THEN
! Positive definite constraint
      DO i=1,im
        a4 = -(3.*(ar(i)+al(i)))
        da1 = ar(i) - al(i)
        IF (da1 .GE. 0.) THEN
          abs0 = da1
        ELSE
          abs0 = -da1
        END IF
        IF (abs0 .LT. -a4) THEN
          fmin = a0(i) + 0.25/a4*da1**2 + a4*r12
          IF (fmin .LT. 0.) THEN
            IF (ar(i) .GT. 0. .AND. al(i) .GT. 0.) THEN
              ar_tl(i) = 0.0_8
              ar(i) = 0.
              al_tl(i) = 0.0_8
              al(i) = 0.
            ELSE IF (da1 .GT. 0.) THEN
              ar_tl(i) = -(2.*al_tl(i))
              ar(i) = -(2.*al(i))
            ELSE
              al_tl(i) = -(2.*ar_tl(i))
              al(i) = -(2.*ar(i))
            END IF
          END IF
        END IF
      END DO
    ELSE
! Standard PPM constraint
      DO i=1,im
        IF (al(i)*ar(i) .LT. 0.) THEN
          da1 = al(i) - ar(i)
          da2 = da1**2
          a6da = 3.*(al(i)+ar(i))*da1
          IF (a6da .LT. -da2) THEN
            ar_tl(i) = -(2.*al_tl(i))
            ar(i) = -(2.*al(i))
          ELSE IF (a6da .GT. da2) THEN
            al_tl(i) = -(2.*ar_tl(i))
            al(i) = -(2.*ar(i))
          END IF
        ELSE
! effect of dm=0 included here
          al_tl(i) = 0.0_8
          al(i) = 0.
          ar_tl(i) = 0.0_8
          ar(i) = 0.
        END IF
      END DO
    END IF
  END SUBROUTINE PERT_PPM_TLM

end module tp_core_tlm_mod
