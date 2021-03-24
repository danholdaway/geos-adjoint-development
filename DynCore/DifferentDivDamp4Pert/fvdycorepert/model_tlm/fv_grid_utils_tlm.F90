 module fv_grid_utils_tlm_mod
 
 use fv_mp_mod,         only: is,js,ie,je, isd,jsd,ied,jed
 use fv_grid_utils_mod, only: a11, a12, a21, a22
 use fv_grid_tools_mod, only: g_type=>grid_type

 implicit none
 private

 public C2L_ORD2_TLM

 contains

  SUBROUTINE C2L_ORD2_TLM(u, u_tl, v, v_tl, ua, ua_tl, va, va_tl, dx, dy&
&   , rdxa, rdya, km)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: km
    REAL*8, INTENT(IN) :: u(isd:ied, jsd:jed+1, km)
    REAL*8, INTENT(IN) :: u_tl(isd:ied, jsd:jed+1, km)
    REAL*8, INTENT(IN) :: v(isd:ied+1, jsd:jed, km)
    REAL*8, INTENT(IN) :: v_tl(isd:ied+1, jsd:jed, km)
    REAL, INTENT(IN) :: dx(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: dy(isd:ied+1, jsd:jed)
    REAL, INTENT(IN) :: rdxa(isd:ied, jsd:jed)
    REAL, INTENT(IN) :: rdya(isd:ied, jsd:jed)
!
    REAL*8, INTENT(OUT) :: ua(isd:ied, jsd:jed, km)
    REAL*8, INTENT(OUT) :: ua_tl(isd:ied, jsd:jed, km)
    REAL*8, INTENT(OUT) :: va(isd:ied, jsd:jed, km)
    REAL*8, INTENT(OUT) :: va_tl(isd:ied, jsd:jed, km)
!--------------------------------------------------------------
! Local 
    REAL :: wu(is:ie, js:je+1)
    REAL :: wu_tl(is:ie, js:je+1)
    REAL :: wv(is:ie+1, js:je)
    REAL :: wv_tl(is:ie+1, js:je)
    REAL :: u1(is:ie), v1(is:ie)
    REAL :: u1_tl(is:ie), v1_tl(is:ie)
    INTEGER :: i, j, k
    ua_tl = 0.0_8
    va_tl = 0.0_8
    wu_tl = 0.0
    v1_tl = 0.0
    wv_tl = 0.0
    u1_tl = 0.0
    DO k=1,km
      IF (g_type .LT. 4) THEN
        DO j=js,je+1
          DO i=is,ie
            wu_tl(i, j) = dx(i, j)*u_tl(i, j, k)
            wu(i, j) = u(i, j, k)*dx(i, j)
          END DO
        END DO
        DO j=js,je
          DO i=is,ie+1
            wv_tl(i, j) = dy(i, j)*v_tl(i, j, k)
            wv(i, j) = v(i, j, k)*dy(i, j)
          END DO
        END DO
        DO j=js,je
          DO i=is,ie
! Co-variant to Co-variant "vorticity-conserving" interpolation
            u1_tl(i) = rdxa(i, j)*(wu_tl(i, j)+wu_tl(i, j+1))
            u1(i) = (wu(i, j)+wu(i, j+1))*rdxa(i, j)
            v1_tl(i) = rdya(i, j)*(wv_tl(i, j)+wv_tl(i+1, j))
            v1(i) = (wv(i, j)+wv(i+1, j))*rdya(i, j)
! Cubed (cell center co-variant winds) to lat-lon:
            ua_tl(i, j, k) = a11(i, j)*u1_tl(i) + a12(i, j)*v1_tl(i)
            ua(i, j, k) = a11(i, j)*u1(i) + a12(i, j)*v1(i)
            va_tl(i, j, k) = a21(i, j)*u1_tl(i) + a22(i, j)*v1_tl(i)
            va(i, j, k) = a21(i, j)*u1(i) + a22(i, j)*v1(i)
          END DO
        END DO
      ELSE
! 2nd order:
        DO j=js,je
          DO i=is,ie
            ua_tl(i, j, k) = 0.5*(u_tl(i, j, k)+u_tl(i, j+1, k))
            ua(i, j, k) = 0.5*(u(i, j, k)+u(i, j+1, k))
            va_tl(i, j, k) = 0.5*(v_tl(i, j, k)+v_tl(i+1, j, k))
            va(i, j, k) = 0.5*(v(i, j, k)+v(i+1, j, k))
          END DO
        END DO
      END IF
    END DO
  END SUBROUTINE C2L_ORD2_TLM

 end module fv_grid_utils_tlm_mod
