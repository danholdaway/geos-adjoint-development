 module fv_grid_utils_adm_mod
 
 use fv_mp_mod,         only: is,js,ie,je, isd,jsd,ied,jed
 use fv_grid_utils_mod, only: a11, a12, a21, a22
 use fv_grid_tools_mod, only: g_type=>grid_type

 implicit none
 private

 public C2L_ORD2, C2L_ORD2_ADM

 contains

  SUBROUTINE C2L_ORD2_ADM(u, u_ad, v, v_ad, ua, ua_ad, va, va_ad, dx, dy&
&   , rdxa, rdya, km)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: km
    REAL*8, INTENT(IN) :: u(isd:ied, jsd:jed+1, km)
    REAL*8 :: u_ad(isd:ied, jsd:jed+1, km)
    REAL*8, INTENT(IN) :: v(isd:ied+1, jsd:jed, km)
    REAL*8 :: v_ad(isd:ied+1, jsd:jed, km)
    REAL, INTENT(IN) :: dx(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: dy(isd:ied+1, jsd:jed)
    REAL, INTENT(IN) :: rdxa(isd:ied, jsd:jed)
    REAL, INTENT(IN) :: rdya(isd:ied, jsd:jed)
!
    REAL*8 :: ua(isd:ied, jsd:jed, km)
    REAL*8 :: ua_ad(isd:ied, jsd:jed, km)
    REAL*8 :: va(isd:ied, jsd:jed, km)
    REAL*8 :: va_ad(isd:ied, jsd:jed, km)
!--------------------------------------------------------------
! Local 
    REAL :: wu(is:ie, js:je+1)
    REAL :: wu_ad(is:ie, js:je+1)
    REAL :: wv(is:ie+1, js:je)
    REAL :: wv_ad(is:ie+1, js:je)
    REAL :: u1(is:ie), v1(is:ie)
    REAL :: u1_ad(is:ie), v1_ad(is:ie)
    INTEGER :: i, j, k
    INTEGER :: branch
    REAL :: temp_ad0
    REAL :: temp_ad
    DO k=1,km
      IF (g_type .LT. 4) THEN
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL PUSHCONTROL1B(0)
      END IF
    END DO
    wu_ad = 0.0_8
    v1_ad = 0.0_8
    wv_ad = 0.0_8
    u1_ad = 0.0_8
    DO k=km,1,-1
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO j=je,js,-1
          DO i=ie,is,-1
            v_ad(i, j, k) = v_ad(i, j, k) + 0.5*va_ad(i, j, k)
            v_ad(i+1, j, k) = v_ad(i+1, j, k) + 0.5*va_ad(i, j, k)
            va_ad(i, j, k) = 0.0_8
            u_ad(i, j, k) = u_ad(i, j, k) + 0.5*ua_ad(i, j, k)
            u_ad(i, j+1, k) = u_ad(i, j+1, k) + 0.5*ua_ad(i, j, k)
            ua_ad(i, j, k) = 0.0_8
          END DO
        END DO
      ELSE
        DO j=je,js,-1
          DO i=ie,is,-1
            u1_ad(i) = u1_ad(i) + a11(i, j)*ua_ad(i, j, k) + a21(i, j)*&
&             va_ad(i, j, k)
            v1_ad(i) = v1_ad(i) + a12(i, j)*ua_ad(i, j, k) + a22(i, j)*&
&             va_ad(i, j, k)
            va_ad(i, j, k) = 0.0_8
            ua_ad(i, j, k) = 0.0_8
            temp_ad = rdya(i, j)*v1_ad(i)
            wv_ad(i, j) = wv_ad(i, j) + temp_ad
            wv_ad(i+1, j) = wv_ad(i+1, j) + temp_ad
            v1_ad(i) = 0.0_8
            temp_ad0 = rdxa(i, j)*u1_ad(i)
            wu_ad(i, j) = wu_ad(i, j) + temp_ad0
            wu_ad(i, j+1) = wu_ad(i, j+1) + temp_ad0
            u1_ad(i) = 0.0_8
          END DO
        END DO
        DO j=je,js,-1
          DO i=ie+1,is,-1
            v_ad(i, j, k) = v_ad(i, j, k) + dy(i, j)*wv_ad(i, j)
            wv_ad(i, j) = 0.0_8
          END DO
        END DO
        DO j=je+1,js,-1
          DO i=ie,is,-1
            u_ad(i, j, k) = u_ad(i, j, k) + dx(i, j)*wu_ad(i, j)
            wu_ad(i, j) = 0.0_8
          END DO
        END DO
      END IF
    END DO
  END SUBROUTINE C2L_ORD2_ADM
  SUBROUTINE C2L_ORD2(u, v, ua, va, dx, dy, rdxa, rdya, km)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: km
    REAL*8, INTENT(IN) :: u(isd:ied, jsd:jed+1, km)
    REAL*8, INTENT(IN) :: v(isd:ied+1, jsd:jed, km)
    REAL, INTENT(IN) :: dx(isd:ied, jsd:jed+1)
    REAL, INTENT(IN) :: dy(isd:ied+1, jsd:jed)
    REAL, INTENT(IN) :: rdxa(isd:ied, jsd:jed)
    REAL, INTENT(IN) :: rdya(isd:ied, jsd:jed)
!
    REAL*8, INTENT(OUT) :: ua(isd:ied, jsd:jed, km)
    REAL*8, INTENT(OUT) :: va(isd:ied, jsd:jed, km)
!--------------------------------------------------------------
! Local 
    REAL :: wu(is:ie, js:je+1)
    REAL :: wv(is:ie+1, js:je)
    REAL :: u1(is:ie), v1(is:ie)
    INTEGER :: i, j, k
    DO k=1,km
      IF (g_type .LT. 4) THEN
        DO j=js,je+1
          DO i=is,ie
            wu(i, j) = u(i, j, k)*dx(i, j)
          END DO
        END DO
        DO j=js,je
          DO i=is,ie+1
            wv(i, j) = v(i, j, k)*dy(i, j)
          END DO
        END DO
        DO j=js,je
          DO i=is,ie
! Co-variant to Co-variant "vorticity-conserving" interpolation
            u1(i) = (wu(i, j)+wu(i, j+1))*rdxa(i, j)
            v1(i) = (wv(i, j)+wv(i+1, j))*rdya(i, j)
! Cubed (cell center co-variant winds) to lat-lon:
            ua(i, j, k) = a11(i, j)*u1(i) + a12(i, j)*v1(i)
            va(i, j, k) = a21(i, j)*u1(i) + a22(i, j)*v1(i)
          END DO
        END DO
      ELSE
! 2nd order:
        DO j=js,je
          DO i=is,ie
            ua(i, j, k) = 0.5*(u(i, j, k)+u(i, j+1, k))
            va(i, j, k) = 0.5*(v(i, j, k)+v(i+1, j, k))
          END DO
        END DO
      END IF
    END DO
  END SUBROUTINE C2L_ORD2

 end module fv_grid_utils_adm_mod
