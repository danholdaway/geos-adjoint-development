module fv_fill_adm_mod

 use fv_arrays_mod, only: p_precision

 implicit none
 public fillz_adm

 contains

  SUBROUTINE FILLZ_ADM(im, km, nq, q, q_ad, dp, dp_ad)
    IMPLICIT NONE
! No. of longitudes
    INTEGER, INTENT(IN) :: im
! No. of levels
    INTEGER, INTENT(IN) :: km
! Total number of tracers
    INTEGER, INTENT(IN) :: nq
! pressure thickness
    REAL(p_precision), INTENT(IN) :: dp(im, km)
    REAL(p_precision) :: dp_ad(im, km)
! tracer mixing ratio
    REAL(p_precision), INTENT(INOUT) :: q(im, km, nq)
    REAL(p_precision) :: q_ad(im, km, nq)
! !LOCAL VARIABLES:
    INTEGER :: i, k, ic
    REAL :: qup, qly, dup
    REAL :: qup_ad, qly_ad, dup_ad
    INTRINSIC MAX
    INTRINSIC MIN
    INTEGER :: branch
    REAL*8 :: temp_ad3
    REAL*8 :: temp_ad2
    REAL*8 :: temp_ad1
    REAL*8 :: temp_ad0
    REAL(p_precision) :: temp_ad
    REAL(p_precision) :: temp
    DO ic=1,nq
! Top layer
      DO i=1,im
        IF (q(i, 1, ic) .LT. 0.) THEN
          CALL PUSHREAL8(q(i, 2, ic))
          q(i, 2, ic) = q(i, 2, ic) + q(i, 1, ic)*dp(i, 1)/dp(i, 2)
          CALL PUSHREAL8(q(i, 1, ic))
          q(i, 1, ic) = 0.
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
      END DO
      CALL PUSHINTEGER4(k)
! Interior
      DO k=2,km-1
        DO i=1,im
          IF (q(i, k, ic) .LT. 0.) THEN
            IF (0. .LT. q(i, k-1, ic)*dp(i, k-1)) THEN
              qup = q(i, k-1, ic)*dp(i, k-1)
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHCONTROL1B(1)
              qup = 0.
            END IF
            CALL PUSHREAL4(qly)
            qly = -(q(i, k, ic)*dp(i, k))
            IF (0.5*qly .GT. 0.99*qup) THEN
              CALL PUSHREAL4(dup)
              dup = 0.99*qup
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHREAL4(dup)
              dup = 0.5*qly
              CALL PUSHCONTROL1B(1)
            END IF
            CALL PUSHREAL8(q(i, k-1, ic))
            q(i, k-1, ic) = q(i, k-1, ic) - dup/dp(i, k-1)
! Borrow from below: q(i,k,ic) is still negative at this stage
            CALL PUSHREAL8(q(i, k+1, ic))
            q(i, k+1, ic) = q(i, k+1, ic) - (qly-dup)/dp(i, k+1)
            CALL PUSHREAL8(q(i, k, ic))
            q(i, k, ic) = 0.
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
      END DO
! Bottom layer
      k = km
      DO i=1,im
        IF (q(i, k, ic) .LT. 0. .AND. q(i, k-1, ic) .GT. 0.) THEN
! Borrow from above
          qup = q(i, k-1, ic)*dp(i, k-1)
          CALL PUSHREAL4(qly)
          qly = -(q(i, k, ic)*dp(i, k))
          IF (qly .GT. qup) THEN
            CALL PUSHREAL4(dup)
            dup = qup
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHREAL4(dup)
            dup = qly
            CALL PUSHCONTROL1B(1)
          END IF
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
      END DO
    END DO
    DO ic=nq,1,-1
      DO i=im,1,-1
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) THEN
          temp_ad3 = -(q_ad(i, k-1, ic)/dp(i, k-1))
          temp_ad2 = q_ad(i, k, ic)/dp(i, k)
          dup_ad = temp_ad3 + temp_ad2
          dp_ad(i, k) = dp_ad(i, k) - dup*temp_ad2/dp(i, k)
          dp_ad(i, k-1) = dp_ad(i, k-1) - dup*temp_ad3/dp(i, k-1)
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL4(dup)
            qup_ad = dup_ad
            qly_ad = 0.0
          ELSE
            CALL POPREAL4(dup)
            qly_ad = dup_ad
            qup_ad = 0.0
          END IF
          CALL POPREAL4(qly)
          q_ad(i, k, ic) = q_ad(i, k, ic) - dp(i, k)*qly_ad
          dp_ad(i, k) = dp_ad(i, k) - q(i, k, ic)*qly_ad
          q_ad(i, k-1, ic) = q_ad(i, k-1, ic) + dp(i, k-1)*qup_ad
          dp_ad(i, k-1) = dp_ad(i, k-1) + q(i, k-1, ic)*qup_ad
        END IF
      END DO
      DO k=km-1,2,-1
        DO i=im,1,-1
          CALL POPCONTROL1B(branch)
          IF (branch .NE. 0) THEN
            CALL POPREAL8(q(i, k, ic))
            q_ad(i, k, ic) = 0.0_8
            CALL POPREAL8(q(i, k+1, ic))
            temp_ad0 = -(q_ad(i, k+1, ic)/dp(i, k+1))
            qly_ad = temp_ad0
            dp_ad(i, k+1) = dp_ad(i, k+1) - (qly-dup)*temp_ad0/dp(i, k+1&
&             )
            CALL POPREAL8(q(i, k-1, ic))
            temp_ad1 = -(q_ad(i, k-1, ic)/dp(i, k-1))
            dup_ad = temp_ad1 - temp_ad0
            dp_ad(i, k-1) = dp_ad(i, k-1) - dup*temp_ad1/dp(i, k-1)
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              CALL POPREAL4(dup)
              qup_ad = 0.99*dup_ad
            ELSE
              CALL POPREAL4(dup)
              qly_ad = qly_ad + 0.5*dup_ad
              qup_ad = 0.0
            END IF
            CALL POPREAL4(qly)
            q_ad(i, k, ic) = q_ad(i, k, ic) - dp(i, k)*qly_ad
            dp_ad(i, k) = dp_ad(i, k) - q(i, k, ic)*qly_ad
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              q_ad(i, k-1, ic) = q_ad(i, k-1, ic) + dp(i, k-1)*qup_ad
              dp_ad(i, k-1) = dp_ad(i, k-1) + q(i, k-1, ic)*qup_ad
            END IF
          END IF
        END DO
      END DO
      CALL POPINTEGER4(k)
      DO i=im,1,-1
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) THEN
          CALL POPREAL8(q(i, 1, ic))
          q_ad(i, 1, ic) = 0.0_8
          CALL POPREAL8(q(i, 2, ic))
          temp = dp(i, 1)/dp(i, 2)
          temp_ad = q(i, 1, ic)*q_ad(i, 2, ic)/dp(i, 2)
          q_ad(i, 1, ic) = q_ad(i, 1, ic) + temp*q_ad(i, 2, ic)
          dp_ad(i, 1) = dp_ad(i, 1) + temp_ad
          dp_ad(i, 2) = dp_ad(i, 2) - temp*temp_ad
        END IF
      END DO
    END DO
  END SUBROUTINE FILLZ_ADM

END MODULE FV_FILL_ADM_MOD
