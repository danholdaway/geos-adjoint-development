!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.12 (r6213) - 13 Oct 2016 10:30
!
MODULE FV_FILL_MOD_D
  USE FV_ARRAYS_MOD, ONLY : real8
  IMPLICIT NONE
  PUBLIC fillz
  PUBLIC fillz_tlm

CONTAINS
!  Differentiation of fillz in forward (tangent) mode:
!   variations   of useful results: q
!   with respect to varying inputs: q
  SUBROUTINE FILLZ_TLM(im, km, nq, q, q_tl, dp)
    IMPLICIT NONE
! No. of longitudes
    INTEGER, INTENT(IN) :: im
! No. of levels
    INTEGER, INTENT(IN) :: km
! Total number of tracers
    INTEGER, INTENT(IN) :: nq
! pressure thickness
    REAL(real8), INTENT(IN) :: dp(im, km)
! tracer mixing ratio
    REAL(real8), INTENT(INOUT) :: q(im, km, nq)
    REAL(real8), INTENT(INOUT) :: q_tl(im, km, nq)
! !LOCAL VARIABLES:
    LOGICAL :: zfix(im)
    REAL(real8) :: dm(km)
    REAL(real8) :: dm_tl(km)
    INTEGER :: i, k, ic
    REAL(real8) :: qup, qly, dup, dq, sum0, sum1, fac
    REAL(real8) :: qup_tl, qly_tl, dup_tl, dq_tl, sum0_tl, sum1_tl, &
&   fac_tl
    INTRINSIC MIN
    INTRINSIC MAX
    REAL*8 :: max1
    REAL*8 :: max1_tl
    dm_tl = 0.0_8
    DO ic=1,nq
! Top layer
      DO i=1,im
        IF (q(i, 1, ic) .LT. 0.) THEN
          q_tl(i, 2, ic) = q_tl(i, 2, ic) + dp(i, 1)*q_tl(i, 1, ic)/dp(i&
&           , 2)
          q(i, 2, ic) = q(i, 2, ic) + q(i, 1, ic)*dp(i, 1)/dp(i, 2)
          q_tl(i, 1, ic) = 0.0_8
          q(i, 1, ic) = 0.
        END IF
      END DO
! Interior
      zfix(:) = .false.
      DO k=2,km-1
        DO i=1,im
          IF (q(i, k, ic) .LT. 0.) THEN
            zfix(i) = .true.
            IF (q(i, k-1, ic) .GT. 0.) THEN
              IF (q(i, k-1, ic)*dp(i, k-1) .GT. -(q(i, k, ic)*dp(i, k))&
&             ) THEN
                dq_tl = -(dp(i, k)*q_tl(i, k, ic))
                dq = -(q(i, k, ic)*dp(i, k))
              ELSE
                dq_tl = dp(i, k-1)*q_tl(i, k-1, ic)
                dq = q(i, k-1, ic)*dp(i, k-1)
              END IF
              q_tl(i, k-1, ic) = q_tl(i, k-1, ic) - dq_tl/dp(i, k-1)
              q(i, k-1, ic) = q(i, k-1, ic) - dq/dp(i, k-1)
              q_tl(i, k, ic) = q_tl(i, k, ic) + dq_tl/dp(i, k)
              q(i, k, ic) = q(i, k, ic) + dq/dp(i, k)
            END IF
            IF (q(i, k, ic) .LT. 0.0 .AND. q(i, k+1, ic) .GT. 0.) THEN
              IF (q(i, k+1, ic)*dp(i, k+1) .GT. -(q(i, k, ic)*dp(i, k))&
&             ) THEN
                dq_tl = -(dp(i, k)*q_tl(i, k, ic))
                dq = -(q(i, k, ic)*dp(i, k))
              ELSE
                dq_tl = dp(i, k+1)*q_tl(i, k+1, ic)
                dq = q(i, k+1, ic)*dp(i, k+1)
              END IF
              q_tl(i, k+1, ic) = q_tl(i, k+1, ic) - dq_tl/dp(i, k+1)
              q(i, k+1, ic) = q(i, k+1, ic) - dq/dp(i, k+1)
              q_tl(i, k, ic) = q_tl(i, k, ic) + dq_tl/dp(i, k)
              q(i, k, ic) = q(i, k, ic) + dq/dp(i, k)
            END IF
          END IF
        END DO
      END DO
! Bottom layer
      k = km
      DO i=1,im
        IF (q(i, k, ic) .LT. 0. .AND. q(i, k-1, ic) .GT. 0.) THEN
          zfix(i) = .true.
! Borrow from above
          qup_tl = dp(i, k-1)*q_tl(i, k-1, ic)
          qup = q(i, k-1, ic)*dp(i, k-1)
          qly_tl = -(dp(i, k)*q_tl(i, k, ic))
          qly = -(q(i, k, ic)*dp(i, k))
          IF (qly .GT. qup) THEN
            dup_tl = qup_tl
            dup = qup
          ELSE
            dup_tl = qly_tl
            dup = qly
          END IF
          q_tl(i, k-1, ic) = q_tl(i, k-1, ic) - dup_tl/dp(i, k-1)
          q(i, k-1, ic) = q(i, k-1, ic) - dup/dp(i, k-1)
          q_tl(i, k, ic) = q_tl(i, k, ic) + dup_tl/dp(i, k)
          q(i, k, ic) = q(i, k, ic) + dup/dp(i, k)
        END IF
      END DO
! Perform final check and non-local fix if needed
      DO i=1,im
        IF (zfix(i)) THEN
          sum0 = 0.
          sum0_tl = 0.0_8
          DO k=2,km
            dm_tl(k) = dp(i, k)*q_tl(i, k, ic)
            dm(k) = q(i, k, ic)*dp(i, k)
            sum0_tl = sum0_tl + dm_tl(k)
            sum0 = sum0 + dm(k)
          END DO
          IF (sum0 .GT. 0.) THEN
            sum1 = 0.
            sum1_tl = 0.0_8
            DO k=2,km
              IF (0. .LT. dm(k)) THEN
                max1_tl = dm_tl(k)
                max1 = dm(k)
              ELSE
                max1 = 0.
                max1_tl = 0.0_8
              END IF
              sum1_tl = sum1_tl + max1_tl
              sum1 = sum1 + max1
            END DO
            fac_tl = (sum0_tl*sum1-sum0*sum1_tl)/sum1**2
            fac = sum0/sum1
            DO k=2,km
              IF (0. .LT. fac*dm(k)/dp(i, k)) THEN
                q_tl(i, k, ic) = (fac_tl*dm(k)+fac*dm_tl(k))/dp(i, k)
                q(i, k, ic) = fac*dm(k)/dp(i, k)
              ELSE
                q_tl(i, k, ic) = 0.0_8
                q(i, k, ic) = 0.
              END IF
            END DO
          END IF
        END IF
      END DO
    END DO
  END SUBROUTINE FILLZ_TLM
  SUBROUTINE FILLZ(im, km, nq, q, dp)
    IMPLICIT NONE
! No. of longitudes
    INTEGER, INTENT(IN) :: im
! No. of levels
    INTEGER, INTENT(IN) :: km
! Total number of tracers
    INTEGER, INTENT(IN) :: nq
! pressure thickness
    REAL(real8), INTENT(IN) :: dp(im, km)
! tracer mixing ratio
    REAL(real8), INTENT(INOUT) :: q(im, km, nq)
! !LOCAL VARIABLES:
    LOGICAL :: zfix(im)
    REAL(real8) :: dm(km)
    INTEGER :: i, k, ic
    REAL(real8) :: qup, qly, dup, dq, sum0, sum1, fac
    INTRINSIC MIN
    INTRINSIC MAX
    REAL*8 :: max1
    DO ic=1,nq
! Top layer
      DO i=1,im
        IF (q(i, 1, ic) .LT. 0.) THEN
          q(i, 2, ic) = q(i, 2, ic) + q(i, 1, ic)*dp(i, 1)/dp(i, 2)
          q(i, 1, ic) = 0.
        END IF
      END DO
! Interior
      zfix(:) = .false.
      DO k=2,km-1
        DO i=1,im
          IF (q(i, k, ic) .LT. 0.) THEN
            zfix(i) = .true.
            IF (q(i, k-1, ic) .GT. 0.) THEN
              IF (q(i, k-1, ic)*dp(i, k-1) .GT. -(q(i, k, ic)*dp(i, k))&
&             ) THEN
                dq = -(q(i, k, ic)*dp(i, k))
              ELSE
                dq = q(i, k-1, ic)*dp(i, k-1)
              END IF
              q(i, k-1, ic) = q(i, k-1, ic) - dq/dp(i, k-1)
              q(i, k, ic) = q(i, k, ic) + dq/dp(i, k)
            END IF
            IF (q(i, k, ic) .LT. 0.0 .AND. q(i, k+1, ic) .GT. 0.) THEN
              IF (q(i, k+1, ic)*dp(i, k+1) .GT. -(q(i, k, ic)*dp(i, k))&
&             ) THEN
                dq = -(q(i, k, ic)*dp(i, k))
              ELSE
                dq = q(i, k+1, ic)*dp(i, k+1)
              END IF
              q(i, k+1, ic) = q(i, k+1, ic) - dq/dp(i, k+1)
              q(i, k, ic) = q(i, k, ic) + dq/dp(i, k)
            END IF
          END IF
        END DO
      END DO
! Bottom layer
      k = km
      DO i=1,im
        IF (q(i, k, ic) .LT. 0. .AND. q(i, k-1, ic) .GT. 0.) THEN
          zfix(i) = .true.
! Borrow from above
          qup = q(i, k-1, ic)*dp(i, k-1)
          qly = -(q(i, k, ic)*dp(i, k))
          IF (qly .GT. qup) THEN
            dup = qup
          ELSE
            dup = qly
          END IF
          q(i, k-1, ic) = q(i, k-1, ic) - dup/dp(i, k-1)
          q(i, k, ic) = q(i, k, ic) + dup/dp(i, k)
        END IF
      END DO
! Perform final check and non-local fix if needed
      DO i=1,im
        IF (zfix(i)) THEN
          sum0 = 0.
          DO k=2,km
            dm(k) = q(i, k, ic)*dp(i, k)
            sum0 = sum0 + dm(k)
          END DO
          IF (sum0 .GT. 0.) THEN
            sum1 = 0.
            DO k=2,km
              IF (0. .LT. dm(k)) THEN
                max1 = dm(k)
              ELSE
                max1 = 0.
              END IF
              sum1 = sum1 + max1
            END DO
            fac = sum0/sum1
            DO k=2,km
              IF (0. .LT. fac*dm(k)/dp(i, k)) THEN
                q(i, k, ic) = fac*dm(k)/dp(i, k)
              ELSE
                q(i, k, ic) = 0.
              END IF
            END DO
          END IF
        END IF
      END DO
    END DO
  END SUBROUTINE FILLZ
END MODULE FV_FILL_MOD_D
