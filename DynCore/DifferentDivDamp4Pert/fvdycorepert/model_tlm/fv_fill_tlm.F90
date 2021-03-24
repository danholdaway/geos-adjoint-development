module fv_fill_tlm_mod

 use fv_arrays_mod, only: p_precision

 implicit none
 public fillz_tlm

 contains

  SUBROUTINE FILLZ_TLM(im, km, nq, q, q_tl, dp, dp_tl)
    IMPLICIT NONE
! No. of longitudes
    INTEGER, INTENT(IN) :: im
! No. of levels
    INTEGER, INTENT(IN) :: km
! Total number of tracers
    INTEGER, INTENT(IN) :: nq
! pressure thickness
    REAL(p_precision), INTENT(IN) :: dp(im, km)
    REAL(p_precision), INTENT(IN) :: dp_tl(im, km)
! tracer mixing ratio
    REAL(p_precision), INTENT(INOUT) :: q(im, km, nq)
    REAL(p_precision), INTENT(INOUT) :: q_tl(im, km, nq)
! !LOCAL VARIABLES:
    INTEGER :: i, k, ic
    REAL :: qup, qly, dup
    REAL :: qup_tl, qly_tl, dup_tl
    INTRINSIC MAX
    INTRINSIC MIN
    DO ic=1,nq
! Top layer
      DO i=1,im
        IF (q(i, 1, ic) .LT. 0.) THEN
          q_tl(i, 2, ic) = q_tl(i, 2, ic) + ((q_tl(i, 1, ic)*dp(i, 1)+q(&
&           i, 1, ic)*dp_tl(i, 1))*dp(i, 2)-q(i, 1, ic)*dp(i, 1)*dp_tl(i&
&           , 2))/dp(i, 2)**2
          q(i, 2, ic) = q(i, 2, ic) + q(i, 1, ic)*dp(i, 1)/dp(i, 2)
          q_tl(i, 1, ic) = 0.0_8
          q(i, 1, ic) = 0.
        END IF
      END DO
! Interior
      DO k=2,km-1
        DO i=1,im
          IF (q(i, k, ic) .LT. 0.) THEN
            IF (0. .LT. q(i, k-1, ic)*dp(i, k-1)) THEN
              qup_tl = q_tl(i, k-1, ic)*dp(i, k-1) + q(i, k-1, ic)*dp_tl&
&               (i, k-1)
              qup = q(i, k-1, ic)*dp(i, k-1)
            ELSE
              qup = 0.
              qup_tl = 0.0
            END IF
            qly_tl = -(q_tl(i, k, ic)*dp(i, k)+q(i, k, ic)*dp_tl(i, k))
            qly = -(q(i, k, ic)*dp(i, k))
            IF (0.5*qly .GT. 0.99*qup) THEN
              dup_tl = 0.99*qup_tl
              dup = 0.99*qup
            ELSE
              dup_tl = 0.5*qly_tl
              dup = 0.5*qly
            END IF
            q_tl(i, k-1, ic) = q_tl(i, k-1, ic) - (dup_tl*dp(i, k-1)-dup&
&             *dp_tl(i, k-1))/dp(i, k-1)**2
            q(i, k-1, ic) = q(i, k-1, ic) - dup/dp(i, k-1)
! Borrow from below: q(i,k,ic) is still negative at this stage
            q_tl(i, k+1, ic) = q_tl(i, k+1, ic) - ((qly_tl-dup_tl)*dp(i&
&             , k+1)-(qly-dup)*dp_tl(i, k+1))/dp(i, k+1)**2
            q(i, k+1, ic) = q(i, k+1, ic) - (qly-dup)/dp(i, k+1)
            q_tl(i, k, ic) = 0.0_8
            q(i, k, ic) = 0.
          END IF
        END DO
      END DO
! Bottom layer
      k = km
      DO i=1,im
        IF (q(i, k, ic) .LT. 0. .AND. q(i, k-1, ic) .GT. 0.) THEN
! Borrow from above
          qup_tl = q_tl(i, k-1, ic)*dp(i, k-1) + q(i, k-1, ic)*dp_tl(i, &
&           k-1)
          qup = q(i, k-1, ic)*dp(i, k-1)
          qly_tl = -(q_tl(i, k, ic)*dp(i, k)+q(i, k, ic)*dp_tl(i, k))
          qly = -(q(i, k, ic)*dp(i, k))
          IF (qly .GT. qup) THEN
            dup_tl = qup_tl
            dup = qup
          ELSE
            dup_tl = qly_tl
            dup = qly
          END IF
          q_tl(i, k-1, ic) = q_tl(i, k-1, ic) - (dup_tl*dp(i, k-1)-dup*&
&           dp_tl(i, k-1))/dp(i, k-1)**2
          q(i, k-1, ic) = q(i, k-1, ic) - dup/dp(i, k-1)
          q_tl(i, k, ic) = q_tl(i, k, ic) + (dup_tl*dp(i, k)-dup*dp_tl(i&
&           , k))/dp(i, k)**2
          q(i, k, ic) = q(i, k, ic) + dup/dp(i, k)
        END IF
      END DO
    END DO
  END SUBROUTINE FILLZ_TLM

end module fv_fill_tlm_mod
