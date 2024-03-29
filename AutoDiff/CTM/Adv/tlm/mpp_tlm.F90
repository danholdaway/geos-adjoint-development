!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.12 (r6213) - 13 Oct 2016 10:30
!
MODULE MPP_MOD_D
  USE FV_ARRAYS_MOD, ONLY : real8
  IMPLICIT NONE
  PRIVATE 
  PUBLIC mpp_sum
  PUBLIC mpp_sum_tlm
  INTERFACE MPP_SUM
      MODULE PROCEDURE MPP_SUM_R1
      MODULE PROCEDURE MPP_SUM_R2
  END INTERFACE

  INTERFACE MPP_SUM_TLM
      MODULE PROCEDURE MPP_SUM_R1_TLM
  END INTERFACE


CONTAINS
!  Differentiation of mpp_sum_r1 in forward (tangent) mode:
!   variations   of useful results: psums
!   with respect to varying inputs: psums
  SUBROUTINE MPP_SUM_R1_TLM(psums, psums_tl, dims)
    IMPLICIT NONE
    REAL(real8), INTENT(INOUT) :: psums(:)
    REAL(real8), INTENT(INOUT) :: psums_tl(:)
    INTEGER, INTENT(IN) :: dims
    INTRINSIC SUM
    psums_tl = SUM(psums_tl)
    psums = SUM(psums)
  END SUBROUTINE MPP_SUM_R1_TLM
  SUBROUTINE MPP_SUM_R1(psums, dims)
    IMPLICIT NONE
    REAL(real8), INTENT(INOUT) :: psums(:)
    INTEGER, INTENT(IN) :: dims
    INTRINSIC SUM
    psums = SUM(psums)
  END SUBROUTINE MPP_SUM_R1
  SUBROUTINE MPP_SUM_R2(psums, dims)
    IMPLICIT NONE
    REAL(real8), INTENT(INOUT) :: psums(:, :)
    INTEGER, INTENT(IN) :: dims
    INTRINSIC SUM
    psums = SUM(psums)
  END SUBROUTINE MPP_SUM_R2
END MODULE MPP_MOD_D
