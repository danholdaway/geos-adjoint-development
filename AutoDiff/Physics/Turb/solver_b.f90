SUBROUTINE SOLVER_B(im, jm, lm, a, b, c, yb)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: im, jm, lm
  REAL, DIMENSION(im, jm, lm), INTENT(IN) :: a, b, c
  REAL, DIMENSION(im, jm, lm), INTENT(INOUT) :: yb
  INTEGER :: i, j, l
  REAL, DIMENSION(im, jm, lm) :: a1, b1, c1
  REAL :: tempb(im, jm)
!Make temporary copies so as not to overwrite
  a1 = a
  b1 = b
  c1 = c
!LU Decomposition
!----------------
  DO i=1,im
    DO j=1,jm
      b1(i, j, 1) = 1./b1(i, j, 1)
      DO l=2,lm
        a1(i, j, l) = a1(i, j, l)*b1(i, j, l-1)
        b1(i, j, l) = 1./(b1(i, j, l)-c1(i, j, l-1)*a1(i, j, l))
      END DO
    END DO
  END DO
  DO l=1,lm-1,1
    tempb = b1(:, :, l)*yb(:, :, l)
    yb(:, :, l+1) = yb(:, :, l+1) - c1(:, :, l)*tempb
    yb(:, :, l) = tempb
  END DO
  yb(:, :, lm) = b1(:, :, lm-1)*yb(:, :, lm)/(b1(:, :, lm-1)-a1(:, :, lm&
&   )*(c1(:, :, lm-1)*b1(:, :, lm-1)+1.0))
  DO l=lm,2,-1
    yb(:, :, l-1) = yb(:, :, l-1) - a1(:, :, l)*yb(:, :, l)
  END DO
END SUBROUTINE SOLVER_B

