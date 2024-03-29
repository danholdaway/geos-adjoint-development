module IRRAD_TL

IMPLICIT NONE

PRIVATE
PUBLIC :: irrad_d

contains

!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.9 (r5096) - 24 Feb 2014 16:53
!
!  Differentiation of irrad in forward (tangent) mode:
!   variations   of useful results: flxu_dev dfdts_dev flxd_dev
!   with respect to varying inputs: ta_dev
!   RW status of diff variables: flxu_dev:out ta_dev:in dfdts_dev:out
!                flxd_dev:out
SUBROUTINE IRRAD_D(m, np, ta_dev, ta_devd, flxu_dev, flxu_devd, flxd_dev&
& , flxd_devd, dfdts_dev, dfdts_devd)
  IMPLICIT NONE
!----- INPUTS -----
  INTEGER, INTENT(IN) :: m, np
  REAL, DIMENSION(m, np), INTENT(IN) :: ta_dev
  REAL, DIMENSION(m, np), INTENT(IN) :: ta_devd
  REAL, DIMENSION(m, np + 1), INTENT(OUT) :: flxu_dev
  REAL, DIMENSION(m, np+1), INTENT(OUT) :: flxu_devd
  REAL, DIMENSION(m, np + 1), INTENT(OUT) :: flxd_dev
  REAL, DIMENSION(m, np+1), INTENT(OUT) :: flxd_devd
  REAL, DIMENSION(m, np + 1), INTENT(OUT) :: dfdts_dev
  REAL, DIMENSION(m, np+1), INTENT(OUT) :: dfdts_devd
  INTEGER :: i, j
  flxu_dev = 0.0
  flxd_dev = 0.0
  dfdts_dev = 0.0
  flxu_devd = 0.0
  dfdts_devd = 0.0
  flxd_devd = 0.0
  DO i=1,m
    DO j=1,np
      flxu_devd(i, j) = 1e-6*ta_devd(i, j) + 1e-8*2*ta_dev(i, j)*ta_devd&
&       (i, j)
      flxu_dev(i, j) = ta_dev(i, j)*1e-6 + 1e-8*ta_dev(i, j)**2
      flxd_devd(i, j) = 1e-4*ta_devd(i, j) + 1e-10*3*ta_dev(i, j)**2*&
&       ta_devd(i, j)
      flxd_dev(i, j) = ta_dev(i, j)*1e-4 + 1e-10*ta_dev(i, j)**3
      dfdts_devd(i, j) = 0.5*(flxu_devd(i, j)+flxd_devd(i, j))
      dfdts_dev(i, j) = 0.5*(flxu_dev(i, j)+flxd_dev(i, j))
    END DO
  END DO
END SUBROUTINE IRRAD_D


end module IRRAD_TL
