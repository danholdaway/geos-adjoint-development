 subroutine irrad ( m, np, ta_dev, flxu_dev, flxd_dev, dfdts_dev)

 IMPLICIT NONE

 !----- INPUTS -----
 INTEGER, INTENT(IN)                     :: m, np

 REAL, DIMENSION(m,np),    INTENT(IN)    :: ta_dev

 REAL, DIMENSION(m,np+1), INTENT(OUT)    :: flxu_dev
 REAL, DIMENSION(m,np+1), INTENT(OUT)    :: flxd_dev
 REAL, DIMENSION(m,np+1), INTENT(OUT)    :: dfdts_dev

 INTEGER :: i,j

 flxu_dev = 0.0
 flxd_dev = 0.0
 dfdts_dev = 0.0

 DO i = 1,m
    DO j = 1,np


       flxu_dev(i,j) = ta_dev(i,j) * 1e-6 + 1e-8*ta_dev(i,j)**2
       flxd_dev(i,j) = ta_dev(i,j) * 1e-4 + 1e-10*ta_dev(i,j)**3
       dfdts_dev(i,j) = 0.5*(flxu_dev(i,j) + flxd_dev(i,j))


    endDO
 endDO

   end subroutine irrad


