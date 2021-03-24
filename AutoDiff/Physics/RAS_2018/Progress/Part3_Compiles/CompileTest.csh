#!/bin/csh

reset

module purge
module load comp/intel-18.0.0.128
module load mpi/impi-18.0.0.128

#Clean up
rm -f *.mod *.o

setenv FORT "ifort"

$FORT -c -I./ RASPARAMS.F90
$FORT -c -I./ MAPL_Constants.F90
$FORT -c -I./ aer_cloud.F90
$FORT -c -I./ GEOS_Utilities.F90
$FORT -c -I./ ras.F90

if (-e ras.mod) then
  echo "COMPILED CORRECTLY"
  rm -f *.mod *.o
else
  echo "COMPILATION ERROR"
endif
