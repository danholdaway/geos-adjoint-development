#!/bin/csh
 
reset

module purge
module load comp/intel-18.0.0.128
module load mpi/impi-18.0.0.128

which fpp

set mydir = `pwd`

set FFLAGS = "-DsysLinux -DESMA64 -DHAS_NETCDF4 -DHAS_NETCDF3 -DH5_HAVE_PARALLEL -DNETCDF_NEED_NF_MPIIO -DHAVE_SHMEM"

foreach modelfile (`ls *0 | grep -v fpp_`)
  echo $modelfile
  fpp -P $FFLAGS $modelfile > fpp_$modelfile
end
