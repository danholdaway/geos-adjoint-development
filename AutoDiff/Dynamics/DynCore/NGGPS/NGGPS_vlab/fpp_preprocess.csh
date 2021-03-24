#!/bin/csh
 
reset

module purge
module load comp/intel-18.0.0.128
module load mpi/impi-18.0.0.128

which fpp

rm -rf PreProcessed
mkdir PreProcessed

set mydir = `pwd`

#If using ifort force with -D__INTEL_COMPILER as not automatic with fpp

set FFLAGS = "-DMAPL_MODE -DSPMD -DTIMING -D__INTEL_COMPILER -DOVERLOAD_R4"

cd ProgressTracking/Part3_Compiles/

foreach modelfile (`ls *0`)
  echo $modelfile
  fpp -P $FFLAGS -I${mydir}/fms $modelfile > ${mydir}/ProgressTracking/Part4_PreProcessed/$modelfile
end

cd $mydir
