#!/bin/csh -f

reset

source $CURMODEL/src/g5_modules

#Select compiler
set COMP = "ifort"

#Select complier options
set OPTIONS  = " -vec-report0 -O3"

# Set up the netCDF libraries
set NETCDF = "/usr/local/other/netcdf/4.1.1/serial_Intel-11.1.069"
set HDF    = "/usr/local/other/hdf5/1.8.4/serial_intel-11.1.069"
set ZLIB   = "/usr/local/other/zlib/1.2.3-intel-11.1.069"
set SZIP   = "/usr/local/other/szip/2.1/Intel-11.1.069"

set INPUTDIR = "/home/drholdaw/Lin_Moist_Physics"

set INCS = "-I${NETCDF}/include -I${INPUTDIR}/Inputs"
set LIBS = "-L${NETCDF}/lib -lnetcdf -L${HDF}/lib -lhdf5_hl -lhdf5 -L${ZLIB}/lib -lz -lm -L${SZIP}/lib -lsz -lmkl_core -lmkl_intel_lp64 -lmkl_sequential"

set LDFLAGS = "-lmkl_core -lmkl_intel_lp64 -lmkl_sequential ${INCS} ${LIBS}"


#Main program
set PROGRAM = "compEOFs.f90"

#Optional executable name
set EXEC = "-o compEOFs.x"

$COMP $OPTIONS $EXEC $PROGRAM $LDFLAGS 

echo Running Program

compEOFs.x


rm compEOFs.x


