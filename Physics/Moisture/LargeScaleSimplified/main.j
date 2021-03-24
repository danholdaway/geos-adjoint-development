#!/bin/csh -f

reset

#Select compiler
set COMP = "ifort"

#Select complier options
set OPTIONS  = " -traceback -r{8} -w -O2 -ftz -align all -fno-alias -fp-model precise"

# Set up the netCDF libraries
set NETCDF = "/usr/local/other/netcdf/4.1.1/serial_Intel-11.1.069"
set HDF    = "/usr/local/other/hdf5/1.8.4/serial_intel-11.1.069"
set ZLIB   = "/usr/local/other/zlib/1.2.3-intel-11.1.069"
set SZIP   = "/usr/local/other/szip/2.1/Intel-11.1.069"

set INCS = "-I${NETCDF}/include"
set LIBS = "-L${NETCDF}/lib -lnetcdf -L${HDF}/lib -lhdf5_hl -lhdf5 -L${ZLIB}/lib -lz -lm -L${SZIP}/lib -lsz -lmkl_core -lmkl_intel_lp64 -lmkl_sequential"

set LDFLAGS = "${INCS} ${LIBS}"

# Required sub-routines and modules
set SUBS = "global.F90"

#Main program
set PROGRAM = "driver.f90"
#set PROGRAM = "main1.2.2.f90"

#Optional executable name
set EXEC = "-o main.x"

$COMP $OPTIONS $EXEC $SUBS $PROGRAM $LDFLAGS 
#$COMP $OPTIONS $EXEC $PROGRAM $LDFLAGS

echo StartingComputations
./main.x
echo Done
