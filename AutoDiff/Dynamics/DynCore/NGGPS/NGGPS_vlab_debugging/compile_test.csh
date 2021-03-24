#!/bin/csh
 
reset

module purge
module load comp/intel-18.0.0.128
module load mpi/impi-18.0.0.128

#Clean up
rm -f *.mod *.o fms/*.mod fms/*.o tools/*.mod tools/*.o MAPL/*.mod MAPL/*.o

#Paths to include
set mydir = `pwd`

set FORT = "ifort"

set FFLAGS = "-DsysLinux -DESMA64 -DHAS_NETCDF4 -DHAS_NETCDF3 -DH5_HAVE_PARALLEL -DNETCDF_NEED_NF_MPIIO -DMAPL_MODE -DSKIP_NO_QUAD_PRECISION -DREAD_GRID -DNO_GRID_G -DSPMD -DTIMING -DOLDMPP -DHAVE_SHMEM"

set incs = "-I ${mydir} -I ${mydir}/fms/ -I ${mydir}/tools/ -I ${mydir}/MAPL/"

cd MAPL/
$FORT $FFLAGS -c ${incs} MAPL.F90
cd ../

cd fms/
$FORT $FFLAGS -c ${incs} constants.F90 platform.F90 field_manager.F90 horiz_interp_type.F90 mpp_domains.F90 mpp_io.F90 mpp_paramter.F90 time_manager.F90 fms_io.F90 tracer_manager.F90 mpp.F90 diag_manager.F90
cd ../

$FORT $FFLAGS -c ${incs} fv_arrays.F90

cd tools
$FORT $FFLAGS -c ${incs} fv_timing.F90 fv_mp_mod.F90 fv_restart.F90 fv_nudge.F90 fv_diagnostics.F90
cd ../

$FORT $FFLAGS -c ${incs} fv_grid_utils.F90

cd tools
$FORT $FFLAGS -c ${incs} init_hydro.F90
cd ../

$FORT $FFLAGS -c ${incs} fv_control.F90 fv_cmp.F90 fv_sg.F90 fv_grid_utils.F90 a2b_edge.F90 tp_core.F90 sw_core.F90 nh_utils.F90 nh_core.F90 fv_fill.F90 fv_mapz.F90 boundary.F90 fv_nesting.F90 fv_tracer2d.F90 dyn_core.F90 fv_dynamics.F90

if (-e fv_dynamics_mod.mod) then
  echo "COMPILED CORRECTLY"
else
  echo "COMPILATION ERROR"
endif

#Clean up again once check complete
rm -f *.mod *.o fms/*.mod fms/*.o tools/*.mod tools/*.o MAPL/*.mod MAPL/*.o

