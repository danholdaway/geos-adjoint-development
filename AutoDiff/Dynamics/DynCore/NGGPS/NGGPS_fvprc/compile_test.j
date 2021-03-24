#!/bin/csh
 
reset

source /gpfsm/dnb31/drholdaw/GEOSadas-5_16_5_p1_DEV5/Linux/bin/g5_modules

#Clean up
rm -f *.mod *.o fms/*.mod fms/*.o tools/*.mod tools/*.o MAPL/*.mod MAPL/*.o

#Paths to include
set mydir = `pwd`

set incs = "-I ${mydir} -I ${mydir}/fms/ -I ${mydir}/tools/ -I ${mydir}/MAPL/"

cd MAPL/
ifort -check -c MAPL.F90 ${incs}
cd ../

cd fms/
ifort -check -c constants.F90 field_manager.F90 horiz_interp_type.F90 mpp_domains.F90 mpp_io.F90 mpp_paramter.F90 time_manager.F90 fms_io.F90 tracer_manager.F90 ${incs}
cd ../

cd tools
ifort -check -c fv_diagnostics.F90 fv_nudge.F90  fv_timing.F90 ${incs}
cd ../

ifort -check -c fv_arrays.F90 ${incs}

cd fms/
ifort -check -c mpp.F90 ${incs}
cd ../

cd tools
ifort -check -c fv_mp_mod.F90 fv_restart.F90 ${incs}
cd ../

ifort -check -c fv_grid_utils.F90 ${incs}

cd tools
ifort -check -c init_hydro.F90 ${incs}
cd ../

ifort -check -c fv_grid_utils.F90 a2b_edge.F90 tp_core.F90 sw_core.F90 nh_core.F90 dyn_core.F90 fv_mapz.F90 fv_tracer2d.F90 boundary.F90 fv_nesting.F90 fv_dynamics.F90 ${incs}


#Clean up again once check complete
rm -f *.mod *.o fms/*.mod fms/*.o tools/*.mod tools/*.o MAPL/*.mod MAPL/*.o

