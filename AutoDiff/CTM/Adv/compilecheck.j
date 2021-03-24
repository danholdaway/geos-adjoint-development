#!/bin/csh
 
reset

source /gpfsm/dnb31/drholdaw/Heracles-5_4_p3_CTM-r1/Linux/bin/g5_modules

#Clean up
rm -f *.mod *.o

ifort -check -c -r8 fv_arrays.F90 fv_timing.F90

ifort -check -c -r8 fv_fill.F90 mpp_domains.F90 mpp.F90 fv_mapz.F90

ifort -check -c -r8 fv_mp_mod.F90

ifort -check -c -r8 fv_grid_tools.F90 fv_grid_utils.F90

ifort -check -c -r8 tp_core.F90

ifort -check -c -r8 fv_tracer2d.F90

#Clean up
rm -f *.mod *.o
