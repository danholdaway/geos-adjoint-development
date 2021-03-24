#!/bin/csh

reset

source /discover/home/drholdaw/.cshrc

#Tapenade options
set opts = "-html -tgtvarname _tl -tgtfuncname _tlm -adjvarname _ad -adjfuncname _adm"


#Generate tangent linear code
#----------------------------
rm -f tlm/*.F90

$TAPENADE_HOME/bin/tapenade ${opts} -d -O tlm/ -head "fv_tracer2d_mod.offline_tracer_advection (q)/(q)" fv_tracer2d.F90 fv_fill.F90 fv_grid_utils.F90 fv_mp_mod.F90 mpp.F90 fv_arrays.F90 fv_grid_tools.F90 fv_mapz.F90 fv_timing.F90 mpp_domains.F90 tp_core.F90

cd tlm/
rename mod_d.f90 tlm.F90 *mod_d.f90
cd ../


#Generate adjoint code
#---------------------
rm -f adm/*.F90

$TAPENADE_HOME/bin/tapenade ${opts} -b -O adm/ -head "fv_tracer2d_mod.offline_tracer_advection (q)/(q)" fv_tracer2d.F90 fv_fill.F90 fv_grid_utils.F90 fv_mp_mod.F90 mpp.F90 fv_arrays.F90 fv_grid_tools.F90 fv_mapz.F90 fv_timing.F90 mpp_domains.F90 tp_core.F90

cd adm/
rename mod_b.f90 adm.F90 *mod_b.f90
cd ../



