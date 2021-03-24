#!/bin/csh

reset

source /discover/home/drholdaw/.cshrc

#Tapenade options
set opts = "-tgtvarname _tl -tgtfuncname _tlm -adjvarname _ad -adjfuncname _adm"


#Generate tangent linear code
#----------------------------
@ start_time = `date +%s`

rm -rf tlm/
mkdir tlm/

$TAPENADE_HOME/bin/tapenade ${opts} -d -O tlm/ -head "ConvectionMod.convection (tc bcnv)/(tc bcnv)" fpp_ConvectionMod.F90 

cd tlm/

#Rename file ending to match the model
rename mod_d.f90 _tlm.F90 *mod_d.f90

cd ../

@ end_time = `date +%s`
@ diff = $end_time - $start_time
echo "TLM generation took $diff seconds"

#Generate adjoint code
#---------------------
@ start_time = `date +%s`

#Generate the adjoint with all the adm/fwd/bwd routines we need
#--------------------------------------------------------------

rm -rf adm/
mkdir adm/

$TAPENADE_HOME/bin/tapenade ${opts} -b -O adm/ -head "ConvectionMod.convection (tc bcnv)/(tc bcnv)" fpp_ConvectionMod.F90 

cd adm/
rename mod_b.f90 _adm.F90 *mod_b.f90

cd ../

@ end_time = `date +%s`
@ diff = $end_time - $start_time
echo "Adjoint generation took $diff seconds"

#Done
