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
rm -rf tlm/tapenadehtml

$TAPENADE_HOME/bin/tapenade ${opts} -d -O tlm/ -head "wind_vt_mod.psichi_to_udvd (u v psi chi)/(u v psi chi)" wind_variables.F90 mpp_domains.F90 fv3jedi_geom_mod.F90
cd tlm/

#Rename file ending to match the model
rename mod_diff.f90 tlm.F90 *mod_diff.f90

cd ../

@ end_time = `date +%s`
@ diff = $end_time - $start_time
echo "TLM generation took $diff seconds"

#Generate tangent linear code
#----------------------------
@ start_time = `date +%s`

rm -rf adm/
mkdir adm/
rm -rf adm/tapenadehtml

$TAPENADE_HOME/bin/tapenade ${opts} -b -O adm/ -head "wind_vt_mod.psichi_to_udvd (u v psi chi)/(u v psi chi)" wind_variables.F90 mpp_domains.F90 fv3jedi_geom_mod.F90
cd adm/

#Rename file ending to match the model
rename mod_diff.f90 adm.F90 *mod_diff.f90

cd ../

@ end_time = `date +%s`
@ diff = $end_time - $start_time
echo "ADM generation took $diff seconds"

