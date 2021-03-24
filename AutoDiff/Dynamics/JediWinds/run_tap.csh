#!/bin/csh

reset

source /discover/home/drholdaw/.cshrc

#Tapenade options
set opts = "-tgtvarname _tl -tgtfuncname _tlm -adjvarname _ad -adjfuncname _adm"


#Generate tangent linear code
#----------------------------
# @ start_time = `date +%s`
# 
# rm -rf tlm/
# mkdir tlm/
# rm -rf tlm/tapenadehtml
# 
# $TAPENADE_HOME/bin/tapenade ${opts} -d -O tlm/ -head "a2d_mod.a2d (ua va ud vd)/(ua va ud vd)" a2d_mod.F90
# #$TAPENADE_HOME/bin/tapenade ${opts} -d -O tlm/ -head "d2a_mod.d2a (ua va u v)/(ua va u v)" d2a_mod.F90
# 
# cd tlm/
# 
# #Rename file ending to match the model
# rename mod_diff.f90 tlm.F90 *mod_diff.f90
# 
# cd ../
# 
# @ end_time = `date +%s`
# @ diff = $end_time - $start_time
# echo "TLM generation took $diff seconds"



#Generate adjoint
#----------------
@ start_time = `date +%s`

rm -rf adm/
mkdir adm/
rm -rf adm/tapenadehtml

$TAPENADE_HOME/bin/tapenade ${opts} -b -O adm/ -head "psichi_to_uava_mod.psichi_to_uava (psi chi ua va)/(psi chi ua va)" psichi_to_uava_mod.F90
#$TAPENADE_HOME/bin/tapenade ${opts} -b -O adm/ -head "a2d_mod.a2d (ua va ud vd)/(ua va ud vd)" a2d_mod.F90
#$TAPENADE_HOME/bin/tapenade ${opts} -b -O adm/ -head "d2a_mod.d2a (ua va u v)/(ua va u v)" d2a_mod.F90

cd adm/

#Rename file ending to match the model
rename mod_diff.f90 adm.F90 *mod_diff.f90

cd ../

@ end_time = `date +%s`
@ diff = $end_time - $start_time
echo "adm generation took $diff seconds"

