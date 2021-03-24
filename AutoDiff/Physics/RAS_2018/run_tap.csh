#!/bin/csh

reset

source /discover/home/drholdaw/.cshrc

#Tapenade options
set opts = "-tgtvarname _tl -tgtfuncname _tlm -adjvarname _ad -adjfuncname _adm -r8"

#Generate tangent linear code
#----------------------------
@ start_time = `date +%s`

rm -rf tlm/
mkdir tlm/
rm -rf tlm/tapenadehtml

$TAPENADE_HOME/bin/tapenade ${opts} -d -O tlm/ -head "RAS.RASE (THO QHO UHO VHO CLW FLXD CNV_PRC3 CNV_UPDFRC)/(THO QHO)" ras.F90 aer_cloud.F90 GEOS_Utilities.F90 MAPL_Constants.F90 RASPARAMS.F90

cd tlm/

#Rename file ending to match the model
rename _diff.f90 _tlm.F90 *_diff.f90

cd ../

#Generate adjoint code
#---------------------

@ start_time = `date +%s`

rm -rf adm/
mkdir adm/
rm -rf adm/tapenadehtml

$TAPENADE_HOME/bin/tapenade ${opts} -b -O adm/ -head "RAS.RASE (THO QHO UHO VHO CLW FLXD CNV_PRC3 CNV_UPDFRC)/(THO QHO UHO VHO CLW FLXD CNV_PRC3 CNV_UPDFRC)" ras.F90 aer_cloud.F90 GEOS_Utilities.F90 MAPL_Constants.F90 RASPARAMS.F90

cd adm/

#Rename file ending to match the model
rename _diff.f90 _adm.F90 *_diff.f90

cd ../
